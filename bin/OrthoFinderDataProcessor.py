#!/usr/bin/env python3

import os
import re
import argparse
import logging
import csv
import sys
from collections import defaultdict
from pathlib import Path

# Assuming helper scripts are in the same directory or in the python path
from OrthoFinderStatsProcessor import OrthoFinderStatsProcessor
from OrthoFinderTreeParser import SpeciesTree
from OrthogroupsStratigraphyClassifier import (
    OrthogroupsStratigraphyClassifier,
    SPECIES_SPECIFIC_CLUSTERED,
    SPECIES_SPECIFIC_UNCLUSTERED,
)
from CLSGClassificationVerifier import CLSGClassificationVerifier
from FastaGeneFilter import GeneSequenceExtractor
from CSVFileMerger import CSVFileMerger
from OrthologyCrossChecker import OrthologyCrossChecker

class PipelineInputError(Exception):
    """Invalid or inconsistent input; the run cannot produce trustworthy output."""


class OrthoFinderDataProcessor:
    def __init__(self, out_dir, orthofinder_dir=None, prt_dir=None,
                stats_filename="Statistics_PerSpecies.tsv",
                hog_filename="N0.tsv",
                tree_filename="SpeciesTree_rooted_node_labels.txt",
                orthogroups_filename="Orthogroups.tsv",
                orthologues_dirname="Orthologues",
                skip_cross_check=False,
                skip_hog_validation=False,
                species_map_file=None,
                allow_missing=False):
        self.out_dir = Path(out_dir)
        # Create the output directory immediately, before setting up logging.
        self.out_dir.mkdir(parents=True, exist_ok=True)
        # SpeciesStats.csv, c_lsg_summary.csv, and global_gene_counts_per_node.tsv
        # are merged into global_summary.csv, which is the single canonical
        # master table at the top level. The three sources are kept, not
        # deleted, since they're useful for debugging one merge step at a
        # time -- but live in intermediate/ so they don't clutter the main
        # output alongside the merged table.
        self.intermediate_dir = self.out_dir / "intermediate"
        self.intermediate_dir.mkdir(parents=True, exist_ok=True)

        self.orthofinder_dir = Path(orthofinder_dir) if orthofinder_dir else None
        self.prt_dir = Path(prt_dir) if prt_dir else None
        self.species_map_file = species_map_file
        self.stats_filename = stats_filename
        self.hog_filename = hog_filename
        self.tree_filename = tree_filename
        self.orthogroups_filename = orthogroups_filename
        self.orthologues_dirname = orthologues_dirname
        # Pairwise Orthologues is a secondary, independent cross-check
        # (reciprocal-best-hit -- a different method entirely from the
        # clustering Orthogroups.tsv is built from) on top of the primary
        # classification. Mandatory by default; --skip-cross-check is the
        # explicit opt-out.
        self.skip_cross_check = skip_cross_check
        # HOG-table cross-validation is optional and non-blocking: it can
        # only ever be a secondary sanity check now (see
        # OrthogroupsStratigraphyClassifier's module docstring for why HOG
        # tables are no longer the primary classification source), so a
        # missing HOG file just means validation is skipped, not an error.
        self.skip_hog_validation = skip_hog_validation
        # Species whose protein FASTA is missing/unreadable are always reported
        # and make the run exit non-zero, unless allow_missing is set.
        self.allow_missing = allow_missing
        self.skipped_species = {}   # species basename -> reason
        self.write_failures = []    # (path, error)
        self.species_map = self._load_species_map()
        self.setup_logging()
        self.stats_processor = None
        self.stats_lookup = None
        self.tree_parser = None
        self.stratigraphy_classifier = None
        self.hog_verifier = None
        self.species_with_hog_data = set()
        self.cross_checker = None
        self.hog_disagreements = []
        self.stats_mismatches = []

    def _load_species_map(self):
        if not self.species_map_file: return {}
        s_map = {}
        try:
            with open(self.species_map_file, 'r', newline='', encoding='utf-8-sig') as f:
                sniffer = csv.Sniffer()
                dialect = sniffer.sniff(f.read(1024)); f.seek(0)
                reader = csv.reader(f, dialect)
                header = [h.lower().strip() for h in next(reader)]
                name_col_idx, base_col_idx = header.index('species'), header.index('basename')
                for row in reader:
                    if len(row) > max(name_col_idx, base_col_idx):
                        basename, species_name = row[base_col_idx].strip(), row[name_col_idx].strip()
                        if basename and species_name: s_map[basename] = species_name
            logging.info(f"Successfully loaded {len(s_map)} entries from species map.")
            return s_map
        except Exception as e:
            logging.error(f"Error reading species map file '{self.species_map_file}': {e}", exc_info=True)
            return {}

    def setup_logging(self):
        log_file = self.out_dir / "OrthoFinderDataProcessor.log"
        logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] - %(message)s", handlers=[logging.FileHandler(log_file, mode='w'), logging.StreamHandler()])
        logging.info("Logging initialized.")

    def _list_top_level(self, directory):
        """Contents of directory, for error messages -- lets a wrong or
        incomplete --orthofinder_dir be spotted at a glance instead of
        just "not found somewhere in this tree". Descends through any
        chain of single-entry wrapper directories (e.g. orthofinder_dir/
        containing only Results_Aug25/) so the message shows the actual
        results layout, not just an uninformative single subfolder name."""
        current = directory
        try:
            entries = sorted(os.listdir(current))
        except OSError as e:
            return f"(could not list {directory}: {e})"
        while len(entries) == 1 and os.path.isdir(os.path.join(current, entries[0])):
            current = os.path.join(current, entries[0])
            try:
                entries = sorted(os.listdir(current))
            except OSError as e:
                return f"(could not list {current}: {e})"
        if not entries:
            return f"{current} is empty"
        return f"{current} contains: {', '.join(entries)}"

    def _find_file(self, directory, filename):
        for root, _, files in os.walk(directory):
            if filename in files:
                found_path = os.path.join(root, filename)
                logging.info(f"Found {filename} at {found_path}")
                return found_path
        raise FileNotFoundError(
            f"{filename} not found anywhere under {directory}. "
            f"{self._list_top_level(directory)}"
        )

    def _find_dir(self, directory, dirname):
        for root, dirs, _ in os.walk(directory):
            if dirname in dirs:
                found_path = os.path.join(root, dirname)
                logging.info(f"Found {dirname}/ at {found_path}")
                return found_path
        raise FileNotFoundError(
            f"{dirname}/ not found anywhere under {directory}. "
            f"{self._list_top_level(directory)}"
        )

    def _find_hog_file(self):
        """
        Like _find_file, but for the (now optional, validation-only) HOG
        table: if the requested file (default N0.tsv) isn't there, checks
        whether any other N#.tsv hierarchical orthogroup file exists and
        says so -- this is the common outgroup situation the README
        documents (an outgroup pushes the true root out one level and
        N0.tsv is never written), so a bare "not found" would hide a fix
        that's usually one flag away, even though HOG validation is
        optional now and a missing file just means it's skipped.
        """
        try:
            return self._find_file(self.orthofinder_dir, self.hog_filename)
        except FileNotFoundError as e:
            available = set()
            for root, _, files in os.walk(self.orthofinder_dir):
                for f in files:
                    if re.fullmatch(r'N\d+\.tsv', f):
                        available.add(f)
            if available:
                ordered = sorted(available, key=lambda n: int(n[1:-4]))
                raise FileNotFoundError(
                    f"{self.hog_filename} not found under {self.orthofinder_dir}, but "
                    f"found: {', '.join(ordered)}. Pass the right one via --hog-filename "
                    f"if you want HOG-based cross-validation (see the README's \"Why "
                    f"--hog-filename might not be N0.tsv\")."
                ) from None
            raise

    def initialize_processors(self):
        logging.info("Initializing data processors...")
        if not self.orthofinder_dir:
            logging.warning("OrthoFinder directory not provided."); return
        try:
            self.stats_converted_path = self.intermediate_dir / "SpeciesStats.csv"
            try:
                stats_file = self._find_file(self.orthofinder_dir, self.stats_filename)
                self.stats_processor = OrthoFinderStatsProcessor(stats_file, self.stats_converted_path, self.species_map_file)
                self.stats_processor.process_stats()
            except FileNotFoundError:
                # Not required for classification itself -- it only feeds one
                # section of the final merged global_summary.csv (and, below,
                # an optional validation of our own computed totals). A
                # missing/renamed stats file shouldn't abort the whole run
                # the way a missing Orthogroups.tsv or species tree must --
                # those two are genuinely required and still raise below.
                logging.warning(
                    f"{self.stats_filename} not found under {self.orthofinder_dir} -- "
                    f"SpeciesStats.csv, the stats cross-validation, and global_summary.csv's "
                    f"stats columns will all be skipped. If this OrthoFinder version renamed "
                    f"or moved it, pass the actual filename via --stats-filename."
                )

            # Species tree: required. Needed for MRCA computation, which the
            # primary classification depends on entirely.
            tree_file = self._find_file(self.orthofinder_dir, self.tree_filename)
            self.tree_parser = SpeciesTree(file_path=tree_file, species_map_file=self.species_map_file)
            n_unlabelled = self.tree_parser.count_unlabelled_internal_nodes()
            if n_unlabelled:
                raise PipelineInputError(
                    f"{n_unlabelled} internal node(s) in {tree_file} have no label. Every internal "
                    f"node must be labelled (N0, N1, ...); use OrthoFinder's "
                    f"SpeciesTree_rooted_node_labels.txt."
                )

            # Orthogroups.tsv: required. This is the primary classification
            # source -- see OrthogroupsStratigraphyClassifier's module
            # docstring for why (HOG tables are structurally blind to any
            # species outside the node they cover; Orthogroups.tsv never is).
            orthogroups_file = self._find_file(self.orthofinder_dir, self.orthogroups_filename)
            self.stratigraphy_classifier = OrthogroupsStratigraphyClassifier(orthogroups_file, self.tree_parser)
            self.stratigraphy_classifier.parse()
            tree_leaves = self.tree_parser.get_all_species()
            not_in_tree = sorted(set(self.stratigraphy_classifier.species_in_og_table) - tree_leaves)
            if not_in_tree:
                raise PipelineInputError(
                    f"{len(not_in_tree)} species in {orthogroups_file} are not leaves of the species "
                    f"tree {tree_file}: {', '.join(not_in_tree)}. Tree leaves: "
                    f"{', '.join(sorted(tree_leaves))}."
                )
            logging.info(f"Parsed {len(self.stratigraphy_classifier.og_species_genes)} orthogroups from {orthogroups_file}")

            self._initialize_hog_verifier()
            self._initialize_cross_checker()
            logging.info("All processors initialized successfully.")
        except FileNotFoundError as e:
            logging.error(f"Failed to initialize processors: {e}"); raise

    def _initialize_hog_verifier(self):
        """
        Optional, secondary: loads a HOG table (default N0.tsv) purely for
        cross-validating the primary Orthogroups.tsv-derived classification
        against it, gene by gene, where a gene happens to be visible in the
        HOG table. Never required -- if it can't be found, HOG validation
        is skipped with a warning, not an error.

        Reuses CLSGClassificationVerifier's independent HOG-table parser
        (gene_to_hog / hog_species_genes) rather than writing a second one,
        since that class's whole job is exactly this: parsing a HOG table
        into a per-gene, per-orthogroup-membership structure, independent
        of any other class in this pipeline.
        """
        if self.skip_hog_validation:
            logging.info("HOG-based cross-validation explicitly skipped (--skip-hog-validation).")
            return
        try:
            hog_file = self._find_hog_file()
        except FileNotFoundError as e:
            logging.warning(f"HOG-based cross-validation will be skipped: {e}")
            return

        self.hog_verifier = CLSGClassificationVerifier(hog_file, self.out_dir)
        try:
            self.hog_verifier.parse_hog_table()
            # A species can be a HEADER column in the HOG table (OrthoFinder
            # writes one column per species regardless) while every cell
            # under it is empty -- exactly the outgroup-blindness case this
            # whole redesign is about. species_in_hog_table (header columns)
            # is therefore NOT the right thing to check when diagnosing a
            # disagreement as "missing species" -- what matters is whether
            # the species has ANY actual gene entries anywhere in the table.
            # Computed once here (one pass over every HOG), not per-gene.
            self.species_with_hog_data = {
                sp for species_genes in self.hog_verifier.hog_species_genes.values()
                for sp, genes in species_genes.items() if genes
            }
            logging.info(
                f"Loaded {len(self.hog_verifier.hog_species_genes)} HOGs from {hog_file} "
                f"for cross-validation ({len(self.species_with_hog_data)} of "
                f"{len(self.hog_verifier.species_in_hog_table)} header species actually "
                f"have any data in it)."
            )
        except Exception as e:
            logging.warning(f"Could not parse {hog_file} for HOG cross-validation: {e}")
            self.hog_verifier = None

    def _initialize_cross_checker(self):
        """
        Locates the Orthologues/ directory and builds the secondary
        pairwise-Orthologues cross-check on the species-specific candidate
        bucket. This used to also check Orthogroups.tsv as a primary check,
        back when the HOG table was the primary classification source and
        Orthogroups.tsv was the downstream safety net; now that
        Orthogroups.tsv IS the primary source (see
        OrthogroupsStratigraphyClassifier), checking it against itself
        would be circular, so only the pairwise check remains. It's kept
        because it's a genuinely independent method (reciprocal-best-hit,
        not clustering) and cheap to run.
        """
        if self.skip_cross_check:
            logging.warning(
                "Cross-check explicitly skipped (--skip-cross-check). Candidate species-specific "
                "genes will NOT be verified against pairwise Orthologues."
            )
            return

        try:
            orthologues_dir = self._find_dir(self.orthofinder_dir, self.orthologues_dirname)
        except FileNotFoundError:
            logging.warning(f"{self.orthologues_dirname}/ not found under {self.orthofinder_dir} -- "
                             f"pairwise cross-check will be skipped.")
            return

        self.cross_checker = OrthologyCrossChecker(orthologues_dir=orthologues_dir)

    def _get_all_gene_ids_from_fasta(self, species_basename):
        if not self.prt_dir:
            self.skipped_species[species_basename] = "no --prt_dir given"
            return set()
        fasta_file = Path(self.prt_dir) / f"{species_basename}.prt"
        if not fasta_file.is_file():
            for ext in ['.fa', '.fasta', '.faa', '.fna']:
                fasta_file = Path(self.prt_dir) / f"{species_basename}{ext}"
                if fasta_file.is_file(): break
            else:
                logging.warning(f"Protein FASTA file not found for {species_basename}, cannot classify this species.")
                self.skipped_species[species_basename] = f"protein FASTA not found in {self.prt_dir}"
                return set()
        all_gene_ids = set()
        try:
            with open(fasta_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        gene_id = line[1:].split()[0].strip()
                        all_gene_ids.add(gene_id)
            logging.info(f"Extracted {len(all_gene_ids)} total gene IDs from {fasta_file.name}.")
            return all_gene_ids
        except Exception as e:
            logging.error(f"Could not read FASTA file {fasta_file}: {e}")
            self.skipped_species[species_basename] = f"could not read {fasta_file}: {e}"
            return set()

    def _write_per_species_summary(self, species_out_dir, filename_csn, summary_data):
        """Writes the summary CSV for a single species."""
        summary_file = species_out_dir / f'{filename_csn}_summary_per_node.csv'
        header = ['Node', 'Number_of_Genes']
        try:
            with open(summary_file, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(header)
                writer.writerows(summary_data)
            logging.info(f"Per-species summary written to {summary_file.name}")
        except IOError as e:
            logging.error(f"Failed to write per-species summary for {filename_csn}: {e}")
            self.write_failures.append((summary_file, e))

    def _write_global_node_summary(self, global_summary_data):
        """Writes the wide-format summary table for all species across all nodes."""
        if not global_summary_data:
            logging.warning("No data available for global node summary. Skipping.")
            return

        all_nodes = sorted(list(set(node for species_data in global_summary_data.values() for node in species_data.keys())))
        header = ['Species'] + all_nodes

        summary_file = self.intermediate_dir / "global_gene_counts_per_node.tsv"
        try:
            with open(summary_file, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(header)
                for species_name, node_data in sorted(global_summary_data.items()):
                    row = [species_name] + [node_data.get(node, 'NA') for node in all_nodes]
                    writer.writerow(row)
            logging.info(f"Global node summary written to {summary_file.name}")
        except IOError as e:
            logging.error(f"Failed to write global node summary: {e}")
            self.write_failures.append((summary_file, e))

        return summary_file

    def _cross_validate_against_hog(self, species_basename, node_buckets):
        """
        For every gene that's visible in the (optional) HOG table, derives
        what node the HOG table alone would place it at -- using the same
        MRCA mechanism as the primary classifier, just applied to the HOG
        table's species membership instead of Orthogroups.tsv's -- and
        compares against the primary classification. Every disagreement is
        logged, tagged with a diagnosed reason where determinable: if the
        gene's Orthogroups.tsv-visible species set includes a species that
        the HOG table's own species columns don't include at all, that
        fully explains the disagreement (the HOG table structurally cannot
        see that species). Otherwise the disagreement is tagged
        "unexplained" -- most likely genuine divergence between similarity
        clustering and tree-reconciliation (see this file's earlier
        discussion / the module docstring), not necessarily an error, but
        not silently assumed to be explained either.
        """
        if not self.hog_verifier:
            return

        for node_label, gene_ids in node_buckets.items():
            for gene_id in gene_ids:
                hog_id = self.hog_verifier.gene_to_hog.get((species_basename, gene_id))
                if hog_id is None:
                    continue  # gene not visible in the HOG table at all -- nothing to compare

                hog_species_genes = self.hog_verifier.hog_species_genes.get(hog_id, {})
                hog_other_species = sorted(
                    sp for sp, genes in hog_species_genes.items() if sp != species_basename and genes
                )
                if not hog_other_species:
                    hog_level = SPECIES_SPECIFIC_CLUSTERED
                else:
                    try:
                        hog_level = self.tree_parser.get_mrca([species_basename] + hog_other_species)
                    except ValueError:
                        continue  # a name in the HOG table isn't a tree leaf -- can't compare, not a disagreement

                if hog_level == node_label:
                    continue

                og_id = self.stratigraphy_classifier.gene_to_og.get((species_basename, gene_id))
                og_species = set(self.stratigraphy_classifier.og_species_genes.get(og_id, {}).keys()) if og_id else set()
                missing_from_hog = sorted(og_species - self.species_with_hog_data)
                reason = "missing_species_from_hog" if missing_from_hog else "unexplained"

                self.hog_disagreements.append({
                    "species": species_basename,
                    "gene_id": gene_id,
                    "orthogroups_level": node_label,
                    "hog_level": hog_level,
                    "reason": reason,
                    "missing_species": ";".join(missing_from_hog),
                })

    def _load_stats_lookup(self):
        """
        Reads back OrthoFinderStatsProcessor's own output (SpeciesStats.csv)
        for the stats cross-validation below, keyed by the same
        (species-map-renamed) species name it was written with. Cached
        after the first call.
        """
        if self.stats_lookup is not None:
            return self.stats_lookup
        self.stats_lookup = {}
        if not self.stats_processor or not self.stats_converted_path.exists():
            return self.stats_lookup
        with open(self.stats_converted_path, "r", newline="") as f:
            for row in csv.DictReader(f):
                self.stats_lookup[row["Species"]] = row
        return self.stats_lookup

    def _validate_against_stats(self, correct_species_name, clustered_count, unclustered_count):
        """
        Cross-checks this pipeline's own computed species-specific counts
        against OrthoFinder's own self-reported per-species statistics
        (Statistics_PerSpecies.tsv, converted to SpeciesStats.csv by
        OrthoFinderStatsProcessor). Both are ultimately derived from
        Orthogroups.tsv, so they should match exactly -- a mismatch means
        this pipeline's own counting logic disagrees with OrthoFinder's,
        which is worth surfacing loudly rather than trusting silently.
        """
        stats_row = self._load_stats_lookup().get(correct_species_name)
        if not stats_row:
            return

        try:
            expected_clustered = int(stats_row["Species_Specific_Clustered"])
            expected_unclustered = int(stats_row["Species_Specific_Unclustered"])
        except (KeyError, ValueError):
            return

        if clustered_count != expected_clustered or unclustered_count != expected_unclustered:
            self.stats_mismatches.append({
                "species": correct_species_name,
                "our_clustered": clustered_count,
                "orthofinder_clustered": expected_clustered,
                "our_unclustered": unclustered_count,
                "orthofinder_unclustered": expected_unclustered,
            })
            logging.warning(
                f"Stats mismatch for {correct_species_name}: this pipeline computed "
                f"{clustered_count} species-specific-clustered / {unclustered_count} unclustered, "
                f"but OrthoFinder's own stats report {expected_clustered} / {expected_unclustered}."
            )

    def process_and_classify_genes_for_species(self, species_basename):
        correct_species_name = self.species_map.get(species_basename, species_basename)
        filename_csn = correct_species_name.replace(" ", "_")
        logging.info(f"--- Processing: {species_basename} -> {correct_species_name} ---")

        all_genes_from_fasta = self._get_all_gene_ids_from_fasta(species_basename)
        if not all_genes_from_fasta:
            logging.warning(f"No genes found for {species_basename} (missing/empty protein FASTA). Skipping.")
            self.skipped_species.setdefault(species_basename, "protein FASTA is empty")
            return [], 0, 0, 0, {}

        species_out_dir = self.out_dir / filename_csn
        gene_ids_dir, protein_dir = species_out_dir / 'gene_IDs_by_node', species_out_dir / 'proteins_by_node'
        gene_ids_dir.mkdir(parents=True, exist_ok=True)
        protein_dir.mkdir(parents=True, exist_ok=True)

        gene_classification_data = []
        node_buckets = defaultdict(set)  # node_label -> set of gene_ids (species-specific buckets included)
        candidate_status = {}
        candidate_lsg_genes = set()

        for gene_id in sorted(all_genes_from_fasta):
            status, node_label = self.stratigraphy_classifier.classify_gene(species_basename, gene_id)
            node_buckets[node_label].add(gene_id)

            og_id = self.stratigraphy_classifier.gene_to_og.get((species_basename, gene_id), 'NA')
            gene_classification_data.append([correct_species_name, gene_id, node_label, og_id])

            if status in (SPECIES_SPECIFIC_CLUSTERED, SPECIES_SPECIFIC_UNCLUSTERED):
                candidate_status[gene_id] = status
                candidate_lsg_genes.add(gene_id)

        self._cross_validate_against_hog(species_basename, node_buckets)

        per_species_summary_data = []
        global_summary_node_data = {}
        for node_name, gene_ids in sorted(node_buckets.items()):
            per_species_summary_data.append([node_name, len(gene_ids)])
            global_summary_node_data[node_name] = len(gene_ids)
            self._save_gene_ids(gene_ids_dir / f'{species_basename}_{node_name}_gene_ids.txt', gene_ids)
            # Per-node protein FASTA, re-enabled: every stratigraphic bucket
            # gets its own FASTA, not just the final candidate set.
            self._filter_fasta(species_basename, "prt", protein_dir, node_name, gene_ids)

        clustered_count = len(node_buckets.get(SPECIES_SPECIFIC_CLUSTERED, set()))
        unclustered_count = len(node_buckets.get(SPECIES_SPECIFIC_UNCLUSTERED, set()))
        self._validate_against_stats(correct_species_name, clustered_count, unclustered_count)

        # Secondary cross-check: a gene that looks species-specific from
        # Orthogroups.tsv alone may still have a real ortholog findable via
        # reciprocal-best-hit that clustering missed. Verify every candidate
        # before it's treated as final.
        if self.cross_checker:
            verified_candidate_genes = self.cross_checker.filter_candidates(
                species_basename, candidate_lsg_genes, original_status_lookup=candidate_status
            )
            n_excluded = len(candidate_lsg_genes) - len(verified_candidate_genes)
            if n_excluded:
                logging.info(f"Cross-check excluded {n_excluded} of {len(candidate_lsg_genes)} "
                             f"candidate genes for {species_basename} (had a pairwise ortholog "
                             f"clustering missed).")
        else:
            verified_candidate_genes = candidate_lsg_genes

        self._save_gene_ids(
            gene_ids_dir / f'{species_basename}_candidate_lsg_gene_ids.txt',
            verified_candidate_genes
        )

        if self.prt_dir and verified_candidate_genes:
            logging.info(f"Writing combined candidate FASTA for {species_basename}")
            self._filter_fasta(
                species_basename, "prt", protein_dir, "candidate_LSGs", verified_candidate_genes
            )

        total_genes_accounted_for = len(all_genes_from_fasta)
        c_lsg_percentage = (clustered_count / total_genes_accounted_for) * 100 if total_genes_accounted_for > 0 else 0
        total_candidate_percentage = (
            (len(candidate_lsg_genes) / total_genes_accounted_for) * 100
            if total_genes_accounted_for > 0 else 0
        )
        verified_candidate_percentage = (
            (len(verified_candidate_genes) / total_genes_accounted_for) * 100
            if total_genes_accounted_for > 0 else 0
        )

        self._write_per_species_summary(species_out_dir, filename_csn, per_species_summary_data)
        return (
            gene_classification_data,
            c_lsg_percentage,
            total_candidate_percentage,
            verified_candidate_percentage,
            global_summary_node_data,
        )

###############################################################################

    def _generate_global_summary(self, c_lsg_summary_path, global_node_summary_path):
        """
        Merges the main summary files into a final global summary using CSVFileMerger.
        """
        logging.info("Generating final global summary file...")

        # Define the paths to the files to be merged
        stats_file = self.stats_converted_path

        if not stats_file.exists() or not c_lsg_summary_path.exists() or not global_node_summary_path.exists():
            logging.warning("One or more summary files are missing. Skipping final global merge.")
            return

        merger = CSVFileMerger(
            output_file=str(self.out_dir / "global_summary.csv"),
            # Define the merge column for each file in order
            merge_cols=["Species", "Species", "Species"],
            input_files=[
                str(stats_file),
                str(c_lsg_summary_path),
                str(global_node_summary_path)
            ],
            species_map_file=self.species_map_file
        )
        merger.merge_files()
        logging.info("Final global summary created successfully.")

    def _write_hog_disagreements_log(self, path):
        """Writes the HOG-vs-Orthogroups.tsv disagreement log (TSV). Always written, even if empty."""
        sorted_rows = sorted(self.hog_disagreements, key=lambda r: (r["species"], r["gene_id"]))
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["species", "gene_id", "orthogroups_level", "hog_level", "reason", "missing_species"]
            )
            writer.writeheader()
            writer.writerows(sorted_rows)

    def _write_stats_mismatches_log(self, path):
        """Writes the stats cross-validation mismatch log (TSV). Always written, even if empty."""
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["species", "our_clustered", "orthofinder_clustered", "our_unclustered", "orthofinder_unclustered"]
            )
            writer.writeheader()
            writer.writerows(sorted(self.stats_mismatches, key=lambda r: r["species"]))

    def _print_final_summary(self, global_summary_data):
        """
        Logs a per-species x per-node placement table directly to
        stdout/the log, so the result of a run is visible immediately
        without opening any CSV -- every column here sums to that
        species' total gene count, which is also worth seeing at a
        glance rather than re-deriving by hand.
        """
        if not global_summary_data:
            return

        all_nodes = sorted({node for node_data in global_summary_data.values() for node in node_data.keys()})
        col_width = max(6, max((len(n) for n in all_nodes), default=0))
        species_width = max(len("Species"), max((len(s) for s in global_summary_data), default=0))

        header = f"{'Species':<{species_width}}  " + "  ".join(f"{n:>{col_width}}" for n in all_nodes) + f"  {'Total':>{col_width}}"
        logging.info("Per-species node placement summary:")
        logging.info(header)
        logging.info("-" * len(header))
        for species in sorted(global_summary_data):
            node_data = global_summary_data[species]
            row = f"{species:<{species_width}}  " + "  ".join(
                f"{node_data.get(n, 0):>{col_width}}" for n in all_nodes
            ) + f"  {sum(node_data.values()):>{col_width}}"
            logging.info(row)

    def _print_validation_status(self):
        """
        A clear pass/fail line per validation check, matching the
        [PASS]/[FAIL] convention CLSGClassificationVerifier already uses --
        so a run's overall health is visible at a glance in the log/stdout,
        not just scattered across individual warning lines.
        """
        logging.info("Validation summary:")

        if self.stats_processor:
            if self.stats_mismatches:
                logging.info(f"  [FAIL] Stats validation: {len(self.stats_mismatches)} species mismatch "
                             f"OrthoFinder's own counts -- see stats_validation_mismatches.tsv")
            else:
                logging.info("  [PASS] Stats validation: computed counts match OrthoFinder's own stats for every species.")
        else:
            logging.info("  [SKIPPED] Stats validation: no stats file was available.")

        if self.hog_verifier:
            unexplained = sum(1 for d in self.hog_disagreements if d["reason"] == "unexplained")
            explained = len(self.hog_disagreements) - unexplained
            if unexplained:
                logging.info(f"  [FAIL] HOG cross-validation: {unexplained} unexplained disagreement(s) "
                             f"(plus {explained} explained by a missing species) -- see hog_disagreements.tsv")
            else:
                logging.info(f"  [PASS] HOG cross-validation: {explained} disagreement(s), all explained by a "
                             f"species missing from the HOG table.")
        else:
            logging.info("  [SKIPPED] HOG cross-validation: no HOG file was available or it was skipped.")

        if self.cross_checker:
            n_excluded = len(self.cross_checker.exclusions)
            n_checked = len(self.cross_checker.checked_species)
            n_total = n_checked + len(self.cross_checker.unchecked_species)
            logging.info(f"  [INFO] Pairwise cross-check: cross-checked {n_checked} of {n_total} species; "
                         f"{n_excluded} candidate gene(s) excluded -- see cross_check_exclusions.tsv")
            if self.cross_checker.unchecked_species:
                logging.warning("  Not cross-checked (no Orthologues_<species>/ directory): "
                                + ", ".join(sorted(self.cross_checker.unchecked_species)))
        else:
            logging.info("  [SKIPPED] Pairwise cross-check: no Orthologues/ directory was available or it was skipped.")

###############################################################################

    def run(self):
        """Runs the pipeline. Returns the process exit code (0 on success)."""
        logging.info("Starting OrthoFinder Data Processor pipeline.")
        # The mkdir line that was here has been moved to __init__
        try:
            self.initialize_processors()
            if not self.stratigraphy_classifier or not self.tree_parser:
                logging.error("Processors not initialized correctly. Cannot proceed."); return 1

            all_species_classifications, c_lsg_stats, global_summary_data = [], {}, {}
            all_species_filenames = sorted(self.stratigraphy_classifier.species_in_og_table)
            logging.info(f"Found {len(all_species_filenames)} species basenames to process.")

            for species_filename in all_species_filenames:
                gene_data, c_lsg_percentage, total_candidate_percentage, verified_candidate_percentage, node_data = \
                    self.process_and_classify_genes_for_species(species_filename)
                if species_filename in self.skipped_species:
                    continue
                all_species_classifications.extend(gene_data)
                c_lsg_stats[species_filename] = [
                    c_lsg_percentage, total_candidate_percentage, verified_candidate_percentage
                ]
                correct_species_name = self.species_map.get(species_filename, species_filename)
                global_summary_data[correct_species_name] = node_data

            self._write_gene_classification_summary(all_species_classifications)
            c_lsg_summary_path = self.generate_c_lsg_summary(c_lsg_stats)
            global_node_summary_path = self._write_global_node_summary(global_summary_data)

            self._generate_global_summary(c_lsg_summary_path, global_node_summary_path)

            if self.cross_checker:
                exclusion_log_path = self.out_dir / "cross_check_exclusions.tsv"
                self.cross_checker.write_exclusion_log(exclusion_log_path)
                logging.info(f"Cross-check excluded {len(self.cross_checker.exclusions)} candidate genes "
                             f"total across all species; log written to {exclusion_log_path}")

            hog_disagreements_path = self.out_dir / "hog_disagreements.tsv"
            self._write_hog_disagreements_log(hog_disagreements_path)
            if self.hog_disagreements:
                logging.info(f"{len(self.hog_disagreements)} gene(s) disagree between Orthogroups.tsv-derived "
                             f"and HOG-derived node placement; log written to {hog_disagreements_path}")

            stats_mismatches_path = self.out_dir / "stats_validation_mismatches.tsv"
            self._write_stats_mismatches_log(stats_mismatches_path)
            if self.stats_mismatches:
                logging.warning(f"{len(self.stats_mismatches)} species' computed counts disagree with "
                                 f"OrthoFinder's own stats; log written to {stats_mismatches_path}")

            self._print_final_summary(global_summary_data)
            self._print_validation_status()
            exit_code = self._report_problems()

        except PipelineInputError as e:
            logging.critical(f"Invalid input, stopping: {e}")
            return 1
        except FileNotFoundError as e:
            # Missing/misnamed input file -- initialize_processors() already
            # logged a specific, actionable message (what was expected, what
            # was actually found). The full stack trace adds nothing for
            # this case and just buries that message under noise, so it's
            # deliberately omitted here (unlike the genuine-bug case below).
            logging.critical(f"Cannot proceed: {e}")
            return 1
        except Exception as e:
            logging.critical(f"A critical error occurred in the main pipeline: {e}", exc_info=True)
            return 1

        logging.info("Pipeline finished.")
        return exit_code

    def _report_problems(self):
        """Logs skipped species and failed writes; returns the exit code (1 if the run must fail)."""
        exit_code = 0
        if self.skipped_species:
            level = logging.WARNING if self.allow_missing else logging.ERROR
            logging.log(level, f"{len(self.skipped_species)} species were skipped and are absent from the output:")
            for sp, reason in sorted(self.skipped_species.items()):
                logging.log(level, f"  {sp}: {reason}")
            if self.allow_missing:
                logging.warning("Continuing because --allow-missing was given.")
            else:
                logging.error("Exiting non-zero. Fix the inputs, or pass --allow-missing to accept the skipped species.")
                exit_code = 1
        if self.write_failures:
            logging.error(f"{len(self.write_failures)} output file(s) could not be written:")
            for path, err in self.write_failures:
                logging.error(f"  {path}: {err}")
            exit_code = 1
        return exit_code

    def _save_gene_ids(self, filepath, ids):
        try:
            with open(filepath, 'w') as f:
                for gene_id in sorted(list(ids)): f.write(f"{gene_id}\n")
        except IOError as e:
            logging.error(f"Error saving gene IDs to {filepath}: {e}")
            self.write_failures.append((filepath, e))

    def _filter_fasta(self, species_filename, file_type, output_dir, out_file_name, gene_ids_to_include):
        if not gene_ids_to_include:
            return

        input_dir = self.prt_dir # Assuming protein for now
        temp_gene_id_file = output_dir / f"temp_genes_to_filter_{out_file_name}.txt"
        self._save_gene_ids(temp_gene_id_file, gene_ids_to_include)

        file_ext = 'faa' if file_type == 'prt' else 'fna'
        fasta_file = Path(input_dir) / f"{species_filename}.{file_ext}"
        if not fasta_file.is_file():
            for ext in ['.fa', '.fasta', '.faa', '.fna']:
                fasta_file = Path(input_dir) / f"{species_filename}{ext}"
                if fasta_file.is_file(): break
            else:
                logging.warning(f"FASTA file not found for {species_filename} in {input_dir}, skipping filtering.")
                os.remove(temp_gene_id_file); return

        output_fasta_name = f"{species_filename}_{out_file_name}_proteins.{file_ext}"
        logging.info(f"Filtering {file_type} sequences for {species_filename} into {output_fasta_name}")

        extractor = GeneSequenceExtractor(str(temp_gene_id_file), str(fasta_file), str(output_dir), include=True, output_filename=output_fasta_name)
        extractor.filter_fasta()

        os.remove(temp_gene_id_file)

    def _write_gene_classification_summary(self, all_species_classifications):
        summary_file = self.out_dir / "C_LSG_Iso.tsv"
        header = ['Species', 'Gene_ID', 'C_LSG_Iso_Node', 'Orthogroup_ID']
        try:
            # Sort by (Species, Gene_ID) for reproducible, diffable output --
            # rows were previously in whatever order species/HOGs happened
            # to be processed in, not sorted.
            sorted_rows = sorted(all_species_classifications, key=lambda row: (row[0], row[1]))
            with open(summary_file, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(header)
                writer.writerows(sorted_rows)
            logging.info(f"Successfully wrote detailed gene summary to {summary_file}")
        except Exception as e:
            logging.error(f"Failed to write gene classification summary: {e}")
            self.write_failures.append((summary_file, e))

    def generate_c_lsg_summary(self, gene_stats):
        """
        gene_stats: {species_filename: [species_specific_percentage,
                     total_candidate_percentage, verified_candidate_percentage]}

        species_specific_percentage is the original, unchanged stat
        (clustered-private bucket only, as a % of all genes). The other two
        are a matched before/after pair over the SAME population (clustered
        + unclustered candidates combined) so the effect of the cross-check
        is actually visible and comparable, instead of silently mixing two
        different denominators under similar-looking column names.
        """
        logging.info("Generating c_lsg percentage summary file...")
        smry_file = self.intermediate_dir / 'c_lsg_summary.csv'
        header = [
            'Species',
            'species_specific_percentage',
            'total_candidate_percentage',
            'verified_candidate_percentage',
        ]
        with open(smry_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            for species_filename in sorted(gene_stats):
                raw_pct, total_candidate_pct, verified_pct = gene_stats[species_filename]
                correct_name = self.species_map.get(species_filename, species_filename)
                writer.writerow([
                    correct_name, f"{raw_pct:.2f}", f"{total_candidate_pct:.2f}", f"{verified_pct:.2f}"
                ])
        logging.info(f"c_lsg_summary.csv created at {smry_file}")
        return smry_file

def main():
    parser = argparse.ArgumentParser(description='Classify all genes of every species by phylogenetic depth of origin, using Orthogroups.tsv and the species tree.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('out_dir', help='Path to the output directory')
    parser.add_argument('--orthofinder_dir', help='Path to the parent OrthoFinder results directory.')
    parser.add_argument('--prt_dir', help='Path to protein FASTA files. Required: this is the full gene universe classification is done against.')
    parser.add_argument('--species-map', help='Path to a CSV/TSV file mapping filenames to correct species names. Must have "Species" and "Basename" headers.')
    parser.add_argument('--stats-filename', default='Statistics_PerSpecies.tsv',
                        help='OrthoFinder\'s own per-species stats, used to validate this pipeline\'s '
                             'computed species-specific counts against OrthoFinder\'s own accounting.')
    parser.add_argument('--tree-filename', default='SpeciesTree_rooted_node_labels.txt')
    parser.add_argument('--orthogroups-filename', default='Orthogroups.tsv',
                        help='OrthoFinder\'s flat, whole-dataset clustering table. This is the primary '
                             'classification source: unlike a Hierarchical Orthogroups (HOG) table, it '
                             'always includes every species OrthoFinder was run on, so it has no blind '
                             'spot for species outside a particular node\'s clade.')
    parser.add_argument('--hog-filename', default='N0.tsv',
                        help='Optional: a Hierarchical Orthogroups file (e.g. N0.tsv, N1.tsv) used ONLY '
                             'for secondary cross-validation against the primary Orthogroups.tsv-derived '
                             'classification, gene by gene. Not required -- if not found, HOG validation '
                             'is skipped with a warning, not an error.')
    parser.add_argument('--skip-hog-validation', action='store_true',
                        help='Skip HOG-based cross-validation entirely, without even trying to find a HOG file.')
    parser.add_argument('--orthologues-dirname', default='Orthologues',
                        help='OrthoFinder\'s pairwise Orthologues directory, used as a secondary cross-check '
                             '(reciprocal-best-hit) on the species-specific candidate bucket.')
    parser.add_argument('--skip-cross-check', action='store_true',
                        help='Skip the pairwise-Orthologues cross-check. Only use this if you understand '
                             'the consequence: candidate species-specific genes will not be verified '
                             'against reciprocal-best-hit orthology, so results may include genes with a '
                             'real ortholog that clustering missed.')
    parser.add_argument('--allow-missing', action='store_true',
                        help='Do not fail when a species in Orthogroups.tsv has no readable protein FASTA. '
                             'Such species are still reported and left out of the output.')
    args = parser.parse_args()
    processor = OrthoFinderDataProcessor(
        args.out_dir,
        args.orthofinder_dir,
        args.prt_dir,
        stats_filename=args.stats_filename,
        tree_filename=args.tree_filename,
        species_map_file=args.species_map,
        hog_filename=args.hog_filename,
        orthogroups_filename=args.orthogroups_filename,
        orthologues_dirname=args.orthologues_dirname,
        skip_cross_check=args.skip_cross_check,
        skip_hog_validation=args.skip_hog_validation,
        allow_missing=args.allow_missing,
    )
    sys.exit(processor.run())

if __name__ == "__main__":
    main()
