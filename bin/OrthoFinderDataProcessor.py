#!/usr/bin/env python3

import os
import re
import argparse
import logging
import csv
from pathlib import Path

# Assuming helper scripts are in the same directory or in the python path
from OrthoFinderStatsProcessor import OrthoFinderStatsProcessor
from OrthoFinderN0GeneSorter import N0OrthoGeneSorter
from OrthoFinderTreeParser import SpeciesTree
from FastaGeneFilter import GeneSequenceExtractor
from CSVFileMerger import CSVFileMerger
from OrthologyCrossChecker import OrthologyCrossChecker

class OrthoFinderDataProcessor:
    def __init__(self, out_dir, orthofinder_dir=None, prt_dir=None,
                stats_filename="Statistics_PerSpecies.tsv",
                hog_filename="N0.tsv",
                tree_filename="SpeciesTree_rooted_node_labels.txt",
                orthogroups_filename="Orthogroups.tsv",
                orthologues_dirname="Orthologues",
                skip_cross_check=False,
                species_map_file=None):
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
        # The Orthogroups.tsv / pairwise Orthologues cross-check is mandatory
        # by default: the HOG table alone (N0.tsv/N1.tsv/etc.) is
        # structurally blind to any species outside the clade covered by
        # that specific node -- not just a nominal "outgroup" -- so a
        # candidate can look species-specific purely because the HOG table
        # can't see the species it actually has a relationship with. This
        # can only be skipped explicitly, not silently.
        self.skip_cross_check = skip_cross_check
        self.species_map = self._load_species_map()
        self.setup_logging()
        self.stats_processor = None
        self.n0_sorter = None
        self.tree_parser = None
        self.cross_checker = None

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
        Wraps _find_file for --hog-filename specifically: if the requested
        file (default N0.tsv) isn't there, check whether any other N#.tsv
        hierarchical orthogroup file exists and suggest it. This is the
        exact, common outgroup situation the README documents -- an
        outgroup pushes the true root out one level and N0.tsv is never
        written -- so a bare "not found" here is worse than useless, it
        actively hides the fix that's usually one flag away.
        """
        try:
            return self._find_file(self.orthofinder_dir, self.hog_filename)
        except FileNotFoundError:
            available = set()
            for root, _, files in os.walk(self.orthofinder_dir):
                for f in files:
                    if re.fullmatch(r'N\d+\.tsv', f):
                        available.add(f)
            if available:
                ordered = sorted(available, key=lambda n: int(n[1:-4]))
                raise FileNotFoundError(
                    f"{self.hog_filename} not found under {self.orthofinder_dir}, but "
                    f"found: {', '.join(ordered)}. If an outgroup was used and this "
                    f"version of OrthoFinder didn't write N0.tsv, pass the right one "
                    f"via --hog-filename (see the README's \"Why --hog-filename might "
                    f"not be N0.tsv\" -- check Species_Tree/SpeciesTree_rooted_node_"
                    f"labels.txt to confirm which node is correct, don't just guess)."
                )
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
                # section of the final merged global_summary.csv. OrthoFinder's
                # output layout has changed across versions before (e.g.
                # N0.tsv itself was dropped in v3.1.4+), so a missing/renamed
                # stats file shouldn't abort the whole run the way a missing
                # HOG table or species tree must -- those two are genuinely
                # required and still raise below.
                logging.warning(
                    f"{self.stats_filename} not found under {self.orthofinder_dir} -- "
                    f"SpeciesStats.csv will be skipped and global_summary.csv will be "
                    f"missing that data. If this OrthoFinder version renamed or moved "
                    f"it, pass the actual filename via --stats-filename."
                )
            hog_file = self._find_hog_file()
            self.n0_sorter = N0OrthoGeneSorter(hog_file, self.species_map_file)
            tree_file = self._find_file(self.orthofinder_dir, self.tree_filename)
            self.tree_parser = SpeciesTree(file_path=tree_file, species_map_file=self.species_map_file)
            self._initialize_cross_checker()
            logging.info("All processors initialized successfully.")
        except FileNotFoundError as e:
            logging.error(f"Failed to initialize processors: {e}"); raise

    def _initialize_cross_checker(self):
        """
        Locates Orthogroups.tsv and the Orthologues/ directory the same way
        every other required file is located (recursive search under
        orthofinder_dir), and builds the mandatory cross-check. Only skipped
        if skip_cross_check was explicitly requested, or if the files
        genuinely can't be found -- either way, this is always visible in
        the log, never silent.
        """
        if self.skip_cross_check:
            logging.warning(
                "Cross-check explicitly skipped (--skip-cross-check). Candidate cLSGs will NOT be "
                "verified against Orthogroups.tsv or pairwise Orthologues -- results may include "
                "genes that have a real ortholog outside the HOG table's visible clade."
            )
            return

        try:
            orthogroups_file = self._find_file(self.orthofinder_dir, self.orthogroups_filename)
        except FileNotFoundError:
            orthogroups_file = None
            logging.warning(f"{self.orthogroups_filename} not found under {self.orthofinder_dir} -- "
                             f"primary cross-check will be skipped.")

        try:
            orthologues_dir = self._find_dir(self.orthofinder_dir, self.orthologues_dirname)
        except FileNotFoundError:
            orthologues_dir = None
            logging.warning(f"{self.orthologues_dirname}/ not found under {self.orthofinder_dir} -- "
                             f"secondary pairwise cross-check will be skipped.")

        if not orthogroups_file and not orthologues_dir:
            logging.error(
                "Neither Orthogroups.tsv nor Orthologues/ could be found. The mandatory cross-check "
                "cannot run at all. Pass --skip-cross-check explicitly if this is intentional, "
                "otherwise candidate cLSGs will not be checked against anything but the HOG table."
            )
            return

        self.cross_checker = OrthologyCrossChecker(orthogroups_file, orthologues_dir)

    def _get_all_gene_ids_from_fasta(self, species_basename):
        if not self.prt_dir: return set()
        fasta_file = Path(self.prt_dir) / f"{species_basename}.prt"
        if not fasta_file.is_file():
            for ext in ['.fa', '.fasta', '.faa', '.fna']:
                fasta_file = Path(self.prt_dir) / f"{species_basename}{ext}"
                if fasta_file.is_file(): break
            else:
                logging.warning(f"Protein FASTA file not found for {species_basename}, cannot determine unclustered genes.")
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
            logging.error(f"Could not read FASTA file {fasta_file}: {e}"); return set()

    def _write_per_species_summary(self, species_out_dir, filename_csn, summary_data):
        """Writes the summary CSV for a single species."""
        summary_file = species_out_dir / f'{filename_csn}_summary_per_node.csv'
        header = ['Node', 'Number_of_HOGs', 'Number_of_Genes']
        try:
            with open(summary_file, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(header)
                writer.writerows(summary_data)
            logging.info(f"Per-species summary written to {summary_file.name}")
        except IOError as e:
            logging.error(f"Failed to write per-species summary for {filename_csn}: {e}")
            
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

        return summary_file 

    def process_and_classify_genes_for_species(self, species_basename):
        correct_species_name = self.species_map.get(species_basename, species_basename)
        filename_csn = correct_species_name.replace(" ", "_")
        logging.info(f"--- Processing: {species_basename} -> {correct_species_name} ---")

        gene_classification_data = []
        per_species_summary_data = []
        global_summary_node_data = {}

        try:
            nodes_up_tree = [species_basename] + self.tree_parser.find_nodes_up_tree(species_basename)
        except ValueError:
            # find_nodes_up_tree raises (not returns []) when species_basename
            # isn't a node in the tree -- e.g. a species-name mismatch between
            # the HOG table and the tree file. Catch it here so one bad name
            # only skips this species, rather than propagating uncaught out
            # of run()'s per-species loop and aborting every other species.
            logging.warning(f"'{species_basename}' not found in the species tree. Skipping this species.")
            return [], 0, 0, 0, {}
            
        HOG_ID_levels = {}
        for node in nodes_up_tree:
            if node == self.tree_parser.root.name: continue
            sister_species = self.tree_parser.find_sister_species_to_node(node)
            parent_node = self.tree_parser.get_parent(node)
            if parent_node: HOG_ID_levels[parent_node] = self.n0_sorter.get_shared_hog_ids(species_basename, sister_species)

        for node in HOG_ID_levels:
            for other_node in HOG_ID_levels:
                if node != other_node and self.tree_parser.is_ancestor(other_node, node):
                    HOG_ID_levels[node] -= HOG_ID_levels[other_node]
        
        HOG_ID_levels['species_specific_clustered'] = self.n0_sorter.get_unique_hog_ids(species_basename)

        species_out_dir = self.out_dir / filename_csn
        gene_ids_dir, hog_ids_dir, protein_dir = species_out_dir / 'gene_IDs_by_node', species_out_dir / 'HOG_IDs_by_node', species_out_dir / 'proteins_by_node'
        gene_ids_dir.mkdir(parents=True, exist_ok=True)
        hog_ids_dir.mkdir(parents=True, exist_ok=True)
        protein_dir.mkdir(parents=True, exist_ok=True)

        all_clustered_genes, c_lsg_count = set(), 0
        candidate_lsg_genes = set()
        candidate_status = {}
        for node_name, hog_set in HOG_ID_levels.items():
            gene_ids = self.n0_sorter.get_genes_for_hog_ids(species_basename, hog_set)
            all_clustered_genes.update(gene_ids)

            if node_name == 'species_specific_clustered':
                candidate_lsg_genes.update(gene_ids)
                # HOG IDs for this bucket are written below in the per-node
                # loop as {species}_species_specific_clustered_HOG_ids.txt --
                # no separate tracking needed here.
                for gene_id in gene_ids:
                    candidate_status[gene_id] = 'species_specific_clustered'

            # Populate data for summary tables
            per_species_summary_data.append([node_name, len(hog_set), len(gene_ids)])
            global_summary_node_data[node_name] = len(gene_ids)

            self._save_gene_ids(gene_ids_dir / f'{species_basename}_{node_name}_gene_ids.txt', gene_ids)
            self._save_gene_ids(hog_ids_dir / f'{species_basename}_{node_name}_HOG_ids.txt', hog_set)

            for hog_id in hog_set:
                genes_in_hog = self.n0_sorter.get_genes_for_hog_ids(species_basename, {hog_id})
                for gene_id in genes_in_hog:
                    gene_classification_data.append([correct_species_name, gene_id, node_name, hog_id])
            
            # Per-node FASTA disabled (collapsed output mode)
            pass

            if node_name == 'species_specific_clustered':
                c_lsg_count = len(gene_ids)

        self._save_gene_ids(gene_ids_dir / f'{species_basename}_all_clustered_ids.txt', all_clustered_genes)

        all_genes_from_fasta = self._get_all_gene_ids_from_fasta(species_basename)
        unclustered_genes = set()
        if all_genes_from_fasta:
            unclustered_genes = all_genes_from_fasta - all_clustered_genes
            logging.info(f"Identified {len(unclustered_genes)} unclustered genes.")
            self._save_gene_ids(gene_ids_dir / f'{species_basename}_species_specific_unclustered_gene_ids.txt', unclustered_genes)
            per_species_summary_data.append(["species_specific_unclustered", "NA", len(unclustered_genes)])
            global_summary_node_data["species_specific_unclustered"] = len(unclustered_genes)
            for gene_id in unclustered_genes:
                gene_classification_data.append([correct_species_name, gene_id, 'species_specific_unclustered', 'NA'])
                candidate_status[gene_id] = 'species_specific_unclustered'
            candidate_lsg_genes.update(unclustered_genes)
        else:
            logging.warning("No protein FASTA provided; cannot determine unclustered genes.")

        # Mandatory cross-check: a gene that looks species-specific from the
        # HOG table alone may still have a real ortholog in a species that
        # table structurally cannot see (see OrthologyCrossChecker). Verify
        # every candidate before it's treated as final.
        if self.cross_checker:
            verified_candidate_genes = self.cross_checker.filter_candidates(
                species_basename, candidate_lsg_genes, original_status_lookup=candidate_status
            )
            n_excluded = len(candidate_lsg_genes) - len(verified_candidate_genes)
            if n_excluded:
                logging.info(f"Cross-check excluded {n_excluded} of {len(candidate_lsg_genes)} "
                             f"candidate genes for {species_basename} (had an ortholog the HOG "
                             f"table couldn't see).")
        else:
            verified_candidate_genes = candidate_lsg_genes

        self._save_gene_ids(
            gene_ids_dir / f'{species_basename}_candidate_lsg_gene_ids.txt',
            verified_candidate_genes
        )

        if self.prt_dir and verified_candidate_genes:
            logging.info(f"Writing combined candidate FASTA for {species_basename}")

            self._filter_fasta(
                species_basename,
                "prt",
                protein_dir,
                "candidate_LSGs",
                verified_candidate_genes
            )

        total_genes_accounted_for = len(all_clustered_genes) + len(unclustered_genes)
        # Kept exactly as before: % of genes in the clustered-private bucket
        # ONLY (does not include unclustered genes). Unchanged meaning, so
        # anything already relying on this column isn't silently redefined.
        c_lsg_percentage = (c_lsg_count / total_genes_accounted_for) * 100 if total_genes_accounted_for > 0 else 0
        # These two are a matched pair -- same population (clustered-private
        # + unclustered candidates combined), before and after the
        # cross-check -- so they're actually comparable, unlike mixing
        # c_lsg_percentage (clustered-only) with a post-filter total would be.
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

###############################################################################

    def run(self):
        logging.info("Starting OrthoFinder Data Processor pipeline.")
        # The mkdir line that was here has been moved to __init__
        try:
            self.initialize_processors()
            if not self.n0_sorter or not self.tree_parser:
                logging.error("Processors not initialized correctly. Cannot proceed."); return
            
            all_species_classifications, c_lsg_stats, global_summary_data = [], {}, {}
            all_species_filenames = self.n0_sorter.all_species
            logging.info(f"Found {len(all_species_filenames)} species basenames to process.")

            for species_filename in all_species_filenames:
                gene_data, c_lsg_percentage, total_candidate_percentage, verified_candidate_percentage, node_data = \
                    self.process_and_classify_genes_for_species(species_filename)
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

        except Exception as e:
            logging.critical(f"A critical error occurred in the main pipeline: {e}", exc_info=True)
        
        logging.info("Pipeline finished.")
    def _save_gene_ids(self, filepath, ids):
        try:
            with open(filepath, 'w') as f:
                for gene_id in sorted(list(ids)): f.write(f"{gene_id}\n")
        except IOError as e:
            logging.error(f"Error saving gene IDs to {filepath}: {e}")
    
    def _filter_fasta(self, species_filename, file_type, output_dir, out_file_name, gene_ids_to_include):
        if not gene_ids_to_include:
            return
        
        input_dir = self.prt_dir # Assuming protein for now
        temp_gene_id_file = output_dir / "temp_genes_to_filter.txt"
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
        header = ['Species', 'Gene_ID', 'C_LSG_Iso_Node', 'HOG_ID']
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
    parser = argparse.ArgumentParser(description='Process OrthoFinder data to classify all genes by phylogenetic origin.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('out_dir', help='Path to the output directory')
    parser.add_argument('--orthofinder_dir', help='Path to the parent OrthoFinder results directory.')
    parser.add_argument('--prt_dir', help='Path to protein FASTA files. Required to identify unclustered genes.')
    parser.add_argument('--species-map', help='Path to a CSV/TSV file mapping filenames to correct species names. Must have "Species" and "Basename" headers.')
    parser.add_argument('--stats-filename', default='Statistics_PerSpecies.tsv')
    parser.add_argument('--hog-filename', default='N0.tsv',
                        help='Which hierarchical orthogroup file to classify from, e.g. N0.tsv, '
                             'N1.tsv, N2.tsv. Use a non-N0 file to analyze a specific clade/node '
                             'rather than the whole species tree -- e.g. if an outgroup was added '
                             'because OrthoFinder did not emit N0.tsv, or to focus on a shallower '
                             'node. Whichever file is used, every species outside that node\'s '
                             'clade is structurally invisible to this classification step, which '
                             'is exactly what --orthogroups-filename/--orthologues-dirname below '
                             'compensate for.')
    parser.add_argument('--tree-filename', default='SpeciesTree_rooted_node_labels.txt')
    parser.add_argument('--orthogroups-filename', default='Orthogroups.tsv',
                        help='OrthoFinder\'s flat, whole-dataset clustering table, used as the '
                             'primary mandatory cross-check against candidate cLSGs.')
    parser.add_argument('--orthologues-dirname', default='Orthologues',
                        help='OrthoFinder\'s pairwise Orthologues directory, used as the '
                             'secondary, higher-precision cross-check for candidates that survive '
                             'the primary one.')
    parser.add_argument('--skip-cross-check', action='store_true',
                        help='Skip the mandatory Orthogroups.tsv/pairwise-Orthologues cross-check. '
                             'Only use this if you understand the consequence: candidate cLSGs will '
                             'not be verified against anything outside the HOG table used by '
                             '--hog-filename, so results may include genes with a real ortholog '
                             'that table cannot see.')
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
    )
    processor.run()

if __name__ == "__main__":
    main()
