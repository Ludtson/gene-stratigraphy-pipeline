#!/usr/bin/env python3

"""
CLSGClassificationVerifier.py

Independently audits the gene classifications produced by
OrthoFinderDataProcessor.py. Ground truth is re-derived from the same inputs
the classification uses -- OrthoFinder's Orthogroups.tsv and the rooted
species tree -- rather than from the pipeline's own bookkeeping.

Orthogroups.tsv is re-parsed here from scratch and the expected node of each
gene is re-derived with the same MRCA rule as the classifier (species tree
parsing is shared via SpeciesTree; no classifier code is reused), so a bug in
the classification logic cannot cancel out and hide a real error.

Primary check, per species: every gene listed in a
gene_IDs_by_node/<species>_<node>_gene_ids.txt file must belong to that node
according to the re-derivation:
  - absent from Orthogroups.tsv                 -> species_specific_unclustered
  - orthogroup contains this species only       -> species_specific_clustered
  - otherwise MRCA of this species and every other species in its orthogroup

Optional secondary check (--hog-tsv): the old HOG-table consistency checks
(a "clustered" gene's HOG is private to the species; an "unclustered" gene is
in no HOG). A HOG table cannot see species outside the clade it covers, so
these findings are reported as notes and never affect the exit status.

Expected input layout (matches OrthoFinderDataProcessor.py's output):
  <results_dir>/<species>/gene_IDs_by_node/<species>_<node>_gene_ids.txt

Usage:
  python CLSGClassificationVerifier.py \
    --orthogroups-tsv Orthogroups.tsv \
    --tree SpeciesTree_rooted_node_labels.txt \
    --results-dir results \
    --report violations_report.tsv

Exit status: 0 if no violations found, 1 if any violation is found (or if
no species could be checked at all), so this can be used as a CI/pipeline
gate, not just a manual report.

Author: Adekola Owoyemi (Casola Lab, Texas A&M University)
Version: 2.0.0
"""

import argparse
import csv
import sys
from pathlib import Path

from OrthoFinderTreeParser import SpeciesTree

# Increase field size limit for large bioinformatics files (HOG rows can have
# very long comma-separated gene lists). sys.maxsize overflows the C `long`
# used internally by csv.field_size_limit on platforms where C long is
# 32-bit (e.g. Windows, even on 64-bit Python), so cap at the largest value
# that's safe everywhere.
csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

__version__ = "2.0.0"
__author__ = "Adekola Owoyemi (Casola Lab, Texas A&M University)"


class CLSGClassificationVerifier:
    def __init__(self, hog_tsv_path, results_dir, species_filter=None,
                 orthogroups_tsv_path=None, tree_path=None):
        # hog_tsv_path is optional (secondary check). OrthoFinderDataProcessor
        # also uses this class purely as a HOG-table parser.
        self.hog_tsv_path = hog_tsv_path
        self.orthogroups_tsv_path = orthogroups_tsv_path
        self.tree_path = tree_path
        self.results_dir = Path(results_dir)
        self.species_filter = set(species_filter) if species_filter else None

        # {(species, gene_id): og_id}, {og_id: {species: set(genes)}}
        self.og_gene_to_og = {}
        self.og_species_genes = {}
        self.species_in_og_table = []
        self.tree = None
        self.notes = []  # secondary HOG findings; never affect exit status

        # {hog_id: {species: set(genes)}}
        self.hog_species_genes = {}
        # {(species, gene_id): hog_id}
        self.gene_to_hog = {}
        self.species_in_hog_table = []
        # Species with at least one actual (non-empty) HOG entry, as opposed
        # to species_in_hog_table which is just the header row -- a species
        # can be a header column with every cell empty (e.g. an outgroup a
        # HOG file is structurally blind to), which is a fundamentally
        # different situation from a single gene inexplicably missing. See
        # _audit_clustered.
        self.species_with_hog_data = set()

        self.violations = []
        self.total_genes_checked = 0
        self.skipped_species = []

    def parse_hog_table(self):
        """
        Parses the HOG table from scratch, independent of any other class in
        this pipeline. Populates self.hog_species_genes, self.gene_to_hog,
        and self.species_in_hog_table.
        """
        with open(self.hog_tsv_path, "r", newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration:
                raise ValueError(f"HOG table is empty: {self.hog_tsv_path}")

            if len(header) < 4:
                raise ValueError(
                    f"HOG table {self.hog_tsv_path} has no species columns "
                    f"(expected 'HOG', 'OG', 'Gene Tree Parent Clade', then one column per species)."
                )
            self.species_in_hog_table = header[3:]

            for row in reader:
                if not row:
                    continue
                hog_id = row[0].strip()
                if not hog_id:
                    continue

                species_genes = {}
                for i, species in enumerate(self.species_in_hog_table):
                    col_idx = i + 3
                    if col_idx >= len(row):
                        continue
                    cell = row[col_idx].strip()
                    if not cell or cell == "-":
                        continue
                    genes = [g.strip() for g in cell.split(",") if g.strip()]
                    if genes:
                        species_genes[species] = set(genes)
                        for gene_id in genes:
                            key = (species, gene_id)
                            if key in self.gene_to_hog and self.gene_to_hog[key] != hog_id:
                                print(
                                    f"[WARNING] Gene '{gene_id}' for species '{species}' appears in "
                                    f"more than one HOG ({self.gene_to_hog[key]} and {hog_id}) -- "
                                    f"keeping the first one seen. This may indicate a duplicate "
                                    f"gene ID or a malformed HOG table."
                                )
                                continue
                            self.gene_to_hog[key] = hog_id

                self.hog_species_genes[hog_id] = species_genes

        self.species_with_hog_data = {
            sp for species_genes in self.hog_species_genes.values()
            for sp, genes in species_genes.items() if genes
        }

    def parse_orthogroups(self):
        """Parses Orthogroups.tsv from scratch (independent of the classifier)."""
        with open(self.orthogroups_tsv_path, "r", newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader, None)
            if not header or len(header) < 2:
                raise ValueError(f"Orthogroups table {self.orthogroups_tsv_path} is empty or has no species columns.")
            self.species_in_og_table = header[1:]
            for row in reader:
                if not row or not row[0].strip():
                    continue
                og_id = row[0].strip()
                species_genes = {}
                for i, species in enumerate(self.species_in_og_table):
                    if i + 1 >= len(row):
                        continue
                    genes = [g.strip() for g in row[i + 1].split(",") if g.strip()]
                    if genes:
                        species_genes[species] = set(genes)
                        for gene_id in genes:
                            self.og_gene_to_og.setdefault((species, gene_id), og_id)
                self.og_species_genes[og_id] = species_genes

    def expected_node(self, species, gene_id):
        """Node the gene should be in, re-derived from Orthogroups.tsv + tree. Raises ValueError if unplaceable."""
        og_id = self.og_gene_to_og.get((species, gene_id))
        if og_id is None:
            return "species_specific_unclustered"
        others = sorted(sp for sp, genes in self.og_species_genes[og_id].items() if sp != species and genes)
        if not others:
            return "species_specific_clustered"
        return self.tree.get_mrca([species] + others)

    @staticmethod
    def _load_gene_id_list(path):
        """Loads a newline-delimited gene ID file. Returns an empty set if the file doesn't exist."""
        ids = set()
        if not path.is_file():
            return ids
        with open(path, "r") as fh:
            for line in fh:
                gene_id = line.strip()
                if gene_id:
                    ids.add(gene_id)
        return ids

    def _find_species_dirs(self):
        """
        Each species subdirectory of the pipeline's output is expected to
        contain a gene_IDs_by_node/ subdirectory. Returns {species_name: Path}.
        """
        if not self.results_dir.is_dir():
            raise FileNotFoundError(f"Results directory not found: {self.results_dir}")

        species_dirs = {}
        for entry in sorted(self.results_dir.iterdir()):
            if entry.is_dir() and (entry / "gene_IDs_by_node").is_dir():
                species_dirs[entry.name] = entry

        if self.species_filter:
            species_dirs = {s: d for s, d in species_dirs.items() if s in self.species_filter}
            missing = self.species_filter - set(species_dirs)
            for s in sorted(missing):
                print(f"[WARNING] Requested species '{s}' has no output subdirectory under {self.results_dir}.")

        return species_dirs

    def _audit_nodes(self, species, species_dir):
        """
        Primary check: every gene in each <species>_<node>_gene_ids.txt must
        have that node as its re-derived expected node. Returns genes checked.
        """
        prefix, suffix = f"{species}_", "_gene_ids.txt"
        checked = 0
        for path in sorted((species_dir / "gene_IDs_by_node").glob(f"{prefix}*{suffix}")):
            claimed_node = path.name[len(prefix):-len(suffix)]
            if claimed_node == "candidate_lsg":
                continue  # union of every node except the root, not a node itself
            for gene_id in sorted(self._load_gene_id_list(path)):
                checked += 1
                try:
                    expected = self.expected_node(species, gene_id)
                except ValueError as e:
                    self.violations.append({
                        "species": species, "gene_id": gene_id, "claimed_node": claimed_node,
                        "expected_node": "", "issue": "node_not_derivable", "detail": str(e),
                    })
                    continue
                if expected != claimed_node:
                    self.violations.append({
                        "species": species, "gene_id": gene_id, "claimed_node": claimed_node,
                        "expected_node": expected, "issue": "node_mismatch",
                        "detail": self.og_gene_to_og.get((species, gene_id), "no_orthogroup"),
                    })
        return checked

    def _hog_notes(self, species, species_dir):
        """
        Secondary, informational: HOG-table consistency of the two
        species-specific buckets. Findings go to self.notes only.
        """
        gdir = species_dir / "gene_IDs_by_node"
        clustered = self._load_gene_id_list(gdir / f"{species}_species_specific_clustered_gene_ids.txt")
        unclustered = self._load_gene_id_list(gdir / f"{species}_species_specific_unclustered_gene_ids.txt")
        if species not in self.species_with_hog_data:
            return  # species invisible to this HOG table; nothing comparable
        for gene_id in sorted(clustered):
            hog_id = self.gene_to_hog.get((species, gene_id))
            if hog_id is None:
                self.notes.append({"species": species, "gene_id": gene_id, "claimed_node": "species_specific_clustered",
                                   "expected_node": "", "issue": "hog_note_clustered_gene_missing_from_hog_table",
                                   "detail": ""})
                continue
            others = sorted(sp for sp, g in self.hog_species_genes.get(hog_id, {}).items() if sp != species and g)
            if others:
                self.notes.append({"species": species, "gene_id": gene_id, "claimed_node": "species_specific_clustered",
                                   "expected_node": "", "issue": "hog_note_clustered_hog_shared_with_other_species",
                                   "detail": f"{hog_id}:{';'.join(others)}"})
        for gene_id in sorted(unclustered):
            hog_id = self.gene_to_hog.get((species, gene_id))
            if hog_id is not None:
                self.notes.append({"species": species, "gene_id": gene_id, "claimed_node": "species_specific_unclustered",
                                   "expected_node": "", "issue": "hog_note_unclustered_gene_present_in_a_hog",
                                   "detail": hog_id})

    def run(self):
        """Re-derives expected nodes and audits every discovered (or requested) species."""
        if not self.orthogroups_tsv_path or not self.tree_path:
            raise ValueError("Both an Orthogroups.tsv and a species tree are required for verification.")
        print(f"[INFO] Parsing Orthogroups table: {self.orthogroups_tsv_path}")
        self.parse_orthogroups()
        self.tree = SpeciesTree(file_path=self.tree_path)
        not_in_tree = sorted(set(self.species_in_og_table) - self.tree.get_all_species())
        if not_in_tree:
            raise ValueError(f"Species in Orthogroups table but not leaves of the tree: {', '.join(not_in_tree)}")
        print(f"[INFO] Loaded {len(self.og_species_genes)} orthogroups across "
              f"{len(self.species_in_og_table)} species.")
        if self.hog_tsv_path:
            self.parse_hog_table()
            print(f"[INFO] Loaded {len(self.hog_species_genes)} HOGs (secondary, informational).")

        species_dirs = self._find_species_dirs()
        if not species_dirs:
            raise ValueError(f"No species subdirectories found under {self.results_dir}.")

        print(f"[INFO] Auditing {len(species_dirs)} species: {', '.join(sorted(species_dirs))}")
        print()

        for species, species_dir in sorted(species_dirs.items()):
            if species not in self.species_in_og_table:
                print(f"  [WARNING] '{species}' is not a column in Orthogroups.tsv -- skipping (name mismatch?).")
                self.skipped_species.append(species)
                continue
            n_checked = self._audit_nodes(species, species_dir)
            self.total_genes_checked += n_checked
            if self.hog_tsv_path:
                self._hog_notes(species, species_dir)
            print(f"  {species}: checked {n_checked} genes")

    def write_report(self, report_path):
        """Writes the report (TSV). Violations first, then HOG notes. Always written, even if empty."""
        with open(report_path, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, delimiter="\t", lineterminator="\n",
                fieldnames=["level", "species", "gene_id", "claimed_node", "expected_node", "issue", "detail"],
            )
            writer.writeheader()
            for level, rows in (("violation", self.violations), ("note", self.notes)):
                for r in rows:
                    writer.writerow({"level": level, **r})

    def print_summary(self, report_path):
        print()
        print(f"[INFO] Total genes checked: {self.total_genes_checked}")
        if self.skipped_species:
            print(f"[WARNING] Species skipped (name mismatch with Orthogroups.tsv): {', '.join(self.skipped_species)}")
        print(f"[INFO] Violations found: {len(self.violations)}")
        if self.hog_tsv_path:
            print(f"[INFO] Secondary HOG notes (informational): {len(self.notes)}")
        print(f"[INFO] Report written to: {report_path}")
        print()

        if self.violations:
            by_issue = {}
            for v in self.violations:
                by_issue[v["issue"]] = by_issue.get(v["issue"], 0) + 1
            for issue, count in sorted(by_issue.items()):
                print(f"  {issue}: {count}")
            print()
            print("[FAIL] Classification does NOT match the Orthogroups.tsv + tree re-derivation. See report for details.")
        else:
            print("[PASS] All checked genes match the Orthogroups.tsv + tree re-derivation.")


def main():
    parser = argparse.ArgumentParser(
        description="Independently verify gene classifications by re-deriving each gene's node "
                     "from Orthogroups.tsv and the species tree.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--orthogroups-tsv", required=True,
        help="Path to OrthoFinder's Orthogroups.tsv (the file the classification was run with).",
    )
    parser.add_argument(
        "--tree", required=True,
        help="Path to the rooted species tree with labelled internal nodes "
             "(SpeciesTree_rooted_node_labels.txt).",
    )
    parser.add_argument(
        "--hog-tsv", default=None,
        help="Optional: a HOG table (e.g. N1.tsv) for a secondary, informational consistency check. "
             "Its findings are reported as notes and never change the exit status.",
    )
    parser.add_argument(
        "--results-dir", required=True,
        help="Path to the pipeline's output directory (contains one subdirectory per species, "
             "each with a gene_IDs_by_node/ folder).",
    )
    parser.add_argument(
        "--species", nargs="*", default=None,
        help="Optional: restrict the audit to these species basenames only. "
             "Default: every species subdirectory found under --results-dir.",
    )
    parser.add_argument(
        "--report", default="clsg_verification_report.tsv",
        help="Path to write the report (TSV). Always written, even if empty.",
    )
    args = parser.parse_args()

    verifier = CLSGClassificationVerifier(
        args.hog_tsv, args.results_dir, species_filter=args.species,
        orthogroups_tsv_path=args.orthogroups_tsv, tree_path=args.tree,
    )

    try:
        verifier.run()
    except (ValueError, FileNotFoundError, OSError) as e:
        print(f"[ERROR] {e}")
        sys.exit(1)

    verifier.write_report(args.report)
    verifier.print_summary(args.report)

    sys.exit(1 if verifier.violations else 0)


if __name__ == "__main__":
    main()
