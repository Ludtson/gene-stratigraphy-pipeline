#!/usr/bin/env python3

"""
CLSGClassificationVerifier.py

Independently audits candidate lineage-specific gene (cLSG) classifications
produced by OrthoFinderDataProcessor.py, by re-deriving ground truth
directly from the OrthoFinder hierarchical orthogroup (HOG) table (e.g.
N0.tsv, N1.tsv, ...), rather than trusting the pipeline's own internal
bookkeeping. This exists because a real run (v7.1) produced "clustered"
cLSGs that were later found to co-occur, in the same HOG, with genes from
another species (L. cruciata) -- i.e. not actually species-specific.

This class does NOT reuse N0OrthoGeneSorter or any other pipeline code. It
re-parses the HOG table from scratch, so a bug shared between this script
and the pipeline cannot silently cancel out and hide a real error.

Two checks, run per species:

1. Clustered cLSGs (genes the pipeline assigned to a HOG it believes is
   private to one species): for every such gene, independently look up its
   HOG in the HOG table and confirm that HOG contains genes from that
   species ONLY -- no other species column has any entry in that row.

2. Unclustered cLSGs (genes the pipeline believes are entirely absent from
   the HOG table -- present in the species' protein FASTA but in no HOG at
   all): for every such gene, confirm it truly does not appear in ANY HOG
   row for that species, anywhere in the HOG table.

Note: this checks internal consistency between the pipeline's output and
the HOG table it was built from. It cannot detect relationships that are
structurally invisible to the HOG table itself -- e.g. a species outside
the clade covered by whichever N?.tsv was used (see
OrthoFinderDataProcessor's --hog-filename). That case is what the
mandatory Orthogroups.tsv / pairwise Orthologues cross-check in the main
pipeline is for, not this script.

Expected input layout (matches OrthoFinderDataProcessor.py's output):
  <results_dir>/<species>/gene_IDs_by_node/<species>_species_specific_clustered_gene_ids.txt
  <results_dir>/<species>/gene_IDs_by_node/<species>_species_specific_unclustered_gene_ids.txt

Usage:
  python CLSGClassificationVerifier.py \
    --hog-tsv N1.tsv \
    --results-dir clsg_v7.1 \
    --report violations_report.tsv

Exit status: 0 if no violations found, 1 if any violation is found (or if
no species could be checked at all), so this can be used as a CI/pipeline
gate, not just a manual report.

Author: Adekola Owoyemi (Casola Lab, Texas A&M University)
Version: 1.0.0
"""

import argparse
import csv
import sys
from pathlib import Path

# Increase field size limit for large bioinformatics files (HOG rows can have
# very long comma-separated gene lists). sys.maxsize overflows the C `long`
# used internally by csv.field_size_limit on platforms where C long is
# 32-bit (e.g. Windows, even on 64-bit Python), so cap at the largest value
# that's safe everywhere.
csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

__version__ = "1.0.0"
__author__ = "Adekola Owoyemi (Casola Lab, Texas A&M University)"


class CLSGClassificationVerifier:
    def __init__(self, hog_tsv_path, results_dir, species_filter=None):
        self.hog_tsv_path = hog_tsv_path
        self.results_dir = Path(results_dir)
        self.species_filter = set(species_filter) if species_filter else None

        # {hog_id: {species: set(genes)}}
        self.hog_species_genes = {}
        # {(species, gene_id): hog_id}
        self.gene_to_hog = {}
        self.species_in_hog_table = []

        self.violations = []
        self.total_clustered_checked = 0
        self.total_unclustered_checked = 0
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

    def _audit_clustered(self, species, species_dir):
        """
        Check 1: every gene the pipeline classified as a clustered cLSG must
        belong to a HOG that, independently, contains genes from this
        species only.
        """
        gene_file = species_dir / "gene_IDs_by_node" / f"{species}_species_specific_clustered_gene_ids.txt"
        claimed_genes = self._load_gene_id_list(gene_file)

        checked = 0
        for gene_id in claimed_genes:
            checked += 1
            hog_id = self.gene_to_hog.get((species, gene_id))

            if hog_id is None:
                self.violations.append({
                    "species": species,
                    "gene_id": gene_id,
                    "hog_id": "NOT_FOUND",
                    "issue": "clustered_cLSG_gene_missing_from_hog_table",
                    "other_species_in_hog": "",
                })
                continue

            species_genes = self.hog_species_genes.get(hog_id, {})
            other_species_present = sorted(
                sp for sp, genes in species_genes.items() if sp != species and genes
            )
            if other_species_present:
                self.violations.append({
                    "species": species,
                    "gene_id": gene_id,
                    "hog_id": hog_id,
                    "issue": "clustered_cLSG_HOG_shared_with_other_species",
                    "other_species_in_hog": ";".join(other_species_present),
                })

        return checked

    def _audit_unclustered(self, species, species_dir):
        """
        Check 2: every gene the pipeline classified as unclustered (absent
        from the HOG table entirely) must genuinely not appear in ANY HOG
        row, for this species, anywhere in the HOG table.
        """
        gene_file = species_dir / "gene_IDs_by_node" / f"{species}_species_specific_unclustered_gene_ids.txt"
        claimed_genes = self._load_gene_id_list(gene_file)

        checked = 0
        for gene_id in claimed_genes:
            checked += 1
            hog_id = self.gene_to_hog.get((species, gene_id))
            if hog_id is not None:
                self.violations.append({
                    "species": species,
                    "gene_id": gene_id,
                    "hog_id": hog_id,
                    "issue": "unclustered_cLSG_actually_present_in_a_HOG",
                    "other_species_in_hog": "",
                })

        return checked

    def run(self):
        """Parses the HOG table and audits every discovered (or requested) species."""
        print(f"[INFO] Parsing HOG table: {self.hog_tsv_path}")
        self.parse_hog_table()
        print(f"[INFO] Loaded {len(self.hog_species_genes)} HOGs across "
              f"{len(self.species_in_hog_table)} species columns.")

        species_dirs = self._find_species_dirs()
        if not species_dirs:
            raise ValueError(f"No species subdirectories found under {self.results_dir}.")

        print(f"[INFO] Auditing {len(species_dirs)} species: {', '.join(sorted(species_dirs))}")
        print()

        for species, species_dir in sorted(species_dirs.items()):
            if species not in self.species_in_hog_table:
                print(f"  [WARNING] '{species}' is not a column in the HOG table -- skipping (name mismatch?).")
                self.skipped_species.append(species)
                continue
            n_clustered = self._audit_clustered(species, species_dir)
            n_unclustered = self._audit_unclustered(species, species_dir)
            self.total_clustered_checked += n_clustered
            self.total_unclustered_checked += n_unclustered
            print(f"  {species}: checked {n_clustered} clustered cLSGs, {n_unclustered} unclustered cLSGs")

    def write_report(self, report_path):
        """Writes the violations report (TSV). Always written, even if empty."""
        with open(report_path, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=["species", "gene_id", "hog_id", "issue", "other_species_in_hog"]
            )
            writer.writeheader()
            writer.writerows(self.violations)

    def print_summary(self, report_path):
        print()
        print(f"[INFO] Total clustered cLSGs checked:   {self.total_clustered_checked}")
        print(f"[INFO] Total unclustered cLSGs checked: {self.total_unclustered_checked}")
        if self.skipped_species:
            print(f"[WARNING] Species skipped (name mismatch with HOG table): {', '.join(self.skipped_species)}")
        print(f"[INFO] Violations found: {len(self.violations)}")
        print(f"[INFO] Report written to: {report_path}")
        print()

        if self.violations:
            by_issue = {}
            for v in self.violations:
                by_issue[v["issue"]] = by_issue.get(v["issue"], 0) + 1
            for issue, count in sorted(by_issue.items()):
                print(f"  {issue}: {count}")
            print()
            print("[FAIL] cLSG classification does NOT match the HOG table. See report for details.")
        else:
            print("[PASS] All checked cLSGs are consistent with the HOG table.")


def main():
    parser = argparse.ArgumentParser(
        description="Independently verify candidate lineage-specific gene (cLSG) "
                     "classifications against the raw OrthoFinder HOG table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--hog-tsv", required=True,
        help="Path to the HOG table actually used for classification (e.g. N0.tsv or N1.tsv -- "
             "must match whatever --hog-filename the main pipeline was run with).",
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
        help="Path to write the violations report (TSV). Always written, even if empty.",
    )
    args = parser.parse_args()

    verifier = CLSGClassificationVerifier(args.hog_tsv, args.results_dir, species_filter=args.species)

    try:
        verifier.run()
    except (ValueError, FileNotFoundError) as e:
        print(f"[ERROR] {e}")
        sys.exit(1)

    verifier.write_report(args.report)
    verifier.print_summary(args.report)

    sys.exit(1 if verifier.violations else 0)


if __name__ == "__main__":
    main()
