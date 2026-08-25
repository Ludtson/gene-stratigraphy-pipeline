#!/usr/bin/env python3

"""
OrthologyCrossChecker.py

Cross-checks candidate lineage-specific genes (cLSGs) produced from the
Hierarchical Orthogroups (HOG) table against two other, independent
OrthoFinder outputs, and reports which candidates have a conflicting
ortholog elsewhere.

Why this exists: the HOG table used by OrthoFinderDataProcessor
(N0.tsv/N1.tsv/etc., selected via --hog-filename) is structurally blind to
any species outside the clade covered by that specific node. If an
outgroup is used to work around OrthoFinder not emitting N0.tsv (see the
main pipeline's README), every species outside the analyzed clade -- not
just a nominal "outgroup" -- is invisible to the HOG-based method, for any
node file chosen. A gene can therefore look like a candidate
lineage-specific gene purely because the HOG table cannot see the species
it actually has a relationship with.

Two independent checks, in order:

1. Primary: Orthogroups.tsv (OrthoFinder's flat, whole-dataset clustering).
   Same species-column table shape as the HOG files, but always includes
   every species run through OrthoFinder, including any outgroup. If a
   candidate gene's Orthogroup row contains any other species, it is not
   truly species-specific.

2. Secondary: pairwise Orthologues_<species>/<species>__v__<other>.tsv
   files (reciprocal-best-hit, tree-refined). Only checked for candidates
   that already survived the primary check, as a higher-precision second
   pass -- these files are per-species-pair, so this check is more
   expensive to run than the primary one.

A candidate flagged by either check is excluded, and the reason is
recorded so the exclusion is auditable, not silent.

Author: Adekola Owoyemi (Casola Lab, Texas A&M University)
Version: 1.0.0
"""

import csv
import sys
from pathlib import Path

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

__version__ = "1.0.0"
__author__ = "Adekola Owoyemi (Casola Lab, Texas A&M University)"


class OrthologyCrossChecker:
    def __init__(self, orthogroups_tsv_path=None, orthologues_dir=None):
        """
        orthogroups_tsv_path: path to OrthoFinder's Orthogroups.tsv.
        orthologues_dir: path to OrthoFinder's Orthologues/ directory,
            containing one Orthologues_<species>/ subdirectory per species,
            each with <species>__v__<other>.tsv pairwise files.
        Either or both may be omitted; whichever checks have their input
        provided are run, and the rest are skipped (with a warning), rather
        than failing outright -- useful for a rerun where only one file is
        available.
        """
        self.orthogroups_tsv_path = Path(orthogroups_tsv_path) if orthogroups_tsv_path else None
        self.orthologues_dir = Path(orthologues_dir) if orthologues_dir else None

        # {orthogroup_id: {species: set(genes)}}
        self.og_species_genes = {}
        # {(species, gene_id): orthogroup_id}
        self.gene_to_og = {}
        self.species_in_og_table = []

        self._orthogroups_loaded = False
        # Cache of parsed pairwise files, keyed by species, to avoid
        # re-reading the same file for every gene of that species.
        self._pairwise_cache = {}

        # Exclusion log entries: {gene_id, species, original_status,
        # check, conflicting_species, conflicting_gene, source_id}
        self.exclusions = []

    def load_orthogroups_table(self):
        """
        Parses Orthogroups.tsv. Same shape as the HOG tables, but with a
        single leading metadata column ("Orthogroup") instead of three
        ("HOG", "OG", "Gene Tree Parent Clade").
        """
        if self._orthogroups_loaded or not self.orthogroups_tsv_path:
            return

        with open(self.orthogroups_tsv_path, "r", newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration:
                raise ValueError(f"Orthogroups table is empty: {self.orthogroups_tsv_path}")

            if len(header) < 2:
                raise ValueError(
                    f"Orthogroups table {self.orthogroups_tsv_path} has no species columns "
                    f"(expected 'Orthogroup', then one column per species)."
                )
            self.species_in_og_table = header[1:]

            for row in reader:
                if not row:
                    continue
                og_id = row[0].strip()
                if not og_id:
                    continue

                species_genes = {}
                for i, species in enumerate(self.species_in_og_table):
                    col_idx = i + 1
                    if col_idx >= len(row):
                        continue
                    cell = row[col_idx].strip()
                    if not cell:
                        continue
                    genes = [g.strip() for g in cell.split(",") if g.strip()]
                    if genes:
                        species_genes[species] = set(genes)
                        for gene_id in genes:
                            self.gene_to_og.setdefault((species, gene_id), og_id)

                self.og_species_genes[og_id] = species_genes

        self._orthogroups_loaded = True

    def check_against_orthogroups(self, species, gene_id):
        """
        Returns a sorted list of other species sharing this gene's flat
        Orthogroup, or an empty list if the gene is unique to `species` (or
        not found at all -- absence here is not itself a conflict).
        """
        if not self.orthogroups_tsv_path:
            return []
        self.load_orthogroups_table()

        og_id = self.gene_to_og.get((species, gene_id))
        if og_id is None:
            return []

        species_genes = self.og_species_genes.get(og_id, {})
        return sorted(sp for sp, genes in species_genes.items() if sp != species and genes)

    def _load_pairwise_files_for_species(self, species):
        """
        Parses every Orthologues_<species>/<species>__v__<other>.tsv file
        for one species. Returns {gene_id: [(other_species, other_gene_id), ...]}.
        Cached per species since each species' pairwise files get consulted
        once per candidate gene of that species.
        """
        if species in self._pairwise_cache:
            return self._pairwise_cache[species]

        gene_to_partners = {}
        species_dir = self.orthologues_dir / f"Orthologues_{species}"

        if not species_dir.is_dir():
            self._pairwise_cache[species] = gene_to_partners
            return gene_to_partners

        for pairwise_file in sorted(species_dir.glob(f"{species}__v__*.tsv")):
            other_species = pairwise_file.stem[len(f"{species}__v__"):]

            with open(pairwise_file, "r", newline="", encoding="utf-8") as fh:
                reader = csv.reader(fh, delimiter="\t")
                header = next(reader, None)
                if not header:
                    continue

                for row in reader:
                    if len(row) < 3:
                        continue
                    focal_genes = [g.strip() for g in row[1].split(",") if g.strip()]
                    other_genes = [g.strip() for g in row[2].split(",") if g.strip()]
                    if not other_genes:
                        continue
                    for gene_id in focal_genes:
                        gene_to_partners.setdefault(gene_id, []).extend(
                            (other_species, other_gene_id) for other_gene_id in other_genes
                        )

        self._pairwise_cache[species] = gene_to_partners
        return gene_to_partners

    def check_against_pairwise_orthologues(self, species, gene_id):
        """
        Returns a list of (other_species, other_gene_id) tuples this gene
        is paired with in OrthoFinder's pairwise Orthologues output, or an
        empty list if none (or if --orthologues-dir wasn't provided).
        """
        if not self.orthologues_dir:
            return []
        gene_to_partners = self._load_pairwise_files_for_species(species)
        return gene_to_partners.get(gene_id, [])

    def filter_candidates(self, species, candidate_genes, original_status_lookup=None):
        """
        Given a species and its set of candidate cLSG gene IDs, returns the
        subset that survives both cross-checks, and records every exclusion
        (with reason) in self.exclusions.

        original_status_lookup: optional {gene_id: status_label} (e.g.
        "species_specific_clustered" / "species_specific_unclustered") used
        only to make the exclusion log more informative.
        """
        survivors = set()

        for gene_id in candidate_genes:
            original_status = (original_status_lookup or {}).get(gene_id, "unknown")

            og_conflicts = self.check_against_orthogroups(species, gene_id)
            if og_conflicts:
                self.exclusions.append({
                    "species": species,
                    "gene_id": gene_id,
                    "original_status": original_status,
                    "check": "orthogroups",
                    "conflicting_species": ";".join(og_conflicts),
                    "conflicting_gene": "",
                })
                continue

            pairwise_hits = self.check_against_pairwise_orthologues(species, gene_id)
            if pairwise_hits:
                self.exclusions.append({
                    "species": species,
                    "gene_id": gene_id,
                    "original_status": original_status,
                    "check": "pairwise_orthologues",
                    "conflicting_species": ";".join(sorted({sp for sp, _ in pairwise_hits})),
                    "conflicting_gene": ";".join(sorted({g for _, g in pairwise_hits})),
                })
                continue

            survivors.add(gene_id)

        return survivors

    def write_exclusion_log(self, path):
        """
        Writes the exclusion log (TSV). Always written, even if empty.
        Sorted by (species, gene_id) for reproducible, diffable output --
        self.exclusions is built while iterating candidate_genes, which is a
        set, so insertion order isn't stable across runs (Python's hash
        randomization) even though the actual set of exclusions is.
        """
        sorted_exclusions = sorted(self.exclusions, key=lambda e: (e["species"], e["gene_id"]))
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "species", "gene_id", "original_status", "check",
                    "conflicting_species", "conflicting_gene",
                ],
            )
            writer.writeheader()
            writer.writerows(sorted_exclusions)
