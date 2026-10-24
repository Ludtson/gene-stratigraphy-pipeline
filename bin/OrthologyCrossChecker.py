#!/usr/bin/env python3

"""
OrthologyCrossChecker.py

Cross-checks candidate species-specific genes against OrthoFinder's
pairwise Orthologues output (reciprocal-best-hit, tree-refined), and
reports which candidates have a conflicting ortholog there.

Why this exists: the primary classification (see
OrthogroupsStratigraphyClassifier.py) determines species-specificity from
Orthogroups.tsv's clustering. Clustering can under-merge -- miss a real
relationship -- due to sequence divergence, domain rearrangement, or
alignment noise, even though Orthogroups.tsv itself has no blind spot for
species coverage the way a Hierarchical Orthogroups table does. Pairwise
Orthologues is a genuinely different, independent method (reciprocal best
hit, not similarity clustering), so it can catch a real ortholog that
clustering missed.

(This class used to also check Orthogroups.tsv itself, back when the
primary classification came from a HOG table and Orthogroups.tsv was the
downstream safety net for the HOG table's outgroup blind spot. Now that
Orthogroups.tsv IS the primary source, checking it against itself would be
circular, so only the pairwise check remains.)

A candidate flagged here is excluded, and the reason is recorded so the
exclusion is auditable, not silent.

Author: Adekola Owoyemi (Casola Lab, Texas A&M University)
Version: 2.0.0
"""

import csv
import logging
import sys
from pathlib import Path

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

__version__ = "2.0.0"
__author__ = "Adekola Owoyemi (Casola Lab, Texas A&M University)"


class OrthologyCrossChecker:
    def __init__(self, orthologues_dir=None):
        """
        orthologues_dir: path to OrthoFinder's Orthologues/ directory,
            containing one Orthologues_<species>/ subdirectory per species,
            each with <species>__v__<other>.tsv pairwise files. If omitted,
            filter_candidates() returns every candidate unchanged (with a
            warning, not silently) -- useful for a rerun where the file
            genuinely isn't available.
        """
        self.orthologues_dir = Path(orthologues_dir) if orthologues_dir else None

        # Cache of parsed pairwise files, keyed by species, to avoid
        # re-reading the same file for every gene of that species.
        self._pairwise_cache = {}

        # Exclusion log entries: {gene_id, species, original_status,
        # check, conflicting_species, conflicting_gene}
        self.exclusions = []

        # Species whose candidates were checked / could not be checked
        # because Orthologues_<species>/ is absent.
        self.checked_species = set()
        self.unchecked_species = set()

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
        if not self.orthologues_dir:
            self._pairwise_cache[species] = gene_to_partners
            return gene_to_partners

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
        empty list if none (or if orthologues_dir wasn't provided).
        """
        if not self.orthologues_dir:
            return []
        gene_to_partners = self._load_pairwise_files_for_species(species)
        return gene_to_partners.get(gene_id, [])

    def filter_candidates(self, species, candidate_genes, original_status_lookup=None):
        """
        Given a species and its set of candidate species-specific gene IDs,
        returns the subset that survives the pairwise-Orthologues check,
        and records every exclusion (with reason) in self.exclusions.

        original_status_lookup: optional {gene_id: status_label} (e.g.
        "species_specific_clustered" / "species_specific_unclustered") used
        only to make the exclusion log more informative.
        """
        if not self.orthologues_dir:
            return set(candidate_genes)

        if (self.orthologues_dir / f"Orthologues_{species}").is_dir():
            self.checked_species.add(species)
        else:
            self.unchecked_species.add(species)
            logging.warning(
                f"No Orthologues_{species}/ under {self.orthologues_dir}: "
                f"pairwise cross-check skipped for {species} "
                f"({len(candidate_genes)} candidate gene(s) left unverified)."
            )

        survivors = set()

        for gene_id in candidate_genes:
            original_status = (original_status_lookup or {}).get(gene_id, "unknown")

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
                fh, delimiter="\t", lineterminator="\n",
                fieldnames=[
                    "species", "gene_id", "original_status", "check",
                    "conflicting_species", "conflicting_gene",
                ],
            )
            writer.writeheader()
            writer.writerows(sorted_exclusions)
