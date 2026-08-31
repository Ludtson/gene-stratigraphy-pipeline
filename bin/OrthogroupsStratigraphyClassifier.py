#!/usr/bin/env python3

"""
OrthogroupsStratigraphyClassifier.py

Classifies each gene of a focal species by phylogenetic depth of origin,
using OrthoFinder's flat Orthogroups.tsv (whole-dataset clustering) plus
the rooted species tree, instead of a Hierarchical Orthogroups (HOG) table.

Why: HOG tables (N0.tsv, N1.tsv, ...) are structurally blind to any species
outside the clade covered by that specific node. If an outgroup is used to
work around OrthoFinder not emitting N0.tsv (see the main pipeline's
README), the outgroup -- and, depending on which node file is used, any
other species outside the analyzed clade -- is invisible to the HOG-based
method entirely, for every focal species, not just the ones directly
compared to it. That's not a rare edge case here: it silently forces every
gene's stratigraphic ceiling to whatever the shallowest available HOG node
is, understating true conservation for genes that would also be shared
with the invisible species if the tool could see it.

Orthogroups.tsv has no such blind spot: it is a single, flat clustering
that always includes every species OrthoFinder was run on. It doesn't
carry per-node depth information directly (unlike a HOG table, it has no
"which tree node does this orthogroup correspond to" label) -- but it
doesn't need to: given the set of species sharing a gene's orthogroup, the
correct stratigraphic node is the most recent common ancestor (MRCA) of
that species set on the actual species tree, computed directly from tree
topology (see SpeciesTree.get_mrca in OrthoFinderTreeParser.py).

Three possible outcomes per gene:
  - species_specific_clustered:   in an orthogroup, but private to this species
  - species_specific_unclustered: absent from Orthogroups.tsv entirely
  - shared: in an orthogroup with other species; node = MRCA of
    (this species + every other species in that orthogroup)

Author: Adekola Owoyemi (Casola Lab, Texas A&M University)
Version: 1.0.0
"""

import csv
import sys
from pathlib import Path

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

__version__ = "1.0.0"
__author__ = "Adekola Owoyemi (Casola Lab, Texas A&M University)"

SPECIES_SPECIFIC_CLUSTERED = "species_specific_clustered"
SPECIES_SPECIFIC_UNCLUSTERED = "species_specific_unclustered"


class OrthogroupsStratigraphyClassifier:
    def __init__(self, orthogroups_tsv_path, species_tree):
        """
        orthogroups_tsv_path: path to OrthoFinder's Orthogroups.tsv.
        species_tree: an already-constructed SpeciesTree instance (see
            OrthoFinderTreeParser.py), used to compute MRCA nodes.
        """
        self.orthogroups_tsv_path = Path(orthogroups_tsv_path)
        self.species_tree = species_tree

        # {og_id: {species: set(genes)}}
        self.og_species_genes = {}
        # {(species, gene_id): og_id}
        self.gene_to_og = {}
        self.species_in_og_table = []
        self._loaded = False

    def parse(self):
        """Parses Orthogroups.tsv, populating og_species_genes and gene_to_og."""
        if self._loaded:
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

        self._loaded = True

    def classify_gene(self, species, gene_id):
        """
        Returns (status, node_label):
          - (SPECIES_SPECIFIC_UNCLUSTERED, SPECIES_SPECIFIC_UNCLUSTERED) if
            the gene is absent from Orthogroups.tsv entirely.
          - (SPECIES_SPECIFIC_CLUSTERED, SPECIES_SPECIFIC_CLUSTERED) if the
            gene's orthogroup is private to this species.
          - ("shared", mrca_node_name) otherwise, where mrca_node_name is
            the MRCA (on the species tree) of this species plus every other
            species present in the gene's orthogroup.
        """
        og_id = self.gene_to_og.get((species, gene_id))
        if og_id is None:
            return (SPECIES_SPECIFIC_UNCLUSTERED, SPECIES_SPECIFIC_UNCLUSTERED)

        species_genes = self.og_species_genes.get(og_id, {})
        other_species = sorted(sp for sp, genes in species_genes.items() if sp != species and genes)
        if not other_species:
            return (SPECIES_SPECIFIC_CLUSTERED, SPECIES_SPECIFIC_CLUSTERED)

        try:
            mrca = self.species_tree.get_mrca([species] + other_species)
        except ValueError as e:
            # get_mrca fails only if a species name in this orthogroup isn't
            # a leaf of the tree at all (e.g. a name mismatch between
            # Orthogroups.tsv and the species tree file) -- shouldn't happen
            # since both come from the same OrthoFinder run, but if it does,
            # this gene genuinely cannot be placed at a specific node. The
            # conservative default is the root (N0), not silently landing on
            # whatever node happened to be the deepest one still resolvable
            # -- an unplaceable gene should look maximally conserved/
            # unresolved, not get miscategorized as narrowly shared.
            root_name = self.species_tree.root.name.split(':')[0]
            print(
                f"[WARNING] Could not compute MRCA for {species}/{gene_id} "
                f"(orthogroup {og_id}, other species {other_species}): {e} "
                f"-- defaulting to root ({root_name})."
            )
            return ("shared", root_name)
        return ("shared", mrca)
