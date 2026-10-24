#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SpeciesTreeBranchLengths.py

Adds substitution-based branch lengths to a fixed rooted species tree, for
tools that need a rooted tree with branch lengths (e.g. whole-genome
aligners). The topology and the root are kept exactly as given.

Steps:
  1. Select orthogroups from Orthogroups.tsv with at most one gene per
     species, present in at least --min-species-frac of the species
     (single-copy orthogroups when 1.0). Sorted by orthogroup ID; the first
     --max-orthogroups are used.
  2. Align each orthogroup's proteins with MAFFT.
  3. Concatenate the alignments per species (missing species = gaps) and
     drop columns with more than --max-gap-frac gaps.
  4. FastTree -nome -mllen -intree: maximum-likelihood branch lengths on the
     fixed topology.
  5. FastTree writes an unrooted tree; it is re-rooted on the edge that
     matches the input root, and its topology is checked against the input.

Sequences are taken from --protein-dir (<species>.fa/.faa/.fasta/.fas, the
species name as in Orthogroups.tsv); FASTA IDs are the first word of the
header. MAFFT and FastTree must be on PATH (both are OrthoFinder
dependencies).

Usage:
  python SpeciesTreeBranchLengths.py \\
    --orthogroups Orthogroups.tsv \\
    --protein-dir protein \\
    --tree SpeciesTree_rooted.txt \\
    --out SpeciesTree_rooted_with_lengths.nwk \\
    --threads 8

Author: Adekola Owoyemi (Casola Lab, Texas A&M University)
Version: 1.0.0
"""

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

__version__ = "1.0.0"
__author__ = "Adekola Owoyemi (Casola Lab, Texas A&M University)"

FASTA_EXTS = (".fa", ".faa", ".fasta", ".fas")


# ------------------------------------------------------------------ Newick
class Node:
    __slots__ = ("name", "length", "children")

    def __init__(self, name="", length=None):
        self.name = name
        self.length = length
        self.children = []

    def leaves(self):
        if not self.children:
            return [self.name]
        out = []
        for c in self.children:
            out.extend(c.leaves())
        return out


def parse_newick(text):
    """Parses Newick into Node objects. Internal labels and support values are kept as names."""
    s = text.strip()
    if not s.endswith(";"):
        raise ValueError("Newick string does not end with ';'")
    pos = 0

    def read_label():
        nonlocal pos
        if pos < len(s) and s[pos] == "'":
            end = s.index("'", pos + 1)
            label = s[pos + 1:end]
            pos = end + 1
            return label
        start = pos
        while pos < len(s) and s[pos] not in ",():;":
            pos += 1
        return s[start:pos].strip()

    def read_length():
        nonlocal pos
        if pos < len(s) and s[pos] == ":":
            pos += 1
            start = pos
            while pos < len(s) and s[pos] not in ",();":
                pos += 1
            return float(s[start:pos])
        return None

    def read_node():
        nonlocal pos
        node = Node()
        if s[pos] == "(":
            pos += 1
            while True:
                node.children.append(read_node())
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"Unexpected character {s[pos]!r} at position {pos}")
        node.name = read_label()
        node.length = read_length()
        return node

    root = read_node()
    if s[pos] != ";":
        raise ValueError(f"Unexpected text after tree at position {pos}")
    return root


def to_newick(node, lengths=True, is_root=True):
    """Writes leaf names and (optionally) branch lengths; internal labels are not written."""
    if node.children:
        out = "(" + ",".join(to_newick(c, lengths, False) for c in node.children) + ")"
    else:
        out = node.name
    if lengths and not is_root and node.length is not None:
        out += f":{node.length:.6g}"
    return out + (";" if is_root else "")


def splits(root):
    """Non-trivial bipartitions of the unrooted tree, each as the side without the first leaf (sorted)."""
    all_leaves = frozenset(root.leaves())
    ref = sorted(all_leaves)[0]
    out = set()

    def walk(n):
        if not n.children:
            return frozenset([n.name])
        ls = frozenset().union(*(walk(c) for c in n.children))
        if 1 < len(ls) < len(all_leaves) - 1:
            out.add(ls if ref not in ls else all_leaves - ls)
        return ls

    walk(root)
    return out


def reroot_like(unrooted, rooted):
    """Re-roots `unrooted` on the edge that separates the two sides of `rooted`'s root. Edge length is split in half."""
    if len(rooted.children) != 2:
        raise ValueError("Input species tree must be rooted (bifurcating root).")
    side = frozenset(rooted.children[0].leaves())
    all_leaves = frozenset(unrooted.leaves())
    if all_leaves != frozenset(rooted.leaves()):
        raise ValueError("Leaf sets of the input tree and the FastTree output differ.")

    # undirected adjacency: node id -> list of (neighbour id, length)
    nodes, adj = {}, {}

    def build(n, parent_id):
        i = len(nodes)
        nodes[i] = n
        adj[i] = []
        if parent_id is not None:
            adj[i].append((parent_id, n.length or 0.0))
            adj[parent_id].append((i, n.length or 0.0))
        for c in n.children:
            build(c, i)
        return i

    build(unrooted, None)

    def leafset(i, came_from):
        if not nodes[i].children and len(adj[i]) <= 1:
            return frozenset([nodes[i].name])
        return frozenset().union(*(leafset(j, i) for j, _ in adj[i] if j != came_from))

    target = None
    for i in adj:
        for j, length in adj[i]:
            if i < j:
                ls = leafset(j, i)
                if ls == side or ls == all_leaves - side:
                    target = (i, j, length)
                    break
        if target:
            break
    if target is None:
        raise ValueError("The input root edge does not exist in the FastTree output (topology changed).")

    i, j, length = target

    def grow(k, came_from, edge_len):
        n = Node(nodes[k].name if not nodes[k].children or len(adj[k]) <= 1 else "", edge_len)
        kids = [(m, l) for m, l in adj[k] if m != came_from]
        for m, l in kids:
            n.children.append(grow(m, k, l))
        if len(n.children) == 1:  # suppress degree-2 node (old bifurcating root)
            only = n.children[0]
            only.length = (only.length or 0.0) + (edge_len or 0.0)
            return only
        return n

    root = Node()
    root.children = [grow(i, j, length / 2), grow(j, i, length / 2)]
    return root


# ------------------------------------------------------------------ data
def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0] if line[1:].split() else ""
                buf = []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def select_orthogroups(path, species, min_frac, max_ogs):
    selected = []
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        cols = header[1:]
        missing = sorted(set(species) - set(cols))
        if missing:
            raise ValueError(f"Tree species not in Orthogroups.tsv: {', '.join(missing)}")
        idx = {sp: cols.index(sp) + 1 for sp in species}
        need = max(2, math.ceil(min_frac * len(species) - 1e-9))
        for row in reader:
            if not row:
                continue
            genes = {}
            ok = True
            for sp, k in idx.items():
                cell = row[k].strip() if k < len(row) else ""
                ids = [g.strip() for g in cell.split(",") if g.strip()]
                if len(ids) > 1:
                    ok = False
                    break
                if ids:
                    genes[sp] = ids[0]
            if ok and len(genes) >= need:
                selected.append((row[0], genes))
    selected.sort(key=lambda x: x[0])
    return selected[:max_ogs], need


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed (exit {r.returncode}): {r.stderr.strip()[:500]}")
    return r.stdout


def main():
    ap = argparse.ArgumentParser(description="Add ML branch lengths to a fixed rooted species tree.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--orthogroups", required=True, help="OrthoFinder Orthogroups.tsv")
    ap.add_argument("--protein-dir", required=True, help="Directory with <species>.faa files (names as in Orthogroups.tsv)")
    ap.add_argument("--tree", required=True, help="Rooted species tree (Newick); topology and root are kept")
    ap.add_argument("--out", required=True, help="Output rooted tree with branch lengths (Newick)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-species-frac", type=float, default=1.0,
                    help="Orthogroup must contain this fraction of species (at most one gene each)")
    ap.add_argument("--max-orthogroups", type=int, default=1000)
    ap.add_argument("--min-orthogroups", type=int, default=50, help="Fail if fewer orthogroups qualify")
    ap.add_argument("--max-gap-frac", type=float, default=0.5, help="Drop alignment columns with more gaps than this")
    ap.add_argument("--keep-work", help="Keep alignments in this directory")
    ap.add_argument("--version", action="version", version=__version__)
    a = ap.parse_args()

    for tool in ("mafft", "FastTree"):
        if not shutil.which(tool):
            sys.exit(f"[ERROR] {tool} not found on PATH")

    rooted = parse_newick(Path(a.tree).read_text(encoding="utf-8"))
    species = sorted(rooted.leaves())
    if len(set(species)) != len(species):
        sys.exit("[ERROR] Duplicate leaf names in the species tree")

    ogs, need = select_orthogroups(a.orthogroups, species, a.min_species_frac, a.max_orthogroups)
    print(f"[INFO] {len(ogs)} orthogroups with <= 1 gene per species in >= {need} of {len(species)} species")
    if len(ogs) < a.min_orthogroups:
        sys.exit(f"[ERROR] Only {len(ogs)} orthogroups qualify (< --min-orthogroups {a.min_orthogroups}); "
                 f"lower --min-species-frac")

    pdir = Path(a.protein_dir)
    seqs = {}
    for sp in species:
        path = next((pdir / f"{sp}{e}" for e in FASTA_EXTS if (pdir / f"{sp}{e}").is_file()), None)
        if path is None:
            sys.exit(f"[ERROR] No protein FASTA for {sp} in {pdir}")
        seqs[sp] = read_fasta(path)

    work = Path(a.keep_work) if a.keep_work else Path(tempfile.mkdtemp(prefix="sptree_bl_"))
    work.mkdir(parents=True, exist_ok=True)
    try:
        concat = {sp: [] for sp in species}
        for n, (og, genes) in enumerate(ogs, 1):
            raw = work / f"{og}.fa"
            with open(raw, "w", encoding="utf-8", newline="\n") as fh:
                for sp, g in sorted(genes.items()):
                    if g not in seqs[sp]:
                        sys.exit(f"[ERROR] {g} ({og}) not found in the {sp} protein FASTA")
                    fh.write(f">{sp}\n{seqs[sp][g].rstrip('*')}\n")
            aln_text = run(["mafft", "--auto", "--quiet", "--thread", str(a.threads), str(raw)])
            aln = {}
            name = None
            for line in aln_text.splitlines():
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    aln[name] = []
                elif name:
                    aln[name].append(line.strip())
            aln = {k: "".join(v).upper() for k, v in aln.items()}
            width = len(next(iter(aln.values())))
            for sp in species:
                concat[sp].append(aln.get(sp, "-" * width))
            if n % 100 == 0:
                print(f"[INFO] aligned {n}/{len(ogs)}", flush=True)

        full = {sp: "".join(v) for sp, v in concat.items()}
        ncol = len(full[species[0]])
        keep = [c for c in range(ncol) if sum(full[sp][c] in "-X?" for sp in species) / len(species) <= a.max_gap_frac]
        sup = work / "supermatrix.fa"
        with open(sup, "w", encoding="utf-8", newline="\n") as fh:
            for sp in species:
                fh.write(f">{sp}\n{''.join(full[sp][c] for c in keep)}\n")
        print(f"[INFO] supermatrix: {len(species)} species x {len(keep)} columns (of {ncol})")

        topo = work / "topology.nwk"
        topo.write_text(to_newick(rooted, lengths=False) + "\n", encoding="utf-8")
        env = dict(os.environ, OMP_NUM_THREADS=str(a.threads))
        ft = run(["FastTree", "-quiet", "-nome", "-mllen", "-intree", str(topo), str(sup)], env=env)
        fitted = parse_newick(ft.strip())

        if splits(fitted) != splits(rooted):
            sys.exit("[ERROR] FastTree output topology differs from the input tree")
        result = reroot_like(fitted, rooted)
        if splits(result) != splits(rooted) or frozenset(result.children[0].leaves()) not in (
                frozenset(rooted.children[0].leaves()), frozenset(rooted.children[1].leaves())):
            sys.exit("[ERROR] Re-rooting did not reproduce the input root")
        Path(a.out).write_text(to_newick(result) + "\n", encoding="utf-8")
        print(f"[INFO] written: {a.out}")
    finally:
        if not a.keep_work:
            shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
