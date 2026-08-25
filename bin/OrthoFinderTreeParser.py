#!/usr/bin/env python3

"""
OrthoFinderTreeParser.py

This script parses a phylogenetic tree file in Newick format (e.g., OrthoFinder's SpeciesTree.txt).
It includes the necessary traversal and utility methods (like find_nodes_up_tree) required by the 
main OrthoFinderDataProcessor.py pipeline script.

Required Arguments:
    --tree-file, -t: Path to the input Newick tree file (e.g., SpeciesTree.txt).

Optional Arguments:
    --output, -o: Path to an output file to save the Node -> Species mapping (TSV format).
    --species-map-file, -m: Path to an optional species short/long name mapping file (TSV format).
"""

import os
import argparse
import csv
from collections import OrderedDict
import sys

# Set a larger field size limit for potentially large gene names. sys.maxsize
# overflows the C `long` used internally by csv.field_size_limit on
# platforms where C long is 32-bit (e.g. Windows, even on 64-bit Python), so
# cap at the largest value that's safe everywhere.
csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

# --- Core Data Structures ---

class TreeNode:
    def __init__(self, name):
        self.name = name
        self.children = []
        self.parent = None
        self.branch_length = None

class SpeciesTree:
    def __init__(self, newick=None, file_path=None, species_map_file=None):
        self.root = None
        self._nodes_cache = {} 
        self.species_map = {} # {short_name: long_name}

        if species_map_file:
            self.species_map = self._load_species_map(species_map_file)
            
        if newick:
            self.root = self.parse_newick(newick)
        elif file_path:
            self.root = self.read_tree_from_file(file_path)
        
        # Populate the cache after parsing for easy lookup
        if self.root:
            self._populate_node_cache(self.root)

    def _load_species_map(self, map_file):
        """Loads the species map file (short_name -> long_name)."""
        mapping = {}
        try:
            with open(map_file, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    if 'short_name' in row and 'Species_name' in row:
                        mapping[row['short_name'].strip()] = row['Species_name'].strip()
        except Exception as e:
            print(f"Warning: Could not load species map file {map_file}. Proceeding without mapping. Error: {e}")
            return {}
        return mapping
    
    def _populate_node_cache(self, node):
        """Helper to recursively populate the node cache."""
        clean_name = node.name.split(':')[0]
        if clean_name:
            self._nodes_cache[clean_name] = node
        
        for child in node.children:
            self._populate_node_cache(child)

    def parse_newick(self, newick):
        """
        Parse the Newick formatted string and construct the tree.
        (Original character iteration logic preserved)
        """
        stack = []
        root = None
        current_node = None
        name = ''
        reading_name = False
        reading_branch_length = False
        branch_length_str = ''
        i = 0
        while i < len(newick) and newick[i] != '(':
            name += newick[i]
            i += 1
        root = TreeNode(name.strip())
        current_node = root
        
        for char in newick[i:]:
            if char == '(':
                if reading_name:
                    current_node.name = name.strip()
                    name = ''
                    reading_name = False
                new_node = TreeNode('')
                stack.append(current_node)
                current_node.children.append(new_node)
                new_node.parent = current_node
                current_node = new_node
            elif char == ')':
                if reading_name:
                    current_node.name = name.strip()
                    name = ''
                    reading_name = False
                if reading_branch_length and branch_length_str:
                    try:
                        current_node.branch_length = float(branch_length_str)
                    except ValueError:
                         pass # Ignore if branch length is not a valid float
                    branch_length_str = ''
                    reading_branch_length = False
                if stack:
                    current_node = stack.pop()
                else:
                    # Reached the end of the main root's children, the next token might be a root label
                    pass 
            elif char == ',':
                if reading_name:
                    current_node.name = name.strip()
                    name = ''
                    reading_name = False
                if reading_branch_length and branch_length_str:
                    try:
                        current_node.branch_length = float(branch_length_str)
                    except ValueError:
                         pass # Ignore if branch length is not a valid float
                    branch_length_str = ''
                    reading_branch_length = False
                if stack: # Should only create a sibling if we are not at the root level
                    new_node = TreeNode('')
                    stack[-1].children.append(new_node)
                    new_node.parent = stack[-1]
                    current_node = new_node
            elif char == ';':
                if reading_name:
                    current_node.name = name.strip()
                    name = ''
                    reading_name = False
                if reading_branch_length and branch_length_str:
                    try:
                        current_node.branch_length = float(branch_length_str)
                    except ValueError:
                         pass # Ignore if branch length is not a valid float
                    branch_length_str = ''
                    reading_branch_length = False
                break
            elif char == ':':
                if reading_name:
                    current_node.name = name.strip()
                    name = ''
                    reading_name = False
                reading_branch_length = True
            elif reading_branch_length:
                branch_length_str += char
            else:
                name += char
                reading_name = True
        
        # After loop, if stack is empty, root is the final node
        return root

    def read_tree_from_file(self, file_path):
        """Read a Newick formatted tree from a file and parse it."""
        try:
            with open(file_path, 'r') as file:
                newick_tree = file.read().strip()
            return self.parse_newick(newick_tree)
        except FileNotFoundError:
            raise FileNotFoundError(f"Tree file not found at {file_path}.")

    def _find_node(self, node_name):
        """Helper function to find a node by name using the cache."""
        return self._nodes_cache.get(node_name)

    # --- REQUIRED UTILITY METHODS ---

    def _get_leaf_species_names(self, node):
        """
        Helper: Recursively collects the names (IDs) of all species (leaves) under a node.
        """
        species_names = []
        if not node.children:
            # It's a leaf node, return its clean name
            return [node.name.split(':')[0]]
        
        for child in node.children:
            species_names.extend(self._get_leaf_species_names(child))
            
        return species_names

    def get_all_species(self):
        """Returns a set of all species names (leaf IDs) in the tree."""
        if not self.root:
            return set()
        return set(self._get_leaf_species_names(self.root))

    def get_parent(self, node_name):
        """
        Returns the name of the parent node of a given node.
        """
        node = self._find_node(node_name)

        if not node:
            raise ValueError(f"Node {node_name} not found.")

        if node.parent:
            return node.parent.name.split(':')[0]  # Return parent's clean name
        else:
            return None  # Node is the root

    def is_ancestor(self, ancestor_name, descendant_name):
        """
        Checks if ancestor_name is an ancestor of descendant_name.
        """
        descendant_node = self._find_node(descendant_name)
        
        if not descendant_node:
            raise ValueError(f"Node {descendant_name} not found.")
            
        current = descendant_node.parent
        while current:
            if current.name.split(':')[0] == ancestor_name:
                return True
            current = current.parent
        return False
        
    def find_nodes_up_tree(self, focal_node_name, is_species=False):
        """
        REQUIRED METHOD
        Finds all nodes up the tree from the focal node (excluding the focal node).
        
        Args:
            focal_node_name (str): The name of the node (species ID or N-label) to start from.
            is_species (bool): Ignored in this implementation, kept for compatibility.
            
        Returns:
            list: A list of clean node names (IDs) moving up towards the root.
        """
        node = self._find_node(focal_node_name)
        if not node:
            raise ValueError(f"Focal node {focal_node_name} not found.")
            
        ancestors = []
        current = node.parent
        while current:
            clean_name = current.name.split(':')[0]
            if clean_name:
                ancestors.append(clean_name)
            current = current.parent
            
        return ancestors

    def find_sister_species_to_node(self, focal_node_name):
        """
        REQUIRED METHOD
        Finds all species (leaves) that are descendants of the sister node(s) 
        to the focal node.
        
        Args:
            focal_node_name (str): The name of the node (species ID or N-label).
            
        Returns:
            set: A set of clean species names (leaf IDs).
        """
        node = self._find_node(focal_node_name)
        if not node:
            raise ValueError(f"Focal node {focal_node_name} not found.")

        if not node.parent:
            # Focal node is the root, it has no sister
            return set()

        sister_species = set()
        parent = node.parent
        
        for sibling in parent.children:
            if sibling != node:
                # Collect all leaf species under the sister node
                sister_species.update(self._get_leaf_species_names(sibling))
                
        return sister_species
    
    def nodes_n_species(self):
        """
        Returns a dictionary mapping each node name to the list of species (leaf names) under that node.
        Translates short species names to long names if a map is loaded.
        """
        # 1. First, collect data using the short names (node IDs, species basenames)
        node_species_map = OrderedDict() # Stores {short_node_name: [short_species_list]}

        def collect_species(node):
            if not node.children:
                return [node.name.split(':')[0]]
            
            species = []
            for child in node.children:
                species.extend(collect_species(child))
                
            clean_name = node.name.split(':')[0]
            # Only map internal nodes that have a label (N0, N1, etc.)
            if clean_name and clean_name.startswith('N'):
                node_species_map[clean_name] = species
            
            # Also map leaf nodes if they were not already added by collect_species
            elif not node.children:
                 node_species_map[clean_name] = species
                 
            return species

        if self.root:
            collect_species(self.root)
        
        # 2. If no mapping is loaded, return the result using short names
        if not self.species_map:
            return node_species_map

        # 3. Apply species name mapping (short -> long)
        mapped_species_map = OrderedDict()
        for node_id, short_names in node_species_map.items():
            
            # Map the contained species names (leaves)
            long_names = [self.species_map.get(short, short) for short in short_names]
            
            # Map the node ID if it is a species name (leaf), otherwise keep the ID (N0, N1, etc.)
            mapped_node_id = self.species_map.get(node_id, node_id)
            
            mapped_species_map[mapped_node_id] = long_names

        return mapped_species_map
    

# --- CLI and Main Function ---

def main():
    parser = argparse.ArgumentParser(description='Parse an OrthoFinder Newick species tree and map internal nodes to contained species.')
    
    parser.add_argument("--tree-file", "-t", type=str, required=True, help="Path to the input Newick tree file (e.g., SpeciesTree.txt).")
    parser.add_argument("--output", "-o", type=str, default=None, help="Path to an output file to save the Node -> Species mapping (TSV format). If not provided, prints to console.")
    parser.add_argument("--species-map-file", "-m", type=str, default=None, help="Path to an optional species short/long name mapping file (TSV format).")
    
    args = parser.parse_args()
    
    print(f"Loading tree from: {args.tree_file}")
    try:
        tree = SpeciesTree(file_path=args.tree_file, species_map_file=args.species_map_file)
    except Exception as e:
        print(f"Failed to parse tree: {e}")
        return

    if not tree.root:
        print("Failed to parse the tree. Aborting.")
        return

    # Use the modified method which applies the map if available
    species_map = tree.nodes_n_species()
    
    if not species_map:
        print("No nodes or species found in the tree.")
        return

    # Output generation logic for file or console
    header = ['Node_ID', 'Contained_Species']
    
    # Writing or printing the mapping
    if args.output:
        try:
            with open(args.output, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(header)
                for node_id, species_list in species_map.items():
                    writer.writerow([node_id, ', '.join(species_list)])
            print(f"Successfully generated Node -> Species map and saved to: {args.output}")
        except Exception as e:
            print(f"Error writing to output file {args.output}: {e}")
    else:
        print("\n--- Node -> Species Mapping (Partial Output) ---")
        print(f"{header[0]}\t{header[1]}")
        print("-" * 40)
        
        # Print a few example nodes to demonstrate functionality
        example_keys = list(species_map.keys())
        
        # Collect nodes: first 3, and any labeled N0, N1, N6, N14 for demonstration
        demo_nodes = list(example_keys[:3])
        for node_name in ['N0', 'N1', 'N6', 'N14']:
            if node_name in example_keys and node_name not in demo_nodes:
                demo_nodes.append(node_name)
        
        for node_id in demo_nodes:
            species_list = species_map[node_id]
            display_list = species_list if len(species_list) < 5 else species_list[:4] + [f'...({len(species_list)} total)']
            print(f"{node_id}\t{', '.join(display_list)}")

if __name__ == "__main__":
    main()

