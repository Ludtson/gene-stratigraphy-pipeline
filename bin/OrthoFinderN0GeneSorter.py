#!/usr/bin/env python3

"""
OrthoFinderN0GeneSorter.py

This script processes the OrthoFinder N0.tsv file to create a species-centric data
structure, enabling quick retrieval and classification of Hierarchical Orthogroups (HOGs) 
and the genes they contain (e.g., species-specific HOGs).

This version has been updated to include methods required by the main orchestrator
script (OrthoFinderDataProcessor.py), specifically:
1. All core classification methods now accept the 'is_large_file' argument.
2. It includes 'generate_summary_data()' and 'write_summary_csv()'.
"""

import argparse
import csv
from collections import defaultdict
import os
import sys

# Increase field size limit for large bioinformatics files. sys.maxsize
# overflows the C `long` used internally by csv.field_size_limit on
# platforms where C long is 32-bit (e.g. Windows, even on 64-bit Python), so
# cap at the largest value that's safe everywhere.
csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

class N0OrthoGeneSorter:
    def __init__(self, file_path, species_map_file=None):
        """
        Initialize the N0OrthoGeneSorter with the file path and parse the TSV file.
        """
        self.file_path = file_path
        # {short_name: long_name}
        self.species_map = self._load_species_map(species_map_file) if species_map_file else {}
        # data[species][hog_id] = {set of genes}
        self.data = self.parse_tsv() 
        # Sort keys to ensure consistent iteration order
        self.all_species = sorted(list(self.data.keys()))

    def _load_species_map(self, map_file):
        """Loads the species map file (short_name -> long_name)."""
        mapping = {}
        try:
            with open(map_file, 'r') as f:
                # Use DictReader for files with a header, assuming two columns: 'short_name' and 'Species_name'
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    # Robustly check for the required keys, using .get() with a default empty string
                    short_name = row.get('short_name', '').strip()
                    long_name = row.get('Species_name', '').strip()
                    if short_name and long_name:
                        mapping[short_name] = long_name
        except Exception as e:
            print(f"Warning: Could not load species map file {map_file}. Proceeding without mapping. Error: {e}")
            return {}
        return mapping

    def get_long_species_name(self, short_name):
        """Translates short name to long name if map is available."""
        return self.species_map.get(short_name, short_name)

    def parse_tsv(self):
        """
        Parse the N0.tsv file and return a nested dictionary.
        
        Returns:
            dict: {species: {hog_id: set(genes)}}
        """
        data = defaultdict(lambda: defaultdict(set))
        
        try:
            with open(self.file_path, 'r', newline='') as file:
                reader = csv.reader(file, delimiter='\t')
                
                try:
                    header = next(reader)
                except StopIteration:
                    print(f"Error: N0.tsv file is empty.")
                    return {}

                # Identify the columns corresponding to species
                species_list = header[3:]
                species_columns_map = {species: i + 3 for i, species in enumerate(species_list)}
                
                if not species_columns_map:
                    print("Error: N0.tsv file appears to be misformatted (no species columns found).")
                    return {}

                for row in reader:
                    if not row or len(row) < 3: continue 
                    hog_id = row[0].strip()
                    
                    for species, col_idx in species_columns_map.items():
                        if col_idx < len(row):
                            genes_str = row[col_idx].strip()
                            if genes_str and genes_str != '-': # Check for non-empty/non-placeholder entries
                                genes = genes_str.split(', ') 
                                data[species][hog_id].update(genes)
                
            print(f"Successfully parsed N0.tsv: Loaded data for {len(data)} species.")
            return dict(data)
            
        except FileNotFoundError:
            print(f"Error: Input TSV file not found at {self.file_path}. Exiting.")
            sys.exit(1)
        except Exception as e:
            print(f"An error occurred while parsing {self.file_path}: {e}")
            sys.exit(1)

    # --- HOG Classification Methods (Updated with is_large_file) ---

    def get_unique_hog_ids(self, focal_species, is_large_file=False):
        """
        Get HOG IDs unique to a species. Signature updated to match orchestrator.
        """
        if focal_species not in self.data: return set()
        unique_hog_ids = set(self.data[focal_species].keys())
        for other_species in self.all_species:
            if other_species != focal_species:
                unique_hog_ids -= set(self.data[other_species].keys())
        return unique_hog_ids

    def get_non_unique_hog_ids(self, focal_species, is_large_file=False):
        """
        Get HOG IDs that are not unique to a species. Signature updated to match orchestrator.
        """
        if focal_species not in self.data: return set()
        all_hog_ids = set(self.data[focal_species].keys())
        # Pass is_large_file through to the helper method
        unique_hog_ids = self.get_unique_hog_ids(focal_species, is_large_file)
        return all_hog_ids - unique_hog_ids

    def get_shared_hog_ids(self, focal_species, other_species_list, is_large_file=False):
        """
        Get HOG IDs shared by the focal species and at least one species in the other_species_list.
        Signature updated to match orchestrator.
        """
        if focal_species not in self.data: return set()
        focal_hogs = set(self.data[focal_species].keys())
        other_hogs = set()
        
        for species in other_species_list:
            if species in self.data:
                other_hogs.update(self.data[species].keys())
                
        # Intersection: HOGs present in focal species AND present in any of the other species
        return focal_hogs.intersection(other_hogs)

    def get_genes_for_hog_ids(self, species, hog_ids, is_large_file=False):
        """
        Get genes for the specified HOG IDs in a species. Signature updated to match orchestrator.
        """
        if species not in self.data: return set()
        genes = set()
        species_data = self.data[species]
        for hog_id in hog_ids:
            if hog_id in species_data:
                genes.update(species_data[hog_id])
        return genes

    # --- New Summary Generation Methods (Updated with is_large_file) ---

    def generate_summary_data(self, is_large_file=False):
        """
        Calculates all required HOG and gene metrics for every species.
        
        Returns:
            list of dicts: List of summary rows, one for each species.
        """
        summary_list = []
        
        for species_short in self.all_species:
            species_long = self.get_long_species_name(species_short)
            
            # 1. Total HOGs
            total_hogs = len(self.data[species_short])
            
            # 2 & 3. Unique HOGs and Genes
            # Pass is_large_file through
            unique_hogs = self.get_unique_hog_ids(species_short, is_large_file=is_large_file)
            unique_hogs_count = len(unique_hogs)
            # Pass is_large_file through
            genes_in_unique_hogs_count = len(self.get_genes_for_hog_ids(species_short, unique_hogs, is_large_file=is_large_file))
            
            # 4 & 5. Non-Unique HOGs and Genes
            # Pass is_large_file through
            non_unique_hogs = self.get_non_unique_hog_ids(species_short, is_large_file=is_large_file)
            non_unique_hogs_count = len(non_unique_hogs)
            # Pass is_large_file through
            genes_in_non_unique_hogs_count = len(self.get_genes_for_hog_ids(species_short, non_unique_hogs, is_large_file=is_large_file))

            summary_list.append({
                'Species': species_long,
                'Total_HOGs': total_hogs,
                'Unique_HOGs': unique_hogs_count,
                'Genes_in_Unique_HOGs': genes_in_unique_hogs_count,
                'Non_Unique_HOGs': non_unique_hogs_count,
                'Genes_in_Non_Unique_HOGs': genes_in_non_unique_hogs_count,
            })
            
        return summary_list

    def write_summary_csv(self, output_path, summary_data):
        """Writes the generated summary data to a CSV file."""
        headers = [
            'Species', 
            'Total_HOGs', 
            'Unique_HOGs', 
            'Genes_in_Unique_HOGs', 
            'Non_Unique_HOGs', 
            'Genes_in_Non_Unique_HOGs'
        ]

        try:
            with open(output_path, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=headers)
                writer.writeheader()
                
                # Sort the data by species name for consistent output
                # We need to sort by the *short* name used as keys in summary_data if possible
                
                # Since summary_data is a list of dicts, we sort the list based on the 'Species' key
                sorted_data = sorted(summary_data, key=lambda x: x['Species'])

                writer.writerows(sorted_data)
            print(f"Successfully generated HOG Summary CSV for all species and saved to: {output_path}")
        except Exception as e:
            print(f"Error writing HOG summary CSV to {output_path}: {e}")
            sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description='Process OrthoFinder N0.tsv to classify HOGs and retrieve genes.')
    
    parser.add_argument("--input-tsv", "-i", type=str, required=True, help="Path to the OrthoFinder N0.tsv file.")
    parser.add_argument("--species", "-s", type=str, default=None, help="Focal species name (short name, e.g., 'Athaliana') for a demonstration test run.")
    parser.add_argument("--species-map-file", "-m", type=str, default=None, help="Path to an optional species short/long name mapping file (TSV format).")
    parser.add_argument("--summary-csv", "-c", type=str, default=None, help="Path to an output CSV file where a summary of HOG statistics for ALL species will be written.")
    
    # Note: The 'large_file' flag is present for compatibility but doesn't change
    # the behavior of the in-memory parsing in this script.
    parser.add_argument('--large_file', '-l', action='store_true', 
                        help='Included for pipeline compatibility, but this script uses in-memory parsing.')
                        
    args = parser.parse_args()
    
    # 1. Load and Parse Data
    sorter = N0OrthoGeneSorter(args.input_tsv, args.species_map_file)
    
    # 2. Summary CSV Generation
    if args.summary_csv:
        # Pass the is_large_file flag to the summary generation method
        summary_data = sorter.generate_summary_data(is_large_file=args.large_file)
        if summary_data:
            sorter.write_summary_csv(args.summary_csv, summary_data)
        
    # 3. Demonstration of Core Functionality (only runs if -s is provided)
    if args.species:
        focal_species_short = args.species
        focal_species_long = sorter.get_long_species_name(focal_species_short)
        
        if focal_species_short not in sorter.all_species:
            print(f"Error: Focal species '{focal_species_short}' not found in the N0.tsv data.")
            return

        print(f"\n--- Analysis Summary for Focal Species: {focal_species_long} ({focal_species_short}) ---")
        
        # Calculate for demonstration output
        total_hogs = len(sorter.data[focal_species_short])
        
        # Pass is_large_file to demonstration methods
        unique_hogs = sorter.get_unique_hog_ids(focal_species_short, is_large_file=args.large_file)
        unique_genes = sorter.get_genes_for_hog_ids(focal_species_short, unique_hogs, is_large_file=args.large_file)
        
        non_unique_hogs = sorter.get_non_unique_hog_ids(focal_species_short, is_large_file=args.large_file)
        non_unique_genes = sorter.get_genes_for_hog_ids(focal_species_short, non_unique_hogs, is_large_file=args.large_file)

        print(f"1. Total HOGs for {focal_species_short}: {total_hogs}")
        print(f"2. Unique HOGs (Species-Specific): {len(unique_hogs)}")
        print(f"3. Genes in Unique HOGs: {len(unique_genes)}")
        print(f"4. Non-Unique HOGs (Shared/Multi-Species): {len(non_unique_hogs)}")
        print(f"5. Genes in Non-Unique HOGs: {len(non_unique_genes)}")
        
        # Example of Shared HOGs
        other_species_list = [s for s in sorter.all_species if s != focal_species_short][:3]
        if other_species_list:
             other_species_long_names = [sorter.get_long_species_name(s) for s in other_species_list]
             # Pass is_large_file to demonstration methods
             shared_hogs = sorter.get_shared_hog_ids(focal_species_short, other_species_list, is_large_file=args.large_file)
             print(f"6. Shared HOGs with {', '.join(other_species_long_names)}: {len(shared_hogs)}")
        
    elif not args.summary_csv:
        # Only print this if neither -s nor -c was used.
        print("\nNote: Please run with the '-s [species_name]' flag for demonstration or '-c [output.csv]' to generate the summary.")

if __name__ == "__main__":
    main()
