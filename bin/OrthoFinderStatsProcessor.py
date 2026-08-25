#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict

class OrthoFinderStatsProcessor:
    def __init__(self, tsv_file, csv_file, species_map_file=None):
        """
        Initialize the OrthoFinderStatsProcessor with the file paths.
        """
        self.tsv_file = tsv_file
        self.csv_file = csv_file
        self.species_map_file = species_map_file
        self.species_map = self.load_species_map()
        
        # Column names standardized to Title_Case_With_Underscores throughout
        # (this dict previously had "Species_Specific-Clustered" with a
        # hyphen, inconsistent with "Species_Specific_Unclustered" right
        # next to it in the same output table).
        self.MAPPINGS = {
            "Total_Proteins": "Number of genes",
            "Clustered_Total": "Number of genes in orthogroups",
            "Species_Specific_Unclustered": "Number of unassigned genes",
            "Species_Specific_Clustered": "Number of genes in species-specific orthogroups"
        }
        
    def load_species_map(self):
        """
        Loads a species map file (CSV or TSV) to map base filenames to correct names.
        This version auto-detects the delimiter and reads the header to find the correct columns.
        It expects headers containing 'Species' and 'Basename' (case-insensitive).
        """
        if not self.species_map_file:
            return {}
        
        s_map = {}
        try:
            with open(self.species_map_file, 'r', newline='', encoding='utf-8-sig') as f:
                sniffer = csv.Sniffer()
                dialect = sniffer.sniff(f.read(1024))
                f.seek(0)
                
                reader = csv.reader(f, dialect)
                
                header = [h.lower().strip() for h in next(reader)]
                
                try:
                    # Map file has 'Species' and 'Basename'
                    name_col_idx = header.index('species')
                    base_col_idx = header.index('basename')
                except ValueError:
                    print("Warning: Map file header must contain 'Species' and 'Basename'. Skipping renaming.")
                    return {}

                # The stats processor needs a map from BASENAME -> SPECIES NAME
                for row in reader:
                    if len(row) > max(name_col_idx, base_col_idx):
                        basename = row[base_col_idx].strip()
                        species_name = row[name_col_idx].strip()
                        if basename and species_name:
                            s_map[basename] = species_name
            
            print(f"Loaded {len(s_map)} species mappings from {self.species_map_file}.")
            return s_map
        except FileNotFoundError:
            print(f"Warning: Species map file not found at {self.species_map_file}. Skipping renaming.")
            return {}
        except Exception as e:
            print(f"Warning: Error loading species map: {e}. Skipping renaming.")
            return {}

    def process_stats(self):
        """
        Process the TSV file, extracting metrics and renaming species if a map is provided.
        """
        try:
            data = {}
            species_basenames = []
            required_metrics = set(self.MAPPINGS.values())
            reverse_mappings = {v: k for k, v in self.MAPPINGS.items()}
            
            with open(self.tsv_file, 'r', newline='') as tsv_f:
                reader = csv.reader(tsv_f, delimiter='\t')
                
                header_row = next(reader, None)
                if not header_row:
                    print(f"Error: The file {self.tsv_file} is empty.")
                    return

                species_basenames = header_row[1:]
                for basename in species_basenames:
                    data[basename] = {}

                for row in reader:
                    if not row or not row[0].strip():
                        continue
                    
                    metric_key_tsv = row[0].strip()
                    if metric_key_tsv in required_metrics:
                        metric_key_csv = reverse_mappings[metric_key_tsv]
                        metric_values = row[1:]
                        
                        for i, basename in enumerate(species_basenames):
                            if i < len(metric_values):
                                data[basename][metric_key_csv] = metric_values[i]
            
            headers = ['Species'] + list(self.MAPPINGS.keys())
            
            with open(self.csv_file, 'w', newline='') as csv_f:
                writer = csv.DictWriter(csv_f, fieldnames=headers)
                writer.writeheader()
                
                for basename in species_basenames:
                    row_data = {'Species': self.species_map.get(basename, basename)}
                    row_data.update(data[basename])
                    writer.writerow(row_data)
            
            print(f"Successfully processed stats and saved to {self.csv_file}")

        except Exception as e:
            print(f"An unexpected error occurred during processing: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert OrthoFinder stats TSV to CSV format.')
    parser.add_argument('tsv_file', type=str, help='Input TSV file path')
    parser.add_argument('csv_file', type=str, help='Output CSV file path')
    parser.add_argument('--species_map', '-s', type=str, default=None, help='Optional CSV/TSV file to map species basenames to full names.')
    args = parser.parse_args()
    
    processor = OrthoFinderStatsProcessor(args.tsv_file, args.csv_file, args.species_map)
    processor.process_stats()
