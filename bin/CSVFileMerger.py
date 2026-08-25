#!/usr/bin/env python3

import argparse
import csv
from collections import OrderedDict
import io 

# Increased field size limit for large bioinformatics files
csv.field_size_limit(1000000) 

class CSVFileMerger:
    def __init__(self, output_file, merge_cols, input_files, species_map_file=None, fill_value='NA'):
        """
        Initialize the CSVFileMerger with the file paths and merge columns.

        fill_value: written for any column a given row has no value for
        (e.g. a key present in one input file but not another -- a common
        case here is a species that appears in raw OrthoFinder stats but
        was never analyzed as a focal species, so it has no per-node data
        to merge in). Defaults to the string "NA" to match the missing-value
        convention used elsewhere in this pipeline's output, rather than an
        empty string that reads as a gap with no explanation.
        """
        self.output_file = output_file
        self.merge_cols = merge_cols
        self.input_files = input_files
        self.species_map_file = species_map_file
        self.species_map = self.load_species_map()
        self.fill_value = fill_value
        self.SNIFF_BYTES = 1024 * 10
        self.primary_merge_key = None # Store the first merge key

    # (load_species_map method remains the same)
    def load_species_map(self):
        """
        Loads the species map file (TSV) to rename short species codes to full names.
        """
        if not self.species_map_file:
            return None
        
        species_map = {}
        try:
            with open(self.species_map_file, 'r', newline='') as f:
                reader = csv.reader(f, delimiter='\t')
                header = next(reader, None)
                
                full_name_idx = 0
                short_name_idx = 1
                
                if header:
                    try:
                        full_name_idx = header.index('Species_name')
                        short_name_idx = header.index('short_name')
                    except ValueError:
                        pass
                
                for row in reader:
                    if len(row) > max(full_name_idx, short_name_idx) and row[full_name_idx].strip() and row[short_name_idx].strip():
                        short_name = row[short_name_idx].strip()
                        full_name = row[full_name_idx].strip()
                        species_map[short_name] = full_name

            print(f"Loaded {len(species_map)} species mappings from {self.species_map_file}.")
            return species_map
            
        except FileNotFoundError:
            print(f"Warning: Species map file not found at {self.species_map_file}. Skipping renaming.")
            return None
        except Exception as e:
            print(f"Warning: Error loading species map: {e}. Skipping renaming.")
            return None

    def _sniff_delimiter(self, filepath):
        """
        Automatically determines the delimiter of a file (prioritizing tab and comma).
        """
        try:
            with open(filepath, 'r', newline='') as f:
                sample = f.read(self.SNIFF_BYTES)
                if not sample:
                    return ','
                
                sniffer = csv.Sniffer()
                try:
                    dialect = sniffer.sniff(sample, delimiters=['\t', ','])
                    return dialect.delimiter
                except csv.Error:
                    if '\t' in sample:
                        return '\t'
                    return ','
        except Exception as e:
            print(f"Error during delimiter sniffing for {filepath}: {e}. Defaulting to comma.")
            return ','


    def merge_files(self):
        """
        Merge the tabular files (CSV/TSV) based on the specified columns, using auto-detected delimiters.
        """
        merged_data = OrderedDict() 
        # Changed from set to list to preserve column order based on file input order
        ordered_fieldnames = [] 
        
        # --- Handle Merge Column Logic ---
        effective_merge_cols = []
        if len(self.merge_cols) == 1:
            effective_merge_cols = self.merge_cols * len(self.input_files)
        elif len(self.merge_cols) == len(self.input_files):
            effective_merge_cols = self.merge_cols
        else:
            print(f"Error: The number of merge columns ({len(self.merge_cols)}) must be 1 or match the number of input files ({len(self.input_files)}).")
            return
        
        # Capture the primary merge key
        if effective_merge_cols:
            self.primary_merge_key = effective_merge_cols[0]


        # 1. Read and Merge Data
        for csv_file, merge_col in zip(self.input_files, effective_merge_cols):
            delim = self._sniff_delimiter(csv_file) 
            print(f"Processing {csv_file} with delimiter: '{delim}'")
            
            try:
                with open(csv_file, 'r', newline='') as f:
                    
                    reader = csv.DictReader(f, delimiter=delim) 
                    
                    if merge_col not in reader.fieldnames:
                        print(f"Error: '{merge_col}' column not found in {csv_file}. Skipping file.")
                        continue
                    
                    # Store new fieldnames in order
                    for field in reader.fieldnames:
                        if field not in ordered_fieldnames:
                            ordered_fieldnames.append(field)

                    for row in reader:
                        key = row[merge_col]
                        
                        if key not in merged_data:
                            merged_data[key] = row
                        else:
                            merged_data[key].update(row)

            except FileNotFoundError:
                print(f"Warning: Input file not found: {csv_file}. Skipping.")
                continue
            except Exception as e:
                print(f"An error occurred while processing {csv_file}: {e}. Skipping.")
                continue

        # 2. Determine final header and Write Output (always as standard CSV)
        if merged_data:
            
            # Use the ordered list
            final_fieldnames = ordered_fieldnames
            
            # Ensure the primary merge key is the first column
            if self.primary_merge_key and self.primary_merge_key in final_fieldnames:
                final_fieldnames.remove(self.primary_merge_key)
                final_fieldnames.insert(0, self.primary_merge_key)
                
            with open(self.output_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=final_fieldnames) 
                writer.writeheader()
                
                for key, row in merged_data.items():
                    # Create the row using the final fieldnames to enforce order
                    output_row = {field: row.get(field, self.fill_value) for field in final_fieldnames}
                    
                    if self.species_map and 'Species' in output_row and self.primary_merge_key == 'Species':
                        # If the key is 'Species', apply renaming
                        output_row['Species'] = self.species_map.get(key, key)
                    elif self.primary_merge_key in output_row:
                        output_row[self.primary_merge_key] = key 
                        
                    writer.writerow(output_row)
            
            print(f"Successfully merged {len(self.input_files)} file(s) into {self.output_file}")
        else:
            print("Error: No data to merge. Output file not created.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge tabular files (CSV/TSV) based on columns, using auto-detected delimiters.')
    
    parser.add_argument('--output_file', '-o', type=str, required=True, help='The output CSV file path.')
    parser.add_argument('--input_files', '-i', type=str, required=True, help='Comma-separated list of input file paths.')
    parser.add_argument('--merge_cols', '-c', type=str, required=True, help='Comma-separated list of merge columns (one per file, or one for all).')
    parser.add_argument('--species_map', '-s', type=str, default=None, help='Optional TSV file to map species short names to full names.')
    
    args = parser.parse_args()
    
    input_files = args.input_files.split(',')
    merge_cols = args.merge_cols.split(',')
    
    merger = CSVFileMerger(args.output_file, merge_cols, input_files, args.species_map)
    merger.merge_files()
    