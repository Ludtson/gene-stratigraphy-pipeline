#!/usr/bin/env python3

"""
FastaGeneFilter.py

This script filters gene sequences from a FASTA file (CDS or protein) based on a list of gene IDs
and writes the filtered genes to a new FASTA file.

The script supports either including genes in the list (--include) or excluding them (default).
It assumes the gene ID in the FASTA header is the part before the first period or space.
"""

import os
import argparse

class GeneSequenceExtractor:
    def __init__(self, gene_ids_file, fasta_file, output_dir, include=False, buffer_size=10000, output_filename=None):
        """
        Initializes the extractor.
        
        Args:
            gene_ids_file (str): Path to the file with gene IDs.
            fasta_file (str): Path to the input FASTA file.
            output_dir (str): Directory to save the output.
            include (bool): If True, include genes from the list. If False, exclude them.
            buffer_size (int): Number of lines to buffer before writing.
            output_filename (str, optional): Desired name for the output file. 
                                             If None, a default name will be generated.
        """
        self.gene_ids_file = gene_ids_file
        self.fasta_file = fasta_file
        self.output_dir = output_dir
        self.include = include
        self.buffer_size = buffer_size
        # Store the custom output filename
        self.output_filename = output_filename

    def _load_gene_ids(self):
        """Loads gene IDs from the file into a set for efficient lookup."""
        try:
            with open(self.gene_ids_file, 'r') as f:
                return {line.strip() for line in f}
        except FileNotFoundError:
            print(f"Error: Gene IDs file not found at {self.gene_ids_file}")
            return set()

    def _extract_gene_id(self, header):
        """Extracts the gene ID from a FASTA header line."""
        # Assumes ID is the first part of the header before a space or period
        header_part = header.lstrip('>').split()[0]
        return header_part

    def filter_fasta(self):
        """
        Reads the FASTA file, filters genes based on the ID list, and writes to a new file.
        """
        gene_ids = self._load_gene_ids()
        if not gene_ids:
            print("No gene IDs to filter by. Exiting.")
            return 0, 0

        # --- Determine the output file path ---
        if self.output_filename:
            # Use the user-provided filename
            output_file_path = os.path.join(self.output_dir, self.output_filename)
        else:
            # Generate a default filename if none was provided
            base_name = os.path.basename(self.fasta_file)
            output_file_path = os.path.join(self.output_dir, f"filtered_{base_name}")
        
        os.makedirs(self.output_dir, exist_ok=True)
        
        print(f"Filtering genes from {os.path.basename(self.fasta_file)}...")
        print(f"Output will be saved to: {output_file_path}")

        total_genes_processed = 0
        genes_written = 0
        buffer = []
        keep_sequence = False

        try:
            with open(self.fasta_file, 'r') as infile, open(output_file_path, 'w') as outfile:
                for line in infile:
                    if line.startswith('>'):
                        total_genes_processed += 1
                        gene_id = self._extract_gene_id(line)
                        
                        # Determine if this sequence should be kept
                        if self.include:
                            keep_sequence = gene_id in gene_ids
                        else:
                            keep_sequence = gene_id not in gene_ids
                        
                        if keep_sequence:
                            genes_written += 1
                            # Write any buffered sequence from the previous kept gene
                            if buffer:
                                outfile.writelines(buffer)
                                buffer = []
                            buffer.append(line)
                    elif keep_sequence:
                        buffer.append(line)
                        if len(buffer) >= self.buffer_size:
                            outfile.writelines(buffer)
                            buffer = []
                
                # Write any remaining buffer content
                if buffer:
                    outfile.writelines(buffer)

        except FileNotFoundError:
            print(f"Error: Input FASTA file not found at {self.fasta_file}")
            return 0, 0
        
        return total_genes_processed, genes_written

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter a FASTA file based on a list of gene IDs.")
    
    parser.add_argument("--gene-ids", "-g", type=str, required=True, help="Path to the gene IDs text file.")
    parser.add_argument("--input-fasta", "-f", type=str, required=True, help="Path to the FASTA file (CDS or protein).")
    parser.add_argument("--output-dir", "-o", type=str, required=True, help="Directory to save the filtered FASTA file.")
    parser.add_argument("--output-filename", "-n", type=str, default=None, help="Optional name for the output FASTA file.")
    parser.add_argument("--include", "-I", action="store_true", help="Flag to INCLUDE genes in the list (default is to EXCLUDE).")
    parser.add_argument("--buffer-size", "-b", type=int, default=10000, help="Buffer size for writing.")

    args = parser.parse_args()
    
    extractor = GeneSequenceExtractor(
        args.gene_ids, 
        args.input_fasta, 
        args.output_dir, 
        args.include, 
        args.buffer_size,
        args.output_filename  # Pass the new argument
    )
    total_genes, written_genes = extractor.filter_fasta()
    
    print("\n--- Filtering Summary ---")
    print(f"Input FASTA File: {os.path.basename(args.input_fasta)}")
    print(f"Filter Mode: {'INCLUDE list' if args.include else 'EXCLUDE list'}")
    print(f"Total Genes Processed: {total_genes}")
    print(f"Genes Written to Output File: {written_genes}")