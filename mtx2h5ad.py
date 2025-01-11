#!/usr/bin/env python3

import scanpy as sc
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description="Convert .mtx file and metadata to .h5ad format.")
parser.add_argument("matrix", help="Path to the .mtx file")
parser.add_argument("genes", help="Path to the genes.tsv file")
parser.add_argument("barcodes", help="Path to the barcodes.tsv file")

# Parse the arguments
args = parser.parse_args()

# Read the exported .mtx file and metadata
adata = sc.read_mtx(args.matrix)  # Load the sparse matrix
adata.var_names = [line.strip() for line in open(args.genes)]  # Gene names
adata.obs_names = [line.strip() for line in open(args.barcodes)]  # Cell barcodes

# Save to .h5ad format
adata.write("matrix.h5ad")

print("Conversion complete. File saved as 'matrix.h5ad'.")

