#!/usr/bin/env python3
import pandas as pd
from scipy.io import mmread
import numpy as np
import os, sys

matrix_file = sys.argv[1]
barcodes_file = sys.argv[2]
features_file = sys.argv[3]
out_file = sys.argv[4]

# Step 1: Load the matrix
sparse_matrix = mmread(matrix_file)

# Step 2: Load the barcodes (cells) and features (genes)
barcodes = pd.read_csv(barcodes_file, header=None, sep='\t')
features = pd.read_csv(features_file, header=None, sep='\t')

# Step 3: Convert the sparse matrix to a dense format if needed
# (Keep it sparse if the matrix is too large to fit into memory)
dense_matrix = sparse_matrix.toarray()  # Optional: use `.todense()` for dense matrix

# Step 4: Create a DataFrame with gene names as columns and cell barcodes as rows
cell_by_gene_matrix = pd.DataFrame(dense_matrix.T, index=barcodes[0], columns=features[1])

# Step 5: Save the matrix to a .npz file along with index and column names
np.savez(out_file, data=cell_by_gene_matrix.values, index=cell_by_gene_matrix.index.values, columns=cell_by_gene_matrix.columns.values)

print(f"Cell-by-gene matrix saved to {out_file}")
