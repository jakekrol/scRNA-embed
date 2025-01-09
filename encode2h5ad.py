#!/usr/bin/env python3
import pandas as pd
from scipy.io import mmread
import numpy as np
import anndata as ad
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

# Step 4: Create an AnnData object
adata = ad.AnnData(X=dense_matrix.T)
adata.obs['barcodes'] = barcodes[0].values
adata.var['genes'] = features[1].values

# Step 5: Save the AnnData object to an .h5ad file
adata.write(out_file)

print(f"AnnData object saved to {out_file}")
