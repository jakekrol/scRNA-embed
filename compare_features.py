#!/usr/bin/env python3
import anndata as ad
import pandas as pd
import sys
import numpy
import scipy
import os

# read two .h5ad files
# and compare feature overlap

# Get the .h5ad files from command line arguments
h5ad_cxg = sys.argv[1]
h5ad_encode = sys.argv[2]
outfile = sys.argv[3]

# Load the .h5ad files
adata_cxg = ad.read_h5ad(h5ad_cxg)
adata_encode = ad.read_h5ad(h5ad_encode)

# Extract the features
features_cxg = set(adata_cxg.var['gene_name'])
# use 'gene_name' for cxg columns
# df = adata_cxg.to_df()
# df.columns = adata_cxg.var['gene_name']

features_encode = set(adata_encode.var['genes'])
print("Features in cxg:", len(features_cxg))
print("Features in encode:", len(features_encode))

# Compare the features
common_features = features_cxg.intersection(features_encode)
print("Common features:", len(common_features))

# write the common features to a file
with open(outfile, 'w') as f:
    for feature in common_features:
        f.write(feature + '\n')
