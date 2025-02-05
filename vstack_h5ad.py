#!/usr/bin/env python3
import scanpy as sc
import anndata as ad
import glob
import pandas as pd
import sys, os

# List of file paths
# files = glob.glob("*.h5ad")
# files = pd.read_csv("query_db.tsv", sep="\t", usecols=[2]).iloc[:,0].tolist()
files = pd.read_csv("query_db.tsv", sep="\t", usecols=[1]).iloc[:,0].tolist()
print(files)

# Load all .h5ad files
adata_list = [sc.read_h5ad(fp) for fp in files]
print("Total observations:", sum([adata.n_obs for adata in adata_list]))
for i, adata in enumerate(adata_list):
    adata.obs_names_make_unique()  # Ensure unique cell names
    # prepend filename
    prefix = '_'.join(files[i].split(".")[:2])
    adata.obs_names = [f"{prefix}_{obs}" for obs in adata.obs_names]
    adata.obs['batch'] = prefix
print("Total observations:", sum([adata.n_obs for adata in adata_list]))


# Concatenate them vertically (stack rows/observations)
adata_combined = ad.concat(adata_list, axis=0, join='inner')  # Use 'inner' to match var_names
adata_combined.obs_names_make_unique()  # Ensure unique cell names, again
print("Combined observations:", adata_combined.n_obs)

# Save the combined AnnData object to a new .h5ad file (optional)
# adata_combined.write("combined.h5ad")
# adata_combined.write("database.h5ad")
adata_combined.write("queries.h5ad")