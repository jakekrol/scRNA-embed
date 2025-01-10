#!/usr/bin/env python3

import scanpy as sc
import scvi

# Load the .h5ad file
file_path = '/data/jake/scrna-embed/data/mat.h5ad'
adata = sc.read_h5ad(file_path)

adata.layers["counts"] = adata.X.copy()  # preserve counts
adata.raw = adata  # freeze the state in `.raw`

# Filter and preprocess the data
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1200,
    subset=True,
    layer="counts",
    flavor="seurat_v3"
)

# Setup AnnData for SCVI
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts"
)
# breakpoint()

# Initialize the SCVI model
vae = scvi.model.SCVI(adata)

# Train the model
vae.train(accelerator="cpu", devices=10, strategy="ddp_find_unused_parameters_true")

# Save the latent space
adata.obsm["X_scVI"] = vae.get_latent_representation()

# Save the trained model
vae.save('/data/jake/scrna-embed/models/scvi_model', overwrite=True)

# Optionally, save the AnnData object with the latent representation
adata.write('/data/jake/scrna-embed/data/mat_with_latent.h5ad')
