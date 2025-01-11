#!/usr/bin/env python3

import scanpy as sc
import scvi
from pytorch_lightning.loggers import TensorBoardLogger
import sys
import datetime

d = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

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
vae = scvi.model.SCVI(
    adata,
    n_latent=10,
    n_layers=2,
    n_hidden=128,
    dropout_rate=0.1
)

logger = TensorBoardLogger(save_dir="logs", name=f"{d}-scvi_training")

# Train the model
vae.train(
    accelerator="cpu",
    devices=20,
    strategy="ddp_find_unused_parameters_true",
    max_epochs=10,                 
    logger=logger,
    batch_size=128
)


# Save the latent space
adata.obsm["X_scVI"] = vae.get_latent_representation()
breakpoint()

# Save the trained model
vae.save(f'/data/jake/scrna-embed/models/{d}-scvi_model', overwrite=True)

# Optionally, save the AnnData object with the latent representation
adata.write('/data/jake/scrna-embed/embeddings/mat_with_latent.h5ad')
