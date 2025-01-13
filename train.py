#!/usr/bin/env python3

import scanpy as sc
import scvi
from pytorch_lightning.loggers import TensorBoardLogger
from pytorch_lightning import Trainer
import sys
import datetime
import gc
import torch.distributed as dist

d = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

# hp
# n_latent = 10
# n_layers = 10
# n_hidden = 128
# max_epochs = 30
# batch_size = 256
hps = {
    "n_latent": 500,
    "n_layers": 2,
    "n_hidden": 256,
    "max_epochs": 50,
    "batch_size": 256,
    "preproc": True,
    'batch_label': True
}
print(hps)
s_hps = '-'.join(map(str, hps.values()))
print(s_hps)

# Load the .h5ad file
# file_path = '/data/jake/scrna-embed/data/mat.h5ad'
file_path = 'database.h5ad'
adata = sc.read_h5ad(file_path)

adata.layers["counts"] = adata.X.copy()  # preserve counts
adata.raw = adata  # freeze the state in `.raw`

# Filter and preprocess the data
if hps['preproc']:
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=1200,
        subset=True,
        layer="counts",
        flavor="seurat_v3"
    )
if hps['batch_label']:
    # already have batch set in the data
    # uncomment this if you want to set batch to a constant
    # adata.obs["batch"] = ['batch_1' for i in range(adata.shape[0])]
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="batch"
    )
else:
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts"
    )

# Initialize the SCVI model
vae = scvi.model.SCVI(
    adata,
    n_latent=hps["n_latent"],
    n_layers=hps["n_layers"],
    n_hidden=hps["n_hidden"],
    dropout_rate=0.1
)

# logger = TensorBoardLogger(save_dir="logs", name=f"{s_hps}-scvi_training")
logger = TensorBoardLogger(save_dir="logs", name=f"database-{s_hps}-scvi_training")

# Train the model
# unfornately, distributing over multiple CPUs
# creates strange behavior post training. I don't know how to "clean it up"
# so I'm just training on a single CPU for now.
# Just train on GPU locally instead.
vae.train(
    accelerator="cpu",
    devices=1, # training
    # strategy="ddp_find_unused_parameters_true",
    # max_epochs=5,
    # max_epochs=50,
    max_epochs=hps["max_epochs"],
    logger=logger,
    batch_size=hps["batch_size"]
)


# Save the latent space
adata.obsm["X_scVI"] = vae.get_latent_representation()

# Save the trained model
# vae.save(f'/data/jake/scrna-embed/models/{s_hps}-scvi_model', overwrite=True)
vae.save(f'models/database-scvi_model', overwrite=True)

# Optionally, save the AnnData object with the latent representation
# adata.write(f'/data/jake/scrna-embed/embeddings/{s_hps}-mat_with_latent.h5ad')
adata.write(f'embeddings/database-mat_with_latent.h5ad')
