#!/usr/bin/env python3
import scvi
import scanpy as sc
import sys
import datetime
import pandas as pd

dir_model='/data/jake/scrna-embed/models/database-scvi_model'
path_model = f'{dir_model}/model.pt'
path_db_adata = 'embeddings/database-mat_with_latent.h5ad'

# load model
adata_db = sc.read_h5ad(path_db_adata)
model = scvi.model.SCVI.load(dir_model, adata=adata_db)
print(model)

# process queries
adata_q = sc.read_h5ad("queries.h5ad")
adata_q.layers["counts"] = adata_q.X.copy()  # preserve counts
adata_q.raw = adata_q  # freeze the state in `.raw`

print("Queries shape: ", adata_q.shape)

hvgs = pd.read_csv("db_hvg.txt",header=None).values.flatten().tolist()
print("Number of HVGs: ", len(hvgs))

# subset and order hvgs for model
adata_q = adata_q[:, hvgs].copy()
adata_q

print("Queries shape after HVGs: ", adata_q.shape)

scvi.model.SCVI.setup_anndata(
    adata_q,
    layer="counts",
    batch_key="batch"
)

# https://chanzuckerberg.github.io/cellxgene-census/notebooks/analysis_demo/comp_bio_scvi_model_use.html

scvi.model.SCVI.prepare_query_anndata(adata_q, dir_model)
vae_q = scvi.model.SCVI.load_query_data(
    adata_q,
    dir_model
)
vae_q.is_trained = True
latent = vae_q.get_latent_representation()
print(latent.shape)
adata_q.obsm["X_scVI"] = latent

# save the latent representation
adata_q.write("queries_with_latent.h5ad")


