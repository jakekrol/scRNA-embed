#!/usr/bin/env python3
import sys, os
import numpy as np
import sklearn
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
import time
import torch

k = 5

# Load the .h5ad file
file = sys.argv[1]
outfile_ptn = sys.argv[2]
adata = sc.read_h5ad(file)

X = adata.X
X_scVI = adata.obsm['X_scVI']

# raw
# adata.raw.X
# gene subset
# adata.layers['counts']

# do knn for subset of rows
np.random.seed(0)
m = 200
sample_indices = np.random.choice(X.shape[0], size=m, replace=False)
queries = X[sample_indices, :]

def knn_search(knn, query_row):
    distances, indices = knn.kneighbors([query_row])
    return distances, indices


nn_X = {}
nn_X_scVI = {}
start = time.time()
for i in sample_indices:
    knn_X = NearestNeighbors(n_neighbors=k, metric='cosine')
    knn_X.fit(np.delete(X, i, axis=0))
    knn_X_scVI = NearestNeighbors(n_neighbors=k, metric='cosine')
    knn_X_scVI.fit(np.delete(X_scVI, i, axis=0))

    q_X = X[i]
    q_X_scVI = X_scVI[i]
    distances, indices = knn_search(knn_X, q_X)
    nn_X[i] = {'dist': distances.reshape(-1), 'idx':indices.reshape(-1)}
    distances, indices = knn_search(knn_X_scVI, q_X_scVI)
    nn_X_scVI[i] = {'dist': distances.reshape(-1), 'idx':indices.reshape(-1)}
end = time.time()
print(f"Time for X_scVI: {time.time() - end}")

def recall_at_topk(nn_truth, nn_pred, k=k):
    max_recall = len(nn_truth.items()) * k
    recall_at_k = 0
    for i in nn_truth.keys():
        truth = set(nn_truth[i]['idx'])
        pred = set(nn_pred[i]['idx'])
        inter = truth & pred
        recall_at_k += len(inter)
    return recall_at_k / max_recall
recall = recall_at_topk(nn_X, nn_X_scVI)
print(f"Recall at k: {recall}")

def hit_rate(nn_truth, nn_pred, k=k):
    hits = 0
    for i in nn_truth.keys():
        top1 = {nn_truth[i]['idx'][0]}
        pred = set(nn_pred[i]['idx'])
        inter = top1 & pred
        if len(inter) > 0:
            hits += 1
    return hits / len(nn_truth.keys())
hit = hit_rate(nn_X, nn_X_scVI)
print(f"Hit rate: {hit}")


torch.save(nn_X, f'{outfile_ptn}_X.pt')
torch.save(nn_X_scVI, f'{outfile_ptn}_X_scVI.pt')

