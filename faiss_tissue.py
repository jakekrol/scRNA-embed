#!/usr/bin/env python3
import numpy as np
import faiss
import pandas as pd
import time

# Goal: get tissue distribution of nearest neighbors in db

# k = 1

# # setup db
# X_db = np.load('db_vectors.npy')
# idx_db = np.load('db_idx.npy',allow_pickle=True)
# t_db = np.load('db_tissue.npy',allow_pickle=True)

# # setup query
# X_q = np.load('query_vectors.npy')
# idx_q = np.load('query_idx.npy',allow_pickle=True)
# t_q = np.load('query_tissue.npy',allow_pickle=True)

# # create a FAISS index
# index = faiss.IndexFlatL2(X_db.shape[1])  # L2 distance (Euclidean distance)
# index.add(X_db)  # Add vectors to the index

# # Perform a nearest neighbor query for all query vectors
# start = time.time()
# distances, faiss_indices = index.search(X_q, k)
# end = time.time()
# print(f"Time to query: {end-start}")
# distances = distances.flatten()
# faiss_indices = faiss_indices.flatten()

# # Map FAISS indices back to original indices
# hit_idx = idx_db[faiss_indices]
# t_hit = t_db[faiss_indices]

# # write query idx, query tissue, hit idx, hit tissue, and hit distance to table
# results = []
# for i in range(len(idx_q)):
#     results.append([idx_q[i], t_q[i], hit_idx[i], t_hit[i], distances[i]])
# df = pd.DataFrame(results, columns=["query_index", "query_tissue", "hit_index", "hit_tissue", "distance"])
# df.to_csv("knn_results2.tsv", index=False,sep='\t')

k = 2

# setup db
X_db = np.load('db_vectors.npy')
idx_db = np.load('db_idx.npy',allow_pickle=True)
t_db = np.load('db_tissue.npy',allow_pickle=True)

# setup query
X_q = np.load('query_vectors.npy')
idx_q = np.load('query_idx.npy',allow_pickle=True)
t_q = np.load('query_tissue.npy',allow_pickle=True)

# create a FAISS index
index = faiss.IndexFlatL2(X_db.shape[1])  # L2 distance (Euclidean distance)
index.add(X_db)  # Add vectors to the index

# Perform a nearest neighbor query for all query vectors
start = time.time()
distances, faiss_indices = index.search(X_q, k)
end = time.time()
print(f"Time to query: {end-start}")
# paste
stacked = np.hstack((distances, faiss_indices))
df_nn = pd.DataFrame(stacked, columns=["distance1", "distance2", "idx_hit1", "idx_hit2"])
results = []
for i in range(len(idx_q)):
    results.append([idx_q[i], t_q[i], idx_db[int(df_nn.loc[i,'idx_hit1'])], idx_db[int(df_nn.loc[i,'idx_hit2'])], 
                    df_nn.loc[i,'distance1'], df_nn.loc[i,'distance2']])
    
df = pd.DataFrame(results, columns=["query_index", "query_tissue", "hit_index1", "hit_index2", "distance1", "distance2"])
df.to_csv("25_01_27-knn_results.tsv", index=False,sep='\t')