#!/usr/bin/env python3

import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv('knn_results.tsv', sep='\t')

# do histogram tissue wise nn distance histj
colors = ['blue', 'red',  'purple', 'orange',  'black']

fig, axes = plt.subplots(nrows=5, ncols=5, figsize=(20, 20))
tissue_idx = {'brain': 0, 'eye': 1, 'heart': 2, 'lung': 3, 'testes': 4}

dfs = {}
for t_q, g in df.groupby('query_tissue'):
    for t_h, g2 in g.groupby('hit_tissue'):
        axes[tissue_idx[t_q], tissue_idx[t_h]].hist(g2['distance'].values, bins=30, color=colors[tissue_idx[t_q]])
        axes[tissue_idx[t_q], tissue_idx[t_h]].set_xlabel('L2 Distance')
        axes[tissue_idx[t_q], tissue_idx[t_h]].set_ylabel('Count')
        axes[tissue_idx[t_q], tissue_idx[t_h]].set_title(f'{t_q} to {t_h}')
        axes[tissue_idx[t_q], tissue_idx[t_h]].grid(False)

plt.tight_layout()
plt.savefig(f'nn_dist_hist.png')
