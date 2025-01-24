#!/usr/bin/env python3

# do nn distance histogram by tissue

import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt

df = pd.read_csv('knn_results.tsv', sep='\t')

# do histogram tissue wise nn distance histj
colors = ['blue', 'red',  'purple', 'orange',  'black']
# for c, d in enumerate(df.groupby('query_tissue')):
#     group, data = d
#     data['distance'].hist(bins=30, color=colors[c])
#     plt.xlabel('L2 Distance')
#     plt.ylabel('Count')
#     plt.title(group)
#     plt.grid(False)
#     plt.savefig(f'nn_distance_hist_{group}.png')
#     plt.close()

# do histogram tissue wise nn distance hist 
# strata by hit tissue
# for c, d in enumerate(df.groupby('query_tissue')):
#     tissue, data = d
#     for d2 in data.groupby('hit_tissue'):
#         group, data2 = d2
#         data2['distance'].hist(bins=30, color=colors[c])
#         plt.xlabel('L2 Distance')
#         plt.ylabel('Count')
#         plt.title(f'{tissue} to {group}')
#         plt.grid(False)
#         plt.savefig(f'nn_distance_hist_{tissue}_to_{group}.png')
#         plt.close()


## subplots
# nrows, ncols = 5, 5

# # Create a figure and a set of subplots
# fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20))
# axes = axes.flatten()

# # Define colors for the plots
# colors = plt.cm.tab20.colors

# # Plot tissue-wise nn distance histograms
# plot_index = 0
# for c, (tissue, data) in enumerate(df.groupby('query_tissue')):
#     for group, data2 in data.groupby('hit_tissue'):
#         ax = axes[plot_index]
#         data2['distance'].hist(bins=30, color=colors[c], ax=ax)
#         ax.set_xlabel('L2 Distance')
#         ax.set_ylabel('Count')
#         ax.set_title(f'{tissue} to {group}')
#         ax.grid(False)
#         plot_index += 1
#         if plot_index >= len(axes):
#             break
#     if plot_index >= len(axes):
#         break


# Define the number of rows and columns for the subplot grid
nrows, ncols = 5, 5

# Create a figure and a set of subplots
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20))
axes = axes.flatten()

# Define colors for the plots
colors = plt.cm.tab20.colors

# Plot tissue-wise nn distance histograms
plot_index = 0
row_titles = []
col_titles = []
for c, (tissue, data) in enumerate(df.groupby('query_tissue')):
    row_titles.append(tissue)
    for group, data2 in data.groupby('hit_tissue'):
        if group not in col_titles:
            col_titles.append(group)
        ax = axes[plot_index]
        data2['distance'].hist(bins=30, color=colors[c], ax=ax)
        ax.set_ylabel('Count')
        ax.grid(False)
        plot_index += 1
        if plot_index >= len(axes):
            break
    if plot_index >= len(axes):
        break

# Add row and column subtitles
for ax, col in zip(axes[:ncols], col_titles):
    ax.annotate(col, xy=(0.5, 1.0), xytext=(0, 5),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')

for ax, row in zip(axes[::ncols], row_titles):
    ax.annotate(row, xy=(-0.3, 0.5), xytext=(-ax.yaxis.labelpad - 5, 0),
                xycoords=ax.yaxis.label, textcoords='offset points',
                size='large', ha='right', va='center')

plt.tight_layout()
# Add a super title at the bottom center of the figure
fig.text(0.5, 0.04, 'L2 Distance', ha='center', va='center', fontsize=16)
fig.subplots_adjust(bottom=0.1)
# Adjust layout and save the figure
plt.savefig('nn_distance_hist_subplots.png')
plt.show()