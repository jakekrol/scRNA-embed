#!/usr/bin/env python3

# heatmap of cell types for mismatching tissues
import pandas as pd
import plotly.express as px
import sys


# read query and nn1 cell type
df = pd.read_csv('25_01_27-nn_cell_types.tsv', sep='\t', usecols=[0,1])
print(df)
# do heatmap of cell type pair counts
df['count'] = 1
df = df.groupby(['query_cell_type', 'hit_cell_type1']).count().reset_index()
df = df.pivot_table(index='query_cell_type', columns='hit_cell_type1', values='count', fill_value=0)
# pad columns and row indices to make heatmap square
df = df.reindex(columns=df.columns.union(df.index))
df = df.reindex(index=df.index.union(df.columns))
# set nan to 0
df = df.fillna(0)
print(df)
fig = px.imshow(df.values,
                labels=dict(x="Hit Cell Type", y="Query Cell Type", color="Count"),
                x=df.columns,
                y=df.index,
                title = 'Cell types of query and nearest neighbor',
                color_continuous_scale='Blues')
# make sure all ticks are visible
# set axes labels to fontsize 16
fig.update_layout(xaxis=dict(tickmode='linear'),
                    yaxis=dict(tickmode='linear'),
                    coloraxis_colorbar=dict(
                        title="Count",
                        x=1.02,  # Position the color bar just to the right of the plot
                        y=0.5,
                        len=0.75,
                        thickness=15
                    ))
# center title
fig.update_layout(title_x=0.5)
fig.show()
fig.write_image('25_01_27-celltype_heatmap-nn1.png')

# do the same for nn2
df = pd.read_csv('25_01_27-nn_cell_types.tsv', sep='\t', usecols=[0,2])
print(df)
# do heatmap of cell type pair counts
df['count'] = 1
df = df.groupby(['query_cell_type', 'hit_cell_type2']).count().reset_index()
df = df.pivot_table(index='query_cell_type', columns='hit_cell_type2', values='count', fill_value=0)
# pad columns and row indices to make heatmap square
df = df.reindex(columns=df.columns.union(df.index))
df = df.reindex(index=df.index.union(df.columns))
# set nan to 0
df = df.fillna(0)
print(df)
fig = px.imshow(df.values,
                labels=dict(x="Hit Cell Type", y="Query Cell Type", color="Count"),
                x=df.columns,
                y=df.index,
                title = 'Cell types of query and second nearest neighbor',
                color_continuous_scale='Blues')
# make sure all ticks are visible
fig.update_layout(xaxis=dict(tickmode='linear'),
                    yaxis=dict(tickmode='linear'),
                    coloraxis_colorbar=dict(
                        title="Count",
                        x=1.02,  # Position the color bar just to the right of the plot
                        y=0.5,
                        len=0.75,
                        thickness=15
                    ))
# center title
fig.update_layout(title_x=0.5)
# make labels bigger
fig.show()
fig.write_image('25_01_27-celltype_heatmap-nn2.png')


