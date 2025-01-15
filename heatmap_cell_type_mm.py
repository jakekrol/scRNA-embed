#!/usr/bin/env python3

# heatmap of cell types for mismatching tissues
import pandas as pd
import plotly.express as px

df = pd.read_csv('tissue_mismatch_cell_types.tsv', sep='\t')
df = df[['query_cell_type', 'hit_cell_type']]

# do heatmap of cell type pair counts
df['count'] = 1
df = df.groupby(['query_cell_type', 'hit_cell_type']).count().reset_index()
df = df.pivot_table(index='query_cell_type', columns='hit_cell_type', values='count', fill_value=0)

# Convert the pivot table to a long-form DataFrame for Plotly
df_long = df.reset_index().melt(id_vars='query_cell_type', var_name='hit_cell_type', value_name='count')

# Create the heatmap using Plotly
fig = px.imshow(df.values, 
                labels=dict(x="Hit Cell Type", y="Query Cell Type", color="Count"),
                x=df.columns,
                y=df.index,
                title = 'Cell types of tissue mismatch',
                color_continuous_scale='Blues')

# Update the layout to ensure all tick labels are shown and move the legend
fig.update_layout(
    xaxis=dict(tickmode='linear'),
    yaxis=dict(tickmode='linear'),
    coloraxis_colorbar=dict(
        title="Count",
        x=1.02,  # Position the color bar just to the right of the plot
        y=0.5,
        len=0.75,
        thickness=15
    )
)

# Show and save the plot
fig.show()
fig.write_image('heatmap_cell_type_mm.png')

