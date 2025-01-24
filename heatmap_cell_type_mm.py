#!/usr/bin/env python3

# heatmap of cell types for mismatching tissues
import pandas as pd
import plotly.express as px

# df = pd.read_csv('tissue_mismatch_cell_types.tsv', sep='\t')
# df = df[['query_cell_type', 'hit_cell_type']]

# # do heatmap of cell type pair counts
# df['count'] = 1
# df = df.groupby(['query_cell_type', 'hit_cell_type']).count().reset_index()
# df = df.pivot_table(index='query_cell_type', columns='hit_cell_type', values='count', fill_value=0)

# # Convert the pivot table to a long-form DataFrame for Plotly
# df_long = df.reset_index().melt(id_vars='query_cell_type', var_name='hit_cell_type', value_name='count')

# # Create the heatmap using Plotly
# fig = px.imshow(df.values, 
#                 labels=dict(x="Hit Cell Type", y="Query Cell Type", color="Count"),
#                 x=df.columns,
#                 y=df.index,
#                 title = 'Cell types of tissue mismatch',
#                 color_continuous_scale='Blues')

# # Update the layout to ensure all tick labels are shown and move the legend
# fig.update_layout(
#     xaxis=dict(tickmode='linear'),
#     yaxis=dict(tickmode='linear'),
#     coloraxis_colorbar=dict(
#         title="Count",
#         x=1.02,  # Position the color bar just to the right of the plot
#         y=0.5,
#         len=0.75,
#         thickness=15
#     )
# )

# # Show and save the plot
# fig.show()
# fig.write_image('heatmap_cell_type_mm.png')


df = pd.read_csv('tissue_mismatch_cell_types.tsv', sep='\t')
df = df[['query_tissue', 'hit_cell_type']]

# do heatmap of cell type pair counts
df['count'] = 1
df = df.groupby(['query_tissue', 'hit_cell_type']).count().reset_index()
df = df.pivot_table(index='query_tissue', columns='hit_cell_type', values='count', fill_value=0)
print(df)
# only keep hit cell types with any values > 10
df = df.loc[:, (df >= 20).any(axis=0)]
print(df)
df = df.drop(columns=['Vascular smooth muscle cell', 'Elongated spermatid'])
# get 10 largest values in dataframe

# Convert the pivot table to a long-form DataFrame for Plotly
# df_long = df.reset_index().melt(id_vars='query_tissue', var_name='hit_cell_type', value_name='count')

# remove "cell" from column names
df.columns = df.columns.str.replace(' cell', '')
df.columns = df.columns.str.replace('blast', 'blst')
df.rename(columns={'T': 'T cell'}, inplace=True)
# uppercase hit_cell_type
df.index = df.index.str.title()


# Create the heatmap using Plotly
fig = px.imshow(df.values, 
                labels=dict(x="Hit Cell Type", y="Query Cell Tissue", color="Count"),
                x=df.columns,
                y=df.index,
                title = 'Predicted cell types of tissue mismatch nearest neighbors',
                color_continuous_scale='Blues',
                text_auto=True)
            

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
# rotate x-ticks  90
# bold all text
size_title=36
size_axis=30
size_colorbar=30
size_ticks=26
fig.update_layout(font=dict(family="Arial", size=18, color='black'))
fig.update_xaxes(tickangle=90)
# make title bigger
fig.update_layout(title_font_size=size_title)
# make x and y axis labels bigger
fig.update_layout(xaxis_title_font_size=size_axis)
fig.update_layout(yaxis_title_font_size=size_axis)
# make ticks bigger
fig.update_layout(xaxis_tickfont_size=size_ticks)
fig.update_layout(yaxis_tickfont_size=size_ticks)
# make scale color scale text bigger
fig.update_layout(coloraxis_colorbar_title_font_size=size_colorbar)
# center title
fig.update_layout(title_x=0.49)

# Show and save the plot
fig.show()
fig.write_image('heatmap_tissue2celltype.png')
