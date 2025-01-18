#!/usr/bin/env python3

# plot heatmap of query cell tissue and hit cell tissue
import pandas as pd
import plotly.express as px

df = pd.read_csv('knn_results.tsv', sep='\t', usecols=[1,3])
print(df)

# do heatmap of tissue pair counts
df['count'] = 1
df = df.groupby(['query_tissue', 'hit_tissue']).count().reset_index()
df = df.pivot_table(index='query_tissue', columns='hit_tissue', values='count', fill_value=0)
# title case index and columns
df.index = df.index.str.title()
df.columns = df.columns.str.title()
print(df)


fig = px.imshow(df.values, 
                labels=dict(x="Nearest Neighbor Cell Tissue", y="Query Cell Tissue", color="Count"),
                x=df.columns,
                y=df.index,
                title = 'Tissue of query and nearest neighbor',
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