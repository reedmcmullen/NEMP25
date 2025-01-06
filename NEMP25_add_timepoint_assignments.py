#! /usr/bin/env python3
#Import required packages and modules.
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import os
import numpy as np

#Define variables and settings.
directory_path = '/wynton/home/pollenlab/reedmcmullen/projects/NEMP25/scanpy_NEMP25'
os.chdir(directory_path)
save_name = 'NEMP25'

#Load the dataset.
print('Loading dataset...')
anndata_preprocessed = directory_path + '/NEMP25_preprocessed.h5ad'
adata = sc.read_h5ad(anndata_preprocessed)

#Add TimePoint metadata.
adata.obs['time_point'] = adata.obs['GEMwell'].apply(lambda x: '6hr' if 1 <= x <= 4 else '54hr' if 5 <= x <= 9 else '7day')

#Plot UMAP colored by time point.
plt.rcParams["figure.figsize"] = (6, 6)
sc.pl.umap(adata, color='time_point', save=f'{save_name}_time_point.png')

# Plot bar plot of time point composition
adata.obs['time_point'] = pd.Categorical(adata.obs['time_point'], categories=['6hr', '54hr', '7day'], ordered=True)
data = adata.obs['time_point'].value_counts(normalize=True).reindex(['6hr', '54hr', '7day'])
data = round(data, 3) * 100
colors = adata.uns['time_point_colors']
plt.rcParams["figure.figsize"] = (5, 5)
ax = data.plot.bar(color=colors)
for container in ax.containers:
    ax.bar_label(container)
ax.grid(False)
plt.ylabel('Percent of Total Cells')
plt.xlabel('Time Point')
plt.title('Overall Time Point Composition')
plt.savefig(os.path.join(directory_path, 'figures',f'{save_name}_time_point_barchart.png'), bbox_inches='tight')

#Save the AnnData object as an .h5ad file.
anndata_preprocessed = directory_path + '/NEMP25_preprocessed.h5ad'
adata.write(anndata_preprocessed, compression='gzip')
