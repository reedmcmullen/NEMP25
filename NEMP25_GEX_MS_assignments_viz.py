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

#Plot a UMAP polot of MULTI-seq droplet type.
plt.rcParams["figure.figsize"] = (6, 6)
sc.pl.umap(adata, color='MS_droplet_type', save=f'_{save_name}_MS_droplet_type.png')

#Plot a bar plot of the MULTI-seq droplet type. 
plt.rcParams["figure.figsize"] = (5, 5)
data = adata.obs["MS_droplet_type"].value_counts(normalize=True, dropna=True)
data = round(data, 3) * 100
colors = adata.uns['MS_droplet_type_colors']
ax = data.plot.bar(color=colors)
for container in ax.containers:
    ax.bar_label(container)
ax.grid(False)  # Remove gridlines
plt.ylim(top=70)
plt.ylabel('Percent of Total Cells')
plt.xlabel('MULTI-seq Droplet Type Assignment')
plt.title('Overall Droplet Type Composition')
plt.savefig(os.path.join(directory_path, 'figures',f'{save_name}_MS_droplet_type_barchart.png'), bbox_inches='tight')
plt.show()
