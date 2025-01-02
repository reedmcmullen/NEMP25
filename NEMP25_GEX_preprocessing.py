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

#Load dataset.
print('Loading the dataset...')
anndata_initial = directory_path + '/NEMP25_initial.h5ad'
adata = sc.read_h5ad(anndata_initial)

#Define QC cutoffs.
n_genes_cutoff = 9000
total_counts_cutoff = 30000
mito_cutoff = 20
ribo_cutoff = 10

#Plot violin plots of QC metrics with cutoffs.
fig, axes = plt.subplots(1, 4, figsize=(16, 6))
sc.pl.violin(adata, 'n_genes_by_counts', ax=axes[0], stripplot=False, show=False)
axes[0].axhline(n_genes_cutoff, color='red', linestyle='--')
sc.pl.violin(adata, 'total_counts', ax=axes[1], stripplot=False, show=False)
axes[1].axhline(total_counts_cutoff, color='red', linestyle='--')
sc.pl.violin(adata, 'pct_counts_mito', ax=axes[2], stripplot=False, show=False)
axes[2].axhline(mito_cutoff, color='red', linestyle='--')
sc.pl.violin(adata, 'pct_counts_ribo', ax=axes[3], stripplot=False, show=False)
axes[3].axhline(ribo_cutoff, color='red', linestyle='--')
plt.savefig(os.path.join(directory_path, 'figures', f'_{save_name}_qc_metrics_cutoffs.png'), bbox_inches='tight')

#Filter outlier cells from each AnnData object based on QC metric plots.
adata = adata[adata.obs.n_genes_by_counts < n_genes_cutoff, :]
adata = adata[adata.obs.total_counts < total_counts_cutoff, :]
adata = adata[adata.obs.pct_counts_mito < mito_cutoff, :]
adata = adata[adata.obs.pct_counts_ribo < ribo_cutoff, :]

#Normalize to median total counts and logarithmize.
print('Normalizing and logarithmizing')
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

#Identify and plot highly variable genes from each AnnData object as an H5AD file.
print('Finding highly variable genes')
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="GEMwell", flavor='seurat')
sc.pl.highly_variable_genes(adata, save=f'{save_name}.png')

#Perform PCA dimensional reduction.
print('Running PCA')
sc.tl.pca(adata, svd_solver='arpack')

#Compute the nearest neighbors using the knn algorithm.
print('Running neighbor finding')
sc.pp.neighbors(adata)

#Compute the umap embedding based on knn-computed neighbors.
print('Running UMAP')
sc.tl.umap(adata)

# Cluster umap embeddings using leiden
print('Leiden clustering')
res = 1
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=res)

#Save the AnnData object as an H5AD file.
print('Saving preprocessed AnnData object')
anndata_preprocessed = directory_path + '/NEMP25_preprocessed.h5ad'
adata.write(anndata_preprocessed, compression='gzip')
