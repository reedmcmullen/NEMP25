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
