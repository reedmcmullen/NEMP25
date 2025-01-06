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

#Subset the AnnData object to single cells.
print(adata.shape)
adata = adata[adata.obs['MS_droplet_type'] == 'S']
adata = adata[adata.obs['species_droplet_type']=='S']
#adata = adata[adata.obs['individual_droplet_type']=='S']
print(adata.shape)

#Save the AnnData object as an .h5ad file.
anndata = directory_path + '/NEMP25.h5ad'
adata.write(anndata, compression='gzip')
