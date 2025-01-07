#! /usr/bin/env python3
#Import required packages and modules.
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import os
import numpy as np
import random

#Define variables and settings.
directory_path = '/wynton/home/pollenlab/reedmcmullen/projects/NEMP25/scanpy_NEMP25'
os.chdir(directory_path)
save_name = 'NEMP25'

#Load the dataset.
print('Loading dataset...')
anndata = directory_path + '/NEMP25.h5ad'
adata = sc.read_h5ad(anndata)

#Subset and save dataset
print('Subsetting and saving dataset...')
random.seed(0)
adata_sub=adata[np.random.choice(adata.obs.index, 100_000, replace=False),:]
anndata_subset = directory_path + '/NEMP25_subset.h5ad'
adata_sub.write(anndata_subset, compression='gzip')
