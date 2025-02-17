#! /usr/bin/env python3
#Read in required packages and modules.
import os
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import pandas as pd
import pertpy as pt

#Define and set the current working directory.
directory_path = '/wynton/group/pollen/reedmcmullen/projects/NEMP25/pertpy_NEMP25'
os.chdir(directory_path)

#Load dataset.
print('Loading dataset...')
anndata = '/wynton/group/pollen/reedmcmullen/projects/NEMP25/scanpy_NEMP25/NEMP25_Tel.h5ad'
adata = sc.read_h5ad(anndata)

#Subset dataset.
print('Subsetting and formatting dataset...')
adata = adata[adata.obs['cell_type'].isin(['Tel_NE', 'Tel_NE_ventral', 'Tel_NE_dorsal'])]
adata = adata[adata.obs['species'].isin(['Human', 'Chimp'])]

#Format dataset.
adata.obs['perturbation'] = adata.obs['perturbation'].str.replace(r'\(.*?\)', '', regex=True) #Remove parentheses and any characters between them.
adata.obs['perturbation'] = adata.obs['perturbation'].str.replace('5%', '', regex=False) #Remove 5% from 5%TreMan value of perturbation column.
adata.obs['group'] = adata.obs['species'].astype(str) +'_'+ adata.obs['time_point'].astype(str) +'_'+ adata.obs['perturbation'].astype(str)
adata.obs['simple_group'] = adata.obs['species'].astype(str) +'_'+ adata.obs['perturbation'].astype(str)
adata.obs['group'] = adata.obs['group'].astype(str)
adata.obs['simple_group'] = adata.obs['simple_group'].astype(str)
adata.obs['cell_type'] = adata.obs['cell_type'].astype(str)

#Calculate energy distance between perturbations.
print('Calculating and plotting energy distance between perturbations...')
distance = pt.tl.Distance("edistance", obsm_key="X_pca")
pert_edist_df = distance.pairwise(adata, groupby='perturbation')
pert_edist_df.to_csv(directory_path + '/pert_edist_df.csv')

#Plot energy distance between perturbations.
clustermap(pert_edist_df, robust=True, figsize=(20, 20))
plt.savefig(os.path.join(directory_path, 'perturbation_edist_clustermap.png'), bbox_inches='tight')

#Calculate energy distance between simple groups.
print('Calculating and plotting energy distance between simple groups...')
distance = pt.tl.Distance("edistance", obsm_key="X_pca")
simple_group_edist_df = distance.pairwise(adata, groupby='simple_group')
simple_group_edist_df.to_csv(directory_path + '/simple_group_edist_df.csv')

#Plot energy distance between groups.
clustermap(simple_group_edist_df, robust=True, figsize=(20, 20))
plt.savefig(os.path.join(directory_path, 'simple_group_edist_clustermap.png'), bbox_inches='tight')

#Calculate energy distance between groups.
print('Calculating and plotting energy distance between groups...')
distance = pt.tl.Distance("edistance", obsm_key="X_pca")
group_edist_df = distance.pairwise(adata, groupby='group')
group_edist_df.to_csv(directory_path + '/group_edist_df.csv')

#Plot energy distance between groups.
clustermap(group_edist_df, robust=True, figsize=(20, 20))
plt.savefig(os.path.join(directory_path, 'group_edist_clustermap.png'), bbox_inches='tight')
