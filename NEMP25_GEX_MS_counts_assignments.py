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

#Read in a csv file defining the sample names and the paths to the MULTI-seq barcode counts. 
print('Loading MULTI-seq count data...')
multiseq_counts_paths_df = pd.read_csv("/wynton/home/pollenlab/reedmcmullen/projects/NEMP25/scanpy_NEMP25/NEMP25_multiseq_counts_paths.csv", index_col="sample")

# Read in the MULTI-seq barcode counts data as a dataframe.
multiseq_counts_df = pd.DataFrame()
for idx, (sample_name, row) in enumerate(multiseq_counts_paths_df.iterrows()):
    sample_path = row['path']
    # Read in MULTI-seq metadata.
    df = pd.read_csv(sample_path, sep='\t', index_col=0)
    # Modify the index to match the cell barcode format of adata.obs.
    df.index = df.index.astype(str) + f'-{idx+1}'
    # Append the current dataframe to the multiseq_df.
    multiseq_counts_df = pd.concat([multiseq_counts_df, df])

# Subset the MULTI-seq count data to only the cells observed in the preprocessed AnnData object and write to a csv file.
print('Subsetting and saving MULTI-seq count data...')
multiseq_counts_df=multiseq_counts_df.loc[multiseq_counts_df.index.isin(adata.obs.index),:]
multiseq_counts_df.to_csv(directory_path + '/NEMP25_multiseq_counts_df.csv')

#Read in a csv file defining the sample names and the paths to the MULTI-seq assignment data.
print('Loading MULTI-seq assignment data...')
multiseq_assignments_paths_df = pd.read_csv("/wynton/home/pollenlab/reedmcmullen/projects/NEMP25/scanpy_NEMP25/NEMP25_multiseq_assignments_paths.csv", index_col="sample")

# Read in the MULTI-seq assignment data as a dataframe.
multiseq_df = pd.DataFrame()
for idx, (sample_name, row) in enumerate(multiseq_assignments_paths_df.iterrows()):
    sample_path = row['path']
    # Read in MULTI-seq metadata.
    df = pd.read_csv(sample_path, sep='\t', index_col=0, header=None, names=['cell_barcode', 'perturbation', 'MS_droplet_type', 'MS_LLR'])
    # Modify the index to match the cell barcode format of adata.obs.
    df.index = df.index.astype(str) + f'-{idx+1}'
    # Append the current dataframe to the multiseq_df.
    multiseq_df = pd.concat([multiseq_df, df])

# Subset the MULTI-seq assignment data to only the cells observed in the preprocessed AnnData object and write to a csv file.
print('Subsetting and saving MULTI-seq assignment data...')
multiseq_df=multiseq_df.loc[multiseq_df.index.isin(adata.obs.index),:]
multiseq_df.to_csv(directory_path + '/NEMP25_multiseq_assignments_df.csv')

# Add the MULTI-seq assignments data to the AnnData object.
adata.obs.loc[multiseq_df.index,'perturbation']=multiseq_df['perturbation']
adata.obs.loc[multiseq_df.index,'MS_droplet_type']=multiseq_df['MS_droplet_type']
adata.obs.loc[multiseq_df.index,'MS_LLR']=multiseq_df['MS_LLR']
adata.obs.head(5)

#Save the results.
anndata_preprocessed = directory_path + '/NEMP25_preprocessed.h5ad'
adata.write(anndata_preprocessed, compression='gzip')




