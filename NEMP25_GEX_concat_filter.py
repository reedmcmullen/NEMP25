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
save_name = "_NEMP25"

#Read in the CSV file defining the sample names and the paths to the filtered feature barcode matrices output from 10X Genomics cellranger count pipeline.
matrices = pd.read_csv("/wynton/home/pollenlab/reedmcmullen/projects/NEMP25/scanpy_NEMP25/NEMP25_matrices_paths.csv", index_col="sample")

#Read in the 10X genomics gene expression data as Anndata objects and store in a dictionary with sample_name, sample_adata key, value pairs
#Add the GEMwell metadata to the 'batch' column and make unique indices by appending '-batch' to them.
adata_dict = {}
for idx, (sample_name, row) in enumerate(matrices.iterrows()):
    sample_path = row['path']
    print(f'Reading in data for sample: {sample_name}')
    adata = sc.read_10x_mtx(sample_path, cache=True)
    # Modify the index of adata.obs to append the GEMwell number.
    adata.obs.index = adata.obs.index.str.replace('-1', '', regex=False)
    adata.obs.index = adata.obs.index + f'-{idx+1}'  # Append '-integer' to each index to represent the GEMwell ID.
    adata.obs['GEMwell'] = idx+1 # Add the GEMwell ID to the column 'GEMwell'
    adata.obs['GEMwell'] = adata.obs['GEMwell'].astype('category') # Change the data type of the GEMwell column to be categorical.
    # Store the modified adata in the dictionary
    adata_dict[sample_name] = adata

#Concatenate AnnData objects, taking the union of variables (i.e. 'outer' join).
print('Concatenating AnnData objects')
adata = ad.concat(list(adata_dict.values()), join='outer')

#Identify highly expressed genes.
print('Finding highest expressed genes')
sc.pl.highest_expr_genes(adata, n_top=20)

#Filter out cells based on a very conservative minimum number of genes and cells.
print('QC filtering of genes and cells')
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=10)

#Annotate mitochondrial genes and ribosomal genes and calculate QC metrics for each.
print('QC filtering by QC metrics')
adata.var['mito'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mito'
adata.var['ribo'] = adata.var_names.str.startswith('RPS' or 'RPL') # annotate the group of ribosomal genes as 'ribo'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mito', 'ribo'], percent_top=None, log1p=False, inplace=True)

#Save initial concatenated AnnData object as a H5AD file.
adata.layers['counts'] = adata.X.copy()
print('Saving initial AnnData object')
anndata_initial = directory_path + '/NEMP25_initial.h5ad'
adata.write(anndata_initial)
