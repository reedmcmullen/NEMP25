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

#Read in a csv file defining the sample names and the paths to the species assignments. 
print('Subsetting and saving species assignments data...')
species_assignments_paths_df = pd.read_csv("/wynton/home/pollenlab/reedmcmullen/projects/NEMP25/scanpy_NEMP25/NEMP25_species_assignments_paths.csv", index_col="sample")

# Read in the species assignment data.
species_assignments_df = pd.DataFrame()
for idx, (sample_name, row) in enumerate(species_assignments_paths_df.iterrows()):
    sample_path = row['path']
    df = pd.read_csv(sample_path, sep='\t', index_col=0, header=None, names=['cell_barcode', 'species', 'species_droplet_type', 'species_LLR'])
    df.index = df.index.str.replace(r'NEMP25_GEX_t[1-3]_[1-6]', '', regex=True) + f'{idx+1}'
    species_assignments_df = pd.concat([species_assignments_df, df])

# Subset the species assignments metadata to only the cell barcodes observed in the adata object after preprocessing. 
species_assignments_df=species_assignments_df.loc[species_assignments_df.index.isin(adata.obs.index),:]

# Add the species assignments metadata to the adata.obs dataframe.
adata.obs.loc[species_assignments_df.index,'species']=species_assignments_df['species']
adata.obs.loc[species_assignments_df.index,'species_droplet_type']=species_assignments_df['species_droplet_type']
adata.obs.loc[species_assignments_df.index,'species_LLR']=species_assignments_df['species_LLR']

#Save the AnnData object as an .h5ad file.
anndata_preprocessed = directory_path + '/NEMP25_preprocessed.h5ad'
adata.write(anndata_preprocessed, compression='gzip')

#Plotting
plt.rcParams["figure.figsize"] = (6, 6)
sc.pl.umap(adata, color='species', save=f'_{save_name}_species.png')
sc.pl.umap(adata, color='species_droplet_type', save=f'_{save_name}_species_droplet_type.png')

#Plot a bar plot of the species composition.
plt.rcParams["figure.figsize"] = (5, 5)
data = adata.obs["species"].value_counts(normalize=True, dropna=False)
data = round(data, 3) * 100
colors = adata.uns['species_colors']
ax = data.plot.bar(color=colors)
for container in ax.containers:
    ax.bar_label(container)
ax.grid(False)
plt.ylabel('Percent of Total Cells')
plt.xlabel('Species')
plt.title('Species Composition')
plt.savefig(os.path.join(directory_path, 'figures', f'{save_name}_species_assignment_barchart.png'), bbox_inches='tight')

#Plot a bar plot of the species droplet type composition.
plt.rcParams["figure.figsize"] = (5, 5)
data = adata.obs["species_droplet_type"].value_counts(normalize=True, dropna=False)
data = round(data, 3) * 100
colors = adata.uns['species_droplet_type_colors']
ax = data.plot.bar(color=colors)
for container in ax.containers:
    ax.bar_label(container)
ax.grid(False)
plt.ylabel('Percent of Total Cells')
plt.xlabel('Species Droplet Type Assignment')
plt.title('Species Droplet Type Composition')
plt.savefig(os.path.join(directory_path, 'figures', f'{save_name}_species_droplet_type_barchart.png'), bbox_inches='tight')

#Plot a stacked barchart of the species assignments by GEMwell.
species_data = adata.obs.groupby(['GEMwell', 'species']).size().unstack(fill_value=0)
species_data =species_data.div(species_data.sum(axis=1), axis=0) * 100  # Normalize to percentages
colors = adata.uns['species_assignment_colors']
plt.rcParams["figure.figsize"] = (15, 5)
ax = species_data.plot(kind='bar', stacked=True, color=colors, width=0.8)
ax.grid(False)
plt.ylim(top=100)
plt.ylabel('Percent of Total Cells')
plt.xlabel('GEM well')
plt.title('Species Composition by GEMwell')
plt.legend(title='Species', bbox_to_anchor=(1.05, 1), loc='upper left')  # Adjust legend position
plt.tight_layout()
plt.savefig(os.path.join(directory_path, 'figures', f'{save_name}_species_assignment_stacked_barchart_by_gemwell.png'), bbox_inches='tight')

#Plot a stacked barchart of the species droplet type by GEMwell.
species_data = adata.obs.groupby(['GEMwell', 'species_droplet_type']).size().unstack(fill_value=0)
species_data =species_data.div(species_data.sum(axis=1), axis=0) * 100  # Normalize to percentages
colors = adata.uns['species_droplet_type_colors']
plt.rcParams["figure.figsize"] = (15, 5)
ax = species_data.plot(kind='bar', stacked=True, color=colors, width=0.8)
ax.grid(False)
plt.ylim(top=100)
plt.ylabel('Percent of Total Cells')
plt.xlabel('GEMwell')
plt.title('Species Droplet Type Composition by GEMwell')
plt.legend(title='Species Droplet Type', bbox_to_anchor=(1.05, 1), loc='upper left')  # Adjust legend position
plt.tight_layout()
plt.savefig(os.path.join(directory_path, 'figures', f'{save_name}_species_droplet_type_stacked_barchart_by_gemwell.png'), bbox_inches='tight')

#Plot a stacked barchart of the species assignments by TimePoint.
species_data = adata_species_singlets.obs.groupby(['TimePoint', 'species_assignment']).size().unstack(fill_value=0)
species_data =species_data.div(species_data.sum(axis=1), axis=0) * 100  # Normalize to percentages
colors = adata.uns['species_assignment_colors']
plt.rcParams["figure.figsize"] = (5, 5)
ax = species_data.plot(kind='bar', stacked=True, color=colors, width=0.8)
ax.grid(False)
plt.ylim(top=100)
plt.ylabel('Percent of Total Cells')
plt.xlabel('Time Point')
plt.title('Species Assignment Composition by Time Point')
plt.legend(title='Species Assignment', bbox_to_anchor=(1.05, 1), loc='upper left')  # Adjust legend position
plt.tight_layout()
plt.savefig(os.path.join(directory_path, 'figures', f'{save_name}_species_assignment_stacked_barchart_by_timepoint_singlets.png'), bbox_inches='tight')

