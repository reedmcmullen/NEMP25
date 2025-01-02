


#Normalize to median total counts and logarithmize.
print('Normalizing and logarithmizing')
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

#Identify and plot highly variable genes from each AnnData object as an H5AD file.
print('Finding highly variable genes')
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="GEMwell", flavor='seurat')
sc.pl.highly_variable_genes(adata, save=save_name+'.png')

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
results_file_preprocessed = directory_path + '/NEMP25_preprocessed.h5ad'
adata.write(results_file_preprocessed)

#Subset the AnnData object to 10% of cells for faster DEG testing and visualization.
print('Subsetting AnnData object')
# Set a seed for reproducibility.
np.random.seed(42)
# Calculate 10% of the total number of cells
n_cells = int(adata.n_obs * 0.1)
# Randomly sample indices without replacement
random_indices = np.random.choice(adata.obs.index, size=n_cells, replace=False)
# Subset the AnnData object
adata_subset = adata[random_indices].copy()
adata_subset

#Save the subset AnnData object.
print('Saving subset of the preprocessed AnnData object')
results_file_preprocessed_subset = directory_path + '/NEMP25_preprocessed_subset.h5ad'
adata_subset.write(results_file_preprocessed_subset)
