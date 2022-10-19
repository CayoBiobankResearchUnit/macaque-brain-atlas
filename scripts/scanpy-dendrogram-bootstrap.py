#!/usr/bin/env python
# 🐛

import sys
import scipy
import scipy.cluster.hierarchy
import anndata as ad
import scanpy as sc
import numpy as np
import os

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())
sc._settings.settings._vector_friendly=True

arguments = sys.argv
# arguments = ['scanpy-marker-genes.py','rna','rna','100','0']

prefix = arguments[1]
analysis = arguments[2]

if len(arguments) > 3:
	chunk_size = int(arguments[3])
else:
	chunk_size = 10

if len(arguments) > 4:
	offset = int(arguments[4])
else:
	offset = 0

# Read in main data
adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_slim_dendrogram.h5ad')

# Function for pulling a scipy linkage object and generating a newick file
def export_newick(anndata, dendrogram_key, groupby, savetree=None):
	# Credit: https://stackoverflow.com/a/31878514
	
	def get_newick(node, parent_dist, leaf_names, newick='') -> str:
		"""
		Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.
		
		:param node: output of sciply.cluster.hierarchy.to_tree()
		:param parent_dist: output of sciply.cluster.hierarchy.to_tree().dist
		:param leaf_names: list of leaf names
		:param newick: leave empty, this variable is used in recursion.
		:returns: tree in Newick format
		"""
		if node.is_leaf():
			return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
		else:
			if len(newick) > 0:
				newick = "):%.2f%s" % (parent_dist - node.dist, newick)
			else:
				newick = ");"
			newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
			newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
			newick = "(%s" % (newick)
			return newick
	
	tree = scipy.cluster.hierarchy.to_tree(anndata.uns[dendrogram_key]['linkage'],False)
	leaf_names = anndata.obs[groupby].cat.categories.to_list()
	
	if savetree is None:
		return get_newick(tree, tree.dist, leaf_names)
	else:
		with open(savetree, 'w') as f:
			f.write(get_newick(tree, tree.dist, leaf_names)+'\n')


os.makedirs('stats/dendrograms',exist_ok=True)

# export_newick(adata,'dendrogram_pca_cell_class_keep','cell_class_keep','stats/dendrograms/'+'dendrogram_pca_cellclasskeep.nwk')
# export_newick(adata,'dendrogram_pca_region','region','stats/dendrograms/'+'dendrogram_pca_region.nwk')
# export_newick(adata,'dendrogram_pca_cell_region','cell_region_keep')

for i in range(offset+1,offset+chunk_size+1):
	print(i)
	bdata = adata[np.random.choice(adata.n_obs, size=adata.n_obs, replace=True)]
	sc.tl.dendrogram(bdata,groupby='cell_class_keep',use_rep='X_pca_harmony',key_added='this_class')
	sc.tl.dendrogram(bdata,groupby='cell_class_keep',use_rep='X_umap10',key_added='this_class_umap')
	sc.tl.dendrogram(bdata,groupby='region',use_rep='X_pca_harmony',key_added='this_region')
	export_newick(bdata,'this_class','cell_class_keep','stats/dendrograms/'+'dendrogram_pca_cellclasskeep_bootstrap_'+str(i).zfill(6)+'.nwk')
	export_newick(bdata,'this_class_umap','cell_class_keep','stats/dendrograms/'+'dendrogram_umap_cellclasskeep_bootstrap_'+str(i).zfill(6)+'.nwk')
	export_newick(bdata,'this_region','region','stats/dendrograms/'+'dendrogram_pca_region_bootstrap_'+str(i).zfill(6)+'.nwk')

