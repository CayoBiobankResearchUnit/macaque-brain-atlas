#!/usr/bin/env python
# 🐛

import sys
import pandas as pd
import scipy as sp, scipy.io as si
import anndata as ad
import scanpy as sc
import scanpy.external as sce
import numpy as np
import gzip as gz
import os
import matplotlib.pyplot as plt
import math
import pyreadr as pyr

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())
sc._settings.settings._vector_friendly=True

arguments = sys.argv
# arguments = ['scanpy-marker-genes.py','atac-marker-gene']

plot_list = [arguments[i] for i in range(1,len(arguments))]

# Read in main data
os.makedirs('figures/final',exist_ok=True)

if bool(np.isin('atac-marker-gene',plot_list)):
	adata = ad.read_h5ad('hdf5/anndata_'+'atac'+'_all_markergenes.h5ad')
	push_status('read_h5ad')
	
	def plot_gene_activity(gene,colorbar='right'):
		p = sc.pl.umap(adata,gene_symbols='gene_short_name',color=gene,frameon=None,color_map='viridis',return_fig=True,vmin=0,vmax=3,colorbar_loc=colorbar)
		p.axes[0].set_aspect(aspect=1)
		p.axes[0].set_axis_off()
		p.axes[0].set_title('')
		p.axes[0].set_rasterized(True)
		if colorbar == 'right':
			p.savefig('figures/final/atac-geneactivity-umap-'+gene+'.pdf',transparent=True)
		elif colorbar is None:
			p.savefig('figures/final/atac-geneactivity-umap-noscale-'+gene+'.pdf',transparent=True)
	
	for g in ['ENTPD1','GAD2','SLC1A2','ATP10A']:
		print(g)
		plot_gene_activity(g)
		plot_gene_activity(g,colorbar=None)
	
	
	



# def plot_umap(gene,cell_type=None):
# 	if not bool(np.isin(gene,adata.var.gene_short_name.to_list())):
# 		print('Cannot plot '+gene)
# 		return
# 	print(gene)
# 	a = sc.pl.umap(adata, gene_symbols='gene_short_name', color=gene,return_fig=True)
# 	a.axes[0].set_aspect(aspect=1)
# 	if cell_type is None:
# 		cell_suffix=''
# 	else:
# 		cell_suffix='_'+cell_type
# 	a.savefig('figures/marker_genes/umap_'+prefix+'_'+gene+cell_suffix+'.pdf')
# 	a.clear()
# 	plt.close('all')
# 
# # scanpy should be 1.9.1 for colorbar=None to work
# def plot_umaps(genes,cell_type=None):
#     good_genes = []
#     for gene in genes:
#             if bool(np.isin(gene,adata.var.gene_short_name.to_list())): good_genes.append(gene)
#     ncolumns = math.ceil(math.sqrt(len(good_genes)))
#     if ncolumns < 3:
#     	ncolumns = 3
#     a = sc.pl.umap(adata, gene_symbols='gene_short_name', color=good_genes, colorbar_loc=None,frameon=None,ncols=ncolumns,wspace=0,return_fig=True)
#     for i in range(0,len(good_genes)):
#             a.axes[i].set_aspect(aspect=1)
#             a.axes[i].set_axis_off()
#             a.axes[i].title.set_fontsize('24')
#     if cell_type is None:
#             raise ValueError('Parameter cell_type cannot be None')
#     a.savefig('figures/marker_genes/umap_'+prefix+'_'+cell_type+'_all.pdf')
#     a.clear()
#     plt.close('all')
