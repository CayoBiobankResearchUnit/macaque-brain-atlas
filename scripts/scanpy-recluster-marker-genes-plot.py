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
import pathlib

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())
sc._settings.settings._vector_friendly=True

arguments = sys.argv
# arguments = ['scanpy-preprocess.py','rna','rna','3','t','SST,PVALB,PAX6,VIP,LAMP5,ADARB2']

prefix = arguments[1]
analysis = arguments[2]
this_cluster = int(arguments[3]) - 1
plot_list = arguments[4]
gene_list = arguments[5]

# Read in main data
adata = ad.read_h5ad('hdf5/recluster/anndata_'+prefix+'_class'+str(this_cluster+1)+'_markergenes.h5ad')
push_status(prefix+' read_h5ad')

# # https://github.com/scverse/scanpy/issues/2239#issuecomment-1104178881
adata.uns['log1p']['base'] = None

# Bring in ortholog symbols (genes with no symbol but which have 1:1 human orthologs with symbols)
ortholog_symbols = pd.read_csv('biomart/'+'rna'+'-orthologs-with-symbols.txt',delimiter='\t',header=None,index_col=None)
# ortholog_symbols = pd.read_csv('biomart/'+prefix+'-orthologs-with-symbols.txt',delimiter='\t',header=None,sep='\t',index_col=None)

for i in range(0,len(ortholog_symbols[0])):
	gene = ortholog_symbols[0][i]
	symbol = ortholog_symbols[1][i] + '*'
	print('Renaming '+gene+' to '+symbol)
	adata.var['gene_short_name'].cat.rename_categories({gene:symbol},inplace=True)


if this_cluster == 2:
	adata.uns['cell_subcluster_colors'] = ['#9bfffa','#ffcce7','#9bd352','#dd82c8','#4242d6','#3ac9a8','#b7f28a','#34d822','#9bff84','#e87896','#2c57ba','#8f6dce','#d83227','#4d9901','#66a3cc','#e26f31','#51e8dd','#a5ef8d','#27f468','#67e573']

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
os.makedirs('figures/marker_genes_recluster',exist_ok=True)

do_plots = plot_list.split(',')

plot_dendrogram = 'dendrogram_cell_subcluster' in list(adata.uns.keys())
v_min = None
v_max = 3

if 'u' in do_plots:
	if not 'scaled' in list(adata.layers.keys()):
		adata.layers['scaled'] = np.load('npy/scaled_'+prefix+'_class'+str(this_cluster+1)+'_markergenes.npy')
	
	# scanpy should be 1.9.1 for colorbar=None to work
	def plot_umaps(genes,cell_type=None):
		good_genes = []
		for gene in genes:
				if bool(np.isin(gene,adata.var.gene_short_name.to_list())): good_genes.append(gene)
		ncolumns = math.ceil(math.sqrt(len(good_genes)))
		if ncolumns < 3:
			ncolumns = 3
		# a = sc.pl.umap(adata, gene_symbols='gene_short_name', color=good_genes, colorbar_loc=None,frameon=None,ncols=ncolumns,wspace=0,return_fig=True)
		a = sc.pl.umap(adata, gene_symbols='gene_short_name', color=good_genes, colorbar_loc=None,frameon=None,ncols=ncolumns,wspace=0,return_fig=True,layer='scaled',vmin=v_min,vmax=v_max)
		for i in range(0,len(good_genes)):
				a.axes[i].set_aspect(aspect=1)
				a.axes[i].set_axis_off()
				a.axes[i].title.set_fontsize('24')
		if cell_type is None:
				raise ValueError('Parameter cell_type cannot be None')
		a.savefig('figures/marker_genes_recluster/umap_'+prefix+'_'+cell_type+'.pdf')
		a.clear()
		plt.close('all')
	
	plot_umaps(gene_list.split(','),'class'+str(this_cluster+1)+'_'+gene_list.replace(',','-').replace('*',''))

if 't' in do_plots:
	sc.pl.tracksplot(adata,var_names=gene_list.split(','),groupby='cell_subcluster',dendrogram=plot_dendrogram,gene_symbols='gene_short_name',layer='raw',save=False,xticklabels=False)
	
	
	# sc.pl.tracksplot(adata,var_names=gene_list.split(','),groupby='cell_subcluster',dendrogram=plot_dendrogram,gene_symbols='gene_short_name',layer='raw',save=False,xticklabels=False,cmap=[(0.607843137254902,1,0.980392156862745),(0.56078431372549,0.427450980392157,0.807843137254902),(0.301960784313725,0.6,0.00392156862745098),(0.4,0.63921568627451,0.8),(0.886274509803922,0.435294117647059,0.192156862745098),(0.317647058823529,0.909803921568627,0.866666666666667),(0.647058823529412,0.937254901960784,0.552941176470588),(0.152941176470588,0.956862745098039,0.407843137254902),(0.403921568627451,0.898039215686275,0.450980392156863),(1,0.8,0.905882352941176),(0.607843137254902,0.827450980392157,0.32156862745098),(0.866666666666667,0.509803921568627,0.784313725490196),(0.258823529411765,0.258823529411765,0.83921568627451),(0.227450980392157,0.788235294117647,0.658823529411765),(0.717647058823529,0.949019607843137,0.541176470588235),(0.203921568627451,0.847058823529412,0.133333333333333),(0.607843137254902,1,0.517647058823529),(0.909803921568627,0.470588235294118,0.588235294117647),(0.172549019607843,0.341176470588235,0.729411764705882),(0.847058823529412,0.196078431372549,0.152941176470588)])
	# sc.pl.tracksplot(adata,var_names=gene_list.split(','),groupby='cell_subcluster',dendrogram=plot_dendrogram,gene_symbols='gene_short_name',layer='raw',save=False,xticklabels=False,cmap=["#9bfffa","#8f6dce","#4d9901","#66a3cc","#e26f31","#51e8dd","#a5ef8d","#27f468","#67e573","#ffcce7","#9bd352","#dd82c8","#4242d6","#3ac9a8","#b7f28a","#34d822","#9bff84","#e87896","#2c57ba","#d83227"])
	# sc.pl.tracksplot(adata,var_names=gene_list.split(','),groupby='cell_subcluster',dendrogram=plot_dendrogram,gene_symbols='gene_short_name',save=False,xticklabels=False)
	# sc.pl.rank_genes_groups_tracksplot(adata,var_names=gene_list.split(','),groupby='cell_subcluster',dendrogram=plot_dendrogram,gene_symbols='gene_short_name',save=False,xticklabels=False)
	
	# plt.savefig('figures/marker_genes_recluster/tracks_'+prefix+'_'+'class'+str(this_cluster+1)+'_'+gene_list.replace(',','-').replace('*','')+'.pdf',bbox_inches='tight')
	# plt.savefig('gaba.pdf',bbox_inches='tight')
	# 
	plt.savefig('figures/marker_genes_recluster/tracks_'+prefix+'_'+'class'+str(this_cluster+1)+'_'+gene_list.replace(',','-').replace('*','')+'.png',bbox_inches='tight',dpi=400)

if 'v' in do_plots:
	sc.pl.stacked_violin(adata,var_names=gene_list.split(','),groupby='cell_subcluster',dendrogram=plot_dendrogram,gene_symbols='gene_short_name',cmap='viridis_r',save=False)
	
	plt.savefig('figures/marker_genes_recluster/stackedviolin_'+prefix+'_'+'class'+str(this_cluster+1)+'_'+gene_list.replace(',','-').replace('*','')+'.pdf',bbox_inches='tight')

if 'h' in do_plots:
	if not 'scaled' in list(adata.layers.keys()):
		adata.layers['scaled'] = np.load('npy/scaled_'+prefix+'_class'+str(this_cluster+1)+'_markergenes.npy')
	
	sc.pl.heatmap(adata,var_names=gene_list.split(','),groupby='cell_subcluster',dendrogram=plot_dendrogram,layer='scaled',swap_axes=True,cmap='RdBu_r',vmin=-v_max,vmax=v_max,gene_symbols='gene_short_name',save=False)
	
	plt.savefig('figures/marker_genes_recluster/heatmap_'+prefix+'_'+'class'+str(this_cluster+1)+'_'+gene_list.replace(',','-').replace('*','')+'.pdf',bbox_inches='tight')

if 'd' in do_plots:
	sc.pl.rank_genes_groups_dotplot(adata,var_names=gene_list.split(','),groupby='cell_subcluster',gene_symbols='gene_short_name',vmin=None,vmax=None,save=False)
	
	plt.savefig('figures/marker_genes_recluster/dotplot_'+prefix+'_'+'class'+str(this_cluster+1)+'_'+gene_list.replace(',','-').replace('*','')+'.pdf',bbox_inches='tight')

if 'U' in do_plots:
	# a = sc.pl.umap(adata,color='cell_subcluster',legend_loc='on data',frameon=False,legend_fontsize=6,legend_fontoutline=1,return_fig=True)
	a = sc.pl.umap(adata,color='n.umi',frameon=False,return_fig=True)
	a.axes[0].set_aspect(aspect=1)
	a.axes[0].set_axis_off()
	a.axes[0].title.set_text(None)
	# a.savefig('figures/marker_genes_recluster/umap_'+prefix+'_'+'class'+str(this_cluster+1)+'.pdf')
	a.savefig('figures/marker_genes_recluster/umap_'+prefix+'_'+'class'+str(this_cluster+1)+'_umi.pdf')

