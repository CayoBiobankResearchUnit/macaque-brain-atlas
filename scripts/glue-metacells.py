#!/usr/bin/env python
# 🐛

import sys
import pandas as pd
import scipy as sp, scipy.io as si
import anndata as ad
import scanpy as sc
import scanpy.external as sce
import itertools
import numpy as np
import gzip as gz
import os
import networkx as nx
import seaborn as sns
import scglue
import math
from matplotlib import rcParams
import sklearn
import faiss
import re
import pybedtools as bed

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['glue-metacells.py','biccn','rna','atac']

prefix = arguments[1]
rn_prefix = arguments[2]
at_prefix = arguments[3]

adata = ad.read_h5ad('hdf5/anndata_glue_'+prefix+'_combined.h5ad')
push_status(prefix+' read_h5ad')

adata.obsm['X_glue'] = np.loadtxt('stats/glue/glue_'+prefix+'_all.txt.gz',dtype='float32')
push_status(prefix+' loadtxt')

rnadata = ad.read_h5ad('hdf5/anndata_glue_graph_'+rn_prefix+'_all_preprocessed.h5ad')
push_status(prefix+' read_h5ad '+rn_prefix)
atadata = ad.read_h5ad('hdf5/anndata_glue_graph_'+at_prefix+'_all_preprocessed.h5ad')
push_status(prefix+' read_h5ad '+at_prefix)

adata.obs['modality'] = 'RNA'
is_atac = [bool(re.search('^A',x)) for x in adata.obs.biccn_id]
adata.obs.loc[is_atac,'modality'] = 'ATAC'

adata.obs.drop('_index',axis=1,inplace=True)
# rnadata.obs.drop('_index',axis=1,inplace=True)
# atadata.obs.drop('_index',axis=1,inplace=True)

adata.obs.index = adata.obs['cell'].to_list()
# rnadata.obs.index = rnadata.obs['cell'].to_list()
atadata.obs.index = atadata.obs['cell'].to_list()

# Backup
adata.write_h5ad('hdf5/anndata_glue_'+prefix+'_all_pre_metacell.h5ad')
push_status(prefix+' write_h5ad')

rnaInt = adata[adata.obs['modality'] == 'RNA']
ataInt = adata[adata.obs['modality'] == 'ATAC']

del adata

rnadata.obsm['X_glue'] = rnaInt.obsm['X_glue']
atadata.obsm['X_glue'] = ataInt.obsm['X_glue']

rnaobs = rnadata.obs.merge(rnaInt.obs[['cell','modality']],on='cell',how='left')
rnaobs.index = rnaobs['cell'].to_list()

ataobs = atadata.obs.merge(ataInt.obs[['cell','modality']],on='cell',how='left')
ataobs.index = ataobs['cell'].to_list()

rnadata.obs = rnaobs
atadata.obs = ataobs

# Make counts matrix
rnadata.X = rnadata.layers['counts']

# Add cell type metadata
####################

at_types = pd.read_csv('stats/clusters/'+at_prefix+'-celltype-predictions.txt',header=0,sep='\t',index_col=0).rename(columns={'glue_type':'cell_class','glue_type_confidence':'cell_class_confidence'})
rn_types = pd.read_csv('stats/clusters/'+rn_prefix+'-final-cellclasses.txt',header=None,sep='\t',index_col=0).rename(columns={1:'cell_class'})
rn_types.index = rn_types.index.to_list()

at_types = at_types[['cell_class','cell_class_confidence']]
rn_types['cell_class_confidence'] = np.inf

at_types['cell'] = at_types.index.to_list()
rn_types['cell'] = rn_types.index.to_list()

tmp = rnadata.obs.merge(rn_types,how='left',on='cell',copy=True)
tmp.index = tmp['cell'].to_list()
rnadata.obs = tmp

tmp = atadata.obs.merge(at_types,how='left',on='cell',copy=True)
tmp.index = tmp['cell'].to_list()
atadata.obs = tmp

# celltypes = celltypes[celltypes['glue_subtype_confidence'] >= 0.75]

cell_types = pd.read_csv('stats/clusters/'+rn_prefix+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()

# adata = adata[adata.obs['glue_subtype_confidence'] >= 0.75]

rnadata.obs['cell_class'] = pd.Categorical(
	values = rnadata.obs['cell_class'],
	categories = cell_types
).remove_unused_categories()

atadata.obs['cell_class'] = pd.Categorical(
	values = atadata.obs['cell_class'],
	categories = rnadata.obs['cell_class'].cat.categories.to_list()
)

# Because the file sizes are getting huge, dump the progress thus far
rnadata.write_h5ad('hdf5/anndata_glue_'+rn_prefix+'_all_pre_metacell.h5ad')
atadata.write_h5ad('hdf5/anndata_glue_'+at_prefix+'_all_pre_metacell.h5ad')

# metacell_mean_size = 100
metacell_mean_size = 250
n_meta = math.floor(rnadata.n_obs / metacell_mean_size)

exec(open('scripts/_function_rewrites.py').read())

metacell_results = get_metacells(rnadata,atadata,use_rep='X_glue',n_meta=n_meta,common=True,seed=0)

meta = metacell_results['meta']
adatas = metacell_results['adatas']
meta_summary = []

for i in adatas[0].obs.index.to_list():
	print(i)
	
	this = meta[meta['metacell'] == int(i)]
	
	# pd.crosstab(this['modality'],this['cell_subcluster'])
	
	this_rn = this[this['modality'] == 'RNA']
	this_at = this[this['modality'] == 'ATAC']
	
	rn_counts = this_rn.groupby(['cell_class']).size()
	at_counts = this_at.groupby(['cell_class']).size()
	
	rn_prop = rn_counts / rn_counts.sum()
	at_prop = at_counts / at_counts.sum()
	
	rn_call = this['cell_class'].cat.categories[np.argmax(rn_prop)]
	rn_pct = rn_prop.max()
	rn_score = rn_prop.max() - rn_prop[set(rn_prop.index) ^ set([rn_call])].max()
	
	at_call = this['cell_class'].cat.categories[np.argmax(at_prop)]
	at_pct = at_prop.max()
	at_score = at_prop.max() - at_prop[set(at_prop.index) ^ set([at_call])].max()
	
	m = pd.DataFrame([{
		'metacell_id':str(i),
		'n_rna':this.groupby('modality').size()['RNA'],
		'n_atac':this.groupby('modality').size()['ATAC'],
		'rn_type':rn_call,
		'rn_pct':rn_pct,
		'rn_score':rn_score,
		'at_type':at_call,
		'at_pct':at_pct,
		'at_score':at_score,
		'call_match':rn_call==at_call
	}])
	m.index = m['metacell_id'].to_list()
	
	meta_summary.append(m)

meta_all = pd.concat(meta_summary)

rnameta = adatas[0]
atameta = adatas[1]

rnameta.obs = meta_all
atameta.obs = meta_all

os.makedirs('stats/metacells',exist_ok=True)
os.makedirs('hdf5/metacells',exist_ok=True)

rnameta.write_h5ad('hdf5/metacells/anndata_metacells_'+rn_prefix+'.h5ad')
atameta.write_h5ad('hdf5/metacells/anndata_metacells_'+at_prefix+'.h5ad')

meta_all.to_csv('stats/metacells/metacells_metadata.txt.gz',sep='\t',header=True,index=True)
meta.to_csv('stats/metacells/metacells_metacell_assignments.txt.gz',sep='\t',header=True,index=True)

si.mmwrite('stats/metacells/metacells_'+rn_prefix+'_counts',rnameta.X,field='integer')
si.mmwrite('stats/metacells/metacells_'+at_prefix+'_counts',atameta.X,field='integer')

rnameta.var.to_csv('stats/metacells/metacells_'+rn_prefix+'_var.txt.gz',sep='\t',header=True,index=True)
atameta.var.to_csv('stats/metacells/metacells_'+at_prefix+'_var.txt.gz',sep='\t',header=True,index=True)


# rnameta = ad.read_h5ad('hdf5/metacells/anndata_metacells_'+rn_prefix+'.h5ad')
# atameta = ad.read_h5ad('hdf5/metacells/anndata_metacells_'+at_prefix+'.h5ad')


# Subset peaks to those within 150k of any genes

# genes = scglue.genomics.Bed(rnameta.var.assign(name=rnameta.var_names).query('highly_variable'))
# peaks = scglue.genomics.Bed(atameta.var.assign(name=atameta.var.index.to_list()))
# 
# tss = genes.strand_specific_start_site()
# promoters = tss.expand(2000, 0)


# Bring in highly variable genes from regulatory analysis
rnameta_var = pd.read_csv('stats/glue/'+prefix+'_'+rn_prefix+'_regulatory_var.txt.gz',header=0,sep='\t',index_col=0)

genes = scglue.genomics.Bed(rnameta_var.assign(name=rnameta_var.index.to_list()).query('highly_variable'))
peaks = scglue.genomics.Bed(atameta.var.assign(name=atameta.var.index.to_list()))

# Bedtools requires sorted coordinates
genes = genes.sort_values(by=['chrom','chromStart','chromEnd'])
peaks = peaks.sort_values(by=['chrom','chromStart','chromEnd'])

window_size = 150000

# Write bed files from gene and peaks data frames
genes.df.to_csv('stats/metacells/metacells_'+prefix+'_genes.bed',sep='\t',header=None,index=None)
peaks.df.to_csv('stats/metacells/metacells_'+prefix+'_peaks.bed',sep='\t',header=None,index=None)

genes = bed.BedTool('stats/metacells/metacells_'+prefix+'_genes.bed')
peaks = bed.BedTool('stats/metacells/metacells_'+prefix+'_peaks.bed')

# # Remove all peaks that overlap genes (skip this)
# intergenic_peaks = peaks.subtract(genes)

# Calculate distance between intergenic peaks and their closest genes
nearby = peaks.closest(genes, d=True)

# Append peaks to list only if they are within window_size of a gene
peaks_pass = []
for peak in nearby:
	if int(peak[-1]) < window_size:
		peaks_pass.append(peak.name)

# Subset ATAC metacells to include only peaks within window_size of genes
atameta_pass = atameta[:,np.isin(atameta.var_names.to_list(),peaks_pass).tolist()]

si.mmwrite('stats/metacells/metacells_'+at_prefix+'_counts_pass',atameta_pass.X,field='integer')
atameta_pass.var.to_csv('stats/metacells/metacells_'+at_prefix+'_var_pass.txt.gz',sep='\t',header=True,index=True)

atameta_pass.write_h5ad('hdf5/metacells/anndata_metacells_'+at_prefix+'_pass.h5ad')


# atameta = ad.read_h5ad('hdf5/metacells/anndata_metacells_'+at_prefix+'.h5ad')
# push_status('read_h5ad')

# Bring in regulatory scores
regulatory_scores = pd.read_csv('stats/glue/'+prefix+'_regulatory_regulatory.txt.gz',sep='\t',index_col=None,header=0)

peaks_regulatory = regulatory_scores['target'].unique()
peaks_metacell = atameta.var_names.to_list()

# Below is super slow
# peaks_keep = [i for i in range(0,len(peaks_metacell)) if peaks_metacell[i] in peaks_regulatory]

# atameta_intersect = atameta[:,np.isin(atameta.var_names.to_list(),peaks_intersect).tolist()]

# peaks_intersect = list(set(atameta.var_names.to_list()) & set(peaks_regulatory))
# [i for i in peaks_metacell if i in peaks_regulatory]

peaks_regulatory_df = pd.DataFrame(peaks_regulatory,columns=['peak'],index=peaks_regulatory)
peaks_regulatory_df['regulatory'] = True

metacell_var = atameta.var.merge(peaks_regulatory_df['regulatory'],how='outer',left_index=True,right_index=True)
atameta_intersect = atameta[:,metacell_var['regulatory'] == True]

si.mmwrite('stats/metacells/metacells_'+at_prefix+'_counts_intersect',atameta_intersect.X,field='integer')
atameta_intersect.var.to_csv('stats/metacells/metacells_'+at_prefix+'_var_intersect.txt.gz',sep='\t',header=True,index=True)


# Split into 4 chunks since the matrix is way too big
# In R:
# x=c(223616942,196197964,185288947,169963040,187317192,179085566,169868564,145679320,134124166,99517758,133066086,130043856,108737130,128056306,113283604,79627064,95433459,74474043,58315233,77137495,153388924,11753682); y = cumsum(x) / sum(x); unlist(lapply(seq(0.25,1,0.25),function(i) which(y >= i)[1]))






# chroms = {
# 	'set1':['1','2','3','4'],
# 	'set2':['5','6','7','8'],
# 	'set3':['9','10','11','12','13','14'],
# 	'set4':['15','16','17','18','19','20','X','Y']
# }
# 
# for c in chroms:
# 	print(c)
# 	atameta_this = atameta[:,np.isin(atameta.var['chrom'].to_list(),chroms[c]).tolist()]
# 	si.mmwrite('stats/metacells/metacells_'+at_prefix+'_counts_pass_'+c,atameta_this.X,field='integer')
# 	atameta_this.var.to_csv('stats/metacells/metacells_'+at_prefix+'_var_pass_'+c+'.txt.gz',sep='\t',header=True,index=True)
# 