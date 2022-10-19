#!/usr/bin/env python
# 🐛

import sys
import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
import os

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['scanpy-prep.py','atac','atac']

prefix = arguments[1]
analysis = arguments[2]

# Read in the ID column from the sample metadata data frame
sample_list = pd.read_csv('data/'+prefix+'_metadata.txt',sep='\t').id.tolist()

# Read in updated scrublet scores
scrublet_scores = pd.read_csv('data/'+prefix+'_scrublet_parameters_final.txt',sep='\t')

all_h5ad = []

for i in sample_list:
	sys.stderr.write('Now appending '+str(i)+'\n')
	this = ad.read_h5ad('hdf5/anndata_'+prefix+'-'+i+'.h5ad')
	if 'meta_features' in locals():
		new_features = this.var.loc[list(set.difference(set(this.var.index.tolist()),set(meta_features.index.tolist()))),:]
		meta_features = pd.concat([meta_features,new_features],axis=0)
	else:
		meta_features = this.var
	this.obs.insert(this.obs.columns.get_loc('predicted_doublet')+1,'manual_doublet',(this.obs.doublet_score > scrublet_scores[scrublet_scores['id']==i]['final_threshold'].squeeze()).tolist()) 
	all_h5ad.append(this)

push_status(prefix+' append')

# for i in range(0,len(all_h5ad)):
# 	all_h5ad[i].X = all_h5ad[i].X.astype(np.int64)

adata = ad.concat(all_h5ad,join='outer')

push_status(prefix+' concat')

adata.var = meta_features.loc[adata.var.index.tolist(),:]

adata.write_h5ad(filename='hdf5/anndata_'+prefix+'_all.h5ad')

push_status(prefix+' write_h5ad')
