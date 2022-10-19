#!/usr/bin/env python

import sys
import pandas as pd
import scipy as sp, scipy.io as si
import numpy as np
import gzip as gz
import os
import math
import joblib

arguments = sys.argv
# arguments = ['get-fragments.py','atac','atac','NSM201','12']

prefix = arguments[1]
analysis = arguments[2]
this_sample = int(arguments[3])
n_cores = arguments[4]

# Read in the ID column from the sample metadata data frame
sample_list = pd.read_csv('data/'+prefix+'_metadata.txt',sep='\t')['id'].tolist()

# Calculate the sample ID (NSM number)
query_sample = sample_list[this_sample-1]

fragments = pd.read_csv('fragment-files/'+query_sample+'-fragments.txt.gz',header=None,sep='\t',index_col=None)

celltypes = pd.read_csv('stats/clusters/'+prefix+'-celltype-predictions.txt',header=0,sep='\t',index_col=0)
celltype_levels = pd.read_csv('stats/clusters/'+'rna'+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()

celltypes = celltypes[['id','modality','glue_type','glue_type_confidence']]

celltypes['cell'] = celltypes.index.to_list()

celltypes = celltypes[celltypes['glue_type_confidence'] >= 0.95]

cells_sample = celltypes[(celltypes['id'] == query_sample)]

celltype_levels_indata = list(cells_sample['glue_type'].unique())

# Reorder
celltype_levels_indata = np.array(celltype_levels)[np.isin(celltype_levels,celltype_levels_indata)].tolist()

fragments_pass = fragments.merge(cells_sample[['glue_type']],how='outer',left_on=3,right_index=True)

os.makedirs('fragment-files-cells',exist_ok=True)
def get_cluster(cell):
	this_cluster = celltype_levels.index(cell)+1
	fragment_in = fragments_pass[fragments_pass['glue_type'] == cell]
	fragment_in[[int(i) for i in range(0,5)]].to_csv('fragment-files-cells/'+prefix+'-fragments-'+query_sample+'-class'+str(this_cluster)+'.txt.gz',sep='\t',header=None,index=None)
	return this_cluster

try:
	dev_null = joblib.Parallel(n_jobs=pd.Series([len(celltype_levels_indata),int(n_cores)-1]).min())(joblib.delayed(get_cluster)(i) for i in celltype_levels_indata)
except:
	for cell in celltype_levels_indata:
		this_cluster = celltype_levels.index(cell)+1
		fragment_in = fragments_pass[fragments_pass['glue_type'] == cell]
		fragment_in[[int(i) for i in range(0,5)]].to_csv('fragment-files-cells/'+prefix+'-fragments-'+query_sample+'-class'+str(this_cluster)+'.txt.gz',sep='\t',header=None,index=None)
