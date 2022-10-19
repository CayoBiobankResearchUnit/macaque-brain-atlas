#!/usr/bin/env python

import sys
import pandas as pd
import scipy as sp, scipy.io as si
import numpy as np
import gzip as gz
import os
import math

arguments = sys.argv
# arguments = ['get-fragments.py','atac','atac','NSM201']

prefix = arguments[1]
analysis = arguments[2]
this_sample = int(arguments[3])

# Read in the ID column from the sample metadata data frame
sample_list = pd.read_csv('data/'+prefix+'_metadata.txt',sep='\t')['id'].tolist()

# Calculate the sample ID (NSM number)
query_sample = sample_list[this_sample-1]

celltypes = pd.read_csv('stats/clusters/'+prefix+'-celltype-predictions.txt',header=0,sep='\t',index_col=0)
celltype_levels = pd.read_csv('stats/clusters/'+'rna'+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()

celltypes = celltypes[['id','modality','glue_type','glue_type_confidence']]

celltypes['cell'] = celltypes.index.to_list()

celltypes = celltypes[celltypes['glue_type_confidence'] >= 0.95]

cells_sample = celltypes[(celltypes['id'] == query_sample)]

celltype_levels_indata = list(cells_sample['glue_type'].unique())

# Reorder
celltype_levels_indata = np.array(celltype_levels)[np.isin(celltype_levels,celltype_levels_indata)].tolist()

os.makedirs('whitelists-cells',exist_ok=True)
for cell in celltype_levels_indata:
	this_cluster = celltype_levels.index(cell)+1
	cells_sample[cells_sample['glue_type'] == cell]['cell'].to_csv('whitelists-cells/'+prefix+'-whitelists-'+query_sample+'-class'+str(this_cluster)+'.txt',sep='\t',header=None,index=None)
