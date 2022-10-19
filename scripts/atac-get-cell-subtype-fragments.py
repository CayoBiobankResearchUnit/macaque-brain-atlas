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
# arguments = ['get-fragments.py','atac','atac','6','12']

prefix = arguments[1]
analysis = arguments[2]
this_cluster = int(arguments[3])
n_cores = arguments[4]

# Read in the ID column from the sample metadata data frame
sample_list = pd.read_csv('data/'+prefix+'_metadata.txt',sep='\t')['id'].tolist()

# # Calculate the sample ID (NSM number)
# query_sample = sample_list[this_sample-1]

cell_subclusters = pd.read_csv('stats/clusters/subintegration/'+'subpeak'+'-'+prefix+'sub'+'-class'+str(this_cluster)+'-cellsubtype-predictions.txt',header=0,index_col=0,sep='\t')

cell_subclusters = cell_subclusters[['id','modality','glue_subtype','glue_subtype_confidence']]

fragments_out = []
for i in [sample_list[x] for x in range(0,110)]:
	print(i)
	if os.path.exists('fragment-files-cells/'+prefix+'-fragments-'+i+'-class'+str(this_cluster)+'.txt.gz'):
		fragments_out.append(pd.read_csv('fragment-files-cells/'+prefix+'-fragments-'+i+'-class'+str(this_cluster)+'.txt.gz',header=None,sep='\t',index_col=None))

fragments = pd.concat(fragments_out)

fragments = fragments.merge(cell_subclusters[['glue_subtype','glue_subtype_confidence']],left_on=3,right_index=True,how='left')

fragments_pass = fragments[fragments['glue_subtype_confidence'] >= 0]

celltype_levels = pd.read_csv('stats/clusters/'+'rna'+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()
this_cell = celltype_levels[this_cluster-1]

# fragments_pass['subtype_number'] = [i.replace(this_cell+' ','') for i in fragments_pass['glue_subtype'].astype('str')]

os.makedirs('fragment-files-clusters',exist_ok=True)

# for i in [this_cell +' '+str(i) for i in range(1,16)]:
for i in list(fragments_pass['glue_subtype'].unique()):
	print(i)
	this_integer = i.replace(this_cell+' ','')
	fragment_in = fragments_pass[fragments_pass['glue_subtype'] == i]
	fragment_in[[int(i) for i in range(0,5)]].to_csv('fragment-files-clusters/'+prefix+'-fragments-class'+str(this_cluster)+'_cluster'+str(this_integer)+'.txt.gz',sep='\t',header=None,index=None)
