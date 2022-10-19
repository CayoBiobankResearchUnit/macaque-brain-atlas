#!/usr/bin/env python

import sys
import pandas as pd
import scipy as sp, scipy.io as si
import numpy as np
import gzip as gz
import os
import math
import joblib
import pysam

arguments = sys.argv
# arguments = ['atac-get-cell-type-transposition-sites.py','atac','atac','1','24']

prefix = arguments[1]
analysis = arguments[2]
this_sample = int(arguments[3])
n_cores = arguments[4]

# Read in the ID column from the sample metadata data frame
sample_list = pd.read_csv('data/'+prefix+'_metadata.txt',sep='\t')['id'].tolist()

# Calculate the sample ID (NSM number)
query_sample = sample_list[this_sample-1]

transposition_sites = pd.read_csv('transposition-sites/'+query_sample+'-transposition_sites.bed.gz',header=None,sep='\t',index_col=None)

celltypes = pd.read_csv('stats/clusters/'+prefix+'-celltype-predictions.txt',header=0,sep='\t',index_col=0)
celltype_levels = pd.read_csv('stats/clusters/'+'rna'+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()

celltypes = celltypes[['id','modality','glue_type','glue_type_confidence']]

celltypes['cell'] = celltypes.index.to_list()

celltypes = celltypes[celltypes['glue_type_confidence'] >= 0.95]

cells_sample = celltypes[(celltypes['id'] == query_sample)]

celltype_levels_indata = list(cells_sample['glue_type'].unique())

# Reorder
celltype_levels_indata = np.array(celltype_levels)[np.isin(celltype_levels,celltype_levels_indata)].tolist()

# Merge annotations and keep only those with annotations
transposition_sites_pass = transposition_sites.merge(cells_sample[['glue_type']],how='outer',left_on=3,right_index=True)
transposition_sites_pass = transposition_sites_pass[~transposition_sites_pass['glue_type'].isnull()]

# Make chromosomes categorical
transposition_sites_pass[0] = pd.Categorical(
	[str(x) for x in transposition_sites_pass[0]],
	categories=[str(x) for x in range(1,21)] + ['X','Y']
)

# Sort bed
transposition_sites_pass = transposition_sites_pass.sort_values(by=[0,1,2]).reset_index(drop=True)

os.makedirs('transposition-sites-cells',exist_ok=True)

for cell in celltype_levels_indata:
	this_cluster = celltype_levels.index(cell)+1
	transposition_site_in = transposition_sites_pass[transposition_sites_pass['glue_type'] == cell]
	transposition_site_in[[int(i) for i in range(0,4)]].to_csv('transposition-sites-cells/'+prefix+'-transposition-sites-'+query_sample+'-class'+str(this_cluster)+'.bed',sep='\t',header=None,index=None)
	pysam.tabix_compress('transposition-sites-cells/'+prefix+'-transposition-sites-'+query_sample+'-class'+str(this_cluster)+'.bed','transposition-sites-cells/'+prefix+'-transposition-sites-'+query_sample+'-class'+str(this_cluster)+'.bed.gz',force=True)
	os.remove('transposition-sites-cells/'+prefix+'-transposition-sites-'+query_sample+'-class'+str(this_cluster)+'.bed')
	pysam.tabix_index('transposition-sites-cells/'+prefix+'-transposition-sites-'+query_sample+'-class'+str(this_cluster)+'.bed.gz',force=True,preset='bed')

# def get_cluster(cell):
# 	this_cluster = celltype_levels.index(cell)+1
# 	transposition_site_in = transposition_sites_pass[transposition_sites_pass['glue_type'] == cell]
# 	transposition_site_in[[int(i) for i in range(0,4)]].to_csv('transposition-sites-cells/'+prefix+'-transposition-sites-'+query_sample+'-class'+str(this_cluster)+'.bed',sep='\t',header=None,index=None)
# 	pysam.tabix_compress('transposition-sites-cells/'+prefix+'-transposition-sites-'+query_sample+'-class'+str(this_cluster)+'.bed','transposition-sites-cells/'+prefix+'-transposition-sites-'+query_sample+'-class'+str(this_cluster)+'.bed.gz',force=True)
# 	pysam.tabix_index('transposition-sites-cells/'+prefix+'-transposition-sites-'+query_sample+'-class'+str(this_cluster)+'.bed.gz',force=True)
# 	return this_cluster
# 
# try:
# 	dev_null = joblib.Parallel(n_jobs=pd.Series([len(celltype_levels_indata),int(n_cores)-1]).min())(joblib.delayed(get_cluster)(i) for i in celltype_levels_indata)
# except:
# 	for cell in celltype_levels_indata:
# 		this_cluster = celltype_levels.index(cell)+1
# 		transposition_site_in = transposition_sites_pass[transposition_sites_pass['glue_type'] == cell]
# 		transposition_site_in[[int(i) for i in range(0,4)]].to_csv('transposition-sites-cells/'+prefix+'-transposition-sites-'+query_sample+'-class'+str(this_cluster)+'.bed',sep='\t',header=None,index=None)
# 		pysam.tabix_compress('transposition-sites-cells/'+prefix+'-transposition-sites-'+query_sample+'-class'+str(this_cluster)+'.bed','transposition-sites-cells/'+prefix+'-transposition-sites-'+query_sample+'-class'+str(this_cluster)+'.bed.gz',force=True)
# 		pysam.tabix_index('transposition-sites-cells/'+prefix+'-transposition-sites-'+query_sample+'-class'+str(this_cluster)+'.bed.gz',force=True,preset='bed')
# 