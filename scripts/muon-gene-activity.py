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
import muon as mu
import gffpandas.gffpandas as gffpd

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['scanpy-gene-activity.py','atac','atac','1']

prefix = arguments[1]
analysis = arguments[2]

# Read in the ID column from the sample metadata data frame
sample_list = pd.read_csv('data/'+prefix+'_metadata.txt',sep='\t').id.tolist()

# Calculate the sample ID (NSM number)
query_sample = sample_list[int(arguments[3])-1]

adata = ad.read_h5ad('hdf5/anndata_'+prefix+'-'+query_sample+'.h5ad')

fragment_file = 'fragment-files/'+query_sample+'-fragments.txt.gz'

adata.obs.index = adata.obs['cell'].to_list()

if pd.Series(adata.obs['cell']).str.startswith('NSM').any():
	sys.stderr.write('Cell names properly formatted.\n')
else:
	sys.stderr.write('Cell names not properly formatted. Appending sample ID.\n')
	adata.obs['cell'] = np.char.add(np.char.add(adata.obs['id'].to_list(),'-'),adata.obs['cell'].to_list())

cells_final = pd.read_csv('stats/umap_post/meta_'+prefix+'_all.txt.gz',header=0,sep='\t',index_col=0)
vars_final = pd.read_csv('stats/var_final/var_'+prefix+'_all.txt.gz',header=0,sep='\t',index_col=0)

cells_in = np.isin(adata.obs['cell'].to_list(),cells_final.index.to_list()).tolist()
genes_in = np.isin(adata.var.index.to_list(),vars_final.index.to_list()).tolist()

adata = adata[cells_in]
# adata = adata[:,genes_in]

annotation = gffpd.read_gff3(gff_file).filter_feature_of_type(['gene']).attributes_to_columns()

annotation['seq_id'] = annotation['seq_id'].astype('str').to_list()

annotation['Chromosome'] = annotation['seq_id']
annotation['Start'] = annotation['start']
annotation['End'] = annotation['end']

# annotation = annotation[np.isin(annotation['Chromosome'].to_list(),list(vars_final.chromosome.unique().astype('str'))).tolist()]

# Get chromosome lengths
fai = pd.read_csv(fai_file,header=None,sep='\t',index_col=None)

mu.atac.tl.locate_fragments(adata,fragments=fragment_file)

# Drop chromosomes that are unknown
chr_in_data = adata.var_names.str.split(r'[_]').map(lambda x: x[0]).unique()

annotation = annotation[np.isin(annotation['Chromosome'].to_list(),list(chr_in_data)).tolist()]

upstream_bp = 2000

# The behavior of mu.atac.tl.count_fragments_features is unexpected. See issue below.
# https://github.com/scverse/muon/issues/59
# Split by strand to work around.

# Start with the easy stuff

fwd_annotation = annotation[annotation['strand'] == '+']
rev_annotation = annotation[annotation['strand'] == '-']

# forward-strand genes are straightforward

# muon.atac.tl.count_fragments_features fails if the Start - upstream_bp becomes negative.
# To solve this, split the annotation into safe set (pass) and problem set (fail) and process separately

fwd_annotation_pass = fwd_annotation[fwd_annotation['Start'] > upstream_bp]
fwd_annotation_fail = fwd_annotation[fwd_annotation['Start'] <= upstream_bp]

adatas = [mu.atac.tl.count_fragments_features(adata,fwd_annotation_pass,extend_upstream=upstream_bp,extend_downstream=0)]

# For the "fail" set, set the upstream length to the start position
for bp in list(fwd_annotation_fail.Start.unique()):
	this = fwd_annotation_fail[fwd_annotation_fail.Start == bp]
	out = mu.atac.tl.count_fragments_features(adata,this,extend_upstream=bp,extend_downstream=0)
	adatas.append(out)

# The negative strand fragments can be processed all together. Use "extend_downstream" because the upstream region is End + upstream_bp
out = mu.atac.tl.count_fragments_features(adata,rev_annotation,extend_upstream=0,extend_downstream=upstream_bp)
adatas.append(out)

# for c in chr_in_data:
# 	print(c)
# 	chr_len = int(fai[fai[0]==c][1])
# 	chr_annotation = rev_annotation[rev_annotation['Chromosome'] == c]
# 	rev_annotation_pass = chr_annotation[chr_annotation['End'] <= (chr_len - upstream_bp)]
# 	rev_annotation_fail = chr_annotation[chr_annotation['End'] > (chr_len - upstream_bp)]
# 	
# 	out = mu.atac.tl.count_fragments_features(adata,chr_annotation,extend_upstream=0,extend_downstream=upstream_bp)
# 	adatas.append(out)
# 	
# 	for e in list(rev_annotation_fail.End.unique()):
# 		bp = chr_len - e
# 		this = rev_annotation_fail[rev_annotation_fail.End == e]
# 		out = mu.atac.tl.count_fragments_features(adata,this,extend_upstream=0,extend_downstream=bp)
# 		adatas.append(out)

bdata = ad.concat(adatas,axis=1,join='inner',merge='same')

bdata.var.index = bdata.var.gene_id.to_list()

# Do the following to solve bug "TypeError: Can't implicitly convert non-string objects to strings"
# bdata.var.seq_id = bdata.var.seq_id.astype('int64').to_list()
# bdata.var.Chromosome = bdata.var.Chromosome.astype('int64').to_list()

bdata.write_h5ad('hdf5/anndata_'+prefix+'-'+query_sample+'_geneactivity.h5ad')

push_status(prefix+' geneactivity ' + query_sample)
