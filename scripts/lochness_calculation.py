import sys
import pandas as pd
import scipy as sp, scipy.io as si
import scanpy as sc
import anndata as ad
import scanpy.external as sce
import numpy as np
import gzip as gz
import scrublet as scr
import os
import math
import episcanpy.api as epi
import sklearn as sk
import sklearn.decomposition  # TruncatedSVD
import sklearn.preprocessing  # normalize
import sklearn.feature_extraction.text  # TfidfTransformer
import csv
import seaborn as sns

cell_meta_all = pd.read_csv('../data/rna-scanpy-final-classified-meta.txt', sep='\t',engine="python")
knn_conn = np.load('../data/rna_knn_conn.npy',allow_pickle=True)
knn_dist = np.load('../data/rna_knn_dist.npy',allow_pickle=True)
knn_conn = knn_conn.item()
knn_dist = knn_dist.item()
region_counts = cell_meta_all['region'].value_counts()
regions = cell_meta_all['region'].unique()
target_score = []
kadj = (knn_conn>0).sum(axis=1)

for r in regions:
    target_mask = ((cell_meta_all['region'] == r) *1).to_numpy()
    target_neighbors = knn_conn.dot(target_mask)
    target_score.append(np.multiply(target_neighbors,np.squeeze(np.reciprocal(kadj*1.0)))/region_counts[r]*len(target_neighbors))

target_score_mtx = np.stack(target_score, axis=0)

means = []
for r in regions:
    tmp = target_score_mtx[:,np.where(cell_meta_all['region'] == r)].reshape(len(regions),-1)
    means.append(np.mean(tmp,axis=1).T)

means_mtx = np.stack(means, axis=0)

ct_dict = {1:"excitatory neurons",
2:"medium spiny neurons",
3:"inhibitory neurons",
4:"dopaminergic neurons",
5:"serotonergic neurons",
6:"cerebellar neurons",
7:"basket cells",
8:"astrocytes",
9:"oligodendrocytes",
10:"oligodendrocyte precursor cells",
11:"vascular cells",
12:"microglia"}

for ct in range(12):
    text_file = open('../data/knn_graphs/cellids_rna_class'+str(ct+1)+'.txt', 'r')
    cell_ids = text_file.read().split()
    cell_meta_sub = cell_meta_all.loc[cell_ids]
    knn_conn_sub = np.load('../data/knn_graphs/connectivities_rna_class'+str(ct+1)+'.npy',allow_pickle=True).item()
    region_counts = cell_meta_sub['region'].value_counts()
    regions = cell_meta_sub['region'].unique()
    target_score = []
    kadj = (knn_conn_sub>0).sum(axis=1)
    for r in regions:
        target_mask = ((cell_meta_sub['region'] == r) * 1).to_numpy()
        target_neighbors = knn_conn_sub.dot(target_mask)
        target_score.append(
            np.multiply(target_neighbors, np.squeeze(np.reciprocal(kadj * 1.0))) / region_counts[r] * len(target_neighbors))
    target_score_mtx = np.stack(target_score, axis=0)
    np.savetxt('../similarity_score/'+ct_dict[ct+1]+'_similarity_scores.txt', target_score_mtx.T, delimiter=',')
    np.savetxt('../similarity_score/'+ct_dict[ct+1]+'_similarity_scores_regions.txt', regions, delimiter=',',fmt='%s')
    means = []
    for r in regions:
        tmp = target_score_mtx[:, np.where(cell_meta_sub['region'] == r)].reshape(len(regions), -1)
        means.append(np.mean(tmp, axis=1).T)
    means_mtx = np.stack(means, axis=0)
    np.savetxt('../similarity_score/'+ct_dict[ct+1]+'_similarity_mtx.txt', means_mtx, delimiter=',')
