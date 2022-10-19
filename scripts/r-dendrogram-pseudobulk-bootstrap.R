#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('rna','0')

prefix = arguments[1]
i = as.integer(arguments[2])

library(dendextend)
library(ape)
library(phangorn)
library(phylogram)
library(parallel)
library(data.table)

in.file = file.path('expr_bootstrap',paste0(prefix,'_cellclass_pseudobulk_',formatC(i,width=4,flag=0),'.txt.gz'))

genes.info = readRDS('checkpoints/mmulatta_genes.rds')

genes = scan(gzfile(file.path('expr_bootstrap',paste0(prefix,'_cellclass_pseudobulk_genes.txt.gz'))),what='',sep='\n',quiet=TRUE)

genes.keep = which(genes %in% subset(genes.info,chromosome_name %in% 1:20)$ensembl_gene_id)

cells = scan(gzfile(file.path('expr_bootstrap',paste0(prefix,'_cellclass_pseudobulk_cells.txt.gz'))),what='',sep='\n',quiet=TRUE)

x = suppressWarnings(t(as.matrix(fread(in.file)))[genes.keep,])
dimnames(x) = list(genes[genes.keep],cells)

x.norm = limma::voom(x)$E

# x.pca = prcomp(cor(x.norm))$rotation

x.tree = as.phylo(hclust(dist(t(x.norm))))

dir.create('expr_bootstrap/nwk',showWarnings=FALSE)
write.tree(x.tree,file=file.path('expr_bootstrap/nwk',paste0(prefix,'_cellclass_pseudobulk_',formatC(i,width=4,flag=0),'.nwk')))
