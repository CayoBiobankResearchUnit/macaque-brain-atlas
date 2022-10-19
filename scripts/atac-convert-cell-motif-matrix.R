#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('atac','atac')

prefix = arguments[1]
analysis = 'atac'

library(Matrix)
library(data.table)

mtx = readMM(file.path('stats/motifs',paste0(prefix,'_cell_by_motif-all.mtx.gz')))

rows = fread(file.path('stats/motifs',paste0(prefix,'_cell_by_motif-all-rows.txt.gz')),header=FALSE)[[1]]
cols = fread(file.path('stats/motifs',paste0(prefix,'_cell_by_motif-all-columns.txt.gz')),header=FALSE)[[1]]

rownames(mtx) = rows
colnames(mtx) = cols

m = as(t(mtx),'dgCMatrix')

saveRDS(m,file=file.path('stats/motifs',paste0(prefix,'_cell_by_motif-all.rds')))
