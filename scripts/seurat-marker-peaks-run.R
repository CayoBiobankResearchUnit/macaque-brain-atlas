#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('atac','atac','1','24')

library(Seurat)
library(future)

prefix = arguments[1]
analysis = arguments[2]
this.i = as.integer(arguments[3])

if (length(arguments) > 3) {
	n.cores = as.integer(arguments[4])
} else {
	n.cores = future::availableCores()
}

cds = readRDS(file.path('rds/seurat',paste0(prefix,'_seurat_cds.rds')))

plan(
	strategy='multicore',
	workers = n.cores
)

i1 = levels(cds@meta.data$cell_class)[this.i]
i2 = setdiff(levels(cds@meta.data$cell_class),i1)

cell.types = scan(file.path('stats/clusters',paste0('rna','-final-cellclasses-levels.txt')),what='',sep='\n',quiet=TRUE)

cell.class = match(i1,cell.types)

da.peaks = FindMarkers(
	object = cds,
	assay = 'ATAC',
	slot = 'counts',
	ident.1 = i1,
	ident.2 = i2,
	logfc.threshold = 0,
	min.pct = 0,
	test.use = 'LR',
	latent.vars = 'n.umi',
	only.pos = FALSE,
	verbose = TRUE
)
da.peaks$cluster = i1
da.peaks$peak = rownames(da.peaks) = gsub('-','_',rownames(da.peaks))

saveRDS(da.peaks,file=file.path('rds',paste0(prefix,'_marker_peaks_seurat_all_class',cell.class,'.rds')))
