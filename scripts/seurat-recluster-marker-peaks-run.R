#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('atacsub','subpeak','1','21','24')

library(Seurat)
library(Matrix)
library(data.table)

prefix = arguments[1]
glue.prefix = arguments[2]
this.cluster = as.integer(arguments[3])
this.i = as.integer(arguments[4])

if (length(arguments) > 4) {
	n.cores = as.integer(arguments[5])
	plan(
		strategy='multicore',
		workers = n.cores
	)
} 

cds = readRDS(file.path('rds/seurat',paste0(prefix,'_class',this.cluster,'_seurat_cds.rds')))

i1 = levels(cds@meta.data$cell_subcluster)[this.i]
i2 = setdiff(levels(cds@meta.data$cell_subcluster),i1)

if (nlevels(cds@meta.data$cell_subcluster) == 1) {
	da.peaks = data.frame(
		p_val = NA,
		avg_log2FC = NA,
		pct.1 = rowMeans(cds@assays$ATAC@counts),
		pct.2 = NA,
		p_val_adj = NA,
		cluster = i1,
		peak = gsub('-','_',rownames(cds@assays$ATAC@counts))
	)
	rownames(da.peaks) = rownames(cds@assays$ATAC@counts)
	
	da.peaks = da.peaks[order(da.peaks$pct.1,decreasing=TRUE),]
} else if (sum(cds@meta.data$cell_subcluster == i1) < 3) {
	da.peaks = data.frame(
		p_val = NA,
		avg_log2FC = NA,
		pct.1 = rowMeans(matrix(cds@assays$ATAC@counts[,cds@meta.data$cell_subcluster == i1],ncol=sum(cds@meta.data$cell_subcluster == i1))),
		pct.2 = rowMeans(cds@assays$ATAC@counts[,cds@meta.data$cell_subcluster != i1]),
		p_val_adj = NA,
		cluster = i1,
		peak = gsub('-','_',rownames(cds@assays$ATAC@counts))
	)
	rownames(da.peaks) = rownames(cds@assays$ATAC@counts)
	
	da.peaks = da.peaks[order(da.peaks$pct.1/da.peaks$pct.2,da.peaks$pct.1,decreasing=TRUE),]
} else {
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
}

saveRDS(da.peaks,file=file.path('rds',paste0(prefix,'_class',this.cluster,'_marker_peaks_seurat_subcluster_i',this.i,'.rds')))
