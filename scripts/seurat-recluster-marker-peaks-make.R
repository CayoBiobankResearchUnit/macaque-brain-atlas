#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('atacsub','subpeak','12')

library(Seurat)
library(Matrix)
library(data.table)

prefix = arguments[1]
glue.prefix = arguments[2]
this.cluster = as.integer(arguments[3])

x = readMM(file.path('mm/recluster/seurat',paste0(prefix,'_downsampled_class',this.cluster,'_counts.mtx')))
cell.metadata = as.data.frame(fread(file.path('mm/recluster/seurat',paste0(prefix,'_downsampled_class',this.cluster,'_cell_metadata.txt.gz')),sep='\t',header=TRUE))
feature.metadata = as.data.frame(fread(file.path('mm/recluster/seurat',paste0(prefix,'_downsampled_class',this.cluster,'_feature_metadata.txt.gz')),sep='\t',header=TRUE))

rownames(x) = rownames(cell.metadata) = cell.metadata$V1
colnames(x) = rownames(feature.metadata) = feature.metadata$V1

cell.metadata$V1 = NULL
feature.metadata$V1 = NULL

cell.subtypes = scan(file.path('stats/subclusters',paste0('rna','-cellsubclusters-levels.txt')),what='',sep='\n',quiet=TRUE)

cell.metadata$cell_subcluster = factor(cell.metadata$cell_subcluster,levels=cell.subtypes)[,drop=TRUE]

cds = CreateSeuratObject(
	counts = t(x),
	meta.data = cell.metadata,
	assay = 'ATAC',
	project = paste0(prefix,'-class',this.cluster)
)

# Write in cell predictions
Idents(cds) = cds@meta.data$cell_subcluster

dir.create('rds/seurat',showWarnings=FALSE)
saveRDS(cds,file=file.path('rds/seurat',paste0(prefix,'_class',this.cluster,'_seurat_cds.rds')))

# da.peaks = FindAllMarkers(
# 	object = cds,
# 	assay = 'ATAC',
# 	slot = 'counts',
# 	logfc.threshold = 0,
# 	min.pct = 0, # Set to 0
# 	test.use = 'LR',
# 	latent.vars = 'n.umi',
# 	only.pos = FALSE,
# 	verbose = TRUE
# )
# 
# rownames(da.peaks) = da.peaks$peak = gsub('-','_',rownames(da.peaks))
# da.peaks$gene = NULL
# 
# saveRDS(da.peaks,file=file.path('rds',paste0(prefix,'_class',this.cluster,'_marker_peaks_seurat.rds')))
# 