#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('atac','atac')

library(Seurat)
library(Matrix)
library(data.table)

prefix = arguments[1]
analysis = arguments[2]

x = readMM(file.path('mm/seurat',paste0(prefix,'_downsampled_counts.mtx')))
cell.metadata = as.data.frame(fread(file.path('mm/seurat',paste0(prefix,'_downsampled_cell_metadata.txt.gz')),sep='\t',header=TRUE))
feature.metadata = as.data.frame(fread(file.path('mm/seurat',paste0(prefix,'_downsampled_feature_metadata.txt.gz')),sep='\t',header=TRUE))

rownames(x) = rownames(cell.metadata) = cell.metadata$V1
colnames(x) = rownames(feature.metadata) = feature.metadata$V1

cell.metadata$V1 = NULL
feature.metadata$V1 = NULL

cell.types = scan(file.path('stats/clusters',paste0('rna','-final-cellclasses-levels.txt')),what='',sep='\n',quiet=TRUE)

cell.metadata$cell_class = factor(cell.metadata$cell_class,levels=cell.types)[,drop=TRUE]

cds = CreateSeuratObject(
	counts = t(x),
	meta.data = cell.metadata,
	assay = 'ATAC',
	project = prefix
)

# Write in cell predictions
Idents(cds) = cds@meta.data$cell_class

dir.create('rds/seurat',showWarnings=FALSE)
saveRDS(cds,file=file.path('rds/seurat',paste0(prefix,'_seurat_cds.rds')))


# Instead of running the below code in serial, use "seurat-marker-peaks-run.R" to run one cell type in parallel

# da.peaks = FindAllMarkers(
# 	object = cds,
# 	assay = 'ATAC',
# 	slot = 'counts',
# 	logfc.threshold = 0,
# 	min.pct = 0.1,
# 	test.use = 'LR',
# 	latent.vars = 'n.umi',
# 	only.pos = FALSE,
# 	verbose = TRUE
# )
# 
# rownames(da.peaks) = da.peaks$peak = gsub('-','_',rownames(da.peaks))
# da.peaks$gene = NULL
# 
# saveRDS(da.peaks,file=file.path('rds',paste0(prefix,'_marker_peaks_seurat.rds')))
# 