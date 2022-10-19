#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('rna','rna','cerebellum')

prefix = arguments[1]
analysis = arguments[2]
this.region = arguments[3]
this.region.column = arguments[4]

suppressMessages(library(monocle3))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))

# This script is for performing the clustering from the scanpy workflow

this.sample = 'all'

meta.file = file.path('stats/umap_recluster',paste0('meta_',prefix,'_region_',this.region,'.txt.gz'))
umap02.file = file.path('stats/umap_recluster',paste0('umap02_',prefix,'_region_',this.region,'.txt.gz'))
umap10.file = file.path('stats/umap_recluster',paste0('umap10_',prefix,'_region_',this.region,'.txt.gz'))

meta = fread(meta.file,na.strings=c('','NA'))
cell.names = rownames(meta) = meta$V1
meta$V1 = NULL

umap02 = fread(umap02.file)
names(umap02) = paste('umap',1:2,sep='.')

umap10 = fread(umap10.file)
names(umap10) = paste('umap',1:10,sep='.')

rownames(umap02) = rownames(umap10) = cell.names

e = t(umap10)
colnames(e) = cell.names
g.df = data.frame(umap=rownames(e))
rownames(g.df) = g.df$umap
g.df$gene_short_name = NA

# Now slot the UMAP into monocle3 (it will go into the expression data slot, which won't be used for downstream work).
cds = new_cell_data_set(
	expression_data=e,
	cell_metadata=meta,
	gene_metadata=g.df)

# We're going to call the 10-dimensional UMAP a "PCA", but it's not really
pca.monocle = as.matrix(umap10)
dimnames(pca.monocle) = list(rownames(umap10),NULL)

umap.monocle = as.matrix(umap02)
dimnames(umap.monocle) = list(rownames(umap02),NULL)

reducedDims(cds) = SimpleList(
	PCA = pca.monocle,
	UMAP = umap.monocle
)

# Cluster based on the PCA slot, but note that the PCA slot contains the 10-D UMAP
broad.clusters = cluster_cells(cds, resolution=1e-6, reduction_method='PCA')@clusters$PCA$clusters
cds = cluster_cells(cds, resolution=1e-5, reduction_method='PCA')
# cds = cluster_cells(cds, resolution=1e-5, reduction_method='UMAP')

push.status(paste0('cluster_cells ',prefix,' region ',this.region))

this.umap = data.frame(meta,umap02,partition2=factor(cds@clusters$PCA$partitions),cluster2=factor(cds@clusters$PCA$clusters),cluster3=broad.clusters)
# rownames(this.umap) = this.umap$cell

this.umap$region_major = factor(c(dmPFC='Cortical',vmPFC='Cortical',dlPFC='Cortical',
vlPFC='Cortical',ACC='Cortical',CC='Subcortical',CN='Subcortical',
NAc='Subcortical',EC='Cortical',PC='Cortical',A1='Cortical',
AMY='Subcortical',HIP='Subcortical',M1='Cortical',mdTN='Subcortical',
vlTN='Subcortical',LGN='Subcortical',S1='Cortical',IPP='Cortical',
SPP='Cortical',STS='Cortical',IT='Cortical',V1='Cortical',
CV='Cerebellum',lCb='Cerebellum',MB='Brainstem',MdO='Brainstem',MdC='Brainstem',
Pons='Brainstem')[this.umap$region],levels=c('Cortical','Subcortical','Cerebellum','Brainstem'))

this.umap$region = factor(this.umap$region,levels=region.levels)[,drop=TRUE]

dir.create('figures/umap-recluster',showWarnings=FALSE)

# prefix = paste0(prefix,'_bbknn')

# Save clusters
saveRDS(this.umap,file=file.path('umap',paste0(prefix,'-scanpy-recluster-region_',this.region,'.rds')))

set.seed(42)
this.umap = this.umap[sample(1:nrow(this.umap)),]

# Adjust size and alpha

if ((1e6/nrow(this.umap)) < 1) {
	size=0.1
	alpha=0.05
} else if ((1e6/nrow(this.umap)) < 10) {
	size=0.25
	alpha=0.2
} else if ((1e6/nrow(this.umap)) < 100) {
	size=0.5
	alpha=0.5
} else {
	size=1
	alpha=0.5
}

plot.umap(this.umap,color='cluster2',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-cluster-region_',this.region,'_presentation.pdf')),color.label='Cluster',legend=FALSE,presentation=TRUE,size=size,alpha=alpha)
plot.umap(this.umap,color='partition2',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-partition-region_',this.region,'_presentation.pdf')),color.label='Partition',legend=FALSE,presentation=TRUE,size=size,alpha=alpha)
plot.umap(this.umap,color='region_major',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-regionmajor-region_',this.region,'_presentation.pdf')),color.label='Class',legend=FALSE,presentation=TRUE,size=size,alpha=alpha)
plot.umap(this.umap,color='region',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-region-region_',this.region,'_presentation.pdf')),color.label='Region',legend=FALSE,presentation=TRUE,size=size,alpha=alpha)
plot.umap(this.umap,color='sequencing_run_id',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-batch-region_',this.region,'_presentation.pdf')),color.label='Region',legend=FALSE,presentation=TRUE,size=size,alpha=alpha)

plot.umap(this.umap,color='cluster2',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-cluster-region_',this.region,'.pdf')),color.label='Cluster',legend=TRUE,presentation=FALSE,size=size,alpha=alpha)
plot.umap(this.umap,color='partition2',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-partition-region_',this.region,'.pdf')),color.label='Partition',legend=TRUE,presentation=FALSE,size=size,alpha=alpha)
plot.umap(this.umap,color='region_major',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-regionmajor-region_',this.region,'.pdf')),color.label='Class',legend=TRUE,presentation=FALSE,size=size,alpha=alpha)
plot.umap(this.umap,color='region',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-region-region_',this.region,'.pdf')),color.label='Region',legend=TRUE,presentation=FALSE,size=size,alpha=alpha)
plot.umap(this.umap,color='sequencing_run_id',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-batch-region_',this.region,'.pdf')),color.label='Region',legend=TRUE,presentation=FALSE,size=size,alpha=alpha)

plot.umap(this.umap,color='cluster3',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-cluster3-region_',this.region,'_presentation.pdf')),color.label='Cluster',legend=FALSE,presentation=TRUE,size=size,alpha=alpha)
plot.umap(this.umap,color='cluster3',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-cluster3-region_',this.region,'.pdf')),color.label='Cluster',legend=TRUE,presentation=FALSE,size=size,alpha=alpha)

push.status(paste0('plot.umap ',prefix,' region ',this.region))
