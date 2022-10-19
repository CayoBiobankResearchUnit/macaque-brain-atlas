#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('biccn','rna','atac','3')

prefix = arguments[1]
rn.prefix = arguments[2]
at.prefix = arguments[3]
this.cluster = arguments[4]

suppressMessages(library(monocle3))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))

this.sample = 'all'

meta.file = file.path('stats/glue/subintegration',paste0('meta_',prefix,'_class',this.cluster,'.txt.gz'))
umap02.file = file.path('stats/glue/subintegration',paste0('umap02_',prefix,'_class',this.cluster,'.txt.gz'))
umap10.file = file.path('stats/glue/subintegration',paste0('umap10_',prefix,'_class',this.cluster,'.txt.gz'))

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
cds = cluster_cells(cds, resolution=1e-5, reduction_method='PCA')
# cds = cluster_cells(cds, resolution=1e-5, reduction_method='UMAP')

push.status(paste('cluster_cells',prefix))

this.umap = data.frame(meta,umap02,partition=cds@clusters$PCA$partitions,cluster=cds@clusters$PCA$clusters)
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

dir.create('figures/glue/subintegration',showWarnings=FALSE)

that.umap = readRDS(file.path('umap',paste0(rn.prefix,'-scanpy-recluster-classified.rds')))
that.umap = that.umap[c('cell','cell_subcluster_manual')]

this.umap$index = 1:nrow(this.umap)

this.umap = merge(this.umap,that.umap,by='cell',all.x=TRUE,all.y=FALSE)

this.umap = this.umap[order(this.umap$index),]

this.umap$index = NULL

# Save clusters
saveRDS(this.umap,file=file.path('umap',paste0(prefix,'-glue-recluster-class',this.cluster,'.rds')))

set.seed(42)
this.umap = this.umap[sample(1:nrow(this.umap)),]

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

plot.umap(this.umap,color='modality',file=file.path('figures/glue/subintegration',paste0('umap-',prefix,'-class',this.cluster,'-modality-all.pdf')),color.label='Modality',size=size,alpha=alpha)
plot.umap(this.umap,color='modality',file=file.path('figures/glue/subintegration',paste0('umap-',prefix,'-class',this.cluster,'-modality-all_nolegend.pdf')),color.label='Modality',legend=FALSE,size=size,alpha=alpha)
plot.umap(this.umap,color='modality',file=file.path('figures/glue/subintegration',paste0('umap-',prefix,'-class',this.cluster,'-modality-facet-all.pdf')),color.label='Modality',facet=TRUE,facet.by='modality',legend=FALSE,size=size,alpha=alpha)
plot.umap(this.umap,color='modality',file=file.path('figures/glue/subintegration',paste0('umap-',prefix,'-class',this.cluster,'-modality-all_presentation.pdf')),color.label='Modality',legend=FALSE,presentation=TRUE,size=size,alpha=alpha)

plot.umap(this.umap,color='cell_subcluster_manual',file=file.path('figures/glue/subintegration',paste0('umap-',prefix,'-class',this.cluster,'-cellsubcluster-all.pdf')),color.label='Cell subtype',size=size,alpha=alpha)
plot.umap(this.umap,color='cell_subcluster_manual',file=file.path('figures/glue/subintegration',paste0('umap-',prefix,'-class',this.cluster,'-cellsubcluster-all_nolegend.pdf')),color.label='Cell subtype',legend=FALSE,size=size,alpha=alpha)
plot.umap(this.umap,color='cell_subcluster_manual',file=file.path('figures/glue/subintegration',paste0('umap-',prefix,'-class',this.cluster,'-cellsubcluster-facet-all.pdf')),color.label='Cell subtype',facet=TRUE,facet.by='region',legend=FALSE,size=size,alpha=alpha)
plot.umap(this.umap,color='cell_subcluster_manual',file=file.path('figures/glue/subintegration',paste0('umap-',prefix,'-class',this.cluster,'-cellsubcluster-all_presentation.pdf')),color.label='Cell subtype',legend=FALSE,presentation=TRUE,size=size,alpha=alpha)

plot.umap(this.umap,color='region',file=file.path('figures/glue/subintegration',paste0('umap-',prefix,'-class',this.cluster,'-region-all.pdf')),color.label='Region',size=size,alpha=alpha)
plot.umap(this.umap,color='region',file=file.path('figures/glue/subintegration',paste0('umap-',prefix,'-class',this.cluster,'-region-all_nolegend.pdf')),color.label='Region',legend=FALSE,size=size,alpha=alpha)
plot.umap(this.umap,color='region',file=file.path('figures/glue/subintegration',paste0('umap-',prefix,'-class',this.cluster,'-region-facet-all.pdf')),color.label='Region',facet=TRUE,facet.by='region',legend=FALSE,size=size,alpha=alpha)
plot.umap(this.umap,color='region',file=file.path('figures/glue/subintegration',paste0('umap-',prefix,'-class',this.cluster,'-region-all_presentation.pdf')),color.label='Region',legend=FALSE,presentation=TRUE,size=size,alpha=alpha)

plot.umap(this.umap,color='region_major',file=file.path('figures/glue/subintegration',paste0('umap-',prefix,'-class',this.cluster,'-regionmajor-all.pdf')),color.label='Class',size=size,alpha=alpha)
plot.umap(this.umap,color='region_major',file=file.path('figures/glue/subintegration',paste0('umap-',prefix,'-class',this.cluster,'-regionmajor-all_nolegend.pdf')),color.label='Class',legend=FALSE,size=size,alpha=alpha)
plot.umap(this.umap,color='region_major',file=file.path('figures/glue/subintegration',paste0('umap-',prefix,'-class',this.cluster,'-regionmajor-facet-all.pdf')),color.label='Class',facet=TRUE,facet.by='region_major',legend=FALSE,size=size,alpha=alpha)
plot.umap(this.umap,color='region_major',file=file.path('figures/glue/subintegration',paste0('umap-',prefix,'-class',this.cluster,'-regionmajor-all_presentation.pdf')),color.label='Class',legend=FALSE,presentation=TRUE,size=size,alpha=alpha)
