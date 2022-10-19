#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('biccn')

prefix = arguments[1]

suppressMessages(library(monocle3))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))

this.sample = 'all'

meta.file = file.path('stats/glue',paste0('meta_',prefix,'_all.txt.gz'))
umap02.file = file.path('stats/glue',paste0('umap02_',prefix,'_all.txt.gz'))
umap10.file = file.path('stats/glue',paste0('umap10_',prefix,'_all.txt.gz'))

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

dir.create('figures/umap-glue',showWarnings=FALSE)

# Save clusters
saveRDS(this.umap,file=file.path('umap',paste0(prefix,'-glue-final.rds')))

set.seed(42)
this.umap = this.umap[sample(1:nrow(this.umap)),]

plot.umap(this.umap,color='modality',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-modality-all.pdf')),color.label='Modality')
plot.umap(this.umap,color='modality',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-modality-all_nolegend.pdf')),color.label='Modality',legend=FALSE)
plot.umap(this.umap,color='modality',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-modality-facet-all.pdf')),color.label='Modality',facet=TRUE,facet.by='modality',legend=FALSE)
plot.umap(this.umap,color='modality',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-modality-all_presentation.pdf')),color.label='Modality',legend=FALSE,presentation=TRUE)

plot.umap(this.umap,color='cluster',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-cluster-all.pdf')),color.label='Cluster')
plot.umap(this.umap,color='cluster',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-cluster-all_nolegend.pdf')),color.label='Cluster',legend=FALSE)
plot.umap(this.umap,color='cluster',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-cluster-facet-all.pdf')),color.label='Cluster',facet=TRUE,facet.by='cluster',legend=FALSE)
plot.umap(this.umap,color='cluster',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-cluster-all_presentation.pdf')),color.label='Cluster',legend=FALSE,presentation=TRUE)

plot.umap(this.umap,color='partition',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-partition-all.pdf')),color.label='Partition')
plot.umap(this.umap,color='partition',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-partition-all_nolegend.pdf')),color.label='Partition',legend=FALSE)
plot.umap(this.umap,color='partition',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-partition-facet-all.pdf')),color.label='Partition',facet=TRUE,facet.by='partition',legend=FALSE)
plot.umap(this.umap,color='partition',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-partition-all_presentation.pdf')),color.label='Partition',legend=FALSE,presentation=TRUE)

plot.umap(this.umap,color='region',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-region-all.pdf')),color.label='Region')
plot.umap(this.umap,color='region',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-region-all_nolegend.pdf')),color.label='Region',legend=FALSE)
plot.umap(this.umap,color='region',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-region-facet-all.pdf')),color.label='Region',facet=TRUE,facet.by='region',legend=FALSE)
plot.umap(this.umap,color='region',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-region-all_presentation.pdf')),color.label='Region',legend=FALSE,presentation=TRUE)

plot.umap(this.umap,color='region_major',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-regionmajor-all.pdf')),color.label='Class')
plot.umap(this.umap,color='region_major',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-regionmajor-all_nolegend.pdf')),color.label='Class',legend=FALSE)
plot.umap(this.umap,color='region_major',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-regionmajor-facet-all.pdf')),color.label='Class',facet=TRUE,facet.by='region_major',legend=FALSE)
plot.umap(this.umap,color='region_major',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-regionmajor-all_presentation.pdf')),color.label='Class',legend=FALSE,presentation=TRUE)


# Bring in RNA and ATAC annotations

# Temp file (this.umap with cell annotations)
rna.umap = readRDS('rna_cells_assigned.rds')
cell.assign = rna.umap$cell_final
names(cell.assign) = rna.umap$cell

this.umap$cell_final = NA
this.umap$cell_final = cell.assign[this.umap$cell]

plot.umap(this.umap,color='cell_final',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-cellfinal-all.pdf')),color.label='Cell Type')
plot.umap(this.umap,color='cell_final',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-cellfinal-all_nolegend.pdf')),color.label='Cell Type',legend=FALSE)
plot.umap(this.umap,color='cell_final',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-cellfinal-facet-all.pdf')),color.label='Cell Type',facet=TRUE,facet.by='cell_final')
plot.umap(this.umap,color='cell_final',file=file.path('figures/umap-glue',paste0('umap-',prefix,'-cellfinal-all_presentation.pdf')),color.label='Cell Type',legend=FALSE,presentation=TRUE)
