#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('rna','rna')

prefix = arguments[1]
analysis = arguments[2]

suppressMessages(library(monocle3))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))

sample.list = as.character(read.delim(file.path('data',paste0(prefix,'_metadata.txt')))$id)

this.sample = 'all'

meta.file = file.path('stats/umap_post',paste0('meta_',prefix,'_all_endogenous.txt.gz'))
umap02.file = file.path('stats/umap_post',paste0('umap02_',prefix,'_all_endogenous.txt.gz'))
umap10.file = file.path('stats/umap_post',paste0('umap10_',prefix,'_all_endogenous.txt.gz'))

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

these.partitions = cds@clusters$PCA$partitions

partitions.join = reshape2::melt(these.partitions)
partitions.join = data.frame(cell=rownames(partitions.join),rna_new_partition=partitions.join$value)



these.clusters = cds@clusters$PCA$clusters

clusters.join = reshape2::melt(these.clusters)
clusters.join = data.frame(cell=rownames(clusters.join),rna_new_cluster=clusters.join$value)



# Bring in metadata and update some things
this.umap = readRDS(file.path('umap',paste0('biccn','-glue-final-all-metadata.rds')))
	
new.umap = fread('stats/umap_post/umap02_rna_all_endogenous.txt.gz')
new.umap.labels = fread('stats/umap_post/meta_rna_all_endogenous.txt.gz')

new.umap = as.data.frame(new.umap)
rownames(new.umap) = new.umap.labels$V1
new.umap$cell = rownames(new.umap)
names(new.umap) = c('cell')

names(this.umap)[names(this.umap) %in% c('rna_umap.1','rna_umap.2')] = c('old_rna_umap.1','old_rna_umap.2')
names(new.umap) = c('rna_umap.1','rna_umap.2','cell')

this.umap$order = 1:nrow(this.umap)
this.umap = merge(this.umap,new.umap,by='cell',all.x=TRUE)

this.umap = this.umap[order(this.umap$order),]

this.umap$region = factor(this.umap$region,levels=region.levels)[,drop=TRUE]

set.seed(42)
this.umap$order = sample(1:nrow(this.umap))

this.umap$rna_region = this.umap$atac_region = this.umap$region
this.umap$rna_region[this.umap$modality == 'ATAC'] = NA
this.umap$atac_region[this.umap$modality == 'RNA'] = NA

this.umap$cell_class = this.umap$rna_cell = this.umap$atac_cell = factor(this.umap$cell_class,levels=cell.levels[!cell.levels %in% c('radial glial cells','mesenchymal stem cells')])
this.umap$rna_cell[this.umap$modality == 'ATAC'] = NA
this.umap$atac_cell[this.umap$modality == 'RNA'] = NA
	
saveRDS(this.umap,file=file.path('umap',paste0('biccn','-glue-all-updated-metadata.rds')))

# this.umap = readRDS(file.path('umap',paste0('biccn','-glue-all-updated-metadata.rds')))

# this.umap$rna_new_partition = NA
# this.umap$rna_new_partition = reshape2::melt(these.partitions)[this.umap$cell,]

this.umap = merge(this.umap,partitions.join,by='cell',all.x=TRUE)
this.umap = merge(this.umap,clusters.join,by='cell',all.x=TRUE)
this.umap = this.umap[order(this.umap$order),]

this.umap$rna_new_cell_class = NA

partition.to.cell = c(
	'1' = 'excitatory neurons',
	'2' = 'oligodendrocytes',
	'3' = 'cerebellar neurons',
	'4' = 'astrocytes',
	'5' = NA,
	'6' = 'oligodendrocyte precursor cells',
	'7' = 'basket cells',
	'8' = 'microglia',
	'9' = 'medium spiny neurons',
	'10' = 'vascular cells',
	'11' = 'ependymal cells',
	'12' = 'unknown excitatory neurons',
	'13' = 'dopaminergic neurons',
	'14' = 'AHSG neurons',
	'15' = 'F5 neurons',
	'16' = 'KIRDL12 neurons'
)

cluster.to.cell = c(
	'11' = 'inhibitory neurons',
	'16' = 'inhibitory neurons',
	'18' = 'inhibitory neurons',
	'22' = 'inhibitory neurons',
	'23' = 'inhibitory neurons',
	'28' = 'inhibitory neurons',
	'31' = 'inhibitory neurons',
	'33' = 'inhibitory neurons',
	'60' = 'serotonergic neurons'
)

this.umap$rna_new_cell_class = partition.to.cell[as.character(this.umap$rna_new_partition)]

this.umap$rna_new_cell_class = with(this.umap,ifelse(!is.na(rna_new_cell_class),rna_new_cell_class,cluster.to.cell[as.character(rna_new_cluster)]))

this.umap$cell_class_old = this.umap$cell_class

this.umap$cell_class = factor(this.umap$cell_class,levels=levels(this.umap$cell_class)[!levels(this.umap$cell_class) %in% c('radial glial cells','mesenchymal stem cells')])

this.umap$keep = with(this.umap,!is.na(rna_cell) | !is.na(atac_region))

# Bring in uncorrected UMAP

nobatch_umap = fread(paste0('stats/umap_post/umap02_','rna','_all_endogenous_uncorrected.txt.gz'))

nobatch_umap$cell = new.umap.labels$V1
names(nobatch_umap)[1:2] = c('rna_nobatch.1','rna_nobatch.2')

this.umap = merge(this.umap,nobatch_umap,by='cell',all.x=TRUE)

this.umap = this.umap[this.umap$order,]
rownames(this.umap) = this.umap$cell

saveRDS(this.umap,file=file.path('umap',paste0('biccn','-metadata-patched.rds')))