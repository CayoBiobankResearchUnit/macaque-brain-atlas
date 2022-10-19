#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('atac','atac')

prefix = arguments[1]
analysis = arguments[2]

suppressMessages(library(monocle3))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))

sample.list = as.character(read.delim(file.path('data',paste0(prefix,'_metadata.txt')))$id)

this.sample = 'all'

meta.file = file.path('stats/umap',paste0('meta_',prefix,'_all.txt.gz'))
umap02.file = file.path('stats/umap',paste0('umap02_',prefix,'_all.txt.gz'))
umap10.file = file.path('stats/umap',paste0('umap10_',prefix,'_all.txt.gz'))

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
cds = cluster_cells(cds, resolution=1e-4, reduction_method='PCA')
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

dir.create('figures/umap-all',showWarnings=FALSE)

this.umap$auto_doublet = this.umap$doublet_score > 0.20

cell.names = row.names(this.umap)

# Calculate distribution
# with(this.umap,tapply(auto_doublet,cluster,mean))

cluster.doublets = data.frame(
	cluster = levels(this.umap$cluster),
	doublet_mean = with(this.umap,tapply(auto_doublet,cluster,mean)),
	doublet_cluster = with(this.umap,tapply(manual_doublet | auto_doublet,cluster,mean)) > 0.15
)

this.umap$doublet_cluster = this.umap$cluster %in% subset(cluster.doublets,doublet_cluster)$cluster

# this.umap = merge(this.umap,cluster.doublets,by='cluster',all.x=TRUE,all.y=TRUE)
# this.umap = this.umap[cell.names,]

this.umap$doublet_call = factor(with(this.umap,ifelse(manual_doublet,'manual doublet cell',ifelse(auto_doublet,'auto doublet cell',ifelse(doublet_cluster,'doublet cluster','singlet')))),levels=c('singlet','manual doublet cell','auto doublet cell','doublet cluster'))

# library(parallel)
# this.by.cluster = do.call(rbind,mclapply(split(this.umap,this.umap$cluster),function(x) {
# 	if (mean(x$auto_doublet) > 0.15) {
# 		x$doublet_cluster = TRUE
# 	} else {
# 		x$doublet_cluster = FALSE
# 	}
# 	x
# },mc.cores=n.cores))[cell.names,]

# Save clusters
saveRDS(this.umap,file=file.path('umap',paste0(prefix,'-scanpy-all.rds')))

# Write clusters
dir.create('stats/clusters',showWarnings=FALSE)

write.table(
	data.frame(x=rownames(this.umap),y=as.integer(as.character(this.umap$cluster))-1),
	file=file.path('stats/clusters',paste0(prefix,'-original-clusters.txt')),
	sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE
)
write.table(
	data.frame(x=rownames(this.umap),y=as.integer(as.character(this.umap$partition))-1),
	file=file.path('stats/clusters',paste0(prefix,'-original-partitions.txt')),
	sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE
)

dir.create('stats/cells',showWarnings=FALSE)
write(
	rownames(subset(this.umap,doublet_call == 'singlet')),
	file=file.path('stats/cells',paste0(prefix,'-singlets.txt')),
	sep='\n'
)
write(
	rownames(subset(this.umap,doublet_call != 'singlet')),
	file=file.path('stats/cells',paste0(prefix,'-doublets.txt')),
	sep='\n'
)


if (analysis == 'atac') {
	plot.umap(this.umap,color='umi_binarized',file=file.path('figures/umap-all',paste0('umap-',prefix,'-umi-all.pdf')),color.label='UMI')
	plot.umap(this.umap,color='FRIP',file=file.path('figures/umap-all',paste0('umap-',prefix,'-frip-all.pdf')),color.label='FRIP')
	plot.umap(this.umap,color='FRIT',file=file.path('figures/umap-all',paste0('umap-',prefix,'-frit-all.pdf')),color.label='FRIT')
} else if (analysis == 'rna') {
	plot.umap(this.umap,color='total_counts',file=file.path('figures/umap-all',paste0('umap-',prefix,'-umi-all.pdf')),color.label='UMI')
	plot.umap(this.umap,color='pct_counts_mt',file=file.path('figures/umap-all',paste0('umap-',prefix,'-mt-all.pdf')),color.label='MT (%)')
}

plot.umap(this.umap,color='doublet_score',file=file.path('figures/umap-all',paste0('umap-',prefix,'-doubletscore-all.pdf')),color.label='Doublet Score')

# plot.umap(this.umap,color='cluster',file=file.path('figures/umap-all',paste0('umap-',prefix,'-cluster-all.pdf')),color.label='Cluster')
plot.umap(this.umap,color='partition',file=file.path('figures/umap-all',paste0('umap-',prefix,'-partition-all.pdf')),color.label='Partition')
plot.umap(this.umap,color='region_major',file=file.path('figures/umap-all',paste0('umap-',prefix,'-regionmajor-all.pdf')),color.label='Class')
plot.umap(this.umap,color='region',file=file.path('figures/umap-all',paste0('umap-',prefix,'-region-all.pdf')),color.label='Region')
plot.umap(this.umap,color='id',file=file.path('figures/umap-all',paste0('umap-',prefix,'-sample-all.pdf')),color.label='Sample')
plot.umap(this.umap,color='sex',file=file.path('figures/umap-all',paste0('umap-',prefix,'-sex-all.pdf')),color.label='Sex')
plot.umap(this.umap,color='hemisphere',file=file.path('figures/umap-all',paste0('umap-',prefix,'-hemisphere-all.pdf')),color.label='Hemisphere')
plot.umap(this.umap,color='sequencing_run_id',file=file.path('figures/umap-all',paste0('umap-',prefix,'-sequencingrunid-all.pdf')),color.label='Sequencing Run')
plot.umap(this.umap,color='doublet_call',file=file.path('figures/umap-all',paste0('umap-',prefix,'-doubletcall-all.pdf')),color.label='Doublet Call')

plot.umap(this.umap,color='cluster',file=file.path('figures/umap-all',paste0('umap-',prefix,'-cluster-all_nolegend.pdf')),color.label='Cluster',legend=FALSE)
plot.umap(this.umap,color='partition',file=file.path('figures/umap-all',paste0('umap-',prefix,'-partition-all_nolegend.pdf')),color.label='Partition',legend=FALSE)
plot.umap(this.umap,color='region_major',file=file.path('figures/umap-all',paste0('umap-',prefix,'-regionmajor-all_nolegend.pdf')),color.label='Class',legend=FALSE)
plot.umap(this.umap,color='region',file=file.path('figures/umap-all',paste0('umap-',prefix,'-region-all_nolegend.pdf')),color.label='Region',legend=FALSE)
plot.umap(this.umap,color='id',file=file.path('figures/umap-all',paste0('umap-',prefix,'-sample-all_nolegend.pdf')),color.label='Sample',legend=FALSE)
plot.umap(this.umap,color='sex',file=file.path('figures/umap-all',paste0('umap-',prefix,'-sex-all_nolegend.pdf')),color.label='Sex',legend=FALSE)
plot.umap(this.umap,color='hemisphere',file=file.path('figures/umap-all',paste0('umap-',prefix,'-hemisphere-all_nolegend.pdf')),color.label='Hemisphere',legend=FALSE)
plot.umap(this.umap,color='sequencing_run_id',file=file.path('figures/umap-all',paste0('umap-',prefix,'-sequencingrunid-all_nolegend.pdf')),color.label='Sequencing Run',legend=FALSE)
plot.umap(this.umap,color='doublet_call',file=file.path('figures/umap-all',paste0('umap-',prefix,'-doubletcall-all_nolegend.pdf')),color.label='Doublet Call',legend=FALSE)

# plot.umap(this.umap,color='cluster',file=file.path('figures/umap-all',paste0('umap-',prefix,'-cluster-facet-all.pdf')),color.label='Cluster',facet=TRUE,facet.by='cluster',legend=FALSE)
plot.umap(this.umap,color='partition',file=file.path('figures/umap-all',paste0('umap-',prefix,'-partition-facet-all.pdf')),color.label='Partition',facet=TRUE,facet.by='partition',legend=FALSE)
plot.umap(this.umap,color='region_major',file=file.path('figures/umap-all',paste0('umap-',prefix,'-regionmajor-facet-all.pdf')),color.label='Class',facet=TRUE,facet.by='region_major',legend=FALSE)
plot.umap(this.umap,color='region',file=file.path('figures/umap-all',paste0('umap-',prefix,'-region-facet-all.pdf')),color.label='Region',facet=TRUE,facet.by='region',legend=FALSE)
plot.umap(this.umap,color='id',file=file.path('figures/umap-all',paste0('umap-',prefix,'-sample-facet-all.pdf')),color.label='Sample',facet=TRUE,facet.by='id',legend=FALSE)
plot.umap(this.umap,color='sex',file=file.path('figures/umap-all',paste0('umap-',prefix,'-sex-facet-all.pdf')),color.label='Sex',facet=TRUE,facet.by='sex',legend=FALSE)
plot.umap(this.umap,color='hemisphere',file=file.path('figures/umap-all',paste0('umap-',prefix,'-hemisphere-facet-all.pdf')),color.label='Hemisphere',facet=TRUE,facet.by='hemisphere',legend=FALSE)
plot.umap(this.umap,color='sequencing_run_id',file=file.path('figures/umap-all',paste0('umap-',prefix,'-sequencingrunid-facet-all.pdf')),color.label='Sequencing Run',facet=TRUE,facet.by='sequencing_run_id',legend=FALSE)


push.status(paste('plot.umap',prefix))

# # Sub in old UMAP coordinates
# umap.old.file = file.path('stats/prenorm',paste0('prenorm_umap_',prefix,'_all.txt.gz'))
# umap.old = fread(umap.old.file)
# names(umap.old) = paste0('umap',1:2,sep='.')
# 
# that.umap = this.umap
# that.umap[,paste('umap',1:2,sep='.')] = umap.old
# 
# plot.umap(that.umap,color='cluster',file=file.path('figures/umap-all',paste0('umap_prenorm-',prefix,'-cluster-all.pdf')),color.label='Cluster')
# plot.umap(that.umap,color='partition',file=file.path('figures/umap-all',paste0('umap_prenorm-',prefix,'-partition-all.pdf')),color.label='Partition')
# plot.umap(that.umap,color='region_major',file=file.path('figures/umap-all',paste0('umap_prenorm-',prefix,'-regionmajor-all.pdf')),color.label='Class')


# this.umap$cell_type = factor(this.umap$CellType_1B,
# levels=c('Neurons_Excitatory',
# 'Neurons_Inhibitory_ADARB2',
# 'Neurons_Inhibitory_LHX6/SST',
# 'Neurons_Inhibitory_LHX6/PVALB',
# 'Neurons_MediumSpiny',
# 'Neurons_Thalamic',
# 'Neurons_Purkinje',
# 'Neurons_Basket',
# 'Neurons_Granule',
# 'Oligodendrocyte Precursors',
# 'Astrocytes',
# 'Endothelial',
# 'Microglia',
# 'Oligodendrocytes',
# 'Immune'),
# labels=c('excitatory neurons',
# 'inhibitory neurons (ADARB2)',
# 'inhibitory neurons (LHX6/SST)',
# 'inhibitory neurons (LHX6/PVALB)',
# 'medium spiny neurons',
# 'thalamic neurons',
# 'Purkinje neurons',
# 'basket neurons',
# 'cerebellar granule cells',
# 'oligodendrocyte precursors',
# 'astrocytes',
# 'endothelial cells',
# 'microglia',
# 'oligodendrocytes',
# 'immune cells'))
# 
# 
# 
# plot.umap(this.umap,color='cell_type',file=file.path('figures/umap-final',paste0('umap-',prefix,'-celltype-all.pdf')),color.label='Cell Type',alpha=0.05)
# plot.umap(this.umap,color='cell_type',file=file.path('figures/umap-final',paste0('umap-',prefix,'-celltype-all_nolegend.pdf')),color.label='Cell Type',legend=FALSE,alpha=0.05)
# plot.umap(this.umap,color='cell_type',file=file.path('figures/umap-final',paste0('umap-',prefix,'-celltype-facet-all.pdf')),color.label='Cell Type',facet=TRUE,facet.by='cell_type',legend=FALSE)
# 