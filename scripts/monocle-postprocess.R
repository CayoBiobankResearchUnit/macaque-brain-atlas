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

meta.file = file.path('stats/umap_post',paste0('meta_',prefix,'_all.txt.gz'))
umap02.file = file.path('stats/umap_post',paste0('umap02_',prefix,'_all.txt.gz'))
umap10.file = file.path('stats/umap_post',paste0('umap10_',prefix,'_all.txt.gz'))

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

dir.create('figures/umap-final',showWarnings=FALSE)

# Save clusters
saveRDS(this.umap,file=file.path('umap',paste0(prefix,'-scanpy-final.rds')))

set.seed(42)
this.umap = this.umap[sample(1:nrow(this.umap)),]

if (analysis == 'atac') {
	plot.umap(this.umap,color='umi_binarized',file=file.path('figures/umap-final',paste0('umap-',prefix,'-umi-all.pdf')),color.label='UMI')
	plot.umap(this.umap,color='FRIP',file=file.path('figures/umap-final',paste0('umap-',prefix,'-frip-all.pdf')),color.label='FRIP')
	plot.umap(this.umap,color='FRIT',file=file.path('figures/umap-final',paste0('umap-',prefix,'-frit-all.pdf')),color.label='FRIT')
} else if (analysis == 'rna') {
	plot.umap(this.umap,color='total_counts',file=file.path('figures/umap-final',paste0('umap-',prefix,'-umi-all.pdf')),color.label='UMI')
	plot.umap(this.umap,color='pct_counts_mt',file=file.path('figures/umap-final',paste0('umap-',prefix,'-mt-all.pdf')),color.label='MT (%)')
}

plot.umap(this.umap,color='doublet_score',file=file.path('figures/umap-final',paste0('umap-',prefix,'-doubletscore-all.pdf')),color.label='Doublet Score')

plot.umap(this.umap,color='cluster',file=file.path('figures/umap-final',paste0('umap-',prefix,'-cluster-all.pdf')),color.label='Cluster')
plot.umap(this.umap,color='partition',file=file.path('figures/umap-final',paste0('umap-',prefix,'-partition-all.pdf')),color.label='Partition')
plot.umap(this.umap,color='region_major',file=file.path('figures/umap-final',paste0('umap-',prefix,'-regionmajor-all.pdf')),color.label='Class')
plot.umap(this.umap,color='region',file=file.path('figures/umap-final',paste0('umap-',prefix,'-region-all.pdf')),color.label='Region')
plot.umap(this.umap,color='id',file=file.path('figures/umap-final',paste0('umap-',prefix,'-sample-all.pdf')),color.label='Sample')
plot.umap(this.umap,color='sex',file=file.path('figures/umap-final',paste0('umap-',prefix,'-sex-all.pdf')),color.label='Sex')
plot.umap(this.umap,color='hemisphere',file=file.path('figures/umap-final',paste0('umap-',prefix,'-hemisphere-all.pdf')),color.label='Hemisphere')
plot.umap(this.umap,color='sequencing_run_id',file=file.path('figures/umap-final',paste0('umap-',prefix,'-sequencingrunid-all.pdf')),color.label='Sequencing Run')

plot.umap(this.umap,color='cluster',file=file.path('figures/umap-final',paste0('umap-',prefix,'-cluster-all_nolegend.pdf')),color.label='Cluster',legend=FALSE)
plot.umap(this.umap,color='partition',file=file.path('figures/umap-final',paste0('umap-',prefix,'-partition-all_nolegend.pdf')),color.label='Partition',legend=FALSE)
plot.umap(this.umap,color='region_major',file=file.path('figures/umap-final',paste0('umap-',prefix,'-regionmajor-all_nolegend.pdf')),color.label='Class',legend=FALSE)
plot.umap(this.umap,color='region',file=file.path('figures/umap-final',paste0('umap-',prefix,'-region-all_nolegend.pdf')),color.label='Region',legend=FALSE)
plot.umap(this.umap,color='id',file=file.path('figures/umap-final',paste0('umap-',prefix,'-sample-all_nolegend.pdf')),color.label='Sample',legend=FALSE)
plot.umap(this.umap,color='sex',file=file.path('figures/umap-final',paste0('umap-',prefix,'-sex-all_nolegend.pdf')),color.label='Sex',legend=FALSE)
plot.umap(this.umap,color='hemisphere',file=file.path('figures/umap-final',paste0('umap-',prefix,'-hemisphere-all_nolegend.pdf')),color.label='Hemisphere',legend=FALSE)
plot.umap(this.umap,color='sequencing_run_id',file=file.path('figures/umap-final',paste0('umap-',prefix,'-sequencingrunid-all_nolegend.pdf')),color.label='Sequencing Run',legend=FALSE)

plot.umap(this.umap,color='cluster',file=file.path('figures/umap-final',paste0('umap-',prefix,'-cluster-facet-all.pdf')),color.label='Cluster',facet=TRUE,facet.by='cluster',legend=FALSE)
plot.umap(this.umap,color='partition',file=file.path('figures/umap-final',paste0('umap-',prefix,'-partition-facet-all.pdf')),color.label='Partition',facet=TRUE,facet.by='partition',legend=FALSE)
plot.umap(this.umap,color='region_major',file=file.path('figures/umap-final',paste0('umap-',prefix,'-regionmajor-facet-all.pdf')),color.label='Class',facet=TRUE,facet.by='region_major',legend=FALSE)
plot.umap(this.umap,color='region',file=file.path('figures/umap-final',paste0('umap-',prefix,'-region-facet-all.pdf')),color.label='Region',facet=TRUE,facet.by='region',legend=FALSE)
plot.umap(this.umap,color='id',file=file.path('figures/umap-final',paste0('umap-',prefix,'-sample-facet-all.pdf')),color.label='Sample',facet=TRUE,facet.by='id',legend=FALSE)
plot.umap(this.umap,color='sex',file=file.path('figures/umap-final',paste0('umap-',prefix,'-sex-facet-all.pdf')),color.label='Sex',facet=TRUE,facet.by='sex',legend=FALSE)
plot.umap(this.umap,color='hemisphere',file=file.path('figures/umap-final',paste0('umap-',prefix,'-hemisphere-facet-all.pdf')),color.label='Hemisphere',facet=TRUE,facet.by='hemisphere',legend=FALSE)
plot.umap(this.umap,color='sequencing_run_id',file=file.path('figures/umap-final',paste0('umap-',prefix,'-sequencingrunid-facet-all.pdf')),color.label='Sequencing Run',facet=TRUE,facet.by='sequencing_run_id',legend=FALSE)

plot.umap(this.umap,color='cluster',file=file.path('figures/umap-final',paste0('umap-',prefix,'-cluster-all_presentation.pdf')),color.label='Cluster',legend=FALSE,presentation=TRUE)
plot.umap(this.umap,color='partition',file=file.path('figures/umap-final',paste0('umap-',prefix,'-partition-all_presentation.pdf')),color.label='Partition',legend=FALSE,presentation=TRUE)
plot.umap(this.umap,color='region_major',file=file.path('figures/umap-final',paste0('umap-',prefix,'-regionmajor-all_presentation.pdf')),color.label='Class',legend=FALSE,presentation=TRUE)
plot.umap(this.umap,color='region',file=file.path('figures/umap-final',paste0('umap-',prefix,'-region-all_presentation.pdf')),color.label='Region',legend=FALSE,presentation=TRUE)
plot.umap(this.umap,color='id',file=file.path('figures/umap-final',paste0('umap-',prefix,'-sample-all_presentation.pdf')),color.label='Sample',legend=FALSE,presentation=TRUE)
plot.umap(this.umap,color='sex',file=file.path('figures/umap-final',paste0('umap-',prefix,'-sex-all_presentation.pdf')),color.label='Sex',legend=FALSE,presentation=TRUE)
plot.umap(this.umap,color='hemisphere',file=file.path('figures/umap-final',paste0('umap-',prefix,'-hemisphere-all_presentation.pdf')),color.label='Hemisphere',legend=FALSE,presentation=TRUE)
plot.umap(this.umap,color='sequencing_run_id',file=file.path('figures/umap-final',paste0('umap-',prefix,'-sequencingrunid-all_presentation.pdf')),color.label='Sequencing Run',legend=FALSE,presentation=TRUE)


push.status(paste('plot.umap',prefix))

# this.umap = readRDS(file.path('umap',paste0(prefix,'-scanpy-final.rds')))
# abcam.markers = fread('int_abcam_markers.txt.gz')
# this.umap = cbind(this.umap,abcam.markers)
# 
# dir.create('figures/umap-final-markers',showWarnings=FALSE)
# for (i in setdiff(colnames(abcam.markers),'V1')) {
# 	message(i)
# 	plot.umap(this.umap,color=i,file=file.path('figures/umap-final-markers',paste0('umap-',prefix,'-cluster-all_',i,'_nolegend.pdf')),color.label=i,legend=FALSE)
# }
# 



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


this.umap=readRDS('umap/rna-scanpy-final_mindist_0.25.rds')
library(monocle3)
cds = readRDS('checkpoints/rna_cds.rds')

x = colData(cds)[c('Sample','cell','CellType_1B')]
x$id = paste0('NSM',formatC(as.integer(gsub('NSM','',x$Sample)),width=3,flag=0))
x$cell_id = x$cell
rownames(x) = x$cell = with(x,paste0(id,'-',cell_id))

y = x[c('cell','CellType_1B')]

this.umap = as.data.frame(merge(this.umap,y,by='cell',all.x=TRUE))

this.umap$cell_type = factor(this.umap$CellType_1B,
levels=c('Neurons_Excitatory',
'Neurons_Inhibitory_ADARB2',
'Neurons_Inhibitory_LHX6/SST',
'Neurons_Inhibitory_LHX6/PVALB',
'Neurons_MediumSpiny',
'Neurons_Thalamic',
'Neurons_Purkinje',
'Neurons_Basket',
'Neurons_Granule',
'Oligodendrocyte Precursors',
'Astrocytes',
'Endothelial',
'Microglia',
'Oligodendrocytes',
'Immune'),
labels=c('excitatory neurons',
'inhibitory neurons (ADARB2)',
'inhibitory neurons (LHX6/SST)',
'inhibitory neurons (LHX6/PVALB)',
'medium spiny neurons',
'thalamic neurons',
'Purkinje neurons',
'basket neurons',
'cerebellar granule cells',
'oligodendrocyte precursors',
'astrocytes',
'endothelial cells',
'microglia',
'oligodendrocytes',
'immune cells'))

plot.umap(this.umap,color='cell_type',file=file.path('figures/umap-final_mindist_0.25',paste0('umap-',prefix,'-celltype-all.pdf')),color.label='Cell Type',alpha=0.05)
plot.umap(this.umap,color='cell_type',file=file.path('figures/umap-final_mindist_0.25',paste0('umap-',prefix,'-celltype-all_nolegend.pdf')),color.label='Cell Type',legend=FALSE,alpha=0.05)
plot.umap(this.umap,color='cell_type',file=file.path('figures/umap-final_mindist_0.25',paste0('umap-',prefix,'-celltype-all_presentation.pdf')),color.label='Cell Type',legend=FALSE,alpha=0.05,presentation=TRUE)
plot.umap(this.umap,color='cell_type',file=file.path('figures/umap-final_mindist_0.25',paste0('umap-',prefix,'-celltype-facet-all.pdf')),color.label='Cell Type',facet=TRUE,facet.by='cell_type',legend=FALSE)

for (i in levels(this.umap$cell_type)) {
	message(i)
	j = gsub('[()]','',gsub('[ /]','_',i))
	this = this.umap
	this$this_cell = this$cell_type
	this$this_cell[this$this_cell != i] = NA
	this$this_cell = this$this_cell[,drop=TRUE]
	plot.umap(this,color='this_cell',file=file.path('figures/umap-final_mindist_0.25',paste0('umap-',prefix,'-showcell-',j,'-all_presentation.pdf')),color.label='Cell Type',legend=FALSE,alpha=0.05,presentation=TRUE)
}

# Split by partition

lapply(split(this.umap,this.umap$cluster),function(x) {
	y = subset(x,!is.na(cell_type))
	sort(table(y$cell_type) / nrow(y),decreasing=TRUE)[1:5] * 100
})

# > colSums(table(this.umap$cluster,this.umap$partition) > 0)
#  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
# 29 16  2 11  8  1  2  1  2  4  1  1  1  1  1

# Step one: classify with as much specificity as possible
partition.to.cell = c(
'1' = NA,
'2' = NA,
'3' = 'astrocytes',
'4' = NA,
'5' = 'oligodendrocytes',
'6' = 'oligodendrocyte precursor cells',
'7' = 'radial glial cells', #provisional (clusters 23 and 81)
'8' = 'basket cells',
'9' = 'microglia',
'10' = NA, # vascular cells (endothelial and pericytes)
'11' = 'mesenchymal stem cells', #provisional
'12' = 'ependymal cells', #provisional
'13' = 'CLIC6 neurons', # CLIC6 neurons
'14' = 'AHSG neurons', # AHSG neurons
'15' = 'F5 neurons' # F5 neurons
)

# this.umap$cell_type_round1 = this.umap$cell_type

this.umap$cell_final = NA

this.umap$cell_final = partition.to.cell[as.character(this.umap$partition)]

# lapply(split(subset(this.umap,is.na(cell_final)),subset(this.umap,is.na(cell_final))$cluster),function(x) {
# 	y = subset(x,!is.na(cell_type))
# 	sort(table(y$cell_type) / nrow(y),decreasing=TRUE)[1:5] * 100
# })

cluster.to.cell = c(
'1' = 'excitatory neurons',
'2' = 'cerebellar granule cells',
'3' = 'excitatory neurons',
'5' = 'excitatory neurons',
'6' = 'inhibitory neurons',
'8' = 'excitatory neurons',
'9' = 'inhibitory neurons',
'10' = 'cerebellar granule cells',
'11' = 'cerebellar granule cells',
'13' = 'excitatory neurons',
'14' = 'excitatory neurons',
'15' = 'inhibitory neurons',
'18' = 'excitatory neurons',
'20' = 'excitatory neurons',
'21' = 'inhibitory neurons',
'22' = 'excitatory neurons',
'24' = 'excitatory neurons',
'25' = 'inhibitory neurons',
'28' = 'medium spiny neurons',
'29' = 'endothelial cells',
'31' = 'excitatory neurons',
'32' = 'thalamic interneurons', # thalamic interneurons, subtype of inhibitory neurons
'33' = 'excitatory neurons',
'34' = 'inhibitory neurons',
'35' = 'pericytes',
'36' = 'excitatory neurons',
'37' = 'excitatory neurons',
'38' = 'cerebellar granule cells',
'39' = 'excitatory neurons',
'40' = 'endothelial cells',
'41' = 'excitatory neurons',
'42' = 'cerebellar granule cells',
'43' = 'excitatory neurons',
'44' = 'excitatory neurons',
'46' = 'excitatory neurons',
'47' = 'excitatory neurons',
'49' = 'excitatory neurons',
# '50' = NA, #partition 13
'51' = 'medium spiny neurons',
'52' = 'Purkinje cells',
'53' = 'excitatory neurons',
'54' = 'cerebellar granule cells',
'55' = 'pericytes',
'57' = 'excitatory neurons',
'58' = 'excitatory neurons',
'59' = 'excitatory neurons',
'60' = 'medium spiny neurons',
'61' = 'dopaminergic neurons',
'62' = 'serotinergic neurons',
# '63' = NA, #partition 14
'64' = 'cerebellar granule cells',
'66' = 'excitatory neurons',
'69' = 'cerebellar granule cells',
# '70' = NA, #partition 15
'71' = 'cerebellar granule cells',
'72' = 'cerebellar granule cells',
'73' = 'cerebellar granule cells',
'74' = 'cerebellar granule cells',
'75' = 'cerebellar granule cells',
'76' = 'inhibitory neurons',
'77' = 'cerebellar granule cells',
'78' = 'cerebellar granule cells',
'79' = 'inhibitory neurons'
)

this.umap$cell_final[is.na(this.umap$cell_final)] = cluster.to.cell[as.character(this.umap$cluster[is.na(this.umap$cell_final)])]

this.umap$cell_final = factor(this.umap$cell_final,levels=c(
'excitatory neurons',
'CLIC6 neurons',
'AHSG neurons',
'F5 neurons',
'medium spiny neurons',
'inhibitory neurons',
'dopaminergic neurons',
'serotinergic neurons',
'thalamic interneurons',
'cerebellar granule cells',
'Purkinje cells',
'basket cells',
'astrocytes',
'oligodendrocytes',
'oligodendrocyte precursor cells',
'endothelial cells',
'pericytes',
'microglia',
'radial glial cells',
'mesenchymal stem cells',
'ependymal cells'
),labels=c(
'excitatory neurons',
'CLIC6 neurons',
'AHSG neurons',
'F5 neurons',
'medium spiny neurons',
'inhibitory neurons',
'dopaminergic neurons',
'serotinergic neurons',
'thalamic interneurons',
'cerebellar granule cells',
'Purkinje cells',
'basket cells',
'astrocytes',
'oligodendrocytes',
'oligodendrocyte precursor cells',
'endothelial cells',
'pericytes',
'microglia',
'radial glial cells',
'mesenchymal stem cells',
'ependymal cells'
))

this.umap$cell_round1_level3 = this.umap$cell_final

this.umap$cell_round1_level2 = dplyr::recode(this.umap$cell_final,
'excitatory neurons' = 'excitatory neurons',
'CLIC6 neurons' = 'excitatory neurons',
'AHSG neurons' = 'excitatory neurons',
'F5 neurons' = 'excitatory neurons',
'medium spiny neurons' = 'medium spiny neurons',
'inhibitory neurons' = 'inhibitory neurons',
'dopaminergic neurons' = 'dopaminergic neurons',
'serotinergic neurons' = 'serotinergic neurons',
'thalamic interneurons' = 'inhibitory neurons',
'cerebellar granule cells' = 'cerebellar granule cells',
'Purkinje cells' = 'Purkinje cells',
'basket cells' = 'basket cells',
'astrocytes' = 'astrocytes',
'oligodendrocytes' = 'oligodendrocytes',
'oligodendrocyte precursor cells' = 'oligodendrocyte precursor cells',
'endothelial cells' = 'vascular cells',
'pericytes' = 'vascular cells',
'microglia' = 'microglia',
'radial glial cells' = 'radial glial cells',
'mesenchymal stem cells' = 'mesenchymal stem cells',
'ependymal cells' = 'ependymal cells'
)

this.umap$cell_round1_level1 = dplyr::recode(this.umap$cell_final,
'excitatory neurons' = 'excitatory neurons',
'CLIC6 neurons' = 'CLIC6 neurons',
'AHSG neurons' = 'AHSG neurons',
'F5 neurons' = 'F5 neurons',
'medium spiny neurons' = 'excitatory neurons',
'inhibitory neurons' = 'inhibitory neurons',
'dopaminergic neurons' = 'inhibitory neurons',
'serotinergic neurons' = 'inhibitory neurons',
'thalamic interneurons' = 'inhibitory neurons',
'cerebellar granule cells' = 'cerebellar cells',
'Purkinje cells' = 'cerebellar cells',
'basket cells' = 'basket cells',
'astrocytes' = 'astrocytes',
'oligodendrocytes' = 'oligodendrocytes',
'oligodendrocyte precursor cells' = 'oligodendrocyte precursor cells',
'endothelial cells' = 'vascular cells',
'pericytes' = 'vascular cells',
'microglia' = 'microglia',
'radial glial cells' = 'radial glial cells',
'mesenchymal stem cells' = 'mesenchymal stem cells',
'ependymal cells' = 'ependymal cells'
)

cell.type.colors = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666', '#ed1c24', '#00aeef', '#86328c', '#00aaad', '#ac7eaf', '#2e3192', '#809ead', '#6a3e14') # (extension of palette Dark2)


plot.umap(this.umap,color='cell_final',file=file.path('figures/umap-final_mindist_0.25',paste0('umap-',prefix,'-cellfinal-all.pdf')),color.label='Cell Type',alpha=0.05)
plot.umap(this.umap,color='cell_final',file=file.path('figures/umap-final_mindist_0.25',paste0('umap-',prefix,'-cellfinal-all_nolegend.pdf')),color.label='Cell Type',legend=FALSE,alpha=0.05)
plot.umap(this.umap,color='cell_final',file=file.path('figures/umap-final_mindist_0.25',paste0('umap-',prefix,'-cellfinal-all_presentation.pdf')),color.label='Cell Type',legend=FALSE,alpha=0.05,presentation=TRUE)
plot.umap(this.umap,color='cell_final',file=file.path('figures/umap-final_mindist_0.25',paste0('umap-',prefix,'-cellfinal-facet-all.pdf')),color.label='Cell Type',facet=TRUE,facet.by='cell_final',legend=FALSE)

rm(cds)
save(list=ls(),file='rna_cell_classifications.RData')

this.umap$cluster_unclassified = as.character(this.umap$cluster)
this.umap$cluster_unclassified[!is.na(this.umap$cell_final)] = '0'

plot.umap(this.umap,color='cluster_unclassified',file=file.path('figures/umap-final_mindist_0.25',paste0('umap-',prefix,'-celltodo-all.pdf')),color.label='Cell Type',alpha=0.05)
plot.umap(this.umap,color='cluster_unclassified',file=file.path('figures/umap-final_mindist_0.25',paste0('umap-',prefix,'-celltodo-all_nolegend.pdf')),color.label='Cell Type',legend=FALSE,alpha=0.05)
plot.umap(this.umap,color='cluster_unclassified',file=file.path('figures/umap-final_mindist_0.25',paste0('umap-',prefix,'-celltodo-all_presentation.pdf')),color.label='Cell Type',legend=FALSE,alpha=0.05,presentation=TRUE)
plot.umap(this.umap,color='cluster_unclassified',file=file.path('figures/umap-final_mindist_0.25',paste0('umap-',prefix,'-celltodo-facet-all.pdf')),color.label='Cell Type',facet=TRUE,facet.by='cluster_unclassified',legend=FALSE)





lapply(split(this.umap,this.umap$cluster_unclassified),function(x) {
	sort(table(x$region),decreasing=TRUE)[1:5]
})



# $`0`
# 
#     V1     CV    lCb     IT     MB
# 313058 207900 204923 168750 165770
# 
# $`23`
# 
#   STS   SPP    CC    M1    S1
# 19737 16283  5508  2006  1570
# 
# $`45`
# 
#  STS  SPP   CC   M1   S1
# 6570 3242  457  179  170
# 
# $`48`
# 
#    MB    CC vlPFC   ACC    IT
#  2360  1956   452   391   378
# 
# $`50`
# 
#  HIP Pons   V1  SPP   IT
# 1905  707  331  306  280
# 
# $`61`
# 
#   MB Pons   IT  SPP   A1
#  631  283  131  125  105
# 
# $`62`
# 
#  MB SPP  V1  A1  IT
# 812 121 121  95  82
# 
# $`81`
# 
#   STS   SPP dmPFC vmPFC dlPFC
#    36    17     0     0     0
   
lapply(split(this.umap,this.umap$cluster_unclassified),function(x) {
	sort(table(x$region)/table(this.umap$region),decreasing=TRUE)[1:5]
})


lapply(split(this.umap,this.umap$cluster_unclassified),function(x) {
	sort(table(x$animal_id)/table(this.umap$animal_id),decreasing=TRUE)[1:5]
})

# $`0`
# 
#        CV       lCb        V1       NAc       LGN
# 0.9993799 0.9986014 0.9981030 0.9971616 0.9969790
# 
# $`23`
# 
#        STS         CC        SPP        AMY         M1
# 0.18747507 0.10974954 0.09693993 0.02532009 0.02466070
# 
# $`45`
# 
#         STS         SPP          CC         MdO         MdC
# 0.062406201 0.019301066 0.009105944 0.003601749 0.002333722
# 
# $`48`
# 
#          CC          MB         ACC          CN         MdO
# 0.038974236 0.013896413 0.005294588 0.005182626 0.003859017
# 
# $`50`
# 
#         HIP        Pons         MdO         IPP       dmPFC
# 0.038323811 0.015721243 0.004888089 0.002497182 0.002462115
# 
# $`61`
# 
#        Pons          MB         MdO         MdC        mdTN
# 0.006292944 0.003715524 0.003601749 0.002167028 0.001441199
# 
# $`62`
# 
#          MB         MdO        mdTN         MdC        Pons
# 0.004781308 0.003601749 0.002594158 0.002000333 0.001601032
# 
# $`81`
# 
#          STS          SPP        dmPFC        vmPFC        dlPFC
# 0.0003419518 0.0001012085 0.0000000000 0.0000000000 0.0000000000
# 

this.umap$cluster_classified = as.character(this.umap$cell_final)
this.umap$cluster_classified[is.na(this.umap$cluster_classified)] = as.character(this.umap$cluster)[is.na(this.umap$cluster_classified)]




write.table(this.umap[,c('cell','cluster')],file=file.path('stats/clusters',paste0(prefix,'-final-clusters.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(this.umap[,c('cell','cluster_unclassified')],file='rna_cell_clusterstodo.txt',sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(this.umap[,c('cell','cluster_classified')],file=file.path('stats/clusters',paste0(prefix,'-final-celltypes.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)



moo = split(this.umap,this.umap$id)
poo = do.call(rbind,lapply(moo,function(x) data.frame(unique(x[,c('id','region','hemisphere','animal_id','social_group','sex','region_major')]),do.call(data.frame,lapply(table(x$cell_final),function(y) y)))))
region.celltype.summary = data.frame(poo[,1:7],n=unlist(lapply(moo,nrow)),t(apply(poo[,8:ncol(poo)],1,function(x) x / sum(x))))



foo = readRDS('~/Downloads/region_celltype_summary.rds')

x = unique(foo[1:8])
y = reshape2::melt(foo,id='id',measure.vars=names(foo)[9:ncol(foo)])

poo = merge(x,y,by='id')
levels(poo$variable) = gsub('\\.',' ',levels(poo$variable))

poo$variable = factor(poo$variable,levels=c(
'excitatory neurons',
'inhibitory neurons',
'medium spiny neurons',
'dopaminergic neurons',
'thalamic neurons',
'Purkinje cells',
'basket neurons',
'cerebellar granule cells',
'neural progenitor cells',
'radial glial cells',
'astrocytes',
'oligodendrocytes',
'oligodendrocyte progenitor cells',
'ependymal cells',
'endothelial cells',
'microglia'
))

p = ggplot(poo,aes(id,value,fill=variable)) +
	geom_bar(stat='identity') +
	theme_classic() +
	scale_fill_manual(values=cell.type.colors) +
	facet_wrap(~region,strip.position='top',scales='free_x') +
	theme(axis.text.x=element_blank())

poo$sample_id = with(poo,paste(animal_id,hemisphere,sep='-'))
p = ggplot(poo,aes(sample_id,value,fill=variable)) +
	geom_bar(stat='identity') +
	theme_classic() +
	scale_fill_manual(values=cell.type.colors) +
	facet_wrap(~region,strip.position='top',scales='free_x') +
	theme(axis.text.x=element_text(size=6))