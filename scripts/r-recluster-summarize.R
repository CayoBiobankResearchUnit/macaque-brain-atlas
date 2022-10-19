#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('rna','rna','9')

prefix = arguments[1]
analysis = arguments[2]
this.cluster = as.integer(arguments[3])

# Read in subclustering results
this.umap = readRDS(file.path('umap',paste0(prefix,'-scanpy-recluster-class',this.cluster,'.rds')))

# Read in unified results (to grab round1 umap)
that.umap = readRDS(file.path('umap',paste0(prefix,'-scanpy-final-classified.rds')))

that.umap.coords = that.umap[rownames(this.umap),c('umap.1','umap.2')]
names(that.umap.coords) = paste('umap_all',1:2,sep='.')

names(this.umap)[names(this.umap) %in% c('umap.1','umap.2')] = paste('umap_sub',1:2,sep='.')

this.umap = data.frame(this.umap,that.umap.coords)

this.umap$subcluster = with(this.umap,paste(this.cluster,cluster2,sep='-'))
this.umap$subcluster2 = with(this.umap,paste(this.cluster,cluster3,sep='-'))
this.umap$subpartition = with(this.umap,paste(this.cluster,partition2,sep='-'))

this.umap$cell_subcluster = with(this.umap,if (nlevels(cluster2) > 1) {
	 paste(cell_class,cluster2)
} else {
	cell_class
})
this.umap$cell_subcluster2 = with(this.umap,if (nlevels(cluster3) > 1) {
	 paste(cell_class,cluster3)
} else {
	cell_class
})

this.umap$cell_subpartition = with(this.umap,if (nlevels(partition2) > 1) {
	 paste(cell_class,LETTERS[as.integer(partition2)])
} else {
	cell_class
})

saveRDS(this.umap,file=file.path('umap',paste0(prefix,'-scanpy-recluster-class',this.cluster,'-summarized.rds')))

dir.create('stats/subclusters',showWarnings=FALSE)
write.table(this.umap[,c('cell','subcluster')],file=file.path('stats/subclusters',paste0(prefix,'-class',this.cluster,'-clusters1.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(this.umap[,c('cell','subcluster2')],file=file.path('stats/subclusters',paste0(prefix,'-class',this.cluster,'-cluster2s.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(this.umap[,c('cell','subpartition')],file=file.path('stats/subclusters',paste0(prefix,'-class',this.cluster,'-partitions.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(this.umap[,c('cell','cell_subcluster')],file=file.path('stats/subclusters',paste0(prefix,'-class',this.cluster,'-cellcluster1s.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(this.umap[,c('cell','cell_subcluster2')],file=file.path('stats/subclusters',paste0(prefix,'-class',this.cluster,'-cellcluster2s.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
