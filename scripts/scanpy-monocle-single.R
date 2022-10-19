#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('atac','atac','1')

prefix = arguments[1]
analysis = arguments[2]
this.i = as.integer(arguments[3])

suppressMessages(library(monocle3))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))

sample.list = as.character(read.delim(file.path('data',paste0(prefix,'_metadata.txt')))$id)

this.sample = sample.list[this.i]

meta.file = file.path('stats/umap',paste0('meta_',prefix,'-',this.sample,'.txt.gz'))
umap02.file = file.path('stats/umap',paste0('umap02_',prefix,'-',this.sample,'.txt.gz'))
umap10.file = file.path('stats/umap',paste0('umap10_',prefix,'-',this.sample,'.txt.gz'))

meta = fread(meta.file)
rownames(meta) = meta$V1
meta$V1 = NULL

umap02 = fread(umap02.file)
names(umap02) = paste('umap',1:2,sep='.')

umap10 = fread(umap10.file)
names(umap10) = paste0('umap',1:10,sep='.')

rownames(umap02) = rownames(umap10) = rownames(meta)

e = t(umap10)
colnames(e) = rownames(meta)
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

cds = cluster_cells(cds, resolution=1e-5, reduction_method='PCA')
cds = cluster_cells(cds, resolution=1e-5, reduction_method='UMAP')

plot.umap = function(x,method='umap',color='sample',color.label=NULL,log.transform=FALSE,file=NULL,width=7,height=7,rasterize=TRUE,legend=TRUE,facet=FALSE,facet.by=NULL,size=0.1,alpha=0.05) {
	require(ggplot2)
	require(ggrastr)
	require(RColorBrewer)
	require(viridis)

	y = x

	color.label = if (is.null(color.label)) color else color.label

	scale.color = if (class(y[[color]]) %in% c('integer','numeric')) {
		scale_color_viridis(option='D',name=color.label,trans=if(log.transform) 'log10' else 'identity')
	} else if (class(y[[color]]) %in% c('character','factor','logical')) {
		if (color=='landmark') {
			scale_color_manual(values=c('#cccccc','#ff0000'),name=color.label)
		} else if (length(unique(y[[color]])) <= 8) {
			scale_color_brewer(palette='Dark2',name=color.label)
		} else {
			scale_color_discrete(name=color.label)
		}
	}

	geom.point = if (rasterize) {
		geom_point_rast(size=size,shape=19,alpha=alpha)
	} else {
		geom_point(size=size,shape=19,alpha=alpha)
	}
	
	p = ggplot(data=y,aes_string('umap.1','umap.2',color=color)) +
		geom.point +
		coord_equal() +
		scale.color +
		theme_classic() +
		theme(
			axis.text=element_blank(),
			axis.ticks=element_blank(),
			legend.position=if (legend) 'right' else 'none'
		) +
		guides(color = if (class(y[[color]]) %in% c('integer','numeric')) guide_colorbar() else guide_legend(override.aes = list(size=2,alpha=1))) +
		xlab(paste(ifelse(method=='umap','UMAP',ifelse(method=='tsne','t-SNE','Dim')),1)) +
		ylab(paste(ifelse(method=='umap','UMAP',ifelse(method=='tsne','t-SNE','Dim')),2))

	if (facet) p = p + facet_wrap(ifelse(is.null(facet.by),color,facet.by))

	if (!is.null(file)) ggsave(p,file=file,useDingbats=FALSE,width=width,height=height) else p
}

this.umap = data.frame(meta,umap02,partition=cds@clusters$PCA$partitions,cluster=cds@clusters$PCA$clusters)

this.umap$region_major = factor(c(dmPFC='Cortical',vmPFC='Cortical',dlPFC='Cortical',
vlPFC='Cortical',ACC='Cortical',CC='Subcortical',CN='Subcortical',
NAc='Subcortical',EC='Cortical',PC='Cortical',A1='Cortical',
AMY='Subcortical',HIP='Subcortical',M1='Cortical',mdTN='Subcortical',
vlTN='Subcortical',LGN='Subcortical',S1='Cortical',IPP='Cortical',
SPP='Cortical',STS='Cortical',IT='Cortical',V1='Cortical',
CV='Cerebellum',lCb='Cerebellum',MdO='Brainstem',MdC='Brainstem',
Pons='Brainstem')[this.umap$region],levels=c('Cortical','Subcortical','Cerebellum','Brainstem'))

this.umap$region = factor(this.umap$region,levels=region.levels)

dir.create('figures/umap-single',showWarnings=FALSE)
plot.umap(this.umap,color='cluster',file=file.path('figures/umap-single',paste0('umap-',prefix,'-cluster-',this.sample,'.pdf')),color.label='Cluster',size=0.5,alpha=0.5)
plot.umap(this.umap,color='partition',file=file.path('figures/umap-single',paste0('umap-',prefix,'-partition-',this.sample,'.pdf')),color.label='Partition',size=0.5,alpha=0.5)

plot.umap(this.umap,color='cluster',file=file.path('figures/umap-single',paste0('umap-',prefix,'-cluster-facet-',this.sample,'.pdf')),color.label='Cluster',facet=TRUE,facet.by='cluster',size=0.5,alpha=0.5,legend=FALSE)
plot.umap(this.umap,color='partition',file=file.path('figures/umap-single',paste0('umap-',prefix,'-partition-facet-',this.sample,'.pdf')),color.label='Partition',facet=TRUE,facet.by='partition',size=0.5,alpha=0.5,legend=FALSE)

saveRDS(this.umap,file=file.path('umap',paste0(prefix,'-scanpy-',this.sample,'.rds')))
