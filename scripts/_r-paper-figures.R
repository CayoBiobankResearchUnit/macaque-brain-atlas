#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('atac-umap-region')

figure.list = arguments

colors = as.matrix(read.table('data/colors.txt',row.names=1,sep='\t',comment.char='',header=FALSE))[,1]
cell.levels = c('excitatory neurons','cerebellar neurons','inhibitory neurons','basket cells','medium spiny neurons','dopaminergic neurons','serotonergic neurons','AHSG neurons','F5 neurons','KIR3DL12 neurons','astrocytes','oligodendrocytes','oligodendrocyte precursor cells','radial glial cells','ependymal cells','microglia','KIR3DL12 microglia','vascular cells','mesenchymal stem cells')
cell.classes = scan(what='',sep='\n',file='stats/clusters/rna-final-cellclasses-levels.txt',quiet=TRUE)

colors = c(colors,c('all cells' = '#000000'))

library(ggplot2)
library(ggrastr)
library(RColorBrewer)
library(viridis)
library(randomcoloR)
library(egg)
library(data.table)

region.levels = c('dmPFC','vmPFC','dlPFC','vlPFC','ACC','M1','EC','PC','A1','STS','MT','IT','S1','IPP','SPP','V1','CC','CN','NAc','AMY','HIP','mdTN','vlTN','LGN','CV','lCb','MB','MdO','MdC','Pons')

# Loads cell.subtype.colors
source('subtype_colors.R')


dir.create('figures/final',showWarnings=FALSE)

if ('' %in% figures.list) {

	this.umap = readRDS(file.path('umap',paste0('biccn','-metadata-patched.rds')))
	
	na.umap = subset(this.umap,is.na(cell_class))
	
	notna.umap = do.call(rbind,mclapply(split(this.umap,this.umap$cell_class),function(x) {
		# RNA subcoordinates
		
		range1 = diff(range(x$rna_sub_umap.1,na.rm=TRUE))
		range2 = diff(range(x$rna_sub_umap.2,na.rm=TRUE))
		new_x = (x$rna_sub_umap.1 - min(x$rna_sub_umap.1,na.rm=TRUE)) / range1
		new_y = (x$rna_sub_umap.2 - min(x$rna_sub_umap.2,na.rm=TRUE)) / range2
		if (range1 > range2) {
			new_x = new_x * (range1/range2)
		} else if (range2 > range1) {
			new_y = new_y * (range2/range1)
		}
		new_x = as.numeric(scale(new_x,scale=FALSE))
		new_y = as.numeric(scale(new_y,scale=FALSE))
		x$rna_scale_sub_umap.1 = new_x
		x$rna_scale_sub_umap.2 = new_y
		
		# GLUE subcoordinates
		range1 = diff(range(x$glue_sub_umap.1,na.rm=TRUE))
		range2 = diff(range(x$glue_sub_umap.2,na.rm=TRUE))
		new_x = (x$glue_sub_umap.1 - min(x$glue_sub_umap.1,na.rm=TRUE)) / range1
		new_y = (x$glue_sub_umap.2 - min(x$glue_sub_umap.2,na.rm=TRUE)) / range2
		if (range1 > range2) {
			new_x = new_x * (range1/range2)
		} else if (range2 > range1) {
			new_y = new_y * (range2/range1)
		}
		new_x = as.numeric(scale(new_x,scale=FALSE))
		new_y = as.numeric(scale(new_y,scale=FALSE))
		x$glue_scale_sub_umap.1 = new_x
		x$glue_scale_sub_umap.2 = new_y
		x
	},mc.cores=4))
	
	this.umap = rbind(notna.umap,data.frame(na.umap,rna_scale_sub_umap.1=NA,rna_scale_sub_umap.2=NA,glue_scale_sub_umap.1=NA,glue_scale_sub_umap.2=NA))
	
	this.umap = this.umap[order(this.umap$order),]
	rownames(this.umap) = this.umap$cell
	
	shuffle.umap = this.umap[this.umap$order,]
	
	if (style == 'ggrastr') {
		p = ggplot(data=shuffle.umap,aes_string('glue_umap.1','glue_umap.2',color='rna_region',alpha='modality',size='modality')) +
			geom_point_rast(shape=19,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
			scale_color_manual(values=colors[levels(shuffle.umap$region)],name='Region',na.translate=TRUE,na.value=rgb(211,211,211,1,maxColorValue=255)) +
			scale_alpha_manual(values=c(0.05,0.025),name='Modality') +
			scale_size_manual(values=c(0.1,0.05),name='Modality') +
			coord_equal() +
			theme_void(base_size=12) +
			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
		ggsave(p,file=file.path('figures/final','glue-umap-rna-region.pdf'),useDingbats=FALSE,width=7,height=7)

		p = ggplot(data=shuffle.umap,aes_string('glue_umap.1','glue_umap.2',color='atac_region',alpha='modality',size='modality')) +
			geom_point_rast(shape=19,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
			scale_color_manual(values=colors[levels(shuffle.umap$region)],name='Region',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
			scale_alpha_manual(values=c(0.025,0.05),name='Modality') +
			scale_size_manual(values=c(0.05,0.1),name='Modality') +
			coord_equal() +
			theme_void(base_size=12) +
			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
		ggsave(p,file=file.path('figures/final','glue-umap-atac-region.pdf'),useDingbats=FALSE,width=7,height=7)

		p = ggplot(data=shuffle.umap,aes_string('glue_umap.1','glue_umap.2',color='rna_cell',alpha='modality',size='modality')) +
			geom_point_rast(shape=19,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
			scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
			scale_alpha_manual(values=c(0.05,0.025),name='Modality') +
			scale_size_manual(values=c(0.1,0.05),name='Modality') +
			coord_equal() +
			theme_void(base_size=12) +
			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
		ggsave(p,file=file.path('figures/final','glue-umap-rna-cell.pdf'),useDingbats=FALSE,width=7,height=7)

		p = ggplot(data=shuffle.umap,aes_string('glue_umap.1','glue_umap.2',color='atac_cell',alpha='modality',size='modality')) +
			geom_point_rast(shape=19,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
			scale_color_manual(values=colors[levels(shuffle.umap$atac_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
			scale_alpha_manual(values=c(0.025,0.05),name='Modality') +
			scale_size_manual(values=c(0.05,0.1),name='Modality') +
			coord_equal() +
			theme_void(base_size=12) +
			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
		ggsave(p,file=file.path('figures/final','glue-umap-atac-cell.pdf'),useDingbats=FALSE,width=7,height=7)
	
		p = ggplot(data=subset(shuffle.umap,modality=='ATAC'),aes_string('atac_umap.1','atac_umap.2',color='atac_cell')) +
			geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
			scale_color_manual(values=colors[levels(shuffle.umap$atac_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
			coord_equal() +
			theme_void(base_size=12) +
			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=2))
		ggsave(p,file=file.path('figures/final','atac-umap-cell.pdf'),useDingbats=FALSE,width=7,height=7)

		p = ggplot(data=shuffle.umap,aes_string('glue_umap.1','glue_umap.2',color='modality')) +
			geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=TRUE) +
			scale_color_manual(values=c('#e41a1c','#377eb8'),name='Modality',na.translate=FALSE,na.value=rgb(1,1,1,0)) +
			coord_equal() +
			theme_void(base_size=12) +
			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank()) +
			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
		ggsave(p,file=file.path('figures/final','glue-umap-modality.pdf'),useDingbats=FALSE,width=7,height=7)

		p = ggplot(data=subset(shuffle.umap,modality=='RNA'),aes_string('rna_umap.1','rna_umap.2',color='rna_cell')) +
			geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
			scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
			coord_equal() +
			theme_void(base_size=12) +
			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=2))
		ggsave(p,file=file.path('figures/final','rna-umap-cell.pdf'),useDingbats=FALSE,width=7,height=7)

		p = ggplot(data=subset(shuffle.umap,modality=='RNA'),aes_string('rna_umap.1','rna_umap.2',color='rna_region')) +
			geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
			scale_color_manual(values=colors[levels(shuffle.umap$rna_region)],name='Region',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
			coord_equal() +
			theme_void(base_size=12) +
			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=2))
		ggsave(p,file=file.path('figures/final','rna-umap-region.pdf'),useDingbats=FALSE,width=7,height=7)

		p = ggplot(data=subset(shuffle.umap,modality=='ATAC'),aes_string('atac_umap.1','atac_umap.2',color='atac_region')) +
			geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
			scale_color_manual(values=colors[levels(shuffle.umap$atac_region)],name='Region',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
			coord_equal() +
			theme_void(base_size=12) +
			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=2))
		ggsave(p,file=file.path('figures/final','atac-umap-region.pdf'),useDingbats=FALSE,width=7,height=7)
		
		vasc = subset(this.umap,cell_class == 'vascular cells' & modality == 'RNA') %>% droplevels
		p = ggplot(data=vasc,aes(cell_subtype,rna_doublet_score)) +
			geom_violin(width=0.5) + geom_boxplot(outlier.shape=NA,width=0.05) +
			theme_classic()
		rna.umap = droplevels(subset(this.umap,keep & modality == 'RNA'))
		
		cells.with.subtypes = names(which(table(unique(rna.umap[,c('cell_class','cell_subtype')])[,'cell_class']) > 1))
		rna.umap = droplevels(subset(rna.umap,cell_class %in% cells.with.subtypes))
		rna.umap$subtype = with(rna.umap,as.factor(as.integer(gsub('.+?([0-9]+)$','\\1',as.character(cell_subtype)))))
		
		cell.subtype.colors.map = cell.subtype.colors$color
		names(cell.subtype.colors.map) = cell.subtype.colors$cell_subtype
		
		p = ggplot(data=rna.umap,aes(subtype,rna_doublet_score)) +
			geom_violin(aes(fill=cell_subtype),width=0.5) +
			geom_boxplot(outlier.shape=NA,width=0.05) +
			facet_wrap(~cell_class,ncol=2,scales='free_x') +
			scale_y_continuous(breaks=c(0,0.2)) +
			scale_fill_manual(values=cell.subtype.colors.map[levels(rna.umap$cell_subtype)]) +
			theme_classic() +
			theme(strip.background=element_blank(),legend.position='none') +
			ylab('doublet score') +
			xlab('cell subtype')
		ggsave(p,file='figures/final/rna-subtype-doublet_scores.pdf',useDingbats=FALSE)
		
		p = ggplot(data=subset(rna.umap,cell_class %in% c('excitatory neurons','inhibitory neurons','cerebellar neurons')),aes(subtype,rna_doublet_score)) +
			geom_violin(aes(fill=cell_subtype),width=0.5) +
			geom_boxplot(outlier.shape=NA,width=0.05) +
			facet_wrap(~cell_class,ncol=1,scales='free_x') +
			scale_y_continuous(breaks=c(0,0.2)) +
			scale_fill_manual(values=cell.subtype.colors.map[levels(rna.umap$cell_subtype)]) +
			theme_classic() +
			theme(strip.background=element_blank(),legend.position='none') +
			ylab('doublet score') +
			xlab('cell subtype')
		ggsave(p,file='figures/final/rna-subtype-doublet_scores_part1.pdf',useDingbats=FALSE,width=14)
		
		p = ggplot(data=subset(rna.umap,!cell_class %in% c('excitatory neurons','inhibitory neurons','cerebellar neurons')),aes(subtype,rna_doublet_score)) +
			geom_violin(aes(fill=cell_subtype),width=0.5) +
			geom_boxplot(outlier.shape=NA,width=0.05) +
			facet_wrap(~cell_class,ncol=1,scales='free_x') +
			scale_y_continuous(breaks=c(0,0.2)) +
			scale_fill_manual(values=cell.subtype.colors.map[levels(rna.umap$cell_subtype)]) +
			theme_classic() +
			theme(strip.background=element_blank(),legend.position='none') +
			ylab('doublet score') +
			xlab('cell subtype')
		ggsave(p,file='figures/final/rna-subtype-doublet_scores_part2.pdf',useDingbats=FALSE)
			
	} else if (style == 'rasterize') {
#		p = ggplot(data=shuffle.umap,aes_string('glue_umap.1','glue_umap.2',color='rna_region',alpha='modality',size='modality')) +
#			geom_point(shape=19,na.rm=FALSE) +
#			scale_color_manual(values=colors[levels(shuffle.umap$region)],name='Region',na.translate=TRUE,na.value=rgb(211,211,211,1,maxColorValue=255)) +
#			scale_alpha_manual(values=c(0.05,0.025),name='Modality') +
#			scale_size_manual(values=c(0.1,0.05),name='Modality') +
#			coord_equal() +
#			theme_void(base_size=12) +
#			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
#			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
#		ggsave(p,file=file.path('figures/final','glue-umap-rna-region_rasterized.png'),dpi=300,width=3.9,height=3.9)
#
#		p = ggplot(data=shuffle.umap,aes_string('glue_umap.1','glue_umap.2',color='atac_region',alpha='modality',size='modality')) +
#			geom_point(shape=19,na.rm=FALSE) +
#			scale_color_manual(values=colors[levels(shuffle.umap$region)],name='Region',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
#			scale_alpha_manual(values=c(0.025,0.05),name='Modality') +
#			scale_size_manual(values=c(0.05,0.1),name='Modality') +
#			coord_equal() +
#			theme_void(base_size=12) +
#			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
#			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
#		ggsave(p,file=file.path('figures/final','glue-umap-atac-region_rasterized.png'),dpi=300,width=3.9,height=3.9)
#
#		p = ggplot(data=shuffle.umap,aes_string('glue_umap.1','glue_umap.2',color='rna_cell',alpha='modality',size='modality')) +
#			geom_point(shape=19,na.rm=FALSE) +
#			scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
#			scale_alpha_manual(values=c(0.05,0.025),name='Modality') +
#			scale_size_manual(values=c(0.1,0.05),name='Modality') +
#			coord_equal() +
#			theme_void(base_size=12) +
#			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
#			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
#		ggsave(p,file=file.path('figures/final','glue-umap-rna-cell_rasterized.png'),dpi=300,width=2,height=2)
#
#		p = ggplot(data=shuffle.umap,aes_string('glue_umap.1','glue_umap.2',color='atac_cell',alpha='modality',size='modality')) +
#			geom_point(shape=19,na.rm=FALSE) +
#			scale_color_manual(values=colors[levels(shuffle.umap$atac_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
#			scale_alpha_manual(values=c(0.025,0.05),name='Modality') +
#			scale_size_manual(values=c(0.05,0.1),name='Modality') +
#			coord_equal() +
#			theme_void(base_size=12) +
#			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
#			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
#		ggsave(p,file=file.path('figures/final','glue-umap-atac-cell_rasterized.png'),dpi=300,width=2,height=2)
#	
#		p = ggplot(data=subset(shuffle.umap,modality=='ATAC'),aes_string('atac_umap.1','atac_umap.2',color='atac_cell')) +
#			geom_point(size=0.1,shape=19,alpha=0.05,na.rm=FALSE) +
#			scale_color_manual(values=colors[levels(shuffle.umap$atac_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
#			coord_equal() +
#			theme_void(base_size=12) +
#			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
#			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=2))
#		ggsave(p,file=file.path('figures/final','atac-umap-cell_rasterized.png'),dpi=300,width=3.2,height=3.2)
#
#		p = ggplot(data=shuffle.umap,aes_string('glue_umap.1','glue_umap.2',color='modality')) +
#			geom_point(size=0.1,shape=19,alpha=0.05,na.rm=TRUE) +
#			scale_color_manual(values=c('#e41a1c','#377eb8'),name='Modality',na.translate=FALSE,na.value=rgb(1,1,1,0)) +
#			coord_equal() +
#			theme_void(base_size=12) +
#			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank()) +
#			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
#		ggsave(p,file=file.path('figures/final','glue-umap-modality_rasterized.png'),dpi=300,width=7,height=7)
#
#		p = ggplot(data=subset(shuffle.umap,modality=='RNA'),aes_string('rna_umap.1','rna_umap.2',color='rna_cell')) +
#			geom_point(size=0.1,shape=19,alpha=0.05,na.rm=FALSE) +
#			scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
#			coord_equal() +
#			theme_void(base_size=12) +
#			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
#			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=2))
#		ggsave(p,file=file.path('figures/final','rna-umap-cell_rasterized.png'),dpi=300,width=7,height=7)
#
#		p = ggplot(data=subset(shuffle.umap,modality=='RNA'),aes_string('rna_umap.1','rna_umap.2',color='rna_region')) +
#			geom_point(size=0.1,shape=19,alpha=0.05,na.rm=FALSE) +
#			scale_color_manual(values=colors[levels(shuffle.umap$rna_region)],name='Region',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
#			coord_equal() +
#			theme_void(base_size=12) +
#			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
#			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=2))
#		ggsave(p,file=file.path('figures/final','rna-umap-region_rasterized.png'),dpi=300,width=7,height=7)
#
#		p = ggplot(data=subset(shuffle.umap,modality=='ATAC'),aes_string('atac_umap.1','atac_umap.2',color='atac_region')) +
#			geom_point(size=0.1,shape=19,alpha=0.05,na.rm=FALSE) +
#			scale_color_manual(values=colors[levels(shuffle.umap$atac_region)],name='Region',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
#			coord_equal() +
#			theme_void(base_size=12) +
#			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
#			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=2))
#		ggsave(p,file=file.path('figures/final','atac-umap-region_rasterized.png'),dpi=300,width=3.8,height=3.8)
#
#
#
#


		this.umap = this.umap[order(this.umap$order),]
		
		rna.umap = droplevels(subset(this.umap,modality=='RNA' & keep))
		
		set.seed(42)
		rna.umap = rna.umap[sample(1:nrow(rna.umap)),]
		
		batch.umap = rbind(
			data.frame(rna.umap[,c('cell','region','sequencing_run_id','cell_class','rna_umap.1','rna_umap.2')],analysis=factor('after',levels=c('before','after'))),
			data.frame(within(rna.umap,{rna_umap.1=rna_nobatch.1;rna_umap.2=rna_nobatch.2})[,c('cell','region','sequencing_run_id','cell_class','rna_umap.1','rna_umap.2')],analysis=factor('before',levels=c('before','after')))
		)
		
		library(RColorBrewer)
		p = ggplot(data=batch.umap,aes_string('rna_umap.1','rna_umap.2',alpha='sequencing_run_id',color='sequencing_run_id')) +
			geom_point_rast(size=0.05,shape=19,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
			facet_wrap(~analysis,ncol=1) +
			scale_color_manual(values=c('#7570b3','#1b9e77','#cccccc'),name='Batch',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
			scale_alpha_manual(values=c(0.1,0.05,0.01),guide='none') +
			# scale_color_brewer(palette='Dark2') +
			coord_equal() +
			theme_void(base_size=12) +
			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='right') +
			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
		ggsave(p,file=file.path('figures/final','rna-umap-batch.pdf'),useDingbats=FALSE,width=7,height=7)
		
		
		
		
		# Plot subclustering
		cell.sizes = c(
			0.0000001,  # exc
			0.100,  # msn
			0.010,  # inh
			1.000,  # dop
			1.000,  # ser
			0.010,  # crb
			0.250,  # bsk
			0.050,  # ast
			0.025,  # oli
			0.050,  # opc
			0.100,  # vas
			0.250,  # mgl
			0.500,  # epn
			0.500,  # ahsg
			0.500,  # f5
			0.500,  # k3n
			0.500   # k3m
			
		)
		cell.alphas = c(
			0.05,  # exc
			0.05,  # msn
			0.05,  # inh
			0.05,  # dop
			0.05,  # ser
			0.05,  # crb
			0.05,  # bsk
			0.05,  # ast
			0.05,  # oli
			0.05,  # opc
			0.05,  # vas
			0.05,  # mgl
			0.05,  # epn
			0.05,  # ahsg
			0.05,  # f5
			0.05,  # k3n
			0.05   # k3m
		)
		
		cell.subtypes = unique(subset(this.umap,select=c('cell_class','cell_subtype')))
		cell.subtypes = cell.subtypes[order(cell.subtypes$cell_class,cell.subtypes$cell_subtype),]
		rownames(cell.subtypes) = NULL
		cell.subtypes = cell.subtypes[complete.cases(cell.subtypes),]
		
		this.umap$cell_subtype = factor(this.umap$cell_subtype,levels=cell.subtypes$cell_subtype)
		
		cell.subtypes = do.call(rbind,lapply(split(cell.subtypes,cell.subtypes$cell_class),function(x) {
			x$color = randomcoloR::randomColor(nrow(x))
			x
		}))
		rownames(cell.subtypes) = NULL
		
		# dput(cell.subtypes)
		# cell.subtypes = structure(list(cell_class=structure(c(1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,2L,2L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,4L,4L,5L,5L,6L,6L,6L,6L,6L,6L,6L,6L,6L,6L,6L,6L,6L,6L,6L,6L,7L,8L,8L,8L,8L,8L,8L,9L,9L,9L,9L,9L,9L,9L,9L,10L,11L,11L,11L,11L,11L,11L,12L,12L,13L,14L,15L,16L,17L),.Label=c('excitatory neurons','medium spiny neurons','inhibitory neurons','dopaminergic neurons','serotonergic neurons','cerebellar neurons','basket cells','astrocytes','oligodendrocytes','oligodendrocyte precursor cells','vascular cells','microglia','ependymal cells','AHSG neurons','F5 neurons','KIR3DL12 neurons','KIR3DL12 microglia'),class='factor'),cell_subtype=structure(c(1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,13L,14L,15L,16L,17L,18L,19L,20L,21L,22L,23L,24L,25L,26L,27L,28L,29L,30L,31L,32L,33L,34L,35L,36L,37L,38L,39L,40L,41L,42L,43L,44L,45L,46L,47L,48L,49L,50L,51L,52L,53L,54L,55L,56L,57L,58L,59L,60L,61L,62L,63L,64L,65L,66L,67L,68L,69L,70L,71L,72L,73L,74L,75L,76L,77L,78L,79L,80L,81L,82L,83L,84L,85L,86L,87L,88L,89L,90L,91L,92L,93L,94L,95L,96L,97L,98L,99L,100L,101L,102L,103L,104L,105L,108L,109L,110L,111L,112L),.Label=c('excitatory neurons 1','excitatory neurons 2','excitatory neurons 3','excitatory neurons 4','excitatory neurons 5','excitatory neurons 6','excitatory neurons 7','excitatory neurons 8','excitatory neurons 9','excitatory neurons 10','excitatory neurons 11','excitatory neurons 12','excitatory neurons 13','excitatory neurons 14','excitatory neurons 15','excitatory neurons 16','excitatory neurons 17','excitatory neurons 18','excitatory neurons 19','excitatory neurons 20','excitatory neurons 21','excitatory neurons 22','excitatory neurons 23','excitatory neurons 24','excitatory neurons 25','excitatory neurons 26','excitatory neurons 27','excitatory neurons 28','excitatory neurons 29','excitatory neurons 30','excitatory neurons 31','excitatory neurons 32','excitatory neurons 33','excitatory neurons 34','excitatory neurons 35','excitatory neurons 36','excitatory neurons 37','excitatory neurons 38','excitatory neurons 39','medium spiny neurons 1','medium spiny neurons 2','inhibitory neurons 1','inhibitory neurons 2','inhibitory neurons 3','inhibitory neurons 4','inhibitory neurons 5','inhibitory neurons 6','inhibitory neurons 7','inhibitory neurons 8','inhibitory neurons 9','inhibitory neurons 10','inhibitory neurons 11','inhibitory neurons 12','inhibitory neurons 13','inhibitory neurons 14','inhibitory neurons 15','inhibitory neurons 16','inhibitory neurons 17','inhibitory neurons 18','inhibitory neurons 19','inhibitory neurons 20','dopaminergic neurons 1','dopaminergic neurons 2','serotonergic neurons 1','serotonergic neurons 2','cerebellar neurons 1','cerebellar neurons 2','cerebellar neurons 3','cerebellar neurons 4','cerebellar neurons 5','cerebellar neurons 6','cerebellar neurons 7','cerebellar neurons 8','cerebellar neurons 9','cerebellar neurons 10','cerebellar neurons 11','cerebellar neurons 12','cerebellar neurons 13','cerebellar neurons 14','cerebellar neurons 15','cerebellar neurons 16','basket cells','astrocytes 1','astrocytes 2','astrocytes 3','astrocytes 4','astrocytes 5','astrocytes 6','oligodendrocytes 1','oligodendrocytes 2','oligodendrocytes 3','oligodendrocytes 4','oligodendrocytes 5','oligodendrocytes 6','oligodendrocytes 7','oligodendrocytes 8','oligodendrocyte precursor cells','vascular cells 1','vascular cells 2','vascular cells 3','vascular cells 4','vascular cells 5','vascular cells 6','microglia 1','microglia 2','radial glial cells','mesenchymal stem cells','ependymal cells','AHSG neurons','F5 neurons','KIR3DL12 neurons','KIR3DL12 microglia','serotinergic neurons 1','serotinergic neurons 2'),class='factor'),color=c('#9aaff9','#2b77ad','#6ba4e0','#d6bf3b','#59ffd8','#163175','#f2e0a7','#96bcdd','#4e04a3','#e8985f','#2eb232','#e5647c','#b1f945','#7723a8','#3eefe3','#77d1e5','#a6f975','#f4e99f','#cce4ff','#b1b8ef','#88f73d','#144bcc','#8deeef','#ef73ad','#f9c2db','#ccb20c','#e542cf','#2b86bf','#7562c4','#9495e0','#abe24a','#59eaba','#46a4c4','#cf8bf9','#809302','#1eb24d','#62e0dc','#efe773','#c0fc92','#f791bd','#8b9bf4','#9bfffa','#8f6dce','#4d9901','#66a3cc','#e26f31','#51e8dd','#a5ef8d','#27f468','#67e573','#ffcce7','#9bd352','#dd82c8','#4242d6','#3ac9a8','#b7f28a','#34d822','#9bff84','#e87896','#2c57ba','#d83227','#f4c529','#e8853a','#aa56e2','#fccabf','#7097ea','#1bf922','#3ce093','#840110','#b172cc','#50f4cb','#7d85f2','#d3b5fc','#9cfcb2','#f9f7a7','#abf791','#51db8d','#98e21f','#c9af0a','#dd0063','#ca74ed','#fc499f','#21c44f','#b5345d','#b73e12','#6287ef','#f4dd6b','#1e65d8','#9de8f2','#e3f7a3','#e835d9','#ea56b1','#55f446','#48ce39','#1be8cc','#7c35a5','#f40c04','#f28aab','#b0caf4','#601d91','#225689','#78c6cc','#e028b8','#9680e5','#f4f7a3','#a3f1f7','#e07081','#95cadb','#95d635','#57d34a')),row.names=c(NA,-110L),class='data.frame')
		
		cell.subtype.colors = cell.subtypes$color
		names(cell.subtype.colors) = cell.subtypes$cell_subtype

		p = ggplot(data=subset(this.umap,modality=='RNA' & keep),aes_string('rna_umap.1','rna_umap.2',color='cell_subtype')) +
			geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
			scale_color_manual(values=cell.subtype.colors[levels(this.umap$cell_subtype)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
			coord_equal() +
			theme_void(base_size=12) +
			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=2))
		ggsave(p,file=file.path('figures/final','rna-umap-cellsubtype.pdf'),useDingbats=FALSE,width=7,height=7)

		
		umap.labels = do.call(rbind,lapply(split(this.umap,this.umap$cell_class),function(x) {
			data.frame(
				cell_class = levels(x$cell_class[,drop=TRUE]),
				cell_subtype = levels(x$cell_subtype[,drop=TRUE]),
				rna_scale_sub_umap.1 = with(x,tapply(rna_scale_sub_umap.1,cell_subtype[,drop=TRUE],mean,na.rm=TRUE)),
				rna_scale_sub_umap.2 = with(x,tapply(rna_scale_sub_umap.2,cell_subtype[,drop=TRUE],mean,na.rm=TRUE)),
				subtype_label = gsub('.+?([0-9]*)$','\\1',levels(x$cell_subtype[,drop=TRUE]))
			)
		}))
		rownames(umap.labels) = NULL
		umap.labels$cell_class = factor(umap.labels$cell_class,levels=levels(this.umap$cell_class))
		umap.labels$cell_subtype = factor(umap.labels$cell_subtype,levels=levels(this.umap$cell_subtype))
		
		p = ggplot() +
			geom_point_rast(
				data=droplevels(subset(this.umap,!is.na(cell_class) & !grepl('^[A-Z0-9]+ ',cell_class))),
				aes(rna_scale_sub_umap.1,rna_scale_sub_umap.2,color=cell_subtype,size=cell_class,alpha=cell_class),
				shape=19,dev='ragg_png',raster.dpi=300,na.rm=TRUE) +
			facet_wrap(~cell_class,ncol=3) +
			scale_color_manual(values=cell.subtypes$color,name='Modality') +
			scale_size_manual(values=cell.sizes,guide='none') +
			scale_alpha_manual(values=cell.alphas,guide='none') +
			coord_equal() +
			theme_void(base_size=12) +
			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
		ggsave(p,file=file.path('figures/final','rna-sub-umap-cellsubtype.pdf'),useDingbats=FALSE,width=7,height=7)
		
		library(ggrepel)
		library(ggforce)
		p = ggplot() +
			geom_point_rast(
				data=droplevels(subset(this.umap,!is.na(cell_class) & !grepl('^[A-Z0-9]+ ',cell_class))),
				aes(rna_scale_sub_umap.1,rna_scale_sub_umap.2,color=cell_subtype,size=cell_class,alpha=cell_class),
				shape=19,dev='ragg_png',raster.dpi=300,na.rm=TRUE) +
			geom_text(
				data=droplevels(subset(umap.labels,!is.na(cell_class) & !grepl('^[A-Z0-9]+ ',cell_class))),
				aes(x=rna_scale_sub_umap.1,y=rna_scale_sub_umap.2,label=subtype_label),
				size=1
			) +
			facet_wrap(~cell_class,ncol=3) +
			scale_color_manual(values=cell.subtypes$color,name='Modality') +
			scale_size_manual(values=cell.sizes,guide='none') +
			scale_alpha_manual(values=cell.alphas,guide='none') +
			coord_equal() +
			theme_void(base_size=12) +
			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
		ggsave(p,file=file.path('figures/final','rna-sub-umap-cellsubtype-labels.pdf'),useDingbats=FALSE,width=7,height=7)
		
		
		r = droplevels(subset(this.umap,keep & !is.na(cell_class) & !grepl('^[A-Z0-9]+ ',cell_class)))
		
		# alphas
		# exc: 0.05
		# inh,cb: 0.1
		# 
		# for (i in c('excitatory neurons','inhibitory neurons','cerebellar neurons','oligodendrocytes','astrocytes','vascular cells','oligodendrocyte precursor cells','microglia','medium spiny neurons','microglia')) {
		for (i in c('serotonergic neurons','dopaminergic neurons','ependymal cells')) {
			# i = 'inhibitory neurons'
			a = if (i %in% c('excitatory neurons')) {
				0.05
			} else if (i %in% c('inhibitory neurons','cerebellar neurons')) {
				0.10
			} else if (i %in% c('oligodendrocytes')) {
				0.15
			} else if (i %in% c('astrocytes')) {
				0.15
			} else if (i %in% c('vascular cells')) {
				0.25
			} else if (i %in% c('oligodendrocyte precursor cells','medium spiny neurons')) {
				0.4
			} else if (i %in% c('microglia','basket cells','serotonergic neurons','dopaminergic neurons','ependymal cells')) {
				0.75
			} else if (i %in% c('serotonergic neurons','dopaminergic neurons','ependymal cells')) {
				1
			}
			sz = if (i %in% c('serotonergic neurons','dopaminergic neurons','ependymal cells')) {
				0.5
			} else {
				0.1
			}
			print(i)
		
			x = subset(r,cell_class == i) %>% droplevels
			y = subset(umap.labels,cell_class == i) %>% droplevels
			z = subset(cell.subtypes,cell_class == i) %>% droplevels
			p = ggplot() +
				geom_point_rast(
					data=x,
					aes(rna_scale_sub_umap.1,rna_scale_sub_umap.2,color=cell_subtype),
					size=sz,alpha=a,shape=19,dev='ragg_png',raster.dpi=300,na.rm=TRUE) +
				geom_text(
					data=y,
					aes(x=rna_scale_sub_umap.1,y=rna_scale_sub_umap.2,label=subtype_label),
					size=8
				) +
				scale_color_manual(values=z$color,name='Modality') +
				coord_equal() +
				theme_void(base_size=12) +
				theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
				guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
			ggsave(p,file=file.path('figures/final',paste0('rna-hires-',gsub(' ','-',i),'.pdf')),useDingbats=FALSE,width=7,height=7)
		}
		
		
		
		
		# p = ggplot() +
		# 	geom_point_rast(
		# 		data=droplevels(subset(this.umap,!is.na(cell_class) & !grepl('^[A-Z0-9]+ ',cell_class))),
		# 		aes(rna_scale_sub_umap.1,rna_scale_sub_umap.2,color=cell_subtype,size=cell_class,alpha=cell_class),
		# 		shape=19,dev='ragg_png',raster.dpi=300,na.rm=TRUE) +
		# 	geom_text(
		# 		data=droplevels(subset(umap.labels,!is.na(cell_class) & !grepl('^[A-Z0-9]+ ',cell_class))),
		# 		aes(x=rna_scale_sub_umap.1,y=rna_scale_sub_umap.2,label=subtype_label),
		# 		size=1
		# 	) +
		# 	facet_wrap_paginate(~cell_class,ncol=3) +
		# 	scale_color_manual(values=cell.subtypes$color,name='Modality') +
		# 	scale_size_manual(values=cell.sizes,guide='none') +
		# 	scale_alpha_manual(values=cell.alphas,guide='none') +
		# 	coord_equal() +
		# 	theme_void(base_size=12) +
		# 	theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='none') +
		# 	guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
		# ggsave(p,file=file.path('figures/final','rna-sub-umap-cellsubtype-facets-labels.pdf'),useDingbats=FALSE,width=7,height=7)
		
		
		subint.umap = droplevels(subset(this.umap,keep & !is.na(cell_class) & cell_class %in% atac_cell & !is.na(glue_sub_umap.1) & !is.infinite(glue_sub_umap.1)))
				
		cell.sizes = c(
			0.0000001,  # exc
			0.100,  # msn
			0.010,  # inh
			0.010,  # crb
			0.250,  # bsk
			0.050,  # ast
			0.025,  # oli
			0.050,  # opc
			0.100,  # vas
			0.250,  # mgl
			0.500   # epn
		)
		cell.alphas = c(
			0.05,  # exc
			0.05,  # msn
			0.05,  # inh
			0.05,  # crb
			0.05,  # bsk
			0.05,  # ast
			0.05,  # oli
			0.05,  # opc
			0.05,  # vas
			0.05,  # mgl
			0.05   # epn
		)
		
		names(cell.sizes) = names(cell.alphas) = c(
			'excitatory neurons','medium spiny neurons','inhibitory neurons','cerebellar neurons','basket cells','astrocytes',
			'oligodendrocytes','oligodendrocyte precursor cells','vascular cells','microglia','ependymal cells'
		)
		p = ggplot(subint.umap,aes(glue_scale_sub_umap.1,glue_scale_sub_umap.2,color=modality,size=cell_class,alpha=cell_class)) +
			geom_point_rast(shape=19,dev='ragg_png',raster.dpi=300,na.rm=TRUE) +
			facet_wrap(~cell_class,ncol=3) +
			scale_color_manual(values=c('#e41a1c','#377eb8'),name='Modality') +
			scale_size_manual(values=cell.sizes,guide='none') +
			scale_alpha_manual(values=cell.alphas,guide='none') +
			coord_equal() +
			theme_void(base_size=12) +
			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank()) +
			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
		ggsave(p,file=file.path('figures/final','glue-sub-umap-modality.pdf'),useDingbats=FALSE,width=7,height=7)
		
		p = ggplot(subint.umap,aes(glue_scale_sub_umap.1,glue_scale_sub_umap.2,color=modality,size=cell_class,alpha=cell_class)) +
			geom_point_rast(shape=19,dev='ragg_png',raster.dpi=300,na.rm=TRUE) +
			facet_grid(cell_class~modality) +
			scale_color_manual(values=c('#e41a1c','#377eb8'),name='Modality') +
			scale_size_manual(values=cell.sizes,guide='none') +
			scale_alpha_manual(values=cell.alphas,guide='none') +
			coord_equal() +
			theme_void(base_size=12) +
			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),strip.text.y=element_text(hjust=0)) +
			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
		ggsave(p,file=file.path('figures/final','glue-sub-umap-modality_faceted.pdf'),useDingbats=FALSE,width=7,height=7)

		# cell.subtypes = unique(subset(subint.umap,select=c('cell_class','cell_subtype')))
		# cell.subtypes = cell.subtypes[order(cell.subtypes$cell_class,cell.subtypes$cell_subtype),]
		# rownames(cell.subtypes) = NULL
		# cell.subtypes = do.call(rbind,lapply(split(cell.subtypes,cell.subtypes$cell_class),function(x) {
		# 	x$color = randomcoloR::randomColor(nrow(x))
		# 	x
		# }))
		# rownames(cell.subtypes) = NULL

		umap.labels = do.call(rbind,mclapply(split(subint.umap,subint.umap$cell_class),function(x) {
			data.frame(
				cell_class = levels(x$cell_class[,drop=TRUE]),
				cell_subtype = levels(x$cell_subtype[,drop=TRUE]),
				glue_scale_sub_umap.1 = with(x,tapply(glue_scale_sub_umap.1,cell_subtype[,drop=TRUE],mean,na.rm=TRUE)),
				glue_scale_sub_umap.2 = with(x,tapply(glue_scale_sub_umap.2,cell_subtype[,drop=TRUE],mean,na.rm=TRUE)),
				subtype_label = gsub('.+?([0-9]*)$','\\1',levels(x$cell_subtype[,drop=TRUE]))
			)
		},mc.cores=4))
		rownames(umap.labels) = NULL
		umap.labels$cell_class = factor(umap.labels$cell_class,levels=levels(this.umap$cell_class))
		umap.labels$cell_subtype = factor(umap.labels$cell_subtype,levels=levels(this.umap$cell_subtype))
		
		
		p = ggplot() +
			geom_point_rast(
				data = subint.umap,
				aes(glue_scale_sub_umap.1,glue_scale_sub_umap.2,color=cell_subtype,size=cell_class,alpha=cell_class),
				shape=19,dev='ragg_png',raster.dpi=300,na.rm=TRUE
			) +
			facet_grid(cell_class~modality) +
			scale_color_manual(values=cell.subtype.colors[levels(subint.umap$cell_subtype)],name='Subtype') +
			scale_size_manual(values=cell.sizes,guide='none') +
			scale_alpha_manual(values=cell.alphas,guide='none') +
			coord_equal() +
			theme_void(base_size=12) +
			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),strip.text.y=element_text(hjust=0),legend.position='none') +
			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
		ggsave(p,file=file.path('figures/final','glue-sub-umap-modality_subtype_faceted.pdf'),useDingbats=FALSE,width=7,height=7)
		
		p = ggplot() +
			geom_point_rast(
				data = subint.umap,
				aes(glue_scale_sub_umap.1,glue_scale_sub_umap.2,color=cell_subtype,size=cell_class,alpha=cell_class),
				shape=19,dev='ragg_png',raster.dpi=300,na.rm=TRUE
			) +
			geom_text(data=umap.labels,aes(glue_scale_sub_umap.1,glue_scale_sub_umap.2,label=subtype_label),size=0.5) +
			facet_grid(cell_class~modality) +
			scale_color_manual(values=cell.subtype.colors[levels(subint.umap$cell_subtype)],name='Subtype') +
			scale_size_manual(values=cell.sizes,guide='none') +
			scale_alpha_manual(values=cell.alphas,guide='none') +
			coord_equal() +
			theme_void(base_size=12) +
			theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),strip.text.y=element_text(hjust=0),legend.position='none') +
			guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
		ggsave(p,file=file.path('figures/final','glue-sub-umap-modality_subtype_faceted-labels.pdf'),useDingbats=FALSE,width=7,height=7)

		
		for (i in cell.classes) {
			class.number = match(i,cell.classes)
			class.name = gsub('neuron$','neurons',gsub('serotinergic','serotonergic',i))
			
			if (class.name %in% levels(subint.umap$cell_class)) {
				message(class.number,': ',class.name)
		
				p = ggplot() +
					geom_point_rast(
						data = subset(subint.umap,cell_class %in% class.name),
						aes(glue_scale_sub_umap.1,glue_scale_sub_umap.2,color=cell_subtype,size=cell_class,alpha=cell_class),
						shape=19,dev='ragg_png',raster.dpi=300,na.rm=TRUE
					) +
					geom_text(data=subset(umap.labels,cell_class %in% class.name),aes(glue_scale_sub_umap.1,glue_scale_sub_umap.2,label=subtype_label),size=2.5) +
					facet_wrap(~modality,nrow=1) +
					scale_color_manual(values=cell.subtype.colors[levels(subint.umap$cell_subtype)],name='Subtype',na.translate=TRUE,na.value=rgb(233,233,233,0.1,maxColorValue=255)) +
					scale_size_manual(values=cell.sizes[class.name],guide='none') +
					scale_alpha_manual(values=cell.alphas[class.name],guide='none') +
					coord_equal() +
					theme_void(base_size=12) +
					theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),strip.text.y=element_text(hjust=0),legend.position='none') +
					guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
				ggsave(p,file=file.path('figures/final',paste0('glue-sub-umap-modality_subtype_class',class.number,'.pdf')),useDingbats=FALSE,width=7,height=7)
			}
		}
	}
	









	this.cell.proportions = with(subset(this.umap,!is.na(rna_cell) | !is.na(atac_region)),rbind(
		data.frame(
			region = factor(region.levels,levels=region.levels),
			counts = as.integer(table(atac_region)[region.levels]),
			type = factor('ATAC',levels=c('RNA','ATAC','total'))
		),
		data.frame(
			region = factor(region.levels,levels=region.levels),
			counts = as.integer(table(rna_region)[region.levels]),
			type = factor('RNA',levels=c('RNA','ATAC','total'))
		),
		data.frame(
			region = factor(region.levels,levels=region.levels),
			counts = as.integer(table(region)[region.levels]),
			type = factor('total',levels=c('RNA','ATAC','total'))
		)
	))
	p = ggplot(droplevels(subset(this.cell.proportions,type %in% c('RNA','ATAC','total'))),aes(region,counts/1000,fill=region)) +
		geom_bar(stat='identity') +
		facet_wrap(~type,nrow=1) +
		coord_flip() +
		scale_x_discrete(limits=rev) +
		scale_y_continuous(breaks=seq(0,400,200)) +
		scale_fill_manual(values=colors[levels(this.cell.proportions$region)],name='Region') +
		theme_classic(base_size=12) +
		theme(
			legend.position = 'none',
			# axis.text.x = element_text(vjust=1,hjust=0,angle=-30),
			strip.background = element_blank(),
			axis.title.y = element_blank()
		) +
		ylab('cell count (thousands)')
	ggsave(p,file=file.path('figures/final','integrated-region-counts.pdf'),useDingbats=FALSE,width=7,height=6)
		
	
	# Terrible name, but now calculate proportions
	these.cell.proportions = with(droplevels(subset(this.umap,(!is.na(rna_cell) | !is.na(atac_region)) & cell_class %in% names(which(table(atac_cell) & table(rna_cell))))),data.frame(
		cell = factor(levels(cell_class),levels=levels(cell_class)),
		rna_proportions = as.integer(table(rna_cell)[levels(cell_class)])/sum(!is.na(rna_cell)),
		atac_proportions = as.integer(table(atac_cell)[levels(cell_class)])/sum(!is.na(atac_cell))
	))
	these.cell.proportions = within(these.cell.proportions,{
		log_rna_proportions = log2(rna_proportions + 1)
		log_atac_proportions = log2(atac_proportions + 1)
	})

	p = ggplot(these.cell.proportions) +
		geom_abline(slope=1,intercept=0,size=0.25,color='#000000',linetype=1) +
		# geom_smooth(aes(log_rna_proportions,log_atac_proportions),method=lm,size=1,color='#000000',fill=rgb(233,233,233,10,maxColorValue=255),linetype=1) +
		geom_point(aes(rna_proportions,atac_proportions,color=cell),size=2) +
		coord_equal() +
		scale_x_continuous(limits=c(0,0.5),breaks=seq(0,0.5,0.1)) +
		scale_y_continuous(limits=c(0,0.5),breaks=seq(0,0.5,0.1)) +
		scale_color_manual(values=colors[levels(these.cell.proportions$cell)],name='Cell type') +
		theme_classic(base_size=12) +
		theme(
			legend.position = 'left',
			legend.spacing.y=unit(0,'cm'),
			legend.text=element_text(margin=margin(t=0)),
			legend.title=element_blank()
			# legend.position = 'none',
			# axis.text.x = element_text(vjust=1,hjust=0,angle=-30),
		) +
		xlab('RNA cell proportion') +
		ylab('ATAC cell proportion')
	ggsave(p,file=file.path('figures/final','integrated-cell-proportions-correlation.pdf'),useDingbats=FALSE,width=7,height=6)

	with(these.cell.proportions,cor.test(rna_proportions,atac_proportions,method='spearman'))
	
	# Now cell proportions by region
	rna.region.cell.proportions = with(droplevels(subset(this.umap,(!is.na(rna_cell) | !is.na(atac_region)) & cell_class %in% names(which(table(atac_cell) & table(rna_cell))))),reshape2::melt(apply(as.matrix(table(rna_cell,rna_region)),2,function(x) x/sum(x))))
	atac.region.cell.proportions = with(droplevels(subset(this.umap,(!is.na(rna_cell) | !is.na(atac_region)) & cell_class %in% names(which(table(atac_cell) & table(rna_cell))))),reshape2::melt(apply(as.matrix(table(atac_cell,atac_region)),2,function(x) x/sum(x))))
	names(rna.region.cell.proportions) = c('cell','region','rna_proportion')
	names(atac.region.cell.proportions) = c('cell','region','atac_proportion')
	region.cell.proportions = merge(rna.region.cell.proportions,atac.region.cell.proportions,by=c('cell','region'))
	region.cell.proportions$cell = factor(region.cell.proportions$cell,levels=intersect(cell.levels,region.cell.proportions$cell))
	region.cell.proportions$region = factor(region.cell.proportions$region,levels=intersect(region.levels,region.cell.proportions$region))
	
	# region.cell.proportions = droplevels(subset(region.cell.proportions,!region %in% c('MB','MT')))
	
	region.cell.correlations = do.call(rbind,lapply(
		split(region.cell.proportions,region.cell.proportions$region),
		function(x) {
			result = cor.test(x$rna_proportion,x$atac_proportion)
			data.frame(region=unique(x$region),rho=result$estimate,p=result$p.value)
		}
	))
	region.cell.correlations$padj = p.adjust(region.cell.correlations$p,'fdr')
	region.cell.correlations$n = as.integer(with(droplevels(subset(this.umap,(!is.na(rna_cell) | !is.na(atac_region)) & cell_class %in% names(which(table(atac_cell) & table(rna_cell))))),table(region))[rownames(region.cell.correlations)])

	p = ggplot(region.cell.correlations) +
		geom_point(aes(region,rho,color=region)) +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.5),minor_breaks=seq(0,1,0.25)) +
		scale_color_manual(values=colors[levels(region.cell.correlations$region)],name='Region') +
		theme_classic(base_size=12) +
		theme(
			legend.position = 'none',
			panel.grid.major.y = element_line(color='#aaaaaa',size=0.25),
			panel.grid.minor.y = element_line(color='#cccccc',size=0.1),
			axis.title.x = element_blank(),
			axis.text.x = element_text(vjust=1,hjust=0,angle=-45)
		) +
		ylab(expression('Spearman\'s'~rho))
	ggsave(p,file=file.path('figures/final','integrated-cell-proportions-correlation-region.pdf'),useDingbats=FALSE,width=7,height=2)
	
#	p = ggplot(data=shuffle.umap,aes_string('glue_umap.1','glue_umap.2',color='region')) +
#		geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=TRUE) +
#		scale_color_manual(values=colors[levels(this.umap$region)],name='Region',na.translate=FALSE,na.value=rgb(1,1,1,0)) +
#		facet_wrap(~modality) +
#		coord_equal() +
#		theme_void(base_size=12) +
#		theme(legend.spacing.y=unit(0,'cm'),legend.text=element_text(margin=margin(t=0)),legend.title=element_blank(),legend.position='top') +
#		guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=15))
#	ggsave(p,file=file.path('figures/final','glue-umap-region-modality.pdf'),useDingbats=FALSE,width=7,height=7)
	
	
	
	
	
	
	
	

} else {

fill.colors = colors
fill.colors['all cells'] = '#ffffff'


library(data.table)

regulatory.results = do.call(rbind,mclapply(c(0:3,6:13,15),function(this.cluster) {

	if (!this.cluster) {
		metacell.lr.results = read.table(
			file.path('stats/metacell_lr',paste0('gene_peak_all.tsv')),
			col.names = c('meta.beta','meta.serr','meta.sbet','meta.pval','peak','ensembl_gene_id'),
			sep='\t'
		)

		metacell.lr.null = read.table(
			file.path('stats/metacell_lr',paste0('gene_peak_null_all.tsv')),
			col.names = c('null.beta','null.serr','null.sbet','null.pval','peak','ensembl_gene_id'),
			sep='\t'
		)

		metacell.lr.results = merge(metacell.lr.results,metacell.lr.null,all.x=TRUE,all.y=TRUE,by=c('ensembl_gene_id','peak'))

		glue.regulatory.results = fread(
			file.path('stats/glue',paste0('biccn','_regulatory_regulatory.txt.gz'))
		)[,c('source','target','qval','dist','pval','score','weight','type')]
		names(glue.regulatory.results) = c('ensembl_gene_id','peak','glue.qval','dist','glue.pval','glue.score','glue.weight','glue.type')

		results.combined = merge(metacell.lr.results,glue.regulatory.results,by=c('ensembl_gene_id','peak'),all.y=TRUE,all.x=TRUE)

		results.combined$meta.qval = p.adjust(results.combined$meta.pval,'fdr')
		results.combined$cell = 'all cells'
	} else {
		this.cell.class = cell.classes[this.cluster]

		metacell.lr.results = read.table(
			file.path('stats/metacell_lr',paste0('gene_peak_class',this.cluster,'.tsv')),
			col.names = c('meta.beta','meta.serr','meta.sbet','meta.pval','peak','ensembl_gene_id'),
			sep='\t'
		)

		metacell.lr.null = read.table(
			file.path('stats/metacell_lr',paste0('gene_peak_null_class',this.cluster,'.tsv')),
			col.names = c('null.beta','null.serr','null.sbet','null.pval','peak','ensembl_gene_id'),
			sep='\t'
		)

		metacell.lr.results = merge(metacell.lr.results,metacell.lr.null,all.x=TRUE,all.y=TRUE,by=c('ensembl_gene_id','peak'))

		glue.regulatory.results = fread(
			file.path('stats/glue/subintegration',paste0('subpeak','_class',this.cluster,'_regulatory.txt.gz'))
		#	col.names = c('ensembl_gene_id','peak','glue.qval','dist','glue.pval','glue.score','glue.weight','glue.type')
		)[,c('source','target','qval','dist','pval','score','weight','type')]
		names(glue.regulatory.results) = c('ensembl_gene_id','peak','glue.qval','dist','glue.pval','glue.score','glue.weight','glue.type')

		results.combined = merge(metacell.lr.results,glue.regulatory.results,by=c('ensembl_gene_id','peak'),all.y=TRUE,all.x=TRUE)

		results.combined$meta.qval = p.adjust(results.combined$meta.pval,'fdr')

		if (any(!complete.cases(results.combined))) {
			warning(sum(!complete.cases(results.combined)),' records are incomplete and will be removed.')
			warning(sum(is.na(results.combined$meta.pval)),' metacell p-values are NA.')
			warning(sum(is.na(results.combined$null.pval)),' permuted metacell p-values are NA.')
			results.combined = results.combined[complete.cases(results.combined),]
		}
		results.combined$cell = this.cell.class
	}
	results.combined
},mc.cores=4))


regulatory.results$dist_category = factor(with(regulatory.results,
	ifelse(dist <= 25000,'0-25 kb',
		ifelse(dist <= 50000,'25-50 kb',
			ifelse(dist <= 75000,'50-75 kb',
				ifelse(dist <= 100000,'75-100 kb',
					ifelse(dist <= 125000,'100-125 kb',
						'125-150 kb')
					)
				)
			)
		)
	),levels=c('0-25 kb','25-50 kb','50-75 kb','75-100 kb','100-125 kb','125-150 kb'))

# regulatory.results$cell[regulatory.results$cell == 'serotinergic neurons'] = 'serotonergic neurons'

regulatory.results$cell = factor(regulatory.results$cell,levels=c('all cells',cell.levels))[,drop=TRUE]



# Recalculate distance
head(regulatory.results)
peak.coordinates = as.data.frame.matrix(do.call(rbind,strsplit(regulatory.results$peak,'_')))

names(peak.coordinates) = c('peak_chromosome','peak_start','peak_end')

peak.coordinates = within(peak.coordinates,{
	peak_start = as.integer(peak_start)
	peak_end = as.integer(peak_end)
})

regulatory.results = data.frame(regulatory.results,peak.coordinates)

ens.genes = readRDS(file.path('checkpoints',paste0(ens.species,'_genes.rds')))

rownames(ens.genes) = ens.genes$ensembl_gene_id

regulatory.results = data.frame(regulatory.results,ens.genes[regulatory.results$ensembl_gene_id,2:ncol(ens.genes)])


regulatory.results$tss = with(regulatory.results,
	ifelse(strand > 0,start_position,end_position)
)

regulatory.results$tss_dist = with(regulatory.results,
	ifelse(
		strand > 0,
		ifelse(	# case if + strand
			peak_end < tss, # if the peak ends before tss,
			peak_end - tss, # the peak is upstream and sign of distance is negative
			ifelse(
				peak_start > tss, # If the peak starts after tss
				peak_start - tss, # the peak is downstream and sign of distance is positive
				0
			)
		),
		ifelse(	# case if - strand
			peak_start > tss, # if the peak starts after tss,
			tss - peak_start, # the peak is upstream and sign of distance is negative
			ifelse(
				peak_end < tss, # if the peak ends before tss
				tss - peak_end, # the peak is downstream and sign of distance is positive
				0
			)
		)
	)	
)


regulatory.results$tss_dist_category = factor(with(regulatory.results,
	ifelse(tss_dist == 0,
		'overlapping',
		ifelse(tss_dist > 0,
			ifelse(tss_dist <= 25000,'0-25 kb downstream',
				ifelse(tss_dist <= 50000,'25-50 kb downstream',
					ifelse(tss_dist <= 75000,'50-75 kb downstream',
						ifelse(tss_dist <= 100000,'75-100 kb downstream',
							ifelse(tss_dist <= 125000,'100-125 kb downstream',
								'125-150 kb downstream')
						)
					)
				)
			),
			ifelse(tss_dist >= -25000,'0-25 kb upstream',
				ifelse(tss_dist >= -50000,'25-50 kb upstream',
					ifelse(tss_dist >= -75000,'50-75 kb upstream',
						ifelse(tss_dist >= -100000,'75-100 kb upstream',
							ifelse(tss_dist >= -125000,'100-125 kb upstream',
								'125-150 kb upstream')
						)
					)
				)
			)
		)
	)
),levels=c('125-150 kb upstream','100-125 kb upstream','75-100 kb upstream','50-75 kb upstream','25-50 kb upstream','0-25 kb upstream','overlapping','0-25 kb downstream','25-50 kb downstream','50-75 kb downstream','75-100 kb downstream','100-125 kb downstream','125-150 kb downstream'))



regulatory.results$tss_dist_label = factor(
	regulatory.results$tss_dist_category,
	levels=c('125-150 kb upstream','100-125 kb upstream','75-100 kb upstream','50-75 kb upstream','25-50 kb upstream','0-25 kb upstream','overlapping','0-25 kb downstream','25-50 kb downstream','50-75 kb downstream','75-100 kb downstream','100-125 kb downstream','125-150 kb downstream'),
	labels=c('[125-150]-','[100-125]-','[75-100]-','[50-75]-','[25-50]-','[0-25]-','0','[0-25]+','[25-50]+','[50-75]+','[75-100]+','[100-125]+','[125-150]+')
)


regulatory.results$tss_dist_label2 = factor(
	regulatory.results$tss_dist_category,
	levels=c('125-150 kb upstream','100-125 kb upstream','75-100 kb upstream','50-75 kb upstream','25-50 kb upstream','0-25 kb upstream','overlapping','0-25 kb downstream','25-50 kb downstream','50-75 kb downstream','75-100 kb downstream','100-125 kb downstream','125-150 kb downstream'),
	labels=c(' 125-150 kb',' 100-125 kb','  75-100 kb','  50-75 kb','  25-50 kb','   0-25 kb','0 kb','0-25 kb','25-50 kb','50-75 kb','75-100 kb','100-125 kb','125-150 kb')
)

regulatory.results$meta.padj = regulatory.results$meta.qval
regulatory.results$meta.padj[regulatory.results$meta.padj == 0] = min(regulatory.results$meta.padj[regulatory.results$meta.padj > 0])

regulatory.results$log.meta.padj = -log10(regulatory.results$meta.padj)




regulatory.results$gene_dist = with(regulatory.results,
	ifelse(
		strand > 0,
		ifelse(	# case if + strand
			peak_end < tss, # if the peak ends before tss,
			tss - peak_end, # the distance from the gene is the distance from tss
			ifelse(
				peak_start > end_position, # If the peak starts after the gene end
				peak_start - end_position, # the peak is downstream of the whole gene
				0
			)
		),
		ifelse(	# case if - strand
			peak_start > tss, # if the peak starts after tss,
			peak_start - tss, # the distance from the gene is the distance from tss
			ifelse(
				peak_end < start_position, # if the peak ends before tss
				start_position - peak_start, # the peak is downstream and sign of distance is positive
				0
			)
		)
	)	
)






de.genes = readRDS(paste0('rds/','rna','_marker_genes.rds'))

de.genes = within(de.genes,{
	cell = factor(cell)
	logfoldchanges = as.numeric(logfoldchanges)
	wilcox_score = as.numeric(wilcox_score)
	ttest_score = as.numeric(ttest_score)
	logreg_score = as.numeric(logreg_score)
	bp1 = as.integer(bp1)
	bp2 = as.integer(bp2)
	gene_strand = factor(gene_strand,levels=c('+','-'))
})




p = ggplot(droplevels(regulatory.results),aes(dist_category,glue.score,fill=cell)) +
	geom_boxplot(outlier.shape=NA) +
	theme_classic() +
	facet_wrap(~cell) +
	scale_fill_manual(values=colors[levels(regulatory.results$cell)]) +
	# xlab('Genomic distance') +
	ylab('GLUE regulatory score') +
	theme(
		axis.title.x=element_blank(),
		axis.text.x=element_text(angle=75, vjust = 1, hjust=1),
		legend.position='none',
		strip.background=element_blank())
ggsave(p,file=file.path('figures/final',paste0('regulatory_dist_glue_scores.pdf')),useDingbats=FALSE,height=7,width=7)




p0 = ggplot(droplevels(subset(regulatory.results,cell == 'all cells')),aes(tss_dist_label2,glue.score,fill=cell)) +
	geom_boxplot(outlier.shape=NA,width=0.25) +
	theme_classic(base_size=16) +
	# facet_wrap(~cell) +
	scale_fill_manual(values=fill.colors[levels(regulatory.results$cell)]) +
	# xlab('Genomic distance') +
	ylab('GLUE regulatory score') +
	theme(
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.text.x=element_text(angle=75, vjust = 1, hjust=1),
		legend.position='none',
		strip.background=element_blank())

p = ggplot(droplevels(subset(regulatory.results,cell != 'all cells')),aes(tss_dist_label2,glue.score,fill=cell)) +
	geom_boxplot(outlier.shape=NA) +
	theme_classic() +
	facet_wrap(~cell) +
	scale_fill_manual(values=fill.colors[levels(regulatory.results$cell)]) +
	# xlab('Genomic distance') +
	ylab('GLUE regulatory score') +
	theme(
		axis.title.x=element_blank(),
		axis.text.x=element_text(angle=75, vjust = 1, hjust=1),
		legend.position='none',
		strip.background=element_blank())
# ggsave(p,file=file.path('figures/final',paste0('regulatory_tssdist_glue_scores.pdf')),useDingbats=FALSE,height=7,width=7)

pdf(file='/dev/null')
p = ggarrange(p0,p,nrow=2,heights=c(1,3))
dev.off()

pdf(file=file.path('figures/final',paste0('regulatory_tssdist_glue_scores.pdf')),useDingbats=FALSE,height=7,width=7)
p
dev.off()


# plot.list = lapply(levels(regulatory.results$cell),function(i) {
# 	x = regulatory.results[regulatory.results$cell == i,]
# 	x$sig = factor(x$glue.qval < 0.05,levels=c('TRUE','FALSE'))
# 	p = ggplot(droplevels(x),aes(dist_category,log.meta.padj,fill=cell,alpha=sig)) +
# 		geom_boxplot(outlier.shape=NA) +
# 		theme_classic() +
# 		coord_cartesian(ylim = boxplot.stats(x$log.meta.padj)$stats[c(1,5)]*1.05) +
# 		scale_fill_manual(values=fill.colors[levels(regulatory.results$cell)]) +
# 		scale_alpha_manual(values=c(1,0.5)) +
# 		ylab(expression(-log[10]~italic(p))) +
# 		theme(
# 			axis.title.x=element_blank(),
# 			axis.text.x = if (i %in% rev(levels(regulatory.results$cell))[1:4]) element_text(angle=75, vjust = 1, hjust=1) else element_blank(),
# 			axis.title.y = if (i %in% levels(regulatory.results$cell)[seq(1,nlevels(regulatory.results$cell),4)]) element_text() else element_blank(),
# 			plot.title=element_text(size=6),
# 			legend.position='none'
# 		) +
# 		ggtitle(i)
# 	p
# })
# 
# plot.list = c(plot.list,list(ncol=4))
# 
# pdf(file='/dev/null')
# p = do.call(ggarrange,plot.list)
# dev.off()
# 
# pdf(file=file.path('figures/final',paste0('regulatory_dist_metacell_padj_bysig.pdf')),useDingbats=FALSE,height=7,width=7)
# p
# dev.off()



fill.colors['all cells'] = '#dddddd'

regulatory.results.stats = do.call(rbind,lapply(split(regulatory.results,list(regulatory.results$cell,regulatory.results$tss_dist_label,factor(regulatory.results$glue.qval < 0.05,levels=c('TRUE','FALSE')))),function(x) {
	if (nrow(x) > 1) {
		this.cell = unique(x$cell)
		this.dist = unique(x$tss_dist_label)
		this.sig = factor(unique(x$glue.qval < 0.05),levels=c('TRUE','FALSE'))
		x.summary = do.call(data.frame,as.list(c(boxplot.stats(x$log.meta.padj)$stats,boxplot.stats(x$meta.sbet)$stats,boxplot.stats(x$meta.beta)$stats)))
		names(x.summary) = c('pval.whisker.min','pval.q1','pval.median','pval.q3','pval.whisker.max','sbet.whisker.min','sbet.q1','sbet.median','sbet.q3','sbet.whisker.max','beta.whisker.min','beta.q1','beta.median','beta.q3','beta.whisker.max')
		data.frame(x.summary,cell=this.cell,dist=this.dist,sig=this.sig)
	} else {
		data.frame(pval.whisker.min = numeric(0), pval.q1 = numeric(0), pval.median = numeric(0), pval.q3 = numeric(0), pval.whisker.max = numeric(0), 
		sbet.whisker.min = numeric(0), sbet.q1 = numeric(0), sbet.median = numeric(0), sbet.q3 = numeric(0), sbet.whisker.max = numeric(0), 
		beta.whisker.min = numeric(0), beta.q1 = numeric(0), beta.median = numeric(0), beta.q3 = numeric(0), beta.whisker.max = numeric(0), 
		cell = character(0), dist = character(0), sig = character(0), row.names = integer(0))
	}
}))
rownames(regulatory.results.stats) = NULL

dodge = position_dodge(width=0.75)
p0 = ggplot(droplevels(subset(regulatory.results.stats,cell == 'all cells')),aes(dist,pval.median)) +
	geom_linerange(aes(ymin=pval.whisker.min,ymax=pval.q1,linetype=sig),position=dodge) +
	geom_linerange(aes(ymin=pval.q3,ymax=pval.whisker.max,linetype=sig),position=dodge) +
	geom_crossbar(aes(ymin=pval.q1,ymax=pval.q3,fill=cell,alpha=sig),position=dodge,width=0.5) +
	theme_classic() +
	scale_fill_manual(values=fill.colors['all cells']) +
	scale_alpha_manual(values=c(1,0.5)) +
	scale_linetype_manual(values=c(1,1)) +
	ylab(expression(-log[10]~italic(p['adj']))) +
	theme(
		axis.text.x = element_text(angle=75, vjust = 1, hjust=1),
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		legend.position='none'
	)

dodge = position_dodge(width=0.75)
p1 = ggplot(droplevels(subset(regulatory.results.stats,cell != 'all cells')),aes(dist,pval.median)) +
	geom_linerange(aes(ymin=pval.whisker.min,ymax=pval.q1,linetype=sig),position=dodge) +
	geom_linerange(aes(ymin=pval.q3,ymax=pval.whisker.max,linetype=sig),position=dodge) +
	geom_crossbar(aes(ymin=pval.q1,ymax=pval.q3,fill=cell,alpha=sig),position=dodge,width=0.5) +
	theme_classic() +
	scale_fill_manual(values=fill.colors[levels(droplevels(subset(regulatory.results.stats,cell != 'all cells'))$cell)]) +
	scale_alpha_manual(values=c(1,0.5)) +
	scale_linetype_manual(values=c(1,1)) +
	facet_wrap(~cell,ncol=3,scales='free_y') +
	ylab(expression(-log[10]~italic(p['adj']))) +
	theme(
		axis.text.x = element_text(angle=75, vjust = 1, hjust=1),
		axis.title.x = element_blank(),
		strip.background = element_blank(),
		legend.position='none'
	)

pdf(file='/dev/null')
p = egg::ggarrange(p0,p1,nrow=2,heights=c(1,4))
dev.off()

pdf(file=file.path('figures/final',paste0('regulatory_tssdist_metacell_padj_bysig.pdf')),useDingbats=FALSE,height=7,width=7)
p
dev.off()



dodge = position_dodge(width=0.75)
p0 = ggplot(droplevels(subset(regulatory.results.stats,cell == 'all cells')),aes(dist,beta.median)) +
	geom_linerange(aes(ymin=beta.whisker.min,ymax=beta.q1,linetype=sig),position=dodge) +
	geom_linerange(aes(ymin=beta.q3,ymax=beta.whisker.max,linetype=sig),position=dodge) +
	geom_crossbar(aes(ymin=beta.q1,ymax=beta.q3,fill=cell,alpha=sig),position=dodge,width=0.5) +
	theme_classic() +
	scale_fill_manual(values=fill.colors['all cells']) +
	scale_alpha_manual(values=c(1,0.5)) +
	scale_linetype_manual(values=c(1,1)) +
	ylab(expression(beta)) +
	theme(
		axis.text.x = element_text(angle=75, vjust = 1, hjust=1),
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		legend.position='none'
	)

dodge = position_dodge(width=0.75)
p1 = ggplot(droplevels(subset(regulatory.results.stats,cell != 'all cells')),aes(dist,beta.median)) +
	geom_linerange(aes(ymin=beta.whisker.min,ymax=beta.q1,linetype=sig),position=dodge) +
	geom_linerange(aes(ymin=beta.q3,ymax=beta.whisker.max,linetype=sig),position=dodge) +
	geom_crossbar(aes(ymin=beta.q1,ymax=beta.q3,fill=cell,alpha=sig),position=dodge,width=0.5) +
	theme_classic() +
	scale_fill_manual(values=fill.colors[levels(droplevels(subset(regulatory.results.stats,cell != 'all cells'))$cell)]) +
	scale_alpha_manual(values=c(1,0.5)) +
	scale_linetype_manual(values=c(1,1)) +
	facet_wrap(~cell,ncol=3,scales='free_y') +
	ylab(expression(beta)) +
	theme(
		axis.text.x = element_text(angle=75, vjust = 1, hjust=1),
		axis.title.x = element_blank(),
		strip.background = element_blank(),
		legend.position='none'
	)

pdf(file='/dev/null')
p = egg::ggarrange(p0,p1,nrow=2,heights=c(1,4))
dev.off()

pdf(file=file.path('figures/final',paste0('regulatory_tssdist_metacell_beta_bysig.pdf')),useDingbats=FALSE,height=7,width=7)
p
dev.off()




dodge = position_dodge(width=0.75)
p0 = ggplot(droplevels(subset(regulatory.results.stats,cell == 'all cells')),aes(dist,sbet.median)) +
	geom_linerange(aes(ymin=sbet.whisker.min,ymax=sbet.q1,linetype=sig),position=dodge) +
	geom_linerange(aes(ymin=sbet.q3,ymax=sbet.whisker.max,linetype=sig),position=dodge) +
	geom_crossbar(aes(ymin=sbet.q1,ymax=sbet.q3,fill=cell,alpha=sig),position=dodge,width=0.5) +
	theme_classic() +
	scale_fill_manual(values=fill.colors['all cells']) +
	scale_alpha_manual(values=c(1,0.5)) +
	scale_linetype_manual(values=c(1,1)) +
	ylab(expression('standardized'~beta)) +
	theme(
		axis.text.x = element_text(angle=75, vjust = 1, hjust=1),
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		legend.position='none'
	)

dodge = position_dodge(width=0.75)
p1 = ggplot(droplevels(subset(regulatory.results.stats,cell != 'all cells')),aes(dist,sbet.median)) +
	geom_linerange(aes(ymin=sbet.whisker.min,ymax=sbet.q1,linetype=sig),position=dodge) +
	geom_linerange(aes(ymin=sbet.q3,ymax=sbet.whisker.max,linetype=sig),position=dodge) +
	geom_crossbar(aes(ymin=sbet.q1,ymax=sbet.q3,fill=cell,alpha=sig),position=dodge,width=0.5) +
	theme_classic() +
	scale_fill_manual(values=fill.colors[levels(droplevels(subset(regulatory.results.stats,cell != 'all cells'))$cell)]) +
	scale_alpha_manual(values=c(1,0.5)) +
	scale_linetype_manual(values=c(1,1)) +
	facet_wrap(~cell,ncol=3,scales='free_y') +
	ylab(expression('standardized'~beta)) +
	theme(
		axis.text.x = element_text(angle=75, vjust = 1, hjust=1),
		axis.title.x = element_blank(),
		strip.background = element_blank(),
		legend.position='none'
	)

pdf(file='/dev/null')
p = egg::ggarrange(p0,p1,nrow=2,heights=c(1,4))
dev.off()

pdf(file=file.path('figures/final',paste0('regulatory_tssdist_metacell_sbet_bysig.pdf')),useDingbats=FALSE,height=7,width=7)
p
dev.off()




# GLUE score against metacell standardized beta
p = ggplot(subset(regulatory.results,cell=='all cells'), aes(glue.score,meta.sbet)) +
	geom_point_rast(aes(color=factor(glue.qval < 0.05 & meta.qval < 0.05,levels=c('FALSE','TRUE'))),size=0.005,alpha=0.025) +
	geom_vline(xintercept=with(subset(regulatory.results,cell=='all cells'),mean(c(min(glue.score[glue.qval < 0.05]),max(glue.score[glue.qval > 0.05]))))) +
	geom_hline(yintercept=with(subset(regulatory.results,cell=='all cells'),mean(c(min(abs(meta.sbet)[meta.qval < 0.05]),max(abs(meta.sbet)[meta.qval > 0.05]))))) +
	geom_hline(yintercept=-with(subset(regulatory.results,cell=='all cells'),mean(c(min(abs(meta.sbet)[meta.qval < 0.05]),max(abs(meta.sbet)[meta.qval > 0.05]))))) +
#	geom_density_2d_filled(bins=10,alpha=0.5) +
#	geom_density_2d(bins=10,color='#ffffff') +
#	geom_smooth(method=lm) +
	scale_color_manual(values=c('#000000','#ff0000')) +
#	scale_alpha_manual(values=c(0.025,0.025)) +
#	scale_fill_manual(values=c('#ffffff00',viridis::viridis_pal(option='D',alpha=0.5)(9))) +
	coord_cartesian() +
	theme_classic(base_size=16) +
	theme(legend.position='none') +
	xlab('GLUE score') + 
	ylab(expression('std.'~beta))
ggsave(p,file=file.path('figures/final',paste0('regulatory_glue_score_metacell_sbet.pdf')),useDingbats=FALSE,width=7,height=4)



p0 = ggplot(droplevels(subset(regulatory.results,cell == 'all cells')),aes(tss_dist_label2,glue.score,fill=cell)) +
	geom_boxplot(outlier.shape=NA,width=0.25) +
	theme_classic(base_size=16) +
	# facet_wrap(~cell) +
	scale_fill_manual(values=fill.colors[levels(regulatory.results$cell)]) +
	# xlab('Genomic distance') +
	ylab('GLUE regulatory score') +
	theme(
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.text.x=element_text(angle=75, vjust = 1, hjust=1),
		legend.position='none',
		strip.background=element_blank())
ggsave(p0,file=file.path('figures/final',paste0('regulatory_tssdist_glue_scores_all.pdf')),useDingbats=FALSE,height=4,width=7)


p0 = ggplot(droplevels(subset(regulatory.results,cell == 'all cells')),aes(tss_dist_label2,glue.score,fill=factor(ifelse(meta.qval < 0.05,'metaTRUE','metaFALSE'),levels=c('metaFALSE','metaTRUE')))) +
	geom_boxplot(outlier.shape=NA,width=0.5,position='dodge') +
	theme_classic(base_size=16) +
	# facet_wrap(~cell) +
	scale_fill_manual(values=c('metaFALSE'='#1b9e77','metaTRUE'='#d95f02')) +
	# xlab('Genomic distance') +
	ylab('GLUE regulatory score') +
	theme(
		axis.title.x=element_blank(),
		#axis.title.y=element_blank(),
		axis.text.x=element_text(angle=45, vjust = 1, hjust=1),
		legend.position='none',
		strip.background=element_blank())
ggsave(p0,file=file.path('figures/final',paste0('regulatory_tssdist_glue_scores_all2.pdf')),useDingbats=FALSE,height=4,width=7)



p0 = ggplot(droplevels(subset(regulatory.results,cell == 'all cells')),aes(glue.score,color=factor(ifelse(meta.qval < 0.05,'metaTRUE','metaFALSE'),levels=c('metaFALSE','metaTRUE')),fill=factor(ifelse(meta.qval < 0.05,'metaTRUE','metaFALSE'),levels=c('metaFALSE','metaTRUE')))) +
	geom_density(alpha=0.5) +
	theme_classic(base_size=16) +
	coord_flip() +
	# facet_wrap(~cell) +
	scale_fill_manual(values=c('metaFALSE'='#1b9e77','metaTRUE'='#d95f02')) +
	scale_color_manual(values=c('metaFALSE'='#1b9e77','metaTRUE'='#d95f02')) +
	# xlab('Genomic distance') +
	ylab('GLUE regulatory score') +
	theme(
		axis.title.x=element_blank(),
		#axis.title.y=element_blank(),
		axis.text.x=element_text(angle=45, vjust = 1, hjust=1),
		legend.position='none',
		strip.background=element_blank())
ggsave(p0,file=file.path('figures/final',paste0('regulatory_tssdist_glue_scores_all_density.pdf')),useDingbats=FALSE,height=4,width=7)


foo = droplevels(subset(regulatory.results,cell == 'all cells'))
foo = do.call(rbind,lapply(split(foo,foo$tss_dist_label2),function(x) {
	data.frame(
		tss_dist_label2 = unique(x$tss_dist_label2),
		proportion = as.numeric(table(x$meta.qval < 0.05) / nrow(x)),
		metaSig = factor(c('metaFALSE','metaTRUE'),levels=c('metaFALSE','metaTRUE'))
	)
}))
p0 = ggplot(foo,aes(tss_dist_label2,proportion,fill=metaSig)) +
	geom_bar(stat='identity',position='stack') +
	theme_classic(base_size=16) +
	# facet_wrap(~cell) +
	scale_fill_manual(values=c('metaFALSE'='#1b9e77','metaTRUE'='#d95f02')) +
	# xlab('Genomic distance') +
	ylab('GLUE regulatory score') +
	theme(
		axis.title.x=element_blank(),
		#axis.title.y=element_blank(),
		axis.text.x=element_text(angle=45, vjust = 1, hjust=1),
		legend.position='none',
		strip.background=element_blank())
ggsave(p0,file=file.path('figures/final',paste0('regulatory_tssdist_glue_scores_all_bar.pdf')),useDingbats=FALSE,height=4,width=7)




# saveRDS(regulatory.results,file='rds/figures/regulatory_results.rds')






}














prefix = 'biccn'
atac.prefix = 'atac'
this.gene = 'ENSMMUG00000011333' # MBP

show.background = FALSE

suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrastr))
suppressMessages(library(ggforce))
suppressMessages(library(lemon))

file.colnames = c('peak','ensembl_gene_id','dist','gene_dir','glue.score','glue.pval','glue.qval','meta.beta','meta.pval','meta.qval')

cell.classes = scan(what='',sep='\n',file='stats/clusters/rna-final-cellclasses-levels.txt',quiet=TRUE)
this.cell.class = 'all cells'

genome.padding = 0.05

inc.file = file.path('stats/regulatory',paste0(prefix,'_',atac.prefix,'_all_peaks_inc.txt'))
dec.file = file.path('stats/regulatory',paste0(prefix,'_',atac.prefix,'_all_peaks_dec.txt'))
background.file = file.path('stats/regulatory/bed',paste0(prefix,'_',atac.prefix,'_all_peaks_genes_background.txt'))


results = rbind(
	read.delim(inc.file,header=FALSE,col.names=file.colnames),
	read.delim(dec.file,header=FALSE,col.names=file.colnames)
)

ens.genes = readRDS(file.path('checkpoints',paste0(ens.species,'_genes.rds')))
ens.exons = readRDS(file.path('checkpoints',paste0(ens.species,'_exons.rds')))
ens.tss = readRDS(file.path('checkpoints',paste0(ens.species,'_transcripts.rds')))
ens.utr = readRDS(file.path('checkpoints',paste0(ens.species,'_utr.rds')))

background.bed = read.delim(background.file,header=FALSE,col.names=c('peak','ensembl_gene_id'))

results.bed = data.frame(
	as.data.frame.matrix(do.call(rbind,strsplit(results$peak,'_'))),
	results
)
names(results.bed)[1:3] = c('chr','start','stop')
results.bed$start = as.numeric(results.bed$start)
results.bed$stop = as.numeric(results.bed$stop)

results.bed = merge(results.bed,ens.genes,by='ensembl_gene_id',all.x=TRUE,all.y=FALSE)

results.this = results.bed[results.bed$ensembl_gene_id == this.gene,]

peaks.this = results.this[,c('chr','start','stop','glue.score','glue.pval','meta.beta','meta.pval')]
background.prefix='_no_background'

genome.max = max(c(peaks.this$start,peaks.this$stop,results.this$start_position,results.this$end_position))
genome.min = min(c(peaks.this$start,peaks.this$stop,results.this$start_position,results.this$end_position))

genome.range = genome.max - genome.min

genome.chr = unique(results.this$chr)
genome.max = genome.max + genome.padding * genome.range
genome.min = genome.min - genome.padding * genome.range

gene.this = ens.genes[ens.genes$ensembl_gene_id == this.gene,]
tss.this = ens.tss[ens.tss$ensembl_gene_id == this.gene,]
exons.this = ens.exons[ens.exons$ensembl_gene_id == this.gene,]

tss.this$transcription_end_site = with(tss.this,ifelse(transcription_start_site == transcript_start, transcript_end, ifelse(transcription_start_site == transcript_end,transcript_start,NA)))

tss.this = tss.this[order(tss.this$transcript_start,tss.this$transcript_end),]

# background.this = subset(background.bed,chr == genome.chr & stop > genome.min & start < genome.max)
# background.this$peak = with(background.this,paste(chr,start,stop,sep='_'))

global.tss = if (gene.this$strand > 0) {
	min(tss.this$transcription_start_site)
} else {
	max(tss.this$transcription_start_site)
}

results.this$tss = global.tss
results.this$midpoint = with(results.this,start+(stop-start)/2)

results.this$direction = factor(results.this$meta.beta > 0,levels=c('TRUE','FALSE'),labels=c('upregulate','downregulate'))

links.this = do.call(rbind,lapply(split(results.this,results.this$peak),function(x) {
	this.peak = unique(x$peak)
	this.direction = unique(x$direction)
	this.start = min(c(x$midpoint,x$tss))
	this.end = max(c(x$midpoint,x$tss))
	data.frame(
		x = c(this.start,mean(c(this.start,this.end)),this.end),
		y = c(0,2*(-log10(x$meta.pval)),0),
		peak = this.peak,
		direction = this.direction
	)
}))

inferred.tss = max(tss.this$transcription_start_site)
results.this$tss2 = inferred.tss

links.that = do.call(rbind,lapply(split(results.this,results.this$peak),function(x) {
	this.peak = unique(x$peak)
	this.direction = unique(x$direction)
	this.start = min(c(x$midpoint,x$tss2))
	this.end = max(c(x$midpoint,x$tss2))
	data.frame(
		x = c(this.start,mean(c(this.start,this.end)),this.end),
		y = c(0,2*(-log10(x$meta.pval)),0),
		peak = this.peak,
		direction = this.direction
	)
}))


# exons.this = ens.exons[ens.exons$ensembl_gene_id == this.gene,]

scaling.factor = 20/max(c(10,links.this$y))

tss.this$y = if (nrow(tss.this) == 1) {
	-4/scaling.factor
} else {
	-13/scaling.factor + seq(0,10/scaling.factor,10/scaling.factor/(nrow(tss.this)-1))
}

exons.this = merge(tss.this,exons.this,by=c('ensembl_gene_id','ensembl_transcript_id'))

max.transcript.length = max(abs(tss.this$transcript_end - tss.this$transcript_start))
arrow.interval = round(max.transcript.length/10)
# Split tss.this into chunks
tss.arrow = do.call(rbind,lapply(1:nrow(tss.this),function(i) {
	x = tss.this[i,]
	meta = x[c('ensembl_gene_id','ensembl_transcript_id','y')]
	if (x$transcription_start_site > x$transcription_end_site) {
		new.start = seq(x$transcription_start_site,x$transcription_end_site,-arrow.interval)
		new.end = new.start - arrow.interval
		new.end[new.end < x$transcription_end_site] = x$transcription_end_site
	} else {
		new.start = seq(x$transcription_start_site,x$transcription_end_site,arrow.interval)
		new.end = new.start + arrow.interval
		new.end[new.end > x$transcription_end_site] = x$transcription_end_site
	}
	suppressWarnings(unique(data.frame(meta,transcription_start_site = new.start, transcription_end_site = new.end)))
}))

open.color        = '#ffffff'
open.line.color   = '#000000'
closed.color      = '#cccccc'
closed.line.color = '#000000'
link.inc.color    = '#377eb8'
link.dec.color    = '#e41a1c'
exon.color        = '#ff7f00'
exon.line.color   = NA
gene.color        = '#000000'

p0 = ggplot() +
	geom_bezier(data=links.this,aes(x=x,y=y,group=peak,color=direction),size=0.25) +
#	coord_flex_fixed(ratio = genome.range/50*scaling.factor,ylim=c(0,max(c(10,links.this$y/1.8))),left=capped_vertical(capped='bottom')) +
	coord_cartesian(ylim=c(0,max(links.this$y)/2)) +
	scale_x_continuous(limits=c(genome.min,genome.max)) +
	scale_color_manual(values=c(link.inc.color,link.dec.color)) +
	theme_classic() +
	theme(
		legend.position='none',
		legend.title=element_blank(),
		axis.text.x = element_blank(),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.line.x = element_blank()
	) +
	xlab(paste('Chromosome',gene.this$chromosome_name)) +
	ylab(expression(-log[10]~italic(p))) +
	ggtitle(ifelse(gene.this$ensembl_gene_id == gene.this$external_gene_name,paste0(this.cell.class,': ',gene.this$ensembl_gene_id),eval(parse(text=paste0('expression("',this.cell.class,':"~italic(',gene.this$external_gene_name,'))')))))

p00 = ggplot() +
	geom_bezier(data=links.that,aes(x=x,y=y,group=peak,color=direction),size=0.25) +
#	coord_flex_fixed(ratio = genome.range/50*scaling.factor,ylim=c(0,max(c(10,links.this$y/1.8))),left=capped_vertical(capped='bottom')) +
	coord_cartesian(ylim=c(0,max(links.this$y)/2)) +
	scale_x_continuous(limits=c(genome.min,genome.max)) +
	scale_color_manual(values=c(link.inc.color,link.dec.color)) +
	theme_classic() +
	theme(
		legend.position='none',
		legend.title=element_blank(),
		axis.text.x = element_blank(),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.line.x = element_blank()
	) +
	xlab(paste('Chromosome',gene.this$chromosome_name)) +
	ylab(expression(-log[10]~italic(p))) +
	ggtitle(ifelse(gene.this$ensembl_gene_id == gene.this$external_gene_name,paste0(this.cell.class,': ',gene.this$ensembl_gene_id),eval(parse(text=paste0('expression("',this.cell.class,':"~italic(',gene.this$external_gene_name,'))')))))




p1 = ggplot() +
	geom_rect(data=subset(peaks.this,!is.na(meta.pval)),aes(xmin=start,xmax=stop,ymin=-0.45/scaling.factor,ymax=-0.05/scaling.factor),fill=open.color,color=open.line.color,size=0.25)

if (any(is.na(peaks.this$meta.pval))) {
	p1 = p1 + geom_rect(data=subset(peaks.this,is.na(meta.pval)),aes(xmin=start,xmax=stop,ymin=-0.45/scaling.factor,ymax=-0.05/scaling.factor),fill=closed.color,color=closed.line.color,size=0.25)
}

p1 = p1 + 
#	geom_rect(data=unique(results.this[c('ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','strand')]),aes(xmin=start_position,xmax=end_position,ymin=-1.5/scaling.factor,ymax=-1.0/scaling.factor),fill=gene.color,color='black',size=0.25) +
	scale_x_continuous(limits=c(genome.min,genome.max)) +
	theme_classic() +
	theme(
		legend.position='none',
		axis.text = element_blank(),
		axis.title = element_blank(),
		axis.ticks = element_blank(),
		axis.line = element_blank()
	)

oli = subset(readRDS(file.path('rds',paste0(atac.prefix,'_all_marker_peaks_combined.rds'))),cluster=='oligodendrocytes' & peak %in% with(peaks.this,paste(chr,start,stop,sep='_')))

oli = data.frame(oli,within(as.data.frame.matrix(do.call(rbind,strsplit(oli$peak,'_'))),{
	chr = V1
	start = as.integer(V2)
	stop = as.integer(V3)
})[,c('chr','start','stop')])

oli = within(oli,{
	logreg_ymin = ifelse(logreg_score < 0,logreg_score,0)
	logreg_ymax = ifelse(logreg_score < 0,0,logreg_score)
	logFC_ymin = ifelse(logfoldchanges < 0,logfoldchanges,0)
	logFC_ymax = ifelse(logfoldchanges < 0,0,logfoldchanges)
	mid = start + (stop - start) / 2
})

p2 = ggplot() +
	geom_segment(data=oli,aes(x=mid,xend=mid,y=0,yend=logfoldchanges,color=pts)) +
	scale_x_continuous(limits=c(genome.min,genome.max)) +
	scale_y_continuous(breaks=seq(-2,2,2)) +
	scale_color_viridis() +
	theme_classic() +
	theme(
		legend.position='none',
		legend.title=element_blank(),
		axis.text.x = element_blank(),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.line.x = element_blank()
	) + 
	ylab(expression(log[2]~'FC'))

p3 = ggplot() +
	geom_segment(data=oli,aes(x=mid,xend=mid,y=0,yend=logreg_score,color=pts)) +
	scale_x_continuous(limits=c(genome.min,genome.max)) +
	scale_y_continuous(breaks=c(0,0.2)) +
	scale_color_viridis() +
	theme_classic() +
	theme(
		legend.position='right',
		axis.text.x = element_blank(),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.line.x = element_blank()
	) + 
	ylab(expression('LR score'))

p4 = ggplot() +
	geom_segment(data=tss.arrow,aes(x=transcription_start_site,xend=transcription_end_site,y=y,yend=y),size=0.25,arrow=arrow(length=unit(0.05,'inches'))) +
	geom_rect(data=exons.this,aes(xmin=exon_chrom_start,xmax=exon_chrom_end,ymin=y-0.25/scaling.factor,ymax=y+0.25/scaling.factor),fill=exon.color,color=exon.line.color,size=0.25) +
#	geom_rect(data=unique(results.this[c('ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','strand')]),aes(xmin=start_position,xmax=end_position,ymin=-1.5/scaling.factor,ymax=-1.0/scaling.factor),fill=gene.color,color='black',size=0.25) +
	geom_segment(data=unique(results.this[c('ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','strand')]),aes(x=start_position,xend=end_position,y=-1.5/scaling.factor,yend=-1.5/scaling.factor),color=gene.color,size=0.5,arrow=arrow(length=unit(0.05,'inches'),ends='both',angle=90)) +
	scale_x_continuous(limits=c(genome.min,genome.max)) +
	theme_classic() +
	theme(
		legend.position='none',
		axis.text.y = element_blank(),
		axis.title.y = element_blank(),
		axis.ticks.y = element_blank(),
		axis.line.y = element_blank()
	) +
	xlab(paste('Chromosome',gene.this$chromosome_name))



# region = '18_2923185_3094144'
oli.coverage = read.delim('coverage/class9_18_2923185_3094144_bins.cov.txt',header=FALSE)


library(data.table)

rna.meta = fread('datasets/processed/rna-pseudobulk-meta.txt.gz')

umi.totals = tapply(rna.meta$umi,rna.meta$cellclass,sum)
umi.totals.join = as.integer(umi.totals)
names(umi.totals.join) = names(umi.totals)

cell.rna.coverage = do.call(rbind,lapply(1:19,function(i) {
	this.coverage = read.delim(file.path('coverage',paste0('rna_class',i,'_18_2923185_3094144_bins.cov.txt')),header=FALSE)
	this.coverage$cell_class = factor(cell.classes[i],levels=cell.classes)
	this.coverage
}))
levels(cell.rna.coverage$cell_class) = gsub('serotinergic','serotonergic',levels(cell.rna.coverage$cell_class))
levels(cell.rna.coverage$cell_class) = gsub('neuron$','neurons',levels(cell.rna.coverage$cell_class))

cell.rna.coverage$cell_class = factor(
	cell.rna.coverage$cell_class,
	levels=cell.levels
)

cell.rna.coverage.plot = droplevels(subset(cell.rna.coverage,cell_class %in% c(
'excitatory neurons',
'cerebellar neurons',
'inhibitory neurons',
'basket cells',
'medium spiny neurons',
'dopaminergic neurons',
'serotonergic neurons',
'astrocytes',
'oligodendrocytes',
'oligodendrocyte precursor cells',
'radial glial cells',
'ependymal cells',
'microglia',
'vascular cells',
'mesenchymal stem cells'
)))
cell.rna.coverage.plot$total_umi = umi.totals.join[as.character(cell.rna.coverage.plot$cell_class)]

cell.rna.coverage.plot$cpm = with(cell.rna.coverage.plot,V4 / total_umi * 1e6)

cell.rna.coverage.plot$oli = factor(cell.rna.coverage.plot$cell_class == 'oligodendrocytes',levels=c('TRUE','FALSE'))

p5 = ggplot() +
	geom_path(data=cell.rna.coverage.plot,aes(x=V2+50,y=cpm,color=cell_class,alpha=oli),size=0.5) +
	scale_x_continuous(limits=c(genome.min,genome.max)) +
	facet_wrap(~oli,ncol=1) +
	scale_color_manual(values=colors[levels(cell.rna.coverage.plot$cell_class)]) +
	scale_alpha_manual(values=c(1,0.5)) +
	scale_y_continuous(breaks=c(3000),limits=c(0,ceiling(max(cell.rna.coverage.plot$cpm)))) +
	# scale_color_viridis() +
	theme_classic() +
	theme(
		legend.position='none',
		axis.text.x = element_blank(),
		axis.title = element_blank(),
		axis.ticks.x = element_blank(),
		axis.line.x = element_blank(),
		# axis.text.y = element_blank(),
		# axis.ticks.y = element_blank(),
		strip.background=element_blank(),
		strip.text = element_blank()
	) + 
	ylab(expression('RNA coverage'))


# ATAC

umi.totals.join = unlist(lapply(c(1:3,6:13,15),function(i) {
	out = sum(scan(paste0('bed_region/umi/class',i,'_umi.txt')))
	names(out) = cell.classes[i]
	out
}))

cell.atac.coverage = do.call(rbind,lapply(c(1:3,6:13,15),function(i) {
	this.coverage = read.delim(file.path('coverage',paste0('atac_class',i,'_18_2923185_3094144_bins.cov.txt')),header=FALSE)
	this.coverage$cell_class = factor(cell.classes[i],levels=cell.classes)
	this.coverage
}))
cell.atac.coverage$cell_class = factor(
	cell.atac.coverage$cell_class,
	levels=cell.levels
)[,drop=TRUE]

cell.atac.coverage.plot = droplevels(subset(cell.atac.coverage,TRUE))

cell.atac.coverage.plot$cell_class = factor(cell.atac.coverage.plot$cell_class,levels=c(
'oligodendrocytes',
'oligodendrocyte precursor cells',
'microglia',
'vascular cells',
'ependymal cells',
'astrocytes',
'radial glial cells',
'medium spiny neurons',
'basket cells',
'cerebellar neurons',
'inhibitory neurons',
'excitatory neurons'
))

cell.atac.coverage.plot$total_umi = umi.totals.join[as.character(cell.atac.coverage.plot$cell_class)]

cell.atac.coverage.plot$cpm = with(cell.atac.coverage.plot,V4 / total_umi * 1e6)

cell.atac.coverage.plot$cell_facet = as.character(cell.atac.coverage.plot$cell_class)
cell.atac.coverage.plot$cell_facet[!cell.atac.coverage.plot$cell_facet %in% c('oligodendrocytes','oligodendrocyte precursor cells','microglia')] = 'other'

cell.atac.coverage.plot$cell_facet = factor(cell.atac.coverage.plot$cell_facet,levels=c(
'oligodendrocytes',
'oligodendrocyte precursor cells',
'microglia',
'other'
))

p6 = ggplot() +
	geom_path(data=cell.atac.coverage.plot,aes(x=V2+50,y=cpm,color=cell_class,alpha=factor(cell_facet=='other',levels=c('FALSE','TRUE'))),size=0.5) +
	scale_x_continuous(limits=c(genome.min,genome.max)) +
	facet_wrap(~cell_facet,ncol=1) +
	scale_color_manual(values=colors[levels(cell.atac.coverage.plot$cell_class)]) +
	scale_alpha_manual(values=c(1,0.5)) +
	scale_y_continuous(breaks=c(50)) +
	# scale_color_viridis() +
	theme_classic() +
	theme(
		legend.position='none',
		axis.text.x = element_blank(),
		axis.title = element_blank(),
		axis.ticks.x = element_blank(),
		axis.line.x = element_blank(),
		# axis.text.y = element_blank(),
		# axis.ticks.y = element_blank(),
		strip.background=element_blank(),
		strip.text = element_blank()
	) + 
	ylab(expression('ATAC coverage'))





# RNA coverage, zoomed in
p7 = ggplot() +
	geom_path(data=within(droplevels(subset(cell.rna.coverage.plot,TRUE | !(cell_class %in% c('oligodendrocytes')))),{cell_facet=factor(ifelse(cell_class == 'oligodendrocytes','oligodendrocytes',ifelse(cell_class == 'microglia','microglia',ifelse(cell_class == 'oligodendrocyte precursor cells','oligodendrocyte precursor cells','other'))),levels=c('oligodendrocytes','oligodendrocyte precursor cells','microglia','other'))}),aes(x=V2+50,y=cpm,color=cell_class,alpha=factor(cell_class %in% c('oligodendrocytes','microglia','oligodendrocyte precursor cells'),levels=c('TRUE','FALSE'))),size=0.5) +
	scale_x_continuous(limits=c(genome.min,genome.max)) +
	facet_wrap(~cell_facet,ncol=1) +
	scale_color_manual(values=colors[levels(droplevels(subset(cell.rna.coverage.plot, TRUE | !(cell_class %in% c('oligodendrocytes'))))$cell_class)]) +
	scale_alpha_manual(values=c(1,0.5)) +
	scale_y_continuous(limits=c(0,max(droplevels(subset(cell.rna.coverage.plot,TRUE | !(cell_class %in% c('oligodendrocytes'))))$cpm)),breaks=200) +
	coord_cartesian(ylim=c(0,200)) +
	# scale_color_viridis() +
	theme_classic() +
	theme(
		legend.position='none',
		# axis.text.x = element_blank(),
		axis.title = element_blank(),
		# axis.ticks.x = element_blank(),
		# axis.line.x = element_blank(),
		# axis.text.y = element_blank(),
		# axis.ticks.y = element_blank(),
		strip.background=element_blank(),
		strip.text = element_blank()
	) + 
	ylab(expression('RNA coverage'))
ggsave(p7,file=file.path('figures/final',paste0('regulatory_mbp_rna_rest.pdf')),useDingbats=FALSE,height=5,width=7)

# ggarrange(p0,p1,p2,p3,p4,nrow=5,heights=c(1,0.1,0.2,0.2,2))

egg::ggarrange(p0,p1,p2,p5,p6,p4,nrow=6,heights=c(1.5,0.1,0.5,0.5,1.5,0.75))

pdf(file='/dev/null')
p = egg::ggarrange(p0,p1,p2,p5,p6,p4,nrow=6,heights=c(1,0.075,0.5,0.5,1,0.75))
dev.off()


pdf(file=file.path('figures/final',paste0('regulatory_mbp_links.pdf')),useDingbats=FALSE,height=5,width=7)
p
dev.off()


pdf(file='/dev/null')
p = egg::ggarrange(p00,p1,p2,p5,p6,p4,nrow=6,heights=c(1,0.075,0.5,0.5,1,0.75))
dev.off()


pdf(file=file.path('figures/final',paste0('regulatory_mbp_links2.pdf')),useDingbats=FALSE,height=5,width=7)
p
dev.off()


dir.create('figures/regulatory/tracks',showWarnings=FALSE)




}




# LD score
summary.files = file.path('stats/ldsc/all_peaks/scores',list.files('stats/ldsc/all_peaks/scores',pattern='_summary.txt'))

summary.stats = do.call(rbind,lapply(summary.files,read.delim,row.names=1))

summary.stats = subset(summary.stats,Cell_class!='radial glial cells')

summary.stats = do.call(rbind,lapply(split(summary.stats,summary.stats$Cell_class),function(x) within(x,{Enrichment_padj = p.adjust(Enrichment_p,'fdr')})))

write.table(summary.stats,file='tables/final/ldsc_summary_stats_allpeaks.txt',sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)


ld.stats = read.delim('tables/final/ldsc_summary_stats_allpeaks.txt')

phenotype.metadata = read.delim('data/phenotype_metadata.txt')

phenotype.metadata$phenotype_label = with(phenotype.metadata,paste0(citation,': ',trait))

ld.stats = merge(ld.stats,phenotype.metadata,by.x='Phenotype',by.y='phenotype',all.x=TRUE)


ld.stats$Phenotype = factor(ld.stats$Phenotype,levels=phenotype.metadata$phenotype)
ld.stats$Cell_class = factor(ld.stats$Cell_class,levels=intersect(cell.levels,unique(ld.stats$Cell_class)))

ld.stats$trait = factor(ld.stats$trait,levels=unique(phenotype.metadata$trait))
ld.stats$phenotype_label = factor(ld.stats$phenotype_label,levels=phenotype.metadata$phenotype_label)

ld.stats$or = ld.stats$Enrichment
ld.stats$or[!is.na(ld.stats$or) & ld.stats$or < 0] = NA
ld.stats$or[!is.na(ld.stats$or) & ld.stats$Enrichment_padj > 0.05] = NA

# p = ggplot(tf.keep,aes(motif_name,cell_class,fill=log2OR)) +
# 	geom_tile() +
# 	coord_fixed() +
# 	# scale_fill_viridis() +
# 	scale_fill_gradient2(
# 		high='#2066ac',
# 		mid='#ffffff',
# 		low='#b2182b',
# 		midpoint=0,
# 		breaks=seq(-1,1,1),
# 		# limits=c(-max(abs(tf.keep$log2OR)),max(abs(tf.keep$log2OR))),
# 		limits=c(-1.5,1.5),
# 		oob=scales::squish,
# 		name=expression(log[2]~'enrichment')
# 	) +
# 	scale_y_discrete(limits=rev) +
# 	theme_classic(base_size=8) +
# 	theme(
# 		axis.title=element_blank(),
# 		axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5)
# 	)
# ggsave(p,file='figures/final/tf_jaspar_unique_peaks_enrichment.pdf',width=7,height=4,useDingbats=FALSE)

p = ggplot(droplevels(subset(ld.stats,!is.na(or))),aes(phenotype_label,Cell_class,fill=log2(or))) +
	geom_tile() +
	coord_equal() +
	scale_fill_viridis(
		limits=with(ld.stats,c(0,mean(log2(or),na.rm=T) + 2 * sd(log2(or),na.rm=T))),
		oob=scales::squish,
		name=expression(log[2]~'OR')
	) +
	scale_x_discrete(position='bottom') +
 	scale_y_discrete(limits=rev,position='right') +
	theme_classic(base_size=8) +
 	theme(
 		axis.title=element_blank(),
 		axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
 		axis.text.y=element_text(angle=0,hjust=0,vjust=0.5),
 		legend.position = 'left'
 	)
ggsave(p,file='figures/final/ld-enrich01.pdf',height=7,width=7)


p1 = ggplot(droplevels(subset(ld.stats,!is.na(or))),aes(phenotype_label,Cell_class,fill=log2(or))) +
	geom_tile() +
	coord_equal() +
	scale_fill_viridis(
		limits=with(ld.stats,c(0,mean(log2(or),na.rm=T) + 2 * sd(log2(or),na.rm=T))),
		oob=scales::squish,
		name=expression(log[2]~'OR')
	) +
	scale_x_discrete(position='top') +
 	scale_y_discrete(limits=rev,position='right') +
	theme_classic(base_size=8) +
 	theme(
 		axis.title=element_blank(),
 		axis.text.x=element_text(angle=45,hjust=0,vjust=0.5),
 		axis.text.y=element_text(angle=0,hjust=0,vjust=0.5),
 		legend.position = 'left'
 	)
ggsave(p1,file='figures/final/ld-enrich02.pdf',height=7,width=7)




# Plot QC stats

# Bring in raw RNA files

# 
library(data.table)

rna.metadata = do.call(rbind,mclapply(file.path('mm',list.files('mm',pattern='rna-NSM[0-9]{3}-metadata.txt.gz')),fread,mc.cores=8))
atac.metadata = do.call(rbind,mclapply(file.path('mm',list.files('mm',pattern='atac-NSM[0-9]{3}-metadata.txt.gz')),fread,mc.cores=8))

saveRDS(rna.metadata,file='checkpoints/rna_raw_qc_metadata.rds')
saveRDS(atac.metadata,file='checkpoints/atac_raw_qc_metadata.rds')

p = ggplot(rna.metadata,aes(perc_mitochondrial_umis)) +
	geom_histogram(bins=100) +
	geom_vline(xintercept=5,color='red') +
	scale_y_continuous(trans='log10',breaks=10^seq(1,5,2),labels=gsub(' ','',seq(1,5,2))) +
	theme_classic(base_size=16) +
	xlab('perc. mitochondrial UMIs') +
	ylab(expression(log[10]~'count'))
ggsave(p,file='figures/final/rna_qc_mt.pdf',height=3,width=7,useDingbats=FALSE)

p = ggplot(rna.metadata,aes(n.umi)) +
	geom_histogram(bins=100) +
	geom_vline(xintercept=100,color='red') +
	scale_x_continuous(trans='log10',breaks=10^seq(2,6,1),labels=gsub(' ','',seq(2,6,1))) +
	scale_y_continuous(trans='log10',breaks=10^seq(1,5,2),labels=gsub(' ','',seq(1,5,2))) +
	theme_classic(base_size=16) +
	xlab(expression(log[10]~'nUMIs')) +
	ylab(expression(log[10]~'count'))
ggsave(p,file='figures/final/rna_qc_umi.pdf',height=3,width=7,useDingbats=FALSE)



p = ggplot(atac.metadata,aes(n.umi)) +
	geom_histogram(bins=100) +
	geom_vline(xintercept=1000,color='red') +
	geom_vline(xintercept=100000,color='red') +
	scale_x_continuous(trans='log10',breaks=10^seq(2,6,1),labels=gsub(' ','',seq(2,6,1))) +
	scale_y_continuous(trans='log10',breaks=10^seq(1,5,2),labels=gsub(' ','',seq(1,5,2))) +
	theme_classic(base_size=16) +
	xlab(expression(log[10]~'nUMIs')) +
	ylab(expression(log[10]~'count'))
ggsave(p,file='figures/final/atac_qc_umi.pdf',height=3,width=7,useDingbats=FALSE)

p = ggplot(atac.metadata,aes(FRIP)) +
	geom_histogram(bins=100) +
	geom_vline(xintercept=0.3,color='red') +
	scale_x_continuous(breaks=seq(0.1,0.7,0.2)) +
	scale_y_continuous(breaks=seq(0,1e5,5e4),labels=gsub(' ','',format(seq(0,1e5,5e4),scientific=FALSE,big.mark=''))) +
	theme_classic(base_size=16) +
	xlab(expression('FRIP')) +
	# ylab(expression(log[10]~'count'))
	ylab(expression('count'))
ggsave(p,file='figures/final/atac_qc_frip.pdf',height=3,width=7,useDingbats=FALSE)



# Bring in doublet-clusters

atac.doublet.umap = readRDS(file.path('umap',paste0('atac','-scanpy-all.rds')))
rna.doublet.umap = readRDS(file.path('umap',paste0('rna','-scanpy-all.rds')))

set.seed(42)
rna.doublet.umap$order = sample(1:nrow(rna.doublet.umap))
set.seed(42)
atac.doublet.umap$order = sample(1:nrow(atac.doublet.umap))

rna.doublet.umap = rna.doublet.umap[order(rna.doublet.umap$order),]
atac.doublet.umap = atac.doublet.umap[order(atac.doublet.umap$order),]

p = ggplot(data=rna.doublet.umap,aes_string('umap.1','umap.2',color='doublet_score')) +
	geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
	# scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
	scale_color_viridis(limits=c(0,with(rna.doublet.umap,mean(doublet_score) + 2*sd(doublet_score))),oob=scales::squish) +
	coord_equal() +
	theme_void(base_size=12) +
	theme(legend.position='right') +
	guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=2))
ggsave(p,file='figures/final/rna_qc_doublet_umap_doubletscore.pdf',height=7,width=7,useDingbats=FALSE)



rna.doublet.umap$doublet_class = factor(with(rna.doublet.umap,ifelse(doublet_cluster,'doublet cluster',ifelse(doublet_call == 'singlet','singlet','doublet cell'))),
	levels = c('singlet','doublet cluster','doublet cell'))
atac.doublet.umap$doublet_class = factor(with(atac.doublet.umap,ifelse(doublet_cluster,'doublet cluster',ifelse(doublet_call == 'singlet','singlet','doublet cell'))),
	levels = c('singlet','doublet cluster','doublet cell'))


p = ggplot(data=rna.doublet.umap,aes_string('umap.1','umap.2',color='doublet_score')) +
	geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
	# scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
	# scale_color_viridis(name='Scrublet score',limits=c(0,with(rna.doublet.umap,mean(doublet_score) + 2*sd(doublet_score))),oob=scales::squish) +
	scale_color_viridis(name='Scrublet score',limits=c(0,0.2),oob=scales::squish) +
	coord_equal() +
	theme_void(base_size=12) +
	theme(legend.position='right') +
	guides()
ggsave(p,file='figures/final/rna_qc_doublet_umap_doubletscore.pdf',height=7,width=7,useDingbats=FALSE)
p = ggplot(data=atac.doublet.umap,aes_string('umap.1','umap.2',color='doublet_score')) +
	geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
	# scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
	# scale_color_viridis(name='Scrublet score',limits=c(0,with(rna.doublet.umap,mean(doublet_score) + 2*sd(doublet_score))),oob=scales::squish) +
	scale_color_viridis(name='Scrublet score',limits=c(0,0.2),oob=scales::squish) +
	coord_equal() +
	theme_void(base_size=12) +
	theme(legend.position='right') +
	guides()
ggsave(p,file='figures/final/atac_qc_doublet_umap_doubletscore.pdf',height=7,width=7,useDingbats=FALSE)



p = ggplot(data=rna.doublet.umap,aes_string('umap.1','umap.2',color='doublet_class')) +
	geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
	# scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
	# scale_color_viridis(name='Scrublet score',limits=c(0,with(rna.doublet.umap,mean(doublet_score) + 2*sd(doublet_score))),oob=scales::squish) +
	scale_color_manual(values=c('#cccccc','#377eb8','#e41a1c')) +
	coord_equal() +
	theme_void(base_size=12) +
	theme(legend.position='right') +
	guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
ggsave(p,file='figures/final/rna_qc_doublet_umap_doubletcalls.pdf',height=7,width=7,useDingbats=FALSE)
p = ggplot(data=atac.doublet.umap,aes_string('umap.1','umap.2',color='doublet_class')) +
	geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
	# scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
	# scale_color_viridis(name='Scrublet score',limits=c(0,with(rna.doublet.umap,mean(doublet_score) + 2*sd(doublet_score))),oob=scales::squish) +
	scale_color_manual(values=c('#cccccc','#377eb8','#e41a1c')) +
	coord_equal() +
	theme_void(base_size=12) +
	theme(legend.position='right') +
	guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
ggsave(p,file='figures/final/atac_qc_doublet_umap_doubletcalls.pdf',height=7,width=7,useDingbats=FALSE)


# Visualize UMI and doublet score on final UMAP

p = ggplot(data=this.umap[with(this.umap,keep & modality == 'RNA'),],aes_string('rna_umap.1','rna_umap.2',color='rna_doublet_score')) +
	geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
	# scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
	# scale_color_viridis(name='Scrublet score',limits=c(0,with(rna.doublet.umap,mean(doublet_score) + 2*sd(doublet_score))),oob=scales::squish) +
	scale_color_viridis(name='Scrublet score',limits=c(0,0.2),oob=scales::squish) +
	coord_equal() +
	theme_void(base_size=12) +
	theme(legend.position='right') +
	guides()
ggsave(p,file='figures/final/rna_qc_doublet_umap_final_doubletscore.pdf',height=7,width=7,useDingbats=FALSE)
p = ggplot(data=this.umap[with(this.umap,keep & modality == 'ATAC'),],aes_string('atac_umap.1','atac_umap.2',color='atac_doublet_score')) +
	geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
	# scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
	# scale_color_viridis(name='Scrublet score',limits=c(0,with(rna.doublet.umap,mean(doublet_score) + 2*sd(doublet_score))),oob=scales::squish) +
	scale_color_viridis(name='Scrublet score',limits=c(0,0.2),oob=scales::squish) +
	coord_equal() +
	theme_void(base_size=12) +
	theme(legend.position='right') +
	guides()
ggsave(p,file='figures/final/atac_qc_doublet_umap_final_doubletscore.pdf',height=7,width=7,useDingbats=FALSE)

p = ggplot(data=this.umap[with(this.umap,keep & modality == 'RNA'),],aes_string('rna_umap.1','rna_umap.2',color='log10(rna_umi)')) +
	geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
	# scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
	# scale_color_viridis(name='Scrublet score',limits=c(0,with(rna.doublet.umap,mean(doublet_score) + 2*sd(doublet_score))),oob=scales::squish) +
	scale_color_viridis(name=expression(log[10]~UMI)) +
	coord_equal() +
	theme_void(base_size=12) +
	theme(legend.position='right') +
	guides()
ggsave(p,file='figures/final/rna_qc_doublet_umap_final_umi.pdf',height=7,width=7,useDingbats=FALSE)
p = ggplot(data=this.umap[with(this.umap,keep & modality == 'ATAC'),],aes_string('atac_umap.1','atac_umap.2',color='log10(atac_umi)')) +
	geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
	# scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
	# scale_color_viridis(name='Scrublet score',limits=c(0,with(rna.doublet.umap,mean(doublet_score) + 2*sd(doublet_score))),oob=scales::squish) +
	scale_color_viridis(name=expression(log[10]~UMI)) +
	coord_equal() +
	theme_void(base_size=12) +
	theme(legend.position='right') +
	guides()
ggsave(p,file='figures/final/atac_qc_doublet_umap_final_umi.pdf',height=7,width=7,useDingbats=FALSE)

p = ggplot(data=this.umap[with(this.umap,keep & modality == 'ATAC'),],aes_string('atac_umap.1','atac_umap.2',color='atac_FRIP')) +
	geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
	# scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
	# scale_color_viridis(name='Scrublet score',limits=c(0,with(rna.doublet.umap,mean(doublet_score) + 2*sd(doublet_score))),oob=scales::squish) +
	scale_color_viridis(name=expression('FRIP')) +
	coord_equal() +
	theme_void(base_size=12) +
	theme(legend.position='right') +
	guides()
ggsave(p,file='figures/final/atac_qc_doublet_umap_final_frip.pdf',height=7,width=7,useDingbats=FALSE)


# Visualize UMAPs split by sex, hemisphere






# p = ggplot(data=this.umap[with(this.umap,keep & modality == 'RNA' & sequencing_run_id != 'Snyder-Mackler_RNA3-036'),],aes_string('rna_umap.1','rna_umap.2',color='cell_class')) +
# 	geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
# 	# scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
# 	# scale_color_viridis(name='Scrublet score',limits=c(0,with(rna.doublet.umap,mean(doublet_score) + 2*sd(doublet_score))),oob=scales::squish) +
# 	# scale_color_manual(name='Sex',values=c('#1b9e77','#7570b3')) +
# 	facet_wrap(~sex,nrow=1) +
# 	scale_color_manual(values=colors[levels(this.umap$cell_class)]) +
# 	coord_equal() +
# 	theme_void(base_size=12) +
# 	theme(legend.position='none',legend.title=element_blank()) +
# 	guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
# ggsave(p,file='figures/final/rna_exploratory_umap_sex.pdf',height=7,width=7,useDingbats=FALSE)
# 
# p = ggplot(data=this.umap[with(this.umap,keep & modality == 'ATAC'),],aes_string('atac_umap.1','atac_umap.2',color='cell_class')) +
# 	geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
# 	# scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
# 	# scale_color_viridis(name='Scrublet score',limits=c(0,with(rna.doublet.umap,mean(doublet_score) + 2*sd(doublet_score))),oob=scales::squish) +
# 	# scale_color_manual(name='Sex',values=c('#1b9e77','#7570b3')) +
# 	facet_wrap(~sex,nrow=1) +
# 	scale_color_manual(values=colors[levels(this.umap$cell_class)]) +
# 	coord_equal() +
# 	theme_void(base_size=12) +
# 	theme(legend.position='none',legend.title=element_blank()) +
# 	guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
# ggsave(p,file='figures/final/atac_exploratory_umap_sex.pdf',height=7,width=7,useDingbats=FALSE)





p = ggplot(data=within(this.umap[with(this.umap,keep & modality == 'RNA' & sequencing_run_id != 'Snyder-Mackler_RNA3-036'),],{sex = factor(sex,levels=c('F','M'),labels=c('female','male'))}),aes_string('rna_umap.1','rna_umap.2',color='cell_class')) +
	geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
	# scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
	# scale_color_viridis(name='Scrublet score',limits=c(0,with(rna.doublet.umap,mean(doublet_score) + 2*sd(doublet_score))),oob=scales::squish) +
	# scale_color_manual(name='Sex',values=c('#1b9e77','#7570b3')) +
	facet_wrap(~sex,nrow=1) +
	scale_color_manual(values=colors[levels(this.umap$cell_class)]) +
	coord_equal() +
	theme_void(base_size=12) +
	theme(legend.position='none',legend.title=element_blank()) +
	guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
ggsave(p,file='figures/final/rna_exploratory_umap_sex.pdf',height=7,width=7,useDingbats=FALSE)

p = ggplot(data=within(this.umap[with(this.umap,keep & modality == 'ATAC'),],{sex = factor(sex,levels=c('F','M'),labels=c('female','male'))}),aes_string('atac_umap.1','atac_umap.2',color='cell_class')) +
	geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
	# scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
	# scale_color_viridis(name='Scrublet score',limits=c(0,with(rna.doublet.umap,mean(doublet_score) + 2*sd(doublet_score))),oob=scales::squish) +
	# scale_color_manual(name='Sex',values=c('#1b9e77','#7570b3')) +
	facet_wrap(~sex,nrow=1) +
	scale_color_manual(values=colors[levels(this.umap$cell_class)]) +
	coord_equal() +
	theme_void(base_size=12) +
	theme(legend.position='none',legend.title=element_blank()) +
	guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
ggsave(p,file='figures/final/atac_exploratory_umap_sex.pdf',height=7,width=7,useDingbats=FALSE)




p = ggplot(data=within(this.umap[with(this.umap,keep & modality == 'RNA' & !is.na(hemisphere) & sequencing_run_id != 'Snyder-Mackler_RNA3-036'),],{hemisphere = factor(hemisphere,levels=c('L','R'),labels=c('left','right'))}),aes_string('rna_umap.1','rna_umap.2',color='cell_class')) +
	geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
	# scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
	# scale_color_viridis(name='Scrublet score',limits=c(0,with(rna.doublet.umap,mean(doublet_score) + 2*sd(doublet_score))),oob=scales::squish) +
	# scale_color_manual(name='Sex',values=c('#1b9e77','#7570b3')) +
	facet_wrap(~hemisphere,nrow=1) +
	scale_color_manual(values=colors[levels(this.umap$cell_class)]) +
	coord_equal() +
	theme_void(base_size=12) +
	theme(legend.position='none',legend.title=element_blank()) +
	guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
ggsave(p,file='figures/final/rna_exploratory_umap_hemisphere.pdf',height=7,width=7,useDingbats=FALSE)

p = ggplot(data=within(this.umap[with(this.umap,keep & modality == 'ATAC' & !is.na(hemisphere)),],{hemisphere = factor(hemisphere,levels=c('L','R'),labels=c('left','right'))}),aes_string('atac_umap.1','atac_umap.2',color='cell_class')) +
	geom_point_rast(size=0.1,shape=19,alpha=0.05,dev='ragg_png',raster.dpi=300,na.rm=FALSE) +
	# scale_color_manual(values=colors[levels(shuffle.umap$rna_cell)],name='Cell',na.translate=TRUE,na.value=rgb(233,233,233,0.01,maxColorValue=255)) +
	# scale_color_viridis(name='Scrublet score',limits=c(0,with(rna.doublet.umap,mean(doublet_score) + 2*sd(doublet_score))),oob=scales::squish) +
	# scale_color_manual(name='Sex',values=c('#1b9e77','#7570b3')) +
	facet_wrap(~hemisphere,nrow=1) +
	scale_color_manual(values=colors[levels(this.umap$cell_class)]) +
	coord_equal() +
	theme_void(base_size=12) +
	theme(legend.position='none',legend.title=element_blank()) +
	guides(color=guide_legend(override.aes = list(size=2,alpha=1),ncol=1))
ggsave(p,file='figures/final/atac_exploratory_umap_hemisphere.pdf',height=7,width=7,useDingbats=FALSE)




# Plot hypothetical logistic regression
p = ggplot(data.frame(ge=c(1,0.45,0.78,0.87,0.92,0,0.14,0.4,0.6,0.1),ca=c(rep(1,5),rep(0,5))),aes(ge,ca)) + geom_point(size=3) + geom_smooth(method='glm',color='blue',se=FALSE, method.args = list(family=binomial)) + theme_classic() + scale_x_continuous() + theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank()) + coord_fixed(ratio=1,xlim=c(0,1),ylim=c(0,1))
ggsave(p,file='figures/final/logistic_regression_cartoon.pdf',height=3,width=3,useDingbats=FALSE)


p = ggplot(data.frame(ge=c(0.9,0.6,0.18,0.39,0.31,0.9,0.96,0.05,0.8,0.37),ca=c(rep(1,5),rep(0,5))),aes(ge,ca)) + geom_point(size=3) +
#	geom_smooth(method='glm',color='blue',se=FALSE, method.args = list(family=binomial)) +
	theme_classic() + scale_x_continuous() + theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank()) + coord_fixed(ratio=1,xlim=c(0,1),ylim=c(0,1))
ggsave(p,file='figures/final/logistic_regression_bad_cartoon.pdf',height=3,width=3,useDingbats=FALSE)


d = as.data.frame.matrix(cbind(
c(1:5,1:5,1:5,1:5),
c(rep(0,5),rep(1,5),rep(2,5),rep(3,5)),
c(0.0,0.0,0.5,0.4,0.5,
  0.0,0.3,0.0,0.3,0.5,
  0.0,0.3,0.5,0.7,1.0,
  0.0,0.4,0.2,0.8,1.0)
))

p = ggplot(d,aes(V1,-V2,fill=V3)) +
	geom_tile() +
	coord_fixed() +
	scale_fill_gradientn(colors=brewer.pal(9,'Blues')) +
	theme_void() +
	theme(legend.position='right')
ggsave(p,file='figures/final/regulatory_heatmap1.pdf',height=5,width=5,useDingbats=FALSE)

p = ggplot(d,aes(V1,-V2,fill=V3)) +
	geom_tile() +
	coord_fixed() +
	scale_fill_gradientn(colors=brewer.pal(9,'Reds')) +
	theme_void() +
	theme(legend.position='right')
ggsave(p,file='figures/final/regulatory_heatmap2.pdf',height=5,width=5,useDingbats=FALSE)





# Zoom in on RNA reads in relation to poly-A site
p5 = ggplot() +
	geom_path(data=subset(cell.rna.coverage.plot,V2 >= 3082000),aes(x=V2+50,y=cpm,color=cell_class,alpha=oli),size=0.5) +
	geom_vline(xintercept=3086373) + geom_vline(xintercept=3084965) + facet_wrap(~oli,ncol=1) +
	scale_color_manual(values=colors[levels(cell.rna.coverage.plot$cell_class)]) +
	scale_alpha_manual(values=c(1,0.5)) +
	# scale_color_viridis() +
	theme_classic() +
	theme(
		legend.position='none',
		# axis.text.x = element_blank(),
		axis.title = element_blank(),
		# axis.ticks.x = element_blank(),
		axis.line.x = element_blank(),
		# axis.text.y = element_blank(),
		# axis.ticks.y = element_blank(),
		strip.background=element_blank(),
		strip.text = element_blank()
	) + 
	ylab(expression('RNA coverage'))



# DE gene stats
de.genes = readRDS(paste0('rds/','rna','_marker_genes.rds'))

de.genes = within(de.genes,{
cell = factor(cell)
logfoldchanges = as.numeric(logfoldchanges)
wilcox_score = as.numeric(wilcox_score)
ttest_score = as.numeric(ttest_score)
logreg_score = as.numeric(logreg_score)
bp1 = as.integer(bp1)
bp2 = as.integer(bp2)
gene_strand = factor(gene_strand,levels=c('+','-'))
})

# mbp

regulatory.results$gene_dist = with(regulatory.results,
	ifelse(
		strand > 0,
		ifelse(	# case if + strand
			peak_end < tss, # if the peak ends before tss,
			tss - peak_end, # the distance from the gene is the distance from tss
			ifelse(
				peak_start > end_position, # If the peak starts after the gene end
				peak_start - end_position, # the peak is downstream of the whole gene
				0
			)
		),
		ifelse(	# case if - strand
			peak_start > tss, # if the peak starts after tss,
			peak_start - tss, # the distance from the gene is the distance from tss
			ifelse(
				peak_end < start_position, # if the peak ends before tss
				start_position - peak_start, # the peak is downstream and sign of distance is positive
				0
			)
		)
	)	
)

mbp = droplevels(subset(regulatory.results,cell=='all cells' & external_gene_name == 'MBP'))

mbp$isoform_dist = with(mbp,ifelse(
	peak_start > 3046976,
	peak_start - 3046976,
	ifelse(peak_end > 3046976,0,peak_end-3046976)
))

























# GLUE quality plots

in_metrics = read.delim('stats/glue/biccn_dx_consistency.txt',row.names=1)

dir.create('figures/glue',showWarnings=FALSE)

p = ggplot(in_metrics,aes(n_meta,consistency)) +
	geom_path(color='#0078b9') +
	geom_hline(yintercept=0.05,color='#961a1d',linetype=2) +
	theme_classic(base_size=24) +
	scale_y_continuous(
		limits=c(0,max(in_metrics$consistency)),
		breaks=seq(0,max(in_metrics$consistency),0.1)) +
	theme(panel.grid.major=element_line(size=0.1,color='#cccccc'))
# ggsave(p,file=file.path('figures/glue',paste0(prefix,'_dx_consistency.pdf')))

in_metrics_all = do.call(rbind,lapply(c(1:3,6:12,15),function(i) {
	x = read.delim(file.path('stats/glue/subintegration',paste0('subpeak_dx_consistency_class',i,'.txt')))[c('n_meta','consistency')]
	data.frame(cell_class=cell.classes[i],x)
}))

in_metrics_all = rbind(
	data.frame(cell_class = 'all cells',in_metrics),
	in_metrics_all
)

in_metrics_all$cell_class = factor(
	in_metrics_all$cell_class,
	levels=c('all cells',cell.levels)
)[,drop=TRUE]

in_metrics_all$cell_class = factor(in_metrics_all$cell_class,
		levels = names(sort(tapply(in_metrics_all$consistency, in_metrics_all$cell_class,mean),decreasing=TRUE))
	)

p = ggplot(in_metrics_all,aes(n_meta,consistency,color=cell_class)) +
	geom_path(size=0.5) +
	# geom_path(color='#0078b9') +
	# geom_hline(yintercept=0.05,color='#b45656',linetype=2) +
	geom_hline(yintercept=0.05,color='#000000',linetype=2) +
	facet_wrap(~cell_class) +
	theme_classic() +
	# scale_color_manual(values=c('#0078b9','#b45656')) +
	scale_color_manual(values=colors[levels(in_metrics_all$cell_class)]) +
	scale_y_continuous(
		limits=c(-0.1,max(in_metrics$consistency)),
		breaks=seq(0,max(in_metrics$consistency),0.1)) +
	theme(
		panel.grid.major=element_line(size=0.1,color='#cccccc'),
		legend.title=element_blank(),
		strip.background=element_blank(),
		legend.position='none'
	) +
	xlab('number of metacells')
ggsave(p,file='figures/final/integration_score.pdf',height=5,width=7,useDingbats=FALSE)





these.predictions = read.delim(file.path('stats/clusters',paste0('rnaext','-celltype-predictions.txt')),row.names=1)

query.prefix = 'rnaext'
validation.stats = read.delim(file.path('stats/glue_predictions',paste0(query.prefix,'-prediction-stats.txt')))

library(ggplot2)
library(egg)

glue.confidence.threshold = 0.95

p1 = ggplot(validation.stats,aes(threshold,log10(1-perc.match))) +
	geom_line() +
	geom_vline(xintercept=glue.confidence.threshold,color='red',size=0.5) +
	geom_hline(yintercept=log10(1-subset(validation.stats,threshold==0.95)$perc.match),color='#000000',size=0.5,linetype=3) +
	geom_text(
		data=data.frame(x=glue.confidence.threshold,y=log10(1-subset(validation.stats,threshold==0.95)$perc.match)),
		aes(x,y,label=paste0('error = ',round((1-subset(validation.stats,threshold==0.95)$perc.match) * 100,2),'%')),
		hjust=1,
		nudge_y=-0.05,
		nudge_x=-0.01
	) +
	theme_classic(base_size=16) +
	xlab('confidence score threshold') +
	ylab(expression('log'[10]~'error'))
ggsave(p1,file='figures/final/integration_prediction_accuracy.pdf',height=3,width=7,useDingbats=FALSE)

p2 = ggplot(validation.stats,aes(threshold,log10(n.cells))) +
	geom_line() +
	geom_vline(xintercept=glue.confidence.threshold,color='red',size=0.5) +
	geom_hline(yintercept=log10(subset(validation.stats,threshold==0.95)$n.cells),color='#000000',size=0.5,linetype=3) +
	geom_text(
		data=data.frame(x=glue.confidence.threshold,y=log10(subset(validation.stats,threshold==0.95)$n.cells)),
		# aes(x,y,label=eval(parse(text=paste0('expression(italic(N)~"="~"',format(subset(validation.stats,threshold==0.95)$n.cells,big.mark=','),'")')))),
		aes(x,y,label=paste0('N = ',format(subset(validation.stats,threshold==0.95)$n.cells,big.mark=','),'')),
		hjust=1,
		nudge_y=-0.05,
		nudge_x=-0.01
	) +
	theme_classic(base_size=16) +
	theme(axis.title.x=element_blank(),axis.text.x=element_blank()) +
#	xlab('confidence score threshold') +
	ylab(expression('log'[10]~'num. cells'))

pdf(file='/dev/null')
p = ggarrange(p2,p1,nrow=2)
dev.off()

pdf(file=file.path('figures/final',paste0('integration-prediction-stats-final.pdf')),useDingbats=FALSE,height=5)
p
dev.off()

validation.stats.all = do.call(rbind,lapply(c(1:3,6:12,15),function(i) {
	x = read.delim(file.path('stats/glue_predictions/subintegration',paste0('subpeak-rna-class',i,'-prediction-stats.txt')))
	data.frame(cell_class=cell.classes[i],x)
}))

these.predictions = read.delim(file.path('stats/glue_predictions/subintegration',paste0('subpeak-rna-class',,'-prediction-stats.txt')),row.names=1)


validation.stats.all = rbind(
	data.frame(cell_class='all cells',validation.stats),
	validation.stats.all
)
validation.stats.all$cell_class = factor(
	validation.stats.all$cell_class,
	levels=c('all cells',cell.levels)
)[,drop=TRUE]

p2 = ggplot(validation.stats.all,aes(threshold,log10(1-perc.match))) +
	geom_line() +
	# geom_vline(xintercept=glue.confidence.threshold,color='red',size=0.5) +
	# geom_hline(yintercept=log10(1-subset(validation.stats,threshold==0.95)$perc.match),color='#000000',size=0.5,linetype=3) +
	# geom_text(
	# 	data=data.frame(x=glue.confidence.threshold,y=log10(1-subset(validation.stats,threshold==0.95)$perc.match)),
	# 	aes(x,y,label=paste0('error = ',round((1-subset(validation.stats,threshold==0.95)$perc.match) * 100,2),'%')),
	# 	hjust=1,
	# 	nudge_y=-0.05,
	# 	nudge_x=-0.01
	# ) +
	facet_wrap(~cell_class) +
	theme_classic(base_size=16) +
	xlab('confidence score threshold') +
	ylab(expression('log'[10]~'error'))


p2 = ggplot(validation.stats,aes(threshold,log10(n.cells))) +
	geom_line() +
	geom_vline(xintercept=glue.confidence.threshold,color='red',size=0.5) +
	theme_classic(base_size=16) +
	theme(axis.title.x=element_blank(),axis.text.x=element_blank()) +
#	xlab('confidence score threshold') +
	ylab(expression('log'[10]~'num. cells'))
# 
# pdf(file='/dev/null')
# p = ggarrange(p2,p1,nrow=2)
# dev.off()
# 
# pdf(file=file.path('figures',paste0(query.prefix,'-glue-prediction-stats.pdf')),useDingbats=FALSE,height=4)
# p
# dev.off()



# Marker genes plot

de.genes = readRDS(paste0('rds/','rna','_marker_genes.rds'))

de.genes = within(de.genes,{
	cell = factor(cell)
	logfoldchanges = as.numeric(logfoldchanges)
	wilcox_score = as.numeric(wilcox_score)
	ttest_score = as.numeric(ttest_score)
	logreg_score = as.numeric(logreg_score)
	bp1 = as.integer(bp1)
	bp2 = as.integer(bp2)
	gene_strand = factor(gene_strand,levels=c('+','-'))
})

levels(de.genes$cell) = gsub('serotinergic','serotonergic',gsub('neuron$','neurons',levels(de.genes$cell)))

de.genes$cell = factor(de.genes$cell,levels=cell.levels)[,drop=TRUE]

de.genes$cell = factor(de.genes$cell,levels=c("excitatory neurons", "cerebellar neurons", "inhibitory neurons", 
"basket cells", "medium spiny neurons", "dopaminergic neurons", 
"serotonergic neurons", 
"astrocytes", "oligodendrocytes", "oligodendrocyte precursor cells", 
"ependymal cells", "microglia", "vascular cells","AHSG neurons", "F5 neurons","KIR3DL12 neurons", "KIR3DL12 microglia"
))

keep.genes = do.call(rbind,lapply(split(de.genes,de.genes$cell),function(x) {
	x = subset(x,pts > 0.10 & ttest_pvals_adj  < 0.05 & !grepl('^ENSMMU',gene_short_name))
	if (unique(x$cell) %in% c('F5 neurons','AHSG neurons','KIR3DL12 neurons','KIR3DL12 microglia')) {
		head(x[order(x$logreg_score,decreasing=TRUE),],1)[,c('cell','gene_short_name')]
	} else {
		head(x[order(x$logreg_score,decreasing=TRUE),],3)[,c('cell','gene_short_name')]
	}
}))

cat(unlist(lapply(split(de.genes,de.genes$cell),function(x) {
	x = subset(x,pts > 0.10 & ttest_pvals_adj  < 0.05)
	paste(head(x[order(x$logreg_score,decreasing=TRUE),],10)$gene_short_name,collapse=',')
})),sep='\n')




de.keep = subset(de.genes,gene_short_name %in% keep.genes$gene_short_name)
de.keep$gene_short_name = factor(de.keep$gene_short_name,levels=unique(keep.genes$gene_short_name))


p = ggplot(de.keep,aes(gene_short_name,cell,color=logfoldchanges,size=pts*100)) +
	geom_point() +
	scale_y_discrete(limits=rev) +
	# scale_color_viridis(limits=c(-10,10),oob=scales::squish) +
	scale_color_gradient2(
		high='#2066ac',
		mid='#ffffff',
		low='#b2182b',
		midpoint=0,
		breaks=seq(-10,10,5),
		# limits=c(-max(abs(tf.keep$log2OR)),max(abs(tf.keep$log2OR))),
		limits=c(-10,10),
		oob=scales::squish,
		name=expression(log[2]~'FC')
	) +
	scale_size_continuous(name='% expressed',range=c(1,3)) +
	coord_fixed() +
	theme_classic(base_size=8) +
	theme(
		axis.text.x = element_text(vjust=1,hjust=0,angle=-45),
		axis.title = element_blank(),
		legend.text=element_text(margin=margin(t=0))
	)
ggsave(p,file='figures/final/marker_genes_plot.pdf',height=7,width=7,useDingbats=FALSE)


# Download MBP isoforms
ens = biomaRt::useEnsembl(biomart='ENSEMBL_MART_ENSEMBL',dataset=paste0('hsapiens','_gene_ensembl'))

ens.orthologs = biomaRt::getBM(
				attributes=c('ensembl_gene_id','mmulatta_homolog_ensembl_gene'),
				mart = ens)
names(ens.orthologs) = c('hsa','mmu')

ens.orthologs = unique(subset(ens.orthologs,nchar(hsa) == 15 & nchar(mmu) == 18))
write.table(ens.orthologs,file='stats/hsa_mmu_orthologs.txt',sep='\t',row.names=FALSE,quote=FALSE)




# Summarize regulatory

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

prefix = 'biccn'
atac.prefix = 'atac'

suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrastr))

cell.classes = scan(what='',sep='\n',file='stats/clusters/rna-final-cellclasses-levels.txt',quiet=TRUE)
this.cell.class = 'all cells'

metacell.lr.results = read.table(
	file.path('stats/metacell_lr',paste0('gene_peak_all.tsv')),
	col.names = c('meta.beta','meta.serr','meta.sbet','meta.pval','peak','ensembl_gene_id'),
	sep='\t'
)

metacell.lr.null = read.table(
	file.path('stats/metacell_lr',paste0('gene_peak_null_all.tsv')),
	col.names = c('null.beta','null.serr','null.sbet','null.pval','peak','ensembl_gene_id'),
	sep='\t'
)

metacell.lr.results = merge(metacell.lr.results,metacell.lr.null,all.x=TRUE,all.y=TRUE,by=c('ensembl_gene_id','peak'))

glue.regulatory.results = fread(
	file.path('stats/glue',paste0(prefix,'_regulatory_regulatory.txt.gz'))
#	col.names = c('ensembl_gene_id','peak','glue.qval','dist','glue.pval','glue.score','glue.weight','glue.type')
)[,c('source','target','qval','dist','pval','score','weight','type')]
names(glue.regulatory.results) = c('ensembl_gene_id','peak','glue.qval','dist','glue.pval','glue.score','glue.weight','glue.type')

results.combined = merge(metacell.lr.results,glue.regulatory.results,by=c('ensembl_gene_id','peak'),all.y=TRUE,all.x=TRUE)

results.combined$meta.qval = p.adjust(results.combined$meta.pval,'fdr')

if (any(!complete.cases(results.combined))) {
	warning(sum(!complete.cases(results.combined)),' records are incomplete and will be removed.')
	warning(sum(is.na(results.combined$meta.pval)),' metacell p-values are NA.')
	warning(sum(is.na(results.combined$null.pval)),' permuted metacell p-values are NA.')
	results.combined = results.combined[complete.cases(results.combined),]
}

dir.create('figures/regulatory',showWarnings=FALSE)

results.sig = subset(results.combined,glue.qval < 0.05 & meta.qval < 0.05)

results.sig$gene_dir = with(results.sig,ifelse(meta.beta>0,'+','-'))

gene.counts = as.data.frame.matrix(as.matrix(table(results.sig$ensembl_gene_id,results.sig$gene_dir)))
# gene.counts$direction = ifelse(gene.counts$`-` > 0 & gene.counts$`+` > 0,'±',ifelse(gene.counts$`+` > 0,'+','-'))
gene.counts$direction = ifelse(gene.counts$`-` > 0 & gene.counts$`+` > 0,'±','.')
gene.counts$ensembl_gene_id = rownames(gene.counts)
# gene.counts$direction = factor(gene.counts$direction,levels=c('-','+','±'))
gene.counts$direction = factor(gene.counts$direction,levels=c('.','±'))

gene.counts.long = subset(reshape2::melt(gene.counts,id.vars=c('ensembl_gene_id','direction'),measure.vars=c('-','+')),value > 0)

gene.counts.long$variable_text = factor(
	gene.counts.long$variable,
	levels=levels(gene.counts.long$variable),
	labels=paste0(levels(gene.counts.long$variable),' regulators (',table(gene.counts.long$variable),' genes, ',tapply(gene.counts.long$value,gene.counts.long$variable,sum),' peaks)')
)

# p = ggplot(gene.counts.long,aes(value,fill=direction)) +
# 	geom_histogram() +
# 	facet_wrap(~variable_text,nrow=2,scales='free_y') +
# 	scale_fill_manual(name='Gene direction',values=c('#91bfdb','#fc8d59','#999999')) +
# 	theme_classic(base_size=16) +
# 	xlab('Number of elements') +
# 	ylab('Count') +
# 	ggtitle(this.cell.class)
# suppressMessages(ggsave(p,file=file.path('figures/regulatory',paste0(prefix,'_',atac.prefix,'_all_genes-histogram.pdf')),useDingbats=FALSE))
# 
# p = ggplot(gene.counts.long,aes(value,fill=direction)) +
# 	geom_histogram() +
# 	facet_wrap(~variable_text,nrow=2,scales='fixed') +
# 	scale_fill_manual(name='Gene direction',values=c('#91bfdb','#fc8d59','#999999')) +
# 	theme_classic(base_size=16) +
# 	xlab('Number of elements') +
# 	ylab('Count') +
# 	ggtitle(this.cell.class)
# suppressMessages(ggsave(p,file=file.path('figures/regulatory',paste0(prefix,'_',atac.prefix,'_all_genes-histogram-fixed.pdf')),useDingbats=FALSE))

gene.counts.long$cell_class = 'all cells'

gene.counts.long2 = do.call(rbind,lapply(c(1:3,6:12,15),function(i) {
	
	metacell.lr.results = read.table(
		file.path('stats/metacell_lr',paste0('gene_peak_class',i,'.tsv')),
		col.names = c('meta.beta','meta.serr','meta.sbet','meta.pval','peak','ensembl_gene_id'),
		sep='\t'
	)
	metacell.lr.null = read.table(
		file.path('stats/metacell_lr',paste0('gene_peak_null_class',i,'.tsv')),
		col.names = c('null.beta','null.serr','null.sbet','null.pval','peak','ensembl_gene_id'),
		sep='\t'
	)
	metacell.lr.results = merge(metacell.lr.results,metacell.lr.null,all.x=TRUE,all.y=TRUE,by=c('ensembl_gene_id','peak'))

	glue.regulatory.results = fread(
		file.path('stats/glue/subintegration',paste0('subpeak','_class',i,'_regulatory.txt.gz'))
	#	col.names = c('ensembl_gene_id','peak','glue.qval','dist','glue.pval','glue.score','glue.weight','glue.type')
	)[,c('source','target','qval','dist','pval','score','weight','type')]
	names(glue.regulatory.results) = c('ensembl_gene_id','peak','glue.qval','dist','glue.pval','glue.score','glue.weight','glue.type')

	results.combined = merge(metacell.lr.results,glue.regulatory.results,by=c('ensembl_gene_id','peak'),all.y=TRUE,all.x=TRUE)

	results.combined$meta.qval = p.adjust(results.combined$meta.pval,'fdr')

	if (any(!complete.cases(results.combined))) {
		warning(sum(!complete.cases(results.combined)),' records are incomplete and will be removed.')
		warning(sum(is.na(results.combined$meta.pval)),' metacell p-values are NA.')
		warning(sum(is.na(results.combined$null.pval)),' permuted metacell p-values are NA.')
		results.combined = results.combined[complete.cases(results.combined),]
	}

	results.sig = subset(results.combined,glue.qval < 0.05 & meta.qval < 0.05)

	results.sig$gene_dir = with(results.sig,ifelse(meta.beta>0,'+','-'))

	gene.counts = as.data.frame.matrix(as.matrix(table(results.sig$ensembl_gene_id,results.sig$gene_dir)))
	# gene.counts$direction = ifelse(gene.counts$`-` > 0 & gene.counts$`+` > 0,'±',ifelse(gene.counts$`+` > 0,'+','-'))
	gene.counts$direction = ifelse(gene.counts$`-` > 0 & gene.counts$`+` > 0,'±','.')
	gene.counts$ensembl_gene_id = rownames(gene.counts)
	# gene.counts$direction = factor(gene.counts$direction,levels=c('-','+','±'))
	gene.counts$direction = factor(gene.counts$direction,levels=c('.','±'))

	gene.counts.long = subset(reshape2::melt(gene.counts,id.vars=c('ensembl_gene_id','direction'),measure.vars=c('-','+')),value > 0)

	gene.counts.long$variable_text = factor(
		gene.counts.long$variable,
		levels=levels(gene.counts.long$variable),
		labels=paste0(levels(gene.counts.long$variable),' regulators (',table(gene.counts.long$variable),' genes, ',tapply(gene.counts.long$value,gene.counts.long$variable,sum),' peaks)')
	)
	gene.counts.long$cell_class = cell.classes[i]
	gene.counts.long
}))

gene.counts.long.combined = rbind(gene.counts.long,gene.counts.long2)

gene.counts.long.combined$cell_class = factor(gene.counts.long.combined$cell_class,levels=c('all cells',cell.levels))[,drop=TRUE]




p = ggplot(droplevels(subset(gene.counts.long.combined,cell_class != 'all cells')),aes(value,fill=direction)) +
	geom_histogram() +
	facet_grid(rows=vars(cell_class),cols=vars(variable),scales='free_y') +
	scale_fill_manual(name='Gene direction',values=c('#999999','#e31a1c')) +
	theme_classic(base_size=12) +
	theme(
		strip.text.x = element_text(size=8),
		strip.text.y = element_text(angle=0,hjust=0,size=8),
		strip.background=element_blank(),
		legend.title=element_blank(),
		legend.position='none'
	) +
	xlab('number of elements') +
	ylab('count')
suppressMessages(ggsave(p,file=file.path('figures/final','regulatory_histogram.pdf'),useDingbats=FALSE,width=7,height=10))




da.peaks = readRDS(paste0('rds/','atac_all_marker_peaks_combined.rds'))


# atacsub_class1_marker_peaks_combined.rds




purkinje.de = readRDS(paste0('rds/','rna','_marker_genes.rds'))

purkinje.de =  within(purkinje.de,{
	cell = factor(cell)
	logfoldchanges = as.numeric(logfoldchanges)
	wilcox_score = as.numeric(wilcox_score)
	ttest_score = as.numeric(ttest_score)
	logreg_score = as.numeric(logreg_score)
})


purkinje.promoter = readRDS('rds/atac_class6_marker_genes_geneactivity.rds')

purkinje.promoter =  within(purkinje.promoter,{
	cell = factor(cell)
	logfoldchanges = as.numeric(logfoldchanges)
	wilcox_score = as.numeric(wilcox_score)
	ttest_score = as.numeric(ttest_score)
	logreg_score = as.numeric(logreg_score)
	Start = as.integer(Start)
	End = as.integer(End)
})

purkinje = subset(purkinje.promoter,cell == 'cerebellar neurons 16')

purkinje = purkinje[order(purkinje$logreg_score,decreasing=TRUE),]

top.markers = c('GRID2','PRKG1','AUTS2','ARHGAP26','ITPR1')

purkinje$in_set = purkinje$gene_short_name %in% top.markers
p = ggplot(purkinje,aes(logreg_score)) + geom_histogram() + facet_wrap(~in_set)


purkinje1 = purkinje.de[purkinje.de$cell %in% 'cerebellar neurons 16',c('cell','gene_short_name','logfoldchanges','pts','logreg_score')]
purkinje2 = purkinje.promoter[purkinje.promoter$cell %in% 'cerebellar neurons 16',c('cell','gene_short_name','logfoldchanges','pts','logreg_score')]

names(purkinje1)[3:5] = paste0('de_',names(purkinje1)[3:5])
names(purkinje2)[3:5] = paste0('pa_',names(purkinje2)[3:5])

purkinje = merge(purkinje1,purkinje2,by=c('cell','gene_short_name'))


p = ggplot(purkinje,aes(de_logreg_score,pa_logreg_score)) + geom_point(size=0.2) + geom_smooth(method=lm) +
	theme_classic()
ggsave(p,file='de_pa_purkinje.pdf',useDingbats=F)


# Subtype DE genes

de.genes = readRDS('rds/rna_class6_marker_genes.rds')
de.genes =  within(de.genes,{
	cell = factor(cell)
	logfoldchanges = as.numeric(logfoldchanges)
	wilcox_score = as.numeric(wilcox_score)
	ttest_score = as.numeric(ttest_score)
	logreg_score = as.numeric(logreg_score)
})

de.genes = de.genes[order(de.genes$logreg_score,decreasing=TRUE),]
purkinje.de = subset(de.genes,cell == 'cerebellar neurons 16')


# Subtype DA peaks

# purkinje.monalisa = readRDS(paste0('rds/monalisa/class',6,'_cluster',16,'_','merged_unique_peaks','_TFenrich_df.rds'))

purkinje.monalisa = readRDS(paste0('rds/monalisa/class',6,'_cluster',16,'_','top_1p_peaks','_TFenrich_df.rds'))


prefix='biccn'
atac.prefix='atac'
inc = read.delim(file.path('stats/regulatory',paste0(prefix,'_',atac.prefix,'_all_peaks_inc.txt')),header=FALSE)
dec = read.delim(file.path('stats/regulatory',paste0(prefix,'_',atac.prefix,'_all_peaks_dec.txt')),header=FALSE)
names(inc) = names(dec) = c('peak','ensembl_gene_id','dist','gene_dir','glue.score','glue.pval','glue.qval','meta.beta','meta.pval','meta.qval')

# da.peaks = readRDS(paste0('rds/','atac_all_marker_peaks_combined.rds'))


da.keep = subset(da.peaks,logreg_score > 0 & logfoldchanges > 0 & ttest_pvals_adj < 0.05)
da.keep = droplevels(da.keep[!da.keep$cluster %in% 'radial glial cells',])


gene.peaks = unique(reg[c('peak')])

 
da.keep = subset(da.keep,peak %in% gene.peaks$peak,select=c('cluster','peak')) 

da.merged = do.call(rbind,lapply(split(da.keep,da.keep$peak),function(x) data.frame(peak = unique(x$peak),cluster=paste(x$cluster,collapse='|'))))

reg = rbind(inc,dec)


reg = merge(reg,da.merged,by='peak',all.x=TRUE)
ens.genes = readRDS(file.path('checkpoints',paste0(ens.species,'_genes.rds')))


reg.display = subset(reg,glue.qval < 0.005 & meta.qval < 0.005)

ens.genes = subset(ens.genes,select=c('ensembl_gene_id','external_gene_name'))

reg.display = merge(reg.display,ens.genes,by='ensembl_gene_id',all.x=TRUE)

reg.display = data.frame(as.data.frame.matrix(do.call(rbind,strsplit(reg.display$peak,'_'))),reg.display)

write.table(reg.display,file='tables/ccre_genes.txt',sep='\t',quote=FALSE)




# Now do again at subtype level

prefix='subpeak'
atac.prefix='atacsub'

ens.species='mmulatta'
ens.genes = readRDS(file.path('checkpoints',paste0(ens.species,'_genes.rds')))
ens.genes = subset(ens.genes,select=c('ensembl_gene_id','external_gene_name'))

out = vector(7,mode='list')
names(out) = as.character(c(1,3,6,8:9,11:12))

#reg = do.call(rbind,lapply(c(1,3,6:9,12),function(i) {
for (i in c(1,3,6,8:9,11:12)) {
print(i)
inc = read.delim(file.path('stats/regulatory',paste0(prefix,'_',atac.prefix,'_class',i,'_peaks_inc.txt')),header=FALSE)
dec = read.delim(file.path('stats/regulatory',paste0(prefix,'_',atac.prefix,'_class',i,'_peaks_dec.txt')),header=FALSE)
names(inc) = names(dec) = c('peak','ensembl_gene_id','dist','gene_dir','glue.score','glue.pval','glue.qval','meta.beta','meta.pval','meta.qval')
reg = subset(rbind(inc,dec),glue.qval < 0.005 & meta.qval < 0.005)

da.peaks = readRDS(file.path('rds',paste0('atacsub','_class',i,'_marker_peaks_combined.rds')))
da.keep = subset(da.peaks,logreg_score > 0 & logfoldchanges > 0 & ttest_pvals_adj < 0.05)

gene.peaks = unique(reg[c('peak')])

da.keep = subset(da.keep,peak %in% gene.peaks$peak,select=c('cluster','peak')) 

da.merged = do.call(rbind,lapply(split(da.keep,da.keep$peak),function(x) data.frame(peak = unique(x$peak),cluster=paste(x$cluster,collapse='|'))))

reg = merge(reg,da.merged,by='peak',all.x=TRUE)

out[[as.character(i)]] = reg
}
# reg
# }))
# 

for (i in c(1,3,6,8:9,11:12)) {
out[[as.character(i)]] = data.frame(cell_class=gsub('excitatory','glutamatergic',gsub('inhibitory','GABAergic',cell.classes))[i],out[[as.character(i)]])
}

reg = do.call(rbind,out)

reg$cluster = gsub('inhibitory','GABAergic',gsub('excitatory','glutamatergic',reg$cluster))


reg.display = subset(reg,glue.qval < 0.001 & meta.qval < 0.001)
reg.display$cluster = gsub('[A-z ]+','',reg.display$cluster)


reg.display = merge(reg.display,ens.genes,by='ensembl_gene_id',all.x=TRUE)

reg.display = data.frame(as.data.frame.matrix(do.call(rbind,strsplit(reg.display$peak,'_'))),reg.display)


cell.levels = c('glutamatergic neurons','cerebellar neurons','GABAergic neurons','basket cells','medium spiny neurons','dopaminergic neurons','serotonergic neurons','AHSG neurons','F5 neurons','KIR3DL12 neurons','astrocytes','oligodendrocytes','oligodendrocyte precursor cells','radial glial cells','ependymal cells','microglia','KIR3DL12 microglia','vascular cells','mesenchymal stem cells')


reg.display$cell_class = factor(reg.display$cell_class,levels=cell.levels)[,drop=TRUE]


reg.display = reg.display[with(reg.display,order(cell_class,V1,V2)),]

reg.display$gene_dir = factor(reg.display$gene_dir,levels=c('-','+'),labels=c('negative','positive'))

write.table(reg.display,file='tables/final/regulatory_subtype.txt',sep='\t',quote=FALSE,row.names=FALSE)




# Put together metadata table

library(data.table)
rna = subset(read.delim('data/rna_metadata_final.txt',row.names=1),select=c('id','region','hemisphere','animal_id','social_group','sex'))
atac = subset(read.delim('data/atac_metadata.txt',row.names=1),select=c('region','hemisphere','animal_id','social_group','sex'))
atac$id = rownames(atac)

rna2 = fread('data/rna_run-data.txt.gz')[,c('run_id','sample')]
names(rna2)[2] = 'id'
atac2 = fread('data/atac_run-data.txt.gz')[,c('run_id','sample')]
names(atac2)[2] = 'id'


rna = merge(rna,rna2,by='id')
atac = merge(atac,atac2,by='id')

rna = subset(rna,id!='NSM349.b')

foo = merge(rna,atac,by=c('animal_id','region','hemisphere'),all.x=T,all.y=T)
subset(foo,region=='STS' & animal_id == '2C0')

foo$region = factor(foo$region,levels=region.levels)
foo$animal_sort = factor(foo$animal_id,levels=c('3I4','4I3','6J2','2C0','3R7'),labels=c('F1','F1','M1','F2','M2'))

foo$sex = with(foo,ifelse(is.na(sex.x),sex.y,sex.x))
foo$social_group = with(foo,ifelse(is.na(social_group.x),social_group.y,social_group.x))

# foo$hemisphere = factor(foo$hemisphere,levels=c('R','L'))

foo2 = foo[with(foo,order(animal_sort,region,hemisphere)),c('animal_id','region','hemisphere','sex','social_group','id.x','run_id.x','id.y','run_id.y')]


write.table(foo2,file='tables/sample_metadata.txt',sep='\t',row.names=FALSE,quote=FALSE)

# 

rna = droplevels(subset(this.umap,modality=='RNA' & keep))

levels(rna$cell_subtype) = gsub('^AHSG','APOA2',gsub('^inhibitory','GABAergic',gsub('^excitatory','glutamatergic',levels(rna$cell_subtype))))
levels(rna$cell_class) = gsub('^AHSG','APOA2',gsub('^inhibitory','GABAergic',gsub('^excitatory','glutamatergic',levels(rna$cell_class))))

rna$cell_subtype = factor(rna$cell_subtype,levels=cell.subtype.levels)

rna$id2 = gsub('\\.[ab]$','',rna$id)
x = table(rna$id2)




atac = droplevels(subset(this.umap,modality=='ATAC' & keep))

levels(atac$cell_subtype) = gsub('^AHSG','APOA2',gsub('^inhibitory','GABAergic',gsub('^excitatory','glutamatergic',levels(atac$cell_subtype))))
levels(atac$cell_class) = gsub('^AHSG','APOA2',gsub('^inhibitory','GABAergic',gsub('^excitatory','glutamatergic',levels(atac$cell_class))))

atac$cell_subtype = factor(atac$cell_subtype,levels=cell.subtype.levels)

atac$id2 = gsub('\\.[ab]$','',atac$id)
x = table(atac$id2)


cell.type.column = 'cell_subtype'
region.column = 'region'

cell.metadata = rna

# Compute Jensen-Shannon divergence
specificity.out = data.frame(cell_type = levels(cell.metadata[[cell.type.column]]))

# Make a table of cell type by region
m = as.matrix(table(cell.metadata[[cell.type.column]],cell.metadata[[region.column]]))

# For each cell type, compare its frequency distribution to the full distribution across all regions
specificity.out$specificity_score = unlist(lapply(1:nrow(m),function(i) {
	out = suppressMessages(JSD(rbind(m[i,],colSums(m)),est.prob='empirical'))
	names(out) = rownames(m)[i]
	out
}))




# 

cell.levels = c('excitatory neurons','cerebellar neurons','inhibitory neurons','basket cells','medium spiny neurons','dopaminergic neurons','serotonergic neurons','AHSG neurons','F5 neurons','KIR3DL12 neurons','astrocytes','oligodendrocytes','oligodendrocyte precursor cells','radial glial cells','ependymal cells','microglia','KIR3DL12 microglia','vascular cells','mesenchymal stem cells')
region.levels = c('dmPFC','vmPFC','dlPFC','vlPFC','ACC','M1','EC','PC','A1','STS','MT','IT','S1','IPP','SPP','V1','CC','CN','NAc','AMY','HIP','mdTN','vlTN','LGN','CV','lCb','MB','MdO','MdC','Pons')


rna$cell_class = factor(rna$cell_class,levels=cell.levels,
labels=c('glutamatergic neurons','cerebellar neurons','GABAergic neurons','basket cells','medium spiny neurons','dopaminergic neurons','serotonergic neurons','APOA2 neurons','F5 neurons','KIR3DL12 neurons','astrocytes','oligodendrocytes','oligodendrocyte precursor cells','radial glial cells','ependymal cells','microglia','KIR3DL12 microglia','vascular cells','mesenchymal stem cells')
)[,drop=TRUE]


with(rna,table(cell_class,region)) %>% apply(2,function(x) x/sum(x)) %>%  as.matrix %>% as.data.frame.matrix %>% write.table(file='tables/cell_composition.txt',sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)


with(rna,table(cell_subtype,region)) %>% apply(2,function(x) x/sum(x)) %>%  as.matrix %>% as.data.frame.matrix %>% write.table(file='tables/cell_subtype_composition.txt',sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)





de.genes.all = do.call(rbind,lapply(file.path('rds',list.files('rds',pattern='rna_class[0-9]+_marker_genes.rds')),readRDS))

de.genes.all$cell = factor(de.genes.all$cell)

levels(de.genes.all$cell) = gsub('neuron$','neurons',gsub('serotinergic','serotonergic',gsub('^AHSG','APOA2',gsub('^inhibitory','GABAergic',gsub('^excitatory','glutamatergic',levels(de.genes.all$cell))))))


de.genes.all = droplevels(subset(de.genes.all,!cell %in% c('radial glial cells','mesenchymal stem cells')))

de.genes.all$cell = factor(de.genes.all$cell,levels=cell.subtype.levels)


unlist(lapply(split(de.genes.all,de.genes.all$cell),function(x) {
x = x[x$pts > 0.1 & x$ttest_pvals_adj < 0.05,]
paste(head(x[order(x$logreg_score,decreasing=TRUE),],5)$gene_short_name,collapse=',')
})) %>% cat(sep='\n')



tf.subtype.results = do.call(rbind,lapply(file.path('rds/monalisa',list.files('rds/monalisa',pattern='class[0-9]+_cluster[0-9]+_top_1p_peaks_TFenrich_df.rds')),function(x) data.frame(x,readRDS(x))))


tf.subtype.results$class = as.integer(gsub('rds/monalisa/class([0-9]+)_cluster([0-9]+)_top_1p_peaks_TFenrich_df.rds','\\1',tf.subtype.results$x))
tf.subtype.results$cluster = as.integer(gsub('rds/monalisa/class([0-9]+)_cluster([0-9]+)_top_1p_peaks_TFenrich_df.rds','\\2',tf.subtype.results$x))

tf.subtype.results$cell_class = (gsub('^excitatory','glutamatergic',cell.classes) %>% { gsub('^inhibitory','GABAergic',.) } %>% { gsub('serotinergic','serotonergic',.) })[tf.subtype.results$class]

tf.subtype.results$cell_subtype = paste(tf.subtype.results$cell_class,tf.subtype.results$cluster)

tf.subtype.results$cell_subtype = factor(tf.subtype.results$cell_subtype,levels=cell.subtype.levels)[,drop=TRUE]



subtype.levels = c("excitatory neurons 1", "excitatory neurons 2", "excitatory neurons 3",
"excitatory neurons 4", "excitatory neurons 5", "excitatory neurons 6",
"excitatory neurons 7", "excitatory neurons 8", "excitatory neurons 9",
"excitatory neurons 10", "excitatory neurons 11", "excitatory neurons 12",
"excitatory neurons 13", "excitatory neurons 14", "excitatory neurons 15",
"excitatory neurons 16", "excitatory neurons 17", "excitatory neurons 18",
"excitatory neurons 19", "excitatory neurons 20", "excitatory neurons 21",
"excitatory neurons 22", "excitatory neurons 23", "excitatory neurons 24",
"excitatory neurons 25", "excitatory neurons 26", "excitatory neurons 27",
"excitatory neurons 28", "excitatory neurons 29", "excitatory neurons 30",
"excitatory neurons 31", "excitatory neurons 32", "excitatory neurons 33",
"excitatory neurons 34", "excitatory neurons 35", "excitatory neurons 36",
"excitatory neurons 37", "excitatory neurons 38", "excitatory neurons 39",
"medium spiny neurons 1", "medium spiny neurons 2", "inhibitory neurons 1",
"inhibitory neurons 2", "inhibitory neurons 3", "inhibitory neurons 4",
"inhibitory neurons 5", "inhibitory neurons 6", "inhibitory neurons 7",
"inhibitory neurons 8", "inhibitory neurons 9", "inhibitory neurons 10",
"inhibitory neurons 11", "inhibitory neurons 12", "inhibitory neurons 13",
"inhibitory neurons 14", "inhibitory neurons 15", "inhibitory neurons 16",
"inhibitory neurons 17", "inhibitory neurons 18", "inhibitory neurons 19",
"inhibitory neurons 20", "dopaminergic neurons 1", "dopaminergic neurons 2",
"serotinergic neurons 1", "serotinergic neurons 2", "cerebellar neurons 1",
"cerebellar neurons 2", "cerebellar neurons 3", "cerebellar neurons 4",
"cerebellar neurons 5", "cerebellar neurons 6", "cerebellar neurons 7",
"cerebellar neurons 8", "cerebellar neurons 9", "cerebellar neurons 10",
"cerebellar neurons 11", "cerebellar neurons 12", "cerebellar neurons 13",
"cerebellar neurons 14", "cerebellar neurons 15", "cerebellar neurons 16",
"astrocytes 1", "astrocytes 2", "astrocytes 3", "astrocytes 4",
"astrocytes 5", "astrocytes 6", "oligodendrocytes 1", "oligodendrocytes 2",
"oligodendrocytes 3", "oligodendrocytes 4", "oligodendrocytes 5",
"oligodendrocytes 6", "oligodendrocytes 7", "oligodendrocytes 8",
"vascular cells 1", "vascular cells 2", "vascular cells 3", "vascular cells 4",
"vascular cells 5", "vascular cells 6", "microglia 1", "microglia 2"
)

subtype.levels.corrected = subtype.levels %>% { gsub('^excitatory','glutamatergic',.) } %>% { gsub('^inhibitory','GABAergic',.) } %>% { gsub('serotinergic','serotonergic',.) }


disease.ft.results = readRDS(paste0('rds/','rna','-recluster-marker-gene-enrichment-disease-ft.rds'))
disease.ft.results$cell = factor(disease.ft.results$cell,levels=subtype.levels,labels=subtype.levels.corrected)

disease.ft.display = subset(disease.ft.results,padj < 0.05,select=c('cell','pathway_id','pathway_name','estimate','padj'))
write.table(disease.ft.display,file='tables/final/disease_marker_gene_enrich.txt',row.names=FALSE,sep='\t',quote=FALSE)






# LD scores
ld.stats = read.delim('tables/final/ldsc_summary_stats_allpeaks.txt')

ld.stats$Cell_class = factor(ld.stats$Cell_class,levels=gsub('GABAergic','inhibitory',gsub('glutamatergic','excitatory',cell.levels)),labels=cell.levels)[,drop=TRUE]

ld.stats$Phenotype=factor(ld.stats$Phenotype,levels=scan(what='',sep='\n'))

ld.stats = ld.stats[with(ld.stats,order(Cell_class,Enrichment_p)),]

ld.display = subset(ld.stats,Enrichment > 0 & Enrichment_padj < 0.05)

write.table(ld.display,file='tables/final/ld_final.txt',sep='\t',quote=FALSE,row.names=FALSE)



# 
x = readRDS('rds/rna_class3_marker_genes.rds')
gaba = subset(x,gene_short_name %in% c('SST','PVALB','PAX6','VIP','LAMP5','ADARB2'))
