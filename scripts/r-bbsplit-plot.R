#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

library(Matrix)
library(tidyverse)
library(data.table)
library(viridis)
library(RColorBrewer)
library(ggrastr)

arguments = commandArgs(trailing=TRUE)
# arguments = c('rna','NSM084','macaque','human','mouse')

prefix = arguments[1]
this_id = arguments[2]
genome1 = arguments[3]
genome2 = arguments[4]
genome3 = arguments[5]

meta = fread(file.path(paste0('stats/contamination/',genome1,'_',genome2,'_',genome3),paste0(this_id,'_assigned.txt.gz')))

colors = as.matrix(read.table('data/colors.txt',row.names=1,sep='\t',comment.char='',header=FALSE))[,1]
cell.levels = c('excitatory neurons','cerebellar neurons','inhibitory neurons','basket cells','medium spiny neurons','dopaminergic neurons','serotonergic neurons','AHSG neurons','F5 neurons','KIR3DL12 neurons','astrocytes','oligodendrocytes','oligodendrocyte precursor cells','radial glial cells','ependymal cells','microglia','KIR3DL12 microglia','vascular cells','mesenchymal stem cells')

colors = c(colors,c('all cells' = '#000000'))

cell.classes = c('excitatory neurons','medium spiny neurons','inhibitory neurons','dopaminergic neurons','serotinergic neurons','cerebellar neurons','basket cells','astrocytes','oligodendrocytes','oligodendrocyte precursor cells','vascular cells','microglia','radial glial cells','mesenchymal stem cells','ependymal cells','AHSG neuron','F5 neuron','KIR3DL12 neuron','KIR3DL12 microglia')
meta$cell_class = factor(meta$cell_class,levels=c(cell.classes,''),labels=c(gsub('neuron$','neurons',gsub('serotinergic','serotonergic',cell.classes)),'unclassified')) 

meta$cell_class = factor(meta$cell_class,levels=c(cell.levels,'unclassified'))

colors = c(colors,c('unclassified' = '#cccccc'))

names(colors)[names(colors) == 'serotinergic neurons'] = 'serotonergic neurons'

names(colors) = gsub('neuron$','neurons',gsub('serotinergic','serotonergic',names(colors)))

dir.create('figures/contamination',showWarnings=FALSE)

p = ggplot(subset(meta,n_assigned >= 10 & cell_class != 'unclassified'),aes_string(paste0('n_',genome1),paste0('n_',genome2),color='cell_class')) +
	geom_point_rast(size=0.5,alpha=0.5) +
	coord_cartesian() +
	scale_color_manual(values=colors[levels(meta$cell_class)]) +
	guides(color=guide_legend(override.aes=list(size=2,alpha=1))) +
	theme_classic() +
	theme(legend.position = 'right',legend.title=element_blank()) +
	xlab(paste0(genome1,' count')) +
	ylab(paste0(genome2,' count'))
ggsave(p,file=file.path('figures/contamination',paste0(this_id,'-',genome1,'_',genome2,'_barnyard.pdf')),useDingbats=FALSE)

p = ggplot(subset(meta,n_assigned >= 10 & cell_class != 'unclassified'),aes_string(paste0('n_',genome1),paste0('n_',genome3),color='cell_class')) +
	geom_point_rast(size=0.5,alpha=0.5) +
	coord_cartesian() +
	scale_color_manual(values=colors[levels(meta$cell_class)]) +
	guides(color=guide_legend(override.aes=list(size=2,alpha=1))) +
	theme_classic() +
	theme(legend.position = 'right',legend.title=element_blank()) +
	xlab(paste0(genome1,' count')) +
	ylab(paste0(genome3,' count'))
ggsave(p,file=file.path('figures/contamination',paste0(this_id,'-',genome1,'_',genome3,'_barnyard.pdf')),useDingbats=FALSE)

p = ggplot(subset(meta,n_assigned >= 10),aes_string(paste0('frac_',genome2,'+frac_',genome3))) +
	geom_histogram(bins=25) +
	facet_wrap(~rt,ncol=2) +
	coord_cartesian(xlim=c(0,1)) +
	theme_classic() +
	theme(strip.background=element_blank()) +
	xlab(paste0('fraction reads ',genome2,'+',genome3)) +
	ylab('count')
ggsave(p,file=file.path('figures/contamination',paste0(this_id,'-',genome1,'_',genome2,'_',genome3,'_rt_wells.pdf')),useDingbats=FALSE)

dir.create('rds/contamination',showWarnings=FALSE)

saveRDS(meta,file=file.path('rds/contamination',paste0(this_id,'-',genome1,'_',genome2,'_',genome3,'.rds')))

