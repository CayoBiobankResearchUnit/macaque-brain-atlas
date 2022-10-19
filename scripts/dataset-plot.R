#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('rna','human-cortex')

prefix = arguments[1]
dataset_prefix = arguments[2]

cor.result = readRDS(file.path('rds/nnls/',paste0('nnls_results-',prefix,'_',dataset_prefix,'.rds')))

query.levels = scan(file='stats/subclusters/rna-cellsubclusters-levels.txt',what='',sep='\n',quiet=TRUE)
ref.levels = sort(unique(cor.result$target))

if (any(!query.levels %in% cor.result$source)) {
	# Temporary fixes
	query.levels = gsub('neuron$','neurons',query.levels)
	query.levels = gsub('serotinergic','serotonergic',query.levels)
}

query.classes = gsub(' [0-9]+$','',query.levels)

if (prefix == 'rna-class') {
	query.levels = unique(query.classes)
	cor.result$source = gsub('neuron$','neurons',cor.result$source)
	cor.result$source = gsub('serotinergic','serotonergic',cor.result$source)
}

# cor.matrix = matrix(0,nrow=length(query.levels),ncol=length(ref.levels),dimnames=list(query.levels,ref.levels))
# 
# cor.matrix[as.matrix(cor.result[,c('source','target')])] = cor.result$beta_1

cor.result$source = factor(cor.result$source,levels=query.levels,labels=gsub('inhibitory','GABAergic',gsub('excitatory','glutamatergic',query.levels)))
cor.result$target = factor(cor.result$target,levels=ref.levels)

query.levels = gsub('inhibitory','GABAergic',gsub('excitatory','glutamatergic',query.levels))
query.classes = gsub('inhibitory','GABAergic',gsub('excitatory','glutamatergic',query.classes))

# cortex-dev
# human-cortex
# dev-brain-regions
# adult-brain-vasc-peri
# adult-brain-vasc-endo
# vascular-dev
# early-brain
# brain-vasc-atlas
# mouse-drg-injury
# dev-inhibitory-neurons-macaque

library(ggplot2)
library(viridis)

dir.create('figures/nnls',showWarnings=FALSE)

x.font.size = if (length(query.levels) > 100) 4 else if (length(query.levels) < 20) 8
y.font.size = if (length(ref.levels) > 75) {
	4
} else if (length(ref.levels) > 50) {
	6
} else if (length(ref.levels) < 15) {
	10
} else {
	8
}

p = ggplot(cor.result,aes(source,target,fill=beta_2)) +
	geom_tile() +
	scale_fill_viridis(name='NNLS',option='D') +
# 	coord_fixed() +
	theme_classic() +
	theme(
		axis.line=element_blank(),
		axis.title=element_blank(),
		axis.text.x = element_text(vjust=1,hjust=0,angle=-30,size=x.font.size),
		axis.text.y = element_text(size=y.font.size)
	)
ggsave(p,file=file.path('figures/nnls',paste0('nnls-',prefix,'_',dataset_prefix,'_all_direction1.pdf')),useDingbats=FALSE)

p = ggplot(cor.result,aes(source,target,fill=beta_1)) +
	geom_tile() +
	scale_fill_viridis(name='NNLS',option='D') +
#	coord_fixed() +
	theme_classic() +
	theme(
		axis.line=element_blank(),
		axis.title=element_blank(),
		axis.text.x = element_text(vjust=1,hjust=0,angle=-30,size=x.font.size),
		axis.text.y = element_text(size=y.font.size)
	)
ggsave(p,file=file.path('figures/nnls',paste0('nnls-',prefix,'_',dataset_prefix,'_all_direction2.pdf')),useDingbats=FALSE)

if (dataset_prefix %in% c('vascular-dev','brain-vasc-atlas','adult-brain-vasc-endo','adult-brain-vasc-peri','myeloid-neuroinflam')) {
	# keep.levels = query.levels[query.classes %in% c('astrocytes','oligodendrocytes','oligodendrocyte precursor cells','vascular cells','microglia','radial glial cells','mesenchymal stem cells','ependymal cells')]
	keep.list = c('vascular cells','microglia','radial glial cells','mesenchymal stem cells','ependymal cells')
	keep.levels = query.levels[query.classes %in% keep.list]
} else if (dataset_prefix %in% c('cortex-dev','dev-brain-regions','early-brain','dev-inhibitory-neurons-macaque')) {
	keep.list = c('oligodendrocyte precursor cells','radial glial cells','mesenchymal stem cells')
	keep.levels = query.levels[query.classes %in% keep.list]
} else if (dataset_prefix %in% c('mouse-drg-injury')) {
	keep.list = c('astrocytes','oligodendrocytes','oligodendrocyte precursor cells','radial glial cells','mesenchymal stem cells')
	keep.levels = query.levels[query.classes %in% keep.list]
} else if (dataset_prefix %in% c('human-cortex')) {
	keep.list = c('glutamatergic neurons','GABAergic neurons','dopaminergic neurons','serotonergic neurons')
	keep.levels = query.levels[query.classes %in% keep.list]
} else if (dataset_prefix %in% c('fang-merfish-l2','fang-merfish-l3')) {
	keep.list = c('glutamatergic neurons','GABAergic neurons','dopaminergic neurons','serotonergic neurons')
	keep.levels = query.levels[query.classes %in% keep.list]
} else {
	keep.levels = query.levels
}

if (prefix != 'rna-class') {
	cor.filtered = droplevels(subset(cor.result,source %in% keep.levels))

	x.font.size = if (nlevels(cor.filtered$source) > 75) {
		4
	} else if (nlevels(cor.filtered$source) > 50) {
		6
	} else if (nlevels(cor.filtered$source) < 15) {
		10
	} else {
		8
	}

	p = ggplot(cor.filtered,aes(source,target,fill=beta_2)) +
		geom_tile() +
		scale_fill_viridis(name='NNLS',option='D') +
	#	coord_fixed() +
		theme_classic() +
		theme(
			axis.line=element_blank(),
			axis.title=element_blank(),
			axis.text.x = element_text(vjust=1,hjust=0,angle=-30,size=x.font.size),
			axis.text.y = element_text(size=y.font.size)
		)
	ggsave(p,file=file.path('figures/nnls',paste0('nnls-',prefix,'_',dataset_prefix,'_filtered_direction1.pdf')),useDingbats=FALSE)
	
	if (dataset_prefix %in% c('fang-merfish-l2','fang-merfish-l3')) {

		p = ggplot(cor.result,aes(source,target,fill=ifelse(beta_2>1,1,beta_2))) +
			geom_tile() +
			scale_fill_viridis(name='NNLS',option='D') +
		# 	coord_fixed() +
			theme_classic() +
			theme(
				axis.line=element_blank(),
				axis.title=element_blank(),
				axis.text.x = element_text(vjust=1,hjust=0,angle=-30,size=x.font.size),
				axis.text.y = element_text(size=y.font.size)
			)
		ggsave(p,file=file.path('figures/nnls',paste0('nnls-',prefix,'_',dataset_prefix,'_all_direction1-clipped.pdf')),useDingbats=FALSE)

		p = ggplot(cor.filtered,aes(source,target,fill=ifelse(beta_2>1,1,beta_2))) +
			geom_tile() +
			scale_fill_viridis(name='NNLS',option='D') +
		#	coord_fixed() +
			theme_classic() +
			theme(
				axis.line=element_blank(),
				axis.title=element_blank(),
				axis.text.x = element_text(vjust=1,hjust=0,angle=-30,size=x.font.size),
				axis.text.y = element_text(size=y.font.size)
			)
		ggsave(p,file=file.path('figures/nnls',paste0('nnls-',prefix,'_',dataset_prefix,'_filtered_direction1-clipped.pdf')),useDingbats=FALSE)
	}

	p = ggplot(cor.filtered,aes(source,target,fill=beta_1)) +
		geom_tile() +
		scale_fill_viridis(name='NNLS',option='D') +
	#	coord_fixed() +
		theme_classic() +
		theme(
			axis.line=element_blank(),
			axis.title=element_blank(),
			axis.text.x = element_text(vjust=1,hjust=0,angle=-30,size=x.font.size),
			axis.text.y = element_text(size=y.font.size)
		)
	ggsave(p,file=file.path('figures/nnls',paste0('nnls-',prefix,'_',dataset_prefix,'_filtered_direction2.pdf')),useDingbats=FALSE)
}
