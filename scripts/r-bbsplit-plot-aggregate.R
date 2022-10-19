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
# arguments = c('rna','macaque','human','mouse')

prefix = arguments[1]
genome1 = arguments[2]
genome2 = arguments[3]
genome3 = arguments[4]


if (!file.exists('rds/contamination/all-macaque_human_mouse.rds')) {
	files = file.path('rds/contamination',list.files('rds/contamination',pattern=paste0(genome1,'_',genome2,'_',genome3,'.rds')))

	meta = do.call(rbind,parallel::mclapply(files,readRDS,mc.cores=n.cores))

	saveRDS(meta,file='rds/contamination/all-macaque_human_mouse.rds')
} else {
	meta = readRDS('rds/contamination/all-macaque_human_mouse.rds')
}

colors = as.matrix(read.table('data/colors.txt',row.names=1,sep='\t',comment.char='',header=FALSE))[,1]
cell.levels = c('excitatory neurons','cerebellar neurons','inhibitory neurons','basket cells','medium spiny neurons','dopaminergic neurons','serotonergic neurons','AHSG neurons','F5 neurons','KIR3DL12 neurons','astrocytes','oligodendrocytes','oligodendrocyte precursor cells','radial glial cells','ependymal cells','microglia','KIR3DL12 microglia','vascular cells','mesenchymal stem cells')

colors = c(colors,c('all cells' = '#000000'))
colors = c(colors,c('unclassified' = '#cccccc'))

names(colors) = gsub('neuron$','neurons',gsub('serotinergic','serotonergic',names(colors)))

dir.create('figures/contamination',showWarnings=FALSE)

# length(table(meta$rt))

meta.filtered = subset(meta,!is.na(n_assigned) & n_assigned >= 10)

meta.filtered$rt_index = meta.filtered$rt

meta.filtered$rt = with(meta.filtered,paste(id,rt,sep='-'))


median.rt.cells = median(table(meta.filtered$rt))

set.seed(42)
average.distribution = sample(with(meta.filtered,frac_human+frac_mouse),median.rt.cells,replace=FALSE)
average.contaminated = with(meta.filtered,sum((frac_human + frac_mouse) >= 0.75)) / nrow(meta.filtered)


rt.results = do.call(rbind,mclapply(split(meta.filtered,meta.filtered$rt),function(x) {
	ks.result = ks.test(with(x,frac_human+frac_mouse),average.distribution,alternative='less')
	bt.result = binom.test(with(x,sum((frac_human+frac_mouse) >= 0.75)),nrow(x),p=average.contaminated,alternative='greater')
	data.frame(
		rt=unique(x$rt),
		n_cells_assigned=nrow(x),
		exo_fraction=with(x,sum((frac_human+frac_mouse) >= 0.75)/nrow(x)),
		exo_median=with(x,median(frac_human+frac_mouse)),
		ks_stat=ks.result$statistic,
		ks_pval=ks.result$p.value,
		bt_coef=bt.result$estimate,
		bt_pval=bt.result$p.value
	)
},mc.cores=4))

rt.results$ks_padj = p.adjust(rt.results$ks_pval,'fdr')
rt.results$bt_padj = p.adjust(rt.results$bt_pval,'fdr')

meta.filtered$rt_ordered_ks = factor(meta.filtered$rt,levels=rt.results[order(rt.results$ks_stat,decreasing=TRUE),]$rt)
meta.filtered$rt_ordered_bt = factor(meta.filtered$rt,levels=rt.results[order(-rt.results$bt_pval,rt.results$bt_coef,decreasing=TRUE),]$rt)

p = ggplot(subset(meta.filtered,cell_class != 'unclassified'),aes_string(paste0('n_',genome1),paste0('n_',genome2),color='cell_class')) +
	geom_point_rast(size=0.5,alpha=0.5) +
	coord_cartesian() +
	facet_wrap(~id,scales='free') +
	scale_color_manual(values=colors[levels(meta$cell_class)]) +
	guides(color=guide_legend(override.aes=list(size=2,alpha=1))) +
	theme_classic() +
	theme(
		legend.position = 'right',
		legend.title=element_blank(),
		strip.background=element_blank(),
		strip.text = element_text(size=4),
		axis.text = element_text(size=4)
	) +
	xlab(paste0(genome1,' count')) +
	ylab(paste0(genome2,' count'))
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_barnyard.pdf')),useDingbats=FALSE,width=15,height=15)

p = ggplot(subset(meta.filtered,cell_class != 'unclassified'),aes_string(paste0('n_',genome1),paste0('n_',genome3),color='cell_class')) +
	geom_point_rast(size=0.5,alpha=0.5) +
	coord_cartesian() +
	facet_wrap(~id,scales='free') +
	scale_color_manual(values=colors[levels(meta$cell_class)]) +
	guides(color=guide_legend(override.aes=list(size=2,alpha=1))) +
	theme_classic() +
	theme(
		legend.position = 'right',
		legend.title=element_blank(),
		strip.background=element_blank(),
		strip.text = element_text(size=4),
		axis.text = element_text(size=4)
	) +
	xlab(paste0(genome1,' count')) +
	ylab(paste0(genome3,' count'))
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome3,'_barnyard.pdf')),useDingbats=FALSE,width=15,height=15)

p = ggplot(subset(meta.filtered,cell_class != 'unclassified'),aes_string(paste0('n_',genome1),paste0('n_',genome2,'+n_',genome3),color='cell_class')) +
	geom_point_rast(size=0.5,alpha=0.5) +
	coord_cartesian() +
	facet_wrap(~id,scales='free') +
	scale_color_manual(values=colors[levels(meta$cell_class)]) +
	guides(color=guide_legend(override.aes=list(size=2,alpha=1))) +
	theme_classic() +
	theme(
		legend.position = 'right',
		legend.title=element_blank(),
		strip.background=element_blank(),
		strip.text = element_text(size=4),
		axis.text = element_text(size=4)
	) +
	xlab(paste0(genome1,' count')) +
	ylab(paste0(genome2,'/',genome3,' count'))
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_barnyard.pdf')),useDingbats=FALSE,width=15,height=15)


p = ggplot(meta.filtered,aes_string(paste0('frac_',genome2,'+frac_',genome3))) +
	geom_histogram(binwidth=0.1) +
	facet_wrap(~rt,scales='free_y') +
	theme_classic() +
	theme(strip.background=element_blank(),strip.text=element_text(size=6),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=4)) +
	xlab(paste0('fraction reads ',genome2,'/',genome3)) +
	ylab('count')
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_rt_wells.pdf')),useDingbats=FALSE,height=15,width=15)

p = ggplot(meta.filtered,aes_string(paste0('frac_',genome2,'+frac_',genome3))) +
	geom_histogram(binwidth=0.1) +
	facet_wrap(~rt_ordered_ks,scales='free_y') +
	theme_classic() +
	theme(strip.background=element_blank(),strip.text=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=4)) +
	xlab(paste0('fraction reads ',genome2,'/',genome3)) +
	ylab('count')
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_rt_wells_ordered_ks.pdf')),useDingbats=FALSE,height=15,width=15)

p = ggplot(meta.filtered,aes_string(paste0('frac_',genome2,'+frac_',genome3))) +
	geom_histogram(binwidth=0.1) +
	facet_wrap(~rt_ordered_bt,scales='free_y') +
	theme_classic() +
	theme(strip.background=element_blank(),strip.text=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=4)) +
	xlab(paste0('fraction reads ',genome2,'/',genome3)) +
	ylab('count')
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_rt_wells_ordered_bt.pdf')),useDingbats=FALSE,height=15,width=15)

rt.bad = subset(rt.results,bt_padj < 0.05 & ks_padj < 0.05)

p = ggplot(subset(meta.filtered,rt %in% rt.bad$rt),aes_string(paste0('frac_',genome2,'+frac_',genome3))) +
	geom_histogram(binwidth=0.1) +
	facet_wrap(~rt_ordered_ks,scales='free_y') +
	theme_classic() +
	theme(strip.background=element_blank(),strip.text=element_text(size=8),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=4)) +
	xlab(paste0('fraction reads ',genome2,'/',genome3)) +
	ylab('count')
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_rt_wells_ordered_ks_bad.pdf')),useDingbats=FALSE,height=15,width=15)

p = ggplot(subset(meta.filtered,rt %in% rt.bad$rt),aes_string(paste0('frac_',genome2,'+frac_',genome3))) +
	geom_histogram(binwidth=0.1) +
	facet_wrap(~rt_ordered_bt,scales='free_y') +
	theme_classic() +
	theme(strip.background=element_blank(),strip.text=element_text(size=8),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=4)) +
	xlab(paste0('fraction reads ',genome2,'/',genome3)) +
	ylab('count')
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_rt_wells_ordered_bt_bad.pdf')),useDingbats=FALSE,height=15,width=15)





p = ggplot(meta.filtered,aes_string(paste0('frac_',genome2,'+frac_',genome3))) +
	geom_histogram(binwidth=0.1) +
	facet_wrap(~cell_class,scales='free_y') +
	theme_classic() +
	theme(strip.background=element_blank(),strip.text=element_text(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text()) +
	xlab(paste0('fraction reads ',genome2,'/',genome3)) +
	ylab('count')
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_by_class.pdf')),useDingbats=FALSE,height=7,width=7)

meta.long = tidyr::pivot_longer(meta.filtered,cols=c(paste0('frac_',c(genome1,genome2,genome3))))

meta.long$name = factor(meta.long$name,levels=paste0('frac_',c(genome1,genome2,genome3)),labels=paste0('fraction ',c(genome1,genome2,genome3)))

p = ggplot(meta.long,aes_string('value')) +
	geom_histogram(binwidth=0.05) +
	facet_grid(cell_class~name,scales='free_y') +
	theme_classic() +
	theme(
		strip.background=element_blank(),
		strip.text=element_text(),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank(),
		axis.text.x=element_text(),
		strip.text.y=element_text(angle=0,hjust=0),
	) +
	xlab(paste0('fraction reads')) +
	ylab('count')
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_by_class_facet.pdf')),useDingbats=FALSE,height=7,width=7)





p = ggplot(meta.filtered,aes_string(paste0('frac_',genome2))) +
	geom_histogram(binwidth=0.1) +
	facet_wrap(~rt_ordered_bt,scales='free_y') +
	theme_classic() +
	theme(strip.background=element_blank(),strip.text=element_text(size=8),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=4)) +
	xlab(paste0('fraction reads ',genome2)) +
	ylab('count')
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_rt_wells.pdf')),useDingbats=FALSE,height=25,width=25)

p = ggplot(meta.filtered,aes_string(paste0('frac_',genome3))) +
	geom_histogram(binwidth=0.1) +
	facet_wrap(~rt_ordered_bt,scales='free_y') +
	theme_classic() +
	theme(strip.background=element_blank(),strip.text=element_text(size=8),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=4)) +
	xlab(paste0('fraction reads ',genome3)) +
	ylab('count')
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome3,'_rt_wells.pdf')),useDingbats=FALSE,height=25,width=25)


# Now redo RT wells plot but with only showing previously assigned cell classes
p = ggplot(subset(meta.filtered,rt %in% rt.bad$rt & !cell_class %in% c('unclassified')),aes_string(paste0('frac_',genome2,'+frac_',genome3))) +
	geom_histogram(binwidth=0.1) +
	facet_wrap(~rt_ordered_ks,scales='free_y') +
	theme_classic() +
	theme(strip.background=element_blank(),strip.text=element_text(size=8),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=4)) +
	xlab(paste0('fraction reads ',genome2,'/',genome3)) +
	ylab('count')
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_rt_wells_ordered_ks_assignedcells_bad.pdf')),useDingbats=FALSE,height=15,width=15)

p = ggplot(subset(meta.filtered,rt %in% rt.bad$rt & !cell_class %in% c('unclassified')),aes_string(paste0('frac_',genome2,'+frac_',genome3))) +
	geom_histogram(binwidth=0.1) +
	facet_wrap(~rt_ordered_bt,scales='free_y') +
	theme_classic() +
	theme(strip.background=element_blank(),strip.text=element_text(size=8),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=4)) +
	xlab(paste0('fraction reads ',genome2,'/',genome3)) +
	ylab('count')
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_rt_wells_ordered_bt_assigned_cells_bad.pdf')),useDingbats=FALSE,height=15,width=15)

# Now redo RT wells plot but excluding mouse/human clusters
p = ggplot(subset(meta.filtered,rt %in% rt.bad$rt & !cell_class %in% c('unclassified','radial glial cells','mesenchymal stem cells')),aes_string(paste0('frac_',genome2,'+frac_',genome3))) +
	geom_histogram(binwidth=0.1) +
	facet_wrap(~rt_ordered_ks,scales='free_y') +
	theme_classic() +
	theme(strip.background=element_blank(),strip.text=element_text(size=8),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=4)) +
	xlab(paste0('fraction reads ',genome2,'/',genome3)) +
	ylab('count')
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_rt_wells_ordered_ks_assignedcells_endogenous_bad.pdf')),useDingbats=FALSE,height=15,width=15)

p = ggplot(subset(meta.filtered,rt %in% rt.bad$rt & !cell_class %in% c('unclassified','radial glial cells','mesenchymal stem cells')),aes_string(paste0('frac_',genome2,'+frac_',genome3))) +
	geom_histogram(binwidth=0.1) +
	facet_wrap(~rt_ordered_bt,scales='free_y') +
	theme_classic() +
	theme(strip.background=element_blank(),strip.text=element_text(size=8),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=4)) +
	xlab(paste0('fraction reads ',genome2,'/',genome3)) +
	ylab('count')
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_rt_wells_ordered_bt_assigned_cells_endogenous_bad.pdf')),useDingbats=FALSE,height=15,width=15)



# Focus on RT wells of one individual
p = ggplot(subset(meta.filtered,id %in% 'NSM345' & !cell_class %in% c('unclassified')),aes_string(paste0('frac_',genome2,'+frac_',genome3))) +
	geom_histogram(binwidth=0.1) +
	facet_wrap(~rt_index,scales='free_y',ncol=1) +
	theme_classic() +
	theme(strip.background=element_blank(),strip.text=element_text(),axis.text.y=element_text(),axis.text.x=element_text()) +
	xlab(paste0('fraction reads ',genome2,'/',genome3)) +
	ylab('count')
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_rt_wells_NSM345_all.pdf')),useDingbats=FALSE,height=7,width=7)

# Focus on RT wells of one individual
p = ggplot(subset(meta.filtered,id %in% 'NSM345' & !cell_class %in% c('unclassified','radial glial cells','mesenchymal stem cells')),aes_string(paste0('frac_',genome2,'+frac_',genome3))) +
	geom_histogram(binwidth=0.1) +
	facet_wrap(~rt_index,scales='free_y',ncol=1) +
	theme_classic() +
	theme(strip.background=element_blank(),strip.text=element_text(),axis.text.y=element_text(),axis.text.x=element_text()) +
	xlab(paste0('fraction reads ',genome2,'/',genome3)) +
	ylab('count')
ggsave(p,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_rt_wells_NSM345_endogenous.pdf')),useDingbats=FALSE,height=7,width=7)


nsm.345 = rbind(
	data.frame(subset(meta.filtered,id %in% 'NSM345' & !cell_class %in% c('unclassified')),stage=factor('before',levels=c('before','after'))),
	data.frame(subset(meta.filtered,id %in% 'NSM345' & !cell_class %in% c('unclassified','radial glial cells','mesenchymal stem cells')),stage=factor('after',levels=c('before','after')))
)


# Focus on RT wells of one individual
p0 = ggplot(nsm.345,aes_string(paste0('frac_',genome2,'+frac_',genome3))) +
	geom_histogram(binwidth=0.05) +
	facet_grid(rt_index~stage,scales='free_y') +
	theme_classic(base_size=16) +
	theme(strip.background=element_blank(),strip.text.x=element_text(),strip.text.y=element_text(angle=0,hjust=0),axis.text.y=element_text(),axis.text.x=element_text()) +
	xlab(paste0('fraction reads ',genome2,'/',genome3)) +
	ylab('count')
ggsave(p0,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_rt_wells_NSM345_comparison.pdf')),useDingbats=FALSE,height=7,width=7)

nsm.346 = rbind(
	data.frame(subset(meta,id %in% 'NSM346' & !cell_class %in% c('unclassified')),stage=factor('before',levels=c('before','after'))),
	data.frame(subset(meta,id %in% 'NSM346' & !cell_class %in% c('unclassified','radial glial cells','mesenchymal stem cells')),stage=factor('after',levels=c('before','after')))
)

# Focus on RT wells of one individual
p00 = ggplot(nsm.346,aes_string(paste0('frac_',genome2,'+frac_',genome3))) +
	geom_histogram(binwidth=0.05) +
	facet_grid(rt~stage,scales='free_y') +
	theme_classic(base_size=16) +
	theme(strip.background=element_blank(),strip.text.x=element_text(),strip.text.y=element_text(angle=0,hjust=0),axis.text.y=element_text(),axis.text.x=element_text()) +
	xlab(paste0('fraction reads ',genome2,'/',genome3)) +
	ylab('count')
ggsave(p00,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_rt_wells_NSM346_comparison.pdf')),useDingbats=FALSE,height=7,width=7)


meta.long = tidyr::pivot_longer(meta.filtered,cols=c(paste0('frac_',c(genome1,genome2,genome3))))

meta.long$name = factor(meta.long$name,levels=paste0('frac_',c(genome1,genome2,genome3)),labels=paste0('fraction ',c(genome1,genome2,genome3)))

meta.long$cell_ordered = factor(meta.long$cell_class,levels=c(
'radial glial cells',
'mesenchymal stem cells',
'excitatory neurons',
'inhibitory neurons',
'cerebellar neurons',
'oligodendrocytes',
'astrocytes',
'vascular cells',
'oligodendrocyte precursor cells',
'medium spiny neurons',
'microglia',
'basket cells',
'ependymal cells',
'dopaminergic neurons',
'serotonergic neurons',
'AHSG neurons',
'F5 neurons',
'KIR3DL12 neurons',
'KIR3DL12 microglia',
'unclassified'
),labels=c(
'unknown cluster 1 (ASPM/CENPF)',
'unknown cluster 2 (COL1A1/COL1A2)',
'excitatory neurons',
'inhibitory neurons',
'cerebellar neurons',
'oligodendrocytes',
'astrocytes',
'vascular cells',
'oligodendrocyte precursor cells',
'medium spiny neurons',
'microglia',
'basket cells',
'ependymal cells',
'dopaminergic neurons',
'serotonergic neurons',
'AHSG neurons',
'F5 neurons',
'KIR3DL12 neurons',
'KIR3DL12 microglia',
'unannotated'
))

p1 = ggplot(meta.long,aes_string('value')) +
	geom_histogram(binwidth=0.05) +
	facet_grid(cell_ordered~name,scales='free_y') +
	scale_x_continuous(breaks=seq(0,1,0.5)) +
	theme_classic(base_size=12) +
	theme(
		strip.background=element_blank(),
		strip.text=element_text(),
		axis.text.y=element_text(),
		axis.text.x=element_text(),
		strip.text.y=element_text(angle=0,hjust=0),
	) +
	xlab(paste0('fraction reads')) +
	ylab('count (thousands)')
ggsave(p1,file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_by_class_facet_reorganized.pdf')),useDingbats=FALSE,height=7,width=7)

library(egg)

pdf(file='/dev/null')
p = egg::ggarrange(p0,p1,nrow=1)
dev.off()

pdf(file=file.path('figures/contamination',paste0('all-',genome1,'_',genome2,'_',genome3,'_qc_final.pdf')),useDingbats=FALSE,height=7,width=7)
p
dev.off()