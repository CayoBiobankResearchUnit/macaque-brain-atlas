#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('atac','marker')

prefix = arguments[1] # atac for parent-levels, atacsub for class level
peak.set = arguments[2]

if (prefix == 'atac' & peak.set == 'marker') {
	files = file.path('stats/tf/homer/markerpeaks/class',list.files('stats/tf/homer/markerpeaks/class',pattern='class[0-9]+_marker_jasparTFs.txt'))
} else if (prefix == 'atac' & peak.set == 'regulatory') {
	files = file.path('stats/tf/homer/regulatory/class',list.files('stats/tf/homer/regulatory/class',pattern='allpeaks_inc_v_dec_jasparTFs.txt'))
} else if (prefix == 'atacsub' & peak.set == 'marker') {
	files = file.path('stats/tf/homer/markerpeaks/cluster',list.files('stats/tf/homer/markerpeaks/cluster',pattern='class[0-9]+_cluster[0-9]+_marker_jasparTFs.txt'))
} else if (prefix == 'atacsub' & peak.set == 'regulatory') {
	files = file.path('stats/tf/homer/regulatory/cluster',list.files('stats/tf/homer/regulatory/cluster',pattern='class[0-9]+_[di][en]c_jasparTFs.txt'))
}

homer.results = lapply(files,function(f) {
	# x = read.delim(f,col.names=c('motif','consensus','pval_homer','log_pval','qval_homer','n_target','pct_target','n_background','pct_background'))
	x = read.delim(f)
	total_target = as.integer(gsub('.+?\\.of\\.([0-9]+)\\.$','\\1',names(x)[6]))
	total_background = as.integer(gsub('.+?\\.of\\.([0-9]+)\\.$','\\1',names(x)[8]))
	colnames(x) = c('motif','consensus','pval_homer','log_pval','qval_homer','n_target','pct_target_homer','n_background','pct_background_homer')
	x$pct_target = x$n_target / total_target
	x$pct_background = x$n_background / total_background
	x$odds_ratio = x$pct_target / x$pct_background
	x$pval = exp(x$log_pval)
	x$qval = p.adjust(x$pval,'fdr')
	x
})
names(homer.results) = files

cell.levels = scan(file='stats/clusters/rna-final-cellclasses-levels.txt',what='',sep='\n')
if (prefix == 'atac' & peak.set == 'marker') {
	homer.results = droplevels(do.call(rbind,lapply(names(homer.results),function(i) {
		x = homer.results[[i]]
		this.class = as.integer(gsub('.+class([0-9]+)_.*','\\1',i))
		x$cell_class = factor(cell.levels[this.class],levels=cell.levels)
		x
	})))
} else if (prefix == 'atac' & peak.set == 'regulatory') {
} else if (prefix == 'atacsub' & peak.set == 'marker') {
} else if (prefix == 'atacsub' & peak.set == 'regulatory') {
}

jaspar.info = read.delim('data/jaspar2018_CORE_vertebrates_non-redundant.txt',header=FALSE,col.names=c('motif','motif_name'))

homer.results = merge(homer.results,jaspar.info,all.x=TRUE,all.y=FALSE)

homer.results = homer.results[order(homer.results$log_pval,homer.results$cell_class),]

tfs.view = unique(subset(do.call(rbind,lapply(split(homer.results,homer.results$cell_class),head,3)),qval_homer < 0.1)$motif)

homer.view = subset(homer.results, motif %in% tfs.view)

homer.view$motif_name = factor(homer.view$motif_name,levels=unique(subset(do.call(rbind,lapply(split(homer.results,homer.results$cell_class),head,3)),qval_homer < 0.1)$motif_name))


homer.view$frac_target = with(homer.view,pct_target / (pct_target + pct_background))
homer.view$frac_background = with(homer.view,pct_background / (pct_target + pct_background))

homer.target = homer.view[,c('pval_homer','qval_homer','pct_target','frac_target','odds_ratio','cell_class','motif_name')]
homer.background = homer.view[,c('pval_homer','qval_homer','pct_background','frac_background','odds_ratio','cell_class','motif_name')]


names(homer.target)[3:4] = names(homer.background)[3:4] = c('pct','frac')
homer.target$type = factor('target',levels=c('target','background'))
homer.background$type = factor('background',levels=c('target','background'))


homer.pie = rbind(homer.target,homer.background)

p = ggplot(homer.pie) +
	geom_bar(aes(x=1,y=frac,fill=type,alpha=factor(qval_homer < 0.05,levels=c('TRUE','FALSE'))),stat='identity') +
	facet_grid(cell_class ~ motif_name) +
	# scale_fill_viridis(option='D') +
	# scale_y_continuous(breaks=seq(-2,2,2)) +
	scale_fill_manual(values=c('#4daf4a','#cccccc')) +
	scale_alpha_manual(values=c(1,0.25)) +
	coord_polar('y',start=0) +
	theme_classic() +
	guides(alpha = 'none') +
	theme(
		legend.position='none',
		axis.text=element_blank(),
		axis.line=element_blank(),
		axis.ticks=element_blank(),
		axis.title=element_blank(),
		strip.background = element_blank(),
		strip.text.x = element_text(angle=45,hjust=0),
		strip.text.y = element_text(angle=0,hjust=0),
	) +
	ylab(expression(log[2]~'odds ratio'))

p = ggplot(homer.view) +
	geom_bar(aes(x=1,y=log2(odds_ratio),fill=-log10(pval_homer),alpha=factor(qval_homer < 0.05,levels=c('TRUE','FALSE'))),stat='identity') +
	facet_grid(cell_class ~ motif_name) +
	scale_fill_viridis(option='D') +
	scale_y_continuous(breaks=seq(-2,2,2)) +
	scale_alpha_manual(values=c(1,0.1)) +
	theme_classic() +
	guides(alpha = 'none') +
	theme(
		axis.text.x=element_blank(),
		axis.line.x=element_blank(),
		axis.ticks.x=element_blank(),
		axis.title.x=element_blank(),
		strip.background = element_blank(),
		strip.text.y = element_text(angle=0,hjust=0),
	) +
	ylab(expression(log[2]~'odds ratio'))
	