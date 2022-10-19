#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('biccn','rna','subpeak')

prefix = arguments[1]
rn.prefix = arguments[2]

if (length(arguments) > 2) {
	prefix2 = arguments[3]
}

suppressMessages(library(Matrix))
suppressMessages(library(data.table))

step1 = file.path('stats/glue',paste0(prefix,'_dx_consistency.txt'))

step2 = file.path('stats/glue/subintegration',list.files(path='stats/glue/subintegration',pattern=paste0(prefix,'_dx_consistency_class[0-9]+\\.txt')))
step2 = step2[order(as.integer(gsub('.+?class([0-9]+)\\.txt','\\1',step2)))]

cell.levels = scan(file.path('stats/clusters',paste0(rn.prefix,'-final-cellclasses-levels.txt')),what='',sep='\n',quiet=TRUE)

step2 = data.frame(
	cell.i = as.integer(gsub('.+?class([0-9]+)\\.txt','\\1',step2)),
	cell.class=cell.levels[as.integer(gsub('.+?class([0-9]+)\\.txt','\\1',step2))],
	file = step2)

library(ggplot2)
library(egg)

if (file.exists(step1)) {
	in_metrics = read.delim(step1,row.names=1)

	dir.create('figures/glue',showWarnings=FALSE)

	p = ggplot(in_metrics,aes(n_meta,consistency)) +
		geom_path(color='#0078b9') +
		geom_hline(yintercept=0.05,color='#961a1d',linetype=2) +
		theme_article(base_size=24) +
		scale_y_continuous(
			limits=c(0,max(in_metrics$consistency)),
			breaks=seq(0,max(in_metrics$consistency),0.1)) +
		theme(panel.grid.major=element_line(size=0.1,color='#cccccc'))
	ggsave(p,file=file.path('figures/glue',paste0(prefix,'_dx_consistency.pdf')))
}

in_metrics_all = do.call(rbind,lapply(1:nrow(step2),function(i) {
	in_metrics = read.delim(step2$file[i],row.names=1)
	data.frame(
		cell.class = factor(step2$cell.class[i],levels=cell.levels),
		in_metrics)
}))

if (exists('prefix2')) {
	step3 = file.path('stats/glue/subintegration',list.files(path='stats/glue/subintegration',pattern=paste0(prefix2,'_dx_consistency_class[0-9]+\\.txt')))
	step3 = step3[order(as.integer(gsub('.+?class([0-9]+)\\.txt','\\1',step3)))]
	
	step3 = data.frame(
		cell.i = as.integer(gsub('.+?class([0-9]+)\\.txt','\\1',step3)),
		cell.class=cell.levels[as.integer(gsub('.+?class([0-9]+)\\.txt','\\1',step3))],
		file = step3)
	
	in_metrics_all = rbind(
		data.frame(in_metrics_all,set=factor('before',levels=c('before','after'))),
		data.frame(do.call(rbind,lapply(1:nrow(step3),function(i) {
			in_metrics = read.delim(step3$file[i],row.names=1)
			data.frame(
				cell.class = factor(step3$cell.class[i],levels=cell.levels),
				in_metrics)
		})),set=factor('after',levels=c('before','after')))
	)
	
	p = ggplot(in_metrics_all,aes(n_meta,consistency,color=set)) +
		geom_path() +
		geom_hline(yintercept=0.05,color='#3c8048',linetype=2) +
		facet_wrap(~cell.class) +
		theme_article() +
		scale_color_manual(values=c('#0078b9','#b45656')) +
		scale_y_continuous(
			limits=c(0,max(in_metrics$consistency)),
			breaks=seq(0,max(in_metrics$consistency),0.1)) +
		theme(panel.grid.major=element_line(size=0.1,color='#cccccc'),legend.title=element_blank())
	ggsave(p,file=file.path('figures/glue',paste0(prefix,'_vs_',prefix2,'_dx_consistency_all.pdf')),width=10)
	
} else {
	p = ggplot(in_metrics_all,aes(n_meta,consistency)) +
		geom_path(color='#0078b9') +
		geom_hline(yintercept=0.05,color='#961a1d',linetype=2) +
		facet_wrap(~cell.class) +
		theme_article() +
		scale_y_continuous(
			limits=c(0,max(in_metrics$consistency)),
			breaks=seq(0,max(in_metrics$consistency),0.1)) +
		theme(panel.grid.major=element_line(size=0.1,color='#cccccc'))
	ggsave(p,file=file.path('figures/glue',paste0(prefix,'_dx_consistency_all.pdf')))
}