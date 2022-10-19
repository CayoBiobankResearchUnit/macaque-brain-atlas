#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('atac','atac','1')
prefix = arguments[1]
analysis = arguments[2]

if (!analysis %in% c('rna','atac')) {
	stop('Argument 2 must be either "rna" or "atac"')
}

suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(parallel))
suppressMessages(library(egg))

sample.list = as.character(read.delim(file.path('data',paste0(prefix,'_metadata.txt')))$id)

n.cores = 4

if (!file.exists(file.path('data',paste0(prefix,'_scrublet_parameters.txt')))) {
	dev.null = mclapply(sample.list,function(i) {
		params = read.delim(paste0('stats/scrublet/scrublet_parameters_',prefix,'-',i,'.txt'))
		sim = fread(paste0('stats/scrublet/scrublet_score_distribution_',prefix,'-',i,'-sim.txt.gz'))
		obs = fread(paste0('stats/scrublet/scrublet_score_distribution_',prefix,'-',i,'-obs.txt.gz'))
		scr = rbind(data.frame(type='Simulated doublets',values=unlist(sim)),data.frame(type='Observed indices',values=unlist(obs)))
		p = ggplot(scr,aes(values)) +
			geom_histogram(bins=50) +
			geom_vline(xintercept=params$threshold) +
			facet_wrap(~type,nrow=1,scales='free_y') +
			scale_x_continuous(breaks=seq(0,1,0.1),minor_breaks=seq(0,1,0.01),limits=c(0,1)) +
	#		coord_cartesian(xlim=c(0,1)) +
			theme_article() +
			theme(
				axis.text.y=element_blank(),
				axis.ticks.y=element_blank(),
				panel.grid.major.x=element_line(size=0.2,color='gray'),
				panel.grid.minor.x=element_line(size=0.1,color='gray')) +
			xlab('Doublet score') +
			ylab('Prob. density')
		ggsave(p,file=file.path('figures/scrublet',paste0('scrublet_score_distribution2_',prefix,'-',i,'.pdf')),useDingbats=FALSE,width=7,height=2)
	},mc.cores=n.cores)
	
	all.parameters = do.call(rbind,lapply(sample.list,function(i) {
		read.delim(paste0('stats/scrublet/scrublet_parameters_',prefix,'-',i,'.txt'))
	}))
	
	# Create a new column for manually set thresholds
	all.parameters$manual_threshold = ''
	
	write.table(all.parameters,file=file.path('data',paste0(prefix,'_scrublet_parameters.txt')),sep='\t',quote=FALSE,row.names=FALSE)
} else {
	all.parameters = read.delim(file.path('data',paste0(prefix,'_scrublet_parameters.txt')))
	all.parameters$final_threshold = with(all.parameters,ifelse(manual_threshold=='auto',threshold,as.numeric(manual_threshold)))

	dev.null = mclapply(sample.list,function(i) {
		params = subset(all.parameters,id == i)
		new.threshold = with(params,ifelse(threshold == final_threshold,FALSE,TRUE))
		sim = fread(paste0('stats/scrublet/scrublet_score_distribution_',prefix,'-',i,'-sim.txt.gz'))
		obs = fread(paste0('stats/scrublet/scrublet_score_distribution_',prefix,'-',i,'-obs.txt.gz'))
		scr = rbind(data.frame(type='Simulated doublets',values=unlist(sim)),data.frame(type='Observed indices',values=unlist(obs)))
		p = ggplot(scr,aes(values)) +
			geom_histogram(bins=50) +
			geom_vline(xintercept=params$threshold) +
			geom_vline(xintercept=params$final_threshold,color=ifelse(new.threshold,'red',NA)) +
			facet_wrap(~type,nrow=1,scales='free_y') +
			scale_x_continuous(breaks=seq(0,1,0.1),minor_breaks=seq(0,1,0.01),limits=c(0,1)) +
	#		coord_cartesian(xlim=c(0,1)) +
			theme_article() +
			theme(
				axis.text.y=element_blank(),
				axis.ticks.y=element_blank(),
				panel.grid.major.x=element_line(size=0.2,color='gray'),
				panel.grid.minor.x=element_line(size=0.1,color='gray')) +
			xlab('Doublet score') +
			ylab('Prob. density')
		ggsave(p,file=file.path('figures/scrublet',paste0('scrublet_score_distribution3_',prefix,'-',i,'.pdf')),useDingbats=FALSE,width=7,height=2)
	},mc.cores=n.cores)
	write.table(all.parameters,file=file.path('data',paste0(prefix,'_scrublet_parameters_final.txt')),sep='\t',quote=FALSE,row.names=FALSE)
}