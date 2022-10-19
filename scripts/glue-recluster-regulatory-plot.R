#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

# arguments = c('biccn','subpeak')
prefix = arguments[1]

if (length(arguments) > 1) {
	prefix2 = arguments[2]
}
library(data.table)
library(parallel)
library(ggplot2)
library(egg)
library(ggrastr)

cell.classes = scan('stats/clusters/rna-final-cellclasses-levels.txt',what='',sep='\n',quiet=TRUE)

files1 = list.files(path='stats/glue/subintegration',pattern=paste0(prefix,'_class[0-9]+_regulatory.txt.gz'))

file.list = data.frame(files=files1,cell_class=factor(cell.classes[as.integer(gsub('.+?class([0-9]+).+','\\1',files1))],levels=cell.classes))
file.list = file.list[order(file.list$cell_class),]

peaks.binned = do.call(rbind,mclapply(1:nrow(file.list),function(i) {
	this.file = file.path('stats/glue/subintegration',file.list[i,'files'])
	this.class = file.list[i,'cell_class']
	x = fread(this.file)
	x$dist_category = factor(with(x,
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
	x$cell_class = this.class
	x$set = factor('before',levels=c('before','after'))
	x
},mc.cores=4))

if (exists('prefix2')) {
	files2 = list.files(path='stats/glue/subintegration',pattern=paste0(prefix2,'_class[0-9]+_regulatory.txt.gz'))

	file.list = data.frame(files=files2,cell_class=factor(cell.classes[as.integer(gsub('.+?class([0-9]+).+','\\1',files2))],levels=cell.classes))
	file.list = file.list[order(file.list$cell_class),]

	peaks.binned2 = do.call(rbind,mclapply(1:nrow(file.list),function(i) {
		this.file = file.path('stats/glue/subintegration',file.list[i,'files'])
		this.class = file.list[i,'cell_class']
		x = fread(this.file)
		x$dist_category = factor(with(x,
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
		x$cell_class = this.class
		x$set = factor('after',levels=c('before','after'))
		x
	},mc.cores=4))
	peaks.binned = rbind(peaks.binned,peaks.binned2)
}

if (exists('prefix2')) {
	p = ggplot(droplevels(peaks.binned),aes(dist_category,score,fill=set)) +
		geom_boxplot_jitter(outlier.size = 0.0001,size=0.5,outlier.jitter.width=0) +
		theme_article() +
		facet_wrap(~cell_class) +
		scale_fill_manual(values=c('#ff7800','#0077bd')) +
		# xlab('Genomic distance') +
		ylab('GLUE regulatory score') +
		theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=75, vjust = 1, hjust=1),legend.title=element_blank())
	ggsave(p,file=file.path('figures/glue',paste0(prefix,'_vs_',prefix2,'_regulatory_score.pdf')),useDingbats=FALSE,height=7,width=7)
} else {
	p = ggplot(droplevels(peaks.binned),aes(dist_category,score)) +
		geom_boxplot_jitter(fill='#ff7800',outlier.size = 0.0001,size=0.5,outlier.jitter.width=0) +
		theme_article() +
		facet_wrap(~cell_class) +
		# xlab('Genomic distance') +
		ylab('GLUE regulatory score') +
		theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=75, vjust = 1, hjust=1),legend.title=element_blank())
	ggsave(p,file=file.path('figures/glue',paste0(prefix,'_regulatory_score.pdf')),useDingbats=FALSE,height=7,width=7)
}