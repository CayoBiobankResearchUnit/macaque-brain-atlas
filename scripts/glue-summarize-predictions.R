#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('biccn','rna','atac','1')

prefix = arguments[1]
rn.prefix = arguments[2]
at.prefix = arguments[3]
do.matrix = as.integer(arguments[4])

glue.confidence.threshold = 0.95

if (do.matrix == 1) {
	query.prefix = at.prefix
	validation.column = 'atactype'
} else if (do.matrix == 2) {
	query.prefix = paste0(rn.prefix,'ext')
	validation.column = 'rnatype'
} else if (do.matrix == 3) {
	query.prefix = paste0(rn.prefix,'int')
	validation.column = 'rnatype'
} else {
	stop('Error')
}

these.predictions = read.delim(file.path('stats/clusters',paste0(query.prefix,'-celltype-predictions.txt')),row.names=1)

these.predictions = these.predictions[,c('id','animal_id','n.umi','region','modality',validation.column,'glue_type','glue_type_confidence')]

cell.levels = scan(file=file.path('stats/clusters',paste0(rn.prefix,'-final-cellclasses-levels.txt')),what='',sep='\n')

these.predictions$region = factor(these.predictions$region,levels=region.levels)
these.predictions$glue_type = factor(these.predictions$glue_type,levels=cell.levels)

if (do.matrix > 1) {
	these.predictions[[validation.column]] = factor(these.predictions[[validation.column]],levels=cell.levels)

	validation.stats = do.call(rbind,parallel::mclapply(seq(0,1,0.01),function(i) {
		x = subset(these.predictions,glue_type_confidence >= i)
		n = nrow(x)
		n.match = sum(x[[validation.column]] == x$glue_type)
		data.frame(threshold = i,n.cells = n, perc.match = n.match/n)
	}))

	dir.create('stats/glue_predictions',showWarnings=FALSE)
	write.table(validation.stats,file=file.path('stats/glue_predictions',paste0(query.prefix,'-prediction-stats.txt')),sep='\t',row.names=FALSE,quote=FALSE)

	library(ggplot2)
	library(egg)

	p1 = ggplot(validation.stats,aes(threshold,log10(1-perc.match))) +
		geom_line() +
		geom_vline(xintercept=glue.confidence.threshold,color='red',size=0.5) +
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

	pdf(file='/dev/null')
	p = ggarrange(p2,p1,nrow=2)
	dev.off()

	pdf(file=file.path('figures',paste0(query.prefix,'-glue-prediction-stats.pdf')),useDingbats=FALSE,height=4)
	p
	dev.off()
} else {
	predictions.pass = subset(these.predictions,glue_type_confidence >= glue.confidence.threshold)
	predictions.split = split(predictions.pass,predictions.pass$glue_type)
	dir.create('stats/subintegration',showWarnings=FALSE)
	dev.null = lapply(seq_len(length(cell.levels)),function(i) {
		x = predictions.split[[i]]
		write.table(x['glue_type'],sep='\t',row.names=TRUE,col.names=FALSE,quote=FALSE,file=file.path('stats/subintegration',paste0(query.prefix,'-class',i,'-cells.txt')))
	})
	# Now do the RNA too
	rna.pass = read.delim(file.path('stats/clusters',paste0(rn.prefix,'-final-cellclasses.txt')),header=FALSE)
	rna.pass$V2 = factor(rna.pass$V2,levels=cell.levels)
	rna.split = split(rna.pass,rna.pass$V2)
	dev.null = lapply(seq_len(length(cell.levels)),function(i) {
		x = rna.split[[i]]
		write.table(x,sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE,file=file.path('stats/subintegration',paste0(rn.prefix,'-class',i,'-cells.txt')))
	})
}
