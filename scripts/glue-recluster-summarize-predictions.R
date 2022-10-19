#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('biccn','rna','atac','9','2')

prefix = arguments[1]
rn.prefix = arguments[2]
at.prefix = arguments[3]
this.cluster = as.integer(arguments[4])
do.matrix = as.integer(arguments[5])

glue.confidence.threshold = 0.75

if (do.matrix == 1) {
	query.prefix = at.prefix
	validation.column = NULL
} else if (do.matrix == 2) {
	query.prefix = rn.prefix
	validation.column = 'rnasubtype'
} else {
	stop('Error')
}

these.predictions = read.delim(file.path('stats/clusters/subintegration',paste0(prefix,'-',query.prefix,'-class',this.cluster,'-celltype-predictions.txt')),row.names=1)

these.predictions = these.predictions[,c('id','animal_id','n.umi','region','modality',validation.column,'glue_subtype','glue_subtype_confidence')]

cell.levels = scan(file=file.path('stats/subclusters',paste0(rn.prefix,'-cellsubclusters-levels.txt')),what='',sep='\n')

these.predictions$region = factor(these.predictions$region,levels=region.levels)
these.predictions$glue_subtype = factor(these.predictions$glue_subtype,levels=cell.levels)[,drop=TRUE]

if (do.matrix > 1) {
	these.predictions[[validation.column]] = factor(these.predictions[[validation.column]],levels=cell.levels)[,drop=TRUE]
	
	if (!identical(levels(these.predictions$glue_subtype),levels(these.predictions[[validation.column]]))) {
		message('Prediction levels missing')
		levels(these.predictions$glue_subtype) = levels(these.predictions[[validation.column]])
	}
	validation.stats = do.call(rbind,parallel::mclapply(seq(0,1,0.01),function(i) {
		x = subset(these.predictions,glue_subtype_confidence >= i)
		n = nrow(x)
		n.match = sum(x[[validation.column]] == x$glue_subtype)
		data.frame(threshold = i,n.cells = n, perc.match = n.match/n)
	}))

	dir.create('stats/glue_predictions/subintegration',showWarnings=FALSE)
	write.table(validation.stats,file=file.path('stats/glue_predictions/subintegration',paste0(prefix,'-',query.prefix,'-class',this.cluster,'-prediction-stats.txt')),sep='\t',row.names=FALSE,quote=FALSE)

	library(ggplot2)
	library(egg)

	p1 = ggplot(validation.stats,aes(threshold,log10((1-perc.match)+0.01))) +
		geom_line() +
		geom_vline(xintercept=glue.confidence.threshold,color='red',size=0.5) +
		theme_classic(base_size=16) +
		xlab('confidence score threshold') +
		ylab(expression('log'[10]~'error + 0.01'))

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

	dir.create('figures/glue/subintegration',showWarnings=FALSE)

	pdf(file=file.path('figures/glue/subintegration',paste0(prefix,'-',query.prefix,'-class',this.cluster,'-glue-prediction-stats.pdf')),useDingbats=FALSE,height=4)
	p
	dev.off()
} else {
	predictions.pass = subset(these.predictions,glue_subtype_confidence >= glue.confidence.threshold)
	
	write.table(these.predictions,file=file.path('stats/clusters/subintegration',paste0(prefix,'-',query.prefix,'-class',this.cluster,'-cellsubtype-predictions.txt')),sep='\t',row.names=TRUE,col.names=NA,quote=FALSE)
	write.table(predictions.pass,file=file.path('stats/clusters/subintegration',paste0(prefix,'-',query.prefix,'-class',this.cluster,'-cellsubtypes.txt')),sep='\t',row.names=TRUE,col.names=NA,quote=FALSE)
}
