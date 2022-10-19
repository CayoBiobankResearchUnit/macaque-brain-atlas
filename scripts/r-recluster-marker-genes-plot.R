#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('rna','rna')

prefix = arguments[1]
analysis = arguments[2]

suppressMessages(library(monocle3))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(egg))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gridExtra))
suppressMessages(library(parallel))

# Cell class levels

dir.create('rds/markers',showWarnings=FALSE)
cell.classes = scan(file.path('stats/clusters',paste0(prefix,'-final-cellclasses-levels.txt')),what='',sep='\n')

marker.genes = parallel::mclapply(seq_along(cell.classes),function(i) {
#	readRDS(file.path('rds',paste0(prefix,'_class',i,'_marker_genes.rds')))
#	readRDS(file.path('rds/resolution1e-5',paste0(prefix,'_class',i,'_marker_genes.rds')))
	readRDS(file.path('rds/resolution1e-6',paste0(prefix,'_class',i,'_marker_genes.rds')))
},mc.cores=n.cores)

names(marker.genes) = cell.classes

cell.types = unlist(mclapply(cell.classes,function(i) {
	x = marker.genes[[i]]
	cell.subtypes = unique(x$cell)
	if (length(cell.subtypes) > 1) {
		cell.subtype.integers = as.integer(gsub(paste0(i,' '),'',cell.subtypes))
		paste(i,sort(cell.subtype.integers))
	} else {
		i
	}
},mc.cores=n.cores))

marker.genes = do.call(rbind,marker.genes)
rownames(marker.genes) = NULL

marker.genes = within(marker.genes,{
	cell = factor(cell,levels=cell.types)
	logfoldchanges = as.numeric(logfoldchanges)
	wilcox_score = as.numeric(wilcox_score)
	ttest_score = as.numeric(ttest_score)
	logreg_score = as.numeric(logreg_score)
	bp1 = as.integer(bp1)
	bp2 = as.integer(bp2)
	gene_strand = factor(gene_strand,levels=c('+','-'))
})

# Classify pts

marker.genes$pts_group = factor(with(marker.genes,ifelse(
	pts >= 0.5,
	'>50%',
	ifelse(
		pts >= 0.25,
		'>25%',
		ifelse(
			pts >= 0.10,
			'>10%',
			'<10%'
		)
	)
)),levels=c('<10%','>10%','>25%','>50%'))


# Gene names for homologs were already crunched in round 1. Just grab this file

add.gene.names = read.delim(file.path('biomart',paste0(prefix,'-orthologs-with-symbols.txt')),header=FALSE)
names(add.gene.names) = c('ensembl_gene_id','hsapiens_homolog_associated_gene_name')
rownames(add.gene.names) = add.gene.names$ensembl_gene_id
# library(biomaRt)
# mmul = useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL', dataset='mmulatta_gene_ensembl',version=101)
# mmul.orthologs = getBM(attributes=c('ensembl_gene_id','external_gene_name','hsapiens_homolog_ensembl_gene','hsapiens_homolog_orthology_type','hsapiens_homolog_associated_gene_name'),mart=mmul)
# 
# add.gene.names = subset(mmul.orthologs,hsapiens_homolog_orthology_type=='ortholog_one2one' & nchar(hsapiens_homolog_associated_gene_name) & !grepl('\\.[0-9]+$',hsapiens_homolog_associated_gene_name) & !nchar(external_gene_name) & ensembl_gene_id %in% marker.genes$gene)
# rownames(add.gene.names) = add.gene.names$ensembl_gene_id
# 
# dir.create('biomart',showWarnings=FALSE)
# write.table(add.gene.names[,c('ensembl_gene_id','hsapiens_homolog_associated_gene_name')],file=file.path('biomart',paste0(prefix,'-orthologs-with-symbols.txt')),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)

for (i in 1:nrow(add.gene.names)) {
	this.gene = add.gene.names[i,'ensembl_gene_id']
	this.symbol = add.gene.names[this.gene,'hsapiens_homolog_associated_gene_name']
	message(paste0('Now replacing ID ',this.gene,' with gene name ',this.symbol,'*.'))
	marker.genes[marker.genes$gene == this.gene,]$gene_short_name = paste0(this.symbol,'*')
}

n.show = 20
expr.frac = 0.1
top.expr.markers = do.call(rbind,lapply(split(marker.genes,marker.genes$cell),function(x) {
	y = subset(x,!is.nan(logreg_score) & pts >= expr.frac)
	y = y[order(y$logreg_score,decreasing=TRUE),]
	out = as.data.frame(y[,c('cell','gene_short_name','logfoldchanges','pts','pts_group','logreg_score')]) %>% head(n.show)
	if (nrow(out)) {
		out$i = 1:min(nrow(out),n.show)
	} else {
		out$i = integer(0)
	}
	out
}))

top.markers = do.call(rbind,lapply(split(marker.genes,marker.genes$cell),function(x) {
	y = subset(x,!is.nan(logreg_score) & TRUE)
	y = y[order(y$logreg_score,decreasing=TRUE),]
	out = as.data.frame(y[,c('cell','gene_short_name','logfoldchanges','pts','pts_group','logreg_score')]) %>% head(n.show)
	if (nrow(out)) {
		out$i = 1:min(nrow(out),n.show)
	} else {
		out$i = integer(0)
	}
	out
}))

saveRDS(top.markers,file=file.path('rds/markers',paste0(prefix,'_subclusters_','logreg','_marker_genes_','all','.rds')))
saveRDS(top.expr.markers,file=file.path('rds/markers',paste0(prefix,'_subclusters_','logreg','_marker_genes_','expr','.rds')))

# Multi page (need gridExtra)

top.expr.markers$parent_cell = with(top.expr.markers,factor(gsub(' [0-9]*$','',cell),levels=unique(gsub(' [0-9]*$','',levels(cell)))))[,drop=TRUE]
top.markers$parent_cell = with(top.markers,factor(gsub(' [0-9]*$','',cell),levels=unique(gsub(' [0-9]*$','',levels(cell)))))[,drop=TRUE]

pl = parallel::mclapply(levels(top.expr.markers$parent_cell),function(this.i) {
	this = subset(top.expr.markers,parent_cell == this.i)
	ncolumns = floor(sqrt(nlevels(this$cell[,drop=TRUE])))
	bsize = ceiling(40/max(ncolumns,3))
	p = ggplot(this) +
		geom_text(aes(i,logreg_score,color=pts_group,label=gene_short_name),hjust='left',vjust='center',angle=90,size=1,show.legend=FALSE) +
		geom_point(aes(i,logreg_score,color=pts_group),alpha=0) +
		facet_wrap(~cell,scales='free_y',ncol=ncolumns) +
		scale_color_manual(values=c('#cccccc','#1b9e77','#d95f02','#7570b3'),name='expr',drop=FALSE) +
		scale_y_continuous(expand = expansion(mult=c(0.1,0.5))) +
		theme_article(base_size=bsize) +
		theme(
			axis.text.x=element_blank(),
			axis.title.x=element_blank(),
			axis.ticks.x=element_blank()
		) +
		guides(color=guide_legend(override.aes=list(shape=15,size=4,alpha=1))) +
		ylab('Score (logistic regression)') +
		ggtitle(paste0('Top markers (expressed in >',round(expr.frac*100,2),'% of cells)'))
		p
},mc.cores=n.cores)
pdf(file='/dev/null'); glist = lapply(pl,ggplotGrob); suppressMessages(dev.off())
ggsave(file.path('figures',paste0(prefix,'-recluster-multipage-marker-genes-logreg-top.pdf')), marrangeGrob(glist, nrow=1, ncol=1, top=NULL))

pl = parallel::mclapply(levels(top.markers$parent_cell),function(this.i) {
	this = subset(top.markers,parent_cell == this.i)
	ncolumns = floor(sqrt(nlevels(this$cell[,drop=TRUE])))
	bsize = ceiling(40/max(ncolumns,3))
	p = ggplot(this) +
		geom_text(aes(i,logreg_score,color=pts_group,label=gene_short_name),hjust='left',vjust='center',angle=90,size=1,show.legend=FALSE) +
		geom_point(aes(i,logreg_score,color=pts_group),alpha=0) +
		facet_wrap(~cell,scales='free_y',ncol=ncolumns) +
		scale_color_manual(values=c('#cccccc','#1b9e77','#d95f02','#7570b3'),name='expr',drop=FALSE) +
		scale_y_continuous(expand = expansion(mult=c(0.1,0.5))) +
		theme_article(base_size=bsize) +
		theme(
			axis.text.x=element_blank(),
			axis.title.x=element_blank(),
			axis.ticks.x=element_blank()
		) +
		guides(color=guide_legend(override.aes=list(shape=15,size=4,alpha=1))) +
		ylab('Score (logistic regression)') +
		ggtitle(paste0('Top markers (expressed in >',round(expr.frac*100,2),'% of cells)'))
		p
},mc.cores=n.cores)
pdf(file='/dev/null'); glist = lapply(pl,ggplotGrob); suppressMessages(dev.off())
ggsave(file.path('figures',paste0(prefix,'-recluster-multipage-marker-genes-logreg-all.pdf')), marrangeGrob(glist, nrow=1, ncol=1, top=NULL))

p = ggplot(top.expr.markers) +
	geom_text(aes(i,logreg_score,color=pts_group,label=gene_short_name),hjust='left',vjust='center',angle=90,size=1,show.legend=FALSE) +
	geom_point(aes(i,logreg_score,color=pts_group),alpha=0) +
	facet_wrap(~cell,scales='free_y',ncol=floor(sqrt(nlevels(top.expr.markers$cell)))) +
	scale_color_manual(values=c('#cccccc','#1b9e77','#d95f02','#7570b3'),name='expr',drop=FALSE) +
	scale_y_continuous(expand = expansion(mult=c(0.1,0.5))) +
	theme_article(base_size=8) +
	theme(
		axis.text.x=element_blank(),
		axis.title.x=element_blank(),
		axis.ticks.x=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(shape=15,size=4,alpha=1))) +
	ylab('Score (logistic regression)') +
	ggtitle(paste0('Top markers (expressed in >',round(expr.frac*100,2),'% of cells)'))
ggsave(p,file=file.path('figures',paste0(prefix,'-recluster-marker-genes-logreg-top.pdf')),width=10,height=10)

p = ggplot(top.markers) +
	geom_text(aes(i,logreg_score,color=pts_group,label=gene_short_name),hjust='left',vjust='center',angle=90,size=1,show.legend=FALSE) +
	geom_point(aes(i,logreg_score,color=pts_group),alpha=0) +
	facet_wrap(~cell,scales='free_y',ncol=floor(sqrt(nlevels(top.expr.markers$cell)))) +
	scale_color_manual(values=c('#cccccc','#1b9e77','#d95f02','#7570b3'),name='expr',drop=FALSE) +
	scale_y_continuous(expand = expansion(mult=c(0.1,0.5))) +
	theme_article(base_size=8) +
	theme(
		axis.text.x=element_blank(),
		axis.title.x=element_blank(),
		axis.ticks.x=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(shape=15,size=4,alpha=1))) +
	ylab('Score (logistic regression)') +
	ggtitle(paste0('Top markers (all)'))
ggsave(p,file=file.path('figures',paste0(prefix,'-recluster-marker-genes-logreg-all.pdf')),width=10,height=10)

# Wilcoxon
top.expr.markers = do.call(rbind,lapply(split(marker.genes,marker.genes$cell),function(x) {
	y = subset(x,!is.nan(wilcox_score) & pts >= expr.frac)
	y = y[order(y$wilcox_score,decreasing=TRUE),]
	out = as.data.frame(y[,c('cell','gene_short_name','logfoldchanges','pts','pts_group','wilcox_score')]) %>% head(n.show)
	if (nrow(out)) {
		out$i = 1:min(nrow(out),n.show)
	} else {
		out$i = integer(0)
	}
	out
}))

top.markers = do.call(rbind,lapply(split(marker.genes,marker.genes$cell),function(x) {
	y = subset(x,!is.nan(wilcox_score) & TRUE)
	y = y[order(y$wilcox_score,decreasing=TRUE),]
	out = as.data.frame(y[,c('cell','gene_short_name','logfoldchanges','pts','pts_group','wilcox_score')]) %>% head(n.show)
	if (nrow(out)) {
		out$i = 1:min(nrow(out),n.show)
	} else {
		out$i = integer(0)
	}
	out
}))

saveRDS(top.markers,file=file.path('rds/markers',paste0(prefix,'_subclusters_','wilcox','_marker_genes_','all','.rds')))
saveRDS(top.expr.markers,file=file.path('rds/markers',paste0(prefix,'_subclusters_','wilcox','_marker_genes_','expr','.rds')))


top.expr.markers$parent_cell = with(top.expr.markers,factor(gsub(' [0-9]*$','',cell),levels=unique(gsub(' [0-9]*$','',levels(cell)))))[,drop=TRUE]
top.markers$parent_cell = with(top.markers,factor(gsub(' [0-9]*$','',cell),levels=unique(gsub(' [0-9]*$','',levels(cell)))))[,drop=TRUE]

pl = parallel::mclapply(levels(top.expr.markers$parent_cell),function(this.i) {
	this = subset(top.expr.markers,parent_cell == this.i)
	ncolumns = floor(sqrt(nlevels(this$cell[,drop=TRUE])))
	bsize = ceiling(40/max(ncolumns,3))
	p = ggplot(this) +
		geom_text(aes(i,wilcox_score,color=pts_group,label=gene_short_name),hjust='left',vjust='center',angle=90,size=1,show.legend=FALSE) +
		geom_point(aes(i,wilcox_score,color=pts_group),alpha=0) +
		facet_wrap(~cell,scales='free_y',ncol=ncolumns) +
		scale_color_manual(values=c('#cccccc','#1b9e77','#d95f02','#7570b3'),name='expr',drop=FALSE) +
		scale_y_continuous(expand = expansion(mult=c(0.1,0.5))) +
		theme_article(base_size=bsize) +
		theme(
			axis.text.x=element_blank(),
			axis.title.x=element_blank(),
			axis.ticks.x=element_blank()
		) +
		guides(color=guide_legend(override.aes=list(shape=15,size=4,alpha=1))) +
		ylab('Score (logistic regression)') +
		ggtitle(paste0('Top markers (expressed in >',round(expr.frac*100,2),'% of cells)'))
		p
},mc.cores=n.cores)
pdf(file='/dev/null'); glist = lapply(pl,ggplotGrob); suppressMessages(dev.off())
ggsave(file.path('figures',paste0(prefix,'-recluster-multipage-marker-genes-wilcox-top.pdf')), marrangeGrob(glist, nrow=1, ncol=1, top=NULL))

pl = parallel::mclapply(levels(top.markers$parent_cell),function(this.i) {
	this = subset(top.markers,parent_cell == this.i)
	ncolumns = floor(sqrt(nlevels(this$cell[,drop=TRUE])))
	bsize = ceiling(40/max(ncolumns,3))
	p = ggplot(this) +
		geom_text(aes(i,wilcox_score,color=pts_group,label=gene_short_name),hjust='left',vjust='center',angle=90,size=1,show.legend=FALSE) +
		geom_point(aes(i,wilcox_score,color=pts_group),alpha=0) +
		facet_wrap(~cell,scales='free_y',ncol=ncolumns) +
		scale_color_manual(values=c('#cccccc','#1b9e77','#d95f02','#7570b3'),name='expr',drop=FALSE) +
		scale_y_continuous(expand = expansion(mult=c(0.1,0.5))) +
		theme_article(base_size=bsize) +
		theme(
			axis.text.x=element_blank(),
			axis.title.x=element_blank(),
			axis.ticks.x=element_blank()
		) +
		guides(color=guide_legend(override.aes=list(shape=15,size=4,alpha=1))) +
		ylab('Score (logistic regression)') +
		ggtitle(paste0('Top markers (expressed in >',round(expr.frac*100,2),'% of cells)'))
		p
},mc.cores=n.cores)
pdf(file='/dev/null'); glist = lapply(pl,ggplotGrob); suppressMessages(dev.off())
ggsave(file.path('figures',paste0(prefix,'-recluster-multipage-marker-genes-wilcox-all.pdf')), marrangeGrob(glist, nrow=1, ncol=1, top=NULL))

p = ggplot(top.expr.markers) +
	geom_text(aes(i,wilcox_score,color=pts_group,label=gene_short_name),hjust='left',vjust='center',angle=90,size=1,show.legend=FALSE) +
	geom_point(aes(i,wilcox_score,color=pts_group),alpha=0) +
	facet_wrap(~cell,scales='free_y',ncol=floor(sqrt(nlevels(top.expr.markers$cell)))) +
	scale_color_manual(values=c('#cccccc','#1b9e77','#d95f02','#7570b3'),name='expr',drop=FALSE) +
	scale_y_continuous(expand = expansion(mult=c(0.1,0.5))) +
	theme_article(base_size=8) +
	theme(
		axis.text.x=element_blank(),
		axis.title.x=element_blank(),
		axis.ticks.x=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(shape=15,size=4,alpha=1))) +
	ylab('Score (Wilcoxon test)') +
	ggtitle(paste0('Top markers (expressed in >',round(expr.frac*100,2),'% of cells)'))
ggsave(p,file=file.path('figures',paste0(prefix,'-recluster-marker-genes-wilcox-top.pdf')),width=10,height=10)

p = ggplot(top.markers) +
	geom_text(aes(i,wilcox_score,color=pts_group,label=gene_short_name),hjust='left',vjust='center',angle=90,size=1,show.legend=FALSE) +
	geom_point(aes(i,wilcox_score,color=pts_group),alpha=0) +
	facet_wrap(~cell,scales='free_y',ncol=floor(sqrt(nlevels(top.expr.markers$cell)))) +
	scale_color_manual(values=c('#cccccc','#1b9e77','#d95f02','#7570b3'),name='expr',drop=FALSE) +
	scale_y_continuous(expand = expansion(mult=c(0.1,0.5))) +
	theme_article(base_size=8) +
	theme(
		axis.text.x=element_blank(),
		axis.title.x=element_blank(),
		axis.ticks.x=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(shape=15,size=4,alpha=1))) +
	ylab('Score (Wilcoxon test)') +
	ggtitle(paste0('Top markers (all)'))
ggsave(p,file=file.path('figures',paste0(prefix,'-recluster-marker-genes-wilcox-all.pdf')),width=10,height=10)

top.expr.markers = do.call(rbind,lapply(split(marker.genes,marker.genes$cell),function(x) {
	y = subset(x,!is.nan(ttest_score) & pts >= expr.frac)
	y = y[order(y$ttest_score,decreasing=TRUE),]
	out = as.data.frame(y[,c('cell','gene_short_name','logfoldchanges','pts','pts_group','ttest_score')]) %>% head(n.show)
	if (nrow(out)) {
		out$i = 1:min(nrow(out),n.show)
	} else {
		out$i = integer(0)
	}
	out
}))

top.markers = do.call(rbind,lapply(split(marker.genes,marker.genes$cell),function(x) {
	y = subset(x,!is.nan(ttest_score) & TRUE)
	y = y[order(y$ttest_score,decreasing=TRUE),]
	out = as.data.frame(y[,c('cell','gene_short_name','logfoldchanges','pts','pts_group','ttest_score')]) %>% head(n.show)
	if (nrow(out)) {
		out$i = 1:min(nrow(out),n.show)
	} else {
		out$i = integer(0)
	}
	out
}))

saveRDS(top.markers,file=file.path('rds/markers',paste0(prefix,'_subclusters_','ttest','_marker_genes_','all','.rds')))
saveRDS(top.expr.markers,file=file.path('rds/markers',paste0(prefix,'_subclusters_','ttest','_marker_genes_','expr','.rds')))

top.expr.markers$parent_cell = with(top.expr.markers,factor(gsub(' [0-9]*$','',cell),levels=unique(gsub(' [0-9]*$','',levels(cell)))))[,drop=TRUE]
top.markers$parent_cell = with(top.markers,factor(gsub(' [0-9]*$','',cell),levels=unique(gsub(' [0-9]*$','',levels(cell)))))[,drop=TRUE]

pl = parallel::mclapply(levels(top.expr.markers$parent_cell),function(this.i) {
	this = subset(top.expr.markers,parent_cell == this.i)
	ncolumns = floor(sqrt(nlevels(this$cell[,drop=TRUE])))
	bsize = ceiling(40/max(ncolumns,3))
	p = ggplot(this) +
		geom_text(aes(i,ttest_score,color=pts_group,label=gene_short_name),hjust='left',vjust='center',angle=90,size=1,show.legend=FALSE) +
		geom_point(aes(i,ttest_score,color=pts_group),alpha=0) +
		facet_wrap(~cell,scales='free_y',ncol=ncolumns) +
		scale_color_manual(values=c('#cccccc','#1b9e77','#d95f02','#7570b3'),name='expr',drop=FALSE) +
		scale_y_continuous(expand = expansion(mult=c(0.1,0.5))) +
		theme_article(base_size=bsize) +
		theme(
			axis.text.x=element_blank(),
			axis.title.x=element_blank(),
			axis.ticks.x=element_blank()
		) +
		guides(color=guide_legend(override.aes=list(shape=15,size=4,alpha=1))) +
		ylab('Score (logistic regression)') +
		ggtitle(paste0('Top markers (expressed in >',round(expr.frac*100,2),'% of cells)'))
		p
},mc.cores=n.cores)
pdf(file='/dev/null'); glist = lapply(pl,ggplotGrob); suppressMessages(dev.off())
ggsave(file.path('figures',paste0(prefix,'-recluster-multipage-marker-genes-ttest-top.pdf')), marrangeGrob(glist, nrow=1, ncol=1, top=NULL))

pl = parallel::mclapply(levels(top.markers$parent_cell),function(this.i) {
	this = subset(top.markers,parent_cell == this.i)
	ncolumns = floor(sqrt(nlevels(this$cell[,drop=TRUE])))
	bsize = ceiling(40/max(ncolumns,3))
	p = ggplot(this) +
		geom_text(aes(i,ttest_score,color=pts_group,label=gene_short_name),hjust='left',vjust='center',angle=90,size=1,show.legend=FALSE) +
		geom_point(aes(i,ttest_score,color=pts_group),alpha=0) +
		facet_wrap(~cell,scales='free_y',ncol=ncolumns) +
		scale_color_manual(values=c('#cccccc','#1b9e77','#d95f02','#7570b3'),name='expr',drop=FALSE) +
		scale_y_continuous(expand = expansion(mult=c(0.1,0.5))) +
		theme_article(base_size=bsize) +
		theme(
			axis.text.x=element_blank(),
			axis.title.x=element_blank(),
			axis.ticks.x=element_blank()
		) +
		guides(color=guide_legend(override.aes=list(shape=15,size=4,alpha=1))) +
		ylab('Score (logistic regression)') +
		ggtitle(paste0('Top markers (expressed in >',round(expr.frac*100,2),'% of cells)'))
		p
},mc.cores=n.cores)
pdf(file='/dev/null'); glist = lapply(pl,ggplotGrob); suppressMessages(dev.off())
ggsave(file.path('figures',paste0(prefix,'-recluster-multipage-marker-genes-ttest-all.pdf')), marrangeGrob(glist, nrow=1, ncol=1, top=NULL))

p = ggplot(top.expr.markers) +
	geom_text(aes(i,ttest_score,color=pts_group,label=gene_short_name),hjust='left',vjust='center',angle=90,size=1,show.legend=FALSE) +
	geom_point(aes(i,ttest_score,color=pts_group),alpha=0) +
	facet_wrap(~cell,scales='free_y',ncol=floor(sqrt(nlevels(top.expr.markers$cell)))) +
	scale_color_manual(values=c('#cccccc','#1b9e77','#d95f02','#7570b3'),name='expr',drop=FALSE) +
	scale_y_continuous(expand = expansion(mult=c(0.1,0.5))) +
	theme_article(base_size=8) +
	theme(
		axis.text.x=element_blank(),
		axis.title.x=element_blank(),
		axis.ticks.x=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(shape=15,size=4,alpha=1))) +
	ylab('Score (t-test)') +
	ggtitle(paste0('Top markers (expressed in >',round(expr.frac*100,2),'% of cells)'))
ggsave(p,file=file.path('figures',paste0(prefix,'-recluster-marker-genes-ttest-top.pdf')),width=10,height=10)

p = ggplot(top.markers) +
	geom_text(aes(i,ttest_score,color=pts_group,label=gene_short_name),hjust='left',vjust='center',angle=90,size=1,show.legend=FALSE) +
	geom_point(aes(i,ttest_score,color=pts_group),alpha=0) +
	facet_wrap(~cell,scales='free_y',ncol=floor(sqrt(nlevels(top.expr.markers$cell)))) +
	scale_color_manual(values=c('#cccccc','#1b9e77','#d95f02','#7570b3'),name='expr',drop=FALSE) +
	scale_y_continuous(expand = expansion(mult=c(0.1,0.5))) +
	theme_article(base_size=8) +
	theme(
		axis.text.x=element_blank(),
		axis.title.x=element_blank(),
		axis.ticks.x=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(shape=15,size=4,alpha=1))) +
	ylab('Score (t-test)') +
	ggtitle(paste0('Top markers (all)'))
ggsave(p,file=file.path('figures',paste0(prefix,'-recluster-marker-genes-ttest-all.pdf')),width=10,height=10)


# Automate naming of cell classes
expr.frac = 0.1
top.expr.markers = do.call(rbind,lapply(split(marker.genes,marker.genes$cell),function(x) {
	y = subset(x,!is.nan(logreg_score) & pts >= expr.frac)
	y = y[order(y$logreg_score,decreasing=TRUE),]
	y$parent_cell = gsub(' [0-9]+$','',y$cell)
	out = as.data.frame(y[,c('cell','parent_cell','gene_short_name','logfoldchanges','pts','pts_group','logreg_score')])
	if (nrow(out)) {
		out$i = 1:nrow(out)
	} else {
		out$i = integer(0)
	}
	out
}))

cells.out = data.frame(within(unique(top.expr.markers[c('cell','parent_cell')]),{cell_subcluster=cell})[c('cell_subcluster','parent_cell')],cell_subcluster_label='',complete=FALSE)
rownames(cells.out) = NULL
for (parent.cell in unique(top.expr.markers$parent_cell)) {
	this.cell = subset(top.expr.markers,parent_cell == parent.cell)
	j = 1
	while(nrow(this.cell)) {
		next.genes = subset(this.cell,i==j)
		for (k in next.genes$cell) {
			cells.out[cells.out$cell_subcluster==k,]$cell_subcluster_label = ifelse(nchar(cells.out[cells.out$cell_subcluster==k,]$cell_subcluster_label),paste(cells.out[cells.out$cell_subcluster==k,]$cell_subcluster_label,subset(next.genes,cell == k)$gene_short_name,sep='-'),subset(next.genes,cell == k)$gene_short_name)
		}
		good.labels = names(which(table(subset(cells.out,parent_cell == parent.cell)$cell_subcluster_label) == 1))
		cells.out[cells.out$parent_cell == parent.cell & cells.out$cell_subcluster_label %in% good.labels,]$complete = TRUE
		
		this.cell = subset(this.cell,!cell %in% subset(cells.out,complete)$cell_subcluster)
		j = j + 1
	}
}

cells.out$cell_subcluster_label = with(cells.out,paste(cell_subcluster_label,parent_cell))
write.table(cells.out[c('cell_subcluster','parent_cell','cell_subcluster_label')],file=file.path('stats',paste0(prefix,'-recluster-auto-names-logreg.txt')),sep='\t',row.names=FALSE,quote=FALSE)



# Automate naming of cell classes
expr.frac = 0.1
top.expr.markers = do.call(rbind,lapply(split(marker.genes,marker.genes$cell),function(x) {
	y = subset(x,!is.nan(wilcox_score) & pts >= expr.frac)
	y = y[order(y$wilcox_score,decreasing=TRUE),]
	y$parent_cell = gsub(' [0-9]+$','',y$cell)
	out = as.data.frame(y[,c('cell','parent_cell','gene_short_name','logfoldchanges','pts','pts_group','wilcox_score')])
	if (nrow(out)) {
		out$i = 1:nrow(out)
	} else {
		out$i = integer(0)
	}
	out
}))

cells.out = data.frame(within(unique(top.expr.markers[c('cell','parent_cell')]),{cell_subcluster=cell})[c('cell_subcluster','parent_cell')],cell_subcluster_label='',complete=FALSE)
rownames(cells.out) = NULL
for (parent.cell in unique(top.expr.markers$parent_cell)) {
	this.cell = subset(top.expr.markers,parent_cell == parent.cell)
	j = 1
	while(nrow(this.cell)) {
		next.genes = subset(this.cell,i==j)
		for (k in next.genes$cell) {
			cells.out[cells.out$cell_subcluster==k,]$cell_subcluster_label = ifelse(nchar(cells.out[cells.out$cell_subcluster==k,]$cell_subcluster_label),paste(cells.out[cells.out$cell_subcluster==k,]$cell_subcluster_label,subset(next.genes,cell == k)$gene_short_name,sep='-'),subset(next.genes,cell == k)$gene_short_name)
		}
		good.labels = names(which(table(subset(cells.out,parent_cell == parent.cell)$cell_subcluster_label) == 1))
		cells.out[cells.out$parent_cell == parent.cell & cells.out$cell_subcluster_label %in% good.labels,]$complete = TRUE
		
		this.cell = subset(this.cell,!cell %in% subset(cells.out,complete)$cell_subcluster)
		j = j + 1
	}
}

cells.out$cell_subcluster_label = with(cells.out,paste(cell_subcluster_label,parent_cell))
write.table(cells.out[c('cell_subcluster','parent_cell','cell_subcluster_label')],file=file.path('stats',paste0(prefix,'-recluster-auto-names-wilcox.txt')),sep='\t',row.names=FALSE,quote=FALSE)
