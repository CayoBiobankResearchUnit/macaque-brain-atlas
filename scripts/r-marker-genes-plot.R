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

marker.genes = readRDS(file.path('rds',paste0(prefix,'_marker_genes.rds')))

cell.types = c(
'excitatory neurons',
'medium spiny neurons',
'inhibitory neurons',
'dopaminergic neurons',
'serotinergic neurons',
'thalamic interneurons',
'cerebellar granule cells',
'Purkinje cells',
'basket cells',
'astrocytes',
'oligodendrocytes',
'oligodendrocyte precursor cells',
'endothelial cells',
'pericytes',
'microglia',
'radial glial cells',
'mesenchymal stem cells',
'ependymal cells',
'partition 14',
'partition 15',
'partition 16',
'partition 17'
)

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


library(biomaRt)
mmul = useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL', dataset='mmulatta_gene_ensembl',version=101)
mmul.orthologs = getBM(attributes=c('ensembl_gene_id','external_gene_name','hsapiens_homolog_ensembl_gene','hsapiens_homolog_orthology_type','hsapiens_homolog_associated_gene_name'),mart=mmul)

add.gene.names = subset(mmul.orthologs,hsapiens_homolog_orthology_type=='ortholog_one2one' & nchar(hsapiens_homolog_associated_gene_name) & !grepl('\\.[0-9]+$',hsapiens_homolog_associated_gene_name) & !nchar(external_gene_name) & ensembl_gene_id %in% marker.genes$gene)
rownames(add.gene.names) = add.gene.names$ensembl_gene_id

dir.create('biomart',showWarnings=FALSE)
write.table(add.gene.names[,c('ensembl_gene_id','hsapiens_homolog_associated_gene_name')],file=file.path('biomart',paste0(prefix,'-orthologs-with-symbols.txt')),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)

for (i in 1:nrow(add.gene.names)) {
	this.gene = add.gene.names[i,'ensembl_gene_id']
	this.symbol = add.gene.names[this.gene,'hsapiens_homolog_associated_gene_name']
	message(paste0('Now replacing ID ',this.gene,' with gene name ',this.symbol,'*.'))
	marker.genes[marker.genes$gene == this.gene,]$gene_short_name = paste0(this.symbol,'*')
}

n.show = 20
expr.frac = 0.1
top.expr.markers = do.call(rbind,lapply(split(marker.genes,marker.genes$cell),function(x) {
	y = subset(x,pts >= expr.frac)
	y = y[order(y$logreg_score,decreasing=TRUE),]
	out = as.data.frame(y[,c('cell','gene_short_name','logfoldchanges','pts','pts_group','logreg_score')]) %>% head(n.show)
	out$i = 1:n.show
	out
}))

top.markers = do.call(rbind,lapply(split(marker.genes,marker.genes$cell),function(x) {
	y = subset(x,TRUE)
	y = y[order(y$logreg_score,decreasing=TRUE),]
	out = as.data.frame(y[,c('cell','gene_short_name','logfoldchanges','pts','pts_group','logreg_score')]) %>% head(n.show)
	out$i = 1:n.show
	out
}))

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
ggsave(p,file=file.path('figures',paste0(prefix,'-marker-genes-logreg-top.pdf')))

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
ggsave(p,file=file.path('figures',paste0(prefix,'-marker-genes-logreg-all.pdf')))

top.expr.markers = do.call(rbind,lapply(split(marker.genes,marker.genes$cell),function(x) {
	y = subset(x,pts >= expr.frac)
	y = y[order(y$wilcox_score,decreasing=TRUE),]
	out = as.data.frame(y[,c('cell','gene_short_name','logfoldchanges','pts','pts_group','wilcox_score')]) %>% head(n.show)
	out$i = 1:n.show
	out
}))

top.markers = do.call(rbind,lapply(split(marker.genes,marker.genes$cell),function(x) {
	y = subset(x,TRUE)
	y = y[order(y$wilcox_score,decreasing=TRUE),]
	out = as.data.frame(y[,c('cell','gene_short_name','logfoldchanges','pts','pts_group','wilcox_score')]) %>% head(n.show)
	out$i = 1:n.show
	out
}))

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
ggsave(p,file=file.path('figures',paste0(prefix,'-marker-genes-wilcox-top.pdf')))

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
ggsave(p,file=file.path('figures',paste0(prefix,'-marker-genes-wilcox-all.pdf')))

top.expr.markers = do.call(rbind,lapply(split(marker.genes,marker.genes$cell),function(x) {
	y = subset(x,pts >= expr.frac)
	y = y[order(y$ttest_score,decreasing=TRUE),]
	out = as.data.frame(y[,c('cell','gene_short_name','logfoldchanges','pts','pts_group','ttest_score')]) %>% head(n.show)
	out$i = 1:n.show
	out
}))

top.markers = do.call(rbind,lapply(split(marker.genes,marker.genes$cell),function(x) {
	y = subset(x,TRUE)
	y = y[order(y$ttest_score,decreasing=TRUE),]
	out = as.data.frame(y[,c('cell','gene_short_name','logfoldchanges','pts','pts_group','ttest_score')]) %>% head(n.show)
	out$i = 1:n.show
	out
}))

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
ggsave(p,file=file.path('figures',paste0(prefix,'-marker-genes-ttest-top.pdf')))

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
ggsave(p,file=file.path('figures',paste0(prefix,'-marker-genes-ttest-all.pdf')))

