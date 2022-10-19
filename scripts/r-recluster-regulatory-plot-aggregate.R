#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('subpeak','atacsub','ENSMMUG00000008778',0)

# PVALB:  ENSMMUG00000008778
# SST:    ENSMMUG00000000850
# ADARB2: ENSMMUG00000017120
prefix = arguments[1]
atac.prefix = arguments[2]
this.gene = arguments[3]

if (length(arguments) > 3) {
	show.background = as.logical(as.integer(arguments[4]))
} else {
	show.background = TRUE
}

suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrastr))
suppressMessages(library(ggforce))
suppressMessages(library(lemon))
suppressMessages(library(egg))
# suppressMessages(library(gtable))
# suppressMessages(library(grid))
# suppressMessages(library(gridExtra))
# suppressMessages(library(ggplotify))

file.colnames = c('peak','ensembl_gene_id','dist','gene_dir','glue.score','glue.pval','glue.qval','meta.beta','meta.pval','meta.qval')

cell.classes = scan(what='',sep='\n',file='stats/clusters/rna-final-cellclasses-levels.txt',quiet=TRUE)

genome.padding = 0.05

all.results = lapply(1:length(cell.classes),function(this.cluster) {
	this.cell.class = cell.classes[this.cluster]
	inc.file = file.path('stats/regulatory',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_peaks_inc.txt'))
	dec.file = file.path('stats/regulatory',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_peaks_dec.txt'))
	background.file = file.path('stats/regulatory/bed',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_peaks_background.bed'))
	if (file.exists(inc.file) && file.exists(dec.file)) {
		results = rbind(
			read.delim(inc.file,header=FALSE,col.names=file.colnames),
			read.delim(dec.file,header=FALSE,col.names=file.colnames)
		)
		results$cell_class = this.cell.class
		results
	} else if (file.exists(inc.file)) {
		results = read.delim(inc.file,header=FALSE,col.names=file.colnames)
		results$cell_class = this.cell.class
		results
	} else if (file.exists(dec.file)) {
		results = read.delim(dec.file,header=FALSE,col.names=file.colnames)
		results$cell_class = this.cell.class
		results
	} else {
		NULL
	}
})

all.results = do.call(rbind,all.results[!unlist(lapply(all.results,is.null))])

# Get gene annotations
if (!file.exists(file.path('checkpoints',paste0(ens.species,'_genes.rds')))) {
	suppressMessages(library(biomaRt))

	ens = biomaRt::useEnsembl(biomart='ENSEMBL_MART_ENSEMBL',dataset=paste0(ens.species,'_gene_ensembl'))

	ens.genes = biomaRt::getBM(
				attributes=c('ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','strand'),
				mart = ens)

	ens.exons = biomaRt::getBM(
				attributes=c('ensembl_gene_id','ensembl_transcript_id','ensembl_exon_id','exon_chrom_start','exon_chrom_end'),
				mart = ens)

	ens.tss = biomaRt::getBM(
				attributes=c('ensembl_gene_id','ensembl_transcript_id','transcript_start','transcript_end','transcription_start_site','transcript_biotype'),
				mart = ens)

	ens.utr = biomaRt::getBM(
				attributes=c('ensembl_gene_id','ensembl_transcript_id','5_utr_start','3_utr_start'),
				mart = ens)

	saveRDS(ens.genes,file=file.path('checkpoints',paste0(ens.species,'_genes.rds')))
	saveRDS(ens.exons,file=file.path('checkpoints',paste0(ens.species,'_exons.rds')))
	saveRDS(ens.tss,file=file.path('checkpoints',paste0(ens.species,'_transcripts.rds')))
	saveRDS(ens.utr,file=file.path('checkpoints',paste0(ens.species,'_utr.rds')))
} else {
	ens.genes = readRDS(file.path('checkpoints',paste0(ens.species,'_genes.rds')))
	ens.exons = readRDS(file.path('checkpoints',paste0(ens.species,'_exons.rds')))
	ens.tss = readRDS(file.path('checkpoints',paste0(ens.species,'_transcripts.rds')))
	ens.utr = readRDS(file.path('checkpoints',paste0(ens.species,'_utr.rds')))
}

background.bed = do.call(rbind,lapply(match(unique(all.results$cell_class),cell.classes),function(this.cluster) {
	this.cell.class = cell.classes[this.cluster]
	background.file = file.path('stats/regulatory/bed',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_peaks_genes_background.txt'))
	background.bed = read.delim(background.file,header=FALSE,col.names=c('peak','ensembl_gene_id'))
	background.bed$cell_class = this.cell.class
	background.bed
}))

all.results.bed = data.frame(
	as.data.frame.matrix(do.call(rbind,strsplit(all.results$peak,'_'))),
	all.results
)
names(all.results.bed)[1:3] = c('chr','start','stop')
all.results.bed$start = as.numeric(all.results.bed$start)
all.results.bed$stop = as.numeric(all.results.bed$stop)

all.results.bed = merge(all.results.bed,ens.genes,by='ensembl_gene_id',all.x=TRUE,all.y=FALSE)

results.this = all.results.bed[all.results.bed$ensembl_gene_id == this.gene,]
background.this = background.bed[background.bed$ensembl_gene_id == this.gene,]

if (any(results.this$meta.pval==0)) {
	warning('p values of 0 found, setting them to minimum nonzero p value instead')
	results.this[results.this$meta.pval==0,'meta.pval'] = min(all.results.bed$meta.pval[all.results.bed$meta.pval > 0])
}

background.this = data.frame(as.data.frame.matrix(do.call(rbind,strsplit(background.this$peak,'_'))),background.this)
names(background.this)[1:3] = c('chr','start','stop')

background.this = within(background.this,{
	glue.score = NA
	glue.pval = NA
	meta.beta = NA
	meta.pval = NA
	start = as.integer(start)
	stop = as.integer(stop)
})

if (show.background) {
	peaks.this = rbind(
		results.this[,c('chr','start','stop','glue.score','glue.pval','meta.beta','meta.pval','cell_class')],
		subset(background.this,!peak %in% results.this$peak,select=c('chr','start','stop','glue.score','glue.pval','meta.beta','meta.pval','cell_class'))
	)
	background.prefix=''
} else {
	peaks.this = results.this[,c('chr','start','stop','glue.score','glue.pval','meta.beta','meta.pval','cell_class')]
	background.prefix='_no_background'
}

genome.max = max(c(peaks.this$start,peaks.this$stop,results.this$start_position,results.this$end_position))
genome.min = min(c(peaks.this$start,peaks.this$stop,results.this$start_position,results.this$end_position))

genome.range = genome.max - genome.min

genome.chr = unique(results.this$chr)
genome.max = genome.max + genome.padding * genome.range
genome.min = genome.min - genome.padding * genome.range

gene.this = ens.genes[ens.genes$ensembl_gene_id == this.gene,]
tss.this = ens.tss[ens.tss$ensembl_gene_id == this.gene,]
exons.this = ens.exons[ens.exons$ensembl_gene_id == this.gene,]

tss.this$transcription_end_site = with(tss.this,ifelse(transcription_start_site == transcript_start, transcript_end, ifelse(transcription_start_site == transcript_end,transcript_start,NA)))

tss.this = tss.this[order(tss.this$transcript_start,tss.this$transcript_end),]

global.tss = if (gene.this$strand > 0) {
	min(tss.this$transcription_start_site)
} else {
	max(tss.this$transcription_start_site)
}

results.this$tss = global.tss
results.this$midpoint = with(results.this,start+(stop-start)/2)

results.this$direction = factor(results.this$meta.beta > 0,levels=c('TRUE','FALSE'),labels=c('upregulate','downregulate'))

links.this = do.call(rbind,lapply(split(results.this,with(results.this,paste0(peak,cell_class))),function(x) {
	this.peak = unique(x$peak)
	this.cell.class = unique(x$cell_class)
	this.direction = unique(x$direction)
	this.start = min(c(x$midpoint,x$tss))
	this.end = max(c(x$midpoint,x$tss))
	data.frame(
		x = c(this.start,mean(c(this.start,this.end)),this.end),
		y = c(0,2*(-log10(x$meta.pval)),0),
		peak = this.peak,
		direction = this.direction,
		cell_class = this.cell.class
	)
}))

# exons.this = ens.exons[ens.exons$ensembl_gene_id == this.gene,]

peaks.this$cell_class = factor(peaks.this$cell_class,levels=cell.classes)[,drop=TRUE]
links.this$cell_class = factor(links.this$cell_class,levels=cell.classes)[,drop=TRUE]

scaling.factor = 20/max(c(10,links.this$y)) * (nlevels(links.this$cell_class) * 0.67)

tss.this$y = if (nrow(tss.this) == 1) {
	-4/scaling.factor
} else {
	-5/scaling.factor + seq(0,2/scaling.factor,2/scaling.factor/(nrow(tss.this)-1))
}

exons.this = merge(tss.this,exons.this,by=c('ensembl_gene_id','ensembl_transcript_id'))

max.transcript.length = max(abs(tss.this$transcript_end - tss.this$transcript_start))
arrow.interval = round(max.transcript.length/10)
# Split tss.this into chunks
tss.arrow = do.call(rbind,lapply(1:nrow(tss.this),function(i) {
	x = tss.this[i,]
	meta = x[c('ensembl_gene_id','ensembl_transcript_id','y')]
	if (x$transcription_start_site > x$transcription_end_site) {
		new.start = seq(x$transcription_start_site,x$transcription_end_site,-arrow.interval)
		new.end = new.start - arrow.interval
		new.end[new.end < x$transcription_end_site] = x$transcription_end_site
	} else {
		new.start = seq(x$transcription_start_site,x$transcription_end_site,arrow.interval)
		new.end = new.start + arrow.interval
		new.end[new.end > x$transcription_end_site] = x$transcription_end_site
	}
	suppressWarnings(unique(data.frame(meta,transcription_start_site = new.start, transcription_end_site = new.end)))
}))

open.color        = '#ffffff'
open.line.color   = '#000000'
closed.color      = '#cccccc'
closed.line.color = '#000000'
link.inc.color    = '#377eb8'
link.dec.color    = '#e41a1c'
exon.color        = '#000000'
exon.line.color   = NA
gene.color        = '#000000'

peak.height = max(c(10,links.this$y/1.8))/(25 * nlevels(links.this$cell_class)) 

peak.start = -0.05 * peak.height
peak.end = -0.95 * peak.height

plot.breaks = seq(0,max(c(10,links.this$y/1.8)),10^round(log10(max(c(10,links.this$y/1.8))/4)))

if (length(plot.breaks) > 5) {
	plot.breaks = plot.breaks[rep(c(TRUE,FALSE),ceiling(length(plot.breaks)/2))[1:length(plot.breaks)]]
}

p = ggplot() +
	geom_rect(data=subset(peaks.this,!is.na(meta.pval)),aes(xmin=start,xmax=stop,ymin=peak.end/scaling.factor,ymax=peak.start/scaling.factor),fill=open.color,color=open.line.color,size=0.25)

if (any(is.na(peaks.this$meta.pval))) {
	p = p + geom_rect(data=subset(peaks.this,is.na(meta.pval)),aes(xmin=start,xmax=stop,ymin=peak.end/scaling.factor,ymax=peak.start/scaling.factor),fill=closed.color,color=closed.line.color,size=0.25)
}

p = p +
	geom_bezier(data=links.this,aes(x=x,y=y,group=peak,color=direction)) +
	facet_rep_grid(cell_class~1,scales='fixed') +
	coord_flex_fixed(ratio = genome.range/50*scaling.factor,ylim=c(peak.end/scaling.factor,max(c(10,links.this$y/1.8))),left=capped_vertical(capped='bottom')) +
	scale_x_continuous(limits=c(genome.min,genome.max)) +
	scale_y_continuous(name=expression(-log[10]~italic('p')),c(peak.end/scaling.factor,max(c(10,links.this$y/1.8))),breaks=plot.breaks) +
	scale_color_manual(values=c(link.inc.color,link.dec.color)) +
	theme_classic() +
	theme(
		legend.position='none',
		strip.background=element_blank(),
		axis.title.x = element_blank(),
		axis.title.y = element_text(),
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		strip.text.x = element_blank()
	) +
	ylab(expression(-log[10]~italic('p'))) +
	ggtitle(ifelse(gene.this$ensembl_gene_id == gene.this$external_gene_name,gene.this$ensembl_gene_id,eval(parse(text=paste0('expression(italic(',gene.this$external_gene_name,'))')))))

p2 = ggplot() +
	geom_segment(data=tss.arrow,aes(x=transcription_start_site,xend=transcription_end_site,y=y,yend=y),size=0.25,arrow=arrow(length=unit(0.05,'inches'))) +
	geom_rect(data=exons.this,aes(xmin=exon_chrom_start,xmax=exon_chrom_end,ymin=y-0.25/scaling.factor,ymax=y+0.25/scaling.factor),fill=exon.color,color=exon.line.color,size=0.25) +
	geom_segment(data=unique(results.this[c('ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','strand')]),aes(x=start_position,xend=end_position,y=-1.5/scaling.factor,yend=-1.5/scaling.factor),color=gene.color,size=0.5,arrow=arrow(length=unit(0.05,'inches'),ends='both',angle=90)) +
	coord_flex_fixed(ratio = genome.range/50*scaling.factor,ylim=c(min(tss.this$y) - 2/scaling.factor,1/scaling.factor),left=capped_vertical(capped='bottom')) +
	scale_y_continuous(c(min(tss.this$y) - 2/scaling.factor,1/scaling.factor)) +
	scale_x_continuous(limits=c(genome.min,genome.max)) +
	theme_classic() +
	theme(
		legend.position='none',
		axis.text.y = element_blank(),
		axis.title.y = element_blank(),
		axis.ticks.y = element_blank(),
		axis.line.y = element_line()
	) +
	xlab(paste('Chromosome',gene.this$chromosome_name))

dir.create('figures/regulatory/tracks',showWarnings=FALSE)
# ggsave(p,file=file.path('figures/regulatory/tracks',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_peaks_tracks_',this.gene,'_',gene.this$external_gene_name,'.pdf')),useDingbats=FALSE)

# # ggplotGrob draws to device. Disable this by sending the output to /dev/null (not ideal but I give up on better solutions)
# pdf(file='/dev/null')
# 	g = ggplotGrob(p)
# dev.off()
# 
# g$grobs[[as.integer(rownames(g$layout)[g$layout$name == 'ylab-r'])]] = textGrob(expression(-log[10]~italic(p)), rot = 90, vjust = 0)
# 
# p1 = as.ggplot(g) 

pdf(file='/dev/null')
p3 = ggarrange(p,p2,ncol=1,heights = c(nlevels(links.this$cell_class)*2,0.5))
dev.off()

pdf(file=file.path('figures/regulatory/tracks',paste0(prefix,'_',atac.prefix,background.prefix,'_all_peaks_tracks_',this.gene,'_',gene.this$external_gene_name,'.pdf')),useDingbats=FALSE,height=7)
p3
dev.off()

