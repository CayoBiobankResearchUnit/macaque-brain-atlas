#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('subpeak','atacsub','7')

prefix = arguments[1]
atac.prefix = arguments[2]
this.cluster = as.integer(arguments[3])

suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrastr))

cell.classes = scan(what='',sep='\n',file='stats/clusters/rna-final-cellclasses-levels.txt',quiet=TRUE)
this.cell.class = cell.classes[this.cluster]

metacell.lr.results = read.table(
	file.path('stats/metacell_lr',paste0('gene_peak_class',this.cluster,'.tsv')),
	col.names = c('meta.beta','meta.serr','meta.sbet','meta.pval','peak','ensembl_gene_id'),
	sep='\t'
)

metacell.lr.null = read.table(
	file.path('stats/metacell_lr',paste0('gene_peak_null_class',this.cluster,'.tsv')),
	col.names = c('null.beta','null.serr','null.sbet','null.pval','peak','ensembl_gene_id'),
	sep='\t'
)

metacell.lr.results = merge(metacell.lr.results,metacell.lr.null,all.x=TRUE,all.y=TRUE,by=c('ensembl_gene_id','peak'))

glue.regulatory.results = fread(
	file.path('stats/glue/subintegration',paste0(prefix,'_class',this.cluster,'_regulatory.txt.gz'))
#	col.names = c('ensembl_gene_id','peak','glue.qval','dist','glue.pval','glue.score','glue.weight','glue.type')
)[,c('source','target','qval','dist','pval','score','weight','type')]
names(glue.regulatory.results) = c('ensembl_gene_id','peak','glue.qval','dist','glue.pval','glue.score','glue.weight','glue.type')

results.combined = merge(metacell.lr.results,glue.regulatory.results,by=c('ensembl_gene_id','peak'),all.y=TRUE,all.x=TRUE)

results.combined$meta.qval = p.adjust(results.combined$meta.pval,'fdr')

if (any(!complete.cases(results.combined))) {
	warning(sum(!complete.cases(results.combined)),' records are incomplete and will be removed.')
	warning(sum(is.na(results.combined$meta.pval)),' metacell p-values are NA.')
	warning(sum(is.na(results.combined$null.pval)),' permuted metacell p-values are NA.')
	results.combined = results.combined[complete.cases(results.combined),]
}

dir.create('figures/regulatory',showWarnings=FALSE)

# GLUE p values against metacell p values
p = ggplot(results.combined, aes(-log10(glue.pval),-log10(meta.pval))) +
	geom_point_rast(size=0.01,alpha=0.1) +
	geom_smooth(method=lm) +
#	coord_cartesian(xlim=c(0,5),ylim=c(0,5)) +
	theme_classic(base_size=16) +
	xlab(expression('-log'[10]~italic(p)[GLUE])) +
	ylab(expression('-log'[10]~italic(p)[metacell])) +
	ggtitle(this.cell.class)
suppressMessages(ggsave(p,file=file.path('figures/regulatory',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_glueP-v-metacellP.pdf')),useDingbats=FALSE))

# # GLUE score against abs(metacell standardized beta)
# p = ggplot(results.combined, aes(glue.score,abs(meta.sbet))) +
# 	geom_point_rast(size=0.01,alpha=0.1) +
# 	geom_smooth(method=lm) +
# 	theme_classic(base_size=16) +
# 	coord_cartesian(ylim=c(0,quantile(abs(results.combined$meta.sbet),0.999))) +
# 	xlab('GLUE score') + 
# 	ylab(expression(paste('|','std.'~beta[observed],'|')))

# GLUE score against metacell standardized beta
p = ggplot(results.combined, aes(glue.score,meta.sbet,color=(glue.qval < 0.05 & meta.qval<0.05))) +
	geom_point_rast(size=0.01,alpha=0.1) +
#	geom_smooth(method=lm) +
	scale_color_manual(values=c('black','red')) +
	coord_cartesian(ylim=c(-quantile(results.combined$meta.sbet,0.999),quantile(results.combined$meta.sbet,0.999))) +
	theme_classic(base_size=16) +
	theme(legend.position='none') +
	xlab('GLUE score') + 
	ylab(expression('std.'~beta[observed])) +
	ggtitle(this.cell.class)
suppressMessages(ggsave(p,file=file.path('figures/regulatory',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_glueScore-v-metacellStdBeta.pdf')),useDingbats=FALSE))

# GLUE score against metacell beta
p = ggplot(results.combined, aes(glue.score,meta.beta,color=(glue.qval < 0.05 & meta.qval<0.05))) +
	geom_point_rast(size=0.01,alpha=0.1) +
#	geom_smooth(method=lm) +
	scale_color_manual(values=c('black','red')) +
	coord_cartesian(ylim=c(-quantile(results.combined$meta.beta,0.999),quantile(results.combined$meta.beta,0.999))) +
	theme_classic(base_size=16) +
	theme(legend.position='none') +
	xlab('GLUE score') + 
	ylab(expression(beta[observed])) +
	ggtitle(this.cell.class)
suppressMessages(ggsave(p,file=file.path('figures/regulatory',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_glueScore-v-metacellBeta.pdf')),useDingbats=FALSE))

# GLUE p-values against GLUE score
p = ggplot(results.combined, aes(glue.score,-log10(glue.pval),color=glue.qval<0.05)) +
	geom_point_rast(size=0.01,alpha=0.1) +
	scale_color_manual(values=c('black','red')) +
	theme_classic(base_size=16) +
	theme(legend.position='none') +
	xlab('GLUE score') + 
	ylab(expression('-log'[10]~italic(p)[GLUE])) +
	ggtitle(this.cell.class)
suppressMessages(ggsave(p,file=file.path('figures/regulatory',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_glueScore-v-glueP.pdf')),useDingbats=FALSE))

# QQ plot of p values
p = ggplot(results.combined, aes(-log10(sort(null.pval)),-log10(sort(meta.pval)))) +
	geom_point_rast(size=0.01,alpha=1) +
	geom_abline(slope=1,intercept=0) +
#	coord_equal() +
	theme_classic(base_size=16) +
	xlab(expression('-log'[10]~italic(p)[permuted])) +
	ylab(expression('-log'[10]~italic(p)[observed])) +
	ggtitle(this.cell.class)
suppressMessages(ggsave(p,file=file.path('figures/regulatory',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_metacellP-qqplot.pdf')),useDingbats=FALSE))

# QQ plot of standardized beta
p = ggplot(results.combined, aes(sort(null.sbet),sort(meta.sbet))) +
	geom_point_rast(size=0.01,alpha=1) +
	geom_abline(slope=1,intercept=0) +
#	coord_equal() +
	theme_classic(base_size=16) +
	xlab(expression('std.'~italic(beta)[permuted])) +
	ylab(expression('std.'~italic(beta)[observed])) +
	ggtitle(this.cell.class)
suppressMessages(ggsave(p,file=file.path('figures/regulatory',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_metacellBeta-qqplot.pdf')),useDingbats=FALSE))

results.sig = subset(results.combined,glue.qval < 0.05 & meta.qval < 0.05)

results.sig$gene_dir = with(results.sig,ifelse(meta.beta>0,'+','-'))

gene.counts = as.data.frame.matrix(as.matrix(table(results.sig$ensembl_gene_id,results.sig$gene_dir)))
gene.counts$direction = ifelse(gene.counts$`-` > 0 & gene.counts$`+` > 0,'±',ifelse(gene.counts$`+` > 0,'+','-'))
gene.counts$ensembl_gene_id = rownames(gene.counts)
gene.counts$direction = factor(gene.counts$direction,levels=c('-','+','±'))

gene.counts.long = subset(reshape2::melt(gene.counts,id.vars=c('ensembl_gene_id','direction'),measure.vars=c('-','+')),value > 0)

gene.counts.long$variable_text = factor(
	gene.counts.long$variable,
	levels=levels(gene.counts.long$variable),
	labels=paste0(levels(gene.counts.long$variable),' regulators (',table(gene.counts.long$variable),' genes, ',tapply(gene.counts.long$value,gene.counts.long$variable,sum),' peaks)')
)

p = ggplot(gene.counts.long,aes(value,fill=direction)) +
	geom_histogram() +
	facet_wrap(~variable_text,nrow=2,scales='free_y') +
	scale_fill_manual(name='Gene direction',values=c('#91bfdb','#fc8d59','#999999')) +
	theme_classic(base_size=16) +
	xlab('Number of elements') +
	ylab('Count') +
	ggtitle(this.cell.class)
suppressMessages(ggsave(p,file=file.path('figures/regulatory',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_genes-histogram.pdf')),useDingbats=FALSE))

p = ggplot(gene.counts.long,aes(value,fill=direction)) +
	geom_histogram() +
	facet_wrap(~variable_text,nrow=2,scales='fixed') +
	scale_fill_manual(name='Gene direction',values=c('#91bfdb','#fc8d59','#999999')) +
	theme_classic(base_size=16) +
	xlab('Number of elements') +
	ylab('Count') +
	ggtitle(this.cell.class)
suppressMessages(ggsave(p,file=file.path('figures/regulatory',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_genes-histogram-fixed.pdf')),useDingbats=FALSE))

peak.counts = as.data.frame.matrix(as.matrix(table(results.sig$peak,results.sig$gene_dir)))
peak.counts$direction = ifelse(peak.counts$`-` > 0 & peak.counts$`+` > 0,'±',ifelse(peak.counts$`+` > 0,'+','-'))
peak.counts$peak = rownames(peak.counts)
peak.counts$direction = factor(peak.counts$direction,levels=c('-','+','±'))

peak.counts.long = subset(reshape2::melt(peak.counts,id.vars=c('peak','direction'),measure.vars=c('-','+')),value > 0)

peak.counts.long$variable_text = factor(
	peak.counts.long$variable,
	levels=levels(peak.counts.long$variable),
	labels=paste0(levels(gene.counts.long$variable),' regulators (',tapply(gene.counts.long$value,gene.counts.long$variable,sum),' peaks, ',table(gene.counts.long$variable),' genes)')
)

p = ggplot(peak.counts.long,aes(value,fill=direction)) +
	geom_histogram() +
	facet_wrap(~variable_text,nrow=2,scales='free_y') +
	scale_fill_manual(name='Gene direction',values=c('#91bfdb','#fc8d59','#999999')) +
	theme_classic(base_size=16) +
	xlab('Number of genes') +
	ylab('Count') +
	ggtitle(this.cell.class)
suppressMessages(ggsave(p,file=file.path('figures/regulatory',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_peaks-histogram.pdf')),useDingbats=FALSE))

p = ggplot(peak.counts.long,aes(value,fill=direction)) +
	geom_histogram() +
	facet_wrap(~variable_text,nrow=2,scales='fixed') +
	scale_fill_manual(name='Gene direction',values=c('#91bfdb','#fc8d59','#999999')) +
	theme_classic(base_size=16) +
	xlab('Number of genes') +
	ylab('Count') +
	ggtitle(this.cell.class)
suppressMessages(ggsave(p,file=file.path('figures/regulatory',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_peaks-histogram-fixed.pdf')),useDingbats=FALSE))


results.sig = results.sig[,c('peak','ensembl_gene_id','dist','gene_dir','glue.score','glue.pval','glue.qval','meta.beta','meta.pval','meta.qval')]

results.inc = subset(results.sig,gene_dir == '+')
results.dec = subset(results.sig,gene_dir == '-')

dir.create('stats/regulatory',showWarnings=FALSE)

if (nrow(results.inc)) {
	write.table(
		results.inc,
		file=file.path('stats/regulatory',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_peaks_inc.txt')),
		sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
}

if (nrow(results.dec)) {
	write.table(
		results.dec,
		file=file.path('stats/regulatory',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_peaks_dec.txt')),
		sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
}

results.sig.bed = unique(results.sig[,c('peak','gene_dir')])
results.background.bed = unique(results.combined[,c('peak')])

results.background.bed = as.data.frame.matrix(do.call(rbind,strsplit(results.background.bed,'_')))
results.background.bed = within(results.background.bed,{
	V1 = factor(V1,levels=c(1:20,'X','Y'))
	V2 = as.numeric(V2)
	V3 = as.numeric(V3)
})
results.background.bed = results.background.bed[with(results.background.bed,order(V1,V2,V3)),]

results.sig.bed = data.frame(
	as.data.frame.matrix(do.call(rbind,strsplit(results.sig.bed$peak,'_'))),
	results.sig.bed
)
results.sig.bed = within(results.sig.bed,{
	V1 = factor(V1,levels=c(1:20,'X','Y'))
	V2 = as.numeric(V2)
	V3 = as.numeric(V3)
})
results.sig.bed = results.sig.bed[with(results.sig.bed,order(V1,V2,V3)),]

rownames(results.sig.bed) = rownames(results.background.bed) = NULL

dir.create('stats/regulatory/bed',showWarnings=FALSE)

results.background.tests = unique(results.combined[,c('peak','ensembl_gene_id')])

write.table(
	results.background.tests,
	file=file.path('stats/regulatory/bed',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_peaks_genes_background.txt')),
	sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)

write.table(
	subset(results.sig.bed,gene_dir == '-',select=paste0('V',1:3)),
	file=file.path('stats/regulatory/bed',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_peaks_significant_dec.bed')),
	sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)

write.table(
	subset(results.sig.bed,gene_dir == '+',select=paste0('V',1:3)),
	file=file.path('stats/regulatory/bed',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_peaks_significant_inc.bed')),
	sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)

write.table(
	unique(subset(results.sig.bed,TRUE,select=paste0('V',1:3))),
	file=file.path('stats/regulatory/bed',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_peaks_significant.bed')),
	sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)

write.table(
	unique(subset(results.background.bed,TRUE,select=paste0('V',1:3))),
	file=file.path('stats/regulatory/bed',paste0(prefix,'_',atac.prefix,'_class',this.cluster,'_peaks_background.bed')),
	sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
