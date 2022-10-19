#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

library(parallel)
library(GenomicRanges)

if (file.exists('rds/peaks_summary.rds')) {

all.peaks.df = readRDS('rds/peaks_summary.rds')

} else {

main.file = 'data/merged_peaks.bed'
class.files = file.path('bed',list.files('bed',pattern='all_peaks_class[0-9]+.bed$'))

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

rownames(ens.genes) = ens.genes$ensembl_gene_id

# Find the TSS
ens.genes$tss = with(ens.genes,ifelse(strand > 0,start_position,end_position))

# Calculate promoter as TSS extended 2000 bp
ens.genes$promoter_start = with(ens.genes,ifelse(strand > 0,tss-2000,tss))
ens.genes$promoter_end = with(ens.genes,ifelse(strand > 0,tss,tss+2000))

## A very slow way to do the same as above
# ens.genes2 = do.call(rbind,mclapply(rownames(ens.genes),function(i) {
# 	gene.this = ens.genes[i,]
# 	tss.this = ens.tss[ens.tss$ensembl_gene_id == i,]
# 	gene.this$tss = if (gene.this$strand > 0) {
# 		min(tss.this$transcription_start_site)
# 	} else {
# 		max(tss.this$transcription_start_site)
# 	}
# 	gene.this
# },mc.cores = 8))

i.genes = with(ens.genes,GRanges(chromosome_name,IRanges(start_position,end_position,names=ensembl_gene_id),'*'))
i.promoters = with(ens.genes,GRanges(chromosome_name,IRanges(promoter_start,promoter_end,names=ensembl_gene_id),'*'))

all.peaks = read.delim(main.file,header=FALSE)
class.peaks = lapply(class.files,read.delim,header=FALSE)

peaks = c(list(all.peaks),class.peaks)
names(peaks) = c(main.file,class.files)

all.peaks = lapply(peaks,function(x) {
	x$V4 = with(x,paste(V1,V2,V3,sep='_'))
	x$i = 1:nrow(x)
	i.this = with(x,GRanges(V1,IRanges(V2,V3,names=V4),'*'))
	gene.overlap = findOverlaps(i.this,i.genes)
	this.gene.overlaps = data.frame(
		peak = x$V4[queryHits(findOverlaps(i.this,i.genes))],
		gene = ens.genes$ensembl_gene_id[subjectHits(findOverlaps(i.this,i.genes))]
	)
	x$gene_overlap = x$V4 %in% this.gene.overlaps$peak
	
	which.genes = with(this.gene.overlaps,tapply(gene,peak,paste,collapse=','))
	x$which_genes = which.genes[x$V4]
	
	x$distance_to_gene = distanceToNearest(i.this,i.genes)@elementMetadata$distance
	
	promoter.overlap = findOverlaps(i.this,i.promoters)
	this.promoter.overlaps = data.frame(
		peak = x$V4[queryHits(findOverlaps(i.this,i.promoters))],
		gene = ens.genes$ensembl_gene_id[subjectHits(findOverlaps(i.this,i.promoters))]
	)
	x$promoter_overlap = x$V4 %in% this.promoter.overlaps$peak
	
	which.promoters = with(this.promoter.overlaps,tapply(gene,peak,paste,collapse=','))
	x$which_promoters = which.promoters[x$V4]
	
	x$distance_to_promoter = distanceToNearest(i.this,i.promoters)@elementMetadata$distance
	
	names(x)[1:4] = c('chr','peak_start','peak_end','peak')
	x
})

cell.classes = scan(what='',sep='\n',file='stats/clusters/rna-final-cellclasses-levels.txt',quiet=TRUE)

all.peaks.df = do.call(rbind,lapply(1:length(all.peaks),function(i) {
	bed.file = names(all.peaks)[i]
	x = all.peaks[[i]]
	cell.type = if (bed.file == 'data/merged_peaks.bed') {
		0
	} else {
		as.integer(gsub('bed/all_peaks_class([0-9]+).bed','\\1',bed.file))
	}
	x$cell_class_integer = cell.type
	x$cell_class = factor(ifelse(cell.type == 0,'all cells',cell.classes[cell.type]),levels=c('all cells',cell.classes))
	x
}))

all.peaks.df$cell_class = all.peaks.df$cell_class[,drop=TRUE]
all.peaks.df$peak_length = with(all.peaks.df,peak_end-peak_start)

all.peaks.df$chr = factor(all.peaks.df$chr,levels=unique(all.peaks.df$chr))

saveRDS(all.peaks.df,'rds/peaks_summary.rds')

}

fai = read.delim('data/Macaca_mulatta.Mmul_10.dna.toplevel.fa.fai',header=FALSE)
fai = fai[fai$V1 %in% levels(all.peaks.df$chr),]
fai$V1 = factor(fai$V1,levels=fai$V1)

genome.size = sum(fai$V2)

with(all.peaks.df,tapply(peak_length,cell_class,length))
with(all.peaks.df,tapply(peak_length,cell_class,sum))

p = ggplot(all.peaks.df,aes(log10(peak_length))) +
	geom_histogram() +
	facet_wrap(~cell_class,scales='free_y') +
	theme_classic()

p = ggplot(all.peaks.df,aes(log10(distance_to_gene + 1))) +
	geom_histogram() +
	facet_wrap(~cell_class,scales='free_y') +
	scale_y_continuous(trans='log10') +
	theme_classic()

p = ggplot(all.peaks.df,aes(log10(distance_to_promoter + 1))) +
	geom_histogram() +
	facet_wrap(~cell_class,scales='free_y') +
	scale_y_continuous(trans='log10') +
	theme_classic()

peaks.summary = do.call(rbind,lapply(split(all.peaks.df,all.peaks.df$cell_class),function(x) {
	data.frame(
		cell_class = unique(x$cell_class),
		n_peaks = nrow(x),
		peak_length = sum(x$peak_length),
		peak_length_pct = sum(x$peak_length)/genome.size,
		pct_overlap_genes = sum(x$gene_overlap) / nrow(x),
		pct_overlap_promoters = sum(x$promoter_overlap) / nrow(x),
		pct_overlap_genes_or_promoters = sum(x$gene_overlap | x$promoter_overlap) / nrow(x),
		pct_overlap_genes_and_promoters = sum(x$gene_overlap & x$promoter_overlap) / nrow(x),
		pct_overlap_genes_not_promoters = sum(x$gene_overlap & !x$promoter_overlap) / nrow(x),
		pct_overlap_promoters_not_genes = sum(!x$gene_overlap & x$promoter_overlap) / nrow(x),
		pct_no_overlap = sum(!(x$gene_overlap | x$promoter_overlap)) / nrow(x),
		n_overlap_genes = sum(x$gene_overlap),
		n_overlap_promoters = sum(x$promoter_overlap),
		n_overlap_genes_or_promoters = sum(x$gene_overlap | x$promoter_overlap),
		n_no_overlap = sum(!(x$gene_overlap | x$promoter_overlap))
	)
}))

peaks.summary.long = reshape2::melt(peaks.summary[,c('cell_class','n_peaks','peak_length_pct','pct_overlap_genes_not_promoters','pct_overlap_genes_and_promoters','pct_overlap_promoters_not_genes','pct_no_overlap')])

colors = as.matrix(read.table('data/colors.txt',row.names=1,sep='\t',comment.char='',header=FALSE))[,1]
cell.levels = c('excitatory neurons','cerebellar neurons','inhibitory neurons','basket cells','medium spiny neurons','dopaminergic neurons','serotonergic neurons','AHSG neurons','F5 neurons','KIR3DL12 neurons','astrocytes','oligodendrocytes','oligodendrocyte precursor cells','radial glial cells','ependymal cells','microglia','KIR3DL12 microglia','vascular cells','mesenchymal stem cells')

colors = c(colors,c('all cells' = '#000000'))
fill.colors = colors
fill.colors['all cells'] = '#ffffff'

peaks.summary.long$cell_class = factor(peaks.summary.long$cell_class,levels=c('all cells',cell.levels))[,drop=TRUE]

peaks.summary.long = droplevels(subset(peaks.summary.long,!cell_class %in% c('radial glial cells')))

peaks.summary.long$variable_facet = factor(peaks.summary.long$variable,levels=c('n_peaks','peak_length_pct','pct_overlap_genes_not_promoters','pct_overlap_genes_and_promoters','pct_overlap_promoters_not_genes','pct_no_overlap'),labels=c('number of peaks','fraction of genome','percent overlapping','percent overlapping','percent overlapping','percent overlapping'))

peaks.summary.long$variable = factor(
	peaks.summary.long$variable,
	levels=c(
		'n_peaks',
		'peak_length_pct',
		'pct_no_overlap',
		'pct_overlap_promoters_not_genes',
		'pct_overlap_genes_and_promoters',
		'pct_overlap_genes_not_promoters'
	)
)

alpha.values = c(
'n_peaks' = 1,
'peak_length_pct' = 1,
'pct_no_overlap' = 0.25,
'pct_overlap_promoters_not_genes' = 0.5,
'pct_overlap_genes_and_promoters' = 0.75,
'pct_overlap_genes_not_promoters' = 1
)

peaks.summary.long$variable2 = factor(
	peaks.summary.long$variable,
	levels=c(
		'n_peaks',
		'peak_length_pct',
		'pct_no_overlap',
		'pct_overlap_promoters_not_genes',
		'pct_overlap_genes_and_promoters',
		'pct_overlap_genes_not_promoters'
	),
	labels=c(
		'n_peaks',
		'peak_length_pct',
		'pct_no_overlap',
		'pct_overlap',
		'pct_overlap',
		'pct_overlap'
	)
)

alpha.values2 = c(
'n_peaks' = 1,
'peak_length_pct' = 1,
'pct_no_overlap' = 0.25,
'pct_overlap' = 1
)


p = ggplot(peaks.summary.long,aes(cell_class,value,fill=cell_class,alpha=variable)) +
	geom_bar(stat='identity',position='stack') +
	facet_wrap(~variable_facet,scales='free_x',nrow=1) +
	scale_fill_manual(values=colors[levels(peaks.summary.long$cell_class)]) +
	scale_alpha_manual(values=alpha.values[levels(peaks.summary.long$variable)]) +
	scale_x_discrete(limits=rev) +
	theme_classic() +
	theme(strip.background=element_blank(),axis.title = element_blank(),legend.position='none',axis.text.x=element_text(vjust=1,hjust=0,angle=-30)) +
	coord_flip()

peaks.summary.long = rbind(
	peaks.summary.long,
	within(subset(peaks.summary.long,variable=='peak_length_pct'),{value=1-value; variable2='no_peak_pct'})
)

peaks.summary.long$variable2 = factor(
	peaks.summary.long$variable2,
	levels=c(
		'n_peaks',
		'no_peak_pct',
		'peak_length_pct',
		'pct_no_overlap',
		'pct_overlap'
	)
)

alpha.values2 = c(
'n_peaks' = 1,
'no_peak_pct' = 0.25,
'peak_length_pct' = 1,
'pct_no_overlap' = 0.25,
'pct_overlap' = 1
)

peaks.summary.long$variable_facet = factor(
	peaks.summary.long$variable_facet,
	levels=c("number of peaks","fraction of genome","percent overlapping"),
	labels=c('number of peaks (thousands)','percent genome in peaks','percent peaks in promoters/genes')
)

peaks.summary.long$value[peaks.summary.long$variable2 == 'n_peaks'] = peaks.summary.long$value[peaks.summary.long$variable2 == 'n_peaks'] / 1000
peaks.summary.long$value[peaks.summary.long$variable2 %in% c('no_peak_pct','peak_length_pct','pct_no_overlap','pct_overlap')] = peaks.summary.long$value[peaks.summary.long$variable2 %in% c('no_peak_pct','peak_length_pct','pct_no_overlap','pct_overlap')] * 100

p2 = ggplot(peaks.summary.long,aes(cell_class,value,fill=cell_class,alpha=variable2)) +
	geom_bar(stat='identity',position='stack') +
	facet_wrap(~variable_facet,scales='free_x',nrow=1) +
	scale_fill_manual(values=colors[levels(peaks.summary.long$cell_class)]) +
	scale_alpha_manual(values=alpha.values2[levels(peaks.summary.long$variable2)]) +
	scale_x_discrete(limits=rev) +
	theme_classic(base_size=16) +
	theme(strip.background=element_blank(),axis.title = element_blank(),legend.position='none',axis.text.x=element_text(vjust=1,hjust=0,angle=-30)) +
	coord_flip()

ggsave(p2,file='figures/final/atac_summary.pdf',width=7,height=4,useDingbats=FALSE)
