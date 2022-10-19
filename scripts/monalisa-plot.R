#!/usr/bin/env Rscript

library(monaLisa)
library(ggplot2)


TFpfms=readRDS('rds-noah/JASPAR2018_pfms.rds')
TF_log_enrich = readRDS('rds-noah/TF_cell_class_log2enr.rds') # log2 OR 
TF_log_padj =readRDS('rds-noah/TF_cell_class_negLog10Padj.rds')#  - log10 padj

jaspar.info = read.delim('data/jaspar2018_CORE_vertebrates_non-redundant.txt',header=FALSE,col.names=c('motif','motif_name'))



colnames(TF_log_enrich) = colnames(TF_log_padj) = c(
'excitatory neurons',
'oligodendrocyte precursor cells',
'vascular cells',
'microglia',
'ependymal cells',
'medium spiny neurons',
'inhibitory neurons',
'cerebellar neurons',
'basket cells',
'astrocytes',
'oligodendrocytes'
)

TF_log_enrich.melt = reshape2::melt(TF_log_enrich)
TF_log_padj.melt = reshape2::melt(TF_log_padj)

tf.stats = TF_log_enrich.melt
names(tf.stats) = c('motif_name','cell_class','log2OR')

tf.stats = data.frame(tf.stats,log10padj = TF_log_padj.melt$value)

tf.stats = merge(tf.stats,jaspar.info,by='motif_name',all.x=TRUE)

tf.motifs.pass = unique(subset(tf.stats,log10padj > -log10(0.05) & log2OR > 0)$motif_name)

tf.pass = droplevels(subset(tf.stats,motif_name %in% tf.motifs.pass))

tf.best = do.call(rbind,lapply(split(tf.pass,tf.pass$cell_class),function(x) head(x[order(x$log10padj,decreasing=TRUE),],10)))

tf.plot = character(5 * nlevels(tf.best$cell_class))

cells.sorted = sort(tapply(tf.best$log2OR, tf.best$cell_class,mean),decreasing=TRUE)
for (i in 1:length(cells.sorted)) {
	this.cell = names(cells.sorted[i])
	this.tf = as.character(head(subset(tf.best,cell_class == this.cell & !motif_name %in% tf.plot)$motif_name,5))
	tf.plot[((i-1)*5 + 1):(((i-1)*5 + 1)+4)] = this.tf
}

tf.keep = subset(tf.stats,motif_name %in% tf.plot)

tf.keep$cell_class = factor(tf.keep$cell_class,levels=names(cells.sorted))
tf.keep$motif_name = factor(tf.keep$motif_name,levels=tf.plot)

tf.keep = tf.keep[order(tf.keep$cell_class,tf.keep$motif_name),]

tf.ids = unique(tf.keep[order(tf.keep$motif_name),]$motif)

p = ggplot(tf.keep,aes(motif_name,cell_class,fill=log2OR)) +
	geom_tile() +
	coord_fixed() +
	# scale_fill_viridis() +
	scale_fill_gradient2(
		high='#2066ac',
		mid='#ffffff',
		low='#b2182b',
		midpoint=0,
		breaks=seq(-1,1,1),
		# limits=c(-max(abs(tf.keep$log2OR)),max(abs(tf.keep$log2OR))),
		limits=c(-1.5,1.5),
		oob=scales::squish,
		name=expression(log[2]~'enrichment')
	) +
	scale_y_discrete(limits=rev) +
	theme_classic(base_size=8) +
	theme(
		axis.title=element_blank(),
		axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5)
	)
ggsave(p,file='figures/final/tf_jaspar_unique_peaks_enrichment.pdf',width=7,height=4,useDingbats=FALSE)


## select PFMs, where TFs is an matrix of the TFs you include
## TFs should be JASPAR names (e.g., MA0623.1)

# TFS = TF_log_enrich[as.character(unique(tf.keep$motif_name)),]
# rownames(TFS) = as.character(unique(tf.keep$motif))
# pfms <- TFpfms[TFS]

pfms = TFpfms

# tf.keep = subset(tf.stats,motif_name %in% c('SPI1', 'ELF1', 'NFATC2', 'NEUROD2', 'FLI1', 'IRF5', 'MAFF', 'ETV5', 'EBF1', 'BHLHE22'))

tf.to.logo = do.call(c,lapply(split(tf.keep,tf.keep$cell_class),function(x) head(x[order(x$log10padj,decreasing=TRUE),],1)$motif))

for (i in unique(tf.keep$motif)) {
# for (i in tf.to.logo) {
	print(i)
	this.pfm = pfms[i]
	this.motif.name = unique(as.character(subset(tf.keep,motif == i)$motif_name))
	this.grob = list(seqLogoGrob(pfms[[i]]))
	this.logo = ComplexHeatmap::HeatmapAnnotation(
		logo = annoSeqlogo(
			grobL = this.grob,
			which = 'row',
			space = unit(0.5, "mm"),
			width = unit(1.5, "inch")
		),
		show_legend=FALSE,show_annotation_name = FALSE,which='row')
	
	pdf(file=paste0('figures/motif/pfm_',this.motif.name,'.pdf'),useDingbats=FALSE,width=7,height=1)
		plot(this.logo)
	dev.off()
}



tf.sig = subset(tf.stats,log10padj > -log10(0.05) & log2OR > 0)

tf.sig$cell_class = factor(tf.sig$cell_class,levels=names(tf.to.logo))

tf.sig = tf.sig[order(tf.sig$cell_class,-tf.sig$log10padj,-tf.sig$log2OR),]

tf.sig$padj = 10^-tf.sig$log10padj
tf.sig$odds_ratio = 2^tf.sig$log2OR

tf.sig.out = tf.sig[,c('cell_class','motif_name','motif','odds_ratio','padj')]
rownames(tf.sig.out) = NULL

names(tf.sig.out) = c('cell_class','motif_name','jaspar_id','odds_ratio','p_adj')

dir.create('tables',showWarnings=FALSE)
dir.create('tables/final',showWarnings=FALSE)

write.table(tf.sig.out,file='tables/final/monalisa_enrichment.txt',quote=FALSE,row.names=FALSE,sep='\t')


# ## identify the longest one
# maxwidth <- max(vapply(TFBSTools::Matrix(pfms), ncol, 0L))
# 
# ## create list of logo grob for TFs from the PFMs above
# grobL <- lapply(pfms, seqLogoGrob, xmax = maxwidth, xjust = "center")
# 
# ## create object with PWM logos
# hmSeqlogo <- ComplexHeatmap::HeatmapAnnotation(
#     logo = annoSeqlogo(grobL = grobL, which = "row",
#     space = unit(0.5, "mm"),
#     width = unit(1.5, "inch")),
#     show_legend = FALSE, show_annotation_name = FALSE, which = "row")
# 
# ## plot heatmap (add hmSeqLogo)
# pdf(file='monalisa.pdf')
# ComplexHeatmap::Heatmap(
#   matrix = TF_log_enrich, ## this is the TF x cell matrix
#   name = "logFC",
#   width = unit(3,"inch"),
#   column_title = "x", show_heatmap_legend = TRUE,
#   heatmap_legend_param = list(color_bar = "continuous"),
#   use_raster = F, show_row_dend = F,show_column_dend = T,
#   left_annotation = hmSeqlogo,row_names_side = "left")
# dev.off()