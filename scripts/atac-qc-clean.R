#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')

arguments = commandArgs(trailing=TRUE)
prefix = arguments[1]

frip.min=0
frit.min=0
umi.max=Inf
if (length(arguments) > 1) frip.min = as.numeric(arguments[2])
if (length(arguments) > 2) frit.min = as.numeric(arguments[3])
if (length(arguments) > 3) umi.max = as.numeric(arguments[4])

library(monocle3)
library(ggplot2)
library(ggrastr)

cds = readRDS(file.path('checkpoints',paste0(prefix,'_cds.rds')))

cells.to.remove = Reduce(intersect,list(
	rownames(subset(colData(cds),TRUE))[subset(colData(cds),TRUE)$FRIP < frip.min],
	rownames(subset(colData(cds),TRUE))[subset(colData(cds),TRUE)$FRIT < frit.min],
	rownames(subset(colData(cds),TRUE))[subset(colData(cds),TRUE)$n.umi > umi.max]
))

keep.cells = setdiff(rownames(colData(cds)),cells.to.remove)

# cds.keep = cds[,keep.cells]

frip.thresholds = read.delim(file.path('data',paste0(prefix,'_frip_thresholds.txt')))
frit.thresholds = read.delim(file.path('data',paste0(prefix,'_frit_thresholds.txt')))
# umi.thresholds = read.delim(file.path('data',paste0(prefix,'_umi_thresholds.txt')))

frip.to.remove = unlist(lapply(1:nrow(frip.thresholds),function(i) {
	this.id = frip.thresholds$id[i]
	this.threshold = frip.thresholds$threshold[i]
	if (!is.na(this.threshold)) {
		rownames(subset(colData(cds),id == this.id))[subset(colData(cds),id == this.id)$FRIP < this.threshold]
	} else {
		character(0L)
	}
}))

p = ggplot() +
	geom_density(data=as.data.frame(colData(cds)),aes(FRIP)) +
	geom_vline(data=frip.thresholds,aes(xintercept=threshold)) +
	facet_wrap(~id,ncol=ceiling(sqrt(length(unique(colData(cds)$id))))) +
	scale_x_continuous(trans='log10') +
	theme_minimal(base_size=16) +
	theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())
ggsave(p,file=paste0('figures/',prefix,'_qc_frip_thresholds.pdf'),useDingbats=FALSE)

frit.to.remove = unlist(lapply(1:nrow(frit.thresholds),function(i) {
	this.id = frit.thresholds$id[i]
	this.threshold = frit.thresholds$threshold[i]
	if (!is.na(this.threshold)) {
		rownames(subset(colData(cds),id == this.id))[subset(colData(cds),id == this.id)$FRIT < this.threshold]
	} else {
		character(0L)
	}
}))

p = ggplot() +
	geom_density(data=as.data.frame(colData(cds)),aes(FRIT)) +
	geom_vline(data=frit.thresholds,aes(xintercept=threshold)) +
	facet_wrap(~id,ncol=ceiling(sqrt(length(unique(colData(cds)$id))))) +
	scale_x_continuous(trans='log10') +
	theme_minimal(base_size=16) +
	theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())
ggsave(p,file=paste0('figures/',prefix,'_qc_frit_thresholds.pdf'),useDingbats=FALSE)

keep.cells = setdiff(rownames(colData(cds)),cells.to.remove)
keep.cells = setdiff(keep.cells,intersect(frip.to.remove,frit.to.remove))

cds.keep = cds[,keep.cells]

cds.keep = preprocess_cds(cds.keep, method = 'LSI')
cds.keep = reduce_dimension(cds.keep,
	reduction_method = 'UMAP', 
	preprocess_method = 'LSI')

cds.keep = cluster_cells(cds.keep, resolution=1e-5)

saveRDS(cds.keep,file=file.path('checkpoints',paste0(prefix,'_cds_clean.rds')))