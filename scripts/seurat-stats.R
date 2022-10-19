#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
prefix = arguments[1]
analysis = arguments[2]
n.cells = if (length(arguments) < 3) '' else as.numeric(arguments[3])

if (!analysis %in% c('rna','atac')) {
	stop('Argument 2 must be either "rna" or "atac"')
}

library(monocle3)
library(Seurat)

cds.atac = readRDS(file.path('checkpoints',paste0(prefix,'_',n.cells,'_cds_seurat.rds')))

# Number of candidate CRE
nrow(cds.atac)

# Stats on features
Annotation(cds.atac)

tss = read.delim('data/tss.txt')

i.tss = with(tss,GRanges(chr,IRanges(tss,width=1),strand=strand))
i.promoters = with(tss,GRanges(chr,IRanges(tss-1000,tss+1000),strand=strand))

atac.features = rownames(cds.atac@assays$ATAC@counts)

atac.features.split = strsplit(atac.features,'-')

atac.features = do.call(rbind,lapply(atac.features.split,function(x) {
	out = as.data.frame(matrix(x,nrow=1))
	names(out) = c('chr','start','stop')
	out
}))

atac.features = within(atac.features,{
	start = as.numeric(start)
	stop = as.numeric(stop)
})

atac.features$size = with(atac.features,stop-start + 1)

i.features = with(atac.features,GRanges(chr,IRanges(start,stop),'*'))

distToNearest = GenomicRanges::distanceToNearest(i.features,Annotation(cds.atac),ignore.strand=TRUE)

distToNearestTSS = GenomicRanges::distanceToNearest(i.features,i.tss,ignore.strand=TRUE)

distToNearestPromoter = GenomicRanges::distanceToNearest(i.features,i.promoters,ignore.strand=TRUE)

library(ggplot2)

atac.features = data.frame(
	atac.features,
	d2gene=distToNearest@elementMetadata$distance,
	d2TSS=distToNearestTSS@elementMetadata$distance,
	d2promoter=distToNearestPromoter@elementMetadata$distance
)

p = ggplot(atac.features,aes(d2TSS+1)) + geom_histogram(bins=50) + scale_x_continuous(trans='log10') + theme_classic()

Annotation(cds.atac)







colData(cds.rna)$CellType_1C <- dplyr::recode(colData(cds.rna)$CellType_1B,
                                                         "Astrocytes"="Astrocytes",
                                                         "Endothelial"="Endothelial",
                                                         "Immune"="Endothelial",
                                                         "Microglia"="Microglia",
                                                         "Neurons_Basket"="Neurons_Basket",
                                                         "Neurons_Excitatory"="Neurons_Excitatory",
                                                         "Neurons_Granule"="Neurons_Granule",
                                                         "Neurons_Inhibitory_ADARB2"="Neurons_Inhibitory",
                                                         "Neurons_Inhibitory_LHX6/PVALB"="Neurons_Inhibitory",
                                                         "Neurons_Inhibitory_LHX6/SST"="Neurons_Inhibitory",
                                                         "Neurons_MediumSpiny"="Neurons_MediumSpiny",
                                                         "Neurons_Purkinje"="Neurons_Purkinje",
                                                         "Neurons_Thalamic"="Neurons_Thalamic",
                                                         "Oligodendrocyte Precursors"="Oligodendrocyte Precursors",
                                                         "Oligodendrocytes"="Oligodendrocytes")

table(cds.rna$CellType_1C)
table(cds.rna$Cell_Class_Oct_NeuronSubClasses)


library(philentropy)
# Make a table of cell type by region
m = as.matrix(table(cds.rna$CellType_1C,cds.rna$Region))

# For each cell type, compare its frequency distribution to the full distribution across all regions
cell.type.specificity = unlist(lapply(1:nrow(m),function(i) {
	out = JSD(rbind(m[i,],colSums(m)),est.prob='empirical')
	names(out) = rownames(m)[i]
	out
}))


# Cell types per region corrected for sampling
m.corrected = apply(m,2,function(x) x/sum(x))

100*(rowSums(m.corrected) / sum(rowSums(m.corrected)))

library('philentropy')