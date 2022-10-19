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

sample.list = as.character(read.delim(file.path('data',paste0(prefix,'_metadata.txt')))$id)

this.sample = 'all'

# this.umap = readRDS(file.path('umap',paste0(prefix,'-scanpy-final.rds')))

this.umap = readRDS(file.path('umap',paste0(prefix,'-scanpy-final-classified.rds')))

region.subclass = region.hierarchy$region.subclass
names(region.subclass) = region.hierarchy$region

this.umap$region_subclass = region.subclass[as.character(this.umap$region)]

n.cells = nrow(this.umap)
# Table of cells recovered by region subclass
cells.by.region = table(this.umap$region_subclass)

region.enrichments = do.call(rbind,parallel::mclapply(levels(this.umap$cell_round1_level3),function(this.cell) {
	this.umap.cell = subset(this.umap,cell_round1_level3 == this.cell)
	# Table of cells recovered by region subclass for this cell type
	cells.by.region.this = table(this.umap.cell$region_subclass)
	# Number of cells for this cell type
	n.cell = nrow(this.umap.cell)
	# given how many cells we sequenced from each region, how enriched are they in a given cell type?
	out = do.call(rbind,lapply(levels(this.umap$region_subclass),function(this.region) {
		n.region = cells.by.region[this.region]
		cell.and.region = cells.by.region.this[this.region]
		cell.not.region = n.cell - cell.and.region
		region.not.cell = n.region - cell.and.region
		not.cell.region = n.cells - cell.and.region - cell.not.region - region.not.cell
		contingency.matrix = matrix(c(cell.and.region,cell.not.region,region.not.cell,not.cell.region),nrow=2,dimnames=list(c('in region','not in region'),c('in cell','not in cell')))
		res = fisher.test(contingency.matrix,alternative='greater')
		data.frame(
			cell = this.cell,
			proportion.of.region = cell.and.region/n.region,
			proportion.of.cell = cell.and.region/n.cell,
			region = this.region,
			odds.ratio = res$estimate,
			p.value = res$p.value
		)
	}))
	out$q.value = p.adjust(out$p.value,method='fdr')
	out
},mc.cores=n.cores))

rownames(region.enrichments) = NULL

region.enrichments$cell = factor(region.enrichments$cell,levels=levels(this.umap$cell_round1_level3))

lapply(split(region.enrichments,region.enrichments$cell),function(x) {
	out = x[order(1-(x$q.value < 0.01),-x$odds.ratio),]
	#out = subset(out,q.value < 0.01)
	rownames(out) = NULL
	out$OR = format(round(out$odds.ratio,2),nsmall=2,width=6)
	out$p = format(round(out$p.value,4),nsmall=4)
	out$q = format(round(out$q.value,4),nsmall=4)
	out$prop.region = format(round(out$proportion.of.region,6),nsmall=6)
	out$prop.cell = format(round(out$proportion.of.cell,6),nsmall=6)
	out[,c('region','OR','p','q','prop.region','prop.cell')]
})


library(philentropy)

# Compute Jensen-Shannon divergence
specificity.out = data.frame(cell_type = levels(this.umap$cell_round1_level3))

# Make a table of cell type by region
m = as.matrix(table(this.umap$cell_round1_level3,this.umap$region_major))

# For each cell type, compare its frequency distribution to the full distribution across all regions
specificity.out$specificity_class = unlist(lapply(1:nrow(m),function(i) {
	out = suppressMessages(JSD(rbind(m[i,],colSums(m)),est.prob='empirical'))
	names(out) = rownames(m)[i]
	out
}))

m = as.matrix(table(this.umap$cell_round1_level3,this.umap$region_subclass))

# For each cell type, compare its frequency distribution to the full distribution across all regions
specificity.out$specificity_subclass = unlist(lapply(1:nrow(m),function(i) {
	out = suppressMessages(JSD(rbind(m[i,],colSums(m)),est.prob='empirical'))
	names(out) = rownames(m)[i]
	out
}))

m = as.matrix(table(this.umap$cell_round1_level3,this.umap$region))

# For each cell type, compare its frequency distribution to the full distribution across all regions
specificity.out$specificity_region = unlist(lapply(1:nrow(m),function(i) {
	out = suppressMessages(JSD(rbind(m[i,],colSums(m)),est.prob='empirical'))
	names(out) = rownames(m)[i]
	out
}))


this.umap$animal_region = paste(this.umap$animal_id,this.umap$region,sep='_')

this.umap$animal_region_hemisphere = paste(this.umap$animal_id,this.umap$region,this.umap$hemisphere,sep='_')


cell.proportions = do.call(rbind,lapply(split(this.umap,this.umap$animal_region_hemisphere),function(x) {
	x$id = paste(unique(x$id),collapse=', ')
	meta = unique(x[c('id','region','hemisphere','animal_id','social_group','sex','region_major')])
	data.frame(
		meta,
		n=nrow(meta),
		do.call(data.frame,as.list(table(x$cell_round1_level3)/nrow(x)))
	)
}))

cell.proportions.long = merge(
	unique(cell.proportions[1:8]),
	reshape2::melt(cell.proportions,id='id',measure.vars=names(cell.proportions)[9:ncol(cell.proportions)]),
	by='id'
)

levels(cell.proportions.long$variable) = gsub('\\.',' ',levels(cell.proportions.long$variable))

cell.type.colors = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666', '#ed1c24', '#00aeef', '#86328c', '#00aaad', '#ac7eaf', '#2e3192', '#809ead', '#6a3e14','#ea9423','#b0a690','#feea3f','#fbc9c3','#c5df94','#000000') # (extension of palette Dark2)

cell.proportions.long$is_2C0 = factor(cell.proportions.long$animal_id == '2C0',levels=c('FALSE','TRUE'))

p = ggplot(cell.proportions.long,aes(id,value,fill=variable,alpha=is_2C0)) +
	geom_bar(stat='identity') +
	theme_classic(base_size=6) +
	scale_fill_manual(values=cell.type.colors) +
	scale_alpha_manual(values=c(0.75,1)) +
	facet_grid(cols=vars(region),space='free',scales='free') +
	theme(
		axis.text.x=element_text(angle=-90,hjust=0),
		axis.title.x=element_blank(),
		legend.title=element_blank(),
		strip.background=element_blank(),
		strip.text=element_text(size=10)
	) +
	guides(alpha='none',fill=guide_legend(ncol=1)) +
	ylab('Proportion')
ggsave(file='cell_proportions.pdf',useDingbats=FALSE,height=6,width=20)


