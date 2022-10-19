#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('rna','rna')

prefix = arguments[1]
analysis = arguments[2]

library(dendextend)
library(ape)
library(phangorn)
library(phylogram)
library(parallel)

cell.tree = read.tree(file.path('stats/dendrograms/','dendrogram_pca_cellclasskeep.nwk'))
cell2.tree = read.tree(file.path('stats/dendrograms/','dendrogram_umap_cellclasskeep.nwk'))
region.tree = read.tree(file.path('stats/dendrograms/','dendrogram_pca_region.nwk'))


# cell.tree = read.tree(file.path('stats/dendrograms/','dendrogram_pca_cellclasskeep_reordered.nwk'))
# cell2.tree = read.tree(file.path('stats/dendrograms/','dendrogram_umap_cellclasskeep_reordered.nwk'))
# region.tree = read.tree(file.path('stats/dendrograms/','dendrogram_pca_region_reordered.nwk'))




# region.newick = read.tree()
# (((lCb:0.02,CV:0.02):1.4,((((Pons:0.4,((MdC:0.19,MdO:0.19):0.13,MB:0.32):0.09):0.1,CC:0.5):0.35,(LGN:0.13,(vlTN:0.1,mdTN:0.1):0.03):0.72):0.36,((NAc:0.08,CN:0.08):0.63,(HIP:0.31,(AMY:0.1,EC:0.1):0.22):0.39):0.51):0.2):0.19,(MT:0.88,((V1:0.53,(SPP:0.2,M1:0.2):0.33):0.2,((A1:0.31,ACC:0.31):0.07,(((S1:0.09,(IPP:0.04,STS:0.04):0.05):0.08,dlPFC:0.17):0.12,((IT:0.08,PC:0.08):0.12,(vlPFC:0.1,(vmPFC:0.08,dmPFC:0.08):0.02):0.1):0.09):0.09):0.35):0.15):0.73);
# region.tree = region.newick

# cell.newick = read.tree()
# (((((mesenchymalstemcells:0.71,radialglialcells:0.71):0.34,((ependymalcells:0.82,astrocytes:0.82):0.14,oligodendrocyteprecursorcells:0.96):0.09):0.01,(vascularcells:0.92,microglia:0.92):0.14):0.04,oligodendrocytes:1.1):0.34,(cerebellarneurons:1.15,((serotinergicneurons:0.94,dopaminergicneurons:0.94):0.05,((basketcells:0.82,inhibitoryneurons:0.82):0.13,(mediumspinyneurons:0.73,excitatoryneurons:0.73):0.23):0.03):0.16):0.29);
# cell.tree = cell.newick

cell.class.levels = setdiff(scan(file='stats/clusters/rna-final-cellclasses-levels.txt',what='',sep='\n',quiet=TRUE),c('radial glial cells','mesenchymal stem cells'))
names(cell.class.levels) = gsub(' ','',cell.class.levels)

# Temporary fix
if (any(grepl('neuron$',cell.class.levels))) cell.class.levels[grepl('neuron$',cell.class.levels)] = gsub('neuron$','neurons',cell.class.levels[grepl('neuron$',cell.class.levels)])
if (any(grepl('serotinergic',cell.class.levels))) cell.class.levels[grepl('serotinergic',cell.class.levels)] = gsub('serotinergic','serotonergic',cell.class.levels[grepl('serotinergic',cell.class.levels)])

cell.tree$tip.label = gsub(' ','_',as.character(cell.class.levels[cell.tree$tip.label]))
cell2.tree$tip.label = gsub(' ','_',as.character(cell.class.levels[cell2.tree$tip.label]))

cell.tree = as.dendrogram(cell.tree)
cell2.tree = as.dendrogram(cell2.tree)
region.tree = as.dendrogram(region.tree)

cell.bootstrap = mclapply(
	file.path('stats/dendrograms',paste0('dendrogram_pca_cellclasskeep_bootstrap_',formatC(1:1000,width=6,flag=0),'.nwk')),
	function(x) {
		y = read.tree(x)
		y$tip.label = gsub(' ','_',as.character(cell.class.levels[y$tip.label]))
		as.phylo(y)
	},
	mc.cores = 52
)

cell2.bootstrap = mclapply(
	file.path('stats/dendrograms',paste0('dendrogram_umap_cellclasskeep_bootstrap_',formatC(1:1000,width=6,flag=0),'.nwk')),
	function(x) {
		y = read.tree(x)
		y$tip.label = gsub(' ','_',as.character(cell.class.levels[y$tip.label]))
		as.phylo(y)
	},
	mc.cores = 52
)

region.bootstrap = mclapply(
	file.path('stats/dendrograms',paste0('dendrogram_pca_region_bootstrap_',formatC(1:1000,width=6,flag=0),'.nwk')),
	function(x) {
		as.phylo(read.tree(x))
	},
	mc.cores = 52
)

colors = as.matrix(read.table('data/colors.txt',row.names=1,sep='\t',comment.char='',header=FALSE))[,1]

names(colors) = gsub(' ','_',names(colors))

region.levels = c('dmPFC','vmPFC','dlPFC','vlPFC','ACC','M1','EC','PC','A1','STS','MT','IT','S1','IPP','SPP','V1','CC','CN','NAc','AMY','HIP','mdTN','vlTN','LGN','CV','lCb','MB','MdO','MdC','Pons')

colReg = function(n) {
    if (is.leaf(n)) {
        a = attributes(n)
        labCol = as.character(colors[a$label])
        attr(n, 'nodePar') = list(lab.cex = 1, lab.col = labCol, col = labCol, pch=19, cex=0.5)
        attr(n, 'edgePar') = list(col = labCol, lwd=2.5)
    }
    n
}

edgeReg = list(lwd=2.5)

cell.tree = dendrapply(cell.tree, colReg)
cell2.tree = dendrapply(cell2.tree, colReg)
region.tree = dendrapply(region.tree, colReg)

cell.bootstrap.multiphylo = do.call(c,cell.bootstrap)
cell2.bootstrap.multiphylo = do.call(c,cell2.bootstrap)
region.bootstrap.multiphylo = do.call(c,region.bootstrap)

dir.create('figures/dendrograms',showWarnings=FALSE)

pdf(file='figures/dendrograms/densitree_pca_cellclasskeep.pdf',useDingbats=FALSE,height=10)
	densiTree(
		cell.bootstrap.multiphylo,
		alpha=0.005,
		consensus=as.phylo(cell.tree),
		direction='rightwards',
		scaleX=TRUE,
		col='#000000',
		width=1, lty=1, cex=1, font=1,
		tip.color=rapply(cell.tree,function(n) {
			if (is.leaf(n)) attr(n, 'edgePar')[['col']]
		}), srt=0, adj=0,
		label.offset=0,
		scale.bar = FALSE,
		jitter=list(amount = 0.15, random = TRUE)
	)
	# Add the consensus manually
	consensus = as.phylo(cell.tree)
	xy = ape::plotPhyloCoor(consensus, tip.height = 1:length(unlist(cell.tree)), direction = 'rightwards')
	xx = xy[,1]; yy = xy[,2]; xx = xx/max(xx); xx = xx + (1 - max(xx))
	ape::cladogram.plot(consensus$edge, xx, yy, edge.color = '#0000ff', edge.width = 1, edge.lty = 1)
dev.off()

pdf(file='figures/dendrograms/densitree_umap_cellclasskeep.pdf',useDingbats=FALSE,height=10)
	densiTree(
		cell2.bootstrap.multiphylo,
		alpha=0.005,
		consensus=as.phylo(cell2.tree),
		direction='rightwards',
		scaleX=TRUE,
		col='#000000',
		width=1, lty=1, cex=1, font=1,
		tip.color=rapply(cell2.tree,function(n) {
			if (is.leaf(n)) attr(n, 'edgePar')[['col']]
		}), srt=0, adj=0,
		label.offset=0,
		scale.bar = FALSE,
		jitter=list(amount = 0.15, random = TRUE)
	)
	# Add the consensus manually
	consensus = as.phylo(cell2.tree)
	xy = ape::plotPhyloCoor(consensus, tip.height = 1:length(unlist(cell2.tree)), direction = 'rightwards')
	xx = xy[,1]; yy = xy[,2]; xx = xx/max(xx); xx = xx + (1 - max(xx))
	ape::cladogram.plot(consensus$edge, xx, yy, edge.color = '#0000ff', edge.width = 1, edge.lty = 1)
dev.off()

pdf(file='figures/dendrograms/densitree_pca_region.pdf',useDingbats=FALSE,height=10)
	densiTree(
		region.bootstrap.multiphylo,
		alpha=0.005,
		consensus=as.phylo(region.tree),
		direction='rightwards',
		scaleX=TRUE,
		col='#000000',
		width=1, lty=1, cex=1, font=1,
		tip.color=rapply(region.tree,function(n) {
			if (is.leaf(n)) attr(n, 'edgePar')[['col']]
		}), srt=0, adj=0,
		label.offset=0,
		scale.bar = FALSE,
		jitter=list(amount = 0.15, random = TRUE)
	)
	# Add the consensus manually
	consensus = as.phylo(region.tree)
	xy = ape::plotPhyloCoor(consensus, tip.height = 1:length(unlist(region.tree)), direction = 'rightwards')
	xx = xy[,1]; yy = xy[,2]; xx = xx/max(xx); xx = xx + (1 - max(xx))
	ape::cladogram.plot(consensus$edge, xx, yy, edge.color = '#0000ff', edge.width = 1, edge.lty = 1)
dev.off()




# region.tree = as.dendrogram(region.newick)

pdf(file='figures/dendrograms/dendrogram_pca_region.pdf',useDingbats=FALSE,height=10)
	plot(region.tree)
dev.off()

pdf(file='figures/dendrograms/dendrogram_umap_cellclasskeep.pdf',useDingbats=FALSE,height=10)
	plot(cell2.tree)
dev.off()

pdf(file='figures/dendrograms/dendrogram_pca_cellclasskeep.pdf',useDingbats=FALSE,height=10)
	plot(cell.tree)
dev.off()

