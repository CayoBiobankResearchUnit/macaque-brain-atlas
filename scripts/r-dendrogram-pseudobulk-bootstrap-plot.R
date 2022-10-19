#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('rna')

prefix = arguments[1]

library(dendextend)
library(ape)
library(phangorn)
library(phylogram)
library(parallel)
library(data.table)

in.files = file.path('expr_bootstrap/nwk',paste0(prefix,'_cellclass_pseudobulk_',formatC(1:1000,width=4,flag=0),'.nwk'))

cell.bootstrap = mclapply(in.files,function(x) {
	out = read.tree(x)
	out$tip.label = gsub('serotinergic','serotonergic',out$tip.label)
	out
},mc.cores=52)

cell.multiphylo = do.call(c,cell.bootstrap)

cell.tree = read.tree(file.path('expr_bootstrap/nwk',paste0(prefix,'_cellclass_pseudobulk_',formatC(0,width=4,flag=0),'.nwk')))

dir.create('figures/dendrograms',showWarnings=FALSE)

colors = as.matrix(read.table('data/colors.txt',row.names=1,sep='\t',comment.char='',header=FALSE))[,1]
names(colors) = gsub(' ','_',names(colors))

# names(colors) = gsub('serotinergic','serotonergic',names(colors))

cell.tree$tip.label = gsub('serotinergic','serotonergic',cell.tree$tip.label)

colReg = function(n) {
    if (is.leaf(n)) {
        a = attributes(n)
        labCol = as.character(colors[a$label])
        attr(n, 'nodePar') = list(lab.cex = 1, lab.col = labCol, col = labCol, pch=19, cex=0.5)
        attr(n, 'edgePar') = list(col = labCol, lwd=2.5)
    }
    n
}
cell.tree = dendrapply(as.dendrogram(cell.tree), colReg)

edgeReg = list(lwd=2.5)

pdf(file='figures/dendrograms/densitree_pseudobulk_cellclass_expression.pdf',useDingbats=FALSE,height=10)
	densiTree(
		cell.multiphylo,
		alpha=0.005,
		# consensus=as.phylo(cell.tree),
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
	# consensus = as.phylo(cell.tree)
	# xy = ape::plotPhyloCoor(consensus, tip.height = 1:length(unlist(cell.tree)), direction = 'rightwards')
	# xx = xy[,1]; yy = xy[,2]; xx = xx/max(xx); xx = xx + (1 - max(xx))
	# ape::cladogram.plot(consensus$edge, xx, yy, edge.color = '#0000ff', edge.width = 1, edge.lty = 1)
dev.off()


# region.cell.tree = as.dendrogram(hclust(dist(t(region.cell.proportions))))

# region.cell.tree = as.dendrogram(region.cell.newick)

pdf(file='figures/dendrograms/dendrogram_pseudobulk_cellclass_expression.pdf',useDingbats=FALSE,height=10)
	plot(cell.tree)
dev.off()

