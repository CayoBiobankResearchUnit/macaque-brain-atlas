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

meta.file = file.path('umap',paste0('biccn','-metadata-patched.rds'))
this.umap = readRDS(meta.file)

this.umap = subset(this.umap,modality == 'RNA' & !is.na(rna_cell))

keep.cells = with(this.umap,intersect(levels(cell_class),names(which(table(cell_class) > 2000))))

this.umap = droplevels(subset(this.umap,cell_class %in% keep.cells))

region.cell.proportions = apply(as.matrix(with(this.umap,table(cell_class,region))),2,function(x) x/sum(x))

region.cell.tree = as.dendrogram(hclust(dist(t(region.cell.proportions))))

# region.cell.tree = as.dendrogram(read.tree('stats/dendrograms/dendrogram_prop_region_reordered.nwk'))
##### (((((((((((dmPFC:0.02634766052,A1:0.02634766052):0.009608057811,vlPFC:0.03595571833):0.02885432896,IT:0.06481004729):0.0187834031,(M1:0.02637947181,S1:0.02637947181):0.05721397858):0.02465171664,(vmPFC:0.03669493761,ACC:0.03669493761):0.07155022942):0.06625133881,SPP:0.1744965058):0.0290631233,((EC:0.05761297173,AMY:0.05761297173):0.06685205166,HIP:0.1244650234):0.07909460574):0.06071288035,((((dlPFC:0.01199051235,IPP:0.01199051235):0.01694709349,STS:0.02893760583):0.03569717647,PC:0.0646347823):0.02693447169,V1:0.09156925399):0.1727032555):0.1050489813,MT:0.3693214907):0.5165565863,((CN:0.1904597643,NAc:0.1904597643):0.5786940285,((mdTN:0.1335550625,LGN:0.1335550625):0.5535253247,(CC:0.4739481821,(((vlTN:0.14890841,MdC:0.14890841):0.05724255461,MdO:0.2061509646):0.1487240421,(MB:0.1926181009,Pons:0.1926181009):0.1622569058):0.1190731754):0.2131322051):0.08207340565):0.1167242842):0.2548519871,(CV:0.04297087383,lCb:0.04297087383):1.09775919);


# region.cell.newick = read.tree()
# ((lCb:0.04297087383,CV:0.04297087383):1.09775919,((((((Pons:0.1926181009,MB:0.1926181009):0.1622569058,(MdO:0.2061509646,(MdC:0.14890841,vlTN:0.14890841):0.05724255461):0.1487240421):0.1190731754,CC:0.4739481821):0.2131322051,(LGN:0.1335550625,mdTN:0.1335550625):0.5535253247):0.08207340565,(NAc:0.1904597643,CN:0.1904597643):0.5786940285):0.1167242842,(MT:0.3693214907,((V1:0.09156925399,(PC:0.0646347823,(STS:0.02893760583,(IPP:0.01199051235,dlPFC:0.01199051235):0.01694709349):0.03569717647):0.02693447169):0.1727032555,((HIP:0.1244650234,(AMY:0.05761297173,EC:0.05761297173):0.06685205166):0.07909460574,(SPP:0.1744965058,((ACC:0.03669493761,vmPFC:0.03669493761):0.07155022942,((S1:0.02637947181,M1:0.02637947181):0.05721397858,(IT:0.06481004729,(vlPFC:0.03595571833,(A1:0.02634766052,dmPFC:0.02634766052):0.009608057811):0.02885432896):0.0187834031):0.02465171664):0.06625133881):0.0290631233):0.06071288035):0.1050489813):0.5165565863):0.2548519871);
# region.cell.tree = as.dendrogram(region.cell.newick)

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

region.cell.tree = dendrapply(region.cell.tree, colReg)

if (!file.exists('checkpoints/region_cell_proportions_bootstrap.rds')) {
	region.cell.bootstrap = mclapply(1:1000,function(i) {
		as.phylo(hclust(dist(t(apply(as.matrix(with(this.umap[sample(1:nrow(this.umap),replace=TRUE),],table(cell_class,region))),2,function(x) x/sum(x))))))
	},mc.cores=52)
	push.status('bootstrap')
	saveRDS(region.cell.bootstrap,file='checkpoints/region_cell_proportions_bootstrap.rds')
} else {
	region.cell.bootstrap = readRDS('checkpoints/region_cell_proportions_bootstrap.rds')
}

region.cell.multiphylo = do.call(c,region.cell.bootstrap)

dir.create('figures/dendrograms',showWarnings=FALSE)

pdf(file='figures/dendrograms/densitree_pca_region_proportions.pdf',useDingbats=FALSE,height=10)
	densiTree(
		region.cell.multiphylo,
		alpha=0.005,
		consensus=as.phylo(region.cell.tree),
		direction='rightwards',
		scaleX=TRUE,
		col='#000000',
		width=1, lty=1, cex=1, font=1,
		tip.color=rapply(region.cell.tree,function(n) {
			if (is.leaf(n)) attr(n, 'edgePar')[['col']]
		}), srt=0, adj=0,
		label.offset=0,
		scale.bar = FALSE,
		jitter=list(amount = 0.15, random = TRUE)
	)
	# Add the consensus manually
	consensus = as.phylo(region.cell.tree)
	xy = ape::plotPhyloCoor(consensus, tip.height = 1:length(unlist(region.cell.tree)), direction = 'rightwards')
	xx = xy[,1]; yy = xy[,2]; xx = xx/max(xx); xx = xx + (1 - max(xx))
	ape::cladogram.plot(consensus$edge, xx, yy, edge.color = '#0000ff', edge.width = 1, edge.lty = 1)
dev.off()

# region.cell.tree = as.dendrogram(hclust(dist(t(region.cell.proportions))))

# region.cell.tree = as.dendrogram(region.cell.newick)

pdf(file='figures/dendrograms/dendrogram_pca_region_proportions.pdf',useDingbats=FALSE,height=10)
	plot(region.cell.tree)
dev.off()



