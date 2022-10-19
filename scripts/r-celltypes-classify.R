#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('rna','rna')

prefix = arguments[1]
analysis = arguments[2]

this.sample = 'all'

this.umap = readRDS(file.path('umap',paste0(prefix,'-scanpy-final.rds')))

if (analysis == 'rna') {
	# this.umap=readRDS('umap/rna-scanpy-final_mindist_0.25.rds')

	# > colSums(table(this.umap$cluster,this.umap$partition) > 0)
	#  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
	# 29 16  2 11  8  1  2  1  2  4  1  1  1  1  1

	# Step one: classify with as much specificity as possible
	partition.to.cell = c(
	'1' = NA,
	'2' = 'oligodendrocytes',
	'3' = NA,
	'4' = 'astrocytes',
	'5' = NA,
	'6' = 'oligodendrocyte precursor cells',
	'7' = 'radial glial cells',
	'8' = 'basket cells',
	'9' = 'microglia',
	'10' = 'medium spiny neurons',
	'11' = NA, # endothelial cells and pericytes
	'12' = 'mesenchymal stem cells',
	'13' = 'ependymal cells',
	'14' = 'partition 14', # AHSG neurons?
	'15' = 'partition 15', # F5 neurons?
	'16' = 'partition 16', # dopaminergic neurons?
	'17' = 'partition 17' # microglia?
	)

	# this.umap$cell_type_round1 = this.umap$cell_type

	this.umap$cell_final = NA

	this.umap$cell_final = partition.to.cell[as.character(this.umap$partition)]

	# lapply(split(subset(this.umap,is.na(cell_final)),subset(this.umap,is.na(cell_final))$cluster),function(x) {
	# 	y = subset(x,!is.na(cell_type))
	# 	sort(table(y$cell_type) / nrow(y),decreasing=TRUE)[1:5] * 100
	# })

	cluster.to.cell = c(
	'1' = 'excitatory neurons',
	'3' = 'cerebellar granule cells',
	'4' = 'cerebellar granule cells',
	'6' = 'excitatory neurons',
	'7' = 'excitatory neurons',
	'8' = 'excitatory neurons',
	'9' = 'excitatory neurons',
	'11' = 'excitatory neurons',
	'12' = 'excitatory neurons',
	'13' = 'inhibitory neurons',
	'14' = 'inhibitory neurons',
	'15' = 'cerebellar granule cells',
	'17' = 'excitatory neurons',
	'18' = 'inhibitory neurons',
	'19' = 'excitatory neurons',
	'20' = 'inhibitory neurons',
	'22' = 'excitatory neurons',
	'24' = 'excitatory neurons',
	'25' = 'excitatory neurons',
	'26' = 'cerebellar granule cells',
	'29' = 'inhibitory neurons',
	'31' = 'excitatory neurons',
	'32' = 'inhibitory neurons',
	'34' = 'endothelial cells',
	'36' = 'excitatory neurons',
	'37' = 'inhibitory neurons',
	'38' = 'thalamic interneurons',
	'39' = 'excitatory neurons',
	'40' = 'pericytes',
	'41' = 'cerebellar granule cells',
	'42' = 'cerebellar granule cells',
	'43' = 'excitatory neurons',
	'44' = 'excitatory neurons',
	'46' = 'excitatory neurons',
	'47' = 'excitatory neurons',
	'48' = 'excitatory neurons',
	'49' = 'excitatory neurons',
	'50' = 'endothelial cells',
	'51' = 'endothelial cells',
	'52' = 'cerebellar granule cells',
	'54' = 'excitatory neurons',
	'55' = 'excitatory neurons',
	'56' = 'cerebellar granule cells',
	'57' = 'excitatory neurons',
	'58' = 'cerebellar granule cells',
	'59' = 'excitatory neurons',
	'61' = 'Purkinje cells',
	'62' = 'excitatory neurons',
	'63' = 'pericytes',
	'64' = 'excitatory neurons',
	'66' = 'dopaminergic neurons',
	'67' = 'serotinergic neurons',
	'68' = 'excitatory neurons',
	'70' = 'excitatory neurons',
	'71' = 'cerebellar granule cells',
	'75' = 'excitatory neurons',
	'77' = 'inhibitory neurons',
	'82' = 'excitatory neurons'
	)

	this.umap$cell_final[is.na(this.umap$cell_final)] = cluster.to.cell[as.character(this.umap$cluster[is.na(this.umap$cell_final)])]

	this.umap$cell_final = factor(this.umap$cell_final,levels=c(
	'excitatory neurons',
	'medium spiny neurons',
	'inhibitory neurons',
	'dopaminergic neurons',
	'serotinergic neurons',
	'thalamic interneurons',
	'cerebellar granule cells',
	'Purkinje cells',
	'basket cells',
	'astrocytes',
	'oligodendrocytes',
	'oligodendrocyte precursor cells',
	'endothelial cells',
	'pericytes',
	'microglia',
	'radial glial cells',
	'mesenchymal stem cells',
	'ependymal cells',
	'partition 14',
	'partition 15',
	'partition 16',
	'partition 17'
	),labels=c(
	'excitatory neurons',
	'medium spiny neurons',
	'inhibitory neurons',
	'dopaminergic neurons',
	'serotinergic neurons',
	'thalamic interneurons',
	'cerebellar granule cells',
	'Purkinje cells',
	'basket cells',
	'astrocytes',
	'oligodendrocytes',
	'oligodendrocyte precursor cells',
	'endothelial cells',
	'pericytes',
	'microglia',
	'radial glial cells',
	'mesenchymal stem cells',
	'ependymal cells',
	'partition 14',
	'partition 15',
	'partition 16',
	'partition 17'
	))

	this.umap$cell_round1_level3 = dplyr::recode(this.umap$cell_final,
	'excitatory neurons' = 'excitatory neurons',
	'medium spiny neurons' = 'medium spiny neurons',
	'inhibitory neurons' = 'inhibitory neurons',
	'dopaminergic neurons' = 'dopaminergic neurons',
	'serotinergic neurons' = 'serotinergic neurons',
	'thalamic interneurons' = 'thalamic interneurons',
	'cerebellar granule cells' = 'cerebellar granule cells',
	'Purkinje cells' = 'Purkinje cells',
	'basket cells' = 'basket cells',
	'astrocytes' = 'astrocytes',
	'oligodendrocytes' = 'oligodendrocytes',
	'oligodendrocyte precursor cells' = 'oligodendrocyte precursor cells',
	'endothelial cells' = 'endothelial cells',
	'pericytes' = 'pericytes',
	'microglia' = 'microglia',
	'radial glial cells' = 'radial glial cells',
	'mesenchymal stem cells' = 'mesenchymal stem cells',
	'ependymal cells' = 'ependymal cells',
	'partition 14' = 'AHSG neuron',
	'partition 15' = 'F5 neuron',
	'partition 16' = 'KIR3DL12 neuron',
	'partition 17' = 'KIR3DL12 microglia'
	)

	this.umap$cell_round1_level2 = dplyr::recode(this.umap$cell_final,
	'excitatory neurons' = 'excitatory neurons',
	'medium spiny neurons' = 'medium spiny neurons',
	'inhibitory neurons' = 'inhibitory neurons',
	'dopaminergic neurons' = 'dopaminergic neurons',
	'serotinergic neurons' = 'serotinergic neurons',
	'thalamic interneurons' = 'inhibitory neurons',
	'cerebellar granule cells' = 'cerebellar granule cells',
	'Purkinje cells' = 'Purkinje cells',
	'basket cells' = 'basket cells',
	'astrocytes' = 'astrocytes',
	'oligodendrocytes' = 'oligodendrocytes',
	'oligodendrocyte precursor cells' = 'oligodendrocyte precursor cells',
	'endothelial cells' = 'vascular cells',
	'pericytes' = 'vascular cells',
	'microglia' = 'microglia',
	'radial glial cells' = 'radial glial cells',
	'mesenchymal stem cells' = 'mesenchymal stem cells',
	'ependymal cells' = 'ependymal cells',
	'partition 14' = 'AHSG neuron',
	'partition 15' = 'F5 neuron',
	'partition 16' = 'KIR3DL12 neuron',
	'partition 17' = 'KIR3DL12 microglia'
	)

	this.umap$cell_round1_level1 = dplyr::recode(this.umap$cell_final,
	'excitatory neurons' = 'excitatory neurons',
	'medium spiny neurons' = 'medium spiny neurons',
	'inhibitory neurons' = 'inhibitory neurons',
#	'dopaminergic neurons' = 'inhibitory neurons',
#	'serotinergic neurons' = 'inhibitory neurons',
	'dopaminergic neurons' = 'dopaminergic neurons',
	'serotinergic neurons' = 'serotinergic neurons',
	'thalamic interneurons' = 'inhibitory neurons',
	'cerebellar granule cells' = 'cerebellar neurons',
	'Purkinje cells' = 'cerebellar neurons',
	'basket cells' = 'basket cells',
	'astrocytes' = 'astrocytes',
	'oligodendrocytes' = 'oligodendrocytes',
	'oligodendrocyte precursor cells' = 'oligodendrocyte precursor cells',
	'endothelial cells' = 'vascular cells',
	'pericytes' = 'vascular cells',
	'microglia' = 'microglia',
	'radial glial cells' = 'radial glial cells',
	'mesenchymal stem cells' = 'mesenchymal stem cells',
	'ependymal cells' = 'ependymal cells',
	'partition 14' = 'AHSG neuron',
	'partition 15' = 'F5 neuron',
	'partition 16' = 'KIR3DL12 neuron',
	'partition 17' = 'KIR3DL12 microglia'
	)

	write.table(this.umap[,c('cell','cluster')],file=file.path('stats/clusters',paste0(prefix,'-final-clusters.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
	write.table(this.umap[,c('cell','partition')],file=file.path('stats/clusters',paste0(prefix,'-final-partitions.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
	write.table(this.umap[,c('cell','cell_round1_level3')],file=file.path('stats/clusters',paste0(prefix,'-final-celltypes.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
	write.table(this.umap[,c('cell','cell_round1_level1')],file=file.path('stats/clusters',paste0(prefix,'-final-cellclasses.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
	# write.table(this.umap[,c('cell','cluster_unclassified')],file='rna_cell_clusterstodo.txt',sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
	# write.table(this.umap[,c('cell','cluster_classified')],file=file.path('stats/clusters',paste0(prefix,'-final-celltypes.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

	write(levels(this.umap[,c('cell_round1_level3')]),file=file.path('stats/clusters',paste0(prefix,'-final-celltypes-levels.txt')),sep='\n')
	write(levels(this.umap[,c('cell_round1_level1')]),file=file.path('stats/clusters',paste0(prefix,'-final-cellclasses-levels.txt')),sep='\n')

	saveRDS(this.umap,file=file.path('umap',paste0(prefix,'-scanpy-final-classified.rds')))

	cell.type.colors = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666', '#ed1c24', '#00aeef', '#86328c', '#00aaad', '#ac7eaf', '#2e3192', '#809ead', '#6a3e14') # (extension of palette Dark2)
	cell.type.colors = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666', '#ed1c24', '#00aeef', '#86328c', '#00aaad', '#ac7eaf', '#2e3192', '#809ead', '#6a3e14','#ea9423','#b0a690','#feea3f','#fbc9c3','#c5df94','#000000') # (extension of palette Dark2)

	set.seed(42)
	this.umap = this.umap[sample(1:nrow(this.umap)),]

	plot.umap(this.umap,color='cell_final',file=file.path('figures/umap-final',paste0('umap-',prefix,'-cellfinal-all.pdf')),color.label='Cell Type',alpha=0.05)
	plot.umap(this.umap,color='cell_final',file=file.path('figures/umap-final',paste0('umap-',prefix,'-cellfinal-all_nolegend.pdf')),color.label='Cell Type',legend=FALSE,alpha=0.05)
	plot.umap(this.umap,color='cell_final',file=file.path('figures/umap-final',paste0('umap-',prefix,'-cellfinal-all_presentation.pdf')),color.label='Cell Type',legend=FALSE,alpha=0.05,presentation=TRUE)
	plot.umap(this.umap,color='cell_final',file=file.path('figures/umap-final',paste0('umap-',prefix,'-cellfinal-facet-all.pdf')),color.label='Cell Type',facet=TRUE,facet.by='cell_final',legend=FALSE)

} else if (analysis == 'atac') {

	partition.to.cell = c(
	'1' = 'excitatory neurons',
	'2' = NA,
	'3' = 'inhibitory neurons',
	'4' = 'oligodendrocytes',
	'5' = 'cerebellar neurons',
	'6' = 'medium spiny neurons',
	'7' = 'oligodendrocyte precursor cells',
	'8' = 'microglia',
	'9' = 'basket cells' # provisional
	)

	this.umap$cell_final = NA

	this.umap$cell_final = partition.to.cell[as.character(this.umap$partition)]

	cluster.to.cell = c(
	'4' = 'astrocytes',
	'8' = 'astrocytes',
	'15' = 'vascular cells',
	'16' = 'astrocytes',
	'34' = 'astrocytes',
	'38' = 'astrocytes'
	)

	this.umap$cell_final[is.na(this.umap$cell_final)] = cluster.to.cell[as.character(this.umap$cluster[is.na(this.umap$cell_final)])]

	this.umap$cell_final = factor(this.umap$cell_final,levels=c(
	'excitatory neurons',
	'medium spiny neurons',
	'inhibitory neurons',
	'cerebellar neurons',
	'basket cells',
	'astrocytes',
	'oligodendrocytes',
	'oligodendrocyte precursor cells',
	'vascular cells',
	'microglia'
	),labels=c(
	'excitatory neurons',
	'medium spiny neurons',
	'inhibitory neurons',
	'cerebellar neurons',
	'basket cells',
	'astrocytes',
	'oligodendrocytes',
	'oligodendrocyte precursor cells',
	'vascular cells',
	'microglia'
	))

	this.umap$cell_round1_level1 = this.umap$cell_round1_level2 = this.umap$cell_round1_level3 = this.umap$cell_final

	write.table(this.umap[,c('cell','cluster')],file=file.path('stats/clusters',paste0(prefix,'-final-clusters.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
	write.table(this.umap[,c('cell','partition')],file=file.path('stats/clusters',paste0(prefix,'-final-partitions.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
	write.table(this.umap[,c('cell','cell_round1_level3')],file=file.path('stats/clusters',paste0(prefix,'-final-celltypes.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
	write.table(this.umap[,c('cell','cell_round1_level1')],file=file.path('stats/clusters',paste0(prefix,'-final-cellclasses.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
	# write.table(this.umap[,c('cell','cluster_unclassified')],file='rna_cell_clusterstodo.txt',sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
	# write.table(this.umap[,c('cell','cluster_classified')],file=file.path('stats/clusters',paste0(prefix,'-final-celltypes.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

	write(levels(this.umap[,c('cell_round1_level3')]),file=file.path('stats/clusters',paste0(prefix,'-final-celltypes-levels.txt')),sep='\n')
	write(levels(this.umap[,c('cell_round1_level1')]),file=file.path('stats/clusters',paste0(prefix,'-final-cellclasses-levels.txt')),sep='\n')

	saveRDS(this.umap,file=file.path('umap',paste0(prefix,'-scanpy-final-classified.rds')))

	cell.type.colors = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666', '#ed1c24', '#00aeef', '#86328c', '#00aaad', '#ac7eaf', '#2e3192', '#809ead', '#6a3e14') # (extension of palette Dark2)
	cell.type.colors = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666', '#ed1c24', '#00aeef', '#86328c', '#00aaad', '#ac7eaf', '#2e3192', '#809ead', '#6a3e14','#ea9423','#b0a690','#feea3f','#fbc9c3','#c5df94','#000000') # (extension of palette Dark2)

	set.seed(42)
	this.umap = this.umap[sample(1:nrow(this.umap)),]

	plot.umap(this.umap,color='cell_final',file=file.path('figures/umap-final',paste0('umap-',prefix,'-cellfinal-all.pdf')),color.label='Cell Type',alpha=0.05)
	plot.umap(this.umap,color='cell_final',file=file.path('figures/umap-final',paste0('umap-',prefix,'-cellfinal-all_nolegend.pdf')),color.label='Cell Type',legend=FALSE,alpha=0.05)
	plot.umap(this.umap,color='cell_final',file=file.path('figures/umap-final',paste0('umap-',prefix,'-cellfinal-all_presentation.pdf')),color.label='Cell Type',legend=FALSE,alpha=0.05,presentation=TRUE)
	plot.umap(this.umap,color='cell_final',file=file.path('figures/umap-final',paste0('umap-',prefix,'-cellfinal-facet-all.pdf')),color.label='Cell Type',facet=TRUE,facet.by='cell_final',legend=FALSE)

}

# plot.umap(subset(this.umap,region==''),color='cell_final',file=file.path('figures/umap-final',paste0('umap-',prefix,'-cellfinal-all.pdf')),color.label='Cell Type',alpha=0.05)
# 
# 
# 
# 
# this.umap$cluster_unclassified = as.character(this.umap$cluster)
# this.umap$cluster_unclassified[!is.na(this.umap$cell_final)] = '0'
# 
# plot.umap(this.umap,color='cluster_unclassified',file=file.path('figures/umap-final_mindist_0.25',paste0('umap-',prefix,'-celltodo-all.pdf')),color.label='Cell Type',alpha=0.05)
# plot.umap(this.umap,color='cluster_unclassified',file=file.path('figures/umap-final_mindist_0.25',paste0('umap-',prefix,'-celltodo-all_nolegend.pdf')),color.label='Cell Type',legend=FALSE,alpha=0.05)
# plot.umap(this.umap,color='cluster_unclassified',file=file.path('figures/umap-final_mindist_0.25',paste0('umap-',prefix,'-celltodo-all_presentation.pdf')),color.label='Cell Type',legend=FALSE,alpha=0.05,presentation=TRUE)
# plot.umap(this.umap,color='cluster_unclassified',file=file.path('figures/umap-final_mindist_0.25',paste0('umap-',prefix,'-celltodo-facet-all.pdf')),color.label='Cell Type',facet=TRUE,facet.by='cluster_unclassified',legend=FALSE)
