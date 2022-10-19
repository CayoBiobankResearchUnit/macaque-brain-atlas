#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('rna','rna')

prefix = arguments[1]
analysis = arguments[2]

suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(primEnrich))
suppressMessages(library(parallel))

# Cell class levels

dir.create('rds/markers',showWarnings=FALSE)
cell.classes = scan(file.path('stats/clusters',paste0(prefix,'-final-cellclasses-levels.txt')),what='',sep='\n')

marker.genes = parallel::mclapply(seq_along(cell.classes)[!cell.classes %in% c('radial glial cells','mesenchymal stem cells','AHSG neuron','F5 neuron','KIR3DL12 neuron','KIR3DL12 microglia')],function(i) {
	data.frame(cell_class=cell.classes[i],readRDS(file.path('rds/',paste0(prefix,'_class',i,'_marker_genes.rds'))))
},mc.cores=n.cores)

cell.classes.keep = cell.classes[!cell.classes %in% c('radial glial cells','mesenchymal stem cells','AHSG neuron','F5 neuron','KIR3DL12 neuron','KIR3DL12 microglia')]
names(marker.genes) = cell.classes.keep

cell.types = unlist(mclapply(cell.classes.keep,function(i) {
	x = marker.genes[[i]]
	cell.subtypes = unique(x$cell)
	if (length(cell.subtypes) > 1) {
		cell.subtype.integers = as.integer(gsub(paste0(i,' '),'',cell.subtypes))
		paste(i,sort(cell.subtype.integers))
	} else {
		i
	}
},mc.cores=n.cores))

marker.genes = do.call(rbind,marker.genes)
rownames(marker.genes) = NULL

marker.genes = within(marker.genes,{
	cell_class = factor(cell_class,levels=cell.classes.keep)
	cell = factor(cell,levels=cell.types)
	logfoldchanges = as.numeric(logfoldchanges)
	wilcox_score = as.numeric(wilcox_score)
	ttest_score = as.numeric(ttest_score)
	logreg_score = as.numeric(logreg_score)
	bp1 = as.integer(bp1)
	bp2 = as.integer(bp2)
	gene_strand = factor(gene_strand,levels=c('+','-'))
})

saveRDS(marker.genes,file='rds/subclass_marker_genes.rds')

marker.genes.sorted = marker.genes[with(marker.genes,order(cell,-logreg_score)),]

# Drop cell subtypes that do not end in a number (cell classes that have only one subtype are pointless to include and will lead to errors)
marker.genes.sorted = droplevels(subset(marker.genes.sorted,cell %in% levels(marker.genes.sorted$cell)[grepl('[0-9]$',levels(marker.genes.sorted$cell))]))

disease.ks.results = vector(mode='list',length=nlevels(marker.genes.sorted$cell))
names(disease.ks.results) = levels(marker.genes.sorted$cell)
for (i in levels(marker.genes.sorted$cell)) {
	message(i)
	x = subset(marker.genes.sorted,cell==i)
	test.values = x$logreg_score
	names(test.values) = x$id
	kst = primEnrich(
		test.values,
		alternative='greater',
		statistic='ks',
		annotation = macaqueDISEASES,
		n.cores = n.cores
	)
	disease.ks.results[[i]] = data.frame(cell = i,kst)
}

disease.ks.results = do.call(rbind,disease.ks.results)
rownames(disease.ks.results) = NULL
disease.ks.results$cell_class = gsub(' [0-9]+$','',disease.ks.results$cell)
saveRDS(disease.ks.results,file=paste0('rds/',prefix,'-recluster-marker-gene-enrichment-disease-ks.rds'))


disease.ft.results = vector(mode='list',length=nlevels(marker.genes.sorted$cell))
names(disease.ft.results) = levels(marker.genes.sorted$cell)
for (i in levels(marker.genes.sorted$cell)) {
	message(i)
	x = subset(marker.genes.sorted,cell==i)
	test.genes = x$id[1:100]
	fet = primEnrich(
		test.genes,
		background=x$id,
		alternative='greater',
		statistic='fisher',
		annotation = macaqueDISEASES,
		n.cores = n.cores
	)
	disease.ft.results[[i]] = data.frame(cell = i,fet)
}
disease.ft.results = do.call(rbind,disease.ft.results)
rownames(disease.ft.results) = NULL
disease.ft.results$cell_class = gsub(' [0-9]+$','',disease.ft.results$cell)

saveRDS(disease.ft.results,file=paste0('rds/',prefix,'-recluster-marker-gene-enrichment-disease-ft.rds'))



go.ks.results = vector(mode='list',length=nlevels(marker.genes.sorted$cell))
names(go.ks.results) = levels(marker.genes.sorted$cell)
for (i in levels(marker.genes.sorted$cell)) {
	message(i)
	x = subset(marker.genes.sorted,cell==i)
	test.values = x$logreg_score
	names(test.values) = x$id
	kst = topgoKS(
		test.values,
		alternative='greater',
		go.ontology='BP',
		annotation = macaqueGO
	)
	go.ks.results[[i]] = data.frame(cell = i,kst)
}
go.ks.results = do.call(rbind,go.ks.results)
rownames(go.ks.results) = NULL
go.ks.results$cell_class = gsub(' [0-9]+$','',go.ks.results$cell)
saveRDS(go.ks.results,file=paste0('rds/',prefix,'-recluster-marker-gene-enrichment-go-ks.rds'))

go.ft.results = vector(mode='list',length=nlevels(marker.genes.sorted$cell))
names(go.ft.results) = levels(marker.genes.sorted$cell)
for (i in levels(marker.genes.sorted$cell)) {
	message(i)
	x = subset(marker.genes.sorted,cell==i)
	test.genes = x$id[1:100]
	fet = topgoFisher(
		test.genes,
		background=x$id,
		alternative='greater',
		go.ontology='BP',
		annotation = macaqueGO
	)
	go.ft.results[[i]] = data.frame(cell = i,fet)
}
go.ft.results = do.call(rbind,go.ft.results)
rownames(go.ft.results) = NULL
go.ft.results$cell_class = gsub(' [0-9]+$','',go.ft.results$cell)
saveRDS(go.ft.results,file=paste0('rds/',prefix,'-recluster-marker-gene-enrichment-go-ft.rds'))





# Focus just on cerebellar neurons
x = split(marker.genes.sorted,marker.genes.sorted$cell)[['cerebellar neurons 16']]
test.values = x$logreg_score
names(test.values) = x$id
kst = primEnrich(
	test.values,
	alternative='greater',
	statistic='ks',
	annotation = macaqueDISEASES,
	n.cores = n.cores
)

# Run Fisher tests over a range of N (top N marker genes)
autism.results = do.call(rbind,lapply(c(10,25,50,100,200,500,1000,2000,5000),function(i) {
	test.genes = x$id[1:i] # x$id[x$ttest_pvals_adj < 0.05 & x$logreg_score > 0]
	out = primEnrich(
		test.genes,
		background=x$id,
		alternative='greater',
		statistic='fisher',
		annotation = macaqueDISEASES
	)
	this.rank = which(out$pathway_name == 'Autistic disorder')
	data.frame(subset(out,pathway_name == 'Autistic disorder'),n = i,rank=this.rank)
}))

