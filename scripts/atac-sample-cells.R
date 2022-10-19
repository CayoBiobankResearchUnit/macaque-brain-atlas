#!/usr/bin/env Rscript

source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
this = arguments[1]
n.cells = if (length(arguments) > 1) as.integer(arguments[2]) else 1000

num.cells = n.cells

this.file = file.path('raw-cds',paste0(this,'-cds.rds'))

library(monocle3)
library(data.table)

if (!is.null(n.cells)) {
	cds = readRDS(this.file)

	n.cells.total = nrow(colData(cds))

	if (n.cells.total < n.cells) n.cells = n.cells.total

	set.seed(42)
	sampled.cells = sample(rownames(colData(cds)),n.cells,replace=FALSE)

	cds.sampled = cds[,sampled.cells]

	dir.create(paste0('smp-cds',num.cells),showWarnings=FALSE)

	saveRDS(cds.sampled,file.path(paste0('smp-cds',num.cells),paste0(this,'-cds.rds')))

	dir.create('smp-cells',showWarnings=FALSE)
	write(rownames(colData(cds.sampled)),sep='\n',file=file.path('smp-cells',paste0(this,'-',num.cells,'-cells.txt')))

	# Write new fragment file
	fragment.file = file.path('fragment-files',paste0(this,'-fragments.txt.gz'))

	fragments = fread(fragment.file)

	fragments.sampled = subset(fragments,V4 %in% rownames(colData(cds.sampled)))

	fwrite(fragments.sampled,file=file.path(paste0('smp-cds',num.cells),paste0(this,'-fragments.txt')),sep='\t',row.names=FALSE,col.names=FALSE)

	# Block compress with bgzip
	system(paste('bgzip -f',file.path(paste0('smp-cds',num.cells),paste0(this,'-fragments.txt'))))
} else {
	dir.create(paste0('smp-cds',num.cells),showWarnings=FALSE)

	system(paste('cp',file.path('raw-cds',paste0(this,'-cds.rds')),paste0('smp-cds',num.cells)))
}