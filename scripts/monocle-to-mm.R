#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('atac','atac','1')
prefix = arguments[1]
analysis = arguments[2]
this.i = as.integer(arguments[3])

if (!analysis %in% c('rna','atac')) {
	stop('Argument 2 must be either "rna" or "atac"')
}

suppressMessages(library(monocle3))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))

sample.list = as.character(read.delim(file.path('data',paste0(prefix,'_metadata.txt')))$id)

this.sample = sample.list[this.i]

if (analysis == 'atac') {
	cds.monocle = readRDS(file.path(paste0('raw-cds'),paste0(this.sample,'-cds.rds')))
} else if (analysis == 'rna') {
	cds.monocle = readRDS(file.path(paste0('checkpoints/',prefix,'_monocle3_cds'),paste0(this.sample,'_cds.RDS')))
}
dir.create('mm',showWarnings=FALSE)

if (analysis == 'rna') {
	# Don't write the matrix for ATAC because the peak set should be unified across all runs
	writeMM(exprs(cds.monocle),file=file.path('mm',paste0(prefix,'-',this.sample,'-expr.mm')))

	fwrite(list(rownames(exprs(cds.monocle))),sep='\n',file=file.path('mm',paste0(prefix,'-',this.sample,'-expr.rows.txt.gz')))
	fwrite(list(paste(this.sample,colnames(exprs(cds.monocle)),sep='-')),sep='\n',file=file.path('mm',paste0(prefix,'-',this.sample,'-expr.cols.txt.gz')))
}
meta = read.delim(file.path('data',paste0(prefix,'_metadata.txt')))

out = as.data.frame(colData(cds.monocle))
out$id = this.sample
out = merge(out,meta,by='id')
rownames(out) = with(out,paste(id,cell,sep='-'))
out = out[paste(this.sample,colnames(exprs(cds.monocle)),sep='-'),]

fwrite(out,file=file.path('mm',paste0(prefix,'-',this.sample,'-metadata.txt.gz')),sep='\t')
fwrite(as.data.frame(rowData(cds.monocle)),file=file.path('mm',paste0(prefix,'-',this.sample,'-features.txt.gz')),sep='\t')
