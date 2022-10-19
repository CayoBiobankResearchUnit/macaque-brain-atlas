#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

# arguments = c('unique_peaks','6','16')

peak.set = arguments[1]
this.class = arguments[2]
this.cluster = arguments[3]

if (!file.exists(paste0('rds/monalisa/class',this.class,'_cluster',this.cluster,'_',peak.set,'_TFenrich_df.rds'))) {
	library(SummarizedExperiment)
	library(TFBSTools)
	library(JASPAR2018)
	library(BSgenome.Mmulatta.UCSC.rheMac10)
	library(monaLisa)

	if (!file.exists('rds/jaspar2018_pfms.rds')) {
		pwms = getMatrixSet(JASPAR2018,
			opts = list(matrixtype = 'PWM',
			tax_group = 'vertebrates'))
		saveRDS(pwms,file='rds/jaspar2018_pfms.rds')
	} else {
		pwms = readRDS('rds/jaspar2018_pfms.rds')
	}

	lmr = rtracklayer::import(
		con = paste0('bed/subsubpeaks/',peak.set,'-class',this.class,'_cluster',this.cluster,'.bed'),
		format = 'bed')

	seqlevelsStyle(lmr) = 'UCSC'

	lmrseqs = getSeq(BSgenome.Mmulatta.UCSC.rheMac10, lmr)

	se = calcBinnedMotifEnrR(seqs = lmrseqs, pwmL = pwms, BPPARAM = BiocParallel::MulticoreParam(20), verbose = T, background = 'genome' , genome=BSgenome.Mmulatta.UCSC.rheMac10)

	dir.create('rds/monalisa',showWarnings=FALSE)

	saveRDS(se,paste0('rds/monalisa/class',this.class,'_cluster',this.cluster,'_',peak.set,'_TFenrich.rds'))

	results.out = data.frame(
		motif = se@elementMetadata$motif.id,
		motif_name = se@elementMetadata$motif.name,
		negLog10P = se@assays@data$negLog10P[,1],
		negLog10Padj = se@assays@data$negLog10Padj[,1],
		log2enr = se@assays@data$log2enr[,1]
	)

	results.out = results.out[with(results.out,order(-negLog10Padj,-negLog10P,-log2enr)),]

	saveRDS(results.out,paste0('rds/monalisa/class',this.class,'_cluster',this.cluster,'_',peak.set,'_TFenrich_df.rds'))
}
# results.out$padj = 10^-results.out$negLog10Padj



