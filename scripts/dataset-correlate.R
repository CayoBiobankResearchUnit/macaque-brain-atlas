#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')
source('scripts/_include_nnls.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('rna','dev-inhibitory-neurons-macaque')

prefix = arguments[1]
dataset_prefix = arguments[2]

if (length(arguments) > 2) {
	n.cores = as.integer(arguments[3])
}

suppressMessages(library(data.table))

# Determine species

organism = if (dataset_prefix %in% c('cortex-dev','human-cortex','dev-brain-regions','adult-brain-vasc-endo','adult-brain-vasc-peri','vascular-dev','early-brain','brain-vasc-atlas','myeloid-neuroinflam','fang-merfish-l2','fang-merfish-l3')) {
	'human'
} else if (dataset_prefix %in% c('mouse-drg-injury')) {
	'mouse'
} else if (dataset_prefix %in% c('dev-inhibitory-neurons-macaque')) {
	'macaque'
}

# first matrix
expr.file = file.path('datasets/processed',paste0(dataset_prefix,'-pseudobulk-expr.txt.gz'))
rows.file = file.path('datasets/processed',paste0(dataset_prefix,'-pseudobulk-features.txt.gz'))
cols.file = file.path('datasets/processed',paste0(dataset_prefix,'-pseudobulk-meta.txt.gz'))

e = as.matrix(fread(expr.file,header=FALSE))
these.genes = rownames(e) = fread(rows.file,header=TRUE)$V1
meta = fread(cols.file,header=TRUE)
these.cells = colnames(e) = meta$V1
meta$V1 = NULL

ref.e = e
ref.m = meta

# second matrix
expr.file = file.path('datasets/processed',paste0(prefix,'-pseudobulk-expr.txt.gz'))
rows.file = file.path('datasets/processed',paste0(prefix,'-pseudobulk-features.txt.gz'))
cols.file = file.path('datasets/processed',paste0(prefix,'-pseudobulk-meta.txt.gz'))

e = as.matrix(fread(expr.file,header=FALSE))
these.genes = rownames(e) = fread(rows.file,header=TRUE)$V1
meta = fread(cols.file,header=TRUE)
these.cells = colnames(e) = meta$V1
meta$V1 = NULL

qry.e = e
qry.m = meta

if (organism == 'human') {

	if (file.exists('rds/mmul_hsap_orthologs.rds')) {
		mmul.orthologs = readRDS('rds/mmul_hsap_orthologs.rds')
	} else {
		suppressMessages(library(biomaRt))
		mmul = useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL', dataset='mmulatta_gene_ensembl',version=101)
		mmul.orthologs = getBM(attributes=c('ensembl_gene_id','external_gene_name','hsapiens_homolog_ensembl_gene','hsapiens_homolog_orthology_type','hsapiens_homolog_associated_gene_name'),mart=mmul)
	}
	mmul.hsap.orthologs = subset(mmul.orthologs,ensembl_gene_id %in% rownames(qry.e) & hsapiens_homolog_associated_gene_name %in% rownames(ref.e) & hsapiens_homolog_orthology_type == 'ortholog_one2one')

	qry.e = qry.e[mmul.hsap.orthologs$ensembl_gene_id,]
	ref.e = ref.e[mmul.hsap.orthologs$hsapiens_homolog_associated_gene_name,]

} else if (organism == 'mouse') {

	if (file.exists('rds/mmul_mmus_orthologs.rds')) {
		mmul.orthologs = readRDS('rds/mmul_mmus_orthologs.rds')
	} else {
		suppressMessages(library(biomaRt))
		mmul = useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL', dataset='mmulatta_gene_ensembl',version=101)
		mmul.orthologs = getBM(attributes=c('ensembl_gene_id','external_gene_name','mmusculus_homolog_ensembl_gene','mmusculus_homolog_orthology_type','mmusculus_homolog_associated_gene_name'),mart=mmul)
	}
	mmul.mmus.orthologs = subset(mmul.orthologs,ensembl_gene_id %in% rownames(qry.e) & mmusculus_homolog_associated_gene_name %in% rownames(ref.e) & mmusculus_homolog_orthology_type == 'ortholog_one2one')

	qry.e = qry.e[mmul.mmus.orthologs$ensembl_gene_id,]
	ref.e = ref.e[mmul.mmus.orthologs$mmusculus_homolog_associated_gene_name,]

} else if (organism == 'macaque') {

	if (file.exists('rds/mmul_mmul_orthologs.rds')) {
		mmul.orthologs = readRDS('rds/mmul_mmul_orthologs.rds')
	} else {
		suppressMessages(library(biomaRt))
		mmul = useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL', dataset='mmulatta_gene_ensembl',version=101)
		mmul.orthologs = getBM(attributes=c('ensembl_gene_id','external_gene_name'),mart=mmul)
	}
	mmul.mmul.orthologs = subset(mmul.orthologs,ensembl_gene_id %in% rownames(qry.e) & external_gene_name %in% rownames(ref.e))

	qry.e = qry.e[mmul.mmul.orthologs$ensembl_gene_id,]
	ref.e = ref.e[mmul.mmul.orthologs$external_gene_name,]

}
rownames(ref.e) = rownames(qry.e)

# Accidental empty column
ref.e = ref.e[,colnames(ref.e) != '']

# Run Jun's NNLS code
cor.result = correlation_analysis_bidirection(qry.e,ref.e)

cor.result = cor.result[order(cor.result$source,-cor.result$beta_1),]

dir.create('rds/nnls',showWarnings=FALSE)
saveRDS(cor.result,file=file.path('rds/nnls',paste0('nnls_results-',prefix,'_',dataset_prefix,'.rds')))