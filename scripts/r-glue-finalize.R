#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('biccn','rna','atac')

prefix = arguments[1]
rn.prefix = arguments[2]
at.prefix = arguments[3]

rn.umap = readRDS(file.path('umap',paste0(rn.prefix,'-scanpy-recluster-classified.rds')))
at.umap = readRDS(file.path('umap',paste0(prefix,'-',at.prefix,'-glue-subintegration-classified.rds')))

rownames(at.umap) = at.umap$cell

glue.umap = readRDS(file.path('umap',paste0(prefix,'-glue-final.rds')))

glue.umap = subset(within(glue.umap,{
	glue_umap.1 = umap.1
	glue_umap.2 = umap.2
	glue_partition = partition
	glue_cluster = cluster
}),select=c('glue_umap.1','glue_umap.2','glue_partition','glue_cluster'))

# Assign the overall cell class 
at.umap$cell_class = at.umap$cell_class_prediction
at.umap$cell_class[at.umap$cell_class_prediction_conf < 0.95] = NA
at.umap$cell_subtype = at.umap$cell_subcluster_prediction

rn.umap$cell_subtype = rn.umap$cell_subcluster_manual

meta.columns = c('id','cell','biccn_id','region','region_label','region_class','region_subclass','hemisphere','animal_id','social_group','age','sex','extractor','isolation_site','isolation_order','extraction_batch','extraction.date','fixation','lysis','sequencing_run_id','cell_class','cell_subtype')

glue.meta = rbind(
	data.frame(rn.umap[,meta.columns],modality=factor('RNA',levels=c('RNA','ATAC'))),
	data.frame(at.umap[,meta.columns],modality=factor('ATAC',levels=c('RNA','ATAC')))
)

glue.meta$biccn_id = gsub('^A','',glue.meta$biccn_id)

# Crunch assay specific columns
rn.meta = subset(within(rn.umap,{
	rna_umi = n.umi
	rna_mt_umi = perc_mitochondrial_umis
	rna_doublet_score = doublet_score
	rna_umap.1 = umap.1
	rna_umap.2 = umap.2
	rna_partition = partition
	rna_cluster = cluster
	rna_cell_round1_type = cell_round1_type
	rna_cell_class = cell_class
	rna_subcluster_hi_res = subcluster_hi_res
	rna_cell_subcluster_hi_res = cell_subcluster_hi_res
	rna_subcluster_lo_res = subcluster_lo_res
	rna_cell_subcluster_lo_res = cell_subcluster_lo_res
	rna_subcluster_manual = subcluster_manual
	rna_cell_subcluster_manual = cell_subcluster_manual
}),select=c('rna_umi','rna_mt_umi','rna_doublet_score','rna_umap.1','rna_umap.2','rna_partition','rna_cluster','rna_cell_round1_type','rna_cell_class','rna_subcluster_hi_res','rna_cell_subcluster_hi_res','rna_subcluster_lo_res','rna_cell_subcluster_lo_res','rna_subcluster_manual','rna_cell_subcluster_manual'))

at.meta = subset(within(at.umap,{
	atac_umi = n.umi
	atac_doublet_score = doublet_score
	atac_total_fragments = total_fragments
	atac_total_dedup_fragments = total_dedup_fragments
	atac_total_binarized_fragments = total_binarized_fragments
	atac_RIP = RIP
	atac_FRIP = FRIP
	atac_FRIT = FRIT
	atac_umap.1 = umap.1
	atac_umap.2 = umap.2
	atac_partition = partition
	atac_cluster = cluster
	atac_cell_class_manual = cell_class_manual
	atac_cell_class_prediction = cell_class_prediction
	atac_cell_class_prediction_conf = cell_class_prediction_conf
	atac_cell_subcluster_prediction = cell_subcluster_prediction
	atac_cell_subcluster_prediction_conf = cell_subcluster_prediction_conf
}),select=c('atac_umi','atac_doublet_score','atac_total_fragments','atac_total_dedup_fragments','atac_total_binarized_fragments','atac_RIP','atac_FRIP','atac_FRIT','atac_umap.1','atac_umap.2','atac_partition','atac_cluster','atac_cell_class_manual','atac_cell_class_prediction','atac_cell_class_prediction_conf','atac_cell_subcluster_prediction','atac_cell_subcluster_prediction_conf'))

# Incorporate cell class-specific umap coordinates

# file.path('umap',paste0(prefix,'-scanpy-recluster-class',this.cluster,'.rds'))
# file.path('umap',paste0(prefix,'-glue-recluster-class',this.cluster,'.rds'))

rn.sub.umap = do.call(rbind,lapply(seq_along(levels(glue.meta$cell_class)),function(i) {
	this.umap = readRDS(file.path('umap',paste0(rn.prefix,'-scanpy-recluster-class',i,'.rds')))[,c('umap.1','umap.2')]
	names(this.umap) = paste0('rna_sub_',names(this.umap))
	this.umap
}))

glue.sub.umap = do.call(rbind,lapply(seq_along(levels(glue.meta$cell_class))[file.exists(file.path('umap',paste0(prefix,'-glue-recluster-class',seq_along(levels(glue.meta$cell_class)),'.rds')))],function(i) {
	this.glue.umap = readRDS(file.path('umap',paste0(prefix,'-glue-recluster-class',i,'.rds')))[,c('cell','umap.1','umap.2','partition','cluster')]
	rownames(this.glue.umap) = this.glue.umap$cell
	this.glue.umap = within(this.glue.umap,{
		glue_sub_umap.1 = umap.1
		glue_sub_umap.2 = umap.2
		glue_subcluster = paste0(i,'-',cluster)
	})[,c('glue_sub_umap.1','glue_sub_umap.2','glue_subcluster')]
	this.glue.umap
}))

glue.umap$cell = rownames(glue.umap)
rn.meta$cell = rownames(rn.meta)
at.meta$cell = rownames(at.meta)
rn.sub.umap$cell = rownames(rn.sub.umap)
glue.sub.umap$cell = rownames(glue.sub.umap)

this.umap = merge(glue.meta,glue.umap,by='cell',all=TRUE)
this.umap = merge(this.umap,rn.meta,by='cell',all=TRUE)
this.umap = merge(this.umap,at.meta,by='cell',all=TRUE)
this.umap = merge(this.umap,rn.sub.umap,by='cell',all=TRUE)
this.umap = merge(this.umap,glue.sub.umap,by='cell',all=TRUE)

this.umap$biccn_id = factor(this.umap$biccn_id)

subcluster.levels = matrix(as.integer(do.call(rbind,strsplit(unique(this.umap$glue_subcluster[!is.na(this.umap$glue_subcluster)]),'-'))),ncol=2)
this.umap$glue_subcluster = factor(this.umap$glue_subcluster,levels=apply(subcluster.levels[order(subcluster.levels[,1],subcluster.levels[,2]),],1,function(x) paste(x[1],x[2],sep='-')))

rownames(this.umap) = this.umap$cell

saveRDS(this.umap,file=file.path('umap',paste0(prefix,'-glue-final-all-metadata.rds')))