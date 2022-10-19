#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
rna.prefix = arguments[1]
atac.prefix = arguments[2]
n.cells = if (length(arguments) < 3) '' else as.numeric(arguments[3])

library(Seurat)
library(monocle3)
library(Signac)

cds.rna = readRDS(file.path('checkpoints',paste0(rna.prefix,'_cds_seurat.rds')))
cds.atac = readRDS(file.path('checkpoints',paste0(atac.prefix,'_',n.cells,'_cds_seurat.rds')))

gene.activities = GeneActivity(cds.atac, features = VariableFeatures(cds.rna))

push.status(paste('GeneActivity',n.cells))

cds.atac[['ACTIVITY']] = CreateAssayObject(counts = gene.activities)

DefaultAssay(cds.atac) = 'ACTIVITY'
cds.atac = NormalizeData(cds.atac)

push.status(paste('NormalizeData',n.cells))

cds.atac = ScaleData(cds.atac, features = rownames(cds.atac))

push.status(paste('ScaleData',n.cells))

transfer.anchors = FindTransferAnchors(reference = cds.rna, query = cds.atac, features = VariableFeatures(object = cds.rna),
    reference.assay = 'RNA', query.assay = 'ACTIVITY', reduction = 'cca')

push.status(paste('FindTransferAnchors',n.cells))

celltype.predictions = TransferData(anchorset = transfer.anchors, refdata = cds.rna$seurat_annotations,
    weight.reduction = cds.atac[['lsi']], dims = 1:75)

push.status(paste('TransferData',n.cells))

cds.atac = AddMetaData(cds.atac, metadata = celltype.predictions)

# Co-embed

push.status(paste('Label transfer',n.cells))

genes.use = VariableFeatures(cds.rna)

push.status(paste('VariableFeatures',n.cells))

refdata = GetAssayData(cds.rna, assay = 'RNA', slot = 'data')[genes.use, ]

imputation = TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = cds.atac[['lsi']],
    dims = 1:75)

push.status(paste('TransferData',n.cells))

cds.atac[['RNA']] = imputation

cds.cmbd = merge(x = cds.rna, y = cds.atac)

cds.cmbd = ScaleData(cds.cmbd, features = genes.use, do.scale = FALSE)

push.status(paste('ScaleData',n.cells))

cds.cmbd = RunPCA(cds.cmbd, features = genes.use, verbose = FALSE)

push.status(paste('RunPCA',n.cells))

cds.cmbd = RunUMAP(cds.cmbd, dims = 1:30)

push.status(paste('RunUMAP',n.cells))

push.status(paste('Co-embedding',n.cells))

if (nchar(n.cells)) {
	save(list=ls(),file=file.path('checkpoints',paste0('coembed_',rna.prefix,'-',atac.prefix,'_',n.cells,'.RData')))
} else {
	save(list=ls(),file=file.path('checkpoints',paste0('coembed_',rna.prefix,'-',atac.prefix,'.RData')))
}
