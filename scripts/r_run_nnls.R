library(plyr)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(stringr)
library(Matrix)
library(nnls)
library(tidyr)
library(FNN)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(viridis)
library(gplots)

source(paste0(work_path, 'scripts/r_dataset_correlations.R'))

run_nnls_seurat = function(seurat_obj, external_seurat_obj, field1, field2, min_cells_per_group = 50, fold.change = 1.5, gene_num = 200) {
  common_genes = intersect(rownames(seurat_obj), rownames(external_seurat_obj))
  tmp =  table(seurat_obj@meta.data[colnames(seurat_obj), field1])
  cells1 = colnames(seurat_obj)[seurat_obj@meta.data[, field1] %in% names(tmp)[tmp>min_cells_per_group]]
  seurat_obj = seurat_obj[,cells1]
  if (length(common_genes)<100) { print('common genes < 100') }
  combined.seurat = CreateSeuratObject(counts = cbind(GetAssayData(seurat_obj)[common_genes, ], GetAssayData(external_seurat_obj)[common_genes, ]), min.cells = 1)
  combined.seurat <- NormalizeData(combined.seurat)
  
  data.use1 = GetAssayData(combined.seurat)[, colnames(seurat_obj)]
  data.use2 = GetAssayData(combined.seurat)[, colnames(external_seurat_obj)]
  
  data.use1.average_profiles = get_average_expression(data.use1, grouping_vector=seurat_obj@meta.data[colnames(data.use1), field1])
  data.use2.average_profiles = get_average_expression(data.use2, grouping_vector=external_seurat_obj@meta.data[colnames(data.use2), field2])
  
  # Calculate correlations using NNLS approach
  nnls_correlations = bidirectional_correlation(data.use1.average_profiles, data.use2.average_profiles, fold.change = fold.change, top_gene_num = gene_num, spec_gene_num = gene_num)
  nnls_correlations = nnls_correlations %>% mutate(beta=2 * (beta_1 + 0.01) * (beta_2 + 0.01))
  
  nnls_correlations = nnls_correlations %>%
    dplyr::select(source, target, beta) %>%
    tidyr::spread(key=source, value=beta) %>%
    as.data.frame()
  
  rownames(nnls_correlations) = nnls_correlations$target
  res = as.matrix(nnls_correlations[, -1])
  return(res)
}

###cortical
mtx = as(t(Matrix::readMM(paste0(work_path, 'data/rna_cortical.mtx'))),'dgCMatrix')
gene_meta = read.csv2(paste0(work_path, 'data/rna_cortical.gene.txt'),sep=',',row.names = 1)
cell_meta = read.csv2(paste0(work_path, 'data/rna_cortical.cell.txt'),sep=',',row.names = 1)   
rownames(gene_meta) = gene_meta$id
rownames(mtx) = gene_meta$gene_short_name
colnames(mtx) = cell_meta$cell
seurat_obj = CreateSeuratObject(mtx,meta.data=cell_meta)

external_mtx = as(as.matrix(read.csv(paste0(work_path, 'data/external/human-cortex/human-cortex-expr.tsv'), sep="\t", row.names=1,header=T)),'dgCMatrix')
external_cell_meta = read.csv2(paste0(work_path, 'data/external/human-cortex/human-cortex-meta.tsv'),sep='\t',row.names = 1)
rownames(external_cell_meta) = gsub('-','.',rownames(external_cell_meta))
external_cell_meta = external_cell_meta[colnames(external_mtx),]
external_seurat_obj = CreateSeuratObject(external_mtx,meta.data=external_cell_meta)
external_seurat_obj = NormalizeData(external_seurat_obj)

field1 = 'cell_subcluster_manual'
field2 = 'subclass_label'
res = run_nnls_seurat(seurat_obj, external_seurat_obj, field1, field2, min_cells_per_group=50)
saveRDS(res,paste0(work_path, 'nnls/rna_cortex_Allen_human_cortex_nnls.rds'))

field1 = 'cell_class'
field2 = 'subclass_label'
res = run_nnls_seurat(seurat_obj, external_seurat_obj, field1, field2, min_cells_per_group=50)
saveRDS(res,paste0(work_path, 'nnls/rna_cortex_class_Allen_human_cortex_nnls.rds'))

####hippo
mtx = as(t(Matrix::readMM(paste0(work_path, 'data/rna_HIP.mtx'))),'dgCMatrix')
gene_meta = read.csv2(paste0(work_path, 'data/rna_HIP.gene.txt'),sep=',',row.names = 1)
cell_meta = read.csv2(paste0(work_path, 'data/rna_HIP.cell.txt'),sep=',',row.names = 1)   
rownames(gene_meta) = gene_meta$id
rownames(mtx) = gene_meta$gene_short_name
colnames(mtx) = cell_meta$cell
seurat_obj = CreateSeuratObject(mtx,meta.data=cell_meta)

external_mtx = as(Matrix::readMM(paste0(work_path, 'data/external/hippocampus/fm_hippo_207785.mtx')),'dgCMatrix')
external_mtx = t(external_mtx)
external_gene_meta = read.csv2(paste0(work_path, 'data/external/hippocampus/fm_hippo_207785_gene.csv'),sep=',',row.names = 1)
external_cell_meta = read.csv2(paste0(work_path, 'data/external/hippocampus/fm_hippo_207785_barcode.csv'),sep=',',row.names = 1)
rownames(external_mtx) = rownames(external_gene_meta)
colnames(external_mtx) = rownames(external_cell_meta)
external_seurat_obj = CreateSeuratObject(external_mtx,meta.data=external_cell_meta)
external_seurat_obj = NormalizeData(external_seurat_obj)

field1 = 'cell_subcluster_manual'
field2 = 'label_20'
res = run_nnls_seurat(seurat_obj, external_seurat_obj, field1, field2)
saveRDS(res,paste0(work_path, 'nnls/rna_HIP_fm_hippo_207785_nnls.rds'))

field1 = 'cell_class'
field2 = 'label_20'
res = run_nnls_seurat(seurat_obj, external_seurat_obj, field1, field2)
saveRDS(res,paste0(work_path, 'nnls/rna_HIP_class_fm_hippo_207785_nnls.rds'))

###cb
mtx = as(t(Matrix::readMM(paste0(work_path, 'data/rna_cerebellum.mtx'))),'dgCMatrix')
gene_meta = read.csv2(paste0(work_path, 'data/rna_cerebellum.gene.txt'),sep=',',row.names = 1)
cell_meta = read.csv2(paste0(work_path, 'data/rna_cerebellum.cell.txt'),sep=',',row.names = 1)   
rownames(gene_meta) = gene_meta$id
rownames(mtx) = gene_meta$gene_short_name
colnames(mtx) = cell_meta$cell
seurat_obj = CreateSeuratObject(mtx,meta.data=cell_meta)

external_seurat_obj = readRDS(paste0(work_path, 'data/external/cerebellum/Seurat4.0_Cerebellum.rds'))
external_seurat_obj = NormalizeData(external_seurat_obj)

field1 = 'cell_subcluster_manual'
field2 = 'Celltype'
res = run_nnls_seurat(seurat_obj, external_seurat_obj, field1, field2, min_cells_per_group = 20)
pheatmap(res, cluster_rows=F, cluster_cols=F)
saveRDS(res,paste0(work_path, 'nnls/rna_CB_BGI_nnls.rds'))

field1 = 'cell_class'
field2 = 'Celltype'
res = run_nnls_seurat(seurat_obj, external_seurat_obj, field1, field2, min_cells_per_group = 20)
pheatmap(res, cluster_rows=F, cluster_cols=F)
saveRDS(res,paste0(work_path, 'nnls/rna_CB_class_BGI_nnls.rds'))
