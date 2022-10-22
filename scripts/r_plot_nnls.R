library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(viridis)
library(gplots)
library(SeuratDisk)
library(loomR)
library(Matrix)
library(dplyr)

work_path = '/Users/fanny/Documents/Shendure/macaque/'
res = readRDS(paste0(work_path, 'nnls/nnls_results-rna_brain-vasc-atlas.rds'))
res = res %>% mutate(beta=2 * (beta_1 + 0.01) * (beta_2 + 0.01))
res = res %>%
  dplyr::select(source, target, beta) %>%
  tidyr::spread(key=source, value=beta) %>%
  as.data.frame()
rownames(res) = res$target
res = as.matrix(res[, -1])
res_to_plot = cbind(res[,colnames(res)[grep('vascular',colnames(res))]],
                res[,colnames(res)[grep('microglia ',colnames(res))]],
                res[,colnames(res)[grep('ependymal ',colnames(res))]])
pheatmap(res_to_plot, filename = paste0(work_path, 'nnls/nnls_results-rna_brain-vasc-atlas.pdf'), 
         cluster_rows=F, cluster_cols=F, color= colorRampPalette(brewer.pal(n = 7, name = "Blues"))(50),)
while (!is.null(dev.list()))  dev.off()

res = readRDS(paste0(work_path, 'nnls/nnls_results-rna_human-cortex.rds'))
res = res %>% mutate(beta=2 * (beta_1 + 0.01) * (beta_2 + 0.01))
res = res %>%
  dplyr::select(source, target, beta) %>%
  tidyr::spread(key=source, value=beta) %>%
  as.data.frame()
rownames(res) = res$target
res = as.matrix(res[, -1])
res_to_plot = res[rownames(res)[grep('Exc',rownames(res))],colnames(res)[grep('excitatory',colnames(res))]]
pheatmap(res_to_plot, filename = paste0(work_path, 'nnls/nnls_results-rna_human-cortex-Exc.pdf'), 
         cluster_rows=F, cluster_cols=F, color= colorRampPalette(brewer.pal(n = 7, name = "Blues"))(50),)
while (!is.null(dev.list()))  dev.off()

res_to_plot = res[rownames(res)[grep('Inh',rownames(res))],colnames(res)[grep('inhibitory',colnames(res))]]
pheatmap(res_to_plot, filename = paste0(work_path, 'nnls/nnls_results-rna_human-cortex-Inh.pdf'), 
         cluster_rows=F, cluster_cols=F, color= colorRampPalette(brewer.pal(n = 7, name = "Blues"))(50),)
while (!is.null(dev.list()))  dev.off()

res_to_plot = res[c(rownames(res)[grep('Astro',rownames(res))],
                        rownames(res)[grep('Oligo',rownames(res))],
                        rownames(res)[grep('Micro',rownames(res))],
                        rownames(res)[grep('OPC',rownames(res))]
                        ),
                  c(colnames(res)[grep('astrocyte',colnames(res))],
                        colnames(res)[grep('microglia ',colnames(res))],
                        colnames(res)[grep('oligodendrocyte',colnames(res))]
                        )]
pheatmap(res_to_plot, filename = paste0(work_path, 'nnls/nnls_results-rna_human-cortex-glia.pdf'), 
         cluster_rows=F, cluster_cols=F, color= colorRampPalette(brewer.pal(n = 7, name = "Blues"))(50),)
while (!is.null(dev.list()))  dev.off()


res = readRDS(paste0(work_path, 'nnls/rna_HIP_fm_hippo_207785_nnls.rds'))

pheatmap(res, filename = paste0(work_path, 'nnls/rna_HIP_fm_hippo_207785_nnls.pdf'), 
         cluster_rows=F, cluster_cols=F, color= colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100))
while (!is.null(dev.list()))  dev.off()


res = readRDS(paste0(work_path, 'nnls/rna_cortex_Allen_human_cortex_nnls.rds'))
pheatmap(res, filename = paste0(work_path, 'nnls/rna_cortex_Allen_human_cortex_nnls.pdf'), 
         cluster_rows=F, cluster_cols=F, color= colorRampPalette(brewer.pal(n = 7, name = "Oranges"))(100),
         width = 15,height=10)
while (!is.null(dev.list()))  dev.off()
temp = scale(res,center=F)
pheatmap(temp, filename = paste0(work_path, 'nnls/rna_cortex_Allen_human_cortex_nnls_norm.pdf'), 
         cluster_rows=F, cluster_cols=T, color= colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
         width = 15,height=10)
while (!is.null(dev.list()))  dev.off()

sub = res[c("L4 IT","L5 ET","L5/6 IT Car3","L5/6 NP","L6 CT","L6b"),grep('excitatory', colnames(res), value=TRUE)]
pheatmap(sub, filename = paste0(work_path, 'nnls/rna_cortex_Allen_human_cortex_nnls_sub.pdf'), 
         cluster_rows=F, cluster_cols=T, color= colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
         width = 15,height=10)
while (!is.null(dev.list()))  dev.off()
temp_sub = scale(sub,center=F)
pheatmap(temp_sub, filename = paste0(work_path, 'nnls/rna_cortex_Allen_human_cortex_nnls_sub_norm.pdf'), 
         cluster_rows=F, cluster_cols=T, color= colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
         width = 15,height=10)
while (!is.null(dev.list()))  dev.off()
res = readRDS(paste0(work_path, 'nnls/rna_cortex_class_Allen_human_cortex_nnls.rds'))
pheatmap(res, filename = paste0(work_path, 'nnls/rna_cortex_class_Allen_human_cortex_nnls.pdf'), 
         cluster_rows=F, cluster_cols=F, color= colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
         width = 15,height=10)
while (!is.null(dev.list()))  dev.off()

res = readRDS(paste0(work_path, 'nnls/rna_HIP_fm_hippo_207785_nnls.rds'))
pheatmap(res, filename = paste0(work_path, 'nnls/rna_HIP_fm_hippo_207785_nnls.pdf'), 
         cluster_rows=F, cluster_cols=F, color= colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100))
while (!is.null(dev.list()))  dev.off()

res = readRDS(paste0(work_path, 'nnls/rna_CB_class_BGI_nnls.rds'))
pheatmap(res, filename = paste0(work_path, 'nnls/rna_CB_class_BGI_nnls.pdf'), 
         cluster_rows=F, cluster_cols=F, color= colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100))
while (!is.null(dev.list()))  dev.off()

res = readRDS(paste0(work_path, 'nnls/rna_CB_BGI_nnls.rds'))
pheatmap(res, filename = paste0(work_path, 'nnls/rna_CB_BGI_nnls.pdf'), 
         cluster_rows=F, cluster_cols=F, color= colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100))
while (!is.null(dev.list()))  dev.off()

res = readRDS(paste0(work_path, 'nnls/rna_CB_class_BGI_nnls.rds'))
pheatmap(res, filename = paste0(work_path, 'nnls/rna_CB_class_BGI_nnls.pdf'), 
         cluster_rows=F, cluster_cols=F, color= colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100))
while (!is.null(dev.list()))  dev.off()
