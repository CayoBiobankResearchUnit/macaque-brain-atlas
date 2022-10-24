work_path = '~/macaque/'

library(RColorBrewer)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(viridis)
library(gplots)
library(monocle3)

regions = c('NAc', 'LGN', 'ACC', 'STS', 'IPP', 'V1', 'CN', 'Pons', 'MdC',
            'HIP', 'S1', 'lCb', 'AMY', 'mdTN', 'dlPFC', 'dmPFC', 'vlPFC',
            'vmPFC', 'CC', 'vlTN', 'M1', 'A1', 'MdO', 'EC', 'PC', 'SPP', 'CV',
            'IT', 'MB')
region_subclasses = c('basal ganglia','thalamus','frontal lobe', 'temporal lobe', 'parietal lobe', 'occipital lobe',
                      'brainstem','hippocampus','cerebellum', 'amygdala', 'corpus callosum')

cts = c("excitatory neurons","medium spiny neurons","inhibitory neurons","dopaminergic neurons",
        "serotinergic neurons","cerebellar neurons","basket cells","astrocytes",
        "oligodendrocytes","oligodendrocyte precursor cells",
        "vascular cells","microglia","ependymal cells","AHSG neuron","F5 neuron","KIR3DL12 neuron","KIR3DL12 microglia")
##
pd_all = readRDS(paste0(work_path, 'data/biccn-metadata-patched.rds'))
pd = pd_all[pd_all$modality=='RNA',]
region_order = unique(pd[,c("region", "region_subclass", "region_class")]) %>%
  as.data.frame() %>%
  arrange(region_class, region_subclass, region)

colors = read.table(paste0(work_path, 'data/colors.txt'),header=F,sep='\t',comment.char = "")
color_plate = as.vector(colors$V2)
names(color_plate) = colors$V1

require(data.table)
require(plotly)
library(ggrastr)
library(viridis)
##region subclass new
ct = "astrocytes"
i = 8
df2 = fread(paste0(work_path, 'similarity_score/',ct,'_similarity_scores_subclass.txt'),header=FALSE)
cell_ids = read.table(paste0(work_path, '/data/knn_graphs/cellids_rna_class',i,'.txt'), header=F)
rownames(pd) = pd$cell
sub_pd = pd[cell_ids$V1,]
df2 = as.data.frame(df2)
rownames(df2) = cell_ids$V1
colnames(df2) = region_subclasses
df = df2
df$max_region = colnames(df2)[max.col(df2,ties.method="random")]
df$max_region = factor(df$max_region,levels = unique(region_order$region_subclass))
df$sum_values = rowSums(df2)
df$max_value = apply(df2, MARGIN =  1, FUN = max, na.rm = T)
df$max_pct = df$max_value/df$sum_values
df = cbind(sub_pd[,c("region","region_label","region_class","region_subclass","rna_umap.1","rna_umap.2",
                     "rna_sub_umap.1","rna_sub_umap.2","rna_cell_class","rna_cell_subcluster_manual")],df)
p1 = ggplot(df, aes(x=rna_sub_umap.1, y=rna_sub_umap.2, alpha = 0.01, color = max_region)) + 
  geom_point() + facet_wrap(~max_region) + scale_color_manual(guide="none", values=color_plate) + theme_classic() 
ggsave(rasterize(p1, dpi=300),filename = paste0(work_path, 'plots/',ct,'_lochness_multi_subclass_facet_new.pdf'), width=10, height=10)
p2 = ggplot(df, aes(x=rna_sub_umap.1, y=rna_sub_umap.2, alpha = 0.01, color = max_region)) + 
  geom_point() + scale_color_manual(guide="none", values=color_plate) + theme_classic() 
ggsave(rasterize(p2, dpi=300),filename = paste0(work_path, 'plots/',ct,'_lochness_multi_subclass_new.pdf'), width=10, height=10)

p3 = ggplot(df, aes(x=rna_sub_umap.1, y=rna_sub_umap.2, color = max_pct)) + 
  geom_point(size=0.1) + theme_classic() 
p3

temp = df[df[['thalamus']]>0,]
p3 = ggplot(temp, aes(x=rna_sub_umap.1, y=rna_sub_umap.2,color = `thalamus`,alpha=0.01)) + 
  geom_point(size=0.01) + scale_colour_gradient(low = "gray88", high = "blue", na.value = NA) + theme_classic() 
ggsave(rasterize(p3, dpi=300),filename = paste0(work_path, 'plots/',ct,'_lochness_thalamus.pdf'), width=10, height=10)

temp = df[df[['brainstem']]>0,]
p3 = ggplot(temp, aes(x=rna_sub_umap.1, y=rna_sub_umap.2,color = `brainstem`,alpha=0.01)) + 
  geom_point(size=0.01) + scale_colour_gradient(low = "gray88", high = "blue", na.value = NA) + theme_classic() 
ggsave(rasterize(p3, dpi=300),filename = paste0(work_path, 'plots/',ct,'_lochness_brainstem.pdf'), width=10, height=10)

temp = df[df[['cerebellum']]>0,]
p3 = ggplot(temp, aes(x=rna_sub_umap.1, y=rna_sub_umap.2,color = `cerebellum`,alpha=0.01)) + 
  geom_point(size=0.01) + scale_colour_gradient(low = "gray88", high = "blue", na.value = NA) + theme_classic() 
ggsave(rasterize(p3, dpi=300),filename = paste0(work_path, 'plots/',ct,'_lochness_cerebellum.pdf'), width=10, height=10)

temp = df[df[['occipital lobe']]>0,]
p3 = ggplot(temp, aes(x=rna_sub_umap.1, y=rna_sub_umap.2,color = `occipital lobe`,alpha=0.01)) + 
  geom_point(size=0.01) + scale_colour_gradient(low = "gray88", high = "blue", na.value = NA) + theme_classic() 
ggsave(rasterize(p3, dpi=300),filename = paste0(work_path, 'plots/',ct,'_lochness_occipital lobe.pdf'), width=10, height=10)

temp = df[df[['basal ganglia']]>0,]
p3 = ggplot(temp, aes(x=rna_sub_umap.1, y=rna_sub_umap.2,color = `basal ganglia`,alpha=0.01)) + 
  geom_point(size=0.01) + scale_colour_gradient(low = "gray88", high = "blue", na.value = NA) + theme_classic() 
ggsave(rasterize(p3, dpi=300),filename = paste0(work_path, 'plots/',ct,'_lochness_basal ganglia.pdf'), width=10, height=10)

temp = df[df[['frontal lobe']]>0,]
p3 = ggplot(temp, aes(x=rna_sub_umap.1, y=rna_sub_umap.2,color = `frontal lobe`,alpha=0.01)) + 
  geom_point(size=0.01) + scale_colour_gradient(low = "gray88", high = "blue", na.value = NA) + theme_classic() 
ggsave(rasterize(p3, dpi=300),filename = paste0(work_path, 'plots/',ct,'_lochness_frontal lobe.pdf'), width=10, height=10)

##region new, just cortical
ct = "excitatory neurons"
i = 1
df2 = fread(paste0(work_path, 'similarity_score/',ct,'_similarity_scores.txt'),header=FALSE)
cell_ids = read.table(paste0(work_path, '/data/knn_graphs/cellids_rna_class',i,'.txt'), header=F)
rownames(pd) = pd$cell
sub_pd = pd[cell_ids$V1,]
df2 = as.data.frame(df2)
rownames(df2) = cell_ids$V1
colnames(df2) = regions[1:ncol(df2)]

df2 = df2[,c('dmPFC','vmPFC','dlPFC','vlPFC','ACC','M1','EC','PC','A1','STS','IT','S1','IPP','SPP','V1')]
df = df2

df$max_region = colnames(df2)[max.col(df2,ties.method="random")]
df$max_region = factor(df$max_region,levels = region_order$region)
df$sum_values = rowSums(df2)
df$max_value = apply(df2, MARGIN =  1, FUN = max, na.rm = T)
df$max_pct = df$max_value/df$sum_values

df = cbind(sub_pd[,c("region","region_label","region_class","region_subclass","rna_umap.1","rna_umap.2",
                     "rna_sub_umap.1","rna_sub_umap.2","rna_cell_class","rna_cell_subcluster_manual")],df)
p1 = ggplot(df, aes(x=rna_sub_umap.1, y=rna_sub_umap.2, alpha = 0.01, color = max_region)) + 
  geom_point() + facet_wrap(~max_region) + scale_color_manual(guide="none", values=color_plate) + theme_classic() 
ggsave(rasterize(p1, dpi=300),filename = paste0(work_path, 'plots/',ct,'_lochness_multi_facet_new_cortical.pdf'), width=10, height=10)
p2 = ggplot(df, aes(x=rna_sub_umap.1, y=rna_sub_umap.2, alpha = 0.01, color = max_region)) + 
  geom_point() + scale_color_manual(guide="none", values=color_plate) + theme_classic() 
ggsave(rasterize(p2, dpi=300),filename = paste0(work_path, 'plots/',ct,'_lochness_multi_new_cortical.pdf'), width=10, height=10)

for (i in c('dmPFC','vmPFC','dlPFC','vlPFC','ACC','M1','EC','PC','A1','STS','IT','S1','IPP','SPP','V1')) {
  temp = df[df[[i]]>0,]
  temp$value = temp[[i]]
  temp = temp[order(temp$value),]
  p3 = ggplot(temp, aes(x=rna_sub_umap.1, y=rna_sub_umap.2, color=value, alpha=0.01)) + 
    geom_point(size=0.01) + scale_colour_gradient(low = "gray88", high = "blue", na.value = NA) + theme_classic() 
  ggsave(rasterize(p3, dpi=300),filename = paste0(work_path, 'plots/ex_lochness_',i,'.pdf'), width=10, height=10)
}


ct = "oligodendrocytes"
i = 9
df2 = fread(paste0(work_path, 'similarity_score/',ct,'_similarity_scores.txt'),header=FALSE)
cell_ids = read.table(paste0(work_path, '/data/knn_graphs/cellids_rna_class',i,'.txt'), header=F)
rownames(pd) = pd$cell
sub_pd = pd[cell_ids$V1,]
df2 = as.data.frame(df2)
rownames(df2) = cell_ids$V1
colnames(df2) = regions[1:ncol(df2)]

df2 = df2[,c('CC','CN','NAc','AMY','HIP','mdTN','vlTN','LGN')]
df = df2

df$max_region = colnames(df2)[max.col(df2,ties.method="random")]
df$max_region = factor(df$max_region,levels = region_order$region)
df$sum_values = rowSums(df2)
df$max_value = apply(df2, MARGIN =  1, FUN = max, na.rm = T)
df$max_pct = df$max_value/df$sum_values

df = cbind(sub_pd[,c("region","region_label","region_class","region_subclass","rna_umap.1","rna_umap.2",
                     "rna_sub_umap.1","rna_sub_umap.2","rna_cell_class","rna_cell_subcluster_manual")],df)
p1 = ggplot(df, aes(x=rna_sub_umap.1, y=rna_sub_umap.2, alpha = 0.01, color = max_region)) + 
  geom_point() + facet_wrap(~max_region) + scale_color_manual(guide="none", values=color_plate) + theme_classic() 
ggsave(rasterize(p1, dpi=300),filename = paste0(work_path, 'plots/',ct,'_lochness_multi_facet_new_subcortical.pdf'), width=10, height=10)
p2 = ggplot(df, aes(x=rna_sub_umap.1, y=rna_sub_umap.2, alpha = 0.01, color = max_region)) + 
  geom_point() + scale_color_manual(guide="none", values=color_plate) + theme_classic() 
ggsave(rasterize(p2, dpi=300),filename = paste0(work_path, 'plots/',ct,'_lochness_multi_new_subcortical.pdf'), width=10, height=10)

for (i in c('CC','CN','NAc','AMY','HIP','mdTN','vlTN','LGN')) {
  temp = df[df[[i]]>0,]
  temp$value = temp[[i]]
  temp = temp[order(temp$value),]
  p3 = ggplot(temp, aes(x=rna_sub_umap.1, y=rna_sub_umap.2, color=value, alpha=0.01)) + 
    geom_point(size=0.01) + scale_colour_gradient(low = "gray88", high = "blue", na.value = NA) + theme_classic() 
  ggsave(rasterize(p3, dpi=300),filename = paste0(work_path, 'plots/oligo_lochness_',i,'.pdf'), width=10, height=10)
}

ct = "inhibitory neurons"
i = 3
df2 = fread(paste0(work_path, 'similarity_score/',ct,'_similarity_scores.txt'),header=FALSE)
cell_ids = read.table(paste0(work_path, '/data/knn_graphs/cellids_rna_class',i,'.txt'), header=F)
rownames(pd) = pd$cell
sub_pd = pd[cell_ids$V1,]
df2 = as.data.frame(df2)
rownames(df2) = cell_ids$V1
colnames(df2) = regions[1:ncol(df2)]

df2 = df2[,c('CC','CN','NAc','AMY','HIP','mdTN','vlTN','LGN')]
df = df2

df$max_region = colnames(df2)[max.col(df2,ties.method="random")]
df$max_region = factor(df$max_region,levels = region_order$region)
df$sum_values = rowSums(df2)
df$max_value = apply(df2, MARGIN =  1, FUN = max, na.rm = T)
df$max_pct = df$max_value/df$sum_values

df = cbind(sub_pd[,c("region","region_label","region_class","region_subclass","rna_umap.1","rna_umap.2",
                     "rna_sub_umap.1","rna_sub_umap.2","rna_cell_class","rna_cell_subcluster_manual")],df)
p1 = ggplot(df, aes(x=rna_sub_umap.1, y=rna_sub_umap.2, alpha = 0.01, color = max_region)) + 
  geom_point() + facet_wrap(~max_region) + scale_color_manual(guide="none", values=color_plate) + theme_classic() 
ggsave(rasterize(p1, dpi=300),filename = paste0(work_path, 'plots/',ct,'_lochness_multi_facet_new_subcortical.pdf'), width=10, height=10)
p2 = ggplot(df, aes(x=rna_sub_umap.1, y=rna_sub_umap.2, alpha = 0.01, color = max_region)) + 
  geom_point() + scale_color_manual(guide="none", values=color_plate) + theme_classic() 
ggsave(rasterize(p2, dpi=300),filename = paste0(work_path, 'plots/',ct,'_lochness_multi_new_subcortical.pdf'), width=10, height=10)

for (i in c('CC','CN','NAc','AMY','HIP','mdTN','vlTN','LGN')) {
  temp = df[df[[i]]>0,]
  temp$value = temp[[i]]
  temp = temp[order(temp$value),]
  p3 = ggplot(temp, aes(x=rna_sub_umap.1, y=rna_sub_umap.2, color=value, alpha=0.01)) + 
    geom_point(size=0.01) + scale_colour_gradient(low = "gray88", high = "blue", na.value = NA) + theme_classic() 
  ggsave(rasterize(p3, dpi=300),filename = paste0(work_path, 'plots/inh_lochness_',i,'.pdf'), width=10, height=10)
}
