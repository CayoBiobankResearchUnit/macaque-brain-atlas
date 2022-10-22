work_path = '~/macaque/'

library(RColorBrewer)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(viridis)
library(gplots)
library(monocle3)
require(plotly)
library(ggrastr)

pd = readRDS(paste0(work_path, 'data/biccn-metadata-patched.rds'))
colors = read.table(paste0(work_path, 'data/colors.txt'),header=F,sep='\t',comment.char = "")
color_plate = as.vector(colors$V2)
names(color_plate) = colors$V1
cts = c("excitatory neurons","cerebellar neurons","inhibitory neurons","basket cells","medium spiny neurons",
        "dopaminergic neurons","serotonergic neurons","AHSG neurons","F5 neurons","KIR3DL12 neurons",
        "astrocytes","oligodendrocytes","oligodendrocyte precursor cells","ependymal cells",
        "microglia","KIR3DL12 microglia","vascular cells")

###astrocytes
ct = "astrocytes"
cds = readRDS(paste0(work_path, 'subclustering/rna_',ct,'.cds.rds'))
scores = read.csv2(paste0(work_path, 'similarity_score/',ct,'_similarity_scores_subclass.txt'),sep=',',header=F)
regions = read.csv2(paste0(work_path, 'similarity_score/',ct,'_similarity_scores_subclass_regions.txt'),sep=',',header=F)
cell_ids = read.csv2(paste0(work_path, 'data/knn_graphs/cellids_rna_class8.txt'),header=F)
rownames(scores) = cell_ids$V1
colnames(scores) = regions$V1
scores = scores %>% mutate_if(is.character,as.numeric)
df = scores
df$max_region = colnames(scores)[max.col(scores,ties.method="random")]
df$sum_values = rowSums(scores)
df$max_value = apply(scores, MARGIN =  1, FUN = max, na.rm = T)
df$max_pct = df$max_value/df$sum_values
regions$V1 = sub(" ", "_", regions$V1)
colnames(scores) = regions$V1
scores = scores %>% mutate_if(is.character,as.numeric)
df = scores
df$max_region = colnames(scores)[max.col(scores,ties.method="random")]
df$sum_values = rowSums(scores)
df$max_value = apply(scores, MARGIN =  1, FUN = max, na.rm = T)
df$max_pct = df$max_value/df$sum_values
sub_pd = pd[cell_ids$V1,]
sub_pd = cbind(sub_pd[,c("region_label","region_class","region_subclass","umap.1","umap.2",
                         "umap_sub.1","umap_sub.2","cell_class","cell_subcluster_manual")],df)
cds = new_cell_data_set(cds@assays@data$counts,cell_metadata = sub_pd, gene_metadata = rowData(cds))
gene_fits <- fit_models(cds, model_formula_str = "~ basal_ganglia + thalamus + frontal_lobe + 
                      temporal_lobe + parietal_lobe + occipital_lobe + brainstem + hippocampus + amygdala + corpus_callosum + cerebellum") ###TODO
fit_coefs <- coefficient_table(gene_fits)
fit_terms <- fit_coefs %>% filter(term  %in% c("basal_ganglia","thalamus","frontal_lobe","temporal_lobe","parietal_lobe",
                                               "occipital_lobe","brainstem","hippocampus","amygdala","corpus_callosum","cerebellum")) ###TODO
sig_fit_terms = fit_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, gene_id, term, q_value, estimate)
saveRDS(sig_fit_terms, paste0(work_path, 'similarity_score/',cts[i],'_lochness_terms.rds'))
write.table(sig_fit_terms, file = paste0(work_path, 'similarity_score/',cts[i],'_lochness_terms.txt'), quote = F, sep = "\t",
            row.names = F,col.names = TRUE)
