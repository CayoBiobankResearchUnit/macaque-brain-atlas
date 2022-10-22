work_path = '~/macaque/'

library(RColorBrewer)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(reshape2)
library(scales)
library(corrplot)
library(ggrastr)

pd_all = readRDS(paste0(work_path, 'data/biccn-metadata-patched.rds'))

cell.subtypes = structure(list(cell_class=structure(c(1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,2L,2L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,3L,4L,4L,5L,5L,6L,6L,6L,6L,6L,6L,6L,6L,6L,6L,6L,6L,6L,6L,6L,6L,7L,8L,8L,8L,8L,8L,8L,9L,9L,9L,9L,9L,9L,9L,9L,10L,11L,11L,11L,11L,11L,11L,12L,12L,13L,14L,15L,16L,17L),
                                                    .Label=c('excitatory neurons','medium spiny neurons','inhibitory neurons','dopaminergic neurons','serotonergic neurons','cerebellar neurons','basket cells','astrocytes','oligodendrocytes','oligodendrocyte precursor cells','vascular cells','microglia','ependymal cells','AHSG neurons','F5 neurons','KIR3DL12 neurons','KIR3DL12 microglia'),class='factor'),
                               cell_subtype=structure(c(1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,13L,14L,15L,16L,17L,18L,19L,20L,21L,22L,23L,24L,25L,26L,27L,28L,29L,30L,31L,32L,33L,34L,35L,36L,37L,38L,39L,40L,41L,42L,43L,44L,45L,46L,47L,48L,49L,50L,51L,52L,53L,54L,55L,56L,57L,58L,59L,60L,61L,62L,63L,64L,65L,66L,67L,68L,69L,70L,71L,72L,73L,74L,75L,76L,77L,78L,79L,80L,81L,82L,83L,84L,85L,86L,87L,88L,89L,90L,91L,92L,93L,94L,95L,96L,97L,98L,99L,100L,101L,102L,103L,104L,105L,108L,109L,110L,111L,112L),
                                                      .Label=c('excitatory neurons 1','excitatory neurons 2','excitatory neurons 3','excitatory neurons 4','excitatory neurons 5','excitatory neurons 6','excitatory neurons 7','excitatory neurons 8','excitatory neurons 9','excitatory neurons 10','excitatory neurons 11','excitatory neurons 12','excitatory neurons 13','excitatory neurons 14','excitatory neurons 15','excitatory neurons 16','excitatory neurons 17','excitatory neurons 18','excitatory neurons 19','excitatory neurons 20','excitatory neurons 21','excitatory neurons 22','excitatory neurons 23','excitatory neurons 24','excitatory neurons 25','excitatory neurons 26','excitatory neurons 27','excitatory neurons 28','excitatory neurons 29','excitatory neurons 30','excitatory neurons 31','excitatory neurons 32','excitatory neurons 33','excitatory neurons 34','excitatory neurons 35','excitatory neurons 36','excitatory neurons 37','excitatory neurons 38','excitatory neurons 39','medium spiny neurons 1','medium spiny neurons 2','inhibitory neurons 1','inhibitory neurons 2','inhibitory neurons 3','inhibitory neurons 4','inhibitory neurons 5','inhibitory neurons 6','inhibitory neurons 7','inhibitory neurons 8','inhibitory neurons 9','inhibitory neurons 10','inhibitory neurons 11','inhibitory neurons 12','inhibitory neurons 13','inhibitory neurons 14','inhibitory neurons 15','inhibitory neurons 16','inhibitory neurons 17','inhibitory neurons 18','inhibitory neurons 19','inhibitory neurons 20','dopaminergic neurons 1','dopaminergic neurons 2','serotonergic neurons 1','serotonergic neurons 2','cerebellar neurons 1','cerebellar neurons 2','cerebellar neurons 3','cerebellar neurons 4','cerebellar neurons 5','cerebellar neurons 6','cerebellar neurons 7','cerebellar neurons 8','cerebellar neurons 9','cerebellar neurons 10','cerebellar neurons 11','cerebellar neurons 12','cerebellar neurons 13','cerebellar neurons 14','cerebellar neurons 15','cerebellar neurons 16','basket cells','astrocytes 1','astrocytes 2','astrocytes 3','astrocytes 4','astrocytes 5','astrocytes 6','oligodendrocytes 1','oligodendrocytes 2','oligodendrocytes 3','oligodendrocytes 4','oligodendrocytes 5','oligodendrocytes 6','oligodendrocytes 7','oligodendrocytes 8','oligodendrocyte precursor cells','vascular cells 1','vascular cells 2','vascular cells 3','vascular cells 4','vascular cells 5','vascular cells 6','microglia 1','microglia 2','radial glial cells','mesenchymal stem cells','ependymal cells','AHSG neurons','F5 neurons','KIR3DL12 neurons','KIR3DL12 microglia','serotinergic neurons 1','serotinergic neurons 2'),class='factor'),
                               color=c('#9aaff9','#2b77ad','#6ba4e0','#d6bf3b','#59ffd8','#163175','#f2e0a7','#96bcdd','#4e04a3','#e8985f','#2eb232','#e5647c','#b1f945','#7723a8','#3eefe3','#77d1e5','#a6f975','#f4e99f','#cce4ff','#b1b8ef','#88f73d','#144bcc','#8deeef','#ef73ad','#f9c2db','#ccb20c','#e542cf','#2b86bf','#7562c4','#9495e0','#abe24a','#59eaba','#46a4c4','#cf8bf9','#809302','#1eb24d','#62e0dc','#efe773','#c0fc92','#f791bd','#8b9bf4','#9bfffa','#8f6dce','#4d9901','#66a3cc','#e26f31','#51e8dd','#a5ef8d','#27f468','#67e573','#ffcce7','#9bd352','#dd82c8','#4242d6','#3ac9a8','#b7f28a','#34d822','#9bff84','#e87896','#2c57ba','#d83227','#f4c529','#e8853a','#aa56e2','#fccabf','#7097ea','#1bf922','#3ce093','#840110','#b172cc','#50f4cb','#7d85f2','#d3b5fc','#9cfcb2','#f9f7a7','#abf791','#51db8d','#98e21f','#c9af0a','#dd0063','#ca74ed','#fc499f','#21c44f','#b5345d','#b73e12','#6287ef','#f4dd6b','#1e65d8','#9de8f2','#e3f7a3','#e835d9','#ea56b1','#55f446','#48ce39','#1be8cc','#7c35a5','#f40c04','#f28aab','#b0caf4','#601d91','#225689','#78c6cc','#e028b8','#9680e5','#f4f7a3','#a3f1f7','#e07081','#95cadb','#95d635','#57d34a')),row.names=c(NA,-110L),class='data.frame')

pd = pd_all[pd_all$modality=='RNA',]

pd$cell_class = as.character(pd$rna_cell_class)
pd = pd[!is.na(pd$cell_class),]

colors = read.table(paste0(work_path, 'data/colors.txt'),header=F,sep='\t',comment.char = "")
color_plate = as.vector(colors$V2)
names(color_plate) = colors$V1

p1 = ggplot(sub_pd, aes(x=rna_sub_umap.1, y=rna_sub_umap.2, color = region_subclass)) + 
  geom_point(alpha=0.1) + scale_color_manual(values=color_plate, limits = force) + theme_classic() 
ggsave(rasterize(p1, dpi=300),filename=paste0(work_path, "subclustering/rna_inh_sub_umap_region_subclass.pdf"),width=15,height=12)


region_order = unique(pd[,c("region", "region_subclass", "region_class")]) %>%
  as.data.frame() %>%
  arrange(region_class, region_subclass, region)

cts = c("excitatory neurons","cerebellar neurons","inhibitory neurons","basket cells","medium spiny neurons",
        "dopaminergic neurons","serotonergic neurons","AHSG neurons","F5 neurons","KIR3DL12 neurons",
        "astrocytes","oligodendrocytes","oligodendrocyte precursor cells","ependymal cells",
        "microglia","KIR3DL12 microglia","vascular cells")

color_plate = as.vector(colors$V2)
names(color_plate) = colors$V1

##
pd$cell_class = factor(pd$cell_class,levels = cts)
pd$region = factor(pd$region,levels = region_order$region)
pd$region_subclass = factor(pd$region_subclass,levels = unique(region_order$region_subclass))
p1 = ggplot(pd, aes(x=rna_umap.1, y=rna_umap.2, color = region)) + 
  geom_point(alpha=0.05) + scale_color_manual(values=color_plate, limits = force) + theme_classic() 
p1
ggsave(rasterize(p1, dpi=300),filename=paste0(work_path, "plots/rna_umap_region_new.pdf"),width=15,height=12)

p2 = ggplot(pd, aes(x=rna_umap.1, y=rna_umap.2, color = cell_class)) + 
  geom_point(alpha=0.05) + scale_color_manual(limits = force, values=color_plate) + theme_classic() 
p2
ggsave(rasterize(p2, dpi=300),filename=paste0(work_path, "plots/rna_umap_cell_class_new.pdf"),width=15,height=12)

##
data_type = 'rna'
pd$cell_class = factor(pd$cell_class,levels = rev(cts))
pd$region = factor(pd$region,levels = rev(region_order$region))
pd$region_subclass = factor(pd$region_subclass,levels = rev(unique(region_order$region_subclass)))

generate_pct_plots(pd,data_type)

generate_pct_plots = function(pd,data_type) {
  ######
  #1d: ct precentages
  ######
  df1 = pd %>%
    group_by(cell_class,region) %>%
    tally() %>%
    rename(cell_num = n)
  tmp = pd %>% 
    group_by(cell_class) %>%
    tally() %>%
    rename(total_cell_num = n) %>%
    mutate(log_total_cell_num = log2(total_cell_num))
  df2 = df1 %>%
    left_join(tmp, by = "cell_class") %>%
    mutate(frac = 100*cell_num/total_cell_num) %>%
    select(-c(total_cell_num))
  df2$region = factor(df2$region,levels = region_order$region)
  p3 = ggplot(data=df2, aes(x=cell_class, y=frac, fill=region)) +
    geom_bar(position="stack", stat="identity", width = 0.6) +
    scale_fill_manual(values=color_plate, limits = force) +
    labs(x = "", y = "cell fraction (%)") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    theme(axis.text.x = element_text(angle = 60, hjust=1)) + coord_flip()
  ggsave(rasterize(p3, dpi=300),filename=paste0(work_path, "plots/",data_type,"_ct_pct_region_new.pdf"),width=8,height=12)
  
  df1 = pd %>%
    group_by(cell_class,region_subclass) %>%
    tally() %>%
    rename(cell_num = n)
  tmp = pd %>% 
    group_by(cell_class) %>%
    tally() %>%
    rename(total_cell_num = n) %>%
    mutate(log_total_cell_num = log2(total_cell_num))
  df2 = df1 %>%
    left_join(tmp, by = "cell_class") %>%
    mutate(frac = 100*cell_num/total_cell_num) %>%
    select(-c(total_cell_num))
  df2$region_subclass = factor(df2$region_subclass,levels = unique(region_order$region_subclass))
  p3 = ggplot(data=df2, aes(x=cell_class, y=frac, fill=region_subclass)) +
    geom_bar(position="stack", stat="identity", width = 0.6) +
    scale_fill_manual(values=color_plate, limits = force) +
    labs(x = "", y = "cell fraction (%)") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    theme(axis.text.x = element_text(angle = 60, hjust=1)) + coord_flip()
  ggsave(rasterize(p3, dpi=300),filename=paste0(work_path, "plots/",data_type,"_ct_pct_region_subclass_new.pdf"),width=8,height=12)
}

##specificity
library(philentropy)
cell.type.column = 'cell_class'
region.column = 'region'
m = as.matrix(table(pd[[cell.type.column]],pd[[region.column]]))
specificity.out = data.frame(cell_type = levels(pd[[cell.type.column]]))
specificity.out$specificity_score = unlist(lapply(1:nrow(m),function(i) {
  out = suppressMessages(JSD(rbind(m[i,],colSums(m)),est.prob='empirical'))
  names(out) = rownames(m)[i]
  out
}))

specificity.out$cell_type = factor(specificity.out$cell_type, levels = rev(cts))
p = ggplot(data=specificity.out, aes(x=cell_type, y=specificity_score)) + #, fill=cell_type)) +
  geom_bar(position="stack", stat="identity", width = 0.6) +
  scale_fill_manual(values=color_plate, limits = force) +
  labs(x = "", y = "regional specificity score") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1)) + coord_flip()
write.table(specificity.out,paste0(work_path, "data/rna_ct_region_specificity.txt"),row.names = F,quote = F,sep=',')
ggsave(rasterize(p, dpi=300),filename=paste0(work_path, "plots/rna_ct_region_specificity_new.pdf"),width=6,height=12)

##counts
tmp = pd %>% 
  group_by(cell_class) %>%
  tally() %>%
  rename(total_cell_num = n) %>%
  mutate(log_total_cell_num = log2(total_cell_num))

tmp$cell_class = factor(tmp$cell_class,levels = rev(cts))
p = ggplot(data=tmp, aes(x=cell_class, y=log_total_cell_num, fill=cell_class)) +
  geom_bar(position="stack", stat="identity", width = 0.6) +
  scale_fill_manual(values=color_plate, limits = force) +
  labs(x = "", y = "log2(cell count)") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) + coord_flip()
ggsave(rasterize(p, dpi=300),filename=paste0(work_path, "plots/",data_type,"_ct_logcount.pdf"),width=8,height=12)

p = ggplot(data=tmp, aes(x=cell_class, y=total_cell_num, fill=cell_class)) +
  geom_bar(position="stack", stat="identity", width = 0.6) +
  scale_fill_manual(values=color_plate, limits = force) +
  labs(x = "", y = "cell count") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) + coord_flip()
ggsave(rasterize(p, dpi=300),filename=paste0(work_path, "plots/",data_type,"_ct_count.pdf"),width=8,height=12)

sub_color_map = cell.subtypes$color
names(sub_color_map) = cell.subtypes$cell_subtype
sub_pd = pd[pd$cell_class=='inhibitory neurons',]
sub_pd$cell_subcluster_manual = gsub('inhibitory neurons','',sub_pd$rna_cell_subcluster_manual)

p2 = ggplot(sub_pd, aes(x=rna_sub_umap.1, y=rna_sub_umap.2, color = rna_cell_subcluster_manual)) + 
  geom_point(alpha=0.1) + scale_color_manual(values=sub_color_map,limits=force) + theme_classic() 
p2
ggsave(rasterize(p2, dpi=300),filename=paste0(work_path, "subclustering/rna_inh_sub_umap_cell_class.pdf"),width=20,height=12)

