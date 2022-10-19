#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)

arguments='rna'

prefix = arguments[1]

library(monocle3)

cds = readRDS(file.path('checkpoints',paste0(prefix,'_cds.rds')))

cds$Region = gsub('([A-Z]+)_([a-z]+)','\\2\\1',cds$Region)
cds$Region[cds$Region == 'Vrm'] = 'CV'
cds$Region = factor(cds$Region,levels=region.levels)

# Plot manual cell type assignments
p = plot_cells(cds, color_cells_by='CellType_1B', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),alpha=0.05,size=0.1,stroke=0.1) +
	coord_equal() +
	scale_color_discrete(name='Cell type',drop=TRUE) +
	guides(color = guide_legend(override.aes = list(size=2,alpha=1))) +
	theme_classic(base_size=16) +
	theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',prefix,'_clustering_manual_cell_type.pdf'),useDingbats=FALSE,height=7)
ggsave(p+theme(legend.position='none'),file=paste0('figures/',prefix,'_clustering_manual_cell_type_nolegend.pdf'),useDingbats=FALSE,height=7)

# Plot manual cell type assignments
p = plot_cells(cds, color_cells_by='Region', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),alpha=0.05,size=0.1,stroke=0.1) +
	coord_equal() +
	scale_color_discrete(name='Region',drop=TRUE) +
	guides(color = guide_legend(override.aes = list(size=2,alpha=1))) +
	theme_classic(base_size=16) +
	theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',prefix,'_clustering_manual_region.pdf'),useDingbats=FALSE,height=7)
ggsave(p+theme(legend.position='none'),file=paste0('figures/',prefix,'_clustering_manual_region_nolegend.pdf'),useDingbats=FALSE,height=7)

# Plot manual cell type assignments
p = plot_cells(cds, color_cells_by='Region_Major', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),alpha=0.05,size=0.1,stroke=0.1) +
	coord_equal() +
	scale_color_discrete(name='Region (class)',drop=TRUE) +
	guides(color = guide_legend(override.aes = list(size=2,alpha=1))) +
	theme_classic(base_size=16) +
	theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',prefix,'_clustering_manual_regionmajor.pdf'),useDingbats=FALSE,height=7)
ggsave(p+theme(legend.position='none'),file=paste0('figures/',prefix,'_clustering_manual_regionmajor_nolegend.pdf'),useDingbats=FALSE,height=7)


# Plot manual cell type assignments
p = plot_cells(cds, color_cells_by='sex', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),alpha=0.05,size=0.1,stroke=0.1) +
	coord_equal() +
	scale_color_discrete(name='Sex',drop=TRUE) +
	guides(color = guide_legend(override.aes = list(size=2,alpha=1))) +
	theme_classic(base_size=16) +
	theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',prefix,'_clustering_manual_sex.pdf'),useDingbats=FALSE,height=7)
ggsave(p+theme(legend.position='none'),file=paste0('figures/',prefix,'_clustering_manual_sex_nolegend.pdf'),useDingbats=FALSE,height=7)



# Plot manual cell type assignments
p = plot_cells(cds, color_cells_by='hemisphere', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),alpha=0.05,size=0.1,stroke=0.1) +
	coord_equal() +
	scale_color_discrete(name='Hemisphere',drop=TRUE) +
	guides(color = guide_legend(override.aes = list(size=2,alpha=1))) +
	theme_classic(base_size=16) +
	theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',prefix,'_clustering_manual_hemisphere.pdf'),useDingbats=FALSE,height=7)
ggsave(p+theme(legend.position='none'),file=paste0('figures/',prefix,'_clustering_manual_hemisphere_nolegend.pdf'),useDingbats=FALSE,height=7)

# Save dimension reduction

cds.meta = colData(cds)
cds.umap = as.data.frame(reducedDims(cds)$UMAP)
names(cds.umap) = c('dim1','dim2')
cds.umap$cluster = clusters(cds)
cds.umap$partition = partitions(cds)

cds.out = data.frame(cds.meta[rownames(cds.umap),],cds.umap)

saveRDS(cds.out,file=file.path('checkpoints',paste0(prefix,'_umap.rds')))