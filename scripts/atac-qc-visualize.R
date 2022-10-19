#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')

arguments = commandArgs(trailing=TRUE)
prefix = arguments[1]

library(monocle3)
library(ggplot2)
library(ggrastr)
library(viridis)
library(egg)

cds = readRDS(file.path('checkpoints',paste0(prefix,'_cds.rds')))

# UMAP colored by partition
p = plot_cells(cds, color_cells_by='partition', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),size=0.2,alpha=0.02) +
	coord_equal() +
 	guides(color = guide_legend(title='Partition',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',prefix,'_clustering_partitions.pdf'),useDingbats=FALSE)

# UMAP colored by cluster
p = plot_cells(cds, color_cells_by='cluster', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),size=0.2,alpha=0.02) +
	coord_equal() +
 	guides(color = guide_legend(title='Cluster',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',prefix,'_clustering_clusters.pdf'),useDingbats=FALSE)

# UMAP colored by region
p = plot_cells(cds, color_cells_by='region', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),size=0.2,alpha=0.02) +
	coord_equal() +
 	guides(color = guide_legend(title='Region',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',prefix,'_clustering_regions.pdf'),useDingbats=FALSE)

# UMAP colored by n.umi
p = plot_cells(cds, color_cells_by='n.umi', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),size=0.2,alpha=0.02) +
	coord_equal() +
 	scale_color_viridis(option='D',name='UMI',trans='log10') +
	guides(color = guide_colorbar()) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',prefix,'_clustering_umi.pdf'),useDingbats=FALSE)

# UMAP colored by total
p = plot_cells(cds, color_cells_by='total', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),size=0.2,alpha=0.02) +
	coord_equal() +
 	scale_color_viridis(option='D',name='Depth',trans='log10') +
	guides(color = guide_colorbar()) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',prefix,'_clustering_total.pdf'),useDingbats=FALSE)

# UMAP colored by frip
p = plot_cells(cds, color_cells_by='FRIP', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),size=0.2,alpha=0.02) +
	coord_equal() +
	scale_color_viridis(option='D',name='FRIP') +
	guides(color = guide_colorbar()) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',prefix,'_clustering_frip.pdf'),useDingbats=FALSE)

# UMAP colored by frit
p = plot_cells(cds, color_cells_by='FRIT', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),size=0.2,alpha=0.02) +
	coord_equal() +
	scale_color_viridis(option='D',name='FRIT',trans='log10') +
	guides(color = guide_colorbar()) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',prefix,'_clustering_frit.pdf'),useDingbats=FALSE)

# Per-sample distribution of FRIP
p = ggplot(as.data.frame(colData(cds)),aes(FRIP)) +
	geom_density() +
	facet_wrap(~id,ncol=ceiling(sqrt(length(unique(colData(cds)$id))))) +
	scale_x_continuous(trans='log10') +
	theme_minimal(base_size=16) +
	theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank()) +
	xlab('FRIP')
ggsave(p,file=paste0('figures/',prefix,'_qc_frip.pdf'),useDingbats=FALSE)

# Per-sample distribution of FRIT
p = ggplot(as.data.frame(colData(cds)),aes(FRIT)) +
	geom_density() +
	facet_wrap(~id,ncol=ceiling(sqrt(length(unique(colData(cds)$id))))) +
	scale_x_continuous(trans='log10') +
	theme_minimal(base_size=16) +
	theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank()) +
	xlab('FRIT')
ggsave(p,file=paste0('figures/',prefix,'_qc_frit.pdf'),useDingbats=FALSE)

# Per-sample distribution of binarized UMI after filtering
p = ggplot(as.data.frame(colData(cds)),aes(n.umi)) +
	geom_density() +
	facet_wrap(~id,ncol=ceiling(sqrt(length(unique(colData(cds)$id))))) +
	scale_x_continuous(trans='log10') +
	theme_minimal(base_size=16) +
	theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank()) +
	xlab('UMI')
ggsave(p,file=paste0('figures/',prefix,'_qc_umi.pdf'),useDingbats=FALSE)

# Per-sample distribution of total before binarizing and filtering
p = ggplot(as.data.frame(colData(cds)),aes(total)) +
	geom_density() +
	facet_wrap(~id,ncol=ceiling(sqrt(length(unique(colData(cds)$id))))) +
	scale_x_continuous(trans='log10') +
	theme_minimal(base_size=16) +
	theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank()) +
	xlab('Depth')
ggsave(p,file=paste0('figures/',prefix,'_qc_total.pdf'),useDingbats=FALSE)

# Global distribution of FRIP
p = ggplot(as.data.frame(colData(cds)),aes(FRIP)) +
	geom_density() +
	scale_x_continuous(trans='log10') +
	theme_minimal(base_size=16) +
	theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank()) +
	xlab('FRIP')
ggsave(p,file=paste0('figures/',prefix,'_qc_global_frip.pdf'),useDingbats=FALSE)

# Global distribution of FRIT
p = ggplot(as.data.frame(colData(cds)),aes(FRIT)) +
	geom_density() +
	scale_x_continuous(trans='log10') +
	theme_minimal(base_size=16) +
	theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank()) +
	xlab('FRIT')
ggsave(p,file=paste0('figures/',prefix,'_qc_global_frit.pdf'),useDingbats=FALSE)

# Global distribution of binarized UMI after filtering
p = ggplot(as.data.frame(colData(cds)),aes(n.umi)) +
	geom_density() +
	scale_x_continuous(trans='log10') +
	theme_minimal(base_size=16) +
	theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank()) +
	xlab('UMI')
ggsave(p,file=paste0('figures/',prefix,'_qc_global_umi.pdf'),useDingbats=FALSE)

# Global distribution of depth before filtering
p = ggplot(as.data.frame(colData(cds)),aes(total)) +
	geom_density() +
	scale_x_continuous(trans='log10') +
	theme_minimal(base_size=16) +
	theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank()) +
	xlab('Depth')
ggsave(p,file=paste0('figures/',prefix,'_qc_global_total.pdf'),useDingbats=FALSE)

# Per-sample depth vs FRIP
p = ggplot(as.data.frame(colData(cds)),aes(total,FRIP)) +
#	stat_density_2d(aes(fill=..level..),geom='polygon',colour='black',size=0.2) +
	geom_point_rast(size=0.005,alpha=0.01) +
#	geom_density_2d(size=0.25,color='#0000ff') +
	facet_wrap(~id,ncol=ceiling(sqrt(length(unique(colData(cds)$id))))) +
	scale_x_continuous(trans='log10') +
	scale_y_continuous() +
	scale_fill_distiller(palette=4, direction=1) +
	theme_minimal(base_size=16) +
	coord_flip() +
	theme() +
	xlab('Depth') + ylab('FRIP')
ggsave(p,file=paste0('figures/',prefix,'_qc_total_v_frip.pdf'),useDingbats=FALSE)

# Per-sample depth vs FRIP
p = ggplot(as.data.frame(colData(cds)),aes(total,FRIP)) +
	stat_density_2d(aes(fill=..level..),geom='polygon',colour='black',size=0.2) +
#	geom_point_rast(size=0.005,alpha=0.001) +
#	geom_density_2d(size=0.25,color='#0000ff') +
	facet_wrap(~id,ncol=ceiling(sqrt(length(unique(colData(cds)$id))))) +
	scale_x_continuous(trans='log10') +
	scale_y_continuous() +
#	scale_fill_viridis(option='D',direction=-1,name='Level') +
	scale_fill_distiller(palette=4, direction=1,name='Level') +
	theme_minimal(base_size=16) +
	coord_flip() +
	theme() +
	xlab('Depth') + ylab('FRIP')
ggsave(p,file=paste0('figures/',prefix,'_qc_total_v_frip_density.pdf'),useDingbats=FALSE)

# Per-sample FRIP vs FRIT
p = ggplot(as.data.frame(colData(cds)),aes(FRIP,FRIT)) +
	geom_point_rast(size=0.01,alpha=0.05) +
	facet_wrap(~id,ncol=ceiling(sqrt(length(unique(colData(cds)$id))))) +
	scale_x_continuous(trans='log10') +
	scale_y_continuous(trans='log10') +
	theme_minimal(base_size=16) +
	theme()
	xlab('FRIP') + ylab('FRIT')
ggsave(p,file=paste0('figures/',prefix,'_qc_frip_v_frit.pdf'),useDingbats=FALSE)

# Global depth vs FRIP
p = ggplot(as.data.frame(colData(cds)),aes(total,FRIP)) +
#	stat_density_2d(aes(fill=..level..),geom='polygon',colour='black',size=0.2) +
	geom_point_rast(size=0.005,alpha=0.01) +
#	geom_density_2d(size=0.25,color='#0000ff') +
	scale_x_continuous(trans='log10') +
	scale_y_continuous() +
#	scale_fill_viridis(option='D',direction=-1,name='Level') +
	scale_fill_distiller(palette=4, direction=1,name='Level') +
	theme_minimal(base_size=16) +
	coord_flip() +
	theme() +
	xlab('Depth') + ylab('FRIP')
ggsave(p,file=paste0('figures/',prefix,'_qc_global_total_v_frip.pdf'),useDingbats=FALSE)

# Global depth vs FRIP
p = ggplot(as.data.frame(colData(cds)),aes(total,FRIP)) +
	stat_density_2d(aes(fill=..level..),geom='polygon',colour='black',size=0.2) +
#	geom_point_rast(size=0.005,alpha=0.001) +
#	geom_density_2d(size=0.25,color='#0000ff') +
	scale_x_continuous(trans='log10') +
	scale_y_continuous() +
#	scale_fill_viridis(option='D',direction=-1,name='Level') +
	scale_fill_distiller(palette=4, direction=1,name='Level') +
	theme_minimal(base_size=16) +
	coord_flip() +
	theme() +
	xlab('Depth') + ylab('FRIP')
ggsave(p,file=paste0('figures/',prefix,'_qc_global_total_v_frip_density.pdf'),useDingbats=FALSE)

# Global UMI vs FRIP
p = ggplot(as.data.frame(colData(cds)),aes(n.umi,FRIP)) +
#	stat_density_2d(aes(fill=..level..),geom='polygon',colour='black',size=0.2) +
	geom_point_rast(size=0.005,alpha=0.01) +
#	geom_density_2d(size=0.25,color='#0000ff') +
	scale_x_continuous(trans='log10') +
	scale_y_continuous() +
#	scale_fill_viridis(option='D',direction=-1,name='Level') +
	scale_fill_distiller(palette=4, direction=1,name='Level') +
	theme_minimal(base_size=16) +
	coord_flip() +
	theme() +
	xlab('UMI') + ylab('FRIP')
ggsave(p,file=paste0('figures/',prefix,'_qc_global_umi_v_frip.pdf'),useDingbats=FALSE)

# Global UMI vs FRIP
p = ggplot(as.data.frame(colData(cds)),aes(n.umi,FRIP)) +
	stat_density_2d(aes(fill=..level..),geom='polygon',colour='black',size=0.2) +
#	geom_point_rast(size=0.005,alpha=0.001) +
#	geom_density_2d(size=0.25,color='#0000ff') +
	scale_x_continuous(trans='log10') +
	scale_y_continuous() +
#	scale_fill_viridis(option='D',direction=-1,name='Level') +
	scale_fill_distiller(palette=4, direction=1,name='Level') +
	theme_minimal(base_size=16) +
	coord_flip() +
	theme() +
	xlab('UMI') + ylab('FRIP')
ggsave(p,file=paste0('figures/',prefix,'_qc_global_umi_v_frip_density.pdf'),useDingbats=FALSE)
