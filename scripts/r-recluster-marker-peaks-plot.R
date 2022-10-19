#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('atacsub','atac','6')

prefix = arguments[1]
analysis = arguments[2]
this.cluster = as.integer(arguments[3])

suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(ggrastr))

cell.levels = scan(what='',sep='\n',file='stats/subclusters/rna-cellsubclusters-levels.txt')
cell.classes = scan(what='',sep='\n',file='stats/clusters/rna-final-cellclasses-levels.txt')

use.checkpoint = TRUE
make.figures = TRUE

checkpoint.file = file.path('rds',paste0(prefix,'_class',this.cluster,'_marker_peaks_combined.rds'))
if (!file.exists(checkpoint.file) || !use.checkpoint) {
	scanpy = readRDS(file.path('rds',paste0(prefix,'_class',this.cluster,'_marker_peaks.rds')))

	scanpy$cell = factor(scanpy$cell,levels=cell.levels)[,drop=TRUE]
	scanpy$ttest_score = as.numeric(scanpy$ttest_score)
	scanpy$logreg_score = as.numeric(scanpy$logreg_score)
	scanpy$logfoldchanges = as.numeric(scanpy$logfoldchanges)

	seurat = do.call(rbind,lapply(file.path('rds',list.files('rds',pattern=paste0(prefix,'_class',this.cluster,'_marker_peaks_seurat_subcluster_i[0-9]+.rds'))),readRDS))
	seurat$cluster = factor(seurat$cluster,levels=cell.levels)[,drop=TRUE]

	names(scanpy)[names(scanpy) == 'cell'] = 'cluster'

	methods.merged = merge(seurat,scanpy,by=c('cluster','peak'),all.x=TRUE)

	methods.merged.pass = droplevels(subset(methods.merged,cluster %in% names(which(tapply(methods.merged$logreg_score,methods.merged$cluster,function(x) !all(is.na(x))) & tapply(methods.merged$p_val,methods.merged$cluster,function(x) !all(is.na(x)))))))

	methods.merged.pass$cluster_label = factor(
		methods.merged.pass$cluster,
		levels=levels(methods.merged.pass$cluster),
		labels=gsub(paste0(cell.classes[this.cluster],' '),'clust. ',levels(methods.merged.pass$cluster))
	)

	methods.merged.pass = do.call(rbind,lapply(split(methods.merged.pass,methods.merged.pass$cluster),function(x) {
		x = x[order(x$p_val),]
		x$p_null = seq(1/(sum(!is.na(x$p_val))),1,1/(sum(!is.na(x$p_val))))
		x
	}))
	
	saveRDS(methods.merged.pass,file=checkpoint.file)
} else {
	methods.merged.pass = readRDS(checkpoint.file)
}

marker.expr.threshold = 0.05

if (!is.null(methods.merged.pass)) {

if (make.figures) {

dir.create('figures/markerpeaks',showWarnings=FALSE)

if (F) {
p = ggplot(methods.merged.pass,aes(log10(pts),log10(pct.1))) +
	geom_point_rast(size=1e-2,alpha=1) +
	geom_abline(slope=1,intercept=0,color='black') +
	geom_smooth(method=lm,color='blue') +
	coord_equal() +
	facet_wrap(~cluster_label) +
	theme_classic(base_size=16) +
	xlab(expression('log'[10]*'(% acc'[Scanpy]*')')) +
	ylab(expression('log'[10]*'(% acc'[Seurat]*')')) +
	ggtitle(cell.classes[this.cluster])
suppressMessages(ggsave(p,file=file.path('figures/markerpeaks',paste0(prefix,'_class',this.cluster,'_seuratPercExpr_v_scanpyPercExpr.pdf')),useDingbats=FALSE))

p = ggplot(methods.merged.pass,aes(log10(pts_rest),log10(pct.2))) +
	geom_point_rast(size=1e-2,alpha=1) +
	geom_abline(slope=1,intercept=0,color='black') +
	geom_smooth(method=lm,color='blue') +
	coord_equal() +
	facet_wrap(~cluster_label) +
	theme_classic(base_size=16) +
	xlab(expression('log'[10]*'(% rest'[Scanpy]*')')) +
	ylab(expression('log'[10]*'(% rest'[Seurat]*')')) +
	ggtitle(cell.classes[this.cluster])
suppressMessages(ggsave(p,file=file.path('figures/markerpeaks',paste0(prefix,'_class',this.cluster,'_seuratPercRest_v_scanpyPercRest.pdf')),useDingbats=FALSE))

p = ggplot() +
	geom_point_rast(data=methods.merged.pass,aes(-log10(p_val),abs(logreg_score),color=log10(pts)),size=1e-2,alpha=1) +
	geom_smooth(data=subset(methods.merged.pass,pct.1>=marker.expr.threshold),aes(-log10(p_val),abs(logreg_score)),method=lm,color='blue') +
	geom_smooth(data=subset(methods.merged.pass,TRUE),aes(-log10(p_val),abs(logreg_score)),method=lm,color='black',size=0.5) +
#	geom_abline(slope=1,intercept=0) +
#	coord_equal() +
	facet_wrap(~cluster_label,scales='free') +
	scale_color_viridis(option='D',name=expression('log'[10]*'(% accessible)')) +
	theme_classic(base_size=16) +
	theme(legend.position='top',axis.text=element_text(size=8),strip.text=element_text(size=10)) +
	xlab(expression('-log'[10]~italic(p)[Seurat])) +
	ylab(expression('|Coefficient'[Scanpy]*'|')) +
	ggtitle(cell.classes[this.cluster])
suppressMessages(ggsave(p,file=file.path('figures/markerpeaks',paste0(prefix,'_class',this.cluster,'_seuratP_v_scanpyAbsLR.pdf')),useDingbats=FALSE))

p = ggplot() +
	geom_point_rast(data=methods.merged.pass,aes(-log10(p_val),logreg_score,color=log10(pts)),size=1e-2,alpha=1) +
#	geom_abline(slope=1,intercept=0) +
#	coord_equal() +
	facet_wrap(~cluster_label,scales='free') +
	scale_color_viridis(option='D',name=expression('log'[10]*'(% accessible)')) +
	theme_classic(base_size=16) +
	theme(legend.position='top',axis.text=element_text(size=8),strip.text=element_text(size=10)) +
	xlab(expression('-log'[10]~italic(p)[Seurat])) +
	ylab(expression('Coefficient'[Scanpy])) +
	ggtitle(cell.classes[this.cluster])
suppressMessages(ggsave(p,file=file.path('figures/markerpeaks',paste0(prefix,'_class',this.cluster,'_seuratP_v_scanpyLR.pdf')),useDingbats=FALSE))

p = ggplot(subset(methods.merged.pass,pct.1>=marker.expr.threshold),aes(-log10(p_val),abs(logreg_score),color=log10(pct.1))) +
	geom_point_rast(size=1e-2,alpha=1) +
	geom_smooth(method=lm,color='blue') +
#	geom_abline(slope=1,intercept=0) +
#	coord_equal() +
	facet_wrap(~cluster_label,scales='free') +
	scale_color_viridis(option='D',name=expression('log'[10]*'(% accessible)')) +
	theme_classic(base_size=16) +
	theme(legend.position='top',axis.text=element_text(size=8),strip.text=element_text(size=10)) +
	xlab(expression('-log'[10]~italic(p)[Seurat])) +
	ylab(expression('|Coefficient'[Scanpy]*'|')) +
	ggtitle(cell.classes[this.cluster])
suppressMessages(ggsave(p,file=file.path('figures/markerpeaks',paste0(prefix,'_class',this.cluster,'_seuratP_v_scanpyAbsLR-filtered.pdf')),useDingbats=FALSE))

p = ggplot(subset(methods.merged.pass,pct.1>=marker.expr.threshold),aes(-log10(p_val),logreg_score,color=log10(pct.1))) +
	geom_point_rast(size=1e-2,alpha=1) +
#	geom_abline(slope=1,intercept=0) +
#	coord_equal() +
	facet_wrap(~cluster_label,scales='free') +
	scale_color_viridis(option='D',name=expression('log'[10]*'(% accessible)')) +
	theme_classic(base_size=16) +
	theme(legend.position='top',axis.text=element_text(size=8),strip.text=element_text(size=10)) +
	xlab(expression('-log'[10]~italic(p)[Seurat])) +
	ylab(expression('Coefficient'[Scanpy])) +
	ggtitle(cell.classes[this.cluster])
suppressMessages(ggsave(p,file=file.path('figures/markerpeaks',paste0(prefix,'_class',this.cluster,'_seuratP_v_scanpyLR-filtered.pdf')),useDingbats=FALSE))

p = ggplot() +
	geom_point_rast(data=methods.merged.pass,aes(avg_log2FC,logfoldchanges,color=log10(pts)),size=1e-2,alpha=1) +
	geom_smooth(data=subset(methods.merged.pass,pts>=marker.expr.threshold),aes(avg_log2FC,logfoldchanges),method=lm,color='blue') +
	geom_smooth(data=subset(methods.merged.pass,pct.1>=marker.expr.threshold),aes(avg_log2FC,logfoldchanges),method=lm,color='red') +
	geom_smooth(data=subset(methods.merged.pass,TRUE),aes(avg_log2FC,logfoldchanges),method=lm,color='black',size=0.5) +
	facet_wrap(~cluster_label,scales='free') +
	scale_color_viridis(option='D',name=expression('log'[10]*'(% accessible)')) +
	theme_classic(base_size=16) +
	theme(legend.position='top',axis.text=element_text(size=8),strip.text=element_text(size=10)) +
	xlab(expression('log'[2]*'FC'[Seurat])) +
	ylab(expression('log'[2]*'FC'[Scanpy])) +
	ggtitle(cell.classes[this.cluster])
suppressMessages(ggsave(p,file=file.path('figures/markerpeaks',paste0(prefix,'_class',this.cluster,'_seuratFC_v_scanpyFC.pdf')),useDingbats=FALSE))

p = ggplot(subset(methods.merged.pass,pct.1>=marker.expr.threshold),aes(avg_log2FC,logfoldchanges,color=log10(pct.1))) +
	geom_point_rast(size=1e-2,alpha=1) +
	geom_smooth(method=lm,color='blue') +
	facet_wrap(~cluster_label,scales='free') +
	scale_color_viridis(option='D',name=expression('log'[10]*'(% accessible)')) +
	theme_classic(base_size=16) +
	theme(legend.position='top',axis.text=element_text(size=8),strip.text=element_text(size=10)) +
	xlab(expression('log'[2]*'FC'[Seurat])) +
	ylab(expression('log'[2]*'FC'[Scanpy])) +
	ggtitle(cell.classes[this.cluster])
suppressMessages(ggsave(p,file=file.path('figures/markerpeaks',paste0(prefix,'_class',this.cluster,'_seuratFC_v_scanpyFC-filtered.pdf')),useDingbats=FALSE))

p = ggplot() +
	geom_point_rast(data=methods.merged.pass,aes(avg_log2FC,logreg_score,color=log10(pct.1)),size=1e-2,alpha=1) +
	geom_smooth(data=subset(methods.merged.pass,pct.1>=marker.expr.threshold),aes(avg_log2FC,logreg_score),method=lm,color='blue') +
	geom_smooth(data=subset(methods.merged.pass,TRUE),aes(avg_log2FC,logreg_score),method=lm,color='black',size=0.5) +
	facet_wrap(~cluster_label,scales='free') +
	scale_color_viridis(option='D',name=expression('log'[10]*'(% accessible)')) +
	theme_classic(base_size=16) +
	theme(legend.position='top',axis.text=element_text(size=8),strip.text=element_text(size=10)) +
	xlab(expression('log'[2]*'FC'[Seurat])) +
	ylab(expression('Coefficient'[Scanpy])) +
	ggtitle(cell.classes[this.cluster])
suppressMessages(ggsave(p,file=file.path('figures/markerpeaks',paste0(prefix,'_class',this.cluster,'_seuratFC_v_scanpyLR.pdf')),useDingbats=FALSE))

p = ggplot(subset(methods.merged.pass,pct.1>=marker.expr.threshold),aes(avg_log2FC,logreg_score,color=log10(pct.1))) +
	geom_point_rast(size=1e-2,alpha=1) +
	geom_smooth(method=lm,color='blue') +
	facet_wrap(~cluster_label,scales='free') +
	scale_color_viridis(option='D',name=expression('log'[10]*'(% accessible)')) +
	theme_classic(base_size=16) +
	theme(legend.position='top',axis.text=element_text(size=8),strip.text=element_text(size=10)) +
	xlab(expression('log'[2]*'FC'[Seurat])) +
	ylab(expression('Coefficient'[Scanpy])) +
	ggtitle(cell.classes[this.cluster])
suppressMessages(ggsave(p,file=file.path('figures/markerpeaks',paste0(prefix,'_class',this.cluster,'_seuratFC_v_scanpyLR-filtered.pdf')),useDingbats=FALSE))

p = ggplot(methods.merged.pass, aes(-log10(p_null),-log10(p_val))) +
	geom_point_rast(size=1e-2,alpha=1) +
	geom_abline(slope=1,intercept=0) +
	facet_wrap(~cluster_label,scales='free_y') +
#	coord_equal() +
	theme_classic(base_size=16) +
	theme(axis.text=element_text(size=8),strip.text=element_text(size=10)) +
	xlab(expression('-log'[10]~italic(p)[expected])) +
	ylab(expression('-log'[10]~italic(p)[observed])) +
	ggtitle(cell.classes[this.cluster])
suppressMessages(ggsave(p,file=file.path('figures/markerpeaks',paste0(prefix,'_class',this.cluster,'_seuratP-qqplot.pdf')),useDingbats=FALSE))
}

fdr.cutoff = 0.2

methods.merged.pass = do.call(rbind,lapply(split(methods.merged.pass,methods.merged.pass$cluster),function(x) {
	logreg.threshold = median(subset(x,pts>=marker.expr.threshold & avg_log2FC > 0 & p_val_adj < fdr.cutoff)$logreg_score)
	x$logreg_threshold = logreg.threshold
	x$sig = factor(
		with(x,ifelse(
			pts<marker.expr.threshold,
			'not accessible',
			ifelse(
				avg_log2FC > 0 & p_val_adj < fdr.cutoff,
				'marker peak (Seurat)',
				ifelse(!is.na(logreg.threshold) & logreg_score > logreg.threshold,
					'marker peak (Scanpy)',
					'not significant')
				)
			)
		),
		levels=c('marker peak (Seurat)','marker peak (Scanpy)','not significant','not accessible')
	)
	x
}))

p = ggplot() +
	geom_point_rast(data=methods.merged.pass,aes(-log10(p_val),logreg_score,color=sig),size=1e-2,alpha=0.5) +
	geom_hline(data=unique(methods.merged.pass[,c('cluster_label','logreg_threshold')]),aes(yintercept=logreg_threshold),size=0.5,linetype=3) +
	facet_wrap(~cluster_label,scales='free') +
	scale_color_manual(values=c('#e41a1c','#377eb8','#000000','#cccccc'),name='Significant') +
	theme_classic(base_size=16) +
	theme(
		legend.position='top',
		axis.text=element_text(size=8),
#		legend.title=element_text(size=12),
		legend.title=element_blank(),
		legend.text=element_text(size=8),
		strip.text=element_text(size=10)
	) +
	guides(color=guide_legend(override.aes = list(size=2,alpha=1),nrow=1)) + 
	xlab(expression('-log'[10]~italic(p)[Seurat])) +
	ylab(expression('Coefficient'[Scanpy])) +
	ggtitle(cell.classes[this.cluster])
suppressMessages(ggsave(p,file=file.path('figures/markerpeaks',paste0(prefix,'_class',this.cluster,'_marker_peaks.pdf')),useDingbats=FALSE))

dir.create('stats/markerpeaks',showWarnings=FALSE)
dir.create('stats/markerpeaks/bed',showWarnings=FALSE)

methods.merged.pass$cluster_number = factor(
	gsub(paste0(cell.classes[this.cluster],' '),'',methods.merged.pass$cluster),
	levels=sort(unique(as.integer(gsub(paste0(cell.classes[this.cluster],' '),'',methods.merged.pass$cluster))))
)

methods.merged.sig = subset(methods.merged.pass,sig %in% c('marker peak (Seurat)','marker peak (Scanpy)'))

methods.merged.sig = data.frame(
	as.data.frame.matrix(do.call(rbind,strsplit(methods.merged.sig$peak,'_'))),
	methods.merged.sig
)
methods.merged.sig = within(methods.merged.sig,{
	V1 = factor(V1,levels=c(1:20,'X','Y'))
	V2 = as.numeric(V2)
	V3 = as.numeric(V3)
})
methods.merged.sig = methods.merged.sig[with(methods.merged.sig,order(V1,V2,V3)),]
rownames(methods.merged.sig) = NULL

for (i in levels(methods.merged.sig$cluster_number)) {
	cat(i,'\n')
	this = subset(methods.merged.sig,cluster_number == i,select=paste0('V',1:3))
	if (!nrow(this)) {
		warning('No significant peaks')
	} else {
		write.table(
			this,
			file=file.path('stats/markerpeaks/bed',paste0(prefix,'_class',this.cluster,'_cluster',i,'_marker_peaks.bed')),
			sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
	}
}

methods.merged = data.frame(
	as.data.frame.matrix(do.call(rbind,strsplit(methods.merged.pass$peak,'_'))),
	methods.merged.pass
)

n.peaks = length(unique(methods.merged$peak))

# Export top 1% of peaks
for (i in levels(methods.merged$cluster_number)) {
	cat(i,'\n')
	this = subset(methods.merged,cluster_number == i & logreg_score > 0,select=c(paste0('V',1:3),'logreg_score'))
	this = this[order(this$logreg_score,decreasing=TRUE),]
	this.out = head(this,round(0.01 * n.peaks))
	this.out = this.out[with(this.out,order(V1,V2,V3)),paste0('V',1:3)]
	if (!nrow(this)) {
		warning('No significant peaks')
	} else {
		write.table(
			this,
			file=file.path('bed/subsubpeaks',paste0('top_1p_peaks','-class',this.cluster,'_cluster',i,'.bed')),
			sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
	}
}



}

} else {
	warning('Not enough cell clusters.')
}

