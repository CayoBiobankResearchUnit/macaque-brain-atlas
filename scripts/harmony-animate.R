#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('atac','atac','2')

prefix = arguments[1]
analysis = arguments[2]
n.iter = as.integer(arguments[3])

if (length(arguments) < 4) {
	color.by = 'cluster'
} else {
	color.by = arguments[4]
}

suppressMessages(library(Matrix))
suppressMessages(library(data.table))

sample.list = as.character(read.delim(file.path('data',paste0(prefix,'_metadata.txt')))$id)

this.sample = 'all'

# Read in the "good" UMAP (mostly for metadata and cluster annotations)
this.umap = readRDS(file.path('umap',paste0(prefix,'-scanpy-all.rds')))

# Read in progressive stages of UMAPs

umap.files = file.path('stats/harmony',c(
	paste0('umap_',c('prenorm','preharm'),'_',prefix,'_n',1,'_all.txt.gz'),
	paste0('umap_',c('postharm'),'_',prefix,'_n',1:n.iter,'_all.txt.gz')
))

library(parallel)
umaps = mclapply(
	1:length(umap.files),
	function(i) {
		out = fread(umap.files[i])
		names(out) = paste('umap',1:2,sep='.')
		out = data.frame(this.umap[,setdiff(colnames(this.umap),paste('umap',1:2,sep='.'))],out)
		out$iter = i
		out
},mc.cores=length(umap.files))

umaps.all = do.call(rbind,umaps)

plot.umap = function(x,method='umap',color='sample',color.label=NULL,log.transform=FALSE,file=NULL,width=7,height=7,rasterize=TRUE,legend=TRUE,facet=FALSE,facet.by=NULL,size=0.1,alpha=0.05) {
	require(ggplot2)
	require(ggrastr)
	require(RColorBrewer)
	require(viridis)
	require(randomcoloR)

	y = x

	color.label = if (is.null(color.label)) color else color.label

	set.seed(42)
	scale.color = if (class(y[[color]]) %in% c('integer','numeric')) {
		scale_color_viridis(option='D',name=color.label,trans=if(log.transform) 'log10' else 'identity')
	} else if (class(y[[color]]) %in% c('character','factor','logical')) {
		if (color=='landmark') {
			scale_color_manual(values=c('#cccccc','#ff0000'),name=color.label)
		} else if (length(unique(y[[color]])) <= 8) {
			scale_color_brewer(palette='Dark2',name=color.label)
		} else {
			# scale_color_manual(
			# 	name=color.label,
			# 	values=scan(file.path('palettes',paste0('randomcoloR-umap-',formatC(length(unique(y[[color]])),width=4,flag=0),'.txt')),quiet=TRUE,what='',sep='\t')
			# )
			scale_color_manual(name=color.label,values=distinctColorPalette(k=length(unique(y[[color]]))))
			# scale_color_discrete(name=color.label)
		}
	}

	geom.point = if (rasterize) {
		geom_point_rast(size=size,shape=19,alpha=alpha,dev='ragg_png')
	} else {
		geom_point(size=size,shape=19,alpha=alpha)
	}
	
	p = ggplot(data=y,aes_string('umap.1','umap.2',color=color)) +
		geom.point +
		coord_equal() +
		scale.color +
		theme_classic() +
		theme(
			axis.text=element_blank(),
			axis.ticks=element_blank(),
			legend.position=if (legend) 'right' else 'none'
		) +
		guides(color = if (class(y[[color]]) %in% c('integer','numeric')) guide_colorbar() else guide_legend(override.aes = list(size=2,alpha=1))) +
		xlab(paste(ifelse(method=='umap','UMAP',ifelse(method=='tsne','t-SNE','Dim')),1)) +
		ylab(paste(ifelse(method=='umap','UMAP',ifelse(method=='tsne','t-SNE','Dim')),2))

	if (facet) p = p + facet_wrap(ifelse(is.null(facet.by),color,facet.by))

	if (!is.null(file)) ggsave(p,file=file,useDingbats=FALSE,width=width,height=height) else p
}

stages = c('Pre-Normalization','Post-Normalization',paste('Harmony Iteration',1:n.iter))

umaps.all$iter = factor(stages[umaps.all$iter],levels=stages)

library(gganimate)
# assignInNamespace('draw_frames',draw_frames.parallel,'gganimate')

dir.create('animations',showWarnings=FALSE)

p = plot.umap(umaps.all,color=color.by,file=NULL,legend=FALSE) +
	ggtitle('{closest_state}') +
	transition_states(iter, transition_length = 4, state_length = 1) +
    ease_aes('cubic-in-out')
push.status('plot.umap')

a = animate(
    plot = p,
#	renderer = gifski_renderer(loop = FALSE),
    nframes = 200, width = 800, height = 640, res = 300
)
push.status('animate')

anim_save(
	file.path('animations',paste0('umap-harmony-',prefix,'-',gsub('_','',color.by),'-all_nolegend.gif')),
	a
)
push.status('anim_save')

