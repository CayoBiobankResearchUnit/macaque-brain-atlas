#!/usr/bin/env Rscript
# 🐛

# Function to send myself push alerts when things take a long time (R version)

# user ID and token must be set (omitted below)
uid = ''
token = ''

push.status = function(topic=NULL,msg=NULL,user=uid,token=token) {
	suppressMessages(require(pushoverr))
	if (is.null(msg)) m.msg = paste0(paste(Sys.info()[c('user','nodename')],collapse='@'),': ',ifelse(is.null(topic),'Action',topic),' completed') else m.msg = msg
	subj = ifelse(is.null(topic),'r-pushover',paste0('r-pushover - ',topic))
	pushover(
		message=m.msg,
		title=subj,
		user=user,
		app=token
	)
}

# Plot umaps
plot.umap = function(x,method='umap',color='sample',color.label=NULL,log.transform=FALSE,file=NULL,width=7,height=7,rasterize=TRUE,legend=TRUE,facet=FALSE,facet.by=NULL,size=0.1,alpha=0.05,presentation=FALSE,dev='ragg_png',dpi=300) {
	require(ggplot2)
	require(ggrastr)
	require(RColorBrewer)
	require(viridis)
	require(randomcoloR)
	require(egg)

	y = x

	color.label = if (is.null(color.label)) color else color.label
	
	# y[[color]][y[[color]] < 0] = NA

	set.seed(42)
	scale.color = if (class(y[[color]]) %in% c('integer','numeric')) {
		scale_color_viridis(option='D',name=color.label,trans=if(log.transform) 'log10' else 'identity',na.value='grey80')
	} else if (class(y[[color]]) %in% c('character','factor','logical')) {
		if (color=='landmark') {
			scale_color_manual(values=c('#cccccc','#ff0000'),name=color.label,na.translate=FALSE,na.value=rgb(1,1,1,0))
		} else if (color=='sex') {
			scale_color_manual(values=sex.colors,name=color.label,na.translate=FALSE,na.value=rgb(1,1,1,0))
		} else if (color=='region') {
			scale_color_manual(values=region.colors[which(table(y[[color]]) > 0)],name=color.label,na.translate=FALSE,na.value=rgb(1,1,1,0))
		} else if (color=='cell_final') {
			scale_color_manual(values=cell.type.colors[which(table(y[[color]]) > 0)],name=color.label,na.translate=TRUE,na.value=rgb(0.9,0.9,0.9,0.25))
		} else if (color=='hemisphere') {
			scale_color_manual(values=hemisphere.colors,name=color.label,na.translate=FALSE,na.value=rgb(1,1,1,0))
		} else if (length(unique(y[[color]])) <= 8) {
			scale_color_brewer(palette='Dark2',name=color.label,na.translate=FALSE,na.value=rgb(1,1,1,0))
		} else {
			# scale_color_manual(
			# 	name=color.label,
			# 	values=scan(file.path('palettes',paste0('randomcoloR-umap-',formatC(length(unique(y[[color]])),width=4,flag=0),'.txt')),quiet=TRUE,what='',sep='\t')
			# )
			scale_color_manual(name=color.label,values=distinctColorPalette(k=length(unique(y[[color]]))),na.translate=FALSE,na.value=rgb(1,1,1,0))
			# scale_color_discrete(name=color.label)
		}
	}

	geom.point = if (rasterize) {
		geom_point_rast(size=size,shape=19,alpha=alpha,dev=dev,raster.dpi=dpi,na.rm=TRUE)
	} else {
		geom_point(size=size,shape=19,alpha=alpha,na.rm=TRUE)
	}

	if (sum(is.na(y[[color]]))) {
		p = ggplot(data=y,aes_string('umap.1','umap.2',color=color)) +
		if (rasterize) {
			geom_point_rast(data=y[is.na(y[[color]]),],aes_string('umap.1','umap.2'),size=size,shape=19,color='#cccccc',alpha=alpha/5,dev=dev,raster.dpi=dpi)
		} else {
			geom_point(data=y[is.na(y[[color]]),],aes_string('umap.1','umap.2'),size=size,shape=19,color='#cccccc',alpha=alpha/5)
		}
	} else {
		p = ggplot(data=y,aes_string('umap.1','umap.2',color=color))
	}
	
	p = p +
		geom.point +
		coord_equal() +
		scale.color +
		theme_void() +
		(if (presentation) theme_article(base_size=24) else theme_article()) +
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

# Source: https://github.com/thomasp85/gganimate/blob/11edd90389fa6b51c86dd2d41a25ca5e1231ac56/R/animate.R
# See also: https://github.com/thomasp85/gganimate/issues/78

draw_frames.parallel <- function(plot, frames, device, ref_frame, ...) {
  require(future.apply)
  require(progress)
  require(grid)
  require(grDevices)
  require(ggplot2)

  stream <- device == 'current'

  dims <- tryCatch(
    plot_dims(plot, ref_frame),
    error = function(e) {
      warning('Cannot get dimensions of plot table. Plot region might not be fixed', call. = FALSE)
      list(widths = NULL, heights = NULL)
    }
  )

  dir <- tempfile(pattern = '')
  dir.create(dir, showWarnings = FALSE)
  files <- file.path(dir, sprintf('gganim_plot%04d', seq_along(frames)))
  files <- switch(
    tolower(device),
    png = paste0(files, '.png'),
    jpg = ,
    jpeg = paste0(files, '.jpg'),
    tif = ,
    tiff = paste0(files, '.tif'),
    bmp = paste0(files, '.bmp'),
    svglite = ,
    svg = paste0(files, '.svg'),
    current = files,
    stop('Unsupported device', call. = FALSE)
  )
  device <- switch(
    device,
    png = png,
    jpg = ,
    jpeg = jpeg,
    tif = ,
    tiff = tiff,
    bmp = bmp,
    svg = svg,
    svglite = svglite::svglite
  )

  pb <- progress_bar$new(
    'Rendering [:bar] at :fps fps ~ eta: :eta',
    total = length(frames)
  )
  start <- Sys.time()
  pb$tick(0)

  void <- future_mapply(frames, files, seq_along(frames), FUN = function(frame, file, i, stream, ..., plot, dims, pb = NULL) {
    if (!stream) {
      device(file, ...)
      on.exit(dev.off())
    }

    tryCatch(
      plot$scene$plot_frame(plot, frame, widths = dims$widths, heights = dims$heights),
      error = function(e) {
        warning(conditionMessage(e), call. = FALSE)
      }
    )

    if (!is.null(pb)) {
      rate <- i/as.double(Sys.time() - start, units = 'secs')
      if (is.nan(rate)) rate <- 0
      rate <- format(rate, digits = 2)
      pb$tick(tokens = list(fps = rate))
    }
  },
  MoreArgs = list(stream = stream, ..., plot = plot, dims = dims, pb = pb),
  SIMPLIFY = FALSE, USE.NAMES = FALSE)

  frame_vars <- plot$scene$frame_vars[frames, , drop = FALSE]
  if (!stream) frame_vars$frame_source <- files
  frame_vars
}

# Source: https://github.com/bbi-lab/bbi-sciatac-analyze/blob/master/src/reduce_dimensions.R

suppressMessages(library(monocle3))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(readr))
suppressMessages(library(Matrix))
suppressMessages(library(argparse))

########################################
# Message output.
########################################
message_log <- function(msg, ...) {
    message(msg, ...)
#    cat(paste(msg, ...), '\n')
}


########################################
# Quit message logging and success code.
########################################
quit_log <- function( msg, ...) {
    message_log(msg, ...)
    quit(save='no', status=0, runLast=TRUE)
}


########################################
# I/O.
########################################
.get_aux_files <- function(mtx_file) {
    base_name <- str_replace(mtx_file, '[.]gz$', '')
    base_name <- str_replace(base_name, '[.]mtx$', '')
    
    features_file <- paste0(base_name, '.rows.txt')
    cells_file <- paste0(base_name, '.columns.txt')
    
    return(list(features=features_file, cells=cells_file))
}

load_mtx_file <- function(mtx_file) {
    mat <- readMM(mtx_file)
    dim_files <- .get_aux_files(mtx_file)
    
    if (! file.exists(dim_files$features)) {
      quit_log(paste0(dim_files$features, ' file not found when loading ', mtx_file))
    }
    
    if (! file.exists(dim_files$cells)) {
      quit_log(paste0(dim_files$cells, ' file not found when loading ', mtx_file))
    }
    
    rownames(mat) <- read.delim(dim_files$features, header=FALSE)$V1
    colnames(mat) <- read.delim(dim_files$cells, header=FALSE)$V1
    return(as(mat, "dgCMatrix"))
}

########################################
# Utility functions.
########################################
filter_features_z <- function(bmat, lower_bound=-1.5, upper_bound=1.5, downsample=NULL) {
    feature_totals <- log10(Matrix::rowSums(bmat) + 1)
    avg <- mean(feature_totals)
    stdev <- sd(feature_totals)
  
    feature_z <- (feature_totals-avg)/stdev
  
    # Deal with corner case for missing or entirely non-zero windows (latter probably not a thing)
    feature_z[feature_totals == 0] <- -Inf
    feature_z[feature_totals == ncol(bmat)] <- Inf
  
    feature_totals <- feature_totals[feature_z > lower_bound & feature_z < upper_bound]
    
    if (!is.null(downsample)) {
      probabilities <- pnorm(feature_totals, mean=avg, sd=stdev)
      sampled_features <- sample(names(feature_totals), prob=probabilities, size=downsample, replace=FALSE)
      feature_totals <- feature_totals[sampled_features]
    }
    return(bmat[names(feature_totals),])
}


# Slightly better binarization from SnapATAC
binarize_matrix <- function(mat, outlier.filter=1e-3) {
    if (max(mat@x) == 1) {
        return(mat)
    }
    # identify the cutoff using outlier.filter
    count_cutoff <- max(1, quantile(mat@x, 1 - outlier.filter))
    mat@x[mat@x > count_cutoff] <- 0
    mat@x[mat@x > 0] <- 1
    return(mat)
}


filter_regions <- function(features, blacklist_df) {
    column_names <- c('chrom', 'start', 'end')
  
    if (! all(column_names %in% colnames(blacklist_df))) {
      quit_log('chrom, start, and end must be columns in blacklist df.')
    }
  
    features_df <- as.data.frame(str_split_fixed(features, '_', n=3))
    colnames(features_df) <- column_names
    features_df$start <- as.numeric(features_df$start)
    features_df$end <- as.numeric(features_df$end)
  
    features_df.gr <- GRanges(
      features_df$chrom,
      IRanges(features_df$start, features_df$end)
    )
  
    black_list.gr <- GRanges(
      blacklist_df$chrom,
      IRanges(blacklist_df$start, blacklist_df$end)
    )
  
    matching_hits <- queryHits(findOverlaps(features_df.gr, black_list.gr))

    return(features[-matching_hits])
}


test_peak_matrix <- function( peak_matrix, min_feature, min_cell, sample_name, umap_plot_file )
{
    num_feature <- dim(peak_matrix)[1]
    num_cell <- dim(peak_matrix)[2]
    bad_matrix_mesg <- character()
    if( num_feature < min_feature )
    {
        bad_matrix_mesg <- paste0(bad_matrix_mesg, '    too few features in peak matrix: ', num_feature, '\n')
    }
    if( num_cell < min_cell )
    {
        bad_matrix_mesg <- paste0(bad_matrix_mesg, '    too few cells in peak matrix: ', num_cell, '\n')
    }
    if(length(bad_matrix_mesg>0))
    {
        plot_message(paste0(umap_plot_file, '.pdf'), sample_name, bad_matrix_mesg)
        plot_message(paste0(umap_plot_file, '.png'), sample_name, bad_matrix_mesg)
        quit_log('ReduceDimensions: error:\n', bad_matrix_mesg, '  Stopping.\n')
    }
}


########################################
# Plotting functions.
########################################
plot_umap_pdf <- function(cds, umap_plot_file, sample_name) {
    umap_plot_file <- str_replace(umap_plot_file, '[.]pdf', '')
    umap_plot_file <- paste0(umap_plot_file, '.pdf')
    monocle3::plot_cells(cds, reduction_method='UMAP', show_trajectory_graph=FALSE) + ggplot2::ggtitle(sample_name)
    ggsave(umap_plot_file)
}

plot_umap_png <- function(cds, umap_plot_file, sample_name) {
    umap_plot_file <- str_replace(umap_plot_file, '[.]png', '')
    umap_plot_file <- paste0(umap_plot_file, '.png')
    monocle3::plot_cells(cds, reduction_method='UMAP', show_trajectory_graph=FALSE) + ggplot2::ggtitle(sample_name)
    ggsave(umap_plot_file)
}

plot_message <- function(umap_plot_file, sample_name, message) {
  ggplot2::ggplot() +
  ggplot2::geom_text(ggplot2::aes(x=1, y=1, label=message)) +
  monocle3:::monocle_theme_opts() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(x="UMAP 1", y="UMAP 2") +
  ggplot2::ggtitle(sample_name)
  ggplot2::ggsave(umap_plot_file)
}


########################################
# Preprocess peak matrix.
########################################
#
# Notes:
#   o  the default umi_cutoff, frip_cutoff, and frit_cutoff values are
#      meant to be absolute minima.
#
preprocess_peak_matrix <- function( mat_file, count_file, sample_name, umi_cutoff=100, frip_cutoff=0.1, frit_cutoff=0.05, black_list_file=NULL, doublet_predict=FALSE, doublet_predict_top_ntile=0.1, cds_file='monocle3_cds.rds', umap_plot_file='umap_null' )
{
    min_feature <- 10
    min_cell <- 10

    message_log('ReduceDimensions: sample: ', sample_name)

    # load peak matrix for sample.
    pMat <- load_mtx_file(mat_file)
  
    # binarize peak matrix
    pMat <- binarize_matrix(pMat)

    num_features_input <- dim(pMat)[1]
    num_cells_input <- dim(pMat)[2]
    message_log('ReduceDimensions: peak matrix dimensions: ', num_features_input, ' x ', num_cells_input)

    test_peak_matrix(pMat, min_feature, min_cell, sample_name, umap_plot_file)

    ######################################################################################
    # filter cells
    ######################################################################################
    # load summary stats for adding to col Data and reorder to match Cell order in matrices
    cDat <- read.table(count_file, header = TRUE)
    cDat_f <- cDat[match(colnames(pMat), cDat$cell),]

    row.names(cDat_f) <- cDat_f$cell
    colnames(cDat_f) <- c("cell", "total", "umi", "RIP", "RIT")
    cDat_f$FRIP <- cDat_f$RIP / cDat_f$umi
    cDat_f$FRIT <- cDat_f$RIT / cDat_f$umi
    cDat_f$umi_binarized <- Matrix::colSums(pMat)
  
    # filter cells based on unique reads, FRIP and FRIT  (cutoffs set above)

    message_log('ReduceDimensions: UMI cutoff: ', umi_cutoff)
    message_log('ReduceDimensions: FRIP cutoff: ', frip_cutoff)
    message_log('ReduceDimensions: FRIT cutoff: ', frit_cutoff)

    qc_cells <- filter(cDat_f, umi_binarized > umi_cutoff, FRIP > frip_cutoff, FRIT > frit_cutoff) %>% select(cell)
    pMat <- pMat[,colnames(pMat) %in% qc_cells$cell]

    num_features_cell_filter <- dim(pMat)[1]
    num_cells_cell_filter <- dim(pMat)[2]
    message_log('ReduceDimensions: peak matrix dimensions post-cell filter: ', num_features_cell_filter, ' x ', num_cells_cell_filter)

    test_peak_matrix(pMat, min_feature, min_cell, sample_name, umap_plot_file)

    ######################################################################################
    # filter features
    ######################################################################################
    # remove outlier features (z-score based) .
    pMat <- filter_features_z(bmat = pMat, lower_bound=-2, upper_bound=4, downsample=NULL)

    num_features_feature_filter <- dim(pMat)[1]
    num_cells_feature_filter <- dim(pMat)[2]
    message_log('ReduceDimensions: peak matrix dimensions post-feature filter: ', num_features_feature_filter, ' x ', num_cells_feature_filter)
  
    test_peak_matrix(pMat, min_feature, min_cell, sample_name, umap_plot_file)

    # the filter_features function from the atac_helper script removes features that overlap blackout regions of the genome.
    # load blacklist
    if(!is.null(black_list_file))
    {
        message_log('ReduceDimensions: black list region file: ', black_list_file)
        blacklist <- read.table(black_list_file, sep='\t')
        colnames(blacklist) <- c('chrom', 'start', 'end')
        features_f <- filter_regions(features=row.names(pMat), blacklist_df=blacklist)
        pMat <- pMat[row.names(pMat) %in% features_f,]
        num_features_feature_black_list_filter <- dim(pMat)[1]
        num_cells_feature_black_list_filter <- dim(pMat)[2]
        message_log('ReduceDimensions: peak matrix dimensions post-black-list filter: ', num_features_feature_black_list_filter, ' x ', num_cells_feature_black_list_filter)

    }

    ######################################################################################
    # filter cells again
    ######################################################################################
    # remove cells that have zero counts after feature filtering
    # Note: this is (probably) essential when doublet filtering
    #       with scrublet in order to prevent (numpy) divide by
    #       zero runtime warning because scrublet normalizes the
    #       cell read counts.
    umi_cells <- Matrix::colSums(pMat)
    pMat <- pMat[,umi_cells>0]

    num_features_cell_filter2 <- dim(pMat)[1]
    num_cells_cell_filter2 <- dim(pMat)[2]
    message_log('ReduceDimensions: peak matrix dimensions secondary post-cell filter: ', num_features_cell_filter2, ' x ', num_cells_cell_filter2)
  
    test_peak_matrix(pMat, min_feature, min_cell, sample_name, umap_plot_file)

    ######################################################################################
    # filter out doublets
    ######################################################################################
    if(doublet_predict)
    {
        message_log('ReduceDimensions: doublet top ntile cutoff: ', sprintf( '%.4f', doublet_predict_top_ntile))
        # this block calculates the doublet score and adds a column to the colData
        # it does not filter out doublets
        message_log('ReduceDimensions: read scrublet data')
        scrub_res <- read.table(paste0(sample_name, '-scrublet_table.csv'),sep=',')
        message_log('ReduceDimensions: read scrublet column names')
        cell_names <- read.table(paste0(sample_name, '-scrublet_columns.txt'),header=FALSE)
        scrub_res <- cbind(cell_names, scrub_res )
        colnames(scrub_res) <- c('cell', 'doublet_score', 'predicted_doublet')

        if(!anyNA(scrub_res$doublet_score) && !anyNA(scrub_res$predicted_doublet))
        {
            # mark top n-tile of cells with the highest doublet scores
            message_log('ReduceDimensions: mark top ntile doublets')
            threshold <- quantile(scrub_res$doublet_score, 1.0 - doublet_predict_top_ntile)
            message_log('ReduceDimensions: doublet top ntile score threshold: ', sprintf('%.4f', threshold))
            scrub_res$ntile_doublet <- sapply(scrub_res$doublet_score, function(x){
              ifelse(x < threshold, "singlet", "doublet")})
        }
        else
        {
            message_log('ReduceDimensions: disable doublet detection: scrublet returned at least one NA')
            scrub_res$ntile_doublet <- rep(NA, nrow(scrub_res))
        }

        num_rows_scrub_res <- nrow(scrub_res)
        num_cols_scrub_res <- ncol(scrub_res)
        message_log('ReduceDimensions: scrub_res dimensions: ', num_rows_scrub_res, ' x ', num_cols_scrub_res)
        num_rows_coldata_1 <- nrow(cDat_f)
        num_cols_coldata_1 <- ncol(cDat_f)
        message_log('ReduceDimensions: colData dimensions pre-inner_join: ', num_rows_coldata_1, ' x ', num_cols_coldata_1)
        if(anyDuplicated(scrub_res$cell))
        {
            message_log('ReduceDimensions: error: duplicate cell names in scrub_res')
            exit(-1)
        }
        if(anyDuplicated(colnames(pMat)))
        {
            message_log('ReduceDimensions: error: duplicate cell names in pMat')
            exit(-1)
        }
        cDat_f <- inner_join(cDat_f, scrub_res, by = "cell")

        num_rows_coldata_2 <- nrow(cDat_f)
        num_cols_coldata_2 <- ncol(cDat_f)
        message_log('ReduceDimensions: colData dimensions post-inner_join: ', num_rows_coldata_2, ' x ', num_cols_coldata_2)
        num_features_doublet_filter <- dim(pMat)[1]
        num_cells_doublet_filter <- dim(pMat)[2]
        num_ntile_doublets <- nrow(scrub_res[scrub_res$ntile_doublet=='doublet',])
        message_log('ReduceDimensions: top ntile doublet count: ', num_ntile_doublets)
        message_log('ReduceDimensions: peak matrix dimensions post-doublet filter: ', num_features_doublet_filter, ' x ', num_cells_doublet_filter)
    }

    cDat_f <- cDat_f[match(colnames(pMat), cDat_f$cell),]
    row.names(cDat_f) <- cDat_f$cell

    num_rows_coldata_3 <- nrow(cDat_f)
    num_cols_coldata_3 <- ncol(cDat_f)
    message_log('ReduceDimensions: colData dimensions end: ', num_rows_coldata_3, ' x ', num_cols_coldata_3)

    return(list(pMat=pMat, cDat_f=cDat_f))
}


######################################################################################
# make Monocle3 cell_data_set
######################################################################################
# Notes:
#  o  the align_cds is used to remove residual effects of differences in fragment
#     counts. Some people drop the first PC after PCA to do this but subtracting
#     out the effects more precisely targets the problem.

# Light version of monocle3 cds.
make_cds <- function(matrix_data, cds_file=NULL)
{
    message_log('ReduceDimensions: new_cell_data_set')
    cds <- monocle3::new_cell_data_set(matrix_data$pMat, cell_metadata=matrix_data$cDat_f)
    colData(cds)$n.umi <- Matrix::colSums(exprs(cds))
    message_log('ReduceDimensions: detect_genes')
    cds <- monocle3::detect_genes(cds, min_expr=0)
    if(!is.null(cds_file))
    {
        message_log('ReduceDimensions: write CDS file: ', cds_file)
        saveRDS(cds, cds_file)
    }
    return(cds)
}

make_monocle3_cds <- function(matrix_data, num_lsi_dimensions=75, cluster_resolution=1.0e-3, cds_file=NULL)
{
    message_log('ReduceDimensions: new_cell_data_set')
    cds <- monocle3::new_cell_data_set(matrix_data$pMat, cell_metadata=matrix_data$cDat_f)
    colData(cds)$n.umi <- Matrix::colSums(exprs(cds))
    message_log('ReduceDimensions: detect_genes')
    cds <- monocle3::detect_genes(cds, min_expr=0)
    message_log('ReduceDimensions: preprocess_cds')
    message_log('ReduceDimensions: number of LSI dimensions to keep: ', num_lsi_dimensions)
    cds <- monocle3::preprocess_cds(cds, method='LSI', num_dim=num_lsi_dimensions)
    message_log('ReduceDimensions: align_cds')
    cds <- monocle3::align_cds(cds, preprocess_method='LSI', residual_model_formula_str='~n.umi')
    message_log('ReduceDimensions: estimate_size_factors')
    cds <- monocle3::estimate_size_factors(cds)
    message_log('ReduceDimensions: reduce_dimension')
    cds <- monocle3::reduce_dimension(cds, preprocess_method='Aligned', reduction_method='UMAP')
    message_log('ReduceDimensions: cluster_cells')
    message_log('ReduceDimensions: cluster resolution: ', sprintf( '%.4e', cluster_resolution))
    cds <- monocle3::cluster_cells(cds, reduction_method='UMAP', resolution=cluster_resolution)
    if(!is.null(cds_file))
    {
        message_log('ReduceDimensions: write CDS file: ', cds_file)
        saveRDS(cds, cds_file)
    }
    return(cds)
}


######################################################################################
# write reduced dimensions to files
######################################################################################
write_reduced_dimensions <- function(cds, lsi_coords_file=NULL, umap_coords_file=NULL)
{
    if(!is.null(lsi_coords_file))
    {
        message_log('ReduceDimensions: write LSI coordinates file: ', lsi_coords_file)
        lsi_coords <- reducedDims(cds)$LSI
        readr::write_delim(data.frame(lsi_coords), file=lsi_coords_file, delim='\t')
    }

    if(!is.null(umap_coords_file))
    {
        message_log('ReduceDimensions: write UMAP coordinates file: ', umap_coords_file)
        umap_coords <- reducedDims(cds)$UMAP
        readr::write_delim(data.frame(umap_coords), file=umap_coords_file, delim='\t')
    }
    return(NULL)
}


######################################################################################
# make UMAP plot
######################################################################################
make_umap_plot <- function(cds, umap_plot_file=NULL, sample_name)
{
    if(!is.null(umap_plot_file))
    {
        message_log('ReduceDimensions: write UMAP plot: ', umap_plot_file)
        plot_umap_png(cds=cds, umap_plot_file, sample_name)
    }
    if(!is.null(umap_plot_file))
    {
        message_log('ReduceDimensions: write UMAP plot: ', umap_plot_file)
        plot_umap_pdf(cds=cds, umap_plot_file, sample_name)
    }
    return(NULL)
}