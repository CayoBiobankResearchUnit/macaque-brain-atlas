#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('rna','rna','seurat-umap-ref','15')

prefix = arguments[1]
analysis = arguments[2]
integration.method = arguments[3] # possible values: seurat-umap-ref, seurat-pca-ref, harmony, harmony-seurat, batchelor, none
this.cluster = as.integer(arguments[4])

# integration.method = 'harmony'

suppressMessages(library(monocle3))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(harmony))
suppressMessages(library(Seurat))
suppressMessages(library(parallel))

# This script is for performing subclustering from scratch (redoing preprocessing)

this.sample = 'all'

meta.file = file.path('umap',paste0(prefix,'-scanpy-final-classified.rds'))
this.umap = readRDS(meta.file)

expr.file = file.path('mm/recluster',paste0(prefix,'_class',this.cluster,'.mtx'))
rows.file = file.path('mm/recluster',paste0(prefix,'_class',this.cluster,'_rows.txt.gz'))
cols.file = file.path('mm/recluster',paste0(prefix,'_class',this.cluster,'_cols.txt.gz'))
features.file = file.path('stats/var_final',paste0('var_',prefix,'_all.txt.gz'))

e = t(readMM(expr.file))
these.genes = rownames(e) = fread(rows.file,header=FALSE)$V1
these.cells = colnames(e) = fread(cols.file,header=FALSE)$V1

meta = this.umap[these.cells,]
names(meta)[names(meta) %in% c('umap.1','umap.2','partition','cluster')] = paste0(names(meta)[names(meta) %in% c('umap.1','umap.2','partition','cluster')],'_round1')

features = as.data.frame(fread(features.file))
rownames(features) = features$id
features$V1 = NULL
features = features[intersect(rownames(features),these.genes),]

e = e[rownames(features),]

# Assemble into new monocle3 cell data set
cds = new_cell_data_set(
	expression_data=e,
	cell_metadata=meta,
	gene_metadata=features)

cds = detect_genes(cds,min_expr=0)
cds = cds[rowData(cds)$num_cells_expressed > 0,]

cds = estimate_size_factors(cds)
cds = preprocess_cds(cds, method='PCA')
push.status(paste0('preprocess_cds cluster_',this.cluster))

if (integration.method %in% c('harmony','harmony-seurat')) {
	# Regress out UMI
	cds = align_cds(cds,preprocess_method='PCA',residual_model_formula_str='~total_counts')
	this.pca = reducedDims(cds)$Aligned
	
	if (integration.method == 'harmony') {
		# Use harmony directly on PCA matrix
		reducedDims(cds)$Aligned = harmony::HarmonyMatrix(
			data_mat = this.pca, ## PCA embedding matrix of cells
			meta_data = as.data.frame(colData(cds)), ## dataframe with cell labels
			theta = NULL, ## cluster diversity enforcement
			lambda = NULL, ## ridge regression penalty parameter
			sigma = 0.1, # width of soft kmeans clusters
			nclust = NULL, # number of clusters in model. nclust=1 equivalent to simple linear regression
			tau = 0, # Protection against overclustering small datasets with large ones. tau is the expected number of cells per cluster.
			block.size = 0.05, # what proportion of cells to update during clustering
			epsilon.cluster = 1e-05, # convergence tolerance for clustering round of Harmony
			epsilon.harmony = 1e-04, # convergence tolerance for Harmony
			vars_use = 'sequencing_run_id', ## variable to integrate out
			max.iter.harmony = 10, ## stop after maximum 10 iterations
			max.iter.cluster = 20, ## stop after maximum 20 iterations
			return_object = FALSE, ## return corrected matrix only
			do_pca = FALSE ## don't recompute PCs
		)
	} else {
		# Use harmony via Seurat object (should be identical)
		
		# Create Seurat object
		cds.seurat = CreateSeuratObject(counts=assays(cds)$counts,meta.data=as.data.frame(colData(cds)),project=prefix)

		# Extract monocle PCA object and slot it into the Seurat object
		pca.monocle = this.pca
		colnames(pca.monocle) = paste0('PC_',1:ncol(pca.monocle))

		cds.seurat$pca = SeuratObject::CreateDimReducObject(
			embeddings=pca.monocle,
			assay='RNA',key='PC_'
		)

		# Run harmony on Seurat-formatted data
		cds.harmony = RunHarmony(cds.seurat,group.by.vars='sequencing_run_id',project.dim=FALSE)
		push.status(paste0('RunHarmony cluster_',this.cluster))

		# Extract harmony cell embeddings and slot it into monocle Aligned slot (where batchelor output normally goes)
		reducedDims(cds)$Aligned = cds.harmony$harmony@cell.embeddings
		colnames(reducedDims(cds)$Aligned) = gsub('^harmony_','PC',colnames(reducedDims(cds)$Aligned))
		
		rm(list=c('cds.seurat','cds.harmony','pca.monocle'))
	}
	prefix = paste0(prefix,'_harmony')
} else if (integration.method %in% c('seurat-umap-ref','seurat-pca-ref')) {
	# Create Seurat object
	cds.seurat = CreateSeuratObject(counts=assays(cds)$counts,meta.data=as.data.frame(colData(cds)),project=prefix)
	
	# Split by batch ID to make a list
	batch.list = SplitObject(cds.seurat, split.by = 'sequencing_run_id')
	
	# # If any batches only have one cell, merge all queries together
	# if (any(unlist(lapply(batch.list,function(x) ncol(x))) == 1)) {
	# 	batch.list = c(
	# 		list(query = Reduce(merge,batch.list[setdiff(names(batch.list),reference.batch)])),
	# 		batch.list[reference.batch]
	# 	)
	# }
	
	# If any batches have fewer than 10 cells, just do a UMAP on everything together
	if (any(unlist(lapply(batch.list,function(x) ncol(x))) < 10)) {
		message('Some batches are too small. Skipping integration')
		
		reference.cds = NormalizeData(cds.seurat,verbose=FALSE)
		reference.cds = FindVariableFeatures(reference.cds,selection.method='vst',nfeatures=2000,verbose = FALSE)
		reference.cds = ScaleData(reference.cds,vars.to.regress='total_counts',verbose = FALSE)
		reference.cds = RunPCA(reference.cds, npcs = min(ncol(reference.cds)-1,50), verbose = FALSE)
				
		reducedDims(cds)$Aligned = reference.cds$pca@cell.embeddings[rownames(reducedDims(cds)$PCA),]
		
		if (integration.method == 'seurat-umap-ref') {
			reference.cds = RunUMAP(reference.cds, dims = 1:min(ncol(reference.cds)-1,50), reduction='pca', n.neighbors = 15, min.dist = 0.25, spread = 1, return.model = TRUE)
			reducedDims(cds)$UMAP = reference.cds$umap@cell.embeddings[rownames(reducedDims(cds)$PCA),]
		}
		rm(list=c('cds.seurat','batch.list'))
	} else {
		# Assemble reference dataset. Integrate reference batches if there are more than one
		if (length(reference.batch) > 1) {
			# reference.anchors = FindIntegrationAnchors(object.list = batch.list[reference.batch], dims = 1:50)
			# reference.cds = IntegrateData(anchorset = reference.anchors, dims = 1:50)
			reference.cds = Reduce(merge,batch.list[reference.batch])
		} else {
			reference.cds = batch.list[[reference.batch]]
		}
		reference.cds = NormalizeData(reference.cds,verbose=FALSE)
		reference.cds = FindVariableFeatures(reference.cds,selection.method='vst',nfeatures=2000,verbose = FALSE)
		
		# Subset query datasets
		query.list = batch.list[setdiff(names(batch.list),reference.batch)]
		rm(list=c('batch.list'))
		
		# Normalize and find variable features within each batch
		query.list = mclapply(query.list,function(x) {
			x = NormalizeData(x,verbose=FALSE)
			x = FindVariableFeatures(x,selection.method='vst',nfeatures=2000,verbose = FALSE)
			x = ScaleData(x,vars.to.regress='total_counts',verbose = FALSE)
			x
		},mc.cores=min(length(query.list),n.cores))
	
		# Finish preprocessing reference. Run PCA.
		reference.cds = ScaleData(reference.cds,vars.to.regress='total_counts',verbose = FALSE)
		reference.cds = RunPCA(reference.cds, npcs = min(ncol(reference.cds)-1,50), verbose = FALSE)
		
		if (integration.method == 'seurat-umap-ref') {
			reference.cds = RunUMAP(reference.cds, dims = 1:min(ncol(reference.cds)-1,50), reduction='pca', n.neighbors = 15, min.dist = 0.25, spread = 1, return.model = TRUE)
			
			# For each query dataset, perform reference-guided integration
			# integrated.list = vector(length=length(query.list),mode='list')
			# names(integrated.list) = names(query.list)
			integrated.list = mclapply(1:length(query.list),function(i) {
				this.anchors = FindTransferAnchors(
					reference = reference.cds,
					query = query.list[[i]],
					dims = 1:min(ncol(reference.cds)-1,50),
					reference.reduction = 'pca',
					k.anchor = min(ncol(reference.cds)-1,ncol(query.list[[i]])-1,5), # Make sure k.anchor < number of query cells
					k.score = min(ncol(reference.cds)-1,ncol(query.list[[i]])-1,30), # Make sure k.score < number of query cells
					k.filter = min(ncol(reference.cds),ncol(query.list[[i]]),200) # Make sure k.filter <= number of query cells
				)
				out = MapQuery(
					anchorset = this.anchors,
					reference = reference.cds,
					query = query.list[[i]],
					refdata = NULL,
					reference.reduction = 'pca',
					reduction.model = 'umap',
					new.reduction.name = 'pca',
					integrateembeddings.args = list(k.weight=min(ncol(reference.cds),ncol(query.list[[i]]),100)) # Make sure k.weight <= number of query cells
				)
				out
			},mc.cores=min(length(query.list),n.cores))
			names(integrated.list) = names(query.list)
			push.status(paste0('Seurat integrate cluster_',this.cluster))
	
			# Concatenate all PCAs back together and reorder them
			reducedDims(cds)$Aligned = do.call(rbind,lapply(c(list(reference=reference.cds),integrated.list),function(x) {
				pc.embeddings = x$pca@cell.embeddings
				colnames(pc.embeddings) = gsub('^PC_','PC',gsub('^pca_','PC',colnames(pc.embeddings)))
				pc.embeddings
			}))[rownames(reducedDims(cds)$PCA),]
	
			reducedDims(cds)$UMAP = do.call(rbind,lapply(c(list(reference=reference.cds),integrated.list),function(x) {
				umap.embeddings = x[[intersect(c('ref.umap','umap'),names(x))]]@cell.embeddings
				colnames(umap.embeddings) = NULL
				umap.embeddings
			}))[rownames(reducedDims(cds)$PCA),]
		} else {

			# For each query dataset, perform reference-guided integration
			# integrated.list = vector(length=length(query.list),mode='list')
			# names(integrated.list) = names(query.list)
			integrated.list = mclapply(1:length(query.list),function(i) {
				this.anchors = FindTransferAnchors(
					reference = reference.cds,
					query = query.list[[i]],
					dims = 1:min(ncol(reference.cds)-1,50),
					reference.reduction = 'pca',
					k.anchor = min(ncol(reference.cds)-1,ncol(query.list[[i]])-1,5), # Make sure k.anchor < number of query cells
					k.score = min(ncol(reference.cds)-1,ncol(query.list[[i]])-1,30), # Make sure k.score < number of query cells
					k.filter = min(ncol(reference.cds),ncol(query.list[[i]]),200) # Make sure k.filter <= number of query cells
				)
				out = IntegrateEmbeddings(
					anchorset = this.anchors,
					reference = reference.cds,
					query = query.list[[i]],
					new.reduction.name = 'pca',
					reductions = 'pcaproject',
					dims.to.integrate = 1:min(ncol(reference.cds)-1,50),
					k.weight=min(ncol(reference.cds),ncol(query.list[[i]]),100) # Make sure k.weight <= number of query cells
				)
				out
			},mc.cores=min(length(query.list),n.cores))
			names(integrated.list) = names(query.list)
			push.status(paste0('Seurat integrate cluster_',this.cluster))
	
			# Concatenate all PCAs back together and reorder them
			reducedDims(cds)$Aligned = do.call(rbind,lapply(c(list(reference=reference.cds),integrated.list),function(x) {
				pc.embeddings = x$pca@cell.embeddings
				colnames(pc.embeddings) = gsub('^PC_','PC',gsub('^pca_','PC',colnames(pc.embeddings)))
				pc.embeddings
			}))[rownames(reducedDims(cds)$PCA),]
		}
		rm(list=c('integrated.list','query.list','reference.cds'))
	}
	if (integration.method == 'seurat-umap-ref') {
		prefix = paste0(prefix,'_seuratref')
	} else if (integration.method =='seurat-pca-ref') {
		prefix = paste0(prefix,'_seuratrefpca')
	}
}  else if (integration.method == 'batchelor') {
	# Use batchelor batch correction
	cds = align_cds(cds, alignment_group='sequencing_run_id',residual_model_formula_str='~total_counts',preprocess_method='PCA')
	push.status(paste0('align_cds cluster_',this.cluster))
	prefix = paste0(prefix,'_batchelor')
} else if (integration.method == 'none') {
	# If no integration, then do the minimum (regress out UMI)
	cds = align_cds(cds,residual_model_formula_str='~total_counts',preprocess_method='PCA')
	# reducedDims(cds)$Aligned = reducedDims(cds)$PCA
	prefix = paste0(prefix,'_nobatch')
}

# Skip UMAP if UMAP was already calculated from integration
if (!integration.method %in% 'seurat-umap-ref') {
	cds = reduce_dimension(cds, reduction_method='UMAP',preprocess_method='Aligned',umap.min_dist=0.25,umap.n_neighbors=15)
	push.status(paste0('reduce_dimension cluster_',this.cluster))
}

cds = cluster_cells(cds,reduction_method='Aligned',cluster_method='leiden',resolution=1e-4)
push.status(paste0('cluster_cells cluster_',this.cluster))

out.umap = reducedDim(cds,'UMAP')
colnames(out.umap) = c('umap.1','umap.2')
this.umap = data.frame(as.data.frame(colData(cds)),out.umap,partition2=factor(cds@clusters$Aligned$partitions),cluster2=factor(cds@clusters$Aligned$clusters))

# temp
# prefix = paste0(prefix,'_nobatch')

dir.create('rds/recluster',showWarnings=FALSE)
saveRDS(cds,file=file.path('rds/recluster',paste0(prefix,'_cds_class',this.cluster,'.rds')))
saveRDS(this.umap,file=file.path('rds/recluster',paste0(prefix,'_df_class',this.cluster,'.rds')))

# Control overplotting based on how many cells
if ((1e6/nrow(this.umap)) < 1) {
	size=0.1
	alpha=0.05
} else if ((1e6/nrow(this.umap)) < 10) {
	size=0.2
	alpha=0.1
} else if ((1e6/nrow(this.umap)) < 100) {
	size=0.5
	alpha=0.5
} else {
	size=1
	alpha=0.5
}

set.seed(42)
this.umap = this.umap[sample(1:nrow(this.umap)),]

dir.create('figures/recluster',showWarnings=FALSE)

plot.umap(this.umap,color='cluster2',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-cluster-class',this.cluster,'_presentation.pdf')),color.label='Cluster',legend=FALSE,presentation=TRUE,size=size,alpha=alpha)
plot.umap(this.umap,color='partition2',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-partition-class',this.cluster,'_presentation.pdf')),color.label='Partition',legend=FALSE,presentation=TRUE,size=size,alpha=alpha)
plot.umap(this.umap,color='region_major',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-regionmajor-class',this.cluster,'_presentation.pdf')),color.label='Class',legend=FALSE,presentation=TRUE,size=size,alpha=alpha)
plot.umap(this.umap,color='region',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-region-class',this.cluster,'_presentation.pdf')),color.label='Region',legend=FALSE,presentation=TRUE,size=size,alpha=alpha)
plot.umap(this.umap,color='sequencing_run_id',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-batch-class',this.cluster,'_presentation.pdf')),color.label='Batch',legend=FALSE,presentation=TRUE,size=size,alpha=alpha)

plot.umap(this.umap,color='cluster2',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-cluster-class',this.cluster,'.pdf')),color.label='Cluster',legend=TRUE,presentation=FALSE,size=size,alpha=alpha)
plot.umap(this.umap,color='partition2',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-partition-class',this.cluster,'.pdf')),color.label='Partition',legend=TRUE,presentation=FALSE,size=size,alpha=alpha)
plot.umap(this.umap,color='region_major',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-regionmajor-class',this.cluster,'.pdf')),color.label='Class',legend=TRUE,presentation=FALSE,size=size,alpha=alpha)
plot.umap(this.umap,color='region',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-region-class',this.cluster,'.pdf')),color.label='Region',legend=TRUE,presentation=FALSE,size=size,alpha=alpha)
plot.umap(this.umap,color='sequencing_run_id',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-batch-class',this.cluster,'.pdf')),color.label='Batch',legend=TRUE,presentation=FALSE,size=size,alpha=alpha)


# plot.umap(this.umap,color='sequencing_run_id',color.label='Batch',legend=TRUE,presentation=FALSE,size=size,alpha=alpha)


# plot.umap(this.umap,color='sequencing_run_id',file='test-cell-type3.pdf',color.label='Batch',legend=TRUE,presentation=FALSE,size=size,alpha=alpha)

# plot.umap(this.umap,color='sequencing_run_id',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-batch-class',this.cluster,'.pdf')),color.label='Cluster',legend=TRUE,presentation=FALSE,size=size,alpha=alpha)
# plot.umap(this.umap,color='animal_id',file=file.path('figures/umap-recluster',paste0('umap-',prefix,'-animal-class',this.cluster,'.pdf')),color.label='Partition',legend=TRUE,presentation=FALSE,size=size,alpha=alpha)


