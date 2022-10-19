#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('rna','rna')

prefix = arguments[1]
analysis = arguments[2]

this.umap = readRDS(file.path('umap',paste0(prefix,'-scanpy-final-classified.rds')))

# Read in UMAPs

umaps = do.call(rbind,mclapply(1:nlevels(this.umap$cell_round1_level1),function(this.cluster) {
	that.umap = readRDS(file.path('umap',paste0(prefix,'-scanpy-recluster-class',this.cluster,'-summarized.rds')))
	that.umap = subset(that.umap,select=c(
		'umap_sub.1','umap_sub.2','subcluster','subcluster2','subpartition','cell_subcluster','cell_subcluster2','cell_subpartition'
	))
	that.umap
},mc.cores=n.cores))

umaps = umaps[row.names(this.umap),]

if (!identical(rownames(umaps),rownames(this.umap))) stop('Error')

this.umap = data.frame(this.umap,umaps)

this.umap$cell_class = this.umap$cell_round1_level1
this.umap$cell_class_integer = as.integer(this.umap$cell_round1_level1)

# Manual resolution

this.umap$cell_type = NA

# Preliminary
this.umap$subcluster_manual = NA

hi.res = c(1,3,6,8,9,11,12) # For these, I decided to subcluster with higher resolution
lo.res = c(2,10) # For these, I decided to subcluster with lower resolution
eq.res = c(4,5,7,14,15,17,18,19) # For these, my decision didn't matter since they're the same
no.res = c(13,16) # For these, clusters seemed spurious and I decided not to subcluster at all

this.umap[this.umap$cell_class_integer %in% hi.res,]$subcluster_manual = this.umap[this.umap$cell_class_integer %in% hi.res,]$subcluster
this.umap[this.umap$cell_class_integer %in% lo.res,]$subcluster_manual = this.umap[this.umap$cell_class_integer %in% lo.res,]$subcluster2
this.umap[this.umap$cell_class_integer %in% eq.res,]$subcluster_manual = this.umap[this.umap$cell_class_integer %in% eq.res,]$subcluster
this.umap[this.umap$cell_class_integer %in% no.res,]$subcluster_manual = this.umap[this.umap$cell_class_integer %in% no.res,]$subpartition

this.umap$cell_subcluster_manual = with(this.umap,paste(cell_class,gsub('^[0-9]+?-','',subcluster_manual)))
this.umap$cell_subcluster_manual[this.umap$cell_class_integer %in% no.res] = as.character(this.umap$cell_class)[this.umap$cell_class_integer %in% no.res]

cell.subcluster.labels = do.call(rbind,lapply(split(this.umap,this.umap$cell_class),function(x) {
	this.class = unique(x$cell_class_integer)
	this.class.label = as.character(unique(x$cell_class))
	subclusters = unique(x$subcluster_manual)
	if (length(subclusters) > 1) {
		subclusters = subclusters[order(unlist(lapply(strsplit(subclusters,'-'),function(x) as.integer(x[2]))))]
		data.frame(subcluster_manual=subclusters,subcluster_label=paste(this.class.label,gsub(paste0('^',this.class,'-'),'',subclusters)))
	} else {
		data.frame(subcluster_manual=subclusters,subcluster_label=this.class.label)
	}
}))

this.umap$cell_subcluster_manual = factor(this.umap$subcluster_manual,levels=cell.subcluster.labels$subcluster_manual,labels=cell.subcluster.labels$subcluster_label)
this.umap$subcluster_manual = factor(this.umap$subcluster_manual,levels=cell.subcluster.labels$subcluster_manual)

# Now output human-readable labels

# Get rid of cell_final column because it's confusing
this.umap$cell_round1_type = this.umap$cell_round1_level3
this.umap$cell_round1_level3 = NULL # Backed up as cell_round1_type
this.umap$cell_round1_level2 = NULL # Unused
this.umap$cell_round1_level1 = NULL # Backed up as cell_round1_id

this.umap$cell_final = NULL
this.umap$cell_round1_id = NULL

this.umap$subcluster_hi_res = this.umap$subcluster
this.umap$subcluster_lo_res = this.umap$subcluster2
this.umap$cell_subcluster_hi_res = this.umap$cell_subcluster
this.umap$cell_subcluster_lo_res = this.umap$cell_subcluster2

this.umap$subcluster = NULL
this.umap$subcluster2 = NULL
this.umap$cell_subcluster = NULL
this.umap$cell_subcluster2 = NULL

this.umap$region_subclass = factor(this.umap$region,levels=region.hierarchy$region,labels=region.hierarchy$region.subclass)
this.umap$region_class = factor(this.umap$region,levels=region.hierarchy$region,labels=region.hierarchy$region.class)
this.umap$region_label = factor(this.umap$region,levels=region.hierarchy$region,labels=region.hierarchy$region.full)

dobs = c('2C0' = '2006-12-22', '3R7' = '2013-08-19', '3I4' = '2009-11-08', '4I3' = '2009-09-17', '6J2' = '2009-09-10')
dods = c('2C0' = '2016-12-20', '3R7' = '2016-12-14', '3I4' = '2019-10-25', '4I3' = '2019-10-24', '6J2' = '2019-10-31')

dobs = do.call(c,lapply(dobs,as.Date))
dods = do.call(c,lapply(dods,as.Date))

this.umap$dob = dobs[this.umap$animal_id]
this.umap$dod = dods[this.umap$animal_id]

this.umap$age = with(this.umap,as.numeric(dod-dob)/365.2425)

# Function for smartly reordering levels of subclusters/subpartitions
reorder.levels = function(x,sep='-',recode=NULL,letters=FALSE) {
	out = apply(do.call(rbind,strsplit(unique(as.character(x)),sep)),2,as.integer)
	if (!is.null(recode)) {
		# If recode is set, use the labels in recode to replace the cell levels
		out = as.data.frame(out[order(out[,1],out[,2]),])
		out$V1 = recode[out$V1]
		apply(out,1,function(x) {
			if (table(out$V1)[x[1]] > 1) {
				paste(x[1],if (letters) LETTERS[as.integer(x[2])] else x[2],sep=' ')
			} else {
				x[1]
			}
		})
	} else {
		apply(out[order(out[,1],out[,2]),],1,paste,collapse=sep)
	}
}

this.umap = within(this.umap,{
	id = factor(id)
	scrublet_call = factor(scrublet_call,levels=c('Singlet','Doublet'))
	biccn_id = factor(biccn_id)
	hemisphere = factor(hemisphere,levels=c('R','L'))
	animal_id = factor(animal_id)
	social_group = factor(social_group)
	sex = factor(sex)
	extractor = factor(extractor)
	isolation_site = factor(isolation_site)
	extraction_batch = factor(extraction_batch)
	fixation = factor(fixation)
	lysis = factor(lysis)
	sequencing_run_id = factor(sequencing_run_id)
	subpartition = factor(subpartition,levels=reorder.levels(subpartition))
	cell_subpartition =  factor(subpartition,levels=reorder.levels(subpartition,recode=levels(cell_class),letters=TRUE))
	subcluster_hi_res = factor(subcluster_hi_res,levels=reorder.levels(subcluster_hi_res))
	cell_subcluster_hi_res = factor(subcluster_hi_res,levels=reorder.levels(subcluster_hi_res),labels=reorder.levels(subcluster_hi_res,recode=levels(cell_class)))
	subcluster_lo_res = factor(subcluster_lo_res,levels=reorder.levels(subcluster_lo_res))
	cell_subcluster_lo_res = factor(subcluster_lo_res,levels=reorder.levels(subcluster_lo_res),labels=reorder.levels(subcluster_lo_res,recode=levels(cell_class)))
})

this.umap = this.umap[,c(
'id',                      # library ID
'cell',                    # cell ID
'Size_Factor',             # size factor
'n.umi',                   # Total library size
'perc_mitochondrial_umis', # percentage of UMI mapping to mt genome
'scrublet_score',          # Scrublet kNN score (BBI preprocessing)
'scrublet_call',           # Scrublet call (BBI preprocessing)
'biccn_id',                # Tissue sample ID
'region',                  # brain region (abbreivation)
'region_label',            # brain region (full name)
'region_class',            # Broadest-level brain anatomical unit
'region_subclass',         # Intermediate-level brain anatomical unit
# 'animal_integer',        # integer assigned to animal (dropped)
'hemisphere',              # brain hemisphere (NAs indicate midline structures)
'animal_id',               # Cayo Santiago donor animal ID
'social_group',            # Cayo Santiago social group ID
'age',                     # Age (years)
'sex',                     # Sex
'extractor',               # Nuclei extracted by
'isolation_site',          # Nuclei isolation location
# 'organization_order',    # Order in which samples were organized (drop)
# 'dissection_order',      # Order in which samples were dissected (drop)
# 'run_id',                # sci-*-seq ID
'isolation_order',         # Order isolated
'extraction_batch',        # Extraction batch ID
'extraction.date',         # Extraction date
'fixation',                # Fixation method
# 'fixation_volume',       # Volume fixation
'lysis',                   # Lysis method
# 'freezing_medium',       # Medium in which isolated nuclei were flash-frozen (drop)
'sequencing_run_id',       # Sequencing batch
'doublet_score',           # Doublet score (Reran Scrublet)
# 'predicted_doublet',     # whether predicted doublet (all FALSE)
# 'manual_doublet',        # whether manual doublet
# 'is_singlet',            # Whether singlet
'n_genes_by_counts',       # Number of genes expressed
'total_counts',            # Library size used for normalization in scanpy (slightly different than BBI preprocessed due to filtering)
'total_counts_mt',         # Total counts mitochondrial calculated in scanpy (slightly different than BBI preprocessed due to filtering)
'pct_counts_mt',           # percentage counts mitochondrial (slightly different than BBI preprocessed due to filtering)
'umap.1',                  # Global UMAP dimension 1
'umap.2',                  # Global UMAP dimension 2
'partition',               # Global UMAP partition
'cluster',                 # Global UMAP cluster (resolution = 1e-5)
'cell_round1_type',        # Cell types called from global analysis
'cell_class_integer',      # Cell class integers
'umap_sub.1',              # Cell-class-specific UMAP dimension 1 (only meaningful for each cell_class level)
'umap_sub.2',              # Cell-class-specific UMAP dimension 2 (only meaningful for each cell_class level)
'subpartition',            # Cell-class subset reclustering partition ID
'cell_subpartition',       # Cell-class subset reclustering partition label (synonymous with subpartition)
'cell_class',              # Cell class (manually curated levels used for reclustering)
'cell_type',               # Cell type (cell types called from global UMAP), placeholder for manual biologically meaningful layers (e.g, breaking down vascular into endothelial and pericytes, see cell_round1_type for similar attempt from the combined global analysis)
'subcluster_hi_res',       # Cell subcluster called with high-resolution settings (resolution=1e-5)
'cell_subcluster_hi_res',  # Synonymous with above, except with more readable labels
'subcluster_lo_res',       # Cell subcluster called with low-resolution settings (resolution=1e-6)
'cell_subcluster_lo_res',  # Synonymous with above, except with more readable labels
'subcluster_manual',       # Manually assigned cell-class-specific clusters (mix of resolution settings by cell class)
'cell_subcluster_manual'   # Synonymous with above, except with more readable labels
)]

# Make some temporary corrections
this.umap = within(this.umap,{
	cell_round1_type = factor(cell_round1_type,levels=levels(cell_round1_type),labels=gsub('neuron$','neurons',gsub('serotinergic','serotonergic',levels(cell_round1_type))))
	cell_subpartition = factor(cell_subpartition,levels=levels(cell_subpartition),labels=gsub('neuron$','neurons',gsub('serotinergic','serotonergic',levels(cell_subpartition))))
	cell_class = factor(cell_class,levels=levels(cell_class),labels=gsub('neuron$','neurons',gsub('serotinergic','serotonergic',levels(cell_class))))
	cell_subcluster_hi_res = factor(cell_subcluster_hi_res,levels=levels(cell_subcluster_hi_res),labels=gsub(' +',' ',gsub('neuron ','neurons ',gsub('neuron$','neurons',gsub('serotinergic','serotonergic',levels(cell_subcluster_hi_res))))))
	cell_subcluster_lo_res = factor(cell_subcluster_lo_res,levels=levels(cell_subcluster_lo_res),labels=gsub(' +',' ',gsub('neuron ','neurons ',gsub('neuron$','neurons',gsub('serotinergic','serotonergic',levels(cell_subcluster_lo_res))))))
	cell_subcluster_manual = factor(cell_subcluster_manual,levels=levels(cell_subcluster_manual),labels=gsub(' +',' ',gsub('neuron ','neurons ',gsub('neuron$','neurons',gsub('serotinergic','serotonergic',levels(cell_subcluster_manual))))))
})

# Correct regions
region.corrections = data.frame(region=c('dmPFC','vmPFC','dlPFC','vlPFC','ACC','CC','CN','NAc','EC','PC','A1','AMY','HIP','M1','mdTN','vlTN','LGN','S1','IPP','SPP','STS','MT','IT','V1','CV','lCb','MB','MdO','MdC','Pons'),region_label=c('dorsomedial prefrontal cortex','ventromedial prefrontal cortex','dorsolateral prefrontal cortex','ventrolateral prefrontal cortex','anterior cingulate cortex','corpus callosum','head of the caudate nucleus','nucleus accumbens','entorhinal cortex','perirhinal cortex','primary auditory cortex','amygdala','hippocampus','primary motor cortex','mediodorsal thalamic nucleus','ventrolateral thalamic nucleus','lateral geniculate nucleus','primary somatosensory cortex','inferior posterior parietal gyrus','superior posterior parietal gyrus','superior temporal sulcus','middle temporal visual area','inferior temporal gyrus','primary visual cortex','cerebellar vermis','lateral cerebellar cortex','midbrain','open medulla','closed medulla','pons'))

this.umap = within(this.umap,{
	region = factor(region,levels=region.corrections$region)
	region_label = factor(region_label,levels=region.corrections$region_label)
})

# Anything from batch 3
this.umap[(this.umap$sequencing_run_id %in% 'Snyder-Mackler_RNA3-036') & (this.umap$region == 'STS'),]$region = 'MT'
this.umap[(this.umap$sequencing_run_id %in% 'Snyder-Mackler_RNA3-036') & (this.umap$region == 'MT'),]$region_label = 'middle temporal visual area'

# Write final clusters to file
this.umap.split = split(this.umap,this.umap$cell_class)
for (i in levels(this.umap$cell_class)) {
	x = this.umap.split[[i]]
	this.class = unique(x$cell_class_integer)
	write.table(x[,c('cell','subcluster_manual')],file=file.path('stats/subclusters',paste0(prefix,'-class',this.class,'-clusters.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
	write.table(x[,c('cell','cell_subcluster_manual')],file=file.path('stats/subclusters',paste0(prefix,'-class',this.class,'-cellclusters.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
}

write.table(this.umap[,c('cell','subcluster_manual')],file=file.path('stats/subclusters',paste0(prefix,'-subclusters.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(this.umap[,c('cell','cell_subcluster_manual')],file=file.path('stats/subclusters',paste0(prefix,'-cellsubclusters.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
write(levels(this.umap$cell_subcluster_manual),file=file.path('stats/subclusters',paste0(prefix,'-cellsubclusters-levels.txt')),sep='\n')

saveRDS(this.umap,file=file.path('umap',paste0(prefix,'-scanpy-recluster-classified.rds')))

# Write finalized sample to file

this.meta = unique(subset(this.umap,select=c('id','biccn_id','region','region_label','region_class','region_subclass','hemisphere','animal_id','social_group','age','sex')))
rownames(this.meta) = this.meta$id

# Make a new column that is programming friendly (region_subclass values have spaces)
this.meta$region_subclass2 = with(this.meta,factor(region_subclass,levels=levels(region_subclass),labels=gsub(' ','_',levels(region_subclass))))

write.table(this.meta,file=file.path('data',paste0(prefix,'_metadata_final.txt')),sep='\t',row.names=TRUE,col.names=NA,quote=FALSE)
