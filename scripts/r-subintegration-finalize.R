#!/usr/bin/env Rscript
# 🐛

source('scripts/_include_options.R')
source('scripts/_include_functions.R')

arguments = commandArgs(trailing=TRUE)
# arguments = c('atac','atac')

prefix = arguments[1]
analysis = arguments[2]
glue.prefix = if (length(arguments) > 2) arguments[3] else 'biccn'

this.umap = readRDS(file.path('umap',paste0(prefix,'-scanpy-final-classified.rds')))

# Read in RNA levels

cell.levels = scan(file.path('stats/clusters','rna-final-cellclasses-levels.txt'),what='',sep='\n')

glue.predictions = read.delim(file.path('stats/clusters',paste0(prefix,'-celltype-predictions.txt')),row.names=1)

glue.predictions$glue_type = factor(glue.predictions$glue_type,levels=cell.levels)
glue.predictions$cell = factor(glue.predictions$cell,levels=rownames(this.umap))
glue.predictions = glue.predictions[order(glue.predictions$cell),]

if (!identical(rownames(glue.predictions),rownames(this.umap))) {
	stop('Error')
}

this.umap = data.frame(this.umap,glue.predictions[,c('glue_type','glue_type_confidence')])

glue.predictions = do.call(rbind,parallel::mclapply(1:length(cell.levels),function(this.cluster) {
	# If cells exist and the files exist
	if (as.logical(table(this.umap$glue_type)[cell.levels[this.cluster]]) & file.exists(file.path('stats/clusters/subintegration',paste0(glue.prefix,'-',prefix,'-class',this.cluster,'-cellsubtype-predictions.txt')))) {
		this.subtypes = read.delim(file.path('stats/clusters/subintegration',paste0(glue.prefix,'-',prefix,'-class',this.cluster,'-cellsubtype-predictions.txt')),row.names=1)[,c('glue_subtype','glue_subtype_confidence')]
		this.subtypes$cell = rownames(this.subtypes)
	} else {
		this.subtypes = data.frame(glue_subtype=character(0),glue_subtype_confidence=numeric(0),cell=character(0))
	}
	this.subtypes
},mc.cores=n.cores))

this.umap$cell = factor(this.umap$cell,levels=rownames(this.umap))

this.umap = merge(this.umap,glue.predictions,by='cell',all.x=TRUE,all.y=FALSE)
this.umap = this.umap[order(this.umap$cell),]

cell.subtype.levels = scan(file.path('stats/subclusters','rna-cellsubclusters-levels.txt'),what='',sep='\n')

this.umap$glue_subtype = factor(this.umap$glue_subtype,levels=cell.subtype.levels)

this.umap$region_subclass = factor(this.umap$region,levels=region.hierarchy$region,labels=region.hierarchy$region.subclass)
this.umap$region_class = factor(this.umap$region,levels=region.hierarchy$region,labels=region.hierarchy$region.class)
this.umap$region_label = factor(this.umap$region,levels=region.hierarchy$region,labels=region.hierarchy$region.full)

# Temporary: correct error
levels(this.umap$glue_type) = gsub('^serotinergic neurons$','serotonergic neurons',levels(this.umap$glue_type))
levels(this.umap$glue_subtype) = gsub('^serotinergic neurons$','serotonergic neurons',levels(this.umap$glue_subtype))

levels(this.umap$glue_type) = gsub('^AHSG neuron$','AHSG neurons',levels(this.umap$glue_type))
levels(this.umap$glue_subtype) = gsub('^AHSG neuron$','AHSG neurons',levels(this.umap$glue_subtype))

levels(this.umap$glue_type) = gsub('^F5 neuron$','F5 neurons',levels(this.umap$glue_type))
levels(this.umap$glue_subtype) = gsub('^F5 neuron$','F5 neurons',levels(this.umap$glue_subtype))

levels(this.umap$glue_type) = gsub('^KIR3DL12 neuron$','KIR3DL12 neurons',levels(this.umap$glue_type))
levels(this.umap$glue_subtype) = gsub('^KIR3DL12 neuron$','KIR3DL12 neurons',levels(this.umap$glue_subtype))

this.umap$extraction.date = this.umap$extraction_bate

dobs = c('2C0' = '2006-12-22', '3R7' = '2013-08-19', '3I4' = '2009-11-08', '4I3' = '2009-09-17', '6J2' = '2009-09-10')
dods = c('2C0' = '2016-12-20', '3R7' = '2016-12-14', '3I4' = '2019-10-25', '4I3' = '2019-10-24', '6J2' = '2019-10-31')

dobs = do.call(c,lapply(dobs,as.Date))
dods = do.call(c,lapply(dods,as.Date))

this.umap$dob = dobs[this.umap$animal_id]
this.umap$dod = dods[this.umap$animal_id]

this.umap$age = with(this.umap,as.numeric(dod-dob)/365.2425)

this.umap$cell_class_manual = this.umap$cell_final
this.umap$cell_class_prediction = this.umap$glue_type
this.umap$cell_class_prediction_conf = this.umap$glue_type_confidence

this.umap$cell_subcluster_prediction = this.umap$glue_subtype
this.umap$cell_subcluster_prediction_conf = this.umap$glue_subtype_confidence

this.umap$total_fragments = this.umap$total
this.umap$total_dedup_fragments = this.umap$umi
this.umap$total_binarized_fragments = this.umap$umi_binarized

this.umap = within(this.umap,{
	id = factor(id)
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
})

this.umap = this.umap[,c(
'id',                               # library ID
'cell',                             # cell ID
'Size_Factor',                      # size factor
'total_fragments',                  # read counts before deduplication  (2x fragment count)
'total_dedup_fragments',            # read counts after deduplication  (2x deduplicated fragment count)
'RIP',                              # Sum reads in peaks
'RIT',                              # Sum reads in TSS
'FRIP',                             # Fraction reads in peaks
'FRIT',                             # Fraction reads in TSS
'total_binarized_fragments',        # number of transposition sites that intersect the feature after binarizing the input (peak) matrix
'n.umi',                            # number of transposition sites that intersect the feature after filtering the binarized input matrix
'biccn_id',                         # Tissue sample ID
'region',                           # brain region (abbreivation)
'region_label',                     # brain region (full name)
'region_class',                     # Broadest-level brain anatomical unit
'region_subclass',                  # Intermediate-level brain anatomical unit
# 'animal_integer',                 # integer assigned to animal (dropped)
'hemisphere',                       # brain hemisphere (NAs indicate midline structures)
'animal_id',                        # Cayo Santiago donor animal ID
'social_group',                     # Cayo Santiago social group ID
'age',                              # Age (years)
'sex',                              # Sex
'extractor',                        # Nuclei extracted by
'isolation_site',                   # Nuclei isolation location
# 'organization_order',             # Order in which samples were organized (drop)
# 'dissection_order',               # Order in which samples were dissected (drop)
# 'run_id',                         # sci-*-seq ID
'isolation_order',                  # Order isolated
'extraction_batch',                 # Extraction batch ID
'extraction.date',                  # Extraction date
'fixation',                         # Fixation method
# 'fixation_volume',                # Volume fixation
'lysis',                            # Lysis method
# 'freezing_medium',                # Medium in which isolated nuclei were flash-frozen (drop)
'sequencing_run_id',                # Sequencing batch
'doublet_score',                    # Doublet score (Reran Scrublet)
# 'predicted_doublet',              # whether predicted doublet (all FALSE)
# 'manual_doublet',                 # whether manual doublet
# 'is_singlet',                     # Whether singlet
# 'nb_features',                      # number of features (after filtering features)
'umap.1',                           # Global UMAP dimension 1
'umap.2',                           # Global UMAP dimension 2
'partition',                        # Global UMAP partition
'cluster',                          # Global UMAP cluster (resolution = 1e-5)
'cell_class_manual',                # Cell types called from visualization of gene activity scores
'cell_class_prediction',            # Cell types predicted from global glue integration
'cell_class_prediction_conf',       # Cell type prediction score from global glue integration
'cell_subcluster_prediction',       # Cell subtypes predicted from global glue integration
'cell_subcluster_prediction_conf'   # Cell subtype prediction score from global glue integration
)]

# write.table(this.umap[,c('cell','subcluster_manual')],file=file.path('stats/subclusters',paste0(prefix,'-subclusters.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
# write.table(this.umap[,c('cell','cell_subcluster_manual')],file=file.path('stats/subclusters',paste0(prefix,'-cellsubclusters.txt')),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
# write(levels(this.umap$cell_subcluster_manual),file=file.path('stats/subclusters',paste0(prefix,'-cellsubclusters-levels.txt')),sep='\n')

# Correct regions
region.corrections = data.frame(region=c('dmPFC','vmPFC','dlPFC','vlPFC','ACC','CC','CN','NAc','EC','PC','A1','AMY','HIP','M1','mdTN','vlTN','LGN','S1','IPP','SPP','STS','MT','IT','V1','CV','lCb','MB','MdO','MdC','Pons'),region_label=c('dorsomedial prefrontal cortex','ventromedial prefrontal cortex','dorsolateral prefrontal cortex','ventrolateral prefrontal cortex','anterior cingulate cortex','corpus callosum','head of the caudate nucleus','nucleus accumbens','entorhinal cortex','perirhinal cortex','primary auditory cortex','amygdala','hippocampus','primary motor cortex','mediodorsal thalamic nucleus','ventrolateral thalamic nucleus','lateral geniculate nucleus','primary somatosensory cortex','inferior posterior parietal gyrus','superior posterior parietal gyrus','superior temporal sulcus','middle temporal visual area','inferior temporal gyrus','primary visual cortex','cerebellar vermis','lateral cerebellar cortex','midbrain','open medulla','closed medulla','pons'))

this.umap = within(this.umap,{
	region = factor(region,levels=region.corrections$region)
	region_label = factor(region_label,levels=region.corrections$region_label)
})

saveRDS(this.umap,file=file.path('umap',paste0(glue.prefix,'-',prefix,'-glue-subintegration-classified.rds')))