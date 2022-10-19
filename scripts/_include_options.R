#!/usr/bin/env Rscript
# 🐛

# Options
options(future.globals.maxSize = 100000 * 1024^2)

# Number of cores to use for parallel computing (set to 0 to use all detected cores)
n.cores = 24

# Set region levels
region.levels = c('dmPFC','vmPFC','dlPFC','vlPFC','ACC','CC','CN','NAc','EC','PC','A1','AMY','HIP','M1','mdTN','vlTN','LGN','S1','IPP','SPP','STS','IT','V1','CV','lCb','MB','MdO','MdC','Pons')

region.colors = c("#82ED3B", "#D7B3E8", "#8A37ED", "#6E5EDA", "#7BD3E7", "#DFE194",
"#DD5664", "#8B9159", "#E38D56", "#E68FC0", "#DCA79C", "#E2BD4E",
"#E3DFBE", "#6D93A3", "#897FD2", "#D147DA", "#5EE375", "#AFDF6A",
"#DAE241", "#D9D3DF", "#73EEE0", "#80B1E7", "#DA499E", "#BEE9E0",
"#D077DC", "#5EE7AE", "#6DB6A3", "#8F6078", "#9EDC9F")

cell.type.colors = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666', '#ed1c24', '#00aeef', '#86328c', '#00aaad', '#ac7eaf', '#2e3192', '#809ead', '#6a3e14')

# Set region full names
region.full = c('dorsomedial prefrontal cortex','ventromedial prefrontal cortex','dorsolateral prefrontal cortex','ventrolateral prefrontal cortex','anterior cingulate cortex','corpus callosum','head of the caudate nucleus','nucleus accumbens','entorhinal cortex','perirhinal cortex','primary auditory cortex','amygdala','hippocampus','primary motor cortex','mediodorsal thalamic nucleus','ventrolateral thalamic nucleus','lateral geniculate nucleus','primary somatosensory cortex','inferior posterior parietal gyrus','superior posterior parietal gyrus','superior temporal sulcus','inferior temporal gyrus','primary visual cortex','cerebellar vermis','lateral cerebellar cortex','midbrain','open medulla','closed medulla','pons')
names(region.full) = region.levels

region.hierarchy = data.frame(
	region = factor(
		region.levels,
		levels = region.levels
	),
	region.full = factor(
		region.full,
		levels = region.full
	),
	region.class = factor(
		c('cortical','cortical','cortical','cortical','cortical','subcortical','subcortical','subcortical','cortical','cortical','cortical','subcortical','subcortical','cortical','subcortical','subcortical','subcortical','cortical','cortical','cortical','cortical','cortical','cortical','cerebellum','cerebellum','brainstem','brainstem','brainstem','brainstem'),
		levels = c('cortical','subcortical','cerebellum','brainstem')
	),
	region.subclass = factor(
		c('frontal lobe','frontal lobe','frontal lobe','frontal lobe','frontal lobe','corpus callosum','basal ganglia','basal ganglia','temporal lobe','temporal lobe','temporal lobe','hippocampus/amygdala','hippocampus/amygdala','frontal lobe','thalamus','thalamus','thalamus','parietal lobe','parietal lobe','parietal lobe','temporal lobe','temporal lobe','occipital lobe','cerebellum','cerebellum','brainstem','brainstem','brainstem','brainstem'),
		levels = c('frontal lobe','temporal lobe','parietal lobe','occipital lobe','corpus callosum','hippocampus/amygdala','thalamus','basal ganglia','cerebellum','brainstem')
	)
)

ens.species = 'mmulatta'

ens.version = 101

full.species = 'Macaca mulatta'

genome.code = 'Mmul_10'

genome.short = 'mmul'

sex.colors = c('#4daf4a','#984ea3') # (colors based on palette PRGn)
hemisphere.colors = c('#bf812d','#35978f','') # (colors based on palette BrBG)

# If using a reference-based mapping, identify the batch (sequencing_run_id level) that is the reference
# reference.batch = c('Snyder-Mackler_RNA3-019','Snyder-Mackler_RNA3-026') # 'Snyder-Mackler_RNA3-036'
reference.batch = 'Snyder-Mackler_RNA3-036'

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
#                             End configurations
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #

# If not set, sets n.cores to detected cores
if (!n.cores) n.cores = ifelse('future' %in% rownames(installed.packages()),future::availableCores(methods='mc.cores') + 1,parallel::detectCores(logical=FALSE))
