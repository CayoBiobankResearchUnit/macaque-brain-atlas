# cayo_sc_brain


# # # # # # # # # # # # # # # # # # # # # # #
#               ATAC pipeline               #
# # # # # # # # # # # # # # # # # # # # # # #

# Merge peaks using unified peak set
sbatch --mem=16G --array=1-110 scripts/atac-merge-peaks.sh

# Convert monocle cds to scanpy objects
sbatch --array=1-110 slurm/monocle-to-scanpy.sh atac atac

# Visualize scrublet distribution
scripts/scrublet-visualize.R atac atac

# Edit data/atac_scrublet_parameters.txt file to set manual thresholds. Then plot again
scripts/scrublet-visualize.R atac atac

# Process each sample individually (for QC)
sbatch --array=1-110 slurm/scanpy-single.sh atac atac

# Concatenate all scanpy objects together
sbatch slurm/scanpy-concat.sh atac atac

# Preprocess
sbatch --mem=800G slurm/scanpy-preprocess.sh atac atac

# The preprocess automatically identifies and labels doublet clusters. Now remove and reprocess
sbatch  --mem=700G slurm/scanpy-postprocess.sh atac atac

# Classify manually (this will be superseded later — just a preliminary exercise)
scripts/r-celltypes-classify.R atac atac

# Calculate promoter accessibility score (file name uses gene activity but this is not the best term)
sbatch --array=1-110 --mem=16G slurm/muon-gene-activity.sh atac atac

# Prep marker genes and plot
# Concatenate gene activity scores into one combined annotated dataset
scripts/muon-geneactivity-marker-genes-concat.py atac atac

# Normalize gene activity scores and prep for plotting
scripts/muon-geneactivity-marker-genes-prep.py atac atac

# Plot marker genes
scripts/scanpy-marker-genes-plot.py atac atac

# Call peak motifs (break in 25 chunks)
sbatch --array=1-25 slurm/atac-call-peak-motifs.sh atac

# Make peak-TF matrix
sbatch slurm/atac-generate-motif-matrix.sh atac

# Make cell-motif matrix
sbatch --array=1-110 --mem=24G slurm/atac-generate-cell-motif-matrix.sh atac

# Concatenate cell-motif matrix
sbatch slurm/atac-concat-cell-motif-matrix.sh atac


# # # # # # # # # # # # # # # # # # # # # # #
#               RNA pipeline                #
# # # # # # # # # # # # # # # # # # # # # # #

#    *    *    *    *    *    *    *    *
# Main pipeline (unified across all data)
#    *    *    *    *    *    *    *    *

sbatch --array=1-139 slurm/monocle-to-scanpy.sh rna rna

scripts/scrublet-visualize.R rna rna

# Edit data/rna_scrublet_parameters.txt file to set manual thresholds. Then plot again
scripts/scrublet-visualize.R rna rna

# Run QC
sbatch --array=1-139 --exclusive --mem=0 slurm/scanpy-single.sh rna rna

sbatch --mem=500G slurm/scanpy-concat.sh rna rna

# Preprocess
sbatch --exclusive --mem=1T slurm/scanpy-preprocess.sh rna rna

# The preprocess automatically identifies and labels doublet clusters. Now remove and reprocess
sbatch --mem=700G --time=4-00:00:00 slurm/scanpy-postprocess.sh rna rna

# Classify manually
scripts/r-celltypes-classify.R rna rna

# Summarize classifications
scripts/r-celltypes-summarize.R rna rna

# Prep and plot marker genes (python script highlights markers in UMAP, R script visualizes top markers)
scripts/scanpy-marker-genes-prep.py rna rna
scripts/scanpy-marker-genes-plot.py rna rna
scripts/r-marker-genes-plot.R rna rna

## Quick aside to detect and remove contamination

# Test for exogenous contamination (mapping to three genomes simultaneously)
sbatch --mem=150G slurm/bbsplit-index.sh macaque human mouse
sbatch --array=1-$(cat data/rna_metadata.txt | tail -n+2 | wc -l) --mem=150G --cpus-per-task=24 slurm/bbsplit-run.sh macaque human mouse

# Parse bbsplit results
sbatch --array=1-$(cat data/rna_metadata.txt | tail -n+2 | wc -l) --mem=8G --time=00:30:00 slurm/py-bbsplit-parse.sh rna macaque human mouse

# Redo marker genes
sbatch --exclusive --mem=0 --time=3-00:00:00 slurm/scanpy-marker-genes.sh rna rna

# Redo dendrograms
sbatch --mem=64G --array=1-20 slurm/scanpy-dendrogram-bootstrap.sh rna rna 50






#    *    *    *    *    *    *    *    *
# Subcluster RNA ("recluster" pipeline)
#    *    *    *    *    *    *    *    *

# Prep subclustering
sbatch --time=1-00:00:00 --mem=200G --array=1-$(wc -l stats/clusters/rna-final-cellclasses-levels.txt | cut -d' ' -f1 | xargs) slurm/scanpy-recluster-prep.sh rna rna

# Recluster
sbatch --time=1-00:00:00 --array=1-17 --mem=750G slurm/scanpy-recluster-process.sh rna rna

# The script below was useful for exploring reclustering/integrating methods but ultimately dropped
# slurm/monocle-recluster-full.sh

# Finalize
scripts/r-recluster-finalize.R rna rna

# Script to visualize marker gene expression, here is just one example run (class 3 = GABAergic neurons)
sbatch --time=1-00:00:00 --mem=200G slurm/scanpy-recluster-plot.sh rna rna 3 u,t,v,h,d SST,PVALB,PAX6,VIP,SNCG,LAMP5,ADARB2,LHX6,GAD1,GAD2


#    *    *    *    *    *    *    *    *
# Cluster identification using reference data
#    *    *    *    *    *    *    *    *

# Download datasets
scripts/dataset-download.sh

# Preprocess (pseudobulk) query dataset (both subtypes as well as parent cell classes)
sbatch slurm/dataset-preprocess.sh rna
sbatch slurm/dataset-preprocess.sh rna-class

# Preprocess (pseudobulk) each reference dataset
sbatch slurm/dataset-preprocess.sh early-brain
sbatch slurm/dataset-preprocess.sh brain-vasc-atlas
sbatch slurm/dataset-preprocess.sh mouse-drg-injury
sbatch slurm/dataset-preprocess.sh vascular-dev
sbatch slurm/dataset-preprocess.sh dev-inhibitory-neurons-macaque
sbatch slurm/dataset-preprocess.sh dev-brain-regions

# Run NNLS on subtypes
sbatch slurm/dataset-correlate.sh rna cortex-dev
sbatch slurm/dataset-correlate.sh rna human-cortex
sbatch slurm/dataset-correlate.sh rna adult-brain-vasc-peri
sbatch slurm/dataset-correlate.sh rna adult-brain-vasc-endo
sbatch slurm/dataset-correlate.sh rna vascular-dev
sbatch slurm/dataset-correlate.sh rna early-brain
sbatch slurm/dataset-correlate.sh rna brain-vasc-atlas
sbatch slurm/dataset-correlate.sh rna mouse-drg-injury
sbatch slurm/dataset-correlate.sh rna dev-inhibitory-neurons-macaque
sbatch slurm/dataset-correlate.sh rna myeloid-neuroinflam
sbatch slurm/dataset-correlate.sh rna fang-merfish-l2
sbatch slurm/dataset-correlate.sh rna fang-merfish-l3

# Run NNLS on cell classes
sbatch slurm/dataset-correlate.sh rna-class cortex-dev
sbatch slurm/dataset-correlate.sh rna-class human-cortex
sbatch slurm/dataset-correlate.sh rna-class adult-brain-vasc-peri
sbatch slurm/dataset-correlate.sh rna-class adult-brain-vasc-endo
sbatch slurm/dataset-correlate.sh rna-class vascular-dev
sbatch slurm/dataset-correlate.sh rna-class early-brain
sbatch slurm/dataset-correlate.sh rna-class brain-vasc-atlas
sbatch slurm/dataset-correlate.sh rna-class brain-vasc-atlas
sbatch slurm/dataset-correlate.sh rna-class mouse-drg-injury
sbatch slurm/dataset-correlate.sh rna-class dev-inhibitory-neurons-macaque
sbatch slurm/dataset-correlate.sh rna-class myeloid-neuroinflam
sbatch slurm/dataset-correlate.sh rna-class fang-merfish-l2
sbatch slurm/dataset-correlate.sh rna-class fang-merfish-l3

# Plot NNLS
for i in human-cortex adult-brain-vasc-peri adult-brain-vasc-endo brain-vasc-atlas myeloid-neuroinflam fang-merfish-l2 fang-merfish-l3; do
scripts/dataset-plot.R rna ${i}
scripts/dataset-plot.R rna-class ${i}
done

for i in human-cortex fang-merfish-l2 fang-merfish-l3; do
scripts/dataset-plot.R rna ${i}
done

# # # # # # # # # # # # # # # # # # # # # # #
#          JOINT multiomic pipeline         #
# # # # # # # # # # # # # # # # # # # # # # #

# GLUE pipeline, preprocess 

sbatch --mem=1T slurm/glue-preprocess.sh rna rna
sbatch --mem=1T slurm/glue-preprocess.sh atac atac

# Run on 1.5 T node
sbatch --gres=gpu:4 --constraint=V100 --exclusive --mem=0 slurm/glue-model.sh biccn rna atac

# Prepare GLUE label transfer pipeline
sbatch --exclusive --mem=0 slurm/glue-transfer-labels-prep.sh biccn rna atac

# Make kNN neighbor graphs, run in parallel instead of in sequence.
sbatch --array=1-10 --exclusive --mem=0 slurm/glue-transfer-labels-knn.sh biccn rna atac

# Transfer labels
sbatch --array=1-3 --exclusive --mem=0 slurm/glue-transfer-labels-map.sh biccn rna atac

# Run GLUE regulatory inference
sbatch slurm/glue-model-regulatory.sh biccn rna atac

# Make metacells using the GLUE cell embeddings
sbatch --mem=500G slurm/glue-metacells.sh biccn rna atac

# Copy over metacell regression results
cat /scratch/nsnyderm/atac_rna/allcells_part{1..60}_gene-peak.tsv | grep -v Error | grep -v subscript | grep -v '^$' | cut -f 1-6 > stats/metacell_lr/gene_peak_all.tsv
cat /scratch/nsnyderm/atac_rna/allcells_part{1..60}_perm_gene-peak.tsv | grep -v Error | grep -v subscript | grep -v '^$' | cut -f 1-6 > stats/metacell_lr/gene_peak_null_all.tsv

# Summarize regulatory results
scripts/r-regulatory-summarize.R biccn atac

# Calculate marker peaks
sbatch --mem=1T slurm/scanpy-marker-peaks.sh atac atac
## Calculate peaks the seurat way
sbatch --mem=500G slurm/seurat-marker-peaks-prep.sh atac atac

# For Seurat, parallelize by calculating marker peaks separately for each cell class
sbatch --array=1-12 --exclusive --mem=0 slurm/seurat-marker-peaks-run.sh atac atac

# Combine Seurat and Scanpy files
sbatch --mem=128G --time=12:00:00 slurm/r-marker-peaks-plot.sh atac atac

# Prep LDSC files
scripts/ldsc-prep-inputs.sh

# # # # # # # # # # # # # # # # # # # # # # #
#              Re-called peaks              #
# # # # # # # # # # # # # # # # # # # # # # #

# For cell-class specific analyses, peak-calling was repeated on cell-class-specific data

# Generate cell-type specific transposition sites, fragments, and whitelists
sbatch --array=1-110 --mem=200G slurm/atac-get-cell-type-fragments.sh
sbatch --array=1-110 --mem=200G slurm/atac-get-cell-type-transposition-sites.sh
sbatch --array=1-110 --mem=4G slurm/atac-get-cell-type-whitelists.sh

# Use MACS3 to call subpeaks on the split fragments files
sbatch --array=1-3,6-13,15 slurm/macs3-call-subpeaks.sh

# Divvy up the jobs by size
paste <(seq 1 $(ls transposition-sites-cells/*.bed.gz | wc -l)) \
<(/usr/bin/ls -lS transposition-sites-cells/*.bed.gz | \
	tr -s ' ' | cut -d ' ' -f 5,9 | \
	sed 's/\([0-9]*\) .*\?\(NSM[0-9]\{3\}\)-class\([0-9]*\)\.bed\.gz/\2'$'\t''\3'$'\t''\1/g') > \
	data/transposition_file_sizes.txt

sbatch slurm/atac-merge-subpeaks-parallel.sh 1 1172

# The script below is effective at running stragglers (jobs that failed)
# sbatch --mem=200G --array=1 slurm/atac-merge-subpeaks-single.sh

sbatch slurm/scanpy-subpeaks-prep-parallel.sh 1 1172
# sbatch --array=1-20 --mem=64G slurm/scanpy-subpeaks-prep-single.sh

# Concat all samples together (within a given cell class) using new peaks
sbatch --array=1-3,6-13,15 slurm/scanpy-subpeaks-concat.sh

# The above generates a new prefix "atacsub". Now run GLUE integration
sbatch --array=1-3,6-13,15 --exclusive --mem=0 slurm/glue-recluster-model.sh subpeak rna atacsub

########
# Predictions

# Make a redundant copy of cell classes and levels
ln stats/clusters/atac-final-cellclasses.txt stats/clusters/atacsub-final-cellclasses.txt
ln stats/clusters/atac-final-cellclasses-levels.txt stats/clusters/atacsub-final-cellclasses-levels.txt
ln data/atac_metadata.txt data/atacsub_metadata.txt

# Transfer labels
sbatch --array=1-3,6-13,15 slurm/glue-recluster-transfer-labels-prep.sh subpeak rna atacsub

for i in {1..3} {6..13} 15; do
sbatch --array=1-7 --mem=50G slurm/glue-recluster-transfer-labels-knn.sh subpeak rna atacsub $i
done

for i in {1..3} {6..13} 15; do
sbatch --array=1-2 --mem=50G slurm/glue-recluster-transfer-labels-map.sh subpeak rna atacsub $i
done

# Make metacells
sbatch --array=1-3,6-13,15 --mem=200G slurm/glue-recluster-metacells.sh subpeak rna atacsub

# Again, pull the metacell LR results from elsewhere
mkdir -p stats/metacell_lr

for i in {1..3} {6..13} 15; do
echo $i
cat /scratch/nsnyderm/atac_rna/class${i}_part{1..12}_gene-peak.tsv | grep -v Error | grep -v subscript | grep -v '^$' | cut -f 1-6 > stats/metacell_lr/gene_peak_class${i}.tsv
cat /scratch/nsnyderm/atac_rna/class${i}_part{1..12}_perm_gene-peak.tsv | grep -v Error | grep -v subscript | grep -v '^$' | cut -f 1-6 > stats/metacell_lr/gene_peak_null_class${i}.tsv
done

# Calculate marker peaks (first in scanpy)
# All done
sbatch --array=1-3,6-13,15 --mem=200G slurm/scanpy-recluster-marker-peaks.sh atacsub subpeak

# Now in seurat
sbatch --array=1-3,6-13,15  --mem=100G slurm/seurat-recluster-marker-peaks.sh atacsub subpeak

# Copy unique peaks (generated by Noah)
for i in {1..3} {6..13} 15; do
	echo class${i}
	cp ../../nsnyderm/bbi/atac_peaks/class${i}_peaks.bed.unique bed/unique_peaks_class${i}.bed
done

# copy all peaks
for i in {1..3} {6..13} 15; do
	echo class${i}
	cp bed/subpeaks/merged_peaks-class${i}.bed bed/all_peaks_class${i}.bed
done

# Create a list of cell class / subtype index combinations as a task list
rm -rf data/atacsub_seurat_cell_subclusters.txt
for i in {1..3} {6..13} 15; do
	column_number=$(zcat mm/recluster/seurat/atacsub_downsampled_class${i}_cell_metadata.txt.gz | head -n1 | sed -n '1 s/cell_subcluster.*//p' | sed 's/[^\t*]//g' | wc -c)
	for j in $(seq 1 $(zcat mm/recluster/seurat/atacsub_downsampled_class${i}_cell_metadata.txt.gz | tail -n+2 | cut -f ${column_number} | sort -u | wc -l)); do
#	for j in $(seq 1 $(zcat mm/recluster/seurat/atacsub_downsampled_class${i}_cell_metadata.txt.gz | tail -n+2 | cut -f 36 | sort -u | wc -l)); do
	echo ${i}$'\t'${j} >> data/atacsub_seurat_cell_subclusters.txt
	done
done

# Calculate marker peaks
sbatch --array=1-83 slurm/seurat-recluster-marker-peaks-run.sh atacsub subpeak

# Plot marker peaks and summarize regulatory inference results
sbatch --array=1-3,6-13,15 --mem=24G --time=12:00:00 slurm/r-recluster-marker-peaks-plot.sh atacsub atac
sbatch --array=1-3,6-13,15 --mem=8G --time=12:00:00 slurm/r-recluster-regulatory-summarize.sh subpeak atacsub


# Grab TF enrichment results

mkdir -p stats/tf/homer/markerpeaks/class
mkdir -p stats/tf/homer/markerpeaks/cluster
mkdir -p stats/tf/homer/regulatory/class
mkdir -p stats/tf/homer/regulatory/cluster

cp /scratch/nsnyderm/bbi/reg_peaks/class{{1..3},{6..13},15}_marker_jasparTFs.txt stats/tf/homer/markerpeaks/class
cp /scratch/nsnyderm/bbi/reg_peaks/class*_cluster*_marker_jasparTFs.txt stats/tf/homer/markerpeaks/cluster
cp /scratch/nsnyderm/bbi/reg_peaks/allpeaks_inc_v_dec_jasparTFs.txt stats/tf/homer/regulatory/class
cp /scratch/nsnyderm/bbi/reg_peaks/class{{1..3},{6..13},15}_*c_jasparTFs.txt stats/tf/homer/regulatory/cluster

# copy JASPAR names (to keep all TFs in sync
cp /scratch/nsnyderm/bbi/reg_peaks/jaspar2018_CORE_vertebrates_non-redundant.txt data/

# Summarize TF enrichment results
scripts/r-tf-summarize.R

# Run LDSC
# Liftover all peaks
sbatch --array=1-3,6-13,15 slurm/liftover-bed.sh all_peaks

# Calculate LD scores
sbatch --array=1-3,6-13,15 slurm/ldsc-ld-score.sh all_peaks

# Calculate partitioned h2
# 17322530
sbatch --array=1-3,6-13,15 slurm/ldsc-partitioned-h2.sh all_peaks



# Purkinje-specific analyses

sbatch --mem=500G slurm/atac-get-cell-subtype-fragments.sh atac atac 6

sbatch --mem=260G slurm/atac-get-cell-subtype-fragments.sh atac atac 16

# Call peaks at the cell-subtype level
sbatch --mem=32G slurm/macs3-call-sub-subpeaks.sh 6

# grep -n '^6'$'\t'16 data/atacsub_markerpeaks_cell_subclusters.txt
# sed -n 40p data/atacsub_markerpeaks_cell_subclusters.txt
sbatch --time=1:00:00 --array=40 slurm/liftover-recluster-bed.sh merged_peaks atacsub
sbatch --time=1:00:00 --array=40 slurm/liftover-recluster-bed.sh unique_peaks atacsub

# After liftover, some cell clusters have very few peaks, remove them and make an updated task list
wc -l bed/marker_peaks_hg19_class*_cluster*.bed | tr -s ' ' | grep -v 'total$' | sed -E 's:([0-9]+) bed/marker_peaks_hg19_class([0-9]+)_cluster([0-9]+).bed:\1\t\2\t\3:g' | awk '$1 > 10' | cut -f 2-3 | sort -n -k1,1 -k2,2 > data/atacsub_markerpeaks_pass_cell_subclusters.txt

# grep -n '^6'$'\t'16 data/atacsub_markerpeaks_pass_cell_subclusters.txt
# sed -n 35p data/atacsub_markerpeaks_pass_cell_subclusters.txt (purkinje is number 35 on the task list)
sbatch --exclusive --mem=0 --requeue --array=35 slurm/ldsc-recluster-ld-score.sh top_1p_peaks # 17600324
sbatch --exclusive --mem=0 --requeue --array=35 slurm/ldsc-recluster-partitioned-h2.sh top_1p_peaks # 17600325

# Make task list for monaLisa
ls bed/subsubpeaks/top_1p_peaks-class*_cluster*.bed | sed -E 's/.*class([0-9]+)_cluster([0-9]+).bed/\1\t\2/g' | sort -k1,1n -k2,2n > data/monalisa_subpeaks_list.txt
sbatch --exclusive --mem=0 --requeue slurm/monalisa-sub-subpeaks.sh top_1p_peaks 6 16  # 17600336

# Now plot regulatory links
scripts/r-recluster-regulatory-plot-aggregate.R subpeak atacsub ENSMMUG00000011333 0 # MBP

# Create a list of cell class / subtype index combinations as a task list FOR MARKER PEAKS
ls -1 stats/markerpeaks/bed/atacsub_class*_cluster*_marker_peaks.bed | sed -E 's:.*?bed/atacsub_class([0-9]+)_cluster([0-9]+)_marker_peaks.bed:\1\t\2:g' | sort -n -k1,1 -k2,2 > data/atacsub_markerpeaks_cell_subclusters.txt



# Some figures in the paper are generated in the (extreme unorganized) scripts below:
# scripts/_r-paper-figures.R
# Scripts/_py-paper-figures.py



