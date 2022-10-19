#!/bin/bash
#SBATCH --job-name="integrate"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err

module load r/4.0.2

rna_prefix=$1
atac_prefix=$2
ncells=$3

Rscript --min-vsize=1040G --min-nsize=40M  scripts/seurat-integrate.R $rna_prefix $atac_prefix $ncells