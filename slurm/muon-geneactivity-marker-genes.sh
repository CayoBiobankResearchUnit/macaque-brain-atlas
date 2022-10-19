#!/bin/bash
#SBATCH --job-name="muon"
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load python/3.7.1

scripts/muon-geneactivity-concat.py atac atac

scripts/muon-geneactivity-marker-genes-prep.py atac atac

scripts/muon-geneactivity-marker-genes-plot.py atac atac
