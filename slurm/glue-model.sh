#!/bin/bash
#SBATCH --job-name="glue"
#SBATCH --mem=0
#SBATCH --gres=gpu:4
#SBATCH --constraint=V100
#SBATCH --time=7-00:00:00
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err

module load python/3.7.1
module load bedtools/2.30.0

scripts/glue-model.py $1 $2 $3

module purge
module load r/4.0.2

scripts/glue-monocle.R $1

module purge
module load python/3.7.1
module load bedtools/2.30.0

scripts/glue-regulatory inference.py $1 $2 $3
