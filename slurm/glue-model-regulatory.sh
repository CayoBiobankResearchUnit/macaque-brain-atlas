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

scripts/glue-model-regulatory.py $1 $2 $3

scripts/glue-model-regulatory-inference.py $1 $2 $3
