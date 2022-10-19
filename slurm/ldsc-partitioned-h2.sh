#!/bin/bash
#SBATCH --job-name="ldsc"
#SBATCH --time=7-00:00:00
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load perl/5.26.0
module load parallel/20180922

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

bed_prefix=$1

module load python/3.7.1

sumstats_file=sumstats_liberal.txt

n_phenotypes=$(tail -n+2 data/1000G/${sumstats_file} | wc -l | xargs)

parallel -j $slots scripts/ldsc-partitioned-h2.sh {1} {2} {3} ::: $bed_prefix ::: $SLURM_ARRAY_TASK_ID ::: $(seq 1 $n_phenotypes)
