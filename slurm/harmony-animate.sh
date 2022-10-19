#!/bin/bash
#SBATCH --job-name="harmony"
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err

module load r/4.0.2

# Third argument is number of iterations
Rscript scripts/harmony-animate.R $1 $2 $3 $4

# sbatch -p $WILDFIRE -q wildfire --exclusive slurm/harmony-animate.sh atac atac 2 cluster
# sbatch -p $WILDFIRE -q wildfire --exclusive slurm/harmony-animate.sh atac atac 2 region_major
# sbatch -p $WILDFIRE -q wildfire --exclusive slurm/harmony-animate.sh atac atac 2 region
# sbatch -p $WILDFIRE -q wildfire --exclusive slurm/harmony-animate.sh atac atac 2 sex
# sbatch -p $WILDFIRE -q wildfire --exclusive slurm/harmony-animate.sh atac atac 2 hemisphere
# sbatch -p $WILDFIRE -q wildfire --exclusive slurm/harmony-animate.sh atac atac 2 id
# sbatch -p $WILDFIRE -q wildfire --exclusive slurm/harmony-animate.sh atac atac 2 partition
# sbatch -p $WILDFIRE -q wildfire --exclusive slurm/harmony-animate.sh atac atac 2 sequencing_run_id
