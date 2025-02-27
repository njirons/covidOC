#!/bin/bash
#SBATCH --job-name npi-bhm              
#SBATCH --partition partition-name             
#SBATCH --ntasks 1                    # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 2-00:00                # Time limit in D-HH:MM
#SBATCH --mem-per-cpu=4000             # Memory limit for each tasks (in MB)
#SBATCH --array=1-100                 # The number of iterations for this job
#SBATCH -o messages/npi-bhm%a.out            # File to which STDOUT will be written
#SBATCH -e messages/npi-bhm%a.err            # File to which STDERR will be written

source /etc/profile.d/z00_lmod.sh
module load R
Rscript scripts/npi-bhm-ar1-nX1-cluster.R $SLURM_ARRAY_TASK_ID
Rscript scripts/npi-bhm-ar1-nX2-cluster.R $SLURM_ARRAY_TASK_ID
Rscript scripts/npi-bhm-ar1-nX3-cluster.R $SLURM_ARRAY_TASK_ID