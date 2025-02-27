#!/bin/bash
#SBATCH --job-name seir              
#SBATCH --partition partition-name             
#SBATCH --ntasks 1                    # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 2-00:00                # Time limit in D-HH:MM
#SBATCH --mem-per-cpu=4000             # Memory limit for each tasks (in MB)
#SBATCH --array=1-51                 # The number of iterations for this job
#SBATCH -o messages/seir%a.out            # File to which STDOUT will be written
#SBATCH -e messages/seir%a.err            # File to which STDERR will be written

source /etc/profile.d/z00_lmod.sh
module load R
Rscript scripts/seir_state_cluster.R $SLURM_ARRAY_TASK_ID