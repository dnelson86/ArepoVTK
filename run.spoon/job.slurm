#!/bin/sh
##SBATCH --mail-user dnelson@cfa.harvard.edu
##SBATCH --mail-type=fail
#SBATCH -p itc_cluster
#SBATCH --array=472-575
#SBATCH -J render
#SBATCH -o frames_final/run_%a.txt
#SBATCH -e frames_final/run_%a.txt
#SBATCH -t 60 # 1 hour in min
#SBATCH --mem=22000
##SBATCH --mem-per-cpu=3900 # MB per cpu
##SBATCH --exclusive
##SBATCH --dependency=afterany:731591 # run after job xx
#SBATCH --ntasks 8
#SBATCH --ntasks-per-node=8 # ptile

./ArepoRT config_movie.txt ${SLURM_ARRAY_TASK_ID} >> frames_final/run_${SLURM_ARRAY_TASK_ID}.txt
