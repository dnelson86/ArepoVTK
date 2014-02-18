#!/bin/sh
##SBATCH --mail-user dnelson@cfa.harvard.edu
##SBATCH --mail-type=fail
#SBATCH -p itc_cluster
#SBATCH --array=3,19,33
#SBATCH -J ren18n
#SBATCH -o output/1820_256_16k_nni/run_%a.txt
#SBATCH -e output/1820_256_16k_nni/run_%a.txt
#SBATCH -t 2- # 2 days
#SBATCH --mem=58000
##SBATCH --exclusive
##SBATCH --dependency=afterany:731591 # run after job xx
#SBATCH --ntasks 16
#SBATCH --ntasks-per-node=16 # ptile

./ArepoRT_nni config_1820_256_16k_nni.txt -j ${SLURM_ARRAY_TASK_ID} >> output/1820_256_16k_nni/run_${SLURM_ARRAY_TASK_ID}.txt
