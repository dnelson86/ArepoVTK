#!/bin/sh
##SBATCH --mail-user dnelson@cfa.harvard.edu
##SBATCH --mail-type=fail
#SBATCH -p hernquist
#SBATCH --array=0-63
#SBATCH -J ren18g
#SBATCH -o output/1820_256_16k_cellgrad/run_%a.txt
#SBATCH -e output/1820_256_16k_cellgrad/run_%a.txt
#SBATCH -t 1- # 1 days
#SBATCH --mem=58000
##SBATCH --exclusive
##SBATCH --dependency=afterany:731591 # run after job xx
#SBATCH --ntasks 16
#SBATCH --ntasks-per-node=16 # ptile

./ArepoRT_grad config_1820_256_16k_cellgrad.txt -j ${SLURM_ARRAY_TASK_ID} >> output/1820_256_16k_cellgrad/run_${SLURM_ARRAY_TASK_ID}.txt
