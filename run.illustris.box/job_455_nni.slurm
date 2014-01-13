#!/bin/sh
##SBATCH --mail-user dnelson@cfa.harvard.edu
##SBATCH --mail-type=fail
#SBATCH -p hernquist
#SBATCH --array=0-15
#SBATCH -J ren455n
#SBATCH -o output/455_2k_nni/run_%a.txt
#SBATCH -e output/455_2k_nni/run_%a.txt
#SBATCH -t 0-12 # 12 hours
#SBATCH --mem=58000
##SBATCH --exclusive
##SBATCH --dependency=afterany:731591 # run after job xx
#SBATCH --ntasks 16
#SBATCH --ntasks-per-node=16 # ptile

./ArepoRT_nni config_455.txt -j ${SLURM_ARRAY_TASK_ID} >> output/455_2k_nni/run_${SLURM_ARRAY_TASK_ID}.txt
