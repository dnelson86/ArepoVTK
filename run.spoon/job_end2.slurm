#!/bin/sh
#SBATCH --mail-user dnelson@cfa.harvard.edu
#SBATCH --mail-type=fail
#SBATCH -p itc_cluster
#SBATCH -J rendEnd2
#SBATCH -o frames_end2/run.txt
#SBATCH -e frames_end2/run.txt
#SBATCH -t 2880 # 2 day in min
#SBATCH --mem=220000
##SBATCH --mem-per-cpu=3900 # MB per cpu
#SBATCH --exclusive
##SBATCH --dependency=afterany:731591 # run after job xx
#SBATCH --ntasks 64
#SBATCH --ntasks-per-node=8 # ptile

./ArepoRT config_end2.txt >> frames_end2/run.txt
