#!/bin/sh
#SBATCH --mail-user dnelson@cfa.harvard.edu
#SBATCH --mail-type=fail
#SBATCH -p itc_cluster
#SBATCH -J tempRot
#SBATCH -o frames_rotTemp/run.txt
#SBATCH -e frames_rotTemp/run.txt
#SBATCH -t 5760 # 4 day in min
#SBATCH --mem=50000
##SBATCH --mem-per-cpu=3900 # MB per cpu
#SBATCH --exclusive
##SBATCH --dependency=afterany:731591 # run after job xx
#SBATCH --ntasks 64
#SBATCH --ntasks-per-node=64 # ptile

./ArepoRT config_rotTemp.txt >> frames_rotTemp/run.txt
