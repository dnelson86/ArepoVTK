#!/bin/sh
#SBATCH --mail-user dnelson@cfa.harvard.edu
#SBATCH --mail-type=fail
#SBATCH -p itc_cluster
#SBATCH -J densRot
#SBATCH -o frames_rotDens/run.txt
#SBATCH -e frames_rotDens/run.txt
#SBATCH -t 5760 # 4 day in min
#SBATCH --mem=50000
##SBATCH --mem-per-cpu=3900 # MB per cpu
#SBATCH --exclusive
##SBATCH --dependency=afterany:731591 # run after job xx
#SBATCH --ntasks 64
#SBATCH --ntasks-per-node=64 # ptile

./ArepoRT config_rotDens.txt >> frames_rotDens/run.txt
