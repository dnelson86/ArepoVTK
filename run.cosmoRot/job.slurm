#!/bin/sh
##SBATCH --mail-user dnelson@cfa.harvard.edu
##SBATCH --mail-type=fail
#SBATCH -p itc_cluster
#SBATCH -J render
#SBATCH -o run.txt
#SBATCH -e run.txt
#SBATCH -t 0-12 # 12 hours
#SBATCH --mem=22000
##SBATCH --mem-per-cpu=3900 # MB per cpu
##SBATCH --exclusive
##SBATCH --dependency=afterany:731591 # run after job xx
#SBATCH --ntasks 16
#SBATCH --ntasks-per-node=16 # ptile

./ArepoRT config.txt >> run.txt
