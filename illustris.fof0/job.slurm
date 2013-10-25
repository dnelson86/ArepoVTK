#!/bin/sh
##SBATCH --mail-user dnelson@cfa.harvard.edu
##SBATCH --mail-type=fail
#SBATCH -p itc_cluster
#SBATCH -J rend1820d
#SBATCH -o run_1820_innerDens.txt
#SBATCH -e run_1820_innerDens.txt
#SBATCH -t 5760 # 4 days in min
#SBATCH --mem=90000
##SBATCH --mem-per-cpu=3900 # MB per cpu
##SBATCH --exclusive
##SBATCH --dependency=afterany:731591 # run after job xx
#SBATCH --ntasks 32
#SBATCH --ntasks-per-node=32 # ptile

./ArepoRT_sphInner config_1820_innerDens.txt >> run_1820_innerDens.txt
