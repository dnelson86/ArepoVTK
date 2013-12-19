#!/bin/sh
##SBATCH --mail-user dnelson@cfa.harvard.edu
##SBATCH --mail-type=fail
#SBATCH -p itc_cluster
#SBATCH -J R1820t60T
#SBATCH -o run_1820_tree60T.txt
#SBATCH -e run_1820_tree60T.txt
#SBATCH -t 5760 # 4 days in min
#SBATCH --mem=240000
##SBATCH --mem-per-cpu=3900 # MB per cpu
#SBATCH --exclusive
##SBATCH --dependency=afterany:731591 # run after job xx
#SBATCH --ntasks 64
#SBATCH --ntasks-per-node=64 # ptile

./ArepoRT_sph config_1820_tree60T.txt >> run_1820_tree60T.txt
