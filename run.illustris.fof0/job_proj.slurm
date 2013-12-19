#!/bin/sh
##SBATCH --mail-user dnelson@cfa.harvard.edu
##SBATCH --mail-type=fail
#SBATCH -p itc_cluster
#SBATCH -J render
#SBATCH -o run.txt
#SBATCH -e run.txt
#SBATCH -t 360 # 6 hour in min
#SBATCH --mem=115200 # 3600MB/core for 32 cores
##SBATCH --mem-per-cpu=3900 # MB per cpu
##SBATCH --exclusive
##SBATCH --dependency=afterany:731591 # run after job xx
#SBATCH --ntasks 32
#SBATCH --ntasks-per-node=32 # ptile

#I455 (boxcen 1860,1860,1860) (rvir 930)
mpirun --mca mpi_leave_pinned 0 -n 32 ./Arepo_proj param_proj.txt 5 123 1440 1440 0 1 2 40 3680 40 3680 40 3680 >> run.txt

# CHANGE IC FILE IN PARAM_PROJ.TXT
# CHANGE BOXSIZE IN PARAM_PROJ.TXT

#I1820 123 (boxcen 1840,1840,1840) (rvir 920)
# SLICEWIDTH=2
#mpirun --mca mpi_leave_pinned 0 -n 32 ./Arepo_proj param_proj.txt 5 123 1440 1440 0 1 2 20 3660 20 3660 20 3660
# SLICEWIDTH=1
#mpirun --mca mpi_leave_pinned 0 -n 32 ./Arepo_proj param_proj.txt 5 123 1440 1440 0 1 2 20 3660 20 3660 920 2760
# SLICEWIDTH=0.5
#mpirun --mca mpi_leave_pinned 0 -n 32 ./Arepo_proj param_proj.txt 5 123 1440 1440 0 1 2 20 3660 20 3660 1380 2300

#I1820 135 sf20 (boxcen 1780,1780,1780) (rvir 890)
# SLICEWIDTH=2
#mpirun --mca mpi_leave_pinned 0 -n 32 ./Arepo_proj param_proj.txt 5 135 2048 2048 0 1 2 20 3540 20 3540 20 3540
# SLICEWIDTH=1
#mpirun --mca mpi_leave_pinned 0 -n 32 ./Arepo_proj param_proj.txt 5 135 2048 2048 0 1 2 20 3540 20 3540 890 2670
# SLICEWIDTH=0.5
#mpirun --mca mpi_leave_pinned 0 -n 32 ./Arepo_proj param_proj.txt 5 135 2048 2048 0 1 2 20 3540 20 3540 1335 2225

#I1820 135 sf20 1rvir cube, highres for composite image
#mpirun --mca mpi_leave_pinned 0 -n 64 ./Arepo_proj param_proj.txt 5 135 3000 3000 0 1 2 890 2670 890 2670 890 2670

#I1820 135 sf150 (boxcen 6675 6675 6675) (boxsize 13350) (rvir 890)
# same 2rvir x,y dimensions as above, but 15rvir deep
#mpirun --mca mpi_leave_pinned 0 -n 64 ./Arepo_proj param_proj.txt 5 135 2048 2048 0 1 2 4915 8435 4915 8435 20 13330
