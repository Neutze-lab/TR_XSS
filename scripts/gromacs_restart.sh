#!/usr/bin/env csh
#SBATCH -A C3SE408-22-1
#SBATCH -p vera
#SBATCH -N 4
#SBATCH -t 25:00:00
set n_cores=120

ml "GCC/8.3.0"
ml OpenMPI/3.1.4
ml GROMACS/2019.4

set run_name=$1

set cnt    = 1
set cntmax = $2

mpirun -np ${n_cores} gmx_mpi mdrun -s $run_name"_"${cnt}.tpr -v -cpi $run_name"_"${cnt}_prev.cpt -deffnm $run_name"_"${cnt} -append
