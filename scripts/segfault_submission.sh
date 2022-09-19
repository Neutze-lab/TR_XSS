#!/usr/bin/env csh
#SBATCH -A C3SE408-22-1
#SBATCH -p vera
#SBATCH -N 20
#SBATCH -t 20:00:00
set n_cores=400

ml "GCC/7.3.0-2.30"
ml OpenMPI/3.1.1
ml GROMACS/2018.2

mpirun -np ${n_cores} gmx_mpi mdrun -s the_name_the_number.tpr -v -cpi the_name_the_number_prev.cpt -deffnm the_name_the_number -append
