#!/usr/bin/env csh
#SBATCH -A C3SE408-22-1
#SBATCH -p vera
#SBATCH -N 6
#SBATCH -c 2
#SBATCH -t 25:00:00
set n_cores=120

ml "GCC/8.3.0"
ml OpenMPI/3.1.4
ml GROMACS/2020

set run_name=$1


###############Production#############

set cnt    = 1
set cntmax = $2

while ( ${cnt} <= ${cntmax} )
        gmx grompp -f step7_production_macro.mdp -o $run_name"_"${cnt}.tpr -c step6.6_equilibration.gro -p topol.top -r target_pdb -maxwarn 2
        mpirun -np $n_cores gmx_mpi mdrun -v -deffnm $run_name"_"${cnt} -x $run_name"_"${cnt}.xtc -ntomp 1
    @ cnt += 1
end
