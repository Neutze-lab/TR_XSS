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

##############Minimization##############

# In the case that there is a problem during minimization using a single precision of GROMACS, please try to use
# a double precision of GROMACS only for the minimization step.

# step6.0
gmx grompp -f step6.0_minimization.mdp -o step6.0_minimization.tpr -c step5_input_aligned.gro -r step5_input_aligned.gro -p topol.top
mpirun -np $n_cores gmx_mpi mdrun -v -deffnm step6.0_minimization



##############Equilibration#############

# Equilibration
set cnt = 1
set cntmax = 6

#if ( ! -f step6.6_equilibration.tpr ) then
echo Running 6 equilibration steps as output is not detected
        while ( ${cnt} <= ${cntmax} )
            @ pcnt = ${cnt} - 1
            if ( ${cnt} == 1 ) then
                gmx grompp -f step6.{$cnt}_equilibration.mdp -o step6.{$cnt}_equilibration.tpr -c step6.{$pcnt}_minimization.gro -r step5_input_aligned.gro -p topol.top -maxwarn -1
                mpirun -np $n_cores gmx_mpi mdrun -v -deffnm step6.{$cnt}_equilibration -tunepme
            else
                gmx grompp -f step6.{$cnt}_equilibration.mdp -o step6.{$cnt}_equilibration.tpr -c step6.{$pcnt}_equilibration.gro -r step5_input_aligned.gro -p topol.top -maxwarn -1
                mpirun -np $n_cores gmx_mpi mdrun -v -deffnm step6.{$cnt}_equilibration -tunepme -ntomp 1
            endif
            @ cnt += 1
        end
#endif


###############Production#############

set cnt    = 1
set cntmax = $2

while ( ${cnt} <= ${cntmax} )
        gmx grompp -f step7_production_macro.mdp -o $run_name"_"${cnt}.tpr -c step6.6_equilibration.gro -p topol.top -r target_pdb -maxwarn 2
        mpirun -np $n_cores gmx_mpi mdrun -v -deffnm $run_name"_"${cnt}
    @ cnt += 1
end
