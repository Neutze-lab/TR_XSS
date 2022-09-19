#!/bin/bash

## SETTING INITIAL VARIABLES

run_name=step7
run_length=10000000
steps_per_readout=0
xtc_steps_per_readout=5000
production_run_number=1
no_stacks=2
no_cores=$((no_stacks*32))
hours=20
target_pdb=target.gro
id=

## Choose to recieve protein specify YES or NO
target=YES
crysol=YES
full_GROMACS=NO
#fake_skip=YES


## Print settings for reference

echo "Name:			${run_name}"
echo "Run Length:		${run_length}"
echo "Steps per PDB readout:	${steps_per_readout}"
echo "Number of runs:		${production_run_number}"
echo "Number of Nodes:	${no_stacks}"
echo "Hours: ${hours}"

###################################################################
###################################################################
####################     MAIN SCRIPT    ###########################
###################################################################
###################################################################

ml "GCC/8.3.0"
ml OpenMPI/3.1.4
ml GROMACS/2020

echo Copying over input files to pwd
rsync -av --exclude='topol.top' /path/bR-190BOG-inputs/* .

echo Creating gro file from MatWAXS pdb using charmm36 forcefield 
gmx pdb2gmx -f $id.pdb -o $id.gro < input_command_pdb2gmx.tmp

rm topol.top
rm posre.itp

cp /path/bR-190BOG-inputs/topol.top .
export id

./es2gs.dms

echo Starting Gromacs

if [ "$full_GROMACS" == "YES" ]
then
	cp /path/step7_production_macro.mdp .
	cp /path/danMac2.sh .
	sed -i "s/target_pdb/${target_pdb}/g" danMac2.sh
	sed -i "s/SBATCH -N 6/SBATCH -N ${no_stacks}/g" danMac2.sh
	sed -i "s/SBATCH -t 25:00:00/SBATCH -t ${hours}:00:00/g" danMac2.sh
	sed -i "s/set n_cores=120/set n_cores=${no_cores}/g" danMac2.sh
	sed -i "s/nsteps                  = /nsteps                  = ${run_length}/g" step7_production_macro.mdp
	sed -i "s/nstxout                 = /nstxout                 = ${steps_per_readout}/g" step7_production_macro.mdp
	sed -i "s/nstvout                 = /nstvout                 = ${steps_per_readout}/g" step7_production_macro.mdp
	sed -i "s/nstfout                 = /nstfout                 = ${steps_per_readout}/g" step7_production_macro.mdp
	sed -i "s/nstxtcout               = /nstxtcout               = ${xtc_steps_per_readout}/g" step7_production_macro.mdp
	mv danMac2.sh dan-${id}.sh
	echo "Full GROMACS Target MD simulation with micelle is initiated"
	sbatch dan-${id}.sh ${run_name} ${production_run_number} > tmp.file
fi

if [ "$full_GROMACS" == "NO" ]
then
	cp /path/step7_production_macro.mdp .
	cp /path/danMac2_production.sh .
	sed -i "s/target_pdb/${target_pdb}/g" danMac2_production.sh
	sed -i "s/SBATCH -N 6/SBATCH -N ${no_stacks}/g" danMac2_production.sh
	sed -i "s/SBATCH -t 25:00:00/SBATCH -t ${hours}:00:00/g" danMac2_production.sh
	sed -i "s/set n_cores=120/set n_cores=${no_cores}/g" danMac2_production.sh
	sed -i "s/nsteps                  = /nsteps                  = ${run_length}/g" step7_production_macro.mdp
	sed -i "s/nstxout                 = /nstxout                 = ${steps_per_readout}/g" step7_production_macro.mdp
	sed -i "s/nstvout                 = /nstvout                 = ${steps_per_readout}/g" step7_production_macro.mdp
	sed -i "s/nstfout                 = /nstfout                 = ${steps_per_readout}/g" step7_production_macro.mdp
	sed -i "s/nstxtcout               = /nstxtcout               = ${xtc_steps_per_readout}/g" step7_production_macro.mdp
	mv danMac2_production.sh dan-${id}.sh
	echo "Production GROMACS Target MD simulation with micelle is initiated"
	sbatch dan-${id}.sh ${run_name} ${production_run_number} > tmp.file
fi

if [ "$full_GROMACS" == "RESTART" ]
then
	cp /path/gromacs_restart.sh .
	cp bR-190BOG-exc-mem-exc-prot/step7_1.* .
	sed -i "s/set n_cores=120/set n_cores=${no_cores}/g" gromacs_restart.sh
	echo "Production GROMACS Target MD simulation with micelle is continued."
	sbatch gromacs_restart.sh ${run_name} ${production_run_number} > tmp.file
fi

run_number=$(head -1 tmp.file | cut -c 20-30)
echo $run_number
squeue -j $run_number > tmp.file
while [[ $(wc -l tmp.file | cut -c 1) != [1] ]];
do
  	echo Gromacs still running
        squeue -j $run_number > tmp.file
        sleep 20
done

a=0
while [ $(grep "Segmentation fault" slurm*out | wc -l) != 0 ];
do
        zip ${run_name}_slurms_segf_${a}.zip slurm*out
	cp /path/segfault_handling.sh .
	cp /path/segfault_submission.sh .

        sed -i "s/SBATCH -N 20/SBATCH -N ${no_stacks}/g" segfault_submission.sh
        sed -i "s/SBATCH -t 20:00:00/SBATCH -t ${hours}:00:00/g" segfault_submission.sh
        sed -i "s/set n_cores=400/set n_cores=${no_cores}/g" segfault_submission.sh

	./segfault_handling.sh ${run_name} ${production_run_number}
	rm segfault_handling.sh
	rm segfault_submission.sh
	a=$((a+1))
done

rm tmp.file

mkdir ${run_name}

#### CONVERTING TRAJ FILES AND GENERATING PDBS ####
###################################################
gmx make_ndx -f $run_name"_1".gro -o ${run_name}"_index.ndx" < /path/input_commands.tmp

for x in `seq 1 ${production_run_number}`
do
	if [ $target == "YES" ]
	then
		gmx_mpi trjconv -s $run_name"_"$x.gro -f ${run_name}"_"$x.xtc -dt 0.002 -pbc whole -nzero 4 -n ${run_name}_index.ndx -sep -o ${run_name}_${x}_protmem.pdb < /path/traj_input_commands_PROT_LYR_BOG.tmp
		gmx_mpi trjconv -s $run_name"_"$x.gro -f ${run_name}"_"$x.xtc -dt 0.002 -pbc whole -nzero 4 -n ${run_name}_index.ndx -sep -o ${run_name}_${x}_protein.pdb < /path/traj_input_commands_PROT_LYR.tmp
		gmx_mpi trjconv -s $run_name"_"$x.gro -f ${run_name}"_"$x.xtc -dt 0.002 -pbc whole -nzero 4 -n ${run_name}_index.ndx -sep -o ${run_name}_${x}_membrane.pdb < /path/traj_input_commands_BOG.tmp
	fi
        mv $run_name"_"${x}.trr $run_name/
        mv $run_name"_"${x}.tpr $run_name/
        mv $run_name"_"${x}.cpt $run_name/
        mv $run_name"_"${x}.edr $run_name/
        mv $run_name"_"${x}.log $run_name/
        mv $run_name"_"${x}.gro $run_name/
	mv $run_name"_"${x}.xtc $run_name/
done

#### CLEANING UP FILES ####
###########################
mv step7_production_macro.mdp ${run_name}/
mv ${run_name}/step7_production_macro.mdp ${run_name}/${run_name}"_production.mdp"
mv step7_production.mdp step7_production_original.mdp
mv ${run_name}"_index.ndx" ${run_name}/
mv step7_production_original.mdp ${run_name}/

#### RUNNING CRYSOL ####
########################
if [ $crysol == "YES" ]
then
	cp /path/crysol .
	cp /path/run_crysol_calc.sh .
	cp /path/submit_parallelized_crysol_macro.sh .

	if [ $target == "YES" ]
	then
		./submit_parallelized_crysol_macro.sh "${run_name}*protmem"
		./submit_parallelized_crysol_macro.sh "${run_name}*protein"
		./submit_parallelized_crysol_macro.sh "${run_name}*membrane"
	fi
fi
rm crysol
rm run_crysol_calc.sh
rm submit_parallelized_crysol_macro.sh

#### GENERATE THE .MAT FILE ####
################################
cp /path/process_intensity_files.m .
cp /path/run_matlab_macro.sh .
sed -i "s/input_one/${production_run_number}/g" run_matlab_macro.sh
sed -i "s/input_two/${run_name}/g" run_matlab_macro.sh
sed -i "s/input_three/protein/g" run_matlab_macro.sh

	sbatch run_matlab_macro.sh ${crysol} > matlab.tmp
	run_number=$(head -1 matlab.tmp | cut -c 20-30)
	echo $run_number
	squeue -j $run_number > matlab.tmp
	while [[ $(wc -l matlab.tmp | cut -c 1) != [1] ]];
	do
  		echo Still waiting for Matlab
        	squeue -j $run_number > matlab.tmp
        	sleep 5
	done

if [ ${target} == "YES" ]
then
	sed -i "s/protein/protmem/g" run_matlab_macro.sh
	sbatch run_matlab_macro.sh ${crysol} > matlab.tmp
	run_number=$(head -1 matlab.tmp | cut -c 20-30)
	echo $run_number
	squeue -j $run_number > matlab.tmp
	while [[ $(wc -l matlab.tmp | cut -c 1) != [1] ]];
	do
					echo Still waiting for Matlab
					squeue -j $run_number > matlab.tmp
					sleep 5
	done
	sed -i "s/protmem/membrane/g" run_matlab_macro.sh
        sbatch run_matlab_macro.sh ${crysol}  > matlab.tmp
        run_number=$(head -1 matlab.tmp | cut -c 20-30)
        echo $run_number
        squeue -j $run_number > matlab.tmp
        while [[ $(wc -l matlab.tmp | cut -c 1) != [1] ]];
        do
                echo Still waiting for Matlab
                squeue -j $run_number > matlab.tmp
                sleep 5
        done
fi



rm matlab.tmp
rm run_matlab_macro.sh
rm process_intensity_files.m
rm process_debye_intensity_files.m

echo "Finished Calculations with matlab"

#### FINAL FILE CLEAN UP ####
#############################

find -maxdepth 1 -type f -name "${run_name}*.pdb" > pdb_list.tmp
find -maxdepth 1 -type f -name "${run_name}*.log" > log_list.tmp
find -maxdepth 1 -type f -name "${run_name}*.int" > int_list.tmp
find -maxdepth 1 -type f -name "${run_name}*.inf" > inf_list.tmp
find -maxdepth 1 -type f -name "${run_name}*.mat" > mat_list.tmp

cat pdb_list.tmp log_list.tmp int_list.tmp inf_list.tmp mat_list.tmp > all_files.tmp
tar -cvz -T all_files.tmp -f ${run_name}.tar.gz

find -maxdepth 1 -type f -name "${run_name}*.pdb" -delete
find -maxdepth 1 -type f -name "${run_name}*.log" -delete
find -maxdepth 1 -type f -name "${run_name}*.int" -delete
find -maxdepth 1 -type f -name "${run_name}*.inf" -delete
find -maxdepth 1 -type f -name "${run_name}*_D.mat" -delete

rm pdb_list.tmp log_list.tmp int_list.tmp inf_list.tmp mat_list.tmp all_files.tmp

zip ${run_name}_slurms.zip slurm*out
mv ${run_name}_slurms.zip ${run_name}/
mv ${run_name}.tar.gz ${run_name}/
mv ${run_name}*.mat ${run_name}/


cp check_files-create-pseudo-protmem.sh step7/ && cd step7 && ./check_files-create-pseudo-protmem.sh 

cd exc-mem-rest-prot && cp /path/macro_pseudo_processing.sh . \
&& ./macro_pseudo_processing.sh && cd ../rest-mem-exc-prot && cp /path/macro_pseudo_processing.sh . \
&& ./macro_pseudo_processing.sh && cd ../../

echo "Job started at" && date +%T

##################### rename & move ########################

name0=bR-190BOG-${id}
name1=bR-190BOG-${id}-exc-mem-exc-prot
name2=bR-190BOG-${id}-exc-mem-rest-prot
name3=bR-190BOG-${id}-rest-mem-exc-prot

mkdir ../${name0}

mkdir input_fit
echo "Created folder input_fitting"

mv step7 ${name1} && cd ${name1} && rename step7 ${name1} *.mat
echo "Renamed dir and files of exc-mem-exc-prot"

cd exc-mem-rest-prot && mv step7 ${name2} && cd ${name2} && rename step7 ${name2} *.mat \
&& cd ../ && mkdir -p ../../input_fit/${name2} && cp ${name2}/*.mat ../../input_fit/${name2}/
echo "Renamed and moved dir and files of exc-mem-rest-prot"

cd ../rest-mem-exc-prot/ && mv step7 ${name3} && cd ${name3} && rename step7 ${name3} *.mat \
&& cd ../ && mkdir -p ../../input_fit/${name3} && cp ${name3}/*.mat ../../input_fit/${name3}/
echo "Renamed and moved dir and files of rest-mem-exc-prot"

cd ../../ && mkdir -p input_fit/${name1} && cp ${name1}/*.mat input_fit/${name1}/
echo "Moved exc-mem-exc-prot"

mv input_fit ../${name0}/



echo "Job finished at" && date +%T

module purge

echo Everything is finished please check it has run successfully.
