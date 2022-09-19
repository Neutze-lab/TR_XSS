#!/bin/bash

## SETTING INITIAL VARIABLES

run_name=step7
production_run_number=1

## Choose to recieve protein specify YES or NO
target=YES
crysol=YES

 ml "GCC/8.3.0"
 ml OpenMPI/3.1.4
 ml GROMACS/2020

mkdir ${run_name}


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


module purge

echo Everything is finished please check it has run successfully.
