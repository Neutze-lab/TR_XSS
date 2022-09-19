#!/bin/bash

function monitor {
        all_lines=1
        lines=$(wc -l ${1} | cut -d" " -f1)
        while [[ $all_lines != $lines ]];
       	do
                all_lines=0
                for y in $(tail -$lines ${1} | cut -c 20-30)
               	do
			squeue -j $y &> ${y}_job.tmp
                        current_lines=$(wc -l ${y}_job.tmp | cut -d" " -f1)
                        all_lines=$((all_lines+current_lines))
			echo $(squeue -j ${y})
			echo ${current_lines} " - current lines"
			echo ${all_lines} " - all lines"
			echo ${lines} " - lines"
               	done
		echo " "
               	sleep 5
       	done
rm ${1}
rm *_job.tmp
}

function check_and_submit {
	pdb_list=$(ls ${1}*pdb)
	int_list=$(ls ${1}*int)
	found=0
	run_no=1
	if [ ! $2 == $(ls ${1}*int | wc -l) ]
	then
		echo Missing crysol files detected
		for x in $pdb_list
		do
			for y in $int_list
			do
				if [ ${x%????} == ${y%????} ]
				then
					found=$((found+1))
				fi
			done

			if [ $found == 0 ]
			then
				echo Submitting $x to crysol again
				sbatch run_crysol_calc.sh $x $run_no >> crysol_catchup_job.tmp	
				run_no=$((run_no+1))
			fi
			found=0
		done
	else
		echo "Crysol didn't miss any"
	fi
}


## Here beginth the real script

total_files=$(ls ${1}*pdb | wc -l)

files_per_cpu=10
to_use_cpus=$((total_files/files_per_cpu))

files_list=$(ls ${1}*pdb) 
files_list=($files_list)
date=$(date +%Y-%m-%d)

echo $total_files > crysol_log_file.log
echo $files_per_cpu >> crysol_log_file.log
echo $to_use_cpus >> crysol_log_file.log
echo >> crysol_log_file.log

echo Submiting Crysol Jobs to queue.
for run in `seq 1 $((to_use_cpus+1))`
do
	if [ $run -eq 1 ]
	then
		first_end_point=$((files_per_cpu))
		sub_list=("${files_list[@]:0:$files_per_cpu}")
		sbatch run_crysol_calc.sh "${sub_list[*]}" $run > crysol_job.tmp
		printf '%s\n' "${sub_list[@]}" >> crysol_log_file.log
		echo >> crysol_log_file.log
	else
		current_list_start=$(($((run-1))*files_per_cpu))
		sub_list=("${files_list[@]:$current_list_start:$files_per_cpu}")
		sbatch run_crysol_calc.sh "${sub_list[*]}" $run >> crysol_job.tmp
		printf '%s\n' "${sub_list[@]}" >> crysol_log_file.log
		echo >> crysol_log_file.log
		
	fi
done

monitor crysol_job.tmp
echo Crysol finished.
rm temp*

echo Checking all files went through crysol
for i in `seq 1 5`
do
	check_and_submit $1 $total_files
	if [ -e crysol_catchup_job.tmp ]
	then
		if [ i == 5 ]
		then
			echo Files still not processing after 4 runs, crysol run time problem
			mv crysol_catchup_job.tmp crysol_failed_jobs.lst
			exit
		fi
		monitor crysol_catchup_job.tmp
	fi
done

if [ $total_files == $(ls ${1}*int | wc -l) ]
then
	echo All files processed through crysol
fi
