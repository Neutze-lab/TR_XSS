#!/bin/bash

sed -i "s/the_name/${1}/g" segfault_submission.sh
last_successful_run=$(ls ${1}*prev.cpt | wc -l)
sed -i "s/the_number/${2}/g" segfault_submission.sh

sbatch segfault_submission.sh > tmp.file

run_number=$(head -1 tmp.file | cut -c 20-30)
echo $run_number
squeue -j $run_number > tmp.file
while [[ $(wc -l tmp.file | cut -c 1) != [1] ]];
do
        echo Gromacs still running
        squeue -j $run_number > tmp.file
        sleep 20
done

