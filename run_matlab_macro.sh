#!/usr/bin/env bash
#SBATCH -A C3SE408-22-1
#SBATCH -p vera
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -t 01:30:00

ml MATLAB/2019a

if [ ${1} == "YES" ]
then
	matlab -nodesktop -nojvm -nosplash -r "process_intensity_files(input_one,'input_two','input_three','nope');exit"
fi


module purge
echo "Intensity files have been converted to .mat files"
