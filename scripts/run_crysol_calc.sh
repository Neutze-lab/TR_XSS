#!/usr/bin/env bash
#SBATCH -A C3SE408-22-1
#SBATCH -p vera
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -t 01:00:00


#if [ -z $1 ]
#then
#	echo "no file extension specified"
#fi

#if [ -z $2 ]
#then
#	echo "Using only one core"
#fi

for f in $1

do
	echo "Processing $f file.";
	head -n -2 $f > temp$2.pdb
	echo "Calculating intensites for $f via crysol";
	singularity exec /path/ubuntu.simg bash -c "LD_LIBRARY_PATH=~/macro_v2/ATSAS-3.0.3-1/lib/atsas/ && ./crysol temp$2.pdb -lm 50 -fb 18 -sm 2 -ns 201 -dro 0 -eh >> temp$2.inf"

	if [ ! -e temp"$2"00.int ]
	then
		echo "couldn't find the $2 '.int' file"
		exit
	fi

	mv temp"$2"00.int ${f%????}.int
	mv temp"$2"00.log ${f%????}.log
	mv temp$2.inf ${f%????}.inf
	rm temp$2.pdb
	rm temp$2.alm
done
