#!/bin/bash

echo "Job started at" && date +"%T"

################## GROUNDSTATE #####################

mv ../step7-ground.tar.gz .


files=501

tar xzf step7-ground.tar.gz --wildcards '*.int' '*.pdb'

echo "These are the pdb files"
membrane_pdb="`ls *membrane*.pdb | wc -l`"
echo  "$membrane_pdb" 
protein_pdb="`ls *protein*.pdb | wc -l`"
echo "$protein_pdb"
protmem_pdb="`ls *protmem*.pdb | wc -l`"
echo  "$protmem_pdb"

echo  "These are the int files"
membrane_int="`ls *membrane*.int | wc -l`"
echo  "$membrane_int"
protein_int="`ls *protein*.int | wc -l`"
echo  "$protein_int"
protmem_int="`ls *protmem*.int | wc -l`"
echo  "$protmem_int"

mkdir ground

# CHECKING INT FILES 
if [ "$membrane_int" == "$files" ]
then
	echo "correct number of membrane int files"
	rm *membrane*.int
elif [ "$membrane_int" -lt "$files" ]
then
	echo "incorrect number of membrane int files. Check run."
fi

if [ "$protein_int" == "$files" ]
then
	echo "correct number of protein int files"
	rm *protein*.int
elif [ "$protein_int" -lt "$files" ]
then
	echo "incorrect number of protein int files. Check run."
fi

if [ "$protmem_int" == "$files" ]
then
	echo "correct number of protmem int files"
	rm *protmem*.int
elif [ "$protmem_int" -lt "$files" ]
then
	 echo "incorrect number of protmem int files. Check run."
fi

# CHECKING PDB FILES

if [ "$membrane_pdb" == "$files" ]
then
	echo "correct number of membrane pdb files"
	mv *membrane*.pdb ground
elif [ "$membrane_pdb" -lt "$files" ]
then
	echo "incorrect number of membrane pdb files. Check run."
fi

if [ "$protein_pdb" == "$files" ]
then
	echo "correct number of protein pdb files"
	mv *protein*.pdb ground
elif [ "$protein_pdb" -lt "$files" ]
then
	echo "incorrect number of protein pdb files. Check run."
fi

if [ "$protmem_pdb" == "$files" ]
then
	echo "correct number of protmem pdb files"
	mv *protmem*.pdb ground
elif [ "$protmem_pdb" -lt "$files" ]
then
	echo "incorrect number of protmem pdb files. Check run."
fi

################## EXCITED STATE #####################

files=501

tar xzf step7.tar.gz --wildcards '*.int' '*.pdb'

echo "These are the pdb files"
membrane_pdb="`ls *membrane*.pdb | wc -l`"
echo  "$membrane_pdb" 
protein_pdb="`ls *protein*.pdb | wc -l`"
echo "$protein_pdb"
protmem_pdb="`ls *protmem*.pdb | wc -l`"
echo  "$protmem_pdb"

echo  "These are the int files"
membrane_int="`ls *membrane*.int | wc -l`"
echo  "$membrane_int"
protein_int="`ls *protein*.int | wc -l`"
echo  "$protein_int"
protmem_int="`ls *protmem*.int | wc -l`"
echo  "$protmem_int"

mkdir excited

# CHECKING INT FILES 
if [ "$membrane_int" == "$files" ]
then
	echo "correct number of membrane int files"
	rm *membrane*.int
elif [ "$membrane_int" -lt "$files" ]
then
	echo "incorrect number of membrane int files. Check run."
fi

if [ "$protein_int" == "$files" ]
then
	echo "correct number of protein int files"
	rm *protein*.int
elif [ "$protein_int" -lt "$files" ]
then
	echo "incorrect number of protein int files. Check run."
fi

if [ "$protmem_int" == "$files" ]
then
	echo "correct number of protmem int files"
	rm *protmem*.int
elif [ "$protmem_int" -lt "$files" ]
then
	 echo "incorrect number of protmem int files. Check run."
fi

# CHECKING PDB FILES

if [ "$membrane_pdb" == "$files" ]
then
	echo "correct number of membrane pdb files"
	mv *membrane*.pdb excited
elif [ "$membrane_pdb" -lt "$files" ]
then
	echo "incorrect number of membrane pdb files. Check run."
fi

if [ "$protein_pdb" == "$files" ]
then
	echo "correct number of protein pdb files"
	mv *protein*.pdb excited
elif [ "$protein_pdb" -lt "$files" ]
then
	echo "incorrect number of protein pdb files. Check run."
fi

if [ "$protmem_pdb" == "$files" ]
then
	echo "correct number of protmem pdb files"
	mv *protmem*.pdb excited
elif [ "$protmem_pdb" -lt "$files" ]
then
	echo "incorrect number of protmem pdb files. Check run."
fi

################## EXC-MEM-REST-PROT #####################

for f in excited/step7_1_membrane*.pdb
do
echo "Processing $f file.";
tail -n +6 $f > temp.pdb
mv temp.pdb ${f%.*}_MOD.pdb
done

for g in ground/step7_1_protein*.pdb
do
echo "Processing $g file.";
head -n -2 $g > temp.pdb
mv temp.pdb ${g%.*}_MOD.pdb
done

mv excited/step7_1_membrane*MOD.pdb .
mv ground/step7_1_protein*MOD.pdb .

membrane_files=$(ls step7_1_membrane*_MOD.pdb)
protein_files=$(ls step7_1_protein*_MOD.pdb)
nr_files=$(wc -w <<< $membrane_files)

for x in `seq -f "%04g" 1 1 $nr_files`
do
echo "Concatenating $x file.";
membrane_file=$(echo $membrane_files | cut -d" " -f${x})
protein_file=$(echo $protein_files | cut -d" " -f${x})
cat $protein_file $membrane_file > protmem${x}.pdb
rename protmem step7_1_protmem protmem*.pdb
done

#echo "Creating folder for modified files and moving them there"
mkdir MOD_files
mv *membrane*MOD.pdb MOD_files
mv *protein*MOD.pdb MOD_files

rm -r MOD_files

#echo "Copying original excited membrane files and ground protein files to exc-mem-rest-prot"
mkdir exc-mem-rest-prot
cp -r -v excited/*step7_1_membrane*.pdb exc-mem-rest-prot
cp -r -v ground/*step7_1_protein*.pdb exc-mem-rest-prot 
mv *protmem*.pdb exc-mem-rest-prot

#tar cvzf excited_ground.tar.gz excited ground
#rm -r excited
#rm -r ground



echo "Finished exc-mem-rest-prot"

################## REST-MEM-EXC-PROT #####################

for h in excited/step7_1_protein*.pdb
do
echo "Processing $h file.";
head -n -2 $h > temp.pdb
mv temp.pdb ${h%.*}_MOD.pdb
done

for q in ground/step7_1_membrane*.pdb
do
echo "Processing $q file.";
tail -n +6 $q > temp.pdb
mv temp.pdb ${q%.*}_MOD.pdb
done

mv excited/*step7_1_protein*MOD.pdb .
mv ground/*step7_1_membrane*MOD.pdb .

membrane_files=$(ls step7_1_membrane*_MOD.pdb)
protein_files=$(ls step7_1_protein*_MOD.pdb)
nr_files=$(wc -w <<< $membrane_files)

for x in `seq -f "%04g" 1 1 $nr_files`
do
echo "Concatenating $x file.";
membrane_file=$(echo $membrane_files | cut -d" " -f${x})
protein_file=$(echo $protein_files | cut -d" " -f${x})
cat $protein_file $membrane_file > protmem${x}.pdb
rename protmem step7_1_protmem protmem*.pdb
done

echo "Creating folder for modified files and moving them there"
mkdir MOD_files
mv *membrane*MOD.pdb MOD_files
mv *protein*MOD.pdb MOD_files

rm -r MOD_files

echo "Copying original excite membrane files and ground protein files to rest-mem-exc-prot"
mkdir rest-mem-exc-prot
cp -r -v excited/*step7_1_protein*.pdb rest-mem-exc-prot
cp -r -v ground/*step7_1_membrane*.pdb rest-mem-exc-prot
mv *protmem*.pdb rest-mem-exc-prot

tar cvzf excited_ground.tar.gz excited ground
rm -r excited
rm -r ground

echo "Finished at" && date +"%T"
