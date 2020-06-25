#!/bin/bash

if [ "$#" -ne 1 ]; then
	echo "Usage: ./addMakefile.sh directory"
	exit 1
fi

cd source

if [ ! -d $1 ]; then
	echo "Directory does not exist in /source/"
	exit 2
fi

cd $1

if [ -f "Makefile" ]; then
	echo "Makefile already exists; enter 'y' to overwrite."
	read conf

	if [ $conf != "y" ]; then
		echo "Exiting."
		exit 0
	fi
fi

cp ../../metamake.py .  

echo "Enter your desired..."
echo "Name of program (_____.exe): "
read fname
echo "FC (Fortran compiler): "
read fcomp
echo "FFLAGS (compiler flags, e.g. -fopenmp): "
read fflags
echo "Extension (e.g.: .F90): "
read ver

./metamake.py $fname $fcomp $fflags $ver
rm metamake.py
