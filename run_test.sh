#!/bin/bash

cd source
for d in */ ; do
    cd $d 
    echo
    echo "$(pwd):"

	# checking for need of metamake
    if [ ! -f ./Makefile ]; then
		cp ../../metamake.py .	

		echo "Makefile does not exist. Enter your desired..."
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
    fi

	# testing to see if generic test is needed
	allbats=$(find *.bats > /dev/null 2>&1)
	if [ $? -ne 0 ]; 
	then
		echo "No .bats test found in $d. You can run ./addTest.sh <directory> to help."
    else 
		bats $allbats
	fi
    cd ..
done
echo
