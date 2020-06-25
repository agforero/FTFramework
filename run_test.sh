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
	allbats=$(find *.bats)
	if [ $? -ne 0 ]; 
	then
		echo 
		echo "No .bats test found in $d. Writing generic make test..."
		echo "Enter preemptive environmental commands to be included in .bats test, e.g. \`module load gcc\`."
		echo "Enter \`ctrl + D\` to conclude your input."
		echo
		echo "---------------------------------------- BEGIN INPUT ----------------------------------------"
		echo "# Environmental commands:" >> test.bats	
		while read comm
		do
			echo "$comm" >> test.bats
		done

		echo "----------------------------------------- END INPUT -----------------------------------------"
		echo
		echo >> test.bats
		echo "# Generic Makefile test" >> test.bats
		echo "@test \"$d make\" {" >> test.bats
		echo "  make" >> test.bats
		echo "}" >> test.bats
		bats test.bats
    else 
		bats $allbats
	fi
    cd ..
done
