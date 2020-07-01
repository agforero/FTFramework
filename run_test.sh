#!/bin/bash

cd source
for d in */ ; do
	cd $d 
	echo
	echo "$(pwd):"
	cp ../../findRelevant.py .
	# checking for need of metamake
	if [ ! -f ./Makefile ]; then
		cp ../../metamake.py .	
		echo "Makefile does not exist. Creating generic Makefile."
		./metamake.py
		rm metamake.py
	fi
	# testing to see if generic test is needed
	../../addTest.sh $d $@ # refreshes and runs test.bats
	bats test.bats
	restOfBats=$(./findRelevant.py -e .bats test.bats) # attempts to find other .bats files
	if [ $? -eq 0 ]; then
		bats $restOfBats
	fi
	cd ..
done
echo