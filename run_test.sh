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
	# running generic test
	# refreshes and runs test.bats
	if [ ! -f tests.bats ]; then
		cp ../../addTests.py .
		./addTests.py
		rm addTests.py
	fi
	bats tests.bats
	restOfBats=$(./findRelevant.py -e .bats tests.bats) # attempts to find other .bats files
	if [ $? -eq 0 ]; then
		bats $restOfBats
	fi
	rm findRelevant.py
	cd ..
done
echo