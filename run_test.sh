#!/bin/bash

# if they specify a compiler in the arg
if [ $# -eq 1 ]; then
	export FC="$1"
fi

cd source
for d in */ ; do
	cd $d 
	echo
	echo "$(pwd):"
	# checking for need of metamake
	if [ ! -f ./Makefile ]; then
		cp ../../metamake.py .	
		echo "Makefile does not exist. Creating generic Makefile."
		./metamake.py
		rm metamake.py
	fi
	# running generic test
	# refreshes and runs test.bats
	if [ ! -f "tests.bats" ]; then
		cp ../../addTests.py .
		./addTests.py
		rm addTests.py
	fi
	bats tests.bats
	cp ../../findRelevant.py .
	restOfBats=$(./findRelevant.py -e .bats tests.bats) # attempts to find other .bats files
	if [ $? -eq 0 ]; then
		rm -f findRelevant.py
		bats $restOfBats
	else
		rm -f findRelevant.py
	fi
	cd ..
done
echo
