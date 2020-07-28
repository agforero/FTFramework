#!/bin/bash

# if they specify a compiler in the arg
if [ $# -eq 1 ]; then
	export FC="$1"
	comp=$1
else
	comp=$(echo $FC)
fi

foo() {
	cd ../source
	echo "COMPILER USED: $FC"
	for d in */ ; do
		cd $d 
		echo
		echo "DIR: $(pwd)"
		# checking for need of metamake
		if [ ! -f ./Makefile ]; then
			../../main/metamake.py
		fi
		if [ ! -f ./"tests.bats" ]; then
			../../main/addTests.py
		fi
		bats tests.bats
		restOfBats=$(../../main/findRelevant.py -e .bats tests.bats) # attempts to find other .bats files
		if [ $? -eq 0 ]; then
			bats $restOfBats
		fi
		cd ..
	done
	echo
}

foo |& tee ../logs/$comp.log
