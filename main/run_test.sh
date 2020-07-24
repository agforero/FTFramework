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
	for d in */ ; do
		cd $d 
		echo
		echo "$(pwd):"
		# checking for need of metamake
		if [ ! -f ./Makefile ]; then
			cp ../../main/metamake.py .
			./metamake.py
			rm metamake.py
		fi
		if [ ! -f "tests.bats" ]; then
			cp ../../main/addTests.py .
			./addTests.py
			rm addTests.py
		fi
		bats tests.bats
		cp ../../main/findRelevant.py .
		restOfBats=$(./findRelevant.py -e .bats tests.bats) # attempts to find other .bats files
		if [ $? -eq 0 ]; then
			bats $restOfBats
		fi
		rm -f findRelevant.py
		cd ..
	done
	echo
}

foo |& tee ../logs/$comp.log
