#!/bin/bash

# if they specify a compiler in the arg
if [ $# -eq 1 ]; then
	export FC="$1"
	comp=$1
else
	comp=$(echo $FC)
fi

# redirecting to logs/<compiler>.log
# mkdir -p logs
# set -o errexit
# readonly LOG_FILE="./logs/$comp.log"
# touch comp.log
# exec 1>$LOG_FILE
# exec 2>&1

foo() {
	cd source
	for d in */ ; do
		cd $d 
		echo
		echo "$(pwd):"
		# checking for need of metamake
		if [ ! -f ./Makefile ]; then
			cp ../../metamake.py .
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
		make clean
		bats tests.bats
		cp ../../findRelevant.py .
		restOfBats=$(./findRelevant.py -e .bats tests.bats) # attempts to find other .bats files
		if [ $? -eq 0 ]; then
			bats $restOfBats
		fi
		rm -f findRelevant.py
		cd ..
	done
	echo
}

foo_parallel() {
	cd source
	for d in */ ; do
		cd $d 
		echo
		echo "$(pwd):"
		# checking for need of metamake
		if [ ! -f ./Makefile ]; then
			cp ../../metamake.py .
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
		make clean
		jobs=$(nproc)
		bats -j $((jobs/2)) tests.bats
		cp ../../findRelevant.py .
		restOfBats=$(./findRelevant.py -e .bats tests.bats) # attempts to find other .bats files
		if [ $? -eq 0 ]; then
			bats -j $((jobs/2)) $restOfBats
		fi
		rm -f findRelevant.py
		cd ..
	done
	echo
}

which parallel # check to see if user has GNU Parallel installed
if [ $? -eq 0 ]; then
	foo_parallel |& tee logs/$comp.log
else
	foo |& tee logs/$comp.log
fi