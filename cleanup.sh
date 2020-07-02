#!/bin/bash

# check to see they entered an argument
if [ $# -ne 1 ]; then
	echo "Usage: ./cleanup.sh <directory>"
	exit 1
fi
# check to see if the directory exists in /source/
cd source
if [ ! -d $1 ]; then
	echo "/$1/ does not exist in /source/"
	exit 2
fi
# check to see if there are even .bats files to interact with
# this logically shouldn't happen if they've run ./run_test at least once
cd $1
if [ ! -f *.bats ]; then
	echo "No .bats file(s) found."
	exit 3
fi
# at this point, we're clear for takeoff
# ...but we should really ask for confirmation
echo 
echo -e "\tEnter:"
echo -e "- b to backup bats test(s) to backup to $$.backup before deleting,"
echo -e "- d to delete without backup, and"
echo -e "- c or anything else to cancel."
echo
echo "Your input:"
read input
if [ $input == "b" ]; then
	echo "Deleting; backing up to $$.backup."
	for f in *.bats; do
		cat "$f" >> "$$.backup"
		rm "$f"
	done
elif [ $input == "d" ]; then
	echo "Deleting without backup."
	for f in *.bats; do
		rm "$f"
	done
else 
	echo "Done."
	exit 0
fi
if [ -f env.txt ]; then
	echo -e "\nWould like to also remove env.txt?"
	echo -e "- d to delete, and"
	echo -e "- c or anything else to cancel."
	echo "Your input:"
	read input
	if [ $input == "d" ]; then
		rm env.txt
	else
		echo "Exiting."
	fi
fi