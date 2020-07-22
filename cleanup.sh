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
echo -e "-\tb to backup bats test(s) to $$.backup before deleting,"
echo -e "-\td to delete without backup, and"
echo -e "-\tc or anything else to cancel.\n"
echo "Your input:"
read input
if [ $input == "b" ]; then
	echo -e "\n\tDeleting; backing up to $$.backup."
	for f in *.bats; do
		cat "$f" >> "$$.backup"
		rm "$f"
	done
elif [ $input == "d" ]; then
	echo -e "\n\tDeleting without backup."
	for f in *.bats; do
		rm "$f"
	done
else 
	echo -e "\n\tCancelling."
	exit 0
fi
if [ -f env.txt ]; then
	echo -e "\nWould like to also remove env.txt?"
	echo -e "-\td to delete, and"
	echo -e "-\tc or anything else to cancel.\n"
	echo "Your input:"
	read input
	if [ $input == "d" ]; then
		echo -e "\n\tDeleting."
		rm env.txt
	else
		echo -e "\n\tExiting."
	fi
fi