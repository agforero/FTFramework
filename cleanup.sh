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
cd $1
if [ ! -f *.bats ]; then
	echo "No .bats file(s) found."
	exit 3
fi
cd ..

# at this point, we're clear for takeoff
# ...but we should really ask for confirmation
echo 
echo -e "\tEnter:"
echo -e "-\t<b> to backup bats test(s) to backup to $$.backup before deleting,"
echo -e "-\t<d> to delete without backup, and"
echo -e "-\t<c> or anything else to cancel."

echo
echo "Your input:"
read input

if [ $input == "b" ]; then
	echo "Deleting; backing up to $$.backup."
	cd $1 
	for f in *.bats; do
		cat "$f" >> "$$.backup"
		rm "$f"
	done

elif [ $input == "d" ]; then
	echo "Deleting without backup."
	cd $1
	for f in *.bats; do
		rm "$f"
	done

else 
	echo "Exiting."
	exit 0
fi
