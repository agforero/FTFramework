#!/bin/bash

if [ "$#" -ne 1 ]; then
	echo "Usage: ./addMakefile.sh directory"
	exit 1
fi
cd ../source
if [ ! -d $1 ]; then
	echo "/$1/ does not exist in /source/"
	exit 2
fi
cd $1
if [ -f "Makefile" ]; then
	echo "Makefile already exists; enter 'y' to overwrite."
	read conf
	if [ $conf != "y" ]; then
		echo "Exiting."
		exit 0
	fi
fi
cp ../../main/metamake.py . 
./metamake.py 
rm metamake.py