#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: ./addNew.sh <path-to-directory>"
    exit 1
fi
cp -r $1 ../source/
cd ../source/$(basename $1)
if [ ! -f ./Makefile ]; then
    echo "Makefile does not exist. Creating generic Makefile."
    ../../main/metamake.py 
fi
../../main/addTests.py