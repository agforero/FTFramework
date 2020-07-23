#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: ./addNew.sh <path-to-directory>"
    exit 1
fi
cp -r $1 source/
cd source/$(basename $1)
if [ ! -f ./Makefile ]; then
    echo "Makefile does not exist. Creating generic Makefile."
    cp ../../metamake.py .
    ./metamake.py 
    rm metamake.py
fi
cp ../../addTests.py .
./addTests.py
rm addTests.py