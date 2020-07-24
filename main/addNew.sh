#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: ./addNew.sh <path-to-directory>"
    exit 1
fi
cp -r $1 ../source/
cd ../source/$(basename $1)
if [ ! -f ./Makefile ]; then
    echo "Makefile does not exist. Creating generic Makefile."
    cp ../../main/metamake.py .
    ./metamake.py 
    rm metamake.py
fi
cp ../../main/addTests.py .
./addTests.py
rm addTests.py