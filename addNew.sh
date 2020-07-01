#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: ./addNew.sh <path-to-directory>"
    exit 1
fi
cp -r $1 source/
cd source/$1
cp ../../findRelevant.py .
if [ ! -f ./Makefile ]; then
    cp ../../metamake.py .
    ./metamake.py
    rm ../../metamake.py 
fi

