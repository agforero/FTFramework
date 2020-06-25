#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: ./addNew.sh <path-to-directory>"
    exit 1
fi

cp -r $1 source/