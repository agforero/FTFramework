#!/bin/bash

cd ../source
for d in */ ; do
    cd $d 
    echo "Resetting $(pwd)"
    rm -f tests.bats
    cd ../
done