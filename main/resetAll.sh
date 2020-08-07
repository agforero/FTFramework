#!/bin/bash

cd ../source
for d in */ ; do
    cd $d 
    echo "deleting $(pwd)/tests.bats"
    rm -f tests.bats
    cd ../
done