#!/bin/bash

cd runs
for d in */ ; do
    rm -r $d
done