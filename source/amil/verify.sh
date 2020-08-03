#!/bin/bash

for f in *.f90 ; do
	echo $f:
	ls "${f%%.*}".o
done
