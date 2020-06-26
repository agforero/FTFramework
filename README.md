A framework to help test compilers against their stable counterparts, e.g. a bleeding-edge Fortran compiler versus `gfortran`.
This framework uses [BATS](https://github.com/bats-core/bats-core) to test the functionality of its `Makefile`s.

To add a directory to the testing framework, cd into this directory, and use `./addNew.sh <path-to-directory>`.
Then, use `./run_test.sh` to test all .bats files in the subdirectories of /source/. 

Other commands include:
- `./addMakefile.sh <directory>`: open a lightweight Makefile generator.
- `./cleanup.sh <directory>`: a script to help delete and/or backup your .bats files.
