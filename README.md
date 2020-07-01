A framework to help test compilers against their stable counterparts, e.g. a bleeding-edge Fortran compiler versus `gfortran`.
This framework uses [BATS](https://github.com/bats-core/bats-core).

To add a directory to the testing framework, cd into this directory, and use `./addNew.sh <path-to-directory>`.
Then, use `./run_test.sh` to test all .bats files in the subdirectories of /source/. 

Other commands include:
- `./addMakefile.sh <directory>`: open a lightweight Makefile generator.
- `./addEnv.sh <directory>`: adds `env.txt` a subdirectory, allowing BATS to use necessary environmental commands before testing, e.g. `module load`.
- `./cleanup.sh <directory>`: helps delete and/or backup your .bats files.
- `./addTest.sh <directory>`: creates a `.bats` test for a directory's `Makefile`, and prompts for necessary environmental commands.

Originally written for Argonne National Laboratory.
