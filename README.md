A framework to help test Fortran compilers against `gfortran`.
This framework uses [BATS](https://github.com/bats-core/bats-core).

To add a directory to the testing framework, `cd` into this directory, and use `./addNew.sh <path-to-directory>`.

Then, use `./run_test.sh (compiler) (flags...)` to test all .bats files in the subdirectories of /source/. If no compiler is specified, it will default to using `gfortran` with no flags; otherwise, it will use the compiler and/or flags specified in the command line.

Other commands include:
- `./addEnv.sh <directory>`: adds `comm.env` a subdirectory, allowing BATS to use necessary environmental commands before testing, e.g. `module load`.
- `./addMakefile.sh <directory>`: creates a generic Makefile.
- `./cleanup.sh <directory>`: helps delete and/or backup your .bats files.

Originally written for Argonne National Laboratory.
