A framework to help test bleeding-edge Fortran compilers.
This framework uses [BATS](https://github.com/bats-core/bats-core).

To add a directory to the testing framework, `cd` into this directory, and use `./addNew.sh <path-to-directory>`.

Then, use `./run_test.sh (compiler)` to test all .bats files in the subdirectories of /source/. If no compiler is specified, it will use the current value of `$FC` (which is set to `f77` by default). You can use `export FC=<compiler>` to use your compiler of choice, e.g. `gfortran` a different compiler you'd like to test. You can also specify your compiler as the argument of `./run_test.sh`, which will `export` it for you.

Other commands include:
<<<<<<< HEAD
- `./addEnv.sh <directory>`: adds `env.txt` a subdirectory, allowing BATS to use necessary environmental commands before testing, e.g. `module load`.
- `./addMakefile.sh <directory>`: creates a generic Makefile. Runs automatically during `./run_test.sh` if no Makefile is present.
=======
- `./addEnv.sh <directory>`: adds `comm.env` a subdirectory, allowing BATS to use necessary environmental commands before testing, e.g. `module load`.
- `./addMakefile.sh <directory>`: creates a generic Makefile.
>>>>>>> ea02dffde10117d0f1b68c5f443dae6a68c02750
- `./cleanup.sh <directory>`: helps delete and/or backup your .bats files.

Originally written for Argonne National Laboratory.
