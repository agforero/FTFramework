# Fortran Testing Framework #
### A framework to help test Fortran compilers. Written by Agustin Forero for Argonne National Laboratory. ###

#### Summary: ####
To add a directory to the testing framework, `cd` into this directory, and use `./addNew.sh <path-to-directory>`.

Then, use `./run_test.sh (compiler)` to compile all files found in the subdirectories of `/source/`. If no compiler is specified, it will use the current value of `$FC` (which is set to `f77` by default within the Makefile). 

Each time `./run_test.sh` is executed, the results will be saved to `/logs/` as `<compiler>.log`, where you can view all compilation errors encountered during runtime. Additionally, `./compare.py <compiler1> <compiler2> (-g / -v)` outputs differences in errors between two given compilers, where:

- no flag outputs raw data without additional formatting,
- `-g` outputs only the first line in each error to a two-column table, and
- `-v` outputs errors in full to a two-column table.

#### Other commands include: ####
- `./addMakefile.sh <directory>`: creates a generic Makefile. Runs automatically during `./run_test.sh` if no Makefile is present.
- `./addEnv.sh <directory>`: adds `comm.env` a subdirectory, allowing BATS to use necessary environmental commands before testing, e.g. `module load`.
- `./cleanup.sh <directory>`: helps delete and/or backup your .bats files. This is important for when you add or delete files from a source directory.

#### Example output: ####

```
$ ./compare.py compiler1 compiler2 -g

                                               compiler1 | compiler2
-------------------------------------------------------- | --------------------------------------------------------
                                                         |
---------------------------------------- /some/directory/on-your-machine: ----------------------------------------
not ok some error 1 not encountered in other file        |
not ok some error 2 not encountered in other file        |
not ok some error 3 not encountered in other file        |
                                                         | not ok some error 4 not encountered in other file       
                                                         | not ok some error 5 not encountered in other file       
                                                         | not ok some error 6 not encountered in other file    
```

This framework uses [BATS](https://github.com/bats-core/bats-core).
