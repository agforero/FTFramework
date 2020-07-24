# Fortran Testing Framework #
### A framework to help test Fortran compilers. Written by Agustin Forero for Argonne National Laboratory. ###

#### A typical test might look like: ####
```
$ ./run_test.sh bleedingedgecompiler
$ ./compare.py gfortran bleedingedgecompiler
```

#### Summary: ####
This is a framework built to help test bleeding-edge FORTRAN compilers. By testing the 
compilation of a wide variety of FORTRAN programs, and cross-checking results with stable 
compilers like `gfortran`, one can find where a compiler might be going wrong.

To add a directory containing FORTRAN programs to use for testing compilation, `cd` into 
this directory, and use `./addNew.sh <path-to-directory>`.

Then, use `./run_test.sh (compiler)` to compile all files found in the subdirectories of 
`/source/`. If no compiler is specified, it will use the current value of `$FC` (which is 
set to `f77` by default within the Makefile). Each time `./run_test.sh` is executed, the 
results will be saved to `/logs/` as `<compiler>.log`, where you can view all compilation 
errors encountered during runtime. 

Additionally, running `./compare.py <compiler1> <compiler2> (-g / -v)` from the base 
directory outputs differences in errors between two given compilers, where:

- no flag outputs differences without additional formatting,
- `-g` outputs only the first line in each error to a two-column table, and
- `-v` outputs errors in full to a two-column table.

This framework uses [BATS](https://github.com/bats-core/bats-core).

#### Other commands include: ####
- `./addMakefile.sh <directory>`: creates a generic Makefile. Runs automatically during `./run_test.sh` if no Makefile is present.
- `./addEnv.sh <directory>`: adds `comm.env` a subdirectory, allowing BATS to use necessary environmental commands before compiling, e.g. `module load`.
- `./cleanup.sh <directory>`: helps delete and/or backup your .bats files. This is important for when you add or delete files from a source directory.

#### Scraping: ####
Scrapers can be immensely helpful with gathering huge swathes of data at once. Within the 
`/scraping/` directory, one can find `scraper.py`, a Python program that can help pull multiple 
FORTRAN files from a single website. It will save these files to a directory named from the 
website. To use, type:

`./scraper.py <website url> <path to files> <extension>`
- `<website url>`: the page on which links to each FORTRAN file are found.
- `<path to files>`: URL path preceding each file link. Example below.
- `<extension>`: the desired extension to search for and download, e.g. `.f90`.

Using the example of [Michel Olagnon's ORDERPACK 2.0](http://www.fortran-2000.com/rank/index.html), one might execute:

`./scraper.py http://www.fortran-2000.com/rank/index.html http://www.fortran-2000.com/rank/ .f90`

In this example, `http://www.fortran-2000.com/rank/` is used as the second argument because 
every FORTRAN file on the page is located at `http://www.fortran-2000.com/rank/<FILE.f90>`.

Though not all of these need a scraper, other websites one might pull files from include:
- http://jean-pierre.moreau.pagesperso-orange.fr/f_function2.html
- http://www.fortran-2000.com/rank/index.html
- https://jblevins.org/mirror/amiller/
- https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html
- https://people.sc.fsu.edu/~jburkardt/f_src/stroud/stroud.html
- https://github.com/agforero/nSTREAM

#### Example outputs: ####

```
$ ./compare.py compiler1 compiler2 -g

                                               compiler1 | compiler2
-------------------------------------------------------- | --------------------------------------------------------
                                                         |
------------------------------- /fortran-testing-framework/source/sample-directory: -------------------------------
not ok error 1 not encountered in compiler2              |
not ok error 2 not encountered in compiler2              |
not ok error 3 not encountered in compiler2              |
                                                         | not ok error 4 not encountered in compiler1             
                                                         | not ok error 5 not encountered in compiler1             
                                                         | not ok error 6 not encountered in compiler1     
                                                         |
```
```
$ ./compare.py gfortran bleeding-edge-compiler -v

                                                gfortran | bleedingedgecompiler
-------------------------------------------------------- | --------------------------------------------------------
                                                         |
------------------------------- /fortran-testing-framework/source/sample-directory: -------------------------------
                                                         | --------------------------------------------------------
                                                         | not ok sample-directory: make -j all                 
                                                         |
                                                         | # (in test file tests.bats, line 5)                     
                                                         | #   `make -j all' failed with status 1                  
                                                         | #     file1.f90(1): here, the compiler might throw      
                                                         | #     some kind of random error. maybe a file is        
                                                         | #     missing, or a module isn't included, or an        
                                                         | #     integer is misplaced. or, if cross-checking       
                                                         | #     against a stable compiler like gfortran, you      
                                                         | #     might have found an error with the compiler       
                                                         | #     itself. if a line in the error message is too     
                                                         | #     long, the verbose table will display it like th...
                                                         |                                                         
                                                         |
                                                         |
------------------------------ /fortran-testing-framework/source/sample-directory-2: ------------------------------
                                                         | --------------------------------------------------------
                                                         | not ok sample-directory-2: make -j all               
                                                         |
                                                         | # (in test file tests.bats, line 5)                     
                                                         | #   `make -j all' failed with status 2                  
                                                         | #     with multiple errors across different   
                                                         | #     directories, you can see how the table uses       
                                                         | #     hyphens to neatly chop things up.                 
                                                         |                                                         
                                                         | --------------------------------------------------------
                                                         | not ok error 2 in sample-directory-2                   
                                                         |
                                                         | # (in test file tests.bats, line 5)                     
                                                         | #   `make -j all' failed with status 3                  
                                                         | #     in the event that there are multiple errors,      
                                                         | #     the table will also use hyphens to organize       
                                                         | #     then nicely while in verbose mode.                
                                                         |                                                         
                                                         |
```
