# Fortran Testing Framework 
### A framework to help test FORTRAN compilers. Written by Agustin Forero for Argonne National Laboratory. 



#### A typical test might look like: 
```
$ ./run_test.sh gfortran
$ ./run_test.sh bleedingedgecompiler
$ ./compare.py gfortran bleedingedgecompiler
```



#### A flowchart describing runtime:
![Diagram](https://docs.google.com/drawings/d/e/2PACX-1vRsiWq6hmns9xMxgzui_GlJXo_F3Glt130am-J2Os0Dtb4hX1VsA9-vQ6hjQ3oMJJmnCq4-IKVzV0Zv/pub?w=1508&h=1110)



#### Requirements: 
- [Bash Automated Testing System (BATS)](https://github.com/bats-core/bats-core)
- [Python 3 or above](https://docs.python-guide.org/starting/install3/linux/)
- [GNU Make 4 or above](http://ftp.gnu.org/gnu/make/)



#### Summary: 
This is a framework built to help test bleeding-edge FORTRAN compilers. By testing the compilation of a wide variety of FORTRAN programs, and cross-checking results with stable compilers like `gfortran`, one can find where a compiler might be going wrong.

To add a directory containing FORTRAN programs to use for testing compilation, `cd` into `/main/`, and use `./addNew.sh <path-to-directory>`.

Then, use `./run_test.sh (compiler)` to compile all files found in the subdirectories of `/source/`. If no compiler is specified, it will use the current value of `$FC` (which is set to `f77` by default within the Makefile). Each time `./run_test.sh` is executed, the results will be saved to `/logs/` as `<compiler>.log`, where you can view all compilation errors encountered during runtime. 

Additionally, running `./compare.py <compiler1> <compiler2> (-g / -v <col>)` from the base directory outputs differences in errors between two given compilers, where:

- no flag outputs differences without additional formatting,
- `-b` outputs a very basic summary,
- `-g` outputs only the first line in each error to a two-column table, and
- `-v <col>` outputs errors in full to a two-column table with `<col>` columns on each side, up to 256.



#### Other commands include: 
- `./addMakefile.sh <directory>`: creates a generic Makefile. Runs automatically during `./run_test.sh` if no Makefile is present.
- `./verify.sh`: checks all subdirectories of `/source/` for possible race conditions.
- `./addEnv.sh <directory>`: adds `comm.env` to a subdirectory, allowing BATS to use necessary environmental commands before compiling, e.g. `module load`.
- `./cleanup.sh <directory>`: helps delete and/or backup your .bats files. This is important for when you add or delete files from a source directory.
- `./resetAll.sh`: deletes all instances of `tests.bats` across all directories -- that is, the BATS file generated during `./run_test.sh`.



#### Scraping: 
Scrapers can be immensely helpful with gathering huge swathes of data at once. Within the `/scraping/` directory, one can find `scraper.py`, a Python program that can help pull multiple FORTRAN files from a single website. It will then save these files to a directory named from the website, at which point the user can use `./addNew.sh ../scraper/<directory>` in `/main/` to add the files. To use, type:

`./scraper.py <website url> <path to files> <extension>`
- `<website url>`: the page on which links to each FORTRAN file are found.
- `<path to files>`: URL path preceding each file location. Example below.
- `<extension>`: the desired extension to search for and download, e.g. `.f90`.

Using the example of [Michel Olagnon's ORDERPACK 2.0](http://www.fortran-2000.com/rank/index.html), one might execute:

`./scraper.py http://www.fortran-2000.com/rank/index.html http://www.fortran-2000.com/rank/ .f90`

In this example, `http://www.fortran-2000.com/rank/` is used as the second argument because every FORTRAN file on the page is located at `http://www.fortran-2000.com/rank/<FILE.f90>`.

Though not all of these need a scraper, other websites one might pull FORTRAN files from include:
- http://jean-pierre.moreau.pagesperso-orange.fr/f_function2.html
- http://www.fortran-2000.com/rank/index.html
- https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html
- https://people.sc.fsu.edu/~jburkardt/f_src/stroud/stroud.html
- https://github.com/agforero/nSTREAM



#### A warning against race conditions: 
Use `./verify.py` from `/main/` to check for possible race conditions!

For the files being compiled, it's helpful that all files involved would already `make` perfectly with the `-j` flag. If there's a problem with the code itself, both compilers will encounter it, which should be alright; but sometimes, race conditions present within directories can cause sporadic errors that do not occur every runtime, thereby rendering the testing inconsistent and therefore invalid.

For example, if a directory has multiple declarations of the same module across different files, multiple instances of `make` might try to write to `<modulename>.mod` at once during `make`, thereby allowing for a race condition. This mod will subsequently be corrupted when referenced later on.

In short, [race conditions](https://www.youtube.com/watch?v=7aF0q7NfwfA) are bad.



#### Sample output:
```
$ ./compare.py gfortran ifx -b
================================================================
DIR: home/agforero/ftframework/source/stroud
DIFF_COUNT: gfortran: 7 ifx: 24
TOTAL_ERRORS: gfortran: 145 ifx: 162

COMPILER: gfortran

genin.o
niederreiter.o
triangle_lyness_rule.o
mxm2_openmp.exe
mxm_openmp.exe
sgefa_openmp.exe
sparse_grid_mixed_size_table.exe

COMPILER: ifx

asa091.o
besselj.o
cdflib.o
corkscrew_plot3d.o
fair_dice_simulation.o
fd1d_heat_implicit.o
fd_predator_prey.o
fem2d_sample.o
file_transpose.o
mxm2_openmp.o
mxm_openmp.o
power_method.o
nms.o
predator_plot3d.o
polpak.o
prob.o
sgefa_openmp.o
shallow_water_1d.o
sparse_grid_mixed_size_table.o
test_nearest.o
subset.o
timer_etime.o
triangulation_node_to_element.o
triangle_lyness_rule.exe
```

```
$ ./compare.py gfortran ifx -g
                        gfortran | ifx                             
7 ----------------------- /source/stroud ----------------------- 24
                         genin.o | asa091.o                        
                  niederreiter.o | besselj.o                       
          triangle_lyness_rule.o | cdflib.o                        
                 mxm2_openmp.exe | corkscrew_plot3d.o              
                  mxm_openmp.exe | fair_dice_simulation.o          
                sgefa_openmp.exe | fd1d_heat_implicit.o            
sparse_grid_mixed_size_table.exe | fd_predator_prey.o              
                                 | fem2d_sample.o
                                 | file_transpose.o
                                 | mxm2_openmp.o
                                 | mxm_openmp.o
                                 | power_method.o
                                 | nms.o
                                 | predator_plot3d.o
                                 | polpak.o
                                 | prob.o
                                 | sgefa_openmp.o
                                 | shallow_water_1d.o
                                 | sparse_grid_mixed_size_table.o
                                 | test_nearest.o
                                 | subset.o
                                 | timer_etime.o
                                 | triangulation_node_to_element.o
                                 | triangle_lyness_rule.exe
```

```
$ ./compare.py gfortran ifx -v 50
                                          gfortran | ifx                                               
7 ----------------------------------------- /source/stroud ----------------------------------------- 24
genin.o                                            | asa091.o                                          
                                                   |                                                   
  gfortran -O2 -fopenmp -c genin.f90               |   ifx -O2 -fopenmp -c asa091.f90                  
  genin.f90:535:21:                                |   asa091.f90(553): error  6404: This name does ...
        |                     1                    |       arg = p * log ( x ) - x - lgamma ( p + 1....
  Error: Actual argument contains too few eleme... |   ------------------------------^                 
  make: *** [Makefile:1189: genin.o] Error 1       |   compilation aborted for asa091.f90 (code 1)     
                                                   |   make: *** [Makefile:73: asa091.o] Error 1       
-------------------------------------------------- |                                                   
niederreiter.o                                     | --------------------------------------------------
                                                   | besselj.o                                         
  gfortran -O2 -fopenmp -c niederreiter.f90        |                                                   
  niederreiter.f90:220:21:                         |   ifx -O2 -fopenmp -c besselj.f90                 
        |                     1                    |   besselj.f90(1860): error  6404: This name doe...
  Error: Actual argument contains too few eleme... |     DLNGAM = LOG (ABS (DGAMMA(X)) )               
  make: *** [Makefile:1897: niederreiter.o] Err... |   ---------------------^                          
                                                   |   besselj.f90(1860): error  6362: The data type...
-------------------------------------------------- |     DLNGAM = LOG (ABS (DGAMMA(X)) )               
triangle_lyness_rule.o                             |   ---------------------^                          
                                                   |   besselj.f90(1860): error  6362: The data type...
  gfortran -O2 -fopenmp -c triangle_lyness_rule... |     DLNGAM = LOG (ABS (DGAMMA(X)) )               
  triangle_lyness_rule.f90:1:2:                    |   ----------------^                               
        |  1                                       |   compilation aborted for besselj.f90 (code 1)    
  Error: Non-numeric character in statement lab... |   make: *** [Makefile:181: besselj.o] Error 1     
  make: *** [Makefile:3364: triangle_lyness_rul... |                                                   
                                                   | --------------------------------------------------
-------------------------------------------------- | cdflib.o                                          
mxm2_openmp.exe                                    |                                                   
                                                   |   ifx -O2 -fopenmp -c cdflib.f90                  
  gfortran -O2 -fopenmp mxm2_openmp.o z_sample_... |   cdflib.f90(7861): catastrophic error: **Inter...
  /soft/packaging/spack-builds/linux-rhel7-x86_... |   compilation aborted for cdflib.f90 (code 1)     
  z_sample_st.f90:(.text.startup+0x0): multiple... |   make: *** [Makefile:334: cdflib.o] Error 1      
  /soft/packaging/spack-builds/linux-rhel7-x86_... |                                                   
  z_sample_st.f90:(.text+0x194b): undefined ref... | --------------------------------------------------
  /soft/packaging/spack-builds/linux-rhel7-x86_... | corkscrew_plot3d.o                                
  /soft/packaging/spack-builds/linux-rhel7-x86_... |                                                   
  collect2: error: ld returned 1 exit status       |   ifx -O2 -fopenmp -c corkscrew_plot3d.f90        
  make: *** [Makefile:1849: mxm2_openmp.exe] Er... |   corkscrew_plot3d.f90(134): error  6353: A RET...
                                                   |     return                                        
-------------------------------------------------- |   --^                                             
mxm_openmp.exe                                     |   compilation aborted for corkscrew_plot3d.f90 ...
                                                   |   make: *** [Makefile:481: corkscrew_plot3d.o] ...
  gfortran -O2 -fopenmp mxm_openmp.o z_sample_s... |                                                   
  /soft/packaging/spack-builds/linux-rhel7-x86_... | --------------------------------------------------
  z_sample_st.f90:(.text.startup+0x0): multiple... | fair_dice_simulation.o                            
  /soft/packaging/spack-builds/linux-rhel7-x86_... |                                                   
  z_sample_st.f90:(.text+0x194b): undefined ref... |   ifx -O2 -fopenmp -c fair_dice_simulation.f90    
  /soft/packaging/spack-builds/linux-rhel7-x86_... |   fair_dice_simulation.f90(147): error  6353: A...
  /soft/packaging/spack-builds/linux-rhel7-x86_... |     return                                        
  collect2: error: ld returned 1 exit status       |   --^                                             
  make: *** [Makefile:1861: mxm_openmp.exe] Err... |   compilation aborted for fair_dice_simulation....
                                                   |   make: *** [Makefile:847: fair_dice_simulation...
-------------------------------------------------- |                                                   
sgefa_openmp.exe                                   | --------------------------------------------------
                                                   | fd1d_heat_implicit.o                              
  gfortran -O2 -fopenmp sgefa_openmp.o z_sample... |                                                   
  /soft/packaging/spack-builds/linux-rhel7-x86_... |   ifx -O2 -fopenmp -c fd1d_heat_implicit.f90      
  z_sample_st.f90:(.text+0x0): multiple definit... |   fd1d_heat_implicit.f90(256): error  6353: A R...
  /soft/packaging/spack-builds/linux-rhel7-x86_... |     return                                        
  z_sample_st.f90:(.text.startup+0x0): multiple... |   --^                                             
  /soft/packaging/spack-builds/linux-rhel7-x86_... |   compilation aborted for fd1d_heat_implicit.f9...
  z_sample_st.f90:(.text+0x194b): undefined ref... |   make: *** [Makefile:901: fd1d_heat_implicit.o...
  /soft/packaging/spack-builds/linux-rhel7-x86_... |                                                   
  /soft/packaging/spack-builds/linux-rhel7-x86_... | --------------------------------------------------
  collect2: error: ld returned 1 exit status       | fd_predator_prey.o                                
  make: *** [Makefile:2614: sgefa_openmp.exe] E... |                                                   
                                                   |   ifx -O2 -fopenmp -c fd_predator_prey.f90        
-------------------------------------------------- |   fd_predator_prey.f90(238): error  6353: A RET...
sparse_grid_mixed_size_table.exe                   |     return                                        
                                                   |   --^                                             
  gfortran -O2 -fopenmp sparse_grid_mixed_size_... |   compilation aborted for fd_predator_prey.f90 ...
  /soft/packaging/spack-builds/linux-rhel7-x86_... |   make: *** [Makefile:928: fd_predator_prey.o] ...
  sparse_grid_mixed_size_table.f90:(.text+0x68e... |                                                   
  /soft/packaging/spack-builds/linux-rhel7-x86_... | --------------------------------------------------
  /soft/packaging/spack-builds/linux-rhel7-x86_... | fem2d_sample.o                                    
  /soft/packaging/spack-builds/linux-rhel7-x86_... |                                                   
  sparse_grid_mixed_size_table.f90:(.text+0x115... |   ifx -O2 -fopenmp -c fem2d_sample.f90            
  /soft/packaging/spack-builds/linux-rhel7-x86_... |   fem2d_sample.f90(203): error  6353: A RETURN ...
  /soft/packaging/spack-builds/linux-rhel7-x86_... |       return                                      
  /soft/packaging/spack-builds/linux-rhel7-x86_... |   ----^                                           
  /soft/packaging/spack-builds/linux-rhel7-x86_... |   compilation aborted for fem2d_sample.f90 (cod...
  /soft/packaging/spack-builds/linux-rhel7-x86_... |   make: *** [Makefile:1021: fem2d_sample.o] Err...
  collect2: error: ld returned 1 exit status       |                                                   
  make: *** [Makefile:2728: sparse_grid_mixed_s... | --------------------------------------------------
                                                   | file_transpose.o                                  
-------------------------------------------------- |                                                   
                                                   |   ifx -O2 -fopenmp -c file_transpose.f90
                                                   |   file_transpose.f90(92): error  6353: A RETURN...
                                                   |     return
                                                   |   --^
                                                   |   compilation aborted for file_transpose.f90 (c...
                                                   |   make: *** [Makefile:1114: file_transpose.o] E...
                                                   | 
                                                   | --------------------------------------------------
                                                   | mxm2_openmp.o
                                                   | 
                                                   |   ifx -O2 -fopenmp -c mxm2_openmp.f90
                                                   |   mxm2_openmp.f90(55): error  6353: A RETURN st...
                                                   |     return
                                                   |   --^
                                                   |   compilation aborted for mxm2_openmp.f90 (code 1)
                                                   |   make: *** [Makefile:1852: mxm2_openmp.o] Error 1
                                                   | 
                                                   | --------------------------------------------------
                                                   | mxm_openmp.o
                                                   | 
                                                   |   ifx -O2 -fopenmp -c mxm_openmp.f90
                                                   |   mxm_openmp.f90(51): error  6353: A RETURN sta...
                                                   |     return
                                                   |   --^
                                                   |   compilation aborted for mxm_openmp.f90 (code 1)
                                                   |   make: *** [Makefile:1864: mxm_openmp.o] Error 1
                                                   | 
                                                   | --------------------------------------------------
                                                   | power_method.o
                                                   | 
                                                   |   ifx -O2 -fopenmp -c power_method.f90
                                                   |   power_method.f90(383): error  6404: This name...
                                                   |       lambda = complex ( lambda_real, lambda_im...
                                                   |   -------------^
                                                   |   compilation aborted for power_method.f90 (cod...
                                                   |   make: *** [Makefile:2074: power_method.o] Err...
                                                   | 
                                                   | --------------------------------------------------
                                                   | nms.o
                                                   | 
                                                   |   ifx -O2 -fopenmp -c nms.f90
                                                   |   nms.f90(11512): catastrophic error: **Interna...
                                                   |   compilation aborted for nms.f90 (code 1)
                                                   |   make: *** [Makefile:1924: nms.o] Error 1
                                                   | 
                                                   | --------------------------------------------------
                                                   | predator_plot3d.o
                                                   | 
                                                   |   ifx -O2 -fopenmp -c predator_plot3d.f90
                                                   |   predator_plot3d.f90(130): error  6353: A RETU...
                                                   |     return
                                                   |   --^
                                                   |   compilation aborted for predator_plot3d.f90 (...
                                                   |   make: *** [Makefile:2101: predator_plot3d.o] ...
                                                   | 
                                                   | --------------------------------------------------
                                                   | polpak.o
                                                   | 
                                                   |   ifx -O2 -fopenmp -c polpak.f90
                                                   |   polpak.f90(3787): error  6404: This name does...
                                                   |     facn = lgamma ( arg )
                                                   |   ---------^
                                                   |   polpak.f90(4235): error  6404: This name does...
                                                   |     prod = complex ( 1.0D+00, 0.0D+00 )
                                                   |   ---------^
                                                   |   polpak.f90(13365): error  6404: This name doe...
                                                   |     r8_beta = exp ( lgamma ( x ) &
                                                   |   ------------------^
                                                   |   polpak.f90(13367): warning  7319: This argume...
                                                   |                   - lgamma ( x + y ) )
                                                   |   ----------------^
                                                   |   polpak.f90(13365): error  6404: This name doe...
                                                   |     r8_beta = exp ( lgamma ( x ) &
                                                   |   ------------^
                                                   |   polpak.f90(13447): error  6404: This name doe...
                                                   |       facn = lgamma ( arg )
                                                   |   -----------^
                                                   |   compilation aborted for polpak.f90 (code 1)
                                                   |   make: *** [Makefile:2047: polpak.o] Error 1
                                                   | 
                                                   | --------------------------------------------------
                                                   | prob.o
                                                   | 
                                                   |   ifx -O2 -fopenmp -c prob.f90
                                                   |   prob.f90(2789): error  6404: This name does n...
                                                   |     beta_log = lgamma ( p ) + lgamma ( q ) - lg...
                                                   |   -------------^
                                                   |   prob.f90(7663): error  6404: This name does n...
                                                   |     g = lgamma ( a / 2.0D+00 )
                                                   |   ------^
                                                   |   prob.f90(9889): error  6404: This name does n...
                                                   |       pdf = real ( cnk * dnmk, kind = 8 ) / lga...
                                                   |   ------------------------------------------^
                                                   |   prob.f90(10936): error  6404: This name does ...
                                                   |       - lgamma ( c_sum + real ( a, kind = 8 ) ) &
                                                   |   ------^
                                                   |   prob.f90(13816): error  6404: This name does ...
                                                   |       pdf = y ** ( c - 1 ) / ( b * lgamma ( rea...
                                                   |   ---------------------------------^
                                                   |   prob.f90(20227): error  6404: This name does ...
                                                   |         lgamma ( real ( n + 1, kind = 8 ) ) &
                                                   |   ------^
                                                   |   prob.f90(24471): error  6404: This name does ...
                                                   |     facn = lgamma ( real ( n + 1, kind = 8 ) )
                                                   |   ---------^
                                                   |   prob.f90(24800): error  6404: This name does ...
                                                   |     pdf_log = lgamma ( real ( a + 1, kind = 8 ) )
                                                   |   ------------^
                                                   |   prob.f90(30779): error  6404: This name does ...
                                                   |     r8_beta = exp ( lgamma ( a ) + lgamma ( b )...
                                                   |   ------------------^
                                                   |   prob.f90(30779): warning  7319: This argument...
                                                   |     r8_beta = exp ( lgamma ( a ) + lgamma ( b )...
                                                   |   ----------------------------------------------^
                                                   |   prob.f90(30779): error  6404: This name does ...
                                                   |     r8_beta = exp ( lgamma ( a ) + lgamma ( b )...
                                                   |   ------------^
                                                   |   prob.f90(31553): error  6404: This name does ...
                                                   |       arg = p * log ( x ) - x - lgamma ( p + 1....
                                                   |   ------------------------------^
                                                   |   prob.f90(31935): error  6404: This name does ...
                                                   |     r8_gamma_log_int = lgamma ( real ( n, kind ...
                                                   |   ---------------------^
                                                   |   prob.f90(36520): error  6404: This name does ...
                                                   |       a = sqrt ( 0.5D+00 * f ) * exp ( lgamma (...
                                                   |   -------------------------------------^
                                                   |   prob.f90(36521): error  6362: The data types ...
                                                   |         - lgamma ( 0.5D+00 * f ) ) * d
                                                   |   ------^
                                                   |   compilation aborted for prob.f90 (code 1)
                                                   |   make: *** [Makefile:2134: prob.o] Error 1
                                                   | 
                                                   | --------------------------------------------------
                                                   | sgefa_openmp.o
                                                   | 
                                                   |   ifx -O2 -fopenmp -c sgefa_openmp.f90
                                                   |   sgefa_openmp.f90(88): error  6353: A RETURN s...
                                                   |     return
                                                   |   --^
                                                   |   sgefa_openmp.f90(153): remark  8291: Recommen...
                                                   |     write ( *, '(a,i8,2x,g10.4,2x,f10.4)' ) &
                                                   |   ------------------------^
                                                   |   sgefa_openmp.f90(226): remark  8291: Recommen...
                                                   |     write ( *, '(a,i8,2x,g10.4,2x,f10.4)' ) &
                                                   |   ------------------------^
                                                   |   sgefa_openmp.f90(299): remark  8291: Recommen...
                                                   |     write ( *, '(a,i8,2x,g10.4,2x,f10.4)' ) &
                                                   |   ------------------------^
                                                   |   sgefa_openmp.f90(372): remark  8291: Recommen...
                                                   |     write ( *, '(a,i8,2x,g10.4,2x,f10.4)' ) &
                                                   |   ------------------------^
                                                   |   sgefa_openmp.f90(445): remark  8291: Recommen...
                                                   |     write ( *, '(a,i8,2x,g10.4,2x,f10.4)' ) &
                                                   |   ------------------------^
                                                   |   compilation aborted for sgefa_openmp.f90 (cod...
                                                   |   make: *** [Makefile:2617: sgefa_openmp.o] Err...
                                                   | 
                                                   | --------------------------------------------------
                                                   | shallow_water_1d.o
                                                   | 
                                                   |   ifx -O2 -fopenmp -c shallow_water_1d.f90
                                                   |   shallow_water_1d.f90(281): error  6353: A RET...
                                                   |     return
                                                   |   --^
                                                   |   compilation aborted for shallow_water_1d.f90 ...
                                                   |   make: *** [Makefile:2635: shallow_water_1d.o]...
                                                   | 
                                                   | --------------------------------------------------
                                                   | sparse_grid_mixed_size_table.o
                                                   | 
                                                   |   ifx -O2 -fopenmp -c sparse_grid_mixed_size_ta...
                                                   |   sparse_grid_mixed_size_table.f90(288): error ...
                                                   |     return
                                                   |   --^
                                                   |   compilation aborted for sparse_grid_mixed_siz...
                                                   |   make: *** [Makefile:2731: sparse_grid_mixed_s...
                                                   | 
                                                   | --------------------------------------------------
                                                   | test_nearest.o
                                                   | 
                                                   |   ifx -O2 -fopenmp -c test_nearest.f90
                                                   |   test_nearest.f90(113): error  6353: A RETURN ...
                                                   |     return
                                                   |   --^
                                                   |   compilation aborted for test_nearest.f90 (cod...
                                                   |   make: *** [Makefile:3043: test_nearest.o] Err...
                                                   | 
                                                   | --------------------------------------------------
                                                   | subset.o
                                                   | 
                                                   |   ifx -O2 -fopenmp -c subset.f90
                                                   |   subset.f90(17400): error  6404: This name doe...
                                                   |     facn = lgamma ( arg )
                                                   |   ---------^
                                                   |   subset.f90(22720): error  6404: This name doe...
                                                   |       facn = lgamma ( arg )
                                                   |   -----------^
                                                   |   compilation aborted for subset.f90 (code 1)
                                                   |   make: *** [Makefile:2935: subset.o] Error 1
                                                   | 
                                                   | --------------------------------------------------
                                                   | timer_etime.o
                                                   | 
                                                   |   ifx -O2 -fopenmp -c timer_etime.f90
                                                   |   timer_etime.f90(134): error  6404: This name ...
                                                   |         time1 = etime ( tarray )
                                                   |   --------------^
                                                   |   compilation aborted for timer_etime.f90 (code 1)
                                                   |   make: *** [Makefile:3199: timer_etime.o] Error 1
                                                   | 
                                                   | --------------------------------------------------
                                                   | triangulation_node_to_element.o
                                                   | 
                                                   |   ifx -O2 -fopenmp -c triangulation_node_to_ele...
                                                   |   triangulation_node_to_element.f90(223): error...
                                                   |     return
                                                   |   --^
                                                   |   compilation aborted for triangulation_node_to...
                                                   |   make: *** [Makefile:3442: triangulation_node_...
                                                   | 
                                                   | --------------------------------------------------
                                                   | triangle_lyness_rule.exe
                                                   | 
                                                   |   ifx -O2 -fopenmp triangle_lyness_rule.o -o tr...
                                                   |   /soft/packaging/spack-builds/linux-rhel7-x86_...
                                                   |   for_main.c:(.text+0x2e): undefined reference ...
                                                   |   make: *** [Makefile:3361: triangle_lyness_rul...
                                                   | 
                                                   | --------------------------------------------------
                                                   |                                                   
```