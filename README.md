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
- [beautifulsoup4 (for scraper)](https://pypi.org/project/beautifulsoup4/)



#### Summary: 
This is a framework built to help test bleeding-edge FORTRAN compilers. By testing the compilation of a wide variety of FORTRAN programs, and cross-checking results with stable compilers like `gfortran`, one can find where a compiler might be going wrong.

To add a directory containing FORTRAN programs to use for testing compilation, `cd` into `/main/`, and use `./addNew.sh <path-to-directory>`.

Then, use `./run_test.sh (compiler)` to compile all files found in the subdirectories of `/source/`. If no compiler is specified, it will use the current value of `$FC` (which is set to `f77` by default within the Makefile). Each time `./run_test.sh` is executed, the results will be saved to `/logs/` as `<compiler>.log`, where you can view all compilation errors encountered during runtime. 

Additionally, running `./compare.py <compiler1> <compiler2> (-b / -g / -v <col>)` from the base directory outputs differences in errors between two given compilers, where:

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



#### Example output:
```bash
$ ./compare.py gfortran flang
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
DIR: home/agustin/Projects/ftframework/source/bu
DIFF_COUNT: gfortran: 0 flang: 1
TOTAL_ERRORS: gfortran: 0 flang: 1

================================================================
COMPILER: flang
ERROR: real.o
# flang  -c real.f90
# F90-S-0081-Illegal selector - KIND parameter has unknown value for data type  (real.f90: 5)
# F90-S-0034-Syntax error at or near identifier q0 (real.f90: 14)
#   0 inform,   0 warnings,   2 severes, 0 fatal for MAIN
# make: *** [Makefile:151: real.o] Error 1

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
DIR: home/agustin/Projects/ftframework/source/nasa_f
DIFF_COUNT: gfortran: 0 flang: 10
TOTAL_ERRORS: gfortran: 0 flang: 10

================================================================
COMPILER: flang
ERROR: prcset.o
# flang  -c prcset.f
# flang-7: error: unable to execute command: Segmentation fault (core dumped)
# flang-7: error: Fortran frontend to LLVM command failed due to signal (use -v to see invocation)
# clang version 7.0.1 
# Target: x86_64-pc-linux-gnu
# Thread model: posix
# InstalledDir: /usr/bin
# flang-7: note: diagnostic msg: PLEASE submit a bug report to  and include the crash backtrace, preprocessed source, and associated run script.
# flang-7: note: diagnostic msg: Error generating preprocessed source(s) - no preprocessable inputs.
# make: *** [Makefile:433: prcset.o] Error 254

ERROR: pool.o
# flang  -c pool.f
# flang-7: error: unable to execute command: Segmentation fault (core dumped)
# flang-7: error: Fortran frontend to LLVM command failed due to signal (use -v to see invocation)
# clang version 7.0.1 
# Target: x86_64-pc-linux-gnu
# Thread model: posix
# InstalledDir: /usr/bin
# flang-7: note: diagnostic msg: PLEASE submit a bug report to  and include the crash backtrace, preprocessed source, and associated run script.
# flang-7: note: diagnostic msg: Error generating preprocessed source(s) - no preprocessable inputs.
# make: *** [Makefile:2179: pool.o] Error 254

ERROR: zzdskbsr.o
# flang  -c zzdskbsr.f
# /tmp/zzdskbsr-0fa576.ll:2936:12: error: '@_master___zzdskbsr__' defined with type 'void (i32, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64, i64*)*'
#         call void @_master___zzdskbsr__(i32 0, i64* %fname, i64* %bodyid, i64* %handle, i64* %cmpfun, i64* %usrctr, i64* %update, i64* %dladsc, i64* %dskdsc, i64* %found, i64* %cmpfun$sd, i64%.U0001.addr)
#                   ^
# 1 error generated.
# make: *** [Makefile:3667: zzdskbsr.o] Error 1

ERROR: zzgfdiu.o
# flang  -c zzgfdiu.f
# /tmp/zzgfdiu-739e1b.ll:313:12: error: '@_master___zzgfdiu__' defined with type 'void (i32, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64, i64, i64, i64*)*'
#         call void @_master___zzgfdiu__(i32 0, i64* %target, i64* %abcorr, i64* %obsrvr, i64* %udfunc, i64* %et, i64* %decres, i64* %dist, i64* %udfunc$sd, i64%.U0003.addr, i64%.U0004.addr, i64%.U0005.addr)
#                   ^
# 1 error generated.
# make: *** [Makefile:4297: zzgfdiu.o] Error 1

ERROR: zzgfcou.o
# flang  -c zzgfcou.f
# /tmp/zzgfcou-02d33d.ll:2110:12: error: '@_master___zzgfcou__' defined with type 'void (i32, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64, i64, i64, i64, i64, i64, i64, i64, i64, i64*)*'
#         call void @_master___zzgfcou__(i32 0, i64* %vecdef, i64* %method, i64* %target, i64* %et, i64* %ref, i64* %abcorr, i64* %obsrvr, i64* %dref, i64* %dvec, i64* %crdsys, i64* %crdnam, i64* %decres, i64* %crdval, i64* %crdfnd, i64* %udfunc, i64* %udfunc$sd, i64%.U0009.addr, i64%.U0010.addr, i64%.U0011.addr, i64%.U0012.addr, i64%.U0013.addr, i64%.U0014.addr, i64%.U0015.addr, i64%.U0016.addr, i64%.U0017.addr)
#                   ^
# 1 error generated.
# make: *** [Makefile:4285: zzgfcou.o] Error 1

ERROR: zzgfilu.o
# flang  -c zzgfilu.f
# /tmp/zzgfilu-0c1f3b.ll:788:12: error: '@_master___zzgfilu__' defined with type 'void (i32, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64, i64, i64, i64, i64, i64, i64, i64*)*'
#         call void @_master___zzgfilu__(i32 0, i64* %method, i64* %angtyp, i64* %target, i64* %illum, i64* %fixref, i64* %abcorr, i64* %obsrvr, i64* %spoint, i64* %et, i64* %udfunc, i64* %decres, i64* %angle, i64* %udfunc$sd, i64%.U0007.addr, i64%.U0008.addr, i64%.U0009.addr, i64%.U0010.addr, i64%.U0011.addr, i64%.U0012.addr, i64%.U0013.addr)
#                   ^
# 1 error generated.
# make: *** [Makefile:4306: zzgfilu.o] Error 1

ERROR: zzgfpau.o
# flang  -c zzgfpau.f
# /tmp/zzgfpau-22f585.ll:552:12: error: '@_master___zzgfpau__' defined with type 'void (i32, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64, i64, i64, i64, i64, i64*)*'
#         call void @_master___zzgfpau__(i32 0, i64* %target, i64* %illmn, i64* %abcorr, i64* %obsrvr, i64* %udfunc, i64* %et, i64* %decres, i64* %rvl, i64* %xtarg, i64* %xillmn, i64* %xabcor, i64* %xobs, i64* %xablk, i64* %udfunc$sd, i64%.U0005.addr, i64%.U0006.addr, i64%.U0007.addr, i64%.U0008.addr, i64%.U0009.addr)
#                   ^
# 1 error generated.
# make: *** [Makefile:4318: zzgfpau.o] Error 1

ERROR: zzgfrru.o
# flang  -c zzgfrru.f
# /tmp/zzgfrru-a1f258.ll:447:12: error: '@_master___zzgfrru__' defined with type 'void (i32, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64, i64, i64, i64, i64*)*'
#         call void @_master___zzgfrru__(i32 0, i64* %target, i64* %abcorr, i64* %obsrvr, i64* %dt, i64* %udfunc, i64* %et, i64* %decres, i64* %rvl, i64* %xtarg, i64* %xabcor, i64* %xobs, i64* %xdt, i64* %udfunc$sd, i64%.U0004.addr, i64%.U0005.addr, i64%.U0006.addr, i64%.U0007.addr)
#                   ^
# 1 error generated.
# make: *** [Makefile:4336: zzgfrru.o] Error 1

ERROR: zzgfspu.o
# flang  -c zzgfspu.f
# /tmp/zzgfspu-f6f7b3.ll:1041:12: error: '@_master___zzgfspu__' defined with type 'void (i32, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64, i64, i64, i64, i64, i64, i64, i64, i64*)*'
#         call void @_master___zzgfspu__(i32 0, i64* %of, i64* %from, i64* %shape, i64* %frame, i64* %et, i64* %udfunc, i64* %abcorr, i64* %decres, i64* %sep, i64* %xabcr, i64* %xbod, i64* %yref, i64* %xref, i64* %xobs, i64* %xrad, i64* %xshp, i64* %udfunc$sd, i64%.U0008.addr, i64%.U0009.addr, i64%.U0010.addr, i64%.U0011.addr, i64%.U0012.addr, i64%.U0013.addr, i64%.U0014.addr, i64%.U0015.addr)
#                   ^
# 1 error generated.
# make: *** [Makefile:4348: zzgfspu.o] Error 1

ERROR: tabrpt.o
# flang  -c tabrpt.f
# /tmp/tabrpt-ca1837.ll:1687:12: error: '@_master___tabrpt__' defined with type 'void (i32, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64*, i64, i64*)*'
#         call void @_master___tabrpt__(i32 0, i64* %nitems, i64* %item, i64* %size, i64* %width, i64* %justr, i64* %presrv, i64* %spcial, i64* %lmarge, i64* %space, i64* %fetch, i64* %fetch$sd, i64%.U0001.addr)
#                   ^
# 1 error generated.
# make: *** [Makefile:5314: tabrpt.o] Error 1
```

```bash
$ ./compare.py gfortran flang -b
================================================================
DIR: home/agustin/Projects/ftframework/source/bu
DIFF_COUNT: gfortran: 0 flang: 1
TOTAL_ERRORS: gfortran: 0 flang: 1

COMPILER: flang

real.o

================================================================
DIR: home/agustin/Projects/ftframework/source/nasa_f
DIFF_COUNT: gfortran: 0 flang: 10
TOTAL_ERRORS: gfortran: 0 flang: 10

COMPILER: flang

prcset.o
pool.o
zzdskbsr.o
zzgfdiu.o
zzgfcou.o
zzgfilu.o
zzgfpau.o
zzgfrru.o
zzgfspu.o
tabrpt.o
```

```bash
$ ./compare.py gfortran flang -g
            gfortran | flang               
0 -------------- /source/bu ------------- 1
                     | real.o
                     |                     
0 ----------- /source/nasa_f ----------- 10
                     | prcset.o
                     | pool.o
                     | zzdskbsr.o
                     | zzgfdiu.o
                     | zzgfcou.o
                     | zzgfilu.o
                     | zzgfpau.o
                     | zzgfrru.o
                     | zzgfspu.o
                     | tabrpt.o
                     |                     
```

```bash
$ ./compare.py gfortran flang -v 64
                                                        gfortran | flang                                                           
0 ---------------------------------------------------------- /source/bu --------------------------------------------------------- 1
                                                                 | real.o
                                                                 | 
                                                                 | # flang  -c real.f90
                                                                 | # F90-S-0081-Illegal selector - KIND parameter has unknown va...
                                                                 | # F90-S-0034-Syntax error at or near identifier q0 (real.f90:...
                                                                 | #   0 inform,   0 warnings,   2 severes, 0 fatal for MAIN
                                                                 | # make: *** [Makefile:151: real.o] Error 1
                                                                 | 
                                                                 | ----------------------------------------------------------------
                                                                 |                                                                 
0 ------------------------------------------------------- /source/nasa_f ------------------------------------------------------- 10
                                                                 | prcset.o
                                                                 | 
                                                                 | # flang  -c prcset.f
                                                                 | # flang-7: error: unable to execute command: Segmentation fau...
                                                                 | # flang-7: error: Fortran frontend to LLVM command failed due...
                                                                 | # clang version 7.0.1 
                                                                 | # Target: x86_64-pc-linux-gnu
                                                                 | # Thread model: posix
                                                                 | # InstalledDir: /usr/bin
                                                                 | # flang-7: note: diagnostic msg: PLEASE submit a bug report t...
                                                                 | # flang-7: note: diagnostic msg: Error generating preprocesse...
                                                                 | # make: *** [Makefile:433: prcset.o] Error 254
                                                                 | 
                                                                 | ----------------------------------------------------------------
                                                                 | pool.o
                                                                 | 
                                                                 | # flang  -c pool.f
                                                                 | # flang-7: error: unable to execute command: Segmentation fau...
                                                                 | # flang-7: error: Fortran frontend to LLVM command failed due...
                                                                 | # clang version 7.0.1 
                                                                 | # Target: x86_64-pc-linux-gnu
                                                                 | # Thread model: posix
                                                                 | # InstalledDir: /usr/bin
                                                                 | # flang-7: note: diagnostic msg: PLEASE submit a bug report t...
                                                                 | # flang-7: note: diagnostic msg: Error generating preprocesse...
                                                                 | # make: *** [Makefile:2179: pool.o] Error 254
                                                                 | 
                                                                 | ----------------------------------------------------------------
                                                                 | zzdskbsr.o
                                                                 | 
                                                                 | # flang  -c zzdskbsr.f
                                                                 | # /tmp/zzdskbsr-0fa576.ll:2936:12: error: '@_master___zzdskbs...
                                                                 | #         call void @_master___zzdskbsr__(i32 0, i64* %fname,...
                                                                 | #                   ^
                                                                 | # 1 error generated.
                                                                 | # make: *** [Makefile:3667: zzdskbsr.o] Error 1
                                                                 | 
                                                                 | ----------------------------------------------------------------
                                                                 | zzgfdiu.o
                                                                 | 
                                                                 | # flang  -c zzgfdiu.f
                                                                 | # /tmp/zzgfdiu-739e1b.ll:313:12: error: '@_master___zzgfdiu__...
                                                                 | #         call void @_master___zzgfdiu__(i32 0, i64* %target,...
                                                                 | #                   ^
                                                                 | # 1 error generated.
                                                                 | # make: *** [Makefile:4297: zzgfdiu.o] Error 1
                                                                 | 
                                                                 | ----------------------------------------------------------------
                                                                 | zzgfcou.o
                                                                 | 
                                                                 | # flang  -c zzgfcou.f
                                                                 | # /tmp/zzgfcou-02d33d.ll:2110:12: error: '@_master___zzgfcou_...
                                                                 | #         call void @_master___zzgfcou__(i32 0, i64* %vecdef,...
                                                                 | #                   ^
                                                                 | # 1 error generated.
                                                                 | # make: *** [Makefile:4285: zzgfcou.o] Error 1
                                                                 | 
                                                                 | ----------------------------------------------------------------
                                                                 | zzgfilu.o
                                                                 | 
                                                                 | # flang  -c zzgfilu.f
                                                                 | # /tmp/zzgfilu-0c1f3b.ll:788:12: error: '@_master___zzgfilu__...
                                                                 | #         call void @_master___zzgfilu__(i32 0, i64* %method,...
                                                                 | #                   ^
                                                                 | # 1 error generated.
                                                                 | # make: *** [Makefile:4306: zzgfilu.o] Error 1
                                                                 | 
                                                                 | ----------------------------------------------------------------
                                                                 | zzgfpau.o
                                                                 | 
                                                                 | # flang  -c zzgfpau.f
                                                                 | # /tmp/zzgfpau-22f585.ll:552:12: error: '@_master___zzgfpau__...
                                                                 | #         call void @_master___zzgfpau__(i32 0, i64* %target,...
                                                                 | #                   ^
                                                                 | # 1 error generated.
                                                                 | # make: *** [Makefile:4318: zzgfpau.o] Error 1
                                                                 | 
                                                                 | ----------------------------------------------------------------
                                                                 | zzgfrru.o
                                                                 | 
                                                                 | # flang  -c zzgfrru.f
                                                                 | # /tmp/zzgfrru-a1f258.ll:447:12: error: '@_master___zzgfrru__...
                                                                 | #         call void @_master___zzgfrru__(i32 0, i64* %target,...
                                                                 | #                   ^
                                                                 | # 1 error generated.
                                                                 | # make: *** [Makefile:4336: zzgfrru.o] Error 1
                                                                 | 
                                                                 | ----------------------------------------------------------------
                                                                 | zzgfspu.o
                                                                 | 
                                                                 | # flang  -c zzgfspu.f
                                                                 | # /tmp/zzgfspu-f6f7b3.ll:1041:12: error: '@_master___zzgfspu_...
                                                                 | #         call void @_master___zzgfspu__(i32 0, i64* %of, i64...
                                                                 | #                   ^
                                                                 | # 1 error generated.
                                                                 | # make: *** [Makefile:4348: zzgfspu.o] Error 1
                                                                 | 
                                                                 | ----------------------------------------------------------------
                                                                 | tabrpt.o
                                                                 | 
                                                                 | # flang  -c tabrpt.f
                                                                 | # /tmp/tabrpt-ca1837.ll:1687:12: error: '@_master___tabrpt__'...
                                                                 | #         call void @_master___tabrpt__(i32 0, i64* %nitems, ...
                                                                 | #                   ^
                                                                 | # 1 error generated.
                                                                 | # make: *** [Makefile:5314: tabrpt.o] Error 1
                                                                 | 
                                                                 | ----------------------------------------------------------------
                                                                 |                                                                 
```