#!/usr/bin/env python
import sys, os, glob
from textwrap import dedent

def main():
    t = open("tests.bats", 'w')
    t.write("# automatically generated .bats tests to check the compilation of each program in this directory\n\n")

    # setting gfortran to the FC flag should be left to ./run_test.sh, if at all.
    # we do, however, need to account for other environmental commands found in comm.env...
    env = glob.glob("comm.env")
    if len(env) == 1: # if we've found our file
        t.write("# preemptive environmental commands, as per comm.env:\n")
        e = open(env[0], 'r')
        for line in e:
            t.write(line + "\n")
        t.write("\n")

    flags = "-j -k -O --trace"
    if os.popen("make -v").read().split()[2][0] != '4': # if we can't use -O (or need -t) because too antiquated Make
        flags = flags[:-len(" -O --trace")]

    t.write(dedent(f"""\
    @test recompiling {os.path.basename(os.getcwd())} using $FC {{
    \tmake clean
    \tmake {flags} all
    \techo 
    }}
    """))
    t.close()

if __name__ == "__main__":
    main()