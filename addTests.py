#!/usr/bin/env python
import sys, os, glob
from textwrap import dedent

def justTheName(st): # helps generate things like <filename>.exe
    if len(st) == 0: return ''
    elif st[0] == '.': return ''
    else: return st[0] + justTheName(st[1:])

def splitAndLower(line): 
    ls = line.split()
    for i in range(len(ls)): ls[i] = ls[i].lower()
    return ls

def main(): # god I love Python
    targets = {}
    os.chdir("./") # necessary?
    l_f90 = glob.glob("*.f90") # so far this only works with .f90 files
    for file in l_f90: # we need to find if the file has its own main
        c = open(file, 'r') # either because it has a <program> call, 
        program = False # or it has lines independent of a module or interface declaration
        inMod = 0
        inInt = 0
        inSub = 0
        inFunc = 0
        for line in c:
            try:
                # the file is relied on by another file
                if line.split()[0].lower() == "module":
                    if not inInt: inMod += 1
                elif splitAndLower(line)[:2] == ["end", "module"]: inMod -= 1

                # interface
                elif line.split()[0].lower() == "interface": inInt += 1
                elif splitAndLower(line)[:2] == ["end", "interface"]: inInt -= 1

                # subroutine
                elif line.split()[0].lower() == "subroutine": inSub += 1
                elif splitAndLower(line)[:2] == ["end", "subroutine"]: inSub -= 1

                # function
                elif line.split()[0].lower() == "function": inFunc += 1
                elif splitAndLower(line)[:2] == ["end", "function"]: inFunc -= 1

                # is the thing an executable?
                elif line.split()[0].lower() == "program" and not (inMod or inSub or inFunc) and not line.isspace():
                    program = True
            except: continue # the line is empty
        if program:
            targets[file] = justTheName(file)
        else:
            targets[file] = justTheName(file) + ".o"
        c.close()

    # then we determine if the Makefile calls us to make <filename>.exe, or just <filename>
    exe = False
    m = open("Makefile", 'r')
    inClean = False # this has to NOT include anything in a <clean> rule, since we otherwise might try to make something like *.exe
    for line in m:
        if line[:5].lower() == "clean":
            inClean = True
        if inClean and len(line.split()) != 0 and line[:1] != "\t" and line[:5].lower() != "clean":
            inClean = False
        if not inClean:
            for word in line.split():
                if word[-4:] == ".exe" or word[-5:] == ".exe:":
                    targets[justTheName(word) + ".f90"] = justTheName(word) + ".exe"
    m.close()

    # then we actually create the tests in our tests.bats
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

    for v in targets.values():
        t.write(dedent(f"""\
        @test make {v} {{
        \trm -f {justTheName(v)}.o {justTheName(v)}.exe {justTheName(v)}.mod
        \tmake {v}
        }}
        """))
    t.close()

if __name__ == "__main__":
    main()