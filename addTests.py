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

def ignoreCommas(term):
    if term[-1] == ',': return term[:-1]
    else: return term    

def allExcept(allUsedMods, fc):
    if fc.range != [0, 0]: return allUsedMods[:fc.range[0]] + allUsedMods[fc.range[1]:]
    else: return allUsedMods

def findl_files():
    os.chdir("./")
    for ext in ("*.f90", "*.f", ".f95", ".f03", ".f08", ".for", ".f77", ".ftn"):
        if len(glob.glob(ext)) != 0: return glob.glob(ext)
        elif len(glob.glob(ext.upper())) != 0: return glob.glob(ext.upper())
    print("No FORTRAN files found.")
    sys.exit(1)

class FFile: # possibly a target, but only if it DOESN'T feed into anything else.
    def __init__(self, name, source):
        self.targetName = name
        self.src = source
        self.definedMods = [] # a list of mod names defined in this file. I should make anything with a main, but only make objects that aren't needed anywhere
        self.usedMods = [] # a list of mods this FFile uses
        self.program = False # initialize as False, since we don't know yet
        self.range = [0, 0] # aka, between x and y, this FFile holds those positions in the allUsedMods array

def main(): # god I love Python
    targets = {}
    allUsedMods = []
    os.chdir("./") # necessary?
    l_files = findl_files()
    l_fc = [] # list of FFile objects
    for file in l_files: # we need to find if the file has its own main
        try: # because sometimes we get a weird UnicodeDecodeError with some of the files.
            c = open(file, 'r') # either because it has a <program> call, 
            program = False # or it has lines independent of a module or interface or whatever
            firstUse = True
            l_fc.append(FFile(justTheName(file), file)) # we can rectify the targetName later depending on if its executable or object
            inMod, inInt, inSub, inFunc = 0, 0, 0, 0
            for line in c:
                try:
                    if line.split()[0][0] == '!' or line.split()[0][0] == 'C': continue
                    # the file is relied on by another file
                    if line.split()[0].lower() == "module":
                        if not inInt: inMod += 1
                        if inMod == 1: 
                            if line.split()[1].lower() != "procedure": l_fc[-1].definedMods.append(line.split()[1]) # add mod name to most recent FFile.definedMods[]

                    elif splitAndLower(line)[:2] == ["end", "module"]: inMod -= 1

                    # the file relies on another file
                    elif line.split()[0].lower() == "use":
                        if firstUse:
                            l_fc[-1].range[0] = len(allUsedMods)
                            l_fc[-1].range[1] = len(allUsedMods) # x and y of the range are set to the index AFTER the last element in allUsedMods
                            firstUse = False
                        for term in line.split()[1:]:
                            allUsedMods.append(ignoreCommas(term))
                            l_fc[-1].usedMods.append(ignoreCommas(term)) 
                            l_fc[-1].range[1] += 1 
                            # ^^^ increase y by 1. for example, if there's only ever one USE call in file where allUsedMods was previously len 10,
                            # range[] = [10, 11]. allUsedMods[10] is the only index inhabited by file, so for the next file, 
                            # allUsedMods[l_fc[-2].range[1]], aka allUsedMods[11], will be vacant, and thus the next x val. professor Olaf would be proud!!

                    # interface
                    elif line.split()[0].lower() == "interface": inInt += 1
                    elif splitAndLower(line)[:2] == ["end", "interface"]: inInt -= 1

                    # subroutine
                    elif line.split()[0] != "end" and "subroutine" in splitAndLower(line): inSub += 1
                    elif splitAndLower(line)[:2] == ["end", "subroutine"]: inSub -= 1

                    # function
                    elif line.split()[0] != "end" and "function" in splitAndLower(line): inFunc += 1
                    elif splitAndLower(line)[:2] == ["end", "function"]: inFunc -= 1

                    # is the thing an executable?
                    elif line.split()[0].lower() == "program" or ((inMod, inInt, inSub, inFunc) == (0, 0, 0, 0) and not line.isspace()):
                        program = True
                        l_fc[-1].program = True
                except: continue # the line is empty
            if program:
                targets[file] = justTheName(file) # for now, only add executables to targets[]
            else:
                l_fc[-1].targetName = justTheName(file) + ".o"
            c.close()
        except UnicodeDecodeError: print(f"UnicodeDecodeError in {file}; not adding to tests.bats.")
        except: print(f"{file} could not be read; not adding to tests.bats.")

    # then, add .o files that are not relied on anywhere else.
    for fc in l_fc:
        if not fc.program:
            bottomfeeder = True # assume innocent of not being at bottom of tree until proven guilty
            for mod in fc.definedMods:
                if mod in allExcept(allUsedMods, fc):
                    bottomfeeder = False # proven guilty.
            if bottomfeeder: # if innocent,
                targets[fc.targetName] = fc.targetName

    # executables should always be both COMPILED and LINKED, regardless of whether or not other files rely on them.
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
        \tmake {v}
        }}
        """))
    t.close()

if __name__ == "__main__":
    main()
