#!/usr/bin/env python
import sys

def justTheName(st): # helps generate things like <filename>.exe
    if len(st) == 0: return ''
    elif st[0] == '.': return ''
    else: return st[0] + justTheName(st[1:])

def main(): # god I love Python
    possibilities = [] # possible make <whatever> 
    allLines = [] # that is, all relevant lines
    m = open("Makefile", 'r')
    for line in m: 
        if line != "\n": allLines.append(line)
    for i in range(len(allLines)):
        if ':' in allLines[i] and not "clean" in allLines[i] and not "all" in allLines[i]: # exclude all and clean
            try: 
                if allLines[i+1][:1] == "\t": # if the line below the possible addition has a tab,
                    possibilities.append(allLines[i]) # then I'm pretty sure we've found a candidate
            except: continue # we've reached end of file
    m.close()

    # then we actually create the tests in our tests.bats
    t = open("tests.bats", 'w')
    t.write("# automatically generated .bats tests to check the compilation of each file in this directory.\n\n")

    # setting gfortran to the FC flag should be left to ./run_test.sh, if at all.

    for e in possibilities:
        t.write("@test make " + str(e.split()[0][:-1]) + " {\n")
        t.write("\trm -f " + justTheName(str(e.split()[0][:-1])) + ".o " + justTheName(str(e.split()[0][:-1])) + ".exe " + justTheName(str(e.split()[0][:-1])) + ".mod\n")
        t.write("\tmake " + str(e.split()[0][:-1]) + "\n")
        t.write("}\n\n")
    t.close()

if __name__ == "__main__":
    main()
    