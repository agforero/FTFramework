#!/usr/bin/env python
import sys, os
import linecache as lc 
from textwrap import dedent

def justTheName(st): # helps generate things like <filename>.exe
    if len(st) == 0: return ''
    elif st[0] == '.': return ''
    else: return st[0] + justTheName(st[1:])

def splitAndLower(line): # split and lower? I hardly even...
    ls = line.split()
    for i in range(len(ls)): ls[i] = ls[i].lower()
    return ls

def correctErrorBody(body, target, compiler):
    ret = ""
    targetFound = False
    for line in body.split("\n"):
        for word in line.split():
            if justTheName(word) == justTheName(target) and compiler in line.split():
                targetFound = True
        if targetFound:
            ret += line + '\n'
    return ret

def buildDirectory(osarg):
    try: f = open(f"../logs/{justTheName(sys.argv[osarg])}.log", 'r')
    except FileNotFoundError:
        print(f"Error: logs/{justTheName(sys.argv[osarg])}.log does not exist.")
        sys.exit(2)

    ret = {} 
    f_compiler = lc.getline(f"../logs/{justTheName(sys.argv[osarg])}.log", 1).split()[2]
    currentDirectory = ""
    inError = False
    lines = f.readlines()
    for i in range(1, len(lines)): # by definition, any execution of <compiler> must start with it as the first term
        if len(lines[i].split()) > 0 and lines[i].split()[0] == "DIR:":
            currentDirectory = lines[i][6:-1] # ignore DIR: and newline
            ret[currentDirectory] = {} # initialize as empty dict
        elif len(lines[i-1].split()) < 2 or len(lines[i].split()) < 2: continue # this should also cover whitespace, right?
        elif inError and splitAndLower(lines[i-1])[:3] == ["#", "make:", "***"]:
            inError = False
            ret[currentDirectory][lines[i-1].split()[-3][:-1]] = correctErrorBody(errorBody, lines[i-1].split()[-3][:-1], f_compiler)
            errorBody = ""
        elif inError:
            errorBody += lines[i]
        elif lines[i-1].split()[1] == f_compiler and lines[i].split()[1] != f_compiler:
            inError = True
            errorBody = lines[i-1] + lines[i]

    return ret

def displayDirectory(osarg, drct): # mainly for testing
    lc.getline(f"../logs/{justTheName(sys.argv[osarg])}.log", 1).split()[2]
    for d in drct.keys():
        print(f"DIR: {d}")
        for e in drct[d].keys():
            print(f"TARGET: {e}\nERROR:")
            print(drct[d][e])

def folderName(path):
    for i in range(len(path)):
        if path[i:i + len("/source")] == "/source":
            return path[i:]
    return path

def rawOutput(master):
    compilers = list(master.keys())
    for dr in master[compilers[0]].keys():
        display1, display2 = False, False # seeing if we even need to display anything
        err1Count, err2Count = 0, 0
        for e in list(master[compilers[0]][dr].keys()):
            if e not in list(master[compilers[1]][dr].keys()):
                display1 = True
                err1Count += 1
        for e in list(master[compilers[1]][dr].keys()):
            if e not in list(master[compilers[0]][dr].keys()):
                display2 = True
                err2Count += 1

        if display1 or display2: 
            print(dedent(f"""\
            {'/' * 64}
            DIR: {dr}
            DIFF_COUNT: {compilers[0]}: {err1Count} {compilers[1]}: {err2Count}
            """))
            
        if display1:
            print(f"{'=' * 32}\nCOMPILER: {compilers[0]}")
            for e in list(master[compilers[0]][dr].keys()):
                if e not in list(master[compilers[1]][dr].keys()):
                    print(f"ERROR: {e}")
                    print(master[compilers[0]][dr][e])
        if display2:
            print(f"{'=' * 32}\nCOMPILER: {compilers[1]}")
            for e in list(master[compilers[1]][dr].keys()):
                if e not in list(master[compilers[0]][dr].keys()):
                    print(f"ERROR: {e}")
                    print(master[compilers[1]][dr][e])

def main():
    version = os.popen("make -v").read().split("\n")[0]
    if version.split()[2][0] != '4':
        print(f"GNU Make 4.x required for compare.py; this shell is currently running {version}.")
        sys.exit(3)

    if len(sys.argv) == 4:
        if sys.argv != "-g" or sys.argv != "-v":
            print(f"{sys.argv[3]}: Invalid flag.")
            sys.exit(1)
    elif len(sys.argv) != 3:
        print(dedent(f"""\
        Usage: ./compare.py <compiler1> <compiler2> (-g / -v)
        no flag outputs data with no additional formatting,
        -g outputs a compact table of differences in error, and
        -v displays errors to a table in full.\
        """))
        sys.exit(1)

    master = {} # T R I P L E N E S T E D D I C T I O N A R Y M A N E U V E R
    master[lc.getline(f"../logs/{justTheName(sys.argv[1])}.log", 1).split()[2]] = buildDirectory(1)
    master[lc.getline(f"../logs/{justTheName(sys.argv[2])}.log", 1).split()[2]] = buildDirectory(2)

    rawOutput(master)

if __name__ == "__main__":
    main()