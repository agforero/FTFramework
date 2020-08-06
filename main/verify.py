#!/usr/bin/env python3

import os, glob, sys
from textwrap import dedent

def findl_files(dr, orig):
    os.chdir(dr)
    files = []
    for ext in ("*.f90", "*.f", "*.f95", "*.f03", "*.f08", "*.for", "*.f77", "*.ftn"):
        if len(glob.glob(ext)) != 0: files += glob.glob(ext)
        elif len(glob.glob(ext.upper())) != 0: files += glob.glob(ext.upper())
    os.chdir(orig)
    return files

def justTheName(st): # helps generate things like <filename>.exe
    if len(st) == 0: return ''
    elif st[0] == '.': return ''
    else: return st[0] + justTheName(st[1:])

def splitAndLower(line):
    ls = line.split()
    for i in range(len(ls)): ls[i] = ls[i].lower()
    return ls

def folderName(path):
    for i in range(len(path)):
        if path[i:i + len("/source")] == "/source":
            return path[i:]
    return path

def main():
    os.chdir("../")
    origDir = f"{os.getcwd()}/source"
    fileList = {}
    first = True
    for dr in os.walk(origDir):
        fileList[dr[0]] = findl_files(dr[0], origDir)
    for dr in fileList.keys():
        allMods = {}
        output = []
        uniErrorFound = False
        for ffile in fileList[dr]:
            f = open(f"{dr}/{ffile}", 'r')
            inInt = 0
            try:
                for line in f:
                    if line.split()[0].lower() == "interface": inInt += 1 
                    elif splitAndLower(line)[:2] == ["end", "interface"]: inInt -= 1 

                    elif line.split()[0].lower() == "module": 
                        if line.split()[1] not in list(allMods.keys()) and not inInt:  
                            allMods[line.split()[1]] = [] 
                            allMods[line.split()[1]].append(ffile) 
                        elif line.split()[1] in list(allMods.keys()) and allMods[line.split()[1]][0] == ffile:  
                            continue # skip it, we're good 
                        elif not inInt: 
                            allMods[line.split()[1]].append(ffile) 
            except UnicodeDecodeError: 
                output.append(f"ERROR: unreadable char in {ffile}; skipping.") 
                uniErrorFound = True
            except IndexError: continue 
            f.close() 
        if uniErrorFound: output[-1] += "\n"

        for mod in allMods.keys():
            if len(allMods[mod]) > 1:
                output.append(f"RACE CONDITION: \"{mod}.mod\" is outputted by:")
                for f in allMods[mod]:
                    output[-1] += f"\n- {f}"
                output[-1] += "\n"

        if len(output) > 0:
            if first:
                first = False
                print(f"{'-' * 64}\nDIR: {folderName(dr)}/")
            else: print(f"{'-' * 64}\nDIR: {folderName(dr)}/") 
            for line in output:
                print(line)

if __name__ == "__main__":
    main()