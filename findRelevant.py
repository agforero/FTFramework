#!/usr/bin/env python
import glob, os, sys

def printUsage():
    print("Usage:\t./findRelevant.py -flag <extension> (exclusions)")

def isIn(e, ls):
    for i in ls:
        if i == e:
            return True
    return False

def check(ls): # checks if any file in l_files is supposed to be excluded
    ret = []
    for i in ls:
        if not isIn(i, sys.argv[3:]): # if the file in question isn't supposed to be excluded,
            ret.append(i) # add it to the return list
    return ret

def main():
    os.chdir("./")

    if len(sys.argv) == 1: # if they used no arguments
        printUsage()
        sys.exit(1)

    elif sys.argv[1][0] != '-': # they didn't use a flag
        printUsage()
        sys.exit(2)

    elif sys.argv[1] == "-help" or sys.argv[1] == "-h":
        printUsage()
        print("Where:\t-help or -h displays options,")
        print("\t-f returns only filenames with no extension,")
        print("\t-o returns filenames with \'.o\' appended,")
        print("\t-x returns filenames with \'.exe\' appended, and ")
        print("\t-e returns filenames with extensions kept.")
        print("\t<extension> designates the extension to search for, and\n\t(exclusions) are files to be excluded from output.")
        sys.exit(0)

    disp = ""
    l_files = check(glob.glob("*"+sys.argv[2]))

    for file in l_files:
        if sys.argv[1] == "-f": # just the filename
            disp = disp + file[:0-len(sys.argv[2])] + " "

        elif sys.argv[1] == "-o": # with .o added
            disp = disp + file[:0-len(sys.argv[2])] + ".o "

        elif sys.argv[1] == "-x": # with .exe appended
            disp = disp + file[:0-len(sys.argv[2])] + ".exe "

        elif sys.argv[1] == "-e": # with extension maintained
            disp = disp + file + " "

        else: # invalid argument
            printUsage()
            sys.exit(1)

    print(disp)
    if disp == "":
        sys.exit(1)

if __name__ == "__main__":
    main()