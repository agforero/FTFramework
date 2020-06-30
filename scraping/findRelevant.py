#!/usr/bin/env python
import glob, os, sys

def isIn(e, ls):
    for i in ls:
        if i == e:
            return True
    return False

def check(ls): # checks if any file in l_files is supposed to be excluded
    ret = []
    for i in ls:
        if not isIn(i, sys.argv[2:]): # if the file in question isn't supposed to be excluded,
            ret.append(i) # add it to the return list
    return ret

def main():
    os.chdir("./")
    disp = ""
    l_files = check(glob.glob("*"+sys.argv[1]))
    for file in l_files:
        disp = disp + file[:0-len(sys.argv[1])] + " "
    print(disp)

if __name__ == "__main__":
    main()