#!/usr/bin/env python
import sys

class Error():
    def __init__(self, startLine, parent):
        self.body = startLine
        self.firstLine = startLine
        self.dir = parent
    
    def addLine(self, line):
        self.body += line

    def getFirstLine(self):
        return self.firstLine
        
def splitAndLower(line):
    ls = line.split()
    for i in range(len(ls)): ls[i] = ls[i].lower()
    return ls

def printNicely(dc): # only the "not ok" part
    for key in dc.keys():
        print(repr(key))
        for e in dc[key]:
            print(e.getFirstLine()[:-1])
        print()

def printEverything(dc): # DUMP TRUCK THE WHOLE ERROR MESSAGE BABYYYYY
    for key in dc.keys():
        print(f"KEY:\n{repr(key)}\n\nERRORS: {len(dc[key])}\n\nBODIES:")
        for e in dc[key]:
            print(e.body + "\n\n")
        print()

def reduceToFirstLines(dc):
    ret = {}
    for key in dc.keys():
        ret[key] = []
        for e in dc[key]:
            ret[key].append(e.getFirstLine())
    return ret

def getDifferences(f1e, f2e): # this should take in only the first lines...in a list, in a dictionary
    ret = {}
    for key in f1e.keys():
        ret[key] = []
        for line in f1e[key]:
            if line not in f2e[key]:
                ret[key].append(f"{sys.argv[1]}: {line}")
    for key in f2e.keys():
        for line in f2e[key]:
            if line not in f1e[key]:
                ret[key].append(f"{sys.argv[2]}: {line}")
    finalRet = {}
    for key in ret.keys():
        if len(ret[key]) != 0:
            finalRet[key] = ret[key]
    return finalRet

def evenOut(st, mx): # help make all the strings uniform size
    for i in range(mx - len(st)):
        st += ' '
    return st

def displaySummary(f1e, f2e, mx):
    print(f"\n{' ' * (mx-len(sys.argv[1]))}{sys.argv[1]} | {sys.argv[2]}")
    print(f"{'-' * mx} | {'-' * mx}")
    for key in f1e.keys():
        if (len(f1e[key]), len(f2e[key])) == (0, 0): continue
        print(f"{evenOut(key[:-1], mx)} | {evenOut(key[:-1], mx)}")
        if len(f1e[key]) >= len(f2e[key]):
            i2 = 0
            for i in range(len(f2e[key])):
                print(f"{evenOut(f1e[key][i][:-1], mx)} | {evenOut(f2e[key][i][:-1], mx)}")
                i2 += 1
            for i in range(len(f1e[key]) - len(f2e[key])):
                print(f"{evenOut(f1e[key][i2][:-1], mx)} |")
                i2 += 1
        else:
            i2 = 0
            for i in range(len(f1e[key])):
                print(f"{evenOut(f1e[key][i][:-1], mx)} | {evenOut(f2e[key][i][:-1], mx)}")
                i2 += 1
            for i in range(len(f2e[key]) - len(f1e[key])):
                print(f"{' ' * mx} | {evenOut(f2e[key][i2][:-1], mx)}")
                i2 += 1

        print(f"{' ' * mx} |")

def main():
    # check if args are ok
    if len(sys.argv) != 3:
        print("Usage: ./compare.py <log1> <log2>")
        sys.exit(1)

    f1 = open(sys.argv[1], 'r')
    f2 = open(sys.argv[2], 'r')
    maxRelLen = 0 # maximum relevant line length
    f1Errors = {}
    f2Errors = {}
    inError = False
    currentDir = ""
    for line in f1:
        try:
            if inError and line.split()[0] == "ok":
                inError = False
            
            elif ((line[0], line[-2]) == ('/', ':')) or ((line[0], line[-1]) == ('/', ':')):
                inError = False
                maxRelLen = max(len(line), maxRelLen)
                currentDir = line
                f1Errors[currentDir] = []
            
            elif splitAndLower(line)[:2] == ["not", "ok"]:
                inError = True
                maxRelLen = max(len(line), maxRelLen)
                f1Errors[currentDir].append(Error(line, currentDir))

            elif inError and line.split()[0] != "ok":
                f1Errors[currentDir][-1].addLine(line)
        
        except: continue

    for line in f2:
        try:
            if inError and line.split()[0] == "ok":
                inError = False

            elif ((line[0], line[-2]) == ('/', ':')) or ((line[0], line[-1]) == ('/', ':')):
                inError = False
                maxRelLen = max(len(line), maxRelLen)
                currentDir = line
                f2Errors[currentDir] = []

            elif splitAndLower(line)[:2] == ["not", "ok"]:
                inError = True
                maxRelLen = max(len(line), maxRelLen)
                f2Errors[currentDir].append(Error(line, currentDir))

            elif inError and line.split()[0] != "ok":
                f2Errors[currentDir][-1].addLine(line)

        except: continue
    
    f1.close()
    f2.close()

    if getDifferences(reduceToFirstLines(f1Errors), reduceToFirstLines(f2Errors)) != {}:
        displaySummary(reduceToFirstLines(f1Errors), reduceToFirstLines(f2Errors), maxRelLen)
    else: print("No differences found. Nice compiler!")

    sys.exit(0)

if __name__ == "__main__":
    main()