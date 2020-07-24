#!/usr/bin/env python
import sys, os

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

def makeEven(mx):
    if mx % 2 != 0: return mx + 1
    else: return mx

def evenOut(st, mx, cha = ' '): # help make all the strings uniform size
    if len(st) > mx-1: return st[:mx-3] + "..."
    else: return st + " " + (cha * (mx - len(st) - 1))

def evenOutRight(st, mx): # nice in theory, looks a bit messy in application
    return (' ' * (mx - len(st))) + st

def center(st, mx):
    return ('-' * (mx - int(len(st)/2))) + ' ' + st + ' ' + ('-' * (mx - int(len(st)/2)))

def matchFirstLine(dc, first): # returns index location of Error object that matches given first line
    for key in dc.keys():
        for i in range(len(dc[key])):
            if dc[key][i].getFirstLine() == first:
                return i
    return -1 # no such first line was found

def justTheName(st): # helps generate things like <filename>.exe
    if len(st) == 0: return ''
    elif st[0] == '.': return ''
    else: return st[0] + justTheName(st[1:])

def displayBasicData(f1e, f2e, mx):
    for key in f1e.keys():
        if (len(f1e[key]), len(f2e[key])) == (0, 0): continue
        cont1, cont2 = False, False
        for e in f1e[key]: # for error in current directory of f1e
            if matchFirstLine(f2e, e.getFirstLine()) == -1:
                cont1 = True
        for e in f2e[key]: # for error in current directory of f2e
            if matchFirstLine(f1e, e.getFirstLine()) == -1:
                cont1 = True
        if cont1 or cont2:
            print(f"{('-' * (mx + 4))}")
            print(f"DIR: {key[:-2]}")
            if cont1: print(f"COMPILER: {sys.argv[1]}")
            for e in f1e[key]: # for error in current directory of f1e
                if matchFirstLine(f2e, e.getFirstLine()) == -1:
                    print(e.getFirstLine()[:-1])
                    for line in e.body.split("\n")[1:]:
                        print(line)
            if cont2: print(f"COMPILER: {sys.argv[2]}")
            for e in f2e[key]: # for error in current directory of f2e
                if matchFirstLine(f1e, e.getFirstLine()) == -1:
                    print(e.getFirstLine()[:-1])
                    for line in e.body.split("\n")[1:]:
                        print(line)

def displaySummary(f1e, f2e, mx):
    print(f"\n{' ' * (mx-len(sys.argv[1]))}{sys.argv[1]} | {sys.argv[2]}")
    print(f"{'-' * mx} | {'-' * mx}")
    for key in f1e.keys():
        if (len(f1e[key]), len(f2e[key])) == (0, 0): continue
        cont = False # determine if the directory actually needs to be displayed or not
        for e in f1e[key]: # for error in current directory of f1e
            if matchFirstLine(f2e, e.getFirstLine()) == -1:
                cont = True
        for e in f2e[key]: # for error in current directory of f2e
            if matchFirstLine(f1e, e.getFirstLine()) == -1:
                cont = True
        if cont:
            print(f"{' ' * mx} |\n{center(key[:-1], mx)}")
            for e in f1e[key]: # for error in current directory of f1e
                if matchFirstLine(f2e, e.getFirstLine()) == -1:
                    print(f"{evenOut(e.getFirstLine()[:-1], mx)} |")
            for e in f2e[key]: # for error in current directory of f2e
                if matchFirstLine(f1e, e.getFirstLine()) == -1:
                    print(f"{' ' * mx} | {evenOut(e.getFirstLine()[:-1], mx)}")
            print(f"{' ' * mx} |")

def displayVerboseSummary(f1e, f2e, mx):
    print(f"\n{' ' * (mx-len(sys.argv[1]))}{sys.argv[1]} | {sys.argv[2]}")
    print(f"{'-' * mx} | {'-' * mx}")
    for key in f1e.keys():
        if (len(f1e[key]), len(f2e[key])) == (0, 0): continue
        cont = False
        for e in f1e[key]: # for error in current directory of f1e
            if matchFirstLine(f2e, e.getFirstLine()) == -1:
                cont = True
        for e in f2e[key]: # for error in current directory of f2e
            if matchFirstLine(f1e, e.getFirstLine()) == -1:
                cont = True
        if cont:
            print(f"{' ' * mx} |\n{center(key[:-1], mx)}")
            for e in f1e[key]: # for error in current directory of f1e
                if matchFirstLine(f2e, e.getFirstLine()) == -1:
                    print(f"{'-' * mx} |")
                    print(f"{evenOut(e.getFirstLine()[:-1] + ':', mx)} |")
                    print(f"{' ' * mx} |")
                    for line in e.body.split("\n")[1:]:
                        print(f"{evenOut(line, mx)} |")
            for e in f2e[key]: # for error in current directory of f2e
                if matchFirstLine(f1e, e.getFirstLine()) == -1:
                    print(f"{' ' * mx} | {'-' * mx}")
                    print(f"{' ' * mx} | {evenOut(e.getFirstLine()[:-1] + ':', mx)}")
                    print(f"{' ' * mx} |")
                    for line in e.body.split("\n")[1:]:
                        print(f"{' ' * mx} | {evenOut(line, mx)}")
            print(f"{' ' * mx} |")

def main():
    # check if args are ok
    if len(sys.argv) == 1 or len(sys.argv) > 4 or (len(sys.argv) == 4 and not (sys.argv[3] == "-v" or sys.argv[3] == "-g")):
        print("Usage: ./compare.py <log1> <log2> (-g / -v)")
        print("no flag outputs to data with no extra formatting,")
        print("-g outputs to a compact table, and")
        print("-v outputs to a verbose table.")
        sys.exit(1)

    os.chdir("../")
    f1 = open(f"logs/{justTheName(sys.argv[1])}.log", 'r')
    f2 = open(f"logs/{justTheName(sys.argv[2])}.log", 'r')
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
                maxRelLen = min(60, max(len(line), maxRelLen))
                currentDir = line
                f1Errors[currentDir] = []
            
            elif splitAndLower(line)[:2] == ["not", "ok"]:
                inError = True
                maxRelLen = min(60, max(len(line), maxRelLen))
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
                maxRelLen = min(60, max(len(line), maxRelLen))
                currentDir = line
                f2Errors[currentDir] = []

            elif splitAndLower(line)[:2] == ["not", "ok"]:
                inError = True
                maxRelLen = min(60, max(len(line), maxRelLen))
                f2Errors[currentDir].append(Error(line, currentDir))

            elif inError and line.split()[0] != "ok":
                f2Errors[currentDir][-1].addLine(line)
        except: continue
    
    f1.close()
    f2.close()

    maxRelLen = makeEven(maxRelLen)
    if len(sys.argv) == 3:
        if getDifferences(reduceToFirstLines(f1Errors), reduceToFirstLines(f2Errors)) != {}:
            displayBasicData(f1Errors, f2Errors, maxRelLen)
        else: print("No differences found. Nice compiler!")
    elif len(sys.argv) == 4:
        if sys.argv[3] == "-g":
            if getDifferences(reduceToFirstLines(f1Errors), reduceToFirstLines(f2Errors)) != {}:
                displaySummary(f1Errors, f2Errors, maxRelLen)
            else: print("No differences found. Nice compiler!")
        elif sys.argv[3] == "-v":
            if getDifferences(reduceToFirstLines(f1Errors), reduceToFirstLines(f2Errors)) != {}:
                displayVerboseSummary(f1Errors, f2Errors, maxRelLen)
            else: print("No differences found. Nice compiler!")

    sys.exit(0)

if __name__ == "__main__":
    main()