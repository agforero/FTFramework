#!/usr/bin/env python3
import sys, os
import linecache as lc 
from textwrap import dedent

def justTheName(st): # helps generate things like <filename>.exe
    if len(st) == 0: return ''
    elif st[0] == '.': return ''
    else: return st[0] + justTheName(st[1:])

def splitAndLower(line):
    ls = line.split()
    for i in range(len(ls)): ls[i] = ls[i].lower()
    return ls

def traceLine(line): # is the line the result of the make --trace flag?
    return splitAndLower(line)[2:4] + splitAndLower(line)[5:7] == ['update', 'target', 'due', 'to:']

def enteredError(f_compiler, lines, i):
    return lines[i-1].split()[1] == f_compiler and lines[i].split()[1] != f_compiler and not traceLine(lines[i])

def errorDeclaration(line): # is the line make telling us that there's been an error
    # this function was made by referring to https://www.gnu.org/software/make/manual/html_node/Error-Messages.html
    banana = line.split()
    if len(banana) < 2: return False # we need at least a make: and something after
    if banana[0] == "#": banana = banana[1:] # shave off the pound sign if it's there
    if banana[0][:5] == "make:": # according to link above, this will always be true unless the Makefile itself is botched
        if banana[1] == "***": return True # should always be enough with my flags
        elif banana[-2] == "Error": 
            try: 
                int(banana[-1])
                return True # if not, see if "Error" is in the line
            except: pass # suppress the error message, even though this case is weird
        elif "Segfault" in banana: return True # maybe it segfaulted? we could go on with elifs at this point
    return False

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
        elif len(lines[i-1].split()) < 2 or len(lines[i].split()) < 2: continue # this should also cover whitespace
        elif inError and traceLine(lines[i]):
            inError = False
            if currentError in ret[currentDirectory].keys(): ret[currentDirectory][currentError] += errorBody
            else: ret[currentDirectory][currentError] = errorBody
        elif inError:
            errorBody += lines[i]
        elif enteredError(f_compiler, lines, i):
            inError = True
            currentError = lines[i-2].split()[4][1:-1]
            errorBody = lines[i-1] + lines[i]

    # at this point, the "error" isn't actually an error unless make declares that it is. thus:
    for dr in ret.keys():
        invalid = []
        for e in ret[dr].keys():
            valid = False
            for line in ret[dr][e].split("\n"):
                if errorDeclaration(line):
                    valid = True
            if not valid:
                invalid.append(e)
        for e in invalid:
            del ret[dr][e] # this should filter out anything that's purely a warning

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

class table:
    def left(st, mx):
        if len(st) > mx: st = st[:(mx - 3)] + "..."
        return f"{' ' * (mx - len(st))}{st} |"

    def leftLeft(st, mx): # yes I know it's a dumb name, leave me alone
        if len(st) > mx: st = st[:(mx - 3)] + "..."
        return f"{st}{' ' * (mx - len(st))} |"

    def right(st, mx):
        if len(st) > mx: st = st[:(mx - 3)] + "..."
        return f"{' ' * (mx)} | {st}"

    def center(st, mx, count1, count2): # expected length is 2 * mx + 3
        expected = (2 * mx) + 3
        totalHyphens = expected - len(st) - 4 - len(count1) - len(count2)
        ret = f"{count1} {'-' * int(totalHyphens / 2)} {st} {'-' * int(totalHyphens / 2)} {count2}"
        if expected - len(ret) != 0:
            ret = f"{count1} {'-' * (int(totalHyphens / 2) + 1)} {st} {'-' * int(totalHyphens / 2)} {count2}"
        return ret

    def empty(mx, ch=' '):
        return f"{ch * mx} | {ch * mx}"

    def bothSides(term1, term2, mx):
        if len(term1) > mx: term1 = term1[:(mx - 3)] + "..."
        if len(term2) > mx: term2 = term2[:(mx - 3)] + "..."
        return f"{' ' * (mx - len(term1))}{term1} | {term2}{' ' * (mx - len(term2))}"

    def bothSidesLeft(term1, term2, mx):
        if len(term1) > mx: term1 = term1[:(mx - 3)] + "..."
        if len(term2) > mx: term2 = term2[:(mx - 3)] + "..."
        return f"{term1}{' ' * (mx - len(term1))} | {term2}{' ' * (mx - len(term2))}"

# outputs [[], [], []], where 0 and 1 are list of keys that differ,
# and [2][0] is max target string length
def getDifferences(master, dr):
    ret = [[],[], [0]]
    compilers = list(master.keys())
    for e in list(master[compilers[0]][dr].keys()):
        if e not in list(master[compilers[1]][dr].keys()):
            ret[2][0] = max(ret[2][0], len(e), len(compilers[0]))
            ret[0].append(e)
    for e in list(master[compilers[1]][dr].keys()):
        if e not in list(master[compilers[0]][dr].keys()):
            ret[2][0] = max(ret[2][0], len(e), len(compilers[1]))
            ret[1].append(e)
    return ret

def rawOutput(master):
    compilers = list(master.keys())
    for dr in master[compilers[0]].keys():
        diffs = getDifferences(master, dr)
        if len(diffs[0]) > 0 or len(diffs[1]) > 0: 
            print(dedent(f"""\
            {'/' * 128}
            DIR: {dr}
            DIFF_COUNT: {compilers[0]}: {len(diffs[0])} {compilers[1]}: {len(diffs[1])}
            """))
        
        if len(diffs[0]) > 0:
            print(f"{'=' * 64}\nCOMPILER: {compilers[0]}")
            for e in diffs[0]:
                print(f"ERROR: {e}")
                print(master[compilers[0]][dr][e])
        
        if len(diffs[1]) > 0:
            print(f"{'=' * 64}\nCOMPILER: {compilers[1]}")
            for e in diffs[1]:
                print(f"ERROR: {e}")
                print(master[compilers[1]][dr][e])

def conciseTable(master):
    compilers = list(master.keys())
    diffs = []
    mx, i = 0, 0
    for dr in master[compilers[0]].keys():
        diffs.append(getDifferences(master, dr))
        mx = max(mx, diffs[-1][2][0])
    print(table.bothSides(compilers[0], compilers[1], mx))
    for dr in master[compilers[0]].keys():
        if len(diffs[i][0]) > 0 or len(diffs[i][1]) > 0: 
            print(f"{table.center(folderName(dr), mx, str(len(diffs[i][0])), str(len(diffs[i][1])))}")
        
        while len(diffs[i][0]) > 0 and len(diffs[i][1]) > 0:
            print(table.bothSides(diffs[i][0][0], diffs[i][1][0], mx))
            try: diffs[i][0], diffs[i][1] = diffs[i][0][1:], diffs[i][1][1:]
            except IndexError: break
        
        if len(diffs[i][0]) > 0:
            for e in diffs[i][0]:
                print(f"{table.left(e, mx)}")
        
        if len(diffs[i][1]) > 0:
            for e in diffs[i][1]:
                print(f"{table.right(e, mx)}")
        
        if len(diffs[i][0]) > 0 or len(diffs[i][1]) > 0: 
            print(table.empty(mx))
        
        i += 1

def verboseTable(master):
    compilers = list(master.keys())
    diffs = []
    mx = min(int(sys.argv[4]), 256)
    i = 0
    for dr in master[compilers[0]].keys():
        diffs.append(getDifferences(master, dr))
    print(table.bothSides(compilers[0], compilers[1], mx))
    for dr in master[compilers[0]].keys(): # let the
        bodies = [[], []] # hit the flooooooor
        
        for j in range(len(bodies)):
            for e in diffs[i][j]:
                bodies[j] += [e, ""] 
                for line in master[compilers[j]][dr][e].split("\n"):
                    bodies[j].append(line)
                bodies[j].append(f"{'-' * mx}")
        
        if len(bodies[0]) > 0 or len(bodies[1]) > 0: 
            print(f"{table.center(folderName(dr), mx, str(len(diffs[i][0])), str(len(diffs[i][1])))}")
        
        while len(bodies[0]) > 0 and len(bodies[1]) > 0:
            print(table.bothSidesLeft(bodies[0][0], bodies[1][0], mx))
            try: bodies[0], bodies[1] = bodies[0][1:], bodies[1][1:]
            except IndexError: break
        
        if len(bodies[0]) > 0:
            for line in bodies[0]:
                print(f"{table.leftLeft(line, mx)}")
        
        if len(bodies[1]) > 0:
            for line in bodies[1]:
                print(f"{table.right(line, mx)}")
        
        if len(bodies[0]) > 0 or len(bodies[1]) > 0: 
            print(table.empty(mx))
        
        i += 1

def printUsage():
    print(dedent(f"""\
    usage: ./compare.py <compiler1> <compiler2> (-g / -v <col>)
    no flag outputs data with no additional formatting,
    -g outputs a compact table of differences in error, and
    -v <col> displays full error messages to a table with <col> columns, up to 128.\
    """))

def main():
    version = os.popen("make -v").read().split("\n")[0]
    if int(version.split()[2][0]) < 4: # we might have a higher version some day!
        print(f"GNU Make 4.x required for compare.py; this shell is currently running {version}.")
        sys.exit(3)

    if len(sys.argv) >= 4:
        if sys.argv[3] != "-g" and sys.argv[3] != "-v":
            print(f"{sys.argv[3]}: invalid flag.")
            sys.exit(1)
        if sys.argv[3] == "-v":
            if len(sys.argv) != 5:
                print("incorrect number of arguments after -v.")
                sys.exit(4)
            try: int(sys.argv[4])
            except:
                print(f"{sys.argv[4]} cannot be converted to int.")
                sys.exit(5)
    elif len(sys.argv) == 1 or (len(sys.argv) >= 2 and (sys.argv[1] == "-h" or sys.argv[1] == "--help")):
        printUsage()
        sys.exit(0)
    elif len(sys.argv) != 3:
        printUsage()
        sys.exit(1)

    master = {} # T R I P L E N E S T E D D I C T I O N A R Y M A N E U V E R
    master[lc.getline(f"../logs/{justTheName(sys.argv[1])}.log", 1).split()[2]] = buildDirectory(1)
    master[lc.getline(f"../logs/{justTheName(sys.argv[2])}.log", 1).split()[2]] = buildDirectory(2)

    if len(sys.argv) == 3: rawOutput(master)
    elif len(sys.argv) >= 4:
        if sys.argv[3] == "-g": conciseTable(master)
        if sys.argv[3] == "-v": verboseTable(master)

if __name__ == "__main__":
    main()