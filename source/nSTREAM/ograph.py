#!/usr/bin/env python

import matplotlib.pyplot as plt

def fix(ls):
    if len(ls) == 2:
        ls[1] = ls[1].strip()
        return int(ls[1])
    ret = []
    for i in range(len(ls)):
        if ls[i] != '':
            ls[i] = ls[i].strip()
            if ls[i][-2:] == "\n":
                ls[i] = ls[i][:-2]
            try:
                ret.append(float(ls[i]))
            except:
                ret.append(0.0) # if fortran outputs something w/ asterisks
    return ret

def checkMax(start, check):
    if start > check: # if the checked number isn't the new maximum
        return start # return the old maximum
    else:   
        return check # otherwise, return the new maximum. the new KING

def main():
    f = open("output.txt", 'r')
    count = 0
    nSize = []
    
    cRate = [] # Copy
    sRate = [] # Scale
    aRate = [] # Add
    tRate = [] # Triad
    
    ticks = []
    subticks = [[],[],[],[]] # 0 = copy, 1 = scale, 2 = add and 3 = triad
    submax = [0,0,0,0] # keeps track of maximum value per operation

    for i in f:
        ls = i.split("      ")
        if count % 5 == 0: 
            ticks.append(fix(ls))
            try:
                nSize.append(fix(ls)[0])
            except:
                nSize.append(fix(ls))

        elif count % 5 == 1: # Copy
            cRate.append(fix(ls)[0])
            submax[0] = checkMax(submax[0], fix(ls)[0])

        elif count % 5 == 2: # Scale
            sRate.append(fix(ls)[0])
            submax[1] = checkMax(submax[1], fix(ls)[0])

        elif count % 5 == 3: # Add
            aRate.append(fix(ls)[0])
            submax[2] = checkMax(submax[2], fix(ls)[0])

        elif count % 5 == 4: # Triad
            tRate.append(fix(ls)[0])
            submax[3] = checkMax(submax[3], fix(ls)[0])
            
        count += 1

    f.close()

    # now, generate the y ticks
    # let's add cushioning to make the values prettier.
    for i in range(len(submax)):
        submax[i] = int(submax[i] - (submax[i] % 1000) + 1000) # rounds UP to nearest 1000

    # then divide y ticks evenly according to maxval
    for i in range(4):
        for j in range(1,5):
            subticks[i].append(int(((j)/4) * submax[i]))

    figure, axes = plt.subplots(nrows=4, ncols=1)

    # Copy
    axes[0].plot(nSize, cRate, 'r')
    axes[0].set_xticks([])
    axes[0].set_ylabel("Copy")
    axes[0].set_yticks(subticks[0])

    # Scale
    axes[1].plot(nSize, sRate, 'g')
    axes[1].set_xticks([])
    axes[1].set_ylabel("Scale")
    axes[1].set_yticks(subticks[1])

    # Add
    axes[2].plot(nSize, aRate, 'b')
    axes[2].set_xticks([])
    axes[2].set_ylabel("Add")
    axes[2].set_yticks(subticks[2])

    # Triad
    axes[3].plot(nSize, tRate, 'm')
    axes[3].set_xlabel("Memory Footprint (MB)")
    axes[3].set_ylabel("Triad")
    axes[3].set_yticks(subticks[3])
  
    figure.suptitle("Rate (MB/s)")
    plt.savefig("graph.png", dpi=960, bbox_inches="tight", pad_inches=0.5)

if __name__ == "__main__":
    main()

