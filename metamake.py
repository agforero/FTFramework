#!/usr/bin/env python
import sys

def main():
	f = open("Makefile", 'w')
	f.write("""FC ?= gfortran
FFLAGS ?= 

SRC = $(shell ./findRelevant.py -e .f90)
EXE = $(shell ./findRelevant.py -x .f90)
OBJS = $(shell ./findRelevant.py -o .f90)

all: $(EXE)

%.exe: %.o
	$(FC) $(FFLAGS) $< -o $@

%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

clean: 
	rm -f $(EXE) $(OBJS)""")

	f.close()
	return 0

if __name__ == "__main__":
	main()