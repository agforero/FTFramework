#!/usr/bin/env python
import sys

def main():
	if len(sys.argv) < 4 or len(sys.argv) > 5:
		print("metamake: expected 4-5 env args, recieved %d." % len(sys.argv))
		return 1

	f = open("Makefile", 'w')
	name = sys.argv[1]
	comp = sys.argv[2]
	if len(sys.argv) == 5:
		flags = "\nFLAGS = " + sys.argv[3]
		ext = sys.argv[4]
	else:
		flags = ""
		ext = sys.argv[3]

	f.write(

"""FC = %s%s

all: %s.exe

%s.exe: %s%s 
	$(FC) $(FFLAGS) -c %s%s
	$(FC) $(FFLAGS) %s.o -o %s.exe

clean:
	rm -f %s.exe *.o""" 

	% (comp, flags, name, name, name, ext, name, ext, name, name, name)
  	)

	f.close()
	return 0

if __name__ == "__main__":
	main()