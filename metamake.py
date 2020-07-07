#!/usr/bin/env python
import os, sys, glob

def find(l_fc, f): # dumb function to make up for the fact that I didn't make l_fc a dictionary
	for i in range(len(l_fc)):
		if l_fc[i].name == f:
			return i
	return -1 # should never get here unless error

class rDepList():
	def __init__(self):
		self.ls = []

	def fix(self, head): # where head is the original file. we don't wanna read it twice
		self.ls.remove(head) # should snip out the parent file

	def printAll(self):
		for e in self.ls:
			print(e)

class FFile:
	def __init__(self, nm):
		self.name = nm
		self.importantLines = [] # should cut down running time...if not, delete. should be a list of lists
		self.independent = True # independent to start with
		self.isProgram = False # is an actual program rather than just some functions
		self.dependencies = {} # no dependencies to start with
		self.explored = False

	def printInfo(self):
		print(f"\nname: {self.name}")
		print(f"importantLines: {self.importantLines}")
		print(f"independent: {self.independent}")
		print(f"dependencies: {self.dependencies}")

	def printEssentials(self):
		print(f"\nname: {self.name}")
		print(f"dependencies: {self.dependencies}")
	
	def depOfDep(self, l_fc): # DFS BABYYYYY
		r_ls = rDepList()
		self.explore(l_fc, r_ls)
		r_ls.fix(self.name)
		for fc in l_fc: fc.explored = False
		return r_ls.ls

	def explore(self, l_fc, r_ls): # in two acts
		self.explored = True
		for key in self.dependencies:
			if not l_fc[find(l_fc, self.dependencies[key])].explored:
				l_fc[find(l_fc, self.dependencies[key])].explore(l_fc, r_ls)
		r_ls.ls.append(self.name)

	def outputInfo(self): # it's like a beaten up pickup truck, held together by duct tape and glue...but it works.
		if self.isProgram: ret = justTheName(self.name) + ".exe: " + self.name # <filename>.exe: <filename>.f90
		else: ret = justTheName(self.name) + ".o: " + self.name # <filename>.o: <filename>.f90
		
		for key in self.dependencies.keys(): # add needed .o files if necessary
			ret += " " + justTheName(self.dependencies[key]) + ".o"
		
		ret += "\n\t$(FC) $(FFLAGS) -c " + self.name # $(FC) $(FFLAGS) -c <filename>.f90
		if self.isProgram: # $(FC) $(FFLAGS) <filename>.o (depedencies).o -o <filename>.exe
			ret += "\n\t$(FC) $(FFLAGS) " + justTheName(self.name) + ".o"
			for key in self.dependencies.keys():
				ret += " " + justTheName(self.dependencies[key]) + ".o"
			ret += " -o " + justTheName(self.name) + ".exe"

		ret += "\n"
		return ret

def justTheName(st): # helps generate things like <filename>.exe
	if len(st) == 0: return ''
	elif st[0] == '.': return ''
	else: return st[0] + justTheName(st[1:])

def allHelper(l_fc):
	ret = ""
	for fc in l_fc:
		if fc.isProgram: ret += " " + justTheName(fc.name) + ".exe" 
		else: ret += " " + justTheName(fc.name) + ".o"
	return ret

def main():
	os.chdir("./") # necessary?
	l_f90 = glob.glob("*.f90") # so far this only works with .f90 files
	l_fc = [] # list containing FFile objects
	allDependencies = {} # dictionary for easy access to dependencies later

	# first, let's find all files that are needed elsewhere (aka needed in allDependencies)
	for file in l_f90: 
		cur = FFile(file)
		f = open(file, 'r')
		for line in f:
			try: 
				# the file relies on another file
				if line.split()[0].lower() == "use":
					cur.independent = False 
					cur.importantLines.append(line.split()) 
					cur.dependencies[line.split()[1]] = "" # we don't know yet where to find this module
					if line.split()[1] not in allDependencies: allDependencies[line.split()[1]] = "" 

				# the file is relied on by another file
				elif line.split()[0].lower() == "module":
					cur.importantLines.append(line.split())
					if line.split()[1] not in allDependencies: allDependencies[line.split()[1]] = "" # we don't yet try to tag the dependencies

				elif line.split()[0].lower() == "program":
					cur.isProgram = True

			except: continue # if the line is empty, or error
		l_fc.append(cur)
		f.close()

	# by now, we should have all dependencies ever needed
	# next, we locate all the filenames containing what we've put into allDependencies
	for fc in l_fc: # we gotta scan through all files 
		for line in fc.importantLines:
			if line[0].lower() == "module": allDependencies[line[1]] = fc.name # this means that the file in question has the module we're looking for

	# at this point, allDependencies should be filled with every module we'll ever need, and where to find them.
	# now, we need to link the .dependencies attribute of FFile.
	for fc in l_fc:
		for key in fc.dependencies.keys():
			fc.dependencies[key] = allDependencies[key] # matches all keys in the current object with its counterpart in allDepdendencies
		for e in fc.depOfDep(l_fc):
			if e not in fc.dependencies.values():
				fc.dependencies["rec_" + justTheName(e)] = e

	# beautiful! it's working smoothly. now to translate the result into a Makefile...
	m = open("Makefile", 'w')
	m.write("FC ?= gfortran\n")
	m.write("FFLAGS ?=\n")
	m.write(f"\nall: {allHelper(l_fc)}\n\n")
	for fc in l_fc:
		m.write(fc.outputInfo()+"\n")
	m.write(f"clean:\n\trm -f *.exe *.o *.mod")
	m.close()

	return 0

if __name__ == "__main__":
	main()