# automatically generated .bats tests to check the compilation of each program in this directory

@test make follow.exe {
	rm -f follow.o follow.exe follow.mod
	make follow.exe
}

@test make givcor.exe {
	rm -f givcor.o givcor.exe givcor.mod
	make givcor.exe
}

@test make sort7.exe {
	rm -f sort7.o sort7.exe sort7.mod
	make sort7.exe
}

@test make tstvalnth.exe {
	rm -f tstvalnth.o tstvalnth.exe tstvalnth.mod
	make tstvalnth.exe
}

