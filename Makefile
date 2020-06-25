
FC = my
FFLAGS = name

all: hello.exe

hello.exe: hello.is 
	$(FC) $(FFLAGS) -c hello.is
	$(FC) $(FFLAGS) hello.o -o hello.exe

clean:
	rm -f hello.exe *.o
