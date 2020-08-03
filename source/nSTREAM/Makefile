CC = gcc
CFLAGS = -O2 -fopenmp

FC = gfortran
FFLAGS = -O2 -fopenmp

all: nstream.exe

nstream.exe: nstream.F90 mysecond.o
	$(FC) $(FFLAGS) -c nstream.F90
	$(FC) $(FFLAGS) nstream.o mysecond.o -o nstream.exe

mysecond.o: mysecond.c
	$(CC) $(CFLAGS) -c mysecond.c

clean:
	rm -f nstream.exe *.o