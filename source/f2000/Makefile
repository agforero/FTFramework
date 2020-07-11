FC ?= gfortran
FFLAGS ?=

.PHONY: all clean
all:  ctrper.o fndnth.o follow.exe givcor.exe indmed.o indnth.o inspar.o inssor.o median.o mrgref.o mrgrnk.o mulcnt.o rapknr.o refpar.o refsor.o rinpar.o rnkpar.o sort7.exe tstvalnth.exe uniinv.o unipar.o unirnk.o unista.o valmed.o valnth.o

ctrper.o: ctrper.f90 mrgrnk.o
	$(FC) $(FFLAGS) -c ctrper.f90

fndnth.o: fndnth.f90
	$(FC) $(FFLAGS) -c fndnth.f90

follow.exe: follow.o
	$(FC) $(FFLAGS) follow.o -o follow.exe

follow.o: follow.f90
	$(FC) $(FFLAGS) -c follow.f90

givcor.exe: givcor.o ctrper.o refsor.o mrgrnk.o
	$(FC) $(FFLAGS) givcor.o ctrper.o refsor.o mrgrnk.o -o givcor.exe

givcor.o: givcor.f90 ctrper.o refsor.o mrgrnk.o
	$(FC) $(FFLAGS) -c givcor.f90

indmed.o: indmed.f90
	$(FC) $(FFLAGS) -c indmed.f90

indnth.o: indnth.f90
	$(FC) $(FFLAGS) -c indnth.f90

inspar.o: inspar.f90
	$(FC) $(FFLAGS) -c inspar.f90

inssor.o: inssor.f90
	$(FC) $(FFLAGS) -c inssor.f90

median.o: median.f90
	$(FC) $(FFLAGS) -c median.f90

mrgref.o: mrgref.f90
	$(FC) $(FFLAGS) -c mrgref.f90

mrgrnk.o: mrgrnk.f90
	$(FC) $(FFLAGS) -c mrgrnk.f90

mulcnt.o: mulcnt.f90 uniinv.o
	$(FC) $(FFLAGS) -c mulcnt.f90

rapknr.o: rapknr.f90
	$(FC) $(FFLAGS) -c rapknr.f90

refpar.o: refpar.f90
	$(FC) $(FFLAGS) -c refpar.f90

refsor.o: refsor.f90
	$(FC) $(FFLAGS) -c refsor.f90

rinpar.o: rinpar.f90
	$(FC) $(FFLAGS) -c rinpar.f90

rnkpar.o: rnkpar.f90
	$(FC) $(FFLAGS) -c rnkpar.f90

sort7.exe: sort7.o
	$(FC) $(FFLAGS) sort7.o -o sort7.exe

sort7.o: sort7.f90
	$(FC) $(FFLAGS) -c sort7.f90

tstvalnth.exe: tstvalnth.o valnth.o mrgrnk.o
	$(FC) $(FFLAGS) tstvalnth.o valnth.o mrgrnk.o -o tstvalnth.exe

tstvalnth.o: tstvalnth.f90 valnth.o mrgrnk.o
	$(FC) $(FFLAGS) -c tstvalnth.f90

uniinv.o: uniinv.f90
	$(FC) $(FFLAGS) -c uniinv.f90

unipar.o: unipar.f90
	$(FC) $(FFLAGS) -c unipar.f90

unirnk.o: unirnk.f90
	$(FC) $(FFLAGS) -c unirnk.f90

unista.o: unista.f90 uniinv.o
	$(FC) $(FFLAGS) -c unista.f90

valmed.o: valmed.f90
	$(FC) $(FFLAGS) -c valmed.f90

valnth.o: valnth.f90
	$(FC) $(FFLAGS) -c valnth.f90

clean:
	rm -f *.exe *.o *.mod