FORTRAN=ifort
FLAGS=-i-static -r8 -C -w95 -implicitnone

all : extract

symd3t.mod : symd3t.f90 
	$(FORTRAN) $(FLAGS) -c $<

symvbb.mod : symvbb.f90 symd3t.mod
	$(FORTRAN) $(FLAGS) -c $<

symvc.mod : symvc.f90 symd3t.mod
	$(FORTRAN) $(FLAGS) -c $<

extract.o : extract.f90 symd3t.mod symvbb.mod symvc.mod
	$(FORTRAN) $(FLAGS) -c $<

extract : extract.o symd3t.o symvbb.o symvc.o
	$(FORTRAN) $(FLAGS) -o $@ $^

clean :
	rm -f *.o *.mod *~ extract
