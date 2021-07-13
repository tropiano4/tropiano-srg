FORTRAN = gfortran-4.6
#FLAGS = -i-dynamic -r8 -O3 -w95 -
DEBUG = -Warray-bounds -Wall -Wunderflow -Warray-temporaries \
       -Wcharacter-truncation -Wtabs -Wintrinsic-shadow -Walign-commons -frange-check \
       -Wconversion -g
OPTIM = -ffree-form -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -O3
FLAGS = $(DEBUG)

all : main

input.mod : input.f90
	$(FORTRAN) $(FLAGS) -c $<

math.mod : math.f90 input.mod
	$(FORTRAN) $(FLAGS) -c $<

formfactor.mod : formfactor.f90 math.mod input.mod
	$(FORTRAN) $(FLAGS) -c $<

main.o : main.f90 math.mod input.mod formfactor.mod
	$(FORTRAN) $(FLAGS) -c $<

main : main.o math.o input.o formfactor.o math.o
	$(FORTRAN) $(FLAGS) -o $@ $^

clean :
	rm -f *.o *.mod *~ main
