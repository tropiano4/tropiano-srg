FORTRAN=ifort
FLAGS=-i-static -r8 -O3 -w95 -implicitnone

all : main_PAV

input.mod : input.f90
	$(FORTRAN) $(FLAGS) -c $<

lgfactor.mod : lgfactor.f90
	$(FORTRAN) $(FLAGS) -c $<

constant.mod : constant.f90
	$(FORTRAN) $(FLAGS) -c $<

complex_eigen.mod : complex_eigen.f90 constant.mod
	$(FORTRAN) $(FLAGS) -c $<

deter.mod : deter.f90
	$(FORTRAN) $(FLAGS) -c $<

indexx.mod : indexx.f90
	$(FORTRAN) $(FLAGS) -c $<

math.mod : math.f90 input.mod lgfactor.mod complex_eigen.mod indexx.mod
	$(FORTRAN) $(FLAGS) -c $<

symd3t_proj.mod : symd3t_proj.f90 input.mod
	$(FORTRAN) $(FLAGS) -c $<

symtalm.mod : symtalm.f90 input.mod lgfactor.mod math.mod
	$(FORTRAN) $(FLAGS) -c $<

gauss.mod : gauss.f90
	$(FORTRAN) $(FLAGS) -c $<

angmom.mod : angmom.f90 input.mod lgfactor.mod math.mod
	$(FORTRAN) $(FLAGS) -c $<

bessik.mod : bessik.f90
	$(FORTRAN) $(FLAGS) -c $<

global.mod : global.f90 bessik.mod input.mod lgfactor.mod symd3t_proj.mod math.mod \
	symtalm.mod gauss.mod angmom.mod
	$(FORTRAN) $(FLAGS) -c $<


GLOBAL_OBJECTS = input.o lgfactor.o indexx.o symd3t_proj.o math.o deter.o complex_eigen.o constant.o symtalm.o gauss.o \
	angmom.o global.o bessik.o

r1r1.mod : r1r1.f90
	$(FORTRAN) $(FLAGS) -c $<

nucleus.mod : nucleus.f90 math.mod r1r1.mod
	$(FORTRAN) $(FLAGS) -c $<

symfield_proj.mod : symfield_proj.f90 input.mod global.mod symd3t_proj.mod symtalm.mod
	$(FORTRAN) $(FLAGS) -c $<

symden_proj.mod : symden_proj.f90 input.mod math.mod global.mod symd3t_proj.mod \
	nucleus.mod symfield_proj.mod deter.mod
	$(FORTRAN) $(FLAGS) -c $<

symgdhf_proj.mod : symgdhf_proj.f90 input.mod math.mod symd3t_proj.mod symfield_proj.mod
	$(FORTRAN) $(FLAGS) -c $<

symgden_proj.mod : symgden_proj.f90 input.mod symd3t_proj.mod nucleus.mod symfield_proj.mod \
	symden_proj.mod symgdhf_proj.mod
	$(FORTRAN) $(FLAGS) -c $<

energy_proj.mod : energy_proj.f90 symgden_proj.mod input.mod symd3t_proj.mod nucleus.mod symfield_proj.mod \
	symden_proj.mod symgdhf_proj.mod
	$(FORTRAN) $(FLAGS) -c $<

eigenval.mod : eigenval.f90
	$(FORTRAN) $(FLAGS) -c $<

jacobi.mod : jacobi.f90
	$(FORTRAN) $(FLAGS) -c $<

diagmeth.mod : diagmeth.f90 input.mod symd3t_proj.mod r1r1.mod nucleus.mod \
	symden_proj.mod symgden_proj.mod symgdhf_proj.mod symke2b.mod symvbb.mod \
	symvc.mod symvls.mod symgdd_proj.mod selfc_proj.mod \
	brent_d1.f90 brent_d2.f90 indexx.mod energy_proj.mod jacobi.mod eigenval.mod global.mod
	$(FORTRAN) $(FLAGS) -c $<

symke2b.mod : symke2b.f90 input.mod global.mod symd3t_proj.mod symfield_proj.mod
	$(FORTRAN) $(FLAGS) -c $<

ibb.mod : ibb.f90 global.mod lgfactor.mod symtalm.mod
	$(FORTRAN) $(FLAGS) -c $<

symvbb.mod : symvbb.f90 input.mod angmom.mod symd3t_proj.mod symfield_proj.mod ibb.mod
	$(FORTRAN) $(FLAGS) -c $<

ic.mod : ic.f90 global.mod lgfactor.mod symtalm.mod
	$(FORTRAN) $(FLAGS) -c $<

symvc.mod : symvc.f90 input.mod angmom.mod symd3t_proj.mod symfield_proj.mod ic.mod
	$(FORTRAN) $(FLAGS) -c $<

ils.mod : ils.f90 symtalm.mod
	$(FORTRAN) $(FLAGS) -c $<

symvls.mod : symvls.f90 input.mod global.mod math.mod symd3t_proj.mod \
	symfield_proj.mod ils.mod
	$(FORTRAN) $(FLAGS) -c $<

symgdd_proj.mod : symgdd_proj.f90 input.mod global.mod symtalm.mod symfield_proj.mod
	$(FORTRAN) $(FLAGS) -c $<

selfc_proj.mod : selfc_proj.f90 input.mod nucleus.mod symden_proj.mod \
	symke2b.mod symvbb.mod symvc.mod symvls.mod symgdd_proj.mod
	$(FORTRAN) $(FLAGS) -c $<

main.o : main.f90 input.mod global.mod nucleus.mod symden_proj.mod \
	diagmeth.mod
	$(FORTRAN) $(FLAGS) -c $<

main_PAV : main.o $(GLOBAL_OBJECTS) \
	r1r1.o symden_proj.o nucleus.o indexx.o diagmeth.o selfc_proj.o \
	symfield_proj.o symke2b.o symvbb.o ibb.o ic.o ils.o symvc.o symvls.o \
	symgdd_proj.o symgden_proj.o energy_proj.o symgdhf_proj.o eigenval.o jacobi.o
	$(FORTRAN) $(FLAGS) -o $@ $^

clean :
	rm -f *.o *.mod *~ main
