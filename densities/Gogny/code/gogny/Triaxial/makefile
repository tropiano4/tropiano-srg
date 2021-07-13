#------------------------------------------------------------------
#|          Makefile to compile the triaxial HFB code             | 
#------------------------------------------------------------------
#| Requires the Linear algebra libraries    Lapack and Blas       |
#|                                                                |
#| An optimized version of Blas and Lapack for PC architectures   |
#| called Atlas can also be used                                  |
#------------------------------------------------------------------
#------------------------------------------------------------------
#|       Atlas + optimized lapack from debian package             |
#------------------------------------------------------------------
#
.f.o:
	$(FOR) -c $(FFLAGS) $*.f
#
# ----------------------------------------------------------------------
#                           Intel Fortran  compiler (Default)
# ----------------------------------------------------------------------
FOR=ftn
FFLAGS=-O3
LIBS= #-L/opt/intel/Compiler/11.1/059/mkl/lib/em64t/ -lmkl_lapack -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -liomp5
#LIBS= -L/opt/intel/Compiler/11.0/069/mkl/lib/em64t/ -lmkl_lapack -lmkl_core -lmkl_em64t -lmkl_intel_thread -lmkl_intel_lp64 -liomp5
#OBJINTEL=dlamch.o xerbla.o
# ----------------------------------------------------------------------
#                           gfortran SSE2 atlas
# ----------------------------------------------------------------------
ifeq ($(ARCH),gfortran) 
FOR=gfortran-4.3
LIBS= -L/usr/lib/sse2/atlas/ -llapack -lblas
FFLAGS=-O3 -malign-double -funroll-loops 
endif
# ----------------------------------------------------------------------
#                           gfortran 3dNOW atlas
# ----------------------------------------------------------------------
ifeq ($(ARCH),Athlon) 
FOR=gfortran-4.3
LIBS= -L/usr/lib/3dnow/atlas/ -llapack -lblas
FFLAGS=-O3 -malign-double -funroll-loops
endif
# ----------------------------------------------------------------------
#                       PowerPC linux + ESSL from IBM
# ----------------------------------------------------------------------
ifeq ($(ARCH),PPC) 
FOR=xlf_r
FFLAGS=-O3
LIBS=-L/gpfs/apps/LAPACK/3.0/32/ -llapack -lessl
endif
#
#------------------------------------------------------------------
#  Object files
#------------------------------------------------------------------
#
OBJECTS0= HFBG2_BETA2.o ITLOC.o ROKAPA.o QUASI.o CALCL.o ONEBODY.o UTIL.o  \
HFBFIELD.o bgdir.o gognyf.o pack.o printd.o second.o GUTS.o binit.o \
denr.o readdata.o RES.o COLLMAS.o YOCCOZ.o RESFLUC.o KIN2BC.o \
BPAISO.o BGECH.o TALMANF.o HOINIT.o
OBJECTS=$(OBJECTS0) $(OBJINTEL)
#
# ----------------------------------------------------------------------
#
#            Targets
#
#  HFBtriaxial:  Standard target 
#
# ----------------------------------------------------------------------
HFBtriaxial: $(OBJECTS)
	$(FOR) $(LFLAGS) $(OBJECTS) $(LIBS) -o HFBtri
#
old2wf: old2wf.o
	$(FOR) $(LFLAGS) old2wf.o $(LIBS) -o old2wf
#	
clean:
	rm COMDIM DIMPERM DIMTRIAX HFBtri *.o
#
distrib:
	tar zcvf HFBtriaxial.tgz HOW-TO LOG *.f makefile *.INP
# ----------------------------------------------------------------------
#                              Dependencies
# ----------------------------------------------------------------------
bgdir.o: bgdir.f COMDIM
BGECHCOU.o: BGECHCOU.f COMDIM
BGECH.o: BGECH.f COMDIM
binit.o: binit.f COMDIM DIMTRIAX
biowf.o: biowf.f DIMPERM
HOINIT.o: HOINIT.f DIMTRIAX MAXDIM
ONEBODY.o: ONEBODY.f MAXDIM
UTIL.o: UTIL.f MAXDIM
KIN2BC.o: KIN2BC.f DIMTRIAX MAXDIM
#
old2wf.o: old2wf.f DIMTRIAX
