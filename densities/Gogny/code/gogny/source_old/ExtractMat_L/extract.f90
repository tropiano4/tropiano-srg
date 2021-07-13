PROGRAM main

 	USE symd3t
 	USE symvbb
 	USE symvc
 
	IMPLICIT NONE
 
	TYPE (SymVBBph) :: vBBph
        TYPE (SymVBBpp) :: vBBpp
        
	TYPE (SymVCph) :: vCph
        TYPE (SymVCpp) :: vCpp
        
	INTEGER :: Lmax, Lsmall, Nmax
	
	DOUBLE PRECISION :: b = -1.0
	
	READ(*,*) Lmax, Lsmall, Nmax

	! BRINK-BOEKER MATRIX ELEMENTS
	
        ! Particle-hole Channel
        
	CALL SymVBBph_new(vBBph, b, Lmax, Nmax)

	CALL SymVBBph_OpenFiles(vBBph, Lmax, Lsmall)
	
	CALL SymVBBph_read(vBBph, Lmax, Nmax)
	
	CALL SymVBBph_write(vBBph, Lsmall, Nmax)
	
	CALL SymVBBph_del(vBBph, Lmax, Nmax)
	
	CALL SymVBBph_CloseFiles()

        ! Particle-particle Channel

        CALL SymVBBpp_new(vBBpp, b, Lmax, Nmax)

	CALL SymVBBpp_OpenFiles(vBBpp, Lmax, Lsmall)
	
	CALL SymVBBpp_read(vBBpp, Lmax, Nmax)
	
	CALL SymVBBpp_write(vBBpp, Lsmall, Nmax)
	
	CALL SymVBBpp_del(vBBpp, Lmax, Nmax)
	
	CALL SymVBBpp_CloseFiles()

	! COULOMB MATRIX ELEMENTS
	
        ! Particle-hole Channel
        
	CALL SymVCph_new(vCph, Lmax, Nmax)

	CALL SymVCph_OpenFiles(Lmax, Lsmall)
	
	CALL SymVCph_read(vCph, Lmax, Nmax)
	
	CALL SymVCph_write(vCph, Lsmall, Nmax)
	
	CALL SymVCph_del(vCph, Lmax, Nmax)
	
	CALL SymVCph_CloseFiles()

        ! Particle-particle Channel

        CALL SymVCpp_new(vCpp, Lmax, Nmax)

	CALL SymVCpp_OpenFiles(Lmax, Lsmall)
	
	CALL SymVCpp_read(vCpp, Lmax, Nmax)
	
	CALL SymVCpp_write(vCpp, Lsmall, Nmax)
	
	CALL SymVCpp_del(vCpp, Lmax, Nmax)
	
	CALL SymVCpp_CloseFiles()
	
	STOP
 
END PROGRAM main

