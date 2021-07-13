PROGRAM main

 	USE symd3t
 	USE symvbb
 	USE symvc
 
	IMPLICIT NONE
 
	TYPE (SymVBBph) :: vBBph
        TYPE (SymVBBpp) :: vBBpp
        
	TYPE (SymVCph) :: vCph
        TYPE (SymVCpp) :: vCpp
        
	INTEGER :: Lmax, Nsmall, Nmax
	
	DOUBLE PRECISION :: b = -1.0
	
	READ(*,*) Nmax, Nsmall, Lmax

	! BRINK-BOEKER MATRIX ELEMENTS
	
        ! Particle-hole Channel
        
	CALL SymVBBph_new(vBBph, b, Lmax, Nmax)

	CALL SymVBBph_OpenFiles(vBBph, Lmax, Nsmall)
	
	CALL SymVBBph_read(vBBph, Lmax, Nmax)
	
	CALL SymVBBph_write(vBBph, Lmax, Nsmall)
	
	CALL SymVBBph_del(vBBph, Lmax, Nmax)
	
	CALL SymVBBph_CloseFiles()

        ! Particle-particle Channel

        CALL SymVBBpp_new(vBBpp, b, Lmax, Nmax)

	CALL SymVBBpp_OpenFiles(vBBpp, Lmax, Nsmall)
	
	CALL SymVBBpp_read(vBBpp, Lmax, Nmax)
	
	CALL SymVBBpp_write(vBBpp, Lmax, Nsmall)
	
	CALL SymVBBpp_del(vBBpp, Lmax, Nmax)
	
	CALL SymVBBpp_CloseFiles()

	! COULOMB MATRIX ELEMENTS
	
        ! Particle-hole Channel
        
	CALL SymVCph_new(vCph, Lmax, Nmax)

	CALL SymVCph_OpenFiles(Lmax, Nsmall)
	
	CALL SymVCph_read(vCph, Lmax, Nmax)
	
	CALL SymVCph_write(vCph, Lmax, Nsmall)
	
	CALL SymVCph_del(vCph, Lmax, Nmax)
	
	CALL SymVCph_CloseFiles()

        ! Particle-particle Channel

        CALL SymVCpp_new(vCpp, Lmax, Nmax)

	CALL SymVCpp_OpenFiles(Lmax, Nsmall)
	
	CALL SymVCpp_read(vCpp, Lmax, Nmax)
	
	CALL SymVCpp_write(vCpp, Lmax, Nsmall)
	
	CALL SymVCpp_del(vCpp, Lmax, Nmax)
	
	CALL SymVCpp_CloseFiles()
	
	STOP
 
END PROGRAM main

