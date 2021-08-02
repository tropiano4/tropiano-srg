!****************************************************************************
!  file: input_file.f90
!
!  Revision: 0.50
!
!  Subroutine that takes an input filename from the command
!   line and reads in parameters in the input file.
!
!  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
!
!  Revision history:
!    0.55  04-Oct-2006 --- upgraded to include Similarity RG option.
!    0.50  06-Aug-2006 --- adapted for standalone program
!    0.10  10-May-2006 --- original version of inputoutput.f90
!
!  Notes:
!    * uses the f2kcli command line routines
!    * f2kcli_ifort works with both Intel and Portland compilers 
!
!****************************************************************************
!***************************************************************************
! Specify the bare NN potential
!  kvnn ---  indicates the potential used
!     kvnn=1 Paris
!         =2 Bonn
!         =3 reid93
!         =4 nij1
!         =5 nij2
!         =6 v18
!         =7 cdbonn
!         =8 idaho(NNLO) 
!         =9 chiral(N3LO epelbaum) [see below for change]
!         10 idaho(N3LO) 500 MeV 
!         11 large cutoff LO+promoted cts
!         12 idaho(N3LO) 600 MeV
!         13 idaho(N3LO) 400 MeV (2.098 fm^{-1})
!  
!     kvnn = 20-24  Epelbaum chiral NLO:  ostat=0; cutnum=1-5 
!          = 25-29  Epelbaum chiral N2LO: ostat=1; cutnum=1-5 
!          = 30-34  Epelbaum chiral N3LO: ostat=2; cutnum=1-5 
!	 **  1) NLO (OSTAT=0):  					**! 
!	 **  kvnn	 CUTNUM 	LS cut-off	 SFR cut-off	**! 
!	 **   20	   1		 400 MeV	   500 MeV	**! 
!	 **   21	   2		 550 MeV	   500 MeV	**! 
!	 **   22	   3		 550 MeV	   600 MeV	**! 
!	 **   23	   4		 400 MeV	   700 MeV	**! 
!	 **   24	   5		 550 MeV	   700 MeV	**! 
!	 **  2) NNLO (OSTAT=1): 					**! 
!	 **  kvnn	 CUTNUM 	LS cut-off	 SFR cut-off	**! 
!	 **   25	   1		 450 MeV	   500 MeV	**! 
!	 **   26	   2		 600 MeV	   500 MeV	**! 
!	 **   27	   3		 550 MeV	   600 MeV	**! 
!	 **   28	   4		 450 MeV	   700 MeV	**! 
!	 **   29	   5		 600 MeV	   700 MeV	**! 
!	 **  3) NNNLO (OSTAT=2):					**! 
!	 **  kvnn	 CUTNUM 	LS cut-off	 SFR cut-off	**! 
!	 **   30	   1		 450 MeV	   500 MeV	**! 
!	 **   31	   2		 600 MeV	   600 MeV	**! 
!	 **   32	   3		 550 MeV	   600 MeV	**! 
!	 **   33	   4		 450 MeV	   700 MeV	**! 
!	 **   34	   5		 600 MeV	   700 MeV	**! 
!
!***************************************************************************
!
!***************************************************************************
!
! Parameters that go into the MethodsType structure "methods" to specify
!  the regulator.
!
!   ireg --- regulator used for vlowk (integer):         
!	      1 => theta function  [sharp --> \theta(\Lambda-k)]      
!	      2 => exponential  [exp(-(k^2/\Lambda^2)^n)]   	       
!	      3 => Woods-Saxon  []  	       
!	      4 => tanh  []		       
!	      5 => power law  []		       
!
!   nsmooth --- integer smoothness parameter for regulator functions 
!                exponential and power law, etc.
!   rsmooth --- real smoothness parameter for regulators Woods-Saxon, tanh               
!
!   imethod --- method for constructing vlowk (integer)
!                1 => Lee-Suzuki (only for ireg=1 at present)
!                2 => 3-step procedure 
!                3 => SRG
!
!   iherm --- Hermitization method (integer)	     
! 	       1 => Gram-Schmidt		    
! 	       2 => Okubo			    
! 	       3 => Cholesky			    
! 	       4 => Kato (really poor!) 	    
!
! Sharp example
!   ireg=1
!   imethod=1
!   iherm=1 

! Smooth exponential example
!   nsmooth=4     
!   rsmooth=0._dp 
!   ireg=2   
!   imethod=2 
!   iherm = 1
!
!***************************************************************************
!
! Parameters that go into the MeshType structure "kmesh" to specify
!  the mesh on which the potential is created.
!
!      ntot --- momentum mesh size for P+Q = full = bare space (integer)
!      nmod --- momentum mesh size for low-k = P space (integer)
!      kmax --- maximum momentum for bare potential (integer)
!    lambda --- momentum cutoff for vlowk potential (real)
!       fac --- actual mesh cutoff is at fac*Lambda (real)
!                Note: this should be 1.0 if using a sharp cutoff (ireg=1)
!
!  Smooth cutoff example (fac = 1.0_dp for sharp cutoff)
!    ntot = 100
!    nmod = 54
!    kmax = 30.0_dp
!    lambda = 2.0_dp
!    fac = 1.2_dp  
!
!***************************************************************************
!


MODULE vlowk_input
  USE f2kcli
  IMPLICIT NONE
  PRIVATE

  INTEGER :: kvnn, smooth, ireg, nsmooth, imethod, iherm, ntot, nmod
  INTEGER :: l, s, jt, itz, it 
  REAL(8) :: Lambda, fac, rsmooth, kmax 

  PUBLIC :: Input_parameters, Input_channel
  PUBLIC :: kvnn, smooth, ireg, nsmooth, imethod, iherm, ntot, nmod, &
            Lambda, fac, rsmooth, kmax, &
	    l, s, jt, itz, it
	     
  CONTAINS
  
!***************************************************************************

  SUBROUTINE Input_channel

    LOGICAL :: input_file_exists
    CHARACTER (len=70) :: input_file
    INTEGER :: input_unit = 43

    ! Get filename of the input data file from command line
    CALL GET_COMMAND_ARGUMENT(2,input_file)  

    ! Check whether the input file exists; if so, open it for reading
    INQUIRE (FILE=input_file, EXIST=input_file_exists)
    IF (input_file_exists) THEN
       WRITE (UNIT=*,FMT=*) "Opening input file ", TRIM(input_file), " ..." 
       OPEN (UNIT=input_unit, FILE=input_file, STATUS='OLD')
     ELSE
       WRITE (UNIT=*,FMT=*) "Input file ", TRIM(input_file), &
         " does not exist!"
       STOP  
    ENDIF

    ! Read from the input file and then close it
    READ (input_unit,*)
    READ (input_unit,*) l
    READ (input_unit,*) s
    READ (input_unit,*) jt
    READ (input_unit,*) itz
  
  CLOSE (input_unit)
  WRITE (UNIT=*,FMT=*) "... channel input complete."
  WRITE (UNIT=*,FMT=*) " "

  END SUBROUTINE Input_channel
!***************************************************************************
!***************************************************************************

  SUBROUTINE Input_parameters
  
    LOGICAL :: input_file_exists
    CHARACTER (len=70) :: input_file
    INTEGER :: input_unit = 10

    ! Get filename of the input data file from command line
    CALL GET_COMMAND_ARGUMENT(1,input_file)  

    ! Check whether the input file exists; if so, open it for reading
    INQUIRE (FILE=input_file, EXIST=input_file_exists)
    IF (input_file_exists) THEN
       WRITE (UNIT=*,FMT=*) "Opening input file ", TRIM(input_file), " ..." 
       OPEN (UNIT=input_unit, FILE=input_file, STATUS='OLD')
     ELSE
       WRITE (UNIT=*,FMT=*) "Input file ", TRIM(input_file), &
         " does not exist!"
       STOP  
    ENDIF

    ! Read from the input file and then close it
    READ (input_unit,*)
    READ (input_unit,*) kvnn
    READ (input_unit,*) Lambda
    READ (input_unit,*) smooth
    READ (input_unit,*) fac
    READ (input_unit,*) ireg
    READ (input_unit,*) nsmooth
    READ (input_unit,*) rsmooth
    READ (input_unit,*) imethod
    READ (input_unit,*) iherm
    READ (input_unit,*)
    READ (input_unit,*)
    READ (input_unit,*) ntot
    READ (input_unit,*) nmod
    READ (input_unit,*) kmax
    READ (input_unit,*)

    IF ((imethod.NE.3).AND.(smooth.EQ.0)) fac = 1.0  ! sharp cutoff needs fac=1.0    
  
    CLOSE (input_unit)
    WRITE (UNIT=*,FMT=*) "... input complete."
    WRITE (UNIT=*,FMT=*) " "

  END SUBROUTINE Input_parameters
!***************************************************************************
  
END MODULE vlowk_input

!***************************************************************************
!***************************************************************************
