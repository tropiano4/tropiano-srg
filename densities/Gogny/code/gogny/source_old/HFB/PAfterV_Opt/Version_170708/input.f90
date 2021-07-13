! Input parameters module
 MODULE input

	IMPLICIT NONE

	!-----------------------------------------------------------------------------------------------------------------------!
	!															!
	!					GLOBAL VARIABLES								!
	!															!
	!    Meaning of some parameters												!
	!      - Nsize	   : Maximum number of single-particle levels taken in the basis. Used for arrays declaration		!
	!      - Npoint	   : Number of points on the lattice									!
	!      - Lmax	   : Maximum orbital angular momentum l (determined automatically in the case of the HO basis)		!
	!      - Nmax	   : Maximum number of radial wave-functions (determined automatically in the case of the HO basis)	!
	!      - IndexWave : At convergence, the program gives the radial wave-functions of the s.p. in the 			!
	!                    canonical basis and of the q.p. in the HFB basis. IndexWave specifies which one should		!
	!                    be written on disk											!
	!      - Nunity	   : Number of radial wave-functions taken into account when the unity is inserted for the 1-body 	!
	!		     kinetic energy											!
	!      - IndexVecNL: Array that gives which basis wave-function corresponds to a given n and l				!
	!      - MeshStep  : Lattice discretization step denoted by h in the notes						!
	!      - RadMesh   : Vector containing the radius on the lattice							!
	!      - Energy	   : Single-particle energies of the basis								!
	!      - WaveFun   : Array containing the single-particle wave-functions of the basis states on the lattice		!
	!		     First index refers to the isospin, second to the lattice and third to the state			!
	!      - WaveDeri  : Same as above but for the derivatives of the basis wave-functions					!
	!      - anneal    : At each iteration, the density is taken as a mixing between the new density and the old one. This	!
	!		     is controlled by the anneal parameter								!
	!															!
	!-----------------------------------------------------------------------------------------------------------------------!
	
	! More variables (input parameters)
	
	DOUBLE PRECISION :: b_0 = 1.6       	! Oscillator length ("b")
	INTEGER :: N_0 = 8                   	! Number of shells
	INTEGER :: protons = 8, neutrons = 8    ! Number of protons and neutrons
	INTEGER :: Gogny = 1        		! 1 = D1S, 2 = D1, 3 = D1prime
	INTEGER :: Basis = 1			! 1 = HO basis, 2 = arbitrary basis 
	INTEGER :: CompHO = 0			! 1 = compatibility with the HO basis in terms of quantum numbers
						! 0 = Lmax and Nmax required
	
	INTEGER :: Nsize, Npoint, Lmax, Nmax, Isospin

        DOUBLE PRECISION :: Ecut = -1.0, Convergence = 1.e-3, anneal = 0.5, facPair = 1.0, ShiftLambda = 0.0, RadPoles = 1.0
	
	INTEGER :: Nlevel(1:2), Nunity = 20, IndexWave = 1, HFOnly=1, Optimization = 0, NITER_MAX = 10
	
	INTEGER, ALLOCATABLE :: LmaxiIso(:)
	INTEGER, ALLOCATABLE :: NmaxOfL(:)
	INTEGER, ALLOCATABLE :: Nmain(:,:), Lmain(:,:), Jmain(:,:)
	INTEGER, ALLOCATABLE :: IndexVecNL(:,:), IndexAux(:,:)
	
	DOUBLE PRECISION :: MeshStep
	DOUBLE PRECISION, ALLOCATABLE :: RadMesh(:)
	DOUBLE PRECISION, ALLOCATABLE :: Energy(:,:), WaveFun(:,:), WaveDeri(:,:)
		
 CONTAINS

	!-----------------------------------------------------------------------!
	!									!
	!       This subroutine reads the file "input.txt" that contains 	!
	!       the standard input to the program.				!
	!									!
	!-----------------------------------------------------------------------!
	
	SUBROUTINE Input_read

		INTEGER, PARAMETER :: file_desc = 17
		INTEGER file_error, N, L, i
		CHARACTER(LEN = 256) :: param
		INTEGER param_len

		OPEN(file_desc, FILE="data/input.txt", ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) STOP "Impossible to read: ./data/input.txt"

		param_len = IO_get_Param(file_desc, param)
		
		DO WHILE (param_len .GT. 0)
			IF (param .EQ. "b_0") THEN
				b_0 = IO_get_RealValue(file_desc)
			ELSE IF (param .EQ. "Slowing-factor") THEN
				anneal = IO_get_RealValue(file_desc)
			ELSE IF (param .EQ. "Iterations") THEN
				NITER_MAX = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "Convergence") THEN
				Convergence = IO_get_RealValue(file_desc)
			ELSE IF (param .EQ. "Optimization") THEN
				Optimization = IO_get_IntegerValue(file_desc)
                        ELSE IF (param .EQ. "Divergence") THEN
                                RadPoles = IO_get_RealValue(file_desc)
			ELSE IF (param .EQ. "IsoBasis") THEN
				Isospin = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "HFOnly") THEN
				HFOnly = IO_get_IntegerValue(file_desc)
                        ELSE IF (param .EQ. "FactorPairing") THEN
				facPair = IO_get_RealValue(file_desc)
                        ELSE IF (param .EQ. "ShiftLambda") THEN
                                ShiftLambda = IO_get_RealValue(file_desc)
			ELSE IF (param .EQ. "N_0") THEN
				N_0 = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "Lmax") THEN
				Lmax = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "Nmax") THEN
				Nmax = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "Unity") THEN
				Nunity = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "IndexWave") THEN
				IndexWave = IO_get_IntegerValue(file_desc)
                        ELSE IF (param .EQ. "Ecut") THEN
                                Ecut = IO_get_RealValue(file_desc)
			ELSE IF (param .EQ. "N") THEN
				neutrons = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "Z") THEN
				protons = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "GOGNY") THEN
				Gogny = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "Basis") THEN
				Basis = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "CompHO") THEN
			
				CompHO = IO_get_IntegerValue(file_desc)
				
				! Forcing Lmax to N_0 if oscillator basis or compatibility required
				IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) THEN
                                        Lmax = N_0
                                        Nmax = N_0
                                END IF
				
				! Forcing N_0 to Lmax if general basis used
				IF (Basis .EQ. 2 .OR. CompHO .EQ. 0) THEN
                                        N_0 = Lmax
                                END IF
				
			ELSE
				PRINT "(A,A)", "Attention! Unknown parameter: ", param(1:param_len)
				! We ignore the remaining lines
				READ (file_desc, "(A)") param
			END IF
			param_len = IO_get_Param(file_desc, param)
		END DO

		CLOSE(file_desc)
		
		! In case of a basis compatible with the harmonic oscillator basis, one forces
		! b (the oscillator length) to take a value proportional to the mass
		
		IF (b_0 .LT. 0.0 .AND. (Basis .EQ. 1 .OR. CompHO .EQ. 1)) b_0 = 1.05*REAL(protons + neutrons)**(1.0/6.0)
		
		WRITE(*,'(5X,"DATA READ FROM INPUT")')
		WRITE(*,'(5X,"====================",/)')
		
		PRINT "(A,F8.5)", "Oscillator Length (b_0)       : ", b_0
		PRINT "(A,F8.5)", "Slowing factor                : ", anneal
		PRINT "(A,I3)",   "Maximum Number of Iterations  : ", NITER_MAX
		PRINT "(A,ES9.2)","Convergence in Energy         : ", Convergence
                PRINT "(A,F8.5)", "Pairing Enhancement Factor    : ", facPair
                PRINT "(A,F8.5)", "Shift of Neutron Number       : ", ShiftLambda
		PRINT "(A,I3)",   "Optimization for BB Mat. Els. : ", Optimization
		PRINT "(A,I3)",   "Isospin of the Basis WF       : ", Isospin
		PRINT "(A,I3)",   "Index of the Wave-function    : ", IndexWave
		PRINT "(A,I3)",   "Pure HF Calculations (=1)     : ", HFOnly
		PRINT "(A,I3)",   "Number of shells (N_0)        : ", N_0
		PRINT "(A,I3)",   "Maximum Ang. Mom. (Lmax)      : ", Lmax
		PRINT "(A,I3)",   "Maximum n (nmax)              : ", Nmax
		PRINT "(A,I3)",   "Maximum n for Unity           : ", Nunity
		PRINT "(A,I3)",   "Index of the Wave-function    : ", IndexWave
		PRINT "(A,F8.5)", "Cut-off energy (if > 0)       : ", Ecut
		PRINT "(A,I3)",   "Number of neutrons (N)        : ", neutrons
		PRINT "(A,I3)",   "Number of protons  (Z)        : ", protons
		SELECT CASE (Gogny)
			CASE (1)
				PRINT "(A)", "Parameters of the Gogny force : D1S"
			CASE (2)
				PRINT "(A)", "Parameters of the Gogny force : D1"
			CASE (3)
				PRINT "(A)", "Parameters of the Gogny force : D1prime"
			CASE DEFAULT
				STOP "Invalid parameters of the Gogny force"
		END SELECT
		PRINT "(A,I3)",   "Type of Basis                 : ", Basis
		PRINT "(A,I3)",   "Compatibility with the HO     : ", CompHO
		
		
		RETURN
		
	END SUBROUTINE Input_read

        !-----------------------------------------------------------------------!
	!									!
        !   Subroutine that reads the spherical basis. It gives the total 	!
	!   number of s.p. levels, the s.p. energies, wave-functions, deri-	!
	!   vatives of the wave-function and quantum numbers n, l and j (the 	!
	!   latter is not used in this version of the code			!
	!									!
        !-----------------------------------------------------------------------!

        SUBROUTINE ReadBasis()
	
		INTEGER, PARAMETER :: file_unit = 12
		
	        INTEGER, ALLOCATABLE :: Ndummy(:)
		
		INTEGER file_error
		INTEGER param_len, stat
		INTEGER i, j, L, N, nb1, nb2, jtemp, kappa, Nmaxi, Lmaxi, Lref, ntheo, Imax
		
		DOUBLE PRECISION :: Emax, E
		
		CHARACTER(LEN = 6)   :: keyword
		
                ! Opening the file to read
		
		OPEN(file_unit, FILE="basis.txt", ACTION="READ", IOSTAT=file_error)
		
		IF (file_error .NE. 0) STOP "Impossible to open: basis.txt"

		READ(file_unit, FMT="(A)", ADVANCE="NO", IOSTAT=stat) keyword
		
                ! We read the input file by using keywords
		
		DO WHILE (keyword .NE. "FINISH")
			 
                        SELECT CASE (keyword)
			  
				CASE ("boxpar") 
			  
					READ(file_unit, *, IOSTAT=stat) Npoint,MeshStep
			      
					ALLOCATE(RadMesh(0:Npoint))
			      
					DO i = 0, Npoint
						RadMesh(i) = i*MeshStep
					END DO
			      
				CASE ("noflev") 
			  
					READ(file_unit, *, IOSTAT=stat) nb1,nb2,Nlevel(1),Nlevel(2)
			      
					Nsize = max(Nlevel(1),Nlevel(2))
			    
				!    
				! NEUTRONS - NEUTRONS - NEUTRONS - NEUTRONS - NEUTRONS - NEUTRONS
				!
			  
				! Neutron single-particle energies
				CASE ("esingn") 
			  
					IF (Nlevel(1) .EQ. 0) STOP "Nlevel = 0 in ReadBasis !"
			  
					ALLOCATE(Energy(Nsize,2))
			      
					READ(file_unit, *, IOSTAT=stat) (Energy(i,1),i=1,Nlevel(1))
			      
				! Neutron angular momentum j (in Jmain)
				CASE("spin_n")
			  
					IF (Nlevel(1) .EQ. 0) STOP "Nlevel = 0 in ReadBasis !"
			  
					ALLOCATE(Jmain(Nsize,2))
					ALLOCATE(Ndummy(Nlevel(1)))
			      
					READ(file_unit, *, IOSTAT=stat) (Ndummy(i),i=1,Nlevel(1))
			      
					DO i=1,Nlevel(1)
						Jmain(i,1) = 2*Ndummy(i) - 1
					END DO
			      
					DEALLOCATE(Ndummy)
			      			      
				! Neutron orbital angular momentum l (in Lmain)
				CASE ("kappan")
			  
					IF (Nlevel(1) .EQ. 0) STOP "Nlevel = 0 in ReadBasis !"
			  
					ALLOCATE(Lmain(Nsize,2))
					ALLOCATE(Ndummy(Nlevel(1)))
			      
					READ(file_unit, *, IOSTAT=stat) (Ndummy(i),i=1,Nlevel(1))
			      
					! Lmaxi is the maximum orbital angular momentum (used to define 
					! a fake number of shells N0)
					Lmaxi = 0
			      
					DO i=1,Nlevel(1)
						kappa = Ndummy(i)
						jtemp = iabs(kappa)
						Lmain(i,1) = jtemp - (1 - kappa/jtemp)/2
						IF (Lmain(i,1) .GT. Lmaxi) Lmaxi = Lmain(i,1)
					END DO
			      
					DEALLOCATE(Ndummy)
			      
					ALLOCATE(LmaxiIso(2))
					LmaxiIso(1) = Lmaxi
			      			      
				CASE ("numbrn")
			  
					IF (Nlevel(1) .EQ. 0) STOP "Nlevel = 0 in ReadBasis !"
			  
					ALLOCATE(Nmain(Nsize,2))
			      
					READ(file_unit, *, IOSTAT=stat) (Nmain(i,1),i=1,Nlevel(1))
			      
				CASE ("wave_n")
			  
					IF (Nlevel(1) .EQ. 0) STOP "Nlevel = 0 in ReadBasis !"
			  
					IF (Basis .EQ. 2 .AND. Isospin .EQ. 1) THEN
				
						ALLOCATE(WaveFun(0:Npoint,Nsize))
					
						DO j=1,Nlevel(1)
							READ(file_unit, *, IOSTAT=stat) (WaveFun(i,j),i=0,Npoint)
						END DO
				
					END IF
			      
			  	CASE ("derivn")
			  
					IF (Nlevel(1) .EQ. 0) STOP "Nlevel = 0 in ReadBasis !"
			  
					IF (Basis .EQ. 2 .AND. Isospin .EQ. 1) THEN
				
						ALLOCATE(WaveDeri(0:Npoint,Nsize))
			      
			      			DO j=1,Nlevel(1)
							READ(file_unit, *, IOSTAT=stat) (WaveDeri(i,j),i=0,Npoint)
						END DO
				
					END IF
			    
				!    
                          	! PROTONS - PROTONS - PROTONS - PROTONS - PROTONS - PROTONS
			  	!
			  
				CASE ("esingp") 
			  
					READ(file_unit, *, IOSTAT=stat) (Energy(i,2),i=1,Nlevel(2))
			      
				CASE( "spin_p")
			  
					ALLOCATE(Ndummy(Nlevel(2)))
			      
					READ(file_unit, *, IOSTAT=stat) (Ndummy(i),i=1,Nlevel(2))
			      
					DO i=1,Nlevel(2)
						Jmain(i,2) = 2*Ndummy(i) - 1
					END DO
			      
					DEALLOCATE(Ndummy)
			      
				CASE ("kappap")
			  
					ALLOCATE(Ndummy(Nlevel(2)))
			      
					READ(file_unit, *, IOSTAT=stat) (Ndummy(i),i=1,Nlevel(2))
			      
					Lmaxi = 0
			      
					DO i=1,Nlevel(2)
						kappa = Ndummy(i)
						jtemp = iabs(kappa)
						Lmain(i,2) = jtemp - (1 - kappa/jtemp)/2
						IF (Lmain(i,2) .GT. Lmaxi) Lmaxi = Lmain(i,2)
					END DO
			      
					DEALLOCATE(Ndummy)
			      
					LmaxiIso(2) = Lmaxi
			      
				CASE ("numbrp")
			  
					READ(file_unit, *, IOSTAT=stat) (Nmain(i,2),i=1,Nlevel(2))
			      
				CASE ("wave_p")
			  
					IF (Basis .EQ. 2 .AND. Isospin .EQ. 2) THEN
				
						DO j=1,Nlevel(2)
							READ(file_unit, *, IOSTAT=stat) (WaveFun(i,j),i=0,Npoint)
						END DO
					
					END IF
			      
			  	CASE ("derivp")
			  
					IF (Basis .EQ. 2 .AND. Isospin .EQ. 2) THEN
				
						DO j=1,Nlevel(2)
							READ(file_unit, *, IOSTAT=stat) (WaveDeri(i,j),i=0,Npoint)
						END DO
					
					END IF
			    
			  END SELECT
			  
                          ! Until we get to the end of the file, we read the next keyword

		          READ(file_unit, FMT="(A)", IOSTAT=stat) keyword

                          ! Before we exit this subroutine, we determine a few useful quantities
			  !	- First of all, we establish the correspondence between the index basis wave-function
			  !	  and the quantum numbers n and l which are used in the calculation of he matrix elements
			  !	- We also determine have the choice between taking all states with E < Ecut and N < Nmax
			  !       or all stated with L < Lmax and N < Nmax. 
			  !	- IndexAux is a special array used for the calculation of 1-body kinetic energy when inserting
			  !	  the unity: We need more wave-functions (more n values) for this.

		          IF (stat .NE. 0) THEN
		
                          	IF (Ecut .LT. 0.0) THEN
                                        
                                        ! IndexVecNL gives the correspondence N, L -> index i of the basis
                		        ALLOCATE(IndexVecNL(0:2*Nmax,0:2*Lmax))

                                        ALLOCATE(IndexAux(0:Nunity,0:2*Lmax))
				
			        	DO N = 0,2*Nmax
					        DO L = 0,2*Lmax
					        	IndexVecNL(N,L) = 0
					        END DO
				        END DO

                                        DO N=0,Nunity
                                                DO L=0,2*Lmax
                                                        IndexAux(N,L) = 0
                                                END DO
                                        END DO
			      
			      	        Emax = -1000.0
                                        
				        DO i=1,Nlevel(Isospin)
				                N = Nmain(i,Isospin)
				                L = Lmain(i,Isospin)
					        IF (N .LE. Nmax .AND. L .LE. Lmax) IndexVecNL(N,L) = i
					        IF (N .LE. Nmax .AND. L .LE. Lmax .AND. Energy(i,Isospin) .GT. Emax) Emax = Energy(i,Isospin)
                                                IF (N .LE. Nunity .AND. L .LE. (Lmax+1)) IndexAux(N,L) = i
				        END DO

                                ELSE
                                        
                                       DO i=1,Nlevel(Isospin)
			                        N = Nmain(i,Isospin)
				                L = Lmain(i,Isospin)
                                                E = Energy(i,Isospin)
					        IF (E .LT. Ecut .AND. L .LE. Lmax .AND. N .GT. Nmax) Nmax = N
				        END DO

                                        ! IndexVecNL gives the correspondence N, L -> index i of the basis
                		        ALLOCATE(IndexVecNL(0:2*Nmax,0:2*Lmax))

                                        ALLOCATE(IndexAux(0:Nunity,0:2*Lmax))

			        	DO N = 0,2*Nmax
					        DO L = 0,2*Lmax
					        	IndexVecNL(N,L) = 0
					        END DO
				        END DO
			      
                                        DO N=0,Nunity
                                                DO L=0,2*Lmax
                                                        IndexAux(N,L) = 0
                                                END DO
                                        END DO
			      
			      	        Emax = Ecut
                                        
				        DO i=1,Nlevel(Isospin)
				                N = Nmain(i,Isospin)
				                L = Lmain(i,Isospin)
                                                E = Energy(i,Isospin)
					        IF (E .LT. Ecut .AND. L .LE. Lmax) IndexVecNL(N,L) = i
                                                IF (N .LE. Nunity .AND. L .LE. (Lmax+1)) IndexAux(N,L) = i
				        END DO

                                        
                                END IF

                		ALLOCATE(NmaxOfL(0:Lmax))
				
			 	 ! Printing the characteristics of the basis with respect to quantum numbers
				 
				WRITE(*,'()')
				WRITE(*,'(5X,"DEFINING THE BASIS - CUT-OFFS")')
				WRITE(*,'(5X,"=============================",/)')
		
                                WRITE(*,'("    Actual Cut-off Energy : Ecut = ",F10.3)') Emax                                
                                WRITE(*,'("    Maximum L value       : lmax = ",I7)') Lmax
                                WRITE(*,'("    Maximum n value       : Nmax = ",I7)') Nmax
				WRITE(*,'("    Maximum n-value of a give l",//,"       L    n  n(HO)")')
				
				DO Lref = 0, Lmax
				
					Nmaxi = 0
					
					DO i=1, Nlevel(Isospin)
						N = Nmain(i,Isospin)
						L = Lmain(i,Isospin)
                                                E = Energy(i,Isospin)
						IF (L .EQ. Lref .AND. E .LE. Emax) THEN
					                IF (N .GE. Nmaxi) Nmaxi = N
						END IF
					END DO
					
					ntheo = (N_0 - Lref)/2 + 1
					
					IF (Ecut .LT. 0.0) NmaxOfL(Lref) = MIN(Nmaxi, Nmax)
					IF (Ecut .GT. 0.0) NmaxOfL(Lref) = Nmaxi
					
					IF (CompHO .EQ. 1) NmaxOfL(Lref) = ntheo
					
					WRITE(*,'(3X,3I5)') Lref,Nmaxi,ntheo
				END DO
			      				
		          	CLOSE(file_unit)
				
		          	RETURN
		   
		          END IF
		
		END DO
	
		RETURN
	END SUBROUTINE ReadBasis
	
	FUNCTION IO_get_Param(fd, param)
		INTEGER IO_get_Param
		INTEGER, INTENT(IN) :: fd ! Descriptor de fichero
		CHARACTER(*), INTENT(INOUT) :: param

		CHARACTER c
		INTEGER num, stat

		READ (fd, FMT="(A1)", ADVANCE="NO", IOSTAT=stat) c
		IF (stat .NE. 0) THEN
			IO_get_Param = 0
			RETURN
		END IF
		! Ignoramos los espacios en blanco iniciales
		DO WHILE (c == " ")
			READ (fd, FMT="(A1)", ADVANCE="NO", IOSTAT=stat) c
			IF (stat .NE. 0) THEN
				IO_get_Param = 0
				RETURN
			END IF
		END DO
		num = 0
		DO WHILE ((c .NE. " ") .AND. (c .NE. "="))
			IF (num .EQ. 0) THEN
				param = c
			ELSE
				param = param(1:num) // c
			END IF
			num = num + 1
			READ (fd, FMT="(A1)", ADVANCE="NO", IOSTAT=stat) c
			IF (stat .NE. 0) THEN
				IO_get_Param = 0
				RETURN
			END IF
		END DO
		! Si no se encontra ningun nombre de parametro, salimos con error
		IF (num .EQ. 0) THEN
			IO_get_Param = 0
			RETURN
		END IF

		DO WHILE (c .EQ. " ")
			READ (fd, FMT="(A1)", ADVANCE="NO", IOSTAT=stat) c
			IF (stat .NE. 0) THEN
				IO_get_Param = 0
				RETURN
			END IF
		END DO

		IF (c .EQ. "=") THEN
			IO_get_Param = num
		ELSE
			IO_get_Param = 0
		END IF
		RETURN
	END FUNCTION IO_get_Param

	FUNCTION IO_get_IntegerValue(fd)
		INTEGER IO_get_IntegerValue
		INTEGER, INTENT(IN) :: fd ! Descriptor de fichero

		INTEGER stat

		READ (fd, FMT="(I)", IOSTAT=stat) IO_get_IntegerValue
		IF (stat .NE. 0) STOP "Imposible leer parametro de entrada"
		RETURN
	END FUNCTION IO_get_IntegerValue

	FUNCTION IO_get_RealValue(fd)
		DOUBLE PRECISION IO_get_RealValue
		INTEGER, INTENT(IN) :: fd ! Descriptor de fichero

		INTEGER stat

		READ (fd, FMT="(E)", IOSTAT=stat) IO_get_RealValue
		IF (stat .NE. 0) STOP "Imposible leer parametro de entrada"
		RETURN
	END FUNCTION IO_get_RealValue

	FUNCTION IO_get_String(fd)
		CHARACTER(LEN=30) IO_get_String
		INTEGER, INTENT(IN) :: fd ! Descriptor de fichero

		INTEGER stat

		READ (fd, FMT="(A)", IOSTAT=stat) IO_get_String
		IF (stat .NE. 0) STOP "Imposible leer parametro de entrada"
		RETURN
	END FUNCTION IO_get_String

END MODULE input
