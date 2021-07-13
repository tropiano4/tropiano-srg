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

	! Input parameters

	LOGICAL, PUBLIC :: regularized_Gaussian = .False.

	INTEGER, PUBLIC :: protons = 8, neutrons = 8    ! Number of protons and neutrons
	INTEGER, PUBLIC :: N_0 = 8                   	! Number of shells
	INTEGER, PUBLIC :: NITER_MAX = 10               ! Maximum number of iterations
	INTEGER, PUBLIC :: HFonly = 0                   ! Switches to HF mode (1) instead of default HFB (0)
	INTEGER, PUBLIC :: Gogny = 1        		! 1 = D1S, 2 = D1, 3 = D1prime, 4 = D1N
	INTEGER, PUBLIC :: Basis = 1			! 1 = HO basis, 2 = arbitrary basis
	INTEGER, PUBLIC :: CompHO = 1			! 1 = compatibility with the HO basis in terms of quantum numbers
						        ! 0 = Lmax and Nmax required
	INTEGER, PUBLIC :: Optimization = 0             ! Optimizes calculation of matrix elements for WS-like bases
	INTEGER, PUBLIC :: Nunity = 20                  ! Maximum number of n values for resolution of the identity
	INTEGER, PUBLIC :: n_reg_max = 0                ! Maximum number of Gaussians used for regularized delta forces
	INTEGER, PUBLIC :: Nsize = 2, Npoint            ! Npoint: number of points of the mesh
	INTEGER, PUBLIC :: Lmax, Nmax                   ! Lmax: maximum number of L values in basis, Nmax: maximum values of n
	INTEGER, PUBLIC :: Isospin = 1                  ! Isospin of the basis functions (1: neutrons, 2: protons)
        INTEGER, PUBLIC :: IsoFacPair = 2               ! multiplies pairing fields for protons (0), neutrons (1) or both (2)
        INTEGER, PUBLIC :: IsoShiftLambda = 2           ! shifts particle number condition for p (0), n(1) or both (2)

        INTEGER, PUBLIC :: switch_Coulomb = 2           ! 0: no Coulomb, 1: direct only, 2: direct + full exchange
        INTEGER, PUBLIC :: switch_LS = 1                ! 0: no spin-orbit, 1: with spin-orbit
        INTEGER, PUBLIC :: switch_CM = 1                ! 0: no 2-body center of mass, 1: with two-body center of mass
        INTEGER, PUBLIC :: switch_DD = 1                ! 0: no density dependence, 1: with density dependence

	DOUBLE PRECISION, PUBLIC :: b_0 = 1.6D0      	! Oscillator length ("b")
        DOUBLE PRECISION, PUBLIC :: Ecut = -1.0D0       ! Cut-off in energy of the basis
        DOUBLE PRECISION, PUBLIC :: Emax = -1.0D0       ! Actual value of the cut-off energy in the basis
        DOUBLE PRECISION, PUBLIC :: Convergence = 1.D-3 ! Criterion to stop convergence
        DOUBLE PRECISION, PUBLIC :: anneal = 0.5D0      ! Convergence
        DOUBLE PRECISION, PUBLIC :: facPair = 1.0D0     ! Value of pairing enhancement if IsoFacPair > 2
        DOUBLE PRECISION, PUBLIC :: ShiftLambda = 0.0D0 ! Value of shifted particle number if IsoShiftLambda>0
        DOUBLE PRECISION, PUBLIC :: RadPoles = 1.0D0    ! ...?
	DOUBLE PRECISION, PUBLIC :: range1 = 0.7D0, range2 = 1.2D0 ! Range of the Gaussians (initialized to D1S values)


	INTEGER, PUBLIC :: Nlevel(1:2), IndexWave = 1

	INTEGER, ALLOCATABLE, PUBLIC :: LmaxiIso(:), NmaxOfL(:)
	INTEGER, ALLOCATABLE, PUBLIC :: Nmain(:,:), Lmain(:,:), Jmain(:,:), IndexVecNL(:,:), IndexAux(:,:)

	DOUBLE PRECISION, PUBLIC :: MeshStep
	DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: RadMesh(:)
	DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: Energy(:,:), WaveFun(:,:), WaveDeri(:,:)

 CONTAINS

	!-----------------------------------------------------------------------!
	!									!
	!       This subroutine reads the file "input.txt" that contains 	!
	!       the standard input to the program.				!
	!									!
	!-----------------------------------------------------------------------!

	SUBROUTINE Input_read

		INTEGER, PARAMETER :: file_desc = 17
		INTEGER :: file_error
		CHARACTER(LEN = 256) :: param
		INTEGER :: i_reg, param_len

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
                        ELSE IF (param .EQ. "IsoFactorPairing") THEN
				IsoFacPair = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "ShiftLambda") THEN
                                ShiftLambda = IO_get_RealValue(file_desc)
                        ELSE IF (param .EQ. "IsoShiftLambda") THEN
                                IsoShiftLambda = IO_get_IntegerValue(file_desc)

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
			ELSE IF (param .EQ. "Regularized") THEN
				i_reg = IO_get_IntegerValue(file_desc)
				IF (i_reg .EQ. 1) THEN
				        regularized_Gaussian = .True.
				ELSE
				        regularized_Gaussian = .False.
				END IF
			ELSE IF (param .EQ. "OrderExpansion") THEN
				n_reg_max = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "Range1") THEN
				range1 = IO_get_RealValue(file_desc)
			ELSE IF (param .EQ. "Range2") THEN
				range2 = IO_get_RealValue(file_desc)

			ELSE IF (param .EQ. "Coulomb") THEN
				switch_Coulomb = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "SpinOrbit") THEN
				switch_LS = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "CenterOfMass") THEN
				switch_CM = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "DensityDependence") THEN
				switch_DD = IO_get_IntegerValue(file_desc)

			ELSE IF (param .EQ. "Basis") THEN
				Basis = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "CompHO") THEN
				CompHO = IO_get_IntegerValue(file_desc)
			ELSE
				WRITE(6,'("Attention! Unknown parameter: ",A)') param(1:param_len)
				! We ignore the remaining lines
				READ (file_desc, "(A)") param
			END IF
			param_len = IO_get_Param(file_desc, param)
		END DO

		CLOSE(file_desc)

		IF (Basis .EQ. 1) CompHO = 1

		! Forcing Lmax to N_0 if oscillator basis or compatibility required
		IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) THEN
                        Lmax = N_0
                        Nmax = N_0
                        Nsize = (2*Lmax+1)*(2*Nmax+1)
                END IF

		! Forcing N_0 to Lmax if general basis used
		IF (Basis .EQ. 2 .OR. CompHO .EQ. 0) THEN
                        N_0 = Lmax
                END IF

		! In case of a basis compatible with the harmonic oscillator basis, one forces
		! b (the oscillator length) to take a value proportional to the mass
		IF (b_0 .LT. 0.0D0 .AND. (Basis .EQ. 1 .OR. CompHO .EQ. 1)) b_0 = 1.05D0*DBLE(protons + neutrons)**(1.0D0/6.0D0)

		! Printing input data
		WRITE(6,'(5X,"DATA READ FROM INPUT")')
		WRITE(6,'(5X,"====================",/)')
		WRITE(6,'("Oscillator Length (b_0) ........: ",F15.12)') b_0
		WRITE(6,'("Slowing factor .................: ",F8.5)') anneal
		WRITE(6,'("Maximum Number of Iterations ...: ",I3)') NITER_MAX
		WRITE(6,'("Convergence in Energy ..........: ",ES9.2)') Convergence
		WRITE(6,'("Optimization for BB Mat. Els. ..: ",I3)') Optimization
		WRITE(6,'("Pairing Enhancement Factor .....: ",F8.5)') facPair
		WRITE(6,'("Shift in Lagrange Parameter ....: ",F8.5)') ShiftLambda
		WRITE(6,'("Isospin of the Basis WF ........: ",I3)') Isospin
		WRITE(6,'("Index of the Wave-function .....: ",I3)') IndexWave
		WRITE(6,'("Pure HF Calculations (=1) ......: ",I3)') HFOnly
		WRITE(6,'("Number of shells (N_0) .........: ",I3)') N_0
		WRITE(6,'("Maximum Ang. Mom. (Lmax) .......: ",I3)') Lmax
		WRITE(6,'("Maximum n (nmax) ...............: ",I3)') Nmax
		WRITE(6,'("Maximum n for Unity ............: ",I3)') Nunity
		WRITE(6,'("Index of the Wave-function .....: ",I3)') IndexWave
		WRITE(6,'("Cut-off energy (if > 0) ........: ",F8.5)') Ecut
		WRITE(6,'("Number of neutrons (N) .........: ",I3)') neutrons
		WRITE(6,'("Number of protons  (Z) .........: ",I3)') protons
		IF (regularized_Gaussian) THEN
		        WRITE(6,'("Regularized Gaussians ...........")')
		        WRITE(6,'("    - order of expansion .......: ",I2)') n_reg_max
		        WRITE(6,'("    - range1 ...................: ",F10.7)') range1
		ELSE
		        WRITE(6,'("Original Gogny interaction.......")')
		        WRITE(6,'("    - range1 ...................: ",F10.7)') range1
		        WRITE(6,'("    - range2 ...................: ",F10.7)') range2
		        SELECT CASE (Gogny)
				CASE (1)
					WRITE(6,'("    - Parametrization ..........: D1S")')
				CASE (2)
					WRITE(6,'("    - Parametrization ..........: D1")')
				CASE (3)
					WRITE(6,'("    - Parametrization ..........: D1p")')
				CASE (4)
					WRITE(6,'("    - Parametrization ..........: D1N")')
				CASE DEFAULT
					STOP "Invalid parameters of the Gogny force"
			END SELECT
		END IF
		IF (switch_CM .EQ. 0) THEN
			WRITE(6,'("    - 2-body center of mass correction turned off")')
		END IF
		IF (switch_LS .EQ. 0) THEN
			WRITE(6,'("    - Spin-orbit term turned off")')
		END IF
		IF (switch_DD .EQ. 0) THEN
			WRITE(6,'("    - Density-dependent term turned off")')
		END IF
		IF (switch_Coulomb .EQ. 0) THEN
			WRITE(6,'("    - Both direct and exchange Coulomb terms switched off")')
		END IF
		IF (switch_Coulomb .EQ. 1) THEN
			WRITE(6,'("    - Exchange Coulomb terms switched off")')
		END IF
		WRITE(6,'("Type of Basis ..................: ",I3)') Basis
		WRITE(6,'("Compatibility with the HO ......: ",I3)') CompHO

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

        SUBROUTINE ReadBasisWS()

		INTEGER, PARAMETER :: file_unit = 12

	        INTEGER, ALLOCATABLE :: Ndummy(:)

		INTEGER :: file_error
		INTEGER :: stat, icount
		INTEGER :: i, j, L, N, nb1, nb2, jtemp, kappa, Nmaxi, Lmaxi, Lref, ntheo

		DOUBLE PRECISION :: E

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

					IF (ALLOCATED(RadMesh)) DEALLOCATE(RadMesh)
					ALLOCATE(RadMesh(0:Npoint))
					DO i = 0, Npoint
						RadMesh(i) = DBLE(i)*MeshStep
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

					ALLOCATE(WaveFun(0:Npoint,Nsize))

					IF (Basis .EQ. 2 .AND. Isospin .EQ. 1) THEN

						DO j=1,Nlevel(1)
							READ(file_unit, *, IOSTAT=stat) (WaveFun(i,j),i=0,Npoint)
						END DO

					END IF

			  	CASE ("derivn")

					IF (Nlevel(1) .EQ. 0) STOP "Nlevel = 0 in ReadBasis !"

					ALLOCATE(WaveDeri(0:Npoint,Nsize))

					IF (Basis .EQ. 2 .AND. Isospin .EQ. 1) THEN

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
			  !	  and the quantum numbers n and l which are used in the calculation of the matrix elements
			  !	- We also determine have the choice between taking all states with E < Ecut and N < Nmax
			  !       or all stated with L < Lmax and N < Nmax.
			  !	- IndexAux is a special array used for the calculation of 1-body kinetic energy when inserting
			  !	  the unity: We need more wave-functions (more n values) for this.

		          IF (stat .NE. 0) THEN

                          	IF (Ecut .LT. 0.0D0) THEN

                                        ! IndexVecNL gives the correspondence N, L -> index i of the basis
                		        IF (.NOT.ALLOCATED(IndexVecNL)) ALLOCATE(IndexVecNL(0:2*Nmax,0:2*Lmax))

                                        ALLOCATE(IndexAux(0:Nunity,0:2*Lmax))

					DO L = 0,2*Lmax
			        		DO N = 0,2*Nmax
					        	IndexVecNL(N,L) = 0
					        END DO
					END DO

                                        DO L=0,2*Lmax
                                        	DO N=0,Nunity
                                                        IndexAux(N,L) = 0
                                                END DO
                                        END DO

			      	        Emax = -1000.0D0

				        DO i=1,Nlevel(Isospin)
				                N = Nmain(i,Isospin)
				                L = Lmain(i,Isospin)
					        IF (N .LE. Nmax .AND. L .LE. Lmax) IndexVecNL(N,L) = i
					        IF (N .LE. Nmax .AND. L .LE. Lmax .AND. Energy(i,Isospin) .GT. Emax) Emax = Energy(i,Isospin)
                                                IF (N .LE. Nunity .AND. L .LE. (Lmax+1)) IndexAux(N,L) = i
                                                write(6,'("i = ",i6," Emax = ",f20.14)') i,Emax
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

		          	CLOSE(file_unit)

		          	RETURN

		          END IF

		END DO

		RETURN
	END SUBROUTINE ReadBasisWS

        !-----------------------------------------------------------------------!
	!									!
        !   Subroutine that defines the correspondence between basis indexes	!
        !   and quantum numbers in the case of the analytical HO basis.	        !
	!									!
        !-----------------------------------------------------------------------!

	SUBROUTINE ReadBasisHO()
		INTEGER :: icount,N,L

                IF (.NOT.ALLOCATED(IndexVecNL)) ALLOCATE(IndexVecNL(0:2*Nmax,0:2*Lmax))

                icount = 0
                DO N = 0,2*Nmax
                	DO L = 0,2*Lmax
                		icount = icount + 1
                		IndexVecNL(N,L) = icount
			END DO
		END DO

		RETURN
	END SUBROUTINE ReadBasisHO

        !-----------------------------------------------------------------------!
	!									!
        !   Subroutine prints basic characteristics of the basis		!
	!									!
        !-----------------------------------------------------------------------!

	SUBROUTINE PrintBasis()
		INTEGER :: Lref,i,N,L,Nmaxi,ntheo
		DOUBLE PRECISION :: E

	        ALLOCATE(NmaxOfL(0:Lmax))

		! Printing the characteristics of the basis with respect to quantum numbers
		WRITE(6,'(/,5X,"DEFINING THE BASIS - CUT-OFFS")')
		WRITE(6,'(5X,"=============================",/)')
		IF (Basis.EQ.1) THEN
			WRITE(6,'("Analytical HO basis ............. ")')
		ELSE
			WRITE(6,'("Numerical basis on a mesh ....... ")')
			WRITE(6,'("Mesh size ......................: ",F10.5)') MeshStep
			WRITE(6,'("Box radius .....................: ",F10.5)') RadMesh(Npoint)
			WRITE(6,'("Number of mesh points ..........: ",I5)') Npoint+1
		END IF
		WRITE(6,'("Compatibility with the HO ......: ",I3)') CompHO
                WRITE(6,'("Actual Cut-off Energy Ecut .....: ",F10.3)') Emax
                WRITE(6,'("Maximum L value lmax ...........: ",I7)') Lmax
                WRITE(6,'("Maximum n value ................: ",I7)') Nmax
		IF (Basis.EQ.1) THEN
			WRITE(6,'("Maximum n-value for given l ....",//,"       L   n(HO)")')
		ELSE
			WRITE(6,'("Maximum n-value for given l ....",//,"       L    n  n(HO)")')
		END IF

		IF (compHO .EQ. 1) THEN

			DO Lref = 0, Lmax
				ntheo = (N_0 - Lref)/2 + 1
				IF (CompHO .EQ. 1) NmaxOfL(Lref) = ntheo
				WRITE(6,'(3X,2I5)') Lref,ntheo
			END DO

		ELSE

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

				IF (Ecut .LT. 0.0) NmaxOfL(Lref) = MIN(Nmaxi, Nmax)
				IF (Ecut .GT. 0.0) NmaxOfL(Lref) = Nmaxi

				WRITE(6,'(3X,3I5)') Lref,Nmaxi,ntheo

			END DO

		END IF

		RETURN
	END SUBROUTINE PrintBasis

	FUNCTION IO_get_Param(fd, param)
		INTEGER IO_get_Param
		INTEGER, INTENT(IN) :: fd ! Descriptor de fichero
		CHARACTER(*), INTENT(INOUT) :: param

		CHARACTER (LEN=1) c
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

		READ (fd, *, IOSTAT=stat) IO_get_IntegerValue
		IF (stat .NE. 0) STOP "Imposible leer parametro de entrada"
		RETURN
	END FUNCTION IO_get_IntegerValue

	FUNCTION IO_get_RealValue(fd)
		DOUBLE PRECISION IO_get_RealValue
		INTEGER, INTENT(IN) :: fd ! Descriptor de fichero

		INTEGER stat

		READ (fd, *, IOSTAT=stat) IO_get_RealValue
		IF (stat .NE. 0) STOP "Imposible leer parametro de entrada"
		RETURN
	END FUNCTION IO_get_RealValue

	FUNCTION IO_get_String(fd)
		CHARACTER(LEN=30) IO_get_String
		INTEGER, INTENT(IN) :: fd ! Descriptor de fichero

		INTEGER stat

		READ (fd, *, IOSTAT=stat) IO_get_String
		IF (stat .NE. 0) STOP "Imposible leer parametro de entrada"
		RETURN
	END FUNCTION IO_get_String

END MODULE input
