 MODULE global

	USE input
	USE lgfactor
	USE symd3t_proj
	USE symtalm
	USE gauss
	USE angmom
	USE bessik

	IMPLICIT NONE

	! PUBLIC VARIABLES
	LOGICAL, PUBLIC :: test_regularization = .False.
	INTEGER, PUBLIC :: Ngaussian = 1
	DOUBLE PRECISION, PUBLIC :: delta_a = 1.D-10
	DOUBLE PRECISION, ALLOCATABLE :: range_gaussian(:)

	! We define the flags that will identify each isospin throughout the program
	INTEGER, PARAMETER, PUBLIC :: PROTON  = 0
	INTEGER, PARAMETER, PUBLIC :: NEUTRON = 1

	! The list of spherical magic numbers
	INTEGER, DIMENSION(0:8), PUBLIC :: MagicNumber
	DATA MagicNumber/ 2, 8, 20, 28, 50, 82, 126, 184, 256/

	! m is equal to mc2/(hbar*c)2

	DOUBLE PRECISION, DIMENSION(0:1), PUBLIC :: m
	!DATA m / 0.02411186800036718652D0, 0.02411186800036718652D0 /
	DATA m / 0.024111868D0, 0.024111868D0 /

	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, PUBLIC :: sq, sq2
	DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, PUBLIC :: BesselTable

	! Parameters of the Gogny force
	DOUBLE PRECISION, DIMENSION(0:1, 4), PUBLIC :: Gogny_W
	DOUBLE PRECISION, DIMENSION(0:1, 4), PUBLIC :: Gogny_B
	DOUBLE PRECISION, DIMENSION(0:1, 4), PUBLIC :: Gogny_H
	DOUBLE PRECISION, DIMENSION(0:1, 4), PUBLIC :: Gogny_M
	DOUBLE PRECISION, PARAMETER, PUBLIC :: x0 = 1.0D0

	DATA Gogny_W / -1720.30D0,  103.64D0, -402.40D0, -21.30D0, -402.40D0, -21.30D0, -2047.61D0,  293.02D0 /
	DATA Gogny_B /  1300.00D0, -163.48D0, -100.00D0, -11.77D0, -100.00D0, -11.77D0,  1700.00D0, -300.78D0 /
	DATA Gogny_H / -1813.53D0,  162.81D0, -496.20D0,  37.27D0, -496.20D0,  37.27D0, -2414.93D0,  414.59D0 /
	DATA Gogny_M /  1397.60D0, -223.93D0,  -23.56D0, -68.81D0,  -23.56D0, -68.81D0,  1519.35D0, -316.84D0 /

	DOUBLE PRECISION, DIMENSION(4), PUBLIC :: Gogny_W0
	DATA Gogny_W0 / 130.0D0, 115.0D0, 130.0D0, 115.0D0 / ! D1S, D1, D1prime, D1N

	DOUBLE PRECISION, DIMENSION(4), PUBLIC :: Gogny_t0
	DATA Gogny_t0 / 1390.6D0, 1350.0D0, 1350.0D0, 1609.46D0 / ! D1S, D1, D1prime, D1N

	! One-body matrix elements of the kinetic energy and Gauss-Laguerre quadrature points
	TYPE (SymD3Tensor), PUBLIC :: EkField, R2Field
	TYPE (GaussLaguerreQuadrature), PUBLIC :: GaussLQ

	! PRIVATE VARIABLES
	INTEGER :: A
	INTEGER :: NLag

 CONTAINS

	SUBROUTINE Global_new

		INTEGER :: max_1, max_2

		! Reading the input parameters
		CALL Input_read

		! Initialize a default mesh
		CALL InitializeMesh()

        	! Reading the basis form WSCOOR (only if general basis is used)
                IF (Basis .EQ. 2 .OR. CompHO .EQ. 0) THEN
                	CALL ReadBasisWS()
                ELSE
                	CALL ReadBasisHO()
                END IF

        	! Printing some basis characteristics of the basis
                CALL PrintBasis()

		! Printing the parameters of the force
		WRITE(*,'()')
		WRITE(*,'(5X,"PARAMETERS OF THE GOGNY FORCE")')
		WRITE(*,'(5X,"=============================",/)')
		WRITE(*,'("Protons : W = ",F8.2," MeV")') Gogny_W(1, Gogny)
		WRITE(*,'("        : B = ",F8.2," MeV")') Gogny_B(1, Gogny)
		WRITE(*,'("        : H = ",F8.2," MeV")') Gogny_H(1, Gogny)
		WRITE(*,'("        : M = ",F8.2," MeV",/)') Gogny_M(1, Gogny)
		WRITE(*,'("Neutrons: W = ",F8.2," MeV")') Gogny_W(0, Gogny)
		WRITE(*,'("        : B = ",F8.2," MeV")') Gogny_B(0, Gogny)
		WRITE(*,'("        : H = ",F8.2," MeV")') Gogny_H(0, Gogny)
		WRITE(*,'("        : M = ",F8.2," MeV",/)') Gogny_M(0, Gogny)
		WRITE(*,'("WLS = ",F8.2," MeV.fm5")') Gogny_W0(Gogny)
		WRITE(*,'("t_0 = ",F8.2," MeV.fm4",/)') Gogny_t0(Gogny)

		NLag = MIN(N_0 * 32, 156) ! NLag = UMIN(N_0 << 5, 156)

		max_1 = 100
		max_2 = 2000

		IF (regularized_Gaussian) THEN
			IF (test_regularization) THEN
			        Ngaussian=3; ALLOCATE(range_gaussian(0:Ngaussian-1))
				range_gaussian(0) = range1 + delta_a
				range_gaussian(1) = range1 - delta_a
				range_gaussian(2) = range1
			ELSE
			        Ngaussian=1; ALLOCATE(range_gaussian(0:Ngaussian-1))
				range_gaussian(0) = range1
			END IF
		ELSE
			Ngaussian=2; ALLOCATE(range_gaussian(0:Ngaussian-1))
			range_gaussian(0) = range1
			range_gaussian(1) = range2
		END IF

!TODO El autor lo implementa, pero no se utiliza en ninguna parte del código
!		CALL LogFactorials_new(max_2)
!		CALL LogSemiFactorials_new(max_2)

                ! Defined in module "lgfactor.f90"
		CALL GammaFunction_new(max_1)

                ! Defined in module "lgfactor.f90"
		CALL DDLogFactorials_new(max_2)
		CALL DDLogSemiFactorials_new(max_2)

                ! Defined here...
		CALL SquareRoot_new(max_1)
		CALL SemiSquareRoot_new(max_1)

                ! Calculate one-body kinetic energy and r.m.s. radius
		CALL Global_start

                ! Defined in module "symtalm.f90" (only if analytical HO basis is used)
		IF (Basis .EQ. 1) THEN

			CALL SymCoefficientB_new
			CALL SymKumar_new

			CALL GaussLaguerreQuadrature_new(GaussLQ, NLag, DBLE(0.5))

		END IF

                ! Defined in module "angmom.f90"
		CALL ThreeJSymbols_new

		RETURN

	CONTAINS

                ! Subroutine initializing a default box with the corresponding mesh

		SUBROUTINE InitializeMesh()
			INTEGER :: i
			DOUBLE PRECISION :: Rbox

			Npoint = 201
			Rbox = 20.0D0
			MeshStep = Rbox / DBLE(Npoint - 1)

			ALLOCATE(RadMesh(0:Npoint))

			DO i = 0, Npoint
				RadMesh(i) = DBLE(i)*MeshStep
			END DO

			RETURN
		END SUBROUTINE InitializeMesh

                ! Subroutine storing in a vector of size imax the square root of
		! all integers from 1 to imax: sqrt(i)

		SUBROUTINE SquareRoot_new(imax)
			INTEGER, INTENT(IN) :: imax

			INTEGER i

			ALLOCATE(sq(0:imax))
			sq(0) = 0.0D0
			DO i = 1, imax
				sq(i) = SQRT(DBLE(i))
			END DO
			RETURN
		END SUBROUTINE SquareRoot_new

                ! Subroutine storing in a vector of size imax the square root of
		! all integers from 1 to imax, plus one half: sqrt(i+0.5)

		SUBROUTINE SemiSquareRoot_new(imax)
			INTEGER, INTENT(IN) :: imax

			INTEGER i

			ALLOCATE(sq2(0:imax))
			sq2(0) = 0.5d0*SQRT(2.0D0)
			DO i = 1, imax
				sq2(i) = SQRT(i + 0.5d0)
			END DO
			RETURN
		END SUBROUTINE SemiSquareRoot_new

	END SUBROUTINE Global_new

        ! Subroutine storing in a vector of size imax the square root of
	! all integers from 1 to imax: sqrt(i)

	SUBROUTINE BesselTabularize(Norder)
		INTEGER, INTENT(IN) :: Norder

		INTEGER :: i_r1, i_r2, k, igauss

#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: xarg, BesselFunc
		REAL(KIND = 16) :: rk, rip, rkp, pi, Order, r1, r2
#else
		DOUBLE PRECISION :: xarg, BesselFunc
		DOUBLE PRECISION :: rk, rip, rkp, pi, Order, r1, r2
#endif
		DOUBLE PRECISION :: range

		ALLOCATE(BesselTable(1:Npoint,1:Npoint,0:Norder,0:Ngaussian-1))

		pi = 4.0D0*ATAN(1.0D0)

		WRITE(*,'("Tabularization of Bessel Functions - k_max = ",i4)') Norder

                DO igauss = 0, Ngaussian-1

                        range = range_gaussian(igauss)

                        DO k = 0, Norder

			        WRITE(*,'("Order k = ",i4," Range = ",f20.14)') k, range

!$OMP Parallel Default(None) &
!$OMP& SHARED(igauss,k,Npoint,BesselTable,range,RadMesh,pi) &
!$OMP& PRIVATE(i_r1,i_r2,Order,r1,r2,xarg,BesselFunc,rk,rip,rkp)
!$OMP DO SCHEDULE(DYNAMIC)
		                Order = DBLE(k) + 0.5D0

			        DO i_r2 = 1, Npoint

				        r2 = RadMesh(i_r2)

		                	DO i_r1 = 1, Npoint

		                		r1 = RadMesh(i_r1)

		                		xarg = 2.0D0*r1*r2/(range**2)
		                		CALL bessel(xarg, Order, BesselFunc, rk, rip, rkp)

		                		BesselTable(i_r1, i_r2, k, igauss) = BesselFunc * EXP(-(r1**2+r2**2)/(range**2)) &
		                		                                                * SQRT(0.5D0*pi/xarg)

		                	END DO
		                END DO
!$OMP End Do
!$OMP End Parallel
                        END DO ! end of order
                END DO ! end of # gaussians

		WRITE(*,'("...DONE")')

		RETURN
	END SUBROUTINE BesselTabularize

	SUBROUTINE BesselFree()

	        DEALLOCATE(BesselTable)

		RETURN
        END SUBROUTINE BesselFree

	!-------------------------------------------------------------------------------!
	!  In this subroutine, we calculate the one-body kinetic energy matrix elements	!
	!  as well as the one-body r.m.s radius matrix elements  			!
	!-------------------------------------------------------------------------------!

	SUBROUTINE Global_start

		INTEGER ta, la, na, nc
		DOUBLE PRECISION factor
		INTEGER nmaxi

		CHARACTER(64) filename
		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

	        !  Creating two new pointers to the 3D tensors EkField and R2Field and allocating
		!  memory for them

		CALL SymD3Tensor_new(EkField)
		CALL SymD3Tensor_new(R2Field)

	        !  Filling in these two pointers with the ONE-BODY kinetic energy and square radius
		!  matrix elements. We don't calculate <na nb | T | nc nd> here, but only things like
		!  <n'|T|n>. Same with the radius

		SELECT CASE(Basis)

		CASE(1)

			DO ta = 0, 1  ! Loop over isospin

				factor = 1.0D0 / (2.0D0 * m(ta))

				DO la = 0, N_0  ! Loop over "bra" orbital angular momentum

					nmaxi = ((N_0 - la) / 2) + 1

					DO na = 1, nmaxi  ! Sum over "bra" main quantum number
						EkField%d3tensor(la)%d2(na, na) = factor * nabla2HO(na - 1, na - 1, la)
						R2Field%d3tensor(la)%d2(na, na) = r2_cutonHO(na - 1, na - 1, la)
					END DO

					DO na = 1, nmaxi - 1

						EkField%d3tensor(la)%d2(na + 1, na) = factor * nabla2HO(na - 1, na, la)
						R2Field%d3tensor(la)%d2(na + 1, na) = r2_cutonHO(na - 1, na, la)

						DO nc = na + 2, nmaxi ! Loop over "ket" main quantum number
							EkField%d3tensor(la)%d2(nc, na) = 0.0D0
							R2Field%d3tensor(la)%d2(nc, na) = 0.0D0
						END DO

					END DO

				END DO
			END DO

		CASE(2)

			DO ta = 0, 1  ! Loop over isospin

				factor = 1.0D0 / (2.0D0 * m(ta))

				DO la = 0, Lmax  ! Loop over "bra" orbital angular momentum

							   nmaxi = Min(Nmax, NmaxOfL(la))
				 	IF (CompHO .EQ. 1) nmaxi = ((Lmax - la) / 2) + 1

					DO na = 1, nmaxi  ! Loop over "bra" main quantum number
						DO nc = 1, na  ! Loop over "bra" main quantum number
							EkField%d3tensor(la)%d2(na, nc) = factor * nabla2(na, nc, la)
							R2Field%d3tensor(la)%d2(na, nc) = r2_cuton(na, nc, la)
						END DO
					END DO

				END DO
			END DO

		END SELECT


                 ! Writing output: the one-body kinetic energy

		SELECT CASE (Basis)

		CASE(1)
			IF (N_0 < 10) THEN
				WRITE(filename, "(A,I1,A)") "data/Ek", N_0, "_HO.txt"
			ELSE
				WRITE(filename, "(A,I2,A)") "data/Ek", N_0, "_HO.txt"
			END IF
		CASE(2)
			IF (N_0 < 10) THEN
				WRITE(filename, "(A,I1,A)") "data/Ek", N_0, "_WS.txt"
			ELSE
				WRITE(filename, "(A,I2,A)") "data/Ek", N_0, "_WS.txt"
			END IF
		END SELECT

		OPEN(file_desc, FILE=filename, ACTION="WRITE", IOSTAT=file_error)

		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to write the results in file ", filename
		ELSE
			DO la = 0, Lmax
								     nmaxi = Min(Nmax, NmaxOfL(la))
				IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) nmaxi = ((N_0 - la) / 2) + 1

				DO na = 1, nmaxi
					DO nc = 1, na
						WRITE (file_desc, "(I3,I3,I3,E24.16)", IOSTAT=file_error) &
							la, na, nc, EkField%d3tensor(la)%d2(na, nc)
					END DO
				END DO
			END DO
			CLOSE(file_desc)
		END IF

                 ! Writing output: the one-body root mean square radius

		SELECT CASE (Basis)

		CASE(1)
			IF (N_0 < 10) THEN
				WRITE(filename, "(A,I1,A)") "data/R2", N_0, "_HO.txt"
			ELSE
				WRITE(filename, "(A,I2,A)") "data/R2", N_0, "_HO.txt"
			END IF
		CASE(2)
			IF (N_0 < 10) THEN
				WRITE(filename, "(A,I1,A)") "data/R2", N_0, "_WS.txt"
			ELSE
				WRITE(filename, "(A,I2,A)") "data/R2", N_0, "_WS.txt"
			END IF
		END SELECT

		OPEN(file_desc, FILE=filename, ACTION="WRITE", IOSTAT=file_error)

		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention:  Impossible to write the results in file ", filename
		ELSE
			DO la = 0, Lmax
								     nmaxi = Min(Nmax, NmaxOfL(la))
				IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) nmaxi = ((N_0 - la) / 2) + 1

				DO na = 1, nmaxi
					DO nc = 1, na
						WRITE (file_desc, "(I3,I3,I3,E24.16)", IOSTAT=file_error) &
							la, na, nc, R2Field%d3tensor(la)%d2(na, nc)
					END DO
				END DO
			END DO
			CLOSE(file_desc)
		END IF
		RETURN

	CONTAINS

                !
		!  Function calculating the matrix elements of the form
		!           < n'l' | nabla^2 | nl >
		!
		!  Refs: Sec. 5.1.1., Page 45,
		!        Appendix C, Page 115
		!

		FUNCTION nabla2HO(na, nc, la)
			DOUBLE PRECISION nabla2HO
			INTEGER, INTENT(IN) :: na, nc, la

			IF ((na + 1) .EQ. nc) THEN
				nabla2HO = sq(nc) * sq2(nc + la) /b_0**2
			ELSE IF (na .EQ. (nc + 1)) THEN
				nabla2HO = sq(na) * sq2(na + la) /b_0**2
			ELSE IF (na .eq. nc) THEN
				nabla2HO = (DBLE(2*na + la) + 1.5D0) /b_0**2
			ELSE
				nabla2HO = 0.0D0
			END IF

			RETURN
		END FUNCTION nabla2HO

                !
		!  Function calculating the matrix elements of the squared radius
		!           < n'l' | r^2 | nl >
		!
		!  Refs: Sec. 5.1.1., Page 45
		!

		FUNCTION r2_cutonHO(na, nc, la)
			DOUBLE PRECISION r2_cutonHO
			INTEGER, INTENT(IN) :: na, nc, la

			IF ((na + 1) .eq. nc) THEN
				r2_cutonHO = - sq(nc) * sq2(nc + la) *(b_0**2)
			ELSE IF (na .eq. (nc + 1)) THEN
				r2_cutonHO = - sq(na) * sq2(na + la) *(b_0**2)
			ELSE IF (na .eq. nc) THEN
				r2_cutonHO = (DBLE(2*na + la) + 1.5D0 ) *(b_0**2)
			ELSE
				r2_cutonHO = 0.0D0
			END IF

			RETURN
		END FUNCTION r2_cutonHO

	END SUBROUTINE Global_start

        !---------------------------------------------------------------------------------!
	!		                       		  			          !
	! Function giving the reduced matrix elements < na la || NABLA || nb lb > 	  !
	! for an arbitrary spherical basis (not necessarily the Harmonic Oscillator)      !
	!		                       		  			          !
        !---------------------------------------------------------------------------------!

	FUNCTION nabla2(na, nc, la)
		DOUBLE PRECISION nabla2
		INTEGER, INTENT(IN) :: na, nc, la
		INTEGER :: IndexBra, IndexKet, IndexLoop, n_loop, l_loop
		DOUBLE PRECISION :: sum

		IndexBra = IndexVecNL(na,la)
		IndexKet = IndexVecNL(nc,la)

		IF (IndexBra .EQ. 0 .OR. IndexKet .EQ. 0) THEN
		    	nabla2 = 0.0D0
		    	RETURN
		END IF

		sum = 0.0D0

		! The summation must go up to N_0 (complete basis). Even though, the matrix element
		! <lmax, nmax || Nabla || lmax, nmax> is wrong as it contains a term proportional to
		! lmax+1

		DO l_loop = 0, Lmax+1
                        DO n_loop = 1, Nunity

                                IndexLoop = IndexAux(n_loop, l_loop)

				IF (IndexLoop .NE. 0) THEN
					sum = sum + SymKineticEnergy2Body_nabla(na, la, n_loop, l_loop) &
						  * SymKineticEnergy2Body_nabla(n_loop, l_loop, nc, la)
				END IF

			END DO
		END DO

		nabla2 = sum / DBLE(2*la + 1)

		RETURN
	END FUNCTION nabla2

	FUNCTION r2_cuton(na, nc, la)
		DOUBLE PRECISION r2_cuton
		INTEGER, INTENT(IN) :: na, nc, la
		INTEGER :: IndexBra, IndexKet

		IndexBra = IndexVecNL(na,la)
		IndexKet = IndexVecNL(nc,la)

		IF (IndexBra .EQ. 0 .OR. IndexKet .EQ. 0) THEN
		    	r2_cuton = 0.0D0
		    	RETURN
		END IF

		r2_cuton = IntegralR(IndexBra,IndexKet)

		RETURN
	END FUNCTION r2_cuton

        !---------------------------------------------------------------------------------!
	!		                       		  			          !
	! Function giving the reduced matrix elements < na la || NABLA || nb lb > 	  !
	! for an arbitrary spherical basis (not necessarily the Harmonic Oscillator)      !
	!		                       		  			          !
        !---------------------------------------------------------------------------------!

	FUNCTION SymKineticEnergy2Body_nabla(na, la, nb, lb)
		DOUBLE PRECISION SymKineticEnergy2Body_nabla
		INTEGER, INTENT(IN) :: na, nb, la, lb
		INTEGER :: IndexBra, IndexKet

                IndexBra = IndexAux(na,la)
                IndexKet = IndexAux(nb,lb)

		IF (IndexBra .EQ. 0 .OR. IndexKet .EQ. 0) THEN
			SymKineticEnergy2Body_nabla = 0.0D0
			RETURN
		END IF

		IF (la .EQ. (lb + 1)) THEN
		   	SymKineticEnergy2Body_nabla = + IntegralA(IndexBra, IndexKet, lb)*sq(la)
		ELSE IF (lb .EQ. (la + 1)) THEN
		  	SymKineticEnergy2Body_nabla = - IntegralB(IndexBra, IndexKet, lb)*sq(lb)
		ELSE
		   	SymKineticEnergy2Body_nabla = 0.0D0
		END IF

		RETURN
	END FUNCTION SymKineticEnergy2Body_nabla

        !---------------------------------------------------------!
	!		                       		          !
	! Function giving the A-integral of Notes, Sec. 3.2.3     !
	!		                       		  	  !
        !---------------------------------------------------------!

	FUNCTION IntegralA(IndexBra, IndexKet, lb)
		DOUBLE PRECISION :: IntegralA, Result
		DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)
		INTEGER, INTENT(IN) :: IndexBra, IndexKet, lb
		INTEGER :: Ipoint

	       	ALLOCATE(Integrand(1:Npoint))

		DO Ipoint = 1, Npoint
			  Integrand(Ipoint) = WaveFun(Ipoint,IndexBra)*WaveDeri(Ipoint,IndexKet) - &
		 	        DBLE(lb + 1)* WaveFun(Ipoint,IndexBra)*WaveFun(Ipoint,IndexKet)/RadMesh(Ipoint)
		END DO

		CALL simps(Integrand,Npoint,MeshStep,Result)

		DEALLOCATE(Integrand)

		IntegralA = Result

		RETURN
	END FUNCTION IntegralA

        !---------------------------------------------------------!
	!		                       		          !
	! Function giving the B-integral of Notes, Sec. 3.2.3     !
	!		                       		  	  !
       	!---------------------------------------------------------!

	FUNCTION IntegralB(IndexBra, IndexKet, lb)
		DOUBLE PRECISION :: IntegralB, Result
		DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)
		INTEGER, INTENT(IN) :: IndexBra, IndexKet, lb
		INTEGER :: Ipoint

	       	ALLOCATE(Integrand(1:Npoint))

		DO Ipoint = 1, Npoint
		 	 Integrand(Ipoint) = WaveFun(Ipoint,IndexBra)*WaveDeri(Ipoint,IndexKet) + &
		     	           DBLE(lb)* WaveFun(Ipoint,IndexBra)*WaveFun(Ipoint,IndexKet)/RadMesh(Ipoint)
		END DO

		CALL simps(Integrand,Npoint,MeshStep,Result)

		DEALLOCATE(Integrand)

		IntegralB = Result

		RETURN
	END FUNCTION IntegralB

        !---------------------------------------------------------!
	!		                       		          !
	! Function giving the B-integral of Notes, Sec. 3.2.3     !
	!		                       		  	  !
       	!---------------------------------------------------------!

	FUNCTION IntegralR(IndexBra, IndexKet)
		DOUBLE PRECISION :: IntegralR, Result
		DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)
		INTEGER, INTENT(IN) :: IndexBra, IndexKet
		INTEGER :: Ipoint

	       	ALLOCATE(Integrand(1:Npoint))

		DO Ipoint = 1, Npoint
		 	 Integrand(Ipoint) = WaveFun(Ipoint,IndexBra)*WaveFun(Ipoint,IndexKet)*(RadMesh(Ipoint)**2)
		END DO

		CALL simps(Integrand,Npoint,MeshStep,Result)

		DEALLOCATE(Integrand)

		IntegralR = Result

		RETURN
	END FUNCTION IntegralR

	SUBROUTINE Global_del
		!TODO
	END SUBROUTINE Global_del

END MODULE global
