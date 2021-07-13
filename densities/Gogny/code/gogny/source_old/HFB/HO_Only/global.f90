 MODULE global

	USE input
	USE lgfactor
	USE symd3t
	USE symtalm
	USE gauss
	USE angmom

	IMPLICIT NONE

	! We define the flags that will identify each isospin throughout the program
	
	INTEGER, PARAMETER :: PROTON  = 0
	INTEGER, PARAMETER :: NEUTRON = 1

	INTEGER A
	INTEGER NLag

	! The list of spherical magic numbers
	
	INTEGER, DIMENSION(0:8) :: MagicNumber
	DATA MagicNumber/ 2, 8, 20, 28, 50, 82, 126, 184, 256/

	! m is equal to mc²/(hbar*c)²
	
	DOUBLE PRECISION, DIMENSION(0:1), PARAMETER :: m = (0.024111868, 0.024111868)
        !DOUBLE PRECISION, DIMENSION(0:1), PARAMETER :: m = (0.0240961393, 0.024129354)

	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: sq, sq2

	DOUBLE PRECISION, PARAMETER :: x0 = 1.0

	! Parameters of the Gogny force
	DOUBLE PRECISION, DIMENSION(0:1, 3) :: Gogny_W
	DOUBLE PRECISION, DIMENSION(0:1, 3) :: Gogny_B
	DOUBLE PRECISION, DIMENSION(0:1, 3) :: Gogny_H
	DOUBLE PRECISION, DIMENSION(0:1, 3) :: Gogny_M

	! GOGNYF = DSA1
	DATA Gogny_W / -1720.30,  103.639, -402.40, -21.30, -402.40, -21.30 /
	DATA Gogny_B /  1300.00, -163.483, -100.00, -11.77, -100.00, -11.77 /
	DATA Gogny_H / -1813.53,  162.812, -496.20,  37.27, -496.20,  37.27 /
	DATA Gogny_M /  1397.60, -223.934,  -23.56, -68.81,  -23.56, -68.81 /

	DOUBLE PRECISION, DIMENSION(3) :: Gogny_W0
	DATA Gogny_W0 / 130.0, 115.0, 130.0 / ! D1S, D1, D1prime

	DOUBLE PRECISION, DIMENSION(3) :: Gogny_t0
	DATA Gogny_t0 / 1390.6, 1350.0, 1350.0 / ! D1S, D1, D1prime

	TYPE (SymD3Tensor) EkField, R2Field
	TYPE (GaussLaguerreQuadrature) GaussLQ

 CONTAINS

	SUBROUTINE Global_new

		INTEGER max_1, max_2

		! Reading the input parameters
		CALL Input_read
		
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

!TODO El autor lo implementa, pero no se utiliza en ninguna parte del cÃ³digo
!		CALL LogFactorials_new(max_2)
!		CALL LogSemiFactorials_new(max_2)

                ! Defined in module "lgfactor.f90"
		CALL DDLogFactorials_new(max_2)
		CALL DDLogSemiFactorials_new(max_2)

                ! Defined here...
		CALL SquareRoot_new(max_1)
		CALL SemiSquareRoot_new(max_1)
		
                ! Calculate one-body kinetic energy and r.m.s. radius
		CALL Global_start
		
                ! Defined in module "symtalm.f90" (only of analytical HO basis is used)
		
		CALL SymCoefficientB_new
		CALL SymKumar_new
	
		CALL GaussLaguerreQuadrature_new(GaussLQ, NLag, DBLE(0.5))
		
                ! Defined in module "angmom.f90"
		CALL ThreeJSymbols_new
		
		RETURN

	CONTAINS

                ! Subroutine storing in a vector of size imax the square root of 
		! all integers from 1 to imax: sqrt(i)

		SUBROUTINE SquareRoot_new(imax)
			INTEGER, INTENT(IN) :: imax

			INTEGER i

			ALLOCATE(sq(0:imax))
			sq(0) = 0.0
			DO i = 1, imax
				sq(i) = SQRT(i + 0.0)
			END DO
			RETURN
		END SUBROUTINE SquareRoot_new

                ! Subroutine storing in a vector of size imax the square root of 
		! all integers from 1 to imax, plus one half: sqrt(i+0.5)

		SUBROUTINE SemiSquareRoot_new(imax)
			INTEGER, INTENT(IN) :: imax

			INTEGER i

			ALLOCATE(sq2(0:imax))
			sq2(0) = 0.5*SQRT(2.0)
			DO i = 1, imax
				sq2(i) = SQRT(i + 0.5)
			END DO
			RETURN
		END SUBROUTINE SemiSquareRoot_new

	END SUBROUTINE Global_new

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
		
		DO ta = 0, 1  ! Loop over isospin
		
			factor = 1.0 / (2.0 * m(ta))
			
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
						EkField%d3tensor(la)%d2(nc, na) = 0.0
						R2Field%d3tensor(la)%d2(nc, na) = 0.0
					END DO
					
				END DO
				
			END DO
		END DO

                ! Writing output: the one-body kinetic energy

		IF (N_0 < 10) THEN
			WRITE(filename, "(A,I1,A)") "data/Ek", N_0, "_HO.txt"
		ELSE
			WRITE(filename, "(A,I2,A)") "data/Ek", N_0, "_HO.txt"			
		END IF

		OPEN(file_desc, FILE=filename, ACTION="WRITE", IOSTAT=file_error)
		
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to write the results in file ", filename
		ELSE
			DO la = 0, N_0

				nmaxi = ((N_0 - la) / 2) + 1
				
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

		IF (N_0 < 10) THEN
			WRITE(filename, "(A,I1,A)") "data/R2", N_0, "_HO.txt"
		ELSE
			WRITE(filename, "(A,I2,A)") "data/R2", N_0, "_HO.txt"			
		END IF

		OPEN(file_desc, FILE=filename, ACTION="WRITE", IOSTAT=file_error)
		
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention:  Impossible to write the results in file ", filename
		ELSE
			DO la = 0, N_0

				nmaxi = ((N_0 - la) / 2) + 1
				
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
				nabla2HO = ((na * 2.0) + la + 1.5) /b_0**2
			ELSE
				nabla2HO = 0.0
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
				r2_cutonHO = ((na * 2.0) + la + 1.5 ) *(b_0**2)
			ELSE
				r2_cutonHO = 0.0
			END IF
			
			RETURN
		END FUNCTION r2_cutonHO

	END SUBROUTINE Global_start

	SUBROUTINE Global_del
		!TODO
	END SUBROUTINE Global_del

END MODULE global
