MODULE global

	USE input
	USE lgfactor
	USE symd3t
	USE symtalm
	USE gauss
	USE angmom

	IMPLICIT NONE

	! Definimos los valores que identificarán la separación entre los
	! parámetros de protones y neutrones almacenados en las tablas del núcleo 
	INTEGER, PARAMETER :: PROTON  = 0
	INTEGER, PARAMETER :: NEUTRON = 1

	INTEGER A
	INTEGER NLag

	INTEGER, DIMENSION(8) :: MagicNumber
	DATA MagicNumber/2, 14, 28, 50, 82, 126, 184, 256/

	DOUBLE PRECISION, DIMENSION(0:1), PARAMETER :: m = (0.024111868, 0.024111868)

	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: sq, sq2

	DOUBLE PRECISION, PARAMETER :: x0 = 1.0

	! Parametros de la fuerza de Gogny
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
	DATA Gogny_W0 / 130.0, 115.0, 130.0 / ! DSA1, D1, D1_PRIMA

	DOUBLE PRECISION, DIMENSION(3) :: Gogny_t0
	DATA Gogny_t0 / 1390.6, 1350.0, 1350.0 / ! DSA1, D1, D1_PRIMA

	TYPE (SymD3Tensor) EkField, R2Field
	TYPE (GaussLaguerreQuadrature) GaussLQ

CONTAINS

	SUBROUTINE Global_new

		INTEGER max_1, max_2

		! Leemos los parámetros de entrada para la ejecución del programa
		CALL Input_read
		PRINT "(/A)", "Valores de los parametros de la fuerza de Gogny:"
		PRINT "(A,F8.2)", "W(0) = ", Gogny_W(0, Gogny)
		PRINT "(A,F8.2)", "W(1) = ", Gogny_W(1, Gogny)
		PRINT "(A,F8.2)", "B(0) = ", Gogny_B(0, Gogny)
		PRINT "(A,F8.2)", "B(1) = ", Gogny_B(1, Gogny)
		PRINT "(A,F8.2)", "H(0) = ", Gogny_H(0, Gogny)
		PRINT "(A,F8.2)", "H(1) = ", Gogny_H(1, Gogny)
		PRINT "(A,F8.2)", "M(0) = ", Gogny_M(0, Gogny)
		PRINT "(A,F8.2)", "M(1) = ", Gogny_M(1, Gogny)
		PRINT "(A,F8.2)", "W0 = ", Gogny_W0(Gogny)
		PRINT "(A,F8.2)", "t0 = ", Gogny_t0(Gogny)

		NLag = MIN(N_0 * 32, 156) ! NLag = UMIN(N_0 << 5, 156)

		max_1 = 4 * N_0
		max_2 = max_1 + 1
!TODO El autor lo implementa, pero no se utiliza en ninguna parte del código
!		CALL LogFactorials_new(max_2)
!		CALL LogSemiFactorials_new(max_2)
		CALL DDLogFactorials_new(max_2)
		CALL DDLogSemiFactorials_new(max_2)

		CALL SquareRoot_new(max_1)
		CALL SemiSquareRoot_new(max_1)
		CALL Global_start
		CALL SymCoefficientB_new
		CALL SymKumar_new
		CALL GaussLaguerreQuadrature_new(GaussLQ, NLag, DBLE(0.5))
		CALL ThreeJSymbols_new
		RETURN

	CONTAINS

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

		SUBROUTINE SemiSquareRoot_new(imax)
			INTEGER, INTENT(IN) :: imax

			INTEGER i

			ALLOCATE(sq2(0:imax))
			sq2(0) = 0.7071067811865476
			DO i = 1, imax
				sq2(i) = SQRT(i + 0.5)
			END DO
			RETURN
		END SUBROUTINE SemiSquareRoot_new

	END SUBROUTINE Global_new

	SUBROUTINE Global_start

		INTEGER ta, la, na, nc
		DOUBLE PRECISION factor
		INTEGER d

		CHARACTER(64) filename
		INTEGER, PARAMETER :: file_desc = 6
		INTEGER file_error

		CALL SymD3Tensor_new(EkField)
		CALL SymD3Tensor_new(R2Field)
		DO ta = 0, 1
			factor = 1.0 / (2.0 * m(ta))
			DO la = 0, N_0
				d = ((N_0 - la) / 2) + 1
				DO na = 1, d
					EkField%d3tensor(la)%d2(na, na) = factor * nabla2(na - 1, na - 1, la)
					R2Field%d3tensor(la)%d2(na, na) = r2_cuton(na - 1, na - 1, la)
				END DO
				DO na = 1, d - 1
					EkField%d3tensor(la)%d2(na + 1, na) = factor * nabla2(na - 1, na, la)
					R2Field%d3tensor(la)%d2(na + 1, na) = r2_cuton(na - 1, na, la)
					DO nc = na + 2, d
						EkField%d3tensor(la)%d2(nc, na) = 0.0
						R2Field%d3tensor(la)%d2(nc, na) = 0.0
					END DO
				END DO
			END DO
		END DO

                 ! Writing output: the one-body kinetic energy

		IF (N_0 < 10) THEN
			WRITE(filename, "(A,I1,A)") "data/Ek", N_0, ".txt"
		ELSE
			WRITE(filename, "(A,I2,A)") "data/Ek", N_0, ".txt"
		END IF

		OPEN(file_desc, FILE=filename, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to write the results in file ", filename
		ELSE
			DO la = 0, N_0
				d = ((N_0 - la) / 2) + 1
				DO na = 1, d
					DO nc = 1, na
						WRITE (file_desc, "(I3,I3,I3,E)", IOSTAT=file_error) &
							la, na, nc, EkField%d3tensor(la)%d2(na, nc)
					END DO
				END DO
			END DO
			CLOSE(file_desc)
		END IF

                 ! Writing output: the one-body root mean square radius

		IF (N_0 < 10) THEN
			WRITE(filename, "(A,I1,A)") "data/R2", N_0, ".txt"
		ELSE
			WRITE(filename, "(A,I2,A)") "data/R2", N_0, ".txt"
		END IF

		OPEN(file_desc, FILE=filename, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention:  Impossible to write the results in file ", filename
		ELSE
			DO la = 0, N_0
				d = ((N_0 - la) / 2) + 1
				DO na = 1, d
					DO nc = 1, na
						WRITE (file_desc, "(I3,I3,I3,E)", IOSTAT=file_error) &
							la, na, nc, R2Field%d3tensor(la)%d2(na, nc)
					END DO
				END DO
			END DO
			CLOSE(file_desc)
		END IF
		RETURN

	CONTAINS

		FUNCTION nabla2(na, nc, la)
			DOUBLE PRECISION nabla2
			INTEGER, INTENT(IN) :: na, nc, la

			IF ((na + 1) .EQ. nc) THEN
				nabla2 = sq(nc) * sq2(nc + la)
			ELSE IF (na .EQ. (nc + 1)) THEN
				nabla2 = sq(na) * sq2(na + la)
			ELSE IF (na .eq. nc) THEN
				nabla2 = (na * 2.0) + la + 1.5
			ELSE
				nabla2 = 0.0
			END IF
			RETURN
		END FUNCTION nabla2

		FUNCTION r2_cuton(na, nc, la)
			DOUBLE PRECISION r2_cuton
			INTEGER, INTENT(IN) :: na, nc, la

			IF ((na + 1) .eq. nc) THEN
				r2_cuton = - sq(nc) * sq2(nc + la)
			ELSE IF (na .eq. (nc + 1)) THEN
				r2_cuton = - sq(na) * sq2(na + la)
			ELSE IF (na .eq. nc) THEN
				r2_cuton = (na * 2.0) + la + 1.5
			ELSE
				r2_cuton = 0.0
			END IF
			RETURN
		END FUNCTION r2_cuton

	END SUBROUTINE Global_start

	SUBROUTINE Global_del
		!TODO
	END SUBROUTINE Global_del

END MODULE global
