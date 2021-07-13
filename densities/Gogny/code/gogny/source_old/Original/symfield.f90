MODULE symfield

	USE input
	USE global
	USE symd3t
!	USE symgden
	USE symtalm

	IMPLICIT NONE

	TYPE SymHartreeFockField
		TYPE (SymD3Tensor), DIMENSION(0:1) :: p
		TYPE (SymD3Tensor), DIMENSION(0:1) :: a
	END TYPE

	TYPE SymHartreeFockBogolField
		TYPE (SymHartreeFockField) rho
		TYPE (SymHartreeFockField) kap
	END TYPE

	INTERFACE ASSIGNMENT(=)
		MODULE PROCEDURE &
			SymHartreeFockField_assign, &
			SymHartreeFockBogolField_assign
	END INTERFACE

!	INTERFACE OPERATOR(+)
!		MODULE PROCEDURE SymHartreeFockField_add
!	END INTERFACE

	INTERFACE OPERATOR(*)
		MODULE PROCEDURE SymHartreeFockField_product2
	END INTERFACE

CONTAINS

	SUBROUTINE SymHartreeFockBogolField_new(HFB)
		TYPE (SymHartreeFockBogolField), INTENT(INOUT) :: HFB

		CALL SymHartreeFockField_new(HFB%rho)
		CALL SymHartreeFockField_new(HFB%kap)
		RETURN
	END SUBROUTINE SymHartreeFockBogolField_new

	FUNCTION SymHartreeFockBogolField_read(HFB, filename)
		LOGICAL SymHartreeFockBogolField_read
		TYPE (SymHartreeFockBogolField), INTENT(INOUT) :: HFB
		CHARACTER(*), INTENT(IN) :: filename

		INTEGER, PARAMETER :: file_desc = 6 ! Descriptor de fichero
		INTEGER file_error, new_N_0

		OPEN (file_desc, FILE=TRIM(filename), ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "No se pudo leer ", filename, " del disco"
			SymHartreeFockBogolField_read = .FALSE.
			RETURN
		END IF

		READ (file_desc, *) new_N_0
		PRINT *, "La funcion de onda es para un espacio de dimension ", new_N_0
		IF (new_N_0 .NE. N_0) THEN
			PRINT *, "Asi que no se puede usar"
			CLOSE (file_desc)
			SymHartreeFockBogolField_read = .FALSE.
			RETURN
		END IF

		IF ((SymHartreeFockField_read(HFB%rho, file_desc, file_error)) .AND. &
		    (SymHartreeFockField_read(HFB%kap, file_desc, file_error))) THEN
			SymHartreeFockBogolField_read = .TRUE.
		ELSE
			SymHartreeFockBogolField_read = .FALSE.
		END IF
		CLOSE (file_desc)
		RETURN
	END FUNCTION SymHartreeFockBogolField_read

	! Almacena los datos de los campos rho y kap en un fichero
	SUBROUTINE SymHartreeFockBogolField_write(HBF, filename)
		TYPE (SymHartreeFockBogolField), INTENT(IN) :: HBF
		CHARACTER(*), INTENT(IN) :: filename

		INTEGER, PARAMETER :: file_desc = 6 ! Descriptor de fichero
		INTEGER file_error

		OPEN (file_desc, FILE=TRIM(filename), ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) RETURN
		! El primer valor almacenado en el fichero es el número de capas: N_0
		WRITE (file_desc, FMT="(I)", IOSTAT=file_error) N_0
		! A continuación se almacenan rho y kap
		CALL SymHartreeFockField_write(HBF%rho, file_desc, file_error)
		CALL SymHartreeFockField_write(HBF%kap, file_desc, file_error)
		CLOSE (file_desc)
		RETURN
	END SUBROUTINE SymHartreeFockBogolField_write

	SUBROUTINE SymHartreeFockBogolField_assign(HBF1, HBF2)
		TYPE (SymHartreeFockBogolField), INTENT(INOUT) :: HBF1
		TYPE (SymHartreeFockBogolField), INTENT(IN) :: HBF2

		HBF1%rho = HBF2%rho
		HBF1%kap = HBF2%kap
		RETURN
	END SUBROUTINE SymHartreeFockBogolField_assign

	SUBROUTINE SymHartreeFockBogolField_del(HFB)
		TYPE (SymHartreeFockBogolField), INTENT(INOUT) :: HFB

		CALL SymHartreeFockField_del(HFB%rho)
		CALL SymHartreeFockField_del(HFB%kap)
		RETURN
	END SUBROUTINE SymHartreeFockBogolField_del

	SUBROUTINE SymHartreeFockField_new(HF)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF

		INTEGER ta

		DO ta = 0, 1
			CALL SymD3Tensor_new(HF%p(ta))
			CALL SymD3Tensor_new(HF%a(ta))
		END DO
		RETURN
	END SUBROUTINE SymHartreeFockField_new

	FUNCTION SymHartreeFockField_read(HF, file_desc, file_error)
		LOGICAL SymHartreeFockField_read
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF
		INTEGER, INTENT(IN) :: file_desc
		INTEGER, INTENT(INOUT) :: file_error

		INTEGER ta, length

		DO ta = 0, 1
			IF (.NOT. SymD3Tensor_read(HF%p(ta), file_desc, file_error)) THEN
				SymHartreeFockField_read = .FALSE.
				RETURN
			END IF
			IF (.NOT. SymD3Tensor_read(HF%a(ta), file_desc, file_error)) THEN
				SymHartreeFockField_read = .FALSE.
				RETURN
			END IF
		END DO
		SymHartreeFockField_read = .TRUE.
		RETURN
	END FUNCTION SymHartreeFockField_read

	SUBROUTINE SymHartreeFockField_write(HF, file_desc, file_error)
		TYPE (SymHartreeFockField), INTENT(IN) :: HF
		INTEGER, INTENT(IN) :: file_desc
		INTEGER, INTENT(INOUT) :: file_error

		INTEGER ta

		DO ta = 0, 1
			CALL SymD3Tensor_write(HF%p(ta), file_desc, file_error)
			CALL SymD3Tensor_write(HF%a(ta), file_desc, file_error)
		END DO
		RETURN
	END SUBROUTINE SymHartreeFockField_write

	SUBROUTINE SymHartreeFockField_assign(HF_out, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		INTEGER ta

		DO ta = 0, 1
			HF_out%p(ta) = HF_in%p(ta)
			HF_out%a(ta) = HF_in%a(ta)
		END DO
		RETURN
	END SUBROUTINE SymHartreeFockField_assign

	SUBROUTINE SymHartreeFockField_add(HF_out, HF1_in, HF2_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymHartreeFockField), INTENT(IN) :: HF1_in, HF2_in

		CALL SymD3Tensor_add(HF_out%p(0), HF1_in%p(0), HF2_in%p(0))
		CALL SymD3Tensor_add(HF_out%p(1), HF1_in%p(1), HF2_in%p(1))
		CALL SymD3Tensor_add(HF_out%a(0), HF1_in%a(0), HF2_in%a(0))
		CALL SymD3Tensor_add(HF_out%a(1), HF1_in%a(1), HF2_in%a(1))
		RETURN
	END SUBROUTINE SymHartreeFockField_add

	SUBROUTINE SymHartreeFockField_add_SymD3Tensor(HF_out, t_in, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymD3Tensor), INTENT(IN) :: t_in
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		CALL SymD3Tensor_add(HF_out%p(0), t_in, HF_in%p(0))
		CALL SymD3Tensor_add(HF_out%p(1), t_in, HF_in%p(1))
		HF_out%a(0) = HF_in%a(0)
		HF_out%a(1) = HF_in%a(1)
		RETURN
	END SUBROUTINE SymHartreeFockField_add_SymD3Tensor

	SUBROUTINE SymHartreeFockField_product(HT_out, R1, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HT_out
		DOUBLE PRECISION, INTENT(IN) :: R1
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		CALL SymD3Tensor_product(HT_out%p(0), R1, HF_in%p(0))
		CALL SymD3Tensor_product(HT_out%p(1), R1, HF_in%p(1))
		CALL SymD3Tensor_product(HT_out%a(0), R1, HF_in%a(0))
		CALL SymD3Tensor_product(HT_out%a(1), R1, HF_in%a(1))
		RETURN
	END SUBROUTINE SymHartreeFockField_product

	FUNCTION SymHartreeFockField_product2(HF1_in, HF2_in)
		DOUBLE PRECISION SymHartreeFockField_product2
		TYPE (SymHartreeFockField), INTENT(IN) :: HF1_in, HF2_in

		DOUBLE PRECISION, DIMENSION(0:1) :: sum2
		INTEGER ta

		DO ta = 0, 1
			sum2(ta) = (HF1_in%p(ta) * HF2_in%p(ta)) &
			         + (HF1_in%a(ta) * HF2_in%a(ta))
		END DO
		SymHartreeFockField_product2 = sum2(0) + sum2(1)
		RETURN
	END FUNCTION SymHartreeFockField_product2

	FUNCTION SymHartreeFockField_distance(HF1_in, HF2_in)
		DOUBLE PRECISION SymHartreeFockField_distance
		TYPE (SymHartreeFockField), INTENT(IN) :: HF1_in, HF2_in

		SymHartreeFockField_distance = MAX( &
			SymD3Tensor_distance(HF1_in%p(0), HF2_in%p(0)), &
			SymD3Tensor_distance(HF1_in%p(1), HF2_in%p(1)))
		RETURN
	END FUNCTION SymHartreeFockField_distance

!	SUBROUTINE SymHartreeFockField_inject(HF_out, HF2)
!		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
!		TYPE (SymHartreeFockField), INTENT(IN) :: HF2
!
!		INTEGER ta
!
!		DO ta=0, 1
!			HF_out%p(ta) = HF2%p(ta)
!			HF_out%a(ta) = HF2%a(ta)
!		END DO
!		IF (HF2%temporal) CALL SymHartreeFockField_del(HF2)
!	END SUBROUTINE SymHartreeFockField_inject

	FUNCTION SymHartreeFockField_ChargeDensity(HF, r)
		DOUBLE PRECISION SymHartreeFockField_ChargeDensity
		TYPE (SymHartreeFockField), INTENT(IN) :: HF
		DOUBLE PRECISION, INTENT(IN) :: r

		DOUBLE PRECISION x, d1, xk, sum1, sumlb
		INTEGER k, lb, nb, nd, nbmax, p2max, p2

		x = r * r
		xk = 1.0
		sum1 = 0.0
		DO k = 0, N_0
 			sumlb = 0.0
			DO lb = 0, k
				nbmax = ((N_0 - lb) / 2) + 1
				DO nb = 1, nbmax
					DO nd = 1, nbmax
						p2max = nb + nd - 2
						p2 = k - lb
						IF (p2 .GT. p2max) CYCLE
						d1 = SymCoefficientB_get(nb - 1, lb, nd - 1, lb, p2)
						sumlb = sumlb + (d1 * HF%p(PROTON)%d3tensor(lb)%d2(nd, nb))
					END DO
				END DO
			END DO
			sum1 = sum1 + (xk * sumlb)
			xk = xk * x
		END DO
		SymHartreeFockField_ChargeDensity = EXP(-x) * sum1
		RETURN
	END FUNCTION SymHartreeFockField_ChargeDensity

	SUBROUTINE SymHartreeFockField_del(HF)
		TYPE (SymHartreeFockField), INTENT(IN) :: HF

		INTEGER ta

		DO ta = 0, 1
			CALL SymD3Tensor_del(HF%p(ta))
			CALL SymD3Tensor_del(HF%a(ta))
		END DO
		RETURN
	END SUBROUTINE SymHartreeFockField_del

END MODULE symfield
