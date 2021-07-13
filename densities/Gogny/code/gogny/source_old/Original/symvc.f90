MODULE symvc

	USE input
	USE angmom
	USE symd3t
	USE symfield
	USE ic

	IMPLICIT NONE

	TYPE SymVCph
		TYPE (SymD3Tensor_SymD3Tensor) v_local, v1_exch, v2_exch, v1
		CHARACTER(LEN = 64) filename
	END TYPE

	TYPE SymVCpp
		TYPE (SymD3Tensor_SymD3Tensor) v1_pair, v2_pair
		CHARACTER(LEN = 64) filename
	END TYPE

!	INTERFACE OPERATOR(*)
!		MODULE PROCEDURE SymVCph_product, SymVCpp_product
!	END INTERFACE

CONTAINS

	SUBROUTINE SymVCph_new(vCph)
		TYPE (SymVCph), INTENT(INOUT) :: vCph

		CALL SymD3Tensor_SymD3Tensor_new(vCph%v_local)
		CALL SymD3Tensor_SymD3Tensor_new(vCph%v1_exch)
		CALL SymD3Tensor_SymD3Tensor_new(vCph%v2_exch)
		CALL SymD3Tensor_SymD3Tensor_new(vCph%v1)

		IF (N_0 < 10) THEN
			WRITE(vCph%filename, "(A,I1,A)") "data/vC", N_0, "ph.txt"
		ELSE
			WRITE(vCph%filename, "(A,I2,A)") "data/vC", N_0, "ph.txt"
		END IF
		RETURN
	END SUBROUTINE SymVCph_new

	SUBROUTINE SymVCph_calculate(vCph)
		TYPE (SymVCph), INTENT(INOUT) :: vCph

		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		DOUBLE PRECISION sumi, total_iC
		INTEGER k, kmin, kmax
		DOUBLE PRECISION tres_j_cuad, cuad
		DOUBLE PRECISION sum1, sum2

		INTEGER, PARAMETER :: file_desc = 6
		INTEGER file_error

		IF (N_0 .LE. 16) THEN
			IF (SymVCph_read(vCph)) THEN
				CALL SymD3Tensor_SymD3Tensor_add(vCph%v1, vCph%v_local, vCph%v1_exch)
				RETURN
			END IF

			OPEN(file_desc, FILE=vCph%filename, ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				PRINT *, "*** AVISO: No se pueden escribir los resultados en ", vCph%filename
			END IF
		END IF

		PRINT *
		PRINT *, "Se van a calcular los elementos de matriz ph de C:"
		DO la = 0, N_0
			namax = ((N_0 - la) / 2) + 1
			DO na = 1, namax
				DO nc = 1, na
					DO lb = 0, la
						nbmax = ((N_0 - lb) / 2) + 1
						DO nb = 1, nbmax
							DO nd = 1, nb
								total_iC = ICoulomb(na - 1, la, nb - 1, lb, nc - 1, la, nd - 1, lb, 0)
								IF (nb .NE. nd) THEN
									total_iC = total_iC + ICoulomb(na - 1, la, nd - 1, lb, nc - 1, la, nb - 1, lb, 0)
									total_iC = total_iC / 2.0
								END IF
								sumi = (VC * I_4PI) * total_iC
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v_local, la, na, nc, lb, nb, nd, sumi)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v_local, lb, nb, nd, la, na, nc, sumi)

								kmin = ABS(la - lb)
								kmax = la + lb
								sum1 = 0.0
								sum2 = 0.0
								DO k = kmin, kmax, 2
									tres_j_cuad = (2 * k + 1) * (ThreeJSymbols_get(2 * la, 2 * k, 2 * lb) ** 2)
									cuad = CUAD2(la, lb, k)
									total_iC = ICoulomb(na - 1, la, nb - 1, lb, nd - 1, lb, nc - 1, la, k)
									IF(nb .NE. nd) THEN
										total_iC = total_iC + ICoulomb(na - 1, la, nd - 1, lb, nb - 1, lb, nc - 1, la, k)
										total_iC = total_iC / 2.0
									END IF
									total_iC = total_iC * (VC * I_4PI)
									sum1 = sum1 + (tres_j_cuad * total_iC)
									sum2 = sum2 + (cuad * tres_j_cuad * total_iC)
								END DO

								sum1 = -0.5 * sum1
								sum2 = -sum2
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v1_exch, la, na, nc, lb, nb, nd, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v1_exch, lb, nb, nd, la, na, nc, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v2_exch, la, na, nc, lb, nb, nd, sum2)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v2_exch, lb, nb, nd, la, na, nc, sum2)
								IF ((N_0 .LE. 16) .AND. (file_error .EQ. 0)) THEN
									WRITE (file_desc, FMT="(I3,I3,I3,I3,I3,I3,E,E,E)", IOSTAT=file_error) &
										la, na, nc, lb, nb, nd, sumi, sum1, sum2
								END IF
							END DO
						END DO
					END DO
				END DO
			END DO
			PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
		END DO
		CALL SymD3Tensor_SymD3Tensor_add(vCph%v1, vCph%v_local, vCph%v1_exch)
		IF (N_0 .LE. 16) CLOSE(file_desc)
		PRINT "(A,A)", "Resultado almacenado en: ", vCph%filename
		RETURN
	END SUBROUTINE SymVCph_calculate

	SUBROUTINE SymVCph_product(HF_out, vCph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVCph), INTENT(IN) :: vCph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		CALL SymVCph_get_Gamma(HF_out, vCph, HF_in)
		RETURN
	END SUBROUTINE SymVCph_product

	SUBROUTINE SymVCph_get_Gamma(HF_out, vCph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVCph), INTENT(IN) :: vCph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(0), vCph%v1, HF_in%p(0))
		HF_out%p(1) = DBLE(0.0)
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(0), vCph%v2_exch, HF_in%a(0))
		HF_out%a(1) = DBLE(0.0)
		RETURN
	END SUBROUTINE SymVCph_get_Gamma

	SUBROUTINE SymVCph_get_LocalGamma(HF_out, vCph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVCph), INTENT(IN) :: vCph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(0), vCph%v_local, HF_in%p(0))
		HF_out%p(1) = DBLE(0.0)
		HF_out%a(0) = DBLE(0.0)
		HF_out%a(1) = DBLE(0.0)
		RETURN
	END SUBROUTINE SymVCph_get_LocalGamma

	SUBROUTINE SymVCph_get_ExchangeGamma(HF_out, vCph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVCph), INTENT(IN) :: vCph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(0), vCph%v1_exch, HF_in%p(0))
		HF_out%p(1) = DBLE(0.0)
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(0), vCph%v2_exch, HF_in%a(0))
		HF_out%a(1) = DBLE(0.0)
		RETURN
	END SUBROUTINE SymVCph_get_ExchangeGamma

	FUNCTION SymVCph_read(vCph)
		LOGICAL SymVCph_read
		TYPE (SymVCph), INTENT(INOUT) :: vCph

		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER i1, i2, i3, i4, i5, i6
		INTEGER, PARAMETER :: file_desc = 6
		INTEGER file_error
		DOUBLE PRECISION sumi, sum1, sum2

		OPEN(file_desc, FILE=vCph%filename, ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "No se pudo leer el archivo: ", vCph%filename
			SymVCph_read = .FALSE.
			RETURN
		END IF

!		PRINT *, "Leyendo:", vCph%filename
		DO la = 0, N_0
			namax = ((N_0 - la) / 2) + 1
			DO na = 1, namax
				DO nc = 1, na
					DO lb = 0, la
						nbmax = ((N_0 - lb) / 2) + 1
						DO nb = 1, nbmax
							DO nd = 1, nb
								READ (file_desc, FMT="(I3,I3,I3,I3,I3,I3,E,E,E)", IOSTAT=file_error) &
									i1, i2, i3, i4, i5, i6, sumi, sum1, sum2
								IF ((file_error .NE. 0) .OR. &
									(la .NE. i1) .OR. (na .NE. i2) .OR. (nc .NE. i3) .OR. &
									(lb .NE. i4) .OR. (nb .NE. i5) .OR. (nd .NE. i6)) THEN
									PRINT *, "Informacion no validad en el archivo: ", vCph%filename
									CLOSE(file_desc)
									SymVCph_read = .FALSE.
									RETURN
								END IF
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v_local, la, na, nc, lb, nb, nd, sumi)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v_local, lb, nb, nd, la, na, nc, sumi)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v1_exch, la, na, nc, lb, nb, nd, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v1_exch, lb, nb, nd, la, na, nc, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v2_exch, la, na, nc, lb, nb, nd, sum2)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v2_exch, lb, nb, nd, la, na, nc, sum2)
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
		CLOSE(file_desc)
		SymVCph_read = .TRUE.
		RETURN
	END FUNCTION SymVCph_read

	SUBROUTINE SymVCph_del(vCph)
		TYPE (SymVCph), INTENT(INOUT) :: vCph

		CALL SymD3Tensor_SymD3Tensor_del(vCph%v_local)
		CALL SymD3Tensor_SymD3Tensor_del(vCph%v1_exch)
		CALL SymD3Tensor_SymD3Tensor_del(vCph%v2_exch)
		CALL SymD3Tensor_SymD3Tensor_del(vCph%v1)
		RETURN
	END SUBROUTINE SymVCph_del

	SUBROUTINE SymVCpp_new(vCpp)
		TYPE (SymVCpp), INTENT(INOUT) :: vCpp

		CALL SymD3Tensor_SymD3Tensor_new(vCpp%v1_pair)
		CALL SymD3Tensor_SymD3Tensor_new(vCpp%v2_pair)

		IF (N_0 < 10) THEN
			WRITE(vCpp%filename, "(A,I1,A)") "data/vC", N_0, "pp.txt"
		ELSE
			WRITE(vCpp%filename, "(A,I2,A)") "data/vC", N_0, "pp.txt"
		END IF
		RETURN
	END SUBROUTINE SymVCpp_new

	SUBROUTINE SymVCpp_calculate(vCpp)
		TYPE (SymVCpp), INTENT(INOUT) :: vCpp

		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER k, kmin, kmax
		DOUBLE PRECISION sum1, sum2, total_iC
		DOUBLE PRECISION tres_j_cuad, cuad

		INTEGER, PARAMETER :: file_desc = 6
		INTEGER file_error

		IF (N_0 .LE. 16) THEN
			IF (SymVCpp_read(vCpp)) RETURN

			OPEN(file_desc, FILE=vCpp%filename, ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				PRINT *, "*** AVISO: No se pueden escribir los resultados en ", vCpp%filename
			END IF
		END IF

		PRINT *, "Se van a calcular los elementos de matriz pp de C:"
		DO la = 0, N_0
			namax = ((N_0 - la) / 2) + 1
			DO na = 1, namax
				DO nc = 1, na
					DO lb = 0, la
						nbmax = ((N_0 - lb) / 2) + 1
						DO nb = 1, nbmax
							DO nd = 1, nb
								kmin = ABS(la - lb)
								kmax = la + lb
								sum1 = 0.0
								sum2 = 0.0
								DO k = kmin, kmax, 2
									tres_j_cuad = (2 * k + 1) * (ThreeJSymbols_get(2 * la, 2 * k, 2 * lb) ** 2)
									cuad = CUAD2(la, lb, k)
									total_iC = (VC * I_4PI) * ICoulomb(na - 1, la, nc - 1, la, nb - 1, lb, nd - 1, lb, k)
									IF(nb .NE. nd) THEN
										total_iC = total_iC + (VC * I_4PI) * &
											ICoulomb(na - 1, la, nc - 1, la, nd - 1, la, nb - 1, lb, 0)
										total_iC = total_iC / 2.0
									END IF
									sum1 = sum1 + (       tres_j_cuad * total_iC)
									sum2 = sum2 + (cuad * tres_j_cuad * total_iC)
								END DO
								sum1 = 0.5 * PAR(la + lb) * sum1
								sum2 =       PAR(la + lb) * sum2
								CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v1_pair, la, na, nc, lb, nb, nd, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v1_pair, lb, nb, nd, la, na, nc, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v2_pair, la, na, nc, lb, nb, nd, sum2)
								CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v2_pair, lb, nb, nd, la, na, nc, sum2)
								IF ((N_0 .LE. 16) .AND. (file_error .EQ. 0)) THEN
									WRITE (file_desc, "(I3,I3,I3,I3,I3,I3,E,E)", IOSTAT=file_error) &
										la, na, nc, lb, nb, nd, sum1, sum2
								END IF
							END DO
						END DO
					END DO
				END DO
			END DO
			PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
		END DO
		IF (N_0 .LE. 16) CLOSE(file_desc)
		RETURN
	END SUBROUTINE SymVCpp_calculate

	SUBROUTINE SymVCpp_product(HF_out, vCpp, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVCpp), INTENT(IN) :: vCpp
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		CALL SymVCpp_get_Delta(HF_out, vCpp, HF_in)
		RETURN
	END SUBROUTINE SymVCpp_product

	SUBROUTINE SymVCpp_get_Delta(HF_out, vCpp, P_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVCpp), INTENT(IN) :: vCpp
		TYPE (SymHartreeFockField), INTENT(IN) :: P_in

		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(0), vCpp%v1_pair, P_in%p(0))
		HF_out%p(1) = DBLE(0.0)
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(0), vCpp%v2_pair, P_in%a(0))
		HF_out%a(1) = DBLE(0.0)
		RETURN
	END SUBROUTINE SymVCpp_get_Delta

	FUNCTION SymVCpp_read(vCpp)
		LOGICAL SymVCpp_read
		TYPE (SymVCpp), INTENT(INOUT) :: vCpp

		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER i1, i2, i3, i4, i5, i6
		INTEGER, PARAMETER :: file_desc = 6
		INTEGER file_error
		DOUBLE PRECISION sum1, sum2

		OPEN(file_desc, FILE=vCpp%filename, ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "No se pudo leer el archivo: ", vCpp%filename
			SymVCpp_read = .FALSE.
			RETURN
		END IF

!		PRINT *, "Leyendo:", vCpp%filename
		DO la = 0, N_0
			namax = ((N_0 - la) / 2) + 1
			DO na = 1, namax
				DO nc = 1, na
					DO lb = 0, la
						nbmax = ((N_0 - lb) / 2) + 1
						DO nb = 1, nbmax
							DO nd = 1, nb
								READ (file_desc, FMT="(I3,I3,I3,I3,I3,I3,E,E)", IOSTAT=file_error) &
									i1, i2, i3, i4, i5, i6, sum1, sum2
								IF ((file_error .NE. 0) .OR. &
									(la .NE. i1) .OR. (na .NE. i2) .OR. (nc .NE. i3) .OR. &
									(lb .NE. i4) .OR. (nb .NE. i5) .OR. (nd .NE. i6)) THEN
									PRINT *, "Informacion no validad en el archivo: ", vCpp%filename
									CLOSE(file_desc)
									SymVCpp_read = .FALSE.
									RETURN
								END IF
								CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v1_pair, la, na, nc, lb, nb, nd, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v1_pair, lb, nb, nd, la, na, nc, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v2_pair, la, na, nc, lb, nb, nd, sum2)
								CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v2_pair, lb, nb, nd, la, na, nc, sum2)
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
		CLOSE(file_desc)

		SymVCpp_read = .TRUE.
		RETURN
	END FUNCTION SymVCpp_read

	SUBROUTINE SymVCpp_del(vCpp)
		TYPE (SymVCpp), INTENT(INOUT) :: vCpp

		CALL SymD3Tensor_SymD3Tensor_del(vCpp%v1_pair)
		CALL SymD3Tensor_SymD3Tensor_del(vCpp%v2_pair)
		RETURN
	END SUBROUTINE SymVCpp_del

END MODULE symvc
