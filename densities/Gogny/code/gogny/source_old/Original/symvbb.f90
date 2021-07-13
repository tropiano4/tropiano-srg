MODULE symvbb

	USE input
	USE angmom
	USE symd3t
	USE symfield
	USE ibb

	IMPLICIT NONE

	TYPE SymVBBph
		DOUBLE PRECISION b1
		TYPE (SymD3Tensor_SymD3Tensor) &
			v_local_same_part, v_local_diff_part, &
			v1_exch_same_part, v1_exch_diff_part, &
			v2_exch_same_part, v2_exch_diff_part
		TYPE (SymD3Tensor_SymD3Tensor) &
			v1_same_part, v1_diff_part, &
			v2_same_part, v2_diff_part
		CHARACTER(LEN = 64) filename
	END TYPE

	TYPE SymVBBpp
		DOUBLE PRECISION b1
		TYPE (SymD3Tensor_SymD3Tensor) v1_pair, v2_pair
		CHARACTER(LEN = 64) filename
	END TYPE

!	INTERFACE OPERATOR(*)
!		MODULE PROCEDURE SymVBBph_product, SymVBBpp_product
!	END INTERFACE

	! Alcance de la fuerza de Brink-Booker
	DOUBLE PRECISION, DIMENSION(0:1) :: mu
	DATA mu /0.7, 1.2/

CONTAINS

	SUBROUTINE SymVBBph_new(vBBph, b)
		TYPE (SymVBBph), INTENT(INOUT) :: vBBph
		DOUBLE PRECISION, INTENT(IN) :: b

		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v_local_same_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v_local_diff_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v1_exch_same_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v1_exch_diff_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v2_exch_same_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v2_exch_diff_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v1_same_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v1_diff_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v2_same_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v2_diff_part)
		vBBph%b1 = b
		IF (N_0 < 10) THEN
			WRITE(vBBph%filename, "(A,I1,A)") "data/vBB", N_0, "ph.txt"
		ELSE
			WRITE(vBBph%filename, "(A,I2,A)") "data/vBB", N_0, "ph.txt"
		END IF
		RETURN
	END SUBROUTINE SymVBBph_new

	SUBROUTINE SymVBBph_calculate(vBBph)
		TYPE (SymVBBph), INTENT(INOUT) :: vBBph

		INTEGER i
		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		DOUBLE PRECISION sumi_sp, sumi_dp, total_ibb
		INTEGER k, kmin, kmax
		DOUBLE PRECISION sum1_sp, sum1_dp, sum2_sp, sum2_dp
		DOUBLE PRECISION tres_j_cuad, cuad

		DOUBLE PRECISION, DIMENSION(0:1) :: x, &
			coef_local_same_part, coef_local_diff_part, &
			coef_exch1_same_part, coef_exch1_diff_part, &
			coef_exch2_same_part, coef_exch2_diff_part

		INTEGER, PARAMETER :: file_desc = 6
		INTEGER file_error

		IF (N_0 .LE. 16) THEN
			IF (SymVBBph_read(vBBph)) THEN
				CALL SymVBBph_update(vBBph)
				RETURN
			END IF

			OPEN(file_desc, FILE=vBBph%filename, ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				PRINT *, "*** AVISO: No se pueden escribir los resultados en ", vBBph%filename
			ELSE
				WRITE (file_desc, FMT="(E)", IOSTAT=file_error) vBBph%b1
			END IF
		END IF

		DO i = 0, 1
			x(i) = mu(i) / vBBph%b1
			coef_local_diff_part(i) = I_4PI * (Gogny_W(i, Gogny) + (0.5 * Gogny_B(i, Gogny)))
			coef_local_same_part(i) = coef_local_diff_part(i) - (I_4PI * (Gogny_H(i, Gogny) + (0.5 * Gogny_M(i, Gogny))))
			coef_exch1_diff_part(i) = - I_4PI * ((0.5 * Gogny_H(i, Gogny)) + Gogny_M(i, Gogny))
			coef_exch1_same_part(i) = coef_exch1_diff_part(i) + (I_4PI * ((0.5 * Gogny_W(i, Gogny)) + Gogny_B(i, Gogny)))
			coef_exch2_diff_part(i) = - I_4PI * Gogny_H(i, Gogny)
			coef_exch2_same_part(i) = coef_exch2_diff_part(i) + (I_4PI * Gogny_W(i, Gogny))
		END DO

		PRINT *, "Se van a calcular los elementos de matriz ph de BB:"
		! Los mas altos son los mas problematicos
		! DO la = N_0, 0, -1
		DO la = 0, N_0
			namax = ((N_0 - la) / 2) + 1
			DO na = 1, namax
				DO nc = 1, na
					DO lb = 0, la
						nbmax = ((N_0 - lb) / 2) + 1
						DO nb = 1, nbmax
							DO nd = 1, nb
								sumi_sp = 0.0
								sumi_dp = 0.0
								DO i = 0, 1
									total_ibb = IBrinkBooker(na - 1, la, nb - 1, lb, nc - 1, la, nd - 1, lb, 0, x(i))
									IF (nb .NE. nd) THEN
										total_ibb = total_ibb + IBrinkBooker(na - 1, la, nd - 1, lb, nc - 1, la, nb - 1, lb, 0, x(i))
										total_ibb = total_ibb / 2.0
									END IF
									sumi_sp = sumi_sp + (coef_local_same_part(i) * total_ibb)
									sumi_dp = sumi_dp + (coef_local_diff_part(i) * total_ibb)
								END DO
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_same_part, la, na, nc, lb, nb, nd, sumi_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_same_part, lb, nb, nd, la, na, nc, sumi_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_diff_part, la, na, nc, lb, nb, nd, sumi_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_diff_part, lb, nb, nd, la, na, nc, sumi_dp)

								kmin = ABS(la - lb)
								kmax = la + lb
								sum1_sp = 0.0
								sum1_dp = 0.0
								sum2_sp = 0.0
								sum2_dp = 0.0
								DO k = kmin, kmax, 2
									tres_j_cuad = (2 * k + 1) * (ThreeJSymbols_get(2 * la, 2 * k, 2 * lb) ** 2)
									cuad = CUAD2(la, lb, k)
									DO i = 0, 1
										total_ibb = IBrinkBooker(na - 1, la, nb - 1, lb, nd - 1, lb, nc - 1, la, k, x(i))
										IF(nb .NE. nd) THEN
											total_ibb = total_ibb + IBrinkBooker(na - 1, la, nd - 1, lb, nb - 1, lb, nc - 1, la, k, x(i))
											total_ibb = total_ibb / 2.0
										END IF
										sum1_sp = sum1_sp + (coef_exch1_same_part(i) * tres_j_cuad * total_ibb)
										sum1_dp = sum1_dp + (coef_exch1_diff_part(i) * tres_j_cuad * total_ibb)
										sum2_sp = sum2_sp + (coef_exch2_same_part(i) * cuad * tres_j_cuad * total_ibb)
										sum2_dp = sum2_dp + (coef_exch2_diff_part(i) * cuad * tres_j_cuad * total_ibb)
									END DO     
								END DO
								sum1_sp = -sum1_sp
								sum1_dp = -sum1_dp
								sum2_sp = -sum2_sp
								sum2_dp = -sum2_dp
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_same_part, la, na, nc, lb, nb, nd, sum1_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_same_part, lb, nb, nd, la, na, nc, sum1_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_diff_part, la, na, nc, lb, nb, nd, sum1_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_diff_part, lb, nb, nd, la, na, nc, sum1_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_same_part, la, na, nc, lb, nb, nd, sum2_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_same_part, lb, nb, nd, la, na, nc, sum2_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_diff_part, la, na, nc, lb, nb, nd, sum2_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_diff_part, lb, nb, nd, la, na, nc, sum2_dp)
								IF ((N_0 .LE. 16) .AND. (file_error .EQ. 0)) THEN
									WRITE (file_desc, FMT="(I3,I3,I3,I3,I3,I3,E,E,E,E,E,E)", IOSTAT=file_error) &
										la, na, nc, lb, nb, nd, sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp
								END IF
							END DO
						END DO
					END DO
				END DO
			END DO
			PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
		END DO
		CALL SymVBBph_update(vBBph)
		IF (N_0 .LE. 16) CLOSE(file_desc)
		RETURN

	CONTAINS

		SUBROUTINE SymVBBph_update(vBBph)
			TYPE (SymVBBph), INTENT(INOUT) :: vBBph

			CALL SymD3Tensor_SymD3Tensor_add(vBBph%v1_same_part, vBBph%v_local_same_part, vBBph%v1_exch_same_part)
			CALL SymD3Tensor_SymD3Tensor_add(vBBph%v1_diff_part, vBBph%v_local_diff_part, vBBph%v1_exch_diff_part)
			CALL SymD3Tensor_SymD3Tensor_add(vBBph%v2_same_part, vBBph%v_local_same_part, vBBph%v2_exch_same_part)
			CALL SymD3Tensor_SymD3Tensor_add(vBBph%v2_diff_part, vBBph%v_local_diff_part, vBBph%v2_exch_diff_part)
			RETURN
		END SUBROUTINE SymVBBph_update

	END SUBROUTINE SymVBBph_calculate

	SUBROUTINE SymVBBph_product(HF_out, vBBph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVBBph), INTENT(IN) :: vBBph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		CALL SymVBBph_get_Gamma(HF_out, vBBph, HF_in)
		RETURN
	END SUBROUTINE SymVBBph_product

	SUBROUTINE SymVBBph_get_Gamma(HF_out, vBBph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVBBph), INTENT(IN) :: vBBph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		TYPE (SymD3Tensor) same, diff

		CALL SymD3Tensor_new(same)
		CALL SymD3Tensor_new(diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v1_same_part, HF_in%p(0))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v1_diff_part, HF_in%p(1))
		CALL SymD3Tensor_add(HF_out%p(0), same, diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v1_same_part, HF_in%p(1))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v1_diff_part, HF_in%p(0))
		CALL SymD3Tensor_add(HF_out%p(1), same, diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v2_exch_same_part, HF_in%a(0))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v2_exch_diff_part, HF_in%a(1))
		CALL SymD3Tensor_add(HF_out%a(0), same, diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v2_exch_same_part, HF_in%a(1))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v2_exch_diff_part, HF_in%a(0))
		CALL SymD3Tensor_add(HF_out%a(1), same, diff)

		CALL SymD3Tensor_del(same)
		CALL SymD3Tensor_del(diff)
		RETURN
	END SUBROUTINE SymVBBph_get_Gamma

	SUBROUTINE SymVBBph_get_LocalGamma(HF_out, vBBph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVBBph), INTENT(IN) :: vBBph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		TYPE (SymD3Tensor) same, diff

		CALL SymD3Tensor_new(same)
		CALL SymD3Tensor_new(diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v_local_same_part, HF_in%p(0))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v_local_diff_part, HF_in%p(1))
		CALL SymD3Tensor_add(HF_out%p(0), same, diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v_local_same_part, HF_in%p(1))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v_local_diff_part, HF_in%p(0))
		CALL SymD3Tensor_add(HF_out%p(1), same, diff)

		HF_out%a(0) = DBLE(0.0)
		HF_out%a(1) = DBLE(0.0)

		CALL SymD3Tensor_del(same)
		CALL SymD3Tensor_del(diff)
		RETURN
	END SUBROUTINE SymVBBph_get_LocalGamma

	SUBROUTINE SymVBBph_get_ExchangeGamma(HF_out, vBBph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVBBph), INTENT(IN) :: vBBph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		TYPE (SymD3Tensor) same, diff

		CALL SymD3Tensor_new(same)
		CALL SymD3Tensor_new(diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v1_exch_same_part, HF_in%p(0))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v1_exch_diff_part, HF_in%p(1))
		CALL SymD3Tensor_add(HF_out%p(0), same, diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v1_exch_same_part, HF_in%p(1))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v1_exch_diff_part, HF_in%p(0))
		CALL SymD3Tensor_add(HF_out%p(1), same, diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v2_exch_same_part, HF_in%a(0))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v2_exch_diff_part, HF_in%a(1))
		CALL SymD3Tensor_add(HF_out%a(0), same, diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v2_exch_same_part, HF_in%a(1))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v2_exch_diff_part, HF_in%a(0))
		CALL SymD3Tensor_add(HF_out%a(1), same, diff)

		CALL SymD3Tensor_del(same)
		CALL SymD3Tensor_del(diff)
		RETURN
	END SUBROUTINE SymVBBph_get_ExchangeGamma

	FUNCTION SymVBBph_read(vBBph)
		LOGICAL SymVBBph_read
		TYPE (SymVBBph), INTENT(INOUT) :: vBBph

		DOUBLE PRECISION b1
		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER i1, i2, i3, i4, i5, i6
		INTEGER, PARAMETER :: file_desc = 6
		INTEGER file_error
		DOUBLE PRECISION sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp

		OPEN(file_desc, FILE=vBBph%filename, ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "No se pudo leer el archivo: ", vBBph%filename
			SymVBBph_read = .FALSE.
			RETURN
		END IF

		READ (file_desc, FMT="(E)", IOSTAT=file_error) b1
		IF ((file_error .NE. 0) .OR. (b1 .NE. vBBph%b1)) THEN
			CLOSE(file_desc)
			SymVBBph_read = .FALSE.
			RETURN
		END IF

!		PRINT *, "Leyendo:", vBBph%filename
		DO la = 0, N_0
			namax = ((N_0 - la) / 2) + 1
			DO na = 1, namax
				DO nc = 1, na
					DO lb = 0, la
						nbmax = ((N_0 - lb) / 2) + 1
						DO nb = 1, nbmax
							DO nd = 1, nb
								READ (file_desc, FMT="(I3,I3,I3,I3,I3,I3,E,E,E,E,E,E)", IOSTAT=file_error) &
									i1, i2, i3, i4, i5, i6, sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp
								IF ((file_error .NE. 0) .OR. &
									(la .NE. i1) .OR. (na .NE. i2) .OR. (nc .NE. i3) .OR. &
									(lb .NE. i4) .OR. (nb .NE. i5) .OR. (nd .NE. i6)) THEN
									PRINT *, "Informacion no validad en el archivo: ", vBBph%filename
									CLOSE(file_desc)
									SymVBBph_read = .FALSE.
									RETURN
								END IF
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_same_part, la, na, nc, lb, nb, nd, sumi_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_same_part, lb, nb, nd, la, na, nc, sumi_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_diff_part, la, na, nc, lb, nb, nd, sumi_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_diff_part, lb, nb, nd, la, na, nc, sumi_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_same_part, la, na, nc, lb, nb, nd, sum1_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_same_part, lb, nb, nd, la, na, nc, sum1_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_diff_part, la, na, nc, lb, nb, nd, sum1_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_diff_part, lb, nb, nd, la, na, nc, sum1_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_same_part, la, na, nc, lb, nb, nd, sum2_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_same_part, lb, nb, nd, la, na, nc, sum2_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_diff_part, la, na, nc, lb, nb, nd, sum2_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_diff_part, lb, nb, nd, la, na, nc, sum2_dp)
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
		CLOSE(file_desc)
		SymVBBph_read = .TRUE.
		RETURN
	END FUNCTION SymVBBph_read

	SUBROUTINE SymVBBph_del(vBBph)
		TYPE (SymVBBph), INTENT(INOUT) :: vBBph

		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v_local_same_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v_local_diff_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v1_exch_same_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v1_exch_diff_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v2_exch_same_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v2_exch_diff_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v1_same_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v1_diff_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v2_same_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v2_diff_part)
		RETURN
	END SUBROUTINE SymVBBph_del

	SUBROUTINE SymVBBpp_new(vBBpp, b)
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp
		DOUBLE PRECISION, INTENT(IN) :: b

		CALL SymD3Tensor_SymD3Tensor_new(vBBpp%v1_pair)
		CALL SymD3Tensor_SymD3Tensor_new(vBBpp%v2_pair)

		vBBpp%b1 = b
		IF (N_0 < 10) THEN
			WRITE(vBBpp%filename, "(A,I1,A)") "data/vBB", N_0, "pp.txt"
		ELSE
			WRITE(vBBpp%filename, "(A,I2,A)") "data/vBB", N_0, "pp.txt"
		END IF
		RETURN
	END SUBROUTINE SymVBBpp_new

	SUBROUTINE SymVBBpp_calculate(vBBpp)
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp

		INTEGER i
		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER k, kmin, kmax
		DOUBLE PRECISION sum1, sum2
		DOUBLE PRECISION tres_j_cuad, cuad, total_ibb

		DOUBLE PRECISION, DIMENSION(0:1) :: x, coef_pair1, coef_pair2

		INTEGER, PARAMETER :: file_desc = 6
		INTEGER file_error

		IF (N_0 .LE. 16) THEN
			IF (SymVBBpp_read(vBBpp)) RETURN

			OPEN(file_desc, FILE=vBBpp%filename, ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				PRINT *, "*** AVISO: No se pueden escribir los resultados en ", vBBpp%filename
			ELSE
				WRITE (file_desc, FMT="(E)", IOSTAT=file_error) vBBpp%b1
			END IF
		END IF

		DO i = 0, 1
			x(i) = mu(i) / vBBpp%b1
			coef_pair1(i) = 0.5 * I_4PI * (Gogny_W(i, Gogny) - Gogny_B(i, Gogny) - Gogny_H(i, Gogny) + Gogny_M(i, Gogny))
			coef_pair2(i) =       I_4PI * (Gogny_W(i, Gogny) + Gogny_B(i, Gogny) - Gogny_H(i, Gogny) - Gogny_M(i, Gogny))
		END DO

		PRINT *, "Se van a calcular los elementos de matriz pp de BB:"
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
									DO i = 0, 1
										total_ibb = IBrinkBooker(na - 1, la, nc - 1, la, nb - 1, lb, nd - 1, lb, k, x(i))
										IF(nb .NE. nd) THEN
											total_ibb = total_ibb + IBrinkBooker(na - 1, la, nc - 1, la, nd - 1, lb, nb - 1, lb, k, x(i))
											total_ibb = total_ibb / 2.0
										END IF
										sum1 = sum1 + (coef_pair1(i) * tres_j_cuad * total_ibb)
										sum2 = sum2 + (coef_pair2(i) * cuad * tres_j_cuad * total_ibb)
									END DO     
								END DO
								sum1 = sum1 * PAR(la + lb)
								sum2 = sum2 * PAR(la + lb)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v1_pair, la, na, nc, lb, nb, nd, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v1_pair, lb, nb, nd, la, na, nc, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v2_pair, la, na, nc, lb, nb, nd, sum2)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v2_pair, lb, nb, nd, la, na, nc, sum2)
								IF ((N_0 .LE. 16) .AND. (file_error .EQ. 0)) THEN
									WRITE (file_desc, FMT="(I3,I3,I3,I3,I3,I3,E,E)", IOSTAT=file_error) &
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
	END SUBROUTINE SymVBBpp_calculate

	SUBROUTINE SymVBBpp_product(HF_out, vBBpp, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVBBpp), INTENT(IN) :: vBBpp
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		CALL SymVBBpp_get_Delta(HF_out, vBBpp, HF_in)
		RETURN
	END SUBROUTINE SymVBBpp_product

	SUBROUTINE SymVBBpp_get_Delta(HF_out, vBBpp, P_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVBBpp), INTENT(IN) :: vBBpp
		TYPE (SymHartreeFockField), INTENT(IN) :: P_in

		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(0), vBBpp%v1_pair, P_in%p(0))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(1), vBBpp%v1_pair, P_in%p(1))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(0), vBBpp%v2_pair, P_in%a(0))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(1), vBBpp%v2_pair, P_in%a(1))
		RETURN
	END SUBROUTINE SymVBBpp_get_Delta

	FUNCTION SymVBBpp_read(vBBpp)
		LOGICAL SymVBBpp_read
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp

		DOUBLE PRECISION b1
		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER i1, i2, i3, i4, i5, i6
		INTEGER, PARAMETER :: file_desc = 6
		INTEGER file_error
		DOUBLE PRECISION sum1, sum2

		OPEN(file_desc, FILE=vBBpp%filename, ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "No se pudo leer el archivo: ", vBBpp%filename
			SymVBBpp_read = .FALSE.
			RETURN
		END IF

		READ (file_desc, FMT="(E)", IOSTAT=file_error) b1
		IF ((file_error .NE. 0) .OR. (b1 .NE. vBBpp%b1)) THEN
			CLOSE(file_desc)
			SymVBBpp_read = .FALSE.
			RETURN
		END IF

!		PRINT *, "Leyendo:", vBBpp%filename
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
									PRINT *, "Informacion no validad en el archivo: ", vBBpp%filename
									CLOSE(file_desc)
									SymVBBpp_read = .FALSE.
									RETURN
								END IF
								CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v1_pair, la, na, nc, lb, nb, nd, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v1_pair, lb, nb, nd, la, na, nc, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v2_pair, la, na, nc, lb, nb, nd, sum2)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v2_pair, lb, nb, nd, la, na, nc, sum2)
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
		CLOSE(file_desc)

		SymVBBpp_read = .TRUE.
		RETURN
	END FUNCTION SymVBBpp_read

	SUBROUTINE SymVBBpp_del(vBBpp)
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp

		CALL SymD3Tensor_SymD3Tensor_del(vBBpp%v1_pair)
		CALL SymD3Tensor_SymD3Tensor_del(vBBpp%v2_pair)
		RETURN
	END SUBROUTINE SymVBBpp_del

END MODULE symvbb
