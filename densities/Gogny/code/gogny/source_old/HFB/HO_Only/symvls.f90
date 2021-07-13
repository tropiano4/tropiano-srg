!----------------------------------------------------------------!
!								 !
!    CALCULATION OF THE SPIN-ORBIT TERM OF THE GOGNY FORCE       !
!								 !
!----------------------------------------------------------------!

 MODULE symvls

	USE input
	USE global
	USE symd3t
	USE symfield
	USE ils

	IMPLICIT NONE

	TYPE SymVLSph
		TYPE (SymD3Tensor_SymD3Tensor) v12, v21
		!CHARACTER(64) filename
	END TYPE

	TYPE SymVLSpp
		TYPE (SymD3Tensor_SymD3Tensor) v22
		!CHARACTER(64) filename
	END TYPE
	
	CHARACTER(64) :: FilePH,  FilePP

 CONTAINS

	!----------------------------------------------------------------!
	!								 !
	!    		PARTICLE-HOLE CHANNEL (MEAN-FIELD)  	         !
	!								 !
	!----------------------------------------------------------------!

	SUBROUTINE SymVLSph_new(vLSph)
		TYPE (SymVLSph), INTENT(INOUT) :: vLSph

		CALL SymD3Tensor_SymD3Tensor_new(vLSph%v12)
		CALL SymD3Tensor_SymD3Tensor_new(vLSph%v21)

		IF (N_0 < 10) THEN
			WRITE(FilePH, "(A,I1,A)") "data/vLS", N_0, "ph_HO.txt"
		ELSE
			WRITE(FilePH, "(A,I2,A)") "data/vLS", N_0, "ph_HO.txt"							
		END IF
			
		RETURN
	END SUBROUTINE SymVLSph_new

	SUBROUTINE SymVLSph_calculate(vLSph)
		TYPE (SymVLSph), INTENT(INOUT) :: vLSph

		INTEGER la, na, namax, nc, lb, nb, nbmax, nd
		DOUBLE PRECISION sum1, sum2

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		IF (SymVLSph_read(vLSph)) RETURN

		OPEN(file_desc, FILE=FilePH, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** AVISO: No se pueden escribir los resultados en ", FilePH
		END IF

		PRINT *, "Calculation of the matrix elements of the spin-orbit term:"
		DO la = 0, N_0
			namax = ((N_0 - la) / 2) + 1
			DO na = 1, namax
				DO nc = 1, na
				
					DO lb = 0, la
						nbmax = ((N_0 - lb) / 2) + 1
						DO nb = 1, nbmax
							DO nd = 1, nb
							
								IF (lb .EQ. 0) THEN
									sum1 = DBLE(0.0)
								ELSE
									sum1 = DBLE(lb * (lb + 1)) * IHFLSHO(nb - 1, nd - 1, lb, na - 1, nc - 1, la)
								END IF
								IF (la .EQ. 0) THEN
									sum2 = DBLE(0.0)
								ELSE
									sum2 = DBLE(la * (la + 1)) * IHFLSHO(na - 1, nc - 1, la, nb - 1, nd - 1, lb)
								END IF
								
								CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v12, la, na, nc, lb, nb, nd, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v21, la, na, nc, lb, nb, nd, sum2)
								
								IF (la .NE. lb) THEN
									! Ojo, estan invertidos los terminos
									CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v12, lb, nb, nd, la, na, nc, sum2)
									CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v21, lb, nb, nd, la, na, nc, sum1)
								END IF
									
								IF (file_error .EQ. 0) THEN
									WRITE (file_desc, FMT="(6I3,2E24.16)", IOSTAT=file_error) &
										la, na, nc, lb, nb, nd, sum1, sum2
								END IF
								
							END DO
						END DO
					END DO
					
				END DO
			END DO
			PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculated"
		END DO
		CLOSE(file_desc)
		
		RETURN
	END SUBROUTINE SymVLSph_calculate

	! Nota: si y0 != 1 este factor cambia la fuerza
	! Se pretende que la fuerza de spin-orbita sea dependiente del isosoin
	SUBROUTINE SymVLSph_get_Gamma(HF_out, vLSph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVLSph), INTENT(IN) :: vLSph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		DOUBLE PRECISION factor
		TYPE (SymD3Tensor) tmp1, tmp2, tmp3

		CALL SymD3Tensor_new(tmp1)
		CALL SymD3Tensor_new(tmp2)
		CALL SymD3Tensor_new(tmp3)

		! Atención: si (y0 .NE. 1), este factor cambia la fuerza.
		! Se pretende que la fuerza de spin-orbita sea dependiente del isosoin.
		factor = Gogny_W0(Gogny) * I_4PI

		CALL SymD3Tensor_product(tmp1, 1.0 + x0, HF_in%a(0))
		CALL SymD3Tensor_product(tmp2, x0, HF_in%a(1))
		CALL SymD3Tensor_add(tmp3, tmp1, tmp2)
		CALL SymD3Tensor_SymD3Tensor_product(tmp1, vLSph%v12, tmp3)
		CALL SymD3Tensor_product(HF_out%p(0), factor, tmp1)

		CALL SymD3Tensor_product(tmp1, 1.0 + x0, HF_in%a(1))
		CALL SymD3Tensor_product(tmp2, x0, HF_in%a(0))
		CALL SymD3Tensor_add(tmp3, tmp1, tmp2)
		CALL SymD3Tensor_SymD3Tensor_product(tmp1, vLSph%v12, tmp3)
		CALL SymD3Tensor_product(HF_out%p(1), factor, tmp1)

		CALL SymD3Tensor_product(tmp1, 1.0 + x0, HF_in%p(0))
		CALL SymD3Tensor_product(tmp2, x0, HF_in%p(1))
		CALL SymD3Tensor_add(tmp3, tmp1, tmp2)
		CALL SymD3Tensor_SymD3Tensor_product(tmp1, vLSph%v21, tmp3)
		CALL SymD3Tensor_product(HF_out%a(0), factor, tmp1)

		CALL SymD3Tensor_product(tmp1, 1.0 + x0, HF_in%p(1))
		CALL SymD3Tensor_product(tmp2, x0, HF_in%p(0))
		CALL SymD3Tensor_add(tmp3, tmp1, tmp2)
		CALL SymD3Tensor_SymD3Tensor_product(tmp1, vLSph%v21, tmp3)
		CALL SymD3Tensor_product(HF_out%a(1), factor, tmp1)

		CALL SymD3Tensor_del(tmp1)
		CALL SymD3Tensor_del(tmp2)
		CALL SymD3Tensor_del(tmp3)
		RETURN
	END SUBROUTINE SymVLSph_get_Gamma

	FUNCTION SymVLSph_read(vLSph)
		LOGICAL SymVLSph_read
		TYPE (SymVLSph), INTENT(INOUT) :: vLSph

		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER i1, i2, i3, i4, i5, i6
		DOUBLE PRECISION sum1, sum2

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		OPEN(file_desc, FILE=FilePH, ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "No se pudo leer el archivo: ", FilePH
			SymVLSph_read = .FALSE.
			RETURN
		END IF

		DO la = 0, N_0

			namax = ((N_0 - la) / 2) + 1
			
			DO na = 1, namax
				DO nc = 1, na
				
					DO lb = 0, la

						nbmax = ((N_0 - lb) / 2) + 1
						
						DO nb = 1, nbmax
							DO nd = 1, nb
							
								READ (file_desc, FMT=*, IOSTAT=file_error) &
									i1, i2, i3, i4, i5, i6, sum1, sum2
								IF ((file_error .NE. 0) .OR. &
									(la .NE. i1) .OR. (na .NE. i2) .OR. (nc .NE. i3) .OR. &
									(lb .NE. i4) .OR. (nb .NE. i5) .OR. (nd .NE. i6)) THEN
									PRINT *, "Informacion no validad en el archivo: ", FilePH
									CLOSE(file_desc)
									SymVLSph_read = .FALSE.
									RETURN
								END IF
								
								CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v12, la, na, nc, lb, nb, nd, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v21, la, na, nc, lb, nb, nd, sum2)
								
								IF (la .NE. lb) THEN
									! Ojo, están invertidos los términos
									CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v12, lb, nb, nd, la, na, nc, sum2)
									CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v21, lb, nb, nd, la, na, nc, sum1)
								END IF
								
							END DO
						END DO
					END DO
					
				END DO
			END DO
		END DO
		CLOSE(file_desc)

		SymVLSph_read = .TRUE.
		RETURN
	END FUNCTION SymVLSph_read

	SUBROUTINE SymVLSph_del(vLSph)
		TYPE (SymVLSph), INTENT(INOUT) :: vLSph

		CALL SymD3Tensor_SymD3Tensor_del(vLSph%v12)
		CALL SymD3Tensor_SymD3Tensor_del(vLSph%v21)
		RETURN
	END SUBROUTINE SymVLSph_del

	!----------------------------------------------------------------!
	!								 !
	!    		PARTICLE-PARTICLE CHANNEL (PAIRING)  	         !
	!								 !
	!----------------------------------------------------------------!

	SUBROUTINE SymVLSpp_new(vLSpp)
		TYPE (SymVLSpp), INTENT(INOUT) :: vLSpp

		CALL SymD3Tensor_SymD3Tensor_new(vLSpp%v22)

		IF (N_0 < 10) THEN
			WRITE(FilePP, "(A,I1,A)") "data/vLS", N_0, "pp_HO.txt"
		ELSE
			WRITE(FilePP, "(A,I2,A)") "data/vLS", N_0, "pp_HO.txt"
		END IF
			
		RETURN
	END SUBROUTINE SymVLSpp_new

	SUBROUTINE SymVLSpp_calculate(vLSpp)
		TYPE (SymVLSpp), INTENT(INOUT) :: vLSpp

		INTEGER la, na, namax, nc, lb, nb, nbmax, nd
		DOUBLE PRECISION sumi

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		IF (SymVLSpp_read(vLSpp)) RETURN

		OPEN(file_desc, FILE=FilePP, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** AVISO: No se pueden escribir los resultados en ", FilePP
		END IF

		PRINT *, "Se van a calcular los elementos de matriz de apar de LS"
		DO la = 0, N_0
			namax = ((N_0 - la) / 2) + 1
			DO na = 1, namax
				DO nc = 1, na
				
					DO lb = 0, la
						nbmax = ((N_0 - lb) / 2) + 1
						DO nb = 1, nbmax
							DO nd = 1, nb
							
								IF ((la .EQ. 0) .OR. (lb .EQ. 0)) THEN
									sumi = 0.0
								ELSE
									sumi = DBLE(PAR(la + lb) * 2.0 * la * (la + 1) * lb * (lb + 1)) * &
										IPLSHO(na - 1, nc - 1, la, nb - 1, nd - 1, lb)
								END IF
									
								CALL SymD3Tensor_SymD3Tensor_assign(vLSpp%v22, la, na, nc, lb, nb, nd, sumi)
									
								IF (la .NE. lb) THEN
									CALL SymD3Tensor_SymD3Tensor_assign(vLSpp%v22, lb, nb, nd, la, na, nc, sumi)
								END IF
									
								IF (file_error .EQ. 0) THEN
									WRITE (file_desc, FMT="(6I3,E24.16)", IOSTAT=file_error) &
											la, na, nc, lb, nb, nd, sumi
								END IF
								
							END DO
						END DO
					END DO
					
				END DO
			END DO
			PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
		END DO
			
	 	CLOSE(file_desc)
	
		RETURN
	END SUBROUTINE SymVLSpp_calculate

	SUBROUTINE SymVLSpp_get_Delta(HF_out, vLSpp, P_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVLSpp), INTENT(IN) :: vLSpp
		TYPE (SymHartreeFockField), INTENT(IN) :: P_in

		DOUBLE PRECISION factor
		TYPE (SymD3Tensor) tmp

		CALL SymD3Tensor_new(tmp)

		HF_out%p(0) = DBLE(0.0)
		HF_out%p(1) = DBLE(0.0)

		factor = Gogny_W0(Gogny) * I_4PI

		CALL SymD3Tensor_SymD3Tensor_product(tmp, vLSpp%v22, P_in%a(0))
		CALL SymD3Tensor_product(HF_out%a(0), factor, tmp)

		CALL SymD3Tensor_SymD3Tensor_product(tmp, vLSpp%v22, P_in%a(1))
		CALL SymD3Tensor_product(HF_out%a(1), factor, tmp)

		CALL SymD3Tensor_new(tmp)
		RETURN
	END SUBROUTINE SymVLSpp_get_Delta

	FUNCTION SymVLSpp_read(vLSpp)
		LOGICAL SymVLSpp_read
		TYPE (SymVLSpp), INTENT(INOUT) :: vLSpp

		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER i1, i2, i3, i4, i5, i6
		DOUBLE PRECISION sumi

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		OPEN(file_desc, FILE=FilePP, ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "No se pudo leer el archivo: ", FilePP
			SymVLSpp_read = .FALSE.
			RETURN
		END IF

		DO la = 0, N_0

			namax = ((N_0 - la) / 2) + 1
			
			DO na = 1, namax
				DO nc = 1, na
				
					DO lb = 0, la

						nbmax = ((N_0 - lb) / 2) + 1
						
						DO nb = 1, nbmax
							DO nd = 1, nb
							
								READ (file_desc, FMT=*, IOSTAT=file_error) &
									i1, i2, i3, i4, i5, i6, sumi
								IF ((file_error .NE. 0) .OR. &
									(la .NE. i1) .OR. (na .NE. i2) .OR. (nc .NE. i3) .OR. &
									(lb .NE. i4) .OR. (nb .NE. i5) .OR. (nd .NE. i6)) THEN
									PRINT *, "Informacion no validad en el archivo: ", FilePP
									CLOSE(file_desc)
									SymVLSpp_read = .FALSE.
									RETURN
								END IF
								CALL SymD3Tensor_SymD3Tensor_assign(vLSpp%v22, la, na, nc, lb, nb, nd, sumi)
								IF (la .NE. lb) THEN
									CALL SymD3Tensor_SymD3Tensor_assign(vLSpp%v22, lb, nb, nd, la, na, nc, sumi)
								END IF
								
							END DO
						END DO
					END DO
					
				END DO
			END DO
		END DO
		CLOSE(file_desc)

		SymVLSpp_read = .TRUE.
		RETURN
	END FUNCTION SymVLSpp_read

	SUBROUTINE SymVLSpp_del(vLSpp)
		TYPE (SymVLSpp), INTENT(INOUT) :: vLSpp

		CALL SymD3Tensor_SymD3Tensor_del(vLSpp%v22)
		RETURN
	END SUBROUTINE SymVLSpp_del

END MODULE symvls
