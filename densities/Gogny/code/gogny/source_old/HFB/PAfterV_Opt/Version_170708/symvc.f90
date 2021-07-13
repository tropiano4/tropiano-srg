!---------------------------------------------------------------------!
!                                                                     !
!     MATRIX ELEMENTS OF THE COULOMB PART OF THE INTERACTION          !                
!                                                                     !
!---------------------------------------------------------------------!

 MODULE symvc

	USE input
	USE angmom
	USE symd3t_proj
	USE symfield_proj
	USE ic

	IMPLICIT NONE

	TYPE SymVCph
		TYPE (SymD3Tensor_SymD3Tensor) v_local, v1_exch, v2_exch, v1
		CHARACTER(64) filename
	END TYPE

	TYPE SymVCpp
		TYPE (SymD3Tensor_SymD3Tensor) v1_pair, v2_pair
		CHARACTER(64) filename
	END TYPE

 CONTAINS

	!-----------------------------------------------------------------------!
	!  Initialization of the tensors and definition of the filenames (p.h.)	!               
	!-----------------------------------------------------------------------!

	SUBROUTINE SymVCph_new(vCph)
		TYPE (SymVCph), INTENT(INOUT) :: vCph

		CALL SymD3Tensor_SymD3Tensor_new(vCph%v_local)
		CALL SymD3Tensor_SymD3Tensor_new(vCph%v1_exch)
		CALL SymD3Tensor_SymD3Tensor_new(vCph%v2_exch)
		CALL SymD3Tensor_SymD3Tensor_new(vCph%v1)

                SELECT CASE (Basis)

                CASE (1)
			IF (N_0 < 10) THEN
				WRITE(vCph%filename, "(A,I1,A)") "data/vC", N_0, "ph_HO.txt"
			ELSE
				WRITE(vCph%filename, "(A,I2,A)") "data/vC", N_0, "ph_HO.txt"
			END IF
                CASE(2)
                	IF (N_0 < 10) THEN
				WRITE(vCph%filename, "(A,I1,A)") "data/vC", N_0, "ph_WS.txt"
			ELSE
				WRITE(vCph%filename, "(A,I2,A)") "data/vC", N_0, "ph_WS.txt"
			END IF
                END SELECT
		
		RETURN
	END SUBROUTINE SymVCph_new

	!-----------------------------------------------------------------------!
	!     Calculation of the matrix elements of the Coulomb force (p.h.)	!               
	!-----------------------------------------------------------------------!

	SUBROUTINE SymVCph_calculate(vCph, Read_Cph, Lmin)
		TYPE (SymVCph), INTENT(INOUT) :: vCph
		INTEGER, INTENT(INOUT) :: Lmin
		LOGICAL, INTENT(IN) :: Read_Cph

		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER :: k, kmin, kmax, Lold
		DOUBLE PRECISION :: total_iC, tres_j_cuad, cuad
		DOUBLE PRECISION :: sumi, sum1, sum2

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		! If the matrix elements were already calculated, we read them and avoid thereby recalculating them
		! The tensors are then updated with the values read from tape.
		! In case of the WS basis, we (smartly) read what was already calculated and eventually calculate
		! only those elements that are missing.

		IF (Lmin .EQ. 0) THEN
		
			IF (Read_Cph) THEN
				CALL SymD3Tensor_SymD3Tensor_add(vCph%v1, vCph%v_local, vCph%v1_exch)
				RETURN
			ELSE
				IF (SymVCph_read(vCph)) THEN
					CALL SymD3Tensor_SymD3Tensor_add(vCph%v1, vCph%v_local, vCph%v1_exch)
					RETURN
				END IF
			END IF
			
			OPEN(file_desc, FILE=vCph%filename, ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				PRINT *, "*** Attention: Impossible to write in ", vCph%filename
			END IF

		ELSE
		         
			Lold = Lmax
			Lmax = Lold - 1
			
			IF (Lmax < 10) THEN
				WRITE(vCph%filename, "(A,I1,A)") "data/vC", Lmax, "ph_WS.txt"
			ELSE
				WRITE(vCph%filename, "(A,I2,A)") "data/vC", Lmax, "ph_WS.txt"
			END IF
		
			IF (SymVCph_read(vCph)) write(*,'("Read file....")')
			
			Lmax = Lold
		
			IF (Lmax < 10) THEN
				WRITE(vCph%filename, "(A,I1,A)") "data/vC", Lmax, "ph_WS.txt"
			ELSE
				WRITE(vCph%filename, "(A,I2,A)") "data/vC", Lmax, "ph_WS.txt"
			END IF
		
		END IF

                SELECT CASE (Basis)

		! Calculation in the case of the spherical harmonic oscillator

                CASE (1)
		
			PRINT *, "Calculation of ph Coulomb matrix elements - Hamonic oscillator Basis"
			
			DO la = 0, Lmax
				namax = ((N_0 - la) / 2) + 1
				DO na = 1, namax
					DO nc = 1, na
					
						DO lb = 0, la
							nbmax = ((N_0 - lb) / 2) + 1
							DO nb = 1, nbmax
								DO nd = 1, nb
							
									total_iC = ICoulombHO(na - 1, la, nb - 1, lb, nc - 1, la, nd - 1, lb, 0)
								
									IF (nb .NE. nd) THEN
										total_iC = total_iC + ICoulombHO(na - 1, la, nd - 1, lb, nc - 1, la, nb - 1, lb, 0)
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
									
										total_iC = ICoulombHO(na - 1, la, nb - 1, lb, nd - 1, lb, nc - 1, la, k)
									
										IF(nb .NE. nd) THEN
											total_iC = total_iC + ICoulombHO(na - 1, la, nd - 1, lb, nb - 1, lb, nc - 1, la, k)
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
								
									IF (file_error .EQ. 0) THEN
										WRITE (file_desc, FMT="(6I3,3E24.16)", IOSTAT=file_error) &
											la, na, nc, lb, nb, nd, sumi, sum1, sum2
									END IF
								
								END DO
							END DO
						END DO
					END DO
				END DO
				PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
			END DO
			
		! Calculation in the case of a general spherical basis. The tensor algebra is the same as in the case of the
		! harmonic oscillator, only the radial integral is different

                CASE (2)
		
			PRINT *, "Calculation of ph Coulomb matrix elements - General Basis"
			
			DO la = Lmin, Lmax
									namax = MIN(Nmax, NmaxOfL(la))
				IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1
				
				DO na = 1, namax
					DO nc = 1, na
					
						DO lb = 0, la
												nbmax = MIN(Nmax, NmaxOfL(lb))
							IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1
							
							DO nb = 1, nbmax
								DO nd = 1, nb
							
									total_iC = ICoulomb(na, la, nb, lb, nc, la, nd, lb, 0)
								
									IF (nb .NE. nd) THEN
										total_iC = total_iC + ICoulomb(na, la, nd, lb, nc, la, nb, lb, 0)
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
									
										total_iC = ICoulomb(na, la, nb, lb, nd, lb, nc, la, k)
									
										IF (nb .NE. nd) THEN
											total_iC = total_iC + ICoulomb(na, la, nd, lb, nb, lb, nc, la, k)
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
								
									IF (file_error .EQ. 0 .AND. Lmin .EQ. 0) THEN
										WRITE (file_desc, FMT="(6I3,3E24.16)", IOSTAT=file_error) &
											la, na, nc, lb, nb, nd, sumi, sum1, sum2
									END IF
								
								END DO
							END DO
						END DO
					END DO
				END DO
				PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
			END DO
			
                END SELECT
		
		! In case we want to optimize our calculation, we write the full matrix elements after 
		! having calculated only the few ones required
		
		IF (Optimization .EQ. 1 .AND. Lmin .NE. 0) THEN
			
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) THEN
				WRITE(*,'("Incompatible options: No optimization for HO case (or compatible)")')
				STOP 'Error in SymVCph_calculate - No optimization pssible'
			END IF
			
			CALL SymVCph_write(vCph)
			
		END IF
			
		CALL SymD3Tensor_SymD3Tensor_add(vCph%v1, vCph%v_local, vCph%v1_exch)
		
		CLOSE(file_desc)
		
		RETURN
	END SUBROUTINE SymVCph_calculate

        !---------------------------------------------------------------------------------------!
	!											!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVCph_write(vCph)
		TYPE (SymVCph), INTENT(INOUT) :: vCph

		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd
		DOUBLE PRECISION :: sumi, sum1, sum2

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		OPEN(file_desc, FILE=vCph%filename, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to write in ", vCph%filename
		END IF

		DO la = 0, Lmax
		
						namax = MIN(Nmax, NmaxOfL(la))
			IF (CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1
                                
			DO na = 1, namax
				DO nc = 1, na
				
					DO lb = 0, la
									nbmax = MIN(Nmax, NmaxOfL(lb))
						IF (CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1
							
						DO nb = 1, nbmax
							DO nd = 1, nb
								
								CALL SymD3Tensor_SymD3Tensor_get(sumi, vCph%v_local, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum1, vCph%v1_exch, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum2, vCph%v2_exch, lb, nb, nd, la, na, nc)
								
								IF (file_error .EQ. 0) THEN
									WRITE (file_desc, FMT="(6I3,3E24.16)", IOSTAT=file_error) &
									la, na, nc, lb, nb, nd, sumi, sum1, sum2
								END IF
							END DO
						END DO
					END DO
				END DO
			END DO
			
		END DO
		
		CLOSE(file_desc)

		RETURN

	END SUBROUTINE SymVCph_write
	
	SUBROUTINE SymVCph_get_Gamma(HF_out, vCph, HF_in, ta, tb, ProjectionOn)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_out
		TYPE (SymVCph), INTENT(IN) :: vCph
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: HF_in
		INTEGER, INTENT(IN) :: ta, tb, ProjectionOn

		IF (ta .EQ. 0) THEN
			IF (tb .EQ. ta) THEN
				CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(ta), vCph%v1, HF_in%p(tb))
				CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(ta), vCph%v2_exch, HF_in%a(tb))
			ELSE
				HF_out%p(ta) = 0.0
				HF_out%a(ta) = 0.0
			END IF
		ELSE
			HF_out%p(ta) = 0.0
			HF_out%a(ta) = 0.0
		END IF
			
		RETURN
	END SUBROUTINE SymVCph_get_Gamma

	SUBROUTINE SymVCph_get_LocalGamma(HF_out, vCph, HF_in, ta, tb, ProjectionOn)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_out
		TYPE (SymVCph), INTENT(IN) :: vCph
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: HF_in
		INTEGER, INTENT(IN) :: ta, tb, ProjectionOn

		IF (ta .EQ. 0) THEN
			IF (tb .EQ. ta) THEN
				CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(ta), vCph%v_local, HF_in%p(tb))
			ELSE
				HF_out%p(ta) = 0.0
			END IF
			HF_out%a(ta) = 0.0
		ELSE
			HF_out%p(ta) = 0.0
			HF_out%a(ta) = 0.0
		END IF
			
		RETURN
	END SUBROUTINE SymVCph_get_LocalGamma

	SUBROUTINE SymVCph_get_ExchangeGamma(HF_out, vCph, HF_in, ta, tb, ProjectionOn)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_out
		TYPE (SymVCph), INTENT(IN) :: vCph
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: HF_in
		INTEGER, INTENT(IN) :: ta, tb, ProjectionOn

		IF (ta .EQ. 0) THEN
			IF (tb .EQ. ta) THEN
				CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(ta), vCph%v1_exch, HF_in%p(tb))
				CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(ta), vCph%v2_exch, HF_in%a(tb))
			ELSE
				HF_out%p(ta) = 0.0
				HF_out%a(ta) = 0.0
			END IF
		ELSE
			HF_out%p(ta) = 0.0
			HF_out%a(ta) = 0.0
		END IF
		
		RETURN
	END SUBROUTINE SymVCph_get_ExchangeGamma

	!-----------------------------------------------------------------------!
	!   Reading the matrix elements in the p.h. channel from a file. If the !
	!   file exists and its format is valid, we fill in the appropriate 	!
	!   tensors with the values read.					!
	!-----------------------------------------------------------------------!

	FUNCTION SymVCph_read(vCph)
		LOGICAL SymVCph_read
		TYPE (SymVCph), INTENT(INOUT) :: vCph

		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER i1, i2, i3, i4, i5, i6
		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error
		DOUBLE PRECISION :: sumi, sum1, sum2

		OPEN(file_desc, FILE=vCph%filename, ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "No se pudo leer el archivo: ", vCph%filename
			SymVCph_read = .FALSE.
			RETURN
		END IF

		DO la = 0, Lmax
								namax = MIN(Nmax, NmaxOfL(la))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1
			
			DO na = 1, namax
				DO nc = 1, na
					DO lb = 0, la
											nbmax = MIN(Nmax, NmaxOfL(lb))
						IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1
						
						DO nb = 1, nbmax
							DO nd = 1, nb
								READ (file_desc, FMT=*, IOSTAT=file_error) &
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

	!-----------------------------------------------------------------------!
	!  Initialization of the tensors and definition of the filenames (p.p.)	!               
	!-----------------------------------------------------------------------!

	SUBROUTINE SymVCpp_new(vCpp)
		TYPE (SymVCpp), INTENT(INOUT) :: vCpp

		CALL SymD3Tensor_SymD3Tensor_new(vCpp%v1_pair)
		CALL SymD3Tensor_SymD3Tensor_new(vCpp%v2_pair)

                SELECT CASE (Basis)
		
                CASE(1)
			IF (N_0 < 10) THEN
				WRITE(vCpp%filename, "(A,I1,A)") "data/vC", N_0, "pp_HO.txt"
			ELSE
				WRITE(vCpp%filename, "(A,I2,A)") "data/vC", N_0, "pp_HO.txt"
			END IF
                CASE(2)
			IF (N_0 < 10) THEN
				WRITE(vCpp%filename, "(A,I1,A)") "data/vC", N_0, "pp_WS.txt"
			ELSE
				WRITE(vCpp%filename, "(A,I2,A)") "data/vC", N_0, "pp_WS.txt"
			END IF
			
                END SELECT
		
		RETURN
	END SUBROUTINE SymVCpp_new

	!-----------------------------------------------------------------------!
	!     Calculation of the matrix elements of the Coulomb force (p.p.)	!               
	!-----------------------------------------------------------------------!

	SUBROUTINE SymVCpp_calculate(vCpp, Read_Cpp, Lmin)
		TYPE (SymVCpp), INTENT(INOUT) :: vCpp
		INTEGER, INTENT(INOUT) :: Lmin
		LOGICAL, INTENT(IN) :: Read_Cpp

		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER :: k, kmin, kmax, Lold
		DOUBLE PRECISION :: tres_j_cuad, cuad, total_iC
		DOUBLE PRECISION :: sum1, sum2

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		! If the matrix elements were already calculated, we read them and avoid thereby recalculating them
		! The tensors are then updated with the values read from tape.
		! In case of the WS basis, we (smartly) read what was already calculated and eventually calculate
		! only those elements that are missing.

		IF (Lmin .EQ. 0) THEN
		
			IF (Read_Cpp) RETURN
			
			IF (SymVCpp_read(vCpp)) RETURN
			
			OPEN(file_desc, FILE=vCpp%filename, ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				PRINT *, "*** Attention: Impossible to write in ", vCpp%filename
			END IF

		ELSE
		         
			Lold = Lmax
			Lmax = Lold - 1
			
			IF (Lmax < 10) THEN
				WRITE(vCpp%filename, "(A,I1,A)") "data/vC", Lmax, "pp_WS.txt"
			ELSE
				WRITE(vCpp%filename, "(A,I2,A)") "data/vC", Lmax, "pp_WS.txt"
			END IF
		
			IF (SymVCpp_read(vCpp)) write(*,'("Read file....")')
			
			Lmax = Lold
		
			IF (Lmax < 10) THEN
				WRITE(vCpp%filename, "(A,I1,A)") "data/vC", Lmax, "pp_WS.txt"
			ELSE
				WRITE(vCpp%filename, "(A,I2,A)") "data/vC", Lmax, "pp_WS.txt"
			END IF
		
		END IF

                SELECT CASE (Basis)

		! Calculation in the case of the spherical harmonic oscillator

                CASE (1)
		
			PRINT *, "Calculation of pp Coulomb matrix elements - Harmonic Oscillator Basis"
			
			DO la = 0, Lmax
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
										total_iC = (VC * I_4PI) * ICoulombHO(na - 1, la, nc - 1, la, nb - 1, lb, nd - 1, lb, k)
										
										IF (nb .NE. nd) THEN
											total_iC = total_iC + (VC * I_4PI) * &
												ICoulombHO(na - 1, la, nc - 1, la, nd - 1, lb, nb - 1, lb, 0)
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
									
									IF (file_error .EQ. 0) THEN
										WRITE (file_desc, "(6I3,2E24.16)", IOSTAT=file_error) &
										la, na, nc, lb, nb, nd, sum1, sum2
									END IF
								END DO
							END DO
						END DO
					END DO
				END DO
				PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
			END DO
			
		! Calculation in the case of a general spherical basis. The tensor algebra is the same as in the case of the
		! harmonic oscillator, only the radial integral is different

                CASE (2)
		
			PRINT *, "Calculation of pp Coulomb matrix elements - General Basis"
			
			DO la = Lmin, Lmax
									namax = MIN(Nmax, NmaxOfL(la))
				IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1
				
				DO na = 1, namax
					DO nc = 1, na
					
						DO lb = 0, la
												nbmax = MIN(Nmax, NmaxOfL(lb))
							IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1
							
							DO nb = 1, nbmax
								DO nd = 1, nb
								
									kmin = ABS(la - lb)
									kmax = la + lb
									sum1 = 0.0
									sum2 = 0.0
									
									DO k = kmin, kmax, 2
									
										tres_j_cuad = (2 * k + 1) * (ThreeJSymbols_get(2 * la, 2 * k, 2 * lb) ** 2)
										cuad = CUAD2(la, lb, k)
										total_iC = (VC * I_4PI) * ICoulomb(na, la, nc, la, nb, lb, nd, lb, k)
										
										IF (nb .NE. nd) THEN
											total_iC = total_iC + (VC * I_4PI) * ICoulomb(na, la, nc, la, nd, lb, nb, lb, 0)
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
									
									IF (file_error .EQ. 0 .AND. Lmin .EQ. 0) THEN
										WRITE (file_desc, "(6I3,2E24.16)", IOSTAT=file_error) &
										la, na, nc, lb, nb, nd, sum1, sum2
									END IF
								END DO
							END DO
						END DO
					END DO
				END DO
				PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
			END DO
			
                END SELECT
		
		! In case we want to optimize our calculation, we write the full matrix elements after 
		! having calculated only the few ones required
		
		IF (Optimization .EQ. 1 .AND. Lmin .NE. 0) THEN
			
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) THEN
				WRITE(*,'("Incompatible options: No optimization for HO case (or compatible)")')
				STOP 'Error in SymVCpp_calculate - No optimization pssible'
			END IF
			
			CALL SymVCpp_write(vCpp)
			
		END IF
			
		CLOSE(file_desc)
		
		RETURN
	END SUBROUTINE SymVCpp_calculate

        !---------------------------------------------------------------------------------------!
	!											!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVCpp_write(vCpp)
		TYPE (SymVCpp), INTENT(INOUT) :: vCpp

		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd
		DOUBLE PRECISION :: sumi, sum1, sum2

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		OPEN(file_desc, FILE=vCpp%filename, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to write in ", vCpp%filename
		END IF

		DO la = 0, Lmax
		
						namax = MIN(Nmax, NmaxOfL(la))
			IF (CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1
                                
			DO na = 1, namax
				DO nc = 1, na
				
					DO lb = 0, la
									nbmax = MIN(Nmax, NmaxOfL(lb))
						IF (CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1
							
						DO nb = 1, nbmax
							DO nd = 1, nb
								
								CALL SymD3Tensor_SymD3Tensor_get(sum1, vCpp%v1_pair, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum2, vCpp%v2_pair, la, na, nc, lb, nb, nd)
								
								IF (file_error .EQ. 0) THEN
									WRITE (file_desc, FMT="(6I3,2E24.16)", IOSTAT=file_error) &
									la, na, nc, lb, nb, nd, sum1, sum2
								END IF
							END DO
						END DO
					END DO
				END DO
			END DO
			
		END DO
		
		CLOSE(file_desc)

		RETURN

	END SUBROUTINE SymVCpp_write
	
	SUBROUTINE SymVCpp_get_Delta(HF_out, vCpp, P_in, ta, tb)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_out
		TYPE (SymVCpp), INTENT(IN) :: vCpp
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: P_in
		INTEGER, INTENT(IN) :: ta, tb

		IF (ta .EQ. 0) THEN
			CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(ta), vCpp%v1_pair, P_in%p(tb))
			CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(ta), vCpp%v2_pair, P_in%a(tb))
		ELSE
			HF_out%p(ta) = 0.0			
			HF_out%a(ta) = 0.0
		END IF
		
		RETURN
	END SUBROUTINE SymVCpp_get_Delta

	!-----------------------------------------------------------------------!
	!   Reading the matrix elements in the p.p. channel from a file. If the !
	!   file exists and its format is valid, we fill in the appropriate 	!
	!   tensors with the values read.					!
	!-----------------------------------------------------------------------!

	FUNCTION SymVCpp_read(vCpp)
		LOGICAL SymVCpp_read
		TYPE (SymVCpp), INTENT(INOUT) :: vCpp

		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER i1, i2, i3, i4, i5, i6
		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error
		DOUBLE PRECISION :: sum1, sum2

		OPEN(file_desc, FILE=vCpp%filename, ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "No se pudo leer el archivo: ", vCpp%filename
			SymVCpp_read = .FALSE.
			RETURN
		END IF

		DO la = 0, Lmax
								namax = MIN(Nmax, NmaxOfL(la))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1
			
			DO na = 1, namax
				DO nc = 1, na
					DO lb = 0, la
											nbmax = MIN(Nmax, NmaxOfL(lb))
						IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1
						
						DO nb = 1, nbmax
							DO nd = 1, nb
								READ (file_desc, FMT=*, IOSTAT=file_error) &
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
