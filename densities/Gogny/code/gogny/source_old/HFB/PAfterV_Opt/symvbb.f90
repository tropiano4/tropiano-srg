!----------------------------------------------------------------!
!								 !
!  CALCULATION OF THE BRINK-BOEKER TERM OF THE GOGNY FORCE       !
!								 !
!----------------------------------------------------------------!

 MODULE symvbb

	USE input
	USE angmom
	USE symd3t_proj
	USE symfield_proj
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

	! Range of the Brink-Boker term
	
	INTEGER :: Lmax_read, Nmax_read
	
	DOUBLE PRECISION, DIMENSION(0:1) :: mu
	DATA mu /0.7, 1.2/

 CONTAINS

        !---------------------------------------------------------------------------------------!
	!											!
        !   Creates a new object of the type SymVBBph, that is in fact a collection of 		!
	!   "3D tensors" arrays SymD3Tensor_SymD3Tensor, as defined in module symd3t.f90	!
	!   Here we create each of the corresponding parts of the type SymVBBph and allocate 	!
	!   the appropriate amount of memory							!
	!   											!
	!   This tensor is for the Brink-Boker matrix elements in the particle-hole channel	!
	!											!
        !---------------------------------------------------------------------------------------!
	
	SUBROUTINE SymVBBph_new(vBBph, b)
		TYPE (SymVBBph), INTENT(INOUT) :: vBBph
		DOUBLE PRECISION, INTENT(IN) :: b

		! Create the tensors

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
		
                SELECT CASE (Basis)
		
		CASE (1)
         
			IF (N_0 < 10) THEN
				WRITE(vBBph%filename, "(A,I1,A)") "data/vBB", N_0, "ph_HO.txt"
			ELSE
				WRITE(vBBph%filename, "(A,I2,A)") "data/vBB", N_0, "ph_HO.txt"
			END IF

		CASE (2)
         
			IF (N_0 < 10) THEN
				WRITE(vBBph%filename, "(A,I1,A)") "data/vBB", N_0, "ph_WS.txt"
			ELSE
				WRITE(vBBph%filename, "(A,I2,A)") "data/vBB", N_0, "ph_WS.txt"
			END IF

		END SELECT
		
		RETURN		
	END SUBROUTINE SymVBBph_new

        !---------------------------------------------------------------------------------------!
	!											!
	!   Subroutine calculating the matrix elements of the Brink-Boker term. We distinguish	!
	!   2 special cases, either the calculation is done in the harmonicoscillator basis 	!
	!   (basis = 1) or in a general spherical basis (basis = 2)				!
	!											!
	!   BEWARE: Here the HO length b is implicitely assumed equal to 1, otherwise,there 	!
	!           should be a factor 1/b							!
	!											!
	!   Refs: Appendix F									!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVBBph_calculate(vBBph, read_BB, Lmin)
		TYPE (SymVBBph), INTENT(INOUT) :: vBBph
		INTEGER, INTENT(INOUT) :: Lmin
		LOGICAL, INTENT(IN) :: read_BB

		INTEGER i, icount
		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		DOUBLE PRECISION sumi_sp, sumi_dp, total_ibb
		INTEGER k, kmin, kmax, Lold
		DOUBLE PRECISION sum1_sp, sum1_dp, sum2_sp, sum2_dp
		DOUBLE PRECISION tres_j_cuad, cuad,xxx

		DOUBLE PRECISION, DIMENSION(0:1) :: x, &
			coef_local_same_part, coef_local_diff_part, &
			coef_exch1_same_part, coef_exch1_diff_part, &
			coef_exch2_same_part, coef_exch2_diff_part

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		! If the matrix elements were already calculated, we read them and avoid thereby recalculating them
		! The tensors are then updated with the values read from tape.
		! In case of the WS basis, we (smartly) read what was already calculated and eventually calculate
		! only those elements that are missing.

		IF (Lmin .EQ. 0) THEN
		
			IF (read_BB) THEN
				CALL SymVBBph_update(vBBph)
				RETURN
			ELSE
				IF (SymVBBph_read(vBBph)) THEN
					CALL SymVBBph_update(vBBph)
					RETURN
				END IF
			END IF
			
			OPEN(file_desc, FILE=vBBph%filename, ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				PRINT *, "*** Attention: Impossible to write in ", vBBph%filename
			ELSE
				WRITE (file_desc, FMT="(E)", IOSTAT=file_error) vBBph%b1
			END IF

		ELSE
		         
			Lold = Lmax
			Lmax = Lold - 1
			
			IF (Lmax < 10) THEN
				WRITE(vBBph%filename, "(A,I1,A)") "data/vBB", Lmax, "ph_WS.txt"
			ELSE
				WRITE(vBBph%filename, "(A,I2,A)") "data/vBB", Lmax, "ph_WS.txt"
			END IF
		
			IF (SymVBBph_read(vBBph)) write(*,'("Read file....")')
			
			Lmax = Lold
		
			IF (Lmax < 10) THEN
				WRITE(vBBph%filename, "(A,I1,A)") "data/vBB", Lmax, "ph_WS.txt"
			ELSE
				WRITE(vBBph%filename, "(A,I2,A)") "data/vBB", Lmax, "ph_WS.txt"
			END IF
		
		END IF

		! For each isospin, we calculate the constant numerical factors in front of the matrix elements 
		! names with the suffix "_same_part" refer to the terms with the delta_(ta,tb), names with the
		! suffix "_diff_part" refer to the terms with the other terms:
		!
		! Refs.: Page 130-131, Definition of v1eBB and v2eBB

		DO i = 0, 1
			x(i) = mu(i) / vBBph%b1
			coef_local_diff_part(i) = I_4PI * (Gogny_W(i, Gogny) + (0.5 * Gogny_B(i, Gogny)))
			coef_local_same_part(i) = coef_local_diff_part(i) - (I_4PI * (Gogny_H(i, Gogny) + (0.5 * Gogny_M(i, Gogny))))
			coef_exch1_diff_part(i) = - I_4PI * ((0.5 * Gogny_H(i, Gogny)) + Gogny_M(i, Gogny))
			coef_exch1_same_part(i) = coef_exch1_diff_part(i) + (I_4PI * ((0.5 * Gogny_W(i, Gogny)) + Gogny_B(i, Gogny)))
			coef_exch2_diff_part(i) = - I_4PI * Gogny_H(i, Gogny)
			coef_exch2_same_part(i) = coef_exch2_diff_part(i) + (I_4PI * Gogny_W(i, Gogny))
		END DO

                SELECT CASE (Basis)
		
		CASE (1)
         
			! Calculation of the matrix elements v1eBB and v2eBB as defined in Page 131  ---  CHECKED AND OK

			PRINT *, "Brink-Boker terms: Particle-Hole Channel - Harmonic Oscillator Basis"
			! Los mas altos son los mas problematicos
			DO la = 0, Lmax
				namax = ((N_0 - la) / 2) + 1
				DO na = 1, namax
					DO nc = 1, na
				
						DO lb = 0, la
							nbmax = ((N_0 - lb) / 2) + 1
							DO nb = 1, nbmax
								DO nd = 1, nb
							
									! Calculation of the term k=0 in the multipole expansion over Wk(r1,r2). This actually 
									! corresponds to the local term

									sumi_sp = 0.0
									sumi_dp = 0.0
								
									DO i = 0, 1
										total_ibb = IBrinkBookerHO(na - 1, la, nb - 1, lb, nc - 1, la, nd - 1, lb, 0, x(i))
										IF (nb .NE. nd) THEN
											total_ibb = total_ibb + IBrinkBookerHO(na - 1, la, nd - 1, lb, nc - 1, la, nb - 1, lb, 0, x(i))
											total_ibb = total_ibb / 2.0
										END IF
										sumi_sp = sumi_sp + (coef_local_same_part(i) * total_ibb)
										sumi_dp = sumi_dp + (coef_local_diff_part(i) * total_ibb)
									END DO
								
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_same_part, la, na, nc, lb, nb, nd, sumi_sp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_same_part, lb, nb, nd, la, na, nc, sumi_sp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_diff_part, la, na, nc, lb, nb, nd, sumi_dp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_diff_part, lb, nb, nd, la, na, nc, sumi_dp)

									! Summations over the multipole Wk(r1, r2) (Appendix B, formula above B2, page 113) are restricted 
									! by the angular part. Only the terms with ABS(la - lb) <= k <= la + lb give non-zero terms
		
									kmin = ABS(la - lb)
									kmax = la + lb
									
									sum1_sp = 0.0
									sum1_dp = 0.0
									sum2_sp = 0.0
									sum2_dp = 0.0
								
									DO k = kmin, kmax, 2
										tres_j_cuad = (2*k + 1) * (ThreeJSymbols_get(2*la, 2*k, 2*lb) ** 2)
										cuad = CUAD2(la, lb, k)
										DO i = 0, 1
											total_ibb = IBrinkBookerHO(na - 1, la, nb - 1, lb, nd - 1, lb, nc - 1, la, k, x(i))
											IF(nb .NE. nd) THEN
												total_ibb = total_ibb + IBrinkBookerHO(na - 1, la, nd - 1, lb, nb - 1, lb, nc - 1, la, k, x(i))
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
								
									IF (file_error .EQ. 0) THEN
										WRITE (file_desc, FMT="(6I3,6E24.16)", IOSTAT=file_error) &
											la, na, nc, lb, nb, nd, sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp
									END IF
								END DO
							END DO
						END DO
					END DO
				END DO
				PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
			END DO

		CASE(2)
		
			PRINT *, "Brink-Boker terms: Particle-Hole Channel - General Basis"

			DO la = Lmin, Lmax
		
			 	icount = 0
							namax = MIN(Nmax, NmaxOfL(la))
				IF (CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1
                                
				DO na = 1, namax
					DO nc = 1, na
				
						DO lb = 0, la
										nbmax = MIN(Nmax, NmaxOfL(lb))
							IF (CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1
							
							DO nb = 1, nbmax
								DO nd = 1, nb
							
									sumi_sp = 0.0
									sumi_dp = 0.0
								
									DO i = 0, 1
										total_ibb = IBrinkBooker(na, la, nb, lb, nc, la, nd, lb, 0, mu(i))
										IF (nb .NE. nd) THEN
											total_ibb = total_ibb + IBrinkBooker(na, la, nd, lb, nc, la, nb, lb, 0, mu(i))
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
										tres_j_cuad = (2*k + 1.0) * (ThreeJSymbols_get(2*la, 2*k, 2*lb) ** 2)
										cuad = CUAD2(la, lb, k)
										DO i = 0, 1
											total_ibb = IBrinkBooker(na, la, nb, lb, nd, lb, nc, la, k, mu(i))
											IF(nb .NE. nd) THEN
												total_ibb = total_ibb + IBrinkBooker(na, la, nd, lb, nb, lb, nc, la, k, mu(i))
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
								
									IF (file_error .EQ. 0) THEN
										WRITE (file_desc, FMT="(6I3,6E24.16)", IOSTAT=file_error) &
											la, na, nc, lb, nb, nd, sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp
									END IF
								END DO
							END DO
						END DO
					END DO
				END DO
				WRITE(*,'("Number of re-used integrals = ",I16)') icount
				WRITE(*,'("Completion rate: ",I3,"%")') INT(100 * (la + 1) / (N_0 + 1))
			END DO
		
		END SELECT
			
		! In case we want to optimize our calculation, we write the full matrix elements after 
		! having calculated only the few ones required
		
		IF (Optimization .EQ. 1 .AND. Lmin .NE. 0) THEN
			
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) THEN
				WRITE(*,'("Incompatible options: No optimization for HO case (or compatible)")')
				STOP 'Error in SymVBBph_calculate - No optimization pssible'
			END IF
			
			CALL SymVBBph_write(vBBph)
			
		END IF
			
		CALL SymVBBph_update(vBBph)
                
		CLOSE(file_desc)

		RETURN

	CONTAINS

		!---------------------------------------------------------------------------------------!
		!   Subroutine updating the tensor fields of the Brink-Boeker term in the particle-	!
		!   hole channel after reading the matrix elements from tape				!
        	!---------------------------------------------------------------------------------------!

		SUBROUTINE SymVBBph_update(vBBph)
			TYPE (SymVBBph), INTENT(INOUT) :: vBBph

			CALL SymD3Tensor_SymD3Tensor_add(vBBph%v1_same_part, vBBph%v_local_same_part, vBBph%v1_exch_same_part)
			CALL SymD3Tensor_SymD3Tensor_add(vBBph%v1_diff_part, vBBph%v_local_diff_part, vBBph%v1_exch_diff_part)
			CALL SymD3Tensor_SymD3Tensor_add(vBBph%v2_same_part, vBBph%v_local_same_part, vBBph%v2_exch_same_part)
			CALL SymD3Tensor_SymD3Tensor_add(vBBph%v2_diff_part, vBBph%v_local_diff_part, vBBph%v2_exch_diff_part)
			
			RETURN
		END SUBROUTINE SymVBBph_update

	END SUBROUTINE SymVBBph_calculate

        !---------------------------------------------------------------------------------------!
	!											!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVBBph_write(vBBph)
		TYPE (SymVBBph), INTENT(INOUT) :: vBBph

		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd
		DOUBLE PRECISION :: sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		OPEN(file_desc, FILE=vBBph%filename, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to write in ", vBBph%filename
		ELSE
			WRITE (file_desc, FMT="(E)", IOSTAT=file_error) vBBph%b1
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
								
								CALL SymD3Tensor_SymD3Tensor_get(sumi_sp, vBBph%v_local_same_part, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sumi_dp, vBBph%v_local_diff_part, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum1_sp, vBBph%v1_exch_same_part, lb, nb, nd, la, na, nc)
								CALL SymD3Tensor_SymD3Tensor_get(sum1_dp, vBBph%v1_exch_diff_part, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum2_sp, vBBph%v2_exch_same_part, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum2_dp, vBBph%v2_exch_diff_part, lb, nb, nd, la, na, nc)
								
								IF (file_error .EQ. 0) THEN
									WRITE (file_desc, FMT="(6I3,6E24.16)", IOSTAT=file_error) &
									la, na, nc, lb, nb, nd, sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp
								END IF
							END DO
						END DO
					END DO
				END DO
			END DO
			
		END DO
		
		CLOSE(file_desc)

		RETURN

	END SUBROUTINE SymVBBph_write
	
	! Produce Gamma(tb) by contracting V(ta) with rho(tb)

	SUBROUTINE SymVBBph_get_Gamma(HF_out, vBBph, HF_in, ta, tb)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_out
		TYPE (SymVBBph), INTENT(IN) :: vBBph
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: HF_in
		INTEGER, INTENT(IN) :: ta, tb

		TYPE (SymD3Tensor) :: same, diff

		CALL SymD3Tensor_new(same)
		CALL SymD3Tensor_new(diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v1_same_part, HF_in%p(ta))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v1_diff_part, HF_in%p(1-ta))
		CALL SymD3Tensor_add(HF_out%p(tb), same, diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v2_exch_same_part, HF_in%a(ta))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v2_exch_diff_part, HF_in%a(1-ta))
		CALL SymD3Tensor_add(HF_out%a(tb), same, diff)

		HF_out%GaugeAngle(tb) = HF_in%GaugeAngle(ta)

		CALL SymD3Tensor_del(same)
		CALL SymD3Tensor_del(diff)
		
		RETURN
	END SUBROUTINE SymVBBph_get_Gamma

	SUBROUTINE SymVBBph_get_LocalGamma(HF_out, vBBph, HF_in, ta, tb, ProjectionOn)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_out
		TYPE (SymVBBph), INTENT(IN) :: vBBph
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: HF_in
		INTEGER, INTENT(IN) :: ta, tb, ProjectionOn

		TYPE (SymD3Tensor) :: same, diff

		CALL SymD3Tensor_new(same)
		CALL SymD3Tensor_new(diff)

		IF (ProjectionOn .EQ. 1) THEN
		
			IF (ta .EQ. tb) THEN
				CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(tb), vBBph%v_local_same_part, HF_in%p(ta))
			ELSE
				CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(tb), vBBph%v_local_diff_part, HF_in%p(ta))
			END IF
		
			HF_out%a(tb) = 0.0

			HF_out%GaugeAngle(tb) = HF_in%GaugeAngle(ta)

		ELSE
		
			CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v_local_same_part, HF_in%p(ta))
			CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v_local_diff_part, HF_in%p(1-ta))
			CALL SymD3Tensor_add(HF_out%p(tb), same, diff)
		
			HF_out%a(tb) = 0.0

			HF_out%GaugeAngle(tb) = HF_in%GaugeAngle(ta)

		END IF

		CALL SymD3Tensor_del(same)
		CALL SymD3Tensor_del(diff)
		
		RETURN
	END SUBROUTINE SymVBBph_get_LocalGamma

	SUBROUTINE SymVBBph_get_ExchangeGamma(HF_out, vBBph, HF_in, ta, tb, ProjectionOn)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_out
		TYPE (SymVBBph), INTENT(IN) :: vBBph
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: HF_in
		INTEGER, INTENT(IN) :: ta, tb, ProjectionOn

		TYPE (SymD3Tensor) same, diff

		CALL SymD3Tensor_new(same)
		CALL SymD3Tensor_new(diff)

		IF (ProjectionOn .EQ. 1) THEN
		
			IF (ta .EQ. tb) THEN			
				CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(tb), vBBph%v1_exch_same_part, HF_in%p(ta))
				CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(tb), vBBph%v2_exch_same_part, HF_in%a(ta))
			ELSE
				CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(tb), vBBph%v1_exch_diff_part, HF_in%p(ta))
				CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(tb), vBBph%v2_exch_diff_part, HF_in%a(ta))
			END IF

			HF_out%GaugeAngle(tb) = HF_in%GaugeAngle(ta)

		ELSE
		
			CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v1_exch_same_part, HF_in%p(ta))
			CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v1_exch_diff_part, HF_in%p(1-ta))
			CALL SymD3Tensor_add(HF_out%p(tb), same, diff)

			CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v2_exch_same_part, HF_in%a(ta))
			CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v2_exch_diff_part, HF_in%a(1-ta))
			CALL SymD3Tensor_add(HF_out%a(tb), same, diff)

			HF_out%GaugeAngle(tb) = HF_in%GaugeAngle(ta)
			
		END IF
		
		CALL SymD3Tensor_del(same)
		CALL SymD3Tensor_del(diff)
		
		RETURN
	END SUBROUTINE SymVBBph_get_ExchangeGamma

        !---------------------------------------------------------------------------------------!
	!											!
	!   Subroutine reading the matrix elements of the Brink-Boker term in the particle-hole	!
	!   channel. This subroutine attempts reading a file containing previously calculated	!
	!   matrix elements. It is used to decide whether or not to recalculate eveything	!
	!											!
        !---------------------------------------------------------------------------------------!

	FUNCTION SymVBBph_read(vBBph)
		LOGICAL SymVBBph_read
		TYPE (SymVBBph), INTENT(INOUT) :: vBBph

		DOUBLE PRECISION b1
		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER i1, i2, i3, i4, i5, i6
		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error
		DOUBLE PRECISION sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp

		OPEN(file_desc, FILE=vBBph%filename, ACTION="READ", IOSTAT=file_error)
		
		IF (file_error .NE. 0) THEN
			PRINT *, "Impossible to read file: ", vBBph%filename
			SymVBBph_read = .FALSE.
			RETURN
		END IF

		READ (file_desc, FMT="(E)", IOSTAT=file_error) b1
		
		IF ((file_error .NE. 0) .OR. (b1 .NE. vBBph%b1)) THEN
			CLOSE(file_desc)
			SymVBBph_read = .FALSE.
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
									i1, i2, i3, i4, i5, i6, sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp
									
								IF ((file_error .NE. 0) .OR. &
									(la .NE. i1) .OR. (na .NE. i2) .OR. (nc .NE. i3) .OR. &
									(lb .NE. i4) .OR. (nb .NE. i5) .OR. (nd .NE. i6)) THEN
									PRINT *, "Invalid information in file: ", vBBph%filename
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

        !---------------------------------------------------------------------------------------!
	!											!
	!   Subroutine deleting the tensors corresponding to the matrix elements of the Brink-	!
	!   Boker in the particle-hole channel.							!
	!											!
        !---------------------------------------------------------------------------------------!

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

        !---------------------------------------------------------------------------------------!
	!											!
        !   Creates a new object of the type SymVBBph, that is in fact a collection of 		!
	!   "3D tensors" arrays SymD3Tensor_SymD3Tensor, as defined in module symd3t.f90	!
	!   Here we create each of the corresponding parts of the type SymVBBph and allocate 	!
	!   the appropriate amount of memory							!
	!   											!
	!   This tensor is for the Brink-Boker matrix elements in the particle-particle channel	!
	!											!
        !---------------------------------------------------------------------------------------!
	
	SUBROUTINE SymVBBpp_new(vBBpp, b)
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp
		DOUBLE PRECISION, INTENT(IN) :: b

		CALL SymD3Tensor_SymD3Tensor_new(vBBpp%v1_pair)
		CALL SymD3Tensor_SymD3Tensor_new(vBBpp%v2_pair)

		vBBpp%b1 = b
		
                SELECT CASE (Basis)
		
		CASE (1)
         
			IF (N_0 < 10) THEN
				WRITE(vBBpp%filename, "(A,I1,A)") "data/vBB", N_0, "pp_HO.txt"
			ELSE
				WRITE(vBBpp%filename, "(A,I2,A)") "data/vBB", N_0, "pp_HO.txt"
			END IF
			
		CASE(2)
		
			IF (N_0 < 10) THEN
				WRITE(vBBpp%filename, "(A,I1,A)") "data/vBB", N_0, "pp_WS.txt"
			ELSE
				WRITE(vBBpp%filename, "(A,I2,A)") "data/vBB", N_0, "pp_WS.txt"
			END IF
		
		END SELECT
			
		RETURN
	END SUBROUTINE SymVBBpp_new

        !---------------------------------------------------------------------------------------!
	!											!
	!   Subroutine calculating the matrix elements of the Brink-Boker term in the particle-	!
	!   particle channel. We distinguish 2 special cases, either the calculation is done 	!
	!   in the harmonic oscillator basis (basis = 1) or in a general spherical basis 	!
	!   (basis = 2)										!	
	!											!
	!   BEWARE: Here the HO length b is implicitely assumed equal to 1, otherwise,there 	!
	!           should be a factor 1/b							!
	!											!
	!   Refs: Appendix F, Page 138								!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVBBpp_calculate(vBBpp, Read_BBpp, Lmin)
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp
		INTEGER, INTENT(INOUT) :: Lmin
		LOGICAL, INTENT(IN) :: Read_BBpp

		INTEGER i, icount
		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER k, kmin, kmax, Lold
		DOUBLE PRECISION sum1, sum2
		DOUBLE PRECISION tres_j_cuad, cuad, total_ibb

		DOUBLE PRECISION, DIMENSION(0:1) :: x, coef_pair1, coef_pair2

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		! If the matrix elements were already calculated, we read them and avoid thereby recalculating them
		! The tensors are then updated with the values read from tape

		IF (Lmin .EQ. 0) THEN
		
			IF (Read_BBpp) RETURN
			
			IF (SymVBBpp_read(vBBpp)) RETURN
			
			OPEN(file_desc, FILE=vBBpp%filename, ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				PRINT *, "*** Attention: Impossible to write in ", vBBpp%filename
			ELSE
				WRITE (file_desc, FMT="(E)", IOSTAT=file_error) vBBpp%b1
			END IF

		ELSE
		         
			Lold = Lmax
			Lmax = Lold - 1
			
			IF (Lmax < 10) THEN
				WRITE(vBBpp%filename, "(A,I1,A)") "data/vBB", Lmax, "pp_WS.txt"
			ELSE
				WRITE(vBBpp%filename, "(A,I2,A)") "data/vBB", Lmax, "pp_WS.txt"
			END IF
			
			IF (SymVBBpp_read(vBBpp)) write(*,'("Read file....")')
			
			Lmax = Lold
		
			IF (Lmax < 10) THEN
				WRITE(vBBpp%filename, "(A,I1,A)") "data/vBB", Lmax, "pp_WS.txt"
			ELSE
				WRITE(vBBpp%filename, "(A,I2,A)") "data/vBB", Lmax, "pp_WS.txt"
			END IF
			
		END IF

		! For each isospin, we calculate the constant numerical factors in front of the matrix elements 
		! names with the suffix "_same_part" refer to the terms with the delta_(ta,tb), names with the
		! suffix "_diff_part" refer to the terms with the other terms:
		!
		! Refs.: Page 132, top of th page, definition of v1pBB and v2pBB  	---  CHECKED AND OK
		!

		DO i = 0, 1
			x(i) = mu(i) / vBBpp%b1
			coef_pair1(i) = 0.5 * I_4PI * (Gogny_W(i, Gogny) - Gogny_B(i, Gogny) - Gogny_H(i, Gogny) + Gogny_M(i, Gogny))
			coef_pair2(i) =       I_4PI * (Gogny_W(i, Gogny) + Gogny_B(i, Gogny) - Gogny_H(i, Gogny) - Gogny_M(i, Gogny))
		END DO

                SELECT CASE (Basis)
		
		CASE (1)
         
			! Calculation of the matrix elements v1pBB and v2pBB as defined in Page 132  ---  CHECKED AND OK

			PRINT *, "Brink-Boker terms: Particle-Particle Channel - Harmonic Oscillator Basis"
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
										tres_j_cuad = (2*k + 1) * (ThreeJSymbols_get(2*la, 2*k, 2*lb) ** 2)
										cuad = CUAD2(la, lb, k)
										DO i = 0, 1
											total_ibb = IBrinkBookerHO(na - 1, la, nc - 1, la, nb - 1, lb, nd - 1, lb, k, x(i))
											IF(nb .NE. nd) THEN
												total_ibb = total_ibb + IBrinkBookerHO(na - 1, la, nc - 1, la, nd - 1, lb, nb - 1, lb, k, x(i))
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
									
									IF (file_error .EQ. 0) THEN
										WRITE (file_desc, FMT="(6I3,2E24.16)", IOSTAT=file_error) &
											la, na, nc, lb, nb, nd, sum1, sum2
									END IF
									
								END DO
							END DO
						END DO
					END DO
				END DO
				PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
			END DO
			
		CASE(2)
		
			PRINT *, "Brink-Boker terms: Particle-Particle Channel - General Basis"
			
			DO la = Lmin, Lmax
			
				icount = 0
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
										
										DO i = 0, 1
											total_ibb = IBrinkBooker(na, la, nc, la, nb, lb, nd, lb, k, mu(i))
											IF(nb .NE. nd) THEN
												total_ibb = total_ibb + IBrinkBooker(na, la, nc, la, nd, lb, nb, lb, k, mu(i))
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
									
									IF (file_error .EQ. 0) THEN
										WRITE (file_desc, FMT="(6I3,2E24.16)", IOSTAT=file_error) &
											la, na, nc, lb, nb, nd, sum1, sum2
									END IF
								END DO
							END DO
						END DO
					END DO
				END DO
				WRITE(*,'("Number of re-used integrals = ",I16)') icount
				WRITE(*,'("Completion rate: ",I3,"%")') INT(100 * (la + 1) / (N_0 + 1))
			END DO
			
		END SELECT
			
		! In case we want to optimize our calculation, we write the full matrix elements after 
		! having calculated only the few ones required
		
		IF (Optimization .EQ. 1 .AND. Lmin .NE. 0) THEN
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) THEN
				WRITE(*,'("Incompatible options: No optimization for HO case (or compatible)")')
				STOP 'Error in SymVBBpp_calculate - No optimization pssible'
			END IF
			CALL SymVBBpp_write(vBBpp)
		END IF
			
		CLOSE(file_desc)
		
		RETURN
	END SUBROUTINE SymVBBpp_calculate

        !---------------------------------------------------------------------------------------!
	!											!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVBBpp_write(vBBpp)
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp

		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd
		DOUBLE PRECISION :: sum1, sum2

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		OPEN(file_desc, FILE=vBBpp%filename, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to write in ", vBBpp%filename
		ELSE
			WRITE (file_desc, FMT="(E)", IOSTAT=file_error) vBBpp%b1
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
								
								CALL SymD3Tensor_SymD3Tensor_get(sum1, vBBpp%v1_pair, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum2, vBBpp%v2_pair, la, na, nc, lb, nb, nd)
								
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

	END SUBROUTINE SymVBBpp_write
	
	SUBROUTINE SymVBBpp_get_Delta(HF_out, vBBpp, P_in, ta, tb)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_out
		TYPE (SymVBBpp), INTENT(IN) :: vBBpp
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: P_in
		INTEGER, INTENT(IN) :: ta, tb

		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(tb), vBBpp%v1_pair, P_in%p(ta))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(tb), vBBpp%v2_pair, P_in%a(ta))
		
		HF_out%GaugeAngle(tb) = P_in%GaugeAngle(ta)

		RETURN
	END SUBROUTINE SymVBBpp_get_Delta

        !---------------------------------------------------------------------------------------!
	!											!
	!   Subroutine reading the matrix elements of the Brink-Boker term in the particle-	!
	!   particle channel. This subroutine attempts reading a file containing previously 	!
	!   calculated matrix elements. It is used to decide whether or not to recalculate 	!
	!   everything										!
	!											!
        !---------------------------------------------------------------------------------------!

	FUNCTION SymVBBpp_read(vBBpp)
		LOGICAL SymVBBpp_read
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp

		DOUBLE PRECISION b1
		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER i1, i2, i3, i4, i5, i6
		INTEGER, PARAMETER :: file_desc = 16
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

        !---------------------------------------------------------------------------------------!
	!											!
	!   Subroutine deleting the tensors corresponding to the matrix elements of the Brink-	!
	!   Boker in the particle-particle channel.						!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVBBpp_del(vBBpp)
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp

		CALL SymD3Tensor_SymD3Tensor_del(vBBpp%v1_pair)
		CALL SymD3Tensor_SymD3Tensor_del(vBBpp%v2_pair)
		
		RETURN
	END SUBROUTINE SymVBBpp_del

END MODULE symvbb
