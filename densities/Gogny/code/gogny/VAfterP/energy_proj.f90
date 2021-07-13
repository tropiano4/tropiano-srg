 MODULE energy_proj

	USE input
	USE symd3t_proj
	USE nucleus
	USE deter
	USE symfield_proj
	USE symden_proj
	USE symgden_proj
	USE symgdhf_proj

	IMPLICIT NONE

	TYPE SymGenDensityGaugeProj
		TYPE (SymGenDensityHFProj) :: rho, kap10, kap01
	END TYPE

	TYPE MatrixType
		COMPLEX, POINTER, DIMENSION(:, :) :: quantum
	END TYPE

	TYPE ProjectionCoeffs
		COMPLEX, POINTER, DIMENSION(:, :) :: x_l
		COMPLEX, POINTER, DIMENSION(:, :) :: y_l
	END TYPE
	
 CONTAINS

	SUBROUTINE ProjectionCoeffs_new(Coefficients, NGauge)
		TYPE (ProjectionCoeffs), INTENT(INOUT) :: Coefficients
		INTEGER, INTENT(IN) :: NGauge
	
		INTEGER :: i, ta, la
	
		ALLOCATE(Coefficients%x_l(1:NGauge, 0:1))
		ALLOCATE(Coefficients%y_l(1:NGauge, 0:1))
		
		DO i = 1, NGauge
		
			DO ta = 0, 1
				Coefficients%x_l(i, ta) = CMPLX(0.0, 0.0)
				Coefficients%y_l(i, ta) = CMPLX(0.0, 0.0)
			END DO
		END DO
		
		RETURN
	END SUBROUTINE ProjectionCoeffs_new

	! Create a new density of the type SymDensityProj from an
	! old SymGenDensityProj
	!
	SUBROUTINE SymDensity_new_GenDensityProj10(density, genden)
		TYPE (SymDensityProj), INTENT(INOUT) :: density
		TYPE (SymGenDensityGaugeProj), INTENT(IN) :: genden

		CALL SymHartreeFockBogolFieldProj_new(density%field)
		
		density%field%rho = genden%rho
		density%field%kap = genden%kap10
		
		RETURN
	END SUBROUTINE SymDensity_new_GenDensityProj10

	! Create a new density of the type SymDensityProj from an
	! old SymGenDensityProj
	!
	SUBROUTINE SymDensity_new_GenDensityProj01(density, genden)
		TYPE (SymDensityProj), INTENT(INOUT) :: density
		TYPE (SymGenDensityGaugeProj), INTENT(IN) :: genden

		CALL SymHartreeFockBogolFieldProj_new(density%field)
		
		density%field%rho = genden%rho
		density%field%kap = genden%kap01
		
		RETURN
	END SUBROUTINE SymDensity_new_GenDensityProj01

	! Create a new density of the type SymDensityProj from an
	! old SymGenDensityProj
	!
	SUBROUTINE SymDensity_new_GenDensityProj_only_rho(density, genden)
		TYPE (SymDensityProj), INTENT(INOUT) :: density
		TYPE (SymGenDensityGaugeProj), INTENT(IN) :: genden

		CALL SymHartreeFockBogolFieldProj_new(density%field)
		
		density%field%rho = genden%rho
		
		RETURN
	END SUBROUTINE SymDensity_new_GenDensityProj_only_rho

	SUBROUTINE SymGenDensityGaugeProj_del(genden)
		TYPE (SymGenDensityGaugeProj), INTENT(INOUT) :: genden

		CALL SymGenDensityHFProj_del(genden%rho)
		CALL SymGenDensityHFProj_del(genden%kap10)
		CALL SymGenDensityHFProj_del(genden%kap01)
		
		RETURN
		
	END SUBROUTINE SymGenDensityGaugeProj_del

	SUBROUTINE ProjectionCoeffs_del(Coefficients)
		TYPE (ProjectionCoeffs), INTENT(INOUT) :: Coefficients

		DEALLOCATE(Coefficients%x_l)
		DEALLOCATE(Coefficients%y_l)
		
		RETURN
		
	END SUBROUTINE ProjectionCoeffs_del

	!-------------------------------------------------------------------------------!
	!										!
	!     This subroutine calculates the projected density genden_proj and the 	!
	!     projection coefficients y_l, stored in CoeffsXY, out of the matrices U 	!
	!     and V. To compute the projected density, it requires the "Gauge" density	!
	!     genden_gauge, which corresponds to rho(phi) (and also kappa01(phi) and	!
	!     kappa10(phi)).								!
	!										!
	!-------------------------------------------------------------------------------!

	SUBROUTINE SymGenDensityProj_make_DensityProj(genden_proj, CoeffsXY, genden_gauge, UV, nucleus, NGauge)
		TYPE (SymGenDensityProj), INTENT(OUT) :: genden_proj
		TYPE (ProjectionCoeffs), INTENT(OUT) :: CoeffsXY
		TYPE (NucleusType), INTENT(IN) :: nucleus
		TYPE (MatrixType), POINTER, DIMENSION(:, :) :: UV
		TYPE (SymGenDensityGaugeProj), DIMENSION(:), POINTER :: genden_gauge
		INTEGER, INTENT(IN) :: NGauge

		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: U, V, N11, rho, Identity
		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: Uphi, Vphi
		COMPLEX :: Trace, TraceBlock, DeterminantBlock, Angle, SumOfx, N0Block, N0
		
		DOUBLE PRECISION :: Gauge, det, Nexpe, pi
		DOUBLE PRECISION :: DBlockModule, DBlockAngle, DeterminantModule, DeterminantAngle
		
		INTEGER :: ta, a, ja, d, ii, jj, Npart, IndexGauge, N, Z

		pi = 4.0*ATAN(1.0)
		
		Z = nucleus%np(0)
		N = nucleus%np(1)
		
		CALL SymGenDensityProj_new_Nucleus(genden_proj, N, Z)
		
		! Main loop over the isospin (0 for protons, 1 for neutrons)
		
		DO ta = 0, 1
		
			SumOfx = 0.0
			
			! First loop over the Gauge angle to get the norm x_l = < phi | PN | phi >
		
			DO IndexGauge = 1, NGauge
		
				Gauge = pi * REAL(IndexGauge)/REAL(NGauge)
		
				Trace = CMPLX(0,0)
				
				N0 = CMPLX(0,0)
				
				DeterminantModule = 1.0
				DeterminantAngle  = 0.0
								
				DO a = 0, 2*Lmax
			
					d = DIM(a)
					
					ja = J(a)

					ALLOCATE(U(d, d))
					ALLOCATE(V(d, d))
		
					! We extract the matrices U and V
					
					DO ii = 1, d
						DO jj = 1, d
							U(ii, jj) = UV(ta, a)%quantum(ii    , jj + d)
							V(ii, jj) = UV(ta, a)%quantum(ii + d, jj + d)
						END DO
					END DO

					! Calculating the expectation value N0 = Tr(rho). Each state "a" is 
					! degenerate 2ja+1 so we must multiply the final trace by this number
					! (ja contains in fact the value 2*ja)
		
					ALLOCATE(rho(d, d))
					
					rho = MATMUL(V, TRANSPOSE(V))
					
					N0Block = CMPLX(0,0)
					
					DO ii = 1, d
						N0Block = N0Block + rho(ii, ii)
					END DO
					
					N0 = N0 + (ja + 1.0)*N0Block
					
					DEALLOCATE(rho)
		
					! Construct Norm and density matrix. 
		
					ALLOCATE(N11(d, d))
		
					N11 = MATMUL(TRANSPOSE(U), U) - MATMUL(TRANSPOSE(V), V)
		
					ALLOCATE(Identity(d, d))
					
					Identity = CMPLX(0,0)
					
					DO ii = 1, d
						Identity(ii,ii) = CMPLX(1,0)
					END DO
							
					ALLOCATE(Uphi(d, d))
		
					Uphi = cos(Gauge)*Identity - CMPLX(0,1)*sin(Gauge)*N11

					! The determinant of the block is the calculated value power
					! the degeneracy of the block (2j+1). We have to be very careful,
					! sicne allt hese numbers are complex
		
					DeterminantBlock = DETERM(Uphi, d)
					
					DBlockModule = ABS(DeterminantBlock)
					DBlockAngle  = ATAN2(AIMAG(DeterminantBlock), REAL(DeterminantBlock))
															
					DeterminantModule = DeterminantModule * (DBlockModule**(ja+1))
					DeterminantAngle  = DeterminantAngle + (ja + 1)*DBlockAngle
		
					! Calculating the trace of the matrix CONJG(N11). Same remark as for the 
					! density matrix: the degeneracy of the block must be taken into
					! account.
					
					TraceBlock = CMPLX(0,0)
					
					DO ii = 1, d
						TraceBlock = TraceBlock + CONJG(N11(ii,ii))
					END DO
					
					Trace = Trace + (ja + 1.0) * TraceBlock
							
					DEALLOCATE(Uphi)
			
					DEALLOCATE(U)
					DEALLOCATE(V)
		
					DEALLOCATE(N11)
		
					DEALLOCATE(Identity)
				
				END DO ! End of loop over a
				
				Nexpe = REAL(nucleus%np(ta))
				
				Angle = Gauge * ((N0 - Nexpe) + 0.5*Trace)
				
				IF (ABS(AIMAG(Angle)) .GT. 1.e-10) THEN
					WRITE(*,'(/,"Angle = ",2F15.10)') Angle
					STOP "Complex Angle in SymGenDensityProj_make_DensityProj"
				END IF
						
				CoeffsXY%x_l(IndexGauge, ta) = CEXP(CMPLX(0,1)*(Angle + 0.5*DeterminantAngle)) * SQRT(DeterminantModule)
							
				SumOfx = SumOfx + CoeffsXY%x_l(IndexGauge, ta)
				
			END DO ! End of loop over IndexGauge
			
			! Normalizing the coefficients
			
			DO IndexGauge = 1, NGauge
				CoeffsXY%y_l(IndexGauge, ta) = CoeffsXY%x_l(IndexGauge, ta)/SumOfx
			END DO
			
			write(*,'("ta = ",i2," (1/L) * Sum x_l = ",2f12.7)') ta, SumOfX/REAL(NGauge)
						
			! Calculating the projected density
			
			DO a = 0, 2*Lmax
			
				d = DIM(a)
		
				ALLOCATE(rho(d, d))
					
				rho = CMPLX(0,0)

				DO IndexGauge = 1, NGauge
					rho = rho + CoeffsXY%y_l(IndexGauge, ta)*genden_gauge(IndexGauge)%rho%rho(ta, a)%store
				END DO
					
				genden_proj%rho%rho(ta, a)%store = rho
				
				DEALLOCATE(rho)

			END DO
		END DO
		
		RETURN
	END SUBROUTINE SymGenDensityProj_make_DensityProj
	
	!-------------------------------------------------------------------------------!
	!										!
	!     This subroutine calculates the Gauge densities rho(phi, kappa10(phi) and	!
	!     kappa01(phi) according to equations (9), (10) and (11) of:		!
	!										!
	!				M. Anguiano, J. L. Egido and L. M. Robledo	!
	!				Nucl. Phys. A696 (2001) 467-493			!
	!										!
	!-------------------------------------------------------------------------------!

	SUBROUTINE SymGenDensityProj_make_DensityGauge(UV, genden_gauge, NGauge)
		TYPE (MatrixType), POINTER, DIMENSION(:, :) :: UV
		INTEGER, INTENT(IN) :: NGauge

		TYPE (SymGenDensityGaugeProj), DIMENSION(:), POINTER :: genden_gauge
		
		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: N11, N20, Identity
		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: Uphi, Vphi, Aphi, UBar, VBar, UphiInv
		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: U, V, TempMat
		COMPLEX :: Phase
		
		DOUBLE PRECISION :: Gauge, det, pi
		
		INTEGER :: ta, a, d, ii, jj, k, IndexGauge

		pi = 4.0*ATAN(1.0)
		
		ALLOCATE(genden_gauge(NGauge))
		
		DO IndexGauge = 1, NGauge
			ALLOCATE(genden_gauge(IndexGauge)%rho%rho(0:1, 0:2*Lmax))
			ALLOCATE(genden_gauge(IndexGauge)%kap01%rho(0:1, 0:2*Lmax))
			ALLOCATE(genden_gauge(IndexGauge)%kap10%rho(0:1, 0:2*Lmax))
		END DO
					
		DO ta = 0, 1
		
			DO IndexGauge = 1, NGauge
						
				Gauge = pi*REAL(IndexGauge)/REAL(NGauge)
				
				DO a = 0, 2*Lmax
			
					d = DIM(a)
					
					ALLOCATE(U(d, d))
					ALLOCATE(V(d, d))
		
					! We extract the matrices U and V
					DO ii = 1, d
						DO jj = 1, d
							U(ii, jj) = UV(ta, a)%quantum(ii    , jj + d)
							V(ii, jj) = UV(ta, a)%quantum(ii + d, jj + d)
						END DO
					END DO

					!-----------------------------------------------------------------------!
					!									!
					!         Construct projected densities					!
					!									!
					!-----------------------------------------------------------------------!
		
					! Construct Norm matrices
		
					ALLOCATE(N11(d, d))
					ALLOCATE(N20(d, d))
		
					! Construct Norm and density matrix. The norm matrix has a different
					! expression as in the generic theory, as we consider (j,l,m) blocks
		
					N11 = MATMUL(TRANSPOSE(U), U) - MATMUL(TRANSPOSE(V), V)
					N20 = MATMUL(TRANSPOSE(U), V) + MATMUL(TRANSPOSE(V), U)
		
					ALLOCATE(Identity(d, d))
					
					Identity = CMPLX(0,0)
		
					DO ii = 1, d
						Identity(ii,ii) = 1.0
					END DO
		
					ALLOCATE(Uphi(d, d))
					ALLOCATE(Vphi(d, d))
		
					Uphi = cos(Gauge)*Identity + CMPLX(0,1)*sin(Gauge)*N11
					Vphi = 	                   - CMPLX(0,1)*sin(Gauge)*N20
		
					! Invert the (complex) Gauge-dependent U matrix
		
					ALLOCATE(UphiInv(d, d))
					
					CALL MatInvCmplx(d, 0, Uphi, UphiInv, det)
					
					! Generate the Thouless matrix A(phi) = [ V(phi) x Inv{U(phi)} ]*
		
					ALLOCATE(Aphi(d, d))
		
					Aphi = MATMUL(Vphi, UphiInv)

					! Calculate UBar and VBar
		
					ALLOCATE(Ubar(d, d))
					ALLOCATE(Vbar(d, d))
		
					Ubar = - CONJG(U) + MATMUL(V, CONJG(Aphi))
					Vbar = + CONJG(V) + MATMUL(U, CONJG(Aphi))
					
					DEALLOCATE(Aphi)
					
					! Calculate the Gauge-dependent densities rho and kappa
			
					ALLOCATE(genden_gauge(IndexGauge)%rho%rho(ta, a)%store(d, d))
					ALLOCATE(genden_gauge(IndexGauge)%kap01%rho(ta, a)%store(d, d))
					ALLOCATE(genden_gauge(IndexGauge)%kap10%rho(ta, a)%store(d, d))
					
					! Density matrix rho(phi)
					
					Phase = CMPLX(0,1)*Gauge
					
					ALLOCATE(TempMat(d, d))
					
					TempMat = MATMUL(CONJG(UPhiInv), TRANSPOSE(V))
					
					genden_gauge(IndexGauge)%rho%rho(ta, a)%store   = CEXP(Phase) * MATMUL(CONJG(V), TempMat)
					
					DEALLOCATE(TempMat)
					
					! Pairing tensors kappa10(phi) and kappa01(phi)
					
					!genden_gauge(IndexGauge)%kap10%rho(ta, a)%store = -MATMUL(Vbar, TRANSPOSE(U))
					!genden_gauge(IndexGauge)%kap01%rho(ta, a)%store = -MATMUL(Ubar, TRANSPOSE(V))
					
					Phase = CMPLX(0,1)*Gauge
					
					ALLOCATE(TempMat(d, d))
					
					TempMat = MATMUL(CONJG(UPhiInv), TRANSPOSE(U))
					
					genden_gauge(IndexGauge)%kap10%rho(ta, a)%store  = CEXP(Phase) * MATMUL(CONJG(V), TempMat)
					
					DEALLOCATE(TempMat)
					
					Phase = -CMPLX(0,1)*Gauge
					
					ALLOCATE(TempMat(d, d))
					
					TempMat = MATMUL(CONJG(UPhiInv), TRANSPOSE(V))
					
					genden_gauge(IndexGauge)%kap01%rho(ta, a)%store  = - CEXP(Phase) * MATMUL(CONJG(U), TempMat)
					
					DEALLOCATE(TempMat)
					
					DEALLOCATE(Uphi)
					DEALLOCATE(Vphi)
		
					DEALLOCATE(UphiInv)
					
					DEALLOCATE(Ubar)
					DEALLOCATE(Vbar)
				
					DEALLOCATE(N11)
					DEALLOCATE(N20)
		
					DEALLOCATE(U)
					DEALLOCATE(V)
		
					DEALLOCATE(Identity)
		
				END DO ! End of loop over a
		
			END DO ! End of loop over IndexGauge
					
		END DO ! End of loop over ta
		
		RETURN
	END SUBROUTINE SymGenDensityProj_make_DensityGauge

	!-------------------------------------------------------------------------------!
	!										!
	!   			CONSTRUCTING THE C-MATRIX				!
	!										!
	!-------------------------------------------------------------------------------!

	SUBROUTINE ProjectionPreparation(UV, CoeffsXY, genden_gauge, C_Matrix, Deriv_Y, NGauge)
		TYPE (MatrixType), POINTER, DIMENSION(:, :) :: UV
		TYPE (ProjectionCoeffs), INTENT(IN) :: CoeffsXY
		INTEGER, INTENT(IN) :: NGauge
		
		TYPE (SymGenDensityGaugeProj), DIMENSION(:), POINTER :: genden_gauge

		TYPE (SymGenDensityHFProj), DIMENSION(:), POINTER :: C_Matrix, Deriv_Y
		TYPE (SymGenDensityHFProj) :: gendenhf
		
		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: U, V, rho, TempMat, Identity
		COMPLEX :: Phase
		
		DOUBLE PRECISION :: Gauge, det, pi
		
		INTEGER :: la, ta, a, d, ii, jj, k, IndexGauge

		pi = 4.0*ATAN(1.0)
		
		ALLOCATE(C_Matrix(NGauge))
		ALLOCATE(Deriv_Y(NGauge))
		
		DO IndexGauge = 1, NGauge
			CALL SymGenDensityHFProj_new(C_Matrix(IndexGauge))
			CALL SymGenDensityHFProj_new(Deriv_Y(IndexGauge))
		END DO
			
		! IN THIS LOOP, WE CALCULATE THE C-MATRIX
				
		CALL SymGenDensityHFProj_new(gendenhf)

		DO ta = 0, 1
		
			! Get the density from the U and V vectors
			DO a = 0, 2*Lmax
			
				d = DIM(a)
					
				ALLOCATE(U(d, d))
				ALLOCATE(V(d, d))
		
				! We extract the matrices U and V
				DO ii = 1, d
					DO jj = 1, d
						U(ii, jj) = UV(ta, a)%quantum(ii    , jj + d)
						V(ii, jj) = UV(ta, a)%quantum(ii + d, jj + d)
					END DO
				END DO
					
				gendenhf%rho(ta, a)%store = MATMUL(V, TRANSPOSE(V))

				DEALLOCATE(U)
				DEALLOCATE(V)

			END DO
				
		END DO
					
			
		DO ta = 0, 1
		
			DO IndexGauge = 1, NGauge
			
				Gauge = pi*REAL(IndexGauge)/REAL(NGauge)
				
				DO a = 0, 2*Lmax
					
					d = DIM(a)
					
					ALLOCATE(Identity(d, d))
					
					Identity = CMPLX(0,0)
		
					DO ii = 1, d
						Identity(ii,ii) = 1.0
					END DO
							
					ALLOCATE(rho(d, d))
					ALLOCATE(TempMat(d, d))
		
					rho = gendenhf%rho(ta, a)%store
					
					Phase = CEXP(CMPLX(0,1)*2.0*Gauge)
					
					TempMat = Identity + rho*(Phase - 1.0)
					
					CALL MatInvCmplx(d, 0, TempMat, rho, det)
					
					C_Matrix(IndexGauge)%rho(ta, a)%store = Phase*rho

					DEALLOCATE(rho)
					DEALLOCATE(TempMat)
					
					DEALLOCATE(Identity)
					
				END DO ! End of loop over a
				
			END DO ! End of loop over IndexGauge
					
		END DO ! End of loop over ta
		
		!DO IndexGauge = 1, NGauge
			
		!	DO ta = 0, 1
		!		DO a = 0, 2*Lmax
		!			genden_gauge(IndexGauge)%rho%rho(ta, a)%store = &
		!					MATMUL(C_Matrix(IndexGauge)%rho(ta, a)%store, gendenhf%rho(ta, a)%store)
		!		END DO
		!	END DO
		!		
		!END DO
					
					
		! IN THIS LOOP, WE CALCULATE THE DERIVATIVE OF THE Y-COEFFICIENTS WITH RESPECT TO RHO
				
		DO ta = 0, 1
		
			DO a = 0, 2*Lmax
					
				d = DIM(a)
					
				ALLOCATE(TempMat(d, d))
				
				TempMat = CMPLX(0,0)
				
				DO IndexGauge = 1, NGauge
						
					Gauge = pi*REAL(IndexGauge)/REAL(NGauge)
					
					Phase = 1.0 - CEXP(-CMPLX(0,1)*2.0*Gauge)
					
					TempMat = TempMat + Phase * CoeffsXY%y_l(IndexGauge, ta) * C_Matrix(IndexGauge)%rho(ta, a)%store
					
				END DO
								
				DO IndexGauge = 1, NGauge
						
					Gauge = pi*REAL(IndexGauge)/REAL(NGauge)
					
					Phase = 1.0 - CEXP(-CMPLX(0,1)*2.0*Gauge)
					
					Deriv_Y(IndexGauge)%rho(ta, a)%store &
						= 0.5*CoeffsXY%y_l(IndexGauge, ta) * Phase * C_Matrix(IndexGauge)%rho(ta, a)%store &
						- 0.5*CoeffsXY%y_l(IndexGauge, ta) * TempMat
										
				END DO
		
				DEALLOCATE(TempMat)
				
			END DO ! End of loop over a
					
		END DO ! End of loop over ta
		
		RETURN
	END SUBROUTINE ProjectionPreparation
	
	!-------------------------------------------------------------------------------!
	!										!
	!   				!
	!										!
	!-------------------------------------------------------------------------------!

	SUBROUTINE SymEkTensorVAP_get_GammaProj(HF_out, EkTensor, genden_gauge, C_Matrix, Deriv_Y, CoeffsXY, NGauge)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_out
		!TYPE (SymD3Tensor), INTENT(INOUT) :: HF_out
		TYPE (SymD3Tensor), INTENT(IN) :: EkTensor
		TYPE (ProjectionCoeffs), INTENT(IN) :: CoeffsXY
		
		TYPE (SymGenDensityGaugeProj), DIMENSION(:), POINTER :: genden_gauge
		TYPE (SymGenDensityHFProj), DIMENSION(:), POINTER :: C_Matrix, Deriv_Y
		
		INTEGER, INTENT(IN) :: NGauge
		
		TYPE (SymHartreeFockFieldProj) :: dY_field, C_field, rho_field, rho_degen
		TYPE (SymHartreeFockFieldProj) :: field1, field2, field3, field4, Ekfield
		TYPE (SymDensityProj) :: density_rho

		COMPLEX, DIMENSION(1:Nmax, 1:Nmax) :: M, MT
		
		COMPLEX :: phase, Trace, factor
		DOUBLE PRECISION :: AngleGauge, pi
		INTEGER :: IndexGauge, ta, la, d
		
		CALL SymHartreeFockFieldProj_new(dY_field)
		CALL SymHartreeFockFieldProj_new(C_field)
		CALL SymHartreeFockFieldProj_new(rho_field)
		CALL SymHartreeFockFieldProj_new(rho_degen)

		CALL SymHartreeFockFieldProj_new(Ekfield)
		
		DO ta = 0,1
			Ekfield%p(ta) = EkTensor
			Ekfield%a(ta) = EkTensor
		END DO
		
		CALL SymHartreeFockFieldProj_new(field1)
		CALL SymHartreeFockFieldProj_new(field2)
		CALL SymHartreeFockFieldProj_new(field3)
		CALL SymHartreeFockFieldProj_new(field4)
		
		pi = 4.0*ATAN(1.0)

		DO IndexGauge = 1, NGauge
		
			AngleGauge = pi*REAL(IndexGauge)/REAL(NGauge)
			
			CALL SymGenDensityHFProj_assign_VAP(dY_field, Deriv_Y(IndexGauge), 0)
			!CALL SymGenDensityHFProj_assign_VAP(C_field, C_Matrix(IndexGauge), 0)
			CALL SymGenDensityHFProj_assign2(C_field, C_Matrix(IndexGauge))

			CALL SymGenDensityHFProj_assign_VAP(rho_field, genden_gauge(IndexGauge)%rho, 0)
			CALL SymGenDensityHFProj_assign_VAP(rho_degen, genden_gauge(IndexGauge)%rho, 1)
			
			DO ta = 0, 1
		
				! Second Term: dy/drho * Tr[t*rho]

				Trace = SymHartreeFockFieldProj_product_iso(Ekfield, rho_degen, ta)
		
				CALL SymHartreeFockFieldProj_scalar_iso(field1, Trace, dY_field, ta)
		
				! First Term
				
				CALL SymHartreeFockFieldProj_full_product_iso(field2, Ekfield, C_field, ta, ta)
				CALL SymHartreeFockFieldProj_full_product_iso(field3, rho_field, field2, ta, ta)
				
				phase = - 1.0 + CEXP(-CMPLX(0,1)*2.0*AngleGauge)

				CALL SymHartreeFockFieldProj_scalar_iso(field4, phase, field3, ta)
		
				CALL SymHartreeFockFieldProjIso_add(field3, field2, field4, ta)
		
				CALL SymHartreeFockFieldProj_scalar_iso(field4, CoeffsXY%y_l(IndexGauge, ta), field3, ta)
		
				CALL SymHartreeFockFieldProjIso_add(field2, field1, field4, ta)
				
				CALL SymHartreeFockFieldProjIso_add(HF_out, HF_out, field2, ta)
				
			END DO
			
		END DO
		
		
		!DO ta = 0, 1
		!	DO la = 0, Lmax
		!				factor = 0.0
		!		IF (la .NE. 0)  factor = 1.0/LS(2*la - 1)
		!		
		!		CALL SymD3Tensor_product(field1%a(ta), factor*CMPLX(1,0), HF_out%a(ta))
		!		CALL SymD3Tensor_assign1(HF_out%a(ta), field1%a(ta))
		!
		!	END DO
		!END DO
		
		!CALL SymHartreeFockFieldProj_add(field3, field2, field1)
		!CALL SymHartreeFockFieldProj_product(HF_out, CMPLX(1,0)*0.5, field3)
		
		CALL SymHartreeFockFieldProj_del(field1)
		CALL SymHartreeFockFieldProj_del(field2)
		CALL SymHartreeFockFieldProj_del(field3)
		CALL SymHartreeFockFieldProj_del(field4)
		
		CALL SymHartreeFockFieldProj_del(dY_field)
		CALL SymHartreeFockFieldProj_del(C_field)
		CALL SymHartreeFockFieldProj_del(rho_field)
		CALL SymHartreeFockFieldProj_del(rho_degen)

		RETURN
	END SUBROUTINE SymEkTensorVAP_get_GammaProj

END MODULE energy_proj
