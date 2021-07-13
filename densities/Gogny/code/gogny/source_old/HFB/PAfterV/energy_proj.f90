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

	TYPE Vector
		COMPLEX, POINTER, DIMENSION(:) :: index
	END TYPE
	
	TYPE ProjectionCoeffs
		TYPE (Vector), POINTER, DIMENSION(:) :: x_l
		TYPE (Vector), POINTER, DIMENSION(:) :: y_l
	END TYPE
	
 CONTAINS

	SUBROUTINE ProjectionCoeffs_new(Coefficients, NGauge, ta)
		TYPE (ProjectionCoeffs), INTENT(INOUT) :: Coefficients
		INTEGER, INTENT(IN) :: NGauge, ta
	
		INTEGER :: i
	
		ALLOCATE(Coefficients%x_l(ta)%index(1:NGauge))
		ALLOCATE(Coefficients%y_l(ta)%index(1:NGauge))
		
		DO i = 1, NGauge
		
			Coefficients%x_l(ta)%index(i) = CMPLX(0,0)
			Coefficients%y_l(ta)%index(i) = CMPLX(0,0)
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

	SUBROUTINE SymGenDensityProj_make_DensityProj(genden_proj, CoeffsXY, genden_gauge, UV, nucleus, Factor, NGauge, ta)
		TYPE (SymGenDensityProj), INTENT(OUT) :: genden_proj
		TYPE (ProjectionCoeffs), INTENT(OUT) :: CoeffsXY
		TYPE (NucleusType), INTENT(IN) :: nucleus
		TYPE (MatrixType), POINTER, DIMENSION(:, :) :: UV
		TYPE (SymGenDensityGaugeProj), DIMENSION(:), POINTER :: genden_gauge
		DOUBLE PRECISION, INTENT(IN) :: Factor
		INTEGER, INTENT(IN) :: NGauge, ta

		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: U, V, N11, rho, Identity
		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: Uphi, Vphi
		COMPLEX :: Trace, TraceBlock, DeterminantBlock, Angle, SumOfx, N0Block, N0
		
		DOUBLE PRECISION :: Gauge, det, Nexpe, pi, Signe, NewAngle
		DOUBLE PRECISION :: DBlockModule, DBlockAngle, DeterminantModule, DeterminantAngle
		
		INTEGER :: a, ja, d, ii, jj, Npart, IndexGauge, N, Z

		pi = 4.0*ATAN(1.0)
		
		! Main loop over the isospin (0 for protons, 1 for neutrons)
		
		SumOfx = 0.0
			
		! First loop over the Gauge angle to get the norm x_l = < phi | PN | phi >
		
		DO IndexGauge = 1, NGauge
		
			Gauge = Factor * pi * REAL(IndexGauge)/REAL(NGauge)
		
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
				! since all these numbers are complex
		
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
				
			Nexpe = nucleus%actual_np(ta)
			
			Angle = Gauge * ((N0 - Nexpe) + 0.5*Trace)
				
			IF (ABS(AIMAG(Angle)) .GT. 1.e-10) THEN
				WRITE(*,'(/,"Angle = ",2F15.10)') Angle
				STOP "Complex Angle in SymGenDensityProj_make_DensityProj"
			END IF
			
			CoeffsXY%x_l(ta)%index(IndexGauge) = CEXP(CMPLX(0,1)*(Angle + 0.5*DeterminantAngle)) * SQRT(DeterminantModule)
!	WRITE(*,'("Trace = ",2F15.10," DET = ",f20.10," Angle = ",2f20.10)') Trace, SQRT(DeterminantModule), DeterminantAngle
							
			SumOfx = SumOfx + CoeffsXY%x_l(ta)%index(IndexGauge)
				
		END DO ! End of loop over IndexGauge
			
		! Normalizing the coefficients
			
		DO IndexGauge = 1, NGauge
			CoeffsXY%y_l(ta)%index(IndexGauge) = CoeffsXY%x_l(ta)%index(IndexGauge)/SumOfx
!		write(*,'("i = ",i2," y_l(i) = ",2f20.10)') IndexGauge, CoeffsXY%x_l(ta)%index(IndexGauge)
		END DO
			
		write(*,'("ta = ",i2," (1/L) * Sum x_l = ",2f12.7," N = ",F10.5)') ta, SumOfX/REAL(NGauge),Nexpe
		
		! Calculating the projected density
			
		DO a = 0, 2*Lmax
			
			d = DIM(a)
		
			ALLOCATE(rho(d, d))
					
			rho = CMPLX(0,0)

			DO IndexGauge = 1, NGauge
				rho = rho + CoeffsXY%y_l(ta)%index(IndexGauge) * genden_gauge(IndexGauge)%rho%rho(ta, a)%store
			END DO
					
			genden_proj%rho%rho(ta, a)%store = rho
				
			DEALLOCATE(rho)

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

	SUBROUTINE SymGenDensityProj_make_DensityGauge(UV, genden_gauge, nucleus, Factor, NGauge, ta)
		TYPE (MatrixType), POINTER, DIMENSION(:, :) :: UV
		TYPE (NucleusType), INTENT(IN) :: nucleus
		DOUBLE PRECISION, INTENT(IN) :: Factor
		INTEGER, INTENT(IN) :: NGauge, ta

		TYPE (SymGenDensityGaugeProj), DIMENSION(:), POINTER :: genden_gauge
		
		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: N11, N20, Identity, S0, S1, S2, S3, S00, S11, S22, S33
		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: Uphi, Vphi, Aphi, UBar, VBar, UphiInv
		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: U, V, TempMat
		COMPLEX :: Phase, U0, V0
		
		DOUBLE PRECISION :: Gauge, det, pi
		
		INTEGER :: a, d, ii, jj, k, IndexGauge, Z, N

		pi = 4.0*ATAN(1.0)
		
		Z = nucleus%np(0)
		N = nucleus%np(1)		
		
		DO IndexGauge = 1, NGauge
						
			Gauge = Factor * pi * REAL(IndexGauge)/REAL(NGauge)
				
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
					
				genden_gauge(IndexGauge)%kap10%rho(ta, a)%store = -MATMUL(Vbar, TRANSPOSE(U))
				genden_gauge(IndexGauge)%kap01%rho(ta, a)%store = -MATMUL(Ubar, TRANSPOSE(V))
				
				genden_gauge(IndexGauge)%rho%rho(1-ta, a)%store   = CMPLX(0,0)
				genden_gauge(IndexGauge)%kap10%rho(1-ta, a)%store = CMPLX(0,0)
				genden_gauge(IndexGauge)%kap01%rho(1-ta, a)%store = CMPLX(0,0)
					
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
					
		RETURN
	END SUBROUTINE SymGenDensityProj_make_DensityGauge

END MODULE energy_proj
