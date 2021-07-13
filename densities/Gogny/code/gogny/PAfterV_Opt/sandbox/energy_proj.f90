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
		DOUBLE PRECISION, POINTER, DIMENSION(:, :) :: quantum
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

			Coefficients%x_l(ta)%index(i) = CMPLX(0.0,0.0)
			Coefficients%y_l(ta)%index(i) = CMPLX(0.0,0.0)
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

	SUBROUTINE SymGenDensityProj_make_DensityProj(genden_proj_real, genden_proj_imag, CoeffsXY, genden_gauge_real, &
	                                              genden_gauge_imag, UV, nucleus, Factor, NGauge, ta)
		TYPE (SymGenDensityProj), INTENT(OUT) :: genden_proj_real, genden_proj_imag
		TYPE (ProjectionCoeffs), INTENT(OUT) :: CoeffsXY
		TYPE (NucleusType), INTENT(IN) :: nucleus
		TYPE (MatrixType), POINTER, DIMENSION(:, :) :: UV
		TYPE (SymGenDensityGaugeProj), DIMENSION(:), POINTER :: genden_gauge_real, genden_gauge_imag
		DOUBLE PRECISION, INTENT(IN) :: Factor
		INTEGER, INTENT(IN) :: NGauge, ta

		COMPLEX, DIMENSION(:, :), ALLOCATABLE ::Identity
		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: Uphi, Vphi
		COMPLEX :: Trace, TraceBlock, DeterminantBlock, Angle, SumOfx

		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: U, V,  N11, rho

		DOUBLE PRECISION :: Gauge, det, Nexpe, pi, Signe, NewAngle, N0Block, N0
		DOUBLE PRECISION :: DBlockModule, DBlockAngle, DeterminantModule, DeterminantAngle

		INTEGER :: a, ja, d, ii, jj, Npart, IndexGauge, N, Z

		pi = 4.0*ATAN(1.0)

		! Main loop over the isospin (0 for protons, 1 for neutrons)

		SumOfx = 0.0

		! First loop over the Gauge angle to get the norm x_l = < phi | PN | phi >

		DO IndexGauge = 1, NGauge

			                        Gauge = Factor * pi * REAL(IndexGauge)/REAL(NGauge)
                        IF (NGauge .EQ. 1)      Gauge = 0.0

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

				N0Block = 0.0

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

				Uphi = RadPoles*cos(Gauge)*Identity - CMPLX(0,1)*RadPoles*sin(Gauge)*N11

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

				TraceBlock = CMPLX(0.0,0.0)

				DO ii = 1, d
					TraceBlock = TraceBlock + N11(ii,ii)
				END DO

				Trace = Trace + (ja + 1.0) * TraceBlock

				DEALLOCATE(Uphi)

				DEALLOCATE(U)
				DEALLOCATE(V)

				DEALLOCATE(N11)

				DEALLOCATE(Identity)

			END DO ! End of loop over a

			Nexpe = DBLE(nucleus%np(ta))

			Angle = Gauge * ((N0 - Nexpe) + 0.5*Trace)

			IF (ABS(AIMAG(Angle)) .GT. 1.e-10) THEN
				WRITE(*,'(/,"Angle = ",2F15.10)') Angle
				STOP "Complex Angle in SymGenDensityProj_make_DensityProj"
			END IF

			CoeffsXY%x_l(ta)%index(IndexGauge) = EXP(CMPLX(0,1)*(Angle + 0.5*DeterminantAngle)) * SQRT(DeterminantModule)
                        WRITE(*,'("Trace = ",2F15.10," DET = ",f20.10," Angle = ",2f20.10)') Trace, SQRT(DeterminantModule), Angle

			SumOfx = SumOfx + CoeffsXY%x_l(ta)%index(IndexGauge)

		END DO ! End of loop over IndexGauge

		! Normalizing the coefficients

		DO IndexGauge = 1, NGauge
			CoeffsXY%y_l(ta)%index(IndexGauge) = CoeffsXY%x_l(ta)%index(IndexGauge)/SumOfx
		        write(*,'("i = ",i2," y_l(i) = ",2f20.10)') IndexGauge, CoeffsXY%x_l(ta)%index(IndexGauge)
		END DO

		!write(*,'("ta = ",i2," (1/L) * Sum x_l = ",2f12.7," N = ",F10.5)') ta, SumOfX/REAL(NGauge),Nexpe

		! Calculating the projected density

		DO a = 0, 2*Lmax

			d = DIM(a)

			ALLOCATE(Identity(d, d))

			Identity = CMPLX(0,0)

			DO IndexGauge = 1, NGauge
				Identity = Identity + CoeffsXY%y_l(ta)%index(IndexGauge) * genden_gauge_real(IndexGauge)%rho%rho(ta, a)%store &
						    + CoeffsXY%y_l(ta)%index(IndexGauge) * genden_gauge_imag(IndexGauge)%rho%rho(ta, a)%store * CMPLX(0, 1)
			END DO

			genden_proj_real%rho%rho(ta, a)%store = REAL(Identity)
			genden_proj_imag%rho%rho(ta, a)%store = AIMAG(Identity)

			DEALLOCATE(Identity)

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

	SUBROUTINE SymGenDensityProj_make_DensityGauge(UV, genden_gauge_real, genden_gauge_imag, nucleus, Factor, NGauge, ta)
		TYPE (MatrixType), POINTER, DIMENSION(:, :) :: UV
		TYPE (NucleusType), INTENT(IN) :: nucleus
		DOUBLE PRECISION, INTENT(IN) :: Factor
		INTEGER, INTENT(IN) :: NGauge, ta

		TYPE (SymGenDensityGaugeProj), DIMENSION(:), POINTER :: genden_gauge_real, genden_gauge_imag

		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: Identity, S00, S11, S22, S33
		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: Uphi, Vphi, Aphi, UBar, VBar, UphiInv, TempMat
		COMPLEX :: Phase, U0, V0

		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: N11, N20, U, V
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: S0, S1, S2, S3

		DOUBLE PRECISION :: Gauge, det, pi

		INTEGER :: a, d, ii, jj, k, IndexGauge, Z, N

		pi = 4.0*ATAN(1.0)

		Z = nucleus%np(0)
		N = nucleus%np(1)

		DO IndexGauge = 1, NGauge

			                        Gauge = Factor * pi * REAL(IndexGauge)/REAL(NGauge)
                        IF (NGauge .EQ. 1)      Gauge = 0.0

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

				Uphi = RadPoles*cos(Gauge)*Identity + CMPLX(0,1)*RadPoles*sin(Gauge)*N11
				Vphi = 	                   - RadPoles*CMPLX(0,1)*sin(Gauge)*N20

				! Invert the (complex) Gauge-dependent U matrix

				ALLOCATE(UphiInv(d, d))

				CALL MatInvCmplx(d, 0, Uphi, UphiInv, det)

				! Generate the Thouless matrix A(phi) = [ V(phi) x Inv{U(phi)} ]*

				ALLOCATE(Aphi(d, d))

				Aphi = MATMUL(Vphi, UphiInv)

				! Calculate UBar and VBar

				ALLOCATE(Ubar(d, d))
				ALLOCATE(Vbar(d, d))

				Ubar = - U + MATMUL(V, CONJG(Aphi))
				Vbar = + V + MATMUL(U, CONJG(Aphi))

				DEALLOCATE(Aphi)

				! Calculate the Gauge-dependent densities rho and kappa

				ALLOCATE(genden_gauge_real(IndexGauge)%rho%rho(ta, a)%store(d, d))
				ALLOCATE(genden_gauge_real(IndexGauge)%kap01%rho(ta, a)%store(d, d))
				ALLOCATE(genden_gauge_real(IndexGauge)%kap10%rho(ta, a)%store(d, d))

				ALLOCATE(genden_gauge_imag(IndexGauge)%rho%rho(ta, a)%store(d, d))
				ALLOCATE(genden_gauge_imag(IndexGauge)%kap01%rho(ta, a)%store(d, d))
				ALLOCATE(genden_gauge_imag(IndexGauge)%kap10%rho(ta, a)%store(d, d))

				! Density matrix rho(phi)

				Phase = CMPLX(0,1)*Gauge

				ALLOCATE(TempMat(d, d))

				TempMat = MATMUL(CONJG(UPhiInv), TRANSPOSE(V))

				genden_gauge_real(IndexGauge)%rho%rho(ta, a)%store   = RadPoles * REAL (CEXP(Phase) * MATMUL(V, TempMat))
				genden_gauge_imag(IndexGauge)%rho%rho(ta, a)%store   = RadPoles * AIMAG(CEXP(Phase) * MATMUL(V, TempMat))

				DEALLOCATE(TempMat)

				! Pairing tensors kappa10(phi) and kappa01(phi)

				genden_gauge_real(IndexGauge)%kap10%rho(ta, a)%store = -REAL(MATMUL(Vbar, TRANSPOSE(U)))
				genden_gauge_real(IndexGauge)%kap01%rho(ta, a)%store = -REAL(MATMUL(Ubar, TRANSPOSE(V)))

				genden_gauge_real(IndexGauge)%rho%rho(1-ta, a)%store   = 0.0
				genden_gauge_real(IndexGauge)%kap10%rho(1-ta, a)%store = 0.0
				genden_gauge_real(IndexGauge)%kap01%rho(1-ta, a)%store = 0.0

				genden_gauge_imag(IndexGauge)%kap10%rho(ta, a)%store = -AIMAG(MATMUL(Vbar, TRANSPOSE(U)))
				genden_gauge_imag(IndexGauge)%kap01%rho(ta, a)%store = -AIMAG(MATMUL(Ubar, TRANSPOSE(V)))

				genden_gauge_imag(IndexGauge)%rho%rho(1-ta, a)%store   = 0.0
				genden_gauge_imag(IndexGauge)%kap10%rho(1-ta, a)%store = 0.0
				genden_gauge_imag(IndexGauge)%kap01%rho(1-ta, a)%store = 0.0

				DEALLOCATE(Uphi)
				DEALLOCATE(Vphi)

				DEALLOCATE(UphiInv)

				DEALLOCATE(Ubar)
				DEALLOCATE(Vbar)

				DEALLOCATE(N11)
				DEALLOCATE(N20)

				IF (NGauge .EQ. 1) THEN

					ALLOCATE(S0(d, d))
					ALLOCATE(S2(d, d))

					S0 = MATMUL(V, TRANSPOSE(V))
					S2 = MATMUL(V, TRANSPOSE(U))

					IF (nucleus%is_blocking(ta) .AND. (a .EQ. nucleus%ia(ta))) THEN

						ALLOCATE(S1(d, d))
						ALLOCATE(S3(d, d))

						DO ii = 1, d
							U0 = U(ii, nucleus%mu0(ta)) ! Column???
							V0 = V(ii, nucleus%mu0(ta))
							U(ii, nucleus%mu0(ta)) =   V0
							V(ii, nucleus%mu0(ta)) = - U0
						END DO

						S1 = MATMUL(V, TRANSPOSE(V))
						S3 = MATMUL(V, TRANSPOSE(U))
						S0 = S0 + (1.0 / (nucleus%ja(ta) + 1.0)) * (S1 - S0)
						S2 = S2 + (1.0 / (nucleus%ja(ta) + 1.0)) * (S3 - S2)

						DEALLOCATE(S1)
						DEALLOCATE(S3)

					END IF

					genden_gauge_real(IndexGauge)%rho%rho(ta, a)%store   = S0
					genden_gauge_imag(IndexGauge)%rho%rho(ta, a)%store   = 0.0

					genden_gauge_real(IndexGauge)%kap10%rho(ta, a)%store = S2
					genden_gauge_real(IndexGauge)%kap01%rho(ta, a)%store = S2

					genden_gauge_imag(IndexGauge)%kap10%rho(ta, a)%store = 0.0
					genden_gauge_imag(IndexGauge)%kap01%rho(ta, a)%store = 0.0

					DEALLOCATE(S0)
					DEALLOCATE(S2)

				END IF

				DEALLOCATE(U)
				DEALLOCATE(V)

				DEALLOCATE(Identity)

			END DO ! End of loop over a

		END DO ! End of loop over IndexGauge

		RETURN
	END SUBROUTINE SymGenDensityProj_make_DensityGauge

END MODULE energy_proj
