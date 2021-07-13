 MODULE diagmeth

	USE input
	USE symd3t_proj
	USE r1r1
	USE nucleus
	USE symden_proj
	USE symgden_proj
	USE symgdhf_proj
	USE symke2b
	USE symvbb
	USE symvc
	USE symvls
	USE symgdd_proj
	USE selfc_proj
	USE wave
	USE indexx
	USE energy_proj

	IMPLICIT NONE

	INCLUDE "brent_d1.f90"

	INTEGER, PARAMETER :: MAX_ITER = 5
	
	INTEGER :: NGauge = 9
	
	DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DensityRadial

	TYPE DiagonalMatrix
		COMPLEX, POINTER, DIMENSION(:) :: value
	END TYPE

	TYPE DiagonalizationMethod
		TYPE (SelfConsistencyMethodProj) consistency
		TYPE (R1R1Function) func
		TYPE (MatrixType), POINTER, DIMENSION(:, :) :: UV
		TYPE (DiagonalMatrix), POINTER, DIMENSION(:, :) :: QuasiParticleEnergies
		TYPE (SymGenDensityProj) S ! Superhamiltoniano sin el multiplicador del radio
		TYPE (SymGenDensityProj) iterated ! Density
		INTEGER ta
	END TYPE
	
        TYPE (ProjectionCoeffs) :: CoeffsXY

 CONTAINS

	INCLUDE "brent_d2.f90"

	SUBROUTINE DiagonalizationMethod_new(diagonal, density)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		TYPE (SymDensityProj), INTENT(INOUT) :: density

		INTEGER :: ta, a, max_a, d, dd

		! Inicializamos las variables y subtipos
		CALL SelfConsistencyMethodProj_new(diagonal%consistency, density)
		
		CALL R1R1Function_new(diagonal%func)
		
		CALL SymGenDensityProj_new_Nucleus(diagonal%S, Nucleus_get_N(density%nucleus), Nucleus_get_Z(density%nucleus))
		CALL SymGenDensityProj_new_SymDensity(diagonal%iterated, density)
		
		max_a = 2*Lmax
		
		! Reservamos memoria para las matrices del tipo
		ALLOCATE(diagonal%QuasiParticleEnergies(0:1, 0:max_a))
		IF (.NOT. ASSOCIATED(diagonal%QuasiParticleEnergies)) STOP "Unable to allocate memory in DiagonalizationMethod_new (1)"
		
		ALLOCATE(diagonal%UV(0:1, 0:max_a))
		IF (.NOT. ASSOCIATED(diagonal%UV)) STOP "Unable to allocate memory in DiagonalizationMethod_new (2)"
		
		DO ta = 0, 1
			DO a = 0, max_a
				
				d = DIM(a)
				dd = d + d
				
				ALLOCATE(diagonal%QuasiParticleEnergies(ta, a)%value(dd))
				IF (.NOT. ASSOCIATED(diagonal%QuasiParticleEnergies(ta, a)%value)) STOP "Unable to allocate memory in DiagonalizationMethod_new (3)"
				
				ALLOCATE(diagonal%UV(ta, a)%quantum(dd, dd))
				IF (.NOT. ASSOCIATED(diagonal%UV(ta, a)%quantum)) STOP "Unable to allocate memory in DiagonalizationMethod_new (4)"
				
			END DO
		END DO
		
		RETURN
	END SUBROUTINE DiagonalizationMethod_new

	SUBROUTINE DiagonalizationMethod_get_SymDensity(diagonal, density)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		TYPE (SymDensityProj), INTENT(INOUT) :: density

		CALL SymDensity_new_GenDensityProj(density, diagonal%iterated)
		
		RETURN
	END SUBROUTINE DiagonalizationMethod_get_SymDensity

	SUBROUTINE DiagonalizationMethod_goto_SelfConsistency(diagonal, tolerance)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		DOUBLE PRECISION, INTENT(IN) :: tolerance

		TYPE (OneDimSolve) :: neutron_constrainer, proton_constrainer
		TYPE (SymDensityProj) :: new_density, new_density_proj
		TYPE (SymHartreeFockFieldProj) :: field1, field2
		
		REAL(KIND = 16) :: TestAccuracy
		DOUBLE PRECISION :: b, diff, R2
		DOUBLE PRECISION :: accuracy, Gauge
		
		COMPLEX :: factor

		INTEGER :: A, N, Z, niter, ta
		INTEGER :: cycles_in, cycles_out, cycles_rate, Lold, Lmin
		INTEGER :: IndexGauge, NumberOfGauge
		
		! Leemos el reloj actual del sistema
		CALL SYSTEM_CLOCK(cycles_in)

		diff = 0.0
		
		b = Nucleus_get_b(diagonal%consistency%density%nucleus)
		
		A = Nucleus_get_A(diagonal%consistency%density%nucleus)
		
		WRITE(*,'(/,"--- Starting the self-consistent Calculation ---",/)')
	        WRITE(*,'("Number of oscillator shells N : ",I6)') N_0
		WRITE(*,'("Maximum angular momentum L    : ",I6)') Lmax
		WRITE(*,'("Maximum n                     : ",I6)') MIN(Nmax, NmaxOfL(0))
		WRITE(*,'("Numerical accuracy in Energy  : ",ES10.2)') tolerance
		WRITE(*,'("Atomic number A               : ",i6)') A
		WRITE(*,'("Oscillator length             : ",F6.3)') b
		
		IF ((A .LE. 1) .OR. (A .GE. 300)) THEN
			PRINT *, "Unexpected A value = ", A
			STOP "DiagonalizationMethod_goto_SelfConsistency"
		END IF
		
		factor = 1.0 - (1.0 / A)

		WRITE(*,'(/)')
		WRITE(*,'("Largest Number  : ",E24.16)') HUGE(TestAccuracy)
		WRITE(*,'("Smallest Number : ",E24.16)') TINY(TestAccuracy)
		WRITE(*,'(/)')

		! Calculation of the matrix elements of the Gogny force.
		! This is the most time- and memory-consuming task

		CALL SymVLSph_calculate(diagonal%consistency%vLSph)
		CALL SymVLSpp_calculate(diagonal%consistency%vLSpp)
		
		! This array takes a lot of space in memory: we don't need it from now on
		
		IF (Basis .EQ. 2) DEALLOCATE(WaveDeri)

        	! Creating and allocating memory for arrays helping accelerate
		! the calculation of radial integrals.
		! Defined in module "ibb.f90"
		
		IF (Optimization .EQ. 1) THEN
		
			Lmin = Lmax

			IF (Basis .EQ. 2) CALL StoreInt_new(Lmin)
		
			CALL SymVBBph_calculate(diagonal%consistency%vBBph, Lmin)
			CALL SymVBBpp_calculate(diagonal%consistency%vBBpp, Lmin)
		
			IF (Basis .EQ. 2) CALL StoreInt_del(Lmin)
			
			CALL SymVCph_calculate(diagonal%consistency%vCph, Lmin)
			CALL SymVCpp_calculate(diagonal%consistency%vCpp, Lmin)
		
		ELSE
		
			Lmin = 0
			
			IF (Basis .EQ. 2) CALL StoreInt_new(Lmin)
		
			CALL SymVBBph_calculate(diagonal%consistency%vBBph, Lmin)
			CALL SymVBBpp_calculate(diagonal%consistency%vBBpp, Lmin)
		
			IF (Basis .EQ. 2) CALL StoreInt_del(Lmin)
		
			CALL SymVCph_calculate(diagonal%consistency%vCph, Lmin)
			CALL SymVCpp_calculate(diagonal%consistency%vCpp, Lmin)
		
		END IF
		
		! Creation of the tensor fields

		CALL SymHartreeFockFieldProj_new(field1)
		CALL SymHartreeFockFieldProj_new(field2)

		niter = 0
gsc:		DO
  			N = Nucleus_get_N(diagonal%consistency%density%nucleus)
 			Z = Nucleus_get_Z(diagonal%consistency%density%nucleus)
			
			A = N + Z
			
			IF (A .LE. 1) THEN
				PRINT *, "Unexpected A value = ", A
				STOP "DiagonalizationMethod_goto_SelfConsistency"
			END IF
			
			Gauge = 4.0*ATAN(1.0)
			
			CALL SymGenDensityProj_setGauge(diagonal%iterated, Gauge, 0)
			CALL SymGenDensityProj_setGauge(diagonal%iterated, Gauge, 1)

			CALL DiagonalizationMethod_get_MeanField(diagonal, diagonal%S, niter)
									
			! We obtain the Fermi energies Lambda_n and Lambda_p by considering the constraint over the number of particle

			IF (HFOnly .EQ. 0) THEN

				CALL DiagonalizationMethod_set_ISOSpin(diagonal, 1)

				CALL OneDimSolve_new(neutron_constrainer, diagonal%func, diagonal)
			
				diagonal%consistency%density%nucleus%lambda_np(1) = OneDimSolve_solve(neutron_constrainer, &
					DBLE(N), diagonal%consistency%density%nucleus%lambda_np(1), &
					0.5 * (diagonal%consistency%density%nucleus%np(1) - diagonal%consistency%density%nucleus%actual_np(1)) &
					+ 0.1)

				CALL DiagonalizationMethod_set_ISOSpin(diagonal, 0)

				CALL OneDimSolve_new(proton_constrainer, diagonal%func, diagonal)
			
				diagonal%consistency%density%nucleus%lambda_np(0) = OneDimSolve_solve(proton_constrainer, &
					DBLE(Z), diagonal%consistency%density%nucleus%lambda_np(0), &
					0.5 * (diagonal%consistency%density%nucleus%np(0) - diagonal%consistency%density%nucleus%actual_np(0)) & 
					+ 0.1)
			
			ELSE
				
				diagonal%consistency%density%nucleus%actual_np(1) = diagonal%consistency%density%nucleus%np(0)
				diagonal%consistency%density%nucleus%actual_np(0) = diagonal%consistency%density%nucleus%np(1)				
				
				write(*,'("Lambda (n) = ",F10.5," Lambda (p) = ",F10.5)') &
				diagonal%consistency%density%nucleus%lambda_np(1),diagonal%consistency%density%nucleus%lambda_np(0)
			
				write(*,'("N = ",F10.5," Z = ",F10.5)') &
				diagonal%consistency%density%nucleus%actual_np(1),diagonal%consistency%density%nucleus%actual_np(0)
			
			END IF

			! Store the energies (for future use) of this step after correction for the particle number (above)
			
			CALL SelfConsistencyMethodProj_store_eHFB(diagonal%consistency)

			!-------------------------------------------------------------------------------------------------------!
			!													!
			!			   	ULTRA-IMPORTANT CALL							!
			!				--------------------							!
			!													!
			!    diagonal%iterated - TYPE: SymGenDensityProj   MODULE: symgden_proj.f90				!
			!    new_density       - TYPE: SymDensityProj      MODULE: symden_proj.f90				!
			!													!
			!  The purpose of the following CALL SymDensity_new_GenDensityProj() is actually to COPY the content of !
			!  the density matrix contained in diagonal%iterated, and obtained after the diagonalization of the HFB	! 
			!  matrix from U and V vectors, into new_density, which will be used to calculate the field Gamma and 	!
			!  Delta.												!
			!			  										!
			!-------------------------------------------------------------------------------------------------------!
			  
			CALL SymDensity_new_GenDensityProj(new_density, diagonal%iterated)
			
                        ! The current step corresponds to the density diagonal%consistency%density, the new one to new_density
			! We compute here the "distance" between the two.
			  
			IF (MOD(niter, MAX_ITER) .EQ. 0) THEN
				diff = SymHartreeFockFieldProj_distance(diagonal%consistency%density%field%rho, new_density%field%rho)
				WRITE(*,'("          k = ",I4," Difference in Energy : ",F15.8)') niter,diff
			END IF
			
			niter = niter + 1

                        ! We reduce the step made into the direction of the new density, new_density, from the old one, 
			! diagonal%consistency%density, with some slowing-factor anneal
			!
			!    rho_{n+1} -> (alpha * rho_{n} + rho_{n+1})/(1 + alpha)
			!
			  
			factor = 1.0 / (1.0 + anneal)
			
			! new step for normal density
			
			CALL SymHartreeFockFieldProj_product(field1, anneal*CMPLX(1,0), diagonal%consistency%density%field%rho)
			CALL SymHartreeFockFieldProj_add(field2, new_density%field%rho, field1)
			CALL SymHartreeFockFieldProj_product(diagonal%consistency%density%field%rho, factor, field2)

			CALL SymHartreeFockFieldProj_product(field1, anneal*CMPLX(1,0), diagonal%consistency%density%field%kap)
			CALL SymHartreeFockFieldProj_add(field2, new_density%field%kap, field1)
			CALL SymHartreeFockFieldProj_product(diagonal%consistency%density%field%kap, factor, field2)

			! We update the numbers of particles with the same rule

			DO ta = 0, 1
				diagonal%consistency%density%nucleus%actual_np(ta) = factor * ((new_density%nucleus%actual_np(ta)) &
					+ (anneal * diagonal%consistency%density%nucleus%actual_np(ta)))
                        
			END DO

			accuracy = SelfConsistencyMethodProj_accuracy(diagonal%consistency)
			WRITE(*,'("Iteration k = ",I4," HF Energy : ",F15.8)') niter,diagonal%consistency%density%nucleus%eHFB
			
			IF ((diff .LE. tolerance) .AND. (MOD(niter, MAX_ITER) .EQ. 0)) THEN
				accuracy = SelfConsistencyMethodProj_accuracy(diagonal%consistency)
                                WRITE(*,'("Energy = ",EN15.5," Precision = ",ES12.5)') diagonal%consistency%density%nucleus%eHFB,accuracy
				IF (accuracy .LE. tolerance) THEN
					EXIT gsc
				END IF
			END IF
			
			IF (niter .GE. NITER_MAX) THEN
				EXIT gsc
		        END IF
			
		END DO gsc

		CALL SymHartreeFockFieldProj_del(field1)
		CALL SymHartreeFockFieldProj_del(field2)

		CALL SymDensityProj_save(diagonal%consistency%density)

		CALL SYSTEM_CLOCK(cycles_out, cycles_rate)
		PRINT "(/A,EN10.2)", "Tiempo empleado (segundos): ", (DBLE(cycles_out - cycles_in) / cycles_rate)

		CALL SymDensityProj_store_actual_R2(diagonal%consistency%density)

		CALL SelfConsistencyMethodProj_store_eHFB(diagonal%consistency)
		CALL SelfConsistencyMethodProj_show_Status(diagonal%consistency)
		
		CALL DiagonalizationMethod_ProjectedEnergy(diagonal)
	
		CALL DiagonalizationMethod_show_ParticleEnergies(diagonal, Gauge)
                CALL DiagonalizationMethod_show_QuasiParticleEnergies(diagonal)
		
		RETURN
	END SUBROUTINE DiagonalizationMethod_goto_SelfConsistency



	SUBROUTINE DiagonalizationMethod_get_MeanField(diagonal, S, Niter)		
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		TYPE (SymGenDensityProj), INTENT(INOUT) :: S

		INTEGER, INTENT(IN) :: Niter

		TYPE (SymHartreeFockFieldProj), DIMENSION(:, :), POINTER :: HF_Gamma_VAP
		TYPE (SymHartreeFockFieldProj), DIMENSION(:), POINTER :: HF_Delta01_VAP, HF_Delta10_VAP
		TYPE (SymHartreeFockFieldProj), DIMENSION(:), POINTER :: Mfactor
		TYPE (SymHartreeFockFieldProj) :: HF_Gamma, HF_Delta, HF_Lambda
		TYPE (SymHartreeFockFieldProj) :: ekcm_field, vbb_field, vc_field, vls_field, gdd_field
		TYPE (SymHartreeFockFieldProj) :: field1, field2
		TYPE (SymD3Tensor) :: ek_tensor
		TYPE (SymGenDensityHFProj) :: gendenhf_gamma, gendenhf_delta

		TYPE (SymGenDensityGaugeProj), DIMENSION(:), POINTER :: genden_gauge
		TYPE (SymGenDensityHFProj), DIMENSION(:), POINTER :: C_Matrix, Deriv_Y
		TYPE (ProjectionCoeffs) :: CoeffsXY
		TYPE (SymGenDensityProj) :: genden_proj, genden_C
		TYPE (SymDensityProj) :: density, density01, density10		
		TYPE (SymGDDph) :: gDDph_phi, gDDph_proj

		COMPLEX :: factor, Trace_VAP
		DOUBLE PRECISION :: b, AngleGauge, Trace_NoProj,pi
		INTEGER :: isoCoupling, ta, tb, A, IndexGauge

!	integer la,u1,u2

		b = Nucleus_get_b(diagonal%consistency%density%nucleus)
		A = Nucleus_get_A(diagonal%consistency%density%nucleus)
		
		factor = (1.0 - (1.0 / A))*CMPLX(1, 0)

		IF ((A .LE. 1) .OR. (A .GE. 300)) THEN
			PRINT *, "Unexpected A value = ", A
			STOP "DiagonalizationMethod::get_MeanField"
		END IF
		
		!***********************************************************************!
		!***********************************************************************!
		!***********************************************************************!
		
		!-----------------------------------------------------------------------!
		!  PROJECTION BLOCK:							!
		!	- CALCULATE THE FOMENKO COEFFICIENTS (DEPEND ON DENSITY MATRIX)	!
		!	- CALCULATE GAUGE-DEPENDENT DENSITIES				!
		!	- CALCULATE PROJECTED DENSITY					!
		!	- CALCULATE C-MATRIX AND DERIVATIVES OF FOMENKO COEFFICIENTS	!
		!-----------------------------------------------------------------------!

		IF (Niter .GE. 1) THEN
		
			! Create new Fomenko coefficients
		
			CALL ProjectionCoeffs_new(CoeffsXY, NGauge)
		
			! Calculate phi-dependent densities rho(phi), kappa10(phi), kappa01(phi) according
			! to formulas (9), (10) and (11) of the Anguiano et al, Nucl. Phys. A696 (2001) 467-493
			! All densities are stored into genden_gauge.
			CALL SymGenDensityProj_make_DensityGauge(diagonal%UV, genden_gauge, NGauge)
		
			! Calculate the projected density and the Fomenko coefficients
			CALL SymGenDensityProj_make_DensityProj(genden_proj, CoeffsXY, genden_gauge, &
					diagonal%UV, diagonal%consistency%density%nucleus, NGauge)
		
			! Convert projected density into the right type, which is SymDensityProj
			CALL SymDensity_new_GenDensityProj(density, genden_proj)
		
			! Calculate the c-matrix and the derivative of the y-coefficients with respect to the density matrix
			CALL ProjectionPreparation(diagonal%UV, CoeffsXY, genden_gauge, C_Matrix, Deriv_Y, NGauge)
		
			!-----------------------------------------------------------------------!
			!     PROJECTION BLOCK: CALCULATION OF THE PROJECTED KINETIC ENERGY	!							!
			!-----------------------------------------------------------------------!
			
			pi = 4.0*ATAN(1.0)
			
			CALL SymHartreeFockFieldProj_new(field1)
			
			CALL SymEkTensorVAP_get_GammaProj(field1, EkField, genden_gauge, C_Matrix, Deriv_Y, CoeffsXY, NGauge)
			
			write(*,'(" Ek (PROJECTION) = ",2f20.5)') factor*(field1 * diagonal%consistency%density%field%rho)
		
			CALL SymHartreeFockFieldProj_del(field1)
			
			!-----------------------------------------------------------------------!
			!     PROJECTION BLOCK: CALCULATION OF THE PROJECTED MEAN-FIELD GAMMA	!							!
			!-----------------------------------------------------------------------!
			
			! Calculate the gauge-dependent mean-field Gamma_{ta, tb}(phi)
	
			ALLOCATE(HF_Gamma_VAP(0:1, 1:NGauge))
			ALLOCATE(HF_Delta10_VAP(1:NGauge))
			ALLOCATE(HF_Delta01_VAP(1:NGauge))
		
			Trace_VAP = CMPLX(0,0)
		
			! HF_Gamma_VAP(0, i) contains the fields Gamma^{tau, tau} for each gauge angle i
			! HF_Gamma_VAP(1, i) contains the fields Gamma^{tau, tau'} for each gauge angle i
			
			DO isoCoupling = 0, 1
			
				DO IndexGauge = 1, NGauge
		
					AngleGauge = pi*REAL(IndexGauge)/REAL(NGauge)			
	
					CALL SymHartreeFockFieldProj_new(HF_Gamma_VAP(isoCoupling, IndexGauge))
				
					CALL SymGammaTensorVAP_get_GammaPhi(HF_Gamma_VAP(isoCoupling, IndexGauge), genden_gauge(IndexGauge), &
											genden_proj, isoCoupling, AngleGauge, diagonal%consistency)
											
			
				END DO
				
			END DO
			
			DO IndexGauge = 1, NGauge
		
				CALL SymHartreeFockFieldProj_new(HF_Delta01_VAP(IndexGauge))
				
				CALL SymGammaTensorVAP_get_DeltaPhi(HF_Delta01_VAP(IndexGauge), genden_gauge(IndexGauge), diagonal%consistency, 0)
			
				CALL SymHartreeFockFieldProj_new(HF_Delta10_VAP(IndexGauge))
					
				CALL SymGammaTensorVAP_get_DeltaPhi(HF_Delta10_VAP(IndexGauge), genden_gauge(IndexGauge), diagonal%consistency, 1)
											
			END DO
			
			!-----------------------------------------------------------------------!
			!    GAUGE-DEPENDENT M-INTEGRALS					!
			!-----------------------------------------------------------------------!
			
			!ALLOCATE(Mfactor(1:NGauge))
			
			!CALL SymGDDph_new(gDDph_phi)
			!CALL SymGDDph_new(gDDph_proj)
			
			DO IndexGauge = 1, NGauge
			
				!CALL SymDensity_new_GenDensityProj10(density10, genden_gauge(IndexGauge))
				
				AngleGauge = pi*REAL(IndexGauge)/REAL(NGauge)
				
				!CALL SymHartreeFockFieldProj_new(Mfactor(IndexGauge))
				
				!CALL SymGDDph_GammaVAP(Mfactor(IndexGauge), gDDph_phi, gDDph_proj, density10%field%rho, density%field%rho, &
				!												AngleGauge, 1)
																
			END DO
			
			DO ta = 0, 1
			
				DO IndexGauge = 1, NGauge
			
					AngleGauge = pi*REAL(IndexGauge)/REAL(NGauge)
				
					!CALL SymHartreeFockFieldProj_new(field1)
					
					!CALL SymGammaTensorVAP_get_GammaDDPhi(field1, Mfactor, CoeffsXY, genden_gauge, C_Matrix, Deriv_Y, &
					!										NGauge, IndexGauge, ta)
											
					!CALL SymHartreeFockFieldProjIso_add(HF_Gamma_VAP(ta, ta, IndexGauge), HF_Gamma_VAP(ta, ta, IndexGauge), &
					!												field1, 0)
																					
					!CALL SymHartreeFockFieldProj_del(field1)

				END DO
				
			END DO

			!-----------------------------------------------------------------------!
			!    CONSTRUCTING THE PROJECTED HFB MATRIX				!
			!-----------------------------------------------------------------------!
			
			CALL SymHartreeFockFieldProj_new(HF_Gamma)
			CALL SymHartreeFockFieldProj_new(HF_Lambda)
			CALL SymHartreeFockFieldProj_new(HF_Delta)
			CALL SymHartreeFockFieldProj_new(field1)
			
			CALL SymGammaTensorVAP_get_GammaP(HF_Gamma, CoeffsXY, genden_gauge, C_Matrix, Deriv_Y, &
							NGauge, HF_Gamma_VAP, diagonal%consistency%density%nucleus)
			
			CALL SymGammaTensorVAP_get_LambdaP(HF_Lambda, CoeffsXY, genden_gauge, C_Matrix, Deriv_Y, &
							NGauge, HF_Delta01_VAP, HF_Delta10_VAP, diagonal%consistency%density%nucleus)
							
			CALL SymHartreeFockFieldProj_add(field1, HF_Gamma, HF_Lambda)
			
			CALL SymGammaTensorVAP_get_DeltaP(HF_Delta, CoeffsXY, C_Matrix, NGauge, HF_Delta01_VAP)
								
			write(*,'("vBB (PROJECTION) = ",2f20.5)') 0.5*(field1 * diagonal%consistency%density%field%rho)
			write(*,'("Hpp (PROJECTION) = ",2f20.5)') 0.5*(diagonal%consistency%density%field%kap * HF_Delta) 
			
			CALL SymHartreeFockFieldProj_del(field1)
			CALL SymHartreeFockFieldProj_del(HF_Gamma)
			CALL SymHartreeFockFieldProj_del(HF_Delta)
			CALL SymHartreeFockFieldProj_del(HF_Lambda)
			
			!-----------------------------------------------------------------------!
			!  PROJECTION BLOCK: DEALLOCATING ALL ARRAYS FOR PROJECTION IN ORDER TO	!
			!  BE READY FOR THE NEXT ITERATION					!
			!-----------------------------------------------------------------------!
		
			DO IndexGauge = 1, NGauge		
				
				CALL SymGenDensityHFProj_del(genden_gauge(IndexGauge)%rho)
				!CALL SymGenDensityHFProj_del(C_Matrix(IndexGauge)%rho)
				!CALL SymGenDensityHFProj_del(Deriv_Y(IndexGauge)%rho)
				
				CALL SymGenDensityHFProj_del(C_Matrix(IndexGauge))
				CALL SymGenDensityHFProj_del(Deriv_Y(IndexGauge))
				
				CALL SymHartreeFockFieldProj_del(HF_Delta01_VAP(IndexGauge))
				CALL SymHartreeFockFieldProj_del(HF_Delta10_VAP(IndexGauge))
				
				DO ta = 0, 1
					CALL SymHartreeFockFieldProj_del(HF_Gamma_VAP(ta, IndexGauge))
				END DO
				
			END DO
			
			CALL SymGenDensityProj_del(genden_proj)
			CALL ProjectionCoeffs_del(CoeffsXY)
			CALL SymDensityProj_del(density)
			
			DEALLOCATE(HF_Gamma_VAP)
			DEALLOCATE(HF_Delta01_VAP)
			DEALLOCATE(HF_Delta10_VAP)
			
			DEALLOCATE(C_Matrix)
			DEALLOCATE(Deriv_Y)
		
		END IF
	
		!***********************************************************************!
		!***********************************************************************!
		!***********************************************************************!
		

		CALL SymHartreeFockFieldProj_new(HF_Gamma)
		CALL SymHartreeFockFieldProj_new(HF_Delta)
		
		CALL SymD3Tensor_new(ek_tensor)
		
		CALL SymHartreeFockFieldProj_new(ekcm_field)
		CALL SymHartreeFockFieldProj_new(vbb_field)
		CALL SymHartreeFockFieldProj_new(vc_field)
		CALL SymHartreeFockFieldProj_new(vls_field)
		CALL SymHartreeFockFieldProj_new(gdd_field)
		
		CALL SymHartreeFockFieldProj_new(field1)
		CALL SymHartreeFockFieldProj_new(field2)

		! Mean-field - Kinetic energy

               	CALL SymD3Tensor_product(ek_tensor, factor, EkField)
		write(*,'("Ek (NO PROJECTION) = ",2f20.5)') (ek_tensor * diagonal%consistency%density%field%rho%p(0) + &
							     ek_tensor * diagonal%consistency%density%field%rho%p(1))
	
		! Mean-field - Two-body center of mass correction

		! Protons
		CALL SymKineticEnergy2Body_product(field1, diagonal%consistency%vEkCMph, diagonal%consistency%density%field%rho, 0, 0, 0)
		! Neutrons
		CALL SymKineticEnergy2Body_product(field1, diagonal%consistency%vEkCMph, diagonal%consistency%density%field%rho, 1, 1, 0)
				
		CALL SymHartreeFockFieldProj_product(ekcm_field, CMPLX(1.0/A, 0), field1)
		
		! Mean-field - Brink-Boker term. We contract a matrix element v_{abcd}(ta) o given isospin ta with a projected density
		! of different isospin tb

		! Protons
		CALL SymVBBph_get_Gamma(vbb_field, diagonal%consistency%vBBph, diagonal%consistency%density%field%rho, 0, 0)
		! Neutrons
		CALL SymVBBph_get_Gamma(vbb_field, diagonal%consistency%vBBph, diagonal%consistency%density%field%rho, 1, 1)

		! Mean-field - Coulomb potential

		! Protons
		CALL SymVCph_get_Gamma(vc_field, diagonal%consistency%vCph, diagonal%consistency%density%field%rho, 0, 0, 0)
		! Neutrons (vanishes of course)
		CALL SymVCph_get_Gamma(vc_field, diagonal%consistency%vCph, diagonal%consistency%density%field%rho, 1, 1, 0)

		! Mean-field - Spin-orbit term

		! Protons
		CALL SymVLSph_get_Gamma(vls_field, diagonal%consistency%vLSph, diagonal%consistency%density%field%rho, 0, 0, 0)
		! Neutrons
		CALL SymVLSph_get_Gamma(vls_field, diagonal%consistency%vLSph, diagonal%consistency%density%field%rho, 1, 1, 0)

		! Mean-field - Density-dependent term. This term needs to be calculated with the projected density, contrary
		!              all the others. Therefore, we need to copy this projected density into a suitable array

		CALL SymGDDph_update(gdd_field, diagonal%consistency%gDDph, diagonal%consistency%density%field%rho, AngleGauge)

		! Total Mean-field = Sum of all the preceding terms

		CALL SymHartreeFockFieldProj_add(field2, vls_field, vbb_field)
		CALL SymHartreeFockFieldProj_add(field1, vc_field, field2)
		CALL SymHartreeFockFieldProj_add(field2, ekcm_field, field1)
		CALL SymHartreeFockFieldProj_add(field1, gdd_field, field2)
		CALL SymHartreeFockFieldProj_add_SymD3Tensor(HF_Gamma, ek_tensor, field1)
			write(*,'("vBB (NO PROJECTION) = ",2f20.5)')  0.5*(vbb_field * diagonal%consistency%density%field%rho)

		! Pairing - Brink-Boker term

		IF (HFOnly .EQ. 0) THEN
			CALL SymVBBpp_get_Delta(vbb_field, diagonal%consistency%vBBpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymVBBpp_get_Delta(vbb_field, diagonal%consistency%vBBpp, diagonal%consistency%density%field%kap, 1, 1)
		ELSE
			CALL SymVBBpp_get_Delta(field1, diagonal%consistency%vBBpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymHartreeFockFieldProj_product(vbb_field, CMPLX(0, 0), field1)
		END IF

		! Pairing - Coulomb term

		IF (HFOnly .EQ. 0) THEN
			CALL SymVCpp_get_Delta(vc_field, diagonal%consistency%vCpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymVCpp_get_Delta(vc_field, diagonal%consistency%vCpp, diagonal%consistency%density%field%kap, 1, 1)
		ELSE
			CALL SymVCpp_get_Delta(field1, diagonal%consistency%vCpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymHartreeFockFieldProj_product(vc_field, CMPLX(0, 0), field1)
		END IF

		! Pairing - Spin-orbit term
	
		IF (HFOnly .EQ. 0) THEN
			CALL SymVLSpp_get_Delta(vls_field, diagonal%consistency%vLSpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymVLSpp_get_Delta(vls_field, diagonal%consistency%vLSpp, diagonal%consistency%density%field%kap, 1, 1)
		ELSE
			CALL SymVLSpp_get_Delta(field1, diagonal%consistency%vLSpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymHartreeFockFieldProj_product(vls_field, CMPLX(0, 0), field1)
		END IF

		! Total Pairing = Sum of all the preceding terms

		CALL SymHartreeFockFieldProj_add(field2, vc_field, vls_field)
		CALL SymHartreeFockFieldProj_add(HF_Delta, vbb_field, field2)
			write(*,'("Hpp (NO PROJECTION) = ",2f20.5)')  0.5*(HF_Delta * diagonal%consistency%density%field%kap)

		CALL SymGenDensityHFProj_new(gendenhf_gamma)
		CALL SymGenDensityHFProj_new(gendenhf_delta)



		CALL SymGenDensityHFProj_copy(gendenhf_gamma, HF_Gamma)
		CALL SymGenDensityHFProj_copy(gendenhf_delta, HF_Delta)

		CALL SymGenDensityProj_new_GammaDelta(S, gendenhf_gamma, gendenhf_delta, b)



		CALL SymGenDensityHFProj_del(gendenhf_gamma)
		CALL SymGenDensityHFProj_del(gendenhf_delta)

		CALL SymD3Tensor_del(ek_tensor)
		
		CALL SymHartreeFockFieldProj_del(HF_Gamma)
		CALL SymHartreeFockFieldProj_del(HF_Delta)

		CALL SymHartreeFockFieldProj_del(ekcm_field)
		CALL SymHartreeFockFieldProj_del(vbb_field)
		CALL SymHartreeFockFieldProj_del(vc_field)
		CALL SymHartreeFockFieldProj_del(vls_field)
		CALL SymHartreeFockFieldProj_del(gdd_field)

		CALL SymHartreeFockFieldProj_del(field1)
		CALL SymHartreeFockFieldProj_del(field2)
				
		RETURN
	END SUBROUTINE DiagonalizationMethod_get_MeanField

	SUBROUTINE DiagonalizationMethod_set_ISOSpin(diagonal, ta)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		INTEGER, INTENT(IN) :: ta

		diagonal%ta = ta
		
		RETURN
	END SUBROUTINE DiagonalizationMethod_set_ISOSpin

	FUNCTION DiagonalizationMethod_get_WaveFunction(diagonal)
		TYPE (WaveFunction) DiagonalizationMethod_get_WaveFunction
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal

		INTEGER :: N, Z, A, u1, u2, ta, a2, d, dd, la, ja
		
		COMPLEX :: factor
		DOUBLE PRECISION :: b, trace, Gauge
		DOUBLE PRECISION, DIMENSION(0:1) :: np, R2

		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: tmp, tmp1, tmp2
		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: U, V, UV, SH, Delta, h, I
		COMPLEX, DIMENSION(:), ALLOCATABLE :: Diag2

		TYPE (SymHartreeFockFieldProj) :: HF_Gamma, HF_Delta
		TYPE (SymHartreeFockFieldProj) :: ekcm_field, vbb_field, vc_field, vls_field, gdd_field
		TYPE (SymHartreeFockFieldProj) :: field1, field2
		TYPE (SymGenDensityProj) :: SuperHamiltonian
		TYPE (SymGenDensityHFProj) :: gendenhf_gamma, gendenhf_delta
		TYPE (SymD3Tensor) :: ek_tensor

                np(0) = 0.0
                np(1) = 0.0

                R2(0) = 0.0
                R2(1) = 0.0
                
		N = Nucleus_get_N(diagonal%consistency%density%nucleus)
		Z = Nucleus_get_Z(diagonal%consistency%density%nucleus)
		
		Gauge = diagonal%consistency%density%field%rho%GaugeAngle(0)
		
		A = N + Z
		
		factor = 1.0 - (1.0 / A)
		
		b = Nucleus_get_b(diagonal%consistency%density%nucleus)
		
		! ATENCION: El siguiente procedimiento crea un nucleo a partir
		! del numero de neutrones y protones (N y Z). Nucleo que es sobreescrito
		! a continuacion... Esto ha sido mantenido por "similitud" con el
		! programa original, pero deberia ser corregido en una version posterior
		
		CALL WaveFunction_new(DiagonalizationMethod_get_WaveFunction, N, Z)
		CALL Nucleus_copy(DiagonalizationMethod_get_WaveFunction%nucleus, diagonal%consistency%density%nucleus)

		CALL SymD3Tensor_new(ek_tensor)

		CALL SymHartreeFockFieldProj_new(HF_Gamma)
		CALL SymHartreeFockFieldProj_new(HF_Delta)
		
		CALL SymHartreeFockFieldProj_new(ekcm_field)
		CALL SymHartreeFockFieldProj_new(vbb_field)
		CALL SymHartreeFockFieldProj_new(vc_field)
		CALL SymHartreeFockFieldProj_new(vls_field)
		CALL SymHartreeFockFieldProj_new(gdd_field)
		
		CALL SymHartreeFockFieldProj_new(field1)
		CALL SymHartreeFockFieldProj_new(field2)
		
		! Mean-field - Kinetic energy

                CALL SymD3Tensor_product(ek_tensor, factor, EkField)

		! Mean-field - Center of mass correction

		! Protons
		CALL SymKineticEnergy2Body_product(field1, diagonal%consistency%vEkCMph, diagonal%consistency%density%field%rho, 0, 0, 0)
		! Neutrons
		CALL SymKineticEnergy2Body_product(field1, diagonal%consistency%vEkCMph, diagonal%consistency%density%field%rho, 1, 1, 0)

		CALL SymHartreeFockFieldProj_product(ekcm_field, CMPLX(1.0/A, 0), field1)

		! Mean-field - Brink-Boker term

		! Protons
		CALL SymVBBph_get_Gamma(vbb_field, diagonal%consistency%vBBph, diagonal%consistency%density%field%rho, 0, 0)
		! Neutrons
		CALL SymVBBph_get_Gamma(vbb_field, diagonal%consistency%vBBph, diagonal%consistency%density%field%rho, 1, 1)

		! Mean-field - Coulomb term

		! Protons
		CALL SymVCph_get_Gamma(vc_field, diagonal%consistency%vCph, diagonal%consistency%density%field%rho, 0, 0, 0)
		! Neutrons (vanishes of course)
		CALL SymVCph_get_Gamma(vc_field, diagonal%consistency%vCph, diagonal%consistency%density%field%rho, 1, 1, 0)

		! Mean-field - Spin-orbit term

		! Protons
		CALL SymVLSph_get_Gamma(vls_field, diagonal%consistency%vLSph, diagonal%consistency%density%field%rho, 0, 0, 0)
		! Neutrons
		CALL SymVLSph_get_Gamma(vls_field, diagonal%consistency%vLSph, diagonal%consistency%density%field%rho, 1, 1, 0)

		! Mean-field - Density-dependent term. This term needs to be calculated with the projected density, contrary
		!              all the others. Therefore, we need to copy this projected density into a suitable array

		CALL SymGDDph_update(gdd_field, diagonal%consistency%gDDph, diagonal%consistency%density%field%rho, Gauge)

		! Total Mean-field = Sum of all the preceding terms

		CALL SymHartreeFockFieldProj_add(field1, vls_field, gdd_field)
		CALL SymHartreeFockFieldProj_add(field2, vc_field, field1)
		CALL SymHartreeFockFieldProj_add(field1, vbb_field, field2)
		CALL SymHartreeFockFieldProj_add(field2, ekcm_field, field1)
		CALL SymHartreeFockFieldProj_add_SymD3Tensor(HF_Gamma, ek_tensor, field2)

		! Pairing - Brink-Boker term

		IF (HFOnly .EQ. 0) THEN
			CALL SymVBBpp_get_Delta(vbb_field, diagonal%consistency%vBBpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymVBBpp_get_Delta(vbb_field, diagonal%consistency%vBBpp, diagonal%consistency%density%field%kap, 1, 1)
		ELSE
			CALL SymVBBpp_get_Delta(field1, diagonal%consistency%vBBpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymHartreeFockFieldProj_product(vbb_field, CMPLX(0, 0), field1)
		END IF

		! Pairing - Coulomb term

		IF (HFOnly .EQ. 0) THEN
			CALL SymVCpp_get_Delta(vc_field, diagonal%consistency%vCpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymVCpp_get_Delta(vc_field, diagonal%consistency%vCpp, diagonal%consistency%density%field%kap, 1, 1)
		ELSE
			CALL SymVCpp_get_Delta(field1, diagonal%consistency%vCpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymHartreeFockFieldProj_product(vc_field, CMPLX(0, 0), field1)
		END IF

		! Pairing - Spin-orbit term
	
		IF (HFOnly .EQ. 0) THEN
			CALL SymVLSpp_get_Delta(vls_field, diagonal%consistency%vLSpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymVLSpp_get_Delta(vls_field, diagonal%consistency%vLSpp, diagonal%consistency%density%field%kap, 1, 1)
		ELSE
			CALL SymVLSpp_get_Delta(field1, diagonal%consistency%vLSpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymHartreeFockFieldProj_product(vls_field, CMPLX(0, 0), field1)
		END IF

		! Total Pairing = Sum of all the preceding terms

		CALL SymHartreeFockFieldProj_add(field2, vc_field, vls_field)
		CALL SymHartreeFockFieldProj_add(HF_Delta, vbb_field, field2)

		CALL SymGenDensityHFProj_new(gendenhf_gamma)
		CALL SymGenDensityHFProj_new(gendenhf_delta)



		CALL SymGenDensityHFProj_copy(gendenhf_gamma, HF_Gamma)
		CALL SymGenDensityHFProj_copy(gendenhf_delta, HF_Delta)

		CALL SymGenDensityProj_new_GammaDelta(SuperHamiltonian, gendenhf_gamma, gendenhf_delta, b)



		CALL SymGenDensityHFProj_del(gendenhf_gamma)
		CALL SymGenDensityHFProj_del(gendenhf_delta)

		DO ta = 0, 1
			DO a2 = 0, 2*Lmax
			
				d = DIM(a2)
				dd = d + d
				la = L(a2)
				ja = J(a2)

				ALLOCATE(I(d, d))
				ALLOCATE(h(d, d))
				ALLOCATE(Delta(d, d))
				ALLOCATE(SH(dd, dd))
				ALLOCATE(Diag2(dd))
				ALLOCATE(UV(dd, dd))
				ALLOCATE(U(d, d))
				ALLOCATE(V(d, d))
				ALLOCATE(tmp(d, d))

				I = CMPLX(0,0)
				DO u1 = 1, d
					I(u1, u1) = CMPLX(1,0)
				END DO

				h = SuperHamiltonian%rho%rho(ta, a2)%store &
							- (diagonal%consistency%density%nucleus%lambda_np(ta) * I) 

				Delta = SuperHamiltonian%kap%rho(ta, a2)%store

				! Componemos la matriz SH de la siguiente forma:
				!    h     -Delta
				!   -Delta -h
				DO u1 = 1, d
					DO u2 = 1, d
						SH(u1    , u2    ) =   h(u1, u2)
						SH(u1    , u2 + d) = - Delta(u1, u2)
						SH(u1 + d, u2    ) = - Delta(u1, u2)
						SH(u1 + d, u2 + d) = - h(u1, u2)
					END DO
				END DO

				CALL MatDiagCmplx(SH, UV, Diag2, dd)

				! We extract the U anv V matrices. The HFB matrix yields 2 spectra: one of positive
				! q.p. energies, the other one of negative q.p. energies. The occupied levels are 
				! found among the levels with positive q.p. energies. Hence the second column:
				!  HFB = | V*  U |
				!	 | U*  V |
				DO u1 = 1, d
					DO u2 = 1, d
						U(u1, u2) = UV(u1    , u2 + d)
						V(u1, u2) = UV(u1 + d, u2 + d)
					END DO
				END DO

				DiagonalizationMethod_get_WaveFunction%genden%rho%rho(ta, a2)%value = U
				DiagonalizationMethod_get_WaveFunction%genden%kap%rho(ta, a2)%value = V

				tmp = MATMUL(V, TRANSPOSE(V))
				
				! Calculamos la traza de la matriz "tmp"
				trace = 0.0
				DO u1 = 1, d
					trace = trace + REAL(tmp(u1, u1))
				END DO
				np(ta) = np(ta) + (DBLE(ja + 1.0) * trace)

				tmp2 = SymD3Tensor_matrix(R2Field, la)
				
				tmp = MATMUL(tmp2, MATMUL(V, TRANSPOSE(V)))

				! Calculamos la traza de la matriz "tmp"
				trace = 0.0
				DO u1 = 1, d
					trace = trace + REAL(tmp(u1, u1))
				END DO
				R2(ta) = R2(ta) + (DBLE(ja + 1.0) * trace) !/ (b**2)

				DEALLOCATE(I)
				DEALLOCATE(h)
				DEALLOCATE(Delta)
				DEALLOCATE(SH)
				DEALLOCATE(Diag2)
				DEALLOCATE(UV)
				DEALLOCATE(U)
				DEALLOCATE(V)
				DEALLOCATE(tmp)
			END DO
		END DO
				
		CALL SymD3Tensor_del(ek_tensor)

		CALL SymHartreeFockFieldProj_del(HF_Gamma)
		CALL SymHartreeFockFieldProj_del(HF_Delta)

		CALL SymHartreeFockFieldProj_del(ekcm_field)
		CALL SymHartreeFockFieldProj_del(vbb_field)
		CALL SymHartreeFockFieldProj_del(vc_field)
		CALL SymHartreeFockFieldProj_del(vls_field)
		CALL SymHartreeFockFieldProj_del(gdd_field)

		CALL SymHartreeFockFieldProj_del(field1)
		CALL SymHartreeFockFieldProj_del(field2)
		
		RETURN
	END FUNCTION DiagonalizationMethod_get_WaveFunction

	FUNCTION DiagonalizationMethod_operator(diagonal)
		DOUBLE PRECISION DiagonalizationMethod_operator
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		
		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: copyR2, copyRho, tmp
		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: xI, h, Delta, SH, U, V
		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: N20
		
		COMPLEX :: trace
		DOUBLE PRECISION :: np, R2, b2

		INTEGER :: a, d, dd, la, ja, i, jj, k
		INTEGER :: u1, u2
		
		np = 0.0
		R2 = 0.0
		b2 = Nucleus_get_b(diagonal%iterated%nucleus) ** 2
		
		DO a = 0, 2*Lmax
		
		
			d = DIM(a)
			dd = d + d
			
			la = L(a)
			ja = J(a)
			
			ALLOCATE(xI(d, d))
			ALLOCATE(h(d, d))
			ALLOCATE(Delta(d, d))
			ALLOCATE(SH(dd, dd))
			ALLOCATE(U(d, d))
			ALLOCATE(V(d, d))
			ALLOCATE(tmp(d, d))
			ALLOCATE(copyR2(d, d))
			ALLOCATE(copyRho(d, d))
			
			xI = 0.0
			DO u1 = 1, d
				xI(u1, u1) = diagonal%func%x
			END DO

			h = diagonal%S%rho%rho(diagonal%ta, a)%store - xI

			Delta = 1.0 * diagonal%S%kap%rho(diagonal%ta, a)%store

			! We construct the HFB matrix
			!    h        Delta
			!   -Delta*    -h*
			DO u1 = 1, d
				DO u2 = 1, d
					SH(u1    , u2    ) =   h(u1, u2)
					SH(u1    , u2 + d) = - Delta(u1, u2)
					SH(u1 + d, u2    ) = - Delta(u1, u2)
					SH(u1 + d, u2 + d) = - h(u1, u2)
				END DO
			END DO

			! ATENCION: Este procedimiento no ha sido traducido, sino
			! extraido a traves de Google.
			! Confiamos en su correcta implementacion

			! We obtain the quasi-particle energies in increasing order (from negative to positive values)

			CALL MatDiagCmplx(SH, diagonal%UV(diagonal%ta, a)%quantum, diagonal%QuasiParticleEnergies(diagonal%ta, a)%value, dd)

			! We extract the matrices U and V
			DO u1 = 1, d
				DO u2 = 1, d
					U(u1, u2) = diagonal%UV(diagonal%ta, a)%quantum(u1    , u2 + d)
					V(u1, u2) = diagonal%UV(diagonal%ta, a)%quantum(u1 + d, u2 + d)
				END DO
			END DO

			! We calculate the new densities rho(phi) and kappa(phi) from the U and V vectors
			! genden_proj contains the projected density rho^{Proj}
			
			CALL SymGenDensityProj_make_Block(diagonal%iterated, diagonal%ta, a, U, V)
			
			! We update the Gauge angle
			CALL SymGenDensityProj_setGauge(diagonal%iterated, diagonal%consistency%density%field%rho%GaugeAngle(0), 0)
			CALL SymGenDensityProj_setGauge(diagonal%iterated, diagonal%consistency%density%field%rho%GaugeAngle(1), 1)
			
			! The particle number is obtained by the trace of the density matrix
			trace = 0.0
			DO u1 = 1, d
				trace = trace + diagonal%iterated%rho%rho(diagonal%ta, a)%store(u1, u1)
			END DO
			np = np + DBLE(ja + 1.0) * trace

			! We need to copy the matrix R2Field because of dimensions inconsistencies when in the
			! case of an arbitrary basis. This inconsistency causes the multiplication below to
			! crash when using -C option at compilation stage
			copyR2 = SymD3Tensor_matrix(R2Field, la)
			copyRho = diagonal%iterated%rho%rho(diagonal%ta, a)%store

			tmp = MATMUL(copyR2, copyRho)
			
			trace = 0.0
			DO u1 = 1, d
				trace = trace + REAL(tmp(u1, u1))
			END DO
			R2 = R2 + DBLE(ja + 1.0) * trace !/ b2

			! Liberamos la memoria de las matrices utilizadas
			DEALLOCATE(xI)
			DEALLOCATE(h)
			DEALLOCATE(Delta)
			DEALLOCATE(SH)
			DEALLOCATE(U)
			DEALLOCATE(V)
			DEALLOCATE(tmp)
			DEALLOCATE(copyR2)
			DEALLOCATE(copyRho)
		END DO
		
		diagonal%iterated%nucleus%actual_np(diagonal%ta) = np
		diagonal%iterated%nucleus%actual_R2(diagonal%ta) = R2 / np
		
		! Resultado final
		DiagonalizationMethod_operator = np
		
		RETURN
	END FUNCTION DiagonalizationMethod_operator

	!-------------------------------------------------------------------------------!
	!										!
	!       Diagonalization of the matrix of the normal density rho. 		!
	!										!
	!  The single-particle energies are defined, in the case of HFB, as the		!
	!  expectation values of the hamiltonian in the basis where rho is diagonal.	!
	!  This basis is referred to as the canonical basis.				!
	!										!
	!-------------------------------------------------------------------------------!
	
	SUBROUTINE DiagonalizationMethod_show_ParticleEnergies(diagonal, Gauge)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		DOUBLE PRECISION, INTENT(IN) :: Gauge

		TYPE (SymGenDensityProj) G, SuperHamiltonian
		
		TYPE (DiagonalMatrix), DIMENSION(0:1) :: D, Occup, SingleParticleE, SingleParticleP, SingleParticleR2
		TYPE (MatrixType), DIMENSION(0:1) :: P, E, R2

		INTEGER, DIMENSION(:), ALLOCATABLE :: an, ap
		INTEGER, DIMENSION(:), ALLOCATABLE :: in, ip
		
		INTEGER, DIMENSION(:,:), ALLOCATABLE :: Count
		
		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: R2Matrix, R2Initial
		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: S, V
		DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: RadialWaveFunction
		
		INTEGER :: ndim, a, num, i, k, dim_cur, dim_acc, la, ja, ta, nrot, file_error
		INTEGER :: Lvalue, Jvalue
		INTEGER :: ipoint, loop1, loop2, IndexBra, IndexKet
		INTEGER :: m, n, Dummy
		
		INTEGER, PARAMETER :: file_unit_1 = 16, file_unit_2 = 17
		
		DOUBLE PRECISION :: b, spn, spp
		
		CHARACTER(8) :: label
		CHARACTER(80) :: file_neut, file_prot
		CHARACTER, DIMENSION(0:29) :: spectr
		
		DATA spectr / "s", "p", "d", "f", "g", "h", "i", "j", "k", "l", &
			      "m", "n", "o", "P", "q", "r", "S", "t", "u", "v", &
			      "w", "x", "y", "z", "?", "?", "?", "?", "?", "?" /

		b = Nucleus_get_b(diagonal%consistency%density%nucleus)

		dim_cur = 0
		DO a = 0, 2*Lmax
			dim_cur = dim_cur + DIM(a)
		END DO
		
		ALLOCATE(Occup(0)%value(dim_cur))
		ALLOCATE(Occup(1)%value(dim_cur))
		
		ALLOCATE(SingleParticleE(0)%value(dim_cur))
		ALLOCATE(SingleParticleE(1)%value(dim_cur))
		ALLOCATE(SingleParticleP(0)%value(dim_cur))
		ALLOCATE(SingleParticleP(1)%value(dim_cur))
		ALLOCATE(SingleParticleR2(0)%value(dim_cur))
		ALLOCATE(SingleParticleR2(1)%value(dim_cur))
		
		ALLOCATE(RadialWaveFunction(0:1, 0:Npoint, 1:dim_cur))

		IF (.NOT. ALLOCATED(DensityRadial)) ALLOCATE(DensityRadial(0:1, 0:Npoint))

		! "SuperHamiltonian" contains the generalized hamiltonian 
		!
		!                 |     h      -Delta  |
		!                 |                    |
		!                 |  -Delta     -h     |
		!
		!  SuperHamiltonian%rho = h (type of SymHartreeFockField)
		!  SuperHamiltonian%kap = Delta (type of SymHartreeFockField)

		Dummy = 2

		CALL DiagonalizationMethod_get_MeanField(diagonal, SuperHamiltonian, Dummy)

							ndim = Nmax * (2*Lmax + 1)
		IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	ndim = (N_0 + 1) * (N_0 + 2) / 2
		
		ALLOCATE(an(ndim + 1))
		ALLOCATE(ap(ndim + 1))

		DO ta = 0, 1
			
			dim_acc = 0
			num = 1
		
			! Initialization of the radial density
		
			DO ipoint = 0, Npoint
				DensityRadial(ta, ipoint) = 0.0
			END DO
		
			DO a = 0, 2*Lmax
		
				dim_cur = DIM(a)
			
				la = L(a)
				ja = J(a)
			
				DO i = 1, dim_cur
					an(num) = a
					ap(num) = a
					num = num + 1
				END DO

				ALLOCATE(S(dim_cur, dim_cur))
				ALLOCATE(V(dim_cur, dim_cur))
				ALLOCATE(R2Matrix(dim_cur, dim_cur))
				ALLOCATE(R2Initial(Nmax, Nmax))

				ALLOCATE(D(ta)%value(dim_cur))
				ALLOCATE(E(ta)%quantum(dim_cur, dim_cur))
				ALLOCATE(P(ta)%quantum(dim_cur, dim_cur))
				ALLOCATE(R2(ta)%quantum(dim_cur, dim_cur))
				
				! S contains the matrix of the density rho for the quantum number a = (la,ja) at the convergence

				S = diagonal%iterated%rho%rho(ta, a)%store
				
				!  After diagonalization, D(ta)% value contains the diagonal values of the p-h density rho 
				!  (the occupation probabilities v2), V the eigenvectors that make pass from the HFB basis
				!  to the canonical basis (where the matrix of rho is diagonal). 
				!
				!  The s.p. energies are by definition the expectation value of the hamiltonian in the 
				!  canonical basis. Therefore they read: (V+)*S*V

				CALL MatDiagCmplx(S, V, D(ta)%value, dim_cur)

				E(ta)%quantum = MATMUL(TRANSPOSE(V), MATMUL(SuperHamiltonian%rho%rho(ta, a)%store, V))
				P(ta)%quantum = MATMUL(TRANSPOSE(V), MATMUL(SuperHamiltonian%kap%rho(ta, a)%store, V))
				
				R2Initial = SymD3Tensor_matrix(R2Field, L(a))
				R2Matrix = R2Initial(1:dim_cur,1:dim_cur)
				
				!R2(ta)%value = MATMUL(TRANSPOSE(V), MATMUL(R2Matrix, V))
				

				DO i = 1, dim_cur
				
					Occup(ta)%value(dim_acc + i) = D(ta)%value(i)
					
					SingleParticleE(ta)%value(dim_acc + i) = E(ta)%quantum(i, i)
					SingleParticleP(ta)%value(dim_acc + i) = P(ta)%quantum(i, i)
					SingleParticleR2(ta)%value(dim_acc + i) = R2(ta)%quantum(i, i)
					
					! We calculate the contribution to the radial density of this (la, ja) state. 
					! The state has the degeneracy 2j + 1, hence the factor ja+1 (ja=2j by definition).

					IF (Basis .EQ. 2) THEN

						DO ipoint = 1, Npoint

							RadialWaveFunction(ta, ipoint, dim_acc + i) = 0.0

							DO m = 1, dim_cur
							
								IndexBra = IndexVecNL(m, la)
							
								RadialWaveFunction(ta, ipoint, dim_acc + i) = &
									 RadialWaveFunction(ta, ipoint, dim_acc + i) &
									+ V(m, i)*WaveFun(ipoint,IndexBra) / RadMesh(ipoint)
							
								DO n = 1, dim_cur
								
									IndexKet = IndexVecNL(n, la)
								
									DensityRadial(ta, ipoint) = DensityRadial(ta, ipoint) &
										+ D(ta)%value(i) * (ja + 1.0) &
										* V(m, i) * V(n, i) &
										* WaveFun(ipoint,IndexBra) &
										* WaveFun(ipoint,IndexKet)
								
								END DO
						
							END DO
						
						END DO

					END IF

				END DO

				DEALLOCATE(D(ta)%value)
				DEALLOCATE(E(ta)%quantum)
				DEALLOCATE(P(ta)%quantum)
				DEALLOCATE(R2(ta)%quantum)

				DEALLOCATE(S)
				DEALLOCATE(V)
				DEALLOCATE(R2Matrix)
				DEALLOCATE(R2Initial)

				dim_acc = dim_acc + dim_cur
			
			END DO ! end of loop over a

			IF (Basis .EQ. 2) THEN

				! The basis wave-functions are in fact y(r) = R(r)/r, so we need to correct for this.
				! The 4*pi comes from the integration over the angles.

				DO ipoint=1, Npoint
					DensityRadial(ta, ipoint) = DensityRadial(ta, ipoint) / (4.0*PI*RadMesh(ipoint)**2)
				END DO

				! Extrapolate to find rho(r=0)

                                DensityRadial(ta, 0) = 3.0*(DensityRadial(ta, 1) - DensityRadial(ta, 2)) + DensityRadial(ta, 3)
	
			END IF

		END DO ! end of loop over ta

		! Sorting the single-particle energies by ascending order
		! in() and ip() are the vectors that contain the information on the ordering

		ALLOCATE(in(dim_acc))
		ALLOCATE(ip(dim_acc))

		DO i = 1, dim_acc
			in(i) = i
			ip(i) = i
		END DO

		CALL indexx_real8(dim_acc, REAL(SingleParticleE(1)%value), in)
		CALL indexx_real8(dim_acc, REAL(SingleParticleE(0)%value), ip)

		IF (Basis .EQ. 2) THEN

			! Total density (proton and neutron)

			OPEN(file_unit_1, FILE='data/DensityCB.dat', ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				WRITE(*,'("Impossible to open the file data/DensityCB.dat")') 
				STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
			END IF
		
			DO ipoint = 0, Npoint
				WRITE(file_unit_1,'(3f20.16)') RadMesh(ipoint),DensityRadial(0, ipoint),DensityRadial(1, ipoint)
			END DO
		
			CLOSE(file_unit_1)
		
			! Radial wavefunctions. The user specifies which one to be stored by the index IndexWave

			OPEN(file_unit_1, FILE='data/WavesCB.dat', ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				WRITE(*,'("Impossible to open the file data/WavesCB.dat")') 
				STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
			END IF
		
			RadialWaveFunction(0, 0, ip(1)) = 3.0*(RadialWaveFunction(0, 1, ip(IndexWave)) & 
								- RadialWaveFunction(0, 2, ip(IndexWave))) &
								+ RadialWaveFunction(0, 3, ip(IndexWave))
			RadialWaveFunction(1, 0, in(1)) = 3.0*(RadialWaveFunction(1, 1, in(IndexWave)) &
								- RadialWaveFunction(1, 2, in(IndexWave))) &
								+ RadialWaveFunction(1, 3, in(IndexWave))
		
			DO ipoint = 0, Npoint
				WRITE(file_unit_1,'(3f20.16)') RadMesh(ipoint), RadialWaveFunction(0, ipoint, ip(IndexWave)), &
									RadialWaveFunction(1, ipoint, in(IndexWave))
			END DO
		
			CLOSE(file_unit_1)

		END IF
		
		! Storing single-particle wave-functions in the canonical basis for the neutrons

		WRITE(file_neut,'("data/HF_sp_n.dat")')

		OPEN(file_unit_1, FILE=file_neut, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			WRITE(*,'("Impossible to open the file ",A)') file_neut
			STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
		END IF
		
		WRITE(*,'("NEUTRON SINGLE-PARTICLE ENERGIES (CANONICAL BASIS)")')

		ALLOCATE(Count(0:29,0:100))
		
		DO Lvalue = 0, 29
			DO Jvalue = 0, 100
				Count(Lvalue, Jvalue) = 0
			END DO
		END DO
                
		i = 1
		
		DO WHILE (REAL(SingleParticleE(1)%value(in(i))) .GE. -100.0 .AND. REAL(SingleParticleE(1)%value(in(i))) .LE. 20.0) 
		
			Lvalue = L(an(in(i)))
			Jvalue = J(an(in(i)))
			
			Count(Lvalue, Jvalue) = Count(Lvalue, Jvalue) + 1
		
			IF (Jvalue .GT. 9) THEN
				WRITE(label,'(2x,I1,A1,I2,"/2")') Count(Lvalue, Jvalue), spectr(Lvalue), Jvalue
			ELSE
				WRITE(label,'(3x,I1,A1,I1,"/2")') Count(Lvalue, Jvalue), spectr(Lvalue), Jvalue
			END IF
			
                        WRITE(*,'(I4,")",2X,A2,I2,"/2",F6.2,F10.3,F8.3,F8.3)') i, spectr(Lvalue), Jvalue, Occup(1)%value(in(i)), SingleParticleE(1)%value(in(i))
                        WRITE(file_unit_1,'(i4,3x,a8,f15.8)') i, label, SingleParticleE(1)%value(in(i))
			
			i = i+1
		
		END DO
		
		CLOSE(file_unit_1)
		
		! Storing single-particle wave-functions in the canonical basis for the protons

		WRITE(file_prot,'("data/HF_sp_p.dat")')

		OPEN(file_unit_2, FILE=file_prot, ACTION="WRITE", IOSTAT=file_error)		
		IF (file_error .NE. 0) THEN
			WRITE(*,'("Impossible to open the file ",A)') file_prot
			STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
		END IF

		WRITE(*,'("PROTON SINGLE-PARTICLE ENERGIES (CANONICAL BASIS)")')
                
		DO Lvalue = 0, 29
			DO Jvalue = 0, 100
				Count(Lvalue, Jvalue) = 0
			END DO
		END DO
                
		i = 1
		
		DO WHILE (REAL(SingleParticleE(0)%value(ip(i))) .GE. -100.0 .AND. REAL(SingleParticleE(0)%value(ip(i))) .LE. 20.0) 
		
			Lvalue = L(ap(in(i)))
			Jvalue = J(ap(in(i)))
			
			Count(Lvalue, Jvalue) = Count(Lvalue, Jvalue) + 1
		
			IF (Jvalue .GT. 9) THEN
				WRITE(label,'(2x,I1,A1,I2,"/2")') Count(Lvalue, Jvalue), spectr(Lvalue), Jvalue
			ELSE
				WRITE(label,'(3x,I1,A1,I1,"/2")') Count(Lvalue, Jvalue), spectr(Lvalue), Jvalue
			END IF
                        WRITE(*,'(I4,")",2X,A2,I2,"/2",F6.2,F10.3,F8.3,F8.3)') i, spectr(Lvalue), Jvalue, Occup(0)%value(ip(i)), SingleParticleE(0)%value(ip(i))
                        WRITE(file_unit_2,'(i4,3x,a8,f15.8)') i, label, SingleParticleE(0)%value(ip(i))
			
			i = i+1
		
		END DO
		
		CLOSE(file_unit_2)
		
		DEALLOCATE(Count)

		DEALLOCATE(in)
		DEALLOCATE(ip)

		DEALLOCATE(Occup(0)%value)
		DEALLOCATE(Occup(1)%value)
		
		DEALLOCATE(SingleParticleE(0)%value)
		DEALLOCATE(SingleParticleE(1)%value)
		
		DEALLOCATE(SingleParticleP(0)%value)
		DEALLOCATE(SingleParticleP(1)%value)
		
		DEALLOCATE(SingleParticleR2(0)%value)
		DEALLOCATE(SingleParticleR2(1)%value)

		DEALLOCATE(an)
		DEALLOCATE(ap)

		DEALLOCATE(DensityRadial)
		DEALLOCATE(RadialWaveFunction)
		
		RETURN
	END SUBROUTINE DiagonalizationMethod_show_ParticleEnergies

	!-------------------------------------------------------------------------------!
	!										!
	!    This subroutine calculates and prints the s.p. quasi-particle energies	!
	!    in the HFB basis. It also gives an option to store the radial component	!
	!    of the q.p. wave-functions (both the upper part U and lower part V).	!
	!										!
	!-------------------------------------------------------------------------------!
	
	SUBROUTINE DiagonalizationMethod_show_QuasiParticleEnergies(diagonal)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal

		INTEGER, DIMENSION(:), ALLOCATABLE :: nIndx, pIndx
		
                INTEGER, DIMENSION(:, :), ALLOCATABLE :: Momentum

                DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: nVV, pVV
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: U, V, QPE
		DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: QPWaveFunctionLower, QPWaveFunctionUpper

                INTEGER :: NumberOfStates, a, na, top
                INTEGER :: num, la, ja, d, s1, ta, sa, an, ap, jn, jp
                INTEGER :: IndexBasis, IndexKet, nradial, ipoint, IndxQP
		INTEGER :: u1, u2, s, loop, file_error
		INTEGER :: m, n, i, state
                
		INTEGER, PARAMETER :: file_unit_1 = 16, file_unit_2 = 17
		
		DOUBLE PRECISION :: sumvv, test
                
		CHARACTER :: cn, cp
                
		CHARACTER, DIMENSION(0:19) :: spectr
                
		DATA spectr / "s", "p", "d", "f", "g", "h", "i", "j", "k", "l", &
			"m", "n", "o", "p", "q", "r", "s", "t", "u", "v" /
		
		NumberOfStates = 0;
		DO a = 0, 2*Lmax
			NumberOfStates = NumberOfStates + DIM(a)
		END DO

		! Allocate memory for the quasi-particle energies, occupation
                ! probabilities and angular momentum
                
		ALLOCATE(nVV(NumberOfStates))
		ALLOCATE(pVV(NumberOfStates))		

		ALLOCATE(nIndx(NumberOfStates))
		ALLOCATE(pIndx(NumberOfStates))

		ALLOCATE(QPE(0:1, 1:NumberOfStates))
		ALLOCATE(Momentum(0:1, 1:NumberOfStates))
		
		ALLOCATE(QPWaveFunctionLower(0:1, 0:Npoint, 1:NumberOfStates))
		ALLOCATE(QPWaveFunctionUpper(0:1, 0:Npoint, 1:NumberOfStates))
		
		! Initialization of the neutron and proton radial density
		
		IF (.NOT. ALLOCATED(DensityRadial)) ALLOCATE(DensityRadial(0:1, 0:Npoint))

		! Alamcenar las energias de cuasiparticulas y sus ocupaciones
		
		DO ta = 0, 1
			
			DO ipoint = 0, Npoint
				DensityRadial(ta, ipoint) = 0.0
			END DO
		
			num = 0
			state = 0
                
			DO a = 0, 2*Lmax
		
				la = L(a)
				ja = J(a)
									d = MIN(Nmax, NmaxOfL(la))
				IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - la) / 2) + 1
			
				! Collecting the q.-p. energies and calculating the occupation factors
                        
				DO s1 = d + 1, d + d
				
                                	num = num + 1
				
					sumvv = 0.0
                                
					QPE(ta, num) = REAL(diagonal%QuasiParticleEnergies(ta, a)%value(s1))
                                
					Momentum(ta, num) = a
                                
					IF (ta .EQ. 1) nIndx(num) = num
					IF (ta .EQ. 0) pIndx(num) = num
                                
					DO sa = 1, d
						sumvv = sumvv + (diagonal%UV(ta, a)%quantum(sa + d, s1) ** 2)
					END DO

					IF ((diagonal%consistency%density%nucleus%is_blocking(ta)) .AND. &
						(a .EQ. diagonal%consistency%density%nucleus%ia(ta)) .AND. &
						(s1 .EQ. (d + diagonal%consistency%density%nucleus%mu0(ta)))) THEN

						DO sa = 1, d
							sumvv = sumvv +((diagonal%UV(ta, a)%quantum(sa, s1) ** 2) &
								      - (diagonal%UV(ta, a)%quantum(sa + d, s1) ** 2)) / (ja + 1.0)
						END DO

					END IF
                                
					IF (ta .EQ. 1) nVV(num) = sumvv
					IF (ta .EQ. 0) pVV(num) = sumvv
					
				END DO

				! Collecting the wave-functions: upper part (corresponding to the vector U) and lower part
				! (corresponding to the vector V)
					
				IF (Basis .EQ. 2) THEN
				
					! First extracting the matrices U and V
			
					ALLOCATE(U(d, d))
					ALLOCATE(V(d, d))

					DO u1 = 1, d
						DO u2 = 1, d
							U(u1, u2) = diagonal%UV(ta, a)%quantum(u1    , u2 + d)
							V(u1, u2) = diagonal%UV(ta, a)%quantum(u1 + d, u2 + d)
						END DO
					END DO
						
					DO i = 1, d

						state = state + 1
				
						DO ipoint = 1, Npoint

							QPWaveFunctionUpper(ta, ipoint, state) = 0.0
							QPWaveFunctionLower(ta, ipoint, state) = 0.0

							DO m = 1, d
                                        
								IndexBasis = IndexVecNL(m, la)

								QPWaveFunctionUpper(ta, ipoint, state) = &
								QPWaveFunctionUpper(ta, ipoint, state) &
								+ U(m, i) * WaveFun(ipoint,IndexBasis)/RadMesh(ipoint)

								QPWaveFunctionLower(ta, ipoint, state) = &
								QPWaveFunctionLower(ta, ipoint, state) &
								+ V(m, i) * WaveFun(ipoint,IndexBasis)/RadMesh(ipoint)

								! Calculating here the single-particle density. We have: 
								!
								!		\rho_{mn} =  { V (V+) }_{mn}
								!
								! and in coordinate representation:
								!
								!          rho(r) = \sum_{mn,l,j} (2j+1) \rho_{mn} \phi_{ml}(r)\phi_{nl}(r)
								!
								! or:
								!     rho(r) = \sum_{i, mn,l,j} (2j+1) V_{ni}V_{mi} \phi_{ml}(r)\phi_{nl}(r)
								!
					
								DO n = 1, d
								
									IndexKet = IndexVecNL(n, la)
								
									DensityRadial(ta, ipoint) = DensityRadial(ta, ipoint) &
												+ (ja + 1.0) * V(m, i) * V(n, i) &
												* WaveFun(ipoint,IndexBasis) &
												* WaveFun(ipoint,IndexKet)
								
								END DO
							
							END DO
                                        
						END DO
                                        
					END DO
                                        
					DEALLOCATE(U)
					DEALLOCATE(V)
					
				END IF
						                                
			END DO ! end of loop over a
			
 1              	CONTINUE   
 
             		! Taking into account the fact that phi(r) = y(r)/r, and extrapolating at r=0 for the density
			
			IF (Basis .EQ. 2) THEN

				! The basis wave-functions are in fact y(r) = R(r)/r, so we need to correct for this.
				! The 4*pi comes from the integration over the angles.

				DO ipoint=1, Npoint
					DensityRadial(ta, ipoint) = DensityRadial(ta, ipoint) / (4.0*PI*RadMesh(ipoint)**2)
				END DO

				! Extrapolate to find rho(r=0)

         			DensityRadial(ta, 0) = 3.0*(DensityRadial(ta, 1) - DensityRadial(ta, 2)) + DensityRadial(ta, 3)

			END IF
		
		END DO ! end of loop over ta
					
		! Ordering of the quasi-particle occupation numbers
		
		CALL indexx_real8(NumberOfStates, nVV, nIndx)
		CALL indexx_real8(NumberOfStates, pVV, pIndx)

		WRITE(*,'(A14,A10,A12,A5,A10)') "EQP NEUT", "V2", "", "V2", "EQP PROT"
		
		top = 100
		
		DO num = top, 1, -1

			an = Momentum(1, nIndx(num))
			ap = Momentum(0, pIndx(num))
			
			cn = spectr(L(an))                       
			cp = spectr(L(ap))
			
			jn = J(an)
			jp = J(ap)

                        WRITE(*,'(2i5,F10.4,F10.4,3X,A1,I2,"/2",A5,I2,"/2",3x,F6.4,F10.4)') &
				num,nIndx(num),QPE(1, nIndx(num)), nVV(nIndx(num)), cn, jn, cp, jp, pVV(pIndx(num)), QPE(0, pIndx(num))
		END DO
		
		IF (Basis .EQ. 2) THEN
				
			! Total density (proton and neutron)

			OPEN(file_unit_1, FILE='data/DensityQP.dat', ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				WRITE(*,'("Impossible to open the file data/DensityQP.dat")') 
				STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
			END IF
		
			DO ipoint = 0, Npoint
				WRITE(file_unit_1,'(3f20.16)') RadMesh(ipoint),DensityRadial(0, ipoint),DensityRadial(1, ipoint)
			END DO
		
			CLOSE(file_unit_1)
		
			! Storing quasi-particle wave-functions (proton and neutrons).

			OPEN(file_unit_1, FILE='data/WavesQP.dat', ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				WRITE(*,'("Impossible to open the file data/WavesQP.dat")') 
				STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
			END IF

			IndxQP = NumberOfStates - IndexWave + 1

			write(*,'("IndexWave = ",i4," IndxQP = ",i4," nIndx = ",i4)') IndexWave, IndxQP, nIndx(IndxQP)

			QPWaveFunctionLower(0, 0, pIndx(IndxQP)) = 3.0*(QPWaveFunctionLower(0, 1, pIndx(IndxQP)) &
								 	- QPWaveFunctionLower(0, 2, pIndx(IndxQP))) &
						 			+ QPWaveFunctionLower(0, 3, pIndx(IndxQP))
 
			QPWaveFunctionUpper(0, 0, pIndx(IndxQP)) = 3.0*(QPWaveFunctionUpper(0, 1, pIndx(IndxQP)) &
								  	- QPWaveFunctionUpper(0, 2, pIndx(IndxQP))) &
								  	+ QPWaveFunctionUpper(0, 3, pIndx(IndxQP))

			QPWaveFunctionLower(1, 0, nIndx(IndxQP)) = 3.0*(QPWaveFunctionLower(1, 1, nIndx(IndxQP)) &
								 	- QPWaveFunctionLower(1, 2, nIndx(IndxQP))) &
						 			+ QPWaveFunctionLower(1, 3, nIndx(IndxQP))
 
			QPWaveFunctionUpper(1, 0, nIndx(IndxQP)) = 3.0*(QPWaveFunctionUpper(1, 1, nIndx(IndxQP)) &
								  	- QPWaveFunctionUpper(1, 2, nIndx(IndxQP))) &
								  	+ QPWaveFunctionUpper(1, 3, nIndx(IndxQP))

			DO ipoint = 0, Npoint
				WRITE(file_unit_1,'(5f20.16)') RadMesh(ipoint), QPWaveFunctionUpper(0, ipoint, pIndx(IndxQP)), &
									QPWaveFunctionLower(0, ipoint, pIndx(IndxQP)), &
									QPWaveFunctionUpper(1, ipoint, nIndx(IndxQP)), &
									QPWaveFunctionLower(1, ipoint, nIndx(IndxQP))
			END DO

			CLOSE(file_unit_1)

		END IF

		DEALLOCATE(nIndx)
		DEALLOCATE(pIndx)

		DEALLOCATE(QPE)

		DEALLOCATE(nVV)
		DEALLOCATE(pVV)

		DEALLOCATE(Momentum)

		DEALLOCATE(QPWaveFunctionUpper)
		DEALLOCATE(QPWaveFunctionLower)
		
		RETURN
                
	END SUBROUTINE DiagonalizationMethod_show_QuasiParticleEnergies

	!-------------------------------------------------------------------------------!
	!										!
	!            CALCULATION OF THE PROJECTED ENERGY (AFTER VARIATION) 		!
	!										!
	!  We calculate here the particle number-projected energies after variation. 	!
	!  All components of the force in the p-h and p-p channel are treated separate-	!
	!  ly). The general formula implemented here corresponds to the following 	!
	!  reference, Eq. 21 in:							!
	!				M. Anguinao, J. L. Egido, L. M. Robledo		!
	!				Nucl. Phys. A696 (2001) 467-493			!
	!										!
	!  For the density-dependent term, the projected density prescription is used.	!
	!  The one-body and two-body kinetic energy terms, as well as the Coulomb and	!
	!  density-dependent term do not mix isospin (= no term in Gamma(nu)rho(pi)).	!
	!										!
	!-------------------------------------------------------------------------------!
	
	SUBROUTINE DiagonalizationMethod_ProjectedEnergy(diagonal)
		TYPE (DiagonalizationMethod), INTENT(IN) :: diagonal
		
		TYPE (SymGenDensityGaugeProj), DIMENSION(:), POINTER :: genden_gauge
		TYPE (SymGenDensityProj) :: genden_proj, SuperHamiltonian
		TYPE (SymGenDensityHFProj), DIMENSION(:), POINTER :: C_Matrix, Deriv_Y
		TYPE (SymDensityProj) :: density, density01, density10
		TYPE (SymHartreeFockFieldProj) :: HF_Gamma_Diff, HF_Gamma_Same, HF_Delta
		TYPE (ProjectionCoeffs) :: CoeffsXY
		TYPE (SelfConsistencyMethodProj) :: consistency
		
		COMPLEX, DIMENSION(0:1) :: P
		COMPLEX :: ProjectedEkCM, ProjectedEk
		COMPLEX :: ProjectedVBB_local, ProjectedVBB_exch, ProjectedVC_local, ProjectedVC_exch, ProjectedVLS, ProjectedDD
		COMPLEX :: ProjectedVBB_Pair, ProjectedVC_Pair, ProjectedVLS_Pair
				
		DOUBLE PRECISION :: A, Gauge, factor
		
		INTEGER :: Gauge_n, Gauge_p, Z, N, IndexGauge
		INTEGER :: ProjectionOn = 1
			
		A = Nucleus_get_A(diagonal%consistency%density%nucleus)
		
		IF (A .LE. 1) STOP "Abortado"
		
		Z = diagonal%consistency%density%nucleus%np(0)
		N = diagonal%consistency%density%nucleus%np(1)
		
		! Create new Fomenko coefficients
		
		CALL ProjectionCoeffs_new(CoeffsXY, NGauge)
		
		! Calculate phi-dependent densities rho(phi), kappa10(phi), kappa01(phi) according
		! to formulas (9), (10) and (11) of the Anguiano et al, Nucl. Phys. A696 (2001) 467-493
		! All densities are stored into genden_gauge.
		
		CALL SymGenDensityProj_make_DensityGauge(diagonal%UV, genden_gauge, NGauge)
		
		! Calculate the projected density and the Fomenko coefficients
		
		CALL SymGenDensityProj_make_DensityProj(genden_proj, CoeffsXY, genden_gauge, diagonal%UV, &
							diagonal%consistency%density%nucleus, NGauge)
		
		!****************************************************************************************************!
		
		! Calculate the c-matrix and the derivative of the y-coefficients with respect to the density matrix
		CALL ProjectionPreparation(diagonal%UV, CoeffsXY, genden_gauge, C_Matrix, Deriv_Y, NGauge)
		
		!****************************************************************************************************!
		
		CALL SymHartreeFockFieldProj_new(HF_Gamma_Same)
		
		! Convert projected density into the right type, which is SymDensityProj
		
		CALL SymDensity_new_GenDensityProj(density, genden_proj)
		
		!---------------------------------------------------------------------------------------!
		! 		Terms depending on the projected density rhoP				!
		!---------------------------------------------------------------------------------------!		
		
                factor = (1.0 - (1.0 / A))
		
		ProjectedEk = factor * (EkField * density%field%rho%p(0) + EkField * density%field%rho%p(1))
		
		! Two-body kinetic energy
		CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma_Same, diagonal%consistency%vEkCMph, density%field%rho, 0, 1, ProjectionOn)
		CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma_Same, diagonal%consistency%vEkCMph, density%field%rho, 1, 0, ProjectionOn)
		
		P(0) = SymHartreeFockFieldProj_product_iso(density%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density%field%rho, HF_Gamma_Same, 1)
		
		ProjectedEkCM = (0.5 / A) * (P(0) + P(1))
				
		! Brink-Boeker direct terms
		CALL SymVBBph_get_LocalGamma(HF_Gamma_Same, diagonal%consistency%vBBph, density%field%rho, 0, 1, ProjectionOn)
		CALL SymVBBph_get_LocalGamma(HF_Gamma_Same, diagonal%consistency%vBBph, density%field%rho, 1, 0, ProjectionOn)
		
		P(0) = SymHartreeFockFieldProj_product_iso(density%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density%field%rho, HF_Gamma_Same, 1)
		
		ProjectedVBB_local = 0.5 * (P(0) + P(1))
		
		! Brink-Boeker exchange terms
		CALL SymVBBph_get_ExchangeGamma(HF_Gamma_Same, diagonal%consistency%vBBph, density%field%rho, 0, 1, ProjectionOn)
		CALL SymVBBph_get_ExchangeGamma(HF_Gamma_Same, diagonal%consistency%vBBph, density%field%rho, 1, 0, ProjectionOn)
		
		P(0) = SymHartreeFockFieldProj_product_iso(density%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density%field%rho, HF_Gamma_Same, 1)
		
		ProjectedVBB_exch = 0.5 * (P(0) + P(1))
		
		! Coulomb direct terms
		CALL SymVCph_get_LocalGamma(HF_Gamma_Same, diagonal%consistency%vCph, density%field%rho, 0, 1, ProjectionOn)
		CALL SymVCph_get_LocalGamma(HF_Gamma_Same, diagonal%consistency%vCph, density%field%rho, 1, 0, ProjectionOn)

		P(0) = SymHartreeFockFieldProj_product_iso(density%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density%field%rho, HF_Gamma_Same, 1)
		
		ProjectedVC_local = 0.5 * (P(0) + P(1))
		
		! Coulomb exchange terms
		CALL SymVCph_get_ExchangeGamma(HF_Gamma_Same, diagonal%consistency%vCph, density%field%rho, 0, 1, ProjectionOn)
		CALL SymVCph_get_ExchangeGamma(HF_Gamma_Same, diagonal%consistency%vCph, density%field%rho, 1, 0, ProjectionOn)

		P(0) = SymHartreeFockFieldProj_product_iso(density%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density%field%rho, HF_Gamma_Same, 1)
		
		ProjectedVC_exch = 0.5 * (P(0) + P(1))
		
		! Spin-orbit term
		CALL SymVLSph_get_Gamma(HF_Gamma_Same, diagonal%consistency%vLSph, density%field%rho, 0, 1, ProjectionOn)
		CALL SymVLSph_get_Gamma(HF_Gamma_Same, diagonal%consistency%vLSph, density%field%rho, 1, 0, ProjectionOn)

		P(0) = SymHartreeFockFieldProj_product_iso(density%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density%field%rho, HF_Gamma_Same, 1)
		
		ProjectedVLS = 0.5 * (P(0) + P(1))
			
		! Density-dependent projected energy
		
		Gauge = 0.0
		
		CALL SelfConsistencyMethodProj_new(consistency, density)
		
		IF (Basis .EQ. 1) THEN
			CALL SymGDDph_make_DD(consistency%gDDph, density%field%rho, Gauge)
		ELSE
			CALL Make_DenGenFun(consistency%gDDph, density%field%rho, Gauge)
		END IF
		
		ProjectedDD = SymGDDph_get_edd(consistency%gDDph)
		
		!---------------------------------------------------------------------------------------!
		! 	Terms depending on the Gauge densities rho(phi), kap10(phi), kap01(phi)		!
		!---------------------------------------------------------------------------------------!		
		
		ProjectedVBB_Pair  = CMPLX(0,0)
		ProjectedVC_Pair   = CMPLX(0,0)
		ProjectedVLS_Pair  = CMPLX(0,0)
				
		DO IndexGauge = 1, NGauge
		
			! Convert projected density into the right type, which is SymDensityProj
		
			CALL SymDensity_new_GenDensityProj10(density10, genden_gauge(IndexGauge))
			
			CALL SymHartreeFockFieldProj_new(HF_Gamma_Diff)
			CALL SymHartreeFockFieldProj_new(HF_Delta)
		
			!
			! Mean-field Channel
			!
			
			! Two-body kinetic energy (center of mass motion)
			CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vEkCMph, density10%field%rho, 0, 0, ProjectionOn)
			CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vEkCMph, density10%field%rho, 1, 1, ProjectionOn)
		
			P(0) = SymHartreeFockFieldProj_product_iso(density10%field%rho, HF_Gamma_Diff, 0)
			P(1) = SymHartreeFockFieldProj_product_iso(density10%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedEkCM = ProjectedEkCM + (0.5 / A) * (CoeffsXY%y_l(IndexGauge, 0)*P(0) + CoeffsXY%y_l(IndexGauge, 1)*P(1))
		
			! Brink-Boeker direct term
			CALL SymVBBph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10%field%rho, 0, 0, ProjectionOn)
			CALL SymVBBph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10%field%rho, 1, 1, ProjectionOn)
		
			P(0) = SymHartreeFockFieldProj_product_iso(density10%field%rho, HF_Gamma_Diff, 0)
			P(1) = SymHartreeFockFieldProj_product_iso(density10%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVBB_local = ProjectedVBB_local + 0.5 * (CoeffsXY%y_l(IndexGauge, 0)*P(0) + CoeffsXY%y_l(IndexGauge, 1)*P(1))
		
			! Brink-Boeker exchange term
			CALL SymVBBph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10%field%rho, 0, 0, ProjectionOn)
			CALL SymVBBph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10%field%rho, 1, 1, ProjectionOn)
		
			P(0) = SymHartreeFockFieldProj_product_iso(density10%field%rho, HF_Gamma_Diff, 0)
			P(1) = SymHartreeFockFieldProj_product_iso(density10%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVBB_exch = ProjectedVBB_exch + 0.5 * (CoeffsXY%y_l(IndexGauge, 0)*P(0) + CoeffsXY%y_l(IndexGauge, 1)*P(1))
			
			! Coulomb direct term
			CALL SymVCph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10%field%rho, 0, 0, ProjectionOn)
			CALL SymVCph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10%field%rho, 1, 1, ProjectionOn)

			P(0) = SymHartreeFockFieldProj_product_iso(density10%field%rho, HF_Gamma_Diff, 0)
			P(1) = SymHartreeFockFieldProj_product_iso(density10%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVC_local = ProjectedVC_local + 0.5 * (CoeffsXY%y_l(IndexGauge, 0)*P(0) + CoeffsXY%y_l(IndexGauge, 1)*P(1))
			
			! Coulomb exchange term
			CALL SymVCph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10%field%rho, 0, 0, ProjectionOn)
			CALL SymVCph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10%field%rho, 1, 1, ProjectionOn)

			P(0) = SymHartreeFockFieldProj_product_iso(density10%field%rho, HF_Gamma_Diff, 0)
			P(1) = SymHartreeFockFieldProj_product_iso(density10%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVC_exch = ProjectedVC_exch + 0.5 * (CoeffsXY%y_l(IndexGauge, 0)*P(0) + CoeffsXY%y_l(IndexGauge, 1)*P(1))
		
			! Spin-orbit term
			CALL SymVLSph_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vLSph, density10%field%rho, 0, 0, ProjectionOn)
			CALL SymVLSph_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vLSph, density10%field%rho, 1, 1, ProjectionOn)

			P(0) = SymHartreeFockFieldProj_product_iso(density10%field%rho, HF_Gamma_Diff, 0)
			P(1) = SymHartreeFockFieldProj_product_iso(density10%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVLS = ProjectedVLS + 0.5 * (CoeffsXY%y_l(IndexGauge, 0)*P(0) + CoeffsXY%y_l(IndexGauge, 1)*P(1))
		
			!
			! Pairing Channel
			!
			
			CALL SymDensity_new_GenDensityProj01(density01, genden_gauge(IndexGauge))
			
			!density01 = diagonal%consistency%density
			
			! Brink-Boeker Pairing
			CALL SymVBBpp_get_Delta(HF_delta, diagonal%consistency%vBBpp, density10%field%kap, 0, 0)
			CALL SymVBBpp_get_Delta(HF_delta, diagonal%consistency%vBBpp, density10%field%kap, 1, 1)
			
			P(0) = SymHartreeFockFieldProj_product_iso(density01%field%kap, HF_delta, 0)
			P(1) = SymHartreeFockFieldProj_product_iso(density01%field%kap, HF_delta, 1)
		
			ProjectedVBB_Pair = ProjectedVBB_Pair - 0.5 * (CoeffsXY%y_l(IndexGauge, 0)*P(0) + CoeffsXY%y_l(IndexGauge, 1)*P(1))
		
			! Coulomb Pairing
			CALL SymVCpp_get_Delta(HF_delta, diagonal%consistency%vCpp, density10%field%kap, 0, 0)
			CALL SymVCpp_get_Delta(HF_delta, diagonal%consistency%vCpp, density10%field%kap, 1, 1)
			
			P(0) = SymHartreeFockFieldProj_product_iso(density01%field%kap, HF_delta, 0)
			P(1) = SymHartreeFockFieldProj_product_iso(density01%field%kap, HF_delta, 1)
		
			ProjectedVC_Pair = ProjectedVC_Pair - 0.5 * (CoeffsXY%y_l(IndexGauge, 0)*P(0) + CoeffsXY%y_l(IndexGauge, 1)*P(1))
		
			! Spin-Orbit Pairing
			CALL SymVLSpp_get_Delta(HF_delta, diagonal%consistency%vLSpp, density10%field%kap, 0, 0)
			CALL SymVLSpp_get_Delta(HF_delta, diagonal%consistency%vLSpp, density10%field%kap, 1, 1)
			
			P(0) = SymHartreeFockFieldProj_product_iso(density01%field%kap, HF_delta, 0)
			P(1) = SymHartreeFockFieldProj_product_iso(density01%field%kap, HF_delta, 1)
		
			ProjectedVLS_Pair = ProjectedVLS_Pair - 0.5 * (CoeffsXY%y_l(IndexGauge, 0)*P(0) + CoeffsXY%y_l(IndexGauge, 1)*P(1))
		
			CALL SymDensityProj_del(density01)
			CALL SymDensityProj_del(density10)
			
			CALL SymHartreeFockFieldProj_del(HF_Gamma_Diff)
			CALL SymHartreeFockFieldProj_del(HF_Delta)

		END DO
		
		WRITE(*,'(/,"ISOSPIN-INDEPENDENT TERMS",/)') 
		WRITE(*,'("KINETIC ENERGY            : ",2f17.10)') ProjectedEk
		WRITE(*,'("CENTER OF MASS CORRECTION : ",2f17.10)') ProjectedEkCM
		WRITE(*,'("BRINK-BOEKER (LOCAL)      : ",2f17.10)') ProjectedVBB_local
		WRITE(*,'("BRINK-BOEKER (EXCHANGE)   : ",2f17.10)') ProjectedVBB_exch
		WRITE(*,'("COULOMB (LOCAL)           : ",2f17.10)') ProjectedVC_local 
		WRITE(*,'("COULOMB (EXCHANGE)        : ",2f17.10)') ProjectedVC_exch
		WRITE(*,'("SPIN-ORBIT                : ",2f17.10)') ProjectedVLS
		WRITE(*,'("DENSITY-DEPENDENT         : ",2f17.10)') ProjectedDD
		WRITE(*,'("TOTAL                     : ",2f17.10,/)') ProjectedDD+ProjectedVLS+ProjectedVC_exch+ProjectedVC_local & 
								    +ProjectedVBB_exch+ProjectedVBB_local+ProjectedEkCM+ProjectedEk
								  
		WRITE(*,'("PAIRING - BRINK-BOEKER    : ",2f17.10)') ProjectedVBB_Pair
		WRITE(*,'("PAIRING - COULOMB         : ",2f17.10)') ProjectedVC_Pair
		WRITE(*,'("PAIRING - SPIN-ORBIT      : ",2f17.10,/)') ProjectedVLS_Pair
		WRITE(*,'("GRAND TOTAL               : ",2f17.10,/)') ProjectedDD+ProjectedVLS+ProjectedVC_exch+ProjectedVC_local & 
								    +ProjectedVBB_exch+ProjectedVBB_local+ProjectedEkCM+ProjectedEk &
								    +ProjectedVBB_Pair+ProjectedVC_Pair+ProjectedVLS_Pair
		CALL SymHartreeFockFieldProj_del(HF_delta)				

		RETURN
	END SUBROUTINE DiagonalizationMethod_ProjectedEnergy

	SUBROUTINE DiagonalizationMethod_del(diagonal)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal

		INTEGER ta, a, max_a

		max_a = 2*Lmax
		DO ta = 0, 1
			DO a = 0, max_a
				DEALLOCATE(diagonal%QuasiParticleEnergies(ta, a)%value)
				DEALLOCATE(diagonal%UV(ta, a)%quantum)
			END DO
		END DO
		DEALLOCATE(diagonal%QuasiParticleEnergies)
		DEALLOCATE(diagonal%UV)

		CALL SelfConsistencyMethodProj_del(diagonal%consistency)
		CALL R1R1Function_del(diagonal%func)
		
		CALL SymGenDensityProj_del(diagonal%S)
		CALL SymGenDensityProj_del(diagonal%iterated)
		
		RETURN
	END SUBROUTINE DiagonalizationMethod_del

	SUBROUTINE WaveFunction_new_SymDensity(wave_func, density, tolerance)
		TYPE (WaveFunction), INTENT(INOUT) :: wave_func
		TYPE (SymDensityProj), INTENT(INOUT) :: density
		DOUBLE PRECISION, INTENT(IN) :: tolerance

		TYPE (DiagonalizationMethod) :: evolv

		CALL Nucleus_new_Nucleus(wave_func%nucleus, density%nucleus)
		
		CALL DiagonalizationMethod_new(evolv, density)
		CALL DiagonalizationMethod_goto_SelfConsistency(evolv, tolerance)
		
		wave_func = DiagonalizationMethod_get_WaveFunction(evolv)
		
		RETURN
	END SUBROUTINE WaveFunction_new_SymDensity

END MODULE diagmeth
