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
	USE indexx
	USE energy_proj
	USE eigenval
	USE jacobi
        USE global

	IMPLICIT NONE

	INCLUDE "brent_d1.f90"

	INTEGER, PARAMETER :: MAX_ITER = 5
	INTEGER :: NGauge = 9
	
	DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DensityRadial, DensityPairing, DerivDens

	TYPE DiagonalMatrix
		DOUBLE PRECISION, POINTER, DIMENSION(:) :: value
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

	!-----------------------------------------------------------------------------------------------!
	!   This subroutine performs the self-consistent HFB calculation. Firstly, it calculates the	!
	!   matrix elements. Eventually, if the latter were already calculated before, it reads the 	!
	!   files that contain them. Then, it enters the interation loops. At each step, it calculates	!
	!   the mean-field h and Delta, solves the particle number equation and saves the density thus 	!
	!   obtained and cycles. The convergence is achieved when the difference between the norms of	!
	!   the density between 2 iterations is lower than "tolerance". At the convergence, a summary 	!
	!   of the calculation is printed, then the single-particle energies in the canonical basis and	!
	!   the densities are calculated and printed, then the quasi-particle energies are printed out.	!
	!-----------------------------------------------------------------------------------------------!

	SUBROUTINE DiagonalizationMethod_goto_SelfConsistency(diagonal, tolerance)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		DOUBLE PRECISION, INTENT(IN) :: tolerance

		TYPE (OneDimSolve) :: neutron_constrainer, proton_constrainer
		TYPE (SymDensityProj) :: new_density, new_density_proj
		TYPE (SymHartreeFockFieldProj) :: field1, field2
		TYPE (SelfConsistencyMethodProj) :: consistency
		
		REAL (KIND=16) :: TestAccuracy
		
		LOGICAL :: Read_BBph, Read_BBpp, Read_Cph, Read_Cpp
		
		DOUBLE PRECISION :: b, diff, R2
		DOUBLE PRECISION :: accuracy, Gauge, factor

		INTEGER :: A, N, Z, niter, ta
		INTEGER :: cycles_in, cycles_out, cycles_rate, Lold, Lmin
		INTEGER :: IndexGauge, NumberOfGauge
		
		! Reading the system's clock
		CALL SYSTEM_CLOCK(cycles_in)

		diff = 0.0
		
		b = Nucleus_get_b(diagonal%consistency%density%nucleus)
		A = Nucleus_get_A(diagonal%consistency%density%nucleus)
		
		IF ((A .LE. 1) .OR. (A .GE. 300)) THEN
			PRINT *, "Unexpected A value = ", A
			STOP "DiagonalizationMethod_goto_SelfConsistency"
		END IF
		
		factor = 1.0 - (1.0 / A)

		WRITE(*,'(5X,"STARTING THE SELF-CONSISTENT CALCULATION")')
		WRITE(*,'(5X,"========================================",/)')
		
		WRITE(*,'("Numerical accuracy in Energy  : ",ES10.2)') tolerance
		WRITE(*,'("Oscillator length             : ",F6.3)') b
		WRITE(*,'("Atomic number A               : ",i6)') A
		WRITE(*,'("Largest Number                : ",E24.16)') HUGE(TestAccuracy)
		WRITE(*,'("Smallest Number               : ",E24.16,/)') TINY(TestAccuracy)

		!-------------------------------------------------------!
		! Calculation of the Matrix elements of the Gogny force	!
		!-------------------------------------------------------!
		
		CALL SymKineticEnergy2Body_calculate(diagonal%consistency%vEkCMph, diagonal%consistency%vEkCMpp)

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

			Read_BBph = .FALSE.
			Read_BBpp = .FALSE.
		
			IF (Basis .EQ. 2 .AND. (.NOT. Read_BBph .OR. .NOT. Read_BBpp)) THEN
                                CALL BesselTabularize(2*Lmax)
                        END IF
                        
			CALL SymVBBph_calculate(diagonal%consistency%vBBph, Read_BBph, Lmin)
			CALL SymVBBpp_calculate(diagonal%consistency%vBBpp, Read_BBpp, Lmin)
		
			IF (Basis .EQ. 2 .AND. (.NOT. Read_BBph .OR. .NOT. Read_BBpp)) THEN
                                CALL BesselFree()
                        END IF

			Read_Cph = .FALSE.
			Read_Cpp = .FALSE.
		
			CALL SymVCph_calculate(diagonal%consistency%vCph, Read_Cph, Lmin)
			CALL SymVCpp_calculate(diagonal%consistency%vCpp, Read_Cpp, Lmin)
		
		ELSE
		
			Lmin = 0
			
			IF (Basis .EQ. 2) THEN
				IF (Lmax < 10) THEN
					WRITE(diagonal%consistency%vBBph%filename, "(A,I1,A)") "data/vBB", Lmax, "ph_WS.txt"
				ELSE
					WRITE(diagonal%consistency%vBBph%filename, "(A,I2,A)") "data/vBB", Lmax, "ph_WS.txt"
				END IF
			ELSE
				IF (Lmax < 10) THEN
					WRITE(diagonal%consistency%vBBph%filename, "(A,I1,A)") "data/vBB", Lmax, "ph_HO.txt"
				ELSE
					WRITE(diagonal%consistency%vBBph%filename, "(A,I2,A)") "data/vBB", Lmax, "ph_HO.txt"
				END IF
			END IF
		
			Read_BBph = SymVBBph_read(diagonal%consistency%vBBph)
			Read_BBpp = SymVBBpp_read(diagonal%consistency%vBBpp)

                        IF (Basis .EQ. 2 .AND. (.NOT. Read_BBph .OR. .NOT. Read_BBpp)) THEN
                                CALL BesselTabularize(2*Lmax)
                        END IF                        
		
			CALL SymVBBph_calculate(diagonal%consistency%vBBph, Read_BBph, Lmin)
			CALL SymVBBpp_calculate(diagonal%consistency%vBBpp, Read_BBpp, Lmin)

                        IF (Basis .EQ. 2 .AND. (.NOT. Read_BBph .OR. .NOT. Read_BBpp)) THEN
                                CALL BesselFree()
                        END IF
		
			Read_Cph = SymVCph_read(diagonal%consistency%vCph)
			Read_Cpp = SymVCpp_read(diagonal%consistency%vCpp)
		
			CALL SymVCph_calculate(diagonal%consistency%vCph, Read_Cph, Lmin)
			CALL SymVCpp_calculate(diagonal%consistency%vCpp, Read_Cpp, Lmin)
		
		END IF
		
		! Creation of the tensor fields

		CALL SymHartreeFockFieldProj_new(field1)
		CALL SymHartreeFockFieldProj_new(field2)

		!-------------------------------------------------------!
		!     Starting the iterative process of convergence	!
		!-------------------------------------------------------!
		
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

			CALL DiagonalizationMethod_get_MeanField(diagonal, diagonal%S, Gauge)
			
			! We obtain the Fermi energies Lambda_n and Lambda_p by considering the constraint over the number of particle

			IF (HFOnly .EQ. 0) THEN

				CALL DiagonalizationMethod_set_ISOSpin(diagonal, 1)

				CALL OneDimSolve_new(neutron_constrainer, diagonal%func, diagonal)
			
				diagonal%consistency%density%nucleus%lambda_np(1) = OneDimSolve_solve(neutron_constrainer, &
					DBLE(N)-ShiftLambda, diagonal%consistency%density%nucleus%lambda_np(1), &
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

			! Store the density for an eventual restart
			CALL SymDensityProj_save(diagonal%consistency%density)

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
			
			CALL SymHartreeFockFieldProj_product(field1, anneal, diagonal%consistency%density%field%rho)
			CALL SymHartreeFockFieldProj_add(field2, new_density%field%rho, field1)
			CALL SymHartreeFockFieldProj_product(diagonal%consistency%density%field%rho, factor, field2)

			CALL SymHartreeFockFieldProj_product(field1, anneal, diagonal%consistency%density%field%kap)
			CALL SymHartreeFockFieldProj_add(field2, new_density%field%kap, field1)
			CALL SymHartreeFockFieldProj_product(diagonal%consistency%density%field%kap, factor, field2)

			! We update the numbers of particles with the same rule

			DO ta = 0, 1
				diagonal%consistency%density%nucleus%actual_np(ta) = factor * ((new_density%nucleus%actual_np(ta)) &
					+ (anneal * diagonal%consistency%density%nucleus%actual_np(ta)))
                        
			END DO

			accuracy = SelfConsistencyMethodProj_accuracy(diagonal%consistency)
			WRITE(*,'("Iteration k = ",I4," HF Energy : ",F15.8)') niter,diagonal%consistency%density%nucleus%eHFB
			
			IF (diff .LE. tolerance) THEN
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

		! Store the density for an eventual restart
		CALL SymDensityProj_save(diagonal%consistency%density)

		CALL SYSTEM_CLOCK(cycles_out, cycles_rate)
		PRINT "(/A,EN10.2)", "Tiempo empleado (segundos): ", (DBLE(cycles_out - cycles_in) / cycles_rate)

		CALL SymDensityProj_store_actual_R2(diagonal%consistency%density)

		CALL SelfConsistencyMethodProj_store_eHFB(diagonal%consistency)

                facPair = 1.0
                
		CALL SelfConsistencyMethodProj_show_Status(diagonal%consistency)
		
		CALL DiagonalizationMethod_ProjectedEnergy(diagonal)
	
		CALL DiagonalizationMethod_show_ParticleEnergies(diagonal, Gauge)
                CALL DiagonalizationMethod_show_QuasiParticleEnergies(diagonal)
		
		RETURN
	END SUBROUTINE DiagonalizationMethod_goto_SelfConsistency

	!-----------------------------------------------------------------------------------------------!
	!   In this subroutine, we calculate the mean-field Gamma and pairing field Delta from the 	!
	!   matrix elements of the force and the densities. Both the matrix elements and the densities	!
	!   at a given iteration of the HFB process are stored in the object "diagonal". Once the fields!
	!   are constructed, we build an object that will contain the super-hamiltonian S.
	!-----------------------------------------------------------------------------------------------!

	SUBROUTINE DiagonalizationMethod_get_MeanField(diagonal, S, Gauge)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		TYPE (SymGenDensityProj), INTENT(INOUT) :: S
		DOUBLE PRECISION, INTENT(IN) :: Gauge

		DOUBLE PRECISION :: b, factor
		INTEGER :: ta, A

		TYPE (SymHartreeFockFieldProj) :: HF_Gamma, HF_Delta
		TYPE (SymHartreeFockFieldProj) :: ekcm_field, vbb_field, vc_field, vls_field, gdd_field
		TYPE (SymHartreeFockFieldProj) :: field1, field2
		TYPE (SymD3Tensor) :: ek_tensor
		TYPE (SymGenDensityHFProj) :: gendenhf_gamma, gendenhf_delta

		b = Nucleus_get_b(diagonal%consistency%density%nucleus)
		A = Nucleus_get_A(diagonal%consistency%density%nucleus)
		
		factor = 1.0 - (1.0 / A)

		IF ((A .LE. 1) .OR. (A .GE. 300)) THEN
			PRINT *, "Unexpected A value = ", A
			STOP "DiagonalizationMethod::get_MeanField"
		END IF

		! Creating all the necessary fields

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

		! Mean-field - Two-body center of mass correction

		! Protons
		CALL SymKineticEnergy2Body_get_Gamma(field1, diagonal%consistency%vEkCMph, diagonal%consistency%density%field%rho, 0, 0, 0)
		! Neutrons
		CALL SymKineticEnergy2Body_get_Gamma(field1, diagonal%consistency%vEkCMph, diagonal%consistency%density%field%rho, 1, 1, 0)
		
		CALL SymHartreeFockFieldProj_product(ekcm_field, 1.0/A, field1)
		
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

		CALL SymGDDph_update(gdd_field, diagonal%consistency%gDDph, diagonal%consistency%density%field%rho, Gauge)

		! Total Mean-field = Sum of all the preceding terms

		CALL SymHartreeFockFieldProj_add(field1, vls_field, gdd_field)
		CALL SymHartreeFockFieldProj_add(field2, vc_field, field1)
		CALL SymHartreeFockFieldProj_add(field1, vbb_field, field2)
		CALL SymHartreeFockFieldProj_add(field2, ekcm_field, field1)
		CALL SymHartreeFockFieldProj_add_SymD3Tensor(HF_Gamma, ek_tensor, field2)

		! Pairing - Center of Mass term

		IF (HFOnly .EQ. 0) THEN
			CALL SymKineticEnergy2Body_get_Delta(field1, diagonal%consistency%vEkCMpp, diagonal%consistency%density%field%kap, 0, 0, 0)
			CALL SymKineticEnergy2Body_get_Delta(field1, diagonal%consistency%vEkCMpp, diagonal%consistency%density%field%kap, 1, 1, 0)
			CALL SymHartreeFockFieldProj_product(ekcm_field, 1.0/A, field1)
		ELSE
			CALL SymKineticEnergy2Body_get_Delta(field1, diagonal%consistency%vEkCMpp, diagonal%consistency%density%field%kap, 0, 0, 0)
			CALL SymKineticEnergy2Body_get_Delta(field1, diagonal%consistency%vEkCMpp, diagonal%consistency%density%field%kap, 1, 1, 0)
			CALL SymHartreeFockFieldProj_product(ekcm_field, 0.0, field1)
		END IF

		! Pairing - Brink-Boker term

		IF (HFOnly .EQ. 0) THEN
			CALL SymVBBpp_get_Delta(vbb_field, diagonal%consistency%vBBpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymVBBpp_get_Delta(vbb_field, diagonal%consistency%vBBpp, diagonal%consistency%density%field%kap, 1, 1)
		ELSE
			CALL SymVBBpp_get_Delta(field1, diagonal%consistency%vBBpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymHartreeFockFieldProj_product(vbb_field, 0.0, field1)
		END IF

		! Pairing - Coulomb term

		IF (HFOnly .EQ. 0) THEN
			CALL SymVCpp_get_Delta(vc_field, diagonal%consistency%vCpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymVCpp_get_Delta(vc_field, diagonal%consistency%vCpp, diagonal%consistency%density%field%kap, 1, 1)
		ELSE
			CALL SymVCpp_get_Delta(field1, diagonal%consistency%vCpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymHartreeFockFieldProj_product(vc_field, 0.0, field1)
		END IF

		! Pairing - Spin-orbit term
	
		IF (HFOnly .EQ. 0) THEN
                        CALL SymVLSpp_get_Delta(vls_field, diagonal%consistency%vLSpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymVLSpp_get_Delta(vls_field, diagonal%consistency%vLSpp, diagonal%consistency%density%field%kap, 1, 1)
		ELSE
			CALL SymVLSpp_get_Delta(field1, diagonal%consistency%vLSpp, diagonal%consistency%density%field%kap, 0, 0)
			CALL SymHartreeFockFieldProj_product(vls_field, 0.0, field1)
		END IF

		! Total Pairing = Sum of all the preceding terms

		CALL SymHartreeFockFieldProj_add(field1, vc_field, vls_field)
		CALL SymHartreeFockFieldProj_add(field2, ekcm_field, field1)
		CALL SymHartreeFockFieldProj_add(HF_Delta, vbb_field, field2)

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

	!-----------------------------------------------------------------------------------------------!
	!   In this subroutine, we calculate the particle number. This is used to adjust the latter to	!
	!   the actual number of particles of the system at each iteration.				!
	!-----------------------------------------------------------------------------------------------!

	FUNCTION DiagonalizationMethod_operator(diagonal)
		DOUBLE PRECISION DiagonalizationMethod_operator
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: copyR2, copyRho, tmp
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: xI, h, Delta, SH, U, V
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: N20
		
		DOUBLE PRECISION :: np, R2, b2, trace

		INTEGER :: a, d, dd, la, ja, i, jj, k
		INTEGER :: u1, u2
		
		np = 0.0
		R2 = 0.0
		b2 = Nucleus_get_b(diagonal%iterated%nucleus) ** 2
		
		! Block structure: a refers to the angular momentum, but remember we deal here with a 2N x 2N matrix,
		! if N is the size of the single-particle basis.

		DO a = 0, 2*Lmax
		
		
			d = DIM(a)
			dd = d + d
			
			la = L(a)
			ja = J(a)
			
			ALLOCATE(xI(d, d))
			ALLOCATE(h(d, d))
			ALLOCATE(Delta(d, d))
			ALLOCATE(U(d, d))
			ALLOCATE(V(d, d))
			ALLOCATE(tmp(d, d))
			ALLOCATE(copyR2(d, d))
			ALLOCATE(copyRho(d, d))
			
			ALLOCATE(SH(dd, dd))
			
			! xI contains the Fermi level lambda, which was determined before
			! from the conservation of particle number
			
			xI = 0.0
			DO u1 = 1, d
				xI(u1, u1) = diagonal%func%x
			END DO

			h = diagonal%S%rho%rho(diagonal%ta, a)%store - xI

			Delta = diagonal%S%kap%rho(diagonal%ta, a)%store

			! In spherical symmetry, the HFB matrix to diagonalize has the form:
			!
			!		    h      -Delta
			!		   -Delta    -h
			!
			DO u1 = 1, d
				DO u2 = 1, d
					SH(u1    , u2    ) =   h(u1, u2)
					SH(u1    , u2 + d) = - facPair * Delta(u1, u2)
					SH(u1 + d, u2    ) = - facPair * Delta(u1, u2)
					SH(u1 + d, u2 + d) = - h(u1, u2)
				END DO
			END DO

			! Diagonalizing the super-hamiltonian SH gives the quasi-particle energies
			! ordered from the lowest to the highest
			CALL EigenValues(dd, dd, SH, &
				diagonal%QuasiParticleEnergies(diagonal%ta, a)%value, &
				diagonal%UV(diagonal%ta, a)%quantum)

			! Extracting the matrices U and V
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
			
			! The number of particles is the trace of the rho matrix. Don't forget about the
			! degeneracy 2j+1 of a spherical shell (here ja = 2*j to avoid fractions)
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
				trace = trace + tmp(u1, u1)
			END DO
			R2 = R2 + DBLE(ja + 1.0) * trace

			! Freeing memory for the matrix we used
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
		
		! When we dealt with all (l,j)-channels, we have the final number of particles
		! and radius.
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
		
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: S, V, R2Matrix, R2Initial
		DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: RadialWaveFunction
		
		INTEGER :: ndim, a, num, i, k, dim_cur, dim_acc, la, ja, ta, nrot, file_error
		INTEGER :: Lvalue, Jvalue
		INTEGER :: ipoint, loop1, loop2, IndexBra, IndexKet
		INTEGER :: m, n
		
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
		!                 (     h      -Delta  )
		!                 (                    )
		!                 (  -Delta     -h     )
		!
		!  SuperHamiltonian%rho = h (type of SymHartreeFockField)
		!  SuperHamiltonian%kap = Delta (type of SymHartreeFockField)

		CALL DiagonalizationMethod_get_MeanField(diagonal, SuperHamiltonian, Gauge)

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

				CALL Jacobi_real8(S, dim_cur, D(ta)%value, V, nrot)

				E(ta)%quantum = MATMUL(TRANSPOSE(V), MATMUL(SuperHamiltonian%rho%rho(ta, a)%store, V))
				P(ta)%quantum = MATMUL(TRANSPOSE(V), MATMUL(SuperHamiltonian%kap%rho(ta, a)%store, V))
				
				R2Initial = SymD3Tensor_matrix(R2Field, L(a))
				R2Matrix = R2Initial(1:dim_cur,1:dim_cur)
				
				R2(ta)%quantum = MATMUL(TRANSPOSE(V), MATMUL(R2Matrix, V))
				
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

		CALL indexx_real8(dim_acc, SingleParticleE(1)%value, in)
		CALL indexx_real8(dim_acc, SingleParticleE(0)%value, ip)

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
		
		WRITE(*,'()')
		WRITE(*,'(5X,"NEUTRON SINGLE-PARTICLE ENERGIES (CANONICAL BASIS)")')
		WRITE(*,'(5X,"==================================================")')

		ALLOCATE(Count(0:29,0:100))
		
		DO Lvalue = 0, 29
			DO Jvalue = 0, 100
				Count(Lvalue, Jvalue) = 0
			END DO
		END DO
                
		i = 1
		
		DO WHILE (SingleParticleE(1)%value(in(i)) .GE. -100.0 .AND. SingleParticleE(1)%value(in(i)) .LE. 40.0) 
		
			Lvalue = L(an(in(i)))
			Jvalue = J(an(in(i)))
			Count(Lvalue, Jvalue) = Count(Lvalue, Jvalue) + 1
		
			IF (Jvalue .GT. 9) THEN
				WRITE(label,'(2x,I1,A1,I2,"/2")') Count(Lvalue, Jvalue), spectr(Lvalue), Jvalue
			ELSE
				WRITE(label,'(3x,I1,A1,I1,"/2")') Count(Lvalue, Jvalue), spectr(Lvalue), Jvalue
			END IF
			
                        WRITE(*,'(I4,")",2X,A2,I2,"/2",3x,F10.5,F10.3)') &
					i, spectr(Lvalue), Jvalue,Occup(1)%value(in(i)), SingleParticleE(1)%value(in(i))
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

		WRITE(*,'()')
		WRITE(*,'(5X,"PROTON SINGLE-PARTICLE ENERGIES (CANONICAL BASIS)")')
		WRITE(*,'(5X,"=================================================")')

		DO Lvalue = 0, 29
			DO Jvalue = 0, 100
				Count(Lvalue, Jvalue) = 0
			END DO
		END DO
                
		i = 1
		
		DO WHILE (SingleParticleE(0)%value(ip(i)) .GE. -100.0 .AND. SingleParticleE(0)%value(ip(i)) .LE. 40.0) 
		
			Lvalue = L(ap(in(i)))
			Jvalue = J(ap(in(i)))
			
			Count(Lvalue, Jvalue) = Count(Lvalue, Jvalue) + 1
		
			IF (Jvalue .GT. 9) THEN
				WRITE(label,'(2x,I1,A1,I2,"/2")') Count(Lvalue, Jvalue), spectr(Lvalue), Jvalue
			ELSE
				WRITE(label,'(3x,I1,A1,I1,"/2")') Count(Lvalue, Jvalue), spectr(Lvalue), Jvalue
			END IF
                        WRITE(*,'(I4,")",2X,A2,I2,"/2",3X,F10.8,F10.3)') &
					i, spectr(Lvalue), Jvalue, Occup(0)%value(ip(i)), SingleParticleE(0)%value(ip(i))
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

		TYPE (SymD3Tensor) :: ek_tensor
		
		INTEGER, DIMENSION(:), ALLOCATABLE :: nIndx, pIndx
		
                INTEGER, DIMENSION(:, :), ALLOCATABLE :: Momentum

                DOUBLE PRECISION, DIMENSION(1:Nmax, 1:Nmax) :: Matrix
                DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: nVV, pVV
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: U, V, QPE, Gamma, Delta
		DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: QPWaveFunctionLower, QPWaveFunctionUpper
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: PairingField, MeanField

                INTEGER :: NumberOfStates, a, na, top, Bra, Ket
                INTEGER :: num, la, ja, d, s1, ta, sa, an, ap, jn, jp
                INTEGER :: IndexBasis, IndexKet, nradial, ipoint, IndxQP
		INTEGER :: u1, u2, s, loop, file_error
		INTEGER :: m, n, i, state, Anucl, Nnucl, Znucl
                
		INTEGER, PARAMETER :: file_unit_1 = 16, file_unit_2 = 17
		
		DOUBLE PRECISION :: sumvv, test, factor
                
		CHARACTER :: cn, cp
                
		CHARACTER, DIMENSION(0:19) :: spectr
                
		DATA spectr / "s", "p", "d", "f", "g", "h", "i", "j", "k", "l", &
				"m", "n", "o", "p", "q", "r", "s", "t", "u", "v" /
		
		NumberOfStates = 0
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
		IF (.NOT. ALLOCATED(DensityPairing)) ALLOCATE(DensityPairing(0:1, 0:Npoint))
		IF (.NOT. ALLOCATED(PairingField)) ALLOCATE(PairingField(0:1, 0:Npoint))
                IF (.NOT. ALLOCATED(DerivDens)) ALLOCATE(DerivDens(0:1, 0:Npoint))
		IF (.NOT. ALLOCATED(MeanField)) ALLOCATE(MeanField(0:1, 0:Npoint))

  		Nnucl = Nucleus_get_N(diagonal%consistency%density%nucleus)
 		Znucl = Nucleus_get_Z(diagonal%consistency%density%nucleus)
			
		Anucl = Nnucl + Znucl
		factor = 1.0 - (1.0 / Anucl)

		! Initialize all densities to 0
		DensityRadial = 0.0
		DensityPairing = 0.0
		PairingField = 0.0
		DerivDens = 0.0
		MeanField = 0.0
		
		CALL SymD3Tensor_new(ek_tensor)
		
		DO ta = 0, 1
			
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
				
					QPE(ta, num) = diagonal%QuasiParticleEnergies(ta, a)%value(s1)
                                
					Momentum(ta, num) = a
                                
					IF (ta .EQ. 1) nIndx(num) = num
					IF (ta .EQ. 0) pIndx(num) = num
                                
					sumvv = 0.0
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
						
					ALLOCATE(Delta(d, d))
					ALLOCATE(Gamma(d, d))
					
					! Define the pairing field
					Delta = diagonal%S%kap%rho(ta, a)%store
					
					! Define the mean-field (we've got to substract the kinetic energy here to get the 
					! mean-field potential)
               				CALL SymD3Tensor_product(ek_tensor, factor, EkField)
					
					Matrix = SymD3Tensor_matrix(ek_tensor, la)
					
					Gamma = diagonal%S%rho%rho(ta, a)%store - Matrix(1:d, 1:d)
					
					DO ipoint = 1, Npoint
					
						DO i = 1, d
						
							Bra = IndexVecNL(i, la)
							
							DO m = 1, d
							
								Ket = IndexVecNL(m, la)
								
								PairingField(ta, ipoint) = PairingField(ta, ipoint) &
						+ (-1)**la * (ja + 1.0) * Delta(i, m) * WaveFun(ipoint, Bra)/RadMesh(ipoint) &
										      * WaveFun(ipoint, Ket)/RadMesh(ipoint)

								MeanField(ta, ipoint) = MeanField(ta, ipoint) &
								+ (ja + 1.0) * Gamma(i, m) * WaveFun(ipoint, Bra)/RadMesh(ipoint) &
											   * WaveFun(ipoint, Ket)/RadMesh(ipoint)

							END DO
							
						END DO
						
					END DO
					
					DEALLOCATE(Delta)
					DEALLOCATE(Gamma)
						
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
												
							                DensityPairing(ta, ipoint) = DensityPairing(ta, ipoint) &
                                                                                     - (-1)**la * (ja + 1.0) * V(m, i) * U(n, i) &
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
					DensityPairing(ta, ipoint) = DensityPairing(ta, ipoint) / (4.0*PI*RadMesh(ipoint)**2)
				END DO

				! Extrapolate to find rho(r=0)

         			DensityRadial(ta, 0) = 3.0*(DensityRadial(ta, 1) - DensityRadial(ta, 2)) + DensityRadial(ta, 3)
         			DensityPairing(ta, 0) = 3.0*(DensityPairing(ta, 1) - DensityPairing(ta, 2)) + DensityPairing(ta, 3)
         			PairingField(ta, 0) = 3.0*(PairingField(ta, 1) - PairingField(ta, 2)) + PairingField(ta, 3)
         			MeanField(ta, 0) = 3.0*(MeanField(ta, 1) - MeanField(ta, 2)) + MeanField(ta, 3)

                                DO ipoint=1, Npoint-1
                                        DerivDens(ta, ipoint) = 0.5*(DensityRadial(ta, ipoint+1) - DensityRadial(ta, ipoint-1)) &
                                                              /(RadMesh(ipoint+1) - RadMesh(ipoint))
                                END DO

                                DerivDens(ta, 0) = 3.0*(DerivDens(ta, 1) - DerivDens(ta, 2)) + DerivDens(ta, 3)
                                DerivDens(ta, Npoint) = 3.0*(DerivDens(ta, Npoint-1) - DerivDens(ta, Npoint-2)) + DerivDens(ta, Npoint-3)

			END IF
		
		END DO ! end of loop over ta
					
		! Ordering of the quasi-particle occupation numbers (by decreasing order!!)
		
		CALL indexx_real8(NumberOfStates, 1.0-nVV, nIndx)
		CALL indexx_real8(NumberOfStates, 1.0-pVV, pIndx)

		WRITE(*,'()')
		WRITE(*,'(5X,"QUASI-PARTICLE-PARTICLE ENERGIES")')
		WRITE(*,'(5X,"================================")')

		WRITE(*,'(A14,A10,A12,A5,A10)') "EQP NEUT", "V2", "", "V2", "EQP PROT"
		
		top = MIN(100, NumberOfStates)
		
		DO num = 1, top

			an = Momentum(1, nIndx(num))
			ap = Momentum(0, pIndx(num))
			
			cn = spectr(L(an))                       
			cp = spectr(L(ap))
			
			jn = J(an)
			jp = J(ap)

                        WRITE(*,'(i4,")",i5,F10.3,F10.5,3X,A1,I2,"/2",A5,I2,"/2",3x,F10.5,F10.3)') &
				num, nIndx(num),QPE(1, nIndx(num)), nVV(nIndx(num)), cn, jn, cp, jp, pVV(pIndx(num)), QPE(0, pIndx(num))
		END DO
		
		! Quasi-particle energies (proton and neutron)

		OPEN(file_unit_1, FILE='data/HF_qp_n.dat', ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			WRITE(*,'("Impossible to open the file data/HF_qp_n.dat")') 
			STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
		END IF
		
		DO num = 1, top
			an = Momentum(1, nIndx(num))
			cn = spectr(L(an))                       
			jn = J(an)
			WRITE(file_unit_1,'(I3,2f25.15,1X,A1,i5)') num, QPE(1, nIndx(num)), nVV(nIndx(num)), cn, jn
		END DO
		
		CLOSE(file_unit_1)
		
		OPEN(file_unit_1, FILE='data/HF_qp_p.dat', ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			WRITE(*,'("Impossible to open the file data/HF_qp_p.dat")') 
			STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
		END IF
		
		DO num = 1, top
			ap = Momentum(1, pIndx(num))
			cp = spectr(L(ap))                       
			jp = J(ap)
			WRITE(file_unit_1,'(I3,2f25.15,1X,A1,i5)') num, QPE(0, pIndx(num)), pVV(pIndx(num)), cp, jp
		END DO
		
		CLOSE(file_unit_1)
		
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
		
			! Total pairing density (proton and neutron)

			OPEN(file_unit_1, FILE='data/DensityPair.dat', ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				WRITE(*,'("Impossible to open the file data/DensityPair.dat")') 
				STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
			END IF
		
			DO ipoint = 0, Npoint
				WRITE(file_unit_1,'(3f20.16)') RadMesh(ipoint),DensityPairing(0, ipoint),DensityPairing(1, ipoint)
			END DO
		
			CLOSE(file_unit_1)

                	! Logarithmic derivative of the density (proton and neutron)

			OPEN(file_unit_1, FILE='data/DensityLogDeriv.dat', ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				WRITE(*,'("Impossible to open the file data/DensityLogDeriv.dat")') 
				STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
			END IF
		
			DO ipoint = 0, Npoint
				WRITE(file_unit_1,'(3f25.16)') RadMesh(ipoint),DerivDens(0, ipoint)/(DensityRadial(0, ipoint) + 1.e-20), &
                                                                               DerivDens(1, ipoint)/(DensityRadial(1, ipoint) + 1.e-20)
			END DO
		
			CLOSE(file_unit_1)
		
			! Pairing and mean-field potentials. They are defined as:
			!
			!   Delta(r) = \sum delta_{ac} phi_{a}(r)*phi_{c}(r)
			!   Gamma(r) = \sum gamma_{ac} phi_{a}(r)*phi_{c}(r)

			OPEN(file_unit_1, FILE='data/PairingField.dat', ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				WRITE(*,'("Impossible to open the file data/PairingField.dat")') 
				STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
			END IF
		
			DO ipoint = 0, Npoint
				WRITE(file_unit_1,'(5f20.14)') RadMesh(ipoint), &
					PairingField(0, ipoint)/(2.0*PI), PairingField(1, ipoint)/(2.0*PI),&
					MeanField(0, ipoint)/(2.0*PI), MeanField(1, ipoint)/(2.0*PI)
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
		TYPE (DiagonalizationMethod), INTENT(OUT) :: diagonal
		
		TYPE (SymGenDensityGaugeProj), DIMENSION(:), POINTER :: genden_gauge_real_n, genden_gauge_imag_n, genden_gauge_real_p, genden_gauge_imag_p
		TYPE (SymGenDensityProj) :: genden_proj_real, genden_proj_imag, SuperHamiltonian
		TYPE (SymDensityProj) :: density_real, density01_real, density10_real
		TYPE (SymDensityProj) :: density_imag, density01_imag, density10_imag
		TYPE (SymHartreeFockFieldProj) :: HF_Gamma_Diff, HF_Gamma_Same, HF_Delta
		TYPE (ProjectionCoeffs) :: CoeffsXY
		TYPE (SelfConsistencyMethodProj) :: consistency
                TYPE (SymGDDph) :: gDDph

		COMPLEX, DIMENSION(:, :), ALLOCATABLE :: copyR2, copyRho, tmp
		
		COMPLEX :: ProjectedEkCM, ProjectedEk
		COMPLEX :: ProjectedVBB_local, ProjectedVBB_exch, ProjectedVC_local, ProjectedVC_exch, ProjectedVLS, ProjectedDD
		COMPLEX :: ProjectedVBB_Pair, ProjectedVC_Pair, ProjectedVLS_Pair, ProjectedEkCM_Pair
		COMPLEX :: trace
		COMPLEX :: eHF, pairing, pairing_BB, pairing_Coulomb, pairing_LS, Pairing_EkCM

		DOUBLE PRECISION, DIMENSION(0:1) :: P

		DOUBLE PRECISION :: A, Gauge, factor, FactorNeut, FactorProt
		DOUBLE PRECISION :: np, R2, b2, NReal, ZReal
		
		INTEGER :: Gauge_n, Gauge_p, Z, N, IndexGauge
		INTEGER :: ProjectionOn = 1
		INTEGER :: ta, aa, d, dd, la, ja, i, jj, k
		INTEGER :: u1, u2, NGauge_n, NGauge_p
		INTEGER :: file_error, file_unit_1 = 21
                INTEGER :: IndexBasis, IndexKet, ipoint, mm, nn                
			
		A = Nucleus_get_A(diagonal%consistency%density%nucleus)
		
		IF (A .LE. 1) STOP "Abortado"
		
		Z = diagonal%consistency%density%nucleus%np(0)
		N = diagonal%consistency%density%nucleus%np(1)
		
		FactorProt = 1.0
		FactorNeut = 1.0
		
		NGauge_n = NGauge
		NGauge_p = NGauge

                IF (MOD(Z,2) .EQ. 1) NGauge_p = 1
                IF (MOD(N,2) .EQ. 1) NGauge_n = 1
		
		ALLOCATE(CoeffsXY%x_l(0:1))
		ALLOCATE(CoeffsXY%y_l(0:1))
		
		! Create new Fomenko coefficients
		
		CALL ProjectionCoeffs_new(CoeffsXY, NGauge_p, 0)
		CALL ProjectionCoeffs_new(CoeffsXY, NGauge_n, 1)
		
		! Calculate phi-dependent densities rho(phi), kappa10(phi), kappa01(phi) according
		! to formulas (9), (10) and (11) of the Anguiano et al, Nucl. Phys. A696 (2001) 467-493
		! All densities are stored into genden_gauge.
		
		ALLOCATE(genden_gauge_real_n(NGauge_n))
		ALLOCATE(genden_gauge_imag_n(NGauge_n))
		
		DO IndexGauge = 1, NGauge_n
		
			CALL SymGenDensityHFProj_new(genden_gauge_real_n(IndexGauge)%rho)
			CALL SymGenDensityHFProj_new(genden_gauge_real_n(IndexGauge)%kap01)
			CALL SymGenDensityHFProj_new(genden_gauge_real_n(IndexGauge)%kap10)
			
			CALL SymGenDensityHFProj_new(genden_gauge_imag_n(IndexGauge)%rho)
			CALL SymGenDensityHFProj_new(genden_gauge_imag_n(IndexGauge)%kap01)
			CALL SymGenDensityHFProj_new(genden_gauge_imag_n(IndexGauge)%kap10)
			
		END DO
		
		ALLOCATE(genden_gauge_real_p(NGauge_p))
		ALLOCATE(genden_gauge_imag_p(NGauge_p))
		
		DO IndexGauge = 1, NGauge_p
		
			CALL SymGenDensityHFProj_new(genden_gauge_real_p(IndexGauge)%rho)
			CALL SymGenDensityHFProj_new(genden_gauge_real_p(IndexGauge)%kap01)
			CALL SymGenDensityHFProj_new(genden_gauge_real_p(IndexGauge)%kap10)
			
			CALL SymGenDensityHFProj_new(genden_gauge_imag_p(IndexGauge)%rho)
			CALL SymGenDensityHFProj_new(genden_gauge_imag_p(IndexGauge)%kap01)
			CALL SymGenDensityHFProj_new(genden_gauge_imag_p(IndexGauge)%kap10)
			
		END DO
	
		! Calculate the gauge-dependent density
		
                CALL DiagonalizationMethod_set_Isospin(diagonal, 0)
		CALL SymGenDensityProj_make_DensityGauge(diagonal%UV, genden_gauge_real_p, genden_gauge_imag_p, diagonal%consistency%density%nucleus, &
																FactorProt, NGauge_p, 0)

                CALL DiagonalizationMethod_set_Isospin(diagonal, 1)
		CALL SymGenDensityProj_make_DensityGauge(diagonal%UV, genden_gauge_real_n, genden_gauge_imag_n, diagonal%consistency%density%nucleus, &
																FactorNeut, NGauge_n, 1)
		
		! Calculate the projected density and the Fomenko coefficients
		
		CALL SymGenDensityProj_new_Nucleus(genden_proj_real, N, Z)
		CALL SymGenDensityProj_new_Nucleus(genden_proj_imag, N, Z)
		
		CALL SymGenDensityProj_make_DensityProj(genden_proj_real, genden_proj_imag, CoeffsXY, genden_gauge_real_p, genden_gauge_imag_p, diagonal%UV, &
												diagonal%consistency%density%nucleus, FactorProt, NGauge_p, 0)
		CALL SymGenDensityProj_make_DensityProj(genden_proj_real, genden_proj_imag, CoeffsXY, genden_gauge_real_n, genden_gauge_imag_n, diagonal%UV, &
												diagonal%consistency%density%nucleus, FactorNeut, NGauge_n, 1)
		
		CALL SymHartreeFockFieldProj_new(HF_Gamma_Same)
		
		! Convert projected density into the right type, which is SymDensityProj
		
		CALL SymDensity_new_GenDensityProj(density_real, genden_proj_real)
		CALL SymDensity_new_GenDensityProj(density_imag, genden_proj_imag)
		
		!---------------------------------------------------------------------------------------!
		! 		Terms depending on the projected density rhoP				!
		!---------------------------------------------------------------------------------------!		
		
                factor = 1.0 - (1.0 / A)
		
		ProjectedEk = factor * (EkField * density_real%field%rho%p(0) + EkField * density_real%field%rho%p(1)) &
			    + factor * (EkField * density_imag%field%rho%p(0) + EkField * density_imag%field%rho%p(1))*CMPLX(0, 1)
		
		! Two-body kinetic energy
		CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma_Same, diagonal%consistency%vEkCMph, density_real%field%rho, 0, 1, ProjectionOn)
		CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma_Same, diagonal%consistency%vEkCMph, density_real%field%rho, 1, 0, ProjectionOn)
		
		P(0) = SymHartreeFockFieldProj_product_iso(density_real%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density_real%field%rho, HF_Gamma_Same, 1)
		
		ProjectedEkCM = (0.5 / A) * (P(0) + P(1))*CMPLX(1, 0)
		
		CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma_Same, diagonal%consistency%vEkCMph, density_imag%field%rho, 0, 1, ProjectionOn)
		CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma_Same, diagonal%consistency%vEkCMph, density_imag%field%rho, 1, 0, ProjectionOn)
		
		P(0) = SymHartreeFockFieldProj_product_iso(density_imag%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density_imag%field%rho, HF_Gamma_Same, 1)
		
		ProjectedEkCM = ProjectedEkCM + (0.5 / A) * (P(0) + P(1))*CMPLX(0, 1)
		
		! Brink-Boeker direct terms
		CALL SymVBBph_get_LocalGamma(HF_Gamma_Same, diagonal%consistency%vBBph, density_real%field%rho, 0, 1, ProjectionOn)
		CALL SymVBBph_get_LocalGamma(HF_Gamma_Same, diagonal%consistency%vBBph, density_real%field%rho, 1, 0, ProjectionOn)
		
		P(0) = SymHartreeFockFieldProj_product_iso(density_real%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density_real%field%rho, HF_Gamma_Same, 1)
		
		ProjectedVBB_local = 0.5 * (P(0) + P(1))*CMPLX(1, 0)
		
		CALL SymVBBph_get_LocalGamma(HF_Gamma_Same, diagonal%consistency%vBBph, density_imag%field%rho, 0, 1, ProjectionOn)
		CALL SymVBBph_get_LocalGamma(HF_Gamma_Same, diagonal%consistency%vBBph, density_imag%field%rho, 1, 0, ProjectionOn)
		
		P(0) = SymHartreeFockFieldProj_product_iso(density_imag%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density_imag%field%rho, HF_Gamma_Same, 1)
		
		ProjectedVBB_local = ProjectedVBB_local + 0.5 * (P(0) + P(1))*CMPLX(0, 1)
		
		! Brink-Boeker exchange terms
		CALL SymVBBph_get_ExchangeGamma(HF_Gamma_Same, diagonal%consistency%vBBph, density_real%field%rho, 0, 1, ProjectionOn)
		CALL SymVBBph_get_ExchangeGamma(HF_Gamma_Same, diagonal%consistency%vBBph, density_real%field%rho, 1, 0, ProjectionOn)
		
		P(0) = SymHartreeFockFieldProj_product_iso(density_real%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density_real%field%rho, HF_Gamma_Same, 1)
		
		ProjectedVBB_exch = 0.5 * (P(0) + P(1))*CMPLX(1, 0)
		
		CALL SymVBBph_get_ExchangeGamma(HF_Gamma_Same, diagonal%consistency%vBBph, density_imag%field%rho, 0, 1, ProjectionOn)
		CALL SymVBBph_get_ExchangeGamma(HF_Gamma_Same, diagonal%consistency%vBBph, density_imag%field%rho, 1, 0, ProjectionOn)
		
		P(0) = SymHartreeFockFieldProj_product_iso(density_imag%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density_imag%field%rho, HF_Gamma_Same, 1)
		
		ProjectedVBB_exch = ProjectedVBB_exch + 0.5 * (P(0) + P(1))*CMPLX(0, 1)
		
		! Coulomb direct terms
		CALL SymVCph_get_LocalGamma(HF_Gamma_Same, diagonal%consistency%vCph, density_real%field%rho, 0, 1, ProjectionOn)
		CALL SymVCph_get_LocalGamma(HF_Gamma_Same, diagonal%consistency%vCph, density_real%field%rho, 1, 0, ProjectionOn)

		P(0) = SymHartreeFockFieldProj_product_iso(density_real%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density_real%field%rho, HF_Gamma_Same, 1)
		
		ProjectedVC_local = 0.5 * (P(0) + P(1))*CMPLX(1, 0)
		
		CALL SymVCph_get_LocalGamma(HF_Gamma_Same, diagonal%consistency%vCph, density_imag%field%rho, 0, 1, ProjectionOn)
		CALL SymVCph_get_LocalGamma(HF_Gamma_Same, diagonal%consistency%vCph, density_imag%field%rho, 1, 0, ProjectionOn)

		P(0) = SymHartreeFockFieldProj_product_iso(density_imag%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density_imag%field%rho, HF_Gamma_Same, 1)
		
		ProjectedVC_local = ProjectedVC_local + 0.5 * (P(0) + P(1))*CMPLX(0, 1)
		
		! Coulomb exchange terms
		CALL SymVCph_get_ExchangeGamma(HF_Gamma_Same, diagonal%consistency%vCph, density_real%field%rho, 0, 1, ProjectionOn)
		CALL SymVCph_get_ExchangeGamma(HF_Gamma_Same, diagonal%consistency%vCph, density_real%field%rho, 1, 0, ProjectionOn)

		P(0) = SymHartreeFockFieldProj_product_iso(density_real%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density_real%field%rho, HF_Gamma_Same, 1)
		
		ProjectedVC_exch = 0.5 * (P(0) + P(1))*CMPLX(1, 0)
		
		CALL SymVCph_get_ExchangeGamma(HF_Gamma_Same, diagonal%consistency%vCph, density_imag%field%rho, 0, 1, ProjectionOn)
		CALL SymVCph_get_ExchangeGamma(HF_Gamma_Same, diagonal%consistency%vCph, density_imag%field%rho, 1, 0, ProjectionOn)

		P(0) = SymHartreeFockFieldProj_product_iso(density_imag%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density_imag%field%rho, HF_Gamma_Same, 1)
		
		ProjectedVC_exch = ProjectedVC_exch + 0.5 * (P(0) + P(1))*CMPLX(0, 1)
		
		! Spin-orbit term
		CALL SymVLSph_get_Gamma(HF_Gamma_Same, diagonal%consistency%vLSph, density_real%field%rho, 0, 1, ProjectionOn)
		CALL SymVLSph_get_Gamma(HF_Gamma_Same, diagonal%consistency%vLSph, density_real%field%rho, 1, 0, ProjectionOn)

		P(0) = SymHartreeFockFieldProj_product_iso(density_real%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density_real%field%rho, HF_Gamma_Same, 1)
		
		ProjectedVLS = 0.5 * (P(0) + P(1))*CMPLX(1, 0)
		
		CALL SymVLSph_get_Gamma(HF_Gamma_Same, diagonal%consistency%vLSph, density_imag%field%rho, 0, 1, ProjectionOn)
		CALL SymVLSph_get_Gamma(HF_Gamma_Same, diagonal%consistency%vLSph, density_imag%field%rho, 1, 0, ProjectionOn)

		P(0) = SymHartreeFockFieldProj_product_iso(density_imag%field%rho, HF_Gamma_Same, 0)
		P(1) = SymHartreeFockFieldProj_product_iso(density_imag%field%rho, HF_Gamma_Same, 1)
		
		ProjectedVLS = ProjectedVLS + 0.5 * (P(0) + P(1))*CMPLX(0, 1)
		
		! Density-dependent projected energy
		
		Gauge = 0.0
		
                CALL SymGDDph_new(gDDph)
		
		IF (Basis .EQ. 1) THEN
			CALL SymGDDph_make_DD(gDDph, density_real%field%rho, Gauge)
		ELSE
			CALL Make_DenGenFun(gDDph, density_real%field%rho, Gauge)
		END IF
		
		ProjectedDD = SymGDDph_get_edd(gDDph)*CMPLX(1, 0)
		
		!IF (Basis .EQ. 1) THEN
		!	CALL SymGDDph_make_DD(gDDph, density_imag%field%rho, Gauge)
		!ELSE
		!	CALL Make_DenGenFun(gDDph, density_imag%field%rho, Gauge)
		!END IF
		
		!ProjectedDD = ProjectedDD + SymGDDph_get_edd(gDDph)*CMPLX(0, 1)
		
		!WRITE(*,'(/,14X,"----------------------------------------")')
		!WRITE(*,'(14X,"| PROJECTION AFTER VARIATION - RESULTS |")')
		!WRITE(*,'(14X,"----------------------------------------",/)')
		
		!WRITE(*,'("                                 Energia cinetica:",2F20.10)') ProjectedEk
		!WRITE(*,'("             Energia cinetica del centro de masas:",2F20.10)') ProjectedEkCM
		!WRITE(*,'("                    Energia local de Brink-Booker:",2F20.10)') ProjectedVBB_local
		!WRITE(*,'("                         Energia local de Coulomb:",2F20.10)') ProjectedVC_local
		!WRITE(*,'("           Energia de intercambio de Brink-Booker:",2F20.10)') ProjectedVBB_exch
		!WRITE(*,'("                Energia de intercambio de Coulomb:",2F20.10)') ProjectedVC_exch
		!WRITE(*,'("               Energia dependiente de la densidad:",2F20.10)') ProjectedDD
		!WRITE(*,'("                           Energia de spin-orbita:",2F20.10)') ProjectedVLS
		!WRITE(*,'(70("-"))')

		eHF = ProjectedEk + ProjectedEkCM &
			+ ProjectedVBB_local + ProjectedVBB_exch &
			+ ProjectedVC_local +  ProjectedVC_exch&
			+ ProjectedDD + ProjectedVLS
			
		!WRITE(*,'("                    Energia total de Brink-Booker:",2F20.10)') ProjectedVBB_local + ProjectedVBB_exch
		!WRITE(*,'("                         Energia total de Coulomb:",2F20.10)') ProjectedVC_local + ProjectedVC_exch
		!WRITE(*,'("                     Energia total (Hartree-Fock):",2F20.10)') eHF
                !WRITE(*,'(70("-"))')

		pairing_BB = ProjectedVBB_Pair
		pairing_Coulomb = ProjectedVC_Pair
		pairing_LS = ProjectedVLS_Pair
		
		!WRITE(*,'("                     Apareamiento de Brink-Booker:",2F20.10)') pairing_BB
		!WRITE(*,'("                          Apareamiento de Coulomb:",2F20.10)') pairing_Coulomb
		!WRITE(*,'("                      Apareamiento de spin-orbita:",2F20.10)') pairing_LS
		!WRITE(*,'(70("-"))')

		pairing = pairing_BB + pairing_Coulomb + pairing_LS
		
		!WRITE(*,'("                               Apareamiento total:",2F20.10)') pairing
		!WRITE(*,'(70("-"))')

		!WRITE(*,'("                                            TOTAL:",2F20.10)') eHF + pairing

								    
		!---------------------------------------------------------------------------------------!
		! 	Terms depending on the Gauge densities rho(phi), kap10(phi), kap01(phi)		!
		!---------------------------------------------------------------------------------------!		
		
		ProjectedVBB_Pair  = CMPLX(0,0)
		ProjectedVC_Pair   = CMPLX(0,0)
		ProjectedVLS_Pair  = CMPLX(0,0)
                ProjectedEkCM_Pair = CMPLX(0,0)
				
		DO IndexGauge = 1, NGauge_p
		
			! Convert projected density into the right type, which is SymDensityProj
		
			CALL SymDensity_new_GenDensityProj10(density10_real, genden_gauge_real_p(IndexGauge))
			CALL SymDensity_new_GenDensityProj10(density10_imag, genden_gauge_imag_p(IndexGauge))
			
			CALL SymHartreeFockFieldProj_new(HF_Gamma_Diff)
			CALL SymHartreeFockFieldProj_new(HF_Delta)
		
			!
			! Mean-field Channel
			!
			
			! Two-body kinetic energy (center of mass motion)
			CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vEkCMph, density10_real%field%rho, 0, 0, ProjectionOn)
		
			P(0) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedEkCM = ProjectedEkCM + (0.5 / A) * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
		
			CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vEkCMph, density10_imag%field%rho, 0, 0, ProjectionOn)
		
			P(0) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedEkCM = ProjectedEkCM - (0.5 / A) * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
		
			CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vEkCMph, density10_real%field%rho, 0, 0, ProjectionOn)
		
			P(0) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedEkCM = ProjectedEkCM + (0.5 / A) * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
		
			CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vEkCMph, density10_imag%field%rho, 0, 0, ProjectionOn)
		
			P(0) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedEkCM = ProjectedEkCM + (0.5 / A) * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
		
			! Brink-Boeker direct term
			CALL SymVBBph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10_real%field%rho, 0, 0, ProjectionOn)
		
			P(0) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVBB_local = ProjectedVBB_local + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
		
			CALL SymVBBph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10_imag%field%rho, 0, 0, ProjectionOn)
		
			P(0) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVBB_local = ProjectedVBB_local - 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
		
			CALL SymVBBph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10_real%field%rho, 0, 0, ProjectionOn)
		
			P(0) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVBB_local = ProjectedVBB_local + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
		
			CALL SymVBBph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10_imag%field%rho, 0, 0, ProjectionOn)
		
			P(0) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVBB_local = ProjectedVBB_local + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
		
			! Brink-Boeker exchange term
			CALL SymVBBph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10_real%field%rho, 0, 0, ProjectionOn)
		
			P(0) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVBB_exch = ProjectedVBB_exch + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
			
			CALL SymVBBph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10_imag%field%rho, 0, 0, ProjectionOn)
		
			P(0) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVBB_exch = ProjectedVBB_exch - 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
			
			CALL SymVBBph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10_real%field%rho, 0, 0, ProjectionOn)
		
			P(0) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVBB_exch = ProjectedVBB_exch + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
			
			CALL SymVBBph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10_imag%field%rho, 0, 0, ProjectionOn)
		
			P(0) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVBB_exch = ProjectedVBB_exch + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
			
			! Coulomb direct term
			CALL SymVCph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10_real%field%rho, 0, 0, ProjectionOn)

			P(0) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVC_local = ProjectedVC_local + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
			
			CALL SymVCph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10_imag%field%rho, 0, 0, ProjectionOn)

			P(0) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVC_local = ProjectedVC_local - 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
			
			CALL SymVCph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10_real%field%rho, 0, 0, ProjectionOn)

			P(0) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVC_local = ProjectedVC_local + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
			
			CALL SymVCph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10_imag%field%rho, 0, 0, ProjectionOn)

			P(0) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVC_local = ProjectedVC_local + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
			
			! Coulomb exchange term
			CALL SymVCph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10_real%field%rho, 0, 0, ProjectionOn)

			P(0) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVC_exch = ProjectedVC_exch + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
		
			CALL SymVCph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10_imag%field%rho, 0, 0, ProjectionOn)

			P(0) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVC_exch = ProjectedVC_exch - 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
		
			CALL SymVCph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10_real%field%rho, 0, 0, ProjectionOn)

			P(0) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVC_exch = ProjectedVC_exch + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
		
			CALL SymVCph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10_imag%field%rho, 0, 0, ProjectionOn)

			P(0) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVC_exch = ProjectedVC_exch + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
		
			! Spin-orbit term
			CALL SymVLSph_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vLSph, density10_real%field%rho, 0, 0, ProjectionOn)

			P(0) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVLS = ProjectedVLS + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
		
			CALL SymVLSph_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vLSph, density10_imag%field%rho, 0, 0, ProjectionOn)

			P(0) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVLS = ProjectedVLS - 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
		
			CALL SymVLSph_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vLSph, density10_real%field%rho, 0, 0, ProjectionOn)

			P(0) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVLS = ProjectedVLS + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
		
			CALL SymVLSph_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vLSph, density10_imag%field%rho, 0, 0, ProjectionOn)

			P(0) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 0)
		
			ProjectedVLS = ProjectedVLS + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
		
			!
			! Pairing Channel
			!
			
			CALL SymDensity_new_GenDensityProj01(density01_real, genden_gauge_real_p(IndexGauge))
			CALL SymDensity_new_GenDensityProj01(density01_imag, genden_gauge_imag_p(IndexGauge))
			
			!density01 = diagonal%consistency%density
			
                        ! Center of Mass Pairing
                        CALL SymKineticEnergy2Body_get_Delta(HF_Delta, diagonal%consistency%vEkCMpp, density10_real%field%kap, 0, 0, ProjectionOn)

                        P(0) =  SymHartreeFockFieldProj_product_iso(density01_real%field%kap, HF_Delta, 0)

                        ProjectedEkCM_Pair = ProjectedEkCM_Pair - ( 0.5/A ) * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)

                        CALL SymKineticEnergy2Body_get_Delta(HF_Delta, diagonal%consistency%vEkCMpp, density10_imag%field%kap, 0, 0, ProjectionOn)

                        P(0) =  SymHartreeFockFieldProj_product_iso(density01_imag%field%kap, HF_Delta, 0)

                        ProjectedEkCM_Pair = ProjectedEkCM_Pair + ( 0.5/A ) * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)

                        CALL SymKineticEnergy2Body_get_Delta(HF_Delta, diagonal%consistency%vEkCMpp, density10_real%field%kap, 0, 0, ProjectionOn)

                        P(0) =  SymHartreeFockFieldProj_product_iso(density01_imag%field%kap, HF_Delta, 0)

                        ProjectedEkCM_Pair = ProjectedEkCM_Pair - ( 0.5/A ) * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))

                        CALL SymKineticEnergy2Body_get_Delta(HF_Delta, diagonal%consistency%vEkCMpp, density10_imag%field%kap, 0, 0, ProjectionOn)

                        P(0) =  SymHartreeFockFieldProj_product_iso(density01_real%field%kap, HF_Delta, 0)

                        ProjectedEkCM_Pair = ProjectedEkCM_Pair - ( 0.5/A ) * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))

			! Brink-Boeker Pairing
			CALL SymVBBpp_get_Delta(HF_Delta, diagonal%consistency%vBBpp, density10_real%field%kap, 0, 0)
			
			P(0) = SymHartreeFockFieldProj_product_iso(density01_real%field%kap, HF_Delta, 0)
		
			ProjectedVBB_Pair = ProjectedVBB_Pair - 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
		
			CALL SymVBBpp_get_Delta(HF_Delta, diagonal%consistency%vBBpp, density10_imag%field%kap, 0, 0)
			
			P(0) = SymHartreeFockFieldProj_product_iso(density01_imag%field%kap, HF_Delta, 0)
		
			ProjectedVBB_Pair = ProjectedVBB_Pair + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
		
			CALL SymVBBpp_get_Delta(HF_Delta, diagonal%consistency%vBBpp, density10_real%field%kap, 0, 0)
			
			P(0) = SymHartreeFockFieldProj_product_iso(density01_imag%field%kap, HF_Delta, 0)
		
			ProjectedVBB_Pair = ProjectedVBB_Pair - 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
		
			CALL SymVBBpp_get_Delta(HF_Delta, diagonal%consistency%vBBpp, density10_imag%field%kap, 0, 0)
			
			P(0) = SymHartreeFockFieldProj_product_iso(density01_real%field%kap, HF_Delta, 0)
		
			ProjectedVBB_Pair = ProjectedVBB_Pair - 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
		
			! Coulomb Pairing
			CALL SymVCpp_get_Delta(HF_Delta, diagonal%consistency%vCpp, density10_real%field%kap, 0, 0)
			
			P(0) = SymHartreeFockFieldProj_product_iso(density01_real%field%kap, HF_Delta, 0)
		
			ProjectedVC_Pair = ProjectedVC_Pair - 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
		
			CALL SymVCpp_get_Delta(HF_Delta, diagonal%consistency%vCpp, density10_imag%field%kap, 0, 0)
			
			P(0) = SymHartreeFockFieldProj_product_iso(density01_imag%field%kap, HF_Delta, 0)
		
			ProjectedVC_Pair = ProjectedVC_Pair + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
		
			CALL SymVCpp_get_Delta(HF_Delta, diagonal%consistency%vCpp, density10_real%field%kap, 0, 0)
			
			P(0) = SymHartreeFockFieldProj_product_iso(density01_imag%field%kap, HF_Delta, 0)
		
			ProjectedVC_Pair = ProjectedVC_Pair - 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
		
			CALL SymVCpp_get_Delta(HF_Delta, diagonal%consistency%vCpp, density10_imag%field%kap, 0, 0)
			
			P(0) = SymHartreeFockFieldProj_product_iso(density01_real%field%kap, HF_Delta, 0)
		
			ProjectedVC_Pair = ProjectedVC_Pair - 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
		
			! Spin-Orbit Pairing
			CALL SymVLSpp_get_Delta(HF_Delta, diagonal%consistency%vLSpp, density10_real%field%kap, 0, 0)
			
			P(0) = SymHartreeFockFieldProj_product_iso(density01_real%field%kap, HF_Delta, 0)
		
			ProjectedVLS_Pair = ProjectedVLS_Pair - 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
		
			CALL SymVLSpp_get_Delta(HF_Delta, diagonal%consistency%vLSpp, density10_imag%field%kap, 0, 0)
			
			P(0) = SymHartreeFockFieldProj_product_iso(density01_imag%field%kap, HF_Delta, 0)
		
			ProjectedVLS_Pair = ProjectedVLS_Pair + 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(P(0), 0)
		
			CALL SymVLSpp_get_Delta(HF_Delta, diagonal%consistency%vLSpp, density10_real%field%kap, 0, 0)
			
			P(0) = SymHartreeFockFieldProj_product_iso(density01_imag%field%kap, HF_Delta, 0)
		
			ProjectedVLS_Pair = ProjectedVLS_Pair - 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
		
			CALL SymVLSpp_get_Delta(HF_Delta, diagonal%consistency%vLSpp, density10_imag%field%kap, 0, 0)
			
			P(0) = SymHartreeFockFieldProj_product_iso(density01_real%field%kap, HF_Delta, 0)
		
			ProjectedVLS_Pair = ProjectedVLS_Pair - 0.5 * CoeffsXY%y_l(0)%index(IndexGauge) * CMPLX(0, P(0))
		
			CALL SymDensityProj_del(density01_real)
			CALL SymDensityProj_del(density10_real)
			
			CALL SymDensityProj_del(density01_imag)
			CALL SymDensityProj_del(density10_imag)
			
			CALL SymHartreeFockFieldProj_del(HF_Gamma_Diff)
			CALL SymHartreeFockFieldProj_del(HF_Delta)

		END DO
		
		DO IndexGauge = 1, NGauge_n
		
			! Convert projected density into the right type, which is SymDensityProj
		
			CALL SymDensity_new_GenDensityProj10(density10_real, genden_gauge_real_n(IndexGauge))
			CALL SymDensity_new_GenDensityProj10(density10_imag, genden_gauge_imag_n(IndexGauge))
			
			CALL SymHartreeFockFieldProj_new(HF_Gamma_Diff)
			CALL SymHartreeFockFieldProj_new(HF_Delta)
		
			!
			! Mean-field Channel
			!
			
			! Two-body kinetic energy (center of mass motion)
			CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vEkCMph, density10_real%field%rho, 1, 1, ProjectionOn)
		
			P(1) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedEkCM = ProjectedEkCM + (0.5 / A) * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
		
			CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vEkCMph, density10_imag%field%rho, 1, 1, ProjectionOn)
		
			P(1) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedEkCM = ProjectedEkCM - (0.5 / A) * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
		
			CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vEkCMph, density10_real%field%rho, 1, 1, ProjectionOn)
		
			P(1) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedEkCM = ProjectedEkCM + (0.5 / A) * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
		
			CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vEkCMph, density10_imag%field%rho, 1, 1, ProjectionOn)
		
			P(1) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedEkCM = ProjectedEkCM + (0.5 / A) * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
		
			! Brink-Boeker direct term
			CALL SymVBBph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10_real%field%rho, 1, 1, ProjectionOn)
		
			P(1) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVBB_local = ProjectedVBB_local + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
		
			CALL SymVBBph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10_imag%field%rho, 1, 1, ProjectionOn)
		
			P(1) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVBB_local = ProjectedVBB_local - 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
		
			CALL SymVBBph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10_real%field%rho, 1, 1, ProjectionOn)
		
			P(1) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVBB_local = ProjectedVBB_local + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
		
			CALL SymVBBph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10_imag%field%rho, 1, 1, ProjectionOn)
		
			P(1) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVBB_local = ProjectedVBB_local + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
		
			! Brink-Boeker exchange term
			CALL SymVBBph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10_real%field%rho, 1, 1, ProjectionOn)
		
			P(1) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVBB_exch = ProjectedVBB_exch + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
			
			CALL SymVBBph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10_imag%field%rho, 1, 1, ProjectionOn)
		
			P(1) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVBB_exch = ProjectedVBB_exch - 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
			
			CALL SymVBBph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10_real%field%rho, 1, 1, ProjectionOn)
		
			P(1) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVBB_exch = ProjectedVBB_exch + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
			
			CALL SymVBBph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vBBph, density10_imag%field%rho, 1, 1, ProjectionOn)
		
			P(1) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVBB_exch = ProjectedVBB_exch + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
			
			! Coulomb direct term
			CALL SymVCph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10_real%field%rho, 1, 1, ProjectionOn)

			P(1) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVC_local = ProjectedVC_local + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
			
			CALL SymVCph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10_imag%field%rho, 1, 1, ProjectionOn)

			P(1) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVC_local = ProjectedVC_local - 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
			
			CALL SymVCph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10_real%field%rho, 1, 1, ProjectionOn)

			P(1) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVC_local = ProjectedVC_local + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
			
			CALL SymVCph_get_LocalGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10_imag%field%rho, 1, 1, ProjectionOn)

			P(1) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVC_local = ProjectedVC_local + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
			
			! Coulomb exchange term
			CALL SymVCph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10_real%field%rho, 1, 1, ProjectionOn)

			P(1) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVC_exch = ProjectedVC_exch + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
		
			CALL SymVCph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10_imag%field%rho, 1, 1, ProjectionOn)

			P(1) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVC_exch = ProjectedVC_exch - 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
		
			CALL SymVCph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10_real%field%rho, 1, 1, ProjectionOn)

			P(1) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVC_exch = ProjectedVC_exch + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
		
			CALL SymVCph_get_ExchangeGamma(HF_Gamma_Diff, diagonal%consistency%vCph, density10_imag%field%rho, 1, 1, ProjectionOn)

			P(1) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVC_exch = ProjectedVC_exch + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
		
			! Spin-orbit term
			CALL SymVLSph_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vLSph, density10_real%field%rho, 1, 1, ProjectionOn)

			P(1) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVLS = ProjectedVLS + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
		
			CALL SymVLSph_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vLSph, density10_imag%field%rho, 1, 1, ProjectionOn)

			P(1) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVLS = ProjectedVLS - 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
		
			CALL SymVLSph_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vLSph, density10_real%field%rho, 1, 1, ProjectionOn)

			P(1) = SymHartreeFockFieldProj_product_iso(density10_imag%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVLS = ProjectedVLS + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
		
			CALL SymVLSph_get_Gamma(HF_Gamma_Diff, diagonal%consistency%vLSph, density10_imag%field%rho, 1, 1, ProjectionOn)

			P(1) = SymHartreeFockFieldProj_product_iso(density10_real%field%rho, HF_Gamma_Diff, 1)
		
			ProjectedVLS = ProjectedVLS + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
		
			!
			! Pairing Channel
			!
			
			CALL SymDensity_new_GenDensityProj01(density01_real, genden_gauge_real_n(IndexGauge))
			CALL SymDensity_new_GenDensityProj01(density01_imag, genden_gauge_imag_n(IndexGauge))
			
			!density01 = diagonal%consistency%density
			
                        ! Center of Mass Pairing
                        CALL SymKineticEnergy2Body_get_Delta(HF_Delta, diagonal%consistency%vEkCMpp, density10_real%field%kap, 1, 1, ProjectionOn)
			
			P(1) = SymHartreeFockFieldProj_product_iso(density01_real%field%kap, HF_Delta, 1)
		
			ProjectedEkCM_Pair = ProjectedEkCM_Pair - ( 0.5/A ) * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)

                        CALL SymKineticEnergy2Body_get_Delta(HF_Delta, diagonal%consistency%vEkCMpp, density10_imag%field%kap, 1, 1, ProjectionOn)
			
			P(1) = SymHartreeFockFieldProj_product_iso(density01_imag%field%kap, HF_Delta, 1)
		
			ProjectedEkCM_Pair = ProjectedEkCM_Pair + ( 0.5/A ) * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)

                        CALL SymKineticEnergy2Body_get_Delta(HF_Delta, diagonal%consistency%vEkCMpp, density10_real%field%kap, 1, 1, ProjectionOn)
			
			P(1) = SymHartreeFockFieldProj_product_iso(density01_imag%field%kap, HF_Delta, 1)
		
			ProjectedEkCM_Pair = ProjectedEkCM_Pair - ( 0.5/A ) * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))

                        CALL SymKineticEnergy2Body_get_Delta(HF_Delta, diagonal%consistency%vEkCMpp, density10_imag%field%kap, 1, 1, ProjectionOn)
			
			P(1) = SymHartreeFockFieldProj_product_iso(density01_real%field%kap, HF_Delta, 1)
		
			ProjectedEkCM_Pair = ProjectedEkCM_Pair - ( 0.5/A ) * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))

			! Brink-Boeker Pairing
			CALL SymVBBpp_get_Delta(HF_Delta, diagonal%consistency%vBBpp, density10_real%field%kap, 1, 1)
			
			P(1) = SymHartreeFockFieldProj_product_iso(density01_real%field%kap, HF_Delta, 1)
		
			ProjectedVBB_Pair = ProjectedVBB_Pair - 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
		
			CALL SymVBBpp_get_Delta(HF_Delta, diagonal%consistency%vBBpp, density10_imag%field%kap, 1, 1)
			
			P(1) = SymHartreeFockFieldProj_product_iso(density01_imag%field%kap, HF_Delta, 1)
		
			ProjectedVBB_Pair = ProjectedVBB_Pair + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
		
			CALL SymVBBpp_get_Delta(HF_Delta, diagonal%consistency%vBBpp, density10_real%field%kap, 1, 1)
			
			P(1) = SymHartreeFockFieldProj_product_iso(density01_imag%field%kap, HF_Delta, 1)
		
			ProjectedVBB_Pair = ProjectedVBB_Pair - 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
		
			CALL SymVBBpp_get_Delta(HF_Delta, diagonal%consistency%vBBpp, density10_imag%field%kap, 1, 1)
			
			P(1) = SymHartreeFockFieldProj_product_iso(density01_real%field%kap, HF_Delta, 1)
		
			ProjectedVBB_Pair = ProjectedVBB_Pair - 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
		
			! Coulomb Pairing
			CALL SymVCpp_get_Delta(HF_Delta, diagonal%consistency%vCpp, density10_real%field%kap, 1, 1)
			
			P(1) = SymHartreeFockFieldProj_product_iso(density01_real%field%kap, HF_Delta, 1)
		
			ProjectedVC_Pair = ProjectedVC_Pair - 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
		
			CALL SymVCpp_get_Delta(HF_Delta, diagonal%consistency%vCpp, density10_imag%field%kap, 1, 1)
			
			P(1) = SymHartreeFockFieldProj_product_iso(density01_imag%field%kap, HF_Delta, 1)
		
			ProjectedVC_Pair = ProjectedVC_Pair + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
		
			CALL SymVCpp_get_Delta(HF_Delta, diagonal%consistency%vCpp, density10_real%field%kap, 1, 1)
			
			P(1) = SymHartreeFockFieldProj_product_iso(density01_imag%field%kap, HF_Delta, 1)
		
			ProjectedVC_Pair = ProjectedVC_Pair - 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
		
			CALL SymVCpp_get_Delta(HF_Delta, diagonal%consistency%vCpp, density10_imag%field%kap, 1, 1)
			
			P(1) = SymHartreeFockFieldProj_product_iso(density01_real%field%kap, HF_Delta, 1)
		
			ProjectedVC_Pair = ProjectedVC_Pair - 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
		
			! Spin-Orbit Pairing
			CALL SymVLSpp_get_Delta(HF_Delta, diagonal%consistency%vLSpp, density10_real%field%kap, 1, 1)
			
			P(1) = SymHartreeFockFieldProj_product_iso(density01_real%field%kap, HF_Delta, 1)
		 
			ProjectedVLS_Pair = ProjectedVLS_Pair - 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
		
			CALL SymVLSpp_get_Delta(HF_Delta, diagonal%consistency%vLSpp, density10_imag%field%kap, 1, 1)
			
			P(1) = SymHartreeFockFieldProj_product_iso(density01_imag%field%kap, HF_Delta, 1)
		
			ProjectedVLS_Pair = ProjectedVLS_Pair + 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(P(1), 0)
		
			CALL SymVLSpp_get_Delta(HF_Delta, diagonal%consistency%vLSpp, density10_real%field%kap, 1, 1)
			
			P(1) = SymHartreeFockFieldProj_product_iso(density01_imag%field%kap, HF_Delta, 1)
		 
			ProjectedVLS_Pair = ProjectedVLS_Pair - 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
		
			CALL SymVLSpp_get_Delta(HF_Delta, diagonal%consistency%vLSpp, density10_imag%field%kap, 1, 1)
			
			P(1) = SymHartreeFockFieldProj_product_iso(density01_real%field%kap, HF_Delta, 1)
		
			ProjectedVLS_Pair = ProjectedVLS_Pair - 0.5 * CoeffsXY%y_l(1)%index(IndexGauge) * CMPLX(0, P(1))
		
			CALL SymDensityProj_del(density01_real)
			CALL SymDensityProj_del(density10_real)
			
			CALL SymDensityProj_del(density01_imag)
			CALL SymDensityProj_del(density10_imag)
			
			CALL SymHartreeFockFieldProj_del(HF_Gamma_Diff)
			CALL SymHartreeFockFieldProj_del(HF_Delta)

		END DO
		
		WRITE(*,'(/,14X,"----------------------------------------")')
		WRITE(*,'(14X,"| PROJECTION AFTER VARIATION - RESULTS |")')
		WRITE(*,'(14X,"----------------------------------------",/)')
		
		WRITE(*,'("                                 Energia cinetica:",2F20.10)') ProjectedEk
		WRITE(*,'("             Energia cinetica del centro de masas:",2F20.10)') ProjectedEkCM
		WRITE(*,'("                    Energia local de Brink-Booker:",2F20.10)') ProjectedVBB_local
		WRITE(*,'("                         Energia local de Coulomb:",2F20.10)') ProjectedVC_local
		WRITE(*,'("           Energia de intercambio de Brink-Booker:",2F20.10)') ProjectedVBB_exch
		WRITE(*,'("                Energia de intercambio de Coulomb:",2F20.10)') ProjectedVC_exch
		WRITE(*,'("               Energia dependiente de la densidad:",2F20.10)') ProjectedDD
		WRITE(*,'("                           Energia de spin-orbita:",2F20.10)') ProjectedVLS
		WRITE(*,'(70("-"))')

		eHF = ProjectedEk + ProjectedEkCM &
			+ ProjectedVBB_local + ProjectedVBB_exch &
			+ ProjectedVC_local +  ProjectedVC_exch&
			+ ProjectedDD + ProjectedVLS
			
		WRITE(*,'("                    Energia total de Brink-Booker:",2F20.10)') ProjectedVBB_local + ProjectedVBB_exch
		WRITE(*,'("                         Energia total de Coulomb:",2F20.10)') ProjectedVC_local + ProjectedVC_exch
		WRITE(*,'("                     Energia total (Hartree-Fock):",2F20.10)') eHF
		WRITE(*,'(70("-"))')

		pairing_BB = ProjectedVBB_Pair
		pairing_Coulomb = ProjectedVC_Pair
		pairing_LS = ProjectedVLS_Pair
                Pairing_EkCM = ProjectedEkCM_Pair
		
                WRITE(*,'("                  Apareamiento del centro de masa:",2F20.10)') pairing_EkCM
		WRITE(*,'("                     Apareamiento de Brink-Booker:",2F20.10)') pairing_BB
		WRITE(*,'("                          Apareamiento de Coulomb:",2F20.10)') pairing_Coulomb
		WRITE(*,'("                      Apareamiento de spin-orbita:",2F20.10)') pairing_LS
		WRITE(*,'(70("-"))')

		pairing = pairing_BB + pairing_Coulomb + pairing_LS + pairing_EkCM
		
		WRITE(*,'("                               Apareamiento total:",2F20.10)') pairing
		WRITE(*,'(70("-"))')

		WRITE(*,'("                                            TOTAL:",2F20.10)') eHF + pairing

		DO ta = 0, 1
		        
                        np = 0.0
                        R2 = 0.0

			DO aa = 0, 2*Lmax
		
				d = DIM(aa)
				dd = d + d
			
				la = L(aa)
				ja = J(aa)
			
				ALLOCATE(tmp(d, d))
				ALLOCATE(copyR2(d, d))
				ALLOCATE(copyRho(d, d))
			
				! The particle number is obtained by the trace of the density matrix
				trace = 0.0
				DO u1 = 1, d
					trace = trace + genden_proj_real%rho%rho(ta, aa)%store(u1, u1) + genden_proj_imag%rho%rho(ta, aa)%store(u1, u1)*CMPLX(0, 1)
				END DO
				np = np + DBLE(ja + 1.0) * REAL(trace)

				! We need to copy the matrix R2Field because of dimensions inconsistencies when in the
				! case of an arbitrary basis. This inconsistency causes the multiplication below to
				! crash when using -C option at compilation stage
				copyR2 = SymD3Tensor_matrix(R2Field, la)
				copyRho = genden_proj_real%rho%rho(ta, aa)%store + genden_proj_imag%rho%rho(ta, aa)%store * CMPLX(0,1)

				tmp = MATMUL(copyR2, copyRho)
			
				trace = 0.0
				DO u1 = 1, d
					trace = trace + REAL(tmp(u1, u1))
				END DO
				R2 = R2 + DBLE(ja + 1.0) * trace !/ b2

				! Liberamos la memoria de las matrices utilizadas
				DEALLOCATE(tmp)
				DEALLOCATE(copyR2)
				DEALLOCATE(copyRho)
				
				diagonal%consistency%density%nucleus%actual_np(ta) = np
				diagonal%consistency%density%nucleus%actual_R2(ta) = R2 / np
		
			END DO
		
		END DO

                NReal = diagonal%consistency%density%nucleus%actual_np(1)
                ZReal = diagonal%consistency%density%nucleus%actual_np(0)

		WRITE(*,'(/,"                               Neutrones       Protones          Total")') 
		
		WRITE(*,'("              Particulas:",3F15.5)') NReal, ZReal, NReal + ZReal
                WRITE(*,'("             Fermi Level:",3F15.5)') diagonal%consistency%density%nucleus%lambda_np(1), &
                                                                diagonal%consistency%density%nucleus%lambda_np(0)
		WRITE(*,'("                   Radio:",3F15.5)') SQRT(diagonal%consistency%density%nucleus%actual_R2(1)), &
								SQRT(diagonal%consistency%density%nucleus%actual_R2(0)), &
								SQRT(Nucleus_get_actual_R2(diagonal%consistency%density%nucleus))
		
		IF (Basis .EQ. 2) THEN

			ALLOCATE(DensityRadial(0:1,0:Npoint))
			
			DO ta = 0, 1
			
				DO aa = 0, 2*Lmax
		
					la = L(aa)
					ja = J(aa)
										d = MIN(Nmax, NmaxOfL(la))
					IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - la) / 2) + 1
				
					! Collecting the wave-functions: upper part (corresponding to the vector U) and lower part
					! (corresponding to the vector V)
					
					DO ipoint = 1, Npoint

						DO mm = 1, d
                                        
							IndexBasis = IndexVecNL(mm, la)
								
							! Calculating here the single-particle density. We have: 
							!
							!		\rho_{mn} =  { V (V+) }_{mn}
							!
							! and in coordinate representation:
							!
							!     rho(r) = \sum_{mn,l,j} (2j+1) \rho_{mn} \phi_{ml}(r)\phi_{nl}(r)
							!
							! or:
							!     rho(r) = \sum_{i, mn,l,j} (2j+1) V_{ni}V_{mi} \phi_{ml}(r)\phi_{nl}(r)
							!
						
							DO nn = 1, d
						
								IndexKet = IndexVecNL(nn, la)
								
								DensityRadial(ta, ipoint) = DensityRadial(ta, ipoint) &
												+ (ja + 1.0) * genden_proj_real%rho%rho(ta, aa)%store(mm, nn) &
												* WaveFun(ipoint,IndexBasis) &
												* WaveFun(ipoint,IndexKet)
								
							END DO
							
						END DO
                                        
					END DO
                                        
				END DO ! end of loop over aa
			
             			! Taking into account the fact that phi(r) = y(r)/r, and extrapolating at r=0 for the density
			
				! The basis wave-functions are in fact y(r) = R(r)/r, so we need to correct for this.
				! The 4*pi comes from the integration over the angles.
	
				DO ipoint=1, Npoint
					DensityRadial(ta, ipoint) = DensityRadial(ta, ipoint) / (4.0*PI*RadMesh(ipoint)**2)
				END DO

				! Extrapolate to find rho(r=0)

         			DensityRadial(ta, 0) = 3.0*(DensityRadial(ta, 1) - DensityRadial(ta, 2)) + DensityRadial(ta, 3)

			END DO ! end of loop over ta
		
			OPEN(file_unit_1, FILE='data/DensityPNP.dat', ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				WRITE(*,'("Impossible to open the file data/DensityPNP.dat")') 
				STOP "In DiagonalizationMethod_ProjectedEnergy - Impossible to open file"
			END IF
		
			DO ipoint = 0, Npoint
				WRITE(file_unit_1,'(3f20.16)') RadMesh(ipoint),DensityRadial(0, ipoint),DensityRadial(1, ipoint)
			END DO
                        
		        DEALLOCATE(DensityRadial)
			CLOSE(file_unit_1)
		
		END IF
					
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

END MODULE diagmeth
