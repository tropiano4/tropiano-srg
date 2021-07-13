 MODULE diagmeth

	USE input
	USE symd3t
	USE r1r1
	USE nucleus
	USE symden
	USE symgden
	USE symgdhf
	USE symke2b
	USE symvbb
	USE symvc
	USE symvls
	USE symgdd
	USE selfc
	USE eigenval
	USE jacobi
	USE indexx

	IMPLICIT NONE

	INCLUDE "brent_d1.f90"

	INTEGER, PARAMETER :: MAX_ITER = 5
	
	TYPE MatrixType
		DOUBLE PRECISION, POINTER, DIMENSION(:, :) :: value
	END TYPE

	TYPE DiagonalMatrix
		DOUBLE PRECISION, POINTER, DIMENSION(:) :: value
	END TYPE

	TYPE DiagonalizationMethod
		TYPE (SelfConsistencyMethod) consistency
		TYPE (R1R1Function) func
		TYPE (MatrixType), POINTER, DIMENSION(:, :) :: UV ! U and V vectors
		TYPE (DiagonalMatrix), POINTER, DIMENSION(:, :) :: QuasiParticleEnergies ! Self-explanatory...
		TYPE (SymGenDensity) S ! Superhamiltonian
		TYPE (SymGenDensity) iterated ! Density at each iteration
		INTEGER ta
	END TYPE

 CONTAINS

	INCLUDE "brent_d2.f90"

	! This subroutine initializes the self-consistent calculation by creating all necessary objects

	SUBROUTINE DiagonalizationMethod_new(diagonal, density)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		TYPE (SymDensity), INTENT(INOUT) :: density

		INTEGER ta, a, max_a, d, dd

		! Inicializamos las variables y subtipos
		CALL SelfConsistencyMethod_new(diagonal%consistency, density)
		CALL R1R1Function_new(diagonal%func)
		CALL SymGenDensity_new_Nucleus(diagonal%S, Nucleus_get_N(density%nucleus), Nucleus_get_Z(density%nucleus))
		CALL SymGenDensity_new_SymDensity(diagonal%iterated, density)

		max_a = 2*N_0
		
		! Reservamos memoria para las matrices del tipo
		ALLOCATE(diagonal%QuasiParticleEnergies(0:1, 0:max_a))
		IF (.NOT. ASSOCIATED(diagonal%QuasiParticleEnergies)) &
                        STOP "Unable to allocate memory in DiagonalizationMethod_new (1)"
		
		ALLOCATE(diagonal%UV(0:1, 0:max_a))
		IF (.NOT. ASSOCIATED(diagonal%UV)) &
                        STOP "Unable to allocate memory in DiagonalizationMethod_new (2)"
		
		DO ta = 0, 1
			DO a = 0, max_a
				d = DIM(a)
				dd = d + d
				ALLOCATE(diagonal%QuasiParticleEnergies(ta, a)%value(dd))
				IF (.NOT. ASSOCIATED(diagonal%QuasiParticleEnergies(ta, a)%value)) &
                                        STOP "Unable to allocate memory in DiagonalizationMethod_new (3)"
				
				ALLOCATE(diagonal%UV(ta, a)%value(dd, dd))
				IF (.NOT. ASSOCIATED(diagonal%UV(ta, a)%value))&
                                        STOP "Unable to allocate memory in DiagonalizationMethod_new (4)"
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
		TYPE (SymDensity) :: new_density
		TYPE (SymD3Tensor) :: ek_tensor
		
                DOUBLE PRECISION :: TestAccuracy
		
		LOGICAL :: Read_BBph, Read_BBpp, Read_Cph, Read_Cpp
		
		DOUBLE PRECISION :: accuracy
		DOUBLE PRECISION :: b, diff, factor, R2
		
		INTEGER :: cycles_in, cycles_out, cycles_rate, Lold, Lmin
		INTEGER :: A, N, Z, niter, ta

		TYPE (SymHartreeFockField) field1, field2

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
		
		Lmin = 0
						
		IF (N_0 < 10) THEN
			WRITE(diagonal%consistency%vBBph%filename, "(A,I1,A)") "data/vBB", N_0, "ph_HO.txt"
		ELSE
			WRITE(diagonal%consistency%vBBph%filename, "(A,I2,A)") "data/vBB", N_0, "ph_HO.txt"
		END IF
		
		Read_BBph = SymVBBph_read(diagonal%consistency%vBBph)
		Read_BBpp = SymVBBpp_read(diagonal%consistency%vBBpp)
		
		CALL SymVBBph_calculate(diagonal%consistency%vBBph, Read_BBph, Lmin)
		CALL SymVBBpp_calculate(diagonal%consistency%vBBpp, Read_BBpp, Lmin)
		
		Read_Cph = SymVCph_read(diagonal%consistency%vCph)
		Read_Cpp = SymVCpp_read(diagonal%consistency%vCpp)
		
		CALL SymVCph_calculate(diagonal%consistency%vCph, Read_Cph, Lmin)
		CALL SymVCpp_calculate(diagonal%consistency%vCpp, Read_Cpp, Lmin)
		
		! Creation of the tensor fields

		CALL SymHartreeFockField_new(field1)
		CALL SymHartreeFockField_new(field2)

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

			CALL DiagonalizationMethod_get_MeanField(diagonal, diagonal%S)

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
				
			END IF

			! Store the density for an eventual restart
			CALL SymDensity_save(diagonal%consistency%density)

			! Store the energies (for future use) of this step after correction for the particle number (above)
			
			CALL SelfConsistencyMethod_store_eHFB(diagonal%consistency)

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
			  
			CALL SymDensity_new_GenDensity(new_density, diagonal%iterated)
			
                        ! The current step corresponds to the density diagonal%consistency%density, the new one to new_density
			! We compute here the "distance" between the two.
			  
			IF (MOD(niter, MAX_ITER) .EQ. 0) THEN
				diff = SymHartreeFockField_distance(diagonal%consistency%density%field%rho, new_density%field%rho)
				WRITE(*,'("          k = ",I4," Difference of the densities : ",F15.8)') niter,diff
			END IF
			
			niter = niter + 1

                        ! We reduce the step made into the direction of the new density, new_density, from the old one, 
			! diagonal%consistency%density, with some slowing-factor anneal
			!
			!    rho_{n+1} -> (alpha * rho_{n} + rho_{n+1})/(1 + alpha)
			!
			  
			factor = 1.0 / (1.0 + anneal)
			
			! new step for normal density
			
			CALL SymHartreeFockField_product(field1, anneal, diagonal%consistency%density%field%rho)
			CALL SymHartreeFockField_add(field2, new_density%field%rho, field1)
			CALL SymHartreeFockField_product(diagonal%consistency%density%field%rho, factor, field2)

			CALL SymHartreeFockField_product(field1, anneal, diagonal%consistency%density%field%kap)
			CALL SymHartreeFockField_add(field2, new_density%field%kap, field1)
			CALL SymHartreeFockField_product(diagonal%consistency%density%field%kap, factor, field2)

			! We update the numbers of particles with the same rule

			DO ta = 0, 1
				diagonal%consistency%density%nucleus%actual_np(ta) = factor * ((new_density%nucleus%actual_np(ta)) &
					+ (anneal * diagonal%consistency%density%nucleus%actual_np(ta)))
                        
			END DO

			accuracy = SelfConsistencyMethod_accuracy(diagonal%consistency)
			WRITE(*,'("Iteration k = ",I4," HF Energy : ",F15.8)') niter,diagonal%consistency%density%nucleus%eHFB
			
			IF (diff .LE. tolerance) THEN
				accuracy = SelfConsistencyMethod_accuracy(diagonal%consistency)
                                WRITE(*,'("Energy = ",EN15.5," Precision = ",ES12.5)') &
                                                diagonal%consistency%density%nucleus%eHFB,accuracy
				IF (accuracy .LE. tolerance) THEN
					EXIT gsc
				END IF
			END IF
			
			IF (niter .GE. NITER_MAX) THEN
				EXIT gsc
		        END IF
			
		END DO gsc

		CALL SymHartreeFockField_del(field1)
		CALL SymHartreeFockField_del(field2)
		
		! Store the density for an eventual restart
		CALL SymDensity_save(diagonal%consistency%density)

		CALL SYSTEM_CLOCK(cycles_out, cycles_rate)
		PRINT "(/A,EN10.2)", "Tiempo empleado (segundos): ", (DBLE(cycles_out - cycles_in) / cycles_rate)

		CALL SymDensity_store_actual_R2(diagonal%consistency%density)

		CALL SelfConsistencyMethod_store_eHFB(diagonal%consistency)
		CALL SelfConsistencyMethod_show_Status(diagonal%consistency)
		
		CALL DiagonalizationMethod_show_ParticleEnergies(diagonal)
                CALL DiagonalizationMethod_show_QuasiParticleEnergies(diagonal)
		
		RETURN
	END SUBROUTINE DiagonalizationMethod_goto_SelfConsistency

	!-----------------------------------------------------------------------------------------------!
	!   In this subroutine, we calculate the mean-field Gamma and pairing field Delta from the 	!
	!   matrix elements of the force and the densities. Both the matrix elements and the densities	!
	!   at a given iteration of the HFB process are stored in the object "diagonal". Once the fields!
	!   are constructed, we build an object that will contain the super-hamiltonian S.
	!-----------------------------------------------------------------------------------------------!

	SUBROUTINE DiagonalizationMethod_get_MeanField(diagonal, S)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		TYPE (SymGenDensity), INTENT(INOUT) :: S

		DOUBLE PRECISION :: b, factor
		INTEGER :: A

		TYPE (SymHartreeFockField) :: HF_Gamma, HF_Delta
		TYPE (SymHartreeFockField) :: ekcm_field, vbb_field, vc_field, vls_field, gdd_field, field1, field2
		TYPE (SymD3Tensor) :: ek_tensor
		TYPE (SymGenDensityHF) :: gendenhf_gamma, gendenhf_delta

		b = Nucleus_get_b(diagonal%consistency%density%nucleus)
		A = Nucleus_get_A(diagonal%consistency%density%nucleus)
		
		IF ((A .LE. 1) .OR. (A .GE. 300)) THEN
			PRINT *, "Unexpected A value = ", A
			STOP "DiagonalizationMethod::get_MeanField"
		END IF

		! Creating all the necessary fields

		CALL SymHartreeFockField_new(HF_Gamma)
		CALL SymHartreeFockField_new(HF_Delta)
		CALL SymD3Tensor_new(ek_tensor)
		CALL SymHartreeFockField_new(ekcm_field)
		CALL SymHartreeFockField_new(vbb_field)
		CALL SymHartreeFockField_new(vc_field)
		CALL SymHartreeFockField_new(vls_field)
		CALL SymHartreeFockField_new(gdd_field)
		CALL SymHartreeFockField_new(field1)
		CALL SymHartreeFockField_new(field2)

		factor = 1.0 - (1.0 / A)

		! Mean-field - Kinetic energy

                CALL SymD3Tensor_product(ek_tensor, DBLE(factor), EkField)

		! Mean-field - Two-body center of mass correction

		CALL SymKineticEnergy2Body_get_Gamma(field1, diagonal%consistency%vEkCMph, diagonal%consistency%density%field%rho)
		CALL SymHartreeFockField_product(ekcm_field, DBLE(1.0 / A), field1)
		
		! Mean-field - Brink-Boker term

		CALL SymVBBph_get_Gamma(vbb_field, diagonal%consistency%vBBph, diagonal%consistency%density%field%rho)

		CALL SymVCph_get_Gamma(vc_field, diagonal%consistency%vCph, diagonal%consistency%density%field%rho)

		! Mean-field - Spin-orbit term

		CALL SymVLSph_get_Gamma(vls_field, diagonal%consistency%vLSph, diagonal%consistency%density%field%rho)

		! Mean-field - Density-dependent term

		CALL SymGDDph_update(gdd_field, diagonal%consistency%gDDph, diagonal%consistency%density%field%rho)

		! Total Mean-field = Sum of all the preceding terms

		CALL SymHartreeFockField_add(field1, vls_field, gdd_field)
		CALL SymHartreeFockField_add(field2, vc_field, field1)
		CALL SymHartreeFockField_add(field1, vbb_field, field2)
		CALL SymHartreeFockField_add(field2, ekcm_field, field1)
		CALL SymHartreeFockField_add_SymD3Tensor(HF_Gamma, ek_tensor, field2)

		! Pairing - Center of Mass term

		IF (HFOnly .EQ. 0) THEN
			CALL SymKineticEnergy2Body_get_Delta(field1, diagonal%consistency%vEkCMpp, diagonal%consistency%density%field%kap)
			CALL SymHartreeFockField_product(ekcm_field, DBLE(1.0 / A), field1)
			!CALL SymHartreeFockField_product(ekcm_field, 0.0, field1)
		ELSE
			CALL SymKineticEnergy2Body_get_Delta(field1, diagonal%consistency%vEkCMpp, diagonal%consistency%density%field%kap)
			CALL SymHartreeFockField_product(ekcm_field, 0.0, field1)
		END IF

		! Pairing - Brink-Boker term

		IF (HFOnly .EQ. 0) THEN
			CALL SymVBBpp_get_Delta(vbb_field, diagonal%consistency%vBBpp, diagonal%consistency%density%field%kap)
		ELSE
			CALL SymVBBpp_get_Delta(field1, diagonal%consistency%vBBpp, diagonal%consistency%density%field%kap)
			CALL SymHartreeFockField_product(vbb_field, 0.0, field1)
		END IF

		! Pairing - Coulomb term

		IF (HFOnly .EQ. 0) THEN
			CALL SymVCpp_get_Delta(vc_field, diagonal%consistency%vCpp, diagonal%consistency%density%field%kap)
		ELSE
			CALL SymVCpp_get_Delta(field1, diagonal%consistency%vCpp, diagonal%consistency%density%field%kap)
			CALL SymHartreeFockField_product(vc_field, 0.0, field1)
		END IF

		! Pairing - Spin-orbit term

		IF (HFOnly .EQ. 0) THEN
			CALL SymVLSpp_get_Delta(vls_field, diagonal%consistency%vLSpp, diagonal%consistency%density%field%kap)
		ELSE
			CALL SymVLSpp_get_Delta(field1, diagonal%consistency%vLSpp, diagonal%consistency%density%field%kap)
			CALL SymHartreeFockField_product(vls_field, 0.0, field1)
		END IF

		! Total Pairing = Sum of all the preceding terms

		CALL SymHartreeFockField_add(field1, vc_field, vls_field)
		CALL SymHartreeFockField_add(field2, ekcm_field, field1)
		CALL SymHartreeFockField_add(HF_Delta, vbb_field, field2)

		! Creating an object that contain the super-hamiltonian to diagonalize

		CALL SymGenDensityHF_new(gendenhf_gamma)
		CALL SymGenDensityHF_new(gendenhf_delta)
		
		CALL SymGenDensityHF_copy(gendenhf_gamma, HF_Gamma)
		CALL SymGenDensityHF_copy(gendenhf_delta, HF_Delta)

		CALL SymGenDensity_new_GammaDelta(S, gendenhf_gamma, gendenhf_delta, b)

		CALL SymGenDensityHF_del(gendenhf_gamma)
		CALL SymGenDensityHF_del(gendenhf_delta)

		CALL SymD3Tensor_del(ek_tensor)
		CALL SymHartreeFockField_del(HF_Gamma)
		CALL SymHartreeFockField_del(HF_Delta)
		CALL SymHartreeFockField_del(ekcm_field)
		CALL SymHartreeFockField_del(vbb_field)
		CALL SymHartreeFockField_del(vc_field)
		CALL SymHartreeFockField_del(vls_field)
		CALL SymHartreeFockField_del(gdd_field)
		CALL SymHartreeFockField_del(field1)
		CALL SymHartreeFockField_del(field2)
		
		RETURN
	END SUBROUTINE DiagonalizationMethod_get_MeanField

	! Changing the isospin under consideration

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

		DOUBLE PRECISION np, R2, b2
		INTEGER a, d, dd, la, ja, i

		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: xI, h, Delta, SH, U, V, tmp, copyR2
		INTEGER u1, u2
		DOUBLE PRECISION trace

		np = 0.0
		R2 = 0.0
		b2 = Nucleus_get_b(diagonal%iterated%nucleus) ** 2

		! Block structure: a refers to the angular momentum, but remember we deal here with a 2N x 2N matrix,
		! if N is the size of the single-particle basis.

		DO a = 0, 2*N_0
		
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
					SH(u1    , u2 + d) = - Delta(u1, u2)
					SH(u1 + d, u2    ) = - Delta(u1, u2)
					SH(u1 + d, u2 + d) = - h(u1, u2)
				END DO
			END DO

			! Diagonalizing the super-hamiltonian SH gives the quasi-particle energies
			! ordered from the lowest to the highest
			CALL EigenValues(dd, dd, SH, &
				diagonal%QuasiParticleEnergies(diagonal%ta, a)%value, &
				diagonal%UV(diagonal%ta, a)%value)

			! Extracting the matrices U and V
			DO u1 = 1, d
				DO u2 = 1, d
					U(u1, u2) = diagonal%UV(diagonal%ta, a)%value(u1    , u2 + d)
					V(u1, u2) = diagonal%UV(diagonal%ta, a)%value(u1 + d, u2 + d)
				END DO
			END DO

			! We calculate the new densities rho and kappa from the U and V matrices
			CALL SymGenDensity_make_Block(diagonal%iterated, diagonal%ta, a, U, V)

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

			tmp = MATMUL(copyR2, diagonal%iterated%rho%rho(diagonal%ta, a)%store)
			
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
	
	SUBROUTINE DiagonalizationMethod_show_ParticleEnergies(diagonal)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal

		TYPE (SymGenDensity) G, SuperHamiltonian
		
		TYPE (DiagonalMatrix), DIMENSION(0:1) :: D, Occup, SingleParticleE, SingleParticleP, SingleParticleR2
		TYPE (MatrixType), DIMENSION(0:1) :: P, E, R2
				
		INTEGER, DIMENSION(:), ALLOCATABLE :: an, ap
		INTEGER, DIMENSION(:), ALLOCATABLE :: in, ip
		
		INTEGER, DIMENSION(:,:), ALLOCATABLE :: Count
		
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: S, V, R2Matrix, R2Initial
		
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
		DO a = 0, 2*N_0
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
		
		! "SuperHamiltonian" contains the generalized hamiltonian 
		!
		!                 (     h      -Delta  )
		!                 (                    )
		!                 (  -Delta     -h     )
		!
		!  SuperHamiltonian%rho = h (type of SymHartreeFockField)
		!  SuperHamiltonian%kap = Delta (type of SymHartreeFockField)

		CALL DiagonalizationMethod_get_MeanField(diagonal, SuperHamiltonian)

		ndim = (N_0 + 1) * (N_0 + 2) / 2
		
		ALLOCATE(an(ndim + 1))
		ALLOCATE(ap(ndim + 1))

		DO ta = 0, 1
			
			dim_acc = 0
			num = 1
		
			DO a = 0, 2*N_0
		
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
				ALLOCATE(R2Initial(((N_0 - la) / 2) + 1, ((N_0 - la) / 2) + 1))

				ALLOCATE(D(ta)%value(dim_cur))
				ALLOCATE(E(ta)%value(dim_cur, dim_cur))
				ALLOCATE(P(ta)%value(dim_cur, dim_cur))
				ALLOCATE(R2(ta)%value(dim_cur, dim_cur))
				
				! S contains the matrix of the density rho for the quantum number a = (la,ja) at the convergence

				S = diagonal%iterated%rho%rho(ta, a)%store
				
				!  After diagonalization, D(ta)% value contains the diagonal values of the p-h density rho 
				!  (the occupation probabilities v2), V the eigenvectors that make pass from the HFB basis
				!  to the canonical basis (where the matrix of rho is diagonal). 
				!
				!  The s.p. energies are by definition the expectation value of the hamiltonian in the 
				!  canonical basis. Therefore they read: (V+)*S*V

				CALL Jacobi_real8(S, dim_cur, D(ta)%value, V, nrot)

				E(ta)%value = MATMUL(TRANSPOSE(V), MATMUL(SuperHamiltonian%rho%rho(ta, a)%store, V))
				P(ta)%value = MATMUL(TRANSPOSE(V), MATMUL(SuperHamiltonian%kap%rho(ta, a)%store, V))
				
				R2Initial = SymD3Tensor_matrix(R2Field, L(a))
				R2Matrix = R2Initial(1:dim_cur,1:dim_cur)
				
				R2(ta)%value = MATMUL(TRANSPOSE(V), MATMUL(R2Matrix, V))

				DO i = 1, dim_cur
				
					Occup(ta)%value(dim_acc + i) = D(ta)%value(i)
					
					SingleParticleE(ta)%value(dim_acc + i) = E(ta)%value(i, i)
					SingleParticleP(ta)%value(dim_acc + i) = P(ta)%value(i, i)
					SingleParticleR2(ta)%value(dim_acc + i) = R2(ta)%value(i, i)

				END DO

				DEALLOCATE(D(ta)%value)
				DEALLOCATE(E(ta)%value)
				DEALLOCATE(P(ta)%value)
				DEALLOCATE(R2(ta)%value)

				DEALLOCATE(S)
				DEALLOCATE(V)
				DEALLOCATE(R2Matrix)
				DEALLOCATE(R2Initial)

				dim_acc = dim_acc + dim_cur
			
			END DO ! end of loop over a

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
			
                        WRITE(*,'(I4,")",2X,A2,I2,"/2",3x,F10.8,F10.3,F8.3,F8.3)') &
					i, spectr(Lvalue), Jvalue, Occup(1)%value(in(i)), SingleParticleE(1)%value(in(i))
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
                        WRITE(*,'(I4,")",2X,A2,I2,"/2",3X,F10.8,F10.3,F8.3,F8.3)') &
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

                DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: nVV, pVV
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: U, V, QPE, Gamma, Delta

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
		DO a = 0, 2*N_0
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
		
  		Nnucl = Nucleus_get_N(diagonal%consistency%density%nucleus)
 		Znucl = Nucleus_get_Z(diagonal%consistency%density%nucleus)
			
		Anucl = Nnucl + Znucl
		factor = 1.0 - (1.0 / Anucl)

		CALL SymD3Tensor_new(ek_tensor)
		
		DO ta = 0, 1
			
			num = 0
			state = 0
                
			DO a = 0, 2*N_0
			
				la = L(a)
				ja = J(a)

				d = ((N_0 - la) / 2) + 1
			
				! Collecting the q.-p. energies and calculating the occupation factors
                        
				DO s1 = d + 1, d + d
				
                                	num = num + 1
				
					QPE(ta, num) = diagonal%QuasiParticleEnergies(ta, a)%value(s1)
                                
					Momentum(ta, num) = a
                                
					IF (ta .EQ. 1) nIndx(num) = num
					IF (ta .EQ. 0) pIndx(num) = num
                                
					sumvv = 0.0
					DO sa = 1, d
						sumvv = sumvv + (diagonal%UV(ta, a)%value(sa + d, s1) ** 2)
					END DO

					IF ((diagonal%consistency%density%nucleus%is_blocking(ta)) .AND. &
						(a .EQ. diagonal%consistency%density%nucleus%ia(ta)) .AND. &
						(s1 .EQ. (d + diagonal%consistency%density%nucleus%mu0(ta)))) THEN

						DO sa = 1, d
							sumvv = sumvv +((diagonal%UV(ta, a)%value(sa, s1) ** 2) &
								      - (diagonal%UV(ta, a)%value(sa + d, s1) ** 2)) / (ja + 1.0)
						END DO

					END IF
                                
					IF (ta .EQ. 1) nVV(num) = sumvv
					IF (ta .EQ. 0) pVV(num) = sumvv
					
				END DO
						                                
			END DO ! end of loop over a
			
 1              	CONTINUE   
 
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

                        WRITE(*,'(i4,")",i5,F10.4,3X,F10.8,3X,A1,I2,"/2",A5,I2,"/2",3x,F10.8,F10.4)') &
				num, nIndx(num), QPE(1, nIndx(num)), nVV(nIndx(num)), cn, jn, cp, jp, pVV(pIndx(num)), QPE(0, pIndx(num))
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
		
		DEALLOCATE(nIndx)
		DEALLOCATE(pIndx)

		DEALLOCATE(QPE)

		DEALLOCATE(nVV)
		DEALLOCATE(pVV)

		DEALLOCATE(Momentum)

		RETURN
                
	END SUBROUTINE DiagonalizationMethod_show_QuasiParticleEnergies

	SUBROUTINE DiagonalizationMethod_del(diagonal)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal

		INTEGER ta, a, max_a

		max_a = 2*N_0
		DO ta = 0, 1
			DO a = 0, max_a
				DEALLOCATE(diagonal%QuasiParticleEnergies(ta, a)%value)
				DEALLOCATE(diagonal%UV(ta, a)%value)
			END DO
		END DO
		DEALLOCATE(diagonal%QuasiParticleEnergies)
		DEALLOCATE(diagonal%UV)

		CALL SelfConsistencyMethod_del(diagonal%consistency)
		CALL R1R1Function_del(diagonal%func)
		CALL SymGenDensity_del(diagonal%S)
		CALL SymGenDensity_del(diagonal%iterated)
		RETURN
	END SUBROUTINE DiagonalizationMethod_del

END MODULE diagmeth
