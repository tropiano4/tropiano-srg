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
	USE wave
	USE eigenval
	USE jacobi
	USE indexx

	IMPLICIT NONE

	INCLUDE "brent_d1.f90"

	INTEGER, PARAMETER :: MAX_ITER = 4

	TYPE MatrixType
		DOUBLE PRECISION, POINTER, DIMENSION(:, :) :: value
	END TYPE

	TYPE DiagonalMatrix
		DOUBLE PRECISION, POINTER, DIMENSION(:) :: value
	END TYPE

	TYPE DiagonalizationMethod
		TYPE (SelfConsistencyMethod) consistency
		TYPE (R1R1Function) func
		TYPE (MatrixType), POINTER, DIMENSION(:, :) :: UV
		TYPE (DiagonalMatrix), POINTER, DIMENSION(:, :) :: QuasiParticleEnergies
		TYPE (SymGenDensity) S ! Superhamiltoniano sin el multiplicador del radio
		TYPE (SymGenDensity) iterated ! Densidad
		INTEGER ta
	END TYPE

CONTAINS

	INCLUDE "brent_d2.f90"

	SUBROUTINE DiagonalizationMethod_new(diagonal, density)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		TYPE (SymDensity), INTENT(INOUT) :: density

		INTEGER ta, a, max_a, d, dd

		! Inicializamos las variables y subtipos
		CALL SelfConsistencyMethod_new(diagonal%consistency, density)
		CALL R1R1Function_new(diagonal%func)
		CALL SymGenDensity_new_Nucleus(diagonal%S, &
			Nucleus_get_N(density%nucleus), &
			Nucleus_get_Z(density%nucleus))
		CALL SymGenDensity_new_SymDensity(diagonal%iterated, density)

		max_a = 2 * N_0
		! Reservamos memoria para las matrices del tipo
		ALLOCATE(diagonal%QuasiParticleEnergies(0:1, 0:max_a))
		IF (.NOT. ASSOCIATED(diagonal%QuasiParticleEnergies)) STOP "Unable to allocate memory"
		ALLOCATE(diagonal%UV(0:1, 0:max_a))
		IF (.NOT. ASSOCIATED(diagonal%UV)) STOP "Unable to allocate memory"
		DO ta = 0, 1
			DO a = 0, max_a
				d = DIM(a)
				dd = d + d
				ALLOCATE(diagonal%QuasiParticleEnergies(ta, a)%value(dd))
				IF (.NOT. ASSOCIATED(diagonal%QuasiParticleEnergies(ta, a)%value)) STOP "Unable to allocate memory"
				ALLOCATE(diagonal%UV(ta, a)%value(dd, dd))
				IF (.NOT. ASSOCIATED(diagonal%UV(ta, a)%value)) STOP "Unable to allocate memory"
			END DO
		END DO
		RETURN
	END SUBROUTINE DiagonalizationMethod_new

	SUBROUTINE DiagonalizationMethod_get_SymDensity(diagonal, density)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		TYPE (SymDensity), INTENT(INOUT) :: density

		CALL SymDensity_new_GenDensity(density, diagonal%iterated)
		RETURN
	END SUBROUTINE DiagonalizationMethod_get_SymDensity

	SUBROUTINE DiagonalizationMethod_goto_SelfConsistency(diagonal, tolerance)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		DOUBLE PRECISION, INTENT(IN) :: tolerance

		DOUBLE PRECISION b, diff, factor, anneal, R2
		INTEGER A, niter, N, Z, ta
		TYPE (OneDimSolve) neutron_constrainer, proton_constrainer
		TYPE (SymDensity) new_density
		DOUBLE PRECISION accuracy
		INTEGER cycles_in, cycles_out, cycles_rate

		TYPE (SymHartreeFockField) field1, field2

		PRINT "(/A)", "Iniciando calculo autoconsistente"
		PRINT "(A,I2)", "Numero de capas: ", N_0
		PRINT "(A,ES8.2)", "Tolerancia: ", tolerance
		! Leemos el reloj actual del sistema
		CALL SYSTEM_CLOCK(cycles_in)

		diff = 0.0
		b = Nucleus_get_b(diagonal%consistency%density%nucleus)
		PRINT "(A,F8.5/)", "Longitud del oscilador (b) = ", b
		A = Nucleus_get_A(diagonal%consistency%density%nucleus)
		IF ((A .LE. 1) .OR. (A .GE. 300)) THEN
			PRINT *, "Unexpected A value = ", A
			STOP "DiagonalizationMethod_goto_SelfConsistency"
		END IF
		factor = 1.0 - (1.0 / A)

		CALL SymVBBph_calculate(diagonal%consistency%vBBph)
		CALL SymVBBpp_calculate(diagonal%consistency%vBBpp)
		CALL SymVCph_calculate(diagonal%consistency%vCph)
		CALL SymVCpp_calculate(diagonal%consistency%vCpp)
		CALL SymVLSph_calculate(diagonal%consistency%vLSph)
		CALL SymVLSpp_calculate(diagonal%consistency%vLSpp)

		CALL SymHartreeFockField_new(field1)
		CALL SymHartreeFockField_new(field2)

		anneal = 0.2
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

			CALL DiagonalizationMethod_set_ISOSpin(diagonal, 1)

			CALL OneDimSolve_new(neutron_constrainer, diagonal%func, diagonal)
			diagonal%consistency%density%nucleus%lambda_np(1) = OneDimSolve_Solve(neutron_constrainer, &
				DBLE(N), diagonal%consistency%density%nucleus%lambda_np(1), &
				0.5 * (diagonal%consistency%density%nucleus%np(1) - diagonal%consistency%density%nucleus%actual_np(1)) + 0.1)

			CALL DiagonalizationMethod_set_ISOSpin(diagonal, 0)

			CALL OneDimSolve_new(proton_constrainer, diagonal%func, diagonal)
			diagonal%consistency%density%nucleus%lambda_np(0) = OneDimSolve_Solve(proton_constrainer, &
				DBLE(Z), diagonal%consistency%density%nucleus%lambda_np(0), &
				0.5 * (diagonal%consistency%density%nucleus%np(0) - diagonal%consistency%density%nucleus%actual_np(0)) + 0.1)

			CALL SymDensity_new_GenDensity(new_density, diagonal%iterated)

			IF (MOD(niter, MAX_ITER) .EQ. 0) THEN
				diff = SymHartreeFockField_distance(diagonal%consistency%density%field%rho, new_density%field%rho)
				PRINT "(A12,I4,A14,ES15.5)", "Iteraciones:", niter, "Diferencia:", diff
			END IF
			niter = niter + 1

			IF (MOD(niter, 50) .EQ. 0) THEN
				anneal = anneal + 0.1
				PRINT "(A,F8.5)", "El annealing se aumento a ", anneal
			END IF

			factor = 1.0 / (1.0 + anneal)
			CALL SymHartreeFockField_product(field1, anneal, diagonal%consistency%density%field%rho)
			CALL SymHartreeFockField_add(field2, new_density%field%rho, field1)
			CALL SymHartreeFockField_product(diagonal%consistency%density%field%rho, factor, field2)

			CALL SymHartreeFockField_product(field1, anneal, diagonal%consistency%density%field%kap)
			CALL SymHartreeFockField_add(field2, new_density%field%kap, field1)
			CALL SymHartreeFockField_product(diagonal%consistency%density%field%kap, factor, field2)

			DO ta = 0, 1
				diagonal%consistency%density%nucleus%actual_np(ta) = factor	* ((new_density%nucleus%actual_np(ta)) &
					+ (anneal * diagonal%consistency%density%nucleus%actual_np(ta)))
			END DO

			IF (MOD(niter, (10 * MAX_ITER)) .EQ. 0) THEN
				CALL SymDensity_save(diagonal%consistency%density)
			END IF

			IF ((diff .LE. tolerance) .AND. (MOD(niter, MAX_ITER) .EQ. 0)) THEN
				accuracy = SelfConsistencyMethod_accuracy(diagonal%consistency)
				PRINT "(A30,EN15.5,A12,ES12.5)", "Energia:", diagonal%consistency%density%nucleus%eHFB, "Precision:", accuracy
				IF (accuracy .LE. tolerance) THEN
					EXIT gsc
				END IF
			END IF
		END DO gsc

		CALL SymHartreeFockField_del(field1)
		CALL SymHartreeFockField_del(field2)

		CALL SYSTEM_CLOCK(cycles_out, cycles_rate)
		PRINT "(/A,EN10.2)", "Tiempo empleado (segundos): ", (DBLE(cycles_out - cycles_in) / cycles_rate)

		CALL SymDensity_store_actual_R2(diagonal%consistency%density)

		CALL SelfConsistencyMethod_store_eHFB(diagonal%consistency)
		CALL SelfConsistencyMethod_show_Status(diagonal%consistency)
		RETURN
	END SUBROUTINE DiagonalizationMethod_goto_SelfConsistency

	SUBROUTINE DiagonalizationMethod_get_MeanField(diagonal, S)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		TYPE (SymGenDensity), INTENT(INOUT) :: S

		DOUBLE PRECISION b, factor
		INTEGER A

		TYPE (SymHartreeFockField) HF_Gamma, HF_Delta
		TYPE (SymD3Tensor) ek_tensor
		TYPE (SymHartreeFockField) ekcm_field, vbb_field, vc_field, vls_field, gdd_field, field1, field2
		TYPE (SymGenDensityHF) gendenhf_gamma, gendenhf_delta

!	integer la,u1,u2

		b = Nucleus_get_b(diagonal%consistency%density%nucleus)
		A = Nucleus_get_A(diagonal%consistency%density%nucleus)
		IF ((A .LE. 1) .OR. (A .GE. 300)) THEN
			PRINT *, "Unexpected A value = ", A
			STOP "DiagonalizationMethod::get_MeanField"
		END IF

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

		CALL SymD3Tensor_product(ek_tensor, DBLE(factor / (b * b)), EkField)
!	print *,"Ek:"
!	do la=0,N_0
!		print *,"la=",la
!		do u1=1,((N_0-la)/2)+1
!			print "(5F10.5)",(ek_tensor%d3tensor(la)%d2(u1,u2),u2=1,((N_0-la)/2)+1)
!		end do
!	end do

		CALL SymKineticEnergy2Body_product(field1, diagonal%consistency%vEkCMph, diagonal%consistency%density%field%rho)
		CALL SymHartreeFockField_product(ekcm_field, DBLE(1.0 / (A * (b * b))), field1)
!	print *,"EkCM:", DBLE(1.0 / (A * (b * b)))
!	do la=0,N_0
!		print *,"la=",la
!		do u1=1,((N_0-la)/2)+1
!			print "(5F10.5)",(ekcm_field%p(0)%d3tensor(la)%d2(u1,u2),u2=1,((N_0-la)/2)+1)
!		end do
!	end do

		CALL SymVBBph_product(vbb_field, diagonal%consistency%vBBph, diagonal%consistency%density%field%rho)
!	print *,"BB:"
!	do la=0,N_0
!		print *,"la=",la
!		do u1=1,((N_0-la)/2)+1
!			print "(5F10.5)",(vbb_field%p(0)%d3tensor(la)%d2(u1,u2),u2=1,((N_0-la)/2)+1)
!		end do
!	end do

		CALL SymVCph_product(field1, diagonal%consistency%vCph, diagonal%consistency%density%field%rho)
		CALL SymHartreeFockField_product(vc_field, DBLE(1.0 / b), field1)
!	print *,"C:"
!	do la=0,N_0
!		print *,"la=",la
!		do u1=1,((N_0-la)/2)+1
!			print "(5F10.5)",(vc_field%p(0)%d3tensor(la)%d2(u1,u2),u2=1,((N_0-la)/2)+1)
!		end do
!	end do

		CALL SymVLSph_product(field1, diagonal%consistency%vLSph, diagonal%consistency%density%field%rho)
		CALL SymHartreeFockField_product(vls_field, DBLE(1.0 / (b ** 5.0)), field1)
!	print *,"LS:"
!	do la=0,N_0
!		print *,"la=",la
!		do u1=1,((N_0-la)/2)+1
!			print "(5F10.5)",(vls_field%p(0)%d3tensor(la)%d2(u1,u2),u2=1,((N_0-la)/2)+1)
!		end do
!	end do

		CALL SymGDDph_update(field1, diagonal%consistency%gDDph, diagonal%consistency%density%field%rho)
		CALL SymHartreeFockField_product(gdd_field, DBLE(1.0 / (b ** 4.0)), field1)
!	print *,"DD:"
!	do la=0,N_0
!		print *,"la=",la
!		do u1=1,((N_0-la)/2)+1
!			print "(5F10.5)",(gdd_field%p(0)%d3tensor(la)%d2(u1,u2),u2=1,((N_0-la)/2)+1)
!		end do
!	end do

		CALL SymHartreeFockField_add(field1, vls_field, gdd_field)
		CALL SymHartreeFockField_add(field2, vc_field, field1)
		CALL SymHartreeFockField_add(field1, vbb_field, field2)
		CALL SymHartreeFockField_add(field2, ekcm_field, field1)
		CALL SymHartreeFockField_add_SymD3Tensor(HF_Gamma, ek_tensor, field2)

!	print *,"Gamma:"
!	do la=0,N_0
!		print *,"la=",la
!		do u1=1,((N_0-la)/2)+1
!			print "(5F10.5)",(HF_Gamma%p(0)%d3tensor(la)%d2(u1,u2),u2=1,((N_0-la)/2)+1)
!		end do
!	end do
!		HF_Gamma = (factor / (b * b)) * EkField + &
!			(1.0 / (A * (b * b))) * (diagonal%consistency%vEkCMph * diagonal%consistency%density%field%rho) + &
!			(diagonal%consistency%vBBph * diagonal%consistency%density%field%rho) + &
!			(1.0 / b) * (diagonal%consistency%vCph * diagonal%consistency%density%field%rho) + &
!			(1.0 / (b ** 5.0)) * (diagonal%consistency%vLSph * diagonal%consistency%density%field%rho) + &
!			(1.0 / (b ** 4.0)) * SymGDDph_update(diagonal%consistency%gDDph, diagonal%consistency%density%field%rho)

		CALL SymVBBpp_product(vbb_field, diagonal%consistency%vBBpp, diagonal%consistency%density%field%kap)

		CALL SymVCpp_product(field1, diagonal%consistency%vCpp, diagonal%consistency%density%field%kap)
		CALL SymHartreeFockField_product(vc_field, DBLE(1.0 / b), field1)

		CALL SymVLSpp_product(field1, diagonal%consistency%vLSpp, diagonal%consistency%density%field%kap)
		CALL SymHartreeFockField_product(vls_field, DBLE(1.0 / (b ** 5.0)), field1)

		CALL SymHartreeFockField_add(field1, vc_field, vls_field)
		CALL SymHartreeFockField_add(HF_Delta, vbb_field, field1)

!		HF_Delta = (diagonal%consistency%vBBpp * diagonal%consistency%density%field%kap) + &
!			(1.0 / b) * (diagonal%consistency%vCpp * diagonal%consistency%density%field%kap) + &
!			(1.0 / (b ** 5.0)) * (diagonal%consistency%vLSpp * diagonal%consistency%density%field%kap)

		CALL SymGenDensityHF_new(gendenhf_gamma)
		CALL SymGenDensityHF_new(gendenhf_delta)
		CALL SymGenDensityHF_copy(gendenhf_gamma, HF_Gamma)
		CALL SymGenDensityHF_copy(gendenhf_delta, HF_Delta)

		CALL SymGenDensity_new_GammaDelta(S, gendenhf_gamma, gendenhf_delta, b)

		CALL SymGenDensityHF_del(gendenhf_gamma)
		CALL SymGenDensityHF_del(gendenhf_delta)

		CALL SymHartreeFockField_del(HF_Gamma)
		CALL SymHartreeFockField_del(HF_Delta)
		CALL SymD3Tensor_del(ek_tensor)
		CALL SymHartreeFockField_del(ekcm_field)
		CALL SymHartreeFockField_del(vbb_field)
		CALL SymHartreeFockField_del(vc_field)
		CALL SymHartreeFockField_del(vls_field)
		CALL SymHartreeFockField_del(gdd_field)
		CALL SymHartreeFockField_del(field1)
		CALL SymHartreeFockField_del(field2)
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

		INTEGER N, Z, A
		DOUBLE PRECISION factor, b
		DOUBLE PRECISION, DIMENSION(0:1) :: np, R2
		DATA np /0.0, 0.0/
		DATA R2 /0.0, 0.0/
		TYPE (SymGenDensity) SuperHamiltonian
		INTEGER ta, a2, d, dd, la, ja

		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: I, h, Delta, SH, UV, U, V, tmp
		DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Diag2
		INTEGER u1, u2
		DOUBLE PRECISION trace

		TYPE (SymHartreeFockField) HF_Gamma, HF_Delta
		TYPE (SymD3Tensor) ek_tensor
		TYPE (SymHartreeFockField) ekcm_field, vbb_field, vc_field, vls_field, gdd_field, field1, field2
		TYPE (SymGenDensityHF) gendenhf_gamma, gendenhf_delta

		N = Nucleus_get_N(diagonal%consistency%density%nucleus)
		Z = Nucleus_get_Z(diagonal%consistency%density%nucleus)
		A = N + Z
		factor = 1.0 - (1.0 / A)
		b = Nucleus_get_b(diagonal%consistency%density%nucleus)

		! ATENCIÓN: El siguiente procedimiento crea un núcleo a partir
		! del número de neutrones y protones (N y Z). Núcleo que es sobreescrito
		! a continuación... Esto ha sido mantenido por "similitud" con el
		! programa original, pero debería ser corregido en una versión posterior
		CALL WaveFunction_new(DiagonalizationMethod_get_WaveFunction, N, Z)
		CALL Nucleus_copy(DiagonalizationMethod_get_WaveFunction%nucleus, diagonal%consistency%density%nucleus)

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

		CALL SymD3Tensor_product(ek_tensor, DBLE(factor / (b * b)), EkField)

		CALL SymKineticEnergy2Body_product(field1, diagonal%consistency%vEkCMph, diagonal%consistency%density%field%rho)
		CALL SymHartreeFockField_product(ekcm_field, DBLE(1.0 / (A * (b * b))), field1)

		CALL SymVBBph_product(vbb_field, diagonal%consistency%vBBph, diagonal%consistency%density%field%rho)

		CALL SymVCph_product(field1, diagonal%consistency%vCph, diagonal%consistency%density%field%rho)
		CALL SymHartreeFockField_product(vc_field, DBLE(1.0 / b), field1)

		CALL SymVLSph_product(field1, diagonal%consistency%vLSph, diagonal%consistency%density%field%rho)
		CALL SymHartreeFockField_product(vls_field, DBLE(1.0 / (b ** 5.0)), field1)

		CALL SymGDDph_update(field1, diagonal%consistency%gDDph, diagonal%consistency%density%field%rho)
		CALL SymHartreeFockField_product(gdd_field, DBLE(1.0 / (b ** 4.0)), field1)

		CALL SymHartreeFockField_add(field1, vls_field, gdd_field)
		CALL SymHartreeFockField_add(field2, vc_field, field1)
		CALL SymHartreeFockField_add(field1, vbb_field, field2)
		CALL SymHartreeFockField_add(field2, ekcm_field, field1)
		CALL SymHartreeFockField_add_SymD3Tensor(HF_Gamma, ek_tensor, field2)

		CALL SymVBBpp_product(vbb_field, diagonal%consistency%vBBpp, diagonal%consistency%density%field%kap)

		CALL SymVCpp_product(field1, diagonal%consistency%vCpp, diagonal%consistency%density%field%kap)
		CALL SymHartreeFockField_product(vc_field, DBLE(1.0 / b), field1)

		CALL SymVLSpp_product(field1, diagonal%consistency%vLSpp, diagonal%consistency%density%field%kap)
		CALL SymHartreeFockField_product(vls_field, DBLE(1.0 / (b ** 5.0)), field1)

		CALL SymHartreeFockField_add(field1, vc_field, vls_field)
		CALL SymHartreeFockField_add(HF_Delta, vbb_field, field1)

!		CALL SymGenDensity_new_GammaDelta(SuperHamiltonian, &
!			(factor / (b * b)) * EkField &
!			+ (1.0 / (A * (b * b))) * (diagonal%consistency%vEkCMph * diagonal%consistency%density%field%rho)  &
!			+ (diagonal%consistency%vBBph * diagonal%consistency%density%field%rho)  &
!			+ (1.0 / b) * (diagonal%consistency%vCph * diagonal%consistency%density%field%rho)  &
!			+ (1.0 / (b ** 5.0)) * (diagonal%consistency%vLSph * diagonal%consistency%density%field%rho)  &
!			+ (1.0 / (b ** 4.0)) * SymGDDph_update(diagonal%consistency%gDDph, diagonal%consistency%density%field%rho), & ! pairing field
!			  (diagonal%consistency%vBBpp * diagonal%consistency%density%field%kap)  &
!			+ (1.0 / b) * (diagonal%consistency%vCpp * diagonal%consistency%density%field%kap)  &
!			+ (1.0 / (b ** 5.0)) * (diagonal%consistency%vLSpp * diagonal%consistency%density%field%kap), &
!			b)

		CALL SymGenDensityHF_new(gendenhf_gamma)
		CALL SymGenDensityHF_new(gendenhf_delta)
		CALL SymGenDensityHF_copy(gendenhf_gamma, HF_Gamma)
		CALL SymGenDensityHF_copy(gendenhf_delta, HF_Delta)

		CALL SymGenDensity_new_GammaDelta(SuperHamiltonian, gendenhf_gamma, gendenhf_delta, b)

		CALL SymGenDensityHF_del(gendenhf_gamma)
		CALL SymGenDensityHF_del(gendenhf_delta)

		CALL SymHartreeFockField_del(HF_Gamma)
		CALL SymHartreeFockField_del(HF_Delta)
		CALL SymD3Tensor_del(ek_tensor)
		CALL SymHartreeFockField_del(ekcm_field)
		CALL SymHartreeFockField_del(vbb_field)
		CALL SymHartreeFockField_del(vc_field)
		CALL SymHartreeFockField_del(vls_field)
		CALL SymHartreeFockField_del(gdd_field)
		CALL SymHartreeFockField_del(field1)
		CALL SymHartreeFockField_del(field2)

		DO ta = 0, 1
			DO a2 = 0, 2 * N_0
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

				I = 0.0
				DO u1 = 1, d
					I(u1, u1) = 1.0
				END DO

				h = SuperHamiltonian%rho%rho(ta, a2)%store &
					- (diagonal%consistency%density%nucleus%lambda_np(ta) * I) &
					- (diagonal%consistency%density%nucleus%lambda_R2(ta) * (b ** 2) &
					* SymD3Tensor_matrix(R2Field, la))

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

				CALL EigenValues(dd, dd, SH, Diag2, UV)

				! Extraemos las matrices U y V
				! A partir de la matriz UV, la posición de las matrices U y V
				! es la siguiente:
				!   X X
				!   U V
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
					trace = trace + tmp(u1, u1)
				END DO
				np(ta) = np(ta) + (DBLE(ja + 1.0) * trace)

				tmp = MATMUL(SymD3Tensor_matrix(R2Field, la), MATMUL(V, TRANSPOSE(V)))
				! Calculamos la traza de la matriz "tmp"
				trace = 0.0
				DO u1 = 1, d
					trace = trace + tmp(u1, u1)
				END DO
				R2(ta) = R2(ta) + (DBLE(ja + 1.0) * trace)

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
		RETURN
	END FUNCTION DiagonalizationMethod_get_WaveFunction

	FUNCTION DiagonalizationMethod_operator(diagonal)
		DOUBLE PRECISION DiagonalizationMethod_operator
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal

		DOUBLE PRECISION np, R2, b2
		INTEGER a, d, dd, la, ja

		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: xI, h, Delta, SH, U, V, tmp
		INTEGER u1, u2
		DOUBLE PRECISION trace

		np = 0.0
		R2 = 0.0
		b2 = Nucleus_get_b(diagonal%iterated%nucleus) ** 2

		DO a = 0, 2 * N_0
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

			xI = 0.0
			DO u1 = 1, d
				xI(u1, u1) = diagonal%func%x
			END DO
			h = diagonal%S%rho%rho(diagonal%ta, a)%store - xI

			Delta = 1.0 * diagonal%S%kap%rho(diagonal%ta, a)%store

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

			! ATENCION: Este procedimiento no ha sido traducido, sino
			! extraído a traves de Google.
			! Confiamos en su correcta implementacion

			! Da las enegias ordenadas de menor (negativas) a mayor (positivas)
			CALL EigenValues(dd, dd, SH, &
				diagonal%QuasiParticleEnergies(diagonal%ta, a)%value, &
				diagonal%UV(diagonal%ta, a)%value)

			! Extraemos las matrices U y V
			DO u1 = 1, d
				DO u2 = 1, d
					U(u1, u2) = diagonal%UV(diagonal%ta, a)%value(u1    , u2 + d)
					V(u1, u2) = diagonal%UV(diagonal%ta, a)%value(u1 + d, u2 + d)
				END DO
			END DO

			CALL SymGenDensity_make_Block(diagonal%iterated, diagonal%ta, a, U, V)

			trace = 0.0
			DO u1 = 1, d
				trace = trace + diagonal%iterated%rho%rho(diagonal%ta, a)%store(u1, u1)
			END DO
			np = np + DBLE(ja + 1.0) * trace

			tmp = MATMUL(SymD3Tensor_matrix(R2Field, la), diagonal%iterated%rho%rho(diagonal%ta, a)%store)
			trace = 0.0
			DO u1 = 1, d
				trace = trace + tmp(u1, u1)
			END DO
			R2 = R2 + DBLE(ja + 1.0) * trace

			! Liberamos la memoria de las matrices utilizadas
			DEALLOCATE(xI)
			DEALLOCATE(h)
			DEALLOCATE(Delta)
			DEALLOCATE(SH)
			DEALLOCATE(U)
			DEALLOCATE(V)
			DEALLOCATE(tmp)
		END DO
		diagonal%iterated%nucleus%actual_np(diagonal%ta) = np
		diagonal%iterated%nucleus%actual_R2(diagonal%ta) = (b2 * R2) / np
		! Resultado final
		DiagonalizationMethod_operator = np
		RETURN
	END FUNCTION DiagonalizationMethod_operator

	! Encontrar la base de particulas
	SUBROUTINE DiagonalizationMethod_show_ParticleEnergies(diagonal)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal

		DOUBLE PRECISION b
		TYPE (SymGenDensity) G, SuperHamiltonian
		INTEGER ndim, a, num, i, dim_cur, dim_acc, la, ja, ta
		INTEGER, DIMENSION(:), ALLOCATABLE :: an, ap
		TYPE (DiagonalMatrix), DIMENSION(0:1) :: D
		TYPE (MatrixType), DIMENSION(0:1) :: P, E, R2
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: S, V
		TYPE (DiagonalMatrix), DIMENSION(0:1) :: Occup, SingleParticleE, SingleParticleP, SingleParticleR2
		INTEGER, DIMENSION(:), ALLOCATABLE :: in, ip
		DOUBLE PRECISION spn, spp
		CHARACTER, DIMENSION(0:19) :: spectr
		DATA spectr / "s", "p", "d", "f", "g", "h", "i", "j", "k", "l", &
			"m", "n", "o", "p", "q", "r", "s", "t", "u", "v" /
		INTEGER nrot

		b = Nucleus_get_b(diagonal%consistency%density%nucleus)

		CALL SymGenDensity_new_SymDensity(G, diagonal%consistency%density)
!A	GenDenHF PB = G.rho.get_particle_basis();

		dim_cur = 0
		DO a = 0, 2 * N_0
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

		CALL DiagonalizationMethod_get_MeanField(diagonal, SuperHamiltonian)
!A	SuperHamiltonian.rho.show_particle_energies(PB);

		ndim = (N_0 + 1) * (N_0 + 2) / 2
		ALLOCATE(an(ndim + 1))
		ALLOCATE(ap(ndim + 1))

		dim_acc = 0
		num = 1
		DO a = 0, 2 * N_0
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

			DO ta = 0, 1
				ALLOCATE(D(ta)%value(dim_cur))
				ALLOCATE(E(ta)%value(dim_cur, dim_cur))
				ALLOCATE(P(ta)%value(dim_cur, dim_cur))
				ALLOCATE(R2(ta)%value(dim_cur, dim_cur))

				S = G%rho%rho(ta, a)%store
				CALL Jacobi_real8(S, dim_cur, D(ta)%value, V, nrot)

				E(ta)%value = MATMUL(TRANSPOSE(V), MATMUL(SuperHamiltonian%rho%rho(ta, a)%store, V))
				P(ta)%value = MATMUL(TRANSPOSE(V), MATMUL(SuperHamiltonian%kap%rho(ta, a)%store, V))
				R2(ta)%value = MATMUL(TRANSPOSE(V), (b * b) * MATMUL(SymD3Tensor_matrix(R2Field, L(a)), V))

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
			END DO

			DEALLOCATE(S)
			DEALLOCATE(V)

			dim_acc = dim_acc + dim_cur
		END DO

		ALLOCATE(in(ndim))
		ALLOCATE(ip(ndim))
		CALL indexx_real8(ndim, SingleParticleE(1)%value, in)
		CALL indexx_real8(ndim, SingleParticleE(0)%value, ip)
		PRINT *
		PRINT "(A32,A36)", "neu", "pro"
spe:		DO i = 1, ndim
			spn = SingleParticleE(1)%value(in(i))
			spp = SingleParticleE(0)%value(ip(i))
			IF ((spn .LT. -20.0) .AND. (spp .LT. -20.0)) CYCLE
			IF ((spn .GT.  10.0) .AND. (spp .GT.  10.0)) EXIT spe

			PRINT "(A2,I2,F6.2,F10.3,F8.3,F8.3,A8,I2,F6.2,F10.3,F10.3,F8.3)", &
				spectr(L(an(in(i)))), J(an(in(i))), Occup(1)%value(in(i)), spn, SingleParticleP(1)%value(in(i)), SingleParticleR2(1)%value(in(i)), &
				spectr(L(ap(ip(i)))), J(ap(ip(i))), Occup(0)%value(ip(i)), spp, SingleParticleP(0)%value(ip(i)), SingleParticleR2(0)%value(ip(i))
		END DO spe

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

	SUBROUTINE DiagonalizationMethod_show_QuasiParticleEnergies(diagonal)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal

		INTEGER NN, NumberOfStates, a, na, top
		INTEGER, DIMENSION(:), ALLOCATABLE :: nIndx, pIndx
		CHARACTER, DIMENSION(0:19) :: spectr
		DATA spectr / "s", "p", "d", "f", "g", "h", "i", "j", "k", "l", &
			"m", "n", "o", "p", "q", "r", "s", "t", "u", "v" /
		DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: nQPE, pQPE, nVV, pVV
		INTEGER, DIMENSION(:), ALLOCATABLE :: nMomentum, pMomentum
		INTEGER num, la, ja, d, s1, sa
		DOUBLE PRECISION sumvvn, sumvvp
		INTEGER an, ap, jn, jp
		CHARACTER cn, cp

		NN = N_0 * 2
		NumberOfStates = 0;
		DO a = 0, NN
			NumberOfStates = NumberOfStates + DIM(a)
		END DO

		! Imprimir las energias de cuasiparticulas
		! y las ocupaciones para neutrones y protones
		top = MIN(100, NumberOfStates)

		ALLOCATE(nIndx(top))
		ALLOCATE(pIndx(top))

		! Asignar memoria a la tabla de energias de cuasiparticulas,
		! sus momentos angulares y sus ocupaciones
		ALLOCATE(nQPE(top))
		ALLOCATE(pQPE(top))
		ALLOCATE(nVV(top))
		ALLOCATE(pVV(top))
		ALLOCATE(nMomentum(top))
		ALLOCATE(pMomentum(top))

		! Alamcenar las energias de cuasiparticulas y sus ocupaciones
		num = 1
		DO a = 0, NN
			la = L(a)
			ja = J(a)
			d = ((N_0 - la) / 2) + 1
			DO s1 = d + 1, d + d
				num = num + 1
				sumvvn = 0.0
				sumvvp = 0.0
				IF (num > top) GOTO 1 ! Fin 
				nQPE(num) = diagonal%QuasiParticleEnergies(NEUTRON, a)%value(s1)
				pQPE(num) = diagonal%QuasiParticleEnergies(PROTON, a)%value(s1)
				nMomentum(num) = a
				pMomentum(num) = a
				nIndx(num) = num
				pIndx(num) = num
				DO sa = 1, d
					sumvvn = sumvvn + (diagonal%UV(NEUTRON, a)%value(sa + d, s1) ** 2)
					sumvvp = sumvvp + (diagonal%UV(PROTON, a)%value(sa + d, s1) ** 2)
				END DO

				IF ((diagonal%consistency%density%nucleus%is_blocking(0)) .AND. &
					(a .EQ. diagonal%consistency%density%nucleus%ia(0)) .AND. &
					(s1 .EQ. (d + diagonal%consistency%density%nucleus%mu0(0)))) THEN
					PRINT *, "Se pone la corrección de bloqueo de protones"
					DO sa = 1, d
						sumvvp = sumvvp + ((diagonal%UV(PROTON, a)%value(sa, s1) ** 2) &
							- (diagonal%UV(PROTON, a)%value(sa + d, s1) ** 2)) / (ja + 1.0)
					END DO
				END IF
!A				IF ( Nbl && a == ibln && s1 == mu0n ) 
!A					DO sa = 1, d
!A						sumvvn = sumvvn + ((diagonal%UV(NEUTRON, a)%value(sa, s1) ** 2) - (diagonal%UV(NEUTRON, a)%value(sa + d, s1) ** 2)) / (ja + 1.0)
				nVV(num) = sumvvn
				pVV(num) = sumvvp
			END DO
		END DO

		! Ordenarlas en orden creciente reordenando
		! los indices de momento angular
1		CALL indexx_real8(top, nVV, nIndx)
		CALL indexx_real8(top, pVV, pIndx)
		PRINT "(A14,A10,A12,A5,A10)", "EQP NEUT", "V2", "", "V2", "EQP PROT"
		DO num = top, 1, -1
			an = nMomentum(nIndx(num))
			ap = pMomentum(pIndx(num))
			cn = spectr((an + 1) / 2)                       
			cp = spectr((ap + 1) / 2)
			jn = J(an)
			jp = J(ap)
			PRINT "(F14.4,F10.4,A1,I2,A5,I2,F6.4,F10.4)", &
				nQPE(nIndx(num)), nVV(nIndx(num)), cn, jn, cp, jp, pVV(pIndx(num)), pQPE(pIndx(num))
		END DO
		PRINT *, "Lineas = ", top

		DEALLOCATE(nIndx)
		DEALLOCATE(pIndx)
		DEALLOCATE(nQPE)
		DEALLOCATE(pQPE)
		DEALLOCATE(nVV)
		DEALLOCATE(pVV)
		DEALLOCATE(nMomentum)
		DEALLOCATE(pMomentum)
		RETURN
	END SUBROUTINE DiagonalizationMethod_show_QuasiParticleEnergies

	SUBROUTINE DiagonalizationMethod_del(diagonal)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal

		INTEGER ta, a, max_a

		max_a = 2 * N_0
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

	SUBROUTINE WaveFunction_new_SymDensity(wave_func, density, tolerance)
		TYPE (WaveFunction), INTENT(INOUT) :: wave_func
		TYPE (SymDensity), INTENT(INOUT) :: density
		DOUBLE PRECISION, INTENT(IN) :: tolerance

		TYPE (DiagonalizationMethod) evolv

		CALL Nucleus_new_Nucleus(wave_func%nucleus, density%nucleus)
		CALL DiagonalizationMethod_new(evolv, density)
		CALL DiagonalizationMethod_goto_SelfConsistency(evolv, tolerance)
		wave_func = DiagonalizationMethod_get_WaveFunction(evolv)
		RETURN
	END SUBROUTINE WaveFunction_new_SymDensity

END MODULE diagmeth
