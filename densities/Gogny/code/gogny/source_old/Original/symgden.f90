MODULE symgden

	USE input
	USE symd3t
	USE nucleus
	USE symfield
	USE symden
	USE symgdhf

	IMPLICIT NONE

	TYPE SymGenDensity
		TYPE (NucleusType) nucleus
		TYPE (SymGenDensityHF) rho, kap
	END TYPE

CONTAINS

	SUBROUTINE SymGenDensity_new(genden)
		TYPE (SymGenDensity), INTENT(INOUT) :: genden

		CALL SymGenDensityHF_new(genden%rho)
		CALL SymGenDensityHF_new(genden%kap)
		RETURN
	END SUBROUTINE SymGenDensity_new

	SUBROUTINE SymGenDensity_new_Nucleus(genden, N, Z)
		TYPE (SymGenDensity), INTENT(INOUT) :: genden
		INTEGER, INTENT(IN) :: N, Z

		CALL Nucleus_new(genden%nucleus, N, Z)
		CALL SymGenDensityHF_new(genden%rho)
		CALL SymGenDensityHF_new(genden%kap)
		RETURN
	END SUBROUTINE SymGenDensity_new_Nucleus

	SUBROUTINE SymGenDensity_new_SymDensity(genden, density)
		TYPE (SymGenDensity), INTENT(INOUT) :: genden
		TYPE (SymDensity), INTENT(INOUT) :: density

		INTEGER ta, la, lla
		DOUBLE PRECISION, DIMENSION(1:((N_0 / 2) + 1), 1:((N_0 / 2) + 1)) :: tmp

		CALL Nucleus_new_Nucleus(genden%nucleus, density%nucleus)
		CALL SymGenDensityHF_new(genden%rho)
		CALL SymGenDensityHF_new(genden%kap)
		DO ta = 0, 1
			DO la = 0, N_0
				lla = 2 * (2 * la + 1)

				genden%rho%rho(ta, 2 * la)%store = &
					(SymD3Tensor_matrix(density%field%rho%p(ta), la) * DBLE(2.0 * la) * &
					 SymD3Tensor_matrix(density%field%rho%a(ta), la)) / lla
				genden%kap%rho(ta, 2 * la)%store = &
					(SymD3Tensor_matrix(density%field%kap%p(ta), la) * DBLE(2.0 * la) * &
					 SymD3Tensor_matrix(density%field%kap%a(ta), la)) / lla

				IF (la .EQ. 0) CYCLE

				genden%rho%rho(ta, (2 * la) - 1)%store = &
					(SymD3Tensor_matrix(density%field%rho%p(ta), la) * DBLE(2.0 * (la + 1)) * &
					 SymD3Tensor_matrix(density%field%rho%a(ta), la)) / lla
				genden%kap%rho(ta, (2 * la) - 1)%store = &
					(SymD3Tensor_matrix(density%field%kap%p(ta), la) * DBLE(2.0 * (la + 1)) * &
					 SymD3Tensor_matrix(density%field%kap%a(ta), la)) / lla
			END DO
		END DO
		RETURN
	END SUBROUTINE SymGenDensity_new_SymDensity

	SUBROUTINE SymDensity_new_GenDensity(density, genden)
		! Recibe como parámetros de entrada el registro de densidad
		! y el número de protones y electrones
		TYPE (SymDensity), INTENT(INOUT) :: density
		TYPE (SymGenDensity), INTENT(IN) :: genden

		CALL Nucleus_new_Nucleus(density%nucleus, genden%nucleus)
		CALL SymHartreeFockBogolField_new(density%field)
		density%field%rho = genden%rho
		density%field%kap = genden%kap
		RETURN
	END SUBROUTINE SymDensity_new_GenDensity

	SUBROUTINE SymGenDensity_new_GammaDelta(genden, HF_gamma, HF_delta, b)
		TYPE (SymGenDensity), INTENT(INOUT) :: genden
		TYPE (SymGenDensityHF), INTENT(IN) :: HF_gamma, HF_delta
!		TYPE (SymHartreeFockField), INTENT(IN) :: HF_gamma, HF_delta
		DOUBLE PRECISION, INTENT(IN) :: b

		CALL Nucleus_new(genden%nucleus, 8, 8, b)
		CALL SymGenDensityHF_new(genden%rho)
		CALL SymGenDensityHF_new(genden%kap)
		genden%rho = HF_gamma
		genden%kap = HF_delta
!		CALL SymGenDensityHF_copy(genden%rho, HF_gamma)
!		CALL SymGenDensityHF_copy(genden%kap, HF_delta)
		RETURN
	END SUBROUTINE SymGenDensity_new_GammaDelta

	SUBROUTINE SymGenDensity_make_Block(genden, ta, a, U, V)
		TYPE (SymGenDensity), INTENT(INOUT) :: genden
		INTEGER, INTENT(IN) :: ta, a
		DOUBLE PRECISION, DIMENSION(:, :), INTENT(INOUT) :: U, V !TODO(INOUT?)

		INTEGER d, i
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: S0, S1, S2, S3
		DOUBLE PRECISION U0, V0

		d = DIM(a)
		ALLOCATE(S0(d, d))
		ALLOCATE(S2(d, d))
		S0 = MATMUL(V, TRANSPOSE(V))
		S2 = MATMUL(V, TRANSPOSE(U))
		IF (genden%nucleus%is_blocking(ta) .AND. (a .EQ. genden%nucleus%ia(ta))) THEN
			ALLOCATE(S1(d, d))
			ALLOCATE(S3(d, d))

			DO i = 1, d
				U0 = U(i, genden%nucleus%mu0(ta)) ! Column???
				V0 = V(i, genden%nucleus%mu0(ta))
				U(i, genden%nucleus%mu0(ta)) =   V0
				V(i, genden%nucleus%mu0(ta)) = - U0
			END DO
			S1 = MATMUL(V, TRANSPOSE(V))
			S3 = MATMUL(V, TRANSPOSE(U))
			S0 = S0 + (DBLE(1.0 / (genden%nucleus%ja(ta) + 1.0)) * (S1 - S0))
			S2 = S2 + (DBLE(1.0 / (genden%nucleus%ja(ta) + 1.0)) * (S3 - S2))

			DEALLOCATE(S1)
			DEALLOCATE(S3)
		END IF
		genden%rho%rho(ta, a)%store = S0
		genden%kap%rho(ta, a)%store = S2
		DEALLOCATE(S0)
		DEALLOCATE(S2)
		RETURN
	END SUBROUTINE SymGenDensity_make_Block

	SUBROUTINE SymGenDensity_reduce(genden_out, genden_in, t_in)
		TYPE (SymGenDensity), INTENT(INOUT) :: genden_out
		TYPE (SymGenDensity), INTENT(IN) :: genden_in
		TYPE (SymD3Tensor), INTENT(IN) :: t_in

		DOUBLE PRECISION b
		TYPE (SymGenDensityHF) tmp

		b = Nucleus_get_b(genden_in%nucleus)
		CALL SymGenDensityHF_new(tmp)
		CALL SymGenDensityHF_reduce(tmp, genden_in%rho, t_in)
		CALL SymGenDensity_new_GammaDelta(genden_out, tmp, genden_in%kap, b)
		CALL SymGenDensityHF_del(tmp)
		RETURN
	END SUBROUTINE SymGenDensity_reduce

	SUBROUTINE SymGenDensity_del(genden)
		TYPE (SymGenDensity), INTENT(INOUT) :: genden

		CALL Nucleus_del(genden%nucleus)
		CALL SymGenDensityHF_del(genden%rho)
		CALL SymGenDensityHF_del(genden%kap)
		RETURN
	END SUBROUTINE SymGenDensity_del

END MODULE symgden
