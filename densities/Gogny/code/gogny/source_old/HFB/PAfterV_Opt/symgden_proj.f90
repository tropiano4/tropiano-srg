 MODULE symgden_proj

	USE input
	USE symd3t_proj
	USE nucleus
	USE symfield_proj
	USE symden_proj
	USE symgdhf_proj

	IMPLICIT NONE

	!---------------------------------------------------------------!
	!    SymGenDensity is a type for objects associated to:		!
	!	- a nucleus 						!
	!       - mean-field densities rho_{na,nc,la}			!
	!       - pairing abnormal tensor kappa_{na,nc,la}		!
	!---------------------------------------------------------------!

	TYPE SymGenDensityProj
		TYPE (NucleusType) nucleus
		TYPE (SymGenDensityHFProj) :: rho, kap
	END TYPE

 CONTAINS

	SUBROUTINE SymGenDensityProj_new(genden)
		TYPE (SymGenDensityProj), INTENT(INOUT) :: genden
	
		CALL SymGenDensityHFProj_new(genden%rho)
		CALL SymGenDensityHFProj_new(genden%kap)
		
		RETURN
	END SUBROUTINE SymGenDensityProj_new

	SUBROUTINE SymGenDensityProj_new_Nucleus(genden, N, Z)
		TYPE (SymGenDensityProj), INTENT(INOUT) :: genden
		INTEGER, INTENT(IN) :: N, Z
		DOUBLE PRECISION :: b
	
		CALL Nucleus_new(genden%nucleus, N, Z, b)
		
		CALL SymGenDensityHFProj_new(genden%rho)
		CALL SymGenDensityHFProj_new(genden%kap)
		
		RETURN
	END SUBROUTINE SymGenDensityProj_new_Nucleus

	SUBROUTINE SymGenDensityProj_new_SymDensity(genden, density)
		TYPE (SymGenDensityProj), INTENT(INOUT) :: genden
		TYPE (SymDensityProj), INTENT(INOUT) :: density

		INTEGER ta, la, lla

		CALL Nucleus_new_Nucleus(genden%nucleus, density%nucleus)
		
		! Creates projected densities rho and kappa (each depends on the Gauge angle phi)
		
		CALL SymGenDensityHFProj_new(genden%rho)
		CALL SymGenDensityHFProj_new(genden%kap)
		
		DO ta = 0, 1
			DO la = 0, Lmax
			
				lla = 2 * (2 * la + 1)

				genden%rho%rho(ta, 2 * la)%store = &
					(SymD3Tensor_matrix(density%field%rho%p(ta), la) * DBLE(2.0 * la) * &
					 SymD3Tensor_matrix(density%field%rho%a(ta), la)) / lla
					 
				genden%kap%rho(ta, 2 * la)%store = &
					(SymD3Tensor_matrix(density%field%kap%p(ta), la) * DBLE(2.0 * la) * &
					 SymD3Tensor_matrix(density%field%kap%a(ta), la)) / lla
					 
				IF (HFOnly .EQ. 1) genden%kap%rho(ta, 2 * la)%store = 0.0

				IF (la .EQ. 0) CYCLE

				genden%rho%rho(ta, (2 * la) - 1)%store = &
					(SymD3Tensor_matrix(density%field%rho%p(ta), la) * DBLE(2.0 * (la + 1)) * &
					 SymD3Tensor_matrix(density%field%rho%a(ta), la)) / lla
					 
				genden%kap%rho(ta, (2 * la) - 1)%store = &
					(SymD3Tensor_matrix(density%field%kap%p(ta), la) * DBLE(2.0 * (la + 1)) * &
					 SymD3Tensor_matrix(density%field%kap%a(ta), la)) / lla
					 
				IF (HFOnly .EQ. 1) genden%kap%rho(ta, (2 * la) - 1)%store = 0.0
				
			END DO
		END DO
		
		RETURN
	END SUBROUTINE SymGenDensityProj_new_SymDensity

	! Create a new density of the type SymDensityProj from an
	! old SymGenDensityProj
	!
	SUBROUTINE SymDensity_new_GenDensityProj(density, genden)
		TYPE (SymDensityProj), INTENT(INOUT) :: density
		TYPE (SymGenDensityProj), INTENT(IN) :: genden

		CALL Nucleus_new_Nucleus(density%nucleus, genden%nucleus)
		
		CALL SymHartreeFockBogolFieldProj_new(density%field)
		
		density%field%rho = genden%rho
		density%field%kap = genden%kap			
		
		RETURN
	END SUBROUTINE SymDensity_new_GenDensityProj

	SUBROUTINE SymGenDensityProj_new_GammaDelta(genden, HF_gamma, HF_delta, b)
		TYPE (SymGenDensityProj), INTENT(INOUT) :: genden
		TYPE (SymGenDensityHFProj), INTENT(IN) :: HF_gamma, HF_delta

		DOUBLE PRECISION, INTENT(IN) :: b

		CALL Nucleus_new(genden%nucleus, 8, 8, b)
		
		CALL SymGenDensityHFProj_new(genden%rho)
		CALL SymGenDensityHFProj_new(genden%kap)
		
		genden%rho = HF_gamma
		genden%kap = HF_delta
		
		RETURN
	END SUBROUTINE SymGenDensityProj_new_GammaDelta

	SUBROUTINE SymGenDensityProj_setGauge(genden, Gauge, ta)
		TYPE (SymGenDensityProj), INTENT(OUT) :: genden
		DOUBLE PRECISION, INTENT(IN) :: Gauge
		INTEGER, INTENT(IN) :: ta

		genden%rho%GaugeAngle(ta) = Gauge
		genden%kap%GaugeAngle(ta) = Gauge			
		
		RETURN
	END SUBROUTINE SymGenDensityProj_setGauge

	!---------------------------------------------------------------!
	!								!
	!		ULTRA-IMPORTANT SUBROUTINE			!
	!		--------------------------			!
	!								!
	!    This subroutine calculates the values of the densities	!
	!    rho and kappa from the U and V vectors of the Bogoliubov	!
	!    transformation. We have:					!
	!								!
	!  		rho = (V*)(VT)					!
	! 		kappa = U(V+)					!
	!								!
	!---------------------------------------------------------------!

	SUBROUTINE SymGenDensityProj_make_Block(genden, ta, a, U, V)
		TYPE (SymGenDensityProj), INTENT(INOUT) :: genden
		DOUBLE PRECISION, DIMENSION(:, :), INTENT(INOUT) :: U, V !TODO(INOUT?)
		INTEGER, INTENT(IN) :: ta, a

		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: S0, S1, S2, S3
		DOUBLE PRECISION :: U0, V0, Gauge
		
		INTEGER :: d, i
		
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
			S0 = S0 + (1.0 / (genden%nucleus%ja(ta) + 1.0)) * (S1 - S0)
			S2 = S2 + (1.0 / (genden%nucleus%ja(ta) + 1.0)) * (S3 - S2)

			DEALLOCATE(S1)
			DEALLOCATE(S3)
			
		END IF		
		
		genden%rho%rho(ta, a)%store = S0
		genden%kap%rho(ta, a)%store = S2
		
		DEALLOCATE(S0)
		DEALLOCATE(S2)
		
		RETURN
	END SUBROUTINE SymGenDensityProj_make_Block

	SUBROUTINE SymGenDensityProj_del(genden)
		TYPE (SymGenDensityProj), INTENT(INOUT) :: genden

		CALL Nucleus_del(genden%nucleus)
		
		CALL SymGenDensityHFProj_del(genden%rho)
		CALL SymGenDensityHFProj_del(genden%kap)
		
		RETURN
		
	END SUBROUTINE SymGenDensityProj_del

END MODULE symgden_proj
