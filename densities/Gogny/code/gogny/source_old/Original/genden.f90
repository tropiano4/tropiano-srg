MODULE genden

	USE input
	USE math
!	USE matrix

	IMPLICIT NONE

	TYPE Matrix2D
		DOUBLE PRECISION, DIMENSION(:, :), POINTER :: value
	END TYPE

	TYPE GenDensityHF
		TYPE (Matrix2D), DIMENSION(:, :), POINTER :: rho
	END TYPE

	TYPE GenDensity
		TYPE (GenDensityHF) rho, kap
	END TYPE

CONTAINS

	SUBROUTINE GenDensity_new(density)
		TYPE (GenDensity), INTENT(INOUT) :: density

		CALL GenDensityHF_new(density%rho)
		CALL GenDensityHF_new(density%kap)
		RETURN
	END SUBROUTINE GenDensity_new

	SUBROUTINE GenDensityHF_new(hf)
		TYPE (GenDensityHF), INTENT(INOUT) :: hf

		INTEGER ta, a, d

		ALLOCATE(hf%rho(0:1, 0:(2 * N_0)))
		DO ta = 0, 1
			DO a = 0, 2 * N_0
				d = DIM(a)
				ALLOCATE(hf%rho(ta, a)%value(d, d))
				! Comprobamos que se ha reservado correctamente la memoria
				IF (.NOT. ASSOCIATED(hf%rho(ta, a)%value)) STOP "Unable to allocate memory"
			END DO
		END DO
		RETURN
	END SUBROUTINE GenDensityHF_new

	SUBROUTINE GenDensityHF_del(hf)
		TYPE (GenDensityHF), INTENT(INOUT) :: hf

		INTEGER ta, a

		DO ta = 0, 1
			DO a = 0, 2 * N_0
				DEALLOCATE(hf%rho(ta, a)%value)
			END DO
		END DO
		DEALLOCATE(hf%rho)
		RETURN
	END SUBROUTINE GenDensityHF_del

END MODULE genden
