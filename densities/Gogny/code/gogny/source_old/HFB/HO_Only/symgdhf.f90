 MODULE symgdhf

	USE input
	USE math
	USE symd3t
	USE symfield

	IMPLICIT NONE

	TYPE MatrixType2
		DOUBLE PRECISION, DIMENSION(:, :), POINTER :: store
	END TYPE

	!---------------------------------------------------------------------------------------!
	!    SymGenDensityHF is a type for objects associated to densities like rho_{na,nc,la}	!
	!---------------------------------------------------------------------------------------!

	TYPE SymGenDensityHF
		TYPE (MatrixType2), DIMENSION(:, :), POINTER :: rho
	END TYPE

	INTERFACE ASSIGNMENT(=)
		MODULE PROCEDURE SymGenDensityHF_assign1, SymGenDensityHF_assign2
	END INTERFACE

 CONTAINS

	SUBROUTINE SymGenDensityHF_new(gendenhf)
		TYPE (SymGenDensityHF), INTENT(INOUT) :: gendenhf

		INTEGER ta, a, d

		ALLOCATE(gendenhf%rho(0:1, 0:(2*N_0)))
		
		DO ta = 0, 1
			DO a = 0, 2*N_0
				d = DIM(a)
				ALLOCATE(gendenhf%rho(ta, a)%store(d, d))
			END DO
		END DO
		RETURN
	END SUBROUTINE SymGenDensityHF_new

	SUBROUTINE SymGenDensityHF_copy(gendenhf, HF)
		TYPE (SymGenDensityHF), INTENT(INOUT) :: gendenhf
		TYPE (SymHartreeFockField), INTENT(IN) :: HF

		INTEGER ta, a, la, d

		INTEGER u1, u2
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: HF_p, HF_a

		DO ta = 0, 1
			DO a = 0, 2*N_0
			
				la = L(a)
				d = DIM(a)

				ALLOCATE(HF_p(d, d), HF_a(d, d))

				HF_p = 0.0
				HF_a = 0.0
				
				DO u1 = 1, d
					DO u2 = 1, u1
						HF_p(u1, u2) = HF%p(ta)%d3tensor(la)%d2(u1, u2)
						HF_a(u1, u2) = HF%a(ta)%d3tensor(la)%d2(u1, u2)
						IF (u1 .NE. u2) THEN
							HF_p(u2, u1) = HF%p(ta)%d3tensor(la)%d2(u1, u2)
							HF_a(u2, u1) = HF%a(ta)%d3tensor(la)%d2(u1, u2)
						END IF
					END DO
				END DO
				gendenhf%rho(ta, a)%store = HF_p + (LS(a) * HF_a)

				DEALLOCATE(HF_p, HF_a)
			END DO
		END DO
		RETURN
	END SUBROUTINE SymGenDensityHF_copy

	SUBROUTINE SymGenDensityHF_assign1(gendenhf_out, gendenhf_in)
		TYPE (SymGenDensityHF), INTENT(INOUT) :: gendenhf_out
		TYPE (SymGenDensityHF), INTENT(IN) :: gendenhf_in

		INTEGER ta, a

		DO ta = 0, 1
			DO a = 0, 2*N_0
				gendenhf_out%rho(ta, a)%store = gendenhf_in%rho(ta, a)%store
			END DO
		END DO
		RETURN
	END SUBROUTINE SymGenDensityHF_assign1

	SUBROUTINE SymGenDensityHF_assign2(HF, gendenhf)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF
		TYPE (SymGenDensityHF), INTENT(IN) :: gendenhf

		INTEGER ta, la, d
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: m

		DO ta = 0, 1
			DO la = 0, N_0
				d = DIM(2*la)
				
				ALLOCATE(m(d, d))

				m = DBLE(2.0 * (la + 1)) * gendenhf%rho(ta, 2 * la)%store
				
				IF (la .NE. 0) THEN
					m = m + (DBLE(2.0 * la) * gendenhf%rho(ta, (2 * la) - 1)%store)
				END IF

				CALL SymD3Tensor_assign_Matrix(HF%p(ta), la, m)

				m = gendenhf%rho(ta, 2 * la)%store
				
				IF (la .NE. 0) THEN
					m = m - gendenhf%rho(ta, (2 * la) - 1)%store
				END IF

				CALL SymD3Tensor_assign_Matrix(HF%a(ta), la, m)

				DEALLOCATE(m)
			END DO
		END DO
		RETURN
	END SUBROUTINE SymGenDensityHF_assign2

	SUBROUTINE SymGenDensityHF_reduce(gendenhf_out, gendenhf_in, t_in)
		TYPE (SymGenDensityHF), INTENT(INOUT) :: gendenhf_out
		TYPE (SymGenDensityHF), INTENT(IN) :: gendenhf_in
		TYPE (SymD3Tensor), INTENT(IN) :: t_in

		INTEGER ta, a, d, la
!		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: S1, S2

		DO ta = 0, 1
			DO a = 0, 2*N_0
				d = DIM(a)
				la = L(a)
!				ALLOCATE(S1(d, d))
!				ALLOCATE(S2(d, d))
!				S1 = gendenhf_in%rho(ta, a)%store
!				S2 = t_in%d3tensor(la)%d2
!				SymGenDensityHF_assign(gendenhf_out, ta, a, S1 - S2)
!				DEALLOCATE(S1)
!				DEALLOCATE(S2)
				gendenhf_out%rho(ta, a)%store = gendenhf_in%rho(ta, a)%store - t_in%d3tensor(la)%d2
			END DO
		END DO
		RETURN
	END SUBROUTINE SymGenDensityHF_reduce

	SUBROUTINE SymGenDensityHF_print(gendenhf)
		TYPE (SymGenDensityHF), INTENT(IN) :: gendenhf

		INTEGER ta, a, d, u1, u2

		DO ta = 0, 1
			DO a = 0, 2*N_0
!TODO				PRINT "()"
				d = DIM(a)
				DO u1 = 1, d
					PRINT "(24F10.3)", (gendenhf%rho(ta, a)%store(u1, u2), u2 = 1, d)
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymGenDensityHF_print

	SUBROUTINE SymGenDensityHF_del(gendenhf)
		TYPE (SymGenDensityHF), INTENT(INOUT) :: gendenhf

		INTEGER ta, a

		DO ta = 0, 1
			DO a = 0, 2*N_0
				DEALLOCATE(gendenhf%rho(ta, a)%store)
			END DO
		END DO
		DEALLOCATE(gendenhf%rho)
		RETURN
	END SUBROUTINE SymGenDensityHF_del

END MODULE symgdhf
