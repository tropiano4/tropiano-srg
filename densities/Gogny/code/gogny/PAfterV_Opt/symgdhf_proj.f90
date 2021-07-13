!---------------------------------------------------------------------------------------!
!											!
!    Module "symgdhf_Proj" is a special module similar in spirit to "symgdhf" that 	!
!    defines a new type, SymGenDensityHFProj, and the related operations. This 	!
!    type aims at dealing with projected density matrices like rho_{ac}(phi), where 	!
!    phi is a Gauge angle. The only difference between the type SymGenDensityHF defined !
!    in module "symgdhf" and the current SymGenDensityHFProj is the addition of 	!
!    a real variable that contains the value of the Gauge Angle				!
!											!
!---------------------------------------------------------------------------------------!

 MODULE symgdhf_proj

	USE input
	USE math
	USE symd3t_proj
	USE symfield_proj

	IMPLICIT NONE

	TYPE MatrixType2
		DOUBLE PRECISION, DIMENSION(:, :), POINTER :: store
	END TYPE

	! Special density type to deal with projected generalized densities. It depends on the Gauge angle.
	! rho is an array of pointers to 2D arrays. The first index refers to the isospin, the second to the
	! angular momentum. Therefore:
	!
	!           rho(ta,a) is a 2D array for isospin ta and angular momentum la
	!
	! We don't use the type SymD3Tensor for this array, although formally the structure is the same.
	! The reason is, the size is double here as we will use this type for generalized densities.
	!
	TYPE SymGenDensityHFProj
		TYPE (MatrixType2), DIMENSION(:, :), POINTER :: rho
		DOUBLE PRECISION, DIMENSION(0:1) :: GaugeAngle
	END TYPE

	INTERFACE ASSIGNMENT(=)
		MODULE PROCEDURE SymGenDensityHFProj_assign1, SymGenDensityHFProj_assign2
	END INTERFACE

 CONTAINS

	SUBROUTINE SymGenDensityHFProj_new(gendenhf)
		TYPE (SymGenDensityHFProj), INTENT(INOUT) :: gendenhf

		INTEGER ta, a, d

		ALLOCATE(gendenhf%rho(0:1, 0:(2*Lmax)))

		DO ta = 0, 1

			DO a = 0, 2*Lmax
				d = DIM(a)
				ALLOCATE(gendenhf%rho(ta, a)%store(d, d))
				gendenhf%rho(ta, a)%store(d, d) = CMPLX(0,0)
			END DO

			gendenhf%GaugeAngle(ta) = -999.9

		END DO

		RETURN
	END SUBROUTINE SymGenDensityHFProj_new

	SUBROUTINE SymGenDensityHFProj_copy(gendenhf, HF)
		TYPE (SymGenDensityHFProj), INTENT(INOUT) :: gendenhf
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: HF

		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: HF_p, HF_a
		INTEGER :: ta, a, la, d, u1, u2

		DO ta = 0, 1
			DO a = 0, 2*Lmax

				la = L(a)
				d = DIM(a)

				ALLOCATE(HF_p(d, d), HF_a(d, d))

				HF_p = CMPLX(0,0)
				HF_a = CMPLX(0,0)

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

			gendenhf%GaugeAngle(ta) =HF%GaugeAngle(ta)

		END DO

		RETURN
	END SUBROUTINE SymGenDensityHFProj_copy

	SUBROUTINE SymGenDensityHFProj_assign1(gendenhf_out, gendenhf_in)
		TYPE (SymGenDensityHFProj), INTENT(INOUT) :: gendenhf_out
		TYPE (SymGenDensityHFProj), INTENT(IN) :: gendenhf_in

		INTEGER :: ta, a

		DO ta = 0, 1

			DO a = 0, 2*Lmax
				gendenhf_out%rho(ta, a)%store = gendenhf_in%rho(ta, a)%store
			END DO

			gendenhf_out%GaugeAngle(ta) = gendenhf_in%GaugeAngle(ta)

		END DO

		RETURN
	END SUBROUTINE SymGenDensityHFProj_assign1

	SUBROUTINE SymGenDensityHFProj_assign2(HF, gendenhf)
		TYPE (SymGenDensityHFProj), INTENT(IN) :: gendenhf
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF

		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: m
		INTEGER :: ta, la, d

		DO ta = 0, 1

			DO la = 0, Lmax
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

			HF%GaugeAngle(ta) = gendenhf%GaugeAngle(ta)

		END DO

		RETURN
	END SUBROUTINE SymGenDensityHFProj_assign2

	SUBROUTINE SymGenDensityHFProj_reduce(gendenhf_out, gendenhf_in, t_in)
		TYPE (SymGenDensityHFProj), INTENT(INOUT) :: gendenhf_out
		TYPE (SymGenDensityHFProj), INTENT(IN) :: gendenhf_in
		TYPE (SymD3Tensor), INTENT(IN) :: t_in

		INTEGER :: ta, a, d, la

		DO ta = 0, 1

			DO a = 0, 2*Lmax
				d = DIM(a)
				la = L(a)
				gendenhf_out%rho(ta, a)%store = gendenhf_in%rho(ta, a)%store - t_in%d3tensor(la)%d2
			END DO

			gendenhf_out%GaugeAngle(ta) = gendenhf_in%GaugeAngle(ta)

		END DO

		RETURN
	END SUBROUTINE SymGenDensityHFProj_reduce

	SUBROUTINE SymGenDensityHFProj_print(gendenhf)
		TYPE (SymGenDensityHFProj), INTENT(IN) :: gendenhf

		INTEGER :: ta, a, d, u1, u2

		DO ta = 0, 1
			DO a = 0, 2*Lmax
				d = DIM(a)
				DO u1 = 1, d
					PRINT "(24F10.3)", (gendenhf%rho(ta, a)%store(u1, u2), u2 = 1, d)
				END DO
			END DO
		END DO

		RETURN
	END SUBROUTINE SymGenDensityHFProj_print

	SUBROUTINE SymGenDensityHFProj_del(gendenhf)
		TYPE (SymGenDensityHFProj), INTENT(INOUT) :: gendenhf

		INTEGER :: ta, a

		DO ta = 0, 1
			DO a = 0, 2*Lmax
				DEALLOCATE(gendenhf%rho(ta, a)%store)
			END DO
		END DO
		DEALLOCATE(gendenhf%rho)

		RETURN
	END SUBROUTINE SymGenDensityHFProj_del

END MODULE symgdhf_proj
