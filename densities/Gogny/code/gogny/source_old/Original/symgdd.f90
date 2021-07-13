MODULE symgdd

	USE input
	USE global
	USE math
	USE symtalm
	USE symfield

	IMPLICIT NONE

	TYPE SymGDDph
		DOUBLE PRECISION, DIMENSION(:, :), POINTER :: dLag
	END TYPE

CONTAINS

	SUBROUTINE SymGDDph_new(gDDph)
		TYPE (SymGDDph), INTENT(INOUT) :: gDDph

		ALLOCATE(gDDph%dLag(0:2, NLag))
		IF (.NOT. ASSOCIATED(gDDph%dLag)) STOP "Unable to allocate memory"
		RETURN
	END SUBROUTINE SymGDDph_new

	SUBROUTINE SymGDDph_update(HF_out, gDDph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymGDDph), INTENT(INOUT) :: gDDph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		INTEGER ta, la, na, namax, nc

		CALL SymGDDph_make_DD(gDDph, HF_in)
		DO ta = 0, 1
			DO la = 0, N_0
				namax = ((N_0 - la) / 2) + 1
				DO na = 1, namax
					DO nc = 1, na
						HF_out%p(ta)%d3tensor(la)%d2(na, nc) = SymGDDph_G1dd(gDDph, na - 1, nc - 1, la, ta)
						HF_out%a(ta)%d3tensor(la)%d2(na, nc) = 0.0
					END DO
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymGDDph_update

	FUNCTION SymGDDph_G1dd(gDDph, na, nc, la, ta)
		DOUBLE PRECISION SymGDDph_G1dd
		TYPE (SymGDDph), INTENT(IN) :: gDDph
		INTEGER, INTENT(IN) :: na, nc, la, ta

		DOUBLE PRECISION d1, d2, d3, d4, d5
		INTEGER p1, p1max, i
		REAL(KIND = 16) sump1, sumi

		p1max = na + nc
		sump1 = 0.0
		DO p1 = 0, p1max
			sumi = 0.0
			DO i = 1, NLag
				d1 = gDDph%dLag(    ta, i)
				d2 = gDDph%dLag(1 - ta, i)
				d3 = d1 + d2
				IF (d3 .LT. 0.0) STOP "ERROR! Fuera de rango"
				d4 = d3 ** (ALPHA - 1.0)
				d5 = (1.0 + 0.5 * x0) * d3 * d3 - (x0 + 0.5) * d1 * (d3 - ALPHA * d2)
! aLag = x/(ALPHA+2.)
!TODO			sumi =  sumi + EXP(LOG(GaussLQ%gauss%w(i)) + LOG(aLag(i)) * (p1 + la) + LOG(d3) * (ALPHA - 1) + LOG(d5))
				sumi =  sumi + EXP(LOG(GaussLQ%gauss%w(i)) &
					+ LOG(GaussLQ%gauss%x(i) / (ALPHA + 2.0)) &
					* DBLE(p1 + la) + LOG(d3) * (ALPHA - 1) + LOG(d5))
			END DO
			sump1 = sump1 + SymCoefficientB_get(na, la, nc, la, p1) * sumi
		END DO
		SymGDDph_G1dd = Gogny_t0(Gogny) * 0.5 * sump1 * I_SALPHA3
		RETURN
	END FUNCTION SymGDDph_G1dd

	SUBROUTINE SymGDDph_make_DD(gDDph, HF)
		TYPE (SymGDDph), INTENT(INOUT) :: gDDph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF

		DOUBLE PRECISION d1
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: wksp
		INTEGER s, i, lb, nb, nbmax, nd, p2, p2max
		DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: pows

		ALLOCATE(wksp(0:1, 0:N_0))
		ALLOCATE(pows(NLag))
		DO s = 0, N_0
			wksp(1, s) = 0.0
			wksp(0, s) = 0.0
			DO lb = 0, s
				nbmax = ((N_0 - lb) / 2) + 1
				p2 = s - lb
				DO nb = 1, nbmax
					DO nd = 1, nb
						p2max = nb + nd - 2
						IF (p2 .GT. p2max) CYCLE
						IF (nb .EQ. nd) THEN
							d1 = DBLE(      I_4PI * SymCoefficientB_get(nb - 1, lb, nd - 1, lb, p2))
						ELSE
							d1 = DBLE(2.0 * I_4PI * SymCoefficientB_get(nb - 1, lb, nd - 1, lb, p2))
						END IF
						wksp(1, s) = wksp(1, s) + (d1 * HF%p(1)%d3tensor(lb)%d2(nb, nd))
						wksp(0, s) = wksp(0, s) + (d1 * HF%p(0)%d3tensor(lb)%d2(nb, nd))
					END DO
				END DO
			END DO
		END DO

		! Inicializamos a 1 la tabla "pows" y a 0 las tablas "dLag"
		DO i = 1, NLag
			pows(i) = 1.0
			gDDph%dLag(1, i) = 0.0
			gDDph%dLag(0, i) = 0.0
		END DO
		DO s = 0, N_0
			DO i = 1, NLag
				gDDph%dLag(1, i) = gDDph%dLag(1, i) + (pows(i) * wksp(1, s))
				gDDph%dLag(0, i) = gDDph%dLag(0, i) + (pows(i) * wksp(0, s))
! aLag = x/(ALPHA+2.)
!TODO			pows(i) = pows(i) * aLag(i)
				pows(i) = pows(i) * (GaussLQ%gauss%x(i) / DBLE(ALPHA + 2.0))
			END DO
		END DO
		DEALLOCATE(wksp)
		DEALLOCATE(pows)
		RETURN
	END SUBROUTINE SymGDDph_make_DD

	FUNCTION SymGDDph_get_edd(gDDph)
		DOUBLE PRECISION SymGDDph_get_edd
		TYPE (SymGDDph), INTENT(IN) :: gDDph

		DOUBLE PRECISION sumi, d1, d2, d3
		INTEGER i

		sumi = 0.0
		DO i = 1, NLag
			d1 = gDDph%dLag(NEUTRON, i)
			d2 = gDDph%dLag(PROTON,  i)
			d3 = d1 + d2
			IF (d3 .LT. 0.0) STOP "ERROR! Fuera de rango"
			sumi = sumi + GaussLQ%gauss%w(i) * (d3 ** ALPHA) * d1 * d2
		END DO
		SymGDDph_get_edd = 3.0 * Gogny_t0(Gogny) * I_SALPHA3 * PI * sumi
		RETURN
	END FUNCTION SymGDDph_get_edd

	SUBROUTINE SymGDDph_del(gDDph)
		TYPE (SymGDDph), INTENT(INOUT) :: gDDph

		DEALLOCATE(gDDph%dLag)
		RETURN
	END SUBROUTINE SymGDDph_del

END MODULE symgdd
