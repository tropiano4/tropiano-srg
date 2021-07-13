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

	!-----------------------------------------------------------------------!
	!  Subroutine that calculates the 
	!-----------------------------------------------------------------------!

	SUBROUTINE SymGDDph_update(HF_out, gDDph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymGDDph), INTENT(INOUT) :: gDDph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		INTEGER ta, la, na, namax, nc

		! Filling in the new density-dependent fields Gamma_{na, nc, la, ta} from the new density

		! Calculate the new density for both isospin channels from the current wave-functions and fields (HF_in)
		
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

	!-----------------------------------------------------------------------!
	!									!
	!  Subroutine that calculates the density-dependent HF field 		!
	!  Gamma_{na, nc, la, ta} of page page 134 (but without the last term 	!
	!  ontaining the pairing tensor (for each isospin). It proceeds to the	!
	!  integration of the auxiliary function Gamma^{DD}(r) of page 133.	!
	!									!
	! 	      CASE OF A SPHERICAL HARMONIC OSCILLATOR BASIS		!
	!									!
	!-----------------------------------------------------------------------!

	FUNCTION SymGDDph_G1dd(gDDph, na, nc, la, ta)
		DOUBLE PRECISION SymGDDph_G1dd
		TYPE (SymGDDph), INTENT(IN) :: gDDph
		INTEGER, INTENT(IN) :: na, nc, la, ta

		DOUBLE PRECISION d1, d2, d3, d4, d5
		INTEGER p1, p1max, i
                DOUBLE PRECISION sump1, sumi

		p1max = na + nc
		sump1 = 0.0
		
		DO p1 = 0, p1max
			sumi = 0.0
			
			! Below, there is a hidden integration (due to the Laguerre points)
			
			DO i = 1, NLag
			
				! d1 = density for isospin  ta  at point ri (ri is a node for the Laguerre integration)
				! d2 = density for isospin 1-ta at point ri (ri is a node for the Laguerre integration)
			
				d1 = gDDph%dLag(    ta, i)
				d2 = gDDph%dLag(1 - ta, i)
				
				! d3 = total density at point ri (ri is a node for the Laguerre integration)
				
				d3 = d1 + d2
				
				IF (d3 .LT. 0.0) STOP "ERROR! Fuera de rango"
				
				! d4 IS UNUSED
				d4 = d3 ** (ALPHA - 1.0)
				
				d5 = (1.0 + 0.5 * x0) * d3 * d3 - (x0 + 0.5) * d1 * (d3 - ALPHA * d2)
! aLag = x/(ALPHA+2.)
!TODO			sumi =  sumi + EXP(LOG(GaussLQ%gauss%w(i)) + LOG(aLag(i)) * (p1 + la) + LOG(d3) * (ALPHA - 1) + LOG(d5))

				! 
				
				sumi =  sumi + EXP(LOG(GaussLQ%gauss%w(i)) &
					+ LOG(GaussLQ%gauss%x(i) / (ALPHA + 2.0)) &
					* DBLE(p1 + la) + LOG(d3) * (ALPHA - 1) + LOG(d5))
					
			END DO
			sump1 = sump1 + SymCoefficientB_get(na, la, nc, la, p1) * sumi
		END DO
		
		SymGDDph_G1dd = Gogny_t0(Gogny) * 0.5 * I_SALPHA3 * sump1
		
		RETURN
	END FUNCTION SymGDDph_G1dd

	!-----------------------------------------------------------------------!
	!   Subroutine that calculates the density for each isospin channel	!
	!									!
	! 	      CASE OF A SPHERICAL HARMONIC OSCILLATOR BASIS		!
	!									!
	!-----------------------------------------------------------------------!

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
			
			! Calculation of rho^tau(r) for each isospin (stored in wksp)
		
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
		
		! Schunck: I added the 1/b_0**3 to put all quantities dependent of the oscillator length apart
		!          (that is to say in the calculation of the radial density). This implies that the 
		!          final result obtained from SymGDDph_get_edd should be multipied by some funny factor
		!          to recover the good result. See my notes for further explanations. The point is: 
		!	   the calculation of the DD FIELD is now independent of the oscillator length, AS IT
		!	   SHOULD BE.
		
		DO s = 0, N_0
			DO i = 1, NLag
				gDDph%dLag(1, i) = gDDph%dLag(1, i) + (pows(i) * wksp(1, s))/(b_0**3)
				gDDph%dLag(0, i) = gDDph%dLag(0, i) + (pows(i) * wksp(0, s))/(b_0**3)
! aLag = x/(ALPHA+2.)
!TODO			pows(i) = pows(i) * aLag(i)

				pows(i) = pows(i) * (GaussLQ%gauss%x(i) / DBLE(ALPHA + 2.0))
			END DO
		END DO
		
		DEALLOCATE(wksp)
		DEALLOCATE(pows)
		
		RETURN
	END SUBROUTINE SymGDDph_make_DD

	!-----------------------------------------------------------------------!
	!   This subroutine gives the density-dependent part of the mean-	!
	!   field energy by calculating 0.5*Tr[ Gamma^{DD}_ {ij} rho) ]. This 	!
	!   comes down to integrating the auxiliary field Gamma^{DD}(r) x 	!
	!   rho(r). Note that this function must give the sum of the proton 	!
	!   and neutron	contributions.						!
	!-----------------------------------------------------------------------!

	FUNCTION SymGDDph_get_edd(gDDph)
		DOUBLE PRECISION SymGDDph_get_edd
		TYPE (SymGDDph), INTENT(IN) :: gDDph

		DOUBLE PRECISION sumi, d1, d2, d3, res
		DOUBLE PRECISION, ALLOCATABLE:: Integrand(:)

		INTEGER i, N
		
		sumi = 0.0
		
		DO i = 1, NLag
			d1 = gDDph%dLag(NEUTRON, i) * b_0**(9./7.)
			d2 = gDDph%dLag(PROTON,  i) * b_0**(9./7.)
			d3 = d1 + d2
			IF (d3 .LT. 0.0) STOP "ERROR - Negative density in SymGDDph_get_edd"
			sumi = sumi + GaussLQ%gauss%w(i) * (d3 ** ALPHA) * d1 * d2
		END DO
			
		res = 3.0*PI * I_SALPHA3 * Gogny_t0(Gogny) * sumi
					
		SymGDDph_get_edd = res
		RETURN
	END FUNCTION SymGDDph_get_edd

	SUBROUTINE SymGDDph_del(gDDph)
		TYPE (SymGDDph), INTENT(INOUT) :: gDDph

		DEALLOCATE(gDDph%dLag)
		RETURN
	END SUBROUTINE SymGDDph_del

END MODULE symgdd
