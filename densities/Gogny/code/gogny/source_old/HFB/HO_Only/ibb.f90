!-----------------------------------------------------------------------!
!									!
!									!
!     RADIAL INTEGRALS: BRINK-BOKER FORCE				!
!									!
!  This module computes the radial integral of tge Brink-Boeker force 	!
!  in the 2 distinct cases of an "analytical" harmonic oscillator basis !
!  and a general Woods-Saxon (or other) basis.				!
!									!
!-----------------------------------------------------------------------!

 MODULE ibb

	USE input
	USE global
	USE lgfactor
	USE symtalm

	IMPLICIT NONE
	
	! PI_BrinkBooker = 2.0 * (PI ** 3.0)
        DOUBLE PRECISION, PARAMETER :: PI_BB = 62.01255336059964035
	
 CONTAINS

        !---------------------------------------------------------------------------------!
	! Function giving the radial integral IBB for the Brink-Boeker term in the case   !
	! of the spherical harmonic oscillator basis	  			          !
	! Ref.:
        !---------------------------------------------------------------------------------!

	FUNCTION IBrinkBookerHO(na, la, nb, lb, nc, lc, nd, ld, k, x)
		DOUBLE PRECISION IBrinkBookerHO
		INTEGER, INTENT(IN) :: na, la, nb, lb, nc, lc, nd, ld, k
		DOUBLE PRECISION, INTENT(IN) :: x ! x = mi(i) / b
		INTEGER N1max, N2max, N1min, N2min
		DOUBLE PRECISION d1, d2
                DOUBLE PRECISION sumN1
		INTEGER N1

		N1max = na + nc + ((la + lc - k) / 2)
		N2max = nb + nd + ((lb + ld - k) / 2)
		N1min = MIN_5N(na, la, nc, lc, k)
		N2min = MIN_5N(nb, lb, nd, ld, k)

		! El sumatorio se deberia realizar en el sentido opuesto
		! para ir de menor a mayor en el orden de magnitud de los numeros
		
		d1 = x * x
		d2 = d1 + 2.0
		
		sumN1 = SymKumar_get(na, la, nc, lc, N1max, k) * SumBB(nb, lb, nd, ld, k, N1max, d2)
		
		DO N1 = N1max - 1, N1min, -1
			sumN1 = (SymKumar_get(na, la, nc, lc, N1, k) * SumBB(nb, lb, nd, ld, k, N1, d2)) &
				- (sumN1 * (N1 + N2min + k + 1.5) / (N1 + 1.0) / (N1 + k + 1.5) / d2)
		END DO

		IBrinkBookerHO = PI_BB * PAR(N1min + N2min) &
			* EXP(DDLogSemiFactorials(N1min + N2min + k) &
		            - DDLogFactorials(N1min) &
			    - DDLogSemiFactorials(N1min+ k) &
			    - DDLogFactorials(N2min) &
			    - DDLogSemiFactorials(N2min + k)) &
			* (x ** 3.0) * sumN1 / (d2 ** (N1min + N2min + k + 1.5))
		RETURN

	CONTAINS

       		!-----------------------------------------------------------------------!
		!  Sum over N2 of the c(N2,k) T(nb,lb,nd,ld,N2,k)*Integral(k,N1,)	1
       		!-----------------------------------------------------------------------!

		FUNCTION SumBB(n1, l1, n2, l2, k, M1, y) ! M1 <- N1
                        DOUBLE PRECISION SumBB
			INTEGER, INTENT(IN) :: n1, l1, n2, l2, k, M1
			DOUBLE PRECISION, INTENT(IN) :: y

			INTEGER N, Nmin, Nmax

			Nmax = n1 + n2 + ((l1 + l2 - k) / 2)
			Nmin = MIN_5N(n1, l1, n2, l2, k)

			SumBB = SymKumar_get(n1, l1, n2, l2, Nmax, k)
			
			DO N = Nmax - 1, Nmin, -1
				SumBB = - (SumBB * (N + M1 + k + 1.5) / (N + 1.0) / (N + k + 1.5) / y) &
					+ SymKumar_get(n1, l1, n2, l2, N, k)
			END DO
			
			RETURN
		END FUNCTION SumBB

	END FUNCTION IBrinkBookerHO

END MODULE ibb
