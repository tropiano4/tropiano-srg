!---------------------------------------------------------------------!
!                                                                     !
!     RADIAL INTEGRALS: COULOMB FORCE                                 !                
!                                                                     !
!---------------------------------------------------------------------!

 MODULE ic

	USE input
	USE global
	USE lgfactor
	USE symtalm
	USE math

	IMPLICIT NONE

		DOUBLE PRECISION, PARAMETER :: VC = 1.44197028
		! PI_COU = SQRT(2.0) * (PI ** (5.0 / 2.0))
                DOUBLE PRECISION, PARAMETER :: PI_COU = 24.73942945119314805

 CONTAINS

	!-----------------------------------------------------------------------!
	!     Radial Integral as defined in Appendix H for the spherical	!
	!     harmonic oscillator basis						!                
	!-----------------------------------------------------------------------!

	FUNCTION ICoulombHO(na, la, nb, lb, nc, lc, nd, ld, k)
		DOUBLE PRECISION ICoulombHO
		INTEGER, INTENT(IN) :: na, la, nb, lb, nc, lc, nd, ld, k

		INTEGER N1, N1min, N1max, N2min, N2max
                DOUBLE PRECISION sumN1

		N1max = na + nc + ((la + lc - k) / 2)
		N2max = nb + nd + ((lb + ld - k) / 2)
		N1min = MIN_5N(na, la, nc, lc, k)
		N2min = MIN_5N(nb, lb, nd, ld, k)

		sumN1 = SymKumar_get(na, la, nc, lc, N1max, k) * SumC(nb, lb, nd, ld, k, N1max, DBLE(2.0))
		DO N1 = N1max - 1, N1min, -1
			sumN1 = (-sumN1 * (N1 + N2min + k + 0.5) / (N1 + 1.0) / (N1 + k + 1.5) / 2.0) &
				+ SymKumar_get(na, la, nc, lc, N1, k) * SumC(nb, lb, nd, ld, k, N1, DBLE(2.0))
		END DO

		ICoulombHO = PI_COU * PAR(N1min + N2min) &
			* EXP(DDLogSemiFactorials(N1min + N2min + k - 1) &
			    - DDLogFactorials(N1min) &
			    - DDLogSemiFactorials(N1min + k) &
			    - DDLogFactorials(N2min) &
			    - DDLogSemiFactorials(N2min + k)) &
			* (sumN1 * (2.0 ** DBLE(-N1min - N2min - k))) / b_0
		RETURN

	CONTAINS

		FUNCTION SumC(n1, l1, n2, l2, k, M1, y)
                        DOUBLE PRECISION SumC
			INTEGER, INTENT(IN) :: n1, l1, n2, l2, k, M1
			DOUBLE PRECISION, INTENT(IN) :: y

			INTEGER N, Nmin, Nmax

			Nmax = n1 + n2 + ((l1 + l2 - k) / 2)
	 		Nmin = MIN_5N(n1, l1, n2, l2, k)

			SumC = SymKumar_get(n1, l1, n2, l2, Nmax, k)
			DO N = Nmax - 1, Nmin, -1
				SumC = - (SumC * (N + M1 + k + 0.5) / (N + 1.0) / (N + k + 1.5) / y) &
					+ SymKumar_get(n1, l1, n2, l2, N, k)
			END DO
			RETURN
		END FUNCTION SumC

	END FUNCTION ICoulombHO

END MODULE ic
