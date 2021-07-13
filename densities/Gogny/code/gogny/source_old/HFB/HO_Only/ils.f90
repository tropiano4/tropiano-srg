!---------------------------------------------------------------------!
!                                                                     !
!     RADIAL INTEGRALS: SPIN-ORBIT FORCE                              !                
!                                                                     !
!---------------------------------------------------------------------!

 MODULE ils

	USE input
	USE global
	USE symtalm

	IMPLICIT NONE

        DOUBLE PRECISION, PARAMETER :: Log2 = 0.6931471805599453094172321214581765680755
	DOUBLE PRECISION, PARAMETER :: SQRT_2 = 1.414213562373095

 CONTAINS

	!---------------------------------------------------------------------!
	! Radial Integral IPLS in the case of the harmonic oscillator basis   ! 
	! Refs.: PhD, Page 137, E18                                           !              
	!---------------------------------------------------------------------!

	FUNCTION IPLSHO(na, nc, la, nb, nd, lb)
		DOUBLE PRECISION IPLSHO
		INTEGER, INTENT(IN) :: na, nc, la, nb, nd, lb

		INTEGER p, pmax
		DOUBLE PRECISION x
                DOUBLE PRECISION s

		! OJO: vale para b1==b2

		pmax = na + nb
		x = la + lb - 0.5

		s = SymCoefficientB_get(na, la, nb, lb, pmax) * Sum1(nc, la, nd, lb, pmax + x)
		DO p = pmax - 1, 0, -1
			s = SymCoefficientB_get(na, la, nb, lb, p) * Sum1(nc, la, nd, lb, p + x) &
				+ 0.5 * (p + x + 1.0) * s
		END DO
		IPLSHO = s * EXP(DDLogSemiFactorials(la + lb - 1) - ((la + lb + 1.0) * Log2)) / SQRT_2 / (b_0 ** 5)
		RETURN

	CONTAINS

		FUNCTION Sum1(n1, l1, n2, l2, x)
                DOUBLE PRECISION Sum1
			INTEGER, INTENT(IN) :: n1, l1, n2, l2
			DOUBLE PRECISION, INTENT(IN) :: x

			INTEGER p, pmax

			pmax = n1 + n2
			Sum1 = SymCoefficientB_get(n1, l1, n2, l2, pmax)
			DO p = pmax - 1, 0, -1
				Sum1 = 0.5 * (p + x + 1.0) * Sum1 + SymCoefficientB_get(n1, l1, n2, l2, p)
			END DO
			RETURN
		END FUNCTION Sum1

	END FUNCTION IPLSHO

	!---------------------------------------------------------------------!
	! Radial Integral IHFLS in the case of the harmonic oscillator basis  ! 
	! Refs.: PhD, Page 137, Sec. E3.3                                     !               
	!---------------------------------------------------------------------!

	FUNCTION IHFLSHO(na, nc, la, nb, nd, lb)
		DOUBLE PRECISION IHFLSHO
		INTEGER, INTENT(IN) :: na, nc, la, nb, nd, lb

                ! OJO: solo vale cuando b1==b2

		INTEGER p, pmax
		DOUBLE PRECISION x, y
                DOUBLE PRECISION s

		pmax = na + nc
		x = la + lb - 0.5
		y = lb - la - 0.5

		s = SymCoefficientB_get(na, la, nc, la, pmax) * Sum2(nb, lb, nd, lb, pmax + x, y - pmax)
		
		!WRITE(*,'("Coef1 = ",ES12.5," Coef2 = ",ES12.5)') SymCoefficientB_get(na, la, nc, la, pmax), Sum2(nb, lb, nd, lb, pmax + x, y - pmax)
		
		DO p = pmax - 1, 0, -1
			s = (SymCoefficientB_get(na, la, nc, la, p) * Sum2(nb, lb, nd, lb, p + x, y - p)) &
				+ 0.5 * (p + x + 1.0) * s
		END DO
		IHFLSHO = s * EXP(DDLogSemiFactorials(la + lb - 1) - ((la + lb + 1.0) * Log2)) / SQRT_2 / (b_0 ** 5)
		RETURN

	CONTAINS

		FUNCTION Sum2(n1, l1, n2, l2, x, y)
                        DOUBLE PRECISION Sum2
			INTEGER, INTENT(IN) :: n1, l1, n2, l2
			DOUBLE PRECISION, INTENT(IN) :: x, y

			INTEGER p, pmax

			pmax = n1 + n2
			Sum2 = SymCoefficientB_get(n1, l1, n2, l2, pmax) * (pmax + y)
			DO p = pmax - 1, 0, -1
				Sum2 = 0.5 * (p + x + 1.0) * Sum2 + (SymCoefficientB_get(n1, l1, n2, l2, p) * (p + y))
			END DO
			RETURN
		END FUNCTION Sum2

	END FUNCTION IHFLSHO
	
END MODULE ils
