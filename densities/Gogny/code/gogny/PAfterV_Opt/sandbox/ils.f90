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
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: Log2
		REAL(KIND = 16) :: s
#else
		DOUBLE PRECISION :: Log2
		DOUBLE PRECISION :: s
#endif

		Log2 = LOG(2.0D0)

		pmax = na + nb
		x = la + lb - 0.5D0

		s = SymCoefficientB_get(na, la, nb, lb, pmax) * Sum1(nc, la, nd, lb, pmax + x)
		DO p = pmax - 1, 0, -1
			s = SymCoefficientB_get(na, la, nb, lb, p) * Sum1(nc, la, nd, lb, p + x) &
				+ 0.5D0 * (p + x + 1.0D0) * s
		END DO
		IPLSHO = s * EXP(DDLogSemiFactorials(la + lb - 1) - ((la + lb + 1.0D0) * Log2)) / SQRT(2.0D0) / (b_0 ** 5)
		RETURN

	CONTAINS

		FUNCTION Sum1(n1, l1, n2, l2, x)
#if(USE_QUADRUPLE==1)
			REAL(KIND = 16) Sum1
#else
			DOUBLE PRECISION Sum1
#endif
			INTEGER, INTENT(IN) :: n1, l1, n2, l2
			DOUBLE PRECISION, INTENT(IN) :: x

			INTEGER p, pmax

			pmax = n1 + n2
			Sum1 = SymCoefficientB_get(n1, l1, n2, l2, pmax)
			DO p = pmax - 1, 0, -1
				Sum1 = 0.5D0 * (p + x + 1.0D0) * Sum1 + SymCoefficientB_get(n1, l1, n2, l2, p)
			END DO
			RETURN
		END FUNCTION Sum1

	END FUNCTION IPLSHO

	!---------------------------------------------------------------------!
	! Radial Integral IPLS in the general case of a spherical basis       !
	!---------------------------------------------------------------------!

	FUNCTION IPLS(na, nc, la, nb, nd, lb)
		DOUBLE PRECISION IPLS
		INTEGER, INTENT(IN) :: na, nc, la, nb, nd, lb

		DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)

		INTEGER IndexBraOne, IndexKetOne, IndexBraTwo, IndexKetTwo, Ipoint

		DOUBLE PRECISION RFourth,Result

		IndexBraOne = IndexVecNL(na,la)
		IndexBraTwo = IndexVecNL(nb,lb)

		IndexKetOne = IndexVecNL(nc,la)
		IndexKetTwo = IndexVecNL(nd,lb)

		IF (IndexBraOne .EQ. 0 .OR. IndexKetOne .EQ. 0 .OR. IndexBraTwo .EQ. 0 .OR. IndexKetTwo .EQ. 0) THEN
			IPLS = 0.0D0
			RETURN
		END IF

	        ALLOCATE(Integrand(1:Npoint))

		DO Ipoint = 1, Npoint
			RFourth = RadMesh(Ipoint)**4
			Integrand(Ipoint) = WaveFun(Ipoint,IndexBraOne)*WaveFun(Ipoint,IndexBraTwo)* &
					    WaveFun(Ipoint,IndexKetOne)*WaveFun(Ipoint,IndexKetTwo)/RFourth
		END DO

		CALL simps(Integrand,Npoint,MeshStep,Result)

		DEALLOCATE(Integrand)

		IPLS = Result

		RETURN

	END FUNCTION IPLS

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
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: Log2
		REAL(KIND = 16) :: s
#else
		DOUBLE PRECISION :: Log2
		DOUBLE PRECISION :: s
#endif

		Log2 = LOG(2.0D0)

		pmax = na + nc
		x = DBLE(la + lb) - 0.5D0
		y = DBLE(lb - la) - 0.5D0

		s = SymCoefficientB_get(na, la, nc, la, pmax) * Sum2(nb, lb, nd, lb, pmax + x, y - pmax)

		DO p = pmax - 1, 0, -1
			s = (SymCoefficientB_get(na, la, nc, la, p) * Sum2(nb, lb, nd, lb, p + x, y - p)) &
				+ 0.5D0 * (DBLE(p+1) + x) * s
		END DO
		IHFLSHO = s * EXP(DDLogSemiFactorials(la + lb - 1) - (DBLE(la + lb + 1) * Log2)) / SQRT(2.0D0) / (b_0 ** 5)
		RETURN

	CONTAINS

		FUNCTION Sum2(n1, l1, n2, l2, x, y)
#if(USE_QUADRUPLE==1)
			REAL(KIND = 16) Sum2
#else
			DOUBLE PRECISION Sum2
#endif
			INTEGER, INTENT(IN) :: n1, l1, n2, l2
			DOUBLE PRECISION, INTENT(IN) :: x, y

			INTEGER p, pmax

			pmax = n1 + n2
			Sum2 = SymCoefficientB_get(n1, l1, n2, l2, pmax) * (pmax + y)
			DO p = pmax - 1, 0, -1
				Sum2 = 0.5D0 * (DBLE(p+1) + x) * Sum2 + (SymCoefficientB_get(n1, l1, n2, l2, p) * (DBLE(p) + y))
			END DO
			RETURN
		END FUNCTION Sum2

	END FUNCTION IHFLSHO

	!---------------------------------------------------------------------!
	! Radial Integral IHFLS in the general case of a spherical basis      !
	!---------------------------------------------------------------------!

	FUNCTION IHFLS(na, nc, la, nb, nd, lb)
		DOUBLE PRECISION IHFLS
		INTEGER, INTENT(IN) :: na, nc, la, nb, nd, lb

		DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)

		INTEGER IndexBraOne, IndexKetOne, IndexBraTwo, IndexKetTwo, Ipoint

		DOUBLE PRECISION Result

		IndexBraOne = IndexVecNL(na,la)
		IndexBraTwo = IndexVecNL(nc,la)

		IndexKetOne = IndexVecNL(nb,lb)
		IndexKetTwo = IndexVecNL(nd,lb)

		IF (IndexBraOne .EQ. 0 .OR. IndexKetOne .EQ. 0 .OR. IndexBraTwo .EQ. 0 .OR. IndexKetTwo .EQ. 0) THEN
			IHFLS = 0.0D0
			RETURN
		END IF

	        ALLOCATE(Integrand(1:Npoint))

		DO Ipoint = 1, Npoint
			Integrand(Ipoint) = WaveFun(Ipoint,IndexBraOne) *WaveFun(Ipoint,IndexBraTwo)* &
					  ( WaveFun(Ipoint,IndexKetOne) *WaveDeri(Ipoint,IndexKetTwo) + &
					    WaveDeri(Ipoint,IndexKetOne)*WaveFun(Ipoint,IndexKetTwo) - &
			              2.0D0*WaveFun(Ipoint,IndexKetOne) *WaveFun(Ipoint,IndexKetTwo)/RadMesh(Ipoint) ) &
					   /RadMesh(Ipoint)**3
		END DO

		CALL simps(Integrand,Npoint,MeshStep,Result)

		DEALLOCATE(Integrand)

		IHFLS = Result

		RETURN

	END FUNCTION IHFLS

END MODULE ils
