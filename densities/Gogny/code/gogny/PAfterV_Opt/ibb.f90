
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
	USE bessik

	IMPLICIT NONE

	! PI_BrinkBooker = 2.0 * (PI ** 3.0)
#if(USE_QUADRUPLE==1)
	REAL(KIND = 16), PARAMETER :: PI_BB = 62.01255336059964035
#else
	DOUBLE PRECISION, PARAMETER :: PI_BB = 62.01255336059964035
#endif

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
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) sumN1
#else
		DOUBLE PRECISION sumN1
#endif
		INTEGER N1

		N1max = na + nc + ((la + lc - k) / 2)
		N2max = nb + nd + ((lb + ld - k) / 2)
		N1min = MIN_5N(na, la, nc, lc, k)
		N2min = MIN_5N(nb, lb, nd, ld, k)

		! El sumatorio se deberia realizar en el sentido opuesto
		! para ir de menor a mayor en el orden de magnitud de los numeros

		d1 = x * x
		d2 = d1 + 2.0d0

		sumN1 = SymKumar_get(na, la, nc, lc, N1max, k) * SumBB(nb, lb, nd, ld, k, N1max, d2)

		DO N1 = N1max - 1, N1min, -1
			sumN1 = (SymKumar_get(na, la, nc, lc, N1, k) * SumBB(nb, lb, nd, ld, k, N1, d2)) &
				- (sumN1 * (N1 + N2min + k + 1.5d0) / (N1 + 1.0d0) / (N1 + k + 1.5d0) / d2)
		END DO

		IBrinkBookerHO = PI_BB * PAR(N1min + N2min) &
			* EXP(DDLogSemiFactorials(N1min + N2min + k) &
		            - DDLogFactorials(N1min) &
			    - DDLogSemiFactorials(N1min+ k) &
			    - DDLogFactorials(N2min) &
			    - DDLogSemiFactorials(N2min + k)) &
			* (x ** 3.0d0) * sumN1 / (d2 ** (N1min + N2min + k + 1.5d0))
		RETURN

	CONTAINS

       		!-----------------------------------------------------------------------!
		!  Sum over N2 of the c(N2,k) T(nb,lb,nd,ld,N2,k)*Integral(k,N1,)	1
       		!-----------------------------------------------------------------------!

		FUNCTION SumBB(n1, l1, n2, l2, k, M1, y) ! M1 <- N1
#if(USE_QUADRUPLE==1)
			REAL(KIND = 16) SumBB
#else
			DOUBLE PRECISION SumBB
#endif
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

        !---------------------------------------------------------------------------------!
	! Function giving the radial integral IBB for the Brink-Boeker term in the case   !
	! of a general spherical basis			  			          !
        !---------------------------------------------------------------------------------!

	FUNCTION IBrinkBooker(na, la, nb, lb, nc, lc, nd, ld, k, x)
		DOUBLE PRECISION IBrinkBooker

		INTEGER, INTENT(IN) :: na, la, nb, lb, nc, lc, nd, ld, k
		DOUBLE PRECISION, INTENT(IN) :: x ! x = mu(i)

		INTEGER :: IndexBra, IndexKet, Index1, Index2, IndexTest, i_r1, Iso

		DOUBLE PRECISION :: VBB, res, Pi, h, mu, Functi
		DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)

		IndexBra = IndexVecNL(na,la)
		IndexKet = IndexVecNL(nc,lc)

                Index1 = IndexVecNL(nb, lb)
                Index2 = IndexVecNL(nd, ld)

                IndexTest = Index1*Index2*IndexBra*IndexKet

		mu = x

		IF (ABS(x - 0.7) .LT. 1.E-14) Iso = 0
		IF (ABS(x - 1.2) .LT. 1.E-14) Iso = 1

		IF (IndexTest .EQ. 0) THEN
			IBrinkBooker = 0.0
			RETURN
		END IF

		Pi = 4.0*ATAN(1.0)
		h = MeshStep

	        ALLOCATE(Integrand(1:Npoint))

		DO i_r1 = 1, Npoint
			VBB = IntegralBessel(Index1, Index2, k, i_r1, mu, Iso)
			Functi = WaveFun(i_r1,IndexBra) * WaveFun(i_r1,IndexKet)
			Integrand(i_r1) = Functi * VBB
		END DO

		res = 0.0
		CALL simps(Integrand, Npoint, h, res)

		DEALLOCATE(Integrand)

		IBrinkBooker = 4.0*Pi * res

		RETURN

	CONTAINS

       		!---------------------------------------------------------------------------------------!
		!											!
		!											!
        	!---------------------------------------------------------------------------------------!

		FUNCTION IntegralBessel(IndexBra, IndexKet, k, i_r1, mu, Iso)
			DOUBLE PRECISION :: IntegralBessel

			INTEGER, INTENT(IN) :: IndexBra, IndexKet, k, i_r1, Iso
			INTEGER :: i_r2

			DOUBLE PRECISION, INTENT(IN) :: mu

			DOUBLE PRECISION :: BesselFunc, res, h, Functi,Order
			DOUBLE PRECISION, ALLOCATABLE :: Integral(:)

	       		ALLOCATE(Integral(1:Npoint))

			BesselFunc = 0.0
			h = MeshStep
			Order = k + 0.5

			DO i_r2 = 1, Npoint

				IF (Iso .EQ. 0) BesselFunc = BesselTableProton(i_r1, i_r2, k)
				IF (Iso .EQ. 1) BesselFunc = BesselTableNeutron(i_r1, i_r2, k)

				Functi = WaveFun(i_r2,IndexBra) * WaveFun(i_r2,IndexKet)

				Integral(i_r2) = Functi * BesselFunc

			END DO

			CALL simps(Integral, Npoint, h, res)

			DEALLOCATE(Integral)

			IntegralBessel = res

			RETURN
		END FUNCTION IntegralBessel


	END FUNCTION IBrinkBooker

END MODULE ibb
