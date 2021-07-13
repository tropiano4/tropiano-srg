
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

 CONTAINS

        !---------------------------------------------------------------------------------!
	! Function giving the radial integral IBB for the Brink-Boeker term in the case   !
	! of the spherical harmonic oscillator basis	  			          !
	! Ref.:
        !---------------------------------------------------------------------------------!

	FUNCTION IBrinkBookerHO(na, la, nb, lb, nc, lc, nd, ld, k, x)
		DOUBLE PRECISION :: IBrinkBookerHO
		INTEGER, INTENT(IN) :: na, la, nb, lb, nc, lc, nd, ld, k
		DOUBLE PRECISION, INTENT(IN) :: x ! x = mi(i) / b
		INTEGER :: N1max, N2max, N1min, N2min
		DOUBLE PRECISION :: d1, d2
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: PI_BB
		REAL(KIND = 16) :: sumN1
#else
		DOUBLE PRECISION :: PI_BB
		DOUBLE PRECISION :: sumN1
#endif
		INTEGER :: N1

		PI_BB = 2.0D0 * ( 4.0D0*ATAN(1.0D0) )**3

		N1max = na + nc + ((la + lc - k) / 2)
		N2max = nb + nd + ((lb + ld - k) / 2)
		N1min = MIN_5N(na, la, nc, lc, k)
		N2min = MIN_5N(nb, lb, nd, ld, k)

		! El sumatorio se deberia realizar en el sentido opuesto
		! para ir de menor a mayor en el orden de magnitud de los numeros

		d1 = x * x
		d2 = d1 + 2.0D0

		sumN1 = SymKumar_get(na, la, nc, lc, N1max, k) * SumBB(nb, lb, nd, ld, k, N1max, d2)

		DO N1 = N1max - 1, N1min, -1
			sumN1 = (SymKumar_get(na, la, nc, lc, N1, k) * SumBB(nb, lb, nd, ld, k, N1, d2)) &
				- (sumN1 * (DBLE(N1 + N2min + k) + 1.5d0) / DBLE(N1 + 1) / (DBLE(N1 + k) + 1.5d0) / d2)
		END DO

		IBrinkBookerHO = PI_BB * PAR(N1min + N2min) &
			* EXP(DDLogSemiFactorials(N1min + N2min + k) &
		            - DDLogFactorials(N1min) &
			    - DDLogSemiFactorials(N1min+ k) &
			    - DDLogFactorials(N2min) &
			    - DDLogSemiFactorials(N2min + k)) &
			* (x ** 3.0D0) * sumN1 / (d2 ** (DBLE(N1min + N2min + k) + 1.5d0))
		RETURN

	CONTAINS

       		!-----------------------------------------------------------------------!
		!  Sum over N2 of the c(N2,k) T(nb,lb,nd,ld,N2,k)*Integral(k,N1,)	1
       		!-----------------------------------------------------------------------!

		FUNCTION SumBB(n1, l1, n2, l2, k, M1, y) ! M1 <- N1
#if(USE_QUADRUPLE==1)
			REAL(KIND = 16) :: SumBB
#else
			DOUBLE PRECISION :: SumBB
#endif
			INTEGER, INTENT(IN) :: n1, l1, n2, l2, k, M1
			DOUBLE PRECISION, INTENT(IN) :: y

			INTEGER :: N, Nmin, Nmax

			Nmax = n1 + n2 + ((l1 + l2 - k) / 2)
			Nmin = MIN_5N(n1, l1, n2, l2, k)

			SumBB = SymKumar_get(n1, l1, n2, l2, Nmax, k)

			DO N = Nmax - 1, Nmin, -1
				SumBB = - (SumBB * (DBLE(N + M1 + k) + 1.5D0) / DBLE(N + 1) / (DBLE(N + k) + 1.5D0) / y) &
					+ SymKumar_get(n1, l1, n2, l2, N, k)
			END DO

			RETURN
		END FUNCTION SumBB

	END FUNCTION IBrinkBookerHO

        !---------------------------------------------------------------------------------!
	! Function giving the radial integral IBB for the Brink-Boeker term in the case   !
	! of a general spherical basis			  			          !
        !---------------------------------------------------------------------------------!

	FUNCTION IBrinkBooker(na, la, nb, lb, nc, lc, nd, ld, k, i, n_reg)
		DOUBLE PRECISION :: IBrinkBooker

		INTEGER, INTENT(IN) :: na, la, nb, lb, nc, lc, nd, ld, k, i, n_reg

		INTEGER :: IndexBra, IndexKet, Index1, Index2, IndexTest, i_r1

		DOUBLE PRECISION :: VBB, res, Pi, h, Functi
		DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)

		IndexBra = IndexVecNL(na,la)
		IndexKet = IndexVecNL(nc,lc)

                Index1 = IndexVecNL(nb, lb)
                Index2 = IndexVecNL(nd, ld)

                IndexTest = Index1*Index2*IndexBra*IndexKet

		IF (IndexTest .EQ. 0) THEN
			IBrinkBooker = 0.0D0
			RETURN
		END IF

		Pi = 4.0D0*ATAN(1.0D0)
		h = MeshStep

	        ALLOCATE(Integrand(1:Npoint))

!$OMP Parallel Default(None) &
!$OMP& SHARED(Npoint,Index1,Index2,k,i,n_reg,WaveFun,IndexBra,IndexKet,Integrand) &
!$OMP& PRIVATE(i_r1,VBB,Functi)
!$OMP DO SCHEDULE(DYNAMIC)
		DO i_r1 = 1, Npoint
			VBB = IntegralBessel(Index1, Index2, k, i_r1, i, n_reg)
			Functi = WaveFun(i_r1,IndexBra) * WaveFun(i_r1,IndexKet)
			Integrand(i_r1) = Functi * VBB
		END DO
!$OMP End Do
!$OMP End Parallel

		res = 0.0D0
		CALL simps(Integrand, Npoint, h, res)

		DEALLOCATE(Integrand)

		IBrinkBooker = 4.0D0*Pi * res

		RETURN

	CONTAINS

       		!---------------------------------------------------------------------------------------!
		!											!
		!											!
        	!---------------------------------------------------------------------------------------!

		FUNCTION IntegralBessel(IndexBra, IndexKet, k, i_r1, i, n_reg)
			DOUBLE PRECISION :: IntegralBessel

			INTEGER, INTENT(IN) :: IndexBra, IndexKet, k, i_r1, i, n_reg
			INTEGER :: i_r2

			DOUBLE PRECISION :: BesselFunc, res, h, Functi, Order, r
			DOUBLE PRECISION, ALLOCATABLE :: Integral(:)

	       		ALLOCATE(Integral(1:Npoint))

			BesselFunc = 0.0D0
			h = MeshStep
			Order = k + 0.5D0

			DO i_r2 = 1, Npoint

				BesselFunc = BesselTable(i_r1, i_r2, k, i)

				Functi = WaveFun(i_r2,IndexBra) * WaveFun(i_r2,IndexKet)

				r = ABS(RadMesh(i_r1) - RadMesh(i_r2))

				Integral(i_r2) = Functi * BesselFunc * (r**n_reg)

			END DO

			CALL simps(Integral, Npoint, h, res)

			DEALLOCATE(Integral)

			IntegralBessel = res

			RETURN
		END FUNCTION IntegralBessel

	END FUNCTION IBrinkBooker

END MODULE ibb
