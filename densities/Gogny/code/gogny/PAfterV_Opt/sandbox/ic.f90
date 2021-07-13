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

	DOUBLE PRECISION, PARAMETER :: VC = 1.4399784085965135298D0

 CONTAINS

	!-----------------------------------------------------------------------!
	!     Radial Integral as defined in Appendix H for the spherical	!
	!     harmonic oscillator basis						!
	!-----------------------------------------------------------------------!

	FUNCTION ICoulombHO(na, la, nb, lb, nc, lc, nd, ld, k)
		DOUBLE PRECISION :: ICoulombHO
		INTEGER, INTENT(IN) :: na, la, nb, lb, nc, lc, nd, ld, k

		INTEGER :: N1, N1min, N1max, N2min, N2max
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: PI_COU
		REAL(KIND = 16) :: sumN1
#else
		DOUBLE PRECISION :: PI_COU
		DOUBLE PRECISION :: sumN1
#endif

		PI_COU = SQRT(2.0D0) * ( 4.0D0 * ATAN(1.0D0) )**(2.5D0)

		N1max = na + nc + ((la + lc - k) / 2)
		N2max = nb + nd + ((lb + ld - k) / 2)
		N1min = MIN_5N(na, la, nc, lc, k)
		N2min = MIN_5N(nb, lb, nd, ld, k)

		sumN1 = SymKumar_get(na, la, nc, lc, N1max, k) * SumC(nb, lb, nd, ld, k, N1max, DBLE(2.0))
		DO N1 = N1max - 1, N1min, -1
			sumN1 = (-sumN1 * (DBLE(N1 + N2min + k) + 0.5D0) / DBLE(N1 + 1) / (DBLE(N1 + k) + 1.5) / 2.0D0) &
				+ SymKumar_get(na, la, nc, lc, N1, k) * SumC(nb, lb, nd, ld, k, N1, DBLE(2.0))
		END DO

		ICoulombHO = PI_COU * PAR(N1min + N2min) &
			* EXP(DDLogSemiFactorials(N1min + N2min + k - 1) &
			    - DDLogFactorials(N1min) &
			    - DDLogSemiFactorials(N1min + k) &
			    - DDLogFactorials(N2min) &
			    - DDLogSemiFactorials(N2min + k)) &
			* (sumN1 * (2.0D0 ** DBLE(-N1min - N2min - k))) / b_0
		RETURN

	CONTAINS

		FUNCTION SumC(n1, l1, n2, l2, k, M1, y)
#if(USE_QUADRUPLE==1)
			REAL(KIND = 16) :: SumC
#else
			DOUBLE PRECISION :: SumC
#endif
			INTEGER, INTENT(IN) :: n1, l1, n2, l2, k, M1
			DOUBLE PRECISION, INTENT(IN) :: y

			INTEGER :: N, Nmin, Nmax

			Nmax = n1 + n2 + ((l1 + l2 - k) / 2)
	 		Nmin = MIN_5N(n1, l1, n2, l2, k)

			SumC = SymKumar_get(n1, l1, n2, l2, Nmax, k)
			DO N = Nmax - 1, Nmin, -1
				SumC = - (SumC * (DBLE(N + M1 + k) + 0.5D0) / DBLE(N + 1) / (DBLE(N + k) + 1.5D0) / y) &
					+ SymKumar_get(n1, l1, n2, l2, N, k)
			END DO
			RETURN
		END FUNCTION SumC

	END FUNCTION ICoulombHO

	!-----------------------------------------------------------------------!
	!     Radial 2-body Coulomb Integral as defined in Appendix H for the 	!
	!     spherical	GENERAL basis. Here we integrate over r1 the result of 	!
	!     the Coulomb multipole k integral over r2.				!
	!-----------------------------------------------------------------------!

	FUNCTION ICoulomb(na, la, nb, lb, nc, lc, nd, ld, k)
		DOUBLE PRECISION :: ICoulomb

		INTEGER, INTENT(IN) :: na, la, nb, lb, nc, lc, nd, ld, k

		INTEGER :: IndexBra, IndexKet, i
		DOUBLE PRECISION :: Vcou, res
		DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)

		IndexBra = IndexVecNL(na,la)
		IndexKet = IndexVecNL(nc,lc)

		IF (IndexBra .EQ. 0 .OR. IndexKet .EQ. 0) THEN
			ICoulomb = 0.0D0
			RETURN
		END IF

	       	ALLOCATE(Integrand(1:Npoint))

!$OMP Parallel Default(None) &
!$OMP& SHARED(Npoint,nb,lb,nd,ld,k,WaveFun,IndexBra,IndexKet,Integrand) &
!$OMP& PRIVATE(i,Vcou)
!$OMP DO SCHEDULE(DYNAMIC)
		DO i=1,Npoint
			Vcou = Integral1(nb, lb, nd, ld, k, i)
			Integrand(i) = WaveFun(i,IndexBra)*WaveFun(i,IndexKet)*Vcou
		END DO
!$OMP End Do
!$OMP End Parallel

		CALL simps(Integrand,Npoint,MeshStep,res)

		DEALLOCATE(Integrand)

		Icoulomb = FOUR_PI * res / DBLE(2*k+1)

		RETURN

	CONTAINS

		!-----------------------------------------------------------------------!
		! 									!
		!   Radial 1-body Coulomb Integral as defined in Appendix H for the 	!
		!   spherical GENERAL basis. Here we integrate the coulomb multipole	!
		!   k over r2. We split the integration in 2 parts, r2 < r1 and r2 > r1.!
		!   The Coulomb multipole is defined as:				!
		!									!
		!          W_k(r1, r2) = 4*pi/(2k+1) * r_min^k/r_max^(k+1)		!
		!									!
		!   with: r_min = min(r1, r2) and r_max = max(r1, r2)			!
		! 									!
		!-----------------------------------------------------------------------!

		FUNCTION Integral1(n1, l1, n2, l2, k, index)
#if(USE_QUADRUPLE==1)
			REAL(KIND = 16) :: Integral1
#else
			DOUBLE PRECISION :: Integral1
#endif

			INTEGER, INTENT(IN) :: n1, n2, l1, l2, k, index

			INTEGER :: IndexBra, IndexKet, i
			DOUBLE PRECISION :: r_min, r_max, Vcou, res
			DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)

			IndexBra = IndexVecNL(n1,l1)
			IndexKet = IndexVecNL(n2,l2)

			! Result equal to zero if the quantum numbers are beyond the limits

			IF (IndexBra .EQ. 0 .OR. IndexKet .EQ. 0) THEN
				Integral1 = 0.0D0
				RETURN
			END IF

	       		ALLOCATE(Integrand(1:Npoint))

			DO i=1,Npoint
				Integrand(i) = 0.0D0
			END DO

			! If r2 < r1 (<=> index < i), we integrate a certain expression...

			IF (index .GT. 1) THEN

				r_max = RadMesh(index)

				DO i=1,index-1
					r_min = RadMesh(i)
					Vcou = r_min**k / r_max**(k+1)
					Integrand(i) = WaveFun(i,IndexBra)*WaveFun(i,IndexKet)*Vcou
				END DO

			END IF

			! If r2 > r1 (<=> index > i), we integrate a slightly different expression...

			r_min = RadMesh(index)

			DO i=index, Npoint
				r_max = RadMesh(i)
				Vcou = r_min**k / r_max**(k+1)
				Integrand(i) = WaveFun(i,IndexBra)*WaveFun(i,IndexKet)*Vcou
			END DO

			CALL simps(Integrand,Npoint,MeshStep,res)

			DEALLOCATE(Integrand)

			Integral1 = res
			RETURN
		END FUNCTION Integral1

	END FUNCTION ICoulomb

END MODULE ic
