!-----------------------------------------------------------------------!
!									!
!									!
!     RADIAL INTEGRALS: BRINK-BOKER FORCE				!
!									!
!  This module computes the radial integral of the Brink-Boeker force 	!
!  in the 2 distinct cases of an "analytical" harmonic oscillator basis !
!  and a general Woods-Saxon (or other) basis.				!
!									!
!-----------------------------------------------------------------------!

 MODULE bessik

	IMPLICIT NONE

 CONTAINS

	!---------------------------------------------------------------------
	!     modified bessel function of fractional order
	!     xnu = order of the Bessel function nu
	!     ri   I_nu(x)  first kind
	!     rk   K_nu(x)  second kind
	!     rip  derivative of I(x)
	!     rkp  derivative of K(x)
	!     from Numerical Recipes.
	!
	!---------------------------------------------------------------------

	SUBROUTINE bessel(x, xnu, ri, rk, rip, rkp)
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16), INTENT(IN) :: x, xnu
		REAL(KIND = 16), INTENT(OUT) :: ri, rip, rk, rkp
#else
		DOUBLE PRECISION, INTENT(IN) :: x, xnu
		DOUBLE PRECISION, INTENT(OUT) :: ri, rip, rk, rkp
#endif

		INTEGER :: MAXIT, i, l, nl

#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: EPS,FPMIN,PI,XMIN
		REAL(KIND = 16) :: a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff
		REAL(KIND = 16) :: gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1
		REAL(KIND = 16) :: ripl,ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
#else
		DOUBLE PRECISION :: EPS,FPMIN,PI,XMIN
		DOUBLE PRECISION :: a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff
		DOUBLE PRECISION :: gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1
		DOUBLE PRECISION :: ripl,ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
#endif

		PARAMETER (EPS=1.D-16, FPMIN=1.D-30 ,MAXIT=10000, XMIN=2.0D0)

		PI = 4.0D0*ATAN(1.0D0)

		IF (x.LE.0.0D0 .OR. xnu.LT.0.0D0) THEN
		   WRITE(*,'("x = ",E24.16," nu = ",F10.5)') x,xnu
		   stop "Bad arguments in bessik"
		END IF

		nl = int(xnu + 0.5D0)
		xmu = xnu - nl
		xmu2 = xmu*xmu
		xi = 1.0D0/x
		xi2 = 2.0D0*xi
		h = xnu*xi

		IF (h .LT. FPMIN) h = FPMIN

		b = xi2*xnu
		d = 0.0D0
		c = h

		DO i = 1,MAXIT

			b = b + xi2
			d = 1.0D0/(b + d)
			c = b + 1.0D0/c
			del = c*d
			h = del*h

			IF (ABS(del - 1.0D0) .LT. EPS) EXIT

		END DO

		IF (i .EQ. MAXIT) STOP 'x too large in bessik; try asymptotic expansion'

		ril  = FPMIN
		ripl = h*ril
		ril1 = ril
		rip1 = ripl
		fact = xnu*xi

		DO l = nl, 1, -1
			ritemp = fact*ril + ripl
			fact = fact - xi
			ripl = fact*ritemp + ril
			ril = ritemp
		END DO

		f = ripl/ril

		IF (x .LT. XMIN) THEN

			x2 = 0.5D0*x
			pimu = PI*xmu

			IF (ABS(pimu) .LT. EPS) THEN
				fact = 1.0D0
			ELSE
				fact = pimu/SIN(pimu)
			END IF

			d = -LOG(x2)
			e = xmu*d

			IF (ABS(e) .LT. EPS) THEN
				fact2 = 1.0D0
			ELSE
				fact2 = sinh(e)/e
			END IF

			CALL beschb(xmu, gam1, gam2, gampl, gammi)

			ff = fact*(gam1*COSH(e) + gam2*fact2*d)
			sum = ff
			e = EXP(e)
			p = 0.5D0*e/gampl
			q = 0.5D0/(e*gammi)
			c = 1.0D0
			d = x2*x2
			sum1 = p

			DO i = 1,MAXIT
	!!!!!write(*,'("x < XMIN - i = ",i5)') i
				ff = (i*ff + p + q)/(i*i - xmu2)
				c = c*d/i
				p = p/(i - xmu)
				q = q/(i + xmu)
				del = c*ff
				sum = sum + del
				del1 = c*(p - i*ff)
				sum1 = sum1 + del1
				IF (abs(del) .LT. ABS(sum)*EPS) EXIT
			END DO

			IF (i .EQ. MAXIT) STOP 'bessk series failed to converge'

			rkmu = sum
			rk1 = sum1*xi2

		ELSE

			b = 2.0D0*(1.0D0 + x)
			d = 1.0D0/b
			delh = d
			h = delh
			q1 = 0.0D0
			q2 = 1.0D0
			a1 = 0.25D0 - xmu2
			c = a1
			q = c
			a = -a1
			s = 1.0D0 + q*delh

			DO i = 2,MAXIT
	!write(*,'("x > XMIN - i = ",i5)') i
				a = a - 2*(i-1)
				c = -a*c/i
				qnew = (q1 - b*q2)/a
				q1 = q2
				q2 = qnew
				q = q + c*qnew
				b = b + 2.0D0
				d = 1.0D0/(b + a*d)
				delh = (b*d - 1.0D0)*delh
				h = h + delh
				dels = q*delh
				s = s + dels
				IF (ABS(dels/s) .LT. EPS) EXIT
			END DO

			IF (i .EQ. MAXIT) STOP 'bessik: failure to converge in cf2'

			h = a1*h
			rkmu = SQRT(PI/(2.0D0*x))*exp(-x)/s
			rk1 = rkmu*(xmu + x + 0.5D0 - h)*xi

		END IF

		rkmup = xmu*xi*rkmu - rk1
		rimu = xi/(f*rkmu - rkmup)
		ri = (rimu*ril1)/ril
		rip = (rimu*rip1)/ril

		DO i=1,nl
			rktemp = (xmu + i)*xi2*rk1 + rkmu
			rkmu = rk1
			rk1 = rktemp
		END DO

		rk = rkmu
		rkp = xnu*xi*rkmu - rk1

		RETURN
	END SUBROUTINE bessel


	SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)

		INTEGER :: NUSE1,NUSE2
		PARAMETER (NUSE1=7,NUSE2=8)

#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: gam1, gam2, gammi, gampl, x, xx, one
		REAL(KIND = 16), POINTER, DIMENSION(:):: c1, c2
#else
		DOUBLE PRECISION :: gam1, gam2, gammi, gampl, x, xx, one
		DOUBLE PRECISION, POINTER, DIMENSION(:):: c1, c2
#endif

		ALLOCATE(c1(7))

		one = 1.0D0

		c1(1) = -1.142022680371168D0
		c1(2) =  6.5165112670737D-3
		c1(3) =  3.087090173086D-4
		c1(4) = -3.4706269649D-6
		c1(5) =  6.9437664D-9
		c1(6) =  3.67795D-11
		c1(7) = -1.356D-13

		ALLOCATE(c2(8))

		c2(1) =  1.843740587300905D0
		c2(2) = -7.68528408447867D-2
		c2(3) =  1.2719271366546D-3
		c2(4) = -4.9717367042D-6
		c2(5) = -3.31261198D-8
		c2(6) =  2.423096D-10
		c2(7) = -1.702D-13
		c2(8) = -1.49D-15

		xx = 8.0D0*x*x - 1.0D0

		gam1 = chebev(-one, one, c1, NUSE1, xx)
		gam2 = chebev(-one, one, c2, NUSE2, xx)

		DEALLOCATE(c1)
		DEALLOCATE(c2)

		gampl = gam2 - x*gam1
		gammi = gam2 + x*gam1

		RETURN
	END SUBROUTINE beschb

	FUNCTION chebev(a, b, c, m, x)
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: chebev
#else
		DOUBLE PRECISION :: chebev
#endif

		INTEGER, INTENT(IN) :: m

#if(USE_QUADRUPLE==1)
		REAL(KIND = 16), INTENT(IN) ::a, b, x
#else
		DOUBLE PRECISION, INTENT(IN) ::a, b, x
#endif

#if(USE_QUADRUPLE==1)
		REAL(KIND = 16), POINTER, DIMENSION(:) :: c
#else
		DOUBLE PRECISION, POINTER, DIMENSION(:) :: c
#endif

		INTEGER :: j

#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: d, dd, sv, y, y2
#else
		DOUBLE PRECISION :: d, dd, sv, y, y2
#endif

		IF ((x - a)*(x - b) .GT. 0.0D0) THEN
			write(*,'("x = ",E24.16," a = ",E24.16," b = ",E24.16)') x,a,b
			STOP 'x not in range in chebev'
		END IF

		d = 0.0D0
		dd = 0.0D0
		y = (2.0D0*x - a - b)/(b - a)
		y2 = 2.0D0*y

		DO j = m, 2, -1
			sv = d
			d = y2*d - dd + c(j)
			dd = sv
		END DO

		chebev = y*d - dd + 0.5D0*c(1)

		RETURN
	END FUNCTION chebev

END MODULE bessik
