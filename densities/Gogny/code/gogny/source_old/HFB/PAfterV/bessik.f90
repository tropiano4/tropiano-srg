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

 MODULE bessik

	IMPLICIT NONE

 CONTAINS

	!---------------------------------------------------------------------
	!     modified besselfunction of fractional order
	!     xnu = order of the Bessel function nu
	!     ri   I_nu(x)  first kind
	!     rk   K_nu(x)  second kind
	!     rip  derivative of I(x)
	!     rkp  derivative of K(x)
	!     from Numerical Recipes.
	!
	!---------------------------------------------------------------------
	
	SUBROUTINE bessel(x, xnu, ri, rk, rip, rkp)
		REAL(KIND = 16), INTENT(IN) :: x, xnu
		REAL(KIND = 16), INTENT(OUT) :: ri, rip, rk, rkp
		!DOUBLE PRECISION, INTENT(IN) :: x
		!DOUBLE PRECISION, INTENT(OUT) :: ri, rip, rk, rkp

		INTEGER :: MAXIT, i, l, nl

		!DOUBLE PRECISION :: EPS,FPMIN,PI,XMIN
		!DOUBLE PRECISION :: a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff
		!DOUBLE PRECISION :: gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1
		!DOUBLE PRECISION :: ripl,ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
	
		REAL(KIND = 16) :: EPS,FPMIN,PI,XMIN
		REAL(KIND = 16) :: a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff
		REAL(KIND = 16) :: gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1
		REAL(KIND = 16) :: ripl,ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
	
		PARAMETER (EPS=1.e-16, FPMIN=1.e-30 ,MAXIT=10000, XMIN=2.0)
	
		PI = 4.0*ATAN(1.0)
	
		IF (x.LE.0.0 .OR. xnu.LT.0.0) THEN
		   WRITE(*,'("x = ",E24.16," nu = ",F10.5)') x,xnu
		   stop "Bad arguments in bessik"
		END IF
	
		nl = int(xnu + 0.5)
		xmu = xnu - nl
		xmu2 = xmu*xmu
		xi = 1.0/x
		xi2 = 2.0*xi
		h = xnu*xi
	
		IF (h .LT. FPMIN) h = FPMIN
	
		b = xi2*xnu
		d = 0.0
		c = h
	
		DO i = 1,MAXIT
	
			b = b + xi2
			d = 1.0/(b + d)
			c = b + 1.0/c
			del = c*d
			h = del*h
		
			IF (ABS(del - 1.0) .LT. EPS) EXIT
			
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
	
			x2 = 0.5*x
			pimu = PI*xmu
		
			IF (ABS(pimu) .LT. EPS) THEN
				fact = 1.0
			ELSE
				fact = pimu/SIN(pimu)
			END IF
		
			d = -LOG(x2)
			e = xmu*d
		
			IF (ABS(e) .LT. EPS) THEN
				fact2 = 1.0
			ELSE
				fact2 = sinh(e)/e
			END IF
		
			CALL beschb(xmu, gam1, gam2, gampl, gammi)
		
			ff = fact*(gam1*COSH(e) + gam2*fact2*d)
			sum = ff
			e = EXP(e)
			p = 0.5*e/gampl
			q = 0.5/(e*gammi)
			c = 1.0
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
	
			b = 2.0*(1.0 + x)
			d = 1.0/b
			delh = d
			h = delh
			q1 = 0.0
			q2 = 1.0
			a1 = 0.25 - xmu2
			c = a1
			q = c
			a = -a1
			s = 1.0 + q*delh
		
			DO i = 2,MAXIT
	!write(*,'("x > XMIN - i = ",i5)') i
				a = a - 2*(i-1)
				c = -a*c/i
				qnew = (q1 - b*q2)/a
				q1 = q2
				q2 = qnew
				q = q + c*qnew
				b = b + 2.0
				d = 1.0/(b + a*d)
				delh = (b*d - 1.0)*delh
				h = h + delh
				dels = q*delh
				s = s + dels
				IF (ABS(dels/s) .LT. EPS) EXIT
			END DO
		
			IF (i .EQ. MAXIT) STOP 'bessik: failure to converge in cf2'

			h = a1*h
			rkmu = SQRT(PI/(2.0*x))*exp(-x)/s
			rk1 = rkmu*(xmu + x + 0.5 - h)*xi
		
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
		
		REAL(KIND = 16) :: gam1, gam2, gammi, gampl, x, xx, one
		REAL(KIND = 16), POINTER, DIMENSION(:):: c1, c2
		!DOUBLE PRECISION :: gam1, gam2, gammi, gampl, x, xx
		!DOUBLE PRECISION, POINTER, DIMENSION(:):: c1, c2

		ALLOCATE(c1(7))

		one = 1.0

		c1(1) = -1.142022680371168e0
		c1(2) =  6.5165112670737e-3
		c1(3) =  3.087090173086e-4
		c1(4) = -3.4706269649e-6
		c1(5) =  6.9437664e-9
		c1(6) =  3.67795e-11
		c1(7) = -1.356e-13
		
		ALLOCATE(c2(8))

		c2(1) =  1.843740587300905e0
		c2(2) = -7.68528408447867e-2
		c2(3) =  1.2719271366546e-3
		c2(4) = -4.9717367042e-6
		c2(5) = -3.31261198e-8
		c2(6) =  2.423096e-10
		c2(7) = -1.702e-13
		c2(8) = -1.49e-15
     
		xx = 8.0*x*x - 1.0
		
		gam1 = chebev(-one, one, c1, NUSE1, xx)
		gam2 = chebev(-one, one, c2, NUSE2, xx)
		
		DEALLOCATE(c1)
		DEALLOCATE(c2)
		
		gampl = gam2 - x*gam1
		gammi = gam2 + x*gam1
	
		RETURN
	END SUBROUTINE beschb

	FUNCTION chebev(a, b, c, m, x)
		REAL(KIND = 16) :: chebev
		!DOUBLE PRECISION :: chebev

		INTEGER, INTENT(IN) :: m
		REAL(KIND = 16), INTENT(IN) ::a, b, x
		!DOUBLE PRECISION, INTENT(IN) ::a, b, x
		
		REAL(KIND = 16), POINTER, DIMENSION(:) :: c
		!DOUBLE PRECISION, POINTER, DIMENSION(:) :: c

		INTEGER :: j
		REAL(KIND = 16) :: d, dd, sv, y, y2
		!DOUBLE PRECISION :: d, dd, sv, y, y2
		
		IF ((x - a)*(x - b) .GT. 0.0) THEN
			write(*,'("x = ",E24.16," a = ",E24.16," b = ",E24.16)') x,a,b
			STOP 'x not in range in chebev'
		END IF
		
		d = 0.0
		dd = 0.0
		y = (2.0*x - a - b)/(b - a)
		y2 = 2.0*y
	
		DO j = m, 2, -1
			sv = d
			d = y2*d - dd + c(j)
			dd = sv
		END DO
	
		chebev = y*d - dd + 0.5*c(1)
		
		RETURN
	END FUNCTION chebev

END MODULE bessik
