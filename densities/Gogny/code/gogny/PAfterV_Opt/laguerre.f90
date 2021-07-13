 MODULE laguerre

 	USE gauss
 	USE lgfactor

 CONTAINS

        SUBROUTINE RadialWaveHO(n, l, x, Rwave, boscil)
        	DOUBLE PRECISION, INTENT(OUT) :: Rwave
        	DOUBLE PRECISION, INTENT(IN) :: x, boscil
        	INTEGER, INTENT(IN) :: n, l

        	INTEGER :: order, ArgGamma

        	DOUBLE PRECISION :: NormSquared, Norm
        	DOUBLE PRECISION :: chi, alpha, Poly, b, argPol

        	b = 1.0/boscil
        	order = n - 1

        	ArgGamma = order + l + 1
        	NormSquared = 2.0d0*b**3*EXP(DDLogFactorials(order)-GammaFunction(ArgGamma))
        	Norm = SQRT(NormSquared)

        	chi = b*x
        	alpha = REAL(l) + 0.5d0

        	argPol = chi**2
        	CALL GeneralizedLaguerre(Poly, argPol, order, alpha)

        	Rwave = Norm*EXP(-0.5*chi**2)*chi**l*Poly

        RETURN
        END SUBROUTINE RadialWaveHO


        SUBROUTINE GeneralizedLaguerre(lagPol, x, order, alpha)
        	DOUBLE PRECISION, INTENT(OUT) :: lagPol
        	INTEGER, INTENT(IN) :: order
        	DOUBLE PRECISION, INTENT(IN) :: x, alpha

        	INTEGER :: i

        	DOUBLE PRECISION :: p0, p1, p2

        	DOUBLE PRECISION :: b(order), c(order)

        	p1 = 1.0D+00

        	IF (order .EQ. 0) THEN
        	        lagPol = p1
        	        RETURN
        	END IF

        	p2 = - x + alpha + 1.0D+00

        	IF (order .EQ. 1) THEN
        	        lagPol = p2
        	        RETURN
        	END IF

        	DO i = 1, order
        	        b(i) = + (2.0d0*REAL(i) + alpha + 1.0d0 - x)/REAL(i + 1)
        	        c(i) = - (REAL(i) + alpha)/(REAL(i + 1))
        	END DO

        	DO i = 2, order
        	        p0 = p1
        	        p1 = p2
        	        p2 = b(i-1)*p1 + c(i-1)*p0
        	END DO

        	lagPol = p2

        RETURN
        END SUBROUTINE GeneralizedLaguerre

END MODULE laguerre

