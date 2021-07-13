!---------------------------------------------------------------------!
!                                                                     !
!          MODULE CONTAINING USEFUL MATHEMATICAL FUNCTIONS            !
!                                                                     !
!---------------------------------------------------------------------!

 MODULE math

	USE input

	IMPLICIT NONE

 CONTAINS

        !---------------------------------------------------------------------!
        !   This routine performs the composite Simpson's rule integration    !
        !   of a function f defined by a table of n equispaced values.        !
        !    								      !
        !                        See: Koonin, Computational Physics, p.9      !
        !    								      !
        !    The parameters are:  					      !
        !     f = Array of values of the function f(x)			      !
        !     n = Number of points x_k			      		      !
        !     h = The uniform spacing between x values: h = x_k+1 - x_k       !
        !    result = Estimate of the integral that is returned to caller.    !
        !---------------------------------------------------------------------!

	SUBROUTINE simps(functi,npoint,step,res)

        DOUBLE PRECISION, INTENT(IN) :: step
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: functi
        DOUBLE PRECISION, INTENT(OUT) :: res

        INTEGER :: npanel, npoint, nhalf, nbegin, nend, i

        DOUBLE PRECISION :: x

        ! Check to see if number of panels is even.  Number of panels is n - 1.

	nbegin = 1
	!nbegin = 0

	npanel = npoint - nbegin
	nhalf  = npanel/2

	res = 0.0D0

	! Number of panels is odd.  Use Simpson's 3/8 rule on first three panels, 1/3 rule on rest of them.

	IF ((npanel-2*nhalf).NE.0) THEN

	     res = 3.0D0*step*(functi(nbegin) + 3.0D0*(functi(nbegin+1)+functi(nbegin+2)) + functi(nbegin+3))/8.0D0

	     IF ((npoint-nbegin).EQ.3) RETURN

	     nbegin=nbegin+3

	END IF

	! Apply 1/3 rule - add in first, second, last values

	res = res + step*(functi(nbegin) + 4.0D0*functi(nbegin+1) + functi(npoint))/3.0D0
	nbegin = nbegin+2

	IF (nbegin.EQ.npoint) THEN
	    RETURN
	ELSE

		x = 0.0D0
		nend = npoint - 1

		DO i = nbegin,nend,2
			x = x + functi(i) + 2.0D0*functi(i+1)
		END DO

		res = res + 2.0D0*step*x/3.0D0

		RETURN

	END IF

        END SUBROUTINE simps

	!
	!                 Newton's method subroutine
	!
	!  This routine calculates the zeros of a function Y(x) by Newton"s method.
	!  The routine requires an initial guess, x0, and a convergence factor, e.
	!  Also required is a limit on the number of iterations, m.
	!

	SUBROUTINE newton(fonction, Nmax, Nstep, accuracy, x0, f0, Maximum)
		DOUBLE PRECISION, EXTERNAL :: fonction
		DOUBLE PRECISION :: accuracy, x0, f0, derivF, f_test, deriv_test
		INTEGER :: Nmax, Nstep, Maximum

		Nstep = 0

		! Get the value of the fonction and its derivativ at initial guess x0
		f0=fonction(x0,derivF)

		DO WHILE ( DABS(f0/derivF) > accuracy )

			f0 = fonction(x0,derivF)

			! Update estimate
			x0 = x0 - (f0/derivF)

                        !IF ( x0 <= 0.0 ) THEN
                        !        x0 = x0 + 2.0*(f0/derivF)
                        !END IF

                        IF (Maximum == 1) THEN
                                f_test = fonction(x0,deriv_test)
                                IF (f_test <= f0) x0 = x0 + 2.0*(f0/derivF)
                        END IF

			Nstep = Nstep+1

			IF ( Nstep >= Nmax ) THEN
                                WRITE(*,'("Warning - Number of iterations exceeded in Newton")')
                                RETURN
                        END IF

		END DO
	RETURN
	END SUBROUTINE newton

END MODULE math
