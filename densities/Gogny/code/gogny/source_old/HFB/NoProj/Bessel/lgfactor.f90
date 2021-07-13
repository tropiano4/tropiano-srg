 MODULE lgfactor

	IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lnf, lng, GammaFunc
	REAL(KIND = 16), DIMENSION(:), ALLOCATABLE :: ddlnf, ddlng

	INTEGER lnf_max, lng_max, ddlnf_max, ddlng_max

	PRIVATE lnf, lng, ddlnf, ddlng
	PRIVATE lnf_max, lng_max, ddlnf_max, ddlng_max

 CONTAINS

	SUBROUTINE LogFactorials_new(imax)
		INTEGER, INTENT(IN) :: imax

		INTEGER i

		ALLOCATE(lnf(0:imax))
		lnf_max = imax
		lnf(0) = 0.0
		DO i = 1, imax
			lnf(i) = LOG(DBLE(i)) + lnf(i - 1)
		END DO
		RETURN
	END SUBROUTINE LogFactorials_new

	SUBROUTINE LogSemiFactorials_new(imax)
		INTEGER, INTENT(IN) :: imax

		INTEGER i

		ALLOCATE(lng(0:imax))
		lng_max = imax
		lng(0) = -0.1207822376352452223455184
		DO i = 1, imax
			lng(i) = LOG(DBLE(i) + 0.5) + lng(i - 1)
		END DO
		RETURN
	END SUBROUTINE LogSemiFactorials_new

	! Functions giving Gamma(n+1/2)

	SUBROUTINE GammaFunction_new(imax)
		INTEGER, INTENT(IN) :: imax

		INTEGER i
                DOUBLE PRECISION :: pi

		ALLOCATE(GammaFunc(0:imax))
                pi = 4.0d0*ATAN(1.0d0)
		GammaFunc(0) = LOG(SQRT(pi))
		DO i = 1, imax
			GammaFunc(i) = LOG(DBLE(2*i-1)) - LOG(2.0d0) + GammaFunc(i - 1)
		END DO
		RETURN
	END SUBROUTINE GammaFunction_new

	FUNCTION GammaFunction(i)
		DOUBLE PRECISION :: GammaFunction
		INTEGER, INTENT(IN) :: i

		IF ((i .GE. 0) .AND. (i .LE. ddlnf_max)) THEN
			GammaFunction = GammaFunc(i)
		ELSE
			WRITE(*,'("Beyond range in GammaFunction - n = ",i8," nmax = ",I8)')  i, ddlnf_max
			STOP "GammaFunction"
		END IF
		RETURN
	END FUNCTION GammaFunction

	! Functions giving Gamma(n+1/2)

	SUBROUTINE DDLogFactorials_new(imax)
		INTEGER, INTENT(IN) :: imax

		INTEGER i

		ALLOCATE(ddlnf(0:imax))
		ddlnf_max = imax
		ddlnf(0) = 0.0
		DO i = 1, imax
			ddlnf(i) = LOG(DBLE(i)) + ddlnf(i - 1)
		END DO
		RETURN
	END SUBROUTINE DDLogFactorials_new

	FUNCTION DDLogFactorials(i)
		REAL(KIND = 16) DDLogFactorials
		INTEGER, INTENT(IN) :: i

		IF ((i .GE. 0) .AND. (i .LE. ddlnf_max)) THEN
			DDLogFactorials = ddlnf(i)
		ELSE
			WRITE(*,'("Beyond range in DDLogFactorials - n = ",i8," nmax = ",I8)')  i, ddlnf_max
			STOP "DDLogFactorials"
		END IF
		RETURN
	END FUNCTION DDLogFactorials

	! Functions giving ln( (n + 3/2)! )

	SUBROUTINE DDLogSemiFactorials_new(imax)
		INTEGER, INTENT(IN) :: imax

		REAL(KIND = 16) pi
		INTEGER i

		ddlng_max = imax
		
		ALLOCATE(ddlng(0:imax))
		
		!ddlng(0) = -0.1207822376352452223455184
		
		! For n=0, we have (3/2)! = SQRT(PI) / 2 (Abramovitz, 6.1.9)
		! We calculate ln ( (3/2)! ) at machine precision
		
		pi = 4.0*ATAN(1.0)
		ddlng(0) = LOG(0.5*SQRT(pi))
		
		DO i = 1, imax
			ddlng(i) = LOG(DBLE(i) + 0.5) + ddlng(i - 1)
		END DO
		
		RETURN
	END SUBROUTINE DDLogSemiFactorials_new

	FUNCTION DDLogSemiFactorials(i)
		REAL(KIND = 16) DDLogSemiFactorials
		REAL(KIND = 16) pi
		INTEGER, INTENT(IN) :: i

		! For n=-1, we have (-1/2)! = SQRT(PI) (Abramovitz, 6.1.8)
		! We calculate ln ( (1/2)! ) at machine precision
		
		IF (i .EQ. -1) THEN
			pi = 4.0*ATAN(1.0)
			DDLogSemiFactorials = LOG(SQRT(pi))		
			!DDLogSemiFactorials = 0.5723649429247000870717137
		ELSE IF ((i .GE. 0) .AND. (i .LE. ddlng_max)) THEN
			DDLogSemiFactorials = ddlng(i)
		ELSE
			WRITE(*,'("Beyond range in DDLogSemiFactorials - n = ",i8," nmax = ",I8)')  i, ddlng_max
			STOP "DDLogSemiFactorials"
		END IF
		RETURN
	END FUNCTION DDLogSemiFactorials
	
END MODULE lgfactor
