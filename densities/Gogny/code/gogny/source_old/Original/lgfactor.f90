MODULE lgfactor

	IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lnf, lng
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

	SUBROUTINE DDLogSemiFactorials_new(imax)
		INTEGER, INTENT(IN) :: imax

		INTEGER i

		ALLOCATE(ddlng(0:imax))
		ddlng_max = imax
!		Neper log of (1/2)! = SQRT(PI) / 2
		ddlng(0) = -0.1207822376352452223455184
		DO i = 1, imax
			ddlng(i) = LOG(DBLE(i) + 0.5) + ddlng(i - 1)
		END DO
		RETURN
	END SUBROUTINE DDLogSemiFactorials_new

	FUNCTION DDLogFactorials(i)
		REAL(KIND = 16) DDLogFactorials
		INTEGER, INTENT(IN) :: i

		IF ((i .GE. 0) .AND. (i .LE. ddlnf_max)) THEN
			DDLogFactorials = ddlnf(i)
		ELSE
!			DDLogFactorials = gammln(i + 1.0)
			PRINT *, "Fuera de rango:", i, ddlnf_max
			STOP "DDLogFactorials"
		END IF
		RETURN
	END FUNCTION DDLogFactorials

	FUNCTION DDLogSemiFactorials(i)
		REAL(KIND = 16) DDLogSemiFactorials
		INTEGER, INTENT(IN) :: i

!		Neper logarithm of (-1/2)! = SQRT(PI)
		IF (i .EQ. -1) THEN
			DDLogSemiFactorials = 0.5723649429247000870717137
		ELSE IF ((i .GE. 0) .AND. (i .LE. ddlng_max)) THEN
			DDLogSemiFactorials = ddlng(i)
		ELSE
!			DDLogSemiFactorials = gammln(i + 1.5)
			PRINT *, "Fuera de rango:", i, ddlng_max
			STOP "DDLogSemiFactorials"
		END IF
		RETURN
	END FUNCTION DDLogSemiFactorials

END MODULE lgfactor
