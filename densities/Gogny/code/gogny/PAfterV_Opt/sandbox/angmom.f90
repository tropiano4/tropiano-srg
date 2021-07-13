!---------------------------------------------------------------------!
!                                                                     !
!     MODULE TO CALCULATE CLEBSCH-GORDAN COEFFICIENTS (3J-SYMBOLS)    !
!                                                                     !
!---------------------------------------------------------------------!

 MODULE angmom

	USE input
	USE lgfactor
	USE math

	IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: TJ

	PRIVATE TJ

 CONTAINS

        ! Function calculating the Delta(abc) function as defined in
	! Varshalovitch, Sec. 8.2, Eq. (1), page 237

	FUNCTION DELTA(a, b, c)
		DOUBLE PRECISION DELTA
		INTEGER, INTENT(IN) :: a, b, c

		DELTA = EXP(0.5D0 * &
			(DDLogFactorials((a + b - c) / 2) + &
			 DDLogFactorials((a + c - b) / 2) + &
			 DDLogFactorials((b + c - a) / 2) - &
			 DDLogFactorials((a + b + c + 2 ) / 2)))
		RETURN
	END FUNCTION DELTA

        !  "Create" a new 3j-symbol by allocating the required memory
	!  and initializing all coefficients to 0

	SUBROUTINE ThreeJSymbols_new

		ALLOCATE(TJ (0:(2*Lmax), 0:(2*Lmax), 0:(2*Lmax)))
		TJ(:,:,:) = 0.0D0

		RETURN
	END SUBROUTINE

        !  Calculating the 3j-symbol for all m equal to 0

	FUNCTION ThreeJSymbols_get(j_1, j2, j_3)
		DOUBLE PRECISION ThreeJSymbols_get
		INTEGER, INTENT(IN) :: j_1, j2, j_3

		INTEGER j1, j3
		INTEGER u1, u2, u3, p

		ThreeJSymbols_get = 0.0D0

		p = (j_1 + j2 + j_3) / 2
		IF (MOD(p, 2) .EQ. 1) RETURN ! is odd?

		IF (j_3 .GT. j_1) THEN
			j1 = j_3
			j3 = j_1
		ELSE
			j1 = j_1
			j3 = j_3
		END IF

		IF (j1/2 .GT. 2*Lmax .OR. j3/2 .GT. 2*Lmax .OR. (j2 - ABS(j1 - j3))/2 .GT. 2*Lmax) THEN
			STOP "Too large j1, j2 or j3 in ThreeJSymbols_get"
		END IF

		ThreeJSymbols_get = TJ(j1 / 2, j3 / 2, (j2 - ABS(j1 - j3)) / 2)

                IF (ABS(ThreeJSymbols_get) .GT. 1.0D0) STOP "3J symbol greater than 1!"

		IF (ABS(ThreeJSymbols_get) .GT. 1.D-14) RETURN

		u1 = ( j1 + j2 - j3) / 2
		u2 = ( j1 - j2 + j3) / 2
		u3 = (-j1 + j2 + j3) / 2

		ThreeJSymbols_get = PAR(p / 2) * DELTA(j1, j2, j3) * &
			      EXP(DDLogFactorials(p / 2) &
				-(DDLogFactorials((p - j1) / 2) &
				+ DDLogFactorials((p - j2) / 2) &
				+ DDLogFactorials((p - j3) / 2)))

                IF (ABS(ThreeJSymbols_get) .GT. 1.0D0) STOP "3J symbol greater than 1!"

		TJ(j1 / 2, j3 / 2, (j2 - ABS(j1 - j3)) / 2) = ThreeJSymbols_get
		RETURN
	END FUNCTION ThreeJSymbols_get

END MODULE angmom
