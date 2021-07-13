MODULE angmom

	USE input
	USE lgfactor
	USE math

	IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: TJ

	PRIVATE TJ

CONTAINS

	FUNCTION DELTA(a, b, c)
		DOUBLE PRECISION DELTA
		INTEGER, INTENT(IN) :: a, b, c

		DELTA = EXP(0.5 * &
			(DDLogFactorials((a + b - c) / 2) + &
			 DDLogFactorials((a + c - b) / 2) + &
			 DDLogFactorials((b + c - a) / 2) - &
			 DDLogFactorials((a + b + c + 2 ) / 2)))
		RETURN
	END FUNCTION DELTA

	SUBROUTINE ThreeJSymbols_new

		INTEGER u1, u2, u3

		ALLOCATE(TJ (0:(2 * N_0), 0:(2 * N_0), 0:(2 * N_0)))
		DO u1 = 0, 2 * N_0
			DO u2 = 0, u1
				DO u3 = 0, u2
					TJ(u1, u2, u3) = 0.0
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE

	FUNCTION ThreeJSymbols_get(j_1, j2, j_3)
		DOUBLE PRECISION ThreeJSymbols_get
		INTEGER, INTENT(IN) :: j_1, j2, j_3

		INTEGER j1, j3
		INTEGER u1, u2, u3, p

		ThreeJSymbols_get = 0.0

		p = (j_1 + j2 + j_3) / 2
		IF (MOD(p, 2) .EQ. 1) RETURN ! is odd?

		IF (j_3 .GT. j_1) THEN
			j1 = j_3
			j3 = j_1
		ELSE
			j1 = j_1
			j3 = j_3
		END IF

		ThreeJSymbols_get = TJ(j1 / 2, j3 / 2, (j2 - ABS(j1 - j3)) / 2)
		IF (ThreeJSymbols_get .NE. 0.0) RETURN

		u1 = ( j1 + j2 - j3) / 2
		u2 = ( j1 - j2 + j3) / 2
		u3 = (-j1 + j2 + j3) / 2
		ThreeJSymbols_get = PAR(p / 2) * DELTA(j1, j2, j3) * &
			EXP(DDLogFactorials(p / 2) &
				- (DDLogFactorials((p - j1) / 2) &
				+ DDLogFactorials((p - j2) / 2) &
				+ DDLogFactorials((p - j3) / 2)))

		TJ(j1 / 2, j3 / 2, (j2 - ABS(j1 - j3)) / 2) = ThreeJSymbols_get
		RETURN
	END FUNCTION ThreeJSymbols_get

END MODULE angmom
