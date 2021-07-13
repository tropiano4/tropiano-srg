MODULE gauss

	USE math

	IMPLICIT NONE

	TYPE GaussQuadrature
		DOUBLE PRECISION, DIMENSION(:), POINTER :: x, w
		INTEGER n
	END TYPE

	TYPE GaussLaguerreQuadrature
		TYPE (GaussQuadrature) gauss
		DOUBLE PRECISION alf
	END TYPE

CONTAINS

	SUBROUTINE GaussQuadrature_new(gauss, n)
		TYPE (GaussQuadrature), INTENT(INOUT) :: gauss
		INTEGER, INTENT(IN) :: n

		gauss%n = n
		ALLOCATE(gauss%x(0:n))
		IF (.NOT. ASSOCIATED(gauss%x)) STOP "Unable to allocate memory"
		ALLOCATE(gauss%w(0:n)) !TODO ¿El tamaño máximo es n?
		IF (.NOT. ASSOCIATED(gauss%w)) STOP "Unable to allocate memory"
		RETURN
	END SUBROUTINE GaussQuadrature_new

	SUBROUTINE GaussLaguerreQuadrature_new(laguerre, n, alf)
		TYPE (GaussLaguerreQuadrature), INTENT(INOUT) :: laguerre
		INTEGER, INTENT(IN) :: n
		DOUBLE PRECISION, INTENT(IN) :: alf

		CALL GaussQuadrature_new(laguerre%gauss, n)
		laguerre%alf = alf
		CALL gaulag2(laguerre%gauss%x, laguerre%gauss%w, n, alf)
		RETURN
	END SUBROUTINE GaussLaguerreQuadrature_new

	! Given alf, the parameter alf of the Laguerre polynomials,
	! this routine returns arrays x(1..n)  and w(1..n)
	! containing the abscissas and weights
	! of the n-point Gauss-Laguerre quadrature formula.
	! The smallest abscissa is returned in x(1), the largest in x(n).
	SUBROUTINE gaulag2(x, w, n, alf)
		DOUBLE PRECISION, DIMENSION(0:), INTENT(INOUT) :: x, w
		INTEGER, INTENT(IN) :: n
		DOUBLE PRECISION, INTENT(IN) :: alf

		INTEGER kind2, kpts
		DOUBLE PRECISION bet
		DOUBLE PRECISION, DIMENSION(1:2) :: endpts
		DATA endpts /0.0, 0.0/
		DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: wksp

		ALLOCATE(wksp(1:n))
		kind2 = 6
		bet  = 0.0
		kpts = 0
		CALL gaussq(kind2, n, alf, bet, kpts, endpts, wksp, x, w)
		RETURN
	END SUBROUTINE

	SUBROUTINE gaussq(kind2, n, alf, bet, kpts, endpts, b, t, w)
		INTEGER, INTENT(IN) :: kind2, n, kpts
		DOUBLE PRECISION, INTENT(IN) :: alf, bet
		DOUBLE PRECISION, DIMENSION(1:2), INTENT(IN) :: endpts
		DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: b
		DOUBLE PRECISION, DIMENSION(0:), INTENT(INOUT) :: t, w

		INTEGER ierr, i
		DOUBLE PRECISION muzero, gam, t1

		CALL clase(kind2, n, alf, bet, b, t, muzero)
		IF (kpts .NE. 0) THEN
			IF (kpts .EQ. 2) THEN
				gam = gbslve(endpts(1), n, t, b)
				t1 = DBLE(endpts(1) - endpts(2)) / (gbslve(endpts(2), n, t, b) - gam)
				b(n - 1) = SQRT(t1)
				t(n) = DBLE(endpts(1)) + gam * t1
			ELSE ! Computing 2nd power
				t(n) = gbslve(endpts(1), n, t, b) * (b(n - 1) ** 2.0) + DBLE(endpts(1))
			END IF
		END IF

		w(1) = 1.0
		DO i = 2, n
			w(i) = 0.0
		END DO

		CALL gbtql2(n, t, b, w, ierr)
		DO i = 1, n
			w(i) = muzero * w(i) * w(i)
		END DO
		RETURN
	END SUBROUTINE gaussq

	SUBROUTINE clase(kind2, n, alf, bet, b, a, muzero)
		INTEGER, INTENT(IN) :: kind2, n
		DOUBLE PRECISION, INTENT(IN) :: alf, bet
		DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: b
		DOUBLE PRECISION, DIMENSION(0:), INTENT(INOUT) :: a
		DOUBLE PRECISION, INTENT(INOUT) :: muzero

		INTEGER i
		DOUBLE PRECISION ab, a2b2, abi

		SELECT CASE (kind2)
		CASE (1)
			muzero = 2.0
			DO i = 1, n - 1
				a(i) = 0.0
				b(i) = i / SQRT(i * 4.0 * i - 1.0)
			END DO
			a(n) = 0.0
		CASE (2)
			muzero = PI
			DO i = 1, n - 1
				a(i) = 0.0
				b(i) = 0.5
			END DO
			b(1) = SQRT(0.5)
			a(n) = 0.0
		CASE (3)
			muzero = PI / 2.0
			DO i = 1, n - 1
				a(i) = 0.0
				b(i) = 0.5
			END DO
			a(n) = 0.0
		CASE (4)
			muzero = SQRT(PI)
			DO i = 1, n - 1
				a(i) = 0.0
				b(i) = SQRT(0.5 * i)
			END DO
			a(n) = 0.0
		CASE (5)
			ab = alf + bet
			abi = ab + 2.0
			muzero = (2.0 ** (ab + 1.0)) * dgamma(alf + 1.0) * dgamma(bet + 1.0) / dgamma(ab + 2.0)
			a(1) = (bet - alf) / abi
			b(1) = SQRT((alf + 1.0) * 4.0 * (bet + 1.0) / ((abi + 1.0) * abi * abi))

			a2b2 = bet * bet - alf * alf
			DO i = 2, n - 1
				abi = i * 2.0 + ab
				a(i) = a2b2 / ((abi - 2.) * abi)
				b(i) = SQRT(i * 4.0 * (i + alf) * (i + bet) * (i + ab) / ((abi * abi - 1.0) * abi * abi))
			END DO
			abi = n * 2.0 + ab
			a(n) = a2b2 / ((abi - 2.0) * abi)
		CASE (6)
			muzero = dgamma(alf + 1.0)
			DO i = 1, n - 1
				a(i) = 2.0 * i - 1.0 + alf
				b(i) = SQRT((alf + i) * i)
			END DO
			a(n) = n * 2.0 - 1.0 + alf
		END SELECT
		RETURN
	END SUBROUTINE clase

	FUNCTION gbslve(shift, n, a, b)
		DOUBLE PRECISION gbslve
		DOUBLE PRECISION, INTENT(IN) :: shift
		INTEGER, INTENT(IN) :: n
		DOUBLE PRECISION, DIMENSION(0:), INTENT(IN) :: a
		DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: b

		INTEGER i
		DOUBLE PRECISION alpha

		alpha = a(1) - shift
		DO i = 2, n - 1 ! Computing 2nd power
			alpha = a(i) - shift - (b(i - 1) ** 2.0) / alpha
		END DO
		gbslve = 1.0 / alpha
		RETURN
	END FUNCTION gbslve

	SUBROUTINE gbtql2(n, d, e, z, ierr)
		INTEGER, INTENT(IN) :: n
		DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: e
		DOUBLE PRECISION, DIMENSION(0:), INTENT(INOUT) :: d, z
		INTEGER, INTENT(INOUT) :: ierr

		DOUBLE PRECISION b, c, f, g, p, r, s
		INTEGER i, ii, j, k, l, m, m2

		ierr = 0
		IF (n .EQ. 1) RETURN
		e(n) = 0.0
		DO l = 1, n
			j = 0
bucle:		DO
				DO m2 = l, n
					IF ((m2 .EQ. n) .OR. (ABS(e(m2)) .LE. (EPS * (ABS(d(m2)) + ABS(d(m2 + 1)))))) THEN
						m = m2
						EXIT
					END IF
				END DO
				p = d(l)
				IF (m .EQ. l) EXIT bucle

				IF (j .EQ. 30) THEN
					ierr = 1
					RETURN
				END IF
				j = j + 1
				g = (d(l + 1) - p) / (e(l) * 2.0)
				r = SQRT(g * g + 1.0)
				g = d(m) - p + e(l) / (g + SIGN(r, g))
				s = 1.0
				c = 1.0
				p = 0.0
				DO ii = 1, m - l
					i = m - ii
					f = s * e(i)
					b = c * e(i)
					IF (ABS(f) .LT. ABS(g)) THEN
						s = f / g
						r = SQRT(s * s + 1.0)
						e(i + 1) = g * r
						c = 1.0 / r
						s = s * c
					ELSE
						c = g / f
						r = SQRT(c * c + 1.0)
						e(i + 1) = f * r
						s = 1.0 / r
						c = c * s
					END IF
					g = d(i + 1) - p
					r = (d(i) - g) * s + c * 2.0 * b
					p = s * r
					d(i + 1) = g + p
					g = c * r - b
					f = z(i + 1)
					z(i + 1) = s * z(i) + c * f
					z(i    ) = c * z(i) - s * f
				END DO
				d(l) = d(l) - p
				e(l) = g
				e(m) = 0.0
			END DO bucle
		END DO

		DO ii = 2, n
			i = ii - 1
			k = i
			p = d(i)
			DO j = ii, n
				IF (d(j) .GE. p) CYCLE
				k = j
				p = d(j)
			END DO
			IF (k .EQ. i) CYCLE
			d(k) = d(i)
			d(i) = p
			p = z(i)
			z(i) = z(k)
			z(k) = p
		END DO
		RETURN
	END SUBROUTINE gbtql2

	FUNCTION dgamma(z)
		DOUBLE PRECISION dgamma
		DOUBLE PRECISION, INTENT(IN) :: z

		DOUBLE PRECISION, DIMENSION(0:17) :: a
		DATA a /1.0, &
			 0.4227843350984678, &
			 0.4118403304263672, &
			 0.0815769192502609, &
			 0.0742490106800904, &
			-2.669810333484e-4, &
			 0.0111540360240344, &
			-0.0028525821446197, &
			 0.0021036287024598, &
			-9.184843690991e-4, &
 			 4.874227944768e-4, &
			-2.347204018919e-4, &
			 1.115339519666e-4, &
 			-4.78747983834e-5, &
			 1.75102727179e-5, &
			-4.9203750904e-6, &
			 9.199156407e-7, &
			-8.39940496e-8/

		INTEGER k
		DOUBLE PRECISION p, t

		IF (z .LE. 1.0) THEN
			t = z
		ELSE IF (z .LE. 2.0) THEN
			t = z - 1.0
		ELSE
			t = z - 2.0
		END IF

		p = a(17)
		DO k = 1, 17
			p = t * p + a(17 - k)
		END DO

		IF (z .GT. 2.0) THEN
			dgamma = p
		ELSE IF (z .GT. 1.0) THEN
			dgamma = p / z
		ELSE
			dgamma = p / (z * (z + 1.0))
		END IF
		RETURN
	END FUNCTION dgamma

END MODULE gauss
