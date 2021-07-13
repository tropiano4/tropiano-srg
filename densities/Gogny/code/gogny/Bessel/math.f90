!---------------------------------------------------------------------!
!                                                                     !
!          MODULE CONTAINING USEFUL MATHEMATICAL FUNCTIONS            !
!                                                                     !
!---------------------------------------------------------------------!

 MODULE math

	USE input
	USE lgfactor

	IMPLICIT NONE

#if(USE_QUADRUPLE==1)
	REAL(KIND = 16), PARAMETER :: EPS = 1.387778780781445700000000D-17
	REAL(KIND = 16), PARAMETER :: PI = 3.1415926535897932384626430D+00
	REAL(KIND = 16), PARAMETER :: PI_C = 0.12698727186848193957000D+00
	REAL(KIND = 16), PARAMETER :: FOUR_PI = 12.5663706143591729539D+00
#else
	DOUBLE PRECISION, PARAMETER :: EPS = 1.387778780781445700000000D-17
	DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626430D+00
	DOUBLE PRECISION, PARAMETER :: PI_C = 0.12698727186848193957000D+00
	DOUBLE PRECISION, PARAMETER :: FOUR_PI = 12.5663706143591729539D+00
#endif

	DOUBLE PRECISION, PARAMETER :: I_4PI = 0.07957747154594766788D+00
	DOUBLE PRECISION, PARAMETER :: ALPHA = 0.333333333333333333333333D+00
	DOUBLE PRECISION, PARAMETER :: I_SALPHA3 = 0.28056585887484734734D+00 ! This is 1/(alpha + 2)^(3/2)

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

		res = 0.0d0

		! Number of panels is odd.  Use Simpson's 3/8 rule on first three panels, 1/3 rule on rest of them.

		IF ((npanel-2*nhalf).NE.0) THEN

		     res = 3.0d0*step*(functi(nbegin) + 3.0d0*(functi(nbegin+1)+functi(nbegin+2)) + functi(nbegin+3))/8.0d0

		     IF ((npoint-nbegin).EQ.3) RETURN

		     nbegin=nbegin+3

		END IF

		! Apply 1/3 rule - add in first, second, last values

		res = res + step*(functi(nbegin) + 4.0D0*functi(nbegin+1) + functi(npoint))/3.0d0
		nbegin = nbegin+2

		IF (nbegin.EQ.npoint) THEN
		    RETURN
		ELSE

			x = 0.0d0
			nend = npoint - 1

			DO i = nbegin,nend,2
				x = x + functi(i) + 2.0d0*functi(i+1)
			END DO

			res = res + 2.0d0*step*x/3.0d0

			RETURN

		END IF

        END SUBROUTINE simps

	SUBROUTINE Simpson_Kind16(functi,npoint,step,res)

#if(USE_QUADRUPLE==1)
        	REAL(KIND = 16), INTENT(IN) :: step
        	REAL(KIND = 16), DIMENSION(:), INTENT(IN) :: functi
        	REAL(KIND = 16), INTENT(OUT) :: res
#else
        	DOUBLE PRECISION, INTENT(IN) :: step
        	DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: functi
        	DOUBLE PRECISION, INTENT(OUT) :: res
#endif

        	INTEGER :: npanel, npoint, nhalf, nbegin, nend, i

        	!REAL (KIND= 16) :: x
        	DOUBLE PRECISION :: x

        	! Check to see if number of panels is even.  Number of panels is n - 1.

		nbegin = 1
		!nbegin = 0

		npanel = npoint - nbegin
		nhalf  = npanel/2

		res = 0.0D0

		! Number of panels is odd.  Use Simpson's 3/8 rule on first three panels, 1/3 rule on rest of them.

		IF ((npanel-2*nhalf).NE.0) THEN

		     res = 3.0d0*step*(functi(nbegin) + 3.0d0*(functi(nbegin+1)+functi(nbegin+2)) + functi(nbegin+3))/8.0d0

		     IF ((npoint-nbegin).EQ.3) RETURN

		     nbegin=nbegin+3

		END IF

		! Apply 1/3 rule - add in first, second, last values

		res = res + step*(functi(nbegin) + 4.0d0*functi(nbegin+1) + functi(npoint))/3.0d0
		nbegin = nbegin+2

		IF (nbegin.EQ.npoint) THEN
		    RETURN
		ELSE

			x = 0.0d0
			nend = npoint - 1

			DO i = nbegin,nend,2
				x = x + functi(i) + 2.0d0*functi(i+1)
			END DO

			res = res + 2.0d0*step*x/3.0d0

			RETURN

		END IF

        END SUBROUTINE Simpson_Kind16


	! Funciones de apoyo
	FUNCTION c(N, L)
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) c
#else
		DOUBLE PRECISION c
#endif
		INTEGER, INTENT(IN) :: N, L

		c = PI_C * PAR(n) * EXP(0.5d0 * (DDLogFactorials(n) + DDLogSemiFactorials(n + l)))
		RETURN
	END FUNCTION c

	FUNCTION L(a)
		INTEGER L
		INTEGER, INTENT(IN) :: a

		L = (a + 1) / 2
		RETURN
	END FUNCTION L

	FUNCTION J(a)
		INTEGER J
		INTEGER, INTENT(IN) :: a

		J = ((a - L(a)) * 2) + 1
		RETURN
	END FUNCTION J

	FUNCTION LS(a)
		DOUBLE PRECISION LS
		INTEGER, INTENT(IN) :: a

		LS = (0.5D0 * PAR(a)) / ((a / 2) + 1)
		RETURN
	END FUNCTION LS

	! Maximum number for n for a given l value

	FUNCTION DIM(a)
		INTEGER DIM
		INTEGER, INTENT(IN) :: a

			                		DIM = MIN(Nmax,NmaxOfL(L(a)))
		IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	DIM = ((N_0 - L(a)) / 2) + 1

		RETURN
	END FUNCTION DIM

        ! Giving the parity of an integer n (+1 or -1)
	FUNCTION PAR(n)
		INTEGER PAR
		INTEGER, INTENT(IN) :: n

		IF (MOD(n, 2) .EQ. 0) THEN
			PAR = 1
		ELSE
			PAR = -1
		END IF
		RETURN
	END FUNCTION PAR

	SUBROUTINE SWAP(a, b)
		INTEGER, INTENT(INOUT) :: a, b

		INTEGER c

		c = a
		a = b
		b = c
		RETURN
	END SUBROUTINE SWAP

	FUNCTION CUAD2(la, lb, k)
		DOUBLE PRECISION CUAD2
		INTEGER, INTENT(IN) :: la, lb, k

		CUAD2 = DBLE((la * (la + 1)) + (lb * (lb + 1)) - (k * (k + 1)))
		RETURN
	END FUNCTION CUAD2

	FUNCTION MIN_3N(N1, N2, k)
		INTEGER MIN_3N
		INTEGER, INTENT(IN) :: N1, N2, k

		MIN_3N = MAX((ABS(N1 - N2) - k), 0) / 2
		RETURN
	END FUNCTION MIN_3N

	FUNCTION MIN_5N(n1, l1, n2, l2, k)
		INTEGER MIN_5N
		INTEGER, INTENT(IN) :: n1, l1, n2, l2, k

		INTEGER N1_arg, N2_arg

		N1_arg = 2 * n1 + l1
		N2_arg = 2 * n2 + l2
		MIN_5N = MAX((ABS(N1_arg - N2_arg) - k) / 2, 0)
		RETURN
	END FUNCTION MIN_5N

	FUNCTION Char2Int(str)
		INTEGER Char2Int
		CHARACTER(*), INTENT(IN) :: str

		INTEGER idx

		idx = 1
		Char2Int = 0
		DO WHILE (idx .LE. LEN(str))
			SELECT CASE (str(idx:idx))
			CASE ('0')
				Char2Int = Char2Int * 10
			CASE ('1')
				Char2Int = (Char2Int * 10) + 1
			CASE ('2')
				Char2Int = (Char2Int * 10) + 2
			CASE ('3')
				Char2Int = (Char2Int * 10) + 3
			CASE ('4')
				Char2Int = (Char2Int * 10) + 4
			CASE ('5')
				Char2Int = (Char2Int * 10) + 5
			CASE ('6')
				Char2Int = (Char2Int * 10) + 6
			CASE ('7')
				Char2Int = (Char2Int * 10) + 7
			CASE ('8')
				Char2Int = (Char2Int * 10) + 8
			CASE ('9')
				Char2Int = (Char2Int * 10) + 9
			CASE DEFAULT
				RETURN
			END SELECT
			idx = idx + 1
		END DO
		RETURN
	END FUNCTION Char2Int

	SUBROUTINE Int2Char(char_out, num_in)
		CHARACTER(*), INTENT(INOUT) :: char_out
		INTEGER, INTENT(IN) :: num_in

		INTEGER max_len, cur_len, num, i
		CHARACTER digit

		max_len = LEN(char_out)
		cur_len = 0
		char_out = ""
		num = num_in
		DO WHILE ((cur_len .LT. max_len) .AND. (num .GT. 0))
			digit = CHAR(MOD(num, 10) + 48)
			num = num / 10
			cur_len = cur_len + 1
			IF (cur_len .GT. 1) THEN
				DO i = cur_len, 2, -1
					char_out(i:i) = char_out(i - 1:i - 1)
				END DO
			END IF
			char_out(1:1) = digit
		END DO
		char_out(cur_len + 1:cur_len + 1) = CHAR(0)
		RETURN
	END SUBROUTINE Int2Char

END MODULE math
