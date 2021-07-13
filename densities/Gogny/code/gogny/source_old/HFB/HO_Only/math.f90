!---------------------------------------------------------------------!
!                                                                     !
!          MODULE CONTAINING USEFUL MATHEMATICAL FUNCTIONS            !                
!                                                                     !
!---------------------------------------------------------------------!

 MODULE math

	USE input
	USE lgfactor

	IMPLICIT NONE

        DOUBLE PRECISION, PARAMETER :: EPS = 1.3877787807814457e-17

        DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793238462643
        DOUBLE PRECISION, PARAMETER :: PI_C = 0.12698727186848193957
	DOUBLE PRECISION, PARAMETER :: I_4PI = 0.0795774715459
        DOUBLE PRECISION, PARAMETER :: FOUR_PI = 12.5663706143591729539

	DOUBLE PRECISION, PARAMETER :: ALPHA = 0.333333333333
	DOUBLE PRECISION, PARAMETER :: I_SALPHA3 = 0.280565858875 ! This is 1/(alpha + 2)^(3/2)

 CONTAINS
	
	! Funciones de apoyo
	FUNCTION c(N, L)
                DOUBLE PRECISION c
		INTEGER, INTENT(IN) :: N, L

		c = PI_C * PAR(n) * EXP(0.5 * (DDLogFactorials(n) + DDLogSemiFactorials(n + l)))
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

		LS = (0.5 * PAR(a)) / ((a / 2) + 1)
		RETURN
	END FUNCTION LS

	! Maximum number for n for a given l value

	FUNCTION DIM(a)
		INTEGER DIM
		INTEGER, INTENT(IN) :: a

		DIM = ((N_0 - L(a)) / 2) + 1
		
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

		CUAD2 = (la * (la + 1.0)) + (lb * (lb + 1.0)) - (k * (k + 1.0))
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
