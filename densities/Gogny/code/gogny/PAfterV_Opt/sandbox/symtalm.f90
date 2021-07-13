! TALMAN

MODULE symtalm

	USE input
	USE lgfactor
	USE math

	IMPLICIT NONE

	TYPE SymCoefficientB5
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16), DIMENSION(:), POINTER :: p
#else
		DOUBLE PRECISION, DIMENSION(:), POINTER :: p
#endif
	END TYPE

	TYPE SymCoefficientB4
		TYPE (SymCoefficientB5), DIMENSION(:), POINTER :: l2
	END TYPE

	TYPE SymCoefficientB3
		TYPE (SymCoefficientB4), DIMENSION(:), POINTER :: n2
	END TYPE

	TYPE SymCoefficientB2
		TYPE (SymCoefficientB3), DIMENSION(:), POINTER :: l1
	END TYPE

	TYPE SymCoefficientB
		TYPE (SymCoefficientB2), DIMENSION(:), POINTER :: n1
	END TYPE

	TYPE (SymCoefficientB) CoeB

	TYPE SymKumar6
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16), DIMENSION(:), POINTER :: L
#else
		DOUBLE PRECISION, DIMENSION(:), POINTER :: L
#endif
	END TYPE

	TYPE SymKumar5
		TYPE (SymKumar6), DIMENSION(:), POINTER :: N
	END TYPE

	TYPE SymKumar4
		TYPE (SymKumar5), DIMENSION(:), POINTER :: l2
	END TYPE

	TYPE SymKumar3
		TYPE (SymKumar4), DIMENSION(:), POINTER :: n2
	END TYPE

	TYPE SymKumar2
		TYPE (SymKumar3), DIMENSION(:), POINTER :: l1
	END TYPE

	TYPE SymKumar
		TYPE (SymKumar2), DIMENSION(:), POINTER :: n1
	END TYPE

	TYPE (SymKumar) Kumar

	DOUBLE PRECISION, PARAMETER :: PI_B = 124.025106721199280701917D0 !  = 4*pi**3

 CONTAINS

	! ATENCION: Este modulo ha de reprogramarse al final del desarrollo
	! En dicho momento, se comprobara cuales son las llamadas mas solicitadas
	! para la evaluacion de "CoefficientB" y "Kumar", preparando asi
	! un registro automatico de las solicitudes con el fin de optimizar
	! el calculo

	!
	!  This subroutine creates and allocates memory for the object of the type SymCoefficientB
	!  which corresponds to the integral B(n1, l1, b1; n2, l2, b2; p, b) whose definition can
	!  be found in the Scary PhD, Appendix D2, page 118, D8.
	!
	SUBROUTINE SymCoefficientB_new

		INTEGER la, na, namax, lc, nc, ncmax, p, pmax

		ALLOCATE(CoeB%n1(0:N_0))
		DO la = 0, N_0
			namax = (N_0 - la) / 2
			ALLOCATE(CoeB%n1(la)%l1(0:namax))
			DO na = 0, namax
				ALLOCATE(CoeB%n1(la)%l1(na)%n2(0:la))
				DO lc = 0, la
					ncmax = (N_0 - lc) / 2
					ALLOCATE(CoeB%n1(la)%l1(na)%n2(lc)%l2(0:ncmax))
					DO nc = 0, ncmax
						pmax = na + nc
						ALLOCATE(CoeB%n1(la)%l1(na)%n2(lc)%l2(nc)%p(0:pmax))
						DO p = 0, pmax
							CoeB%n1(la)%l1(na)%n2(lc)%l2(nc)%p(p) = 0.0D0
						END DO
					END DO
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymCoefficientB_new

	!-------------------------------------------------------------------------------------!
	!										      !
	!  This function is used to evaluate the integral B(n1, l1, b1; n2, l2, b2; p, b).    !
	!  Here, b = b1 = b2 is assumed for some reason and the resulting coefficient is      !
	!  multiplied by the normalization coefficients of the wave-functions c(n1,l1) x      !
	!  c(n2,l2).				   					      !
	!                                                         FORMULA CHECKED AND OK      !
	!										      !
	!  Inputs.  n1, l1, n2, l2, p							      !
	!										      !
	!  Output:  B(n1,l1,n2,l2,p) * c(n1,l1)c(n2,l2)					      !
	!										      !
	!  Ref.: Appendix D2, page 118, D8						      !
	!										      !
	!-------------------------------------------------------------------------------------!

	FUNCTION SymCoefficientB_eval(n1, l1, n2, l2, p)
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) SymCoefficientB_eval
#else
		DOUBLE PRECISION SymCoefficientB_eval
#endif
		INTEGER, INTENT(IN) :: n1, l1, n2, l2, p

		INTEGER k, kmin, kmax
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) sumk
#else
		DOUBLE PRECISION sumk
#endif

		kmin = MAX(p - n2, 0)
		kmax = MIN(p, n1)

		sumk = 0.0D0
		DO k = kmin, kmax
			sumk = sumk + EXP(- DDLogFactorials(k) - DDLogFactorials(n1 - k) &
				          - DDLogSemiFactorials(k + l1) - DDLogFactorials(p - k) &
				          - DDLogFactorials(n2 + k - p) - DDLogSemiFactorials(p + l2 - k))
		END DO

		SymCoefficientB_eval = PI_B * PAR(p + n1 + n2) * sumk * c(n1, l1) * c(n2, l2)

		RETURN
	END FUNCTION SymCoefficientB_eval

	!-------------------------------------------------------------------------------------!
	!										      !
	!  This function is used to:							      !
	!    1) fill in the typed object with the result of the evaluation		      !
	!    2) access the content of the typed object					      !
	!  										      !
	!  Inputs.  na, la, nb, lb, p							      !
	!										      !
	!  Output:  B(n1,l1,n2,l2,p) * c(n1,l1)c(n2,l2)					      !
	!										      !
	!  Ref.: Appendix D2, page 118, D8						      !
	!										      !
	!-------------------------------------------------------------------------------------!

	FUNCTION SymCoefficientB_get(na, la, nb, lb, p)
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) SymCoefficientB_get
#else
		DOUBLE PRECISION SymCoefficientB_get
#endif
		INTEGER, INTENT(IN) :: na, la, nb, lb, p

		INTEGER n1, l1, n2, l2

		IF (lb .GT. la) THEN
			n1 = nb
			l1 = lb
			n2 = na
			l2 = la
		ELSE
			n1 = na
			l1 = la
			n2 = nb
			l2 = lb
		END IF

		! Comprobamos si los parametros se salen del rango de valores registrados
		IF ((((2 * n1) + l1) .GT. N_0) .OR. (((2 * n2) + l2) .GT. N_0)) THEN
			! Si los parametros estan fuera de rango, calculamos el valor
			SymCoefficientB_get = SymCoefficientB_eval(n1, l1, n2, l2, p)
		ELSE
			! Si los parametros estan dentro del rango almacenado,
			! comprobamos si el valor ha sido previamente calculado
			IF (CoeB%n1(l1)%l1(n1)%n2(l2)%l2(n2)%p(p) .EQ. 0.0) THEN
				CoeB%n1(l1)%l1(n1)%n2(l2)%l2(n2)%p(p) = SymCoefficientB_eval(n1, l1, n2, l2, p)
			END IF
			SymCoefficientB_get = CoeB%n1(l1)%l1(n1)%n2(l2)%l2(n2)%p(p)
		END IF
		RETURN
	END FUNCTION SymCoefficientB_get

	!
	!  This subroutine creates an object of the type SymKumar, which definition can be found
	!  at the top of this module: it is a pointer to a pointer.... to a pointer equivalent to
	!  a 6D array with indexes referring to the quantum numbers n1, l1, n2, l2, N, L.
	!  This subroutine allocates the memory for the newly-created object
	!
	SUBROUTINE SymKumar_new

		INTEGER la, na, lc, nc, k
		INTEGER namax, ncmax, kmin, kmax, N1min, N1max, N1
		INTEGER Na2, Nc2
		INTEGER i2, i6, i7

		ALLOCATE(Kumar%n1(0:N_0))
		DO la = 0, N_0
			namax = (N_0 - la) / 2
			ALLOCATE(Kumar%n1(la)%l1(0:namax))
			DO na = 0, namax
				ALLOCATE(Kumar%n1(la)%l1(na)%n2(0:la))
				Na2 = (2 * na) + la
				DO lc = 0, la
					ncmax = (N_0 - lc) / 2
					ALLOCATE(Kumar%n1(la)%l1(na)%n2(lc)%l2(0:ncmax))
					kmin = ABS(la - lc)
					kmax = la + lc
					i2 = ((kmax - kmin) / 2) ! + 1
					DO nc = 0, ncmax
						Nc2 = 2 * nc + lc
						ALLOCATE(Kumar%n1(la)%l1(na)%n2(lc)%l2(nc)%N(0:i2))
						DO k = kmin, kmax, +2
							N1min = MIN_3N(Na2, Nc2, k)
							N1max = na + nc + ((la + lc - k) / 2)
							i6 = (k - kmin) / 2
							i7 = N1max - N1min ! + 1
							ALLOCATE(Kumar%n1(la)%l1(na)%n2(lc)%l2(nc)%N(i6)%L(0:i7))
							DO N1 = N1min, N1max
								Kumar%n1(la)%l1(na)%n2(lc)%l2(nc)%N(i6)%L(N1 - N1min) = 0.0D0
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymKumar_new

	! ojo: es c(N,L) * Kumar(n1,l1,n2,l2,N,L)

	!
	!  This function evaluates the so-called Kumar coefficient times the coefficient C(N,L).
	!  Definitions can be found in the Horrible PhD, Appendix F, Page 142.
	!                                                         !!! FORMULA not CHECKED !!!
	!
	!  Inputs.  n1, l1, n2, l2, N, L
	!
	!  Output:  c(N,L) * TK(n1,l1,n2,l2,N,L)
	!
	!  Ref.: Appendix F, page 142
	!
	FUNCTION SymKumar_eval(n1, l1, n2, l2, N, L)
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) SymKumar_eval
#else
		DOUBLE PRECISION SymKumar_eval
#endif
		INTEGER, INTENT(IN) :: n1, l1, n2, l2, N, L

		INTEGER i1, i2, p, pmin, pmax
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) sump,fac
#else
		DOUBLE PRECISION sump,fac
#endif
		fac = 1.0D0 / ( 2.0D0 * (4.0D0*ATAN(1.0D0))**(0.75D0) )

		i1 = (l1 + l2 + L) / 2
		i2 = (l1 + l2 - L) / 2

		pmin = MAX(N - i2, 0)
		pmax = n1 + n2 ! Upper limit of the summation

		sump = SymCoefficientB_get(n1, l1, n2, l2, pmax)

		DO p = pmax - 1, pmin, -1
			sump = SymCoefficientB_get(n1, l1, n2, l2, p) &
				+ ((sump * (DBLE(p + i1) + 1.5D0) * DBLE(p + i2 + 1)) / DBLE(p + i2 + 1 - N))
		END DO

		SymKumar_eval = fac * EXP(DDLogSemiFactorials(i1 + pmin) + DDLogFactorials(i2 + pmin)) * sump

		IF ((i2 - N) .GT. 0) THEN
			SymKumar_eval = SymKumar_eval / EXP(DDLogFactorials(i2 - N))
		END IF

		RETURN
	END FUNCTION SymKumar_eval

	! Atencion: es c(N,L) * Kumar(n1, l1, n2, l2, N, L)

	!
	!  This function is used to:
	!    1) fill in the typed object Kumar with the result of the evaluation
	!    2) access the content of the typed object
	!
	!  Inputs.  na, la, nb, lb, N, L
	!
	!  Output:  c(N,L) * TK(na,la,nb,lb,N,L)
	!
	!  Ref.: Appendix F, page 142
	!
	FUNCTION SymKumar_get(na, la, nb, lb, N, L)
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) SymKumar_get
#else
		DOUBLE PRECISION SymKumar_get
#endif
		INTEGER, INTENT(IN) :: na, la, nb, lb, N, L

		INTEGER n1, l1, n2, l2
		INTEGER M1, M2
		INTEGER i2, i8, i9, Lmi, Lmx, Nmi, Nmx

		M1 = (na * 2) + la
		M2 = (nb * 2) + lb

		! Formula antigua
		i2 = (la + lb - L) / 2
		Lmi = ABS(la - lb)
		Lmx = la + lb
		i8 = (L - Lmi) / 2

		Nmx = na + nb + i2
		Nmi = MIN_3N(M1, M2, L)
		i9 = N - Nmi

		IF ((M1 .GT. N_0) .OR. (M2 .GT. N_0)) THEN
			SymKumar_get = SymKumar_eval(na, la, nb, lb, N, L)
			RETURN
		END IF

		IF ((L .LT. Lmi) .OR. (L .GT. Lmx) .OR. &
		    (N .LT. Nmi) .OR. (N .GT. Nmx)) THEN
			SymKumar_get = 0.0D0
			RETURN
		END IF

		IF (lb .GT. la) THEN
			n1 = nb
			l1 = lb
			n2 = na
			l2 = la
		ELSE
			n1 = na
			l1 = la
			n2 = nb
			l2 = lb
		END IF

		IF (Kumar%n1(l1)%l1(n1)%n2(l2)%l2(n2)%N(i8)%L(i9) .EQ. 0.0) THEN
			Kumar%n1(l1)%l1(n1)%n2(l2)%l2(n2)%N(i8)%L(i9) = SymKumar_eval(n1, l1, n2, l2, N, L)
		END IF

		SymKumar_get = Kumar%n1(l1)%l1(n1)%n2(l2)%l2(n2)%N(i8)%L(i9)

		RETURN
	END FUNCTION SymKumar_get

END MODULE symtalm
