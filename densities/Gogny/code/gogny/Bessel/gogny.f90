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
!-----------------------------------------------------------------------!
!									!
!									!
!     RADIAL INTEGRALS: BRINK-BOKER FORCE				!
!									!
!  This module computes the radial integral of the Brink-Boeker force 	!
!  in the 2 distinct cases of an "analytical" harmonic oscillator basis !
!  and a general Woods-Saxon (or other) basis.				!
!									!
!-----------------------------------------------------------------------!

 MODULE bessik

	IMPLICIT NONE

 CONTAINS

	!---------------------------------------------------------------------
	!     modified bessel function of fractional order
	!     xnu = order of the Bessel function nu
	!     ri   I_nu(x)  first kind
	!     rk   K_nu(x)  second kind
	!     rip  derivative of I(x)
	!     rkp  derivative of K(x)
	!     from Numerical Recipes.
	!
	!---------------------------------------------------------------------

	SUBROUTINE bessel(x, xnu, ri, rk, rip, rkp)
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16), INTENT(IN) :: x, xnu
		REAL(KIND = 16), INTENT(OUT) :: ri, rip, rk, rkp
#else
		DOUBLE PRECISION, INTENT(IN) :: x, xnu
		DOUBLE PRECISION, INTENT(OUT) :: ri, rip, rk, rkp
#endif

		INTEGER :: MAXIT, i, l, nl

#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: EPS,FPMIN,PI,XMIN
		REAL(KIND = 16) :: a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff
		REAL(KIND = 16) :: gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1
		REAL(KIND = 16) :: ripl,ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
#else
		DOUBLE PRECISION :: EPS,FPMIN,PI,XMIN
		DOUBLE PRECISION :: a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff
		DOUBLE PRECISION :: gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1
		DOUBLE PRECISION :: ripl,ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
#endif

		PARAMETER (EPS=1.D-16, FPMIN=1.D-30 ,MAXIT=10000, XMIN=2.0D0)

		PI = 4.0D0*ATAN(1.0D0)

		IF (x.LE.0.0D0 .OR. xnu.LT.0.0D0) THEN
		   WRITE(*,'("x = ",E24.16," nu = ",F10.5)') x,xnu
		   stop "Bad arguments in bessik"
		END IF

		nl = int(xnu + 0.5D0)
		xmu = xnu - nl
		xmu2 = xmu*xmu
		xi = 1.0D0/x
		xi2 = 2.0D0*xi
		h = xnu*xi

		IF (h .LT. FPMIN) h = FPMIN

		b = xi2*xnu
		d = 0.0D0
		c = h

		DO i = 1,MAXIT

			b = b + xi2
			d = 1.0D0/(b + d)
			c = b + 1.0D0/c
			del = c*d
			h = del*h

			IF (ABS(del - 1.0D0) .LT. EPS) EXIT

		END DO

		IF (i .EQ. MAXIT) STOP 'x too large in bessik; try asymptotic expansion'

		ril  = FPMIN
		ripl = h*ril
		ril1 = ril
		rip1 = ripl
		fact = xnu*xi

		DO l = nl, 1, -1
			ritemp = fact*ril + ripl
			fact = fact - xi
			ripl = fact*ritemp + ril
			ril = ritemp
		END DO

		f = ripl/ril

		IF (x .LT. XMIN) THEN

			x2 = 0.5D0*x
			pimu = PI*xmu

			IF (ABS(pimu) .LT. EPS) THEN
				fact = 1.0D0
			ELSE
				fact = pimu/SIN(pimu)
			END IF

			d = -LOG(x2)
			e = xmu*d

			IF (ABS(e) .LT. EPS) THEN
				fact2 = 1.0D0
			ELSE
				fact2 = sinh(e)/e
			END IF

			CALL beschb(xmu, gam1, gam2, gampl, gammi)

			ff = fact*(gam1*COSH(e) + gam2*fact2*d)
			sum = ff
			e = EXP(e)
			p = 0.5D0*e/gampl
			q = 0.5D0/(e*gammi)
			c = 1.0D0
			d = x2*x2
			sum1 = p

			DO i = 1,MAXIT
	!!!!!write(*,'("x < XMIN - i = ",i5)') i
				ff = (i*ff + p + q)/(i*i - xmu2)
				c = c*d/i
				p = p/(i - xmu)
				q = q/(i + xmu)
				del = c*ff
				sum = sum + del
				del1 = c*(p - i*ff)
				sum1 = sum1 + del1
				IF (abs(del) .LT. ABS(sum)*EPS) EXIT
			END DO

			IF (i .EQ. MAXIT) STOP 'bessk series failed to converge'

			rkmu = sum
			rk1 = sum1*xi2

		ELSE

			b = 2.0D0*(1.0D0 + x)
			d = 1.0D0/b
			delh = d
			h = delh
			q1 = 0.0D0
			q2 = 1.0D0
			a1 = 0.25D0 - xmu2
			c = a1
			q = c
			a = -a1
			s = 1.0D0 + q*delh

			DO i = 2,MAXIT
	!write(*,'("x > XMIN - i = ",i5)') i
				a = a - 2*(i-1)
				c = -a*c/i
				qnew = (q1 - b*q2)/a
				q1 = q2
				q2 = qnew
				q = q + c*qnew
				b = b + 2.0D0
				d = 1.0D0/(b + a*d)
				delh = (b*d - 1.0D0)*delh
				h = h + delh
				dels = q*delh
				s = s + dels
				IF (ABS(dels/s) .LT. EPS) EXIT
			END DO

			IF (i .EQ. MAXIT) STOP 'bessik: failure to converge in cf2'

			h = a1*h
			rkmu = SQRT(PI/(2.0D0*x))*exp(-x)/s
			rk1 = rkmu*(xmu + x + 0.5D0 - h)*xi

		END IF

		rkmup = xmu*xi*rkmu - rk1
		rimu = xi/(f*rkmu - rkmup)
		ri = (rimu*ril1)/ril
		rip = (rimu*rip1)/ril

		DO i=1,nl
			rktemp = (xmu + i)*xi2*rk1 + rkmu
			rkmu = rk1
			rk1 = rktemp
		END DO

		rk = rkmu
		rkp = xnu*xi*rkmu - rk1

		RETURN
	END SUBROUTINE bessel


	SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)

		INTEGER :: NUSE1,NUSE2
		PARAMETER (NUSE1=7,NUSE2=8)

#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: gam1, gam2, gammi, gampl, x, xx, one
		REAL(KIND = 16), POINTER, DIMENSION(:):: c1, c2
#else
		DOUBLE PRECISION :: gam1, gam2, gammi, gampl, x, xx, one
		DOUBLE PRECISION, POINTER, DIMENSION(:):: c1, c2
#endif

		ALLOCATE(c1(7))

		one = 1.0D0

		c1(1) = -1.142022680371168D0
		c1(2) =  6.5165112670737D-3
		c1(3) =  3.087090173086D-4
		c1(4) = -3.4706269649D-6
		c1(5) =  6.9437664D-9
		c1(6) =  3.67795D-11
		c1(7) = -1.356D-13

		ALLOCATE(c2(8))

		c2(1) =  1.843740587300905D0
		c2(2) = -7.68528408447867D-2
		c2(3) =  1.2719271366546D-3
		c2(4) = -4.9717367042D-6
		c2(5) = -3.31261198D-8
		c2(6) =  2.423096D-10
		c2(7) = -1.702D-13
		c2(8) = -1.49D-15

		xx = 8.0D0*x*x - 1.0D0

		gam1 = chebev(-one, one, c1, NUSE1, xx)
		gam2 = chebev(-one, one, c2, NUSE2, xx)

		DEALLOCATE(c1)
		DEALLOCATE(c2)

		gampl = gam2 - x*gam1
		gammi = gam2 + x*gam1

		RETURN
	END SUBROUTINE beschb

	FUNCTION chebev(a, b, c, m, x)
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: chebev
#else
		DOUBLE PRECISION :: chebev
#endif

		INTEGER, INTENT(IN) :: m

#if(USE_QUADRUPLE==1)
		REAL(KIND = 16), INTENT(IN) ::a, b, x
#else
		DOUBLE PRECISION, INTENT(IN) ::a, b, x
#endif

#if(USE_QUADRUPLE==1)
		REAL(KIND = 16), POINTER, DIMENSION(:) :: c
#else
		DOUBLE PRECISION, POINTER, DIMENSION(:) :: c
#endif

		INTEGER :: j

#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: d, dd, sv, y, y2
#else
		DOUBLE PRECISION :: d, dd, sv, y, y2
#endif

		IF ((x - a)*(x - b) .GT. 0.0D0) THEN
			write(*,'("x = ",E24.16," a = ",E24.16," b = ",E24.16)') x,a,b
			STOP 'x not in range in chebev'
		END IF

		d = 0.0D0
		dd = 0.0D0
		y = (2.0D0*x - a - b)/(b - a)
		y2 = 2.0D0*y

		DO j = m, 2, -1
			sv = d
			d = y2*d - dd + c(j)
			dd = sv
		END DO

		chebev = y*d - dd + 0.5D0*c(1)

		RETURN
	END FUNCTION chebev

END MODULE bessik
	! This is only used to define a type "OneDimSolve" that will 
	! be used in the determination of the particle number. This file 
	! is not used as module, but directly in an INCLUDE statement.

	TYPE OneDimSolve
		TYPE (R1R1Function), POINTER :: func
		DOUBLE PRECISION AccX ! Accuracy in X direction
		DOUBLE PRECISION AccY ! Accuracy in Y direction
		INTEGER lim ! Maximum number of iterations
		DOUBLE PRECISION, DIMENSION(0:2) :: x, y ! Trial values of X and Y
		INTEGER L, C, U, Last ! Locations of trial values
		INTEGER vpol, hpol ! Polarities
		DOUBLE PRECISION YY ! Target value
		LOGICAL Finish ! .TRUE. if LookAt finds conv.
		LOGICAL Captured ! .TRUE. when target surrounded
		TYPE (DiagonalizationMethod), POINTER :: diagonal
	END TYPE
	! This is the set of subroutines that determine the particle number
	! at each iteration. Look at subroutine OneDimSolve_LookAt(solve, V)
	! where there is a call to DiagonalizationMethod_operator (which
	! calculates the particle number).

	SUBROUTINE OneDimSolve_new(solve, func, diagonal)
		TYPE (OneDimSolve), INTENT(INOUT) :: solve
		TYPE (R1R1Function), TARGET, INTENT(IN) :: func
		TYPE (DiagonalizationMethod), TARGET, INTENT(IN) :: diagonal

		solve%func => func
		solve%AccX = 1.D-12; solve%AccY = 1.D-12
		solve%lim = 100
		solve%x(0) = 0.0D0; solve%x(1) = 0.0D0; solve%x(2) = 0.0D0
		solve%y(0) = 0.0D0; solve%y(1) = 0.0D0; solve%y(2) = 0.0D0
		solve%L = 1; solve%C = 1; solve%U = 1; solve%Last = 1
		solve%vpol = 1; solve%hpol = 1
		solve%YY = 0.0D0
		solve%Finish = .FALSE.
		solve%Captured = .FALSE.
		solve%diagonal => diagonal

		RETURN
	END SUBROUTINE OneDimSolve_new

	! This subroutine solves X = Y (Y target value) using some kind of mysterious method
	! Dev_in is the initial step for the change in X

	FUNCTION OneDimSolve_solve(solve, Y, X, Dev_in)
		DOUBLE PRECISION OneDimSolve_solve
		TYPE (OneDimSolve), INTENT(INOUT) :: solve
		DOUBLE PRECISION, INTENT(IN) :: Y, X, Dev_in

		INTEGER, PARAMETER :: Lim = 100

		INTEGER State ! 0 = start, 1 = captured1, 2 = captured2, 3 = binary, 4 = finish
		INTEGER i
		DOUBLE PRECISION Dev

		Dev = Dev_in

		solve%lim = Lim
		solve%Captured = .FALSE.

		IF (Dev .EQ. 0.0D0) STOP "Dev es nulo"

		solve%L = 0
		solve%C = 1
		solve%U = 2
		solve%vpol = 1
		solve%hpol = 1
		solve%y(solve%C) = 0.0D0
		solve%y(solve%U) = 0.0D0

		IF (Dev < 0.0D0) THEN
			solve%hpol = - 1
			Dev = - Dev
		END IF

		solve%YY = Y ! Target value
		solve%x(solve%L) = X ! Initial trial value

		IF (.NOT. R1R1Function_isValid(solve%func, X)) STOP "Valor inicial no valido"

		State = 0
bucle:	DO
			SELECT CASE (State)
			CASE (0)

				! Left point: Get value of the function
				CALL OneDimSolve_LookAt(Solve, solve%L)

				IF (solve%Finish) THEN
					State = 4
					CYCLE
				END IF

				! Force a negative sign if we encounter a positive value
				IF (solve%y(solve%L) .GT. 0.0D0) THEN
					CALL OneDimSolve_VFlip(Solve) ! So Y(L) < 0
				END IF

				! Right point
				solve%x(solve%U) = X + (Dev * solve%hpol)

				! Check boundaries for right point
				IF ((.NOT. solve%func%infinite_x_max) .AND. (solve%x(solve%U) .GT. solve%func%x_max)) THEN
					solve%x(solve%U) = (solve%func%x_max + X) / 2.0D0
				END IF

				IF ((.NOT. solve%func%infinite_x_min) .AND. (solve%x(solve%U) .LT. solve%func%x_min)) THEN
					solve%x(solve%U) = (solve%func%x_min + X) / 2.0D0
				END IF

				! Right point: Get value of the function
				CALL OneDimSolve_LookAt(Solve, solve%U)

				IF (solve%Finish) THEN
					State = 4
					CYCLE
				END IF

				! If we have the right point positive and the left point negative, we know where the
				! crossing is (Case State = 1 below)
				IF (solve%y(solve%U) .GT. 0.0D0) THEN
					State = 1
					solve%Captured = .TRUE.
					CYCLE
				END IF

				! If both the left and right point have the same value, we have a "plane" function
				IF (solve%y(solve%U) .EQ. solve%y(solve%L)) STOP "Constant function..."

				! If the value of the function on the right is lower than the value on the left, we must change the
				! way we increment X (decrease instead of increase or vice-versa).
				IF (solve%y(solve%U) .LT. solve%y(solve%L)) THEN
					CALL OneDimSolve_HFlip(Solve) ! Change direction
				END IF

				CALL OneDimSolve_State(Solve, solve%L, solve%U, solve%C)

				DO i = 0, 19
					! Have L C then crossing point, Y[L] < Y[C] <0
					solve%x(solve%U) = solve%x(solve%C) + DBLE(Dev * solve%hpol)

					IF ((.NOT. solve%func%infinite_x_max) .AND. (solve%x(solve%U) .GT. solve%func%x_max)) THEN
						solve%x(solve%U) = (solve%func%x_max + solve%x(solve%C)) / 2.0D0
					END IF
					IF ((.NOT. solve%func%infinite_x_min) .AND. (solve%x(solve%U) .LT. solve%func%x_min)) THEN
						solve%x(solve%U) = (solve%func%x_min + solve%x(solve%C)) / 2.0D0
					END IF

					CALL OneDimSolve_LookAt(Solve, solve%U)

					IF (solve%Finish) THEN
						State = 4
						CYCLE bucle
					END IF

					IF (solve%y(solve%U) .GT. 0) THEN
						State = 2
						solve%Captured = .TRUE.
						CYCLE bucle
					END IF

					Dev = Dev * 2.0D0

					CALL OneDimSolve_State(Solve, solve%C, solve%U, solve%L)

				END DO
				IF (State .NE. 0) EXIT bucle
				STOP "Could not find a crossing point"

			CASE (1)
				! We have 2 points L and U with crossing between them
				CALL OneDimSolve_Linear(Solve, solve%L, solve%C, solve%U) ! Linear interpolation

				! Result to C
				CALL OneDimSolve_LookAt(Solve, solve%C)

				IF (solve%Finish) THEN
					State = 4
					CYCLE
				END IF

				IF (solve%y(solve%C) .GT. 0.0D0) THEN
					CALL OneDimSolve_Flip(Solve) ! Want y[C] < 0
				END IF

				IF (solve%y(solve%C) .LT. (0.5D0 * solve%y(solve%L))) THEN
					CALL OneDimSolve_State(Solve, solve%C, solve%L, solve%U)
					State = 3 ! Binary
					CYCLE
				END IF

				State = 2
				CYCLE

			CASE (2)
				! We have L,C before crossing, U after crossing
				CALL OneDimSolve_Quadratic(Solve, solve%L, solve%C, solve%U) ! Quad interpolation
				! Result to L
				CALL OneDimSolve_State(Solve, solve%C, solve%L, solve%U)
				IF ((((solve%x(solve%C) - solve%x(solve%L)) * solve%hpol) .LE. 0.0D0) .OR. &
				    (((solve%x(solve%C) - solve%x(solve%U)) * solve%hpol) .GE. 0.0D0)) THEN
					State = 1
					CYCLE
				END IF

				CALL OneDimSolve_LookAt(Solve, solve%C)
				IF (solve%Finish) THEN
					State = 4
					CYCLE
				END IF

				IF (solve%y(solve%C) .GT. 0.0D0) THEN
					CALL OneDimSolve_Flip(Solve) ! Want y[C] < 0
				END IF

				IF (solve%y(solve%C) .GT. (0.5D0 * solve%y(solve%L))) THEN
					State = 2
					CYCLE
				ELSE
					CALL OneDimSolve_State(Solve, solve%C, solve%L, solve%U)
					State = 1
					CYCLE
				END IF

				State = 3
				CYCLE

			CASE (3)
				! We have L, U around crossing - do binary search
				DO i = 3, 1, -1

					solve%x(solve%C) = 0.5D0 * (solve%x(solve%L) + solve%x(solve%U))

					CALL OneDimSolve_LookAt(Solve, solve%C)

					IF (solve%Finish) THEN
						State = 4
						CYCLE bucle
					END IF

					IF (solve%y(solve%C) .GT. 0.0D0) THEN
						CALL OneDimSolve_State(Solve, solve%L, solve%U, solve%C)
					ELSE
						CALL OneDimSolve_State(Solve, solve%C, solve%L, solve%U)
					END IF
				END DO

				IF (State .NE. 3) CYCLE

				State = 1

				CYCLE

			CASE (4)
				OneDimSolve_solve = solve%x(solve%Last)
				RETURN
			END SELECT
		END DO bucle
		RETURN
	END FUNCTION OneDimSolve_solve

	SUBROUTINE OneDimSolve_LookAt(solve, V)
		TYPE (OneDimSolve), INTENT(INOUT) :: solve
		INTEGER, INTENT(IN) :: V

		DOUBLE PRECISION yy

		solve%lim = solve%lim - 1

		IF (solve%lim .EQ. 0) STOP "No convergence"

		solve%Last = V

!TODO yy = function(solve%x(V)) - solve%YY

		! Set the function to its initial value "properly" (chacking the bounds)
		CALL R1R1Function_set(solve%func, solve%x(V))

		! yy gets the difference between the current number of particles DiagonalizationMethod_operator(solve%diagonal)
		! and the target value
		yy = DiagonalizationMethod_operator(solve%diagonal) - solve%YY

		! If precision criteria are met, we prepare the exit, if not, we will continue by storing the value of the last
		! try
		IF ((ABS(yy) .LE. solve%AccY) .OR. &
		    (solve%Captured .AND. (ABS(solve%x(solve%L) - solve%x(solve%U)) .LE. solve%AccX))) THEN
			solve%Finish = .TRUE.
		ELSE
			solve%Finish = .FALSE.
		END IF

		solve%y(V) = solve%vpol * yy

		RETURN
	END SUBROUTINE OneDimSolve_LookAt

	SUBROUTINE OneDimSolve_HFlip(solve)
		TYPE (OneDimSolve), INTENT(INOUT) :: solve

		solve%hpol = - solve%hpol
		CALL OneDimSolve_State(Solve, solve%U, solve%C, solve%L)
		RETURN
	END SUBROUTINE OneDimSolve_HFlip

	SUBROUTINE OneDimSolve_VFlip(solve)
		TYPE (OneDimSolve), INTENT(INOUT) :: solve

		solve%vpol = - solve%vpol
		solve%y(0) = - solve%y(0) !TODO solve%y = (-1) * solve%y
		solve%y(1) = - solve%y(1)
		solve%y(2) = - solve%y(2)

		RETURN
	END SUBROUTINE OneDimSolve_VFlip

	SUBROUTINE OneDimSolve_Flip(solve)
		TYPE (OneDimSolve), INTENT(INOUT) :: solve

		solve%hpol = - solve%hpol
		solve%vpol = - solve%vpol

		CALL OneDimSolve_State(Solve, solve%U, solve%C, solve%L)

		solve%y(0) = - solve%y(0) !TODO solve%y = (-1) * solve%y
		solve%y(1) = - solve%y(1)
		solve%y(2) = - solve%y(2)

		RETURN
	END SUBROUTINE OneDimSolve_Flip

	SUBROUTINE OneDimSolve_State(solve, L, C, U)
		TYPE (OneDimSolve), INTENT(INOUT) :: solve
		INTEGER, INTENT(IN) :: L, C, U

		! ATENCION: Hacemos una doble asignacion, por ejemplo: L => I => solve%L
		! debido a que a pesar de que los parametros de entrada no son
		! "modificables", la llamada a esta funcion se suele realizar de la
		! forma: OneDimSolve_State(solve, solve%L, solve%C, solve%U)
		! que puede producir una incosistencia por la secuenciación de las
		! operaciones. Por ejemplo:
		! # OneDimSolve_State(solve, solve%C, solve%L, solve%U)
		! # solve%L = solve%C
		! # solve%C = solve%L ERROR!
		INTEGER I, J, K

		I = L
		J = C
		K = U

		solve%L = I
		solve%C = J
		solve%U = K

		RETURN
	END SUBROUTINE OneDimSolve_State

	SUBROUTINE OneDimSolve_Linear(solve, L, C, U)
		TYPE (OneDimSolve), INTENT(INOUT) :: solve
		INTEGER, INTENT(IN) :: L, C, U

		INTEGER I, J, K

		I = L
		J = C
		K = U

		IF (ABS(solve%y(K) - solve%y(I)) .LT. TINY(DBLE(1.0))) THEN
			WRITE(*,'("Module brent_d2.f90 - Subroutine OneDimSolve_Linear")')
			WRITE(*,'("Linear Interpolation attempted with denominator equal to 0")')
			WRITE(*,'("Denominator : ",E24.16)') solve%y(K) - solve%y(I)
			STOP "Error in OneDimSolve_Linear"
		END IF

		solve%x(J) = ((solve%x(I) * solve%y(K)) - (solve%x(K) * solve%y(I))) / (solve%y(K) - solve%y(I))

		RETURN
	END SUBROUTINE OneDimSolve_Linear

	SUBROUTINE OneDimSolve_Quadratic(solve, L, C, U)
		TYPE (OneDimSolve), INTENT(INOUT) :: solve
		INTEGER, INTENT(IN) :: L, C, U
		! Result to overwrite I
		DOUBLE PRECISION YJK, YIK, YIJ, XKI, XKJ

		INTEGER I, J, K

		I = L
		J = C
		K = U

		YJK = solve%y(J) - solve%y(K)
		YIK = solve%y(I) - solve%y(K)
		YIJ = solve%y(I) - solve%y(J)
		XKI = solve%x(K) - solve%x(I)

		XKJ = ((solve%x(K) * solve%y(J)) - (solve%x(J) * solve%y(K))) / YJK

		IF (ABS(XKI) .LT. TINY(DBLE(1.0)) .OR. ABS(YIK) .LT. TINY(DBLE(1.0)) .OR. ABS(YIJ) .LT. TINY(DBLE(1.0))) THEN
			WRITE(*,'("Module brent_d2.f90 - Subroutine OneDimSolve_Quadratic")')
			WRITE(*,'("Quadratic Interpolation attempted with denominator equal to 0")')
			WRITE(*,'("Denominator XKI : ",E24.16)') XKI
			WRITE(*,'("Denominator YIK : ",E24.16)') YIK
			WRITE(*,'("Denominator YIJ : ",E24.16)') YIJ
			STOP "Error in OneDimSolve_Quadratic"
		END IF

		IF ((((YJK / YIK) ** 2) .GT. ((solve%x(K) - solve%x(J)) / XKI)) .OR. &
		    (((YIJ / YIK) ** 2) .GT. ((solve%x(J) - solve%x(I)) / XKI))) THEN
			solve%x(I) = XKJ
		ELSE
			XKI = (solve%x(K) * solve%y(I) - solve%x(I) * solve%y(K)) / YIK
			solve%x(I) = (XKJ * solve%y(I) - XKI * solve%y(J)) / YIJ
		END IF

		RETURN
	END SUBROUTINE OneDimSolve_Quadratic
 MODULE diagmeth

	USE input
	USE symd3t
	USE r1r1
	USE nucleus
	USE symden
	USE symgden
	USE symgdhf
	USE symke2b
	USE symvbb
	USE symvc
	USE symvls
	USE symgdd
	USE selfc
	USE eigenval
	USE jacobi
	USE indexx
        USE global
        USE laguerre

	IMPLICIT NONE

	INCLUDE "brent_d1.f90"

	INTEGER, PARAMETER :: MAX_ITER = 5

	DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DensityRadial, DensityPairing, DerivDens

	TYPE MatrixType
		DOUBLE PRECISION, POINTER, DIMENSION(:, :) :: quantum
	END TYPE

	TYPE DiagonalMatrix
		DOUBLE PRECISION, POINTER, DIMENSION(:) :: value
	END TYPE

	TYPE DiagonalizationMethod
		TYPE (SelfConsistencyMethod) consistency
		TYPE (R1R1Function) func
		TYPE (MatrixType), POINTER, DIMENSION(:, :) :: UV ! U and V vectors
		TYPE (DiagonalMatrix), POINTER, DIMENSION(:, :) :: QuasiParticleEnergies ! Self-explanatory...
		TYPE (SymGenDensity) S ! Superhamiltonian
		TYPE (SymGenDensity) iterated ! Density at each iteration
		INTEGER ta
	END TYPE

 CONTAINS

	INCLUDE "brent_d2.f90"

	! This subroutine initializes the self-consistent calculation by creating all necessary objects

	SUBROUTINE DiagonalizationMethod_new(diagonal, density)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		TYPE (SymDensity), INTENT(INOUT) :: density

		INTEGER :: ta, a, max_a, d, dd

		! Inicializamos las variables y subtipos
		CALL SelfConsistencyMethod_new(diagonal%consistency, density)

		CALL R1R1Function_new(diagonal%func)

		CALL SymGenDensity_new_Nucleus(diagonal%S, Nucleus_get_N(density%nucleus), Nucleus_get_Z(density%nucleus))
		CALL SymGenDensity_new_SymDensity(diagonal%iterated, density)

		max_a = 2*Lmax

		! Reservamos memoria para las matrices del tipo
		ALLOCATE(diagonal%QuasiParticleEnergies(0:1, 0:max_a))
		IF (.NOT. ASSOCIATED(diagonal%QuasiParticleEnergies)) STOP "Unable to allocate memory in DiagonalizationMethod_new (1)"

		ALLOCATE(diagonal%UV(0:1, 0:max_a))
		IF (.NOT. ASSOCIATED(diagonal%UV)) STOP "Unable to allocate memory in DiagonalizationMethod_new (2)"

		DO ta = 0, 1
			DO a = 0, max_a

				d = DIM(a)
				dd = d + d

				ALLOCATE(diagonal%QuasiParticleEnergies(ta, a)%value(dd))
				IF (.NOT. ASSOCIATED(diagonal%QuasiParticleEnergies(ta, a)%value)) &
				     STOP "Unable to allocate memory in DiagonalizationMethod_new (3)"

				ALLOCATE(diagonal%UV(ta, a)%quantum(dd, dd))
				IF (.NOT. ASSOCIATED(diagonal%UV(ta, a)%quantum)) &
				    STOP "Unable to allocate memory in DiagonalizationMethod_new (4)"

			END DO
		END DO

		RETURN
	END SUBROUTINE DiagonalizationMethod_new

	!-----------------------------------------------------------------------------------------------!
	!   This subroutine performs the self-consistent HFB calculation. Firstly, it calculates the	!
	!   matrix elements. Eventually, if the latter were already calculated before, it reads the 	!
	!   files that contain them. Then, it enters the interation loops. At each step, it calculates	!
	!   the mean-field h and Delta, solves the particle number equation and saves the density thus 	!
	!   obtained and cycles. The convergence is achieved when the difference between the norms of	!
	!   the density between 2 iterations is lower than "tolerance". At the convergence, a summary 	!
	!   of the calculation is printed, then the single-particle energies in the canonical basis and	!
	!   the densities are calculated and printed, then the quasi-particle energies are printed out.	!
	!-----------------------------------------------------------------------------------------------!

	SUBROUTINE DiagonalizationMethod_goto_SelfConsistency(diagonal, tolerance)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		DOUBLE PRECISION, INTENT(IN) :: tolerance

		TYPE (OneDimSolve) :: neutron_constrainer, proton_constrainer
		TYPE (SymDensity) :: new_density
		TYPE (SymHartreeFockField) :: field1, field2

		LOGICAL :: Read_BBph, Read_BBpp, Read_Cph, Read_Cpp

		DOUBLE PRECISION :: accuracy, b, diff, factor
		DOUBLE PRECISION :: TestAccuracy = 1.0D0

		INTEGER :: A, N, Z, niter, ta
		INTEGER :: cycles_in, cycles_out, cycles_rate, Lmin

		! Reading the system's clock
		CALL SYSTEM_CLOCK(cycles_in)

		diff = 0.0D0

		b = Nucleus_get_b(diagonal%consistency%density%nucleus)
		A = Nucleus_get_A(diagonal%consistency%density%nucleus)

		IF ((A .LE. 1) .OR. (A .GE. 300)) THEN
			PRINT *, "Unexpected A value = ", A
			STOP "DiagonalizationMethod_goto_SelfConsistency"
		END IF

		factor = 1.0D0 - (1.0D0 / A)

		WRITE(6,'(5X,"STARTING THE SELF-CONSISTENT CALCULATION")')
		WRITE(6,'(5X,"========================================",/)')

		WRITE(6,'("Numerical precision in Energy ..: ",ES10.2)') tolerance
		WRITE(6,'("Oscillator length ..............: ",F14.11)') b
		WRITE(6,'("Atomic number A ................: ",i6)') A
		WRITE(6,'("Largest Number .................: ",E24.16)') HUGE(TestAccuracy)
		WRITE(6,'("Smallest Number ................: ",E24.16,/)') TINY(TestAccuracy)

		!-------------------------------------------------------!
		! Calculation of the Matrix elements of the Gogny force	!
		!-------------------------------------------------------!

		CALL SymKineticEnergy2Body_calculate(diagonal%consistency%vEkCMph, diagonal%consistency%vEkCMpp)

		! Calculation of the matrix elements of the Gogny force.
		! This is the most time- and memory-consuming task

		CALL SymVLSph_calculate(diagonal%consistency%vLSph)
		CALL SymVLSpp_calculate(diagonal%consistency%vLSpp)

		! This array takes a lot of space in memory: we don't need it from now on

		IF (Basis .EQ. 2) THEN
			IF (ALLOCATED(WaveDeri)) DEALLOCATE(WaveDeri)
		END IF

        	! Creating and allocating memory for arrays helping accelerate
		! the calculation of radial integrals.
		! Defined in module "ibb.f90"

		IF (Optimization .EQ. 1) THEN

			Lmin = Lmax

			Read_BBph = .FALSE.
			Read_BBpp = .FALSE.

		        IF (Basis .EQ. 2 .AND. (.NOT. Read_BBph .OR. .NOT. Read_BBpp)) THEN
                                CALL BesselTabularize(2*Lmax)
                        END IF

			CALL SymVBBph_calculate(diagonal%consistency%vBBph, Read_BBph, Lmin)
			CALL SymVBBpp_calculate(diagonal%consistency%vBBpp, Read_BBpp, Lmin)

		        IF (Basis .EQ. 2 .AND. (.NOT. Read_BBph .OR. .NOT. Read_BBpp)) THEN
                                CALL BesselFree()
                        END IF

			Read_Cph = .FALSE.
			Read_Cpp = .FALSE.

			CALL SymVCph_calculate(diagonal%consistency%vCph, Read_Cph, Lmin)
			CALL SymVCpp_calculate(diagonal%consistency%vCpp, Read_Cpp, Lmin)

		ELSE

			Lmin = 0

			IF (Basis .EQ. 2) THEN
				IF (Lmax < 10) THEN
					WRITE(diagonal%consistency%vBBph%filename, "(A,I1,A)") "data/vBB", Lmax, "ph_WS.txt"
				ELSE
					WRITE(diagonal%consistency%vBBph%filename, "(A,I2,A)") "data/vBB", Lmax, "ph_WS.txt"
				END IF
			ELSE
				IF (Lmax < 10) THEN
					WRITE(diagonal%consistency%vBBph%filename, "(A,I1,A)") "data/vBB", Lmax, "ph_HO.txt"
				ELSE
					WRITE(diagonal%consistency%vBBph%filename, "(A,I2,A)") "data/vBB", Lmax, "ph_HO.txt"
				END IF
			END IF

			Read_BBph = SymVBBph_read(diagonal%consistency%vBBph)
			Read_BBpp = SymVBBpp_read(diagonal%consistency%vBBpp)

        		IF (Basis .EQ. 2 .AND. (.NOT. Read_BBph .OR. .NOT. Read_BBpp)) THEN
                                CALL BesselTabularize(2*Lmax)
                        END IF

			CALL SymVBBph_calculate(diagonal%consistency%vBBph, Read_BBph, Lmin)
			CALL SymVBBpp_calculate(diagonal%consistency%vBBpp, Read_BBpp, Lmin)

		        IF (Basis .EQ. 2 .AND. (.NOT. Read_BBph .OR. .NOT. Read_BBpp)) THEN
                                CALL BesselFree()
                        END IF

			Read_Cph = SymVCph_read(diagonal%consistency%vCph)
			Read_Cpp = SymVCpp_read(diagonal%consistency%vCpp)

			CALL SymVCph_calculate(diagonal%consistency%vCph, Read_Cph, Lmin)
			CALL SymVCpp_calculate(diagonal%consistency%vCpp, Read_Cpp, Lmin)

		END IF

		! Creation of the tensor fields

		CALL SymHartreeFockField_new(field1)
		CALL SymHartreeFockField_new(field2)

		!-------------------------------------------------------!
		!     Starting the iterative process of convergence	!
		!-------------------------------------------------------!

		niter = 0

gsc:		DO
  			N = Nucleus_get_N(diagonal%consistency%density%nucleus)
 			Z = Nucleus_get_Z(diagonal%consistency%density%nucleus)

			A = N + Z

			IF (A .LE. 1) THEN
				PRINT *, "Unexpected A value = ", A
				STOP "DiagonalizationMethod_goto_SelfConsistency"
			END IF

			CALL DiagonalizationMethod_get_MeanField(diagonal, diagonal%S)

			! We obtain the Fermi energies Lambda_n and Lambda_p by considering the constraint over the number of particle

			IF (HFOnly .EQ. 0) THEN

				CALL DiagonalizationMethod_set_ISOSpin(diagonal, 1)

				CALL OneDimSolve_new(neutron_constrainer, diagonal%func, diagonal)

				diagonal%consistency%density%nucleus%lambda_np(1) = OneDimSolve_solve(neutron_constrainer, &
					DBLE(N), diagonal%consistency%density%nucleus%lambda_np(1), &
					0.5D0 * (diagonal%consistency%density%nucleus%np(1) - diagonal%consistency%density%nucleus%actual_np(1)) &
					+ 0.1D0)

				CALL DiagonalizationMethod_set_ISOSpin(diagonal, 0)

				CALL OneDimSolve_new(proton_constrainer, diagonal%func, diagonal)

				diagonal%consistency%density%nucleus%lambda_np(0) = OneDimSolve_solve(proton_constrainer, &
					DBLE(Z), diagonal%consistency%density%nucleus%lambda_np(0), &
					0.5D0 * (diagonal%consistency%density%nucleus%np(0) - diagonal%consistency%density%nucleus%actual_np(0)) &
					+ 0.1D0)

			ELSE

				CALL DiagonalizationMethod_set_ISOSpin(diagonal, 1)
				diagonal%consistency%density%nucleus%actual_np(1) = DiagonalizationMethod_operator(diagonal)

				CALL DiagonalizationMethod_set_ISOSpin(diagonal, 0)
				diagonal%consistency%density%nucleus%actual_np(0) = DiagonalizationMethod_operator(diagonal)

			END IF

			! Store the density for an eventual restart
			CALL SymDensity_save(diagonal%consistency%density)

			! Store the energies (for future use) of this step after correction for the particle number (above)
			CALL SelfConsistencyMethod_store_eHFB(diagonal%consistency)

			!-------------------------------------------------------------------------------------------------------!
			!													!
			!			   	ULTRA-IMPORTANT CALL							!
			!				--------------------							!
			!													!
			!    diagonal%iterated - TYPE: SymGenDensityProj   MODULE: symgden_proj.f90				!
			!    new_density       - TYPE: SymDensityProj      MODULE: symden_proj.f90				!
			!													!
			!  The purpose of the following CALL SymDensity_new_GenDensityProj() is actually to COPY the content of !
			!  the density matrix contained in diagonal%iterated, and obtained after the diagonalization of the HFB	!
			!  matrix from U and V vectors, into new_density, which will be used to calculate the field Gamma and 	!
			!  Delta.												!
			!			  										!
			!-------------------------------------------------------------------------------------------------------!

			CALL SymDensity_new_GenDensity(new_density, diagonal%iterated)

                        ! The current step corresponds to the density diagonal%consistency%density, the new one to new_density
			! We compute here the "distance" between the two.

			IF (MOD(niter, MAX_ITER) .EQ. 0) THEN
				diff = SymHartreeFockField_distance(diagonal%consistency%density%field%rho, new_density%field%rho)
				WRITE(*,'("          k = ",I4," Difference of the densities : ",F15.8)') niter,diff
			END IF

			niter = niter + 1

                        ! We reduce the step made into the direction of the new density, new_density, from the old one,
			! diagonal%consistency%density, with some slowing-factor anneal
			!
			!    rho_{n+1} -> (alpha * rho_{n} + rho_{n+1})/(1 + alpha)
			!

			factor = 1.0D0 / (1.0D0 + anneal)

			! new step for normal density

			CALL SymHartreeFockField_product(field1, anneal, diagonal%consistency%density%field%rho)
			CALL SymHartreeFockField_product(field2, 1.0D0-anneal, new_density%field%rho)
			CALL SymHartreeFockField_add(diagonal%consistency%density%field%rho, field2, field1)

			IF (HFonly.EQ.0) THEN

				CALL SymHartreeFockField_product(field1, anneal, diagonal%consistency%density%field%kap)
				CALL SymHartreeFockField_add(field2, new_density%field%kap, field1)
				CALL SymHartreeFockField_product(diagonal%consistency%density%field%kap, factor, field2)

				! We update the numbers of particles with the same rule
				DO ta = 0, 1
					diagonal%consistency%density%nucleus%actual_np(ta) = factor * ((new_density%nucleus%actual_np(ta)) &
						+ (anneal * diagonal%consistency%density%nucleus%actual_np(ta)))

				END DO

			END IF

			accuracy = SelfConsistencyMethod_accuracy(diagonal%consistency)
			WRITE(*,'("Iteration k = ",I4," HF Energy : ",F15.8)') niter,diagonal%consistency%density%nucleus%eHFB

			IF (diff .LE. tolerance) THEN
				accuracy = SelfConsistencyMethod_accuracy(diagonal%consistency)
                                WRITE(*,'("Energy = ",EN15.5," Precision = ",ES12.5)') &
                                          diagonal%consistency%density%nucleus%eHFB,accuracy
				IF (accuracy .LE. tolerance) THEN
					EXIT gsc
				END IF
			END IF

			IF (niter .GE. NITER_MAX) THEN
				EXIT gsc
		        END IF

		END DO gsc

		CALL SymHartreeFockField_del(field1)
		CALL SymHartreeFockField_del(field2)

		! Store the density for an eventual restart
		CALL SymDensity_save(diagonal%consistency%density)

		CALL SYSTEM_CLOCK(cycles_out, cycles_rate)
		PRINT "(/A,EN10.2)", "Elapsed time (seconds): ", (DBLE(cycles_out - cycles_in) / cycles_rate)

		CALL SymDensity_store_actual_R2(diagonal%consistency%density)

		CALL SelfConsistencyMethod_store_eHFB(diagonal%consistency)
		CALL SelfConsistencyMethod_show_Status(diagonal%consistency)

		CALL DiagonalizationMethod_show_ParticleEnergies(diagonal)
                CALL DiagonalizationMethod_show_QuasiParticleEnergies(diagonal)

		RETURN
	END SUBROUTINE DiagonalizationMethod_goto_SelfConsistency

	!-----------------------------------------------------------------------------------------------!
	!   In this subroutine, we calculate the mean-field Gamma and pairing field Delta from the 	!
	!   matrix elements of the force and the densities. Both the matrix elements and the densities	!
	!   at a given iteration of the HFB process are stored in the object "diagonal". Once the fields!
	!   are constructed, we build an object that will contain the super-hamiltonian S.
	!-----------------------------------------------------------------------------------------------!

	SUBROUTINE DiagonalizationMethod_get_MeanField(diagonal, S)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		TYPE (SymGenDensity), INTENT(INOUT) :: S

		DOUBLE PRECISION :: b, factor
		INTEGER :: A

		TYPE (SymHartreeFockField) :: HF_Gamma, HF_Delta
		TYPE (SymHartreeFockField) :: ekcm_field, vbb_field, vc_field, vls_field, gdd_field
		TYPE (SymHartreeFockField) :: field1, field2
		TYPE (SymD3Tensor) :: ek_tensor
		TYPE (SymGenDensityHF) :: gendenhf_gamma, gendenhf_delta

		b = Nucleus_get_b(diagonal%consistency%density%nucleus)
		A = Nucleus_get_A(diagonal%consistency%density%nucleus)

		factor = 1.0D0 - (1.0D0 / A)

		IF ((A .LE. 1) .OR. (A .GE. 300)) THEN
			PRINT *, "Unexpected A value = ", A
			STOP "DiagonalizationMethod::get_MeanField"
		END IF

		! Creating all the necessary fields

		CALL SymHartreeFockField_new(HF_Gamma)
		CALL SymHartreeFockField_new(HF_Delta)

		CALL SymD3Tensor_new(ek_tensor)

		CALL SymHartreeFockField_new(ekcm_field)
		CALL SymHartreeFockField_new(vbb_field)
		CALL SymHartreeFockField_new(vc_field)
		CALL SymHartreeFockField_new(vls_field)
		CALL SymHartreeFockField_new(gdd_field)

		CALL SymHartreeFockField_new(field1)
		CALL SymHartreeFockField_new(field2)

		! Mean-field - Kinetic energy

                CALL SymD3Tensor_product(ek_tensor, DBLE(factor), EkField)

		! Mean-field - Two-body center of mass correction

		CALL SymKineticEnergy2Body_get_Gamma(field1, diagonal%consistency%vEkCMph, diagonal%consistency%density%field%rho)
		IF (switch_CM .GE. 1) THEN
			CALL SymHartreeFockField_product(ekcm_field, DBLE(1.0 / A), field1)
		ELSE
			CALL SymHartreeFockField_product(ekcm_field, DBLE(0.0), field1)
		END IF

		! Mean-field - Brink-Boker term

		CALL SymVBBph_get_Gamma(vbb_field, diagonal%consistency%vBBph, diagonal%consistency%density%field%rho)

		! Mean-field - Coulomb potential

		CALL SymVCph_get_Gamma(field1, diagonal%consistency%vCph, diagonal%consistency%density%field%rho)
		IF (switch_Coulomb .GE. 1) THEN
			CALL SymHartreeFockField_product(vc_field, DBLE(1.0), field1)
		ELSE
			CALL SymHartreeFockField_product(vc_field, DBLE(0.0), field1)
		END IF

		! Mean-field - Spin-orbit term

		CALL SymVLSph_get_Gamma(field1, diagonal%consistency%vLSph, diagonal%consistency%density%field%rho)
		IF (switch_LS .GE. 1) THEN
			CALL SymHartreeFockField_product(vls_field, DBLE(1.0), field1)
		ELSE
			CALL SymHartreeFockField_product(vls_field, DBLE(0.0), field1)
		END IF

		! Mean-field - Density-dependent term

		CALL SymGDDph_update(field1, diagonal%consistency%gDDph, diagonal%consistency%density%field%rho)
		IF (switch_DD .GE. 1) THEN
			CALL SymHartreeFockField_product(gdd_field, DBLE(1.0), field1)
		ELSE
			CALL SymHartreeFockField_product(gdd_field, DBLE(0.0), field1)
		END IF

		! Total Mean-field = Sum of all the preceding terms

		CALL SymHartreeFockField_add(field1, vls_field, gdd_field)
		CALL SymHartreeFockField_add(field2, vc_field, field1)
		CALL SymHartreeFockField_add(field1, vbb_field, field2)
		CALL SymHartreeFockField_add(field2, ekcm_field, field1)
		CALL SymHartreeFockField_add_SymD3Tensor(HF_Gamma, ek_tensor, field2)

		! Pairing - Center of Mass term

		IF (HFOnly .EQ. 0 .AND. switch_CM .GE. 1) THEN
			CALL SymKineticEnergy2Body_get_Delta(field1, diagonal%consistency%vEkCMpp, diagonal%consistency%density%field%kap)
			CALL SymHartreeFockField_product(ekcm_field, DBLE(1.0 / A), field1)
		ELSE
			CALL SymKineticEnergy2Body_get_Delta(field1, diagonal%consistency%vEkCMpp, diagonal%consistency%density%field%kap)
			CALL SymHartreeFockField_product(ekcm_field, DBLE(0.0), field1)
		END IF

		! Pairing - Brink-Boker term

		CALL SymVBBpp_get_Delta(field1, diagonal%consistency%vBBpp, diagonal%consistency%density%field%kap)
		IF (HFOnly .EQ. 0) THEN
		        CALL SymHartreeFockField_product(vbb_field, DBLE(1.0), field1)
		ELSE
			CALL SymHartreeFockField_product(vbb_field, DBLE(0.0), field1)
		END IF

		! Pairing - Coulomb term

		CALL SymVCpp_get_Delta(field1, diagonal%consistency%vCpp, diagonal%consistency%density%field%kap)
		IF (HFOnly .EQ. 0 .AND. switch_Coulomb .GE. 1) THEN
		        CALL SymHartreeFockField_product(vc_field, DBLE(1.0), field1)
		ELSE
			CALL SymHartreeFockField_product(vc_field, DBLE(0.0), field1)
		END IF

		! Pairing - Spin-orbit term

		CALL SymVLSpp_get_Delta(field1, diagonal%consistency%vLSpp, diagonal%consistency%density%field%kap)
		IF (HFOnly .EQ. 0 .AND. switch_LS .GE. 1) THEN
		        CALL SymHartreeFockField_product(vls_field, DBLE(1.0), field1)
		ELSE
			CALL SymHartreeFockField_product(vls_field, DBLE(0.0), field1)
		END IF

		! Total Pairing = Sum of all the preceding terms

		CALL SymHartreeFockField_add(field1, vc_field, vls_field)
		CALL SymHartreeFockField_add(field2, ekcm_field, field1)
		CALL SymHartreeFockField_add(HF_Delta, vbb_field, field2)

		! Creating an object that contain the super-hamiltonian to diagonalize

		CALL SymGenDensityHF_new(gendenhf_gamma)
		CALL SymGenDensityHF_new(gendenhf_delta)

		CALL SymGenDensityHF_copy(gendenhf_gamma, HF_Gamma)
		CALL SymGenDensityHF_copy(gendenhf_delta, HF_Delta)

		CALL SymGenDensity_new_GammaDelta(S, gendenhf_gamma, gendenhf_delta, b)

		CALL SymGenDensityHF_del(gendenhf_gamma)
		CALL SymGenDensityHF_del(gendenhf_delta)

		CALL SymD3Tensor_del(ek_tensor)

		CALL SymHartreeFockField_del(HF_Gamma)
		CALL SymHartreeFockField_del(HF_Delta)

		CALL SymHartreeFockField_del(ekcm_field)
		CALL SymHartreeFockField_del(vbb_field)
		CALL SymHartreeFockField_del(vc_field)
		CALL SymHartreeFockField_del(vls_field)
		CALL SymHartreeFockField_del(gdd_field)

		CALL SymHartreeFockField_del(field1)
		CALL SymHartreeFockField_del(field2)

		RETURN
	END SUBROUTINE DiagonalizationMethod_get_MeanField

	! Changing the isospin under consideration

	SUBROUTINE DiagonalizationMethod_set_ISOSpin(diagonal, ta)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal
		INTEGER, INTENT(IN) :: ta

		diagonal%ta = ta

		RETURN
	END SUBROUTINE DiagonalizationMethod_set_ISOSpin

	!-----------------------------------------------------------------------------------------------!
	!   In this subroutine, we calculate the particle number. This is used to adjust the latter to	!
	!   the actual number of particles of the system at each iteration.				!
	!-----------------------------------------------------------------------------------------------!

	FUNCTION DiagonalizationMethod_operator(diagonal)
		DOUBLE PRECISION DiagonalizationMethod_operator
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal

                DOUBLE PRECISION, DIMENSION(1:Nmax,1:Nmax) :: copyR2
		DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: sp_energy, esp, occup
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: copyRho, tmp
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: xI, h, Delta, SH, U, V

		DOUBLE PRECISION :: np, R2, b2, trace, Factor, E_curr, E_next

		INTEGER :: a, d, dd, la, ja, i, icount, ncount, found
		INTEGER :: u1, u2
		INTEGER :: Isospin, dim_cur, Npart, diffPart

		INTEGER, DIMENSION(:), ALLOCATABLE :: in
		INTEGER, DIMENSION(:), ALLOCATABLE :: sp_block
		INTEGER, DIMENSION(:,:), ALLOCATABLE :: sp_state

		np = 0.0D0
		R2 = 0.0D0
		b2 = Nucleus_get_b(diagonal%iterated%nucleus) ** 2

                Isospin = diagonal%ta

                IF (IsoFacPair .EQ. Isospin .OR. IsoFacPair .EQ. 2) THEN
                        Factor = facPair
                ELSE
                        Factor = 1.0D0
                END IF

                IF (HFonly.EQ.1) THEN
			dim_cur = 0
			DO a = 0, 2*Lmax
				dim_cur = dim_cur + DIM(a)
			END DO
			ALLOCATE(sp_energy(dim_cur))
			ALLOCATE(sp_block(dim_cur))
			ALLOCATE(sp_state(0:2*Lmax,1:Nmax))
			ALLOCATE(occup(dim_cur))
		END IF

		! Block structure: a refers to the angular momentum, but remember we deal here with a 2N x 2N matrix,
		! if N is the size of the single-particle basis.
		icount = 0

		DO a = 0, 2*Lmax

			d = DIM(a)
			dd = d + d

			la = L(a)
			ja = J(a)

			ALLOCATE(xI(d, d))
			ALLOCATE(h(d, d))
			ALLOCATE(Delta(d, d))
			ALLOCATE(U(d, d))
			ALLOCATE(V(d, d))
			ALLOCATE(tmp(d, d))
			ALLOCATE(copyRho(d, d))

			! xI contains the Fermi level lambda, which was determined before
			! from the conservation of particle number

			IF (HFonly.EQ.1) THEN

				ALLOCATE(SH(d, d))
				ALLOCATE(esp(d))

				SH = diagonal%S%rho%rho(diagonal%ta, a)%store

				! Diagonalizing the hamiltonian SH gives the single-particle energies
				! ordered from the lowest to the highest (in this block)
				CALL EigenValues(d, d, SH, esp, V)

				! Storing the s.p. energies, block identification and state number
				DO u1 = 1, d
					icount = icount + 1
					sp_energy(icount) = esp(u1)
					sp_block(icount)  = a
					sp_state(a, u1)   = icount
				END DO

				! Extracting the matrices U and V
				DO u1 = 1, d
					diagonal%QuasiParticleEnergies(diagonal%ta, a)%value(u1) = esp(u1)
					DO u2 = 1, d
						diagonal%UV(diagonal%ta, a)%quantum(u1,  u2  ) = V(u1, u2)
						diagonal%UV(diagonal%ta, a)%quantum(u1+d,u2  ) = 0.0D0
						diagonal%UV(diagonal%ta, a)%quantum(u1,  u2+d) = 0.0D0
						diagonal%UV(diagonal%ta, a)%quantum(u1+d,u2+d) = 0.0D0
					END DO
				END DO

				DEALLOCATE(esp)

			ELSE

				ALLOCATE(SH(dd, dd))

				xI = 0.0D0
				DO u1 = 1, d
					xI(u1, u1) = diagonal%func%x
				END DO

				h     = diagonal%S%rho%rho(diagonal%ta, a)%store - xI
				Delta = diagonal%S%kap%rho(diagonal%ta, a)%store

				! In spherical symmetry, the HFB matrix to diagonalize has the form:
				!
				!		    h      -Delta
				!		   -Delta    -h
				!
				DO u1 = 1, d
					DO u2 = 1, d
						SH(u1    , u2    ) =   h(u1, u2)
						SH(u1    , u2 + d) = - Factor * Delta(u1, u2)
						SH(u1 + d, u2    ) = - Factor * Delta(u1, u2)
						SH(u1 + d, u2 + d) = - h(u1, u2)
					END DO
				END DO

				! Diagonalizing the super-hamiltonian SH gives the quasi-particle energies
				! ordered from the lowest to the highest
				CALL EigenValues(dd, dd, SH, &
					diagonal%QuasiParticleEnergies(diagonal%ta, a)%value, &
					diagonal%UV(diagonal%ta, a)%quantum)

				! Extracting the matrices U and V
				DO u1 = 1, d
					DO u2 = 1, d
						U(u1, u2) = diagonal%UV(diagonal%ta, a)%quantum(u1    , u2 + d)
						V(u1, u2) = diagonal%UV(diagonal%ta, a)%quantum(u1 + d, u2 + d)
					END DO
				END DO

				! We calculate the new densities rho and kappa from the U and V matrices
				CALL SymGenDensity_make_Block(diagonal%iterated, diagonal%ta, a, U, V)

			END IF

			! The number of particles is the trace of the rho matrix. Don't forget about the
			! degeneracy 2j+1 of a spherical shell (here ja = 2*j to avoid fractions)
			trace = 0.0D0
			DO u1 = 1, d
				trace = trace + diagonal%iterated%rho%rho(diagonal%ta, a)%store(u1, u1)
			END DO
			np = np + DBLE(ja + 1) * trace

			! We need to copy the matrix R2Field because of dimensions inconsistencies in the
			! case of an arbitrary basis. This inconsistency causes the multiplication below to
			! crash when using -C option at compilation stage
			copyR2 = SymD3Tensor_matrix(R2Field, la)
			copyRho = diagonal%iterated%rho%rho(diagonal%ta, a)%store

			tmp = MATMUL(copyR2(1:d,1:d), copyRho)

			trace = 0.0D00
			DO u1 = 1, d
				trace = trace + tmp(u1, u1)
			END DO
			R2 = R2 + DBLE(ja + 1) * trace

			! Freeing memory for the matrix we used
			DEALLOCATE(xI)
			DEALLOCATE(h)
			DEALLOCATE(Delta)
			DEALLOCATE(SH)
			DEALLOCATE(U)
			DEALLOCATE(V)
			DEALLOCATE(tmp)
			DEALLOCATE(copyRho)

		END DO

		IF (HFOnly .EQ. 1) THEN

			! COMPUTING THE FERMI LEVEL

			ALLOCATE(in(dim_cur))

			DO i = 1, dim_cur
				in(i) = i
                	END DO

			CALL indexx_real8(dim_cur, sp_energy, in)

			Npart = diagonal%iterated%nucleus%np(Isospin)

			icount = 0; ncount = 0; found = 0

			DO i=1, dim_cur

				a = sp_block(in(i))
				ja = J(a)
				occup(in(i)) = 1.0D0

				icount = icount + ja + 1

				IF (icount.GT.Npart) THEN
					IF (found.EQ.0) THEN
						found = 1
						diffPart = Npart - icount
						E_curr = sp_energy(in(i))
						E_next = sp_energy(in(i+1))
						diagonal%consistency%density%nucleus%lambda_np(diagonal%ta) = E_curr &
						         + (E_next - E_curr)*DBLE(diffPart)*occup(in(i)) / DBLE(ja + 1)
						occup(in(i)) = 0.0D0
						IF (diffPart.NE.0) THEN
							occup(in(i)) = ABS(DBLE(ncount - Npart)/DBLE(ja + 1))
						END IF
					ELSE
						occup(in(i)) = 0.0D0
					END IF
				END IF
				ncount = icount
			END DO

			! Setting the occupation probabilities (per block) based on simplest vacuum configuration

			DO a = 0, 2*Lmax

				ja = J(a)
				d = DIM(a)
				ALLOCATE(V(d, d))

				! Filters the eigenvectors of the Hamiltonian by their occupation
				DO u1 = 1, d
					DO u2 = 1, d
						i = sp_state(a, u2)
						V(u1,u2) = diagonal%UV(diagonal%ta, a)%quantum(u1,u2)*SQRT(occup(i))
						diagonal%UV(diagonal%ta, a)%quantum(u1,u2) = V(u1,u2)
					END DO
				END DO

				! We calculate the new densities rho and kappa from the eigenvectors of hf with proper occuations
				CALL SymGenDensity_make_HF(diagonal%iterated, diagonal%ta, a, V)

				trace=0.0D0
				DO u1 = 1, d
					trace = trace + diagonal%iterated%rho%rho(diagonal%ta, a)%store(u1, u1)*DBLE(ja+1)
				END DO

				DEALLOCATE(V)

			END DO

			DEALLOCATE(in,sp_energy)
			DEALLOCATE(occup,sp_block,sp_state)

		END IF

		! When we dealt with all (l,j)-channels, we have the final number of particles and radius.
		diagonal%iterated%nucleus%actual_np(diagonal%ta) = np
		diagonal%iterated%nucleus%actual_R2(diagonal%ta) = R2 / np
		diagonal%consistency%density%nucleus%actual_np(diagonal%ta) = np
		diagonal%consistency%density%nucleus%actual_R2(diagonal%ta) = R2 / np

		! Resultado final
		DiagonalizationMethod_operator = np

		RETURN
	END FUNCTION DiagonalizationMethod_operator

	!-------------------------------------------------------------------------------!
	!										!
	!       Diagonalization of the matrix of the normal density rho. 		!
	!										!
	!  The single-particle energies are defined, in the case of HFB, as the		!
	!  expectation values of the hamiltonian in the basis where rho is diagonal.	!
	!  This basis is referred to as the canonical basis.				!
	!										!
	!-------------------------------------------------------------------------------!

	SUBROUTINE DiagonalizationMethod_show_ParticleEnergies(diagonal)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal

		TYPE (SymGenDensity) :: SuperHamiltonian

		TYPE (DiagonalMatrix), DIMENSION(0:1) :: D, Occup, SingleParticleE, SingleParticleP, SingleParticleR2, EqpCanonical
		TYPE (MatrixType), DIMENSION(0:1) :: P, E, R2

		INTEGER :: ndim, a, num, i, dim_cur, dim_acc, la, ja, ta, nrot, file_error
		INTEGER :: Lvalue, Jvalue, ipoint, IndexBra, IndexKet, m, n
		INTEGER, PARAMETER :: file_unit_1 = 16, file_unit_2 = 17
		INTEGER, DIMENSION(0:1) :: dim_mat
		INTEGER, ALLOCATABLE :: an(:), ap(:), in(:), ip(:)
		INTEGER, ALLOCATABLE :: Count(:,:)

		DOUBLE PRECISION :: b, e_cano, d_cano
		DOUBLE PRECISION, DIMENSION(0:1) :: normDelta, averageDelta, eFermi
		DOUBLE PRECISION, ALLOCATABLE :: S(:, :), V(:, :), R2Matrix(:, :), R2Initial(:, :)
		DOUBLE PRECISION, ALLOCATABLE :: RadialWaveFunction(:, :, :)

		CHARACTER(8) :: label
		CHARACTER(80) :: file_neut, file_prot
		CHARACTER, DIMENSION(0:29) :: spectr

		DATA spectr / "s", "p", "d", "f", "g", "h", "i", "j", "k", "l", &
			      "m", "n", "o", "P", "q", "r", "S", "t", "u", "v", &
			      "w", "x", "y", "z", "?", "?", "?", "?", "?", "?" /

		b = Nucleus_get_b(diagonal%consistency%density%nucleus)

		dim_cur = 0
		DO a = 0, 2*Lmax
			dim_cur = dim_cur + DIM(a)
		END DO

		ALLOCATE(Occup(0)%value(dim_cur))
		ALLOCATE(Occup(1)%value(dim_cur))

		ALLOCATE(SingleParticleE(0)%value(dim_cur))
		ALLOCATE(SingleParticleE(1)%value(dim_cur))
		ALLOCATE(SingleParticleP(0)%value(dim_cur))
		ALLOCATE(SingleParticleP(1)%value(dim_cur))
		ALLOCATE(SingleParticleR2(0)%value(dim_cur))
		ALLOCATE(SingleParticleR2(1)%value(dim_cur))

		ALLOCATE(EqpCanonical(0)%value(dim_cur))
		ALLOCATE(EqpCanonical(1)%value(dim_cur))

		ALLOCATE(RadialWaveFunction(0:1, 0:Npoint, 1:dim_cur))

		IF (.NOT. ALLOCATED(DensityRadial)) ALLOCATE(DensityRadial(0:1, 0:Npoint))

		! "SuperHamiltonian" contains the generalized hamiltonian
		!
		!                 (     h      -Delta  )
		!                 (                    )
		!                 (  -Delta     -h     )
		!
		!  SuperHamiltonian%rho = h (type of SymHartreeFockField)
		!  SuperHamiltonian%kap = Delta (type of SymHartreeFockField)

		CALL DiagonalizationMethod_get_MeanField(diagonal, SuperHamiltonian)

							ndim = Nmax * (2*Lmax + 1)
		IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	ndim = (N_0 + 1) * (N_0 + 2) / 2

		ALLOCATE(an(ndim + 1))
		ALLOCATE(ap(ndim + 1))

		DO ta = 0, 1

			eFermi(ta) = diagonal%consistency%density%nucleus%lambda_np(ta)

			dim_acc = 0
			num = 1

			normDelta(ta) = 0.0D0
			averageDelta(ta) = 0.0D0

			! Initialization of the radial density

			DO ipoint = 0, Npoint
				DensityRadial(ta, ipoint) = 0.0D0
			END DO

			DO a = 0, 2*Lmax

				dim_cur = DIM(a)

				la = L(a)
				ja = J(a)

				DO i = 1, dim_cur
					an(num) = a
					ap(num) = a
					num = num + 1
				END DO

				ALLOCATE(S(dim_cur, dim_cur))
				ALLOCATE(V(dim_cur, dim_cur))
				ALLOCATE(R2Matrix(dim_cur, dim_cur))
				ALLOCATE(R2Initial(Nmax, Nmax))

				ALLOCATE(D(ta)%value(dim_cur))
				ALLOCATE(E(ta)%quantum(dim_cur, dim_cur))
				ALLOCATE(P(ta)%quantum(dim_cur, dim_cur))
				ALLOCATE(R2(ta)%quantum(dim_cur, dim_cur))

				! S contains the matrix of the density rho for the quantum number a = (la,ja) at the convergence

				S = diagonal%iterated%rho%rho(ta, a)%store

				!  After diagonalization, D(ta)% value contains the diagonal values of the p-h density rho
				!  (the occupation probabilities v2), V the eigenvectors that make pass from the HFB basis
				!  to the canonical basis (where the matrix of rho is diagonal).
				!
				!  The s.p. energies are by definition the expectation value of the hamiltonian in the
				!  canonical basis. Therefore they read: (V+)*S*V

				CALL Jacobi_real8(S, dim_cur, D(ta)%value, V, nrot)

				E(ta)%quantum = MATMUL(TRANSPOSE(V), MATMUL(SuperHamiltonian%rho%rho(ta, a)%store, V))
				P(ta)%quantum = MATMUL(TRANSPOSE(V), MATMUL(SuperHamiltonian%kap%rho(ta, a)%store, V))

				R2Initial = SymD3Tensor_matrix(R2Field, L(a))
				R2Matrix = R2Initial(1:dim_cur,1:dim_cur)

				R2(ta)%quantum = MATMUL(TRANSPOSE(V), MATMUL(R2Matrix, V))

				DO i = 1, dim_cur

					Occup(ta)%value(dim_acc + i) = D(ta)%value(i)

					SingleParticleE(ta)%value(dim_acc + i) = E(ta)%quantum(i, i)
					SingleParticleP(ta)%value(dim_acc + i) = P(ta)%quantum(i, i)
					SingleParticleR2(ta)%value(dim_acc + i) = R2(ta)%quantum(i, i)

					! Energy of QP in canonical basis (BCS-likei) - Sign '-' comes only to sort by increasing order
					e_cano = SingleParticleE(ta)%value(dim_acc + i) - eFermi(ta)
					d_cano = SingleParticleP(ta)%value(dim_acc + i)

					EqpCanonical(ta)%value(dim_acc + i) = sqrt(e_cano*e_cano + d_cano*d_cano)

					averageDelta(ta) = averageDelta(ta) - P(ta)%quantum(i, i) * D(ta)%value(i) * (1.0D0 - D(ta)%value(i))
					normDelta(ta) = normDelta(ta) + D(ta)%value(i) * (1.0D0 - D(ta)%value(i))

					! We calculate the contribution to the radial density of this (la, ja) state.
					! The state has the degeneracy 2j + 1, hence the factor ja+1 (ja=2j by definition).

					IF (Basis .EQ. 2) THEN

						DO ipoint = 1, Npoint

							RadialWaveFunction(ta, ipoint, dim_acc + i) = 0.0D0

							DO m = 1, dim_cur

								IndexBra = IndexVecNL(m, la)

								RadialWaveFunction(ta, ipoint, dim_acc + i) = &
								 	 RadialWaveFunction(ta, ipoint, dim_acc + i) &
									+ V(m, i)*WaveFun(ipoint,IndexBra) / RadMesh(ipoint)

								DO n = 1, dim_cur

									IndexKet = IndexVecNL(n, la)

									DensityRadial(ta, ipoint) = DensityRadial(ta, ipoint) &
										+ D(ta)%value(i) * (ja + 1.0D0) &
										* V(m, i) * V(n, i) &
										* WaveFun(ipoint,IndexBra) &
										* WaveFun(ipoint,IndexKet)

								END DO

							END DO

						END DO

					END IF

				END DO

				DEALLOCATE(D(ta)%value)
				DEALLOCATE(E(ta)%quantum)
				DEALLOCATE(P(ta)%quantum)
				DEALLOCATE(R2(ta)%quantum)

				DEALLOCATE(S)
				DEALLOCATE(V)
				DEALLOCATE(R2Matrix)
				DEALLOCATE(R2Initial)

				dim_acc = dim_acc + dim_cur

			END DO ! end of loop over a

			IF (Basis .EQ. 2) THEN

				! The basis wave-functions are in fact y(r) = R(r)/r, so we need to correct for this.
				! The 4*pi comes from the integration over the angles.

				DO ipoint=1, Npoint
					DensityRadial(ta, ipoint) = DensityRadial(ta, ipoint) / (4.0D0*PI*RadMesh(ipoint)**2)
				END DO

				! Extrapolate to find rho(r=0)

         			DensityRadial(ta, 0) = 3.0D0*(DensityRadial(ta, 1) - DensityRadial(ta, 2)) + DensityRadial(ta, 3)

			END IF

                        dim_mat(ta) = dim_acc

		END DO ! end of loop over ta

		! Sorting the single-particle energies by ascending order
		! in() and ip() are the vectors that contain the information on the ordering

		ALLOCATE(in(dim_mat(1)))
		ALLOCATE(ip(dim_mat(0)))

		DO i = 1, dim_mat(1)
			in(i) = i
                END DO
		DO i = 1, dim_mat(0)
                        ip(i) = i
		END DO

		CALL indexx_real8(dim_mat(1), SingleParticleE(1)%value, in)
		CALL indexx_real8(dim_mat(0), SingleParticleE(0)%value, ip)

		IF (Basis .EQ. 2) THEN

			! Total density (proton and neutron)

			OPEN(file_unit_1, FILE='data/DensityCB.dat', ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				WRITE(*,'("Impossible to open the file data/DensityCB.dat")')
				STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
			END IF

			DO ipoint = 0, Npoint
				WRITE(file_unit_1,'(3f20.16)') RadMesh(ipoint),DensityRadial(0, ipoint),DensityRadial(1, ipoint)
			END DO

			CLOSE(file_unit_1)

			! Radial wavefunctions. The user specifies which one to be stored by the index IndexWave

			OPEN(file_unit_1, FILE='data/WavesCB.dat', ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				WRITE(*,'("Impossible to open the file data/WavesCB.dat")')
				STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
			END IF

			RadialWaveFunction(0, 0, ip(1)) = 3.0D0*( RadialWaveFunction(0, 1, ip(IndexWave)) &
								- RadialWaveFunction(0, 2, ip(IndexWave))) &
								+ RadialWaveFunction(0, 3, ip(IndexWave))
			RadialWaveFunction(1, 0, in(1)) = 3.0D0*( RadialWaveFunction(1, 1, in(IndexWave)) &
								- RadialWaveFunction(1, 2, in(IndexWave))) &
								+ RadialWaveFunction(1, 3, in(IndexWave))

			DO ipoint = 0, Npoint
				WRITE(file_unit_1,'(3f20.16)') RadMesh(ipoint), RadialWaveFunction(0, ipoint, ip(IndexWave)), &
									RadialWaveFunction(1, ipoint, in(IndexWave))
			END DO

			CLOSE(file_unit_1)

		END IF

		! Storing single-particle wave-functions in the canonical basis for the neutrons

		WRITE(file_neut,'("data/HF_sp_n.dat")')

		OPEN(file_unit_1, FILE=file_neut, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			WRITE(*,'("Impossible to open the file ",A)') file_neut
			STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
		END IF

		write(*,'("Average Delta for neutrons: ",f15.10)') averageDelta(1)/normDelta(1)

		WRITE(*,'()')
		WRITE(*,'(5X,"NEUTRON SINGLE-PARTICLE ENERGIES (CANONICAL BASIS)")')
		WRITE(*,'(5X,"==================================================")')

		ALLOCATE(Count(0:29,0:100))

		DO Lvalue = 0, 29
			DO Jvalue = 0, 100
				Count(Lvalue, Jvalue) = 0
			END DO
		END DO

		i = 1

		DO WHILE (SingleParticleE(1)%value(in(i)) .GE. -100.0D0 .AND. SingleParticleE(1)%value(in(i)) .LE. 40.0D0 .AND. i .LT. dim_mat(1))

			Lvalue = L(an(in(i)))
			Jvalue = J(an(in(i)))

			Count(Lvalue, Jvalue) = Count(Lvalue, Jvalue) + 1

			IF (Jvalue .GT. 9) THEN
				WRITE(label,'(2x,I1,A1,I2,"/2")') Count(Lvalue, Jvalue), spectr(Lvalue), Jvalue
			ELSE
				WRITE(label,'(3x,I1,A1,I1,"/2")') Count(Lvalue, Jvalue), spectr(Lvalue), Jvalue
			END IF

                        WRITE(*,'(I4,")",2X,A2,I2,"/2",3x,F10.7,F10.3,F8.3,3F8.3)') &
					i, spectr(Lvalue), Jvalue, Occup(1)%value(in(i)), SingleParticleE(1)%value(in(i)), &
											  SingleParticleP(1)%value(in(i)), &
											  EqpCanonical(1)%value(in(i))
                        WRITE(file_unit_1,'(i4,3x,a8,f15.8)') i, label, SingleParticleE(1)%value(in(i))

			i = i+1

		END DO

		CLOSE(file_unit_1)

		! Storing single-particle wave-functions in the canonical basis for the protons

		WRITE(file_prot,'("data/HF_sp_p.dat")')

		OPEN(file_unit_2, FILE=file_prot, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			WRITE(*,'("Impossible to open the file ",A)') file_prot
			STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
		END IF

		write(*,'("Average Delta for protons: ",f15.10)') averageDelta(0)/normDelta(0)

		WRITE(*,'()')
		WRITE(*,'(5X,"PROTON SINGLE-PARTICLE ENERGIES (CANONICAL BASIS)")')
		WRITE(*,'(5X,"=================================================")')

		DO Lvalue = 0, 29
			DO Jvalue = 0, 100
				Count(Lvalue, Jvalue) = 0
			END DO
		END DO

		i = 1

		DO WHILE (SingleParticleE(0)%value(ip(i)) .GE. -100.0D0 .AND. SingleParticleE(0)%value(ip(i)) .LE. 40.0D0 .AND. i .LT. dim_mat(0))

			Lvalue = L(ap(ip(i)))
			Jvalue = J(ap(ip(i)))

			Count(Lvalue, Jvalue) = Count(Lvalue, Jvalue) + 1

			IF (Jvalue .GT. 9) THEN
				WRITE(label,'(2x,I1,A1,I2,"/2")') Count(Lvalue, Jvalue), spectr(Lvalue), Jvalue
			ELSE
				WRITE(label,'(3x,I1,A1,I1,"/2")') Count(Lvalue, Jvalue), spectr(Lvalue), Jvalue
			END IF
                        WRITE(*,'(I4,")",2X,A2,I2,"/2",3X,F10.7,F10.3,F8.3,2F8.3)') &
					i, spectr(Lvalue), Jvalue, Occup(0)%value(ip(i)), SingleParticleE(0)%value(ip(i)), &
											  SingleParticleP(0)%value(ip(i))
                        WRITE(file_unit_2,'(i4,3x,a8,f15.8)') i, label, SingleParticleE(0)%value(ip(i))

			i = i+1

		END DO

		CLOSE(file_unit_2)

                ! Obtaining the pairing gap of the lowest q.p (in the
                ! canonical basis)

		DO i = 1, dim_mat(1)
			in(i) = i
                END DO
		DO i = 1, dim_mat(0)
                        ip(i) = i
		END DO

		CALL indexx_real8(dim_mat(1), EqpCanonical(1)%value, in)
		CALL indexx_real8(dim_mat(0), EqpCanonical(0)%value, ip)

                WRITE(*,'("Neutrons - Lowest QP (canonical basis) = ",f14.8," Delta = ",f10.5," in = ",i4)') &
                                     EqpCanonical(1)%value(in(1)), SingleParticleP(1)%value(in(1)), in(1)
                WRITE(*,'("Protons  - Lowest QP (canonical basis) = ",f14.8," Delta = ",f10.5," ip = ",i4)') &
                                     EqpCanonical(0)%value(ip(1)), SingleParticleP(0)%value(ip(1)), ip(1)

		DEALLOCATE(in)
		DEALLOCATE(ip)

		DEALLOCATE(Count)

		DEALLOCATE(Occup(0)%value)
		DEALLOCATE(Occup(1)%value)

		DEALLOCATE(SingleParticleE(0)%value)
		DEALLOCATE(SingleParticleE(1)%value)

		DEALLOCATE(SingleParticleP(0)%value)
		DEALLOCATE(SingleParticleP(1)%value)

		DEALLOCATE(SingleParticleR2(0)%value)
		DEALLOCATE(SingleParticleR2(1)%value)

		DEALLOCATE(EqpCanonical(0)%value)
		DEALLOCATE(EqpCanonical(1)%value)

		DEALLOCATE(an)
		DEALLOCATE(ap)

		DEALLOCATE(DensityRadial)
		DEALLOCATE(RadialWaveFunction)

		RETURN
	END SUBROUTINE DiagonalizationMethod_show_ParticleEnergies

	!-------------------------------------------------------------------------------!
	!										!
	!    This subroutine calculates and prints the s.p. quasi-particle energies	!
	!    in the HFB basis. It also gives an option to store the radial component	!
	!    of the q.p. wave-functions (both the upper part U and lower part V).	!
	!										!
	!-------------------------------------------------------------------------------!

	SUBROUTINE DiagonalizationMethod_show_QuasiParticleEnergies(diagonal)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal

		TYPE (SymD3Tensor) :: ek_tensor

		INTEGER, DIMENSION(:), ALLOCATABLE :: nIndx, pIndx

                INTEGER, DIMENSION(:, :), ALLOCATABLE :: Momentum

                DOUBLE PRECISION, DIMENSION(1:Nmax, 1:Nmax) :: Matrix
                DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: nVV, pVV
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: U, V, QPE, Gamma, Delta
		DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: QPWaveFunctionLower, QPWaveFunctionUpper
		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: PairingField, MeanField

                INTEGER :: NumberOfStates, a, top, Bra, Ket
                INTEGER :: num, la, ja, d, s1, ta, sa, an, ap, jn, jp
                INTEGER :: IndexBasis, IndexKet, ipoint, IndxQP
		INTEGER :: u1, u2, file_error
		INTEGER :: m, n, i, state, Anucl, Nnucl, Znucl

		INTEGER, PARAMETER :: file_unit_1 = 16, file_unit_2 = 17

		DOUBLE PRECISION :: sumvv, factor, x, Rwave, boscil

		CHARACTER :: cn, cp

		CHARACTER, DIMENSION(0:19) :: spectr

		DATA spectr / "s", "p", "d", "f", "g", "h", "i", "j", "k", "l", &
				"m", "n", "o", "p", "q", "r", "s", "t", "u", "v" /

		NumberOfStates = 0
		DO a = 0, 2*Lmax
			NumberOfStates = NumberOfStates + DIM(a)
		END DO

		! Allocate memory for the quasi-particle energies, occupation
                ! probabilities and angular momentum

		ALLOCATE(nVV(NumberOfStates))
		ALLOCATE(pVV(NumberOfStates))

		ALLOCATE(nIndx(NumberOfStates))
		ALLOCATE(pIndx(NumberOfStates))

		ALLOCATE(QPE(0:1, 1:NumberOfStates))
		ALLOCATE(Momentum(0:1, 1:NumberOfStates))

		ALLOCATE(QPWaveFunctionLower(0:1, 0:Npoint, 1:NumberOfStates))
		ALLOCATE(QPWaveFunctionUpper(0:1, 0:Npoint, 1:NumberOfStates))

		! Initialization of the neutron and proton radial density

		IF (.NOT. ALLOCATED(DensityRadial)) ALLOCATE(DensityRadial(0:1, 0:Npoint))
		IF (.NOT. ALLOCATED(DensityPairing)) ALLOCATE(DensityPairing(0:1, 0:Npoint))
		IF (.NOT. ALLOCATED(PairingField)) ALLOCATE(PairingField(0:1, 0:Npoint))
                IF (.NOT. ALLOCATED(DerivDens)) ALLOCATE(DerivDens(0:1, 0:Npoint))
		IF (.NOT. ALLOCATED(MeanField)) ALLOCATE(MeanField(0:1, 0:Npoint))

  		Nnucl = Nucleus_get_N(diagonal%consistency%density%nucleus)
 		Znucl = Nucleus_get_Z(diagonal%consistency%density%nucleus)

                boscil = Nucleus_get_b(diagonal%consistency%density%nucleus)

		Anucl = Nnucl + Znucl
		factor = 1.0D0 - (1.0D0 / Anucl)

		! Initialize all densities to 0
		DensityRadial = 0.0D0
		DensityPairing = 0.0D0
		PairingField = 0.0D0
		DerivDens = 0.0D0
		MeanField = 0.0D0

		CALL SymD3Tensor_new(ek_tensor)

		DO ta = 0, 1

			num = 0
			state = 0

			DO a = 0, 2*Lmax

				la = L(a)
				ja = J(a)
									d = MIN(Nmax, NmaxOfL(la))
				IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - la) / 2) + 1

				! Collecting the q.-p. energies and calculating the occupation factors

				DO s1 = d + 1, d + d

                                	num = num + 1

					QPE(ta, num) = diagonal%QuasiParticleEnergies(ta, a)%value(s1)

					Momentum(ta, num) = a

					IF (ta .EQ. 1) nIndx(num) = num
					IF (ta .EQ. 0) pIndx(num) = num

					sumvv = 0.0D0
					DO sa = 1, d
						sumvv = sumvv + (diagonal%UV(ta, a)%quantum(sa + d, s1) ** 2)
					END DO

					IF ((diagonal%consistency%density%nucleus%is_blocking(ta)) .AND. &
						(a .EQ. diagonal%consistency%density%nucleus%ia(ta)) .AND. &
						(s1 .EQ. (d + diagonal%consistency%density%nucleus%mu0(ta)))) THEN

						DO sa = 1, d
							sumvv = sumvv +((diagonal%UV(ta, a)%quantum(sa, s1) ** 2) &
								      - (diagonal%UV(ta, a)%quantum(sa + d, s1) ** 2)) / (ja + 1.0D0)
						END DO

					END IF

					IF (ta .EQ. 1) nVV(num) = sumvv
					IF (ta .EQ. 0) pVV(num) = sumvv

				END DO

				! First extracting the matrices U and V

				ALLOCATE(U(d, d))
				ALLOCATE(V(d, d))

				DO u1 = 1, d
					DO u2 = 1, d
						U(u1, u2) = diagonal%UV(ta, a)%quantum(u1    , u2 + d)
						V(u1, u2) = diagonal%UV(ta, a)%quantum(u1 + d, u2 + d)
					END DO
				END DO

				ALLOCATE(Delta(d, d))
				ALLOCATE(Gamma(d, d))

				! Define the pairing field
				Delta = diagonal%S%kap%rho(ta, a)%store

				! Define the mean-field (we've got to substract the kinetic energy here to get the
				! mean-field potential)
               			CALL SymD3Tensor_product(ek_tensor, factor, EkField)

				Matrix = SymD3Tensor_matrix(ek_tensor, la)

				Gamma = diagonal%S%rho%rho(ta, a)%store - Matrix(1:d, 1:d)

				! Collecting the wave-functions: upper part (corresponding to the vector U) and lower part
				! (corresponding to the vector V)

				IF (Basis .EQ. 1) THEN

					IF (.NOT. ALLOCATED(WaveFun)) ALLOCATE(WaveFun(0:Npoint,Nsize))

                                        DO ipoint = 1, Npoint
                                                x = RadMesh(ipoint)
                                                DO i = 1, d
						        Bra = IndexVecNL(i, la)
                                                        CALL RadialWaveHO(i, la, x, Rwave, boscil)
							WaveFun(ipoint, Bra) = Rwave*RadMesh(ipoint)
                                                END DO
                                        END DO
                                END IF

				DO ipoint = 1, Npoint

					DO i = 1, d

						Bra = IndexVecNL(i, la)

						DO m = 1, d

							Ket = IndexVecNL(m, la)

							PairingField(ta, ipoint) = PairingField(ta, ipoint) &
						+ (-1)**la * (ja + 1.0D0) * Delta(i, m) * WaveFun(ipoint, Bra)/RadMesh(ipoint) &
										        * WaveFun(ipoint, Ket)/RadMesh(ipoint)

							MeanField(ta, ipoint) = MeanField(ta, ipoint) &
						+ (ja + 1.0D0) * Gamma(i, m) * WaveFun(ipoint, Bra)/RadMesh(ipoint) &
									     * WaveFun(ipoint, Ket)/RadMesh(ipoint)

						END DO

					END DO

				END DO

				DEALLOCATE(Delta)
				DEALLOCATE(Gamma)

				DO i = 1, d

					state = state + 1

					DO ipoint = 1, Npoint

						QPWaveFunctionUpper(ta, ipoint, state) = 0.0D0
						QPWaveFunctionLower(ta, ipoint, state) = 0.0D0

						DO m = 1, d

							IndexBasis = IndexVecNL(m, la)

							QPWaveFunctionUpper(ta, ipoint, state) = &
							QPWaveFunctionUpper(ta, ipoint, state) &
							+ U(m, i) * WaveFun(ipoint,IndexBasis)/RadMesh(ipoint)

							QPWaveFunctionLower(ta, ipoint, state) = &
							QPWaveFunctionLower(ta, ipoint, state) &
							+ V(m, i) * WaveFun(ipoint,IndexBasis)/RadMesh(ipoint)

							! Calculating here the single-particle density. We have:
							!
							!		\rho_{mn} =  { V (V+) }_{mn}
							!
							! and in coordinate representation:
							!
							!          rho(r) = \sum_{mn,l,j} (2j+1) \rho_{mn} \phi_{ml}(r)\phi_{nl}(r)
							!
							! or:
							!     rho(r) = \sum_{i, mn,l,j} (2j+1) V_{ni}V_{mi} \phi_{ml}(r)\phi_{nl}(r)
							!

							DO n = 1, d

								IndexKet = IndexVecNL(n, la)

								DensityRadial(ta, ipoint) = DensityRadial(ta, ipoint) &
											+ (ja + 1.0D0) * V(m, i) * V(n, i) &
											* WaveFun(ipoint,IndexBasis) &
											* WaveFun(ipoint,IndexKet)

                                                                DensityPairing(ta, ipoint) = DensityPairing(ta, ipoint) &
                                                                                    - (-1)**la * (ja + 1.0D0) * V(m, i) * U(n, i) &
											* WaveFun(ipoint,IndexBasis) &
											* WaveFun(ipoint,IndexKet)

							END DO

						END DO

					END DO

				END DO

				DEALLOCATE(U)
				DEALLOCATE(V)

			END DO ! end of loop over a

			! The basis wave-functions are in fact y(r) = R(r)/r, so we need to correct for this.
			! The 4*pi comes from the integration over the angles.

			DO ipoint=1, Npoint
				DensityRadial(ta, ipoint) = DensityRadial(ta, ipoint) / (4.0D0*PI*RadMesh(ipoint)**2)
				DensityPairing(ta, ipoint) = DensityPairing(ta, ipoint) / (4.0D0*PI*RadMesh(ipoint)**2)
			END DO

                        ! Extrapolate to find rho(r=0)

         		DensityRadial(ta, 0) = 3.0D0*(DensityRadial(ta, 1) - DensityRadial(ta, 2)) + DensityRadial(ta, 3)
         		DensityPairing(ta, 0) = 3.0D0*(DensityPairing(ta, 1) - DensityPairing(ta, 2)) + DensityPairing(ta, 3)
         		PairingField(ta, 0) = 3.0D0*(PairingField(ta, 1) - PairingField(ta, 2)) + PairingField(ta, 3)
         		MeanField(ta, 0) = 3.0D0*(MeanField(ta, 1) - MeanField(ta, 2)) + MeanField(ta, 3)

                        DO ipoint=1, Npoint-1
                                DerivDens(ta, ipoint) = 0.5d0*(DensityRadial(ta, ipoint+1) - DensityRadial(ta, ipoint-1)) &
                                                              /(RadMesh(ipoint+1) - RadMesh(ipoint))
                        END DO

                        DerivDens(ta, 0) = 3.0D0*(DerivDens(ta, 1) - DerivDens(ta, 2)) + DerivDens(ta, 3)
                        DerivDens(ta, Npoint) = 3.0D0*(DerivDens(ta, Npoint-1) - DerivDens(ta, Npoint-2)) + DerivDens(ta, Npoint-3)

		END DO ! end of loop over ta

		! Ordering of the quasi-particle occupation numbers (by decreasing order!!)

		CALL indexx_real8(NumberOfStates, 1.0D0-nVV, nIndx)
		CALL indexx_real8(NumberOfStates, 1.0D0-pVV, pIndx)

		WRITE(*,'()')
		WRITE(*,'(5X,"QUASI-PARTICLE-PARTICLE ENERGIES")')
		WRITE(*,'(5X,"================================")')

		WRITE(*,'(A14,A10,A12,A5,A10)') "EQP NEUT", "V2", "", "V2", "EQP PROT"

		top = MIN(100, NumberOfStates)

		DO num = 1, top

			an = Momentum(1, nIndx(num))
			ap = Momentum(0, pIndx(num))

			cn = spectr(L(an))
			cp = spectr(L(ap))

			jn = J(an)
			jp = J(ap)

                        WRITE(*,'(i4,")",i5,F10.4,3X,F10.7,3X,A1,I2,"/2",A5,I2,"/2",3x,F10.7,F10.4)') &
				num, nIndx(num), QPE(1, nIndx(num)), nVV(nIndx(num)), cn, jn, cp, jp, pVV(pIndx(num)), QPE(0, pIndx(num))

			! Getting v2 of QPQP of lowest energy

		END DO

		! Quasi-particle energies (proton and neutron)

		OPEN(file_unit_1, FILE='data/HF_qp_n.dat', ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			WRITE(*,'("Impossible to open the file data/HF_qp_n.dat")')
			STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
		END IF

		DO num = 1, top
			an = Momentum(1, nIndx(num))
			cn = spectr(L(an))
			jn = J(an)
			WRITE(file_unit_1,'(I3,2f25.15,1X,A1,i5)') num, QPE(1, nIndx(num)), nVV(nIndx(num)), cn, jn
		END DO

		CLOSE(file_unit_1)

		OPEN(file_unit_1, FILE='data/HF_qp_p.dat', ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			WRITE(*,'("Impossible to open the file data/HF_qp_p.dat")')
			STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
		END IF

		DO num = 1, top
			ap = Momentum(1, pIndx(num))
			cp = spectr(L(ap))
			jp = J(ap)
			WRITE(file_unit_1,'(I3,2f25.15,1X,A1,i5)') num, QPE(0, pIndx(num)), pVV(pIndx(num)), cp, jp
		END DO

		CLOSE(file_unit_1)

                ! Total density (proton and neutron)

                OPEN(file_unit_1, FILE='data/DensityQP.dat', ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			WRITE(*,'("Impossible to open the file data/DensityQP.dat")')
			STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
		END IF

		DO ipoint = 0, Npoint
			WRITE(file_unit_1,'(3f20.16)') RadMesh(ipoint),DensityRadial(0, ipoint),DensityRadial(1, ipoint)
		END DO

		CLOSE(file_unit_1)

		! Total pairing density (proton and neutron)

		OPEN(file_unit_1, FILE='data/DensityPair.dat', ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			WRITE(*,'("Impossible to open the file data/DensityPair.dat")')
			STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
		END IF

		DO ipoint = 0, Npoint
			WRITE(file_unit_1,'(3f20.16)') RadMesh(ipoint),DensityPairing(0, ipoint),DensityPairing(1, ipoint)
		END DO

		CLOSE(file_unit_1)

        	! Logarithmic derivative of the density (proton and neutron)

		OPEN(file_unit_1, FILE='data/DensityLogDeriv.dat', ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			WRITE(*,'("Impossible to open the file data/DensityLogDeriv.dat")')
			STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
		END IF

		DO ipoint = 0, Npoint
			WRITE(file_unit_1,'(3f25.16)') RadMesh(ipoint),DerivDens(0, ipoint)/(DensityRadial(0, ipoint) + 1.e-20), &
                                                                       DerivDens(1, ipoint)/(DensityRadial(1, ipoint) + 1.e-20)
		END DO

		CLOSE(file_unit_1)

		! Pairing and mean-field potentials. They are defined as:
		!
		!   Delta(r) = \sum delta_{ac} phi_{a}(r)*phi_{c}(r)
		!   Gamma(r) = \sum gamma_{ac} phi_{a}(r)*phi_{c}(r)

		OPEN(file_unit_1, FILE='data/PairingField.dat', ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			WRITE(*,'("Impossible to open the file data/PairingField.dat")')
			STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
		END IF

		DO ipoint = 0, Npoint
			WRITE(file_unit_1,'(5f20.14)') RadMesh(ipoint), &
				PairingField(0, ipoint)/(2.0D0*PI), PairingField(1, ipoint)/(2.0D0*PI),&
				MeanField(0, ipoint)/(2.0D0*PI), MeanField(1, ipoint)/(2.0D0*PI)
		END DO

		CLOSE(file_unit_1)

		! Storing quasi-particle wave-functions (proton and neutrons).

		OPEN(file_unit_1, FILE='data/WavesQP.dat', ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			WRITE(*,'("Impossible to open the file data/WavesQP.dat")')
			STOP "In DiagonalizationMethod_show_ParticleEnergies - Impossible to open file"
		END IF

		IndxQP = NumberOfStates - IndexWave + 1

		write(*,'("IndexWave = ",i4," IndxQP = ",i4," nIndx = ",i4)') IndexWave, IndxQP, nIndx(IndxQP)

		QPWaveFunctionLower(0, 0, pIndx(IndxQP)) = 3.0D0*(QPWaveFunctionLower(0, 1, pIndx(IndxQP)) &
							 	- QPWaveFunctionLower(0, 2, pIndx(IndxQP))) &
					 			+ QPWaveFunctionLower(0, 3, pIndx(IndxQP))

		QPWaveFunctionUpper(0, 0, pIndx(IndxQP)) = 3.0D0*(QPWaveFunctionUpper(0, 1, pIndx(IndxQP)) &
							  	- QPWaveFunctionUpper(0, 2, pIndx(IndxQP))) &
							  	+ QPWaveFunctionUpper(0, 3, pIndx(IndxQP))

		QPWaveFunctionLower(1, 0, nIndx(IndxQP)) = 3.0D0*(QPWaveFunctionLower(1, 1, nIndx(IndxQP)) &
							 	- QPWaveFunctionLower(1, 2, nIndx(IndxQP))) &
					 			+ QPWaveFunctionLower(1, 3, nIndx(IndxQP))

		QPWaveFunctionUpper(1, 0, nIndx(IndxQP)) = 3.0D0*(QPWaveFunctionUpper(1, 1, nIndx(IndxQP)) &
							  	- QPWaveFunctionUpper(1, 2, nIndx(IndxQP))) &
							  	+ QPWaveFunctionUpper(1, 3, nIndx(IndxQP))

		DO ipoint = 0, Npoint
			WRITE(file_unit_1,'(5f20.16)') RadMesh(ipoint), QPWaveFunctionUpper(0, ipoint, pIndx(IndxQP)), &
								QPWaveFunctionLower(0, ipoint, pIndx(IndxQP)), &
								QPWaveFunctionUpper(1, ipoint, nIndx(IndxQP)), &
								QPWaveFunctionLower(1, ipoint, nIndx(IndxQP))
		END DO

		CLOSE(file_unit_1)

		DEALLOCATE(nIndx)
		DEALLOCATE(pIndx)

		DEALLOCATE(QPE)

		DEALLOCATE(nVV)
		DEALLOCATE(pVV)

		DEALLOCATE(Momentum)

		DEALLOCATE(QPWaveFunctionUpper)
		DEALLOCATE(QPWaveFunctionLower)

		RETURN

	END SUBROUTINE DiagonalizationMethod_show_QuasiParticleEnergies

	SUBROUTINE DiagonalizationMethod_del(diagonal)
		TYPE (DiagonalizationMethod), INTENT(INOUT) :: diagonal

		INTEGER :: ta, a, max_a

		max_a = 2*Lmax
		DO ta = 0, 1
			DO a = 0, max_a
				DEALLOCATE(diagonal%QuasiParticleEnergies(ta, a)%value)
				DEALLOCATE(diagonal%UV(ta, a)%quantum)
			END DO
		END DO
		DEALLOCATE(diagonal%QuasiParticleEnergies)
		DEALLOCATE(diagonal%UV)

		CALL SelfConsistencyMethod_del(diagonal%consistency)
		CALL R1R1Function_del(diagonal%func)

		CALL SymGenDensity_del(diagonal%S)
		CALL SymGenDensity_del(diagonal%iterated)
		RETURN
	END SUBROUTINE DiagonalizationMethod_del

END MODULE diagmeth
 MODULE eigenval

 CONTAINS

	! A = Matriz simetrica
	! D = Matriz diagonal
	! Z = Matriz
	SUBROUTINE EigenValues(NM, N, A, D, Z)
		INTEGER, INTENT(IN) :: NM, N
		DOUBLE PRECISION, DIMENSION(:, :), INTENT(IN) :: A
		DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: D
		DOUBLE PRECISION, DIMENSION(:, :), INTENT(INOUT) :: Z

		! E = Matriz diagonal
		DOUBLE PRECISION, DIMENSION(NM, N) :: E
		INTEGER IERR

		CALL TRED2(NM, N, A, D, E, Z)
		
		CALL TQL2(NM, N, D, E, Z, IERR)
		
		IF (IERR .NE. 0) STOP "ERROR: EigenValues"
		RETURN
	END SUBROUTINE EigenValues

      SUBROUTINE TQL2(NM,N,D,E,Z,IER)
!-------------------------------------------------------------------------
!     QL METHOD TO DETERMINE THE EIGENVALUES AND EIGENVECTORS OF:
!
!       1)  A SYMMETRIC TRIDIAGONAL MATRIX.
!       2)  A FULL SYMMETRIC MATRIX AFTER A PREVIOUS CALL TO TRED2.
!
!     CALLING MODE:
!               CALL TQL2(NM,N,D,E,Z,IER)
!     INPUTSS:
!     NM  (I4)  1ST DIMENSION OF MATRICES A AND Z IN CALLING PROGRAM
!     N   (I4)  SIZE OF Z
!     D  (R*8)  MAIN DIAGONAL (N) OF THE TRIDIAGONAL MATRIX
!     E  (R*8)  SUB-DIAGONAL (N) OF THE TRIDIAGONAL MATRIX
!     Z  (R*8)  TABLE (NM,N) STORING THE UNITY MATRIX IF THE TRIDIAGONAL
!               MATRIX IS DEFINED BY D AND E, CASE #1.
!               FOR CASE #2, IT CONTAINS THE ELEMENTS OF THE TRANSFORMATION
!               MATRIX AFTER A CALL TO TRED2.
!     OUTPUTS:
!     D  (R*8)  EIGENVALUES
!     Z  (R*8)  EIGENVECTORS
!     IER (I4)  ERROR CODE = 0,  CONVERGENCE OK.
!                          = L,  NO CONVERGENCE FOR THE Lth EIGENVALUE
!
!     REFERENCE:
!     J.H.WILKINSON,-C.REINSCH,R.S.MARTIN
!     HANDBOOK FOR AUTOMATIC COMPUTATION, VOL.2, LINEAR ALGEBRA
!     SPRINGER-VERLAG 1971.
!-------------------------------------------------------------------------
      INTEGER I,J,K,L,M,N,NM,JM,IER
      !REAL *8 D(N),E(N),Z(NM,N),B,C,F,G,H,P,R,S,EPS,EPS1
      DOUBLE PRECISION :: D(N),E(N),Z(NM,N),B,C,F,G,H,P,R,S,EPS,EPS1
      DATA EPS /0.D0/,JM /30/
      IER = 0
      IF (N.EQ.1) GO TO 38
!
!     MACHINE EPSILON
!
      IF (EPS.NE.0.D0) GO TO 12
      EPS = 1.D0
   10 EPS = EPS/2.D0
      EPS1 = 1.D0+EPS
      IF (EPS1.GT.1.D0) GO TO 10
!
   12 DO 14 I = 2,N
   14 E(I-1) = E(I)
      E(N) = 0.D0
      F = 0.D0
      B = 0.D0
!
      DO 28 L = 1,N
      J = 0
      H = EPS*(ABS(D(L))+ABS(E(L)))
      IF (B.LT.H) B = H
!
!     SEEK SMALLEST ELEMENT OF SUBDIAGONAL
!
      DO 16 M = L,N
      IF (ABS(E(M)).LE.B) GO TO 18
   16 CONTINUE
   18 IF (M.EQ.L) GO TO 26

!     START ITERATION

   20 IF (J.EQ.JM) GO TO 36
      J = J+1

!     SHIFT

      G = D(L)
      P = (D(L+1)-G)/(2.D0*E(L))
      R = SQRT(P*P+1.D0)
      D(L) = E(L)/(P+SIGN(R,P))
      H = G-D(L)
      DO 22 I = L+1,N
   22 D(I) = D(I)-H
      F = F+H

!     QL TRANSFORMATION

      P = D(M)
      C = 1.D0
      S = 0.D0
      DO 24 I = M-1,L,-1
      G = C*E(I)
      H = C*P
      IF (ABS(P).GE.ABS(E(I))) THEN
      C = E(I)/P
      R = SQRT(C*C+1.D0)
      E(I+1) = S*P*R
      S = C/R
      C = 1.D0/R
      ELSE
      C = P/E(I)
      R = SQRT(C*C+1.D0)
      E(I+1) = S*E(I)*R
      S = 1.D0/R
      C = C*S
      ENDIF
      P = C*D(I)-S*G
      D(I+1) = H+S*(C*G+S*D(I))

!     ELEMENTS OF EIGENVECTORS

      DO 24 K = 1,N
      H = Z(K,I+1)
      Z(K,I+1) = S*Z(K,I)+C*H
      Z(K,I) = Z(K,I)*C-S*H
   24 CONTINUE
      E(L) = S*P
      D(L) = C*P
      IF (ABS(E(L)).GT.B) GO TO 20

!     CONVERGENCE

   26 D(L) = D(L)+F
   28 CONTINUE

!     SORT EIGENVALUES AND EIGENVECTORS
!     IN ASVENDING ORDER

      DO 34 L = 2,N
      I = L-1
      K = I
      P = D(I)
      DO 30 J = L,N
      IF (D(J).GE.P) GO TO 30
      K = J
      P = D(J)
   30 CONTINUE
      IF (K.EQ.I) GO TO 34
      D(K) = D(I)
      D(I) = P
      DO 32 J = 1,N
      P = Z(J,I)
      Z(J,I) = Z(J,K)
   32 Z(J,K) = P
   34 CONTINUE
      GO TO 38

!     NO CONVERGENCE

   36 IER = L
   38 RETURN
      END SUBROUTINE TQL2

      SUBROUTINE TRED2(NM,N,A,D,E,Z)
!---------------------------------------------------------------------------
!     TRIDIAGONALIZATION OF A SYMMETRIC MATRIX BY ORTHOGONAL TRANSFORMATIONS
!     (ALGORITHM OF HOUSEHOLDER)
!     CALLING MODE:
!               CALL TRED2(NM,N,A,D,E,Z)
!     INPUTS:
!     NM  (I4)  1ST DIMENSION OF MATRICES A AND Z IN CALLING PROGRAM
!     N   (I4)  SIZE OF A
!     A  (R*8)  TABLE(NM,N) STORING THE COEFFICIENTS OF SYMMETRIC A MATRIX
!               (LOWER HALF), A IS NOT DESTROYED DURING THE PROCESS
!               IF Z MATRIX HAS NOT THE SAME ADDRESS.
!     OUTPUTS:
!     D  (R*8)  MAIN DIAGONAL (N) OF REDUCED TRIDIAGONAL MATRIX
!     E  (R*8)  SUB-DIAGONAL (N) OF REDUCED TRIDIAGONAL MATRIX
!     Z  (R*8)  TABLE (NM,N) STORING THE ELEMENTS OF THE ORTHOGONAL 
!               TRANSFORMATION MATRIX.
!     REFERENCE:
!     J.H.WILKINSON,-C.REINSCH,R.S.MARTIN
!     HANDBOOK FOR AUTOMATIC COMPUTATION, VOL.2, LINEAR ALGEBRA
!     SPRINGER-VERLAG 1971.
!-----------------------------------------------------------------------
      INTEGER I,J,K,L,N,NM
      !REAL *8 A(NM,N),D(N),E(N),Z(NM,N),F,G,H,HH,SCALE
      DOUBLE PRECISION :: A(NM,N),D(N),E(N),Z(NM,N),F,G,H,HH,SCALE

!     LOWER HALF OF A PUT INTO Z

      DO 10 I = 1,N
      DO 10 J = 1,I
   10 Z(I,J) = A(I,J)
      IF (N.EQ.1) GO TO 32

!     N-2 STAGE OF TRANSFORMATION

      DO 30 I = N,2,-1
      L = I-1
      H = 0.

!     CONDITIONNING BY NORM OF A

      SCALE = 0.
      IF (L.LT.2) GO TO 14
      DO 12 K = 1,L
   12 SCALE = SCALE+ABS(Z(I,K))
      IF (SCALE.NE.0.) GO TO 16

   14 E(I) = Z(I,L)
      GO TO 28

   16 DO 18 K = 1,L
      Z(I,K) = Z(I,K)/SCALE
      H = H+Z(I,K)*Z(I,K)
   18 CONTINUE

      F = Z(I,L)
      G = -SIGN(SQRT(H),F)
      E(I) = SCALE*G
      H = H-F*G
      Z(I,L) = F-G
      F = 0.
      DO 24 J = 1,L
      Z(J,I) = Z(I,J)/H
      G = 0.

!     ELEMENT OF A*U
      DO 20 K = 1,J
   20 G = G+Z(J,K)*Z(I,K)
      IF (L.GE.J+1) THEN
      DO 22 K = J+1,L
   22 G = G+Z(K,J)*Z(I,K)

!     ELEMENT OF P = A*U/H

      END IF
      E(J) = G/H
      F = F+E(J)*Z(I,J)
   24 CONTINUE

!     ELEMENT OF K

      HH = F/(H+H)

!     REDUCED FORM OF A

      DO 26 J = 1,L
      F = Z(I,J)
      G = E(J)-HH*F
      E(J) = G
      DO 26 K = 1,J
      Z(J,K) = Z(J,K)-F*E(K)-G*Z(I,K)
   26 CONTINUE
!
   28 D(I) = H
   30 CONTINUE

!     END OF TRANSFORMATION

   32 D(1) = 0.
      E(1) = 0.

!     ACCUMULATE TRANSFORMATION MATRICES IN Z

      DO 40 I = 1,N
      L = I-1
      IF (D(I).NE.0.) THEN
      DO 36 J = 1,L
      G = 0.
      DO 34 K = 1,L
   34 G = G+Z(I,K)*Z(K,J)
      DO 36 K = 1,L
      Z(K,J) = Z(K,J)-G*Z(K,I)
   36 CONTINUE
      END IF
      D(I) = Z(I,I)
      Z(I,I) = 1.
      IF (L.LT.1) GO TO 40
      DO 38 J = 1,L
      Z(I,J) = 0.
      Z(J,I) = 0.
   38 CONTINUE
   40 CONTINUE

      RETURN
      END SUBROUTINE TRED2

END MODULE eigenval
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
		ALLOCATE(gauss%w(0:n)) !TODO
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
		DATA endpts /0.0D0, 0.0D0/
		DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: wksp

		ALLOCATE(wksp(1:n))
		kind2 = 6
		bet  = 0.0D0
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
				t(n) = gbslve(endpts(1), n, t, b) * (b(n - 1) ** 2) + DBLE(endpts(1))
			END IF
		END IF

		w(1) = 1.0D0
		DO i = 2, n
			w(i) = 0.0D0
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
			muzero = 2.0D0
			DO i = 1, n - 1
				a(i) = 0.0D0
				b(i) = i / SQRT(i * 4.0D0 * i - 1.0D0)
			END DO
			a(n) = 0.0D0
		CASE (2)
			muzero = PI
			DO i = 1, n - 1
				a(i) = 0.0D0
				b(i) = 0.5D0
			END DO
			b(1) = SQRT(0.5D0)
			a(n) = 0.0D0
		CASE (3)
			muzero = PI / 2.0D0
			DO i = 1, n - 1
				a(i) = 0.0D0
				b(i) = 0.5D0
			END DO
			a(n) = 0.0D0
		CASE (4)
			muzero = SQRT(PI)
			DO i = 1, n - 1
				a(i) = 0.0D0
				b(i) = SQRT(0.5D0 * i)
			END DO
			a(n) = 0.0
		CASE (5)
			ab = alf + bet
			abi = ab + 2.0D0
			muzero = (2.0D0 ** (ab + 1.0D0)) * dgamma(alf + 1.0D0) * dgamma(bet + 1.0D0) / dgamma(ab + 2.0D0)
			a(1) = (bet - alf) / abi
			b(1) = SQRT((alf + 1.0D0) * 4.0D0 * (bet + 1.0D0) / ((abi + 1.0D0) * abi * abi))

			a2b2 = bet * bet - alf * alf
			DO i = 2, n - 1
				abi = i * 2.0D0 + ab
				a(i) = a2b2 / ((abi - 2.) * abi)
				b(i) = SQRT(i * 4.0D0 * (i + alf) * (i + bet) * (i + ab) / ((abi * abi - 1.0D0) * abi * abi))
			END DO
			abi = n * 2.0D0 + ab
			a(n) = a2b2 / ((abi - 2.0D0) * abi)
		CASE (6)
			muzero = dgamma(alf + 1.0D0)
			DO i = 1, n - 1
				a(i) = 2.0D0 * i - 1.0D0 + alf
				b(i) = SQRT((alf + i) * i)
			END DO
			a(n) = n * 2.0D0 - 1.0D0 + alf
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
			alpha = a(i) - shift - (b(i - 1) ** 2.0D0) / alpha
		END DO
		gbslve = 1.0D0 / alpha
		RETURN
	END FUNCTION gbslve

	SUBROUTINE gbtql2(n, d, e, z, ierr)
		INTEGER, INTENT(IN) :: n
		DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: e
		DOUBLE PRECISION, DIMENSION(0:), INTENT(INOUT) :: d, z
		INTEGER, INTENT(INOUT) :: ierr

		DOUBLE PRECISION b, c, f, g, p, r, s
		INTEGER i, ii, j, k, l, m, m2, flag

		ierr = 0
		IF (n .EQ. 1) RETURN
		e(n) = 0.0D0
		DO l = 1, n
			j = 0
bucle:			DO
                                flag = 0
				DO m2 = l, n-1
                                   IF (ABS(e(m2)) .LE. EPS * (ABS(d(m2)) + ABS(d(m2 + 1)))) THEN
                                           m = m2
                                           flag=1
                                                EXIT
                                   END IF
				END DO
                                IF (flag .EQ. 0) m=n
				p = d(l)
				IF (m .EQ. l) EXIT bucle

				IF (j .EQ. 30) THEN
					ierr = 1
					RETURN
				END IF
				j = j + 1
				g = (d(l + 1) - p) / (e(l) * 2.0D0)
				r = SQRT(g * g + 1.0D0)
				g = d(m) - p + e(l) / (g + SIGN(r, g))
				s = 1.0D0
				c = 1.0D0
				p = 0.0D0
				DO ii = 1, m - l
					i = m - ii
					f = s * e(i)
					b = c * e(i)
					IF (ABS(f) .LT. ABS(g)) THEN
						s = f / g
						r = SQRT(s * s + 1.0D0)
						e(i + 1) = g * r
						c = 1.0D0 / r
						s = s * c
					ELSE
						c = g / f
						r = SQRT(c * c + 1.0D0)
						e(i + 1) = f * r
						s = 1.0D0 / r
						c = c * s
					END IF
					g = d(i + 1) - p
					r = (d(i) - g) * s + c * 2.0D0 * b
					p = s * r
					d(i + 1) = g + p
					g = c * r - b
					f = z(i + 1)
					z(i + 1) = s * z(i) + c * f
					z(i    ) = c * z(i) - s * f
				END DO
				d(l) = d(l) - p
				e(l) = g
				e(m) = 0.0D0
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
		DATA a / 1.0000000000000000D0, &
			 0.4227843350984678D0, &
			 0.4118403304263672D0, &
			 0.0815769192502609D0, &
			 0.0742490106800904D0, &
			-2.669810333484000D-4, &
			 0.0111540360240344D0, &
			-0.0028525821446197D0, &
			 0.0021036287024598D0, &
			-9.184843690991000D-4, &
 			 4.874227944768000D-4, &
			-2.347204018919000D-4, &
			 1.115339519666000D-4, &
 			-4.787479838340000D-5, &
			 1.751027271790000D-5, &
			-4.920375090400000D-6, &
			 9.199156407000000D-7, &
			-8.399404960000000D-8/

		INTEGER k
		DOUBLE PRECISION p, t

		IF (z .LE. 1.0D0) THEN
			t = z
		ELSE IF (z .LE. 2.0D0) THEN
			t = z - 1.0D0
		ELSE
			t = z - 2.0D0
		END IF

		p = a(17)
		DO k = 1, 17
			p = t * p + a(17 - k)
		END DO

		IF (z .GT. 2.0D0) THEN
			dgamma = p
		ELSE IF (z .GT. 1.0D0) THEN
			dgamma = p / z
		ELSE
			dgamma = p / (z * (z + 1.0D0))
		END IF
		RETURN
	END FUNCTION dgamma

END MODULE gauss
 MODULE global

	USE input
	USE lgfactor
	USE symd3t
	USE symtalm
	USE gauss
	USE angmom
	USE bessik

	IMPLICIT NONE

	! PUBLIC VARIABLES
	LOGICAL, PUBLIC :: test_regularization = .False.
	INTEGER, PUBLIC :: Ngaussian = 1
	DOUBLE PRECISION, PUBLIC :: delta_a = 1.D-10
	DOUBLE PRECISION, ALLOCATABLE :: range_gaussian(:)

	! We define the flags that will identify each isospin throughout the program
	INTEGER, PARAMETER, PUBLIC :: PROTON  = 0
	INTEGER, PARAMETER, PUBLIC :: NEUTRON = 1

	! The list of spherical magic numbers
	INTEGER, DIMENSION(0:8), PUBLIC :: MagicNumber
	DATA MagicNumber/ 2, 8, 20, 28, 50, 82, 126, 184, 256/

	! m is equal to mc2/(hbar*c)2

	DOUBLE PRECISION, DIMENSION(0:1), PUBLIC :: m
	!DATA m / 0.02411186800036718652D0, 0.02411186800036718652D0 /
	DATA m / 0.024111868D0, 0.024111868D0 /

	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, PUBLIC :: sq, sq2
	DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, PUBLIC :: BesselTable

	! Parameters of the Gogny force
	DOUBLE PRECISION, DIMENSION(0:1, 4), PUBLIC :: Gogny_W
	DOUBLE PRECISION, DIMENSION(0:1, 4), PUBLIC :: Gogny_B
	DOUBLE PRECISION, DIMENSION(0:1, 4), PUBLIC :: Gogny_H
	DOUBLE PRECISION, DIMENSION(0:1, 4), PUBLIC :: Gogny_M
	DOUBLE PRECISION, PARAMETER, PUBLIC :: x0 = 1.0D0

	DATA Gogny_W / -1720.30D0,  103.64D0, -402.40D0, -21.30D0, -402.40D0, -21.30D0, -2047.61D0,  293.02D0 /
	DATA Gogny_B /  1300.00D0, -163.48D0, -100.00D0, -11.77D0, -100.00D0, -11.77D0,  1700.00D0, -300.78D0 /
	DATA Gogny_H / -1813.53D0,  162.81D0, -496.20D0,  37.27D0, -496.20D0,  37.27D0, -2414.93D0,  414.59D0 /
	DATA Gogny_M /  1397.60D0, -223.93D0,  -23.56D0, -68.81D0,  -23.56D0, -68.81D0,  1519.35D0, -316.84D0 /

	DOUBLE PRECISION, DIMENSION(4), PUBLIC :: Gogny_W0
	DATA Gogny_W0 / 130.0D0, 115.0D0, 130.0D0, 115.0D0 / ! D1S, D1, D1prime, D1N

	DOUBLE PRECISION, DIMENSION(4), PUBLIC :: Gogny_t0
	DATA Gogny_t0 / 1390.6D0, 1350.0D0, 1350.0D0, 1609.46D0 / ! D1S, D1, D1prime, D1N

	! One-body matrix elements of the kinetic energy and Gauss-Laguerre quadrature points
	TYPE (SymD3Tensor), PUBLIC :: EkField, R2Field
	TYPE (GaussLaguerreQuadrature), PUBLIC :: GaussLQ

	! PRIVATE VARIABLES
	INTEGER :: A
	INTEGER :: NLag

 CONTAINS

	SUBROUTINE Global_new

		INTEGER :: max_1, max_2

		! Reading the input parameters
		CALL Input_read

		! Initialize a default mesh
		CALL InitializeMesh()

        	! Reading the basis form WSCOOR
                CALL ReadBasis

		! Printing the parameters of the force
		WRITE(*,'()')
		WRITE(*,'(5X,"PARAMETERS OF THE GOGNY FORCE")')
		WRITE(*,'(5X,"=============================",/)')
		WRITE(*,'("Protons : W = ",F8.2," MeV")') Gogny_W(1, Gogny)
		WRITE(*,'("        : B = ",F8.2," MeV")') Gogny_B(1, Gogny)
		WRITE(*,'("        : H = ",F8.2," MeV")') Gogny_H(1, Gogny)
		WRITE(*,'("        : M = ",F8.2," MeV",/)') Gogny_M(1, Gogny)
		WRITE(*,'("Neutrons: W = ",F8.2," MeV")') Gogny_W(0, Gogny)
		WRITE(*,'("        : B = ",F8.2," MeV")') Gogny_B(0, Gogny)
		WRITE(*,'("        : H = ",F8.2," MeV")') Gogny_H(0, Gogny)
		WRITE(*,'("        : M = ",F8.2," MeV",/)') Gogny_M(0, Gogny)
		WRITE(*,'("WLS = ",F8.2," MeV.fm5")') Gogny_W0(Gogny)
		WRITE(*,'("t_0 = ",F8.2," MeV.fm4",/)') Gogny_t0(Gogny)

		NLag = MIN(N_0 * 32, 156) ! NLag = UMIN(N_0 << 5, 156)

		max_1 = 100
		max_2 = 2000

		IF (regularized_Gaussian) THEN
			IF (test_regularization) THEN
			        Ngaussian=3; ALLOCATE(range_gaussian(0:Ngaussian-1))
				range_gaussian(0) = range1 + delta_a
				range_gaussian(1) = range1 - delta_a
				range_gaussian(2) = range1
			ELSE
			        Ngaussian=1; ALLOCATE(range_gaussian(0:Ngaussian-1))
				range_gaussian(0) = range1
			END IF
		ELSE
			Ngaussian=2; ALLOCATE(range_gaussian(0:Ngaussian-1))
			range_gaussian(0) = range1
			range_gaussian(1) = range2
		END IF

!TODO El autor lo implementa, pero no se utiliza en ninguna parte del código
!		CALL LogFactorials_new(max_2)
!		CALL LogSemiFactorials_new(max_2)

                ! Defined in module "lgfactor.f90"
		CALL GammaFunction_new(max_1)

                ! Defined in module "lgfactor.f90"
		CALL DDLogFactorials_new(max_2)
		CALL DDLogSemiFactorials_new(max_2)

                ! Defined here...
		CALL SquareRoot_new(max_1)
		CALL SemiSquareRoot_new(max_1)

                ! Calculate one-body kinetic energy and r.m.s. radius
		CALL Global_start

                ! Defined in module "symtalm.f90" (only if analytical HO basis is used)
		IF (Basis .EQ. 1) THEN

			CALL SymCoefficientB_new
			CALL SymKumar_new

			CALL GaussLaguerreQuadrature_new(GaussLQ, NLag, DBLE(0.5))

		END IF

                ! Defined in module "angmom.f90"
		CALL ThreeJSymbols_new

		RETURN

	CONTAINS

                ! Subroutine initializing a default box with the corresponding mesh

		SUBROUTINE InitializeMesh()
			INTEGER :: i
			DOUBLE PRECISION :: Rbox

			Npoint = 201
			Rbox = 20.0D0
			MeshStep = Rbox / DBLE(Npoint - 1)

			ALLOCATE(RadMesh(0:Npoint))

			DO i = 0, Npoint
				RadMesh(i) = DBLE(i)*MeshStep
			END DO

			RETURN
		END SUBROUTINE InitializeMesh

                ! Subroutine storing in a vector of size imax the square root of
		! all integers from 1 to imax: sqrt(i)

		SUBROUTINE SquareRoot_new(imax)
			INTEGER, INTENT(IN) :: imax

			INTEGER i

			ALLOCATE(sq(0:imax))
			sq(0) = 0.0D0
			DO i = 1, imax
				sq(i) = SQRT(DBLE(i))
			END DO
			RETURN
		END SUBROUTINE SquareRoot_new

                ! Subroutine storing in a vector of size imax the square root of
		! all integers from 1 to imax, plus one half: sqrt(i+0.5)

		SUBROUTINE SemiSquareRoot_new(imax)
			INTEGER, INTENT(IN) :: imax

			INTEGER i

			ALLOCATE(sq2(0:imax))
			sq2(0) = 0.5d0*SQRT(2.0D0)
			DO i = 1, imax
				sq2(i) = SQRT(i + 0.5d0)
			END DO
			RETURN
		END SUBROUTINE SemiSquareRoot_new

	END SUBROUTINE Global_new

        ! Subroutine storing in a vector of size imax the square root of
	! all integers from 1 to imax: sqrt(i)

	SUBROUTINE BesselTabularize(Norder)
		INTEGER, INTENT(IN) :: Norder

		INTEGER :: i_r1, i_r2, k, igauss

#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: xarg, BesselFunc
		REAL(KIND = 16) :: rk, rip, rkp, pi, Order, r1, r2
#else
		DOUBLE PRECISION :: xarg, BesselFunc
		DOUBLE PRECISION :: rk, rip, rkp, pi, Order, r1, r2
#endif
		DOUBLE PRECISION :: range

		ALLOCATE(BesselTable(1:Npoint,1:Npoint,0:Norder,0:Ngaussian-1))

		pi = 4.0D0*ATAN(1.0D0)

		WRITE(*,'("Tabularization of Bessel Functions - k_max = ",i4)') Norder

                DO igauss = 0, Ngaussian-1

                        range = range_gaussian(igauss)

                        DO k = 0, Norder

			        WRITE(*,'("Order k = ",i4," Range = ",f20.14)') k, range

!$OMP Parallel Default(None) &
!$OMP& SHARED(igauss,k,Npoint,BesselTable,range,RadMesh,pi) &
!$OMP& PRIVATE(i_r1,i_r2,Order,r1,r2,xarg,BesselFunc,rk,rip,rkp)
!$OMP DO SCHEDULE(DYNAMIC)
		                Order = DBLE(k) + 0.5D0

			        DO i_r2 = 1, Npoint

				        r2 = RadMesh(i_r2)

		                	DO i_r1 = 1, Npoint

		                		r1 = RadMesh(i_r1)

		                		xarg = 2.0D0*r1*r2/(range**2)
		                		CALL bessel(xarg, Order, BesselFunc, rk, rip, rkp)

		                		BesselTable(i_r1, i_r2, k, igauss) = BesselFunc * EXP(-(r1**2+r2**2)/(range**2)) &
		                		                                                * SQRT(0.5D0*pi/xarg)

		                	END DO
		                END DO
!$OMP End Do
!$OMP End Parallel
                        END DO ! end of order
                END DO ! end of # gaussians

		WRITE(*,'("...DONE")')

		RETURN
	END SUBROUTINE BesselTabularize

	SUBROUTINE BesselFree()

	        DEALLOCATE(BesselTable)

		RETURN
        END SUBROUTINE BesselFree

        !-------------------------------------------------------------------------------!
	!  In this subroutine, we calculate the one-body kinetic energy matrix elements	!
	!  as well as the one-body r.m.s radius matrix elements  			!
	!-------------------------------------------------------------------------------!

	SUBROUTINE Global_start

		INTEGER ta, la, na, nc
		DOUBLE PRECISION factor
		INTEGER nmaxi

		CHARACTER(64) filename
		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

	        !  Creating two new pointers to the 3D tensors EkField and R2Field and allocating
		!  memory for them

		CALL SymD3Tensor_new(EkField)
		CALL SymD3Tensor_new(R2Field)

	        !  Filling in these two pointers with the ONE-BODY kinetic energy and square radius
		!  matrix elements. We don't calculate <na nb | T | nc nd> here, but only things like
		!  <n'|T|n>. Same with the radius

		SELECT CASE(Basis)

		CASE(1)

			DO ta = 0, 1  ! Loop over isospin

				factor = 1.0D0 / (2.0D0 * m(ta))

				DO la = 0, N_0  ! Loop over "bra" orbital angular momentum

					nmaxi = ((N_0 - la) / 2) + 1

					DO na = 1, nmaxi  ! Sum over "bra" main quantum number
						EkField%d3tensor(la)%d2(na, na) = factor * nabla2HO(na - 1, na - 1, la)
						R2Field%d3tensor(la)%d2(na, na) = r2_cutonHO(na - 1, na - 1, la)
					END DO

					DO na = 1, nmaxi - 1

						EkField%d3tensor(la)%d2(na + 1, na) = factor * nabla2HO(na - 1, na, la)
						R2Field%d3tensor(la)%d2(na + 1, na) = r2_cutonHO(na - 1, na, la)

						DO nc = na + 2, nmaxi ! Loop over "ket" main quantum number
							EkField%d3tensor(la)%d2(nc, na) = 0.0D0
							R2Field%d3tensor(la)%d2(nc, na) = 0.0D0
						END DO

					END DO

				END DO
			END DO

		CASE(2)

			DO ta = 0, 1  ! Loop over isospin

				factor = 1.0D0 / (2.0D0 * m(ta))

				DO la = 0, Lmax  ! Loop over "bra" orbital angular momentum

							   nmaxi = Min(Nmax, NmaxOfL(la))
				 	IF (CompHO .EQ. 1) nmaxi = ((Lmax - la) / 2) + 1

					DO na = 1, nmaxi  ! Loop over "bra" main quantum number
						DO nc = 1, na  ! Loop over "bra" main quantum number
							EkField%d3tensor(la)%d2(na, nc) = factor * nabla2(na, nc, la)
							R2Field%d3tensor(la)%d2(na, nc) = r2_cuton(na, nc, la)
						END DO
					END DO

				END DO
			END DO

		END SELECT


                 ! Writing output: the one-body kinetic energy

		SELECT CASE (Basis)

		CASE(1)
			IF (N_0 < 10) THEN
				WRITE(filename, "(A,I1,A)") "data/Ek", N_0, "_HO.txt"
			ELSE
				WRITE(filename, "(A,I2,A)") "data/Ek", N_0, "_HO.txt"
			END IF
		CASE(2)
			IF (N_0 < 10) THEN
				WRITE(filename, "(A,I1,A)") "data/Ek", N_0, "_WS.txt"
			ELSE
				WRITE(filename, "(A,I2,A)") "data/Ek", N_0, "_WS.txt"
			END IF
		END SELECT

		OPEN(file_desc, FILE=filename, ACTION="WRITE", IOSTAT=file_error)

		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to write the results in file ", filename
		ELSE
			DO la = 0, Lmax
								     nmaxi = Min(Nmax, NmaxOfL(la))
				IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) nmaxi = ((N_0 - la) / 2) + 1

				DO na = 1, nmaxi
					DO nc = 1, na
						WRITE (file_desc, "(I3,I3,I3,E24.16)", IOSTAT=file_error) &
							la, na, nc, EkField%d3tensor(la)%d2(na, nc)
					END DO
				END DO
			END DO
			CLOSE(file_desc)
		END IF

                 ! Writing output: the one-body root mean square radius

		SELECT CASE (Basis)

		CASE(1)
			IF (N_0 < 10) THEN
				WRITE(filename, "(A,I1,A)") "data/R2", N_0, "_HO.txt"
			ELSE
				WRITE(filename, "(A,I2,A)") "data/R2", N_0, "_HO.txt"
			END IF
		CASE(2)
			IF (N_0 < 10) THEN
				WRITE(filename, "(A,I1,A)") "data/R2", N_0, "_WS.txt"
			ELSE
				WRITE(filename, "(A,I2,A)") "data/R2", N_0, "_WS.txt"
			END IF
		END SELECT

		OPEN(file_desc, FILE=filename, ACTION="WRITE", IOSTAT=file_error)

		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention:  Impossible to write the results in file ", filename
		ELSE
			DO la = 0, Lmax
								     nmaxi = Min(Nmax, NmaxOfL(la))
				IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) nmaxi = ((N_0 - la) / 2) + 1

				DO na = 1, nmaxi
					DO nc = 1, na
						WRITE (file_desc, "(I3,I3,I3,E24.16)", IOSTAT=file_error) &
							la, na, nc, R2Field%d3tensor(la)%d2(na, nc)
					END DO
				END DO
			END DO
			CLOSE(file_desc)
		END IF
		RETURN

	CONTAINS

                !
		!  Function calculating the matrix elements of the form
		!           < n'l' | nabla^2 | nl >
		!
		!  Refs: Sec. 5.1.1., Page 45,
		!        Appendix C, Page 115
		!

		FUNCTION nabla2HO(na, nc, la)
			DOUBLE PRECISION nabla2HO
			INTEGER, INTENT(IN) :: na, nc, la

			IF ((na + 1) .EQ. nc) THEN
				nabla2HO = sq(nc) * sq2(nc + la) /b_0**2
			ELSE IF (na .EQ. (nc + 1)) THEN
				nabla2HO = sq(na) * sq2(na + la) /b_0**2
			ELSE IF (na .eq. nc) THEN
				nabla2HO = (DBLE(2*na + la) + 1.5D0) /b_0**2
			ELSE
				nabla2HO = 0.0D0
			END IF

			RETURN
		END FUNCTION nabla2HO

                !
		!  Function calculating the matrix elements of the squared radius
		!           < n'l' | r^2 | nl >
		!
		!  Refs: Sec. 5.1.1., Page 45
		!

		FUNCTION r2_cutonHO(na, nc, la)
			DOUBLE PRECISION r2_cutonHO
			INTEGER, INTENT(IN) :: na, nc, la

			IF ((na + 1) .eq. nc) THEN
				r2_cutonHO = - sq(nc) * sq2(nc + la) *(b_0**2)
			ELSE IF (na .eq. (nc + 1)) THEN
				r2_cutonHO = - sq(na) * sq2(na + la) *(b_0**2)
			ELSE IF (na .eq. nc) THEN
				r2_cutonHO = (DBLE(2*na + la) + 1.5D0 ) *(b_0**2)
			ELSE
				r2_cutonHO = 0.0D0
			END IF

			RETURN
		END FUNCTION r2_cutonHO

	END SUBROUTINE Global_start

        !---------------------------------------------------------------------------------!
	!		                       		  			          !
	! Function giving the reduced matrix elements < na la || NABLA || nb lb > 	  !
	! for an arbitrary spherical basis (not necessarily the Harmonic Oscillator)      !
	!		                       		  			          !
        !---------------------------------------------------------------------------------!

	FUNCTION nabla2(na, nc, la)
		DOUBLE PRECISION nabla2
		INTEGER, INTENT(IN) :: na, nc, la
		INTEGER :: IndexBra, IndexKet, IndexLoop, n_loop, l_loop
		DOUBLE PRECISION :: sum

		IndexBra = IndexVecNL(na,la)
		IndexKet = IndexVecNL(nc,la)

		IF (IndexBra .EQ. 0 .OR. IndexKet .EQ. 0) THEN
		    	nabla2 = 0.0D0
		    	RETURN
		END IF

		sum = 0.0D0

		! The summation must go up to N_0 (complete basis). Even though, the matrix element
		! <lmax, nmax || Nabla || lmax, nmax> is wrong as it contains a term proportional to
		! lmax+1

		DO l_loop = 0, Lmax+1
                        DO n_loop = 1, Nunity

                                IndexLoop = IndexAux(n_loop, l_loop)

				IF (IndexLoop .NE. 0) THEN
					sum = sum + SymKineticEnergy2Body_nabla(na, la, n_loop, l_loop) &
						  * SymKineticEnergy2Body_nabla(n_loop, l_loop, nc, la)
				END IF

			END DO
		END DO

		nabla2 = sum / DBLE(2*la + 1)

		RETURN
	END FUNCTION nabla2

	FUNCTION r2_cuton(na, nc, la)
		DOUBLE PRECISION r2_cuton
		INTEGER, INTENT(IN) :: na, nc, la
		INTEGER :: IndexBra, IndexKet

		IndexBra = IndexVecNL(na,la)
		IndexKet = IndexVecNL(nc,la)

		IF (IndexBra .EQ. 0 .OR. IndexKet .EQ. 0) THEN
		    	r2_cuton = 0.0D0
		    	RETURN
		END IF

		r2_cuton = IntegralR(IndexBra,IndexKet)

		RETURN
	END FUNCTION r2_cuton

        !---------------------------------------------------------------------------------!
	!		                       		  			          !
	! Function giving the reduced matrix elements < na la || NABLA || nb lb > 	  !
	! for an arbitrary spherical basis (not necessarily the Harmonic Oscillator)      !
	!		                       		  			          !
        !---------------------------------------------------------------------------------!

	FUNCTION SymKineticEnergy2Body_nabla(na, la, nb, lb)
		DOUBLE PRECISION SymKineticEnergy2Body_nabla
		INTEGER, INTENT(IN) :: na, nb, la, lb
		INTEGER :: IndexBra, IndexKet

                IndexBra = IndexAux(na,la)
                IndexKet = IndexAux(nb,lb)

		IF (IndexBra .EQ. 0 .OR. IndexKet .EQ. 0) THEN
			SymKineticEnergy2Body_nabla = 0.0D0
			RETURN
		END IF

		IF (la .EQ. (lb + 1)) THEN
		   	SymKineticEnergy2Body_nabla = + IntegralA(IndexBra, IndexKet, lb)*sq(la)
		ELSE IF (lb .EQ. (la + 1)) THEN
		  	SymKineticEnergy2Body_nabla = - IntegralB(IndexBra, IndexKet, lb)*sq(lb)
		ELSE
		   	SymKineticEnergy2Body_nabla = 0.0D0
		END IF

		RETURN
	END FUNCTION SymKineticEnergy2Body_nabla

        !---------------------------------------------------------!
	!		                       		          !
	! Function giving the A-integral of Notes, Sec. 3.2.3     !
	!		                       		  	  !
        !---------------------------------------------------------!

	FUNCTION IntegralA(IndexBra, IndexKet, lb)
		DOUBLE PRECISION :: IntegralA, Result
		DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)
		INTEGER, INTENT(IN) :: IndexBra, IndexKet, lb
		INTEGER :: Ipoint

	       	ALLOCATE(Integrand(1:Npoint))

		DO Ipoint = 1, Npoint
			  Integrand(Ipoint) = WaveFun(Ipoint,IndexBra)*WaveDeri(Ipoint,IndexKet) - &
		 	        DBLE(lb + 1)* WaveFun(Ipoint,IndexBra)*WaveFun(Ipoint,IndexKet)/RadMesh(Ipoint)
		END DO

		CALL simps(Integrand,Npoint,MeshStep,Result)

		DEALLOCATE(Integrand)

		IntegralA = Result

		RETURN
	END FUNCTION IntegralA

        !---------------------------------------------------------!
	!		                       		          !
	! Function giving the B-integral of Notes, Sec. 3.2.3     !
	!		                       		  	  !
       	!---------------------------------------------------------!

	FUNCTION IntegralB(IndexBra, IndexKet, lb)
		DOUBLE PRECISION :: IntegralB, Result
		DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)
		INTEGER, INTENT(IN) :: IndexBra, IndexKet, lb
		INTEGER :: Ipoint

	       	ALLOCATE(Integrand(1:Npoint))

		DO Ipoint = 1, Npoint
		 	 Integrand(Ipoint) = WaveFun(Ipoint,IndexBra)*WaveDeri(Ipoint,IndexKet) + &
		     	           DBLE(lb)* WaveFun(Ipoint,IndexBra)*WaveFun(Ipoint,IndexKet)/RadMesh(Ipoint)
		END DO

		CALL simps(Integrand,Npoint,MeshStep,Result)

		DEALLOCATE(Integrand)

		IntegralB = Result

		RETURN
	END FUNCTION IntegralB

        !---------------------------------------------------------!
	!		                       		          !
	! Function giving the B-integral of Notes, Sec. 3.2.3     !
	!		                       		  	  !
       	!---------------------------------------------------------!

	FUNCTION IntegralR(IndexBra, IndexKet)
		DOUBLE PRECISION :: IntegralR, Result
		DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)
		INTEGER, INTENT(IN) :: IndexBra, IndexKet
		INTEGER :: Ipoint

	       	ALLOCATE(Integrand(1:Npoint))

		DO Ipoint = 1, Npoint
		 	 Integrand(Ipoint) = WaveFun(Ipoint,IndexBra)*WaveFun(Ipoint,IndexKet)*(RadMesh(Ipoint)**2)
		END DO

		CALL simps(Integrand,Npoint,MeshStep,Result)

		DEALLOCATE(Integrand)

		IntegralR = Result

		RETURN
	END FUNCTION IntegralR

	SUBROUTINE Global_del
		!TODO
	END SUBROUTINE Global_del

END MODULE global

!-----------------------------------------------------------------------!
!									!
!									!
!     RADIAL INTEGRALS: BRINK-BOKER FORCE				!
!									!
!  This module computes the radial integral of tge Brink-Boeker force 	!
!  in the 2 distinct cases of an "analytical" harmonic oscillator basis !
!  and a general Woods-Saxon (or other) basis.				!
!									!
!-----------------------------------------------------------------------!

 MODULE ibb

	USE input
	USE global
	USE lgfactor
	USE symtalm
	USE bessik

	IMPLICIT NONE

 CONTAINS

        !---------------------------------------------------------------------------------!
	! Function giving the radial integral IBB for the Brink-Boeker term in the case   !
	! of the spherical harmonic oscillator basis	  			          !
	! Ref.:
        !---------------------------------------------------------------------------------!

	FUNCTION IBrinkBookerHO(na, la, nb, lb, nc, lc, nd, ld, k, x)
		DOUBLE PRECISION :: IBrinkBookerHO
		INTEGER, INTENT(IN) :: na, la, nb, lb, nc, lc, nd, ld, k
		DOUBLE PRECISION, INTENT(IN) :: x ! x = mi(i) / b
		INTEGER :: N1max, N2max, N1min, N2min
		DOUBLE PRECISION :: d1, d2
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: PI_BB
		REAL(KIND = 16) :: sumN1
#else
		DOUBLE PRECISION :: PI_BB
		DOUBLE PRECISION :: sumN1
#endif
		INTEGER :: N1

		PI_BB = 2.0D0 * ( 4.0D0*ATAN(1.0D0) )**3

		N1max = na + nc + ((la + lc - k) / 2)
		N2max = nb + nd + ((lb + ld - k) / 2)
		N1min = MIN_5N(na, la, nc, lc, k)
		N2min = MIN_5N(nb, lb, nd, ld, k)

		! El sumatorio se deberia realizar en el sentido opuesto
		! para ir de menor a mayor en el orden de magnitud de los numeros

		d1 = x * x
		d2 = d1 + 2.0D0

		sumN1 = SymKumar_get(na, la, nc, lc, N1max, k) * SumBB(nb, lb, nd, ld, k, N1max, d2)

		DO N1 = N1max - 1, N1min, -1
			sumN1 = (SymKumar_get(na, la, nc, lc, N1, k) * SumBB(nb, lb, nd, ld, k, N1, d2)) &
				- (sumN1 * (DBLE(N1 + N2min + k) + 1.5d0) / DBLE(N1 + 1) / (DBLE(N1 + k) + 1.5d0) / d2)
		END DO

		IBrinkBookerHO = PI_BB * PAR(N1min + N2min) &
			* EXP(DDLogSemiFactorials(N1min + N2min + k) &
		            - DDLogFactorials(N1min) &
			    - DDLogSemiFactorials(N1min+ k) &
			    - DDLogFactorials(N2min) &
			    - DDLogSemiFactorials(N2min + k)) &
			* (x ** 3.0D0) * sumN1 / (d2 ** (DBLE(N1min + N2min + k) + 1.5d0))
		RETURN

	CONTAINS

       		!-----------------------------------------------------------------------!
		!  Sum over N2 of the c(N2,k) T(nb,lb,nd,ld,N2,k)*Integral(k,N1,)	1
       		!-----------------------------------------------------------------------!

		FUNCTION SumBB(n1, l1, n2, l2, k, M1, y) ! M1 <- N1
#if(USE_QUADRUPLE==1)
			REAL(KIND = 16) :: SumBB
#else
			DOUBLE PRECISION :: SumBB
#endif
			INTEGER, INTENT(IN) :: n1, l1, n2, l2, k, M1
			DOUBLE PRECISION, INTENT(IN) :: y

			INTEGER :: N, Nmin, Nmax

			Nmax = n1 + n2 + ((l1 + l2 - k) / 2)
			Nmin = MIN_5N(n1, l1, n2, l2, k)

			SumBB = SymKumar_get(n1, l1, n2, l2, Nmax, k)

			DO N = Nmax - 1, Nmin, -1
				SumBB = - (SumBB * (DBLE(N + M1 + k) + 1.5D0) / DBLE(N + 1) / (DBLE(N + k) + 1.5D0) / y) &
					+ SymKumar_get(n1, l1, n2, l2, N, k)
			END DO

			RETURN
		END FUNCTION SumBB

	END FUNCTION IBrinkBookerHO

        !---------------------------------------------------------------------------------!
	! Function giving the radial integral IBB for the Brink-Boeker term in the case   !
	! of a general spherical basis			  			          !
        !---------------------------------------------------------------------------------!

	FUNCTION IBrinkBooker(na, la, nb, lb, nc, lc, nd, ld, k, i, n_reg)
		DOUBLE PRECISION :: IBrinkBooker

		INTEGER, INTENT(IN) :: na, la, nb, lb, nc, lc, nd, ld, k, i, n_reg

		INTEGER :: IndexBra, IndexKet, Index1, Index2, IndexTest, i_r1

		DOUBLE PRECISION :: VBB, res, Pi, h, Functi
		DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)

		IndexBra = IndexVecNL(na,la)
		IndexKet = IndexVecNL(nc,lc)

                Index1 = IndexVecNL(nb, lb)
                Index2 = IndexVecNL(nd, ld)

                IndexTest = Index1*Index2*IndexBra*IndexKet

		IF (IndexTest .EQ. 0) THEN
			IBrinkBooker = 0.0D0
			RETURN
		END IF

		Pi = 4.0D0*ATAN(1.0D0)
		h = MeshStep

	        ALLOCATE(Integrand(1:Npoint))

!$OMP Parallel Default(None) &
!$OMP& SHARED(Npoint,Index1,Index2,k,i,n_reg,WaveFun,IndexBra,IndexKet,Integrand) &
!$OMP& PRIVATE(i_r1,VBB,Functi)
!$OMP DO SCHEDULE(DYNAMIC)
		DO i_r1 = 1, Npoint
			VBB = IntegralBessel(Index1, Index2, k, i_r1, i, n_reg)
			Functi = WaveFun(i_r1,IndexBra) * WaveFun(i_r1,IndexKet)
			Integrand(i_r1) = Functi * VBB
		END DO
!$OMP End Do
!$OMP End Parallel

		res = 0.0D0
		CALL simps(Integrand, Npoint, h, res)

		DEALLOCATE(Integrand)

		IBrinkBooker = 4.0D0*Pi * res

		RETURN

	CONTAINS

       		!---------------------------------------------------------------------------------------!
		!											!
		!											!
        	!---------------------------------------------------------------------------------------!

		FUNCTION IntegralBessel(IndexBra, IndexKet, k, i_r1, i, n_reg)
			DOUBLE PRECISION :: IntegralBessel

			INTEGER, INTENT(IN) :: IndexBra, IndexKet, k, i_r1, i, n_reg
			INTEGER :: i_r2

			DOUBLE PRECISION :: BesselFunc, res, h, Functi, Order, r
			DOUBLE PRECISION, ALLOCATABLE :: Integral(:)

	       		ALLOCATE(Integral(1:Npoint))

			BesselFunc = 0.0D0
			h = MeshStep
			Order = k + 0.5D0

			DO i_r2 = 1, Npoint

				BesselFunc = BesselTable(i_r1, i_r2, k, i)

				Functi = WaveFun(i_r2,IndexBra) * WaveFun(i_r2,IndexKet)

				r = ABS(RadMesh(i_r1) - RadMesh(i_r2))

				Integral(i_r2) = Functi * BesselFunc * (r**n_reg)

			END DO

			CALL simps(Integral, Npoint, h, res)

			DEALLOCATE(Integral)

			IntegralBessel = res

			RETURN
		END FUNCTION IntegralBessel

	END FUNCTION IBrinkBooker

END MODULE ibb
!---------------------------------------------------------------------!
!                                                                     !
!     RADIAL INTEGRALS: COULOMB FORCE                                 !
!                                                                     !
!---------------------------------------------------------------------!

 MODULE ic

	USE input
	USE global
	USE lgfactor
	USE symtalm
	USE math

	IMPLICIT NONE

	DOUBLE PRECISION, PARAMETER :: VC = 1.4399784085965135298D0

 CONTAINS

	!-----------------------------------------------------------------------!
	!     Radial Integral as defined in Appendix H for the spherical	!
	!     harmonic oscillator basis						!
	!-----------------------------------------------------------------------!

	FUNCTION ICoulombHO(na, la, nb, lb, nc, lc, nd, ld, k)
		DOUBLE PRECISION :: ICoulombHO
		INTEGER, INTENT(IN) :: na, la, nb, lb, nc, lc, nd, ld, k

		INTEGER :: N1, N1min, N1max, N2min, N2max
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: PI_COU
		REAL(KIND = 16) :: sumN1
#else
		DOUBLE PRECISION :: PI_COU
		DOUBLE PRECISION :: sumN1
#endif

		PI_COU = SQRT(2.0D0) * ( 4.0D0 * ATAN(1.0D0) )**(2.5D0)

		N1max = na + nc + ((la + lc - k) / 2)
		N2max = nb + nd + ((lb + ld - k) / 2)
		N1min = MIN_5N(na, la, nc, lc, k)
		N2min = MIN_5N(nb, lb, nd, ld, k)

		sumN1 = SymKumar_get(na, la, nc, lc, N1max, k) * SumC(nb, lb, nd, ld, k, N1max, DBLE(2.0))
		DO N1 = N1max - 1, N1min, -1
			sumN1 = (-sumN1 * (DBLE(N1 + N2min + k) + 0.5D0) / DBLE(N1 + 1) / (DBLE(N1 + k) + 1.5) / 2.0D0) &
				+ SymKumar_get(na, la, nc, lc, N1, k) * SumC(nb, lb, nd, ld, k, N1, DBLE(2.0))
		END DO

		ICoulombHO = PI_COU * PAR(N1min + N2min) &
			* EXP(DDLogSemiFactorials(N1min + N2min + k - 1) &
			    - DDLogFactorials(N1min) &
			    - DDLogSemiFactorials(N1min + k) &
			    - DDLogFactorials(N2min) &
			    - DDLogSemiFactorials(N2min + k)) &
			* (sumN1 * (2.0D0 ** DBLE(-N1min - N2min - k))) / b_0
		RETURN

	CONTAINS

		FUNCTION SumC(n1, l1, n2, l2, k, M1, y)
#if(USE_QUADRUPLE==1)
			REAL(KIND = 16) :: SumC
#else
			DOUBLE PRECISION :: SumC
#endif
			INTEGER, INTENT(IN) :: n1, l1, n2, l2, k, M1
			DOUBLE PRECISION, INTENT(IN) :: y

			INTEGER :: N, Nmin, Nmax

			Nmax = n1 + n2 + ((l1 + l2 - k) / 2)
	 		Nmin = MIN_5N(n1, l1, n2, l2, k)

			SumC = SymKumar_get(n1, l1, n2, l2, Nmax, k)
			DO N = Nmax - 1, Nmin, -1
				SumC = - (SumC * (DBLE(N + M1 + k) + 0.5D0) / DBLE(N + 1) / (DBLE(N + k) + 1.5D0) / y) &
					+ SymKumar_get(n1, l1, n2, l2, N, k)
			END DO
			RETURN
		END FUNCTION SumC

	END FUNCTION ICoulombHO

	!-----------------------------------------------------------------------!
	!     Radial 2-body Coulomb Integral as defined in Appendix H for the 	!
	!     spherical	GENERAL basis. Here we integrate over r1 the result of 	!
	!     the Coulomb multipole k integral over r2.				!
	!-----------------------------------------------------------------------!

	FUNCTION ICoulomb(na, la, nb, lb, nc, lc, nd, ld, k)
		DOUBLE PRECISION :: ICoulomb

		INTEGER, INTENT(IN) :: na, la, nb, lb, nc, lc, nd, ld, k

		INTEGER :: IndexBra, IndexKet, i
		DOUBLE PRECISION :: Vcou, res
		DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)

		IndexBra = IndexVecNL(na,la)
		IndexKet = IndexVecNL(nc,lc)

		IF (IndexBra .EQ. 0 .OR. IndexKet .EQ. 0) THEN
			ICoulomb = 0.0D0
			RETURN
		END IF

	       	ALLOCATE(Integrand(1:Npoint))

!$OMP Parallel Default(None) &
!$OMP& SHARED(Npoint,nb,lb,nd,ld,k,WaveFun,IndexBra,IndexKet,Integrand) &
!$OMP& PRIVATE(i,Vcou)
!$OMP DO SCHEDULE(DYNAMIC)
		DO i=1,Npoint
			Vcou = Integral1(nb, lb, nd, ld, k, i)
			Integrand(i) = WaveFun(i,IndexBra)*WaveFun(i,IndexKet)*Vcou
		END DO
!$OMP End Do
!$OMP End Parallel

		CALL simps(Integrand,Npoint,MeshStep,res)

		DEALLOCATE(Integrand)

		Icoulomb = FOUR_PI * res / DBLE(2*k+1)

		RETURN

	CONTAINS

		!-----------------------------------------------------------------------!
		! 									!
		!   Radial 1-body Coulomb Integral as defined in Appendix H for the 	!
		!   spherical GENERAL basis. Here we integrate the coulomb multipole	!
		!   k over r2. We split the integration in 2 parts, r2 < r1 and r2 > r1.!
		!   The Coulomb multipole is defined as:				!
		!									!
		!          W_k(r1, r2) = 4*pi/(2k+1) * r_min^k/r_max^(k+1)		!
		!									!
		!   with: r_min = min(r1, r2) and r_max = max(r1, r2)			!
		! 									!
		!-----------------------------------------------------------------------!

		FUNCTION Integral1(n1, l1, n2, l2, k, index)
#if(USE_QUADRUPLE==1)
			REAL(KIND = 16) :: Integral1
#else
			DOUBLE PRECISION :: Integral1
#endif

			INTEGER, INTENT(IN) :: n1, n2, l1, l2, k, index

			INTEGER :: IndexBra, IndexKet, i
			DOUBLE PRECISION :: r_min, r_max, Vcou, res
			DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)

			IndexBra = IndexVecNL(n1,l1)
			IndexKet = IndexVecNL(n2,l2)

			! Result equal to zero if the quantum numbers are beyond the limits

			IF (IndexBra .EQ. 0 .OR. IndexKet .EQ. 0) THEN
				Integral1 = 0.0D0
				RETURN
			END IF

	       		ALLOCATE(Integrand(1:Npoint))

			DO i=1,Npoint
				Integrand(i) = 0.0D0
			END DO

			! If r2 < r1 (<=> index < i), we integrate a certain expression...

			IF (index .GT. 1) THEN

				r_max = RadMesh(index)

				DO i=1,index-1
					r_min = RadMesh(i)
					Vcou = r_min**k / r_max**(k+1)
					Integrand(i) = WaveFun(i,IndexBra)*WaveFun(i,IndexKet)*Vcou
				END DO

			END IF

			! If r2 > r1 (<=> index > i), we integrate a slightly different expression...

			r_min = RadMesh(index)

			DO i=index, Npoint
				r_max = RadMesh(i)
				Vcou = r_min**k / r_max**(k+1)
				Integrand(i) = WaveFun(i,IndexBra)*WaveFun(i,IndexKet)*Vcou
			END DO

			CALL simps(Integrand,Npoint,MeshStep,res)

			DEALLOCATE(Integrand)

			Integral1 = res
			RETURN
		END FUNCTION Integral1

	END FUNCTION ICoulomb

END MODULE ic
!---------------------------------------------------------------------!
!                                                                     !
!     RADIAL INTEGRALS: SPIN-ORBIT FORCE                              !
!                                                                     !
!---------------------------------------------------------------------!

 MODULE ils

	USE input
	USE global
	USE symtalm

	IMPLICIT NONE

 CONTAINS

	!---------------------------------------------------------------------!
	! Radial Integral IPLS in the case of the harmonic oscillator basis   !
	! Refs.: PhD, Page 137, E18                                           !
	!---------------------------------------------------------------------!

	FUNCTION IPLSHO(na, nc, la, nb, nd, lb)
		DOUBLE PRECISION IPLSHO
		INTEGER, INTENT(IN) :: na, nc, la, nb, nd, lb

		INTEGER p, pmax
		DOUBLE PRECISION x
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: Log2
		REAL(KIND = 16) :: s
#else
		DOUBLE PRECISION :: Log2
		DOUBLE PRECISION :: s
#endif

		Log2 = LOG(2.0D0)

		pmax = na + nb
		x = la + lb - 0.5D0

		s = SymCoefficientB_get(na, la, nb, lb, pmax) * Sum1(nc, la, nd, lb, pmax + x)
		DO p = pmax - 1, 0, -1
			s = SymCoefficientB_get(na, la, nb, lb, p) * Sum1(nc, la, nd, lb, p + x) &
				+ 0.5D0 * (p + x + 1.0D0) * s
		END DO
		IPLSHO = s * EXP(DDLogSemiFactorials(la + lb - 1) - ((la + lb + 1.0D0) * Log2)) / SQRT(2.0D0) / (b_0 ** 5)
		RETURN

	CONTAINS

		FUNCTION Sum1(n1, l1, n2, l2, x)
#if(USE_QUADRUPLE==1)
			REAL(KIND = 16) Sum1
#else
			DOUBLE PRECISION Sum1
#endif
			INTEGER, INTENT(IN) :: n1, l1, n2, l2
			DOUBLE PRECISION, INTENT(IN) :: x

			INTEGER p, pmax

			pmax = n1 + n2
			Sum1 = SymCoefficientB_get(n1, l1, n2, l2, pmax)
			DO p = pmax - 1, 0, -1
				Sum1 = 0.5D0 * (p + x + 1.0D0) * Sum1 + SymCoefficientB_get(n1, l1, n2, l2, p)
			END DO
			RETURN
		END FUNCTION Sum1

	END FUNCTION IPLSHO

	!---------------------------------------------------------------------!
	! Radial Integral IPLS in the general case of a spherical basis       !
	!---------------------------------------------------------------------!

	FUNCTION IPLS(na, nc, la, nb, nd, lb)
		DOUBLE PRECISION IPLS
		INTEGER, INTENT(IN) :: na, nc, la, nb, nd, lb

		DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)

		INTEGER IndexBraOne, IndexKetOne, IndexBraTwo, IndexKetTwo, Ipoint

		DOUBLE PRECISION RFourth,Result

		IndexBraOne = IndexVecNL(na,la)
		IndexBraTwo = IndexVecNL(nb,lb)

		IndexKetOne = IndexVecNL(nc,la)
		IndexKetTwo = IndexVecNL(nd,lb)

		IF (IndexBraOne .EQ. 0 .OR. IndexKetOne .EQ. 0 .OR. IndexBraTwo .EQ. 0 .OR. IndexKetTwo .EQ. 0) THEN
			IPLS = 0.0D0
			RETURN
		END IF

	        ALLOCATE(Integrand(1:Npoint))

		DO Ipoint = 1, Npoint
			RFourth = RadMesh(Ipoint)**4
			Integrand(Ipoint) = WaveFun(Ipoint,IndexBraOne)*WaveFun(Ipoint,IndexBraTwo)* &
					    WaveFun(Ipoint,IndexKetOne)*WaveFun(Ipoint,IndexKetTwo)/RFourth
		END DO

		CALL simps(Integrand,Npoint,MeshStep,Result)

		DEALLOCATE(Integrand)

		IPLS = Result

		RETURN

	END FUNCTION IPLS

	!---------------------------------------------------------------------!
	! Radial Integral IHFLS in the case of the harmonic oscillator basis  !
	! Refs.: PhD, Page 137, Sec. E3.3                                     !
	!---------------------------------------------------------------------!

	FUNCTION IHFLSHO(na, nc, la, nb, nd, lb)
		DOUBLE PRECISION IHFLSHO
		INTEGER, INTENT(IN) :: na, nc, la, nb, nd, lb

                ! OJO: solo vale cuando b1==b2

		INTEGER p, pmax
		DOUBLE PRECISION x, y
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: Log2
		REAL(KIND = 16) :: s
#else
		DOUBLE PRECISION :: Log2
		DOUBLE PRECISION :: s
#endif

		Log2 = LOG(2.0D0)

		pmax = na + nc
		x = DBLE(la + lb) - 0.5D0
		y = DBLE(lb - la) - 0.5D0

		s = SymCoefficientB_get(na, la, nc, la, pmax) * Sum2(nb, lb, nd, lb, pmax + x, y - pmax)

		DO p = pmax - 1, 0, -1
			s = (SymCoefficientB_get(na, la, nc, la, p) * Sum2(nb, lb, nd, lb, p + x, y - p)) &
				+ 0.5D0 * (DBLE(p+1) + x) * s
		END DO
		IHFLSHO = s * EXP(DDLogSemiFactorials(la + lb - 1) - (DBLE(la + lb + 1) * Log2)) / SQRT(2.0D0) / (b_0 ** 5)
		RETURN

	CONTAINS

		FUNCTION Sum2(n1, l1, n2, l2, x, y)
#if(USE_QUADRUPLE==1)
			REAL(KIND = 16) Sum2
#else
			DOUBLE PRECISION Sum2
#endif
			INTEGER, INTENT(IN) :: n1, l1, n2, l2
			DOUBLE PRECISION, INTENT(IN) :: x, y

			INTEGER p, pmax

			pmax = n1 + n2
			Sum2 = SymCoefficientB_get(n1, l1, n2, l2, pmax) * (pmax + y)
			DO p = pmax - 1, 0, -1
				Sum2 = 0.5D0 * (DBLE(p+1) + x) * Sum2 + (SymCoefficientB_get(n1, l1, n2, l2, p) * (DBLE(p) + y))
			END DO
			RETURN
		END FUNCTION Sum2

	END FUNCTION IHFLSHO

	!---------------------------------------------------------------------!
	! Radial Integral IHFLS in the general case of a spherical basis      !
	!---------------------------------------------------------------------!

	FUNCTION IHFLS(na, nc, la, nb, nd, lb)
		DOUBLE PRECISION IHFLS
		INTEGER, INTENT(IN) :: na, nc, la, nb, nd, lb

		DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)

		INTEGER IndexBraOne, IndexKetOne, IndexBraTwo, IndexKetTwo, Ipoint

		DOUBLE PRECISION Result

		IndexBraOne = IndexVecNL(na,la)
		IndexBraTwo = IndexVecNL(nc,la)

		IndexKetOne = IndexVecNL(nb,lb)
		IndexKetTwo = IndexVecNL(nd,lb)

		IF (IndexBraOne .EQ. 0 .OR. IndexKetOne .EQ. 0 .OR. IndexBraTwo .EQ. 0 .OR. IndexKetTwo .EQ. 0) THEN
			IHFLS = 0.0D0
			RETURN
		END IF

	        ALLOCATE(Integrand(1:Npoint))

		DO Ipoint = 1, Npoint
			Integrand(Ipoint) = WaveFun(Ipoint,IndexBraOne) *WaveFun(Ipoint,IndexBraTwo)* &
					  ( WaveFun(Ipoint,IndexKetOne) *WaveDeri(Ipoint,IndexKetTwo) + &
					    WaveDeri(Ipoint,IndexKetOne)*WaveFun(Ipoint,IndexKetTwo) - &
			              2.0D0*WaveFun(Ipoint,IndexKetOne) *WaveFun(Ipoint,IndexKetTwo)/RadMesh(Ipoint) ) &
					   /RadMesh(Ipoint)**3
		END DO

		CALL simps(Integrand,Npoint,MeshStep,Result)

		DEALLOCATE(Integrand)

		IHFLS = Result

		RETURN

	END FUNCTION IHFLS

END MODULE ils
 MODULE indexx

 CONTAINS

	SUBROUTINE indexx_real8(n, arrin, indx)
		INTEGER, INTENT(IN) :: n
		DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: arrin
		INTEGER, DIMENSION(*), INTENT(INOUT) :: indx

		INTEGER l, j, ir, indxt, i
		DOUBLE PRECISION q

		DO j = 1, n
			indx(j) = j
		END DO

		IF (n .EQ. 1) RETURN

		l = (n / 2) + 1
		ir = n
		DO
			IF (l .GT. 1) THEN
				l = l - 1
				indxt = indx(l)
				q = arrin(indxt)
			ELSE
				indxt = indx(ir)
				q = arrin(indxt)
				indx(ir) = indx(1)
				ir = ir - 1
				IF (ir .EQ. 1) THEN
					indx(1) = indxt
					RETURN
				END IF
			END IF
			i = l
			j = l * 2
			DO WHILE (j .LE. ir)
			
				IF (j .LT. ir) THEN
					IF (arrin(indx(j)) .LT. arrin(indx(j + 1))) THEN
						j = j + 1
					END IF
				END IF
				
				IF (q .LT. arrin(indx(j))) THEN
					indx(i) = indx(j)
					i = j
					j = j + i
				ELSE
					j = ir + 1
				END IF
			END DO
			indx(i) = indxt
		END DO
	END SUBROUTINE indexx_real8

 END MODULE indexx
! Input parameters module
 MODULE input

	IMPLICIT NONE

	!-----------------------------------------------------------------------------------------------------------------------!
	!															!
	!					GLOBAL VARIABLES								!
	!															!
	!    Meaning of some parameters												!
	!      - Nsize	   : Maximum number of single-particle levels taken in the basis. Used for arrays declaration		!
	!      - Npoint	   : Number of points on the lattice									!
	!      - Lmax	   : Maximum orbital angular momentum l (determined automatically in the case of the HO basis)		!
	!      - Nmax	   : Maximum number of radial wave-functions (determined automatically in the case of the HO basis)	!
	!      - IndexWave : At convergence, the program gives the radial wave-functions of the s.p. in the 			!
	!                    canonical basis and of the q.p. in the HFB basis. IndexWave specifies which one should		!
	!                    be written on disk											!
	!      - Nunity	   : Number of radial wave-functions taken into account when the unity is inserted for the 1-body 	!
	!		     kinetic energy											!
	!      - IndexVecNL: Array that gives which basis wave-function corresponds to a given n and l				!
	!      - MeshStep  : Lattice discretization step denoted by h in the notes						!
	!      - RadMesh   : Vector containing the radius on the lattice							!
	!      - Energy	   : Single-particle energies of the basis								!
	!      - WaveFun   : Array containing the single-particle wave-functions of the basis states on the lattice		!
	!		     First index refers to the isospin, second to the lattice and third to the state			!
	!      - WaveDeri  : Same as above but for the derivatives of the basis wave-functions					!
	!      - anneal    : At each iteration, the density is taken as a mixing between the new density and the old one. This	!
	!		     is controlled by the anneal parameter								!
	!															!
	!-----------------------------------------------------------------------------------------------------------------------!

	! Input parameters

	LOGICAL, PUBLIC :: regularized_Gaussian = .False.

	INTEGER, PUBLIC :: protons = 8, neutrons = 8    ! Number of protons and neutrons
	INTEGER, PUBLIC :: N_0 = 8                   	! Number of shells
	INTEGER, PUBLIC :: NITER_MAX = 10               ! Maximum number of iterations
	INTEGER, PUBLIC :: HFonly = 0                   ! Switches to HF mode (1) instead of default HFB (0)
	INTEGER, PUBLIC :: Gogny = 1        		! 1 = D1S, 2 = D1, 3 = D1prime, 4 = D1N
	INTEGER, PUBLIC :: Basis = 1			! 1 = HO basis, 2 = arbitrary basis
	INTEGER, PUBLIC :: CompHO = 1			! 1 = compatibility with the HO basis in terms of quantum numbers
						        ! 0 = Lmax and Nmax required
	INTEGER, PUBLIC :: Optimization = 0             ! Optimizes calculation of matrix elements for WS-like bases
	INTEGER, PUBLIC :: Nunity = 20                  ! Maximum number of n values for resolution of the identity
	INTEGER, PUBLIC :: n_reg_max = 0                ! Maximum number of Gaussians used for regularized delta forces
	INTEGER, PUBLIC :: Nsize = 2, Npoint            ! Npoint: number of points of the mesh
	INTEGER, PUBLIC :: Lmax, Nmax                   ! Lmax: maximum number of L values in basis, Nmax: maximum values of n
	INTEGER, PUBLIC :: Isospin = 1                  ! Isospin of the basis functions (1: neutrons, 2: protons)
        INTEGER, PUBLIC :: IsoFacPair = 2               ! multiplies pairing fields for protons (0), neutrons (1) or both (2)
        INTEGER, PUBLIC :: IsoShiftLambda = 2           ! shifts particle number condition for p (0), n(1) or both (2)

        INTEGER, PUBLIC :: switch_Coulomb = 2           ! 0: no Coulomb, 1: direct only, 2: direct + full exchange
        INTEGER, PUBLIC :: switch_LS = 1                ! 0: no spin-orbit, 1: with spin-orbit
        INTEGER, PUBLIC :: switch_CM = 1                ! 0: no 2-body center of mass, 1: with two-body center of mass
        INTEGER, PUBLIC :: switch_DD = 1                ! 0: no density dependence, 1: with density dependence

	DOUBLE PRECISION, PUBLIC :: b_0 = 1.6D0      	! Oscillator length ("b")
        DOUBLE PRECISION, PUBLIC :: Ecut = -1.0D0       ! Cut-off in energy of the basis
        DOUBLE PRECISION, PUBLIC :: Convergence = 1.D-3 ! Criterion to stop convergence
        DOUBLE PRECISION, PUBLIC :: anneal = 0.5D0      ! Convergence
        DOUBLE PRECISION, PUBLIC :: facPair = 1.0D0     ! Value of pairing enhancement if IsoFacPair > 2
        DOUBLE PRECISION, PUBLIC :: ShiftLambda = 0.0D0 ! Value of shifted particle number if IsoShiftLambda>0
        DOUBLE PRECISION, PUBLIC :: RadPoles = 1.0D0    ! ...?
	DOUBLE PRECISION, PUBLIC :: range1 = 0.7D0, range2 = 1.2D0 ! Range of the Gaussians (initialized to D1S values)


	INTEGER, PUBLIC :: Nlevel(1:2), IndexWave = 1

	INTEGER, ALLOCATABLE, PUBLIC :: LmaxiIso(:), NmaxOfL(:)
	INTEGER, ALLOCATABLE, PUBLIC :: Nmain(:,:), Lmain(:,:), Jmain(:,:), IndexVecNL(:,:), IndexAux(:,:)

	DOUBLE PRECISION, PUBLIC :: MeshStep
	DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: RadMesh(:)
	DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: Energy(:,:), WaveFun(:,:), WaveDeri(:,:)

 CONTAINS

	!-----------------------------------------------------------------------!
	!									!
	!       This subroutine reads the file "input.txt" that contains 	!
	!       the standard input to the program.				!
	!									!
	!-----------------------------------------------------------------------!

	SUBROUTINE Input_read

		INTEGER, PARAMETER :: file_desc = 17
		INTEGER :: file_error
		CHARACTER(LEN = 256) :: param
		INTEGER :: i_reg, param_len

		OPEN(file_desc, FILE="data/input.txt", ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) STOP "Impossible to read: ./data/input.txt"

		param_len = IO_get_Param(file_desc, param)

		DO WHILE (param_len .GT. 0)
			IF (param .EQ. "b_0") THEN
				b_0 = IO_get_RealValue(file_desc)
			ELSE IF (param .EQ. "Slowing-factor") THEN
				anneal = IO_get_RealValue(file_desc)
			ELSE IF (param .EQ. "Iterations") THEN
				NITER_MAX = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "Convergence") THEN
				Convergence = IO_get_RealValue(file_desc)

			ELSE IF (param .EQ. "Optimization") THEN
				Optimization = IO_get_IntegerValue(file_desc)
                        ELSE IF (param .EQ. "Divergence") THEN
                                RadPoles = IO_get_RealValue(file_desc)
			ELSE IF (param .EQ. "IsoBasis") THEN
				Isospin = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "HFOnly") THEN
				HFOnly = IO_get_IntegerValue(file_desc)
                        ELSE IF (param .EQ. "FactorPairing") THEN
				facPair = IO_get_RealValue(file_desc)
                        ELSE IF (param .EQ. "IsoFactorPairing") THEN
				IsoFacPair = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "ShiftLambda") THEN
                                ShiftLambda = IO_get_RealValue(file_desc)
                        ELSE IF (param .EQ. "IsoShiftLambda") THEN
                                IsoShiftLambda = IO_get_IntegerValue(file_desc)

			ELSE IF (param .EQ. "N_0") THEN
				N_0 = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "Lmax") THEN
				Lmax = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "Nmax") THEN
				Nmax = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "Unity") THEN
				Nunity = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "IndexWave") THEN
				IndexWave = IO_get_IntegerValue(file_desc)
                        ELSE IF (param .EQ. "Ecut") THEN
                                Ecut = IO_get_RealValue(file_desc)
			ELSE IF (param .EQ. "N") THEN
				neutrons = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "Z") THEN
				protons = IO_get_IntegerValue(file_desc)

			ELSE IF (param .EQ. "GOGNY") THEN
				Gogny = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "Regularized") THEN
				i_reg = IO_get_IntegerValue(file_desc)
				IF (i_reg .EQ. 1) THEN
				        regularized_Gaussian = .True.
				ELSE
				        regularized_Gaussian = .False.
				END IF
			ELSE IF (param .EQ. "OrderExpansion") THEN
				n_reg_max = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "Range1") THEN
				range1 = IO_get_RealValue(file_desc)
			ELSE IF (param .EQ. "Range2") THEN
				range2 = IO_get_RealValue(file_desc)

			ELSE IF (param .EQ. "Coulomb") THEN
				switch_Coulomb = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "SpinOrbit") THEN
				switch_LS = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "CenterOfMass") THEN
				switch_CM = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "DensityDependence") THEN
				switch_DD = IO_get_IntegerValue(file_desc)

			ELSE IF (param .EQ. "Basis") THEN
				Basis = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "CompHO") THEN
				CompHO = IO_get_IntegerValue(file_desc)
			ELSE
				WRITE(6,'("Attention! Unknown parameter: ",A)') param(1:param_len)
				! We ignore the remaining lines
				READ (file_desc, "(A)") param
			END IF
			param_len = IO_get_Param(file_desc, param)
		END DO

		CLOSE(file_desc)

		IF (Basis .EQ. 1) CompHO = 1

		! Forcing Lmax to N_0 if oscillator basis or compatibility required
		IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) THEN
                        Lmax = N_0
                        Nmax = N_0
                        Nsize = (2*Lmax+1)*(2*Nmax+1)
                END IF

		! Forcing N_0 to Lmax if general basis used
		IF (Basis .EQ. 2 .OR. CompHO .EQ. 0) THEN
                        N_0 = Lmax
                END IF

		! In case of a basis compatible with the harmonic oscillator basis, one forces
		! b (the oscillator length) to take a value proportional to the mass
		IF (b_0 .LT. 0.0D0 .AND. (Basis .EQ. 1 .OR. CompHO .EQ. 1)) b_0 = 1.05D0*DBLE(protons + neutrons)**(1.0D0/6.0D0)

		! Printing input data
		WRITE(6,'(5X,"DATA READ FROM INPUT")')
		WRITE(6,'(5X,"====================",/)')
		WRITE(6,'("Oscillator Length (b_0) ........: ",F15.12)') b_0
		WRITE(6,'("Slowing factor .................: ",F8.5)') anneal
		WRITE(6,'("Maximum Number of Iterations ...: ",I3)') NITER_MAX
		WRITE(6,'("Convergence in Energy ..........: ",ES9.2)') Convergence
		WRITE(6,'("Optimization for BB Mat. Els. ..: ",I3)') Optimization
		WRITE(6,'("Pairing Enhancement Factor .....: ",F8.5)') facPair
		WRITE(6,'("Shift in Lagrange Parameter ....: ",F8.5)') ShiftLambda
		WRITE(6,'("Isospin of the Basis WF ........: ",I3)') Isospin
		WRITE(6,'("Index of the Wave-function .....: ",I3)') IndexWave
		WRITE(6,'("Pure HF Calculations (=1) ......: ",I3)') HFOnly
		WRITE(6,'("Number of shells (N_0) .........: ",I3)') N_0
		WRITE(6,'("Maximum Ang. Mom. (Lmax) .......: ",I3)') Lmax
		WRITE(6,'("Maximum n (nmax) ...............: ",I3)') Nmax
		WRITE(6,'("Maximum n for Unity ............: ",I3)') Nunity
		WRITE(6,'("Index of the Wave-function .....: ",I3)') IndexWave
		WRITE(6,'("Cut-off energy (if > 0) ........: ",F8.5)') Ecut
		WRITE(6,'("Number of neutrons (N) .........: ",I3)') neutrons
		WRITE(6,'("Number of protons  (Z) .........: ",I3)') protons
		IF (regularized_Gaussian) THEN
		        WRITE(6,'("Regularized Gaussians ...........")')
		        WRITE(6,'("    - order of expansion .......: ",I2)') n_reg_max
		        WRITE(6,'("    - range1 ...................: ",F10.7)') range1
		ELSE
		        WRITE(6,'("Original Gogny interaction.......")')
		        WRITE(6,'("    - range1 ...................: ",F10.7)') range1
		        WRITE(6,'("    - range2 ...................: ",F10.7)') range2
		        SELECT CASE (Gogny)
				CASE (1)
					WRITE(6,'("    - Parametrization ..........: D1S")')
				CASE (2)
					WRITE(6,'("    - Parametrization ..........: D1")')
				CASE (3)
					WRITE(6,'("    - Parametrization ..........: D1p")')
				CASE (4)
					WRITE(6,'("    - Parametrization ..........: D1N")')
				CASE DEFAULT
					STOP "Invalid parameters of the Gogny force"
			END SELECT
		END IF
		IF (switch_CM .EQ. 0) THEN
			WRITE(6,'("    - 2-body center of mass correction turned off")')
		END IF
		IF (switch_LS .EQ. 0) THEN
			WRITE(6,'("    - Spin-orbit term turned off")')
		END IF
		IF (switch_DD .EQ. 0) THEN
			WRITE(6,'("    - Density-dependent term turned off")')
		END IF
		IF (switch_Coulomb .EQ. 0) THEN
			WRITE(6,'("    - Both direct and exchange Coulomb terms switched off")')
		END IF
		IF (switch_Coulomb .EQ. 1) THEN
			WRITE(6,'("    - Exchange Coulomb terms switched off")')
		END IF
		WRITE(6,'("Type of Basis ..................: ",I3)') Basis
		WRITE(6,'("Compatibility with the HO ......: ",I3)') CompHO

		RETURN

	END SUBROUTINE Input_read

        !-----------------------------------------------------------------------!
	!									!
        !   Subroutine that reads the spherical basis. It gives the total 	!
	!   number of s.p. levels, the s.p. energies, wave-functions, deri-	!
	!   vatives of the wave-function and quantum numbers n, l and j (the 	!
	!   latter is not used in this version of the code			!
	!									!
        !-----------------------------------------------------------------------!

        SUBROUTINE ReadBasis()

		INTEGER, PARAMETER :: file_unit = 12

	        INTEGER, ALLOCATABLE :: Ndummy(:)

		INTEGER :: file_error
		INTEGER :: stat, icount
		INTEGER :: i, j, L, N, nb1, nb2, jtemp, kappa, Nmaxi, Lmaxi, Lref, ntheo

		DOUBLE PRECISION :: Emax, E

		CHARACTER(LEN = 6)   :: keyword

                ! Opening the file to read
		OPEN(file_unit, FILE="basis.txt", ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) STOP "Impossible to open: basis.txt"

		READ(file_unit, FMT="(A)", ADVANCE="NO", IOSTAT=stat) keyword

                ! We read the input file by using keywords
		DO WHILE (keyword .NE. "FINISH")

                        SELECT CASE (keyword)

				CASE ("boxpar")

					READ(file_unit, *, IOSTAT=stat) Npoint,MeshStep

					IF (ALLOCATED(RadMesh)) DEALLOCATE(RadMesh)
					ALLOCATE(RadMesh(0:Npoint))
					DO i = 0, Npoint
						RadMesh(i) = DBLE(i)*MeshStep
					END DO

				CASE ("noflev")

					READ(file_unit, *, IOSTAT=stat) nb1,nb2,Nlevel(1),Nlevel(2)

					Nsize = max(Nlevel(1),Nlevel(2))

				!
				! NEUTRONS - NEUTRONS - NEUTRONS - NEUTRONS - NEUTRONS - NEUTRONS
				!

				! Neutron single-particle energies
				CASE ("esingn")

					IF (Nlevel(1) .EQ. 0) STOP "Nlevel = 0 in ReadBasis !"

					ALLOCATE(Energy(Nsize,2))

					READ(file_unit, *, IOSTAT=stat) (Energy(i,1),i=1,Nlevel(1))

				! Neutron angular momentum j (in Jmain)
				CASE("spin_n")

					IF (Nlevel(1) .EQ. 0) STOP "Nlevel = 0 in ReadBasis !"

					ALLOCATE(Jmain(Nsize,2))
					ALLOCATE(Ndummy(Nlevel(1)))

					READ(file_unit, *, IOSTAT=stat) (Ndummy(i),i=1,Nlevel(1))

					DO i=1,Nlevel(1)
						Jmain(i,1) = 2*Ndummy(i) - 1
					END DO

					DEALLOCATE(Ndummy)

				! Neutron orbital angular momentum l (in Lmain)
				CASE ("kappan")

					IF (Nlevel(1) .EQ. 0) STOP "Nlevel = 0 in ReadBasis !"

					ALLOCATE(Lmain(Nsize,2))
					ALLOCATE(Ndummy(Nlevel(1)))

					READ(file_unit, *, IOSTAT=stat) (Ndummy(i),i=1,Nlevel(1))

					! Lmaxi is the maximum orbital angular momentum (used to define
					! a fake number of shells N0)
					Lmaxi = 0

					DO i=1,Nlevel(1)
						kappa = Ndummy(i)
						jtemp = iabs(kappa)
						Lmain(i,1) = jtemp - (1 - kappa/jtemp)/2
						IF (Lmain(i,1) .GT. Lmaxi) Lmaxi = Lmain(i,1)
					END DO

					DEALLOCATE(Ndummy)

					ALLOCATE(LmaxiIso(2))
					LmaxiIso(1) = Lmaxi

				CASE ("numbrn")

					IF (Nlevel(1) .EQ. 0) STOP "Nlevel = 0 in ReadBasis !"

					ALLOCATE(Nmain(Nsize,2))

					READ(file_unit, *, IOSTAT=stat) (Nmain(i,1),i=1,Nlevel(1))

				CASE ("wave_n")

					IF (Nlevel(1) .EQ. 0) STOP "Nlevel = 0 in ReadBasis !"

					ALLOCATE(WaveFun(0:Npoint,Nsize))

					IF (Basis .EQ. 2 .AND. Isospin .EQ. 1) THEN

						DO j=1,Nlevel(1)
							READ(file_unit, *, IOSTAT=stat) (WaveFun(i,j),i=0,Npoint)
						END DO

					END IF

			  	CASE ("derivn")

					IF (Nlevel(1) .EQ. 0) STOP "Nlevel = 0 in ReadBasis !"

					ALLOCATE(WaveDeri(0:Npoint,Nsize))

					IF (Basis .EQ. 2 .AND. Isospin .EQ. 1) THEN

			      			DO j=1,Nlevel(1)
							READ(file_unit, *, IOSTAT=stat) (WaveDeri(i,j),i=0,Npoint)
						END DO

					END IF

				!
                          	! PROTONS - PROTONS - PROTONS - PROTONS - PROTONS - PROTONS
			  	!

				CASE ("esingp")

					READ(file_unit, *, IOSTAT=stat) (Energy(i,2),i=1,Nlevel(2))

				CASE( "spin_p")

					ALLOCATE(Ndummy(Nlevel(2)))

					READ(file_unit, *, IOSTAT=stat) (Ndummy(i),i=1,Nlevel(2))

					DO i=1,Nlevel(2)
						Jmain(i,2) = 2*Ndummy(i) - 1
					END DO

					DEALLOCATE(Ndummy)

				CASE ("kappap")

					ALLOCATE(Ndummy(Nlevel(2)))

					READ(file_unit, *, IOSTAT=stat) (Ndummy(i),i=1,Nlevel(2))

					Lmaxi = 0

					DO i=1,Nlevel(2)
						kappa = Ndummy(i)
						jtemp = iabs(kappa)
						Lmain(i,2) = jtemp - (1 - kappa/jtemp)/2
						IF (Lmain(i,2) .GT. Lmaxi) Lmaxi = Lmain(i,2)
					END DO

					DEALLOCATE(Ndummy)

					LmaxiIso(2) = Lmaxi

				CASE ("numbrp")

					READ(file_unit, *, IOSTAT=stat) (Nmain(i,2),i=1,Nlevel(2))

				CASE ("wave_p")

					IF (Basis .EQ. 2 .AND. Isospin .EQ. 2) THEN

						DO j=1,Nlevel(2)
							READ(file_unit, *, IOSTAT=stat) (WaveFun(i,j),i=0,Npoint)
						END DO

					END IF

			  	CASE ("derivp")

					IF (Basis .EQ. 2 .AND. Isospin .EQ. 2) THEN

						DO j=1,Nlevel(2)
							READ(file_unit, *, IOSTAT=stat) (WaveDeri(i,j),i=0,Npoint)
						END DO

					END IF

			  END SELECT

                          ! Until we get to the end of the file, we read the next keyword
		          READ(file_unit, FMT="(A)", IOSTAT=stat) keyword

                          ! Before we exit this subroutine, we determine a few useful quantities
			  !	- First of all, we establish the correspondence between the index basis wave-function
			  !	  and the quantum numbers n and l which are used in the calculation of the matrix elements
			  !	- We also determine have the choice between taking all states with E < Ecut and N < Nmax
			  !       or all stated with L < Lmax and N < Nmax.
			  !	- IndexAux is a special array used for the calculation of 1-body kinetic energy when inserting
			  !	  the unity: We need more wave-functions (more n values) for this.

		          IF (stat .NE. 0) THEN

                          	IF (Ecut .LT. 0.0D0) THEN

                                        ! IndexVecNL gives the correspondence N, L -> index i of the basis
                		        IF (.NOT.ALLOCATED(IndexVecNL)) ALLOCATE(IndexVecNL(0:2*Nmax,0:2*Lmax))

                                        ALLOCATE(IndexAux(0:Nunity,0:2*Lmax))

                                        IF (basis .EQ. 1) THEN
                                        	icount = 0
                                        	DO N = 0,2*Nmax
                                        		DO L = 0,2*Lmax
						        	icount = icount + 1
						        	IndexVecNL(N,L) = icount
						        END DO
				        	END DO
				        ELSE
			        		DO N = 0,2*Nmax
						        DO L = 0,2*Lmax
					        		IndexVecNL(N,L) = 0
					        	END DO
					        END DO
				        END IF

                                        DO N=0,Nunity
                                                DO L=0,2*Lmax
                                                        IndexAux(N,L) = 0
                                                END DO
                                        END DO

			      	        Emax = -1000.0D0

				        DO i=1,Nlevel(Isospin)
				                N = Nmain(i,Isospin)
				                L = Lmain(i,Isospin)
					        IF (N .LE. Nmax .AND. L .LE. Lmax) IndexVecNL(N,L) = i
					        IF (N .LE. Nmax .AND. L .LE. Lmax .AND. Energy(i,Isospin) .GT. Emax) Emax = Energy(i,Isospin)
                                                IF (N .LE. Nunity .AND. L .LE. (Lmax+1)) IndexAux(N,L) = i
				        END DO

                                ELSE

                                       DO i=1,Nlevel(Isospin)
			                        N = Nmain(i,Isospin)
				                L = Lmain(i,Isospin)
                                                E = Energy(i,Isospin)
					        IF (E .LT. Ecut .AND. L .LE. Lmax .AND. N .GT. Nmax) Nmax = N
				        END DO

                                        ! IndexVecNL gives the correspondence N, L -> index i of the basis
                		        ALLOCATE(IndexVecNL(0:2*Nmax,0:2*Lmax))

                                        ALLOCATE(IndexAux(0:Nunity,0:2*Lmax))

                                        IF (basis .EQ. 1) THEN
                                        	icount = 0
			        		DO N = 0,2*Nmax
						        DO L = 0,2*Lmax
						        	icount = icount + 1
						        	IndexVecNL(N,L) = icount
						        END DO
				        	END DO
				        ELSE
			        		DO N = 0,2*Nmax
						        DO L = 0,2*Lmax
						        	IndexVecNL(N,L) = 0
						        END DO
				        	END DO
				        END IF

                                        DO N=0,Nunity
                                                DO L=0,2*Lmax
                                                        IndexAux(N,L) = 0
                                                END DO
                                        END DO

			      	        Emax = Ecut

				        DO i=1,Nlevel(Isospin)
				                N = Nmain(i,Isospin)
				                L = Lmain(i,Isospin)
                                                E = Energy(i,Isospin)
					        IF (E .LT. Ecut .AND. L .LE. Lmax) IndexVecNL(N,L) = i
                                                IF (N .LE. Nunity .AND. L .LE. (Lmax+1)) IndexAux(N,L) = i
				        END DO


                                END IF

                		ALLOCATE(NmaxOfL(0:Lmax))

			 	 ! Printing the characteristics of the basis with respect to quantum numbers
				WRITE(6,'(/,5X,"DEFINING THE BASIS - CUT-OFFS")')
				WRITE(6,'(5X,"=============================",/)')
				IF (Basis.EQ.1) THEN
					WRITE(6,'("Analytical HO basis ............. ")')
				ELSE
					WRITE(6,'("Numerical basis on a mesh ....... ")')
					WRITE(6,'("Mesh size ......................: ",F10.5)') MeshStep
					WRITE(6,'("Box radius .....................: ",F10.5)') RadMesh(Npoint)
					WRITE(6,'("Number of mesh points ..........: ",I5)') Npoint+1
				END IF
				WRITE(6,'("Compatibility with the HO ......: ",I3)') CompHO
                                WRITE(6,'("Actual Cut-off Energy Ecut .....: ",F10.3)') Emax
                                WRITE(6,'("Maximum L value lmax ...........: ",I7)') Lmax
                                WRITE(6,'("Maximum n value ................: ",I7)') Nmax
				WRITE(6,'("Maximum n-value of a give l ....",//,"       L    n  n(HO)")')

				DO Lref = 0, Lmax

					Nmaxi = 0

					DO i=1, Nlevel(Isospin)
						N = Nmain(i,Isospin)
						L = Lmain(i,Isospin)
                                                E = Energy(i,Isospin)
						IF (L .EQ. Lref .AND. E .LE. Emax) THEN
					                IF (N .GE. Nmaxi) Nmaxi = N
						END IF
					END DO

					ntheo = (N_0 - Lref)/2 + 1

					IF (Ecut .LT. 0.0) NmaxOfL(Lref) = MIN(Nmaxi, Nmax)
					IF (Ecut .GT. 0.0) NmaxOfL(Lref) = Nmaxi

					IF (CompHO .EQ. 1) NmaxOfL(Lref) = ntheo

					WRITE(6,'(3X,3I5)') Lref,Nmaxi,ntheo
				END DO

		          	CLOSE(file_unit)

		          	RETURN

		          END IF

		END DO

		RETURN
	END SUBROUTINE ReadBasis

	FUNCTION IO_get_Param(fd, param)
		INTEGER IO_get_Param
		INTEGER, INTENT(IN) :: fd ! Descriptor de fichero
		CHARACTER(*), INTENT(INOUT) :: param

		CHARACTER (LEN=1) c
		INTEGER num, stat

		READ (fd, FMT="(A1)", ADVANCE="NO", IOSTAT=stat) c
		IF (stat .NE. 0) THEN
			IO_get_Param = 0
			RETURN
		END IF
		! Ignoramos los espacios en blanco iniciales
		DO WHILE (c == " ")
			READ (fd, FMT="(A1)", ADVANCE="NO", IOSTAT=stat) c
			IF (stat .NE. 0) THEN
				IO_get_Param = 0
				RETURN
			END IF
		END DO
		num = 0
		DO WHILE ((c .NE. " ") .AND. (c .NE. "="))
			IF (num .EQ. 0) THEN
				param = c
			ELSE
				param = param(1:num) // c
			END IF
			num = num + 1
			READ (fd, FMT="(A1)", ADVANCE="NO", IOSTAT=stat) c
			IF (stat .NE. 0) THEN
				IO_get_Param = 0
				RETURN
			END IF
		END DO
		! Si no se encontra ningun nombre de parametro, salimos con error
		IF (num .EQ. 0) THEN
			IO_get_Param = 0
			RETURN
		END IF

		DO WHILE (c .EQ. " ")
			READ (fd, FMT="(A1)", ADVANCE="NO", IOSTAT=stat) c
			IF (stat .NE. 0) THEN
				IO_get_Param = 0
				RETURN
			END IF
		END DO

		IF (c .EQ. "=") THEN
			IO_get_Param = num
		ELSE
			IO_get_Param = 0
		END IF
		RETURN
	END FUNCTION IO_get_Param

	FUNCTION IO_get_IntegerValue(fd)
		INTEGER IO_get_IntegerValue
		INTEGER, INTENT(IN) :: fd ! Descriptor de fichero

		INTEGER stat

		READ (fd, *, IOSTAT=stat) IO_get_IntegerValue
		IF (stat .NE. 0) STOP "Imposible leer parametro de entrada"
		RETURN
	END FUNCTION IO_get_IntegerValue

	FUNCTION IO_get_RealValue(fd)
		DOUBLE PRECISION IO_get_RealValue
		INTEGER, INTENT(IN) :: fd ! Descriptor de fichero

		INTEGER stat

		READ (fd, *, IOSTAT=stat) IO_get_RealValue
		IF (stat .NE. 0) STOP "Imposible leer parametro de entrada"
		RETURN
	END FUNCTION IO_get_RealValue

	FUNCTION IO_get_String(fd)
		CHARACTER(LEN=30) IO_get_String
		INTEGER, INTENT(IN) :: fd ! Descriptor de fichero

		INTEGER stat

		READ (fd, *, IOSTAT=stat) IO_get_String
		IF (stat .NE. 0) STOP "Imposible leer parametro de entrada"
		RETURN
	END FUNCTION IO_get_String

END MODULE input
 MODULE jacobi

 CONTAINS

	!*************************************************************
	!* This subroutine computes all eigenvalues and eigenvectors *
	!* of a real symmetric square matrix A(N,N). On output, ele- *
	!* ments of A above the diagonal are destroyed. D(N) returns *
	!* the eigenvalues of matrix A. V(N,N) contains, on output,  *
	!* the eigenvectors of A by columns. THe normalization to    *
	!* unity is made by main program before printing results.    *
	!* NROT returns the number of Jacobi matrix rotations which  *
	!* were required.                                            *
	!* --------------------------------------------------------- *
	!* Ref.:"NUMERICAL RECIPES, Cambridge University Press,1986, *
	!*       chap. 11, pages 346-348".                           *
	!*************************************************************
	SUBROUTINE Jacobi_real8(A, N, D, V, NROT)
		INTEGER, INTENT(IN) :: N
		DOUBLE PRECISION, DIMENSION(N,  N), INTENT(INOUT) :: A, V
		DOUBLE PRECISION, DIMENSION(N), INTENT(INOUT) :: D
		INTEGER, INTENT(INOUT) :: NROT
		DOUBLE PRECISION, POINTER :: B(:), Z(:)
		DOUBLE PRECISION c, g, h, s, sm, t, tau, theta, tresh

		INTEGER ialloc, ip, iq, i, j

		ALLOCATE(B(100), stat = ialloc)
		ALLOCATE(Z(100), stat = ialloc)

		DO ip = 1,  N    !initialize V to identity matrix
			DO iq = 1,  N
				V(ip, iq) = 0.d0
			END DO
			V(ip, ip) = 1.d0
		END DO
	  DO ip = 1,  N
	    B(ip) = A(ip, ip)
	    D(ip) = B(ip)
	    Z(ip) = 0.d0
	  END DO
	  NROT = 0
	  DO i = 1,  50
	    sm = 0.d0
	    DO ip = 1,  N-1     !sum off-diagonal elements
	      DO iq = ip + 1,  N
	        sm = sm + DABS(A(ip, iq))
	      END DO
	    END DO
	    IF(sm .EQ. 0.d0) RETURN  !normal return
	    IF(i .LT. 4) THEN
	      tresh = 0.2d0 * sm ** 2
	    ELSE
	      tresh = 0.d0
	    END IF
	    DO ip = 1,  N-1
	      DO iq = ip + 1,  N
	        g = 100.d0 * DABS(A(ip, iq))
	! after 4 sweeps,  skip the rotation IF the off-diagonal element is small
	        IF((i .GT. 4) .AND. (DABS(D(ip)) + g .EQ. DABS(D(ip))) &
			.AND. (DABS(D(iq)) + g .EQ. DABS(D(iq)))) THEN
			  A(ip, iq) = 0.d0
	        ELSE IF(DABS(A(ip, iq)) .GT. tresh) THEN
		  h = D(iq)-D(ip)
		  IF ((DABS(h) + g) .EQ. DABS(h)) THEN
		    t = A(ip, iq)/h
	          ELSE
		    theta = 0.5d0*h/A(ip, iq)
	            t = 1.d0/(DABS(theta) + DSQRT(1.d0 + theta**2))
		    IF (theta .LT. 0.d0) t = -t
	          END IF
		  c = 1.d0/DSQRT(1.d0 + t**2)
		  s = t*c
	          tau = s/(1.d0 + c)
		  h = t*A(ip, iq)
		  Z(ip) = Z(ip)-h
		  Z(iq) = Z(iq) + h
		  D(ip) = D(ip)-h
		  D(iq) = D(iq) + h
		  A(ip, iq) = 0.d0
		  DO j = 1,  ip-1
		    g = A(j, ip)
		    h = A(j, iq)
		    A(j, ip) = g-s*(h + g*tau)
		    A(j, iq) = h + s*(g-h*tau)
	          END DO
		  DO j = ip + 1,  iq-1
		    g = A(ip, j)
		    h = A(j, iq)
		    A(ip, j) = g-s*(h + g*tau)
		    A(j, iq) = h + s*(g-h*tau)
	          END DO
		  DO j = iq + 1,  N
		    g = A(ip, j)
		    h = A(iq, j)
		    A(ip, j) = g-s*(h + g*tau)
		    A(iq, j) = h + s*(g-h*tau)
	          END DO
		  DO j = 1,  N
		    g = V(j, ip)
		    h = V(j, iq)
		    V(j, ip) = g-s*(h + g*tau)
		    V(j, iq) = h + s*(g-h*tau)
	          END DO
	          NROT = NROT + 1
	        END IF !IF ((i.gt.4)...
	      END DO !main iq loop
	    END DO !main ip loop
	    DO ip = 1,  N
	      B(ip) = B(ip) + Z(ip)
	      D(ip) = B(ip)
	      Z(ip) = 0.d0
	    END DO
	  END DO !main i loop
		WRITE(*,'('' 50 iterations !'')')
		RETURN
	END SUBROUTINE Jacobi_real8

END MODULE jacobi
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

        	b = 1.0D0/boscil
        	order = n - 1

        	ArgGamma = order + l + 1
        	NormSquared = 2.0d0*b**3*EXP(DDLogFactorials(order)-GammaFunction(ArgGamma))
        	Norm = SQRT(NormSquared)

        	chi = b*x
        	alpha = DBLE(l) + 0.5d0

        	argPol = chi**2
        	CALL GeneralizedLaguerre(Poly, argPol, order, alpha)

        	Rwave = Norm*EXP(-0.5D0*chi**2)*chi**l*Poly

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
        	        b(i) = + (2.0d0*DBLE(i) + alpha + 1.0d0 - x)/DBLE(i + 1)
        	        c(i) = - (DBLE(i) + alpha)/(DBLE(i + 1))
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

 MODULE lgfactor

	IMPLICIT NONE

#if(USE_QUADRUPLE==1)
	REAL(KIND = 16), DIMENSION(:), ALLOCATABLE :: ddlnf, ddlng
#else
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ddlnf, ddlng
#endif
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lnf, lng, GammaFunc

	INTEGER lnf_max, lng_max, ddlnf_max, ddlng_max

	PRIVATE lnf, lng, ddlnf, ddlng
	PRIVATE lnf_max, lng_max, ddlnf_max, ddlng_max

 CONTAINS

	SUBROUTINE LogFactorials_new(imax)
		INTEGER, INTENT(IN) :: imax

		INTEGER i

		ALLOCATE(lnf(0:imax))
		lnf_max = imax
		lnf(0) = 0.0D0
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
		lng(0) = -0.1207822376352452223455184D0
		DO i = 1, imax
			lng(i) = LOG(DBLE(i) + 0.5D0) + lng(i - 1)
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
		ddlnf(0) = 0.0D0
		DO i = 1, imax
			ddlnf(i) = LOG(DBLE(i)) + ddlnf(i - 1)
		END DO
		RETURN
	END SUBROUTINE DDLogFactorials_new

	FUNCTION DDLogFactorials(i)
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) DDLogFactorials
#else
		DOUBLE PRECISION DDLogFactorials
#endif
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

#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) pi
#else
		DOUBLE PRECISION pi
#endif
		INTEGER i

		ddlng_max = imax

		ALLOCATE(ddlng(0:imax))

		! For n=0, we have (3/2)! = SQRT(PI) / 2 (Abramovitz, 6.1.9)
		! We calculate ln ( (3/2)! ) at machine precision

		pi = 4.0D0*ATAN(1.0D0)
		ddlng(0) = LOG(0.5D0*SQRT(pi))

		DO i = 1, imax
			ddlng(i) = LOG(DBLE(i) + 0.5D0) + ddlng(i - 1)
		END DO

		RETURN
	END SUBROUTINE DDLogSemiFactorials_new

	FUNCTION DDLogSemiFactorials(i)
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) DDLogSemiFactorials
		REAL(KIND = 16) pi
#else
		DOUBLE PRECISION DDLogSemiFactorials
		DOUBLE PRECISION pi
#endif
		INTEGER, INTENT(IN) :: i

		! For n=-1, we have (-1/2)! = SQRT(PI) (Abramovitz, 6.1.8)
		! We calculate ln ( (1/2)! ) at machine precision

		IF (i .EQ. -1) THEN
			pi = 4.0D0*ATAN(1.0D0)
			DDLogSemiFactorials = LOG(SQRT(pi))
		ELSE IF ((i .GE. 0) .AND. (i .LE. ddlng_max)) THEN
			DDLogSemiFactorials = ddlng(i)
		ELSE
			WRITE(*,'("Beyond range in DDLogSemiFactorials - n = ",i8," nmax = ",I8)')  i, ddlng_max
			STOP "DDLogSemiFactorials"
		END IF
		RETURN
	END FUNCTION DDLogSemiFactorials

END MODULE lgfactor
PROGRAM main

	USE input    ! Parametros de entrada
	USE global
	USE nucleus  ! Nucleo
	USE symden
	USE diagmeth ! Metodo de diagonalizacion

	IMPLICIT NONE

        TYPE (SymDensity) :: density
	TYPE (DiagonalizationMethod) :: diagonal

	INTEGER :: cycles_in, cycles_out, cycles_rate

	! Initialize clock
	CALL SYSTEM_CLOCK(cycles_in)

	! Define important quantities and computes one-body kinetic energy and
	! Gauss-Laguerre quadratures (global.f90)
	CALL Global_new

	! Define a new set of densities and a new nucleus
	CALL SymDensity_new(density, neutrons, protons)
	CALL Nucleus_set_b(density%nucleus, b_0)

	! Self-consistent loop
	CALL DiagonalizationMethod_new(diagonal, density)
	CALL DiagonalizationMethod_goto_SelfConsistency(diagonal, Convergence)

	! Nothing happens here
	CALL Global_del

	! Stops the clock and print times
	CALL SYSTEM_CLOCK(cycles_out, cycles_rate)
	WRITE(6,'("Elapsed Time")')
	WRITE(6,'("Seconds .............:",F20.6)') DBLE(cycles_out - cycles_in) / cycles_rate
	WRITE(6,'("Clock cycles ........:",I12)') cycles_out - cycles_in

	STOP

END PROGRAM main
	DOUBLE PRECISION, DIMENSION(0:1, 0:1658) :: Nucleus_MassExcess
	DATA Nucleus_MassExcess / &
		 008071.380, 0000.008,  007289.029, 0000.007, &
		 013135.824, 0000.014,  014949.913, 0000.021, &
		 014931.314, 0000.021,  025840.000, 0380.000, &
		 002424.910, 0000.030,  025120.000, 0300.000, &
		 011390.000, 0050.000,  011680.000, 0050.000, &
		 017592.600, 0000.900,  014085.700, 0000.700, &
		 018374.000, 0005.000,  026110.000, 0030.000, &
		 014907.000, 0000.800,  015768.900, 0000.800, &
		 027870.000, 0070.000,  031598.000, 0007.000, &
		 020945.600, 0000.800,  004941.710, 0000.080, &
		 022920.400, 0001.200,  035094.000, 0024.000, &
		 040810.000, 0120.000,  024954.100, 0002.000, &
		 011347.700, 0000.400,  012415.900, 0001.100, &
		 028913.900, 0002.300,  033840.000, 0250.000, &
		 012607.100, 0000.400,  012050.990, 0000.270, &
		 015699.100, 0000.300,  040900.000, 0110.000, &
		 020174.000, 0006.000,  008668.200, 0000.300, &
		 010650.400, 0000.900,  024890.000, 0140.000, &
		 025077.000, 0015.000,  013369.500, 0001.300, &
		 000000.000, 0000.000,  017338.100, 0001.000, &
		 032060.000, 0040.000,  035000.000, 0500.000, &
		 016562.500, 0001.100,  003125.032, 0000.010, &
		 005345.520, 0000.270,  023113.000, 0010.000, &
		 040100.000, 0130.000,  023664.000, 0021.000, &
		 003019.908, 0000.016,  002863.433, 0000.015, &
		 008006.540, 0000.080,  028970.000, 0022.000, &
		 009873.200, 0000.800,  000101.496, 0000.024, &
		 002855.400, 0000.500,  016770.000, 0130.000, &
		 013694.000, 0004.000,  005682.100, 0002.300, &
		-004737.037, 0000.030,  010680.000, 0008.000, &
		 023989.000, 0020.000,  021035.000, 0017.000, &
		 007871.000, 0015.000, -000809.080, 0000.260, &
		 001951.780, 0000.250,  016480.000, 0050.000, &
		 024920.000, 0030.000,  013117.000, 0020.000, &
		-000782.200, 0000.800,  000873.400, 0000.600, &
		 005319.000, 0005.000,  015871.000, 0019.000, &
		 003332.100, 0002.900, -001487.430, 0000.090, &
		 001751.000, 0000.600,  012928.000, 0012.000, &
		 003796.900, 0001.200, -000017.350, 0000.170, &
		-007047.800, 0001.800,  006839.000, 0007.000, &
		 017570.000, 0027.000,  008066.000, 0015.000, &
		-000047.500, 0001.800, -005737.400, 0001.700, &
		-002189.800, 0002.000,  010913.000, 0016.000, &
		 009440.000, 0090.000,  002830.000, 0030.000, &
		-008027.200, 0001.500, -005185.200, 0001.600, &
		-000397.000, 0001.500,  018090.000, 0070.000, &
		 003350.000, 0170.000, -005156.000, 0002.500, &
		-009532.300, 0000.800, -005473.700, 0001.400, &
		 006767.000, 0025.000, -005950.000, 0010.000, &
		-008420.400, 0000.800, -013933.500, 0000.600, &
		-000055.000, 0004.000,  010755.000, 0019.000, &
		-002060.000, 0040.000, -009360.300, 0001.400, &
		-013192.800, 0000.600, -008915.800, 0000.900, &
		 003827.000, 0010.000,  000440.000, 0070.000, &
		-006904.000, 0016.000, -016214.100, 0000.600, &
		-012210.400, 0000.600, -007145.000, 0003.000, &
		-005600.000, 0040.000, -014586.300, 0000.800, &
		-017197.200, 0000.600, -012385.400, 0000.600, &
		-000750.000, 0040.000, -001140.000, 0140.000, &
		-015019.200, 0002.100, -016851.000, 0000.600, &
		-021492.900, 0000.500, -007161.000, 0004.000, &
		 004130.000, 0160.000,  002650.000, 0150.000, &
		-010661.000, 0029.000, -018215.800, 0001.300, &
		-021895.400, 0000.500, -016951.900, 0000.900, &
		-003160.000, 0050.000,  008210.000, 0250.000, &
		-009100.000, 0210.000, -015890.000, 0040.000, &
		-024433.600, 0000.500, -020200.900, 0000.700, &
		-014063.000, 0003.000,  011830.000, 0580.000, &
		-015050.000, 0070.000, -022950.600, 0000.600, &
		-024441.200, 0000.500, -019045.400, 0001.500, &
		-007060.000, 0050.000,  016550.000, 0740.000, &
		-001750.000, 1580.000, -024081.100, 0001.000, &
		-024305.500, 0000.500, -026016.370, 0000.160, &
		-013330.000, 0008.000, -002180.000, 0050.000, &
		 021470.000, 1140.000, -020492.000, 0016.000, &
		-026338.100, 0001.100, -026586.630, 0000.140, &
		-021003.800, 0000.500, -009380.000, 0030.000, &
		 026650.000, 3570.000, -019958.000, 0015.000, &
		-024557.700, 0001.200, -029932.370, 0000.130, &
		-024440.070, 0000.180, -018379.000, 0003.000, &
		-014320.000, 0050.000, -024857.900, 0001.900, &
		-028846.910, 0000.130, -029013.740, 0000.040, &
		-023048.900, 0001.300, -011167.000, 0020.000, &
		 004450.000, 0060.000, -020251.000, 0013.000, &
		-030664.250, 0000.240, -029522.100, 0000.080, &
		-030230.710, 0000.260, -017425.000, 0008.000, &
		-006440.000, 0040.000, -026896.480, 0000.260, &
		-031761.780, 0000.070, -030948.700, 0000.500, &
		-024798.700, 0001.400, -013159.000, 0022.000, &
		-026861.000, 0007.000, -029798.260, 0000.120, &
		-034715.100, 0000.800, -028802.000, 0000.900, &
		-022059.000, 0005.000, -029802.800, 0002.500, &
		-033241.000, 0005.000, -033806.300, 0001.000, &
		-027275.600, 0002.000, -014300.000, 0050.000, &
		-022520.000, 0040.000, -027530.000, 0040.000, &
		-035039.100, 0001.100, -033534.500, 0001.000, &
		-034846.100, 0001.000, -020526.000, 0004.000, &
		-009063.000, 0011.000, -027400.000, 0150.000, &
		-033066.500, 0001.200, -035558.400, 0001.000, &
		-035137.500, 0001.000, -028643.000, 0001.400, &
		-015690.000, 0013.000, -034420.000, 0040.000, &
		-035020.900, 0001.000, -038547.000, 0001.100, &
		-032121.900, 0001.300, -025121.000, 0006.000, &
		-023130.000, 0060.000, -031980.000, 0070.000, &
		-036593.000, 0009.000, -038408.600, 0001.100, &
		-036187.800, 0002.200, -029320.000, 0007.000, &
		-032260.000, 0020.000, -035810.000, 0040.000, &
		-041469.200, 0001.200, -037815.900, 0002.100, &
		-037548.200, 0001.200, -013450.000, 0030.000, &
		-029720.000, 0060.000, -036614.000, 0010.000, &
		-040812.800, 0001.200, -041069.500, 0001.300, &
		-039006.400, 0001.400, -031875.000, 0017.000, &
		-019410.000, 0100.000, -029720.000, 0040.000, &
		-035418.000, 0016.000, -043140.700, 0002.600, &
		-041758.800, 0001.300, -044125.400, 0001.200, &
		-037075.000, 0001.400, -029472.000, 0020.000, &
		-025910.000, 0100.000, -035696.000, 0008.000, &
		-042345.500, 0002.500, -044330.400, 0002.100, &
		-044931.800, 0001.000, -042004.000, 0001.200, &
		-034553.000, 0014.000, -032122.000, 0024.000, &
		-044214.000, 0004.000, -044492.000, 0005.000, &
		-048487.100, 0001.000, -044474.700, 0002.600, &
		-042818.000, 0007.000, -029211.000, 0021.000, &
		-018130.000, 0110.000, -030770.000, 0300.000, &
		-041289.000, 0004.000, -046558.000, 0004.000, &
		-048558.100, 0001.000, -047956.200, 0001.300, &
		-045328.300, 0002.700, -037611.000, 0024.000, &
		-024580.000, 0160.000, -039570.000, 0009.000, &
		-044537.000, 0016.000, -051426.000, 0001.000, &
		-049219.000, 0001.500, -050257.300, 0001.500, &
		-042625.400, 0001.500, -034470.000, 0060.000, &
		-035010.000, 0080.000, -043218.000, 0020.000, &
		-049726.900, 0001.400, -052199.000, 0001.500, &
		-051447.700, 0001.500, -048238.900, 0001.500, &
		-040217.000, 0015.000, -049464.000, 0007.000, &
		-051438.800, 0001.600, -055414.400, 0001.500, &
		-050702.900, 0002.400, -048331.000, 0010.000, &
		-034287.000, 0022.000, -022640.000, 0070.000, &
		-046830.000, 0100.000, -051846.000, 0003.000, &
		-055282.500, 0001.500, -054686.800, 0001.600, &
		-050943.100, 0002.100, -042639.000, 0018.000, &
		-029380.000, 0160.000, -049889.000, 0015.000, &
		-056930.100, 0001.500, -055553.100, 0001.800, &
		-056250.300, 0001.400, -048007.900, 0001.400, &
		-039210.000, 0050.000, -049150.000, 0100.000, &
		-055105.100, 0001.600, -057708.300, 0001.500, &
		-057476.900, 0001.400, -054025.600, 0001.400, &
		-045330.000, 0011.000, -055290.000, 0010.000, &
		-056907.500, 0001.500, -060603.700, 0001.500, &
		-056037.700, 0002.500, -053901.000, 0011.000, &
		-038584.000, 0017.000, -026130.000, 0080.000, &
		-057487.000, 0003.000, -060178.500, 0001.500, &
		-059342.600, 0001.500, -056077.400, 0003.000, &
		-047350.000, 0050.000, -032700.000, 0120.000, &
		-055830.000, 0030.000, -062151.800, 0001.500, &
		-059844.100, 0001.800, -060225.000, 0001.500, &
		-051662.300, 0002.500, -042210.000, 0100.000, &
		-055476.000, 0029.000, -060661.400, 0001.500, &
		-062226.200, 0001.500, -061153.600, 0001.500, &
		-056353.500, 0001.700, -047260.000, 0040.000, &
		-052900.000, 0100.000, -061406.000, 0004.000, &
		-061646.800, 0001.500, -064470.800, 0001.500, &
		-058344.000, 0002.600, -054185.000, 0011.000, &
		-058919.000, 0020.000, -062897.100, 0001.700, &
		-064219.600, 0001.500, -061982.000, 0001.800, &
		-056343.000, 0016.000, -058896.000, 0015.000, &
		-061423.000, 0019.000, -066745.500, 0001.500, &
		-062797.000, 0004.000, -061170.000, 0010.000, &
		-051999.000, 0028.000, -055190.000, 0060.000, &
		-061839.000, 0020.000, -065512.800, 0001.500, &
		-065578.800, 0001.500, -062211.900, 0002.200, &
		-056690.000, 0100.000, -059791.000, 0020.000, &
		-067098.000, 0001.500, -065423.600, 0001.500, &
		-066002.200, 0001.700, -058837.000, 0004.000, &
		-054430.000, 0250.000, -059160.000, 0050.000, &
		-065124.800, 0001.600, -067262.100, 0001.800, &
		-065910.400, 0001.800, -062654.900, 0002.000, &
		-056410.000, 0100.000, -066029.000, 0016.000, &
		-066256.700, 0001.800, -068898.800, 0001.600, &
		-063724.000, 0003.000, -061620.000, 0030.000, &
		-052070.000, 0060.000, -063743.000, 0019.000, &
		-067303.000, 0008.000, -067879.600, 0001.600, &
		-066878.600, 0001.700, -062656.000, 0005.000, &
		-056650.000, 0100.000, -063483.000, 0017.000, &
		-065540.000, 0050.000, -070006.500, 0001.600, &
		-067085.400, 0002.000, -066978.000, 0006.000, &
		-058880.000, 0100.000, -060460.000, 0150.000, &
		-065741.000, 0008.000, -068417.300, 0001.700, &
		-069322.600, 0002.900, -067097.000, 0004.000, &
		-063080.000, 0030.000, -056300.000, 0030.000, &
		-063390.000, 0110.000, -069561.000, 0003.000, &
		-068905.900, 0003.000, -070561.800, 0001.500, &
		-064340.000, 0050.000, -067323.000, 0011.000, &
		-070139.400, 0001.900, -069906.200, 0001.800, &
		-067894.000, 0004.000, -068131.000, 0006.000, &
		-068589.300, 0002.100, -072583.600, 0001.400, &
		-068228.000, 0004.000, -067897.000, 0012.000, &
		-065410.000, 0040.000, -069705.000, 0006.000, &
		-071295.200, 0001.400, -070955.000, 0004.000, &
		-068215.000, 0011.000, -063600.000, 0230.000, &
		-056890.000, 0140.000, -065708.000, 0019.000, &
		-068060.000, 0070.000, -073423.600, 0001.400, &
		-070861.400, 0002.200, -072215.300, 0001.500, &
		-065301.000, 0015.000, -062130.000, 0060.000, &
		-051670.000, 0460.000, -062530.000, 0090.000, &
		-068466.000, 0007.000, -071858.100, 0001.700, &
		-073035.400, 0001.600, -072171.500, 0001.500, &
		-069142.000, 0014.000, -064214.000, 0020.000, &
		-057210.000, 0100.000, -062290.000, 0170.000, &
		-066440.000, 0150.000, -073214.800, 0001.500, &
		-072290.500, 0001.800, -075254.400, 0001.500, &
		-070291.000, 0009.000, -068965.000, 0012.000, &
		-060530.000, 0060.000, -071216.000, 0001.900, &
		-073918.800, 0002.200, -074601.900, 0001.500, &
		-073237.000, 0003.000, -070194.000, 0017.000, &
		-064917.000, 0030.000, -057880.000, 0150.000, &
		-071863.000, 0004.000, -072819.000, 0010.000, &
		-077028.400, 0001.500, -073455.000, 0004.000, &
		-074147.000, 0008.000, -066980.000, 0030.000, &
		-062720.000, 0120.000, -069490.000, 0090.000, &
		-073639.000, 0006.000, -075920.000, 0001.600, &
		-076070.600, 0002.400, -074445.000, 0004.000, &
		-070839.000, 0023.000, -051890.000, 0360.000, &
		-059380.000, 0300.000, -069380.000, 0030.000, &
		-072165.000, 0024.000, -077762.500, 0001.800, &
		-075891.600, 0002.300, -077894.000, 0006.000, &
		-072176.000, 0018.000, -070190.000, 0030.000, &
		-057990.000, 0190.000, -066310.000, 0120.000, &
		-072536.000, 0006.000, -076392.200, 0001.900, &
		-077978.000, 0005.000, -077697.000, 0005.000, &
		-075459.000, 0021.000, -071470.000, 0040.000, &
		-065950.000, 0070.000, -058790.000, 0300.000, &
		-065380.000, 0140.000, -070078.000, 0025.000, &
		-077596.500, 0002.100, -077499.000, 0005.000, &
		-080592.000, 0005.000, -076203.000, 0017.000, &
		-075998.000, 0008.000, -068180.000, 0100.000, &
		-064180.000, 0510.000, -069880.000, 0220.000, &
		-075343.000, 0004.000, -079010.000, 0004.000, &
		-079982.000, 0003.000, -079049.000, 0021.000, &
		-076781.000, 0021.000, -072370.000, 0060.000, &
		-066350.000, 0100.000, -075952.000, 0015.000, &
		-077776.000, 0025.000, -082430.000, 0003.000, &
		-079748.000, 0003.000, -080641.000, 0004.000, &
		-074230.000, 0170.000, -072420.000, 0100.000, &
		-078607.000, 0019.000, -081477.000, 0003.000, &
		-082164.400, 0002.800, -081099.000, 0004.000, &
		-077845.000, 0025.000, -073150.000, 0100.000, &
		-070540.000, 0130.000, -075640.000, 0060.000, &
		-083262.000, 0005.000, -082744.200, 0002.700, &
		-084518.800, 0002.500, -079279.000, 0014.000, &
		-073856.000, 0025.000, -080706.000, 0005.000, &
		-084593.100, 0002.900, -084875.500, 0002.500, &
		-083014.200, 0002.800, -079348.000, 0008.000, &
		-074180.000, 0060.000, -067440.000, 0310.000, &
		-070720.000, 0130.000, -079688.000, 0014.000, &
		-082601.000, 0004.000, -087916.800, 0002.500, &
		-084294.200, 0002.900, -083626.000, 0010.000, &
		-076720.000, 0050.000, -081709.000, 0007.000, &
		-086211.000, 0004.000, -087703.000, 0002.500, &
		-084871.000, 0003.000, -080580.000, 0040.000, &
		-075005.000, 0015.000, -064650.000, 0120.000, &
		-074947.000, 0027.000, -079350.000, 0013.000, &
		-085942.800, 0002.900, -086488.900, 0002.500, &
		-088770.500, 0002.300, -082659.000, 0005.000, &
		-080170.000, 0006.000, -071370.000, 0080.000, &
		-077786.000, 0010.000, -083652.000, 0009.000, &
		-086349.000, 0003.000, -087893.700, 0002.300, &
		-086640.000, 0003.000, -082208.000, 0012.000, &
		-075990.000, 0200.000, -068650.000, 0080.000, &
		-074811.000, 0012.000, -082923.000, 0013.000, &
		-084833.000, 0010.000, -088457.400, 0002.300, &
		-086451.600, 0002.800, -086809.000, 0004.000, &
		-078939.000, 0026.000, -064160.000, 0120.000, &
		-072688.000, 0015.000, -080160.000, 0016.000, &
		-084245.000, 0011.000, -087120.200, 0002.300, &
		-087210.700, 0002.400, -086805.000, 0004.000, &
		-083607.000, 0004.000, -077270.000, 0090.000, &
		-068518.000, 0019.000, -078836.000, 0007.000, &
		-082348.000, 0006.000, -087268.300, 0002.400, &
		-086368.500, 0002.400, -088413.300, 0002.300, &
		-084158.000, 0005.000, -082569.000, 0013.000, &
		-072940.000, 0450.000, -065813.000, 0024.000, &
		-075050.000, 0050.000, -081214.000, 0006.000, &
		-085659.600, 0002.400, -086783.900, 0002.000, &
		-087709.500, 0002.000, -086018.000, 0006.000, &
		-083451.000, 0012.000, -078340.000, 0150.000, &
		-061150.000, 0030.000, -072880.000, 0040.000, &
		-078300.000, 0040.000, -085442.000, 0003.000, &
		-085606.000, 0004.000, -088792.400, 0002.000, &
		-085819.000, 0006.000, -086073.000, 0008.000, &
		-079626.000, 0013.000, -076180.000, 0150.000, &
		-058290.000, 0040.000, -068810.000, 0070.000, &
		-076270.000, 0060.000, -082950.000, 0003.000, &
		-085608.200, 0002.700, -087542.200, 0002.000, &
		-087222.000, 0005.000, -086113.000, 0008.000, &
		-082590.000, 0040.000, -077800.000, 0300.000, &
		-054090.000, 0060.000, -066380.000, 0060.000, &
		-072520.000, 0160.000, -081283.000, 0020.000, &
		-083528.000, 0006.000, -088113.400, 0002.000, &
		-086429.000, 0004.000, -088225.000, 0006.000, &
		-083168.000, 0012.000, -081301.000, 0022.000, &
		-050860.000, 0110.000, -062150.000, 0100.000, &
		-070170.000, 0080.000, -077790.000, 0040.000, &
		-082328.000, 0013.000, -085967.200, 0002.000, &
		-087324.300, 0002.000, -087617.700, 0002.100, &
		-085519.000, 0010.000, -082193.000, 0016.000, &
		-076760.000, 0150.000, -060200.000, 0210.000, &
		-067290.000, 0140.000, -076590.000, 0040.000, &
		-079929.000, 0028.000, -086186.000, 0006.000, &
		-086017.300, 0002.300, -089219.700, 0002.100, &
		-085590.000, 0020.000, -085221.000, 0013.000, &
		-078170.000, 0090.000, -073380.000, 0070.000, &
		-078950.000, 0040.000, -083513.000, 0006.000, &
		-086337.000, 0024.000, -087950.400, 0002.100, &
		-087410.000, 0017.000, -085430.000, 0018.000, &
		-081190.000, 0120.000, -075660.000, 0180.000, &
		-071770.000, 0060.000, -076350.000, 0050.000, &
		-083559.000, 0021.000, -084569.000, 0009.000, &
		-089099.400, 0002.200, -086821.000, 0017.000, &
		-087918.000, 0008.000, -082080.000, 0050.000, &
		-068290.000, 0130.000, -075240.000, 0090.000, &
		-080760.000, 0080.000, -084601.000, 0010.000, &
		-087260.400, 0002.200, -088024.300, 0003.000, &
		-087471.000, 0008.000, -084787.000, 0017.000, &
		-080650.000, 0016.000, -074607.000, 0025.000, &
		-072260.000, 0120.000, -080370.000, 0060.000, &
		-082490.000, 0050.000, -088093.000, 0004.000, &
		-086952.000, 0003.000, -089393.000, 0005.000, &
		-085114.000, 0006.000, -083977.000, 0010.000, &
		-070940.000, 0120.000, -077360.000, 0090.000, &
		-082350.000, 0060.000, -085932.000, 0004.000, &
		-087849.000, 0005.000, -088416.000, 0005.000, &
		-087078.000, 0009.000, -084339.000, 0010.000, &
		-079493.000, 0016.000, -073240.000, 0080.000, &
		-076270.000, 0090.000, -079790.000, 0060.000, &
		-086326.000, 0008.000, -086365.000, 0008.000, &
		-089907.000, 0005.000, -086941.000, 0005.000, &
		-087135.000, 0006.000, -080617.000, 0013.000, &
		-077450.000, 0060.000, -083710.000, 0300.000, &
		-086862.000, 0018.000, -088374.000, 0006.000, &
		-088407.000, 0005.000, -086990.000, 0007.000, &
		-083568.000, 0013.000, -083760.000, 0230.000, &
		-085080.000, 0200.000, -089523.000, 0004.000, &
		-087605.000, 0005.000, -089253.000, 0005.000, &
		-084112.000, 0029.000, -082050.000, 0040.000, &
		-085021.000, 0020.000, -087605.000, 0004.000, &
		-088721.000, 0003.000, -088507.000, 0004.000, &
		-086487.000, 0006.000, -082633.000, 0010.000, &
		-076253.000, 0019.000, -067620.000, 0070.000, &
		-082940.000, 0220.000, -088345.000, 0015.000, &
		-087459.000, 0003.000, -090351.000, 0003.000, &
		-086410.000, 0200.000, -085834.000, 0016.000, &
		-072300.000, 0060.000, -086030.000, 0040.000, &
		-088217.000, 0004.000, -089254.000, 0003.000, &
		-088391.000, 0006.000, -085943.000, 0007.000, &
		-073470.000, 0070.000, -086333.000, 0019.000, &
		-086624.000, 0017.000, -090581.200, 0002.900, &
		-087995.000, 0005.000, -088658.000, 0004.000, &
		-081603.000, 0023.000, -077270.000, 0170.000, &
		-083680.000, 0150.000, -087040.000, 0020.000, &
		-089050.000, 0002.900, -089368.000, 0003.000, &
		-088330.000, 0004.000, -084424.000, 0024.000, &
		-071120.000, 0050.000, -062090.000, 0090.000, &
		-083460.000, 0029.000, -084960.000, 0070.000, &
		-090021.800, 0002.900, -088571.000, 0003.000, &
		-090560.000, 0003.000, -084680.000, 0200.000, &
		-084950.000, 0070.000, -088091.300, 0002.900, &
		-089539.000, 0004.000, -090034.100, 0003.000, &
		-087004.000, 0020.000, -082360.000, 0230.000, &
		-080140.000, 0150.000, -082760.000, 0110.000, &
		-088720.000, 0003.000, -088252.000, 0005.000, &
		-091526.100, 0003.000, -086819.000, 0006.000, &
		-085290.000, 0100.000, -077550.000, 0140.000, &
		-062290.000, 0290.000, -082250.000, 0050.000, &
		-086416.000, 0013.000, -088945.000, 0005.000, &
		-090399.300, 0002.900, -088644.000, 0009.000, &
		-085110.000, 0019.000, -066260.000, 0180.000, &
		-079580.000, 0100.000, -086709.000, 0020.000, &
		-087232.000, 0008.000, -091654.200, 0002.900, &
		-087998.000, 0004.000, -087653.000, 0023.000, &
		-068270.000, 0130.000, -078590.000, 0070.000, &
		-083940.000, 0060.000, -087733.000, 0008.000, &
		-090068.200, 0002.800, -089475.000, 0008.000, &
		-087182.000, 0008.000, -083780.000, 0100.000, &
		-078750.000, 0140.000, -072240.000, 0100.000, &
		-075770.000, 0100.000, -083973.000, 0019.000, &
		-085800.000, 0170.000, -091103.800, 0002.600, &
		-088423.000, 0008.000, -089386.000, 0019.000, &
		-083771.000, 0024.000, -081810.000, 0050.000, &
		-073820.000, 0080.000, -074550.000, 0190.000, &
		-080950.000, 0150.000, -085841.000, 0028.000, &
		-089203.900, 0002.500, -089591.600, 0002.400, &
		-088551.000, 0026.000, -086270.000, 0020.000, &
		-082510.000, 0060.000, -077110.000, 0060.000, &
		-083580.000, 0050.000, -089946.000, 0002.600, &
		-088327.000, 0002.400, -090307.200, 0002.700, &
		-086073.000, 0006.000, -085050.000, 0140.000, &
		-078140.000, 0060.000, -083420.000, 0030.000, &
		-087820.500, 0002.600, -089223.800, 0002.000, &
		-089171.300, 0001.900, -087937.000, 0004.000, &
		-085258.000, 0016.000, -081070.000, 0040.000, &
		-081060.000, 0050.000, -088237.300, 0001.400, &
		-087619.900, 0002.000, -090525.100, 0001.500, &
		-087368.000, 0004.000, -087659.600, 0002.000, &
		-081740.000, 0040.000, -080420.000, 0080.000, &
		-085898.400, 0002.000, -088258.000, 0003.000, &
		-089024.800, 0002.300, -088846.100, 0002.500, &
		-087191.500, 0002.000, -084113.000, 0017.000, &
		-079550.000, 0250.000, -077810.000, 0080.000, &
		-086021.000, 0011.000, -086400.000, 0030.000, &
		-090067.100, 0002.300, -087916.000, 0005.000, &
		-089174.000, 0007.000, -084347.000, 0024.000, &
		-077010.000, 0070.000, -083504.000, 0025.000, &
		-086705.000, 0006.000, -088286.000, 0004.000, &
		-088982.000, 0004.000, -088319.000, 0005.000, &
		-086243.000, 0012.000, -082790.000, 0100.000, &
		-074020.000, 0170.000, -083330.000, 0050.000, &
		-084610.000, 0040.000, -088992.000, 0003.000, &
		-087736.000, 0004.000, -089860.800, 0001.600, &
		-085928.000, 0006.000, -085470.000, 0018.000, &
		-078820.000, 0400.000, -073020.000, 0170.000, &
		-080620.000, 0120.000, -084624.000, 0022.000, &
		-087006.000, 0004.000, -088507.000, 0004.000, &
		-088698.100, 0001.800, -087506.000, 0005.000, &
		-085080.000, 0011.000, -081360.000, 0050.000, &
		-070010.000, 0200.000, -080130.000, 0080.000, &
		-082330.000, 0070.000, -087348.000, 0004.000, &
		-086897.000, 0010.000, -089881.500, 0001.500, &
		-086853.000, 0008.000, -087291.000, 0007.000, &
		-068490.000, 0140.000, -077380.000, 0070.000, &
		-082020.000, 0070.000, -085206.000, 0004.000, &
		-087457.000, 0004.000, -088428.000, 0004.000, &
		-088076.000, 0006.000, -086714.000, 0007.000, &
		-083750.000, 0100.000, -079730.000, 0410.000, &
		-076610.000, 0080.000, -079730.000, 0080.000, &
		-085222.000, 0012.000, -085715.000, 0011.000, &
		-089292.000, 0004.000, -087171.000, 0005.000, &
		-088447.000, 0008.000, -083740.000, 0050.000, &
		-071190.000, 0220.000, -079020.000, 0210.000, &
		-082970.000, 0080.000, -085888.000, 0026.000, &
		-087659.000, 0005.000, -088086.000, 0005.000, &
		-087570.000, 0005.000, -074020.000, 0160.000, &
		-082430.000, 0110.000, -083990.000, 0060.000, &
		-088125.000, 0007.000, -086906.000, 0005.000, &
		-088965.000, 0005.000, -085252.000, 0026.000, &
		-084750.000, 0200.000, -077870.000, 0090.000, &
		-083821.000, 0023.000, -086506.000, 0011.000, &
		-087662.000, 0007.000, -087867.000, 0005.000, &
		-086667.000, 0011.000, -084641.000, 0012.000, &
		-080920.000, 0150.000, -074460.000, 0050.000, &
		-079550.000, 0040.000, -086429.000, 0007.000, &
		-086354.000, 0005.000, -088903.000, 0005.000, &
		-086030.000, 0070.000, -086500.000, 0050.000, &
		-081370.000, 0050.000, -079160.000, 0060.000, &
		-069480.000, 0300.000, -076507.000, 0029.000, &
		-082383.000, 0007.000, -086556.000, 0005.000, &
		-087732.000, 0005.000, -087130.000, 0050.000, &
		-085910.000, 0050.000, -083200.000, 0050.000, &
		-079700.000, 0070.000, -074020.000, 0150.000, &
		-072290.000, 0080.000, -080110.000, 0040.000, &
		-082896.000, 0022.000, -088272.000, 0005.000, &
		-086531.000, 0005.000, -087574.000, 0011.000, &
		-083137.000, 0015.000, -068880.000, 0120.000, &
		-075690.000, 0060.000, -080710.000, 0007.000, &
		-084924.000, 0005.000, -087238.000, 0004.000, &
		-086973.000, 0013.000, -084844.000, 0013.000, &
		-082060.000, 0040.000, -077540.000, 0050.000, &
		-072080.000, 0120.000, -072990.000, 0060.000, &
		-077053.000, 0016.000, -083273.000, 0012.000, &
		-084327.000, 0004.000, -088088.000, 0004.000, &
		-084700.000, 0007.000, -084471.000, 0020.000, &
		-078380.000, 0040.000, -068320.000, 0090.000, &
		-074472.000, 0016.000, -079732.000, 0021.000, &
		-082983.000, 0025.000, -085445.000, 0004.000, &
		-086026.000, 0003.000, -084203.000, 0004.000, &
		-080472.000, 0029.000, -075943.000, 0013.000, &
		-069980.000, 0040.000, -065500.000, 0100.000, &
		-070538.000, 0020.000, -077847.000, 0020.000, &
		-080027.000, 0007.000, -084542.000, 0004.000, &
		-083798.000, 0003.000, -085960.000, 0003.000, &
		-081090.000, 0040.000, -078986.000, 0015.000, &
		-071590.000, 0100.000, -067745.000, 0028.000, &
		-073979.000, 0028.000, -078200.000, 0017.000, &
		-081616.000, 0004.000, -083078.000, 0003.000, &
		-084012.000, 0003.000, -082970.000, 0004.000, &
		-079526.000, 0005.000, -074380.000, 0030.000, &
		-063370.000, 0040.000, -071840.000, 0050.000, &
		-074940.000, 0060.000, -080441.000, 0004.000, &
		-080760.000, 0004.000, -083758.000, 0003.000, &
		-081425.000, 0004.000, -081975.000, 0004.000, &
		-075646.000, 0022.000, -060210.000, 0050.000, &
		-068120.000, 0070.000, -073020.000, 0070.000, &
		-077110.000, 0040.000, -079636.000, 0008.000, &
		-081442.000, 0003.000, -081278.000, 0004.000, &
		-080660.000, 0004.000, -078000.000, 0005.000, &
		-072950.000, 0040.000, -055700.000, 0080.000, &
		-065060.000, 0080.000, -069200.000, 0070.000, &
		-075730.000, 0070.000, -076760.000, 0060.000, &
		-080935.000, 0003.000, -079458.000, 0006.000, &
		-081000.000, 0005.000, -077125.000, 0008.000, &
		-076099.000, 0011.000, -067860.000, 0150.000, &
		-052300.000, 0110.000, -061500.000, 0100.000, &
		-067250.000, 0080.000, -072190.000, 0060.000, &
		-075470.000, 0040.000, -078156.000, 0003.000, &
		-079052.000, 0003.000, -079276.000, 0003.000, &
		-077555.000, 0004.000, -075367.000, 0004.000, &
		-070880.000, 0060.000, -064330.000, 0110.000, &
		-047580.000, 0360.000, -063810.000, 0150.000, &
		-070430.000, 0120.000, -072490.000, 0220.000, &
		-077418.000, 0004.000, -076874.000, 0011.000, &
		-079346.000, 0003.000, -076239.000, 0018.000, &
		-076278.000, 0004.000, -070680.000, 0050.000, &
		-068000.000, 0060.000, -066800.000, 0080.000, &
		-070988.000, 0011.000, -074385.000, 0004.000, &
		-076073.000, 0005.000, -077146.000, 0003.000, &
		-076455.000, 0005.000, -075135.000, 0005.000, &
		-071499.000, 0005.000, -054950.000, 0900.000, &
		-064990.000, 0120.000, -068000.000, 0080.000, &
		-073693.000, 0004.000, -073606.000, 0020.000, &
		-077060.000, 0003.000, -074800.000, 0007.000, &
		-075771.000, 0007.000, -071113.000, 0009.000, &
		-069324.000, 0012.000, -062210.000, 0160.000, &
		-070956.000, 0004.000, -073398.000, 0006.000, &
		-074587.000, 0003.000, -074663.000, 0003.000, &
		-074199.000, 0004.000, -071633.000, 0005.000, &
		-068764.000, 0005.000, -063720.000, 0060.000, &
		-070160.000, 0030.000, -071270.000, 0070.000, &
		-074773.000, 0003.000, -072899.000, 0003.000, &
		-074718.000, 0003.000, -070770.000, 0070.000, &
		-070127.000, 0006.000, -063750.000, 0060.000, &
		-060640.000, 0060.000, -070669.000, 0016.000, &
		-072569.000, 0003.000, -073378.000, 0003.000, &
		-072893.000, 0003.000, -071322.000, 0005.000, &
		-069152.000, 0005.000, -065023.000, 0007.000, &
		-068410.000, 0110.000, -072465.000, 0003.000, &
		-071748.000, 0003.000, -073717.000, 0003.000, &
		-070150.000, 0050.000, -070399.000, 0009.000, &
		-064647.000, 0009.000, -062622.000, 0012.000, &
		-054700.000, 0160.000, -070201.000, 0003.000, &
		-071829.000, 0003.000, -072081.000, 0003.000, &
		-071261.000, 0012.000, -069166.000, 0012.000, &
		-066064.000, 0023.000, -062220.000, 0050.000, &
		-056730.000, 0060.000, -069374.000, 0010.000, &
		-070096.000, 0007.000, -072546.000, 0003.000, &
		-070102.000, 0005.000, -070536.000, 0007.000, &
		-056980.000, 0070.000, -053410.000, 0060.000, &
		-066870.000, 0200.000, -069472.000, 0006.000, &
		-070834.000, 0003.000, -070772.000, 0003.000, &
		-069434.000, 0007.000, -066890.000, 0050.000, &
		-063420.000, 0090.000, -067220.000, 0080.000, &
		-070701.000, 0003.000, -069480.000, 0003.000, &
		-070418.000, 0004.000, -066200.000, 0030.000, &
		-056022.000, 0016.000, -047490.000, 0170.000, &
		-066058.000, 0009.000, -068572.000, 0003.000, &
		-069542.000, 0003.000, -069176.000, 0003.000, &
		-067338.000, 0004.000, -064570.000, 0005.000, &
		-049770.000, 0080.000, -067953.000, 0003.000, &
		-067846.000, 0003.000, -069682.000, 0003.000, &
		-066391.000, 0011.000, -066063.000, 0028.000, &
		-060460.000, 0300.000, -046080.000, 0060.000, &
		-065517.000, 0003.000, -067471.000, 0003.000, &
		-068064.000, 0003.000, -067207.000, 0004.000, &
		-065203.000, 0010.000, -062100.000, 0200.000, &
		-064240.000, 0110.000, -065680.000, 0040.000, &
		-068189.000, 0003.000, -066050.000, 0005.000, &
		-066346.000, 0004.000, -061550.000, 0040.000, &
		-049178.000, 0017.000, -040060.000, 0170.000, &
		-064700.000, 0040.000, -066389.000, 0003.000, &
		-066386.000, 0003.000, -065177.000, 0006.000, &
		-062738.000, 0006.000, -059370.000, 0100.000, &
		-054770.000, 0220.000, -042600.000, 0090.000, &
		-062090.000, 0100.000, -065976.000, 0003.000, &
		-064990.000, 0003.000, -065952.000, 0004.000, &
		-061990.000, 0020.000, -038380.000, 0060.000, &
		-063621.000, 0003.000, -064907.000, 0003.000, &
		-064530.000, 0004.000, -062938.000, 0004.000, &
		-060175.000, 0020.000, -056260.000, 0080.000, &
		-062593.000, 0003.000, -063079.000, 0003.000, &
		-064933.000, 0003.000, -061894.000, 0012.000, &
		-061589.000, 0008.000, -056110.000, 0160.000, &
		-041898.000, 0018.000, -032130.000, 0180.000, &
		-059940.000, 0060.000, -062291.000, 0006.000, &
		-063298.000, 0003.000, -062550.000, 0004.000, &
		-060596.000, 0005.000, -057470.000, 0100.000, &
		-034910.000, 0100.000, -060260.000, 0100.000, &
		-062998.000, 0003.000, -061319.000, 0004.000, &
		-061575.000, 0004.000, -057090.000, 0060.000, &
		-030130.000, 0060.000, -058805.000, 0020.000, &
		-060930.000, 0003.000, -061280.000, 0003.000, &
		-060371.000, 0004.000, -058078.000, 0005.000, &
		-054810.000, 0080.000, -056250.000, 0050.000, &
		-060117.000, 0003.000, -059802.000, 0003.000, &
		-060770.000, 0003.000, -057311.000, 0019.000, &
		-033933.000, 0018.000, -023530.000, 0190.000, &
		-057727.000, 0003.000, -059217.000, 0003.000, &
		-059314.000, 0003.000, -057834.000, 0004.000, &
		-026420.000, 0110.000, -056491.000, 0005.000, &
		-057382.000, 0006.000, -059262.000, 0003.000, &
		-056741.000, 0004.000, -056390.000, 0050.000, &
		-051470.000, 0190.000, -021240.000, 0060.000, &
		-056265.000, 0005.000, -057558.000, 0003.000, &
		-056886.000, 0004.000, -053870.000, 0040.000, &
		-056951.000, 0003.000, -055575.000, 0004.000, &
		-055851.000, 0003.000, -025324.000, 0019.000, &
		-014330.000, 0190.000, -052300.000, 0050.000, &
		-054702.000, 0003.000, -055171.900, 0002.900, &
		-054488.000, 0003.000, -017210.000, 0110.000, &
		-053501.000, 0003.000, -053394.100, 0002.800, &
		-054582.600, 0003.000, -051470.000, 0100.000, &
		-011890.000, 0060.000, -050996.000, 0003.000, &
		-052394.200, 0002.800, -052892.200, 0002.800, &
		-051726.000, 0004.000, -049705.000, 0010.000, &
		-050338.000, 0024.000, -052446.500, 0002.700, &
		-050530.000, 0100.000, -050440.000, 0100.000, &
		-045780.000, 0210.000, -016321.000, 0019.000, &
		-049110.000, 0040.000, -050475.100, 0002.700, &
		-050365.000, 0006.000, -049306.000, 0016.000, &
		-046620.000, 0050.000, -046690.000, 0070.000, &
		-049791.900, 0002.700, -048939.000, 0003.000, &
		-049647.000, 0005.000, -045840.000, 0030.000, &
		-047416.200, 0002.800, -048444.000, 0003.000, &
		-048256.000, 0005.000, -046062.000, 0007.000, &
		-046436.000, 0003.000, -048250.000, 0003.000, &
		-045450.000, 0100.000, -044542.000, 0025.000, &
		-006874.000, 0028.000, -043290.000, 0030.000, &
		-045299.000, 0003.000, -046369.100, 0003.000, &
		-045813.000, 0009.000, -041500.000, 0040.000, &
		-042844.000, 0026.000, -045709.500, 0003.000, &
		-044220.000, 0005.000, -044259.000, 0003.000, &
		-039540.000, 0250.000, -041402.000, 0014.000, &
		-043393.000, 0003.000, -043826.000, 0003.000, &
		-042813.000, 0003.000, -038620.000, 0060.000, &
		-042515.000, 0003.000, -041933.000, 0003.000, &
		-043003.000, 0003.000, -039172.000, 0020.000, &
		-037790.000, 0030.000, -039910.000, 0003.000, &
		-041222.000, 0003.000, -041224.000, 0003.000, &
		-038673.000, 0004.000, -039022.000, 0003.000, &
		-041142.000, 0003.000, -038333.000, 0007.000, &
		-037827.000, 0006.000, -035480.000, 0200.000, &
		-037985.000, 0010.000, -038993.000, 0003.000, &
		-038462.000, 0013.000, -036491.000, 0011.000, &
		-034310.000, 0160.000, -035580.000, 0150.000, &
		-038714.000, 0003.000, -036710.000, 0200.000, &
		-037331.000, 0006.000, -032889.000, 0016.000, &
		-034360.000, 0011.000, -036401.000, 0003.000, &
		-036715.000, 0004.000, -035701.000, 0006.000, &
		-033870.000, 0050.000, -030690.000, 0090.000, &
		-035892.000, 0004.000, -034843.000, 0004.000, &
		-036303.000, 0004.000, -032787.000, 0016.000, &
		-033405.000, 0004.000, -034544.000, 0004.000, &
		-034487.000, 0004.000, -027450.000, 0200.000, &
		-032442.000, 0004.000, -032539.000, 0004.000, &
		-034787.000, 0004.000, -032295.000, 0012.000, &
		-032255.000, 0023.000, -029700.000, 0500.000, &
		-031700.000, 0004.000, -032821.000, 0004.000, &
		-032594.000, 0004.000, -031070.000, 0050.000, &
		-028270.000, 0140.000, -028300.000, 0040.000, &
		-029460.000, 0040.000, -032671.000, 0004.000, &
		-031166.000, 0005.000, -031852.000, 0005.000, &
		-017970.000, 0700.000, -028292.000, 0021.000, &
		-030446.000, 0004.000, -031165.000, 0004.000, &
		-030566.000, 0005.000, -028400.000, 0050.000, &
		-019640.000, 0170.000, -029932.000, 0005.000, &
		-029606.000, 0004.000, -030979.000, 0004.000, &
		-027520.000, 0080.000, -019540.000, 0150.000, &
		-027432.000, 0019.000, -029119.000, 0004.000, &
		-029572.000, 0004.000, -028140.000, 0100.000, &
		-025270.000, 0070.000, -020920.000, 0100.000, &
		-026627.000, 0021.000, -027280.000, 0050.000, &
		-029529.000, 0004.000, -027073.000, 0007.000, &
		-020400.000, 0080.000, -008940.000, 0690.000, &
		-023750.000, 0050.000, -026413.000, 0015.000, &
		-027688.000, 0004.000, -027205.000, 0016.000, &
		-025300.000, 0030.000, -021470.000, 0050.000, &
		-010740.000, 0170.000, -024420.000, 0170.000, &
		-027370.000, 0004.000, -026006.000, 0015.000, &
		-025957.000, 0010.000, -020800.000, 0060.000, &
		-010770.000, 0150.000, -023153.000, 0016.000, &
		-025292.000, 0004.000, -025784.000, 0004.000, &
		-024810.000, 0007.000, -021580.000, 0040.000, &
		-017350.000, 0070.000, -012290.000, 0100.000, &
		-024716.000, 0004.000, -024369.000, 0004.000, &
		-025132.000, 0004.000, -020730.000, 0040.000, &
		-011900.000, 0070.000,  000650.000, 0690.000, &
		-022312.000, 0006.000, -023846.000, 0004.000, &
		-023793.000, 0004.000, -021084.000, 0008.000, &
		-017555.000, 0030.000, -013030.000, 0050.000, &
		-001270.000, 0170.000, -020969.000, 0020.000, &
		-022278.000, 0004.000, -023809.000, 0004.000, &
		-020052.000, 0009.000, -018205.000, 0010.000, &
		-012490.000, 0060.000, -001420.000, 0150.000, &
		-016270.000, 0150.000, -021049.000, 0006.000, &
		-022476.000, 0004.000, -020079.000, 0004.000, &
		-017169.000, 0007.000, -013290.000, 0040.000, &
		-008670.000, 0070.000, -002960.000, 0100.000, &
		-016774.000, 0004.000, -021772.000, 0004.000, &
		-018894.000, 0004.000, -017492.000, 0004.000, &
		-012560.000, 0040.000, -002710.000, 0070.000, &
		-013652.000, 0010.000, -017638.000, 0004.000, &
		-018282.000, 0004.000, -016390.000, 0004.000, &
		-012902.000, 0008.000, -008973.000, 0029.000, &
		-003830.000, 0050.000,  008890.000, 0170.000, &
		-009262.000, 0012.000, -014752.000, 0004.000, &
		-014815.000, 0004.000, -015977.000, 0004.000, &
		-011995.000, 0009.000, -009623.000, 0011.000, &
		-003400.000, 0050.000,  008620.000, 0160.000, &
		-010494.000, 0003.000, -011873.000, 0006.000, &
		-012457.000, 0004.000, -011674.000, 0005.000, &
		-008780.000, 0008.000, -004200.000, 0040.000, &
		 000800.000, 0070.000,  007080.000, 0110.000, &
		-007571.000, 0004.000, -008142.000, 0004.000, &
		-010394.000, 0004.000, -008640.000, 0005.000, &
		-008682.000, 0005.000, -003600.000, 0040.000, &
		 007240.000, 0080.000, -005244.000, 0008.000, &
		-006676.000, 0005.000, -006603.000, 0013.000, &
		-005722.000, 0008.000, -003572.000, 0009.000, &
		 000311.000, 0030.000,  006100.000, 0070.000, &
		-000188.100, 0002.700, -001218.000, 0012.000, &
		-004493.000, 0004.000, -003403.000, 0006.000, &
		-004343.000, 0010.000, -000983.000, 0010.000, &
		 000075.000, 0011.000,  006380.000, 0070.000, &
		 001710.000, 0090.000, -000542.900, 0003.000, &
		-001269.000, 0007.000, -001193.000, 0008.000, &
		 000292.000, 0008.000,  002509.000, 0008.000, &
		 005970.000, 0060.000,  010890.000, 0070.000, &
		 017680.000, 0120.000,  001760.000, 0004.000, &
		 002231.000, 0005.000,  000231.000, 0008.000, &
		 002960.000, 0013.000,  003269.000, 0009.000, &
		 008060.000, 0040.000,  017680.000, 0100.000, &
		 004383.000, 0008.000,  003634.000, 0005.000, &
		 004293.000, 0015.000,  005864.000, 0010.000, &
		 008685.000, 0013.000,  012160.000, 0030.000, &
		 017020.000, 0090.000,  008351.700, 0002.700, &
		 008090.000, 0013.000,  005199.000, 0004.000, &
		 007036.000, 0006.000,  006627.000, 0012.000, &
		 010820.000, 0050.000,  012348.000, 0014.000, &
		 018600.000, 0090.000,  010520.000, 0080.000, &
		 008828.300, 0003.000,  008609.000, 0008.000, &
		 009363.000, 0012.000,  011540.000, 0050.000, &
		 014450.000, 0050.000,  010590.000, 0004.000, &
		 011456.000, 0005.000,  010250.000, 0011.000, &
		 013730.000, 0050.000,  014647.000, 0022.000, &
		 013266.000, 0008.000,  012938.000, 0007.000, &
		 014500.000, 0050.000,  016917.000, 0011.000, &
		 016367.000, 0002.700,  016380.000, 0040.000, &
		 014303.000, 0006.000,  016603.000, 0006.000, &
		 017182.000, 0013.000,  021940.000, 0070.000, &
		 018381.000, 0003.000,  017232.400, 0002.900, &
		 017817.000, 0008.000,  019357.000, 0028.000, &
		 022310.000, 0070.000,  021620.000, 0050.000, &
		 018804.000, 0004.000,  020204.000, 0005.000, &
		 019980.000, 0013.000,  023780.000, 0070.000, &
		 023840.000, 0090.000,  021988.000, 0003.000, &
		 021626.000, 0008.000,  022283.000, 0009.000, &
		 024310.000, 0060.000,  027210.000, 0140.000, &
		 023662.700, 0002.700,  024303.000, 0004.000, &
		 023183.000, 0006.000,  026015.000, 0012.000, &
		 027170.000, 0030.000,  029590.000, 0090.000, &
		 027172.600, 0002.700,  025848.700, 0002.800, &
		 025803.900, 0002.900,  026824.000, 0010.000, &
		 028936.000, 0004.000,  028890.000, 0004.000, &
		 026749.000, 0004.000,  028856.000, 0006.000, &
		 029209.000, 0016.000,  032660.000, 0110.000, &
		 030900.000, 0110.000,  029581.000, 0003.000, &
		 029887.000, 0009.000,  031181.000, 0009.000, &
		 033740.000, 0080.000,  030858.700, 0002.300, &
		 032168.000, 0004.000,  031600.000, 0006.000, &
		 035220.000, 0050.000,  035910.000, 0100.000, &
		 033812.100, 0002.300,  033422.400, 0002.700, &
		 033780.000, 0050.000,  035620.000, 0050.000, &
		 035444.400, 0002.100,  035924.000, 0009.000, &
		 034587.000, 0004.000,  038349.000, 0019.000, &
		 038729.400, 0002.100,  037485.800, 0002.300, &
		 036915.000, 0003.000,  040020.000, 0050.000, &
		 040607.000, 0004.000,  040334.000, 0005.000, &
		 038141.900, 0002.100,  039952.000, 0009.000, &
		 040335.000, 0008.000,  044250.000, 0050.000, &
		 042330.000, 0050.000,  040915.500, 0002.100, &
		 041039.200, 0002.300,  042160.000, 0050.000, &
		 045340.000, 0200.000,  042441.700, 0002.100, &
		 043370.000, 0050.000,  042879.000, 0004.000, &
		 047640.000, 0100.000,  045387.300, 0002.100, &
		 044868.300, 0002.100,  045090.000, 0006.000, &
		 050910.000, 0300.000,  047305.900, 0002.100, &
		 047451.600, 0002.100,  046160.100, 0002.100, &
		 048420.000, 0050.000,  049380.000, 0040.000, &
		 050570.800, 0002.100,  049306.600, 0002.300, &
		 048584.800, 0002.100,  049385.000, 0003.000, &
		 052711.000, 0005.000,  052321.000, 0014.000, &
		 050122.500, 0002.100,  051498.000, 0013.000, &
		 051702.000, 0005.000,  054260.000, 0070.000, &
		 052952.000, 0002.100,  052931.100, 0002.100, &
		 053700.000, 0006.000,  057410.000, 0200.000, &
		 054713.800, 0002.100,  055463.400, 0002.200, &
		 054800.700, 0002.200,  059320.000, 0040.000, &
		 059922.000, 0011.000,  057751.000, 0003.000, &
		 057169.500, 0002.300,  057177.200, 0002.400, &
		 058683.000, 0005.000,  059802.000, 0005.000, &
		 059877.200, 0002.300,  058449.200, 0002.100, &
		 060700.000, 0050.000,  061460.000, 0005.000, &
		 063175.000, 0014.000,  061891.900, 0002.800, &
		 060998.000, 0002.200,  061809.600, 0002.300, &
		 063380.000, 0006.000,  065391.000, 0015.000, &
		 064990.000, 0018.000,  062614.000, 0003.000, &
		 064087.400, 0002.400,  070120.000, 0040.000, &
		 065528.000, 0004.000,  065484.000, 0006.000, &
		 066130.000, 0008.000,  068550.000, 0040.000, &
		 067388.000, 0005.000,  068107.000, 0020.000, &
		 067237.000, 0004.000,  070290.000, 0060.000, &
		 071888.000, 0014.000,  070746.000, 0005.000, &
		 069842.800, 0002.600,  069717.900, 0002.300, &
		 071110.000, 0050.000,  072985.000, 0011.000, &
		 072951.000, 0006.000,  071167.000, 0003.000, &
		 074060.000, 0021.000,  076642.000, 0023.000, &
		 075222.000, 0011.000,  074129.000, 0005.000, &
		 074506.000, 0006.000,  075978.000, 0009.000, &
		 076030.000, 0005.000,  077290.000, 0050.000, &
		 076814.000, 0005.000,  082857.000, 0018.000, &
		 079296.000, 0007.000,  079007.000, 0002.600, &
		 079339.000, 0005.000,  081338.000, 0012.000, &
		 081994.000, 0006.000,  080900.000, 0004.000, &
		 084711.000, 0024.000,  084083.000, 0011.000, &
		 083788.000, 0005.000,  084835.000, 0007.000, &
		 086848.000, 0013.000,  085482.000, 0007.000, &
		 087550.000, 0050.000,  087793.000, 0017.000, &
		 094234.000, 0029.000,  088585.000, 0007.000, &
		 090220.000, 0030.000,  094018.000, 0011.000, &
		 095840.000, 0050.000,  098280.000, 0040.000, &
		 098130.000, 0070.000,  106580.000, 0040.000, &
		 110090.000, 0060.000 /
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
MODULE nucleus

	USE math
	USE r1r1

	IMPLICIT NONE

	CHARACTER(LEN = 10), DIMENSION(0:1) :: Nucleus_type
	DATA Nucleus_type / "protones", "neutrones" /

	TYPE NucleusType
		TYPE (R1R1Function) func
		INTEGER, DIMENSION(0:2) :: np
		DOUBLE PRECISION, DIMENSION(0:1) :: actual_np, lambda_np, lambda_R2
		DOUBLE PRECISION, DIMENSION(0:2) :: actual_R2
		DOUBLE PRECISION eHFB, R2, norma
		DOUBLE PRECISION pairing ! Apareamiento
		LOGICAL is_blocking(0:1)
		INTEGER num
		INTEGER, DIMENSION(0:1) :: ia, la, ja, mu0
		CHARACTER(LEN = 64) :: filename
	END TYPE

	INTERFACE ASSIGNMENT(=)
		MODULE PROCEDURE Nucleus_copy
	END INTERFACE

	CHARACTER(5), DIMENSION(0:102) :: Nucleus_specie
	DATA Nucleus_specie / "n", &
		"H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", &
		"Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",        &
		"Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",       &
		"Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",       &
		"Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",       &
		"Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", &
		"La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", &
		"Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", &
		"Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", &
		"Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu",             &
		"Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No"/
 CONTAINS

	!---------------------------------------------------------------------------------------!
	! This subroutine initializes an object of type nucleus. It predefines quantities such	!
	! as Z and N, the Fermi levels, the oscillator length, etc. It also defines blocked 	!
	! level as function of the number of particles.						!
	!---------------------------------------------------------------------------------------!

	SUBROUTINE Nucleus_new(nuc, N, Z, b)
		TYPE (NucleusType), INTENT(INOUT) :: nuc
		INTEGER, INTENT(IN) :: N, Z
		DOUBLE PRECISION, INTENT(IN), OPTIONAL:: b

		CHARACTER(LEN = 32) :: specie, N_0_str, A_str
		INTEGER, DIMENSION(1:200) :: ShellModel_L, ShellModel_J

		CALL R1R1Function_new(nuc%func)

		nuc%eHFB    = 0.0D0
		nuc%pairing = 0.0D0
		nuc%R2      = 0.0D0
		nuc%norma   = 1.0D0
		nuc%num     = 1
		nuc%np(0)   = Z
		nuc%np(1)   = N
		nuc%np(2)   = N + Z
		nuc%actual_np = (/ 0.0D0,  0.0D0     /)
		nuc%lambda_np = (/-8.0D0, -4.0D0     /)
		nuc%actual_R2 = (/ 0.0D0,  0.0D0, 0.0D0/)
		nuc%lambda_R2 = (/ 0.0D0,  0.0D0     /)

		IF (PRESENT(b)) THEN
			nuc%func%x = b
		ELSE
			nuc%func%x = Nucleus_get_InitialOscillatorLength(nuc)
		END IF

		! Get the shell model spectrum to define "blocked" levels
		CALL Nucleus_ShellModel(ShellModel_L, ShellModel_J)

		! Defining blocked levels for neutrons or protons (depending)
		CALL Nucleus_set_NeutronBlocking(nuc, ShellModel_L, ShellModel_J)
		CALL Nucleus_set_ProtonBlocking(nuc, ShellModel_L, ShellModel_J)

		! Name of the file that contains the density
		IF (nuc%np(0) .LE. 102) THEN
			specie = Nucleus_specie(nuc%np(0))
		ELSE
			CALL Int2Char(specie, nuc%np(0))
		END IF

		CALL Int2Char(A_str, nuc%np(2))
		CALL Int2Char(N_0_str, N_0)

		nuc%filename = "data/" // TRIM(specie) // TRIM(A_str) // "_" // TRIM(N_0_str)

		RETURN
	END SUBROUTINE Nucleus_new

	!---------------------------------------------------------------------------------------!
	! This subroutine copies the content of one "nucleus" into another			!
	!---------------------------------------------------------------------------------------!

	SUBROUTINE Nucleus_new_Nucleus(n1, n2)
		TYPE (NucleusType), INTENT(INOUT) :: n1
		TYPE (NucleusType), INTENT(IN) :: n2

		CALL R1R1Function_new(n1%func)
		CALL Nucleus_copy(n1, n2)

		RETURN
	END SUBROUTINE Nucleus_new_Nucleus

	!-----------------------------------------------------------------------!
	!  Copy the data of nucleus n2 into nucleus n1. Both n1 and n2 have 	!
	!  the derived type NucleusType						!
	!-----------------------------------------------------------------------!

	SUBROUTINE Nucleus_copy(n1, n2)
		TYPE (NucleusType), INTENT(INOUT) :: n1
		TYPE (NucleusType), INTENT(IN) :: n2

		INTEGER ta

		n1%eHFB = n2%eHFB
		n1%pairing = n2%pairing
		n1%R2 = n2%R2
		n1%num = n2%num

		DO ta = 0, 1

			n1%np(ta) = n2%np(ta)
			n1%actual_np(ta) = n2%actual_np(ta)
			n1%lambda_np(ta) = n2%lambda_np(ta)
			n1%actual_R2(ta) = n2%actual_R2(ta)
			n1%lambda_R2(ta) = n2%lambda_R2(ta)
			n1%is_blocking(ta) = n2%is_blocking(ta)

			IF (n1%is_blocking(ta) .EQV. .TRUE.) THEN
				n1%ia(ta) = n2%ia(ta)
				n1%la(ta) = n2%la(ta)
				n1%ja(ta) = n2%ja(ta)
				n1%mu0(ta) = n2%mu0(ta)
			END IF

		END DO

		n1%np(2) = n1%np(0) + n1%np(1)
		n1%actual_R2(2) = ((n1%np(0) * n1%actual_R2(0)) + (n1%np(1) * n1%actual_R2(1))) / (n1%np(0) + n1%np(1))

		CALL Nucleus_set_b(n1, Nucleus_get_b(n2))

		n1%norma = n2%norma
		n1%filename = n2%filename

		RETURN
	END SUBROUTINE Nucleus_copy

	!-----------------------------------------------------------------------!
	!  In the case of a harmonic oscillator basis, we predefine the oscil-	!
	!  lator length	as: b = 1.01*A**(1/6)					!
	!-----------------------------------------------------------------------!

	FUNCTION Nucleus_get_InitialOscillatorLength(nuc)
		DOUBLE PRECISION Nucleus_get_InitialOscillatorLength
		TYPE (NucleusType), INTENT(IN) :: nuc

		Nucleus_get_InitialOscillatorLength = 1.01D0 * (DBLE(nuc%np(0) + nuc%np(1)) ** (1.0D0 / 6.0D0))

		RETURN
	END FUNCTION Nucleus_get_InitialOscillatorLength

	!-----------------------------------------------------------------------!
	!  We define which level is blocked for the protons and fill in all the	!
	!  relevant quantities of that level					!
	!-----------------------------------------------------------------------!

	SUBROUTINE Nucleus_set_ProtonBlocking(nuc, ShellModel_L, ShellModel_J)
		TYPE (NucleusType), INTENT(INOUT) :: nuc
		INTEGER, DIMENSION(1:200), INTENT(IN) :: ShellModel_L, ShellModel_J

		INTEGER :: Z, N, A
		INTEGER :: lp, jp, ip, mu0p

		IF (nuc%num .EQ. 0) THEN
			nuc%is_blocking(0) = .FALSE.
			RETURN
		END IF

		Z = nuc%np(0)
		N = nuc%np(1)
		A = nuc%np(0) + nuc%np(1)

		nuc%is_blocking(0) = IS_ODD(Z)

		IF (nuc%is_blocking(0)) THEN

			lp = ShellModel_L(Z)
			jp = ShellModel_J(Z)

			ip = lp * 2
			IF (jp .EQ. ((2 * lp) - 1)) THEN
				ip = ip - 1
			END IF
			mu0p = 1

			! Se bloquea el estado mas proximo al nivel de Fermi
			WRITE(*,'("Proton Blocking: l = ",I2," j = ",I2,"/2 - mu = ",i2)') lp,jp,mu0p

			nuc%la(0) = lp
			nuc%ja(0) = jp
			nuc%ia(0) = ip
			nuc%mu0(0) = mu0p

		END IF

		RETURN
	END SUBROUTINE Nucleus_set_ProtonBlocking

	!-----------------------------------------------------------------------!
	!  We define which level is blocked for the neutrons and fill in all the!
	!  relevant quantities of that level					!
	!-----------------------------------------------------------------------!

	SUBROUTINE Nucleus_set_NeutronBlocking(nuc, ShellModel_L, ShellModel_J)
		TYPE (NucleusType), INTENT(INOUT) :: nuc
		INTEGER, DIMENSION(1:200), INTENT(IN) :: ShellModel_L, ShellModel_J

		INTEGER Z, N, A
		INTEGER ln, jn, in, mu0n

		IF (nuc%num .EQ. 0) THEN
			nuc%is_blocking(1) = .FALSE.
			RETURN
		END IF

		Z = nuc%np(0)
		N = nuc%np(1)
		A = nuc%np(0) + nuc%np(1)

		nuc%is_blocking(1) = IS_ODD(N)

		IF (nuc%is_blocking(1)) THEN

			ln = ShellModel_L(N)
			jn = ShellModel_J(N)

			in = ln * 2
			IF (jn .EQ. ((2 * ln) - 1)) THEN
				in = in - 1
			END IF
			mu0n = 1

			! Se bloquea el estado mas proximo al nivel de Fermi
			WRITE(*,'("Neutron Blocking: l = ",I2," j = ",I2,"/2 - mu = ",i2)') ln,jn,mu0n

			nuc%la(1) = ln
			nuc%ja(1) = jn
			nuc%ia(1) = in
			nuc%mu0(1) = mu0n
		END IF
		RETURN
	END SUBROUTINE Nucleus_set_NeutronBlocking

	!-----------------------------------------------------------------------------------------------!
	! This subroutine defines the shell-model orbitals that we will use in the blocking		!
	! For each orbital, we give the orbital angular momentum L and the total angular momentum J	!
	!-----------------------------------------------------------------------------------------------!

	SUBROUTINE Nucleus_ShellModel(ShellModel_L, ShellModel_J)
		INTEGER, DIMENSION(1:200), INTENT(OUT) :: ShellModel_L, ShellModel_J

		INTEGER :: i

		! 1s1/2
		DO i=1,2
			ShellModel_L(i) = 0
			ShellModel_J(i) = 1
		END DO

		! 1p3/2
		DO i=3,6
			ShellModel_L(i) = 1
			ShellModel_J(i) = 3
		END DO

		! 1p1/2
		DO i=7,8
			ShellModel_L(i) = 1
			ShellModel_J(i) = 1
		END DO

		! 1d5/2
		DO i=9,14
			ShellModel_L(i) = 2
			ShellModel_J(i) = 5
		END DO

		! 2s1/2
		DO i=15,16
			ShellModel_L(i) = 0
			ShellModel_J(i) = 1
		END DO

		! 1d3/2
		DO i=17,20
			ShellModel_L(i) = 2
			ShellModel_J(i) = 3
		END DO

		! 1f7/2
		DO i=21,28
			ShellModel_L(i) = 3
			ShellModel_J(i) = 7
		END DO

		! 2p3/2
		DO i=29,32
			ShellModel_L(i) = 1
			ShellModel_J(i) = 3
		END DO

		! 1f5/2
		DO i=33,38
			ShellModel_L(i) = 3
			ShellModel_J(i) = 5
		END DO

		! 2p1/2
		DO i=39,40
			ShellModel_L(i) = 1
			ShellModel_J(i) = 1
		END DO

		! 1g9/2
		DO i=41,50
			ShellModel_L(i) = 4
			ShellModel_J(i) = 9
		END DO

		! 1g7/2
		DO i=51,58
			ShellModel_L(i) = 4
			ShellModel_J(i) = 7
		END DO

		! 2d5/2
		DO i=59,64
			ShellModel_L(i) = 2
			ShellModel_J(i) = 5
		END DO

		! 2d3/2
		DO i=65,68
			ShellModel_L(i) = 2
			ShellModel_J(i) = 3
		END DO

		! 3s1/2
		DO i=69,70
			ShellModel_L(i) = 0
			ShellModel_J(i) = 1
		END DO

		! 1h11/2
		DO i=71,82
			ShellModel_L(i) = 5
			ShellModel_J(i) = 11
		END DO

		! 1h9/2
		DO i=83,92
			ShellModel_L(i) = 5
			ShellModel_J(i) = 9
		END DO

		! 2f7/2
		DO i=93,100
			ShellModel_L(i) = 3
			ShellModel_J(i) = 7
		END DO

		RETURN
	END SUBROUTINE Nucleus_ShellModel

	!-------------------------------------------------------------------------------!
	! Various utilities follow. Taking as input an object of type "Nucleus", it:	!
	!   - returns the mass number (Nucleus_get_A)					!
	!   - returns the oscillator length (Nucleus_get_b)				!
	!   - returns the neutron number (Nucleus_get_N)				!
	!   - returns the proton number (Nucleus_get_Z)					!
	!   - returns the current rms mass radius (Nucleus_get_actual_R2)		!
	!   - returns the square of the current neutron r.m.s. radius (Nucleus_get_R2n)	!
	!   - returns the square of the current proton r.m.s. radius (Nucleus_get_R2p)	!
	!   - sets the neutron number (Nucleus_set_N)					!
	!   - sets the proton number (Nucleus_set_Z)					!
	!   - sets the oscillator length (Nucleus_set_b)				!
	!-------------------------------------------------------------------------------!

	FUNCTION IS_ODD(n)
		LOGICAL IS_ODD
		INTEGER, INTENT(IN) :: n

		IS_ODD = BTEST(n, 0)
		RETURN
	END FUNCTION IS_ODD

	FUNCTION Nucleus_get_A(nuc)
		INTEGER Nucleus_get_A
		TYPE (NucleusType), INTENT(IN) :: nuc

		Nucleus_get_A = nuc%np(0) + nuc%np(1)
		RETURN
	END FUNCTION Nucleus_get_A

	FUNCTION Nucleus_get_b(nuc)
		DOUBLE PRECISION Nucleus_get_b
		TYPE (NucleusType), INTENT(IN) :: nuc

		Nucleus_get_b = nuc%func%x
		RETURN
	END FUNCTION Nucleus_get_b

	! Devuelve el numero de neutrones del nucleo
	FUNCTION Nucleus_get_N(nuc)
		INTEGER Nucleus_get_N
		TYPE (NucleusType), INTENT(IN) :: nuc

		Nucleus_get_N = nuc%np(1)
		RETURN
	END FUNCTION Nucleus_get_N

	! Devuelve el número de protones del núcleo
	FUNCTION Nucleus_get_Z(nuc)
		INTEGER Nucleus_get_Z
		TYPE (NucleusType), INTENT(IN) :: nuc

		Nucleus_get_Z = nuc%np(0)
		RETURN
	END FUNCTION Nucleus_get_Z

	FUNCTION Nucleus_get_actual_R2(nuc)
		DOUBLE PRECISION Nucleus_get_actual_R2
		TYPE (NucleusType), INTENT(IN) :: nuc

		Nucleus_get_actual_R2 = ((nuc%np(0) * nuc%actual_R2(0))  &
		                       + (nuc%np(1) * nuc%actual_R2(1))) &
		                       / DBLE(nuc%np(0) + nuc%np(1))
		RETURN
	END FUNCTION Nucleus_get_actual_R2

	FUNCTION Nucleus_get_R2n(nuc)
		DOUBLE PRECISION Nucleus_get_R2n
		TYPE (NucleusType), INTENT(IN) :: nuc

		Nucleus_get_R2n = nuc%actual_R2(1)
		RETURN
	END FUNCTION Nucleus_get_R2n

	FUNCTION Nucleus_get_R2p(nuc)
		DOUBLE PRECISION Nucleus_get_R2p
		TYPE (NucleusType), INTENT(IN) :: nuc

		Nucleus_get_R2p = nuc%actual_R2(0)
		RETURN
	END FUNCTION Nucleus_get_R2p

	SUBROUTINE Nucleus_set_N(nuc, N)
		TYPE (NucleusType), INTENT(INOUT) :: nuc
		INTEGER, INTENT(IN) :: N

		nuc%np(1) = N
		nuc%np(2) = nuc%np(0) + nuc%np(1)
		RETURN
	END SUBROUTINE Nucleus_set_N

	SUBROUTINE Nucleus_set_Z(nuc, Z)
		TYPE (NucleusType), INTENT(INOUT) :: nuc
		INTEGER, INTENT(IN) :: Z

		nuc%np(0) = Z
		nuc%np(2) = nuc%np(0) + nuc%np(1)
		RETURN
	END SUBROUTINE Nucleus_set_Z

	SUBROUTINE Nucleus_set_b(nuc, b)
		TYPE (NucleusType), INTENT(INOUT) :: nuc
		DOUBLE PRECISION, INTENT(IN) :: b
		CALL R1R1Function_set(nuc%func, b)
		RETURN
	END SUBROUTINE Nucleus_set_b

	!-------------------------------------------------------------------------------!
	!  Subroutine that displays all what is contained in a given "Nucleus" object	!
	!-------------------------------------------------------------------------------!

	SUBROUTINE Nucleus_show_Status(nuc)
		TYPE (NucleusType), INTENT(IN) :: nuc

		INTEGER i

		PRINT *
		PRINT "(A)", "Nucleus:"
		PRINT "(70A1)", ("-", i = 1, 70) ! Muestra una línea de separación

		PRINT "(A30,F15.5)", "Oscillator length (b):", nuc%func%x
		PRINT "(A30,I15,I15)", "Number of particles (np):", nuc%np(0), nuc%np(1)
		PRINT "(A30,F15.5,F15.5)", "(actual_np):", nuc%actual_np(0), nuc%actual_np(1)
		PRINT "(A30,F15.5,F15.5)", "(lambda_np):", nuc%lambda_np(0), nuc%lambda_np(1)
		PRINT "(A30,F15.5,F15.5,F15.5)", "(actual_R2):", nuc%actual_R2(0), nuc%actual_R2(1), nuc%actual_R2(2)
		PRINT "(A30,E15.5,E15.5)", "(lambda_R2):", nuc%lambda_R2(0), nuc%lambda_R2(1)
		PRINT "(A30,F15.5)", "Total Energy (eHFB):", nuc%eHFB
		PRINT "(A30,E15.5)", "Pairing:", nuc%pairing
		PRINT "(A30,E15.5)", "R2:", nuc%R2
		PRINT "(A30,F15.5)", "Norm:", nuc%norma

		PRINT "(70A1)", ("-", i = 1, 70)

		IF (nuc%is_blocking(0)) THEN
			PRINT "(A,I5,A,I5,A,I5)", "Hay bloqueo de neutrones en: la=", nuc%la(0), "ja=", nuc%ja(0), "mu0=", nuc%mu0(0)
		END IF
		IF (nuc%is_blocking(1)) THEN
			PRINT "(A,I5,A,I5,A,I5)", "Hay bloqueo de neutrones en: la=", nuc%la(1), "ja=", nuc%ja(1), "mu0=", nuc%mu0(1)
		END IF
	END SUBROUTINE Nucleus_show_Status

	SUBROUTINE Nucleus_del(nuc)
		TYPE (NucleusType), INTENT(INOUT) :: nuc

		CALL R1R1Function_del(nuc%func)
	END SUBROUTINE Nucleus_del

END MODULE nucleus
 MODULE r1r1

	IMPLICIT NONE

	TYPE R1R1Function
		DOUBLE PRECISION x
		LOGICAL x_set ! Tiene "x" un valor asignado?
		DOUBLE PRECISION x_min, x_max ! Rango de "x"
		LOGICAL infinite_x_min, infinite_x_max ! Es el rango infinito?
		LOGICAL is_even ! is_odd
!		DOUBLE PRECISION func
	END TYPE

 CONTAINS

	SUBROUTINE R1R1Function_new(func)
		TYPE (R1R1Function), INTENT(INOUT) :: func

		func%x_set = .FALSE.
		func%infinite_x_min = .TRUE.
		func%infinite_x_max = .TRUE.
		func%is_even = .FALSE. ! is_odd?
		RETURN
	END SUBROUTINE R1R1Function_new

	SUBROUTINE R1R1Function_new_R1R1Function(func_out, func_in)
		TYPE (R1R1Function), INTENT(INOUT) :: func_out
		TYPE (R1R1Function), INTENT(IN) :: func_in

		func_out%infinite_x_min = func_in%infinite_x_min
		func_out%infinite_x_max = func_in%infinite_x_max
		func_out%x_set = func_in%x_set
		func_out%x = 0.0D0 !TODO func_in%x
		func_out%is_even = func_in%is_even
		RETURN
	END SUBROUTINE R1R1Function_new_R1R1Function

	SUBROUTINE R1R1Function_set(func, X)
		TYPE (R1R1Function), INTENT(INOUT) :: func
		DOUBLE PRECISION, INTENT(IN) :: X

!		The function is assumed constant outside its range
		IF (.NOT. func%infinite_x_min .AND. (X .LE. func%x_min)) THEN
			func%x = func%x_min
		ELSE IF (.NOT. func%infinite_x_max .AND. (X .GE. func%x_max)) THEN
			func%x = func%x_max
		ELSE
			func%x = X
		END IF
		func%x_set = .TRUE.
		RETURN
	END SUBROUTINE R1R1Function_set

	FUNCTION R1R1Function_isValid(func, X)
		LOGICAL R1R1Function_isValid
		TYPE (R1R1Function), INTENT(INOUT) :: func
		DOUBLE PRECISION, INTENT(IN) :: X

		CALL R1R1Function_set(func, X)
		IF (((func%infinite_x_min) .OR. (func%x .GT. func%x_min)) .AND. &
		    ((func%infinite_x_max) .OR. (func%x .LT. func%x_max))) THEN
			R1R1Function_isValid = .TRUE.
		ELSE
			R1R1Function_isValid = .FALSE.
		END IF
		RETURN
	END FUNCTION R1R1Function_isValid

	SUBROUTINE R1R1Function_del(func)
		TYPE (R1R1Function), INTENT(INOUT) :: func

		func%x_set = .FALSE.
		RETURN
	END SUBROUTINE R1R1Function_del

END MODULE r1r1
 !----------------------------------------------------------------------------------------------!
 !    This module defines an object called "selfconsistencymethodproj" and contains a set of	!
 !    functions that calculate various energies (brink-boeker, spin-orbit, coulomb, etc.).	!
 !    It also contains a few utilities to display the status of a fiven hfb iteration on screen.!
 !----------------------------------------------------------------------------------------------!

 MODULE selfc

	USE input
	USE nucleus
	USE symden
	USE symke2b
	USE symvbb
	USE symvc
	USE symvls
	USE symgdd

	IMPLICIT NONE

	! This type is used to contain at once all matrix elements as well as the density

	TYPE SelfConsistencyMethod
		TYPE (SymKineticEnergy2Body) vEkCMph
		TYPE (SymEk2pp) vEkCMpp
		TYPE (SymVBBph) vBBph
		TYPE (SymVBBpp) vBBpp
		TYPE (SymVCph) vCph
		TYPE (SymVCpp) vCpp
		TYPE (SymVLSph) vLSph
		TYPE (SymVLSpp) vLSpp
		TYPE (SymGDDph) gDDph
		TYPE (SymDensity), POINTER :: density
	END TYPE

 CONTAINS

	SUBROUTINE SelfConsistencyMethod_new(consistency, density)
		TYPE (SelfConsistencyMethod), INTENT(INOUT) :: consistency
		TYPE (SymDensity), TARGET, INTENT(IN) :: density

		DOUBLE PRECISION b

		b = Nucleus_get_b(density%nucleus)

		CALL SymKineticEnergy2Body_new(consistency%vEkCMph, consistency%vEkCMpp)
		CALL SymVBBph_new(consistency%vBBph, b)
		CALL SymVBBpp_new(consistency%vBBpp, b)
		CALL SymVCph_new(consistency%vCph)
		CALL SymVCpp_new(consistency%vCpp)
		CALL SymVLSph_new(consistency%vLSph)
		CALL SymVLSpp_new(consistency%vLSpp)
		CALL SymGDDph_new(consistency%gDDph)

		consistency%density => density

		RETURN
	END SUBROUTINE SelfConsistencyMethod_new

	SUBROUTINE SelfConsistencyMethod_store_eHFB(consistency)
		TYPE (SelfConsistencyMethod), INTENT(INOUT) :: consistency

		DOUBLE PRECISION eHF
		DOUBLE PRECISION kinetic_energy, kinetic_CM_energy
		DOUBLE PRECISION local_energy_BB, exchange_energy_BB ! Brink-Booker
		DOUBLE PRECISION local_energy_Coulomb, exchange_energy_Coulomb
		DOUBLE PRECISION dd_energy ! Density Dependent
		DOUBLE PRECISION ls_energy

		DOUBLE PRECISION pairing
		DOUBLE PRECISION pairing_Ek
		DOUBLE PRECISION pairing_BB ! BrinkBooker
		DOUBLE PRECISION pairing_Coulomb
		DOUBLE PRECISION pairing_LS !

		kinetic_energy          = SelfConsistencyMethod_get_Ek(consistency)
		kinetic_CM_energy       = SelfConsistencyMethod_get_EkCM(consistency)
		local_energy_BB         = SelfConsistencyMethod_get_LocalBBEnergy(consistency)
 		local_energy_Coulomb    = SelfConsistencyMethod_get_LocalCoulombEnergy(consistency)
		exchange_energy_BB      = SelfConsistencyMethod_get_ExchangeBBEnergy(consistency)
		exchange_energy_Coulomb = SelfConsistencyMethod_get_ExchangeCoulombEnergy(consistency)
		dd_energy               = SelfConsistencyMethod_get_DDEnergy(consistency)
		ls_energy               = SelfConsistencyMethod_get_LSEnergy(consistency)

		eHF = kinetic_energy + kinetic_CM_energy + local_energy_BB + exchange_energy_BB &
                + local_energy_Coulomb + exchange_energy_Coulomb + dd_energy + ls_energy

		pairing_Ek      = SelfConsistencyMethod_get_EkPairing(consistency)
		pairing_BB      = SelfConsistencyMethod_get_BBPairing(consistency)
		pairing_Coulomb = SelfConsistencyMethod_get_CoulombPairing(consistency)
		pairing_LS      = SelfConsistencyMethod_get_LSPairing(consistency)

		pairing = pairing_BB + pairing_Coulomb + pairing_LS + pairing_Ek

		consistency%density%nucleus%eHFB = eHF + pairing
		consistency%density%nucleus%pairing = pairing

		RETURN
	END SUBROUTINE SelfConsistencyMethod_store_eHFB

	! One-body Kinetic Energy
	FUNCTION SelfConsistencyMethod_get_Ek(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_Ek
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		INTEGER :: A, ta
		DOUBLE PRECISION :: factor, b
		DOUBLE PRECISION, DIMENSION(0:1) :: Ek

		A = Nucleus_get_A(consistency%density%nucleus)
                b = Nucleus_get_b(consistency%density%nucleus)

                factor = (1.0D0 - (1.0D0 / DBLE(A)))

		DO ta = 0, 1
			Ek(ta) = factor * (EkField * consistency%density%field%rho%p(ta))
		END DO

		SelfConsistencyMethod_get_Ek = Ek(0) + Ek(1)

		RETURN
	END FUNCTION SelfConsistencyMethod_get_Ek

	! Two-body Kinetic Energy
	FUNCTION SelfConsistencyMethod_get_EkCM(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_EkCM
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		INTEGER :: A
		TYPE (SymHartreeFockField) :: HF_Gamma

		CALL SymHartreeFockField_new(HF_Gamma)
		CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma, consistency%vEkCMph, consistency%density%field%rho)

		A = Nucleus_get_A(consistency%density%nucleus)

		IF (A .LE. 1) STOP "Abortado"

		IF (switch_CM .GE. 1) THEN
			SelfConsistencyMethod_get_EkCM = (0.5D0 / DBLE(A)) * (consistency%density%field%rho * HF_Gamma)
		ELSE
			SelfConsistencyMethod_get_EkCM = 0.0D0
		END IF

		CALL SymHartreeFockField_del(HF_Gamma)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_EkCM

	! Energia local de Brink-Booker
	FUNCTION SelfConsistencyMethod_get_LocalBBEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_LocalBBEnergy
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		TYPE (SymHartreeFockField) local_gamma

		CALL SymHartreeFockField_new(local_gamma)
		CALL SymVBBph_get_LocalGamma(local_gamma, consistency%vBBph, consistency%density%field%rho)

		SelfConsistencyMethod_get_LocalBBEnergy = 0.5D0 * (consistency%density%field%rho * local_gamma)

		CALL SymHartreeFockField_del(local_gamma)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_LocalBBEnergy

	! Energia de intercambio de Brink-Booker
	FUNCTION SelfConsistencyMethod_get_ExchangeBBEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_ExchangeBBEnergy
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		TYPE (SymHartreeFockField) exchange_gamma

		CALL SymHartreeFockField_new(exchange_gamma)
		CALL SymVBBph_get_ExchangeGamma(exchange_gamma, consistency%vBBph, consistency%density%field%rho)

		SelfConsistencyMethod_get_ExchangeBBEnergy = 0.5D0 * (consistency%density%field%rho * exchange_gamma)

		CALL SymHartreeFockField_del(exchange_gamma)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_ExchangeBBEnergy

	! Coulomb local energy
	FUNCTION SelfConsistencyMethod_get_LocalCoulombEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_LocalCoulombEnergy
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		TYPE (SymHartreeFockField) local_gamma

		CALL SymHartreeFockField_new(local_gamma)
		CALL SymVCph_get_LocalGamma(local_gamma, consistency%vCph, consistency%density%field%rho)

		IF (switch_Coulomb .GE. 1) THEN
			SelfConsistencyMethod_get_LocalCoulombEnergy = 0.5D0 * (consistency%density%field%rho * local_gamma)
		ELSE
			SelfConsistencyMethod_get_LocalCoulombEnergy = 0.0D0
		END IF

		CALL SymHartreeFockField_del(local_gamma)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_LocalCoulombEnergy

	! Exchange Coulomb Energy
	FUNCTION SelfConsistencyMethod_get_ExchangeCoulombEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_ExchangeCoulombEnergy
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		TYPE (SymHartreeFockField) exchange_gamma

		CALL SymHartreeFockField_new(exchange_gamma)
		CALL SymVCph_get_ExchangeGamma(exchange_gamma, consistency%vCph, consistency%density%field%rho)

		IF (switch_Coulomb .GE. 2) THEN
			SelfConsistencyMethod_get_ExchangeCoulombEnergy = 0.5D0 * (consistency%density%field%rho * exchange_gamma)
		ELSE
			SelfConsistencyMethod_get_ExchangeCoulombEnergy = 0.0D0
		END IF

		CALL SymHartreeFockField_del(exchange_gamma)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_ExchangeCoulombEnergy

	! Energia dependiente de la densidad
	FUNCTION SelfConsistencyMethod_get_DDEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_DDEnergy
		TYPE (SelfConsistencyMethod), INTENT(INOUT) :: consistency

		IF (Basis .EQ. 1) THEN
			CALL SymGDDph_make_DD(consistency%gDDph, consistency%density%field%rho)
		ELSE
			CALL Make_DenGenFun(consistency%gDDph, consistency%density%field%rho)
		END IF

		IF (switch_DD .GE. 1) THEN
			SelfConsistencyMethod_get_DDEnergy = SymGDDph_get_edd(consistency%gDDph)
		ELSE
			SelfConsistencyMethod_get_DDEnergy = 0.0D0
		END IF

		RETURN
	END FUNCTION SelfConsistencyMethod_get_DDEnergy

	! Spin-Orbit Energy (direct and exchange)
	FUNCTION SelfConsistencyMethod_get_LSEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_LSEnergy
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		TYPE (SymHartreeFockField) HF_Gamma

		CALL SymHartreeFockField_new(HF_Gamma)
		CALL SymVLSph_get_Gamma(HF_Gamma, consistency%vLSph, consistency%density%field%rho)

		IF (switch_LS .GE. 1) THEN
			SelfConsistencyMethod_get_LSEnergy = 0.5D0* (consistency%density%field%rho * HF_Gamma)
		ELSE
			SelfConsistencyMethod_get_LSEnergy = 0.0D0
		END IF

		CALL SymHartreeFockField_del(HF_Gamma)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_LSEnergy

	! Center of Mass pairing
	FUNCTION SelfConsistencyMethod_get_EkPairing(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_EkPairing
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		TYPE (SymHartreeFockField) :: EkCM_field

		INTEGER :: A
		DOUBLE PRECISION :: factor

		A = Nucleus_get_A(consistency%density%nucleus)

		CALL SymHartreeFockField_new(EkCM_field)
		CALL SymKineticEnergy2Body_get_Delta(EkCM_field, consistency%vEkCMpp, consistency%density%field%kap)

                factor = (1.0D0 / DBLE(A))

		IF (HFOnly .EQ. 0 .AND. switch_CM .EQ. 1) THEN
			SelfConsistencyMethod_get_EkPairing = 0.5D0 * factor * (consistency%density%field%kap * EkCM_field)
		ELSE
			SelfConsistencyMethod_get_EkPairing = 0.0D0
		END IF

		CALL SymHartreeFockField_del(EkCM_field)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_EkPairing

	! Center of Mass pairing (isospin separate)
	FUNCTION SelfConsistencyMethod_get_EkPairingIso(consistency,ta)
		DOUBLE PRECISION SelfConsistencyMethod_get_EkPairingIso
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency
		INTEGER, INTENT(IN) :: ta

		TYPE (SymHartreeFockField) :: EkCM_field

		INTEGER :: A
		DOUBLE PRECISION :: factor

		A = Nucleus_get_A(consistency%density%nucleus)

		CALL SymHartreeFockField_new(EkCM_field)
		CALL SymKineticEnergy2Body_get_DeltaIso(EkCM_field, consistency%vEkCMpp, consistency%density%field%kap, ta)

                factor = (1.0D0 / DBLE(A))

		IF (HFOnly .EQ. 0 .AND. switch_CM .EQ. 1) THEN
			SelfConsistencyMethod_get_EkPairingIso = 0.5D0 * factor * (consistency%density%field%kap * EkCM_field)
		ELSE
			SelfConsistencyMethod_get_EkPairingIso = 0.0D0
		END IF

		CALL SymHartreeFockField_del(EkCM_field)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_EkPairingIso

	! Brink-Boeker pairing
	FUNCTION SelfConsistencyMethod_get_BBPairing(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_BBPairing
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		TYPE (SymHartreeFockField) vbb_field

		CALL SymHartreeFockField_new(vbb_field)
		CALL SymVBBpp_get_Delta(vbb_field, consistency%vBBpp, consistency%density%field%kap)

		IF (HFOnly .EQ. 0) THEN
			SelfConsistencyMethod_get_BBPairing = 0.5D0 * (consistency%density%field%kap * vbb_field)
		ELSE
			SelfConsistencyMethod_get_BBPairing = 0.0D0
		END IF

		CALL SymHartreeFockField_del(vbb_field)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_BBPairing

	! Brink-Boeker pairing (isospin separate)
	FUNCTION SelfConsistencyMethod_get_BBPairingIso(consistency, ta)
		DOUBLE PRECISION SelfConsistencyMethod_get_BBPairingIso
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency
		INTEGER, INTENT(IN) :: ta

		TYPE (SymHartreeFockField) :: vbb_field

		CALL SymHartreeFockField_new(vbb_field)
		CALL SymVBBpp_get_DeltaIso(vbb_field, consistency%vBBpp, consistency%density%field%kap, ta)

		IF (HFOnly .EQ. 0) THEN
			SelfConsistencyMethod_get_BBPairingIso = 0.5D0 * (consistency%density%field%kap * vbb_field)
		ELSE
			SelfConsistencyMethod_get_BBPairingIso = 0.0D0
		END IF

		CALL SymHartreeFockField_del(vbb_field)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_BBPairingIso

	! Coulomb pairing
	FUNCTION SelfConsistencyMethod_get_CoulombPairing(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_CoulombPairing
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		TYPE (SymHartreeFockField) delta

		CALL SymHartreeFockField_new(delta)
		CALL SymVCpp_get_Delta(delta, consistency%vCpp, consistency%density%field%kap)

		IF (HFOnly .EQ. 0 .AND. switch_Coulomb .GE. 1) THEN
			SelfConsistencyMethod_get_CoulombPairing = 0.5D0 * (consistency%density%field%kap * delta)
		ELSE
			SelfConsistencyMethod_get_CoulombPairing = 0.0D0
		END IF

		CALL SymHartreeFockField_del(delta)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_CoulombPairing

	! Spin-orbit pairing
	FUNCTION SelfConsistencyMethod_get_LSPairing(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_LSPairing
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		TYPE (SymHartreeFockField) delta

		CALL SymHartreeFockField_new(delta)
		CALL SymVLSpp_get_Delta(delta, consistency%vLSpp, consistency%density%field%kap)

		IF (HFOnly .EQ. 0 .AND. switch_LS .GE. 1) THEN
			SelfConsistencyMethod_get_LSPairing = 0.5D0 * (consistency%density%field%kap * delta)
		ELSE
			SelfConsistencyMethod_get_LSPairing = 0.0D0
		END IF

		CALL SymHartreeFockField_del(delta)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_LSPairing

	! Spin-orbit pairing (isospin separate)
	FUNCTION SelfConsistencyMethod_get_LSPairingIso(consistency, ta)
		DOUBLE PRECISION SelfConsistencyMethod_get_LSPairingIso
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency
		INTEGER, INTENT(IN) :: ta

		TYPE (SymHartreeFockField) :: delta

		CALL SymHartreeFockField_new(delta)
		CALL SymVLSpp_get_DeltaIso(delta, consistency%vLSpp, consistency%density%field%kap, ta)

		IF (HFOnly .EQ. 0 .AND. switch_LS .GE. 1) THEN
			SelfConsistencyMethod_get_LSPairingIso = 0.5D0 * (consistency%density%field%kap * delta)
		ELSE
			SelfConsistencyMethod_get_LSPairingIso = 0.0D0
		END IF

		CALL SymHartreeFockField_del(delta)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_LSPairingIso

	FUNCTION SelfConsistencyMethod_accuracy(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_accuracy
		TYPE (SelfConsistencyMethod), INTENT(INOUT) :: consistency

		DOUBLE PRECISION old_eHFB

		old_eHFB = consistency%density%nucleus%eHFB

		CALL SelfConsistencyMethod_store_eHFB(consistency)

		SelfConsistencyMethod_accuracy = ABS(old_eHFB - consistency%density%nucleus%eHFB)
		RETURN
	END FUNCTION SelfConsistencyMethod_accuracy

	SUBROUTINE SelfConsistencyMethod_show_Status(consistency)
		TYPE (SelfConsistencyMethod), INTENT(INOUT) :: consistency

		INTEGER :: i
		DOUBLE PRECISION :: b

		DOUBLE PRECISION :: eHF
		DOUBLE PRECISION :: kinetic_energy, kinetic_CM_energy
		DOUBLE PRECISION :: local_energy_BB, exchange_energy_BB ! Brink-Booker
		DOUBLE PRECISION :: local_energy_Coulomb, exchange_energy_Coulomb
		DOUBLE PRECISION :: dd_energy ! Density Dependent
		DOUBLE PRECISION :: ls_energy

		DOUBLE PRECISION :: pairing, pairing_n, pairing_p
		DOUBLE PRECISION :: pairing_Ekn, pairing_BBn, pairing_LSn, pairing_Ekp, pairing_BBp, pairing_LSp
		DOUBLE PRECISION :: pairing_Ek
		DOUBLE PRECISION :: pairing_BB ! BrinkBooker
		DOUBLE PRECISION :: pairing_Coulomb
		DOUBLE PRECISION :: pairing_LS !

		DOUBLE PRECISION :: N, Z

		b = Nucleus_get_b(consistency%density%nucleus)

		WRITE(*,'(/,5X,"SUMMARY OF THE RUN - ENERGIES")')
		WRITE(*,'(5X,"=============================")')

		PRINT *
		PRINT "(A50,F20.10)", "Longitud del oscilador:", b
		PRINT "(70A1)", ("-", i = 1, 70)

		kinetic_energy = SelfConsistencyMethod_get_Ek(consistency)
		kinetic_CM_energy = SelfConsistencyMethod_get_EkCM(consistency)
		local_energy_BB = SelfConsistencyMethod_get_LocalBBEnergy(consistency)
 		local_energy_Coulomb = SelfConsistencyMethod_get_LocalCoulombEnergy(consistency)
		exchange_energy_BB = SelfConsistencyMethod_get_ExchangeBBEnergy(consistency)
		exchange_energy_Coulomb = SelfConsistencyMethod_get_ExchangeCoulombEnergy(consistency)
		dd_energy = SelfConsistencyMethod_get_DDEnergy(consistency)
		ls_energy = SelfConsistencyMethod_get_LSEnergy(consistency)

		PRINT "(A50,F20.10)", "Energia cinetica:", kinetic_energy
		PRINT "(A50,F20.10)", "Energia cinetica del centro de masas:", kinetic_CM_energy
		PRINT "(A50,F20.10)", "Energia local de Brink-Booker:", local_energy_BB
		PRINT "(A50,F20.10)", "Energia local de Coulomb:", local_energy_Coulomb
		PRINT "(A50,F20.10)", "Energia de intercambio de Brink-Booker:", exchange_energy_BB
		PRINT "(A50,F20.10)", "Energia de intercambio de Coulomb:", exchange_energy_Coulomb
		PRINT "(A50,F20.10)", "Energia dependiente de la densidad:", dd_energy
		PRINT "(A50,F20.10)", "Energia de spin-orbita:", ls_energy
		PRINT "(70A1)", ("-", i = 1, 70)

		eHF = kinetic_energy + kinetic_CM_energy &
			+ local_energy_BB + exchange_energy_BB &
			+ local_energy_Coulomb + exchange_energy_Coulomb &
			+ dd_energy + ls_energy

		PRINT "(A50,F20.10)", "Energia total de Brink-Booker:", local_energy_BB + exchange_energy_BB
		PRINT "(A50,F20.10)", "Energia total de Coulomb:", local_energy_Coulomb + exchange_energy_Coulomb
		PRINT "(A50,F20.10)", "Energia total (Hartree-Fock):", eHF
		PRINT "(70A1)", ("-", i = 1, 70)

		pairing_Ek = SelfConsistencyMethod_get_EkPairing(consistency)
		pairing_BB = SelfConsistencyMethod_get_BBPairing(consistency)
		pairing_Coulomb = SelfConsistencyMethod_get_CoulombPairing(consistency)
		pairing_LS = SelfConsistencyMethod_get_LSPairing(consistency)

		pairing_Ekn = SelfConsistencyMethod_get_EkPairingIso(consistency, 1)
		pairing_Ekp = SelfConsistencyMethod_get_EkPairingIso(consistency, 0)
		pairing_BBn = SelfConsistencyMethod_get_BBPairingIso(consistency, 1)
		pairing_BBp = SelfConsistencyMethod_get_BBPairingIso(consistency, 0)
		pairing_LSn = SelfConsistencyMethod_get_LSPairingIso(consistency, 1)
		pairing_LSp = SelfConsistencyMethod_get_LSPairingIso(consistency, 0)

                pairing_n = pairing_Ekn + pairing_BBn + pairing_LSn
                pairing_p = pairing_Ekp + pairing_BBp + pairing_LSp + pairing_Coulomb

		PRINT "(A50,E20.10)", "Apareamiento del centro de masas:", pairing_Ek
		PRINT "(A50,E20.10)", "Apareamiento de Brink-Booker:", pairing_BB
		PRINT "(A50,E20.10)", "Apareamiento de Coulomb:", pairing_Coulomb
		PRINT "(A50,E20.10)", "Apareamiento de spin-orbita:", pairing_LS
		PRINT "(70A1)", ("-", i = 1, 70)

		pairing = pairing_BB + pairing_Coulomb + pairing_LS + pairing_Ek

		PRINT "(A50,E20.10)", "Apareamiento total:", pairing
        	PRINT "(A50,4E20.10)", "Apareamiento neutrons:",  pairing_n, pairing_Ekn, pairing_BBn, pairing_LSn
        	PRINT "(A50,4E20.10)", "Apareamiento protons:",  pairing_p, pairing_Ekp, pairing_BBp, pairing_LSp
		PRINT "(70A1)", ("-", i = 1, 70)

		PRINT "(A50,F20.10)", "TOTAL:", eHF + pairing

		PRINT *
		PRINT "(A40,A15,A15)", "Neutrones", "Protones", "Total"

		N = consistency%density%nucleus%actual_np(1)
		Z = consistency%density%nucleus%actual_np(0)

		PRINT "(A25,F15.5,F15.5,F15.5)", "Particulas:", N, Z, N + Z
		PRINT "(A25,F15.5,F15.5)", "Potenciales quimicos:", &
			consistency%density%nucleus%lambda_np(1), &
			consistency%density%nucleus%lambda_np(0)
		PRINT "(A25,F15.5,F15.5,F15.5)", "Radio:", &
			SQRT(consistency%density%nucleus%actual_R2(1)), &
			SQRT(consistency%density%nucleus%actual_R2(0)), &
			SQRT(Nucleus_get_actual_R2(consistency%density%nucleus))
		PRINT "(A25,F15.5,F15.5,F15.5)", "Multiplicadores y error:", &
			consistency%density%nucleus%lambda_R2(1), &
			consistency%density%nucleus%lambda_R2(0), &
			Nucleus_get_actual_R2(consistency%density%nucleus) - consistency%density%nucleus%R2
! D.showSpatialDistribution(NEU)
!TODO		CALL Nucleus_show_ExperimentalData(consistency%density%nucleus)
		RETURN
	END SUBROUTINE SelfConsistencyMethod_show_Status

	SUBROUTINE SelfConsistencyMethod_del(consistency)
		TYPE (SelfConsistencyMethod), INTENT(INOUT) :: consistency

		CALL SymKineticEnergy2Body_del(consistency%vEkCMph, consistency%vEkCMpp)

		CALL SymVBBph_del(consistency%vBBph)
		CALL SymVBBpp_del(consistency%vBBpp)
		CALL SymVCph_del(consistency%vCph)
		CALL SymVCpp_del(consistency%vCpp)
		CALL SymVLSph_del(consistency%vLSph)
		CALL SymVLSpp_del(consistency%vLSpp)
		CALL SymGDDph_del(consistency%gDDph)

		NULLIFY(consistency%density)

		RETURN
	END SUBROUTINE SelfConsistencyMethod_del

END MODULE selfc
!------------------------------------------------------------------------------!
!                                                                              !
!                       DEFINITION OF SPECIFIC TYPES                           !
!                                                                              !
!                                                                              !
!  In spherical symmetry, the calculation of the matrix elements often makes   !
!  use of tensors of the form:                                                 !
!                           T_{la, na, nb}                                     !
!                                                                              !
!  which correspond to the reduced form of a more general object               !
!                                                                              !
!         T_{na, nb, la} = delta_{ac} * T_{a,c}, with a = (na, la ja)          !
!                                                                              !
!  See infamous PhD, page 45, top of the page. THIS IS SPECIFIC TO SPHERICAL   !
!  SYMMETRY.                                                                   !
!                                                                              !
!  The definitions of the types below and of their related operations provide  !
!  an appropriate storage facility for all these objects. Each of these "3D    !
!  tensors" has a total size:                                                  !
!                                                                              !
!        (N+1) * n_max(l) * n_max(l),     with: n_max(l) = (N - l)/ 2 + 1      !
!                                                                              !
!  SymD2Tensor: n_max(l) * n_max(l) array T_{la, na, nb} for a given la        !
!  SymD3Tensor: Full (N+1) * n_max(l) * n_max(l) array T_{la, na, nb}          !
!  SymD2Tensor_SymD3Tensor:                                                    !
!  SymD3Tensor_SymD3Tensor: "Matrix" of SymD3Tensor's                          !
!                                                                              !
!                                                                              !
!  NOTA BENE: FOR ADAPTATION TO THE CASE OF A RANDOM SPHERICAL BASIS (NOT THE  !
!             HARMONIC OSCILLATOR) A METHOD MUST BE FOUND TO DEFINE N0	       !
!                                                                              !
!------------------------------------------------------------------------------!

 MODULE symd3t

	USE input

	IMPLICIT NONE

        ! SymD2Tensor points to a 2D array or real numbers
	TYPE SymD2Tensor
		DOUBLE PRECISION, DIMENSION(:, :), POINTER :: d2
	END TYPE

        ! SymD3Tensor points to a vector or SymD2Tensor, so it points to a vector of
	! matrices.
	! It is the sophisticated version of a 3D array
	TYPE SymD3Tensor
		TYPE (SymD2Tensor), DIMENSION(:), POINTER :: d3tensor
	END TYPE

        ! SymD2Tensor_ SymD3Tensor points to an array or SymD3Tensor.
	! This is the sophisticated version of a 5D array
	TYPE SymD2Tensor_SymD3Tensor
		TYPE (SymD3Tensor), DIMENSION(:, :), POINTER :: d2
	END TYPE

        ! SymD3Tensor_ SymD3Tensor points to a vector or SymD2Tensor_SymD3Tensor.
	! This is the sophisticated version of a 6D array
	TYPE SymD3Tensor_SymD3Tensor
		TYPE (SymD2Tensor_SymD3Tensor), DIMENSION(:), POINTER :: d3tensor
	END TYPE

	INTERFACE ASSIGNMENT(=)
		MODULE PROCEDURE &
			SymD3Tensor_assign1, &
			SymD3Tensor_assign2, &
			SymD3Tensor_SymD3Tensor_assign1, &
			SymD3Tensor_SymD3Tensor_assign2
	END INTERFACE

!	INTERFACE OPERATOR(+)
!		MODULE PROCEDURE &
!			SymD3Tensor_add, &
!			SymD3Tensor_SymD3Tensor_add
!	END INTERFACE

	INTERFACE OPERATOR(*)
		MODULE PROCEDURE SymD3Tensor_product2
	END INTERFACE

 CONTAINS

	! Creating an object SymD3Tensor by allocating memory for it
	SUBROUTINE SymD3Tensor_new(t)
		TYPE (SymD3Tensor), INTENT(INOUT) :: t

		INTEGER :: u1, d

                ! Allocating memory for the vector of SymD2Tensor tensors
		ALLOCATE(t%d3tensor(0:Lmax))

                ! For each SymD2Tensor of the aforementioned vector, allocating memory for the
		! corresponding array
		DO u1 = 0, Lmax
								d = Min(Nmax, NmaxOfL(u1))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

			ALLOCATE(t%d3tensor(u1)%d2(d, d))

			t%d3tensor(u1)%d2(:,:) = 0.0D0

		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_new

	FUNCTION SymD3Tensor_read(t, file_desc, file_error)
		LOGICAL SymD3Tensor_read
		TYPE (SymD3Tensor), INTENT(INOUT) :: t
		INTEGER, INTENT(IN) :: file_desc
		INTEGER, INTENT(INOUT) :: file_error

		DOUBLE PRECISION :: dummy
		INTEGER :: u1, u2, u3, d
		INTEGER :: v1, v2, v3

		DO u1 = 0, Lmax
								d = Min(Nmax, NmaxOfL(u1))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

			DO u2 = 1, d
				DO u3 = 1, u2
					READ (file_desc, FMT=*, IOSTAT=file_error) &
						v1, v2, v3, dummy
					t%d3tensor(u1)%d2(u2, u3) = dummy
					IF ((file_error .NE. 0) .OR. (u1 .NE. v1) .OR. (u2 .NE. v2) .OR. (u3 .NE. v3)) THEN
						SymD3Tensor_read = .FALSE.
						RETURN
					END IF
				END DO
			END DO
		END DO
		SymD3Tensor_read = .TRUE.
		RETURN
	END FUNCTION SymD3Tensor_read

	SUBROUTINE SymD3Tensor_write(t, file_desc, file_error)
		TYPE (SymD3Tensor), INTENT(IN) :: t
		INTEGER, INTENT(IN) :: file_desc
		INTEGER, INTENT(INOUT) :: file_error

		INTEGER :: u1, u2, u3, d

		DO u1 = 0, Lmax
								d = Min(Nmax, NmaxOfL(u1))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

			DO u2 = 1, d
				DO u3 = 1, u2
					WRITE (file_desc, FMT="(3I3,E24.16)", IOSTAT=file_error) &
						u1, u2, u3, t%d3tensor(u1)%d2(u2, u3)
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_write

	! Calculating the trace of a 3D tensor. We consider only hermitian matrices in the following,
	! whose trace is a real number.
	FUNCTION SymD3Tensor_trace(t_in)
		DOUBLE PRECISION SymD3Tensor_trace
		TYPE (SymD3Tensor), INTENT(IN) :: t_in

		INTEGER :: la, na, d

		SymD3Tensor_trace = 0.0D0

		DO la = 0, Lmax
								d = Min(Nmax, NmaxOfL(la))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - la) / 2) + 1

			DO na = 1, d
				SymD3Tensor_trace = SymD3Tensor_trace + t_in%d3tensor(la)%d2(na, na)
			END DO
		END DO
		RETURN
	END FUNCTION SymD3Tensor_trace

	! Functions which give as output a square matrix containing the matrix of a given tensor
	! operator in the harmonic oscillator basis

	! Case of a random basis on a mesh (Basis=0)
	FUNCTION SymD3Tensor_matrix(t_in, la)
		DOUBLE PRECISION, DIMENSION(Nmax, Nmax) :: SymD3Tensor_matrix

		TYPE (SymD3Tensor), INTENT(INOUT) :: t_in
		INTEGER, INTENT(IN) :: la

		INTEGER :: u1, u2, d

		! Reconstruction of the matrix of the tensor from the lower triangular part

							d = Min(Nmax, NmaxOfL(la))
		IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - la) / 2) + 1

		DO u1 = 1, d
			DO u2 = 1, u1
				SymD3Tensor_matrix(u1, u2) = t_in%d3tensor(la)%d2(u1, u2)
				IF (u1 .NE. u2) THEN
					SymD3Tensor_matrix(u2, u1) = t_in%d3tensor(la)%d2(u1, u2)
				END IF
			END DO
		END DO

		! In case of the harmonic oscillator basis (analytical), the matrix is slightly oversized.
		! We make sure that all "extra" elements are null

		IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) THEN

			DO u1 = d+1, Nmax
				DO u2 = 1, u1
					SymD3Tensor_matrix(u1, u2) = 0.0D0
				END DO
			END DO

		END IF

		RETURN
	END FUNCTION SymD3Tensor_matrix

	! Case of the HO basis (Basis=1)
	FUNCTION SymD3Tensor_matrix_1(t_in, la)
		DOUBLE PRECISION, DIMENSION(((N_0 - la) / 2) + 1, ((N_0 - la) / 2) + 1) :: SymD3Tensor_matrix_1

		TYPE (SymD3Tensor), INTENT(INOUT) :: t_in
		INTEGER, INTENT(IN) :: la

		INTEGER :: u1, u2, d

		! Reconstruction of the matrix of the tensor from the lower triangular part

		d = ((N_0 - la) / 2) + 1

		DO u1 = 1, d
			DO u2 = 1, u1
				SymD3Tensor_matrix_1(u1, u2) = t_in%d3tensor(la)%d2(u1, u2)
				IF (u1 .NE. u2) THEN
					SymD3Tensor_matrix_1(u2, u1) = t_in%d3tensor(la)%d2(u1, u2)
				END IF
			END DO
		END DO

		RETURN
	END FUNCTION SymD3Tensor_matrix_1

	! Copy the content of t_in into t_out
	SUBROUTINE SymD3Tensor_assign1(t_out, t_in)
		TYPE (SymD3Tensor), INTENT(INOUT) :: t_out
		TYPE (SymD3Tensor), INTENT(IN) :: t_in

		INTEGER :: u1, u2, u3, d

		DO u1 = 0, Lmax
								d = Min(Nmax, NmaxOfL(u1))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

			DO u2 = 1, d
				DO u3 = 1, u2
					t_out%d3tensor(u1)%d2(u2, u3) = t_in%d3tensor(u1)%d2(u2, u3)
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_assign1

	! Fill in all elements of t_out with a single number R1
	SUBROUTINE SymD3Tensor_assign2(t_out, R1)
		TYPE (SymD3Tensor), INTENT(INOUT) :: t_out
		DOUBLE PRECISION, INTENT(IN) :: R1

		INTEGER :: u1, u2, u3, d

		DO u1 = 0, Lmax
								d = Min(Nmax, NmaxOfL(u1))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

			DO u2 = 1, d
				DO u3 = 1, u2
					t_out%d3tensor(u1)%d2(u2, u3) = R1
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_assign2

	! Fill in all elements of t_out with same l-value with a given square matrix of number
	SUBROUTINE SymD3Tensor_assign_Matrix(t_out, u1, M)
		TYPE (SymD3Tensor), INTENT(INOUT) :: t_out
		INTEGER, INTENT(IN) :: u1
		DOUBLE PRECISION, DIMENSION(:, :), INTENT(IN) :: M

		INTEGER :: u2, u3, d

							d = Min(Nmax, NmaxOfL(u1))
		IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

		DO u2 = 1, d
			DO u3 = 1, u2
				t_out%d3tensor(u1)%d2(u2, u3) = M(u2, u3)
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_assign_Matrix

	SUBROUTINE SymD3Tensor_add(t_out, t1_in, t2_in)
		TYPE (SymD3Tensor), INTENT(INOUT) :: t_out
		TYPE (SymD3Tensor), INTENT(IN) :: t1_in, t2_in

		INTEGER :: u1, u2, u3, d

		DO u1 = 0, Lmax
								d = Min(Nmax, NmaxOfL(u1))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

			DO u2 = 1, d
				DO u3 = 1, u2
					t_out%d3tensor(u1)%d2(u2, u3) = &
						t1_in%d3tensor(u1)%d2(u2, u3) + &
						t2_in%d3tensor(u1)%d2(u2, u3)
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_add

	SUBROUTINE SymD3Tensor_product(t_out, R1, t_in)
		TYPE (SymD3Tensor), INTENT(INOUT) :: t_out
		DOUBLE PRECISION, INTENT(IN) :: R1
		TYPE (SymD3Tensor), INTENT(IN) :: t_in

		INTEGER :: u1, u2, u3, d

		DO u1 = 0, Lmax
								d = Min(Nmax, NmaxOfL(u1))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

			DO u2 = 1, d
				DO u3 = 1, u2
					t_out%d3tensor(u1)%d2(u2, u3) = R1 * t_in%d3tensor(u1)%d2(u2, u3)
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_product

	FUNCTION SymD3Tensor_product2(t1_in, t2_in)
		DOUBLE PRECISION :: SymD3Tensor_product2
		TYPE (SymD3Tensor), INTENT(IN) :: t1_in, t2_in

		DOUBLE PRECISION :: sum_la, sum1
		INTEGER :: u1, u2, u3, d

		SymD3Tensor_product2 = 0.0D0

		DO u1 = 0, Lmax
								d = Min(Nmax, NmaxOfL(u1))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

			sum_la = 0.0D0

			DO u2 = 1, d
				DO u3 = 1, u2
					sum1 = t1_in%d3tensor(u1)%d2(u2, u3) * t2_in%d3tensor(u1)%d2(u2, u3)
					IF (u2 .EQ. u3) THEN
						sum_la = sum_la + sum1
					ELSE
						sum_la = sum_la + (2.0D0 * sum1)
					END IF
				END DO
			END DO

			SymD3Tensor_product2 = SymD3Tensor_product2 + sum_la
		END DO
		RETURN
	END FUNCTION SymD3Tensor_product2

	FUNCTION SymD3Tensor_distance(t1_in, t2_in)
		DOUBLE PRECISION SymD3Tensor_distance
		TYPE (SymD3Tensor), INTENT(IN) :: t1_in, t2_in

		INTEGER :: u1, u2, u3, d
		DOUBLE PRECISION :: distance

		SymD3Tensor_distance = 0.0D0

		DO u1 = 0, Lmax
								d = Min(Nmax, NmaxOfL(u1))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

			DO u2 = 1, d
				DO u3 = 1, u2
					distance = ABS(t1_in%d3tensor(u1)%d2(u2, u3) - t2_in%d3tensor(u1)%d2(u2, u3))
					IF (distance .GT. SymD3Tensor_distance) THEN
						SymD3Tensor_distance = distance
					END IF
				END DO
			END DO
		END DO
		RETURN
	END FUNCTION SymD3Tensor_distance

	SUBROUTINE SymD3Tensor_print(t)
		TYPE (SymD3Tensor), INTENT(IN) :: t

		INTEGER :: u1, u2, u3, d

		DO u1 = 0, Lmax
			PRINT "(I2,A)", u1, ":"

								d = Min(Nmax, NmaxOfL(u1))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

			DO u2 = 1, d
				PRINT "(13F10.5)", (t%d3tensor(u1)%d2(u2, u3), u3=1, u2)
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_print

	! Empty the memory reserved for a 3D tensor

	SUBROUTINE SymD3Tensor_del(t)
		TYPE (SymD3Tensor), INTENT(INOUT) :: t

		INTEGER :: u1

		DO u1 = 0, Lmax
			DEALLOCATE(t%d3tensor(u1)%d2)
		END DO

		DEALLOCATE(t%d3tensor)

		RETURN
	END SUBROUTINE SymD3Tensor_del

	! Reserves (allocates) the memory for a "matrix" of 3D tensors

	SUBROUTINE SymD3Tensor_SymD3Tensor_new(t)
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(INOUT) :: t

		INTEGER :: u1, u2, u3, d

		ALLOCATE(t%d3tensor(0:Lmax))

		DO u1 = 0, Lmax
								d = Min(Nmax, NmaxOfL(u1))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

			ALLOCATE(t%d3tensor(u1)%d2(d, d))

			DO u2 = 1, d
				DO u3 = 1, u2
					CALL SymD3Tensor_new(t%d3tensor(u1)%d2(u2, u3))
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_SymD3Tensor_new

	SUBROUTINE SymD3Tensor_SymD3Tensor_assign(t_out, la, na, nc, lb, nb, nd, R1)
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(INOUT) :: t_out
		INTEGER, INTENT(IN) :: la, na, nc, lb, nb, nd
		DOUBLE PRECISION, INTENT(IN) :: R1

		t_out%d3tensor(la)%d2(na, nc)%d3tensor(lb)%d2(nb, nd) = R1

		RETURN
	END SUBROUTINE SymD3Tensor_SymD3Tensor_assign

	! NEW - NEW
	SUBROUTINE SymD3Tensor_SymD3Tensor_get(R1, t_in, la, na, nc, lb, nb, nd)
		DOUBLE PRECISION, INTENT(OUT) :: R1
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(IN) :: t_in
		INTEGER, INTENT(IN) :: la, na, nc, lb, nb, nd

		R1 = t_in%d3tensor(la)%d2(na, nc)%d3tensor(lb)%d2(nb, nd)
		RETURN
	END SUBROUTINE SymD3Tensor_SymD3Tensor_get

	SUBROUTINE SymD3Tensor_SymD3Tensor_assign1(t_out, t_in)
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(INOUT) :: t_out
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(IN) :: t_in

		INTEGER :: u1, u2, u3, d

		DO u1 = 0, Lmax
								d = Min(Nmax, NmaxOfL(u1))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

			DO u2 = 1, d
				DO u3 = 1, u2
					t_out%d3tensor(u1)%d2(u2, u3) = t_in%d3tensor(u1)%d2(u2, u3)
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_SymD3Tensor_assign1

	SUBROUTINE SymD3Tensor_SymD3Tensor_assign2(t_out, R1)
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(INOUT) :: t_out
		DOUBLE PRECISION, INTENT(IN) :: R1

		INTEGER :: u1, u2, u3, d

		DO u1 = 0, Lmax
								d = Min(Nmax, NmaxOfL(u1))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

			DO u2 = 1, d
				DO u3 = 1, u2
					t_out%d3tensor(u1)%d2(u2, u3) = R1
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_SymD3Tensor_assign2

	! Defines the addition of 2 "matrices" of 3D tensors

	SUBROUTINE SymD3Tensor_SymD3Tensor_add(t_out, t1_in, t2_in)
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(INOUT) :: t_out
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(IN) :: t1_in, t2_in

		INTEGER :: u1, u2, u3, d

		DO u1 = 0, Lmax
								d = Min(Nmax, NmaxOfL(u1))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

			DO u2 = 1, d
				DO u3 = 1, u2
					CALL SymD3Tensor_add(t_out%d3tensor(u1)%d2(u2, u3), &
						t1_in%d3tensor(u1)%d2(u2, u3), &
						t2_in%d3tensor(u1)%d2(u2, u3))
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_SymD3Tensor_add

	! Defines the product of 2 "matrices" of 3D tensors

	SUBROUTINE SymD3Tensor_SymD3Tensor_product(t_out, t1_in, t2_in)
		TYPE (SymD3Tensor), INTENT(INOUT) :: t_out
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(IN) :: t1_in
		TYPE (SymD3Tensor), INTENT(IN) :: t2_in

		INTEGER :: u1, u2, u3, d

		DO u1 = 0, Lmax
								d = Min(Nmax, NmaxOfL(u1))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

			DO u2 = 1, d
				DO u3 = 1, u2
					t_out%d3tensor(u1)%d2(u2, u3) = t1_in%d3tensor(u1)%d2(u2, u3) * t2_in
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_SymD3Tensor_product

	! Empty the memory reserved for a "matrix" of 3D tensors

	SUBROUTINE SymD3Tensor_SymD3Tensor_del(t)
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(INOUT) :: t

		INTEGER :: u1, u2, u3, d

		DO u1 = 0, Lmax
								d = Min(Nmax, NmaxOfL(u1))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

			DO u2 = 1, d
				DO u3 = 1, u2
					CALL SymD3Tensor_del(t%d3tensor(u1)%d2(u2, u3))
				END DO
			END DO

			DEALLOCATE(t%d3tensor(u1)%d2)

		END DO

		DEALLOCATE(t%d3tensor)

		RETURN
	END SUBROUTINE SymD3Tensor_SymD3Tensor_del

END MODULE symd3t
 MODULE symden

	USE input
	USE math
	USE global
	USE symd3t
	USE nucleus
	USE symfield

	IMPLICIT NONE

	!---------------------------------------------------------------!
	!    SymDensity is a type for objects associated to:		!
	!	- a nucleus 						!
	!       - a couple of mean-field and pairing field Gamma and 	!
	!         Delta							!
	!---------------------------------------------------------------!

	TYPE SymDensity
		TYPE (NucleusType) nucleus
		TYPE (SymHartreeFockBogolField) field
	END TYPE

	INTERFACE ASSIGNMENT(=)
		MODULE PROCEDURE SymDensity_assign
	END INTERFACE

 CONTAINS

	SUBROUTINE SymDensity_new(density, N, Z)
		! Recibe como parametros de entrada el registro de densidad
		! y el numero de protones y electrones
		TYPE (SymDensity), INTENT(INOUT) :: density
		INTEGER, INTENT(IN) :: N, Z
		DOUBLE PRECISION :: b

		CALL SymHartreeFockBogolField_new(density%field)
		CALL Nucleus_new(density%nucleus, N, Z, b)

		CALL SymDensity_read(density)

		RETURN
	END SUBROUTINE SymDensity_new

	SUBROUTINE SymDensity_new_SymDensity(density_out, density_in)
		TYPE (SymDensity), INTENT(INOUT) :: density_out
		TYPE (SymDensity), INTENT(IN) :: density_in

		CALL SymHartreeFockBogolField_new(density_out%field)
		CALL Nucleus_new_Nucleus(density_out%nucleus, density_in%nucleus)

		density_out%field%rho = density_in%field%rho
		density_out%field%kap = density_in%field%kap

		RETURN
	END SUBROUTINE SymDensity_new_SymDensity

	SUBROUTINE SymDensity_initialize(density)
		TYPE (SymDensity), INTENT(INOUT) :: density

		INTEGER :: ta, la, na, i
		INTEGER :: nosc, fin_apa, num, den, lup
		DOUBLE PRECISION :: vv, vu

		! Calculo de la densidad inicial
		DO ta = 0, 1

			! Initialization of density matrices
			density%field%rho%p(ta) = 0.0D0
			density%field%rho%a(ta) = 0.0D0
			density%field%kap%p(ta) = 0.0D0
			density%field%kap%a(ta) = 0.0D0

			! Actual number of particles for isospin ta
			num = density%nucleus%np(ta)
			i = 1
			IF (N_0 .LE. 1) THEN
				IF (num .GE. MagicNumber(N_0)) STOP "Abortado: SymDensity_initialize (num >= MagicNumber(N_0))"
			END IF

			DO WHILE (MagicNumber(i) .LE. num)
				i = i + 1
			END DO
			fin_apa = i

			IF (ta .EQ. 1) THEN
				PRINT "(/A,I2,A,I3)", "Hay apareamiento de los ", density%nucleus%np(ta), " neutrones hasta el estado mas bajo de la capa ", i
			ELSE
				PRINT "(/A,I2,A,I3)", "Hay apareamiento de los ", density%nucleus%np(ta), " protones hasta el estado mas bajo de la capa ", i
			END IF

			IF (HFonly.EQ.0) THEN
 				den = MagicNumber(i)
				vv = DBLE(num) / DBLE(den)
   				vu = SQRT(vv) * SQRT(1.0D0 - vv)
				DO nosc = i - 1, 0, -1
					! Convenio u = (-) ^ l * |u|
					vu = vu * PAR(nosc)
                        	        lup = MIN(nosc, Lmax)
					DO la = lup, 0, -2
						na = ((nosc - la) / 2) + 1
						density%field%rho%p(ta)%d3tensor(la)%d2(na, na) = 2.0D0 * DBLE(2*la + 1) * vv
						density%field%kap%p(ta)%d3tensor(la)%d2(na, na) = 2.0D0 * DBLE(2*la + 1) * vu
					END DO
				END DO

				la = MIN(fin_apa, Lmax)
				na = 1
				density%field%rho%p(ta)%d3tensor(la)%d2(na, na) = 2.0D0 * DBLE(la + 1) * vv
				density%field%rho%a(ta)%d3tensor(la)%d2(na, na) = vv
				density%field%kap%p(ta)%d3tensor(la)%d2(na, na) = 2.0D0 * DBLE(la + 1) * vu
				density%field%kap%a(ta)%d3tensor(la)%d2(na, na) = vu
			ELSE
 				den = MagicNumber(fin_apa)
				vv = DBLE(num) / DBLE(den)
   				vu = SQRT(vv) * SQRT(1.0D0 - vv)
				DO nosc = i - 1, 0, -1
					! Convenio u = (-) ^ l * |u|
					vu = vu * PAR(nosc)
                        	        lup = MIN(nosc, Lmax)
					DO la = lup, 0, -2
						na = ((nosc - la) / 2) + 1
						density%field%rho%p(ta)%d3tensor(la)%d2(na, na) = 2.0D0 * DBLE(2*la + 1) * vv
						density%field%kap%p(ta)%d3tensor(la)%d2(na, na) = 2.0D0 * DBLE(2*la + 1) * vu
					END DO
				END DO

				la = MIN(fin_apa, Lmax)
				na = 1
				density%field%rho%p(ta)%d3tensor(la)%d2(na, na) = 2.0D0 * DBLE(la + 1) * vv
				density%field%rho%a(ta)%d3tensor(la)%d2(na, na) = vv
				density%field%kap%p(ta)%d3tensor(la)%d2(na, na) = 2.0D0 * DBLE(la + 1) * vu
				density%field%kap%a(ta)%d3tensor(la)%d2(na, na) = vu
			END IF

			density%nucleus%actual_np(ta) = SymD3Tensor_trace(density%field%rho%p(ta))
		END DO

		RETURN
	END SUBROUTINE SymDensity_initialize

	SUBROUTINE SymDensity_read(density)
		TYPE (SymDensity), INTENT(INOUT) :: density

		IF(SymHartreeFockBogolField_read(density%field, density%nucleus%filename)) THEN
			density%nucleus%actual_np(0) = SymD3Tensor_trace(density%field%rho%p(0))
			density%nucleus%actual_np(1) = SymD3Tensor_trace(density%field%rho%p(1))
		ELSE
			CALL SymDensity_initialize(density)
		END IF
		RETURN
	END SUBROUTINE SymDensity_read

	SUBROUTINE SymDensity_save(density)
		TYPE (SymDensity), INTENT(INOUT) :: density

		CALL SymHartreeFockBogolField_write(density%field, density%nucleus%filename)

		RETURN
	END SUBROUTINE SymDensity_save

	SUBROUTINE SymDensity_assign(density_out, density_in)
		TYPE (SymDensity), INTENT(INOUT) :: density_out
		TYPE (SymDensity), INTENT(IN) :: density_in

		density_out%field = density_in%field

		RETURN
	END SUBROUTINE SymDensity_assign

	SUBROUTINE SymDensity_store_actual_R2(density)
		TYPE (SymDensity), INTENT(INOUT) :: density

		DOUBLE PRECISION :: b2
		INTEGER :: ta

		b2 = Nucleus_get_b(density%nucleus) ** 2

		DO ta=0, 1
			density%nucleus%actual_R2(ta) = R2Field * density%field%rho%p(ta) / density%nucleus%np(ta)
		END DO

		density%nucleus%actual_R2(2) = &
			 ((density%nucleus%np(0) * density%nucleus%actual_R2(0))  &
			+ (density%nucleus%np(1) * density%nucleus%actual_R2(1))) &
			/ (density%nucleus%np(0) + density%nucleus%np(1))

		RETURN
	END SUBROUTINE SymDensity_store_actual_R2

	SUBROUTINE SymDensity_shuffle(density)
		TYPE (SymDensity), INTENT(INOUT) :: density

		INTEGER :: ta, la, nosc
		INTEGER :: i, fin_apa, num, den, na
		DOUBLE PRECISION :: vv, vu

		DO ta = 0, 1

			density%field%kap%p(ta) = 0.0D0
			density%field%kap%a(ta) = 0.0D0
			num = density%nucleus%np(ta)

			IF (N_0 .LE. 8) THEN
				IF (num .GE. MagicNumber(N_0)) STOP "Numero de particulas no valido"
			END IF
			i = 1
			DO WHILE (MagicNumber(i) .LE. num)
				i = i + 1
			END DO
			fin_apa = i

			den = MagicNumber(i)
			vv = DBLE(num) / den
			vu = SQRT(vv) * SQRT(1.0D0 - vv)

			DO nosc = i - 1, 0, -1
				vu = vu * PAR(nosc) ! convenio u = (-)^l * |u|
				DO la = nosc, 0, -2
					na = ((nosc - la) / 2) + 1
					density%field%kap%p(ta)%d3tensor(la)%d2(na, na) = 2.0D0 * DBLE(2*la + 1) * vu
				END DO
			END DO
			la = fin_apa
			na = 1

			density%field%kap%p(ta)%d3tensor(la)%d2(na, na) = 2.0D0 * DBLE(la + 1) * vu
			density%field%kap%a(ta)%d3tensor(la)%d2(na, na) = vu
		END DO
		RETURN
	END SUBROUTINE SymDensity_shuffle

	!SUBROUTINE SymDensity_show_SpatialDistribution(density, ta)
	!	TYPE (SymDensity), INTENT(INOUT) :: density
	!	INTEGER, INTENT(IN) :: ta
        !
	!	RETURN
	!END SUBROUTINE SymDensity_show_SpatialDistribution

	!SUBROUTINE SymDensity_show_ParticleDensity(density)
	!	TYPE (SymDensity), INTENT(INOUT) :: density
        !
	!	RETURN
	!END SUBROUTINE SymDensity_show_ParticleDensity

	SUBROUTINE SymDensity_del(density)
		TYPE (SymDensity), INTENT(INOUT) :: density

		CALL Nucleus_del(density%nucleus)
		CALL SymHartreeFockBogolField_del(density%field)

		RETURN
	END SUBROUTINE SymDensity_del

END MODULE symden
 MODULE symfield

	USE input
	USE global
	USE symd3t
	USE symtalm

	IMPLICIT NONE

	!-------------------------------------------------------------------------------!
	!  SymHartreeFockField : Made of 4 tensors (2 for protons, 2 for neutrons)	!
	!    - p refers to the local part of the HF field: Gamma_local			!
	!    - a refers to exchange part of the HF field: Gamma_exchange		!
	!-------------------------------------------------------------------------------!

	TYPE SymHartreeFockField
		TYPE (SymD3Tensor), DIMENSION(0:1) :: p
		TYPE (SymD3Tensor), DIMENSION(0:1) :: a
	END TYPE

	!-------------------------------------------------------------------------------!
	!  SymHartreeFockBogoField: Made of 2 SymHartreeFockField tensors 		!
	!     - rho refers to the p-h mean-field density: Gamma (local and exchange)	!
	!     - kap refers to p-p pairing field: Delta	 (local and exchange)		!
	!-------------------------------------------------------------------------------!

	TYPE SymHartreeFockBogolField
		TYPE (SymHartreeFockField) rho
		TYPE (SymHartreeFockField) kap
	END TYPE

	INTERFACE ASSIGNMENT(=)
		MODULE PROCEDURE &
			SymHartreeFockField_assign, &
			SymHartreeFockBogolField_assign
	END INTERFACE

	INTERFACE OPERATOR(*)
		MODULE PROCEDURE SymHartreeFockField_product2
	END INTERFACE

 CONTAINS

	SUBROUTINE SymHartreeFockBogolField_new(HFB)
		TYPE (SymHartreeFockBogolField), INTENT(INOUT) :: HFB

		CALL SymHartreeFockField_new(HFB%rho)
		CALL SymHartreeFockField_new(HFB%kap)

		RETURN
	END SUBROUTINE SymHartreeFockBogolField_new

	FUNCTION SymHartreeFockBogolField_read(HFB, filename)
		LOGICAL SymHartreeFockBogolField_read
		TYPE (SymHartreeFockBogolField), INTENT(INOUT) :: HFB
		CHARACTER(*), INTENT(IN) :: filename

		INTEGER, PARAMETER :: file_desc = 16 ! Descriptor de fichero
		INTEGER :: file_error, new_N_0, new_Lmax, new_Nmax, new_Basis, new_CompHO

		OPEN (file_desc, FILE=TRIM(filename), STATUS="OLD", ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			SymHartreeFockBogolField_read = .FALSE.
			RETURN
		END IF

		READ (file_desc, *) new_N_0, new_Lmax, new_Nmax, new_Basis, new_CompHO

		IF (new_Basis .NE. Basis) THEN
			CLOSE (file_desc)
			SymHartreeFockBogolField_read = .FALSE.
			RETURN
		ELSE
			IF (Basis .EQ. 1) THEN ! Case HO analytical
				IF (new_N_0 .NE. N_0) THEN
					CLOSE (file_desc)
					SymHartreeFockBogolField_read = .FALSE.
					RETURN
				END IF
			ELSE ! Case WS or HO numerical
				IF (new_compHO .NE. CompHO) THEN
					CLOSE (file_desc)
					SymHartreeFockBogolField_read = .FALSE.
					RETURN
				ELSE
					IF (compHO .EQ. 1 .AND. new_N_0 .NE. N_0) THEN
						CLOSE (file_desc)
						SymHartreeFockBogolField_read = .FALSE.
						RETURN
					END IF
					IF (compHO .EQ. 0 .AND. (new_Lmax .NE. Lmax .OR. new_Nmax .NE. Nmax)) THEN
						CLOSE (file_desc)
						SymHartreeFockBogolField_read = .FALSE.
						RETURN
					END IF
				END IF
			END IF
		END IF

		IF ((SymHartreeFockField_read(HFB%rho, file_desc, file_error)) .AND. &
		    (SymHartreeFockField_read(HFB%kap, file_desc, file_error))) THEN
			SymHartreeFockBogolField_read = .TRUE.
		ELSE
			SymHartreeFockBogolField_read = .FALSE.
		END IF
		CLOSE (file_desc)

		RETURN
	END FUNCTION SymHartreeFockBogolField_read

	! Almacena los datos de los campos rho y kap en un fichero
	SUBROUTINE SymHartreeFockBogolField_write(HBF, filename)
		TYPE (SymHartreeFockBogolField), INTENT(IN) :: HBF
		CHARACTER(*), INTENT(IN) :: filename

		INTEGER, PARAMETER :: file_desc = 16 ! Descriptor de fichero
		INTEGER file_error

		OPEN (file_desc, FILE=TRIM(filename), ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) RETURN

		! El primer valor almacenado en el fichero es el numero de capas: N_0
		WRITE (file_desc, FMT="(5I3)", IOSTAT=file_error) N_0, Lmax, Nmax, Basis, CompHO

		! A continuacion se almacenan rho y kap
		CALL SymHartreeFockField_write(HBF%rho, file_desc, file_error)
		CALL SymHartreeFockField_write(HBF%kap, file_desc, file_error)

		CLOSE (file_desc)
		RETURN
	END SUBROUTINE SymHartreeFockBogolField_write

	SUBROUTINE SymHartreeFockBogolField_assign(HBF1, HBF2)
		TYPE (SymHartreeFockBogolField), INTENT(INOUT) :: HBF1
		TYPE (SymHartreeFockBogolField), INTENT(IN) :: HBF2

		HBF1%rho = HBF2%rho
		HBF1%kap = HBF2%kap

		RETURN
	END SUBROUTINE SymHartreeFockBogolField_assign

	SUBROUTINE SymHartreeFockBogolField_del(HFB)
		TYPE (SymHartreeFockBogolField), INTENT(INOUT) :: HFB

		CALL SymHartreeFockField_del(HFB%rho)
		CALL SymHartreeFockField_del(HFB%kap)

		RETURN
	END SUBROUTINE SymHartreeFockBogolField_del

	SUBROUTINE SymHartreeFockField_new(HF)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF

		INTEGER ta

		DO ta = 0, 1
			CALL SymD3Tensor_new(HF%p(ta))
			CALL SymD3Tensor_new(HF%a(ta))
		END DO

		RETURN
	END SUBROUTINE SymHartreeFockField_new

	FUNCTION SymHartreeFockField_read(HF, file_desc, file_error)
		LOGICAL SymHartreeFockField_read
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF
		INTEGER, INTENT(IN) :: file_desc
		INTEGER, INTENT(INOUT) :: file_error

		INTEGER ta

		DO ta = 0, 1
			IF (.NOT. SymD3Tensor_read(HF%p(ta), file_desc, file_error)) THEN
				SymHartreeFockField_read = .FALSE.
				RETURN
			END IF
			IF (.NOT. SymD3Tensor_read(HF%a(ta), file_desc, file_error)) THEN
				SymHartreeFockField_read = .FALSE.
				RETURN
			END IF
		END DO
		SymHartreeFockField_read = .TRUE.
		RETURN
	END FUNCTION SymHartreeFockField_read

	SUBROUTINE SymHartreeFockField_write(HF, file_desc, file_error)
		TYPE (SymHartreeFockField), INTENT(IN) :: HF
		INTEGER, INTENT(IN) :: file_desc
		INTEGER, INTENT(INOUT) :: file_error

		INTEGER ta

		DO ta = 0, 1
			CALL SymD3Tensor_write(HF%p(ta), file_desc, file_error)
			CALL SymD3Tensor_write(HF%a(ta), file_desc, file_error)
		END DO

		RETURN
	END SUBROUTINE SymHartreeFockField_write

	SUBROUTINE SymHartreeFockField_assign(HF_out, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		INTEGER ta

		DO ta = 0, 1

			HF_out%p(ta) = HF_in%p(ta)
			HF_out%a(ta) = HF_in%a(ta)

		END DO

		RETURN
	END SUBROUTINE SymHartreeFockField_assign

	SUBROUTINE SymHartreeFockField_add(HF_out, HF1_in, HF2_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymHartreeFockField), INTENT(IN) :: HF1_in, HF2_in

		CALL SymD3Tensor_add(HF_out%p(0), HF1_in%p(0), HF2_in%p(0))
		CALL SymD3Tensor_add(HF_out%p(1), HF1_in%p(1), HF2_in%p(1))

		CALL SymD3Tensor_add(HF_out%a(0), HF1_in%a(0), HF2_in%a(0))
		CALL SymD3Tensor_add(HF_out%a(1), HF1_in%a(1), HF2_in%a(1))

		RETURN
	END SUBROUTINE SymHartreeFockField_add

	SUBROUTINE SymHartreeFockField_add_SymD3Tensor(HF_out, t_in, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymD3Tensor), INTENT(IN) :: t_in
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		CALL SymD3Tensor_add(HF_out%p(0), t_in, HF_in%p(0))
		CALL SymD3Tensor_add(HF_out%p(1), t_in, HF_in%p(1))

		HF_out%a(0) = HF_in%a(0)
		HF_out%a(1) = HF_in%a(1)

		RETURN
	END SUBROUTINE SymHartreeFockField_add_SymD3Tensor

	SUBROUTINE SymHartreeFockField_product(HT_out, R1, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HT_out
		DOUBLE PRECISION, INTENT(IN) :: R1
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		CALL SymD3Tensor_product(HT_out%p(0), R1, HF_in%p(0))
		CALL SymD3Tensor_product(HT_out%p(1), R1, HF_in%p(1))

		CALL SymD3Tensor_product(HT_out%a(0), R1, HF_in%a(0))
		CALL SymD3Tensor_product(HT_out%a(1), R1, HF_in%a(1))

		RETURN
	END SUBROUTINE SymHartreeFockField_product

	FUNCTION SymHartreeFockField_product2(HF1_in, HF2_in)
		DOUBLE PRECISION :: SymHartreeFockField_product2
		TYPE (SymHartreeFockField), INTENT(IN) :: HF1_in, HF2_in

		DOUBLE PRECISION, DIMENSION(0:1) :: sum2
		INTEGER :: ta

		DO ta = 0, 1
			sum2(ta) = (HF1_in%p(ta) * HF2_in%p(ta)) &
			         + (HF1_in%a(ta) * HF2_in%a(ta))
		END DO

		SymHartreeFockField_product2 = sum2(0) + sum2(1)

		RETURN
	END FUNCTION SymHartreeFockField_product2

	! Used to make the contraction of a field of given isospin with
	! a density of the same isospin. The result is a real number
	! (which depends on the isospin).

	FUNCTION SymHartreeFockField_product_iso(HF1_in, HF2_in, ta)
		DOUBLE PRECISION :: SymHartreeFockField_product_iso
		TYPE (SymHartreeFockField), INTENT(IN) :: HF1_in, HF2_in

		DOUBLE PRECISION :: sum2
		INTEGER :: ta

		sum2 = (HF1_in%p(ta) * HF2_in%p(ta)) &
		     + (HF1_in%a(ta) * HF2_in%a(ta))

		SymHartreeFockField_product_iso = sum2

		RETURN
	END FUNCTION SymHartreeFockField_product_iso

	FUNCTION SymHartreeFockField_distance(HF1_in, HF2_in)
		DOUBLE PRECISION SymHartreeFockField_distance
		TYPE (SymHartreeFockField), INTENT(IN) :: HF1_in, HF2_in

		SymHartreeFockField_distance = MAX( &
			SymD3Tensor_distance(HF1_in%p(0), HF2_in%p(0)), &
			SymD3Tensor_distance(HF1_in%p(1), HF2_in%p(1)))

		RETURN
	END FUNCTION SymHartreeFockField_distance

	FUNCTION SymHartreeFockField_ChargeDensity(HF, r)
		DOUBLE PRECISION SymHartreeFockField_ChargeDensity
		TYPE (SymHartreeFockField), INTENT(IN) :: HF
		DOUBLE PRECISION, INTENT(IN) :: r

		DOUBLE PRECISION x, d1, xk, sum1, sumlb
		INTEGER k, lb, nb, nd, nbmax, p2max, p2

		x = r * r
		xk = 1.0D0
		sum1 = 0.0D0

		DO k = 0, Lmax

			sumlb = 0.0D0

			DO lb = 0, k
									nbmax = MIN(Nmax, NmaxOfL(lb))
				IF (basis .EQ. 1 .OR. CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

				DO nb = 1, nbmax
					DO nd = 1, nbmax
						p2max = nb + nd - 2
						p2 = k - lb
						IF (p2 .GT. p2max) CYCLE
						d1 = SymCoefficientB_get(nb - 1, lb, nd - 1, lb, p2)
						sumlb = sumlb + (d1 * HF%p(PROTON)%d3tensor(lb)%d2(nd, nb))
					END DO
				END DO
			END DO
			sum1 = sum1 + (xk * sumlb)
			xk = xk * x
		END DO
		SymHartreeFockField_ChargeDensity = EXP(-x) * sum1
		RETURN
	END FUNCTION SymHartreeFockField_ChargeDensity

	SUBROUTINE SymHartreeFockField_del(HF)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF

		INTEGER ta

		DO ta = 0, 1
			CALL SymD3Tensor_del(HF%p(ta))
			CALL SymD3Tensor_del(HF%a(ta))
		END DO

		RETURN
	END SUBROUTINE SymHartreeFockField_del

END MODULE symfield
 MODULE symgdd

	USE input
	USE global
	USE math
	USE symtalm
	USE symfield

	IMPLICIT NONE

	TYPE SymGDDph
		DOUBLE PRECISION, DIMENSION(:, :), POINTER :: dLag
	END TYPE

 CONTAINS

	SUBROUTINE SymGDDph_new(gDDph)
		TYPE (SymGDDph), INTENT(INOUT) :: gDDph

		SELECT CASE (Basis)

		CASE(1)
			ALLOCATE(gDDph%dLag(0:2, NLag))
			IF (.NOT. ASSOCIATED(gDDph%dLag)) STOP "Unable to allocate memory"
		CASE(2)
			ALLOCATE(gDDph%dLag(0:2, Npoint))
			IF (.NOT. ASSOCIATED(gDDph%dLag)) STOP "Unable to allocate memory"
		END SELECT

		RETURN
	END SUBROUTINE SymGDDph_new

	!-----------------------------------------------------------------------!
	!  Subroutine that calculates the
	!-----------------------------------------------------------------------!

	SUBROUTINE SymGDDph_update(HF_out, gDDph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymGDDph), INTENT(INOUT) :: gDDph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		INTEGER :: ta, la, na, namax, nc

		! Filling in the new density-dependent fields Gamma_{na, nc, la, ta} from the new density

		SELECT CASE (Basis)

		CASE(1)

			! Calculate the new density for both isospin channels from the current wave-functions and fields (HF_in)

                        CALL SymGDDph_make_DD(gDDph, HF_in)

			DO ta = 0, 1

				DO la = 0, Lmax
					namax = ((N_0 - la) / 2) + 1
					DO na = 1, namax
						DO nc = 1, na
							HF_out%p(ta)%d3tensor(la)%d2(na, nc) = SymGDDph_G1dd(gDDph, na - 1, nc - 1, la, ta)
							HF_out%a(ta)%d3tensor(la)%d2(na, nc) = DBLE(0.0)
						END DO
					END DO
				END DO

			END DO

		CASE(2)

			! Calculate the new density for both isospin channels from the current wave-functions and fields (HF_in)

                        CALL Make_DenGenFun(gDDph, HF_in)

			DO ta = 0, 1

				DO la = 0, Lmax
								namax = MIN(Nmax,NmaxOfL(la))
					IF (CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

					DO na = 1, namax
						DO nc = 1, na
							HF_out%p(ta)%d3tensor(la)%d2(na, nc) = SymGDDph_Gendd(gDDph, na, nc, la, ta)
							HF_out%a(ta)%d3tensor(la)%d2(na, nc) = DBLE(0.0)
						END DO
					END DO
				END DO

			END DO

		END SELECT
		RETURN
	END SUBROUTINE SymGDDph_update

	!-----------------------------------------------------------------------!
	!									!
	!  Subroutine that calculates the density-dependent HF field 		!
	!  Gamma_{na, nc, la, ta} of page page 134 (but without the last term 	!
	!  ontaining the pairing tensor (for each isospin). It proceeds to the	!
	!  integration of the auxiliary function Gamma^{DD}(r) of page 133.	!
	!									!
	! 	      CASE OF A SPHERICAL HARMONIC OSCILLATOR BASIS		!
	!									!
	!-----------------------------------------------------------------------!

	FUNCTION SymGDDph_G1dd(gDDph, na, nc, la, ta)
		DOUBLE PRECISION :: SymGDDph_G1dd
		TYPE (SymGDDph), INTENT(IN) :: gDDph
		INTEGER, INTENT(IN) :: na, nc, la, ta

		DOUBLE PRECISION :: d1, d2, d3, d4, d5
		INTEGER :: p1, p1max, i
#if(USE_QUADRUPLE==1)
		REAL(KIND = 16) :: sump1, sumi
#else
		DOUBLE PRECISION :: sump1, sumi
#endif

		p1max = na + nc
		sump1 = 0.0D0

		DO p1 = 0, p1max
			sumi = 0.0D0

			! Below, there is a hidden integration (due to the Laguerre points)

			DO i = 1, NLag

				! d1 = density for isospin  ta  at point ri (ri is a node for the Laguerre integration)
				! d2 = density for isospin 1-ta at point ri (ri is a node for the Laguerre integration)
				d1 = gDDph%dLag(    ta, i)
				d2 = gDDph%dLag(1 - ta, i)

				! d3 = total density at point ri (ri is a node for the Laguerre integration)
				d3 = d1 + d2

				IF (d3 .LT. 0.0D0) STOP "ERROR! Fuera de rango"

				! d4 IS UNUSED
				d4 = d3 ** (ALPHA - 1.0D0)

				d5 = (1.0D0 + 0.5D0 * x0) * d3 * d3 - (x0 + 0.5D0) * d1 * (d3 - ALPHA * d2)
! aLag = x/(ALPHA+2.)
!TODO			sumi =  sumi + EXP(LOG(GaussLQ%gauss%w(i)) + LOG(aLag(i)) * (p1 + la) + LOG(d3) * (ALPHA - 1) + LOG(d5))

				!

				sumi =  sumi + EXP(LOG(GaussLQ%gauss%w(i)) &
					+ LOG(GaussLQ%gauss%x(i) / (ALPHA + 2.0D0)) &
					* DBLE(p1 + la) + LOG(d3) * (ALPHA - 1.0D0) + LOG(d5))

			END DO
			sump1 = sump1 + SymCoefficientB_get(na, la, nc, la, p1) * sumi
		END DO

		SymGDDph_G1dd = Gogny_t0(Gogny) * 0.5D0 * I_SALPHA3 * sump1

		RETURN
	END FUNCTION SymGDDph_G1dd

	!-----------------------------------------------------------------------!
	!									!
	!   Subroutine that calculates the density for each isospin channel	!
	!									!
	! 	      CASE OF A SPHERICAL HARMONIC OSCILLATOR BASIS		!
	!									!
	!-----------------------------------------------------------------------!

	SUBROUTINE SymGDDph_make_DD(gDDph, HF)
		TYPE (SymGDDph), INTENT(INOUT) :: gDDph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF

		INTEGER :: s, i, lb, nb, nbmax, nd, p2, p2max
		DOUBLE PRECISION :: d1
		DOUBLE PRECISION, ALLOCATABLE :: pows(:)
		DOUBLE PRECISION, ALLOCATABLE :: wksp(:, :)

		ALLOCATE(wksp(0:1, 0:N_0))
		ALLOCATE(pows(NLag))

		DO s = 0, N_0

			wksp(1, s) = 0.0D0
			wksp(0, s) = 0.0D0

			! Calculation of rho^tau(r) for each isospin (stored in wksp)

			DO lb = 0, s
				nbmax = ((N_0 - lb) / 2) + 1
				p2 = s - lb
				DO nb = 1, nbmax
					DO nd = 1, nb

						p2max = nb + nd - 2

						IF (p2 .GT. p2max) CYCLE

						IF (nb .EQ. nd) THEN
							d1 =             I_4PI * SymCoefficientB_get(nb - 1, lb, nd - 1, lb, p2)
						ELSE
							d1 = DBLE(2.0) * I_4PI * SymCoefficientB_get(nb - 1, lb, nd - 1, lb, p2)
						END IF

						wksp(1, s) = wksp(1, s) + (d1 * HF%p(1)%d3tensor(lb)%d2(nb, nd))
						wksp(0, s) = wksp(0, s) + (d1 * HF%p(0)%d3tensor(lb)%d2(nb, nd))

					END DO
				END DO
			END DO
		END DO

		! Inicializamos a 1 la tabla "pows" y a 0 las tablas "dLag"
		DO i = 1, NLag
			pows(i) = 1.0D0
			gDDph%dLag(1, i) = 0.0D0
			gDDph%dLag(0, i) = 0.0D0
		END DO

		! Schunck: I added the 1/b_0**3 to put all quantities dependent of the oscillator length apart
		!          (that is to say in the calculation of the radial density). This implies that the
		!          final result obtained from SymGDDph_get_edd should be multipied by some funny factor
		!          to recover the good result. See my notes for further explanations. The point is:
		!	   the calculation of the DD FIELD is now independent of the oscillator length, AS IT
		!	   SHOULD BE.

		DO s = 0, N_0
			DO i = 1, NLag
				gDDph%dLag(1, i) = gDDph%dLag(1, i) + (pows(i) * wksp(1, s))/(b_0**3)
				gDDph%dLag(0, i) = gDDph%dLag(0, i) + (pows(i) * wksp(0, s))/(b_0**3)
! aLag = x/(ALPHA+2.)
!TODO			pows(i) = pows(i) * aLag(i)

				pows(i) = pows(i) * (GaussLQ%gauss%x(i) / DBLE(ALPHA + 2.0))
			END DO
		END DO

		DEALLOCATE(wksp)
		DEALLOCATE(pows)

		RETURN
	END SUBROUTINE SymGDDph_make_DD

	!-----------------------------------------------------------------------!
	!									!
	!   Subroutine that calculates the density rho(r) for each isospin 	!
	!   channel taking as radial wave-functions arbitrary functions		!
	!									!
	!             CASE OF AN ARBITRARY SPHERICAL BASIS			!
	!									!
	!-----------------------------------------------------------------------!

	SUBROUTINE Make_DenGenFun(gDDph, HF)
		TYPE (SymGDDph), INTENT(INOUT) :: gDDph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF

		INTEGER :: IndexBra, IndexKet, nb, nd, lb, nbmax, i

		DOUBLE PRECISION :: d1, Radius
		DOUBLE PRECISION, ALLOCATABLE :: wksp(:, :)

		ALLOCATE(wksp(0:1, Npoint))

		DO i = 1, Npoint

			Radius = RadMesh(i)

			wksp(1, i) = 0.0D0
			wksp(0, i) = 0.0D0

			DO lb = 0, Lmax
							nbmax = MIN(Nmax,NmaxOfL(lb))
				IF (CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

				DO nb = 1, nbmax
					DO nd = 1, nb

						IndexBra = IndexVecNL(nb,lb)
						IndexKet = IndexVecNL(nd,lb)

						IF (IndexBra .EQ. 0 .OR. IndexKet .EQ. 0) CYCLE

						IF (nb .EQ. nd) THEN
							d1 = WaveFun(i,IndexBra)*WaveFun(i,IndexKet) / Radius**2
						ELSE
							d1 = 2.0D0*WaveFun(i,IndexBra)*WaveFun(i,IndexKet) / Radius**2
						END IF

						wksp(1, i) = wksp(1, i) + (d1 * HF%p(1)%d3tensor(lb)%d2(nb, nd))
						wksp(0, i) = wksp(0, i) + (d1 * HF%p(0)%d3tensor(lb)%d2(nb, nd))

					END DO
				END DO
			END DO

		END DO ! end of Ipoint loop

		! Filling in the pointers containing the radial densities for each isospin

		DO i = 1, Npoint
			Radius = RadMesh(i)
			gDDph%dLag(1, i) = I_4PI * wksp(1, i)
			gDDph%dLag(0, i) = I_4PI * wksp(0, i)
		END DO

		DEALLOCATE(wksp)

		RETURN
	END SUBROUTINE Make_DenGenFun

	!-----------------------------------------------------------------------!
	!									!
	!  Subroutine that calculates the density-dependent HF field 		!
	!  Gamma_{na, nc, la, ta} of page page 134 (but without the last term 	!
	!  containing the pairing tensor (for each isospin). It proceeds to the	!
	!  integration of the auxiliary function Gamma^{DD}(r) of page 133.	!
	!									!
	!             CASE OF AN ARBITRARY SPHERICAL BASIS			!
	!									!
	!-----------------------------------------------------------------------!

	FUNCTION SymGDDph_Gendd(gDDph, na, nc, la, ta)
		DOUBLE PRECISION SymGDDph_Gendd
		TYPE (SymGDDph), INTENT(IN) :: gDDph

		INTEGER, INTENT(IN) :: na, nc, la, ta

		INTEGER :: IndexBra, IndexKet, i
		DOUBLE PRECISION :: d1, d2, d3, d5, res
		DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)

		IndexBra = IndexVecNL(na,la)
		IndexKet = IndexVecNL(nc,la)

                IF (IndexBra .EQ. 0 .OR. IndexKet .EQ. 0) THEN
			SymGDDph_Gendd = 0.0D0
			RETURN
                END IF

		! Calculating the Gamma^{DD}(r) field at each point i

	        ALLOCATE(Integrand(1:Npoint))

		DO i = 1, Npoint

			! d1 = density for isospin  ta  at point ri
			! d2 = density for isospin 1-ta at point ri
			d1 = gDDph%dLag(    ta, i)
			d2 = gDDph%dLag(1 - ta, i)

			! d3 = total density at point ri (ri is a node for the Laguerre integration)
			d3 = d1 + d2

			IF (d3 .LT. 0.0D0) THEN
				WRITE(*,'("i = ",I4," d1 = ",F20.15," d2 = ",F20.15," d3 = ",F20.15)') i,d1,d2,d3
				STOP "ERROR - Negative density in SymGDDph_Gendd"
			END IF

			d5 = (1.0D0 + 0.5D0 * x0) * d3 * d3 - (x0 + 0.5D0) * d1 * (d3 - ALPHA * d2)

			! The function to integrate is R_a(r) x Gamma^{DD}(r) x R_c(r)
			! Gamma^{DD}(r) is calculated above with the proton and neutron densities
			! contained in gDDph%dLag
			IF (ABS(d3) .LT. 1.D-14) THEN
				Integrand(i) = 0.0D0
			ELSE
				Integrand(i) = EXP(LOG(d3) * (ALPHA - 1.0D0) + LOG(d5)) &
				     * WaveFun(i,IndexBra)*WaveFun(i,IndexKet)
			END IF

		END DO

		! Calculating the Gamma_{na, nc, la, ta}^{DD} field by integration of Gamma^{DD}(r) over r

		CALL simps(Integrand,Npoint,MeshStep,res)

		DEALLOCATE(Integrand)

		SymGDDph_Gendd = Gogny_t0(Gogny) * res

		RETURN
	END FUNCTION SymGDDph_Gendd

	!-----------------------------------------------------------------------!
	!   This subroutine gives the density-dependent part of the mean-	!
	!   field energy by calculating 0.5*Tr[ Gamma^{DD}_ {ij} rho) ]. This 	!
	!   comes down to integrating the auxiliary field Gamma^{DD}(r) x 	!
	!   rho(r). Note that this function must give the sum of the proton 	!
	!   and neutron	contributions.						!
	!-----------------------------------------------------------------------!

	FUNCTION SymGDDph_get_edd(gDDph)
		DOUBLE PRECISION :: SymGDDph_get_edd
		TYPE (SymGDDph), INTENT(IN) :: gDDph

		INTEGER :: i
		DOUBLE PRECISION :: sumi, d1, d2, d3, res, Radius
		DOUBLE PRECISION, ALLOCATABLE :: Integrand(:)

		SELECT CASE (Basis)

		CASE (1)

			sumi = 0.0D0

			DO i = 1, NLag
				d1 = gDDph%dLag(NEUTRON, i) * b_0**(9.0D0/7.0D0)
				d2 = gDDph%dLag(PROTON,  i) * b_0**(9.0D0/7.0D0)
				d3 = d1 + d2
				IF (d3 .LT. 0.0D0) STOP "ERROR - Negative density in SymGDDph_get_edd"
				sumi = sumi + GaussLQ%gauss%w(i) * (d3 ** ALPHA) * d1 * d2
			END DO

			res = 3.0D0*PI * I_SALPHA3 * Gogny_t0(Gogny) * sumi

		CASE(2)

	        	ALLOCATE(Integrand(1:Npoint))

			DO i = 1, Npoint

				Integrand(i) = 0.0D0

				Radius = RadMesh(i)

				d1 = gDDph%dLag(NEUTRON, i)
				d2 = gDDph%dLag(PROTON,  i)
				d3 = d1 + d2

				Integrand(i) = EXP(LOG(d3) * ALPHA + LOG(d1) + LOG(d2) + 2.0D0*LOG(Radius))

			END DO

			! Calculating rho^(alpha) *rho_n * rho_p * r^2

			CALL simps(Integrand,Npoint,MeshStep,sumi)

			res = 4.0D0 * PI * Gogny_t0(Gogny) * sumi * (0.5D0 +  x0)

			DEALLOCATE(Integrand)

		END SELECT

		SymGDDph_get_edd = res
		RETURN
	END FUNCTION SymGDDph_get_edd

	SUBROUTINE SymGDDph_del(gDDph)
		TYPE (SymGDDph), INTENT(INOUT) :: gDDph

		DEALLOCATE(gDDph%dLag)
		RETURN
	END SUBROUTINE SymGDDph_del

END MODULE symgdd
 MODULE symgden

	USE input
	USE symd3t
	USE nucleus
	USE symfield
	USE symden
	USE symgdhf

	IMPLICIT NONE

	!---------------------------------------------------------------!
	!    SymGenDensity is a type for objects associated to:		!
	!	- a nucleus 						!
	!       - mean-field densities rho_{na,nc,la}			!
	!       - pairing abnormal tensor kappa_{na,nc,la}		!
	!---------------------------------------------------------------!

	TYPE SymGenDensity
		TYPE (NucleusType) nucleus
		TYPE (SymGenDensityHF) rho, kap
	END TYPE

 CONTAINS

	SUBROUTINE SymGenDensity_new(genden)
		TYPE (SymGenDensity), INTENT(INOUT) :: genden

		CALL SymGenDensityHF_new(genden%rho)
		CALL SymGenDensityHF_new(genden%kap)

		RETURN
	END SUBROUTINE SymGenDensity_new

	SUBROUTINE SymGenDensity_new_Nucleus(genden, N, Z)
		TYPE (SymGenDensity), INTENT(INOUT) :: genden
		INTEGER, INTENT(IN) :: N, Z
		DOUBLE PRECISION :: b

		CALL Nucleus_new(genden%nucleus, N, Z, b)

		CALL SymGenDensityHF_new(genden%rho)
		CALL SymGenDensityHF_new(genden%kap)

		RETURN
	END SUBROUTINE SymGenDensity_new_Nucleus

	SUBROUTINE SymGenDensity_new_SymDensity(genden, density)
		TYPE (SymGenDensity), INTENT(INOUT) :: genden
		TYPE (SymDensity), INTENT(INOUT) :: density

		INTEGER :: ta, la, lla

		CALL Nucleus_new_Nucleus(genden%nucleus, density%nucleus)

		! Creates densities rho and kappa

		CALL SymGenDensityHF_new(genden%rho)
		CALL SymGenDensityHF_new(genden%kap)

		IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) THEN
			DO ta = 0, 1
				DO la = 0, Lmax

					lla = 2 * (2 * la + 1)

					! Matrix genden%rho%rho(ta,a)%store has size (d,d), d = DIM(la)
					genden%rho%rho(ta, 2 * la)%store = &
						(SymD3Tensor_matrix_1(density%field%rho%p(ta), la) * DBLE(2*la) * &
						 SymD3Tensor_matrix_1(density%field%rho%a(ta), la)) / lla

					genden%kap%rho(ta, 2 * la)%store = &
						(SymD3Tensor_matrix_1(density%field%kap%p(ta), la) * DBLE(2*la) * &
						 SymD3Tensor_matrix_1(density%field%kap%a(ta), la)) / lla

					IF (HFOnly .EQ. 1) genden%kap%rho(ta, 2 * la)%store = 0.0D0

					IF (la .EQ. 0) CYCLE

					genden%rho%rho(ta, (2 * la) - 1)%store = &
						(SymD3Tensor_matrix_1(density%field%rho%p(ta), la) * DBLE(2*(la + 1)) * &
						 SymD3Tensor_matrix_1(density%field%rho%a(ta), la)) / lla

					genden%kap%rho(ta, (2 * la) - 1)%store = &
						(SymD3Tensor_matrix_1(density%field%kap%p(ta), la) * DBLE(2*(la + 1)) * &
						 SymD3Tensor_matrix_1(density%field%kap%a(ta), la)) / lla

					IF (HFOnly .EQ. 1) genden%kap%rho(ta, (2 * la) - 1)%store = 0.0D0

				END DO
			END DO
		ELSE
			DO ta = 0, 1
				DO la = 0, Lmax

					lla = 2 * (2 * la + 1)

					! Matrix genden%rho%rho(ta,a)%store has size (d,d), d = DIM(la)
					genden%rho%rho(ta, 2 * la)%store = &
						(SymD3Tensor_matrix(density%field%rho%p(ta), la) * DBLE(2*la) * &
						 SymD3Tensor_matrix(density%field%rho%a(ta), la)) / lla

					genden%kap%rho(ta, 2 * la)%store = &
						(SymD3Tensor_matrix(density%field%kap%p(ta), la) * DBLE(2*la) * &
						 SymD3Tensor_matrix(density%field%kap%a(ta), la)) / lla

					IF (HFOnly .EQ. 1) genden%kap%rho(ta, 2 * la)%store = 0.0D0

					IF (la .EQ. 0) CYCLE

					genden%rho%rho(ta, (2 * la) - 1)%store = &
						(SymD3Tensor_matrix(density%field%rho%p(ta), la) * DBLE(2*(la + 1)) * &
						 SymD3Tensor_matrix(density%field%rho%a(ta), la)) / lla

					genden%kap%rho(ta, (2 * la) - 1)%store = &
						(SymD3Tensor_matrix(density%field%kap%p(ta), la) * DBLE(2*(la + 1)) * &
						 SymD3Tensor_matrix(density%field%kap%a(ta), la)) / lla

					IF (HFOnly .EQ. 1) genden%kap%rho(ta, (2 * la) - 1)%store = 0.0D0

				END DO
			END DO
		END IF

		RETURN
	END SUBROUTINE SymGenDensity_new_SymDensity

	! Create a new density of the type SymDensity from an
	! old SymGenDensity
	!
	SUBROUTINE SymDensity_new_GenDensity(density, genden)
		TYPE (SymDensity), INTENT(INOUT) :: density
		TYPE (SymGenDensity), INTENT(IN) :: genden

		CALL Nucleus_new_Nucleus(density%nucleus, genden%nucleus)

		CALL SymHartreeFockBogolField_new(density%field)

		density%field%rho = genden%rho
		density%field%kap = genden%kap

		RETURN
	END SUBROUTINE SymDensity_new_GenDensity

	SUBROUTINE SymGenDensity_new_GammaDelta(genden, HF_gamma, HF_delta, b)
		TYPE (SymGenDensity), INTENT(INOUT) :: genden
		TYPE (SymGenDensityHF), INTENT(IN) :: HF_gamma, HF_delta

		DOUBLE PRECISION, INTENT(IN) :: b

		CALL Nucleus_new(genden%nucleus, 8, 8, b)

		CALL SymGenDensityHF_new(genden%rho)
		CALL SymGenDensityHF_new(genden%kap)

		genden%rho = HF_gamma
		genden%kap = HF_delta

		RETURN
	END SUBROUTINE SymGenDensity_new_GammaDelta

	!---------------------------------------------------------------!
	!								!
	!		ULTRA-IMPORTANT SUBROUTINE			!
	!		--------------------------			!
	!								!
	!    This subroutine calculates the values of the densities	!
	!    rho and kappa from the U and V vectors of the Bogoliubov	!
	!    transformation. We have:					!
	!								!
	!  		rho = (V*)(VT)					!
	! 		kappa = U(V+)					!
	!								!
	!---------------------------------------------------------------!

	SUBROUTINE SymGenDensity_make_Block(genden, ta, a, U, V)
		TYPE (SymGenDensity), INTENT(INOUT) :: genden
		DOUBLE PRECISION, DIMENSION(:, :), INTENT(INOUT) :: U, V !TODO(INOUT?)
		INTEGER, INTENT(IN) :: ta, a

		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: S0, S1, S2, S3
		DOUBLE PRECISION :: U0, V0

		INTEGER :: d, i

		d = DIM(a)

		ALLOCATE(S0(d, d))
		ALLOCATE(S2(d, d))

		S0 = MATMUL(V, TRANSPOSE(V))
		S2 = MATMUL(V, TRANSPOSE(U))

		IF (genden%nucleus%is_blocking(ta) .AND. (a .EQ. genden%nucleus%ia(ta))) THEN

			ALLOCATE(S1(d, d))
			ALLOCATE(S3(d, d))

			DO i = 1, d
				U0 = U(i, genden%nucleus%mu0(ta)) ! Column???
				V0 = V(i, genden%nucleus%mu0(ta))
				U(i, genden%nucleus%mu0(ta)) =   V0
				V(i, genden%nucleus%mu0(ta)) = - U0
			END DO

			S1 = MATMUL(V, TRANSPOSE(V))
			S3 = MATMUL(V, TRANSPOSE(U))
			S0 = S0 + (1.0D0 / DBLE(genden%nucleus%ja(ta) + 1)) * (S1 - S0)
			S2 = S2 + (1.0D0 / DBLE(genden%nucleus%ja(ta) + 1)) * (S3 - S2)

			DEALLOCATE(S1)
			DEALLOCATE(S3)

		END IF

		genden%rho%rho(ta, a)%store = S0
		genden%kap%rho(ta, a)%store = S2

		DEALLOCATE(S0)
		DEALLOCATE(S2)

		RETURN
	END SUBROUTINE SymGenDensity_make_Block

	SUBROUTINE SymGenDensity_make_HF(genden, ta, a, V)
		TYPE (SymGenDensity), INTENT(INOUT) :: genden
		DOUBLE PRECISION, DIMENSION(:, :), INTENT(INOUT) :: V
		INTEGER, INTENT(IN) :: ta, a

		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: S0

		INTEGER :: d

		d = DIM(a)

		ALLOCATE(S0(d, d))

		S0 = MATMUL(V, TRANSPOSE(V))

		genden%rho%rho(ta, a)%store = S0
		genden%kap%rho(ta, a)%store = 0.0d0

		DEALLOCATE(S0)

		RETURN
	END SUBROUTINE SymGenDensity_make_HF

	SUBROUTINE SymGenDensity_del(genden)
		TYPE (SymGenDensity), INTENT(INOUT) :: genden

		CALL Nucleus_del(genden%nucleus)
		CALL SymGenDensityHF_del(genden%rho)
		CALL SymGenDensityHF_del(genden%kap)
		RETURN
	END SUBROUTINE SymGenDensity_del

END MODULE symgden
 MODULE symgdhf

	USE input
	USE math
	USE symd3t
	USE symfield

	IMPLICIT NONE

	TYPE MatrixType2
		DOUBLE PRECISION, DIMENSION(:, :), POINTER :: store
	END TYPE

	!---------------------------------------------------------------------------------------!
	!    SymGenDensityHF is a type for objects associated to densities like rho_{na,nc,la}	!
	!---------------------------------------------------------------------------------------!

	TYPE SymGenDensityHF
		TYPE (MatrixType2), DIMENSION(:, :), POINTER :: rho
	END TYPE

	INTERFACE ASSIGNMENT(=)
		MODULE PROCEDURE SymGenDensityHF_assign1, SymGenDensityHF_assign2
	END INTERFACE

 CONTAINS

	SUBROUTINE SymGenDensityHF_new(gendenhf)
		TYPE (SymGenDensityHF), INTENT(INOUT) :: gendenhf

		INTEGER :: ta, a, d

		ALLOCATE(gendenhf%rho(0:1, 0:(2*Lmax)))

		DO ta = 0, 1

			DO a = 0, 2*Lmax
				d = DIM(a)
				ALLOCATE(gendenhf%rho(ta, a)%store(d, d))
				gendenhf%rho(ta, a)%store(d, d) = 0.0D0
			END DO

		END DO

		RETURN
	END SUBROUTINE SymGenDensityHF_new

	SUBROUTINE SymGenDensityHF_copy(gendenhf, HF)
		TYPE (SymGenDensityHF), INTENT(INOUT) :: gendenhf
		TYPE (SymHartreeFockField), INTENT(IN) :: HF

		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: HF_p, HF_a
		INTEGER :: ta, a, la, d, u1, u2

		DO ta = 0, 1
			DO a = 0, 2*Lmax

				la = L(a)
				d = DIM(a)

				ALLOCATE(HF_p(d, d), HF_a(d, d))

				HF_p = 0.0D0
				HF_a = 0.0D0

				DO u1 = 1, d
					DO u2 = 1, u1
						HF_p(u1, u2) = HF%p(ta)%d3tensor(la)%d2(u1, u2)
						HF_a(u1, u2) = HF%a(ta)%d3tensor(la)%d2(u1, u2)
						IF (u1 .NE. u2) THEN
							HF_p(u2, u1) = HF%p(ta)%d3tensor(la)%d2(u1, u2)
							HF_a(u2, u1) = HF%a(ta)%d3tensor(la)%d2(u1, u2)
						END IF
					END DO
				END DO

				gendenhf%rho(ta, a)%store = HF_p + (LS(a) * HF_a)

				DEALLOCATE(HF_p, HF_a)
			END DO
		END DO

		RETURN
	END SUBROUTINE SymGenDensityHF_copy

	SUBROUTINE SymGenDensityHF_assign1(gendenhf_out, gendenhf_in)
		TYPE (SymGenDensityHF), INTENT(INOUT) :: gendenhf_out
		TYPE (SymGenDensityHF), INTENT(IN) :: gendenhf_in

		INTEGER :: ta, a

		DO ta = 0, 1

			DO a = 0, 2*Lmax
				gendenhf_out%rho(ta, a)%store = gendenhf_in%rho(ta, a)%store
			END DO

		END DO

		RETURN
	END SUBROUTINE SymGenDensityHF_assign1

	SUBROUTINE SymGenDensityHF_assign2(HF, gendenhf)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF
		TYPE (SymGenDensityHF), INTENT(IN) :: gendenhf

		DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: m
		INTEGER :: ta, la, d

		DO ta = 0, 1

			DO la = 0, Lmax
				d = DIM(2*la)

				ALLOCATE(m(d, d))

				m = DBLE(2*(la + 1)) * gendenhf%rho(ta, 2 * la)%store

				IF (la .NE. 0) THEN
					m = m + (DBLE(2*la) * gendenhf%rho(ta, (2 * la) - 1)%store)
				END IF

				CALL SymD3Tensor_assign_Matrix(HF%p(ta), la, m)

				m = gendenhf%rho(ta, 2 * la)%store

				IF (la .NE. 0) THEN
					m = m - gendenhf%rho(ta, (2 * la) - 1)%store
				END IF

				CALL SymD3Tensor_assign_Matrix(HF%a(ta), la, m)

				DEALLOCATE(m)
			END DO

		END DO

		RETURN
	END SUBROUTINE SymGenDensityHF_assign2

	SUBROUTINE SymGenDensityHF_reduce(gendenhf_out, gendenhf_in, t_in)
		TYPE (SymGenDensityHF), INTENT(INOUT) :: gendenhf_out
		TYPE (SymGenDensityHF), INTENT(IN) :: gendenhf_in
		TYPE (SymD3Tensor), INTENT(IN) :: t_in

		INTEGER :: ta, a, d, la

		DO ta = 0, 1

			DO a = 0, 2*Lmax
				d = DIM(a)
				la = L(a)
				gendenhf_out%rho(ta, a)%store = gendenhf_in%rho(ta, a)%store - t_in%d3tensor(la)%d2
			END DO

		END DO

		RETURN
	END SUBROUTINE SymGenDensityHF_reduce

	SUBROUTINE SymGenDensityHF_print(gendenhf)
		TYPE (SymGenDensityHF), INTENT(IN) :: gendenhf

		INTEGER :: ta, a, d, u1, u2

		DO ta = 0, 1
			DO a = 0, 2*Lmax
				d = DIM(a)
				DO u1 = 1, d
					PRINT "(24F10.3)", (gendenhf%rho(ta, a)%store(u1, u2), u2 = 1, d)
				END DO
			END DO
		END DO

		RETURN
	END SUBROUTINE SymGenDensityHF_print

	SUBROUTINE SymGenDensityHF_del(gendenhf)
		TYPE (SymGenDensityHF), INTENT(INOUT) :: gendenhf

		INTEGER :: ta, a

		DO ta = 0, 1
			DO a = 0, 2*Lmax
				DEALLOCATE(gendenhf%rho(ta, a)%store)
			END DO
		END DO
		DEALLOCATE(gendenhf%rho)

		RETURN
	END SUBROUTINE SymGenDensityHF_del

END MODULE symgdhf
!------------------------------------------------------------------!
!								   !
!  CALCULATION OF THE KINETIC ENERGY TERM OF THE GOGNY FORCE       !
!								   !
!------------------------------------------------------------------!

 MODULE symke2b

	USE input
	USE global
	USE symd3t
	USE symfield

	IMPLICIT NONE

	TYPE SymKineticEnergy2Body
		TYPE (SymD3Tensor_SymD3Tensor) :: v11, v22
	END TYPE

	TYPE SymEk2pp
		TYPE (SymD3Tensor_SymD3Tensor) :: v11_pair, v22_pair
	END TYPE

	CHARACTER(64), PRIVATE :: filename
        INTEGER, PRIVATE :: file_desc = 16

	PRIVATE SymKineticEnergy2Body_nablaHO

 CONTAINS

	SUBROUTINE SymKineticEnergy2Body_new(vEkCMph, vEkCMpp)
		TYPE (SymKineticEnergy2Body), INTENT(INOUT) :: vEkCMph
		TYPE (SymEk2pp), INTENT(INOUT) :: vEkCMpp

		!  Create and initializes the tensors
		CALL SymD3Tensor_SymD3Tensor_new(vEkCMph%v11)
		CALL SymD3Tensor_SymD3Tensor_new(vEkCMph%v22)
		CALL SymD3Tensor_SymD3Tensor_new(vEkCMpp%v11_pair)
		CALL SymD3Tensor_SymD3Tensor_new(vEkCMpp%v22_pair)

                SELECT CASE (Basis)

		CASE (1)

			IF (N_0 < 10) THEN
				WRITE(filename, "(A,I1,A)") "data/vEkCM", N_0, "ph_HO.txt"
			ELSE
				WRITE(filename, "(A,I2,A)") "data/vEkCM", N_0, "ph_HO.txt"
			END IF

		CASE (2)

			IF (N_0 < 10) THEN
				WRITE(filename, "(A,I1,A)") "data/vEkCM", N_0, "ph_WS.txt"
			ELSE
				WRITE(filename, "(A,I2,A)") "data/vEkCM", N_0, "ph_WS.txt"
			END IF

		END SELECT

	END SUBROUTINE SymKineticEnergy2Body_new

        !---------------------------------------------------------------------------------------!
	!											!
	!  Sum over na, nb, nc, nd and la, lb to calculate the matrix elements of the type:	!
        !											!
	!                   < na la, nb lb | T | nd lb, nc la >					!
        !											!
	!  where T is the kinetic energy. Since T is a one-body operator, this can actually	!
	!  be separated into:									!
	!                        < na la | T | nd lb >< nb lb | T | nc la >			!
	!											!
	!  WARNING: FOR L, THE SUMMATION IS NOT DONE UP TO A GIVEN L_MAX, BUT OVER N_0		!
	!           FOR N, N_MAX IS EXPRESSED AS FUNCTION OF N_0 TOO				!
        !											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymKineticEnergy2Body_calculate(vEkCMph, vEkCMpp)
		TYPE (SymKineticEnergy2Body), INTENT(INOUT) :: vEkCMph
		TYPE (SymEk2pp), INTENT(INOUT) :: vEkCMpp

		INTEGER :: la, lla, na, namax, nc, lb, llb, nb, nbmax, nd
		DOUBLE PRECISION :: cuad, d1

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER :: file_error

		OPEN(file_desc, FILE=filename, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to write the results in ", filename
                END IF

                SELECT CASE (Basis)

		CASE (1)

			DO la = 0, Lmax

				lla = (2 * la) + 1 ! Degeneracy of a shell with l_a
				namax = ((N_0 - la) / 2) + 1 ! Maximum n value for the "bra"

				DO na = 1, namax
					DO nc = 1, na

			        		DO lb = 0, la

							llb = (2 * lb) + 1 ! Degeneracy of a shell with l_b
							nbmax = ((N_0 - lb) / 2) + 1 ! Maximum n value for the "ket"

							cuad = CUAD2(la, lb, 1) ! Defined in Module math.f90

							DO nb = 1, nbmax
								DO nd = 1, nb

									! Product: < na la || NABLA || nd la >< nb lb || NABLA || nc lb >
									d1 = SymKineticEnergy2Body_nablaHO(na - 1, la, nd - 1, lb) * &
									SymKineticEnergy2Body_nablaHO(nb - 1, lb, nc - 1, la)

									! Anti-symmetrization of the matrix elements by permutation:
									!          < na la || NABLA || nb la >< nd lb || NABLA || nc lb >
									IF (nb .NE. nd) THEN
										d1 = d1 + (SymKineticEnergy2Body_nablaHO(na - 1, la, nb - 1, lb) * &
												SymKineticEnergy2Body_nablaHO(nd - 1, lb, nc - 1, la))
									        d1 = d1 / 2.0D0
									END IF

									! Application of the Wigner-Eckart theorem and multiplication by the mass factor
									! to get the proper energy
									d1 = d1 / (m(1) * DBLE(lla * llb))

									! Filling in the required pointers
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v11, la, na, nc, lb, nb, nd, 0.5D0 * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v11, lb, nb, nd, la, na, nc, 0.5D0 * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v22, la, na, nc, lb, nb, nd, cuad * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v22, lb, nb, nd, la, na, nc, cuad * d1)

									! Writing the results on disk
									IF (file_error .EQ. 0) THEN
										WRITE (file_desc, "(6I3,2E24.16)", IOSTAT=file_error) &
										la, na, nc, lb, nb, nd, 0.5 * d1, cuad * d1
									END IF

									! PAIRING: < na la || NABLA || nb lb >< nc la || NABLA || nd lb >
									d1 = SymKineticEnergy2Body_nablaHO(na - 1, la, nb - 1, lb) * &
										SymKineticEnergy2Body_nablaHO(nc - 1, la, nd - 1, lb)

									! Anti-symmetrization of the matrix elements by permutation:
									!          < nb lb || NABLA || na la >< nd lb || NABLA || nc la >
									IF (nb .NE. nd) THEN
										d1 = d1 + (SymKineticEnergy2Body_nablaHO(na - 1, la, nd - 1, lb) * &
												SymKineticEnergy2Body_nablaHO(nc - 1, la, nb - 1, lb))
										d1 = d1 / 2.0D0
									END IF

									! Application of the Wigner-Eckart theorem and multiplication by the mass factor
									! to get the proper energy
									d1 = - d1 / (m(1) * DBLE(lla * llb))

									! Filling in the required pointers
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMpp%v11_pair, la, na, nc, lb, nb, nd, 0.5D0 * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMpp%v11_pair, lb, nb, nd, la, na, nc, 0.5D0 * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMpp%v22_pair, la, na, nc, lb, nb, nd, cuad * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMpp%v22_pair, lb, nb, nd, la, na, nc, cuad * d1)

								END DO
							END DO

						END DO

					END DO
				END DO

			END DO

		CASE (2)

			DO la = 0, Lmax

				lla = (2 * la) + 1 ! Degeneracy of a shell with l_a

							namax = MIN(Nmax, NmaxOfL(la))
				IF (CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

				DO na = 1, namax
					DO nc = 1, na

						DO lb = 0, la

							llb = (2 * lb) + 1 ! Degeneracy of a shell with l_b

							cuad = CUAD2(la, lb, 1) ! Defined in Module math.f90

										nbmax = MIN(Nmax, NmaxOfL(lb))
							IF (CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

							DO nb = 1, nbmax
								DO nd = 1, nb

									! Product: < na la || NABLA || nd la >< nb lb || NABLA || nc lb >
									d1 = SymKineticEnergy2Body_nabla(na, la, nd, lb) * &
										SymKineticEnergy2Body_nabla(nb, lb, nc, la)

									! Anti-symmetrization of the matrix elements by permutation:
									!          < na la || NABLA || nb la >< nd lb || NABLA || nc lb >
									IF (nb .NE. nd) THEN
										d1 = d1 + (SymKineticEnergy2Body_nabla(na, la, nb, lb) * &
												SymKineticEnergy2Body_nabla(nd, lb, nc, la))
										d1 = d1 / 2.0D0
									END IF

									! Application of the Wigner-Eckart theorem and multiplication by the mass factor
									! to get the proper energy
									d1 = d1 / (m(1) * DBLE(lla * llb))

									! Filling in the required pointers
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v11, la, na, nc, lb, nb, nd, 0.5D0 * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v11, lb, nb, nd, la, na, nc, 0.5D0 * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v22, la, na, nc, lb, nb, nd, cuad * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v22, lb, nb, nd, la, na, nc, cuad * d1)

									! Writing the results on disk
									IF (file_error .EQ. 0) THEN
										WRITE (file_desc, "(6I3,2E24.16)", IOSTAT=file_error) &
												la, na, nc, lb, nb, nd, 0.5 * d1, cuad * d1
									END IF

									! PAIRING: < na la || NABLA || nb lb >< nc la || NABLA || nd lb >
									d1 = SymKineticEnergy2Body_nabla(na, la, nb, lb) * &
										SymKineticEnergy2Body_nabla(nc, la, nd, lb)

									! Anti-symmetrization of the matrix elements by permutation:
									!          < nb lb || NABLA || na la >< nd lb || NABLA || nc la >
									IF (nb .NE. nd) THEN
										d1 = d1 + (SymKineticEnergy2Body_nabla(na, la, nd, lb) * &
												SymKineticEnergy2Body_nabla(nc, la, nb, lb))
										d1 = d1 / 2.0D0
									END IF

									! Application of the Wigner-Eckart theorem and multiplication by the mass factor
									! to get the proper energy
									d1 = - d1 / (m(1) * DBLE(lla * llb))

									! Filling in the required pointers
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMpp%v11_pair, la, na, nc, lb, nb, nd, 0.5D0 * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMpp%v11_pair, lb, nb, nd, la, na, nc, 0.5D0 * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMpp%v22_pair, la, na, nc, lb, nb, nd, cuad * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMpp%v22_pair, lb, nb, nd, la, na, nc, cuad * d1)

								END DO
							END DO

						END DO

					END DO
				END DO

			END DO


		END SELECT

		CLOSE(file_desc)
		RETURN

	END SUBROUTINE SymKineticEnergy2Body_calculate

        !---------------------------------------------------------------------------------!
	!		                       		  			          !
	! Function giving the reduced matrix elements < n_a, l_a || NABLA || n_b, l_b >   !
	! for the special case of the harmonic oscillator basis                           !
	!		                       		  			          !
	! BEWARE: Here the HO length b is implicitely assumed equal to 1, otherwise,	  !
	!         there should be a factor 1/b						  !
	!		                       		  			          !
	! Refs: Appendix C, Page 115, Appendix D, Page 127, D24				  !
	!		                       		  			          !
        !---------------------------------------------------------------------------------!

	FUNCTION SymKineticEnergy2Body_nablaHO(na, la, nb, lb)
		DOUBLE PRECISION :: SymKineticEnergy2Body_nablaHO
		INTEGER, INTENT(IN) :: na, la, nb, lb

		IF (la .EQ. (lb + 1)) THEN
			IF (nb .EQ. (na + 1)) THEN
				SymKineticEnergy2Body_nablaHO = - sq(la) * sq(nb)/b_0
			ELSE IF (nb .EQ. na) THEN
				SymKineticEnergy2Body_nablaHO = - sq(la) * sq2(na + la)/b_0
			ELSE
				SymKineticEnergy2Body_nablaHO = 0.0D0
			END IF
		ELSE IF (lb .EQ. (la + 1)) THEN
			IF (na .EQ. (nb + 1)) THEN
				SymKineticEnergy2Body_nablaHO = - sq(lb) * sq(na)/b_0
			ELSE IF (na .EQ. nb) THEN
				SymKineticEnergy2Body_nablaHO = - sq(lb) * sq2(nb + lb)/b_0
			ELSE
				SymKineticEnergy2Body_nablaHO = 0.0D0
			END IF
		ELSE
			SymKineticEnergy2Body_nablaHO = 0.0D0
		END IF

		RETURN
	END FUNCTION SymKineticEnergy2Body_nablaHO

	!
	!  Calculating the kinetic energy fields by contracting the
	!  matrix elements with the density matrix
	!
	SUBROUTINE SymKineticEnergy2Body_get_Gamma(HF_out, vEkCMph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymKineticEnergy2Body), INTENT(IN) :: vEkCMph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(0), vEkCMph%v11, HF_in%p(0))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(1), vEkCMph%v11, HF_in%p(1))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(0), vEkCMph%v22, HF_in%a(0))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(1), vEkCMph%v22, HF_in%a(1))

		RETURN
	END SUBROUTINE SymKineticEnergy2Body_get_Gamma

	SUBROUTINE SymKineticEnergy2Body_get_Delta(HF_out, vEkCMpp, P_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymEk2pp), INTENT(IN) :: vEkCMpp
		TYPE (SymHartreeFockField), INTENT(IN) :: P_in

		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(0), vEkCMpp%v11_pair, P_in%p(0))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(1), vEkCMpp%v11_pair, P_in%p(1))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(0), vEkCMpp%v22_pair, P_in%a(0))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(1), vEkCMpp%v22_pair, P_in%a(1))

		RETURN
	END SUBROUTINE SymKineticEnergy2Body_get_Delta

	SUBROUTINE SymKineticEnergy2Body_get_DeltaIso(HF_out, vEkCMpp, P_in, ta)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymEk2pp), INTENT(IN) :: vEkCMpp
		TYPE (SymHartreeFockField), INTENT(IN) :: P_in
		INTEGER, INTENT(IN) :: ta

		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(ta), vEkCMpp%v11_pair, P_in%p(ta))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(ta), vEkCMpp%v22_pair, P_in%a(ta))

        	HF_out%p(1-ta)=0.0D0
        	HF_out%a(1-ta)=0.0D0

		RETURN
	END SUBROUTINE SymKineticEnergy2Body_get_DeltaIso

	SUBROUTINE SymKineticEnergy2Body_del(vEkCMph, vEkCMpp)
		TYPE (SymKineticEnergy2Body), INTENT(INOUT) :: vEkCMph
		TYPE (SymEk2pp), INTENT(INOUT) :: vEkCMpp

		CALL SymD3Tensor_SymD3Tensor_del(vEkCMph%v11)
		CALL SymD3Tensor_SymD3Tensor_del(vEkCMph%v22)

		CALL SymD3Tensor_SymD3Tensor_del(vEkCMpp%v11_pair)
		CALL SymD3Tensor_SymD3Tensor_del(vEkCMpp%v22_pair)

		RETURN
	END SUBROUTINE SymKineticEnergy2Body_del

END MODULE symke2b
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
!----------------------------------------------------------------!
!								 !
!  CALCULATION OF THE BRINK-BOEKER TERM OF THE GOGNY FORCE       !
!								 !
!----------------------------------------------------------------!

 MODULE symvbb

	USE input
	USE angmom
	USE symd3t
	USE symfield
	USE ibb

	IMPLICIT NONE

	TYPE SymVBBph
		DOUBLE PRECISION :: b1
		TYPE (SymD3Tensor_SymD3Tensor) :: v_local_same_part, v_local_diff_part, &
						  v1_exch_same_part, v1_exch_diff_part, &
						  v2_exch_same_part, v2_exch_diff_part
		TYPE (SymD3Tensor_SymD3Tensor) :: v1_same_part, v1_diff_part, &
						  v2_same_part, v2_diff_part
		CHARACTER(LEN = 64) filename
	END TYPE

	TYPE SymVBBpp
		DOUBLE PRECISION :: b1
		TYPE (SymD3Tensor_SymD3Tensor) :: v1_pair, v2_pair
		CHARACTER(LEN = 64) :: filename
	END TYPE

	! Range of the Brink-Boker term

	INTEGER :: Lmax_read, Nmax_read, i_gaus_max = 1

	DOUBLE PRECISION, DIMENSION(0:1) :: mu

 CONTAINS

        !---------------------------------------------------------------------------------------!
	!											!
        !   Creates a new object of the type SymVBBph, that is in fact a collection of 		!
	!   "3D tensors" arrays SymD3Tensor_SymD3Tensor, as defined in module symd3t.f90	!
	!   Here we create each of the corresponding parts of the type SymVBBph and allocate 	!
	!   the appropriate amount of memory							!
	!   											!
	!   This tensor is for the Brink-Boker matrix elements in the particle-hole channel	!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVBBph_new(vBBph, b)
		TYPE (SymVBBph), INTENT(INOUT) :: vBBph
		DOUBLE PRECISION, INTENT(IN) :: b

		! Create the tensors

		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v_local_same_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v_local_diff_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v1_exch_same_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v1_exch_diff_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v2_exch_same_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v2_exch_diff_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v1_same_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v1_diff_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v2_same_part)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v2_diff_part)

		vBBph%b1 = b

                SELECT CASE (Basis)

		CASE (1)

			IF (N_0 < 10) THEN
				WRITE(vBBph%filename, "(A,I1,A)") "data/vBB", N_0, "ph_HO.txt"
			ELSE
				WRITE(vBBph%filename, "(A,I2,A)") "data/vBB", N_0, "ph_HO.txt"
			END IF

		CASE (2)

			IF (N_0 < 10) THEN
				WRITE(vBBph%filename, "(A,I1,A)") "data/vBB", N_0, "ph_WS.txt"
			ELSE
				WRITE(vBBph%filename, "(A,I2,A)") "data/vBB", N_0, "ph_WS.txt"
			END IF

		END SELECT

		RETURN
	END SUBROUTINE SymVBBph_new

        !---------------------------------------------------------------------------------------!
	!											!
	!   Subroutine calculating the matrix elements of the Brink-Boker term. We distinguish	!
	!   2 special cases, either the calculation is done in the harmonicoscillator basis 	!
	!   (basis = 1) or in a general spherical basis (basis = 2)				!
	!											!
	!   BEWARE: Here the HO length b is implicitely assumed equal to 1, otherwise,there 	!
	!           should be a factor 1/b							!
	!											!
	!   Refs: Appendix F									!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVBBph_calculate(vBBph, read_BB, Lmin)
		TYPE (SymVBBph), INTENT(INOUT) :: vBBph
		INTEGER, INTENT(INOUT) :: Lmin
		LOGICAL, INTENT(IN) :: read_BB

		INTEGER :: i, icount, la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER :: k, kmin, kmax, Lold, n_reg
		DOUBLE PRECISION :: sumi_sp, sumi_dp, total_ibb
		DOUBLE PRECISION :: sum1_sp, sum1_dp, sum2_sp, sum2_dp
		DOUBLE PRECISION :: dp, sp, dp_direct, sp_direct, dp_plus, sp_plus, dp_minu, sp_minu, range
		DOUBLE PRECISION :: tres_j_cuad, cuad, pi

		DOUBLE PRECISION, DIMENSION(0:1) :: x
		DOUBLE PRECISION, ALLOCATABLE :: coeff_general(:,:,:), coeff_aux(:,:,:)

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER :: file_error

		! If the matrix elements were already calculated, we read them and avoid thereby recalculating them
		! The tensors are then updated with the values read from tape.
		! In case of the WS basis, we (smartly) read what was already calculated and eventually calculate
		! only those elements that are missing.

		IF (Lmin .EQ. 0) THEN

			IF (read_BB) THEN
				CALL SymVBBph_update(vBBph)
				RETURN
			ELSE
				IF (SymVBBph_read(vBBph)) THEN
					CALL SymVBBph_update(vBBph)
					RETURN
				END IF
			END IF

			OPEN(file_desc, FILE=vBBph%filename, ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				PRINT *, "*** Attention: Impossible to write in ", vBBph%filename
			ELSE
				WRITE (file_desc, FMT="(E24.16)", IOSTAT=file_error) vBBph%b1
			END IF

		ELSE

			Lold = Lmax
			Lmax = Lold - 1

			IF (Lmax < 10) THEN
				WRITE(vBBph%filename, "(A,I1,A)") "data/vBB", Lmax, "ph_WS.txt"
			ELSE
				WRITE(vBBph%filename, "(A,I2,A)") "data/vBB", Lmax, "ph_WS.txt"
			END IF

			IF (SymVBBph_read(vBBph)) write(*,'("Read file....")')

			Lmax = Lold

			IF (Lmax < 10) THEN
				WRITE(vBBph%filename, "(A,I1,A)") "data/vBB", Lmax, "ph_WS.txt"
			ELSE
				WRITE(vBBph%filename, "(A,I2,A)") "data/vBB", Lmax, "ph_WS.txt"
			END IF

		END IF

		! For each isospin, we calculate the constant numerical factors in front of the matrix elements
		! names with the suffix "_same_part" refer to the terms with the delta_(ta,tb), names with the
		! suffix "_diff_part" refer to the terms with the other terms:
		!
		! Refs.: Page 130-131, Definition of v1eBB and v2eBB

		IF (regularized_Gaussian) THEN

			WRITE(6,'("Subtituting Gaussian with regularized deltas")')

			! Index 1: isospin, index 2: Gogny component, index 3: regularization index (power n of r^n)
			n_reg_max = 2; i_gaus_max = 0
			ALLOCATE(coeff_general(0:i_gaus_max,1:4,0:n_reg_max)); coeff_general = 0.0D0
			IF (test_regularization) THEN
				ALLOCATE(coeff_aux(0:i_gaus_max,1:4,0:n_reg_max)); coeff_aux = 0.0D0
			END IF

			pi = 4.0D0*ATAN(1.0D0)

			mu(0) = range1
			mu(1) = range1

			DO i = 0, i_gaus_max

			        x(i) = mu(i) / vBBph%b1

			        coeff_general(i, 1, 0) = -2358.431874D0 * I_4PI
			        coeff_general(i, 2, 0) =  2030.262376D0 * I_4PI
			        coeff_general(i, 3, 0) = -2868.158544D0 * I_4PI
			        coeff_general(i, 4, 0) =  1874.817737D0 * I_4PI

			        coeff_general(i, 1, 2) =  1359.897114D0 * I_4PI
			        coeff_general(i, 2, 2) = -1339.802973D0 * I_4PI
			        coeff_general(i, 3, 2) =  1856.769960D0 * I_4PI
			        coeff_general(i, 4, 2) = -1373.442870D0 * I_4PI

				!IF (test_regularization) THEN
                                !
				!	coeff_aux(i, 1, 0) = (mu(0)*SQRT(pi))**3 * (coeff_general(0, 1, 0) + 1.5D0*mu(0)*mu(0)*coeff_general(0, 1, 2))
				!	coeff_aux(i, 2, 0) = (mu(0)*SQRT(pi))**3 * (coeff_general(0, 2, 0) + 1.5D0*mu(0)*mu(0)*coeff_general(0, 2, 2))
				!	coeff_aux(i, 3, 0) = (mu(0)*SQRT(pi))**3 * (coeff_general(0, 3, 0) + 1.5D0*mu(0)*mu(0)*coeff_general(0, 3, 2))
				!	coeff_aux(i, 4, 0) = (mu(0)*SQRT(pi))**3 * (coeff_general(0, 4, 0) + 1.5D0*mu(0)*mu(0)*coeff_general(0, 4, 2))
                                !
				!	coeff_aux(i, 1, 2) = 0.5D0*mu(0)**3 * (mu(0)*SQRT(pi))**3 * coeff_general(0, 1, 2)
				!	coeff_aux(i, 2, 2) = 0.5D0*mu(0)**3 * (mu(0)*SQRT(pi))**3 * coeff_general(0, 2, 2)
				!	coeff_aux(i, 3, 2) = 0.5D0*mu(0)**3 * (mu(0)*SQRT(pi))**3 * coeff_general(0, 3, 2)
				!	coeff_aux(i, 4, 2) = 0.5D0*mu(0)**3 * (mu(0)*SQRT(pi))**3 * coeff_general(0, 4, 2)
                                !
				!END IF

 			END DO
		ELSE
			! Index 1: isospin, index 2: Gogny component, index 3: regularization index (power n of r^n)
			n_reg_max = 0; i_gaus_max = 1
			ALLOCATE(coeff_general(0:i_gaus_max,1:4,0:n_reg_max)); coeff_general = 0.0D0

			mu(0) = range_gaussian(0)
			mu(1) = range_gaussian(1)

			DO i = 0, i_gaus_max
				x(i) = mu(i) / vBBph%b1
				coeff_general(i, 1, 0) = I_4PI * Gogny_W(i, Gogny)
				coeff_general(i, 2, 0) = I_4PI * Gogny_B(i, Gogny)
				coeff_general(i, 3, 0) = I_4PI * Gogny_H(i, Gogny)
				coeff_general(i, 4, 0) = I_4PI * Gogny_M(i, Gogny)
			END DO
		END IF

                SELECT CASE (Basis)

		CASE (1)

			! Calculation of the matrix elements v1eBB and v2eBB as defined in Page 131  ---  CHECKED AND OK

			PRINT *, "Brink-Boker terms: Particle-Hole Channel - Harmonic Oscillator Basis"
			! Los mas altos son los mas problematicos
			DO la = Lmin, Lmax
				namax = ((N_0 - la) / 2) + 1
				DO na = 1, namax
					DO nc = 1, na

						DO lb = 0, la
							nbmax = ((N_0 - lb) / 2) + 1
							DO nb = 1, nbmax
								DO nd = 1, nb

									! Calculation of the term k=0 in the multipole expansion over Wk(r1,r2). This actually
									! corresponds to the local term

									sumi_sp = 0.0D0
									sumi_dp = 0.0D0

									DO n_reg = 0, n_reg_max, 2
									        DO i = 0, i_gaus_max
									        	total_ibb = IBrinkBookerHO(na - 1, la, nb - 1, lb, nc - 1, la, nd - 1, lb, 0, x(i))
									        	IF (nb .NE. nd) THEN
									        		total_ibb = total_ibb + IBrinkBookerHO(na - 1, la, nd - 1, lb, nc - 1, la, nb - 1, lb, 0, x(i))
									        		total_ibb = total_ibb / 2.0D0
									        	END IF
									        	sumi_dp = sumi_dp + (coeff_general(i, 1, n_reg) + 0.5D0*coeff_general(i, 2, n_reg)) * total_ibb
									        	sumi_sp = sumi_sp + (coeff_general(i, 1, n_reg) + 0.5D0*coeff_general(i, 2, n_reg) &
									        			  -  coeff_general(i, 3, n_reg) - 0.5D0*coeff_general(i, 4, n_reg))* total_ibb
									        END DO
									END DO

									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_same_part, la, na, nc, lb, nb, nd, sumi_sp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_same_part, lb, nb, nd, la, na, nc, sumi_sp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_diff_part, la, na, nc, lb, nb, nd, sumi_dp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_diff_part, lb, nb, nd, la, na, nc, sumi_dp)

									! Summations over the multipole Wk(r1, r2) (Appendix B, formula above B2, page 113) are restricted
									! by the angular part. Only the terms with ABS(la - lb) <= k <= la + lb give non-zero terms

									kmin = ABS(la - lb)
									kmax = la + lb

									sum1_sp = 0.0D0
									sum1_dp = 0.0D0
									sum2_sp = 0.0D0
									sum2_dp = 0.0D0

									DO k = kmin, kmax, 2
										tres_j_cuad = (2*k + 1) * (ThreeJSymbols_get(2*la, 2*k, 2*lb) ** 2)
										cuad = CUAD2(la, lb, k)
										DO n_reg = 0, n_reg_max, 2
										        DO i = 0, i_gaus_max
										        	total_ibb = IBrinkBookerHO(na - 1, la, nb - 1, lb, nd - 1, lb, nc - 1, la, k, x(i))
										        	IF(nb .NE. nd) THEN
										        		total_ibb = total_ibb + IBrinkBookerHO(na - 1, la, nd - 1, lb, nb - 1, lb, nc - 1, la, k, x(i))
										        		total_ibb = total_ibb / 2.0D0
										        	END IF
										        	sum1_dp = sum1_dp + (coeff_general(i, 4, n_reg) + 0.5D0*coeff_general(i, 3, n_reg))* tres_j_cuad * total_ibb
										        	sum1_sp = sum1_sp + (coeff_general(i, 4, n_reg) + 0.5D0*coeff_general(i, 3, n_reg) &
										        			  -  coeff_general(i, 2, n_reg) - 0.5D0*coeff_general(i, 1, n_reg))* tres_j_cuad * total_ibb
										        	sum2_dp = sum2_dp + (coeff_general(i, 3, n_reg))* cuad * tres_j_cuad * total_ibb
										        	sum2_sp = sum2_sp + (coeff_general(i, 3, n_reg) - coeff_general(i, 1, n_reg))* cuad * tres_j_cuad * total_ibb
										        END DO
										END DO
									END DO

									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_same_part, la, na, nc, lb, nb, nd, sum1_sp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_same_part, lb, nb, nd, la, na, nc, sum1_sp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_diff_part, la, na, nc, lb, nb, nd, sum1_dp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_diff_part, lb, nb, nd, la, na, nc, sum1_dp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_same_part, la, na, nc, lb, nb, nd, sum2_sp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_same_part, lb, nb, nd, la, na, nc, sum2_sp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_diff_part, la, na, nc, lb, nb, nd, sum2_dp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_diff_part, lb, nb, nd, la, na, nc, sum2_dp)

									IF (file_error .EQ. 0 .AND. Lmin .EQ. 0) THEN
										WRITE (file_desc, FMT="(6I3,6E24.16)", IOSTAT=file_error) &
											la, na, nc, lb, nb, nd, sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp
									END IF
								END DO
							END DO
						END DO
					END DO
				END DO
				PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
			END DO

		CASE(2)

			PRINT *, "Brink-Boker terms: Particle-Hole Channel - General Basis"

			DO la = Lmin, Lmax

			 	icount = 0
							namax = MIN(Nmax, NmaxOfL(la))
				IF (CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

				DO na = 1, namax
					DO nc = 1, na

						DO lb = 0, la
										nbmax = MIN(Nmax, NmaxOfL(lb))
							IF (CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

							DO nb = 1, nbmax
								DO nd = 1, nb

									sumi_sp = 0.0D0
									sumi_dp = 0.0D0
									DO n_reg = 0, n_reg_max, 2
										DO i = 0, i_gaus_max
											total_ibb = IBrinkBooker(na, la, nb, lb, nc, la, nd, lb, 0, i, n_reg)
											IF (nb .NE. nd) THEN
												total_ibb = total_ibb + IBrinkBooker(na, la, nd, lb, nc, la, nb, lb, 0, i, n_reg)
												total_ibb = total_ibb / 2.0D0
											END IF
											sumi_dp = sumi_dp + (coeff_general(i, 1, n_reg) + 0.5D0*coeff_general(i, 2, n_reg)) * total_ibb
											sumi_sp = sumi_sp + (coeff_general(i, 1, n_reg) + 0.5D0*coeff_general(i, 2, n_reg) &
													  -  coeff_general(i, 3, n_reg) - 0.5D0*coeff_general(i, 4, n_reg))* total_ibb
										END DO
									END DO

									IF (test_regularization) THEN
										! Compute derivatives with respect to regularization width a (=mu(i))
										dp_plus = IBrinkBooker(na, la, nb, lb, nc, la, nd, lb, 0, 0, 0)
										dp_minu = IBrinkBooker(na, la, nb, lb, nc, la, nd, lb, 0, 1, 0)
										dp = 0.5D0*(dp_plus - dp_minu)/delta_a * (2.0D0/mu(0)**3)
										! Compute second-order term directly
										dp_direct = IBrinkBooker(na, la, nb, lb, nc, la, nd, lb, 0, 2, 2)
										WRITE(6,'("dp = ",e24.16," dpd = ", e24.16," ratio = ",e24.16)') &
											           dp, dp_direct, dp/dp_direct
									END IF

									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_same_part, la, na, nc, lb, nb, nd, sumi_sp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_same_part, lb, nb, nd, la, na, nc, sumi_sp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_diff_part, la, na, nc, lb, nb, nd, sumi_dp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_diff_part, lb, nb, nd, la, na, nc, sumi_dp)

									kmin = ABS(la - lb)
									kmax = la + lb
									sum1_sp = 0.0D0
									sum1_dp = 0.0D0
									sum2_sp = 0.0D0
									sum2_dp = 0.0D0

									DO n_reg = 0, n_reg_max, 2
										DO k = kmin, kmax, 2
											tres_j_cuad = (2*k + 1.0) * (ThreeJSymbols_get(2*la, 2*k, 2*lb) ** 2)
											cuad = CUAD2(la, lb, k)
											DO i = 0, i_gaus_max
												total_ibb = IBrinkBooker(na, la, nb, lb, nd, lb, nc, la, k, i, n_reg)
												IF (nb .NE. nd) THEN
													total_ibb = total_ibb + IBrinkBooker(na, la, nd, lb, nb, lb, nc, la, k, i, n_reg)
													total_ibb = total_ibb / 2.0D0
												END IF
												sum1_dp = sum1_dp + (coeff_general(i, 4, n_reg) + 0.5D0*coeff_general(i, 3, n_reg))* tres_j_cuad * total_ibb
												sum1_sp = sum1_sp + (coeff_general(i, 4, n_reg) + 0.5D0*coeff_general(i, 3, n_reg) &
														  -  coeff_general(i, 2, n_reg) - 0.5D0*coeff_general(i, 1, n_reg))* tres_j_cuad * total_ibb
												sum2_dp = sum2_dp + (coeff_general(i, 3, n_reg))* cuad * tres_j_cuad * total_ibb
												sum2_sp = sum2_sp + (coeff_general(i, 3, n_reg) - coeff_general(i, 1, n_reg))* cuad * tres_j_cuad * total_ibb
											END DO
										END DO
									END DO

									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_same_part, la, na, nc, lb, nb, nd, sum1_sp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_same_part, lb, nb, nd, la, na, nc, sum1_sp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_diff_part, la, na, nc, lb, nb, nd, sum1_dp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_diff_part, lb, nb, nd, la, na, nc, sum1_dp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_same_part, la, na, nc, lb, nb, nd, sum2_sp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_same_part, lb, nb, nd, la, na, nc, sum2_sp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_diff_part, la, na, nc, lb, nb, nd, sum2_dp)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_diff_part, lb, nb, nd, la, na, nc, sum2_dp)

									IF (file_error .EQ. 0 .AND. Lmin .EQ. 0) THEN
										WRITE (file_desc, FMT="(6I3,6E24.16)", IOSTAT=file_error) &
											la, na, nc, lb, nb, nd, sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp
									END IF
								END DO
							END DO
						END DO
					END DO
				END DO
				WRITE(*,'("Number of re-used integrals = ",I16)') icount
				WRITE(*,'("Completion rate: ",I3,"%")') INT(100 * (la + 1) / (N_0 + 1))
			END DO

		END SELECT

		DEALLOCATE(coeff_general)

		! In case we want to optimize our calculation, we write the full matrix elements after
		! having calculated only the few ones required

		IF (Optimization .EQ. 1 .AND. Lmin .NE. 0) THEN

			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) THEN
				WRITE(*,'("Incompatible options: No optimization for HO case (or compatible)")')
				STOP 'Error in SymVBBph_calculate - No optimization pssible'
			END IF

			CALL SymVBBph_write(vBBph)

		END IF

		CALL SymVBBph_update(vBBph)

		CLOSE(file_desc)

		RETURN

	CONTAINS

		!---------------------------------------------------------------------------------------!
		!   Subroutine updating the tensor fields of the Brink-Boeker term in the particle-	!
		!   hole channel after reading the matrix elements from tape				!
        	!---------------------------------------------------------------------------------------!

		SUBROUTINE SymVBBph_update(vBBph)
			TYPE (SymVBBph), INTENT(INOUT) :: vBBph

			CALL SymD3Tensor_SymD3Tensor_add(vBBph%v1_same_part, vBBph%v_local_same_part, vBBph%v1_exch_same_part)
			CALL SymD3Tensor_SymD3Tensor_add(vBBph%v1_diff_part, vBBph%v_local_diff_part, vBBph%v1_exch_diff_part)
			CALL SymD3Tensor_SymD3Tensor_add(vBBph%v2_same_part, vBBph%v_local_same_part, vBBph%v2_exch_same_part)
			CALL SymD3Tensor_SymD3Tensor_add(vBBph%v2_diff_part, vBBph%v_local_diff_part, vBBph%v2_exch_diff_part)

			RETURN
		END SUBROUTINE SymVBBph_update

	END SUBROUTINE SymVBBph_calculate

        !---------------------------------------------------------------------------------------!
	!											!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVBBph_write(vBBph)
		TYPE (SymVBBph), INTENT(INOUT) :: vBBph

		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd
		DOUBLE PRECISION :: sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		OPEN(file_desc, FILE=vBBph%filename, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to write in ", vBBph%filename
		ELSE
			WRITE (file_desc, FMT="(E24.16)", IOSTAT=file_error) vBBph%b1
		END IF

		DO la = 0, Lmax

						namax = MIN(Nmax, NmaxOfL(la))
			IF (CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

			DO na = 1, namax
				DO nc = 1, na

					DO lb = 0, la
									nbmax = MIN(Nmax, NmaxOfL(lb))
						IF (CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

						DO nb = 1, nbmax
							DO nd = 1, nb

								CALL SymD3Tensor_SymD3Tensor_get(sumi_sp, vBBph%v_local_same_part, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sumi_dp, vBBph%v_local_diff_part, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum1_sp, vBBph%v1_exch_same_part, lb, nb, nd, la, na, nc)
								CALL SymD3Tensor_SymD3Tensor_get(sum1_dp, vBBph%v1_exch_diff_part, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum2_sp, vBBph%v2_exch_same_part, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum2_dp, vBBph%v2_exch_diff_part, lb, nb, nd, la, na, nc)

								IF (file_error .EQ. 0) THEN
									WRITE (file_desc, FMT="(6I3,6E24.16)", IOSTAT=file_error) &
									la, na, nc, lb, nb, nd, sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp
								END IF
							END DO
						END DO
					END DO
				END DO
			END DO

		END DO

		CLOSE(file_desc)

		RETURN

	END SUBROUTINE SymVBBph_write

	SUBROUTINE SymVBBph_get_Gamma(HF_out, vBBph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVBBph), INTENT(IN) :: vBBph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		TYPE (SymD3Tensor) same, diff

		CALL SymD3Tensor_new(same)
		CALL SymD3Tensor_new(diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v1_same_part, HF_in%p(0))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v1_diff_part, HF_in%p(1))
		CALL SymD3Tensor_add(HF_out%p(0), same, diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v1_same_part, HF_in%p(1))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v1_diff_part, HF_in%p(0))
		CALL SymD3Tensor_add(HF_out%p(1), same, diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v2_exch_same_part, HF_in%a(0))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v2_exch_diff_part, HF_in%a(1))
		CALL SymD3Tensor_add(HF_out%a(0), same, diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v2_exch_same_part, HF_in%a(1))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v2_exch_diff_part, HF_in%a(0))
		CALL SymD3Tensor_add(HF_out%a(1), same, diff)

		CALL SymD3Tensor_del(same)
		CALL SymD3Tensor_del(diff)

		RETURN
	END SUBROUTINE SymVBBph_get_Gamma

	SUBROUTINE SymVBBph_get_LocalGamma(HF_out, vBBph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVBBph), INTENT(IN) :: vBBph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		TYPE (SymD3Tensor) same, diff

		CALL SymD3Tensor_new(same)
		CALL SymD3Tensor_new(diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v_local_same_part, HF_in%p(0))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v_local_diff_part, HF_in%p(1))
		CALL SymD3Tensor_add(HF_out%p(0), same, diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v_local_same_part, HF_in%p(1))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v_local_diff_part, HF_in%p(0))
		CALL SymD3Tensor_add(HF_out%p(1), same, diff)

		HF_out%a(0) = 0.0D0
		HF_out%a(1) = 0.0D0

		CALL SymD3Tensor_del(same)
		CALL SymD3Tensor_del(diff)

		RETURN
	END SUBROUTINE SymVBBph_get_LocalGamma

	SUBROUTINE SymVBBph_get_ExchangeGamma(HF_out, vBBph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVBBph), INTENT(IN) :: vBBph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		TYPE (SymD3Tensor) same, diff

		CALL SymD3Tensor_new(same)
		CALL SymD3Tensor_new(diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v1_exch_same_part, HF_in%p(0))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v1_exch_diff_part, HF_in%p(1))
		CALL SymD3Tensor_add(HF_out%p(0), same, diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v1_exch_same_part, HF_in%p(1))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v1_exch_diff_part, HF_in%p(0))
		CALL SymD3Tensor_add(HF_out%p(1), same, diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v2_exch_same_part, HF_in%a(0))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v2_exch_diff_part, HF_in%a(1))
		CALL SymD3Tensor_add(HF_out%a(0), same, diff)

		CALL SymD3Tensor_SymD3Tensor_product(same, vBBph%v2_exch_same_part, HF_in%a(1))
		CALL SymD3Tensor_SymD3Tensor_product(diff, vBBph%v2_exch_diff_part, HF_in%a(0))
		CALL SymD3Tensor_add(HF_out%a(1), same, diff)

		CALL SymD3Tensor_del(same)
		CALL SymD3Tensor_del(diff)

		RETURN
	END SUBROUTINE SymVBBph_get_ExchangeGamma

        !---------------------------------------------------------------------------------------!
	!											!
	!   Subroutine reading the matrix elements of the Brink-Boker term in the particle-hole	!
	!   channel. This subroutine attempts reading a file containing previously calculated	!
	!   matrix elements. It is used to decide whether or not to recalculate eveything	!
	!											!
        !---------------------------------------------------------------------------------------!

	FUNCTION SymVBBph_read(vBBph)
		LOGICAL SymVBBph_read
		TYPE (SymVBBph), INTENT(INOUT) :: vBBph

		DOUBLE PRECISION b1
		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER i1, i2, i3, i4, i5, i6
		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error
		DOUBLE PRECISION sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp

		OPEN(file_desc, FILE=vBBph%filename, ACTION="READ", IOSTAT=file_error)

		IF (file_error .NE. 0) THEN
			PRINT *, "Impossible to read file: ", vBBph%filename
			SymVBBph_read = .FALSE.
			RETURN
		END IF

		READ (file_desc, FMT=*, IOSTAT=file_error) b1

		IF ((file_error .NE. 0) .OR. (ABS(b1 - vBBph%b1) .GE. 1.E-14)) THEN
			CLOSE(file_desc)
			SymVBBph_read = .FALSE.
			RETURN
		END IF

		DO la = 0, Lmax
								namax = MIN(Nmax, NmaxOfL(la))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

			DO na = 1, namax
				DO nc = 1, na

					DO lb = 0, la
											nbmax = MIN(Nmax, NmaxOfL(lb))
						IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

						DO nb = 1, nbmax
							DO nd = 1, nb

								READ (file_desc, FMT=*, IOSTAT=file_error) &
									i1, i2, i3, i4, i5, i6, sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp

								IF ((file_error .NE. 0) .OR. &
									(la .NE. i1) .OR. (na .NE. i2) .OR. (nc .NE. i3) .OR. &
									(lb .NE. i4) .OR. (nb .NE. i5) .OR. (nd .NE. i6)) THEN
									PRINT *, "Invalid information in file: ", vBBph%filename
									CLOSE(file_desc)
									SymVBBph_read = .FALSE.
									RETURN
								END IF

								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_same_part, la, na, nc, lb, nb, nd, sumi_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_same_part, lb, nb, nd, la, na, nc, sumi_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_diff_part, la, na, nc, lb, nb, nd, sumi_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_diff_part, lb, nb, nd, la, na, nc, sumi_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_same_part, la, na, nc, lb, nb, nd, sum1_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_same_part, lb, nb, nd, la, na, nc, sum1_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_diff_part, la, na, nc, lb, nb, nd, sum1_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_diff_part, lb, nb, nd, la, na, nc, sum1_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_same_part, la, na, nc, lb, nb, nd, sum2_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_same_part, lb, nb, nd, la, na, nc, sum2_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_diff_part, la, na, nc, lb, nb, nd, sum2_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_diff_part, lb, nb, nd, la, na, nc, sum2_dp)

							END DO
						END DO
					END DO
				END DO
			END DO
		END DO

		CLOSE(file_desc)

		SymVBBph_read = .TRUE.

		RETURN
	END FUNCTION SymVBBph_read

       !---------------------------------------------------------------------------------------!
	!											!
	!   Subroutine deleting the tensors corresponding to the matrix elements of the Brink-	!
	!   Boker in the particle-hole channel.							!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVBBph_del(vBBph)
		TYPE (SymVBBph), INTENT(INOUT) :: vBBph

		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v_local_same_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v_local_diff_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v1_exch_same_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v1_exch_diff_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v2_exch_same_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v2_exch_diff_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v1_same_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v1_diff_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v2_same_part)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v2_diff_part)

		RETURN
	END SUBROUTINE SymVBBph_del

        !---------------------------------------------------------------------------------------!
	!											!
        !   Creates a new object of the type SymVBBph, that is in fact a collection of 		!
	!   "3D tensors" arrays SymD3Tensor_SymD3Tensor, as defined in module symd3t.f90	!
	!   Here we create each of the corresponding parts of the type SymVBBph and allocate 	!
	!   the appropriate amount of memory							!
	!   											!
	!   This tensor is for the Brink-Boker matrix elements in the particle-particle channel	!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVBBpp_new(vBBpp, b)
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp
		DOUBLE PRECISION, INTENT(IN) :: b

		CALL SymD3Tensor_SymD3Tensor_new(vBBpp%v1_pair)
		CALL SymD3Tensor_SymD3Tensor_new(vBBpp%v2_pair)

		vBBpp%b1 = b

                SELECT CASE (Basis)

		CASE (1)

			IF (N_0 < 10) THEN
				WRITE(vBBpp%filename, "(A,I1,A)") "data/vBB", N_0, "pp_HO.txt"
			ELSE
				WRITE(vBBpp%filename, "(A,I2,A)") "data/vBB", N_0, "pp_HO.txt"
			END IF

		CASE(2)

			IF (N_0 < 10) THEN
				WRITE(vBBpp%filename, "(A,I1,A)") "data/vBB", N_0, "pp_WS.txt"
			ELSE
				WRITE(vBBpp%filename, "(A,I2,A)") "data/vBB", N_0, "pp_WS.txt"
			END IF

		END SELECT

		RETURN
	END SUBROUTINE SymVBBpp_new

        !---------------------------------------------------------------------------------------!
	!											!
	!   Subroutine calculating the matrix elements of the Brink-Boker term in the particle-	!
	!   particle channel. We distinguish 2 special cases, either the calculation is done 	!
	!   in the harmonic oscillator basis (basis = 1) or in a general spherical basis 	!
	!   (basis = 2)										!
	!											!
	!   BEWARE: Here the HO length b is implicitely assumed equal to 1, otherwise,there 	!
	!           should be a factor 1/b							!
	!											!
	!   Refs: Appendix F, Page 138								!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVBBpp_calculate(vBBpp, Read_BBpp, Lmin)
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp
		INTEGER, INTENT(INOUT) :: Lmin
		LOGICAL, INTENT(IN) :: Read_BBpp

		INTEGER i, icount
		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER k, kmin, kmax, Lold
		DOUBLE PRECISION sum1, sum2
		DOUBLE PRECISION tres_j_cuad, cuad, total_ibb

		DOUBLE PRECISION, DIMENSION(0:1) :: x, coef_pair1, coef_pair2

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		! If the matrix elements were already calculated, we read them and avoid thereby recalculating them
		! The tensors are then updated with the values read from tape

		IF (Lmin .EQ. 0) THEN

			IF (Read_BBpp) RETURN

			IF (SymVBBpp_read(vBBpp)) RETURN

			OPEN(file_desc, FILE=vBBpp%filename, ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				PRINT *, "*** Attention: Impossible to write in ", vBBpp%filename
			ELSE
				WRITE (file_desc, FMT="(E24.16)", IOSTAT=file_error) vBBpp%b1
			END IF

		ELSE

			Lold = Lmax
			Lmax = Lold - 1

			IF (Lmax < 10) THEN
				WRITE(vBBpp%filename, "(A,I1,A)") "data/vBB", Lmax, "pp_WS.txt"
			ELSE
				WRITE(vBBpp%filename, "(A,I2,A)") "data/vBB", Lmax, "pp_WS.txt"
			END IF

			IF (SymVBBpp_read(vBBpp)) write(*,'("Read file....")')

			Lmax = Lold

			IF (Lmax < 10) THEN
				WRITE(vBBpp%filename, "(A,I1,A)") "data/vBB", Lmax, "pp_WS.txt"
			ELSE
				WRITE(vBBpp%filename, "(A,I2,A)") "data/vBB", Lmax, "pp_WS.txt"
			END IF

		END IF

		! For each isospin, we calculate the constant numerical factors in front of the matrix elements
		! names with the suffix "_same_part" refer to the terms with the delta_(ta,tb), names with the
		! suffix "_diff_part" refer to the terms with the other terms:
		!
		! Refs.: Page 132, top of th page, definition of v1pBB and v2pBB  	---  CHECKED AND OK
		!

		DO i = 0, 1
			x(i) = mu(i) / vBBpp%b1
			coef_pair1(i) = 0.5D0 * I_4PI * (Gogny_W(i, Gogny) - Gogny_B(i, Gogny) - Gogny_H(i, Gogny) + Gogny_M(i, Gogny))
			coef_pair2(i) =         I_4PI * (Gogny_W(i, Gogny) + Gogny_B(i, Gogny) - Gogny_H(i, Gogny) - Gogny_M(i, Gogny))
		END DO

                SELECT CASE (Basis)

		CASE (1)

			! Calculation of the matrix elements v1pBB and v2pBB as defined in Page 132  ---  CHECKED AND OK

			PRINT *, "Brink-Boker terms: Particle-Particle Channel - Harmonic Oscillator Basis"
			DO la = Lmin, Lmax
				namax = ((N_0 - la) / 2) + 1
				DO na = 1, namax
					DO nc = 1, na
						DO lb = 0, la
							nbmax = ((N_0 - lb) / 2) + 1
							DO nb = 1, nbmax
								DO nd = 1, nb

									kmin = ABS(la - lb)
									kmax = la + lb

									sum1 = 0.0D0
									sum2 = 0.0D0

									DO k = kmin, kmax, 2
										tres_j_cuad = (2*k + 1) * (ThreeJSymbols_get(2*la, 2*k, 2*lb) ** 2)
										cuad = CUAD2(la, lb, k)
										DO i = 0, i_gaus_max
											total_ibb = IBrinkBookerHO(na - 1, la, nc - 1, la, nb - 1, lb, nd - 1, lb, k, x(i))
											IF(nb .NE. nd) THEN
												total_ibb = total_ibb + IBrinkBookerHO(na - 1, la, nc - 1, la, nd - 1, lb, nb - 1, lb, k, x(i))
												total_ibb = total_ibb / 2.0D0
											END IF
											sum1 = sum1 + (coef_pair1(i) * tres_j_cuad * total_ibb)
											sum2 = sum2 + (coef_pair2(i) * cuad * tres_j_cuad * total_ibb)
										END DO
									END DO

									sum1 = sum1 * PAR(la + lb)
									sum2 = sum2 * PAR(la + lb)

									CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v1_pair, la, na, nc, lb, nb, nd, sum1)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v1_pair, lb, nb, nd, la, na, nc, sum1)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v2_pair, la, na, nc, lb, nb, nd, sum2)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v2_pair, lb, nb, nd, la, na, nc, sum2)

									IF (file_error .EQ. 0 .AND. Lmin .EQ. 0) THEN
										WRITE (file_desc, FMT="(6I3,2E24.16)", IOSTAT=file_error) &
											la, na, nc, lb, nb, nd, sum1, sum2
									END IF

								END DO
							END DO
						END DO
					END DO
				END DO
				PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
			END DO

		CASE(2)

			PRINT *, "Brink-Boker terms: Particle-Particle Channel - General Basis"

			DO la = Lmin, Lmax

				icount = 0
									namax = MIN(Nmax, NmaxOfL(la))
				IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

				DO na = 1, namax
					DO nc = 1, na

						DO lb = 0, la
												nbmax = MIN(Nmax, NmaxOfL(lb))
							IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

							DO nb = 1, nbmax
								DO nd = 1, nb

									kmin = ABS(la - lb)
									kmax = la + lb
									sum1 = 0.0D0
									sum2 = 0.0D0

									DO k = kmin, kmax, 2

										tres_j_cuad = (2 * k + 1) * (ThreeJSymbols_get(2 * la, 2 * k, 2 * lb) ** 2)
										cuad = CUAD2(la, lb, k)

										DO i = 0, i_gaus_max
									        	total_ibb = IBrinkBooker(na, la, nc, la, nb, lb, nd, lb, k, i, 0)
											IF(nb .NE. nd) THEN
												total_ibb = total_ibb + IBrinkBooker(na, la, nc, la, nd, lb, nb, lb, k, i, 0)
												total_ibb = total_ibb / 2.0D0
											END IF
											sum1 = sum1 + (coef_pair1(i) * tres_j_cuad * total_ibb)
											sum2 = sum2 + (coef_pair2(i) * cuad * tres_j_cuad * total_ibb)
										END DO

									END DO

									sum1 = sum1 * PAR(la + lb)
									sum2 = sum2 * PAR(la + lb)

									CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v1_pair, la, na, nc, lb, nb, nd, sum1)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v1_pair, lb, nb, nd, la, na, nc, sum1)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v2_pair, la, na, nc, lb, nb, nd, sum2)
									CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v2_pair, lb, nb, nd, la, na, nc, sum2)

									IF (file_error .EQ. 0 .AND. Lmin .EQ. 0) THEN
										WRITE (file_desc, FMT="(6I3,2E24.16)", IOSTAT=file_error) &
											la, na, nc, lb, nb, nd, sum1, sum2
									END IF
								END DO
							END DO
						END DO
					END DO
				END DO
				WRITE(*,'("Number of re-used integrals = ",I16)') icount
				WRITE(*,'("Completion rate: ",I3,"%")') INT(100 * (la + 1) / (N_0 + 1))
			END DO

		END SELECT

		! In case we want to optimize our calculation, we write the full matrix elements after
		! having calculated only the few ones required

		IF (Optimization .EQ. 1 .AND. Lmin .NE. 0) THEN
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) THEN
				WRITE(*,'("Incompatible options: No optimization for HO case (or compatible)")')
				STOP 'Error in SymVBBpp_calculate - No optimization pssible'
			END IF
			CALL SymVBBpp_write(vBBpp)
		END IF

		CLOSE(file_desc)

		RETURN
	END SUBROUTINE SymVBBpp_calculate

        !---------------------------------------------------------------------------------------!
	!											!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVBBpp_write(vBBpp)
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp

		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd
		DOUBLE PRECISION :: sum1, sum2

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		OPEN(file_desc, FILE=vBBpp%filename, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to write in ", vBBpp%filename
		ELSE
			WRITE (file_desc, FMT="(E24.16)", IOSTAT=file_error) vBBpp%b1
		END IF

		DO la = 0, Lmax

						namax = MIN(Nmax, NmaxOfL(la))
			IF (CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

			DO na = 1, namax
				DO nc = 1, na

					DO lb = 0, la
									nbmax = MIN(Nmax, NmaxOfL(lb))
						IF (CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

						DO nb = 1, nbmax
							DO nd = 1, nb

								CALL SymD3Tensor_SymD3Tensor_get(sum1, vBBpp%v1_pair, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum2, vBBpp%v2_pair, la, na, nc, lb, nb, nd)

								IF (file_error .EQ. 0) THEN
									WRITE (file_desc, FMT="(6I3,2E24.16)", IOSTAT=file_error) &
										la, na, nc, lb, nb, nd, sum1, sum2
								END IF
							END DO
						END DO
					END DO
				END DO
			END DO

		END DO

		CLOSE(file_desc)

		RETURN

	END SUBROUTINE SymVBBpp_write

	SUBROUTINE SymVBBpp_get_Delta(HF_out, vBBpp, P_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVBBpp), INTENT(IN) :: vBBpp
		TYPE (SymHartreeFockField), INTENT(IN) :: P_in

		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(0), vBBpp%v1_pair, P_in%p(0))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(1), vBBpp%v1_pair, P_in%p(1))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(0), vBBpp%v2_pair, P_in%a(0))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(1), vBBpp%v2_pair, P_in%a(1))

		RETURN
	END SUBROUTINE SymVBBpp_get_Delta

 	SUBROUTINE SymVBBpp_get_DeltaIso(HF_out, vBBpp, P_in, ta)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVBBpp), INTENT(IN) :: vBBpp
		TYPE (SymHartreeFockField), INTENT(IN) :: P_in
		INTEGER, INTENT(IN) :: ta

		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(ta), vBBpp%v1_pair, P_in%p(ta))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(ta), vBBpp%v2_pair, P_in%a(ta))
        	HF_out%p(1-ta)=0.0D0
        	HF_out%a(1-ta)=0.0D0

		RETURN
	END SUBROUTINE SymVBBpp_get_DeltaIso

        !---------------------------------------------------------------------------------------!
	!											!
	!   Subroutine reading the matrix elements of the Brink-Boker term in the particle-	!
	!   particle channel. This subroutine attempts reading a file containing previously 	!
	!   calculated matrix elements. It is used to decide whether or not to recalculate 	!
	!   everything										!
	!											!
        !---------------------------------------------------------------------------------------!

	FUNCTION SymVBBpp_read(vBBpp)
		LOGICAL SymVBBpp_read
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp

		DOUBLE PRECISION b1
		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER i1, i2, i3, i4, i5, i6
		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error
		DOUBLE PRECISION sum1, sum2

		OPEN(file_desc, FILE=vBBpp%filename, ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "No se pudo leer el archivo: ", vBBpp%filename
			SymVBBpp_read = .FALSE.
			RETURN
		END IF

		READ (file_desc, FMT=*, IOSTAT=file_error) b1
		IF ((file_error .NE. 0) .OR. (ABS(b1 - vBBpp%b1) .GE. 1.E-14)) THEN
			CLOSE(file_desc)
			SymVBBpp_read = .FALSE.
			RETURN
		END IF

		DO la = 0, Lmax
								namax = MIN(Nmax, NmaxOfL(la))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1
			DO na = 1, namax
				DO nc = 1, na

					DO lb = 0, la
											nbmax = MIN(Nmax, NmaxOfL(lb))
						IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

						DO nb = 1, nbmax
							DO nd = 1, nb

								READ (file_desc, FMT=*, IOSTAT=file_error) &
									i1, i2, i3, i4, i5, i6, sum1, sum2

								IF ((file_error .NE. 0) .OR. &
									(la .NE. i1) .OR. (na .NE. i2) .OR. (nc .NE. i3) .OR. &
									(lb .NE. i4) .OR. (nb .NE. i5) .OR. (nd .NE. i6)) THEN
									PRINT *, "Informacion no validad en el archivo: ", vBBpp%filename
									CLOSE(file_desc)
									SymVBBpp_read = .FALSE.
									RETURN
								END IF

								CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v1_pair, la, na, nc, lb, nb, nd, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v1_pair, lb, nb, nd, la, na, nc, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v2_pair, la, na, nc, lb, nb, nd, sum2)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v2_pair, lb, nb, nd, la, na, nc, sum2)

							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
		CLOSE(file_desc)

		SymVBBpp_read = .TRUE.
		RETURN
	END FUNCTION SymVBBpp_read

        !---------------------------------------------------------------------------------------!
	!											!
	!   Subroutine deleting the tensors corresponding to the matrix elements of the Brink-	!
	!   Boker in the particle-particle channel.						!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVBBpp_del(vBBpp)
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp

		CALL SymD3Tensor_SymD3Tensor_del(vBBpp%v1_pair)
		CALL SymD3Tensor_SymD3Tensor_del(vBBpp%v2_pair)

		RETURN
	END SUBROUTINE SymVBBpp_del

END MODULE symvbb
!---------------------------------------------------------------------!
!                                                                     !
!     MATRIX ELEMENTS OF THE COULOMB PART OF THE INTERACTION          !
!                                                                     !
!---------------------------------------------------------------------!

 MODULE symvc

	USE input
	USE angmom
	USE symd3t
	USE symfield
	USE ic

	IMPLICIT NONE

	TYPE SymVCph
		TYPE (SymD3Tensor_SymD3Tensor) :: v_local, v1_exch, v2_exch, v1
		CHARACTER(64) :: filename
	END TYPE

	TYPE SymVCpp
		TYPE (SymD3Tensor_SymD3Tensor) :: v1_pair, v2_pair
		CHARACTER(64) :: filename
	END TYPE

 CONTAINS

	!-----------------------------------------------------------------------!
	!  Initialization of the tensors and definition of the filenames (p.h.)	!
	!-----------------------------------------------------------------------!

	SUBROUTINE SymVCph_new(vCph)
		TYPE (SymVCph), INTENT(INOUT) :: vCph

		CALL SymD3Tensor_SymD3Tensor_new(vCph%v_local)
		CALL SymD3Tensor_SymD3Tensor_new(vCph%v1_exch)
		CALL SymD3Tensor_SymD3Tensor_new(vCph%v2_exch)
		CALL SymD3Tensor_SymD3Tensor_new(vCph%v1)

                SELECT CASE (Basis)

                CASE (1)
			IF (N_0 < 10) THEN
				WRITE(vCph%filename, "(A,I1,A)") "data/vC", N_0, "ph_HO.txt"
			ELSE
				WRITE(vCph%filename, "(A,I2,A)") "data/vC", N_0, "ph_HO.txt"
			END IF
                CASE(2)
                	IF (N_0 < 10) THEN
				WRITE(vCph%filename, "(A,I1,A)") "data/vC", N_0, "ph_WS.txt"
			ELSE
				WRITE(vCph%filename, "(A,I2,A)") "data/vC", N_0, "ph_WS.txt"
			END IF
                END SELECT

		RETURN
	END SUBROUTINE SymVCph_new

	!-----------------------------------------------------------------------!
	!     Calculation of the matrix elements of the Coulomb force (p.h.)	!
	!-----------------------------------------------------------------------!

	SUBROUTINE SymVCph_calculate(vCph, Read_Cph, Lmin)
		TYPE (SymVCph), INTENT(INOUT) :: vCph
		INTEGER, INTENT(INOUT) :: Lmin
		LOGICAL, INTENT(IN) :: Read_Cph

		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER :: k, kmin, kmax, Lold
		DOUBLE PRECISION :: total_iC, tres_j_cuad, cuad
		DOUBLE PRECISION :: sumi, sum1, sum2

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER :: file_error

		! If the matrix elements were already calculated, we read them and avoid thereby recalculating them
		! The tensors are then updated with the values read from tape.
		! In case of the WS basis, we (smartly) read what was already calculated and eventually calculate
		! only those elements that are missing.

		IF (Lmin .EQ. 0) THEN

			IF (Read_Cph) THEN
				CALL SymD3Tensor_SymD3Tensor_add(vCph%v1, vCph%v_local, vCph%v1_exch)
				RETURN
			ELSE
				IF (SymVCph_read(vCph)) THEN
					CALL SymD3Tensor_SymD3Tensor_add(vCph%v1, vCph%v_local, vCph%v1_exch)
					RETURN
				END IF
			END IF

			OPEN(file_desc, FILE=vCph%filename, ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				PRINT *, "*** Attention: Impossible to write in ", vCph%filename
			END IF

		ELSE

			Lold = Lmax
			Lmax = Lold - 1

			IF (Lmax < 10) THEN
				WRITE(vCph%filename, "(A,I1,A)") "data/vC", Lmax, "ph_WS.txt"
			ELSE
				WRITE(vCph%filename, "(A,I2,A)") "data/vC", Lmax, "ph_WS.txt"
			END IF

			IF (SymVCph_read(vCph)) write(*,'("Read file....")')

			Lmax = Lold

			IF (Lmax < 10) THEN
				WRITE(vCph%filename, "(A,I1,A)") "data/vC", Lmax, "ph_WS.txt"
			ELSE
				WRITE(vCph%filename, "(A,I2,A)") "data/vC", Lmax, "ph_WS.txt"
			END IF

		END IF

                SELECT CASE (Basis)

		! Calculation in the case of the spherical harmonic oscillator

                CASE (1)

			PRINT *, "Calculation of ph Coulomb matrix elements - Hamonic oscillator Basis"

			DO la = 0, Lmax
				namax = ((N_0 - la) / 2) + 1
				DO na = 1, namax
					DO nc = 1, na

						DO lb = 0, la
							nbmax = ((N_0 - lb) / 2) + 1
							DO nb = 1, nbmax
								DO nd = 1, nb

									total_iC = ICoulombHO(na - 1, la, nb - 1, lb, nc - 1, la, nd - 1, lb, 0)

									IF (nb .NE. nd) THEN
										total_iC = total_iC + ICoulombHO(na - 1, la, nd - 1, lb, nc - 1, la, nb - 1, lb, 0)
										total_iC = total_iC / 2.0D0
									END IF

									sumi = (VC * I_4PI) * total_iC

									CALL SymD3Tensor_SymD3Tensor_assign(vCph%v_local, la, na, nc, lb, nb, nd, sumi)
									CALL SymD3Tensor_SymD3Tensor_assign(vCph%v_local, lb, nb, nd, la, na, nc, sumi)

									kmin = ABS(la - lb)
									kmax = la + lb

									sum1 = 0.0D0
									sum2 = 0.0D0

									DO k = kmin, kmax, 2

										tres_j_cuad = (2 * k + 1) * (ThreeJSymbols_get(2 * la, 2 * k, 2 * lb) ** 2)

										cuad = CUAD2(la, lb, k)

										total_iC = ICoulombHO(na - 1, la, nb - 1, lb, nd - 1, lb, nc - 1, la, k)

										IF(nb .NE. nd) THEN
											total_iC = total_iC + ICoulombHO(na - 1, la, nd - 1, lb, nb - 1, lb, nc - 1, la, k)
											total_iC = total_iC / 2.0D0
										END IF

										total_iC = total_iC * (VC * I_4PI)

										sum1 = sum1 + (tres_j_cuad * total_iC)
										sum2 = sum2 + (cuad * tres_j_cuad * total_iC)
									END DO

									sum1 = -0.5D0 * sum1
									sum2 = -sum2

									CALL SymD3Tensor_SymD3Tensor_assign(vCph%v1_exch, la, na, nc, lb, nb, nd, sum1)
									CALL SymD3Tensor_SymD3Tensor_assign(vCph%v1_exch, lb, nb, nd, la, na, nc, sum1)
									CALL SymD3Tensor_SymD3Tensor_assign(vCph%v2_exch, la, na, nc, lb, nb, nd, sum2)
									CALL SymD3Tensor_SymD3Tensor_assign(vCph%v2_exch, lb, nb, nd, la, na, nc, sum2)

									IF (file_error .EQ. 0) THEN
										WRITE (file_desc, FMT="(6I3,3E24.16)", IOSTAT=file_error) &
											la, na, nc, lb, nb, nd, sumi, sum1, sum2
									END IF

								END DO
							END DO
						END DO
					END DO
				END DO
				PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
			END DO

		! Calculation in the case of a general spherical basis. The tensor algebra is the same as in the case of the
		! harmonic oscillator, only the radial integral is different

                CASE (2)

			PRINT *, "Calculation of ph Coulomb matrix elements - General Basis"

			DO la = Lmin, Lmax
									namax = MIN(Nmax, NmaxOfL(la))
				IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

				DO na = 1, namax
					DO nc = 1, na

						DO lb = 0, la
												nbmax = MIN(Nmax, NmaxOfL(lb))
							IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

							DO nb = 1, nbmax
								DO nd = 1, nb

									total_iC = ICoulomb(na, la, nb, lb, nc, la, nd, lb, 0)

									IF (nb .NE. nd) THEN
										total_iC = total_iC + ICoulomb(na, la, nd, lb, nc, la, nb, lb, 0)
										total_iC = total_iC / 2.0D0
									END IF

									sumi = (VC * I_4PI) * total_iC

									CALL SymD3Tensor_SymD3Tensor_assign(vCph%v_local, la, na, nc, lb, nb, nd, sumi)
									CALL SymD3Tensor_SymD3Tensor_assign(vCph%v_local, lb, nb, nd, la, na, nc, sumi)

									kmin = ABS(la - lb)
									kmax = la + lb

									sum1 = 0.0D0
									sum2 = 0.0D0

									DO k = kmin, kmax, 2

										tres_j_cuad = (2 * k + 1) * (ThreeJSymbols_get(2 * la, 2 * k, 2 * lb) ** 2)

										cuad = CUAD2(la, lb, k)

										total_iC = ICoulomb(na, la, nb, lb, nd, lb, nc, la, k)

										IF (nb .NE. nd) THEN
											total_iC = total_iC + ICoulomb(na, la, nd, lb, nb, lb, nc, la, k)
											total_iC = total_iC / 2.0D0
										END IF

										total_iC = total_iC * (VC * I_4PI)

										sum1 = sum1 + (tres_j_cuad * total_iC)
										sum2 = sum2 + (cuad * tres_j_cuad * total_iC)
									END DO

									sum1 = -0.5D0 * sum1
									sum2 = -sum2

									CALL SymD3Tensor_SymD3Tensor_assign(vCph%v1_exch, la, na, nc, lb, nb, nd, sum1)
									CALL SymD3Tensor_SymD3Tensor_assign(vCph%v1_exch, lb, nb, nd, la, na, nc, sum1)
									CALL SymD3Tensor_SymD3Tensor_assign(vCph%v2_exch, la, na, nc, lb, nb, nd, sum2)
									CALL SymD3Tensor_SymD3Tensor_assign(vCph%v2_exch, lb, nb, nd, la, na, nc, sum2)

									IF (file_error .EQ. 0 .AND. Lmin .EQ. 0) THEN
										WRITE (file_desc, FMT="(6I3,3E24.16)", IOSTAT=file_error) &
											la, na, nc, lb, nb, nd, sumi, sum1, sum2
									END IF

								END DO
							END DO
						END DO
					END DO
				END DO
				PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
			END DO

                END SELECT

		! In case we want to optimize our calculation, we write the full matrix elements after
		! having calculated only the few ones required

		IF (Optimization .EQ. 1 .AND. Lmin .NE. 0) THEN

			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) THEN
				WRITE(*,'("Incompatible options: No optimization for HO case (or compatible)")')
				STOP 'Error in SymVCph_calculate - No optimization pssible'
			END IF

			CALL SymVCph_write(vCph)

		END IF

		CALL SymD3Tensor_SymD3Tensor_add(vCph%v1, vCph%v_local, vCph%v1_exch)

		CLOSE(file_desc)

		RETURN
	END SUBROUTINE SymVCph_calculate

        !---------------------------------------------------------------------------------------!
	!											!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVCph_write(vCph)
		TYPE (SymVCph), INTENT(INOUT) :: vCph

		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd
		DOUBLE PRECISION :: sumi, sum1, sum2

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER :: file_error

		OPEN(file_desc, FILE=vCph%filename, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to write in ", vCph%filename
		END IF

		DO la = 0, Lmax

						namax = MIN(Nmax, NmaxOfL(la))
			IF (CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

			DO na = 1, namax
				DO nc = 1, na

					DO lb = 0, la
									nbmax = MIN(Nmax, NmaxOfL(lb))
						IF (CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

						DO nb = 1, nbmax
							DO nd = 1, nb

								CALL SymD3Tensor_SymD3Tensor_get(sumi, vCph%v_local, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum1, vCph%v1_exch, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum2, vCph%v2_exch, la, na, nc, lb, nb, nd)

								IF (file_error .EQ. 0) THEN
									WRITE (file_desc, FMT="(6I3,3E24.16)", IOSTAT=file_error) &
									       la, na, nc, lb, nb, nd, sumi, sum1, sum2
								END IF
							END DO
						END DO
					END DO
				END DO
			END DO

		END DO

		CLOSE(file_desc)

		RETURN

	END SUBROUTINE SymVCph_write

	SUBROUTINE SymVCph_get_Gamma(HF_out, vCph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVCph), INTENT(IN) :: vCph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(0), vCph%v1, HF_in%p(0))
		HF_out%p(1) = 0.0D0
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(0), vCph%v2_exch, HF_in%a(0))
		HF_out%a(1) = 0.0D0

		RETURN
	END SUBROUTINE SymVCph_get_Gamma

	SUBROUTINE SymVCph_get_LocalGamma(HF_out, vCph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVCph), INTENT(IN) :: vCph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(0), vCph%v_local, HF_in%p(0))
		HF_out%p(1) = 0.0D0
		HF_out%a(0) = 0.0D0
		HF_out%a(1) = 0.0D0

		RETURN
	END SUBROUTINE SymVCph_get_LocalGamma

	SUBROUTINE SymVCph_get_ExchangeGamma(HF_out, vCph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVCph), INTENT(IN) :: vCph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(0), vCph%v1_exch, HF_in%p(0))
		HF_out%p(1) = 0.0D0
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(0), vCph%v2_exch, HF_in%a(0))
		HF_out%a(1) = 0.0D0

		RETURN
	END SUBROUTINE SymVCph_get_ExchangeGamma

	!-----------------------------------------------------------------------!
	!   Reading the matrix elements in the p.h. channel from a file. If the !
	!   file exists and its format is valid, we fill in the appropriate 	!
	!   tensors with the values read.					!
	!-----------------------------------------------------------------------!

	FUNCTION SymVCph_read(vCph)
		LOGICAL SymVCph_read
		TYPE (SymVCph), INTENT(INOUT) :: vCph

		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER :: i1, i2, i3, i4, i5, i6
		INTEGER, PARAMETER :: file_desc = 16
		INTEGER :: file_error
		DOUBLE PRECISION :: sumi, sum1, sum2

		OPEN(file_desc, FILE=vCph%filename, ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "No se pudo leer el archivo: ", vCph%filename
			SymVCph_read = .FALSE.
			RETURN
		END IF

		DO la = 0, Lmax
								namax = MIN(Nmax, NmaxOfL(la))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

			DO na = 1, namax
				DO nc = 1, na
					DO lb = 0, la
											nbmax = MIN(Nmax, NmaxOfL(lb))
						IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

						DO nb = 1, nbmax
							DO nd = 1, nb
								READ (file_desc, FMT=*, IOSTAT=file_error) &
									i1, i2, i3, i4, i5, i6, sumi, sum1, sum2
								IF ((file_error .NE. 0) .OR. &
									(la .NE. i1) .OR. (na .NE. i2) .OR. (nc .NE. i3) .OR. &
									(lb .NE. i4) .OR. (nb .NE. i5) .OR. (nd .NE. i6)) THEN
									PRINT *, "Informacion no validad en el archivo: ", vCph%filename
									CLOSE(file_desc)
									SymVCph_read = .FALSE.
									RETURN
								END IF
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v_local, la, na, nc, lb, nb, nd, sumi)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v_local, lb, nb, nd, la, na, nc, sumi)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v1_exch, la, na, nc, lb, nb, nd, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v1_exch, lb, nb, nd, la, na, nc, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v2_exch, la, na, nc, lb, nb, nd, sum2)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v2_exch, lb, nb, nd, la, na, nc, sum2)
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
		CLOSE(file_desc)
		SymVCph_read = .TRUE.
		RETURN
	END FUNCTION SymVCph_read

	SUBROUTINE SymVCph_del(vCph)
		TYPE (SymVCph), INTENT(INOUT) :: vCph

		CALL SymD3Tensor_SymD3Tensor_del(vCph%v_local)
		CALL SymD3Tensor_SymD3Tensor_del(vCph%v1_exch)
		CALL SymD3Tensor_SymD3Tensor_del(vCph%v2_exch)
		CALL SymD3Tensor_SymD3Tensor_del(vCph%v1)
		RETURN
	END SUBROUTINE SymVCph_del

	!-----------------------------------------------------------------------!
	!  Initialization of the tensors and definition of the filenames (p.p.)	!
	!-----------------------------------------------------------------------!

	SUBROUTINE SymVCpp_new(vCpp)
		TYPE (SymVCpp), INTENT(INOUT) :: vCpp

		CALL SymD3Tensor_SymD3Tensor_new(vCpp%v1_pair)
		CALL SymD3Tensor_SymD3Tensor_new(vCpp%v2_pair)

                SELECT CASE (Basis)

                CASE(1)
			IF (N_0 < 10) THEN
				WRITE(vCpp%filename, "(A,I1,A)") "data/vC", N_0, "pp_HO.txt"
			ELSE
				WRITE(vCpp%filename, "(A,I2,A)") "data/vC", N_0, "pp_HO.txt"
			END IF
                CASE(2)
			IF (N_0 < 10) THEN
				WRITE(vCpp%filename, "(A,I1,A)") "data/vC", N_0, "pp_WS.txt"
			ELSE
				WRITE(vCpp%filename, "(A,I2,A)") "data/vC", N_0, "pp_WS.txt"
			END IF

                END SELECT

		RETURN
	END SUBROUTINE SymVCpp_new

	!-----------------------------------------------------------------------!
	!     Calculation of the matrix elements of the Coulomb force (p.p.)	!
	!-----------------------------------------------------------------------!

	SUBROUTINE SymVCpp_calculate(vCpp, Read_Cpp, Lmin)
		TYPE (SymVCpp), INTENT(INOUT) :: vCpp
		INTEGER, INTENT(INOUT) :: Lmin
		LOGICAL, INTENT(IN) :: Read_Cpp

		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER :: k, kmin, kmax, Lold
		DOUBLE PRECISION :: tres_j_cuad, cuad, total_iC
		DOUBLE PRECISION :: sum1, sum2

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		! If the matrix elements were already calculated, we read them and avoid thereby recalculating them
		! The tensors are then updated with the values read from tape.
		! In case of the WS basis, we (smartly) read what was already calculated and eventually calculate
		! only those elements that are missing.

		IF (Lmin .EQ. 0) THEN

			IF (Read_Cpp) RETURN

			IF (SymVCpp_read(vCpp)) RETURN

			OPEN(file_desc, FILE=vCpp%filename, ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				PRINT *, "*** Attention: Impossible to write in ", vCpp%filename
			END IF

		ELSE

			Lold = Lmax
			Lmax = Lold - 1

			IF (Lmax < 10) THEN
				WRITE(vCpp%filename, "(A,I1,A)") "data/vC", Lmax, "pp_WS.txt"
			ELSE
				WRITE(vCpp%filename, "(A,I2,A)") "data/vC", Lmax, "pp_WS.txt"
			END IF

			IF (SymVCpp_read(vCpp)) write(*,'("Read file....")')

			Lmax = Lold

			IF (Lmax < 10) THEN
				WRITE(vCpp%filename, "(A,I1,A)") "data/vC", Lmax, "pp_WS.txt"
			ELSE
				WRITE(vCpp%filename, "(A,I2,A)") "data/vC", Lmax, "pp_WS.txt"
			END IF

		END IF

                SELECT CASE (Basis)

		! Calculation in the case of the spherical harmonic oscillator

                CASE (1)

			PRINT *, "Calculation of pp Coulomb matrix elements - Harmonic Oscillator Basis"

			DO la = 0, Lmax
				namax = ((N_0 - la) / 2) + 1
				DO na = 1, namax
					DO nc = 1, na

						DO lb = 0, la
							nbmax = ((N_0 - lb) / 2) + 1
							DO nb = 1, nbmax
								DO nd = 1, nb

									kmin = ABS(la - lb)
									kmax = la + lb
									sum1 = 0.0D0
									sum2 = 0.0D0

									DO k = kmin, kmax, 2

										tres_j_cuad = (2 * k + 1) * (ThreeJSymbols_get(2 * la, 2 * k, 2 * lb) ** 2)
										cuad = CUAD2(la, lb, k)
										total_iC = (VC * I_4PI) * ICoulombHO(na - 1, la, nc - 1, la, nb - 1, lb, nd - 1, lb, k)

										IF (nb .NE. nd) THEN
											total_iC = total_iC + (VC * I_4PI) * &
												ICoulombHO(na - 1, la, nc - 1, la, nd - 1, lb, nb - 1, lb, 0)
											total_iC = total_iC / 2.0D0
										END IF

										sum1 = sum1 + (       tres_j_cuad * total_iC)
										sum2 = sum2 + (cuad * tres_j_cuad * total_iC)

									END DO

									sum1 = 0.5D0 * PAR(la + lb) * sum1
									sum2 =         PAR(la + lb) * sum2

									CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v1_pair, la, na, nc, lb, nb, nd, sum1)
									CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v1_pair, lb, nb, nd, la, na, nc, sum1)
									CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v2_pair, la, na, nc, lb, nb, nd, sum2)
									CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v2_pair, lb, nb, nd, la, na, nc, sum2)

									IF (file_error .EQ. 0) THEN
										WRITE (file_desc, "(6I3,2E24.16)", IOSTAT=file_error) &
										la, na, nc, lb, nb, nd, sum1, sum2
									END IF
								END DO
							END DO
						END DO
					END DO
				END DO
				PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
			END DO

		! Calculation in the case of a general spherical basis. The tensor algebra is the same as in the case of the
		! harmonic oscillator, only the radial integral is different

                CASE (2)

			PRINT *, "Calculation of pp Coulomb matrix elements - General Basis"

			DO la = Lmin, Lmax
									namax = MIN(Nmax, NmaxOfL(la))
				IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

				DO na = 1, namax
					DO nc = 1, na

						DO lb = 0, la
												nbmax = MIN(Nmax, NmaxOfL(lb))
							IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

							DO nb = 1, nbmax
								DO nd = 1, nb

									kmin = ABS(la - lb)
									kmax = la + lb
									sum1 = 0.0D0
									sum2 = 0.0D0

									DO k = kmin, kmax, 2

										tres_j_cuad = (2 * k + 1) * (ThreeJSymbols_get(2 * la, 2 * k, 2 * lb) ** 2)
										cuad = CUAD2(la, lb, k)
										total_iC = (VC * I_4PI) * ICoulomb(na, la, nc, la, nb, lb, nd, lb, k)

										IF (nb .NE. nd) THEN
											total_iC = total_iC + (VC * I_4PI) * ICoulomb(na, la, nc, la, nd, lb, nb, lb, 0)
											total_iC = total_iC / 2.0D0
										END IF

										sum1 = sum1 + (       tres_j_cuad * total_iC)
										sum2 = sum2 + (cuad * tres_j_cuad * total_iC)

									END DO

									sum1 = 0.5D0 * PAR(la + lb) * sum1
									sum2 =         PAR(la + lb) * sum2

									CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v1_pair, la, na, nc, lb, nb, nd, sum1)
									CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v1_pair, lb, nb, nd, la, na, nc, sum1)
									CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v2_pair, la, na, nc, lb, nb, nd, sum2)
									CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v2_pair, lb, nb, nd, la, na, nc, sum2)

									IF (file_error .EQ. 0 .AND. Lmin .EQ. 0) THEN
										WRITE (file_desc, "(6I3,2E24.16)", IOSTAT=file_error) &
										la, na, nc, lb, nb, nd, sum1, sum2
									END IF
								END DO
							END DO
						END DO
					END DO
				END DO
				PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
			END DO

                END SELECT

		! In case we want to optimize our calculation, we write the full matrix elements after
		! having calculated only the few ones required

		IF (Optimization .EQ. 1 .AND. Lmin .NE. 0) THEN

			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) THEN
				WRITE(*,'("Incompatible options: No optimization for HO case (or compatible)")')
				STOP 'Error in SymVCpp_calculate - No optimization pssible'
			END IF

			CALL SymVCpp_write(vCpp)

		END IF

		CLOSE(file_desc)

		RETURN
	END SUBROUTINE SymVCpp_calculate

        !---------------------------------------------------------------------------------------!
	!											!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVCpp_write(vCpp)
		TYPE (SymVCpp), INTENT(INOUT) :: vCpp

		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd
		DOUBLE PRECISION :: sum1, sum2

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER :: file_error

		OPEN(file_desc, FILE=vCpp%filename, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to write in ", vCpp%filename
		END IF

		DO la = 0, Lmax

						namax = MIN(Nmax, NmaxOfL(la))
			IF (CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

			DO na = 1, namax
				DO nc = 1, na

					DO lb = 0, la
									nbmax = MIN(Nmax, NmaxOfL(lb))
						IF (CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

						DO nb = 1, nbmax
							DO nd = 1, nb

								CALL SymD3Tensor_SymD3Tensor_get(sum1, vCpp%v1_pair, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum2, vCpp%v2_pair, la, na, nc, lb, nb, nd)

								IF (file_error .EQ. 0) THEN
									WRITE (file_desc, FMT="(6I3,2E24.16)", IOSTAT=file_error) &
									       la, na, nc, lb, nb, nd, sum1, sum2
								END IF
							END DO
						END DO
					END DO
				END DO
			END DO

		END DO

		CLOSE(file_desc)

		RETURN

	END SUBROUTINE SymVCpp_write

	SUBROUTINE SymVCpp_get_Delta(HF_out, vCpp, P_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVCpp), INTENT(IN) :: vCpp
		TYPE (SymHartreeFockField), INTENT(IN) :: P_in

		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(0), vCpp%v1_pair, P_in%p(0))
		HF_out%p(1) = 0.0D0
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(0), vCpp%v2_pair, P_in%a(0))
		HF_out%a(1) = 0.0D0

		RETURN
	END SUBROUTINE SymVCpp_get_Delta

	!-----------------------------------------------------------------------!
	!   Reading the matrix elements in the p.p. channel from a file. If the !
	!   file exists and its format is valid, we fill in the appropriate 	!
	!   tensors with the values read.					!
	!-----------------------------------------------------------------------!

	FUNCTION SymVCpp_read(vCpp)
		LOGICAL SymVCpp_read
		TYPE (SymVCpp), INTENT(INOUT) :: vCpp

		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER :: i1, i2, i3, i4, i5, i6
		INTEGER, PARAMETER :: file_desc = 16
		INTEGER :: file_error
		DOUBLE PRECISION :: sum1, sum2

		OPEN(file_desc, FILE=vCpp%filename, ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "No se pudo leer el archivo: ", vCpp%filename
			SymVCpp_read = .FALSE.
			RETURN
		END IF

		DO la = 0, Lmax
								namax = MIN(Nmax, NmaxOfL(la))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

			DO na = 1, namax
				DO nc = 1, na
					DO lb = 0, la
											nbmax = MIN(Nmax, NmaxOfL(lb))
						IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

						DO nb = 1, nbmax
							DO nd = 1, nb
								READ (file_desc, FMT=*, IOSTAT=file_error) &
									i1, i2, i3, i4, i5, i6, sum1, sum2
								IF ((file_error .NE. 0) .OR. &
									(la .NE. i1) .OR. (na .NE. i2) .OR. (nc .NE. i3) .OR. &
									(lb .NE. i4) .OR. (nb .NE. i5) .OR. (nd .NE. i6)) THEN
									PRINT *, "Informacion no validad en el archivo: ", vCpp%filename
									CLOSE(file_desc)
									SymVCpp_read = .FALSE.
									RETURN
								END IF
								CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v1_pair, la, na, nc, lb, nb, nd, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v1_pair, lb, nb, nd, la, na, nc, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v2_pair, la, na, nc, lb, nb, nd, sum2)
								CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v2_pair, lb, nb, nd, la, na, nc, sum2)
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
		CLOSE(file_desc)

		SymVCpp_read = .TRUE.
		RETURN
	END FUNCTION SymVCpp_read

	SUBROUTINE SymVCpp_del(vCpp)
		TYPE (SymVCpp), INTENT(INOUT) :: vCpp

		CALL SymD3Tensor_SymD3Tensor_del(vCpp%v1_pair)
		CALL SymD3Tensor_SymD3Tensor_del(vCpp%v2_pair)
		RETURN
	END SUBROUTINE SymVCpp_del

END MODULE symvc
!----------------------------------------------------------------!
!								 !
!    CALCULATION OF THE SPIN-ORBIT TERM OF THE GOGNY FORCE       !
!								 !
!----------------------------------------------------------------!

 MODULE symvls

	USE input
	USE global
	USE symd3t
	USE symfield
	USE ils

	IMPLICIT NONE

	TYPE SymVLSph
		TYPE (SymD3Tensor_SymD3Tensor) v12, v21
		!CHARACTER(64) filename
	END TYPE

	TYPE SymVLSpp
		TYPE (SymD3Tensor_SymD3Tensor) v22
		!CHARACTER(64) filename
	END TYPE

	CHARACTER(64) :: FilePH,  FilePP

 CONTAINS

	!----------------------------------------------------------------!
	!								 !
	!    		PARTICLE-HOLE CHANNEL (MEAN-FIELD)  	         !
	!								 !
	!----------------------------------------------------------------!

	SUBROUTINE SymVLSph_new(vLSph)
		TYPE (SymVLSph), INTENT(INOUT) :: vLSph

		CALL SymD3Tensor_SymD3Tensor_new(vLSph%v12)
		CALL SymD3Tensor_SymD3Tensor_new(vLSph%v21)

                SELECT CASE (Basis)

		CASE (1)

			IF (N_0 < 10) THEN
				WRITE(FilePH, "(A,I1,A)") "data/vLS", N_0, "ph_HO.txt"
			ELSE
				WRITE(FilePH, "(A,I2,A)") "data/vLS", N_0, "ph_HO.txt"
			END IF

		CASE (2)

			IF (N_0 < 10) THEN
				WRITE(FilePH, "(A,I1,A)") "data/vLS", N_0, "ph_WS.txt"
			ELSE
				WRITE(FilePH, "(A,I2,A)") "data/vLS", N_0, "ph_WS.txt"
			END IF

		END SELECT

		RETURN
	END SUBROUTINE SymVLSph_new

	SUBROUTINE SymVLSph_calculate(vLSph)
		TYPE (SymVLSph), INTENT(INOUT) :: vLSph

		INTEGER la, na, namax, nc, lb, nb, nbmax, nd
		DOUBLE PRECISION sum1, sum2

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		IF (SymVLSph_read(vLSph)) RETURN

		OPEN(file_desc, FILE=FilePH, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** AVISO: No se pueden escribir los resultados en ", FilePH
		END IF

                SELECT CASE (Basis)

		CASE (1)

			PRINT *, "Calculation of the matrix elements of the spin-orbit term:"
			DO la = 0, Lmax
				namax = ((N_0 - la) / 2) + 1
				DO na = 1, namax
					DO nc = 1, na

						DO lb = 0, la
							nbmax = ((N_0 - lb) / 2) + 1
							DO nb = 1, nbmax
								DO nd = 1, nb

									IF (lb .EQ. 0) THEN
										sum1 = 0.0D0
									ELSE
										sum1 = DBLE(lb * (lb + 1)) * IHFLSHO(nb - 1, nd - 1, lb, na - 1, nc - 1, la)
									END IF
									IF (la .EQ. 0) THEN
										sum2 = 0.0D0
									ELSE
										sum2 = DBLE(la * (la + 1)) * IHFLSHO(na - 1, nc - 1, la, nb - 1, nd - 1, lb)
									END IF

									CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v12, la, na, nc, lb, nb, nd, sum1)
									CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v21, la, na, nc, lb, nb, nd, sum2)

									IF (la .NE. lb) THEN
										! Ojo, estan invertidos los terminos
										CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v12, lb, nb, nd, la, na, nc, sum2)
										CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v21, lb, nb, nd, la, na, nc, sum1)
									END IF

									IF (file_error .EQ. 0) THEN
										WRITE (file_desc, FMT="(6I3,2E24.16)", IOSTAT=file_error) &
											la, na, nc, lb, nb, nd, sum1, sum2
									END IF

								END DO
							END DO
						END DO

					END DO
				END DO
				PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculated"
			END DO
			CLOSE(file_desc)

		CASE (2)

			PRINT *, "Calculation of the matrix elements of the spin-orbit term:"
			DO la = 0, Lmax
									namax = MIN(Nmax, NmaxOfL(la))
				IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

				DO na = 1, namax
					DO nc = 1, na

						DO lb = 0, la
												nbmax = MIN(Nmax, NmaxOfL(lb))
							IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

							DO nb = 1, nbmax
								DO nd = 1, nb

									IF (lb .EQ. 0) THEN
										sum1 = 0.0D0
									ELSE
										sum1 = DBLE(lb * (lb + 1)) * IHFLS(nb, nd, lb, na, nc, la)
									END IF
									IF (la .EQ. 0) THEN
										sum2 = 0.0D0
									ELSE
										sum2 = DBLE(la * (la + 1)) * IHFLS(na, nc, la, nb, nd, lb)
									END IF

									CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v12, la, na, nc, lb, nb, nd, sum1)
									CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v21, la, na, nc, lb, nb, nd, sum2)

									IF (la .NE. lb) THEN
										! Ojo, estan invertidos los terminos
										CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v12, lb, nb, nd, la, na, nc, sum2)
										CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v21, lb, nb, nd, la, na, nc, sum1)
									END IF

									IF (file_error .EQ. 0) THEN
										WRITE (file_desc, FMT="(6I3,2E24.16)", IOSTAT=file_error) &
											la, na, nc, lb, nb, nd, sum1, sum2
									END IF

								END DO
							END DO
						END DO

					END DO
				END DO
				PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculated"
			END DO
			CLOSE(file_desc)

		END SELECT

		RETURN
	END SUBROUTINE SymVLSph_calculate

	! Nota: si y0 != 1 este factor cambia la fuerza
	! Se pretende que la fuerza de spin-orbita sea dependiente del isosoin
	SUBROUTINE SymVLSph_get_Gamma(HF_out, vLSph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVLSph), INTENT(IN) :: vLSph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		DOUBLE PRECISION factor
		TYPE (SymD3Tensor) tmp1, tmp2, tmp3

		CALL SymD3Tensor_new(tmp1)
		CALL SymD3Tensor_new(tmp2)
		CALL SymD3Tensor_new(tmp3)

		! Atención: si (y0 .NE. 1), este factor cambia la fuerza.
		! Se pretende que la fuerza de spin-orbita sea dependiente del isosoin.
		factor = Gogny_W0(Gogny) * I_4PI

		CALL SymD3Tensor_product(tmp1, 1.0D0 + x0, HF_in%a(0))
		CALL SymD3Tensor_product(tmp2, x0, HF_in%a(1))
		CALL SymD3Tensor_add(tmp3, tmp1, tmp2)

		CALL SymD3Tensor_SymD3Tensor_product(tmp1, vLSph%v12, tmp3)
		CALL SymD3Tensor_product(HF_out%p(0), factor, tmp1)

		CALL SymD3Tensor_product(tmp1, 1.0D0 + x0, HF_in%a(1))
		CALL SymD3Tensor_product(tmp2, x0, HF_in%a(0))
		CALL SymD3Tensor_add(tmp3, tmp1, tmp2)

		CALL SymD3Tensor_SymD3Tensor_product(tmp1, vLSph%v12, tmp3)
		CALL SymD3Tensor_product(HF_out%p(1), factor, tmp1)

		CALL SymD3Tensor_product(tmp1, 1.0D0 + x0, HF_in%p(0))
		CALL SymD3Tensor_product(tmp2, x0, HF_in%p(1))
		CALL SymD3Tensor_add(tmp3, tmp1, tmp2)

		CALL SymD3Tensor_SymD3Tensor_product(tmp1, vLSph%v21, tmp3)
		CALL SymD3Tensor_product(HF_out%a(0), factor, tmp1)

		CALL SymD3Tensor_product(tmp1, 1.0D0 + x0, HF_in%p(1))
		CALL SymD3Tensor_product(tmp2, x0, HF_in%p(0))
		CALL SymD3Tensor_add(tmp3, tmp1, tmp2)

		CALL SymD3Tensor_SymD3Tensor_product(tmp1, vLSph%v21, tmp3)
		CALL SymD3Tensor_product(HF_out%a(1), factor, tmp1)

		CALL SymD3Tensor_del(tmp1)
		CALL SymD3Tensor_del(tmp2)
		CALL SymD3Tensor_del(tmp3)

		RETURN
	END SUBROUTINE SymVLSph_get_Gamma

	FUNCTION SymVLSph_read(vLSph)
		LOGICAL SymVLSph_read
		TYPE (SymVLSph), INTENT(INOUT) :: vLSph

		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER i1, i2, i3, i4, i5, i6
		DOUBLE PRECISION sum1, sum2

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		OPEN(file_desc, FILE=FilePH, ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "No se pudo leer el archivo: ", FilePH
			SymVLSph_read = .FALSE.
			RETURN
		END IF

		DO la = 0, Lmax
								namax = MIN(Nmax, NmaxOfL(la))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

			DO na = 1, namax
				DO nc = 1, na

					DO lb = 0, la
											nbmax = MIN(Nmax, NmaxOfL(lb))
						IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

						DO nb = 1, nbmax
							DO nd = 1, nb

								READ (file_desc, FMT=*, IOSTAT=file_error) &
									i1, i2, i3, i4, i5, i6, sum1, sum2
								IF ((file_error .NE. 0) .OR. &
									(la .NE. i1) .OR. (na .NE. i2) .OR. (nc .NE. i3) .OR. &
									(lb .NE. i4) .OR. (nb .NE. i5) .OR. (nd .NE. i6)) THEN
									PRINT *, "Informacion no validad en el archivo: ", FilePH
									CLOSE(file_desc)
									SymVLSph_read = .FALSE.
									RETURN
								END IF

								CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v12, la, na, nc, lb, nb, nd, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v21, la, na, nc, lb, nb, nd, sum2)

								IF (la .NE. lb) THEN
									! Ojo, están invertidos los términos
									CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v12, lb, nb, nd, la, na, nc, sum2)
									CALL SymD3Tensor_SymD3Tensor_assign(vLSph%v21, lb, nb, nd, la, na, nc, sum1)
								END IF

							END DO
						END DO
					END DO

				END DO
			END DO
		END DO
		CLOSE(file_desc)

		SymVLSph_read = .TRUE.
		RETURN
	END FUNCTION SymVLSph_read

	SUBROUTINE SymVLSph_del(vLSph)
		TYPE (SymVLSph), INTENT(INOUT) :: vLSph

		CALL SymD3Tensor_SymD3Tensor_del(vLSph%v12)
		CALL SymD3Tensor_SymD3Tensor_del(vLSph%v21)

		RETURN
	END SUBROUTINE SymVLSph_del

	!----------------------------------------------------------------!
	!								 !
	!    		PARTICLE-PARTICLE CHANNEL (PAIRING)  	         !
	!								 !
	!----------------------------------------------------------------!

	SUBROUTINE SymVLSpp_new(vLSpp)
		TYPE (SymVLSpp), INTENT(INOUT) :: vLSpp

		CALL SymD3Tensor_SymD3Tensor_new(vLSpp%v22)

                SELECT CASE (Basis)

		CASE (1)

			IF (N_0 < 10) THEN
				WRITE(FilePP, "(A,I1,A)") "data/vLS", N_0, "pp_HO.txt"
			ELSE
				WRITE(FilePP, "(A,I2,A)") "data/vLS", N_0, "pp_HO.txt"
			END IF
		CASE (2)

			IF (N_0 < 10) THEN
				WRITE(FilePP, "(A,I1,A)") "data/vLS", N_0, "pp_WS.txt"
			ELSE
				WRITE(FilePP, "(A,I2,A)") "data/vLS", N_0, "pp_WS.txt"
			END IF

		END SELECT

		RETURN
	END SUBROUTINE SymVLSpp_new

	SUBROUTINE SymVLSpp_calculate(vLSpp)
		TYPE (SymVLSpp), INTENT(INOUT) :: vLSpp

		INTEGER la, na, namax, nc, lb, nb, nbmax, nd
		DOUBLE PRECISION sumi

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		IF (SymVLSpp_read(vLSpp)) RETURN

		OPEN(file_desc, FILE=FilePP, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** AVISO: No se pueden escribir los resultados en ", FilePP
		END IF

                SELECT CASE (Basis)

		CASE (1)

			PRINT *, "Se van a calcular los elementos de matriz de apar de LS"
			DO la = 0, Lmax
				namax = ((N_0 - la) / 2) + 1
				DO na = 1, namax
					DO nc = 1, na

						DO lb = 0, la
							nbmax = ((N_0 - lb) / 2) + 1
							DO nb = 1, nbmax
								DO nd = 1, nb

									IF ((la .EQ. 0) .OR. (lb .EQ. 0)) THEN
										sumi = 0.0D0
									ELSE
										sumi = DBLE(PAR(la + lb) * 2 * la * (la + 1) * lb * (lb + 1)) * &
											IPLSHO(na - 1, nc - 1, la, nb - 1, nd - 1, lb)
									END IF

									CALL SymD3Tensor_SymD3Tensor_assign(vLSpp%v22, la, na, nc, lb, nb, nd, sumi)

									IF (la .NE. lb) THEN
										CALL SymD3Tensor_SymD3Tensor_assign(vLSpp%v22, lb, nb, nd, la, na, nc, sumi)
									END IF

									IF (file_error .EQ. 0) THEN
										WRITE (file_desc, FMT="(6I3,E24.16)", IOSTAT=file_error) &
											la, na, nc, lb, nb, nd, sumi
									END IF

								END DO
							END DO
						END DO

					END DO
				END DO
				PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
			END DO

		CASE(2)

			PRINT *, "Se van a calcular los elementos de matriz de apar de LS"
			DO la = 0, Lmax
									namax = MIN(Nmax, NmaxOfL(la))
				IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

				DO na = 1, namax
					DO nc = 1, na

						DO lb = 0, la
												nbmax = MIN(Nmax, NmaxOfL(lb))
							IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

							DO nb = 1, nbmax
								DO nd = 1, nb

									IF ((la .EQ. 0) .OR. (lb .EQ. 0)) THEN
										sumi = 0.0D0
									ELSE
										sumi = DBLE(PAR(la + lb) * 2 * la * (la + 1) * lb * (lb + 1)) * &
											IPLS(na, nc, la, nb, nd, lb)
									END IF

									CALL SymD3Tensor_SymD3Tensor_assign(vLSpp%v22, la, na, nc, lb, nb, nd, sumi)

									IF (la .NE. lb) THEN
										CALL SymD3Tensor_SymD3Tensor_assign(vLSpp%v22, lb, nb, nd, la, na, nc, sumi)
									END IF

									IF (file_error .EQ. 0) THEN
										WRITE (file_desc, FMT="(6I3,E24.16)", IOSTAT=file_error) &
											la, na, nc, lb, nb, nd, sumi
									END IF

								END DO
							END DO
						END DO

					END DO
				END DO
				PRINT "(I3,A)", INT(100 * (la + 1) / (N_0 + 1)), "% calculado"
			END DO

		END SELECT

		 CLOSE(file_desc)
		RETURN
	END SUBROUTINE SymVLSpp_calculate

	SUBROUTINE SymVLSpp_get_Delta(HF_out, vLSpp, P_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVLSpp), INTENT(IN) :: vLSpp
		TYPE (SymHartreeFockField), INTENT(IN) :: P_in

		TYPE (SymD3Tensor) :: tmp
		DOUBLE PRECISION :: factor

		CALL SymD3Tensor_new(tmp)

		HF_out%p(0) = 0.0D0
		HF_out%p(1) = 0.0D0

		factor = Gogny_W0(Gogny) * I_4PI

		CALL SymD3Tensor_SymD3Tensor_product(tmp, vLSpp%v22, P_in%a(0))
		CALL SymD3Tensor_product(HF_out%a(0), factor, tmp)

		CALL SymD3Tensor_SymD3Tensor_product(tmp, vLSpp%v22, P_in%a(1))
		CALL SymD3Tensor_product(HF_out%a(1), factor, tmp)

		CALL SymD3Tensor_new(tmp)

		RETURN
	END SUBROUTINE SymVLSpp_get_Delta

	SUBROUTINE SymVLSpp_get_DeltaIso(HF_out, vLSpp, P_in, ta)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymVLSpp), INTENT(IN) :: vLSpp
		TYPE (SymHartreeFockField), INTENT(IN) :: P_in
		INTEGER, INTENT(IN) :: ta

		DOUBLE PRECISION :: factor
		TYPE (SymD3Tensor) :: tmp

		CALL SymD3Tensor_new(tmp)

		HF_out%p(ta) = 0.0D0

		factor = Gogny_W0(Gogny) * I_4PI

		CALL SymD3Tensor_SymD3Tensor_product(tmp, vLSpp%v22, P_in%a(ta))
		CALL SymD3Tensor_product(HF_out%a(ta), factor, tmp)
        	HF_out%p(1-ta)=0.0D0
        	HF_out%a(1-ta)=0.0D0

		CALL SymD3Tensor_new(tmp)

		RETURN
	END SUBROUTINE SymVLSpp_get_DeltaIso

	FUNCTION SymVLSpp_read(vLSpp)
		LOGICAL SymVLSpp_read
		TYPE (SymVLSpp), INTENT(INOUT) :: vLSpp

		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER i1, i2, i3, i4, i5, i6
		DOUBLE PRECISION sumi

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER file_error

		OPEN(file_desc, FILE=FilePP, ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "No se pudo leer el archivo: ", FilePP
			SymVLSpp_read = .FALSE.
			RETURN
		END IF

		DO la = 0, Lmax
								namax = MIN(Nmax, NmaxOfL(la))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

			DO na = 1, namax
				DO nc = 1, na

					DO lb = 0, la
											nbmax = MIN(Nmax, NmaxOfL(lb))
						IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

						DO nb = 1, nbmax
							DO nd = 1, nb

								READ (file_desc, FMT=*, IOSTAT=file_error) &
									i1, i2, i3, i4, i5, i6, sumi
								IF ((file_error .NE. 0) .OR. &
									(la .NE. i1) .OR. (na .NE. i2) .OR. (nc .NE. i3) .OR. &
									(lb .NE. i4) .OR. (nb .NE. i5) .OR. (nd .NE. i6)) THEN
									PRINT *, "Informacion no validad en el archivo: ", FilePP
									CLOSE(file_desc)
									SymVLSpp_read = .FALSE.
									RETURN
								END IF
								CALL SymD3Tensor_SymD3Tensor_assign(vLSpp%v22, la, na, nc, lb, nb, nd, sumi)
								IF (la .NE. lb) THEN
									CALL SymD3Tensor_SymD3Tensor_assign(vLSpp%v22, lb, nb, nd, la, na, nc, sumi)
								END IF

							END DO
						END DO
					END DO

				END DO
			END DO
		END DO
		CLOSE(file_desc)

		SymVLSpp_read = .TRUE.
		RETURN
	END FUNCTION SymVLSpp_read

	SUBROUTINE SymVLSpp_del(vLSpp)
		TYPE (SymVLSpp), INTENT(INOUT) :: vLSpp

		CALL SymD3Tensor_SymD3Tensor_del(vLSpp%v22)
		RETURN
	END SUBROUTINE SymVLSpp_del

END MODULE symvls
