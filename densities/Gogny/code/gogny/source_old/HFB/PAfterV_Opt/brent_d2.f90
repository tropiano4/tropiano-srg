	! This is the set of subroutines that determine the particle number
	! at each iteration. Look at subroutine OneDimSolve_LookAt(solve, V)
	! where there is a call to DiagonalizationMethod_operator (which 
	! calculates the particle number).

	SUBROUTINE OneDimSolve_new(solve, func, diagonal)
		TYPE (OneDimSolve), INTENT(INOUT) :: solve
		TYPE (R1R1Function), TARGET, INTENT(IN) :: func
		TYPE (DiagonalizationMethod), TARGET, INTENT(IN) :: diagonal

		solve%func => func
		solve%diagonal => diagonal
		solve%AccX = 1.e-12
		solve%AccY = 1.e-12
		
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
		
		IF (Dev .EQ. 0.0) STOP "Dev es nulo"
		
		solve%L = 0
		solve%C = 1
		solve%U = 2
		solve%vpol = 1
		solve%hpol = 1
		solve%y(solve%C) = 0.0
		solve%y(solve%U) = 0.0

		IF (Dev < 0.0) THEN
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
				IF (solve%y(solve%L) .GT. 0.0) THEN
					CALL OneDimSolve_VFlip(Solve) ! So Y(L) < 0
				END IF

				! Right point
				solve%x(solve%U) = X + (Dev * solve%hpol)
				
				! Check boundaries for right point
				IF ((.NOT. solve%func%infinite_x_max) .AND. (solve%x(solve%U) .GT. solve%func%x_max)) THEN
					solve%x(solve%U) = (solve%func%x_max + X) / 2.0
				END IF
				
				IF ((.NOT. solve%func%infinite_x_min) .AND. (solve%x(solve%U) .LT. solve%func%x_min)) THEN
					solve%x(solve%U) = (solve%func%x_min + X) / 2.0
				END IF

				! Right point: Get value of the function 
				CALL OneDimSolve_LookAt(Solve, solve%U)
				
				IF (solve%Finish) THEN
					State = 4
					CYCLE
				END IF
				
				! If we have the right point positive and the left point negative, we know where the 
				! crossing is (Case State = 1 below)
				IF (solve%y(solve%U) .GT. 0.0) THEN
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
						solve%x(solve%U) = (solve%func%x_max + solve%x(solve%C)) / 2.0
					END IF
					IF ((.NOT. solve%func%infinite_x_min) .AND. (solve%x(solve%U) .LT. solve%func%x_min)) THEN
						solve%x(solve%U) = (solve%func%x_min + solve%x(solve%C)) / 2.0
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
					
					Dev = Dev * 2.0
					
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

				IF (solve%y(solve%C) .GT. 0.0) THEN
					CALL OneDimSolve_Flip(Solve) ! Want y[C] < 0
				END IF

				IF (solve%y(solve%C) .LT. (0.5 * solve%y(solve%L))) THEN
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
				IF ((((solve%x(solve%C) - solve%x(solve%L)) * solve%hpol) .LE. 0.0) .OR. &
				    (((solve%x(solve%C) - solve%x(solve%U)) * solve%hpol) .GE. 0.0)) THEN
					State = 1
					CYCLE
				END IF

				CALL OneDimSolve_LookAt(Solve, solve%C)
				IF (solve%Finish) THEN
					State = 4
					CYCLE
				END IF

				IF (solve%y(solve%C) .GT. 0.0) THEN
					CALL OneDimSolve_Flip(Solve) ! Want y[C] < 0
				END IF

				IF (solve%y(solve%C) .GT. (0.5 * solve%y(solve%L))) THEN
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
				
					solve%x(solve%C) = 0.5 * (solve%x(solve%L) + solve%x(solve%U))
					
					CALL OneDimSolve_LookAt(Solve, solve%C)
					
					IF (solve%Finish) THEN
						State = 4
						CYCLE bucle
					END IF

					IF (solve%y(solve%C) .GT. 0.0) THEN
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
