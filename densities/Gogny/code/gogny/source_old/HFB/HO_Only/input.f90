! Input parameters module
 MODULE input

	IMPLICIT NONE

	! More variables (input parameters)
	
	DOUBLE PRECISION :: b_0 = 1.6       	! Oscillator length ("b")
	INTEGER :: N_0 = 8                   	! Number of shells
	INTEGER :: protons = 8, neutrons = 8    ! Number of protons and neutrons
	INTEGER :: Gogny = 1        		! 1 = D1S, 2 = D1, 3 = D1prime
						! 0 = Lmax and Nmax required
	
        DOUBLE PRECISION :: Convergence = 1.e-3, anneal = 0.5
	
	INTEGER :: HFOnly = 0, NITER_MAX = 10
			
 CONTAINS

	!-----------------------------------------------------------------------!
	!									!
	!       This subroutine reads the file "input.txt" that contains 	!
	!       the standard input to the program.				!
	!									!
	!-----------------------------------------------------------------------!
	
	SUBROUTINE Input_read

		INTEGER, PARAMETER :: file_desc = 17
		INTEGER file_error, N, L, i
		CHARACTER(LEN = 256) :: param
		INTEGER param_len

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
			ELSE IF (param .EQ. "N_0") THEN
				N_0 = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "N") THEN
				neutrons = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "Z") THEN
				protons = IO_get_IntegerValue(file_desc)
			ELSE IF (param .EQ. "GOGNY") THEN
				Gogny = IO_get_IntegerValue(file_desc)
			ELSE
				PRINT "(A,A)", "Attention! Unknown parameter: ", param(1:param_len)
				! We ignore the remaining lines
				READ (file_desc, "(A)") param
			END IF
			param_len = IO_get_Param(file_desc, param)
		END DO

		CLOSE(file_desc)
		
		! In case of a basis compatible with the harmonic oscillator basis, one forces
		! b (the oscillator length) to take a value proportional to the mass
		
		IF (b_0 .LT. 0.0) b_0 = 1.05*REAL(protons + neutrons)**(1.0/6.0)
		
		WRITE(*,'(5X,"DATA READ FROM INPUT")')
		WRITE(*,'(5X,"====================",/)')
		
		PRINT "(A,F8.5)", "Oscillator Length (b_0)       : ", b_0
		PRINT "(A,F8.5)", "Slowing factor                : ", anneal
		PRINT "(A,I3)",   "Maximum Number of Iterations  : ", NITER_MAX
		PRINT "(A,ES9.2)","Convergence in Energy         : ", Convergence
		PRINT "(A,I3)",   "Number of shells (N_0)        : ", N_0
		PRINT "(A,I3)",   "Number of neutrons (N)        : ", neutrons
		PRINT "(A,I3)",   "Number of protons  (Z)        : ", protons
		SELECT CASE (Gogny)
			CASE (1)
				PRINT "(A)", "Parameters of the Gogny force : D1S"
			CASE (2)
				PRINT "(A)", "Parameters of the Gogny force : D1"
			CASE (3)
				PRINT "(A)", "Parameters of the Gogny force : D1prime"
			CASE DEFAULT
				STOP "Invalid parameters of the Gogny force"
		END SELECT		
		
		RETURN
		
	END SUBROUTINE Input_read
	
	FUNCTION IO_get_Param(fd, param)
		INTEGER IO_get_Param
		INTEGER, INTENT(IN) :: fd ! Descriptor de fichero
		CHARACTER(*), INTENT(INOUT) :: param

		CHARACTER c
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

		READ (fd, FMT=*, IOSTAT=stat) IO_get_IntegerValue
		IF (stat .NE. 0) STOP "Imposible leer parametro de entrada"
		RETURN
	END FUNCTION IO_get_IntegerValue

	FUNCTION IO_get_RealValue(fd)
		DOUBLE PRECISION IO_get_RealValue
		INTEGER, INTENT(IN) :: fd ! Descriptor de fichero

		INTEGER stat

		READ (fd, FMT=*, IOSTAT=stat) IO_get_RealValue
		IF (stat .NE. 0) STOP "Imposible leer parametro de entrada"
		RETURN
	END FUNCTION IO_get_RealValue

	FUNCTION IO_get_String(fd)
		CHARACTER(LEN=30) IO_get_String
		INTEGER, INTENT(IN) :: fd ! Descriptor de fichero

		INTEGER stat

		READ (fd, FMT=*, IOSTAT=stat) IO_get_String
		IF (stat .NE. 0) STOP "Imposible leer parametro de entrada"
		RETURN
	END FUNCTION IO_get_String

END MODULE input
