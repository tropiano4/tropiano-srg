! Gestiona los parámetros de entrada para la ejecución del programa
MODULE input

	USE io

	IMPLICIT NONE

	! Parámetros de entrada
	DOUBLE PRECISION :: b_0 = 1.66       ! Longitud del oscilador
	INTEGER :: N_0 = 8                   ! Número de capas
	INTEGER :: protons = 8, neutrons = 8 ! Número de protoles y neutrones
	INTEGER :: Gogny = 1                 ! 1 = DSA1, 2 = D1, 3 = D1_PRIMA

CONTAINS

	SUBROUTINE Input_read

		INTEGER, PARAMETER :: file_desc = 6
		INTEGER file_error
		CHARACTER(LEN = 256) :: param
		INTEGER param_len

		OPEN(file_desc, FILE="data/input.txt", ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) STOP "Imposible leer: input.txt"

		PRINT "(/A)", "Leyendo los parametros de entrada..."
		param_len = IO_get_Param(file_desc, param)
		DO WHILE (param_len .GT. 0)
			IF (param .EQ. "b_0") THEN
				b_0 = IO_get_RealValue(file_desc)
				PRINT "(A,F8.5)", "Longitud del oscilador (b_0) = ", b_0
			ELSE IF (param .EQ. "N_0") THEN
				N_0 = IO_get_IntegerValue(file_desc)
				PRINT "(A,I2)", "Numero de capas (N_0) = ", N_0
			ELSE IF (param .EQ. "N") THEN
				protons = IO_get_IntegerValue(file_desc)
				PRINT "(A,I3)", "Numero de neutrones (N) = ", protons
			ELSE IF (param .EQ. "Z") THEN
				neutrons = IO_get_IntegerValue(file_desc)
				PRINT "(A,I3)", "Numero de protones (Z) = ", neutrons
			ELSE IF (param .EQ. "GOGNY") THEN
				Gogny = IO_get_IntegerValue(file_desc)
				PRINT "(A,I1)", "Parametro de la fuerza de Gogny = ", Gogny
				SELECT CASE (Gogny)
					CASE (1)
						PRINT "(A)", "Gogny = DSA1"
					CASE (2)
						PRINT "(A)", "Gogny = D1"
					CASE (3)
						PRINT "(A)", "Gogny = D1_PRIMA"
					CASE DEFAULT
						STOP "Parametro de la fuerza de Gogny no valido"
				END SELECT
			ELSE
				PRINT "(A,A)", "AVISO! Parametro no reconocido: ", param(1:param_len)
				! Ignoramos el resto de línea
				READ (file_desc, "(A)") param
			END IF
			param_len = IO_get_Param(file_desc, param)
		END DO

		CLOSE(file_desc)
		RETURN
	END SUBROUTINE Input_read

END MODULE input
