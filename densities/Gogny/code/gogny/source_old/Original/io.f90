MODULE io

	IMPLICIT NONE

CONTAINS

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
		! Si no se encontró ningún nombre de parámetro, salimos con error
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

		READ (fd, FMT="(I)", IOSTAT=stat) IO_get_IntegerValue
		IF (stat .NE. 0) STOP "Imposible leer parametro de entrada"
		RETURN
	END FUNCTION IO_get_IntegerValue

	FUNCTION IO_get_RealValue(fd)
		DOUBLE PRECISION IO_get_RealValue
		INTEGER, INTENT(IN) :: fd ! Descriptor de fichero

		INTEGER stat

		READ (fd, FMT="(E)", IOSTAT=stat) IO_get_RealValue
		IF (stat .NE. 0) STOP "Imposible leer parametro de entrada"
		RETURN
	END FUNCTION IO_get_RealValue

END MODULE io
