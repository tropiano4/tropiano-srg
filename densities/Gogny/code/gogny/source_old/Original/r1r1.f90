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
		func_out%x = 0 !TODO func_in%x
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
