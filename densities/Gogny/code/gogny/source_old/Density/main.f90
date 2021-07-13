PROGRAM main

 	USE symd3t
 
	IMPLICIT NONE
 
	TYPE (SymD3Tensor), DIMENSION(:, :), POINTER :: t

	INTEGER :: file_desc = 16, file_error

	READ(*,*) N_0_new, L_new, N_new
	
	CALL Density_read(t, file_desc, file_error)

	CALL Density_write(t, file_desc, file_error)
	
	STOP
 
END PROGRAM main

