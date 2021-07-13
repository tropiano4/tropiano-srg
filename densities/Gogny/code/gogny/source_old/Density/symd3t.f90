MODULE symd3t

	TYPE SymD2Tensor
		DOUBLE PRECISION, DIMENSION(:, :), POINTER :: d2
	END TYPE

	TYPE SymD3Tensor
		TYPE (SymD2Tensor), DIMENSION(:), POINTER :: d3tensor
	END TYPE
	
	INTEGER :: N_0, Lmax, Nmax, Basis, CompHO
	INTEGER :: N_0_new, L_new, N_new

 CONTAINS
 
	SUBROUTINE SymD3Tensor_new(t)
		TYPE (SymD3Tensor), INTENT(INOUT) :: t

		INTEGER u1, d

                ! Allocating memory for the vector of SymD2Tensor tensors
		ALLOCATE(t%d3tensor(0:Lmax))
		
                ! For each SymD2Tensor of the aforementioned vector, allocating memory for the 
		! corresponding array
		DO u1 = 0, Lmax
								d = Nmax
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

			ALLOCATE(t%d3tensor(u1)%d2(d, d))
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_new

	SUBROUTINE Density_read(t, file_desc, file_error)
		TYPE (SymD3Tensor), DIMENSION(:, :), POINTER :: t
		INTEGER, INTENT(IN) :: file_desc
		INTEGER, INTENT(INOUT) :: file_error

		INTEGER :: u1, u2, u3, d
		INTEGER :: v1, v2, v3
		INTEGER :: ta, TypeDens
		
		OPEN(file_desc, FILE="density.txt", ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			RETURN
		END IF
		
		READ (file_desc, FMT=*, IOSTAT=file_error) N_0, Lmax, Nmax, Basis, CompHO
		
		ALLOCATE(t(0:1, 0:3))
		
		DO TypeDens = 0, 3
			
			DO ta = 0, 1
		
				CALL SymD3Tensor_new(t(ta, TypeDens))
		
				DO u1 = 0, Lmax
										d = Nmax
					IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

					DO u2 = 1, d
						DO u3 = 1, u2
							READ (file_desc, FMT=*, IOSTAT=file_error) &
								v1, v2, v3, t(ta, TypeDens)%d3tensor(u1)%d2(u2, u3)
							IF ((file_error .NE. 0) .OR. (u1 .NE. v1) .OR. (u2 .NE. v2) .OR. (u3 .NE. v3)) THEN
								RETURN
							END IF
						END DO
					END DO
				END DO
				
			END DO
		END DO
		
		RETURN
	END SUBROUTINE Density_read


	SUBROUTINE Density_write(t, file_desc, file_error)
		TYPE (SymD3Tensor), DIMENSION(:, :), POINTER :: t
		INTEGER, INTENT(IN) :: file_desc
		INTEGER, INTENT(INOUT) :: file_error

		DOUBLE PRECISION :: Initialize = 1.e-4, WhatToWrite
		INTEGER :: u1, u2, u3, d
		INTEGER :: ta, TypeDens

		OPEN(file_desc, FILE="density_extended.txt", ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			WRITE(*,'("Error in Opening the file density_extended.txt")')
			RETURN
		END IF
		
		WRITE (file_desc, FMT="(5I3)", IOSTAT=file_error) N_0_new, L_new, N_new, Basis, CompHO
		
		DO TypeDens = 0, 3
			
			DO ta = 0, 1
		
				DO u1 = 0, L_new
										d = N_new
					IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

					DO u2 = 1, d
						DO u3 = 1, u2
						
							IF (u1 .GT. Lmax) THEN
								WhatToWrite = Initialize
							ELSE
								IF (u2 .GT. Nmax) THEN
									WhatToWrite = Initialize
								ELSE
									WhatToWrite = t(ta, TypeDens)%d3tensor(u1)%d2(u2, u3)
								END IF								
							END IF
						
							WRITE (file_desc, FMT="(3I3,E24.16)", IOSTAT=file_error) &
								u1, u2, u3, WhatToWrite
						END DO
					END DO
				END DO
				
			END DO
		END DO
		
		RETURN
	END SUBROUTINE Density_write

	
	
END MODULE symd3t
