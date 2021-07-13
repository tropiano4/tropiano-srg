MODULE symd3t

	USE input

	IMPLICIT NONE

	TYPE SymD2Tensor
		DOUBLE PRECISION, DIMENSION(:, :), POINTER :: d2
	END TYPE

	TYPE SymD3Tensor
		TYPE (SymD2Tensor), DIMENSION(:), POINTER :: d3tensor
	END TYPE

	TYPE SymD2Tensor_SymD3Tensor
		TYPE (SymD3Tensor), DIMENSION(:, :), POINTER :: d2
	END TYPE

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

	! Recibe como parámetro de entrada el tensor que vamos a inicializar
	SUBROUTINE SymD3Tensor_new(t)
		TYPE (SymD3Tensor), INTENT(INOUT) :: t

		INTEGER u1, d

		ALLOCATE(t%d3tensor(0:N_0))
		DO u1 = 0, N_0
			d = ((N_0 - u1) / 2) + 1
			ALLOCATE(t%d3tensor(u1)%d2(d, d))
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_new

	FUNCTION SymD3Tensor_read(t, file_desc, file_error)
		LOGICAL SymD3Tensor_read
		TYPE (SymD3Tensor), INTENT(INOUT) :: t
		INTEGER, INTENT(IN) :: file_desc
		INTEGER, INTENT(INOUT) :: file_error

		INTEGER u1, u2, u3, d
		INTEGER v1, v2, v3

		DO u1 = 0, N_0
			d = ((N_0 - u1) / 2) + 1
			DO u2 = 1, d
				DO u3 = 1, u2
					READ (file_desc, FMT="(I3,I3,I3,E)", IOSTAT=file_error) &
						v1, v2, v3, t%d3tensor(u1)%d2(u2, u3)
					IF ((file_error .NE. 0) .OR. &
					    (u1 .NE. v1) .OR. (u2 .NE. v2) .OR. (u3 .NE. v3)) THEN
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

		INTEGER u1, u2, u3, d

		DO u1 = 0, N_0
			d = ((N_0 - u1) / 2) + 1
			DO u2 = 1, d
				DO u3 = 1, u2
					WRITE (file_desc, FMT="(I3,I3,I3,E)", IOSTAT=file_error) &
						u1, u2, u3, t%d3tensor(u1)%d2(u2, u3)
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_write

	FUNCTION SymD3Tensor_trace(t_in)
		DOUBLE PRECISION SymD3Tensor_trace
		TYPE (SymD3Tensor), INTENT(IN) :: t_in

		INTEGER la, na, d

		SymD3Tensor_trace = 0.0
		DO la = 0, N_0
			d = ((N_0 - la) / 2) + 1
			DO na = 1, d
				SymD3Tensor_trace = SymD3Tensor_trace + t_in%d3tensor(la)%d2(na, na)
			END DO
		END DO
		RETURN
	END FUNCTION SymD3Tensor_trace

	FUNCTION SymD3Tensor_matrix(t_in, la)
		DOUBLE PRECISION, DIMENSION(((N_0 - la) / 2) + 1, ((N_0 - la) / 2) + 1) :: SymD3Tensor_matrix
		TYPE (SymD3Tensor), INTENT(INOUT) :: t_in
		INTEGER, INTENT(IN) :: la

		INTEGER u1, u2, d

		! Reconstruimos la matriz del tensor a partir de
		! la matriz triangular inferior
		d = ((N_0 - la) / 2) + 1
		DO u1 = 1, d
			DO u2 = 1, u1
				SymD3Tensor_matrix(u1, u2) = t_in%d3tensor(la)%d2(u1, u2)
				IF (u1 .NE. u2) THEN
					SymD3Tensor_matrix(u2, u1) = t_in%d3tensor(la)%d2(u1, u2)
				END IF
			END DO
		END DO
		RETURN
	END FUNCTION SymD3Tensor_matrix

	SUBROUTINE SymD3Tensor_assign1(t_out, t_in)
		TYPE (SymD3Tensor), INTENT(INOUT) :: t_out
		TYPE (SymD3Tensor), INTENT(IN) :: t_in

		INTEGER u1, u2, u3, d

		DO u1 = 0, N_0
			d = ((N_0 - u1) / 2) + 1
			DO u2 = 1, d
				DO u3 = 1, u2
					t_out%d3tensor(u1)%d2(u2, u3) = t_in%d3tensor(u1)%d2(u2, u3)
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_assign1

	SUBROUTINE SymD3Tensor_assign2(t_out, R1)
		TYPE (SymD3Tensor), INTENT(INOUT) :: t_out
		DOUBLE PRECISION, INTENT(IN) :: R1

		INTEGER u1, u2, u3, d

		DO u1 = 0, N_0
			d = ((N_0 - u1) / 2) + 1
			DO u2 = 1, d
				DO u3 = 1, u2
					t_out%d3tensor(u1)%d2(u2, u3) = R1
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_assign2

	SUBROUTINE SymD3Tensor_assign_Matrix(t_out, u1, M)
		TYPE (SymD3Tensor), INTENT(INOUT) :: t_out
		INTEGER, INTENT(IN) :: u1
		DOUBLE PRECISION, DIMENSION(:, :), INTENT(IN) :: M

		INTEGER u2, u3, d

		d = ((N_0 - u1) / 2) + 1
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

		INTEGER u1, u2, u3, d

		DO u1 = 0, N_0
			d = ((N_0 - u1) / 2) + 1
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

		INTEGER u1, u2, u3, d

		DO u1 = 0, N_0
			d = ((N_0 - u1) / 2) + 1
			DO u2 = 1, d
				DO u3 = 1, u2
					t_out%d3tensor(u1)%d2(u2, u3) = R1 * t_in%d3tensor(u1)%d2(u2, u3)
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_product

	FUNCTION SymD3Tensor_product2(t1_in, t2_in)
		DOUBLE PRECISION SymD3Tensor_product2
		TYPE (SymD3Tensor), INTENT(IN) :: t1_in, t2_in

		INTEGER u1, u2, u3, d
		DOUBLE PRECISION sum_la, sum1

		SymD3Tensor_product2 = 0.0
		DO u1 = 0, N_0
			d = ((N_0 - u1) / 2) + 1
			sum_la = 0.0
			DO u2 = 1, d
				DO u3 = 1, u2
					sum1 = t1_in%d3tensor(u1)%d2(u2, u3) * t2_in%d3tensor(u1)%d2(u2, u3)
					IF (u2 .EQ. u3) THEN
						sum_la = sum_la + sum1
					ELSE
						sum_la = sum_la + (2.0 * sum1)
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

		INTEGER u1, u2, u3, d
		DOUBLE PRECISION distance

		SymD3Tensor_distance = 0.0
		DO u1 = 0, N_0
			d = ((N_0 - u1) / 2) + 1
			DO u2 = 1, d
				DO u3 = 1, u2
!					SymD3Tensor_distance = MAX(SymD3Tensor_distance, ABS(t1_in%d3tensor(u1)%d2(u2, u3) - t2_in%d3tensor(u1)%d2(u2, u3)))
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

		INTEGER u1, u2, u3, d

		DO u1 = 0, N_0
			PRINT "(I2,A)", u1, ":"
			d = ((N_0 - u1) / 2) + 1
			DO u2 = 1, d
				PRINT "(13F10.5)", (t%d3tensor(u1)%d2(u2, u3), u3=1, u2)
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_print

	SUBROUTINE SymD3Tensor_del(t)
		TYPE (SymD3Tensor), INTENT(IN) :: t

		INTEGER u1

		DO u1 = 0, N_0
			DEALLOCATE(t%d3tensor(u1)%d2)
		END DO
		DEALLOCATE(t%d3tensor)
		RETURN
	END SUBROUTINE SymD3Tensor_del

	! Reserva la memoria utilizada por un tensor 3d
	! e inicializa el mismo a 0
	SUBROUTINE SymD3Tensor_SymD3Tensor_new(t)
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(INOUT) :: t

		INTEGER u1, u2, u3, d

		ALLOCATE(t%d3tensor(0:N_0))
		DO u1 = 0, N_0
			d = ((N_0 - u1) / 2) + 1
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

	SUBROUTINE SymD3Tensor_SymD3Tensor_assign1(t_out, t_in)
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(INOUT) :: t_out
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(IN) :: t_in

		INTEGER u1, u2, u3, d

		DO u1 = 0, N_0
			d = ((N_0 - u1) / 2) + 1
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

		INTEGER u1, u2, u3, d

		DO u1 = 0, N_0
			d = ((N_0 - u1) / 2) + 1
			DO u2 = 1, d
				DO u3 = 1, u2
					t_out%d3tensor(u1)%d2(u2, u3) = R1
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_SymD3Tensor_assign2

	SUBROUTINE SymD3Tensor_SymD3Tensor_add(t_out, t1_in, t2_in)
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(INOUT) :: t_out
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(IN) :: t1_in, t2_in

		INTEGER u1, u2, u3, d

		DO u1 = 0, N_0
			d = ((N_0 - u1) / 2) + 1
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

	SUBROUTINE SymD3Tensor_SymD3Tensor_product(t_out, t1_in, t2_in)
		TYPE (SymD3Tensor), INTENT(INOUT) :: t_out
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(IN) :: t1_in
		TYPE (SymD3Tensor), INTENT(IN) :: t2_in

		INTEGER u1, u2, u3, d

		DO u1 = 0, N_0
			d = ((N_0 - u1) / 2) + 1
			DO u2 = 1, d
				DO u3 = 1, u2
					t_out%d3tensor(u1)%d2(u2, u3) = t1_in%d3tensor(u1)%d2(u2, u3) * t2_in
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_SymD3Tensor_product

	SUBROUTINE SymD3Tensor_SymD3Tensor_del(t)
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(IN) :: t

		INTEGER u1, u2, u3, d

		DO u1 = 0, N_0
			d = ((N_0 - u1) / 2) + 1
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
