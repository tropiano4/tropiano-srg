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

 MODULE symd3t_proj

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

			t%d3tensor(u1)%d2 = CMPLX(0,0)

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
					t%d3tensor(u1)%d2(u2, u3) = dummy*CMPLX(1,0)
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

		SymD3Tensor_trace = 0.0

		DO la = 0, Lmax
								d = Min(Nmax, NmaxOfL(la))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - la) / 2) + 1

			DO na = 1, d
				SymD3Tensor_trace = SymD3Tensor_trace + t_in%d3tensor(la)%d2(na, na)
			END DO
		END DO
		RETURN
	END FUNCTION SymD3Tensor_trace

	! Function which gives as output a square matrix containing the matrix of a given tensor
	! operator in the harmonic oscillator basis

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
					SymD3Tensor_matrix(u1, u2) = CMPLX(0.0, 0.0)
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
		DOUBLE PRECISION SymD3Tensor_product2
		TYPE (SymD3Tensor), INTENT(IN) :: t1_in, t2_in

		DOUBLE PRECISION :: sum_la, sum1
		INTEGER :: u1, u2, u3, d

		SymD3Tensor_product2 = 0.0

		DO u1 = 0, Lmax
								d = Min(Nmax, NmaxOfL(u1))
			IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	d = ((N_0 - u1) / 2) + 1

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

		INTEGER :: u1, u2, u3, d
		DOUBLE PRECISION :: distance

		SymD3Tensor_distance = 0.0

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

END MODULE symd3t_proj
