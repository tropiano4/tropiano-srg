MODULE symd3t

	TYPE SymD2Tensor
		DOUBLE PRECISION, DIMENSION(:, :), POINTER :: d2
	END TYPE

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

 CONTAINS
 
	SUBROUTINE SymD3Tensor_new(t, Lmax, Nmax)
		TYPE (SymD3Tensor), INTENT(INOUT) :: t
                INTEGER, INTENT(IN) :: Lmax, Nmax

		INTEGER u1, d

                ! Allocating memory for the vector of SymD2Tensor tensors
		ALLOCATE(t%d3tensor(0:Lmax))
		
                ! For each SymD2Tensor of the aforementioned vector, allocating memory for the 
		! corresponding array
		DO u1 = 0, Lmax

			d = Nmax

			ALLOCATE(t%d3tensor(u1)%d2(d, d))
		END DO
		RETURN
	END SUBROUTINE SymD3Tensor_new

	SUBROUTINE SymD3Tensor_del(t, Lmax, Nmax)
		TYPE (SymD3Tensor), INTENT(IN) :: t
                INTEGER, INTENT(IN) :: Lmax, Nmax

		INTEGER u1

		DO u1 = 0, Lmax
			DEALLOCATE(t%d3tensor(u1)%d2)
		END DO
		DEALLOCATE(t%d3tensor)
		RETURN
	END SUBROUTINE SymD3Tensor_del

	! Reserves (allocates) the memory for a "matrix" of 3D tensors
	
	SUBROUTINE SymD3Tensor_SymD3Tensor_new(t, Lmax, Nmax)
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(INOUT) :: t
                INTEGER, INTENT(IN) :: Lmax, Nmax

		INTEGER u1, u2, u3, d

		ALLOCATE(t%d3tensor(0:Lmax))
		
		DO u1 = 0, Lmax
                
			d = Nmax

			ALLOCATE(t%d3tensor(u1)%d2(d, d))
			
			DO u2 = 1, d
				DO u3 = 1, u2
					CALL SymD3Tensor_new(t%d3tensor(u1)%d2(u2, u3), Lmax, Nmax)
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
	
	SUBROUTINE SymD3Tensor_SymD3Tensor_get(R1, t_in, la, na, nc, lb, nb, nd)
		DOUBLE PRECISION, INTENT(OUT) :: R1
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(IN) :: t_in
		INTEGER, INTENT(IN) :: la, na, nc, lb, nb, nd

		R1 = t_in%d3tensor(la)%d2(na, nc)%d3tensor(lb)%d2(nb, nd)
		RETURN
	END SUBROUTINE SymD3Tensor_SymD3Tensor_get

	SUBROUTINE SymD3Tensor_SymD3Tensor_del(t, Lmax, Nmax)
		TYPE (SymD3Tensor_SymD3Tensor), INTENT(IN) :: t
                INTEGER, INTENT(IN) :: Lmax, Nmax

		INTEGER u1, u2, u3, d

		DO u1 = 0, Lmax
                
			d = Nmax

			DO u2 = 1, d
				DO u3 = 1, u2
					CALL SymD3Tensor_del(t%d3tensor(u1)%d2(u2, u3), Lmax, Nmax)
				END DO
			END DO
			DEALLOCATE(t%d3tensor(u1)%d2)
		END DO
		DEALLOCATE(t%d3tensor)
		RETURN
	END SUBROUTINE SymD3Tensor_SymD3Tensor_del

	
END MODULE symd3t
