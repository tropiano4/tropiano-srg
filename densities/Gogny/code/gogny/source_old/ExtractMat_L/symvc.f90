!---------------------------------------------------------------------!
!                                                                     !
!     MATRIX ELEMENTS OF THE COULOMB PART OF THE INTERACTION          !                
!                                                                     !
!---------------------------------------------------------------------!

 MODULE symvc

	USE symd3t

	IMPLICIT NONE

	TYPE SymVCph
		TYPE (SymD3Tensor_SymD3Tensor) v_local, v1_exch, v2_exch, v1
		CHARACTER(64) filename
	END TYPE

	TYPE SymVCpp
		TYPE (SymD3Tensor_SymD3Tensor) v1_pair, v2_pair
		CHARACTER(64) filename
	END TYPE

 	! File units
	
	INTEGER :: file_in = 16, file_out = 17
	
CONTAINS

	!-----------------------------------------------------------------------!
	!  Initialization of the tensors and definition of the filenames (p.h.)	!               
	!-----------------------------------------------------------------------!

	SUBROUTINE SymVCph_new(vCph, Lvalue, Nmax)
		TYPE (SymVCph), INTENT(INOUT) :: vCph
                INTEGER, INTENT(IN) :: Lvalue, Nmax
                
		CALL SymD3Tensor_SymD3Tensor_new(vCph%v_local, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_new(vCph%v1_exch, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_new(vCph%v2_exch, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_new(vCph%v1, Lvalue, Nmax)
		
		RETURN
	END SUBROUTINE SymVCph_new
	
	!-----------------------------------------------------------------------!
	!  Initialization of the tensors and definition of the filenames (p.h.)	!               
	!-----------------------------------------------------------------------!

	SUBROUTINE SymVCph_OpenFiles(Lmax, Lsmall)
		INTEGER, INTENT(IN) :: Lmax, Lsmall
		INTEGER :: file_error
		
		CHARACTER(LEN = 64) :: filein, fileout
      
		IF (Lmax < 10) THEN
			WRITE(filein, "(A,I1,A)") "data/vC", Lmax, "ph_WS.txt"
		ELSE
			WRITE(filein, "(A,I2,A)") "data/vC", Lmax, "ph_WS.txt"
		END IF
		
		IF (Lsmall < 10) THEN
			WRITE(fileout, "(A,I1,A)") "data/vC", Lsmall, "ph_WS.txt"
		ELSE
			WRITE(fileout, "(A,I2,A)") "data/vC", Lsmall, "ph_WS.txt"
		END IF
                write(*,'("filein = ",a64)') filein
                write(*,'("fileout = ",a64)') fileout
                        		
		OPEN(file_in, FILE=filein, ACTION="READ", IOSTAT=file_error)
		OPEN(file_out, FILE=fileout, ACTION="WRITE", IOSTAT=file_error)	
			
		RETURN
		
		RETURN
	END SUBROUTINE SymVCph_OpenFiles

	SUBROUTINE SymVCph_CloseFiles()

		CLOSE(file_in)
		CLOSE(file_out)
		
		RETURN

	END SUBROUTINE SymVCph_CloseFiles

	SUBROUTINE SymVCph_write(vCph, Lvalue, Nmax)
		TYPE (SymVCph), INTENT(IN) :: vCph
		INTEGER, INTENT(IN) :: Lvalue, Nmax

		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd, file_error
		DOUBLE PRECISION :: sumi, sum1, sum2

		DO la = 0, Lvalue
		
			namax = Nmax
                                
			DO na = 1, namax
				DO nc = 1, na
				
					DO lb = 0, la

						nbmax = Nmax
							
						DO nb = 1, nbmax
							DO nd = 1, nb
								
								CALL SymD3Tensor_SymD3Tensor_get(sumi, vCph%v_local, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum1, vCph%v1_exch, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum2, vCph%v2_exch, la, na, nc, lb, nb, nd)
								
								!IF (file_error .EQ. 0) THEN
									WRITE (file_out, FMT="(6I3,3E24.16)", IOSTAT=file_error) &
									la, na, nc, lb, nb, nd, sumi, sum1, sum2
								!END IF
							END DO
						END DO
					END DO
				END DO
			END DO
			
		END DO
		
		RETURN

	END SUBROUTINE SymVCph_write
	
	SUBROUTINE SymVCph_read(vCph, Lvalue, Nmax)
		TYPE (SymVCph), INTENT(INOUT) :: vCph
		INTEGER, INTENT(IN) :: Lvalue, Nmax
		
		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER i1, i2, i3, i4, i5, i6, file_error
		DOUBLE PRECISION sumi, sum1, sum2

		DO la = 0, Lvalue

			namax = Nmax
			
			DO na = 1, namax
				DO nc = 1, na
				
					DO lb = 0, la

						nbmax = Nmax
						
						DO nb = 1, nbmax
							DO nd = 1, nb
							
								READ (file_in, FMT=*, IOSTAT=file_error) &
									i1, i2, i3, i4, i5, i6, sumi, sum1, sum2

								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v_local, la, na, nc, lb, nb, nd, sumi)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v_local, lb, nb, nd, la, na, nc, sumi)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v1_exch, la, na, nc, lb, nb, nd, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v1_exch, lb, nb, nd, la, na, nc, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v2_exch, la, na, nc, lb, nb, nd, sum2)
								CALL SymD3Tensor_SymD3Tensor_assign(vCph%v2_exch, lb, nb, nd, la, na, nc, sum2)
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
		
		RETURN
	END SUBROUTINE SymVCph_read

	SUBROUTINE SymVCph_del(vCph, Lvalue, Nmax)
		TYPE (SymVCph), INTENT(INOUT) :: vCph
		INTEGER, INTENT(IN) :: Lvalue, Nmax

		CALL SymD3Tensor_SymD3Tensor_del(vCph%v_local, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_del(vCph%v1_exch, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_del(vCph%v2_exch, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_del(vCph%v1, Lvalue, Nmax)
		
		RETURN
	END SUBROUTINE SymVCph_del

	!-----------------------------------------------------------------------!
	!  Initialization of the tensors and definition of the filenames (p.p.)	!               
	!-----------------------------------------------------------------------!

	SUBROUTINE SymVCpp_new(vCpp, Lvalue, Nmax)
		TYPE (SymVCpp), INTENT(INOUT) :: vCpp
                INTEGER, INTENT(IN) :: Lvalue, Nmax
                
		CALL SymD3Tensor_SymD3Tensor_new(vCpp%v1_pair, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_new(vCpp%v2_pair, Lvalue, Nmax)
		
		RETURN
	END SUBROUTINE SymVCpp_new

	SUBROUTINE SymVCpp_OpenFiles(Lmax, Lsmall)
		INTEGER, INTENT(IN) :: Lmax, Lsmall
		INTEGER :: file_error
		
		CHARACTER(LEN = 64) :: filein, fileout
      
		IF (Lmax < 10) THEN
			WRITE(filein, "(A,I1,A)") "data/vC", Lmax, "pp_WS.txt"
		ELSE
			WRITE(filein, "(A,I2,A)") "data/vC", Lmax, "pp_WS.txt"
		END IF
		
		IF (Lsmall < 10) THEN
			WRITE(fileout, "(A,I1,A)") "data/vC", Lsmall, "pp_WS.txt"
		ELSE
			WRITE(fileout, "(A,I2,A)") "data/vC", Lsmall, "pp_WS.txt"
		END IF
                        
		OPEN(file_in, FILE=filein, ACTION="READ", IOSTAT=file_error)
		OPEN(file_out, FILE=fileout, ACTION="WRITE", IOSTAT=file_error)	
			
		RETURN

	END SUBROUTINE SymVCpp_OpenFiles
	
	SUBROUTINE SymVCpp_CloseFiles()

		CLOSE(file_in)
		CLOSE(file_out)		
		
		RETURN

	END SUBROUTINE SymVCpp_CloseFiles

	SUBROUTINE SymVCpp_write(vCpp, Lvalue, Nmax)
		TYPE (SymVCpp), INTENT(INOUT) :: vCpp
		INTEGER, INTENT(IN) :: Lvalue, Nmax

		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd, file_error
		DOUBLE PRECISION :: sum1, sum2

		DO la = 0, Lvalue
		
			namax = Nmax
                                
			DO na = 1, namax
				DO nc = 1, na
				
					DO lb = 0, la

						nbmax = Nmax
							
						DO nb = 1, nbmax
							DO nd = 1, nb
								
								CALL SymD3Tensor_SymD3Tensor_get(sum1, vCpp%v1_pair, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum2, vCpp%v2_pair, la, na, nc, lb, nb, nd)
								
								!IF (file_error .EQ. 0) THEN
									WRITE (file_out, FMT="(6I3,2E24.16)", IOSTAT=file_error) &
									la, na, nc, lb, nb, nd, sum1, sum2
								!END IF
							END DO
						END DO
					END DO
				END DO
			END DO
			
		END DO
		
		RETURN

	END SUBROUTINE SymVCpp_write
	
	SUBROUTINE SymVCpp_read(vCpp, Lvalue, Nmax)
		TYPE (SymVCpp), INTENT(INOUT) :: vCpp
		INTEGER, INTENT(IN) :: Lvalue, Nmax

		INTEGER la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER i1, i2, i3, i4, i5, i6, file_error
		DOUBLE PRECISION sum1, sum2

		DO la = 0, Lvalue

			namax = Nmax
			
			DO na = 1, namax
				DO nc = 1, na
				
					DO lb = 0, la

						nbmax = Nmax
						
						DO nb = 1, nbmax
							DO nd = 1, nb
								READ (file_in, FMT=*, IOSTAT=file_error) &
									i1, i2, i3, i4, i5, i6, sum1, sum2

																	CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v1_pair, la, na, nc, lb, nb, nd, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v1_pair, lb, nb, nd, la, na, nc, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v2_pair, la, na, nc, lb, nb, nd, sum2)
								CALL SymD3Tensor_SymD3Tensor_assign(vCpp%v2_pair, lb, nb, nd, la, na, nc, sum2)
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO

				RETURN
	END SUBROUTINE SymVCpp_read

	SUBROUTINE SymVCpp_del(vCpp, Lvalue, Nmax)
		TYPE (SymVCpp), INTENT(INOUT) :: vCpp
		INTEGER, INTENT(IN) :: Lvalue, Nmax

		CALL SymD3Tensor_SymD3Tensor_del(vCpp%v1_pair, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_del(vCpp%v2_pair, Lvalue, Nmax)
		RETURN
	END SUBROUTINE SymVCpp_del

END MODULE symvc
