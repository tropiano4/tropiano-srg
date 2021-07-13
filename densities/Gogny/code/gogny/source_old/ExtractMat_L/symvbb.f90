!----------------------------------------------------------------!
!								 !
!  CALCULATION OF THE BRINK-BOEKER TERM OF THE GOGNY FORCE       !
!								 !
!----------------------------------------------------------------!

 MODULE symvbb

	USE symd3t

	IMPLICIT NONE

	TYPE SymVBBph
		DOUBLE PRECISION b1
		TYPE (SymD3Tensor_SymD3Tensor) &
			v_local_same_part, v_local_diff_part, &
			v1_exch_same_part, v1_exch_diff_part, &
			v2_exch_same_part, v2_exch_diff_part
		TYPE (SymD3Tensor_SymD3Tensor) &
			v1_same_part, v1_diff_part, &
			v2_same_part, v2_diff_part
		CHARACTER(LEN = 64) filename
	END TYPE

	TYPE SymVBBpp
		DOUBLE PRECISION b1
		TYPE (SymD3Tensor_SymD3Tensor) v1_pair, v2_pair
		CHARACTER(LEN = 64) filename
	END TYPE

	! File units
	
	INTEGER :: file_in = 16, file_out = 17
	
 CONTAINS

        !---------------------------------------------------------------------------------------!
	!											!
        !   Creates a new object of the type SymVBBph, that is in fact a collection of 		!
	!   "3D tensors" arrays SymD3Tensor_SymD3Tensor, as defined in module symd3t.f90	!
	!   Here we create each of the corresponding parts of the type SymVBBph and allocate 	!
	!   the appropriate amount of memory							!
	!   											!
	!   This tensor is for the Brink-Boker matrix elements in the particle-hole channel	!
	!											!
        !---------------------------------------------------------------------------------------!
	
	SUBROUTINE SymVBBph_new(vBBph, b, Lvalue, Nmax)
		TYPE (SymVBBph), INTENT(INOUT) :: vBBph
                INTEGER, INTENT(IN) :: Lvalue, Nmax
		DOUBLE PRECISION, INTENT(IN) :: b

		! Create the tensors

		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v_local_same_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v_local_diff_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v1_exch_same_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v1_exch_diff_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v2_exch_same_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v2_exch_diff_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v1_same_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v1_diff_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v2_same_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_new(vBBph%v2_diff_part, Lvalue, Nmax)
		
		vBBph%b1 = b
		
		RETURN		
	END SUBROUTINE SymVBBph_new

	SUBROUTINE SymVBBph_OpenFiles(vBBph, Lmax, Lsmall)
		TYPE (SymVBBph), INTENT(INOUT) :: vBBph
		
		INTEGER, INTENT(IN) :: Lmax, Lsmall
		INTEGER :: file_error
		
		CHARACTER(LEN = 64) :: filein, fileout
      
		IF (Lmax < 10) THEN
			WRITE(filein, "(A,I1,A)") "data/vBB", Lmax, "ph_WS.txt"
		ELSE
			WRITE(filein, "(A,I2,A)") "data/vBB", Lmax, "ph_WS.txt"
		END IF
		
		IF (Lsmall < 10) THEN
			WRITE(fileout, "(A,I1,A)") "data/vBB", Lsmall, "ph_WS.txt"
		ELSE
			WRITE(fileout, "(A,I2,A)") "data/vBB", Lsmall, "ph_WS.txt"
		END IF
                write(*,'("filein = ",a64)') filein
                write(*,'("fileout = ",a64)') fileout
                        
		OPEN(file_in, FILE=filein, ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to read in ", filein
		ELSE
			READ (file_in, FMT=*, IOSTAT=file_error) vBBph%b1
		END IF

		OPEN(file_out, FILE=fileout, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to write in ", fileout
		ELSE
			WRITE(file_out, FMT="(E24.16)", IOSTAT=file_error) vBBph%b1
		END IF
		
		RETURN

	END SUBROUTINE SymVBBph_OpenFiles
	
	SUBROUTINE SymVBBph_CloseFiles()

		CLOSE(file_in)
		CLOSE(file_out)
		
		RETURN

	END SUBROUTINE SymVBBph_CloseFiles

	SUBROUTINE SymVBBph_write(vBBph, Lvalue, Nmax)
		TYPE (SymVBBph), INTENT(IN) :: vBBph
		INTEGER, INTENT(IN) :: Lvalue, Nmax

		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd, file_error
		DOUBLE PRECISION :: sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp

		DO la = 0, Lvalue

			namax = Nmax
                                
			DO na = 1, namax
				DO nc = 1, na
				
					DO lb = 0, la

						nbmax = Nmax
							
						DO nb = 1, nbmax
							DO nd = 1, nb
								
								CALL SymD3Tensor_SymD3Tensor_get(sumi_sp, vBBph%v_local_same_part, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sumi_dp, vBBph%v_local_diff_part, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum1_sp, vBBph%v1_exch_same_part, lb, nb, nd, la, na, nc)
								CALL SymD3Tensor_SymD3Tensor_get(sum1_dp, vBBph%v1_exch_diff_part, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum2_sp, vBBph%v2_exch_same_part, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum2_dp, vBBph%v2_exch_diff_part, lb, nb, nd, la, na, nc)
								
								!IF (file_error .EQ. 0) THEN
									WRITE (file_out, FMT="(6I3,6E24.16)", IOSTAT=file_error) &
									la, na, nc, lb, nb, nd, sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp
								!END IF
							END DO
						END DO
					END DO
				END DO
			END DO
			
		END DO
		
		RETURN

	END SUBROUTINE SymVBBph_write
	
	SUBROUTINE SymVBBph_read(vBBph, Lvalue, Nmax)
		TYPE (SymVBBph), INTENT(INOUT) :: vBBph
		INTEGER, INTENT(IN) :: Lvalue, Nmax

		DOUBLE PRECISION :: b1
		INTEGER :: la, namax, na, nc, lb, nbmax, nb, nd
		INTEGER :: i1, i2, i3, i4, i5, i6, file_error
		DOUBLE PRECISION :: sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp

		DO la = 0, Lvalue

			namax = Nmax
			
			DO na = 1, namax
				DO nc = 1, na
				
					DO lb = 0, la

						nbmax = Nmax
						
						DO nb = 1, nbmax
							DO nd = 1, nb
							
								READ (file_in, FMT=*, IOSTAT=file_error) &
									i1, i2, i3, i4, i5, i6, sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp
								
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_same_part, la, na, nc, lb, nb, nd, sumi_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_same_part, lb, nb, nd, la, na, nc, sumi_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_diff_part, la, na, nc, lb, nb, nd, sumi_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v_local_diff_part, lb, nb, nd, la, na, nc, sumi_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_same_part, la, na, nc, lb, nb, nd, sum1_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_same_part, lb, nb, nd, la, na, nc, sum1_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_diff_part, la, na, nc, lb, nb, nd, sum1_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v1_exch_diff_part, lb, nb, nd, la, na, nc, sum1_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_same_part, la, na, nc, lb, nb, nd, sum2_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_same_part, lb, nb, nd, la, na, nc, sum2_sp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_diff_part, la, na, nc, lb, nb, nd, sum2_dp)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBph%v2_exch_diff_part, lb, nb, nd, la, na, nc, sum2_dp)
								
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
		
		RETURN
	END SUBROUTINE SymVBBph_read

       !---------------------------------------------------------------------------------------!
	!											!
	!   Subroutine deleting the tensors corresponding to the matrix elements of the Brink-	!
	!   Boker in the particle-hole channel.							!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVBBph_del(vBBph, Lvalue, Nmax)
		TYPE (SymVBBph), INTENT(INOUT) :: vBBph
                INTEGER, INTENT(IN) :: Lvalue, Nmax

		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v_local_same_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v_local_diff_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v1_exch_same_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v1_exch_diff_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v2_exch_same_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v2_exch_diff_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v1_same_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v1_diff_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v2_same_part, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_del(vBBph%v2_diff_part, Lvalue, Nmax)
		
		RETURN
	END SUBROUTINE SymVBBph_del

        !---------------------------------------------------------------------------------------!
	!											!
        !   Creates a new object of the type SymVBBph, that is in fact a collection of 		!
	!   "3D tensors" arrays SymD3Tensor_SymD3Tensor, as defined in module symd3t.f90	!
	!   Here we create each of the corresponding parts of the type SymVBBph and allocate 	!
	!   the appropriate amount of memory							!
	!   											!
	!   This tensor is for the Brink-Boker matrix elements in the particle-particle channel	!
	!											!
        !---------------------------------------------------------------------------------------!
	
	SUBROUTINE SymVBBpp_new(vBBpp, b, Lvalue, Nmax)
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp
                INTEGER, INTENT(IN) :: Lvalue, Nmax
		DOUBLE PRECISION, INTENT(IN) :: b

		CALL SymD3Tensor_SymD3Tensor_new(vBBpp%v1_pair, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_new(vBBpp%v2_pair, Lvalue, Nmax)

		vBBpp%b1 = b
					
		RETURN
	END SUBROUTINE SymVBBpp_new

 	SUBROUTINE SymVBBpp_OpenFiles(vBBpp, Lmax, Lsmall)
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp
		
		INTEGER, INTENT(IN) :: Lmax, Lsmall
		INTEGER :: file_error
		
		CHARACTER(LEN = 64) :: filein, fileout
      
		IF (Lmax < 10) THEN
			WRITE(filein, "(A,I1,A)") "data/vBB", Lmax, "pp_WS.txt"
		ELSE
			WRITE(filein, "(A,I2,A)") "data/vBB", Lmax, "pp_WS.txt"
		END IF
		
		IF (Lsmall < 10) THEN
			WRITE(fileout, "(A,I1,A)") "data/vBB", Lsmall, "pp_WS.txt"
		ELSE
			WRITE(fileout, "(A,I2,A)") "data/vBB", Lsmall, "pp_WS.txt"
		END IF
                        
		OPEN(file_in, FILE=filein, ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to read in ", filein
		ELSE
			READ (file_in, FMT=*, IOSTAT=file_error) vBBpp%b1
		END IF

		OPEN(file_out, FILE=fileout, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to write in ", fileout
		ELSE
			WRITE(file_out, FMT="(E24.16)", IOSTAT=file_error) vBBpp%b1
		END IF
		
		RETURN

	END SUBROUTINE SymVBBpp_OpenFiles
	
	SUBROUTINE SymVBBpp_CloseFiles()

		CLOSE(file_in)
		CLOSE(file_out)		
		
		RETURN

	END SUBROUTINE SymVBBpp_CloseFiles

        !---------------------------------------------------------------------------------------!
	!											!
	!											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymVBBpp_write(vBBpp, Lvalue, Nmax)
		TYPE (SymVBBpp), INTENT(IN) :: vBBpp
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
								
								CALL SymD3Tensor_SymD3Tensor_get(sum1, vBBpp%v1_pair, la, na, nc, lb, nb, nd)
								CALL SymD3Tensor_SymD3Tensor_get(sum2, vBBpp%v2_pair, la, na, nc, lb, nb, nd)
								
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

	END SUBROUTINE SymVBBpp_write
	
	SUBROUTINE SymVBBpp_read(vBBpp, Lvalue, Nmax)
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp
		INTEGER, INTENT(IN) :: Lvalue, Nmax

		DOUBLE PRECISION b1
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
									
								CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v1_pair, la, na, nc, lb, nb, nd, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v1_pair, lb, nb, nd, la, na, nc, sum1)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v2_pair, la, na, nc, lb, nb, nd, sum2)
								CALL SymD3Tensor_SymD3Tensor_assign(vBBpp%v2_pair, lb, nb, nd, la, na, nc, sum2)
								
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO

		RETURN
	END SUBROUTINE SymVBBpp_read

 	SUBROUTINE SymVBBpp_del(vBBpp, Lvalue, Nmax)
		TYPE (SymVBBpp), INTENT(INOUT) :: vBBpp
                INTEGER, INTENT(IN) :: Lvalue, Nmax

		CALL SymD3Tensor_SymD3Tensor_del(vBBpp%v1_pair, Lvalue, Nmax)
		CALL SymD3Tensor_SymD3Tensor_del(vBBpp%v2_pair, Lvalue, Nmax)
		
		RETURN
	END SUBROUTINE SymVBBpp_del

END MODULE symvbb
