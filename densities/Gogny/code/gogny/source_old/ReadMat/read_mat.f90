PROGRAM main

	IMPLICIT NONE
 
	INTEGER :: file_desc = 16
	INTEGER :: file_error, la, na, nc, lb, nb, nd, i, Lmax
	INTEGER, DIMENSION(0:20, 0:16) :: size_mat
	DOUBLE PRECISION :: sumi_sp, sumi_dp, sum1_sp, sum1_dp, sum2_sp, sum2_dp, b
	DOUBLE PRECISION, DIMENSION(0:20, 0:16) :: y_axis

	OPEN(file_desc, FILE="input.txt", ACTION="READ", IOSTAT=file_error)
	
	DO la = 0, 20		
		DO i = 0, 16
			size_mat(la, i) = 0
			y_axis(la, i) = 0.0
		END DO
	END DO
	
	Lmax = 0
	
	READ(file_desc, *, IOSTAT=file_error) b
	
	DO WHILE(file_error .EQ. 0)
	
		READ(file_desc, *, IOSTAT=file_error) la, na, nb, lb, nc, nd, sumi_sp, sumi_dp !, sum1_sp, sum1_dp, sum2_sp, sum2_dp
		IF (file_error .NE. 0) EXIT
	
		IF (la .GE. Lmax) Lmax = la
		
		IF (ABS(sumi_sp) .GT. 1.e-1) THEN
			y_axis(la, 0) = -1.0
			size_mat(la, 0) = size_mat(la, 0) + 1
		ELSE IF (ABS(sumi_sp) .GT. 0.5*1.e-2) THEN
			y_axis(la, 1) = LOG(0.5*1.e-2)/LOG(10.0)
			size_mat(la, 1) = size_mat(la, 1) + 1
		ELSE IF (ABS(sumi_sp) .GT. 1.e-2) THEN
			y_axis(la, 2) = -2.0
			size_mat(la, 2) = size_mat(la, 2) + 1
		ELSE IF (ABS(sumi_sp) .GT. 0.5*1.e-3) THEN
			y_axis(la, 3) = LOG(0.5*1.e-3)/LOG(10.0)
			size_mat(la, 3) = size_mat(la, 3) + 1
		ELSE IF (ABS(sumi_sp) .GT. 1.e-3) THEN
			y_axis(la, 4) = -3.0
			size_mat(la, 4) = size_mat(la, 4) + 1
		ELSE IF (ABS(sumi_sp) .GT. 0.5*1.e-4) THEN
			y_axis(la, 5) = LOG(0.5*1.e-4)/LOG(10.0)
			size_mat(la, 5) = size_mat(la, 5) + 1
		ELSE IF (ABS(sumi_sp) .GT. 1.e-4) THEN
			y_axis(la, 6) = -4.0
			size_mat(la, 6) = size_mat(la, 6) + 1
		ELSE IF (ABS(sumi_sp) .GT. 0.5*1.e-5) THEN
			y_axis(la, 7) = LOG(0.5*1.e-5)/LOG(10.0)
			size_mat(la, 7) = size_mat(la, 7) + 1
		ELSE IF (ABS(sumi_sp) .GT. 1.e-5) THEN
			y_axis(la, 8) = -5.0
			size_mat(la, 8) = size_mat(la, 8) + 1
		ELSE IF (ABS(sumi_sp) .GT. 0.5*1.e-6) THEN
			y_axis(la, 9) = LOG(0.5*1.e-6)/LOG(10.0)
			size_mat(la, 9) = size_mat(la, 9) + 1
		ELSE IF (ABS(sumi_sp) .GT. 1.e-6) THEN
			y_axis(la, 10) = -6.0
			size_mat(la, 10) = size_mat(la, 10) + 1
		ELSE IF (ABS(sumi_sp) .GT. 0.5*1.e-7) THEN
			y_axis(la, 11) = LOG(0.5*1.e-7)/LOG(10.0)
			size_mat(la, 11) = size_mat(la, 11) + 1
		ELSE IF (ABS(sumi_sp) .GT. 1.e-7) THEN
			y_axis(la, 12) = -7.0
			size_mat(la, 12) = size_mat(la, 12) + 1
		ELSE
			y_axis(la, 13) = LOG(0.5*1.e-8)/LOG(10.0)
			size_mat(la, 13) = size_mat(la, 13) + 1
		END IF

	END DO
	
	DO la = 0, Lmax		
		DO i = 0, 13
			WRITE(*,'(i4,f6.2,i12)') la, y_axis(la, i), size_mat(la, i)
		END DO
		write(*,'()')
	END DO
	
	STOP
 
END PROGRAM main

