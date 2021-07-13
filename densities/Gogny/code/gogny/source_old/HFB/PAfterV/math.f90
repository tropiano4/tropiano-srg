!---------------------------------------------------------------------!
!                                                                     !
!          MODULE CONTAINING USEFUL MATHEMATICAL FUNCTIONS            !                
!                                                                     !
!---------------------------------------------------------------------!

 MODULE math

	USE complex_eigen
	USE input
	USE indexx
	USE lgfactor

	IMPLICIT NONE

	REAL(KIND = 16), PARAMETER :: EPS = 1.3877787807814457e-17

	REAL(KIND = 16), PARAMETER :: PI = 3.141592653589793238462643
	REAL(KIND = 16), PARAMETER :: PI_C = 0.12698727186848193957
	DOUBLE PRECISION, PARAMETER :: I_4PI = 0.0795774715459
	REAL(KIND = 16), PARAMETER :: FOUR_PI = 12.5663706143591729539

	DOUBLE PRECISION, PARAMETER :: ALPHA = 0.333333333333
	DOUBLE PRECISION, PARAMETER :: I_SALPHA3 = 0.280565858875 ! This is 1/(alpha + 2)^(3/2)

 CONTAINS

        !---------------------------------------------------------------------!
        !   This routine performs the composite Simpson's rule integration    !
        !   of a function f defined by a table of n equispaced values.        !
        !    								      !
        !                        See: Koonin, Computational Physics, p.9      !
        !    								      !
        !    The parameters are:  					      !
        !     f = Array of values of the function f(x)			      !
        !     n = Number of points x_k			      		      !
        !     h = The uniform spacing between x values: h = x_k+1 - x_k       !
        !    result = Estimate of the integral that is returned to caller.    !
        !---------------------------------------------------------------------!

	SUBROUTINE simps(functi,npoint,step,res)

        DOUBLE PRECISION, INTENT(IN) :: step
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: functi
        DOUBLE PRECISION, INTENT(OUT) :: res
	
        INTEGER :: npanel, npoint, nhalf, nbegin, nend, i

        DOUBLE PRECISION :: x

        ! Check to see if number of panels is even.  Number of panels is n - 1.

	nbegin = 1
	!nbegin = 0
	
	npanel = npoint - nbegin
	nhalf  = npanel/2
	
	res = 0.0
      
	! Number of panels is odd.  Use Simpson's 3/8 rule on first three panels, 1/3 rule on rest of them.

	IF ((npanel-2*nhalf).NE.0) THEN

	     res = 3.0*step*(functi(nbegin) + 3.0*(functi(nbegin+1)+functi(nbegin+2)) + functi(nbegin+3))/8.0
	 
	     IF ((npoint-nbegin).EQ.3) RETURN
	 
	     nbegin=nbegin+3
	 
	END IF

	! Apply 1/3 rule - add in first, second, last values

	res = res + step*(functi(nbegin) + 4.0*functi(nbegin+1) + functi(npoint))/3.0
	nbegin = nbegin+2
      
	IF (nbegin.EQ.npoint) THEN
	    RETURN
	ELSE
      
		x = 0.0
		nend = npoint - 1
		
		DO i = nbegin,nend,2
			x = x + functi(i) + 2.0*functi(i+1) 
		END DO
	  
		res = res + 2.0*step*x/3.0

		RETURN
	  
	END IF

        END SUBROUTINE simps
	
	SUBROUTINE Simpson_Kind16(functi,npoint,step,res)

        REAL(KIND=16), INTENT(IN) :: step
        REAL(KIND=16), DIMENSION(:), INTENT(IN) :: functi
        REAL(KIND=16), INTENT(OUT) :: res
	
        INTEGER :: npanel, npoint, nhalf, nbegin, nend, i

        REAL(KIND=16) :: x

        ! Check to see if number of panels is even.  Number of panels is n - 1.

	nbegin = 1
	!nbegin = 0
	
	npanel = npoint - nbegin
	nhalf  = npanel/2
	
	res = 0.0
      
	! Number of panels is odd.  Use Simpson's 3/8 rule on first three panels, 1/3 rule on rest of them.

	IF ((npanel-2*nhalf).NE.0) THEN

	     res = 3.0*step*(functi(nbegin) + 3.0*(functi(nbegin+1)+functi(nbegin+2)) + functi(nbegin+3))/8.0
	 
	     IF ((npoint-nbegin).EQ.3) RETURN
	 
	     nbegin=nbegin+3
	 
	END IF

	! Apply 1/3 rule - add in first, second, last values

	res = res + step*(functi(nbegin) + 4.0*functi(nbegin+1) + functi(npoint))/3.0
	nbegin = nbegin+2
      
	IF (nbegin.EQ.npoint) THEN
	    RETURN
	ELSE
      
		x = 0.0
		nend = npoint - 1
		
		DO i = nbegin,nend,2
			x = x + functi(i) + 2.0*functi(i+1) 
		END DO
	  
		res = res + 2.0*step*x/3.0

		RETURN
	  
	END IF

        END SUBROUTINE Simpson_Kind16
	
	
	! Funciones de apoyo
	FUNCTION c(N, L)
		REAL(KIND = 16) c
		INTEGER, INTENT(IN) :: N, L

		c = PI_C * PAR(n) * EXP(0.5 * (DDLogFactorials(n) + DDLogSemiFactorials(n + l)))
		RETURN
	END FUNCTION c

	FUNCTION L(a)
		INTEGER L
		INTEGER, INTENT(IN) :: a

		L = (a + 1) / 2
		RETURN
	END FUNCTION L

	FUNCTION J(a)
		INTEGER J
		INTEGER, INTENT(IN) :: a

		J = ((a - L(a)) * 2) + 1
		RETURN
	END FUNCTION J

	FUNCTION LS(a)
		DOUBLE PRECISION LS
		INTEGER, INTENT(IN) :: a

		LS = (0.5 * PAR(a)) / ((a / 2) + 1)
		RETURN
	END FUNCTION LS

	! Maximum number for n for a given l value

	FUNCTION DIM(a)
		INTEGER DIM
		INTEGER, INTENT(IN) :: a

			                		DIM = MIN(Nmax,NmaxOfL(L(a)))
		IF (Basis .EQ. 1 .OR. CompHO .EQ. 1) 	DIM = ((N_0 - L(a)) / 2) + 1
		
		RETURN
	END FUNCTION DIM

        ! Giving the parity of an integer n (+1 or -1)
	FUNCTION PAR(n)
		INTEGER PAR
		INTEGER, INTENT(IN) :: n

		IF (MOD(n, 2) .EQ. 0) THEN
			PAR = 1
		ELSE
			PAR = -1
		END IF
		RETURN
	END FUNCTION PAR

	SUBROUTINE SWAP(a, b)
		INTEGER, INTENT(INOUT) :: a, b

		INTEGER c

		c = a
		a = b
		b = c
		RETURN
	END SUBROUTINE SWAP

	FUNCTION CUAD2(la, lb, k)
		DOUBLE PRECISION CUAD2
		INTEGER, INTENT(IN) :: la, lb, k

		CUAD2 = (la * (la + 1.0)) + (lb * (lb + 1.0)) - (k * (k + 1.0))
		RETURN
	END FUNCTION CUAD2

	FUNCTION MIN_3N(N1, N2, k)
		INTEGER MIN_3N
		INTEGER, INTENT(IN) :: N1, N2, k

		MIN_3N = MAX((ABS(N1 - N2) - k), 0) / 2
		RETURN
	END FUNCTION MIN_3N

	FUNCTION MIN_5N(n1, l1, n2, l2, k)
		INTEGER MIN_5N
		INTEGER, INTENT(IN) :: n1, l1, n2, l2, k

		INTEGER N1_arg, N2_arg

		N1_arg = 2 * n1 + l1
		N2_arg = 2 * n2 + l2
		MIN_5N = MAX((ABS(N1_arg - N2_arg) - k) / 2, 0)
		RETURN
	END FUNCTION MIN_5N

	FUNCTION Char2Int(str)
		INTEGER Char2Int
		CHARACTER(*), INTENT(IN) :: str

		INTEGER idx

		idx = 1
		Char2Int = 0
		DO WHILE (idx .LE. LEN(str))
			SELECT CASE (str(idx:idx))
			CASE ('0')
				Char2Int = Char2Int * 10
			CASE ('1')
				Char2Int = (Char2Int * 10) + 1
			CASE ('2')
				Char2Int = (Char2Int * 10) + 2
			CASE ('3')
				Char2Int = (Char2Int * 10) + 3
			CASE ('4')
				Char2Int = (Char2Int * 10) + 4
			CASE ('5')
				Char2Int = (Char2Int * 10) + 5
			CASE ('6')
				Char2Int = (Char2Int * 10) + 6
			CASE ('7')
				Char2Int = (Char2Int * 10) + 7
			CASE ('8')
				Char2Int = (Char2Int * 10) + 8
			CASE ('9')
				Char2Int = (Char2Int * 10) + 9
			CASE DEFAULT
				RETURN
			END SELECT
			idx = idx + 1
		END DO
		RETURN
	END FUNCTION Char2Int

	SUBROUTINE Int2Char(char_out, num_in)
		CHARACTER(*), INTENT(INOUT) :: char_out
		INTEGER, INTENT(IN) :: num_in

		INTEGER max_len, cur_len, num, i
		CHARACTER digit

		max_len = LEN(char_out)
		cur_len = 0
		char_out = ""
		num = num_in
		DO WHILE ((cur_len .LT. max_len) .AND. (num .GT. 0))
			digit = CHAR(MOD(num, 10) + 48)
			num = num / 10
			cur_len = cur_len + 1
			IF (cur_len .GT. 1) THEN
				DO i = cur_len, 2, -1
					char_out(i:i) = char_out(i - 1:i - 1)
				END DO
			END IF
			char_out(1:1) = digit
		END DO
		char_out(cur_len + 1:cur_len + 1) = CHAR(0)
		RETURN
	END SUBROUTINE Int2Char

!        EIGENVALUES AND EIGENVECTORS OF COMPLEX MATRICES

	SUBROUTINE MatDiagCmplx(Matrix, Eigenvectors, Eigenvalues, N)
	
	COMPLEX, DIMENSION(N), INTENT(OUT) :: Eigenvalues
	COMPLEX, DIMENSION(N, N), INTENT(OUT) :: Eigenvectors
	COMPLEX, DIMENSION(N, N), INTENT(IN) :: Matrix
	INTEGER, INTENT(IN) :: N
	
	LOGICAL :: ibal
	
	INTEGER :: ierr, i, j
	INTEGER, DIMENSION(:), ALLOCATABLE :: indx
	
	DOUBLE PRECISION, DIMENSION(N, N) :: ar, ai, zr, zi
	DOUBLE PRECISION, DIMENSION(N) :: wr, wi
	
	ar = REAL(Matrix)
	ai = AIMAG(Matrix)
	
	CALL dceigv(ibal, ar, ai, n, wr, wi, zr, zi, ierr)

	! We sort the eigenvalues (which are complex) by 
	! ascending order of their real part.

	ALLOCATE(indx(N))

	DO i = 1, N
		indx(i) = i
	END DO

	CALL indexx_real8(N, wr, indx)
		
	DO j = 1, N
		Eigenvalues(j) = wr(indx(j)) + CMPLX(0,1)*wi(indx(j))
		!Eigenvalues(j) = wr(j) + CMPLX(0,1)*wi(j)
		DO i = 1, N
			Eigenvectors(i,j) = zr(i,indx(j)) + CMPLX(0,1)*zi(i,indx(j))
			!Eigenvectors(i,j) = zr(i,j) + CMPLX(0,1)*zi(i,j)
		END DO
	END DO
	
	DEALLOCATE(indx)

	RETURN
	END SUBROUTINE MatDiagCmplx

	! Complex Matrix Inversion routine
	!
	SUBROUTINE MatInvCmplx(N, M, AA, Inverse, DET)
	
	COMPLEX, DIMENSION(N,N), INTENT(OUT) :: Inverse
	DOUBLE PRECISION, INTENT(OUT) :: DET
	COMPLEX, DIMENSION(N,N), INTENT(IN) :: AA
	INTEGER, INTENT(IN) :: N, M

	DOUBLE PRECISION, DIMENSION(2*N,2*N) :: MATRIX
	DOUBLE PRECISION, DIMENSION(2*N,2*N) :: InverseReal
	
	INTEGER :: i, j
      
	DO i = 1, N
	
		DO j=1, N
			MATRIX(i,j) = REAL(AA(i,j))
		END DO
		
		DO j = N+1, 2*N
			MATRIX(i,j) = AIMAG(AA(i,j-N))
		END DO
	
	END DO
      
	DO i = N+1, 2*N
	
		DO j=1, N
			MATRIX(i,j) = -AIMAG(AA(i-N,j))
		END DO
		
		DO j = N+1, 2*N
			MATRIX(i,j) = REAL(AA(i-N,j-N))
		END DO
	
	END DO
      
	CALL MatInvReal(2*N,0,MATRIX,InverseReal,DET)    

	DO i = 1, N
		DO j=1, N
			Inverse(i,j) = InverseReal(i,j) + CMPLX(0,1)*InverseReal(i,j+N)
		END DO
	END DO
      
	RETURN
	END SUBROUTINE MatInvCmplx

!*******************************************
!*  SOLVING A LINEAR MATRIX SYSTEM AX = B  *
!*  with Gauss-Jordan method using full    *
!*  pivoting at each step. During the pro- *
!* cess, original A and B matrices are     *
!* destroyed to spare storage location.    *
!* --------------------------------------- *
!* INPUTS:    A   MATRIX N*N               *
!*            B   MATRIX N*M               *
!* --------------------------------------- *
!* OUTPUS:    A   INVERSE OF A N*N         *
!*            DET  DETERMINANT OF A        *
!*            B   SOLUTION MATRIX N*M      *
!* --------------------------------------- *
!* NOTA - If M=0 inversion of A matrix     *
!*        only.                            *
!*******************************************

      SUBROUTINE MatInvReal(N,M,AA,Inverse,DET)    

      DOUBLE PRECISION, DIMENSION(N,N), INTENT(OUT) :: Inverse
      DOUBLE PRECISION, INTENT(OUT) :: DET
      DOUBLE PRECISION, DIMENSION(N,N), INTENT(IN) :: AA
      INTEGER, INTENT(IN) :: N, M

      DOUBLE PRECISION, DIMENSION(N,N) :: MATRIX
      DOUBLE PRECISION, DIMENSION(N,M) :: BB
      DOUBLE PRECISION, POINTER ::  PC(:),PL(:),CS(:)

      DOUBLE PRECISION :: PV, PAV, TT, EPSMACH
      
      INTEGER :: ialloc, I, J, K, IK, JK

!Initializations :                       
	EPSMACH = 1.E-20

	DO i = 1, N
		DO j = 1, N
			MATRIX(I,J) = AA(i,j)
		END DO
	END DO
	
      allocate(PC(1:N),stat=ialloc)
      allocate(PL(1:N),stat=ialloc)
      allocate(CS(1:N),stat=ialloc)   

      DET=1.D0

      DO I=1,N                            

        PC(I)=0.D0

        PL(I)=0.D0

        CS(I)=0.D0

      END DO

!main loop                                

      DO K=1,N                            

!Searching greatest pivot :               

        PV=MATRIX(K,K)                        

        IK=K                              

        JK=K                              

        PAV=DABS(PV)                      

        DO I=K,N                          

          DO J=K,N                        

            IF (DABS(MATRIX(I,J)).GT.PAV) THEN

              PV=MATRIX(I,J)                  

              PAV=DABS(PV)                

              IK=I                        

              JK=J                        

            ENDIF                         

          ENDDO                           

        ENDDO                             

!Search terminated, the pivot is in location I=IK, J=JK.

!Memorizing pivot location: :                  

        PC(K)=JK                               

        PL(K)=IK                               

!Determinant DET is actualised

!If DET=0, ERROR MESSAGE and STOP

!Machine dependant EPSMACH equals here 1.DE-20 

        IF (IK.NE.K) DET=-DET                  

        IF (JK.NE.K) DET=-DET                  

        DET=DET*PV                             

        IF (DABS(DET).LT.EPSMACH) THEN         

!Error message and Stop                        

          PRINT 10                             

          STOP                                 

        ENDIF                                  

!POSITIONNING PIVOT IN K,K:                    

        IF(IK.NE.K) THEN                       

          DO I=1,N                             

!EXCHANGE LINES IK and K:                      

            TT=MATRIX(IK,I)                        

            MATRIX(IK,I)=MATRIX(K,I)                   

            MATRIX(K,I)=TT                         

          ENDDO                                

        ENDIF                                  

      IF (M.NE.0) THEN                         

        DO I=1,M                               

          TT=BB(IK,I)                          

          BB(IK,I)=BB(K,I)                     

          BB(K,I)=TT                           

        ENDDO                                  

      ENDIF                                    

!Pivot is at correct line                      

        IF(JK.NE.K) THEN                       

          DO I=1,N                             

!Exchange columns JK and K of matrix MATRIX        

            TT=MATRIX(I,JK)                        

            MATRIX(I,JK)=MATRIX(I,K)                   

            MATRIX(I,K)=TT                         

          ENDDO                                

        ENDIF                                  

!Pivot is at correct column and located in K,K 


!Store column K in vector CS                   

!then set column K to zero                     

        DO I=1,N                               

          CS(I)=MATRIX(I,K)                        

          MATRIX(I,K)=0.D0                         

        ENDDO                                  

        CS(K)=0.                               

        MATRIX(K,K)=1.                             

!Modify line K :                               

        IF(DABS(PV).LT.EPSMACH) THEN           

          WRITE(*,*) '  PIVOT TOO SMALL - STOP'

          STOP                                 

        ENDIF                                  

        DO I=1,N                               

          MATRIX(K,I)=MATRIX(K,I)/PV                   

        ENDDO                                  

        IF (M.NE.0) THEN                       

          DO I=1,M                             

            BB(K,I)=BB(K,I)/PV                 

          ENDDO                                

        ENDIF                                  

!Modify other lines of matrix MATRIX:              

        DO J=1,N                               

          IF (J.EQ.K) CONTINUE                 

          DO I=1,N                             

!Modify line J of matrix MATRIX :                  

            MATRIX(J,I)=MATRIX(J,I)-CS(J)*MATRIX(K,I)      

          ENDDO                                

          IF (M.NE.0) THEN                     

            DO I=1,M                           

              BB(J,I)=BB(J,I)-CS(J)*BB(K,I)    

            ENDDO                              

          ENDIF                                

        ENDDO                                  

!Line K is ready.                              

      ENDDO                                    

!End of K loop                                 

!The matrix MATRIX is inverted - Rearrange MATRIX      

!Exchange lines                                

      DO I=N,1,-1                              

        IK=PC(I)                               

        IF (IK.EQ.I) CONTINUE                  

!EXCHANGE LINES I AND PC(I) OF MATRIX:             

        DO J=1,N                               

          TT=MATRIX(I,J)                           

          MATRIX(I,J)=MATRIX(IK,J)                     

          MATRIX(IK,J)=TT                          

        ENDDO                                  

        IF (M.NE.0) THEN                       

          DO J=1,M                             

            TT=BB(I,J)                         

            BB(I,J)=BB(IK,J)                   

            BB(IK,J)=TT                        

          ENDDO                                

        ENDIF                                  

!NO MORE EXCHANGE NEEDED                       

!GO TO NEXT LINE                               

      ENDDO                                    

!EXCHANGE COLUMNS                              

      DO J=N,1,-1                              

        JK=PL(J)                               

        IF (JK.EQ.J) CONTINUE                  

!EXCHANGE COLUMNS J AND PL(J) OF MATRIX :          

        DO I=1,N                               

          TT=MATRIX(I,J)                           

          MATRIX(I,J)=MATRIX(I,JK)                     

          MATRIX(I,JK)=TT                          

        ENDDO                                  

!NO MORE EXCHANGE NEEDED                       

!GO TO NEXT COLUMN   

      ENDDO                                    

!REARRANGEMENT TERMINATED.

	DO i = 1, N
		DO j = 1, N
			Inverse(I,J) = MATRIX(i,j)
		END DO
	END DO
	


      RETURN                                   

   10 FORMAT(///'  DETERMINANT EQUALS ZERO, NO SOLUTION!')                    

      END SUBROUTINE MatInvReal

END MODULE math
