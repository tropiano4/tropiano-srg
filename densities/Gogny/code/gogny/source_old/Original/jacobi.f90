MODULE jacobi

CONTAINS

	!*************************************************************
	!* This subroutine computes all eigenvalues and eigenvectors *
	!* of a real symmetric square matrix A(N,N). On output, ele- *
	!* ments of A above the diagonal are destroyed. D(N) returns *
	!* the eigenvalues of matrix A. V(N,N) contains, on output,  *
	!* the eigenvectors of A by columns. THe normalization to    *
	!* unity is made by main program before printing results.    *
	!* NROT returns the number of Jacobi matrix rotations which  *
	!* were required.                                            *
	!* --------------------------------------------------------- *
	!* Ref.:"NUMERICAL RECIPES, Cambridge University Press,1986, *
	!*       chap. 11, pages 346-348".                           *
	!*************************************************************
	SUBROUTINE Jacobi_real8(A, N, D, V, NROT)
		INTEGER, INTENT(IN) :: N
		DOUBLE PRECISION, DIMENSION(N,  N), INTENT(INOUT) :: A, V
		DOUBLE PRECISION, DIMENSION(N), INTENT(INOUT) :: D
		INTEGER, INTENT(INOUT) :: NROT
		DOUBLE PRECISION, POINTER :: B(:), Z(:)
		DOUBLE PRECISION c, g, h, s, sm, t, tau, theta, tresh

		INTEGER ialloc, ip, iq, i, j

		ALLOCATE(B(100), stat = ialloc)
		ALLOCATE(Z(100), stat = ialloc)

		DO ip = 1,  N    !initialize V to identity matrix
			DO iq = 1,  N
				V(ip, iq) = 0.d0
			END DO
			V(ip, ip) = 1.d0
		END DO
	  DO ip = 1,  N
	    B(ip) = A(ip, ip)
	    D(ip) = B(ip)
	    Z(ip) = 0.d0    
	  END DO
	  NROT = 0
	  DO i = 1,  50
	    sm = 0.d0
	    DO ip = 1,  N-1     !sum off-diagonal elements
	      DO iq = ip + 1,  N
	        sm = sm + DABS(A(ip, iq))
	      END DO
	    END DO
	    IF(sm .EQ. 0.d0) RETURN  !normal return
	    IF(i .LT. 4) THEN
	      tresh = 0.2d0 * sm ** 2
	    ELSE
	      tresh = 0.d0
	    END IF
	    DO ip = 1,  N-1
	      DO iq = ip + 1,  N
	        g = 100.d0 * DABS(A(ip, iq))
	! after 4 sweeps,  skip the rotation IF the off-diagonal element is small
	        IF((i .GT. 4) .AND. (DABS(D(ip)) + g .EQ. DABS(D(ip))) &
			.AND. (DABS(D(iq)) + g .EQ. DABS(D(iq)))) THEN
			  A(ip, iq) = 0.d0
	        ELSE IF(DABS(A(ip, iq)) .GT. tresh) THEN
		  h = D(iq)-D(ip)
		  IF ((DABS(h) + g) .EQ. DABS(h)) THEN
		    t = A(ip, iq)/h
	          ELSE
		    theta = 0.5d0*h/A(ip, iq)  
	            t = 1.d0/(DABS(theta) + DSQRT(1.d0 + theta**2))
		    IF (theta .LT. 0.d0) t = -t
	          END IF
		  c = 1.d0/DSQRT(1.d0 + t**2)
		  s = t*c
	          tau = s/(1.d0 + c)
		  h = t*A(ip, iq)
		  Z(ip) = Z(ip)-h
		  Z(iq) = Z(iq) + h
		  D(ip) = D(ip)-h
		  D(iq) = D(iq) + h
		  A(ip, iq) = 0.d0
		  DO j = 1,  ip-1
		    g = A(j, ip)
		    h = A(j, iq)
		    A(j, ip) = g-s*(h + g*tau)
		    A(j, iq) = h + s*(g-h*tau)
	          END DO
		  DO j = ip + 1,  iq-1
		    g = A(ip, j)
		    h = A(j, iq)
		    A(ip, j) = g-s*(h + g*tau)
		    A(j, iq) = h + s*(g-h*tau)
	          END DO		      
		  DO j = iq + 1,  N
		    g = A(ip, j)
		    h = A(iq, j)
		    A(ip, j) = g-s*(h + g*tau)
		    A(iq, j) = h + s*(g-h*tau)
	          END DO		  
		  DO j = 1,  N
		    g = V(j, ip)
		    h = V(j, iq)
		    V(j, ip) = g-s*(h + g*tau)
		    V(j, iq) = h + s*(g-h*tau)
	          END DO		  
	          NROT = NROT + 1
	        END IF !IF ((i.gt.4)...
	      END DO !main iq loop
	    END DO !main ip loop
	    DO ip = 1,  N
	      B(ip) = B(ip) + Z(ip)
	      D(ip) = B(ip)
	      Z(ip) = 0.d0
	    END DO
	  END DO !main i loop
		PAUSE ' 50 iterations !'
		RETURN
	END SUBROUTINE Jacobi_real8

END MODULE jacobi
