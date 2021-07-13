 MODULE eigenval

 CONTAINS

	! A = Matriz simetrica
	! D = Matriz diagonal
	! Z = Matriz
	SUBROUTINE EigenValues(NM, N, A, D, Z)
		INTEGER, INTENT(IN) :: NM, N
		DOUBLE PRECISION, DIMENSION(:, :), INTENT(IN) :: A
		DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: D
		DOUBLE PRECISION, DIMENSION(:, :), INTENT(INOUT) :: Z

		! E = Matriz diagonal
		DOUBLE PRECISION, DIMENSION(NM, N) :: E
		INTEGER IERR

		CALL TRED2(NM, N, A, D, E, Z)
		
		CALL TQL2(NM, N, D, E, Z, IERR)
		
		IF (IERR .NE. 0) STOP "ERROR: EigenValues"
		RETURN
	END SUBROUTINE EigenValues

      SUBROUTINE TQL2(NM,N,D,E,Z,IER)
!-------------------------------------------------------------------------
!     QL METHOD TO DETERMINE THE EIGENVALUES AND EIGENVECTORS OF:
!
!       1)  A SYMMETRIC TRIDIAGONAL MATRIX.
!       2)  A FULL SYMMETRIC MATRIX AFTER A PREVIOUS CALL TO TRED2.
!
!     CALLING MODE:
!               CALL TQL2(NM,N,D,E,Z,IER)
!     INPUTSS:
!     NM  (I4)  1ST DIMENSION OF MATRICES A AND Z IN CALLING PROGRAM
!     N   (I4)  SIZE OF Z
!     D  (R*8)  MAIN DIAGONAL (N) OF THE TRIDIAGONAL MATRIX
!     E  (R*8)  SUB-DIAGONAL (N) OF THE TRIDIAGONAL MATRIX
!     Z  (R*8)  TABLE (NM,N) STORING THE UNITY MATRIX IF THE TRIDIAGONAL
!               MATRIX IS DEFINED BY D AND E, CASE #1.
!               FOR CASE #2, IT CONTAINS THE ELEMENTS OF THE TRANSFORMATION
!               MATRIX AFTER A CALL TO TRED2.
!     OUTPUTS:
!     D  (R*8)  EIGENVALUES
!     Z  (R*8)  EIGENVECTORS
!     IER (I4)  ERROR CODE = 0,  CONVERGENCE OK.
!                          = L,  NO CONVERGENCE FOR THE Lth EIGENVALUE
!
!     REFERENCE:
!     J.H.WILKINSON,-C.REINSCH,R.S.MARTIN
!     HANDBOOK FOR AUTOMATIC COMPUTATION, VOL.2, LINEAR ALGEBRA
!     SPRINGER-VERLAG 1971.
!-------------------------------------------------------------------------
      INTEGER I,J,K,L,M,N,NM,JM,IER
      !REAL *8 D(N),E(N),Z(NM,N),B,C,F,G,H,P,R,S,EPS,EPS1
      DOUBLE PRECISION :: D(N),E(N),Z(NM,N),B,C,F,G,H,P,R,S,EPS,EPS1
      DATA EPS /0.D0/,JM /30/
      IER = 0
      IF (N.EQ.1) GO TO 38
!
!     MACHINE EPSILON
!
      IF (EPS.NE.0.D0) GO TO 12
      EPS = 1.D0
   10 EPS = EPS/2.D0
      EPS1 = 1.D0+EPS
      IF (EPS1.GT.1.D0) GO TO 10
!
   12 DO 14 I = 2,N
   14 E(I-1) = E(I)
      E(N) = 0.D0
      F = 0.D0
      B = 0.D0
!
      DO 28 L = 1,N
      J = 0
      H = EPS*(ABS(D(L))+ABS(E(L)))
      IF (B.LT.H) B = H
!
!     SEEK SMALLEST ELEMENT OF SUBDIAGONAL
!
      DO 16 M = L,N
      IF (ABS(E(M)).LE.B) GO TO 18
   16 CONTINUE
   18 IF (M.EQ.L) GO TO 26

!     START ITERATION

   20 IF (J.EQ.JM) GO TO 36
      J = J+1

!     SHIFT

      G = D(L)
      P = (D(L+1)-G)/(2.D0*E(L))
      R = SQRT(P*P+1.D0)
      D(L) = E(L)/(P+SIGN(R,P))
      H = G-D(L)
      DO 22 I = L+1,N
   22 D(I) = D(I)-H
      F = F+H

!     QL TRANSFORMATION

      P = D(M)
      C = 1.D0
      S = 0.D0
      DO 24 I = M-1,L,-1
      G = C*E(I)
      H = C*P
      IF (ABS(P).GE.ABS(E(I))) THEN
      C = E(I)/P
      R = SQRT(C*C+1.D0)
      E(I+1) = S*P*R
      S = C/R
      C = 1.D0/R
      ELSE
      C = P/E(I)
      R = SQRT(C*C+1.D0)
      E(I+1) = S*E(I)*R
      S = 1.D0/R
      C = C*S
      ENDIF
      P = C*D(I)-S*G
      D(I+1) = H+S*(C*G+S*D(I))

!     ELEMENTS OF EIGENVECTORS

      DO 24 K = 1,N
      H = Z(K,I+1)
      Z(K,I+1) = S*Z(K,I)+C*H
      Z(K,I) = Z(K,I)*C-S*H
   24 CONTINUE
      E(L) = S*P
      D(L) = C*P
      IF (ABS(E(L)).GT.B) GO TO 20

!     CONVERGENCE

   26 D(L) = D(L)+F
   28 CONTINUE

!     SORT EIGENVALUES AND EIGENVECTORS
!     IN ASVENDING ORDER

      DO 34 L = 2,N
      I = L-1
      K = I
      P = D(I)
      DO 30 J = L,N
      IF (D(J).GE.P) GO TO 30
      K = J
      P = D(J)
   30 CONTINUE
      IF (K.EQ.I) GO TO 34
      D(K) = D(I)
      D(I) = P
      DO 32 J = 1,N
      P = Z(J,I)
      Z(J,I) = Z(J,K)
   32 Z(J,K) = P
   34 CONTINUE
      GO TO 38

!     NO CONVERGENCE

   36 IER = L
   38 RETURN
      END SUBROUTINE TQL2

      SUBROUTINE TRED2(NM,N,A,D,E,Z)
!---------------------------------------------------------------------------
!     TRIDIAGONALIZATION OF A SYMMETRIC MATRIX BY ORTHOGONAL TRANSFORMATIONS
!     (ALGORITHM OF HOUSEHOLDER)
!     CALLING MODE:
!               CALL TRED2(NM,N,A,D,E,Z)
!     INPUTS:
!     NM  (I4)  1ST DIMENSION OF MATRICES A AND Z IN CALLING PROGRAM
!     N   (I4)  SIZE OF A
!     A  (R*8)  TABLE(NM,N) STORING THE COEFFICIENTS OF SYMMETRIC A MATRIX
!               (LOWER HALF), A IS NOT DESTROYED DURING THE PROCESS
!               IF Z MATRIX HAS NOT THE SAME ADDRESS.
!     OUTPUTS:
!     D  (R*8)  MAIN DIAGONAL (N) OF REDUCED TRIDIAGONAL MATRIX
!     E  (R*8)  SUB-DIAGONAL (N) OF REDUCED TRIDIAGONAL MATRIX
!     Z  (R*8)  TABLE (NM,N) STORING THE ELEMENTS OF THE ORTHOGONAL 
!               TRANSFORMATION MATRIX.
!     REFERENCE:
!     J.H.WILKINSON,-C.REINSCH,R.S.MARTIN
!     HANDBOOK FOR AUTOMATIC COMPUTATION, VOL.2, LINEAR ALGEBRA
!     SPRINGER-VERLAG 1971.
!-----------------------------------------------------------------------
      INTEGER I,J,K,L,N,NM
      !REAL *8 A(NM,N),D(N),E(N),Z(NM,N),F,G,H,HH,SCALE
      DOUBLE PRECISION :: A(NM,N),D(N),E(N),Z(NM,N),F,G,H,HH,SCALE

!     LOWER HALF OF A PUT INTO Z

      DO 10 I = 1,N
      DO 10 J = 1,I
   10 Z(I,J) = A(I,J)
      IF (N.EQ.1) GO TO 32

!     N-2 STAGE OF TRANSFORMATION

      DO 30 I = N,2,-1
      L = I-1
      H = 0.

!     CONDITIONNING BY NORM OF A

      SCALE = 0.
      IF (L.LT.2) GO TO 14
      DO 12 K = 1,L
   12 SCALE = SCALE+ABS(Z(I,K))
      IF (SCALE.NE.0.) GO TO 16

   14 E(I) = Z(I,L)
      GO TO 28

   16 DO 18 K = 1,L
      Z(I,K) = Z(I,K)/SCALE
      H = H+Z(I,K)*Z(I,K)
   18 CONTINUE

      F = Z(I,L)
      G = -SIGN(SQRT(H),F)
      E(I) = SCALE*G
      H = H-F*G
      Z(I,L) = F-G
      F = 0.
      DO 24 J = 1,L
      Z(J,I) = Z(I,J)/H
      G = 0.

!     ELEMENT OF A*U
      DO 20 K = 1,J
   20 G = G+Z(J,K)*Z(I,K)
      IF (L.GE.J+1) THEN
      DO 22 K = J+1,L
   22 G = G+Z(K,J)*Z(I,K)

!     ELEMENT OF P = A*U/H

      END IF
      E(J) = G/H
      F = F+E(J)*Z(I,J)
   24 CONTINUE

!     ELEMENT OF K

      HH = F/(H+H)

!     REDUCED FORM OF A

      DO 26 J = 1,L
      F = Z(I,J)
      G = E(J)-HH*F
      E(J) = G
      DO 26 K = 1,J
      Z(J,K) = Z(J,K)-F*E(K)-G*Z(I,K)
   26 CONTINUE
!
   28 D(I) = H
   30 CONTINUE

!     END OF TRANSFORMATION

   32 D(1) = 0.
      E(1) = 0.

!     ACCUMULATE TRANSFORMATION MATRICES IN Z

      DO 40 I = 1,N
      L = I-1
      IF (D(I).NE.0.) THEN
      DO 36 J = 1,L
      G = 0.
      DO 34 K = 1,L
   34 G = G+Z(I,K)*Z(K,J)
      DO 36 K = 1,L
      Z(K,J) = Z(K,J)-G*Z(K,I)
   36 CONTINUE
      END IF
      D(I) = Z(I,I)
      Z(I,I) = 1.
      IF (L.LT.1) GO TO 40
      DO 38 J = 1,L
      Z(I,J) = 0.
      Z(J,I) = 0.
   38 CONTINUE
   40 CONTINUE

      RETURN
      END SUBROUTINE TRED2

END MODULE eigenval
