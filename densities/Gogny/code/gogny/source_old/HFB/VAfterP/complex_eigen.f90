!                     EIGENVALUES AND EIGENVECTORS
!                 OF DOUBLE PRECISION COMPLEX MATRICES

! Adapted from the NSWC Mathematics Library, based upon Eispack routines.
! This version by Alan Miller
! Alan.Miller @ vic.cmis.csiro.au
! http://www.ozemail.com.au/~milleraj

! Latest revision - 26 August 1998

 MODULE complex_eigen
	
	USE constant

	IMPLICIT NONE

 CONTAINS

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	FUNCTION dcpabs(x, y) RESULT(fn_val)

	REAL (dp), INTENT(IN) :: x
	REAL (dp), INTENT(IN) :: y
	REAL (dp)             :: fn_val

	!     --------------------------------------
	!     EVALUATION OF SQRT(X*X + Y*Y)
	!     --------------------------------------
	REAL (dp) :: a

	IF (ABS(x) > ABS(y)) THEN
 		a = y / x
		fn_val = ABS(x)*SQRT(1.0_dp + a*a)
  		RETURN
	END IF

	IF (y == 0.0_dp) THEN
		fn_val = 0.0_dp
		RETURN
	END IF

	a = x / y
	fn_val = ABS(y)*SQRT(1.0_dp + a*a)
	RETURN

	END FUNCTION dcpabs

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE dcsqrt (z, w)

	REAL (dp), INTENT(IN)   :: z(:)
	REAL (dp), INTENT(OUT)  :: w(:)

	! ----------------------------------------------------------------------
	!          W = SQRT(Z) FOR THE COMPLEX NUMBER Z

	!                       -----------

	!   Z AND W ARE INTERPRETED AS COMPLEX NUMBERS.
	!   IT IS ASSUMED THAT Z(1) AND Z(2) ARE THE REAL AND IMAGINARY PARTS OF
	!   THE COMPLEX NUMBER Z, AND THAT W(1) AND W(2) ARE THE REAL AND IMAGINARY
	!   PARTS OF W.
	! ----------------------------------------------------------------------
	REAL (dp)            :: x, y, r
	REAL (dp), PARAMETER :: half = 0.5_dp

	x = z(1)
	y = z(2)
	IF (x < 0.0_dp) THEN
		GO TO 30
	ELSE IF (x > 0.0_dp) THEN
		GO TO 20
	END IF

	IF (y /= 0.0_dp) GO TO 11
	w(1) = 0.0_dp
	w(2) = 0.0_dp
	RETURN

 11	 w(1) = SQRT(half * ABS(y))
	w(2) = SIGN(w(1), y)
	RETURN

 20 	IF (y /= 0.0_dp) GO TO 21
	w(1) = SQRT(x)
	w(2) = 0.0_dp
	RETURN

 21	r = dcpabs(x,y)
	w(1) = SQRT(half * (r + x))
	w(2) = half * y / w(1)
	RETURN

 30	IF (y /= 0.0_dp) GO TO 31
	w(1) = 0.0_dp
	w(2) = SQRT(ABS(x))
	RETURN

 31 	r = dcpabs(x, y)
	w(2) = SQRT(half * (r - x))
	w(2) = SIGN(w(2), y)
	w(1) = half * y / w(2)
	RETURN
	END SUBROUTINE dcsqrt

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE dceigv (ibal, ar, ai, n, wr, wi, zr, zi, ierr)
	!-----------------------------------------------------------------------
	!           EIGENVALUES AND EIGENVECTORS OF COMPLEX MATRICES
	!
	!  ibal = .TRUE. if balancing is required
	!  ar, ai = the real and imaginary parts of the matrix (destroyed)
	!  wr, wi = the real and imaginary parts of the eigenvalues
	!  zr, zi = the real and imaginary parts of the eigenvectors
	!
	!-----------------------------------------------------------------------

	LOGICAL, INTENT(IN)        :: ibal
	REAL (dp), INTENT(IN OUT)  :: ar(:,:)     ! ar(ka,n)
	REAL (dp), INTENT(IN OUT)  :: ai(:,:)     ! ai(ka,n)
	INTEGER, INTENT(IN)        :: n
	REAL (dp), INTENT(OUT)     :: wr(:)
	REAL (dp), INTENT(OUT)     :: wi(:)
	REAL (dp), INTENT(OUT)     :: zr(:,:)     ! zr(ka,n)
	REAL (dp), INTENT(OUT)     :: zi(:,:)     ! zi(ka,n)
	INTEGER, INTENT(OUT)       :: ierr

	! Local variables
	REAL (dp) :: ortr(n), orti(n), scale(n)
	INTEGER   :: high, low

	low = 1
	high = n
	IF (ibal) CALL dcbal(n, ar, ai, low, high, scale)
	CALL dcorth(n, low, high, ar, ai, ortr, orti)
	CALL dcmqr2(n, low, high, ortr, orti, ar, ai, wr, wi, zr, zi, ierr)
	IF (ierr /= 0) RETURN
	IF (ibal) CALL dcbabk(n, low, high, scale, n, zr, zi)
	
	RETURN
	END SUBROUTINE dceigv



SUBROUTINE dcbal(n, ar, ai, low, igh, scale)

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: ar(:,:)      ! ar(nm,n)
REAL (dp), INTENT(IN OUT)  :: ai(:,:)      ! ai(nm,n)
INTEGER, INTENT(OUT)       :: low
INTEGER, INTENT(OUT)       :: igh
REAL (dp), INTENT(OUT)     :: scale(:)

!-----------------------------------------------------------------------

!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE
!     CBALANCE, WHICH IS A COMPLEX VERSION OF BALANCE,
!     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).

!     DCBAL BALANCES A COMPLEX MATRIX AND
!     ISOLATES EIGENVALUES WHENEVER POSSIBLE.

!     ON INPUT-

!        N IS THE ORDER OF THE MATRIX,

!        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
!          RESPECTIVELY, OF THE COMPLEX MATRIX TO BE BALANCED.

!     ON OUTPUT-

!        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
!          RESPECTIVELY, OF THE BALANCED MATRIX,

!        LOW AND IGH ARE TWO INTEGERS SUCH THAT AR(I,J) AND AI(I,J)
!          ARE EQUAL TO ZERO IF
!           (1) I IS GREATER THAN J AND
!           (2) J=1,...,LOW-1 OR I=IGH+1,...,N,

!        SCALE CONTAINS INFORMATION DETERMINING THE
!           PERMUTATIONS AND SCALING FACTORS USED.

!     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH
!     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED
!     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS
!     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN
!        SCALE(J) = P(J),    FOR J = 1,...,LOW-1
!                 = D(J,J)       J = LOW,...,IGH
!                 = P(J)         J = IGH+1,...,N.
!     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,
!     THEN 1 TO LOW-1.

!     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.

!     THE ALGOL PROCEDURE EXC CONTAINED IN CBALANCE APPEARS IN DCBAL IN LINE.
!     (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS K, L HAVE BEEN REVERSED.)

!-----------------------------------------------------------------------
INTEGER   :: i, j, k, l, m, iexc
REAL (dp) :: c, f, g, r, s, b2, radix
LOGICAL   :: noconv

!     ********** RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!                THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION.

radix = ipmpar(4)

!                **********

b2 = radix * radix
k = 1
l = n
GO TO 100
!     ********** IN-LINE PROCEDURE FOR ROW AND COLUMN EXCHANGE **********
20 scale(m) = j
IF (j == m) GO TO 50

DO i = 1, l
  f = ar(i,j)
  ar(i,j) = ar(i,m)
  ar(i,m) = f
  f = ai(i,j)
  ai(i,j) = ai(i,m)
  ai(i,m) = f
END DO

DO i = k, n
  f = ar(j,i)
  ar(j,i) = ar(m,i)
  ar(m,i) = f
  f = ai(j,i)
  ai(j,i) = ai(m,i)
  ai(m,i) = f
END DO

50 IF (iexc == 2) GO TO 130
!  ******** SEARCH FOR ROWS ISOLATING AN EIGENVALUE AND PUSH THEM DOWN *******
IF (l == 1) GO TO 280
l = l - 1
!     ********** FOR J=L STEP -1 UNTIL 1 DO -- **********
100 loop120: DO j = l, 1, -1
  
  DO i = 1, l
    IF (i == j) CYCLE
    IF (ar(j,i) /= 0.0_dp .OR. ai(j,i) /= 0.0_dp) CYCLE loop120
  END DO
  
  m = l
  iexc = 1
  GO TO 20
END DO loop120

GO TO 140
!     ********** SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
!                AND PUSH THEM LEFT **********
130 k = k + 1

140 loop170: DO j = k, l
  
  DO i = k, l
    IF (i == j) CYCLE
    IF (ar(i,j) /= 0.0_dp .OR. ai(i,j) /= 0.0_dp) CYCLE loop170
  END DO
  
  m = k
  iexc = 2
  GO TO 20
END DO loop170
!     ********** NOW BALANCE THE SUBMATRIX IN ROWS K TO L **********
scale(k:l) = 1.0_dp
!     ********** ITERATIVE LOOP FOR NORM REDUCTION **********
190 noconv = .false.

DO i = k, l
  c = 0.0_dp
  r = 0.0_dp
  
  DO j = k, l
    IF (j == i) CYCLE
    c = c + ABS(ar(j,i)) + ABS(ai(j,i))
    r = r + ABS(ar(i,j)) + ABS(ai(i,j))
  END DO
!     ********** GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW **********
  IF (c == 0.0_dp .OR. r == 0.0_dp) CYCLE
  g = r / radix
  f = 1.0_dp
  s = c + r
  210 IF (c >= g) GO TO 220
  f = f * radix
  c = c * b2
  GO TO 210
  220 g = r * radix
  230 IF (c < g) GO TO 240
  f = f / radix
  c = c / b2
  GO TO 230
!     ********** NOW BALANCE **********
  240 IF ((c + r) / f >= 0.95D0 * s) CYCLE
  g = 1.0_dp / f
  scale(i) = scale(i) * f
  noconv = .true.
  
  DO j = k, n
    ar(i,j) = ar(i,j) * g
    ai(i,j) = ai(i,j) * g
  END DO
  
  DO j = 1, l
    ar(j,i) = ar(j,i) * f
    ai(j,i) = ai(j,i) * f
  END DO
  
END DO

IF (noconv) GO TO 190

280 low = k
igh = l
RETURN
!     ********** LAST CARD OF DCBAL **********
END SUBROUTINE dcbal



SUBROUTINE dcorth(n, low, igh, ar, ai, ortr, orti)

INTEGER, INTENT(IN)       :: n
INTEGER, INTENT(IN)       :: low
INTEGER, INTENT(IN)       :: igh
REAL (dp), INTENT(IN OUT) :: ar(:,:)      ! ar(nm,n)
REAL (dp), INTENT(IN OUT) :: ai(:,:)      ! ai(nm,n)
REAL (dp), INTENT(OUT)    :: ortr(:)      ! ortr(igh)
REAL (dp), INTENT(OUT)    :: orti(:)      ! orti(igh)
!-----------------------------------------------------------------------

!     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF THE ALGOL
!     PROCEDURE ORTHES, NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).

!     GIVEN A COMPLEX MATRIX, DCORTH REDUCES A SUBMATRIX SITUATED IN ROWS
!     AND COLUMNS LOW THROUGH IGH TO UPPER HESSENBERG FORM BY UNITARY
!     SIMILARITY TRANSFORMATIONS.

!     ON INPUT-

!        N IS THE ORDER OF THE MATRIX,

!        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING SUBROUTINE DCBAL.
!          IF DCBAL HAS NOT BEEN USED, SET LOW=1, IGH=N,

!        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
!          RESPECTIVELY, OF THE COMPLEX INPUT MATRIX.

!     ON OUTPUT-

!        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY, OF THE
!          HESSENBERG MATRIX.  INFORMATION ABOUT THE UNITARY TRANSFORMATIONS
!          USED IN THE REDUCTION IS STORED IN THE REMAINING TRIANGLES UNDER
!          THE HESSENBERG MATRIX,

!        ORTR AND ORTI CONTAIN FURTHER INFORMATION ABOUT THE
!          TRANSFORMATIONS.  ONLY ELEMENTS LOW THROUGH IGH ARE USED.

!-----------------------------------------------------------------------

INTEGER   :: i, j, m, la, kp1
REAL (dp) :: f, g, h, fi, fr, scale

la = igh - 1
kp1 = low + 1
IF (la < kp1) GO TO 200

DO m = kp1, la
  h = 0.0_dp
  ortr(m) = 0.0_dp
  orti(m) = 0.0_dp
  scale = 0.0_dp
!     ********** SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) **********
  DO i = m, igh
    scale = scale + ABS(ar(i,m-1)) + ABS(ai(i,m-1))
  END DO
  
  IF (scale == 0.0_dp) CYCLE
!     ********** FOR I=IGH STEP -1 UNTIL M DO -- **********
  DO i = igh, m, -1
    ortr(i) = ar(i,m-1) / scale
    orti(i) = ai(i,m-1) / scale
    h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
  END DO
  
  g = SQRT(h)
  f = dcpabs(ortr(m), orti(m))
  IF (f == 0.0_dp) GO TO 103
  h = h + f * g
  g = g / f
  ortr(m) = (1.0_dp + g) * ortr(m)
  orti(m) = (1.0_dp + g) * orti(m)
  GO TO 105
  
  103 ortr(m) = g
  ar(m,m-1) = scale
!     ********** FORM (I-(U*UT)/H) * A **********
  105 DO j = m, n
    fr = 0.0_dp
    fi = 0.0_dp
!     ********** FOR I=IGH STEP -1 UNTIL M DO -- **********
    DO i = igh, m, -1
      fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
      fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
    END DO
    
    fr = fr / h
    fi = fi / h
    
    DO i = m, igh
      ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
      ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
    END DO
    
  END DO
!     ********** FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) **********
  DO i = 1, igh
    fr = 0.0_dp
    fi = 0.0_dp
!     ********** FOR J=IGH STEP -1 UNTIL M DO -- **********
    DO j = igh, m, -1
      fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
      fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
    END DO
    
    fr = fr / h
    fi = fi / h
    
    DO j = m, igh
      ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
      ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
    END DO
    
  END DO
  
  ortr(m) = scale * ortr(m)
  orti(m) = scale * orti(m)
  ar(m,m-1) = -g * ar(m,m-1)
  ai(m,m-1) = -g * ai(m,m-1)
END DO

200 RETURN
!     ********** LAST CARD OF DCORTH **********
END SUBROUTINE dcorth



SUBROUTINE dcbabk (n, low, igh, scale, m, zr, zi)

INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: low
INTEGER, INTENT(IN)     :: igh
REAL (dp), INTENT(IN)   :: scale(:)
INTEGER, INTENT(IN)     :: m
REAL (dp), INTENT(OUT)  :: zr(:,:)      ! zr(nm,m)
REAL (dp), INTENT(OUT)  :: zi(:,:)      ! zi(nm,m)

!-----------------------------------------------------------------------

!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE
!     CBABK2, WHICH IS A COMPLEX VERSION OF BALBAK,
!     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).

!     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL (dp)
!     COMPLEX MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
!     BALANCED MATRIX DETERMINED BY DCBAL.

!     ON INPUT-

!        N IS THE ORDER OF THE MATRIX,

!        LOW AND IGH ARE INTEGERS DETERMINED BY DCBAL,

!        SCALE CONTAINS INFORMATION DETERMINING THE PERMUTATIONS
!          AND SCALING FACTORS USED BY DCBAL,

!        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,

!        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY, OF THE
!          EIGENVECTORS TO BE BACK TRANSFORMED IN THEIR FIRST M COLUMNS.

!     ON OUTPUT-

!        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY, OF THE
!          TRANSFORMED EIGENVECTORS IN THEIR FIRST M COLUMNS.

!-----------------------------------------------------------------------
INTEGER   :: i, j, k, ii
REAL (dp) :: s

IF (m == 0) GO TO 200
IF (igh == low) GO TO 120

DO i = low, igh
  s = scale(i)
!     ********** LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
!                IF THE FOREGOING STATEMENT IS REPLACED BY
!                S=1.0/SCALE(I). **********
  DO j = 1, m
    zr(i,j) = zr(i,j) * s
    zi(i,j) = zi(i,j) * s
  END DO
  
END DO

!     ********** FOR I=LOW-1 STEP -1 UNTIL 1,
!                IGH+1 STEP 1 UNTIL N DO -- **********
120 DO ii = 1, n
  i = ii
  IF (i >= low .AND. i <= igh) CYCLE
  IF (i < low) i = low - ii
  k = scale(i)
  IF (k == i) CYCLE
  
  DO j = 1, m
    s = zr(i,j)
    zr(i,j) = zr(k,j)
    zr(k,j) = s
    s = zi(i,j)
    zi(i,j) = zi(k,j)
    zi(k,j) = s
  END DO
  
END DO

200 RETURN
!     ********** LAST CARD OF DCBABK **********
END SUBROUTINE dcbabk



SUBROUTINE dcmqr2 (n, low, igh, ortr, orti, hr, hi, wr, wi, zr, zi, ierr)

INTEGER, INTENT(IN)       :: n
INTEGER, INTENT(IN)       :: low
INTEGER, INTENT(IN)       :: igh
REAL (dp), INTENT(IN OUT) :: ortr(:)
REAL (dp), INTENT(IN OUT) :: orti(:)
REAL (dp), INTENT(IN OUT) :: hr(:,:)      ! hr(nm,n)
REAL (dp), INTENT(IN OUT) :: hi(:,:)      ! hi(nm,n)
REAL (dp), INTENT(OUT)    :: wr(:)
REAL (dp), INTENT(OUT)    :: wi(:)
REAL (dp), INTENT(OUT)    :: zr(:,:)      ! zr(nm,n)
REAL (dp), INTENT(OUT)    :: zi(:,:)      ! zi(nm,n)
INTEGER, INTENT(OUT)      :: ierr

!-----------------------------------------------------------------------

!   THIS SUBROUTINE IS A TRANSLATION OF A UNITARY ANALOGUE OF THE ALGOL
!   PROCEDURE  COMLR2, NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
!   HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
!   THE UNITARY ANALOGUE SUBSTITUTES THE QR ALGORITHM OF FRANCIS
!   (COMP. JOUR. 4, 332-345(1962)) FOR THE LR ALGORITHM.

!   THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS OF A REAL (dp)
!   COMPLEX UPPER HESSENBERG MATRIX BY THE QR METHOD.  THE EIGENVECTORS OF A
!   COMPLEX GENERAL MATRIX CAN ALSO BE FOUND IF DCORTH HAS BEEN USED TO
!   REDUCE THIS GENERAL MATRIX TO HESSENBERG FORM.

!   ON INPUT-

!      N IS THE ORDER OF THE MATRIX,

!      LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING SUBROUTINE DCBAL.
!        IF DCBAL HAS NOT BEEN USED, SET LOW=1, IGH=N,

!      ORTR AND ORTI CONTAIN INFORMATION ABOUT THE UNITARY TRANSFORMATIONS
!        USED IN THE REDUCTION BY DCORTH, IF PERFORMED.
!        ONLY ELEMENTS LOW THROUGH IGH ARE USED.  IF THE EIGENVECTORS
!        OF THE HESSENBERG MATRIX ARE DESIRED, SET ORTR(J) AND
!        ORTI(J) TO 0.0 FOR THESE ELEMENTS,

!      HR AND HI CONTAIN THE REAL AND IMAGINARY PARTS,
!        RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX.
!        THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAIN FURTHER
!        INFORMATION ABOUT THE TRANSFORMATIONS WHICH WERE USED IN THE
!        REDUCTION BY DCORTH, IF PERFORMED.  IF THE EIGENVECTORS OF THE
!        HESSENBERG MATRIX ARE DESIRED, THESE ELEMENTS MAY BE ARBITRARY.

!   ON OUTPUT-

!      ORTR, ORTI, AND THE UPPER HESSENBERG PORTIONS OF HR AND HI
!        HAVE BEEN DESTROYED,

!      WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY, OF THE
!        EIGENVALUES.  IF AN ERROR EXIT IS MADE, THE EIGENVALUES SHOULD
!        BE CORRECT FOR INDICES IERR+1,...,N,

!      ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY, OF THE
!        EIGENVECTORS.  THE EIGENVECTORS ARE UNNORMALIZED.
!        IF AN ERROR EXIT IS MADE, NONE OF THE EIGENVECTORS HAS BEEN FOUND,

!      IERR IS SET TO
!        ZERO   FOR NORMAL RETURN,
!        J      IF THE J-TH EIGENVALUE HAS NOT BEEN DETERMINED AFTER
!               50 ITERATIONS.

!-----------------------------------------------------------------------
INTEGER   :: en, enm1, i, j, k, l, m, ll, ip1, its, lp1, iend
REAL (dp) :: si, sr, ti, tr, xi, xr, yi, yr, zzi, zzr, norm, machep
REAL (dp) :: r2, w(2), z(2)

!     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.

machep = dpmpar(1)

!                **********

ierr = 0
!     ********** INITIALIZE EIGENVECTOR MATRIX **********
DO i = 1, n
  
  DO j = 1, n
    zr(i,j) = 0.0_dp
    zi(i,j) = 0.0_dp
    IF (i == j) zr(i,j) = 1.0_dp
  END DO
END DO
!     ********** FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS
!                FROM THE INFORMATION LEFT BY DCORTH **********
iend = igh - low - 1
IF (iend < 0) THEN
  GO TO 180
ELSE IF (iend == 0) THEN
  GO TO 150
END IF
!     ********** FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- **********
DO i = igh-1, low+1, -1
  IF (ortr(i) == 0.0_dp .AND. orti(i) == 0.0_dp) CYCLE
  IF (hr(i,i-1) == 0.0_dp .AND. hi(i,i-1) == 0.0_dp) CYCLE
!     ********** NORM BELOW IS NEGATIVE OF H FORMED IN DCORTH **********
  norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)
  ip1 = i + 1
  
  DO k = ip1, igh
    ortr(k) = hr(k,i-1)
    orti(k) = hi(k,i-1)
  END DO
  
  DO j = i, igh
    sr = 0.0_dp
    si = 0.0_dp
    
    DO k = i, igh
      sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
      si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
    END DO
    
    sr = sr / norm
    si = si / norm
    
    DO k = i, igh
      zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
      zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
    END DO
    
  END DO
  
END DO
!     ********** CREATE REAL SUBDIAGONAL ELEMENTS **********
150 l = low + 1

DO i = l, igh
  ll = MIN(i+1, igh)
  IF (hi(i,i-1) == 0.0_dp) CYCLE
  norm = dcpabs(hr(i,i-1), hi(i,i-1))
  yr = hr(i,i-1) / norm
  yi = hi(i,i-1) / norm
  hr(i,i-1) = norm
  hi(i,i-1) = 0.0_dp
  
  DO j = i, n
    si = yr * hi(i,j) - yi * hr(i,j)
    hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
    hi(i,j) = si
  END DO
  
  DO j = 1, ll
    si = yr * hi(j,i) + yi * hr(j,i)
    hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
    hi(j,i) = si
  END DO
  
  DO j = low, igh
    si = yr * zi(j,i) + yi * zr(j,i)
    zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
    zi(j,i) = si
  END DO
  
END DO
!     ********** STORE ROOTS ISOLATED BY DCBAL **********
180 DO i = 1, n
  IF (i >= low .AND. i <= igh) CYCLE
  wr(i) = hr(i,i)
  wi(i) = hi(i,i)
END DO

en = igh
tr = 0.0_dp
ti = 0.0_dp
!     ********** SEARCH FOR NEXT EIGENVALUE **********
220 IF (en < low) GO TO 680
its = 0
enm1 = en - 1
!     ********** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
!                FOR L=EN STEP -1 UNTIL LOW DO -- **********
240 DO ll = low, en
  l = en + low - ll
  IF (l == low) EXIT
  IF (ABS(hr(l,l-1)) <= machep * (ABS(hr(l-1,l-1)) + ABS(hi(l-1,l-1))  &
      + ABS(hr(l,l)) + ABS(hi(l,l)))) EXIT
END DO
!     ********** FORM SHIFT **********
IF (l == en) GO TO 660
IF (its == 50) GO TO 1000
IF (its == 10 .OR. its == 20 .OR. its == 30) GO TO 320
sr = hr(en,en)
si = hi(en,en)
xr = hr(enm1,en) * hr(en,enm1)
xi = hi(enm1,en) * hr(en,enm1)
IF (xr == 0.0_dp .AND. xi == 0.0_dp) GO TO 340
yr = (hr(enm1,enm1) - sr) / 2.0_dp
yi = (hi(enm1,enm1) - si) / 2.0_dp
z(1) = yr*yr - yi*yi + xr
z(2) = 2.0_dp*yr*yi + xi
CALL dcsqrt(z, w)
zzr = w(1)
zzi = w(2)
IF (yr * zzr + yi * zzi >= 0.0_dp) GO TO 310
zzr = -zzr
zzi = -zzi
310 z(1) = yr + zzr
z(2) = yi + zzi
r2 = z(1)**2 + z(2)**2
sr = sr - (xr*z(1) + xi*z(2))/r2
si = si - (xi*z(1) - xr*z(2))/r2
GO TO 340
!     ********** FORM EXCEPTIONAL SHIFT **********
320 sr = ABS(hr(en,enm1)) + ABS(hr(enm1,en-2))
si = 0.0_dp

340 DO i = low, en
  hr(i,i) = hr(i,i) - sr
  hi(i,i) = hi(i,i) - si
END DO

tr = tr + sr
ti = ti + si
its = its + 1
!     ********** REDUCE TO TRIANGLE (ROWS) **********
lp1 = l + 1

DO i = lp1, en
  sr = hr(i,i-1)
  hr(i,i-1) = 0.0_dp
  norm = SQRT(hr(i-1,i-1)*hr(i-1,i-1) + hi(i-1,i-1)*hi(i-1,i-1) + sr*sr)
  xr = hr(i-1,i-1) / norm
  wr(i-1) = xr
  xi = hi(i-1,i-1) / norm
  wi(i-1) = xi
  hr(i-1,i-1) = norm
  hi(i-1,i-1) = 0.0_dp
  hi(i,i-1) = sr / norm
  
  DO j = i, n
    yr = hr(i-1,j)
    yi = hi(i-1,j)
    zzr = hr(i,j)
    zzi = hi(i,j)
    hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
    hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
    hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
    hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  END DO
  
END DO

si = hi(en,en)
IF (si == 0.0_dp) GO TO 540
norm = dcpabs(hr(en, en), si)
sr = hr(en,en) / norm
si = si / norm
hr(en,en) = norm
hi(en,en) = 0.0_dp
IF (en == n) GO TO 540
ip1 = en + 1

DO j = ip1, n
  yr = hr(en,j)
  yi = hi(en,j)
  hr(en,j) = sr * yr + si * yi
  hi(en,j) = sr * yi - si * yr
END DO
!     ********** INVERSE OPERATION (COLUMNS) **********
540 DO j = lp1, en
  xr = wr(j-1)
  xi = wi(j-1)
  
  DO i = 1, j
    yr = hr(i,j-1)
    yi = 0.0_dp
    zzr = hr(i,j)
    zzi = hi(i,j)
    IF (i == j) GO TO 560
    yi = hi(i,j-1)
    hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
    560 hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
    hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
    hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  END DO
  
  DO i = low, igh
    yr = zr(i,j-1)
    yi = zi(i,j-1)
    zzr = zr(i,j)
    zzi = zi(i,j)
    zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
    zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
    zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
    zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  END DO
  
END DO

IF (si == 0.0_dp) GO TO 240

DO i = 1, en
  yr = hr(i,en)
  yi = hi(i,en)
  hr(i,en) = sr * yr - si * yi
  hi(i,en) = sr * yi + si * yr
END DO

DO i = low, igh
  yr = zr(i,en)
  yi = zi(i,en)
  zr(i,en) = sr * yr - si * yi
  zi(i,en) = sr * yi + si * yr
END DO

GO TO 240
!     ********** A ROOT FOUND **********
660 hr(en,en) = hr(en,en) + tr
wr(en) = hr(en,en)
hi(en,en) = hi(en,en) + ti
wi(en) = hi(en,en)
en = enm1
GO TO 220
!     ********** ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
!                VECTORS OF UPPER TRIANGULAR FORM **********
680 norm = 0.0_dp

DO i = 1, n
  DO j = i, n
    norm = norm + ABS(hr(i,j)) + ABS(hi(i,j))
  END DO
END DO

IF (n == 1 .OR. norm == 0.0_dp) GO TO 1001
!     ********** FOR EN=N STEP -1 UNTIL 2 DO -- **********
DO en = n, 2, -1
  xr = wr(en)
  xi = wi(en)
  enm1 = en - 1
!     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- **********
  DO i = enm1, 1, -1
    zzr = hr(i,en)
    zzi = hi(i,en)
    IF (i == enm1) GO TO 760
    ip1 = i + 1
    
    DO j = ip1, enm1
      zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
      zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
    END DO
    
    760 yr = xr - wr(i)
    yi = xi - wi(i)
    IF (yr == 0.0_dp .AND. yi == 0.0_dp) yr = machep * norm
    r2 = yr*yr + yi*yi
    hr(i,en) = (zzr*yr + zzi*yi) / r2
    hi(i,en) = (zzi*yr - zzr*yi) / r2
  END DO
  
END DO
!     ********** END BACKSUBSTITUTION **********

!     ********** VECTORS OF ISOLATED ROOTS **********
DO i = 1, n-1
  IF (i >= low .AND. i <= igh) CYCLE
  ip1 = i + 1
  
  DO j = ip1, n
    zr(i,j) = hr(i,j)
    zi(i,j) = hi(i,j)
  END DO
  
END DO
!     ********** MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
!                VECTORS OF ORIGINAL FULL MATRIX.
!                FOR J=N STEP -1 UNTIL LOW+1 DO -- **********
DO j = n, low+1, -1
  m = MIN(j-1,igh)
  
  DO i = low, igh
    zzr = zr(i,j)
    zzi = zi(i,j)
    
    DO k = low, m
      zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
      zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
    END DO
    
    zr(i,j) = zzr
    zi(i,j) = zzi
  END DO
END DO

GO TO 1001
!     ********** SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 50 ITERATIONS **********
1000 ierr = en
1001 RETURN
!     ********** LAST CARD OF DCMQR2 **********
END SUBROUTINE dcmqr2

END MODULE complex_eigen
