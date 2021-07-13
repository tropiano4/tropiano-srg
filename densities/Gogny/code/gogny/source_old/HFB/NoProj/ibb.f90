!-----------------------------------------------------------------------!
!									!
!									!
!     RADIAL INTEGRALS: BRINK-BOKER FORCE				!
!									!
!  This module computes the radial integral of tge Brink-Boeker force 	!
!  in the 2 distinct cases of an "analytical" harmonic oscillator basis !
!  and a general Woods-Saxon (or other) basis.				!
!									!
!-----------------------------------------------------------------------!

 MODULE ibb

	USE input
	USE global
	USE lgfactor
	USE symtalm

	IMPLICIT NONE

	INTEGER :: Ndmax = 700
	
	! PI_BrinkBooker = 2.0 * (PI ** 3.0)
	REAL(KIND = 16), PARAMETER :: PI_BB = 62.01255336059964035
	
	REAL(KIND = 16), DIMENSION(:, :, :, :, :), ALLOCATABLE :: StoreIntNeu
	REAL(KIND = 16), DIMENSION(:, :, :, :, :), ALLOCATABLE :: StoreIntPro

	PRIVATE StoreIntNeu, StoreIntPro

 CONTAINS

        !---------------------------------------------------------------------------------!
	! Function giving the radial integral IBB for the Brink-Boeker term in the case   !
	! of the spherical harmonic oscillator basis	  			          !
	! Ref.:
        !---------------------------------------------------------------------------------!

	FUNCTION IBrinkBookerHO(na, la, nb, lb, nc, lc, nd, ld, k, x)
		DOUBLE PRECISION IBrinkBookerHO
		INTEGER, INTENT(IN) :: na, la, nb, lb, nc, lc, nd, ld, k
		DOUBLE PRECISION, INTENT(IN) :: x ! x = mi(i) / b
		INTEGER N1max, N2max, N1min, N2min
		DOUBLE PRECISION d1, d2
		REAL(KIND = 16) sumN1
		INTEGER N1

		N1max = na + nc + ((la + lc - k) / 2)
		N2max = nb + nd + ((lb + ld - k) / 2)
		N1min = MIN_5N(na, la, nc, lc, k)
		N2min = MIN_5N(nb, lb, nd, ld, k)

		! El sumatorio se deberia realizar en el sentido opuesto
		! para ir de menor a mayor en el orden de magnitud de los numeros
		
		d1 = x * x
		d2 = d1 + 2.0
		
		sumN1 = SymKumar_get(na, la, nc, lc, N1max, k) * SumBB(nb, lb, nd, ld, k, N1max, d2)
		
		DO N1 = N1max - 1, N1min, -1
			sumN1 = (SymKumar_get(na, la, nc, lc, N1, k) * SumBB(nb, lb, nd, ld, k, N1, d2)) &
				- (sumN1 * (N1 + N2min + k + 1.5) / (N1 + 1.0) / (N1 + k + 1.5) / d2)
		END DO

		IBrinkBookerHO = PI_BB * PAR(N1min + N2min) &
			* EXP(DDLogSemiFactorials(N1min + N2min + k) &
		            - DDLogFactorials(N1min) &
			    - DDLogSemiFactorials(N1min+ k) &
			    - DDLogFactorials(N2min) &
			    - DDLogSemiFactorials(N2min + k)) &
			* (x ** 3.0) * sumN1 / (d2 ** (N1min + N2min + k + 1.5))
		RETURN

	CONTAINS

       		!-----------------------------------------------------------------------!
		!  Sum over N2 of the c(N2,k) T(nb,lb,nd,ld,N2,k)*Integral(k,N1,)	1
       		!-----------------------------------------------------------------------!

		FUNCTION SumBB(n1, l1, n2, l2, k, M1, y) ! M1 <- N1
			REAL(KIND = 16) SumBB
			INTEGER, INTENT(IN) :: n1, l1, n2, l2, k, M1
			DOUBLE PRECISION, INTENT(IN) :: y

			INTEGER N, Nmin, Nmax

			Nmax = n1 + n2 + ((l1 + l2 - k) / 2)
			Nmin = MIN_5N(n1, l1, n2, l2, k)

			SumBB = SymKumar_get(n1, l1, n2, l2, Nmax, k)
			
			DO N = Nmax - 1, Nmin, -1
				SumBB = - (SumBB * (N + M1 + k + 1.5) / (N + 1.0) / (N + k + 1.5) / y) &
					+ SymKumar_get(n1, l1, n2, l2, N, k)
			END DO
			
			RETURN
		END FUNCTION SumBB

	END FUNCTION IBrinkBookerHO

        !-------------------------------------------------------------------------------!
	!										!
	!   Subroutine that initializes an array that will contain the one-body 	!
	!   integrals required in the calculation of the Brink-Boker radial integral	!
	!										!
	!   Note: The two arrays StoreIntNeu and StoreIntPro occupy a VERY large space	!
	!         in memory. However, it is absolutely compulsory to use them if one	!
	!	  wants to obtain results in a human timescale... Typically, we can	!
	!	  avoid the calculation of several dozens of millions of integrals by	!
	!	  checking if they are stored.						!
	!										!
        !-------------------------------------------------------------------------------!

	SUBROUTINE StoreInt_new(Lmin)
		INTEGER, INTENT(IN) :: Lmin
		
		INTEGER ::la, na, lb, nb, n
		
		IF (Lmin .EQ.0) THEN
		
			ALLOCATE(StoreIntNeu(0:Lmax,1:NmaxOfL(0),0:Lmax,1:NmaxOfL(0),0:Ndmax))
			ALLOCATE(StoreIntPro(0:Lmax,1:NmaxOfL(0),0:Lmax,1:NmaxOfL(0),0:Ndmax))

                	DO la = 0, Lmax
                       		DO na = 1,NmaxOfL(0)
                                	DO lb = 0, Lmax
                                        	DO nb = 1, NmaxOfL(0)
                                                	DO n=0, Ndmax
                                                      		StoreIntNeu(la, na, lb, nb, n) = 999999.9
								StoreIntPro(la, na, lb, nb, n) = 999999.9
                                                	END DO
                                        	END DO
                                	END DO
                        	END DO
                	END DO
			
		ELSE
	
			ALLOCATE(StoreIntNeu(0:1,1:NmaxOfL(0),0:Lmax,1:NmaxOfL(0),0:Ndmax))
			ALLOCATE(StoreIntPro(0:1,1:NmaxOfL(0),0:Lmax,1:NmaxOfL(0),0:Ndmax))

                       	DO na = 1,NmaxOfL(0)
                                DO lb = 0, Lmax
                                        DO nb = 1, NmaxOfL(0)
                                                DO n=0, Ndmax
                                                      	  StoreIntNeu(0, na, lb, nb, n) = 999999.9
                                                      	  StoreIntPro(0, na, lb, nb, n) = 999999.9
                                                      	  StoreIntNeu(1, na, lb, nb, n) = 999999.9
                                                      	  StoreIntPro(1, na, lb, nb, n) = 999999.9
                                                END DO
                                        END DO
                                END DO
                        END DO
			
		END IF

		RETURN
	END SUBROUTINE StoreInt_new

        !-------------------------------------------------------------------------------!
	! Subroutine that deallocates memory for the array that contains the one-body 	!
	! integrals required in the calculation of the Brink-Boker radial integral.	!
        !-------------------------------------------------------------------------------!

	SUBROUTINE StoreInt_del(Lmin)
		INTEGER, INTENT(IN) :: Lmin
	
		DEALLOCATE(StoreIntNeu)
		DEALLOCATE(StoreIntPro)

		RETURN
	END SUBROUTINE StoreInt_del

        !---------------------------------------------------------------------------------!
	! Function giving the radial integral IBB for the Brink-Boeker term in the case   !
	! of a general spherical basis			  			          !
        !---------------------------------------------------------------------------------!

	FUNCTION IBrinkBooker(na, la, nb, lb, nc, lc, nd, ld, k, x, icount)
		DOUBLE PRECISION IBrinkBooker
		
		INTEGER, INTENT(IN) :: na, la, nb, lb, nc, lc, nd, ld
		INTEGER, INTENT(OUT) :: icount
		INTEGER, INTENT(IN) ::  k
		DOUBLE PRECISION, INTENT(IN) :: x ! x = mu(i)
		
		INTEGER :: s, n, Index1, Index2, Index3, Index4, Ipoint, MinLoops
		
		REAL(KIND = 16), ALLOCATABLE :: Functi1(:), Functi2(:)
		REAL(KIND = 16) :: epsilo, Pi, SumNew, SumOld, Delta, SumTot
		REAL(KIND = 16) :: Term1, Term2, Term3, Fact1, Fact2, FuncMax
		
		Index1 = IndexVecNL(na, la)
		Index2 = IndexVecNL(nb, lb)
		Index3 = IndexVecNL(nc, lc)
		Index4 = IndexVecNL(nd, ld)
		
		IF (Index1 .EQ. 0 .OR. Index2 .EQ. 0 .OR. Index3 .EQ. 0 .OR. Index4 .EQ. 0) THEN
			IBrinkBooker = 0.0
			RETURN
		END IF						
		
		ALLOCATE(Functi1(1:Npoint))	
		ALLOCATE(Functi2(1:Npoint))	
		
		DO Ipoint=1,Npoint
			Functi1(Ipoint) = WaveFun(Ipoint,Index1)* WaveFun(Ipoint,Index3)
			Functi2(Ipoint) = WaveFun(Ipoint,Index2)* WaveFun(Ipoint,Index4)
		END DO			
			
		Pi = 4.0*ATAN(1.0)
		
		SumOld = 10.0
		SumNew = 9.0
                MinLoops = 50
		Delta = ABS(SumNew - SumOld)
		
		epsilo = 1.e-18
		Sumtot = 0.0
		FuncMax = 0.0

		Term1 = 0.0
		Term2 = 0.0
		
		Fact1 = 0.0
		Fact2 = 0.0
		
		s = 0
		
                DO WHILE (ABS(Delta) .GT. epsilo .OR. s .LE. MinLoops) 
		
			n = k + 2*s
		
			! FuncMax is a totally arbitrary normalization coefficients used to keep the integrals within  
			! the bound of sensible values. Naturally, at the end we must cancel off these spurious terms.
		
			Term1 = OneBodyInt(Functi1, na, la, nc, lc, n, x, FuncMax, icount) 
			Fact1 = EXP( -DDLogFactorials(s) - FuncMax - n*LOG(x))
		   
			Term2 = OneBodyInt(Functi2, nb, lb, nd, ld, n, x, FuncMax, icount)
			Fact2 = EXP( -DDLogSemiFactorials(k+s) - FuncMax - n*LOG(x))
		   
			SumNew = Fact1 * Term1 * Fact2 * Term2
			
			SumTot = SumTot + SumNew
			
			Delta = ABS(SumNew - SumOld)
			SumOld = SumNew
			FuncMax = 0.0
		   
			s = s + 1
			
		END DO
		
		IBrinkBooker = 2.0*Pi*SQRT(Pi) * SumTot

		RETURN

	CONTAINS

       		!---------------------------------------------------------------------------------------!
		!											!
		!    This function gives a one-body integral which we need in the series expansion	!
		!    for the radial Brink-Boker integral. Two remarks must be done:			!
		!     - Since the calculation of the series is very costly, we implement a trick:	!
		!       We check if the integral was already calculated and was stored in StoreIntNeu 	!
		!	or StoreIntPro. In the case labelled "optimization", we calculate the matrix	!
		!	elements only for la = Lmax (see loops in subroutines SymVBBph_calculate and 	!
		!	SymVBBpp_calculate). In this case, we can reduce the size of the matrices 	!
		!	needed to store the integrals significantly (by a factor (Lmax+1)/2.		!
		!     - Since this one-body integral can be very large, and it has to be multiplied	!
		!       by very small numbers, we artifically added a normalization constant FuncMax	!
		!       Its role is simply to make sure we don't overflow.				!
		!											!
        	!---------------------------------------------------------------------------------------!

		FUNCTION OneBodyInt(Functi, na, la, nb, lb, n, x, FuncMax, icount)
			DOUBLE PRECISION OneBodyInt
			
			INTEGER, INTENT(IN) :: n, na, la, nb, lb
			INTEGER, INTENT(OUT) :: icount
			DOUBLE PRECISION, INTENT(IN) :: x
		
			REAL(KIND = 16), INTENT(IN) :: Functi(:)
			REAL(KIND = 16), INTENT(OUT) :: FuncMax
			
			REAL(KIND = 16) :: Result, mu, h
			REAL(KIND = 16), ALLOCATABLE :: Integrand(:)
			
			IF (n .GT. 0) THEN
				FuncMax = -DDLogFactorials(n) - n*LOG(0.5*x*x) + n*LOG(2.0)
			ELSE
				FuncMax = 0.0
			END IF
			
			! Finding out if the integral was not already calculated. We consider 2 cases, depending on
			! whether we "optimize" the calculation of the matrix elements.
			! Case no optimization: the integral I(la, na, lb, nb, n) is stored in StoreIntNeu (n) and 
			!			StoreIntPro (p)
			! Case optimization: Here, we calculate the matrix element corresponding to la = Lmax (in symvbb.f90)
			!		     <la | lb >, la = Lmax, any lb -> StoreIntPro(0, na, la, nb, n)
			!		     <la | lb >, lb = Lmax, any la -> StoreIntPro(0, nb, la, na, n)
			!		     <lb | lb >, any lb  	   -> StoreIntPro(1, nb, la, na, n)
			
			IF (n .LE. Ndmax) THEN
			
				IF (Optimization .EQ. 1) THEN
					IF (ABS(x - DBLE(0.7)) .LT. 1.e-10) THEN
						IF (la .EQ. Lmax) OneBodyInt = StoreIntPro(0, na, lb, nb, n)
						IF (la .NE. Lmax .AND. lb .EQ. Lmax) OneBodyInt = StoreIntPro(0, nb, la, na, n)
						IF (la .NE. Lmax .AND. lb .NE. Lmax) OneBodyInt = StoreIntPro(1, nb, la, na, n)
					END IF
					IF (ABS(x - DBLE(1.2)) .LT. 1.e-10) THEN
						IF (la .EQ. Lmax) OneBodyInt = StoreIntNeu(0, na, lb, nb, n)
						IF (la .NE. Lmax .AND. lb .EQ. Lmax) OneBodyInt = StoreIntNeu(0, nb, la, na, n)
						IF (la .NE. Lmax .AND. lb .NE. Lmax) OneBodyInt = StoreIntNeu(1, nb, la, na, n)
					END IF
				ELSE
					IF (ABS(x - DBLE(0.7)) .LT. 1.e-10) OneBodyInt = StoreIntPro(la, na, lb, nb, n)
					IF (ABS(x - DBLE(1.2)) .LT. 1.e-10) OneBodyInt = StoreIntNeu(la, na, lb, nb, n)
				END IF
				
			END IF
			
			IF (ABS(OneBodyInt - DBLE(999999.9)) .GT. 1.e-10) THEN
				icount = icount + 1
				RETURN
			END IF
			
	       	 	ALLOCATE(Integrand(1:Npoint))		

			! Defining the function to integrate
			
			mu = x
			CALL DefineIntegrand(Functi, Integrand, n, mu, FuncMax)
					
			! Calculating the integral by the Simpson method
			
			h = MeshStep
			CALL Simpson_Kind16(Integrand,Npoint,h,Result)
		
			DEALLOCATE(Integrand)
		
			OneBodyInt = Result
			
			! Filling in the arrays that will give us the radial integrals. See above for a quick comment about
			! the different cases
			
			IF (n .LE. Ndmax) THEN
			
				IF (Optimization .EQ. 1) THEN
					IF (ABS(x - DBLE(0.7)) .LT. 1.e-10) THEN
                       				IF (la .EQ. Lmax) StoreIntPro(0, na, lb, nb, n) = OneBodyInt
                       				IF (la .NE. Lmax .AND. lb .EQ. Lmax) StoreIntPro(0, nb, la, na, n) = OneBodyInt					
                       				IF (la .NE. Lmax .AND. lb .NE. Lmax) StoreIntPro(1, nb, la, na, n) = OneBodyInt					
					END IF
					IF (ABS(x - DBLE(1.2)) .LT. 1.e-10) THEN
                       				IF (la .EQ. Lmax) StoreIntNeu(0, na, lb, nb, n) = OneBodyInt
                       				IF (la .NE. Lmax .AND. lb .EQ. Lmax) StoreIntNeu(0, nb, la, na, n) = OneBodyInt					
                       				IF (la .NE. Lmax .AND. lb .NE. Lmax) StoreIntNeu(1, nb, la, na, n) = OneBodyInt					
					END IF
				ELSE
                       			IF (ABS(x - DBLE(0.7)) .LT. 1.e-10) THEN
						StoreIntPro(la, na, lb, nb, n) = OneBodyInt
						StoreIntPro(lb, nb, la, na, n) = OneBodyInt
					END IF
					IF (ABS(x - DBLE(1.2)) .LT. 1.e-10) THEN
						StoreIntNeu(la, na, lb, nb, n) = OneBodyInt
						StoreIntNeu(lb, nb, la, na, n) = OneBodyInt
					END IF
				END IF
				
			END IF
			
			RETURN
		END FUNCTION OneBodyInt
		
		SUBROUTINE DefineIntegrand(Functi, Integrand, n, x, FuncMax)
			REAL(KIND = 16), INTENT(OUT) :: Integrand(:)
			INTEGER, INTENT(IN) :: n
			REAL(KIND = 16), INTENT(IN) :: FuncMax
			REAL(KIND = 16), INTENT(IN) :: x
			REAL(KIND = 16), INTENT(IN) :: Functi(:)
			
			REAL(KIND = 16) :: Radius, ArgMax
			INTEGER :: Ipoint
			
			DO Ipoint = 1, Npoint
				Radius = RadMesh(Ipoint)
				ArgMax = -(Radius/x)**2 + n*LOG(Radius) + FuncMax
		  		Integrand(Ipoint) = Functi(Ipoint) * EXP(ArgMax)
			END DO
			
			RETURN
		END SUBROUTINE DefineIntegrand
			
		
	END FUNCTION IBrinkBooker

END MODULE ibb
