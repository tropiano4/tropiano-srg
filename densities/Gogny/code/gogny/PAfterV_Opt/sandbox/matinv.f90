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

      SUBROUTINE MATINV(N,M,AA,BB,DET)                                

      PARAMETER(EPSMACH=1.D-20)                                    

      REAL*8 AA(N,N),BB(N,M)                                                   

      REAL*8,POINTER ::  PC(:),PL(:),CS(:)                                                  

      REAL*8 PV,PAV,DET,TT    

                                                            

!Initializations :                       

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

        PV=AA(K,K)                                                              

        IK=K                                                                    

        JK=K                                                                    

        PAV=DABS(PV)                                                            

        DO I=K,N                                                                

          DO J=K,N                                                              

            IF (DABS(AA(I,J)).GT.PAV) THEN                                      

              PV=AA(I,J)                                                        

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

            TT=AA(IK,I)                                                         

            AA(IK,I)=AA(K,I)                                                    

            AA(K,I)=TT                                                          

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

!Exchange columns JK and K of matrix AA                                         

            TT=AA(I,JK)                                                         

            AA(I,JK)=AA(I,K)                                                    

            AA(I,K)=TT                                                          

          ENDDO                                                                 

        ENDIF                                                                   

!Pivot is at correct column and located in K,K                                              

                                                                               

!Store column K in vector CS                             

!then set column K to zero                                             

        DO I=1,N                                                                

          CS(I)=AA(I,K)                                                         

          AA(I,K)=0.D0                                                          

        ENDDO                                                                   

!                                                                               

        CS(K)=0.                                                                

        AA(K,K)=1.                                                              

!Modify line K :                                            

        IF(DABS(PV).LT.EPSMACH) THEN                                            

          WRITE(*,*) '  PIVOT TOO SMALL - STOP'                               

          STOP                                                                  

        ENDIF                                                                   

        DO I=1,N                                                                

          AA(K,I)=AA(K,I)/PV                                                    

        ENDDO                                                                   

        IF (M.NE.0) THEN                                                        

          DO I=1,M                                                             

            BB(K,I)=BB(K,I)/PV                                                  

          ENDDO                                                                 

        ENDIF                                                                   

!Modify other lines of matrix AA:                                        

        DO J=1,N                                                                

          IF (J.EQ.K) CONTINUE                                                  

          DO I=1,N                                                              

!Modify line J of matrix AA :                                            

            AA(J,I)=AA(J,I)-CS(J)*AA(K,I)                                       

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

                                                                               

!The matrix AA is inverted - Rearrange AA                         

                                                                               

!Exchange lines                                                            

      DO I=N,1,-1                                                               

        IK=PC(I)                                                                

        IF (IK.EQ.I) CONTINUE                                                   

!EXCHANGE LINES I AND PC(I) OF AA:                                         

        DO J=1,N                                                                

          TT=AA(I,J)                                                            

          AA(I,J)=AA(IK,J)                                                      

          AA(IK,J)=TT                                                           

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

!EXCHANGE COLUMNS J AND PL(J) OF AA :                                       

        DO I=1,N                                                                

          TT=AA(I,J)                                                            

          AA(I,J)=AA(I,JK)                                                      

          AA(I,JK)=TT                                                           

        ENDDO                                                                   

!NO MORE EXCHANGE NEEDED                                                      

!GO TO NEXT COLUMN   

      ENDDO                                                                     

!REARRANGEMENT TERMINATED.                                                        

      RETURN                                                                    

   10 FORMAT(///'  DETERMINANT EQUALS ZERO, NO SOLUTION!')                    

      END                                                                       

 
