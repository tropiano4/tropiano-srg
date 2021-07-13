      SUBROUTINE PRINTD(GM,MG,N1,N2,TITLE)
      IMPLICIT REAL*8 (A-H,O-U)          
      CHARACTER* 8 TITLE                
      DIMENSION GM(MG,1)               
C                                     
      WRITE(6,102) TITLE,N1,N2       
      NPAR = N2/14                  
      NRES = N2 - NPAR*14          
      IF(NPAR.EQ.0) GOTO 10       
C                                
      DO 1 II= 1, NPAR          
      WRITE(6,101) II          
C                             
      IK = (II-1)*14 + 1     
      IKE= IK + 13          
C                          
      DO 2 J = 1,N1                  
    2 WRITE(6,100) (GM(J,L),L=IK,IKE)  
    1 CONTINUE                        
C                                   
   10 IF(NRES.EQ.0) GOTO 20        
      NPAR1 = NPAR + 1            
      WRITE(6,101) NPAR1         
      IK= NPAR*14 + 1           
      IKE = N2                 
      DO 3 JJ = 1,N1          
    3 WRITE(6,100) (GM(JJ,L),L=IK,IKE)     
C                                         
   20 RETURN                             
  102 FORMAT ( /,2X,' MATRIZ ',A8,' N1 N2 ',2I6) 
  101 FORMAT ( 2X,' PARTE ',I4,' DE LA MATRIZ') 
  100 FORMAT ( 2X, 14(D8.2,1X))                
C                                             
      END                                    
