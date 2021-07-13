c
c  Vectores y definiciones utiles para manejar las rutinas de QUASI.f
c
c   UV => U1PP U2PP U1PM U2PM  V1PP V2PP V1PM V2PM     
c         U1NP U2NP U1NM U2NM  V1NP V2NP V1NM V2NM     
c
c   IT = 1   Protones positivos
c        2   Protones negativos
c        3   Neutrones positivos
c        4   Neutrones negativos 
c
      Subroutine ITLOC(NP,NM,IBLOCK)
      
      Parameter (NCM=6)
      Dimension IBLOCK(4)
      
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC1/ NNU(4),NNV(4),NNRO(4),NNKA(4)
      Common /ITLOC2/ ND2PACK(4),NNGPACK(4),NNDPACK(4)
      Common /ITLOC3/ NNC20(4,NCM),NNC11(4,NCM)
      Common /ITLOC4/ NNH20(4),NNH11(4),NNEQP(4)
      Common /DIMS  / INP,INM,INROP,INROM
            
      INP = NP
      INM = NM
      
      NP2 = NP*NP
      NM2 = NM*NM
      
      NGP = (NP*(NP+1))/2
      NGM = (NM*(NM+1))/2
      
      INROP = NGP
      INROM = NGM
c
c     Dimensiones para cada valor de IT
c      
      ND(1) = NP
      ND(2) = NM
      ND(3) = NP
      ND(4) = NM
c
c     Dimensiones al cuadrado para cada valor de IT
c      
      ND2(1) = NP2
      ND2(2) = NM2
      ND2(3) = NP2
      ND2(4) = NM2
c
c     NBLOCK = nb(+i) - nb(-i) 
c
c     where nb(s) is the number of blocked levels of signature s
c
c     N1(it) = ND(it) - NBLOCK(it)
c     N2(it) = ND(it) + NBLOCK(it)
c
      NBLOCK(1) = IBLOCK(1)
      NBLOCK(2) = IBLOCK(2)
      NBLOCK(3) = IBLOCK(3)
      NBLOCK(4) = IBLOCK(4)
      
c           Posiciones de U y V en UV(8*NP*NP + 8*NM*NM) 
c    
C          UV => U1PP U2PP U1PM U2PM  V1PP V2PP V1PM V2PM     
C                U1NP U2NP U1NM U2NM  V1NP V2NP V1NM V2NM     

      NNU(1) = 1
      NNU(2) = 1+2*NP2
      NNU(3) = 1      +4*NP2+4*NM2
      NNU(4) = 1+2*NP2+4*NP2+4*NM2
      
      NNV(1) = NNU(1)+2*NP2+2*NM2
      NNV(2) = NNU(2)+2*NP2+2*NM2
      NNV(3) = NNU(3)+2*NP2+2*NM2
      NNV(4) = NNU(4)+2*NP2+2*NM2
      
C   RO => RO1PP RO2PP RO1PM RO2PM RO1NP RO2NP RO1NM RO2NM 
    
      NNRO(1) = 1
      NNRO(2) = 1+2*NP2
      NNRO(3) = 1      +2*NP2+2*NM2
      NNRO(4) = 1+2*NP2+2*NP2+2*NM2
      
C   KAPPA => KA1PP KA1PM KA1NP KA1NM
    
      NNKA(1) = 1
      NNKA(2) = 1+NP2
      NNKA(3) = 1    +NP2+NM2
      NNKA(4) = 1+NP2+NP2+NM2
      
c 20 parts of the constraints. It is assumed that the constrained operators
c are positive signature ones and therefore the dimension of the 20 part
c is N1xN2
c
c   C20 --> C20PP, C20PM, C20NP, C20NM

      
      JSHC20 = 0
      do it=1,4
         JSHC20 = JSHC20 + ND(it)**2-NBLOCK(it)**2 ! N1xN2
      end do
      
      do i=1,NCM
         NNC20(1,i) = 1 + (i-1)*JSHC20
         do it=2,4
            NNC20(it,i) = NNC20(it-1,i)+ND(it-1)**2-NBLOCK(it-1)**2
	 end do 
      end do
      
c   H20 --> H20PP, H20PM, H20NP, H20NM

       NNH20(1) = 1 
       do it=2,4
          NNH20(it) = NNH20(it-1) + ND(it-1)**2-NBLOCK(it-1)**2
       end do

c 11 parts of the constraints. It is assumed that the constrained operators
c are positive signature ones and therefore the dimension of the 11 parts
c are N1xN1 and N2xN2
c
c   C11 --> C11PP_1,C11PP_2, C11PM_1, C11PM_2,
c           C11NP_1,C11NP_2, C11NM_1, C11NM_2

      
      JSHC11 = 0
      do it=1,4
         JSHC11 = JSHC11 + 2*(ND(it)**2 + NBLOCK(it)**2) ! N1xN1+N2xN2
      end do 
      
      do i=1,NCM
         NNC11(1,i) = 1 + (i-1)*JSHC11
	 do it=2,4
            NNC11(it,i)=NNC11(it-1,i)+2*(ND(it-1)**2+NBLOCK(it-1)**2)
	 end do
      end do
      
c   H11 --> H11PP_1,H11PP_2, H11PM_1, H11PM_2,
c           H11NP_1,H11NP_2, H11NM_1, H11NM_2 

       NNH11(1) = 1 
       do it=2,4
          NNH11(it) = NNH11(it-1) + 2*(ND(it-1)**2+NBLOCK(it-1)**2)
       end do
       
c
c    EQP
c
       NNEQP(1) = 1
       NNEQP(2) = 1 + 2*NP
       NNEQP(3) = 1 + 2*NP + 2*NM
       NNEQP(4) = 1 + 4*NP + 2*NM

c
c   GFIELD   G1PP G2PP G1NP G2NP D1PP D2PP D1NP D2NP 
c
       ND2PACK(1) = NGP
       ND2PACK(2) = NGM
       ND2PACK(3) = NGP
       ND2PACK(4) = NGM
       
       NNGPACK(1) = 1
       NNGPACK(2) = 1 + 8*NGP 
       NNGPACK(3) = 1 + 2*NGP
       NNGPACK(4) = 1 + 8*NGP + 2*NGM
       
       NNDPACK(1) = NNGPACK(1)+4*NGP
       NNDPACK(2) = NNGPACK(2)+4*NGM
       NNDPACK(3) = NNGPACK(3)+4*NGP
       NNDPACK(4) = NNGPACK(4)+4*NGM
       
       return
       end
