      Implicit real*8 (A-H,O-Z)
      Include 'DIMTRIAX'
c      Parameter (NP=50,NM=34)
c      Parameter (NUV=8*(NP*NP+NM*NM))
      Dimension ND2(4),NNU(4),NNV(4),NNUOLD(4)
      Dimension UVOLD (NUV)
      Dimension UV    (NUV)
      
      NP2 = NP*NP
      NM2 = NM*NM
c
c     Dimensiones al cuadrado para cada valor de IT
c      
      ND2(1) = NP2
      ND2(2) = NM2
      ND2(3) = NP2
      ND2(4) = NM2
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
      
      NNUOLD(1) = 1
      NNUOLD(2) = 1 + 8*NP*NP
      NNUOLD(3) = 1 + 4*NP*NP
      NNUOLD(4) = 1 + 8*NP*NP+4*NM*NM
      
      NLR=10
      NLW=11
      
      read (NLR) bx,by,bz,UVOLD
 
  	  do it=1,4
  	    NC = 2*ND2(it)
  	    call dcopy(NC,UVOLD(NNUOLD(it))   ,1,UV(NNU(it)),1)
  	    call dcopy(NC,UVOLD(NNUOLD(it)+NC),1,UV(NNV(it)),1)
  	  end do
	  
          write(6,'( //,"  READING WITH OLD FORMAT ",//,40("="))')
	  write(6,'( 3x,3f7.4)') bx,by,bz
c
          write(6,'( //,"  WRITTING WITH NEW FORMAT ",//,40("="))')
      write (NLW) bx,by,bz,UV
      stop
      end
