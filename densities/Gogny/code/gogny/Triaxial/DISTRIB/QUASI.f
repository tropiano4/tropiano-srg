c +-------------------------------------------------------------------+
c |                                                                   |
c | Computes the 20 part of an operator with the form respect to      |
c | signature                                                         |
c |                                                                   |
c |      /  O   0  \               /   0    O20 \                     |
c |     |          |              |      T      |                     |
c |      \  0 t O  /               \ -O20    0 /                      |
c |                                                                   |
c |   IT = t                                                          |
c |   IH = hermiticity of O                                           |
c |                                                                   |
c |   U and V are stored as matrices of dimension N x 2N of the form  |
c |                                                                   |
c |      U=(U1,U2)      V=(V1,V2)                                     |
c |                                                                   |
c |   with U1 a matrix N x N1 and U2 a matrix N x N2  (N1+N2=2N)      |
c |                                                                   |
c +-------------------------------------------------------------------+
      Subroutine SO20P(U,V,N1,N2,N,O,AUX,O20,IT,IH)
c
      Implicit Real*8 (A-H,O-Z)
c
      Dimension U(N,*),V(N,*)
      Dimension O(N,N),AUX(N,N2)
      Dimension O20(N1,N2)
C
      if(iabs(it*ih).ne.1) then
        write(6,'(" in SO20P IT IH ",2I3, " ******** ")') it,ih
        stop
      end if
c      
      a  = 1.0d+00 
      b  = 0.0d+00
      call dgemm('n','n',N,N2,N,a,O,N,V(1,N1+1),N,b,AUX,N) ! O V2
c
c      T
c    U1   O V2
c
      call dgemm('t','n',N1,N2,N,a,U,N,AUX,n,b,O20,N1)
c
      call dgemm('n','n',N,N2,N,a,O,N,U(1,N1+1),N,b,AUX,N) ! O U2
c
c      T             T
c    U1  O V2 - ht V1  O  U2
c
      a = -Dfloat(it*ih)
      b = 1.0d+00
      call dgemm('t','n',N1,N2,N,a,V,N,AUX,n,b,O20,N1)
c
      RETURN
      END
c +-------------------------------------------------------------------+
c |                                         _                         |
c |             _20              20   /  0  N \                       |
c |    Computes N     such that N  = !   _T   !                       |
c |                                   \ -N  0 /                       |
c +-------------------------------------------------------------------+
      Subroutine SN20(U,V,N1,N2,N,AN20)
c
      Implicit Real*8 (A-H,O-Z)
c
      Dimension U(N,*),V(N,*)
      Dimension AN20(N1,N2)
C      
c      T
c    U1   V2
c
      a   = 1.0d+00
      b   = 0.0d+00
      call dgemm('t','n',N1,N2,N,a,U,N,V(1,N1+1),n,b,AN20,N1)
c
c
c      T        T
c    U1  V2 - V1  U2
c
      a = -1.0d+00
      b =  1.0d+00
      call dgemm('t','n',N1,N2,N,a,V,N,U(1,N1+1),n,b,AN20,N1)
c
      RETURN
      END
c +-------------------------------------------------------------------+
c |                                                                   |
c | Computes the 20 part of an operator with the form respect to      |
c | signature                                                         |
c |                                                                   |
c |      /  0   O  \               /   O20(1)   0    \                |
c |     |          |              |                  |                |
c |      \ -tO  0 /                \     0    O20(2) /                |
c |                                                                   |
c |   IT = t                     T                                    |
c |   IH = hermiticity of O   ( O  = -ht T )                          |
c |                                                                   |
c |   U and V are stored as matrices of dimension N x 2N of the form  |
c |                                                                   |
c |      U=(U1,U2)      V=(V1,V2)                                     |
c |                                                                   |
c |   with U1 a matrix N x N1 and U2 a matrix N x N2  (N1+N2=2N)      |
c |                                                                   |
c +-------------------------------------------------------------------+
      Subroutine SO20M(U,V,N1,N2,N,O,AUX,O20_1,O20_2,IT,IH)
c
      Implicit Real*8 (A-H,O-Z)
c
      Dimension U(N,*),V(N,*)
      Dimension O(N,N),AUX(N,*)
      Dimension O20_1(N1,N1),O20_2(N2,N2)
C
      if(iabs(it*ih).ne.1) then
        write(6,'(" in SO20M IT IH ",2I3, " ******** ")') it,ih
        stop
      end if
c      
      a = 1.0d+00
      b = 0.0d+00
      call dgemm('n','n',N,N1,N,a,O,N,V,N,b,AUX,N) ! O V1
      a = -dfloat(IT)
      call dgemm('n','n',N,N2,N,a,O,N,V(1,N1+1),N,b,AUX(1,N1+1),N) ! O V2
c
c          T                  T
c        U1   O V1       -t U2  O V2
c
      a = 1.0d+00
      call dgemm('t','n',N1,N1,N,a,U,N,AUX,n,b,O20_1,N1)
      call dgemm('t','n',N2,N2,N,a,U(1,N1+1),N,AUX(1,N1+1),n,b,O20_2,N2)
c ------------------------------------------------------------------------
      a = dfloat(IT*IH)
      call dgemm('n','n',N,N1,N,a,O,N,U,N,b,AUX,N) ! O U1
      a = -dfloat(IH)
      call dgemm('n','n',N,N2,N,a,O,N,U(1,N1+1),N,b,AUX(1,N1+1),N) ! O U2
c
c      T	       T                 T  T	      T  T
c    U1  O  V1 + h t V1  O  U1	    -t U2  O  V2 -h V2  O  U2
c
      a = 1.0d+00
      b = 1.0d+00
      call dgemm('t','n',N1,N1,N,a,V,N,AUX,n,b,O20_1,N1)
      call dgemm('t','n',N2,N2,N,a,V(1,N1+1),N,AUX(1,N1+1),n,b,O20_2,N2)
c
      RETURN
      END
c +-------------------------------------------------------------------+
c |                                                                   |
c | Computes the 11 part of an operator with the form respect to      |
c | signature                                                         |
c |                                                                   |
c |                                                                   |
c |      /  O   0  \               / O11(1)     0   \                 |
c |     |          |              |                 |                 |
c |      \  0 t O  /               \   0     O11(2) /                 |
c |                                                                   |
c |                                                                   |
c |   IT = t                                                          |
c |   IH = hermiticity of O                                           |
c |                                                                   |
c |   U and V are stored as matrices of dimension N x 2N of the form  |
c |                                                                   |
c |      U=(U1,U2)      V=(V1,V2)                                     |
c |                                                                   |
c |   with U1 a matrix N x N1 and U2 a matrix N x N2  (N1+N2=2N)      |
c |                                                                   |
c +-------------------------------------------------------------------+
      Subroutine SO11P(U,V,N1,N2,N,O,AUX,O11_1,O11_2,IT,IH)
c
      Implicit Real*8 (A-H,O-Z)
c
      Dimension U(N,*),V(N,*)
      Dimension O(N,N),AUX(N,*)
      Dimension O11_1(N1,N1),O11_2(N2,N2)
C
      if(iabs(it*ih).ne.1) then
        write(6,'(" in SO11P IT IH ",2I3, " ******** ")') it,ih
        stop
      end if
c      
      a = -dfloat(IT*IH)
      b = 0.0d+00
      call dgemm('n','n',N,N1,N,a,O,N,V,N,b,AUX,N) ! O V1
      a = -1.0d+00
      call dgemm('t','n',N,N2,N,a,O,N,V(1,N1+1),N,b,AUX(1,N1+1),N) ! OT V2
c
c          T                  T
c  -h t  V1   O V1       -h V2  O V2
c
      a = 1.0d+00
      call dgemm('t','n',N1,N1,N,a,V,N,AUX,n,b,O11_1,N1)
      call dgemm('t','n',N2,N2,N,a,V(1,N1+1),N,AUX(1,N1+1),n,b,O11_2,N2)
c ------------------------------------------------------------------------
      a = 1.0d+00
      call dgemm('n','n',N,N1,N,a,O,N,U,N,b,AUX,N) ! O U1
      a = dfloat(ih*IT)
      call dgemm('t','n',N,N2,N,a,O,N,U(1,N1+1),N,b,AUX(1,N1+1),N) ! O U2
c
c      T	       T                 T  T	      T  T
c    U1  O  U1 - h t V1  O  V1	    ht U2  O  U2 -  V2  O  V2
c
      a = 1.0d+00
      b = 1.0d+00
      call dgemm('t','n',N1,N1,N,a,U,N,AUX,n,b,O11_1,N1)
      call dgemm('t','n',N2,N2,N,a,U(1,N1+1),N,AUX(1,N1+1),n,b,O11_2,N2)
c
      RETURN
      END
c +-------------------------------------------------------------------+
c |             11                                                    |
c |   Computes N                                                      |
c |                                                                   |
c +-------------------------------------------------------------------+
      Subroutine SN11(U,V,N1,N2,N,AN11_1,AN11_2)
c
      Implicit Real*8 (A-H,O-Z)
c
      Dimension U(N,*),V(N,*)
      Dimension AN11_1(N1,N1),AN11_2(N2,N2)
C      
c
c       T               T
c  -  V1   V1       - V2  V2
c
      a = -1.0d+00
      b =  0.0d+00
      call dgemm('t','n',N1,N1,N,a,V,N,V,n,b,AN11_1,N1)
      call dgemm('t','n',N2,N2,N,a,V(1,N1+1),N,V(1,N1+1),n,b,AN11_2,N2)
c ------------------------------------------------------------------------
c
c      T	 T             T        T
c    U1  U1 -  V1 V1	     U2  U2 - V2  V2
c
      a = 1.0d+00
      b = 1.0d+00
      call dgemm('t','n',N1,N1,N,a,U,N,U,n,b,AN11_1,N1)
      call dgemm('t','n',N2,N2,N,a,U(1,N1+1),N,U(1,N1+1),n,b,AN11_2,N2)
c
      RETURN
      END
c +-------------------------------------------------------------------+
c |                                                                   |
c | Computes the 11 part of an operator with the form respect to      |
c | signature                                                         |
c |                                                                   |
c |      /  0   O  \               /     0     O11   \                |
c |     |          |              |                  |                |
c |      \ -tO  0 /                \  h O11     0    /                |
c |                                                                   |
c |   IT = t                     T                                    |
c |   IH = hermiticity of O   ( O  = -ht T )                          |
c |                                                                   |
c |   U and V are stored as matrices of dimension N x 2N of the form  |
c |                                                                   |
c |      U=(U1,U2)      V=(V1,V2)                                     |
c |                                                                   |
c |   with U1 a matrix N x N1 and U2 a matrix N x N2  (N1+N2=2N)      |
c |                                                                   |
c +-------------------------------------------------------------------+
      Subroutine SO11M(U,V,N1,N2,N,O,AUX,O11,IT,IH)
c
      Implicit Real*8 (A-H,O-Z)
c
      Dimension U(N,*),V(N,*)
      Dimension O(N,N),AUX(N,N2)
      Dimension O11(N1,N2)
C
      if(iabs(it*ih).ne.1) then
        write(6,'(" in SO11M IT IH ",2I3, " ******** ")') it,ih
        stop
      end if
c      
      a  = dfloat(IT*IH)
      b  = 0.0d+00
      call dgemm('n','n',N,N2,N,a,O,N,V(1,N1+1),N,b,AUX,N) ! O V2
c
c      T
c  htV1   O V2
c
      a = 1.0d+00
      call dgemm('t','n',N1,N2,N,a,V,N,AUX,n,b,O11,N1)
c
      call dgemm('n','n',N,N2,N,a,O,N,U(1,N1+1),N,b,AUX,N) ! O U2
c
c      T             T
c ht V1  O V2 +    U1  O  U2
c
      a = 1.0d+00        
      b = 1.0d+00
      call dgemm('t','n',N1,N2,N,a,U,N,AUX,n,b,O11,N1)
c
      RETURN
      END
c +-------------------------------------------------------------------+
c |                                                                   |
c | Computes the 20 part of the hamiltonian                           |
c |                 					              |
c |                   /   0    H20 \		     		      |
c |                  |      T	   |		     		      |
c |                   \ -H20	0 /		     		      |
c |                                     	     	              |
c |                                                                   |
c |   U and V are stored as matrices of dimension N x 2N of the form  |
c |                                                                   |
c |      U=(U1,U2)      V=(V1,V2)                                     |
c |                                                                   |
c |   with U1 a matrix N x N1 and U2 a matrix N x N2  (N1+N2=2N)      |
c |                                                                   |
c |   AUX(N,2N2)                                                      |
c +-------------------------------------------------------------------+
      Subroutine SH20(U,V,N1,N2,N,G,D,AUX,H20)
c
      Implicit Real*8 (A-H,O-Z)
c      
      Dimension U(N,*),V(N,*)
      Dimension AUX(N,*)
      Dimension G(N,*),D(N,*)
      Dimension H20(N1,N2)
c
c ---------------------------------------------------------------------
c
      I1 = N1+1
      I2 = N2+1
      a = 1.0d+00
      b  = 0.0d+00
      call dgemm('n','n',N,N2,N,a,G,N,V(1,I1),N,b,AUX,N) ! h1 V2
c
      call dgemm('n','n',N,N2,N,a,G(1,N+1),N,U(1,I1),N,b,AUX(1,I2),N) ! h2 U2
c
c
      a = 1.0d+00
      b  = 0.0d+00
c
c      T
c    U1   h1 V2
c
      call dgemm('t','n',N1,N2,N,a,U,N,AUX,N,b,H20,N1)
c
c      T
c    V1   h2 U2
c
      a =-1.0d+00
      b = 1.0d+00
c
      call dgemm('t','n',N1,N2,N,a,V,N,AUX(1,I2),N,b,H20,N1)
c
      a = 1.0d+00
      b  = 0.0d+00
      call dgemm('n','n',N,N2,N,a,D,N,U(1,I1),N,b,AUX,N) ! d1 u2
c
      call dgemm('t','n',N,N2,N,a,D,N,V(1,I1),N,b,AUX(1,I2),N) ! d1T v2
c
      a = 1.0d+00
      b  = 1.0d+00
c
c      T
c    U1   d1 U2
c
      call dgemm('t','n',N1,N2,N,a,U,N,AUX,N,b,H20,N1)
c
c      T    T 
c    V1   D1 V2
c
      a = 1.0d+00
      b = 1.0d+00
c
      call dgemm('t','n',N1,N2,N,a,V,N,AUX(1,I2),N,b,H20,N1)
c
      return
      end
c +-------------------------------------------------------------------+
c |                                                                   |
c | Computes the 11 part of the hamiltonian                           |
c |                 					              |
c |                   /  H11_1    0    \	     		      |
c |                  |       	       |	     		      |
c |                   \   0  	H11_2 /		     		      |
c |                                     	     	              |
c |                                                                   |
c |   U and V are stored as matrices of dimension N x 2N of the form  |
c |                                                                   |
c |      U=(U1,U2)      V=(V1,V2)                                     |
c |                                                                   |
c |   with U1 a matrix N x N1 and U2 a matrix N x N2  (N1+N2=2N)      |
c |                                                                   |
c |  H11_1(N1,N1)   H11_2(N2,N2)                                      |
c |                                                                   |
c |   AUX(N,2N2)                                                      |
c +-------------------------------------------------------------------+
      Subroutine SH11(U,V,N1,N2,N,G,D,AUX,H11_1,H11_2)
c
      Implicit Real*8 (A-H,O-Z)
c      
      Dimension U(N,*),V(N,*)
      Dimension AUX(N,*)
      Dimension G(N,*),D(N,*)
      Dimension H11_1(N1,N1),H11_2(N2,N2)
c
c +-------------------------------------------------------------------+
c |                         11                                        |
c |                       H                                           |
c |                         1                                         | 
c +-------------------------------------------------------------------+
c
      I1 = N1+1
      I2 = N2+1
      
      a = 1.0d+00
      b  = 0.0d+00
      call dgemm('t','n',N,N1,N,a,G(1,N+1),N,V,N,b,AUX,N) ! h2T V1
c
      call dgemm('n','n',N,N1,N,a,G,N,U,N,b,AUX(1,I1),N) ! h1 U1
c
c
      a = 1.0d+00
      b  = 0.0d+00
c
c      T    
c    U1   h1 U1
c
      call dgemm('t','n',N1,N1,N,a,U,N,AUX(1,I1),N,b,H11_1,N1)
c
c      T    T
c   -V1   h2 V1
c
      a =-1.0d+00
      b = 1.0d+00
c
      call dgemm('t','n',N1,N1,N,a,V,N,AUX,N,b,H11_1,N1)
c
      a = 1.0d+00
      b  = 0.0d+00
      call dgemm('n','n',N,N1,N,a,D,N,V,N,b,AUX,N) ! d1 V1
c
      call dgemm('t','n',N,N1,N,a,D,N,U,N,b,AUX(1,I1),N) ! d1T U1
c
      a = 1.0d+00
      b  = 1.0d+00
c
c      T
c    U1   d1 V1
c
      call dgemm('t','n',N1,N1,N,a,U,N,AUX,N,b,H11_1,N1)
c
c      T    T
c    V1   D1 U1
c
      a = 1.0d+00
      b = 1.0d+00
c
      call dgemm('t','n',N1,N1,N,a,V,N,AUX(1,I1),N,b,H11_1,N1)
c +-------------------------------------------------------------------+
c |                         11                                        |
c |                       H                                           |
c |                         2                                         | 
c +-------------------------------------------------------------------+
c
      a = 1.0d+00
      b  = 0.0d+00
      call dgemm('t','n',N,N2,N,a,G,N,V(1,I1),N,b,AUX,N) ! h1T V2
c
      call dgemm('n','n',N,N2,N,a,G(1,N+1),N,U(1,I1),N,b,AUX(1,I2),N) ! h2 U2
c
c
      a = 1.0d+00
      b  = 0.0d+00
c
c      T    
c    U2   h2 U2
c
      call dgemm('t','n',N2,N2,N,a,U(1,I1),N,AUX(1,I2),N,b,H11_2,N2)
c
c      T    T
c   -V2   h1 V2
c
      a =-1.0d+00
      b = 1.0d+00
c
      call dgemm('t','n',N2,N2,N,a,V(1,I1),N,AUX,N,b,H11_2,N2)
c
      a  = -1.0d+00
      b  =  0.0d+00
      call dgemm('t','n',N,N2,N,a,D,N,V(1,I1),N,b,AUX,N) ! -d1T V2
c
      call dgemm('n','n',N,N2,N,a,D,N,U(1,I1),N,b,AUX(1,I2),N) ! -d1 U2
c
      a = 1.0d+00
      b  = 1.0d+00
c
c      T    T  
c  - U2   d1 V2
c
      call dgemm('t','n',N2,N2,N,a,U(1,I1),N,AUX,N,b,H11_2,N2)
c
c      T    
c  - V2   D1 U2
c
      a = 1.0d+00
      b = 1.0d+00
c
      call dgemm('t','n',N2,N2,N,a,V(1,I1),N,AUX(1,I2),N,b,H11_2,N2)
c
      return
      end
c +-------------------------------------------------------------------+
c |                                                                   |
c |    Subroutine          N E W U V                                  |
c |                                                                   |
c +-------------------------------------------------------------------+
      Subroutine NEWUV(U,V,N1,N2,N,Z,AUX)
c
      Implicit Real*8 (A-H,O-Z)
c      
      Dimension U(N,*),V(N,*)
      Dimension AUX(*)
      Dimension Z(N1,N2)
      
      NT12 = N*(N1+N2)
      call dcopy (NT12,U,1,AUX,1) ! U-> AUX
c
c   The U and V matrices are stored in:
c
c   U:   U1 U2
c   V:   V1 V2
c
c
c          / U1   0  \        / 0    V2 \
c      U = |         |    V = |         |
c          \ 0    U2 /        \ V1   0  /
c
c +-------------------------------------------------------------------+
c                                             T
c                               U1 = U1 - V2*Z
c +-------------------------------------------------------------------+
      a= -1.0D+00
      b=  1.0D+00
      call dgemm('N','T',N,N1,N2,a,V(1,N1+1),N,Z,N1,b,U,N)
c
c
c +-------------------------------------------------------------------+
c                               U2 = U2 + V1*Z
c +-------------------------------------------------------------------+
      a=  1.0D+00
      b=  1.0D+00
      call dgemm('N','N',N,N2,N1,a,V,N,Z,N1,b,U(1,N1+1),N)
c
c +-------------------------------------------------------------------+
c                               V2 = V2 + U1*Z
c +-------------------------------------------------------------------+
      a=  1.0D+00
      b=  1.0D+00
      call dgemm('N','N',N,N2,N1,a,AUX,N,Z,N1,b,V(1,N1+1),N)
c +-------------------------------------------------------------------+
c                                             T
c                               V1 = V1 - U2*Z
c +-------------------------------------------------------------------+
      a= -1.0D+00
      b=  1.0D+00
      call dgemm('N','T',N,N1,N2,a,AUX(N*N1+1),N,Z,N1,b,V,N)
      
      one = 1.0d+00
c+---------------------------------------------------------------------+
c|                                  T                                  |
c|                             1+Z*Z                                   |
c+---------------------------------------------------------------------+
      do j=1,N1
         do i=1,N1
	    AUX(i+(j-1)*N1) = 0.0d+00
         end do
	 AUX(j+(j-1)*N1) = 1.0d+00
      end do
      
      call DSYRK('L','N',N1,N2,one,Z,N1,one,AUX,N1)
      
c+---------------------------------------------------------------------+
c|              C H O L E S K Y    F A C T O R I Z A T I O N           |
c|                                                                     |
c|                    T      T                                         |
c|               1+Z*Z  = L L      L lower triangular matrix           |
c|                         1 1                                         |
c+---------------------------------------------------------------------+
         call DPOTRF('L',N1,AUX,N1,info)
c
         if(info.ne.0) then
           write(6,*) ' ****** NEWUV  N1   INFO NE 0 in DPOTRF ',info
           stop
         end if
c
c+---------------------------------------------------------------------+
c|                    -1 T                     -1 T                    |
c|           U -> U  L                V -> V  L                        |
c|            1    1  1                1    1  1                       |
c+---------------------------------------------------------------------+
c
         call DTRSM('R','L','T','N',N,N1,one,AUX,N1,U,N)
         call DTRSM('R','L','T','N',N,N1,one,AUX,N1,V,N)
c+---------------------------------------------------------------------+
c|                                T                                    |
c|                             1+Z*Z                                   |
c+---------------------------------------------------------------------+
      
      do j=1,N2
         do i=1,N2
	    AUX(i+(j-1)*N2) = 0.0d+00
         end do
	 AUX(j+(j-1)*N2) = 1.0d+00
      end do
      
      call DSYRK('L','T',N2,N1,one,Z,N1,one,AUX,N2)
      
c+---------------------------------------------------------------------+
c|              C H O L E S K Y    F A C T O R I Z A T I O N           |
c|                                                                     |
c|                  T        T                                         |
c|               1+Z*Z  = L L      L lower triangular matrix           |
c|                         2 2                                         |
c+---------------------------------------------------------------------+
         call DPOTRF('L',N2,AUX,N2,info)
c
         if(info.ne.0) then
           write(6,*) ' ****** NEWUV  N2   INFO NE 0 in DPOTRF ',info
           stop
         end if
c
c+---------------------------------------------------------------------+
c|                    -1 T                     -1 T                    |
c|           U -> U  L                V -> V  L                        |
c|            2    2  2                2    2  2                       |
c+---------------------------------------------------------------------+
c
         call DTRSM('R','L','T','N',N,N2,one,AUX,N2,U(1,N1+1),N)
         call DTRSM('R','L','T','N',N,N2,one,AUX,N2,V(1,N1+1),N)
      	 
      return
      end

c +-------------------------------------------------------------------+
c |                                                                   |
c |      		    B L O C K I N G			      |
c +-------------------------------------------------------------------+
c |								      |
c |   IB1 (NMAXB)   Cuasiparticulas a bloquear de signatura +i	      |
c |   IB2 (NMAXB)   Cuasiparticulas a bloquear de signatura -i        |
c |								      |
c +-------------------------------------------------------------------+
      Subroutine BLOCKING(U,V,N,UB,VB,N1,N2,IB1,IB2)
c
      Implicit Real*8 (A-H,O-Z)
c
      Parameter (NMAXB=4)
      
      Dimension U (N,*),V (N,*)
      Dimension UB(N,*),VB(N,*)
      Dimension IB1(NMAXB),IB2(NMAXB)
      
c
c Signatura + i     
c  
      N1 = N
      N2 = N
      do i=1,NMAXB
         if(IB1(i).ne.0) then
	    N1 = N1 - 1
	    N2 = N2 + 1
	 end if
         if(IB2(i).ne.0) then
	    N1 = N1 + 1
	    N2 = N2 - 1
	 end if
      end do
      
      i1  = 1
      i2  = 1
      
      do ib=1,NMAXB      
      
         i=IB1(ib)
	 if(i.ne.0) then
	 
            call dcopy(N,U(1,i),1,VB(1,N1+i2),1)
            call dcopy(N,V(1,i),1,UB(1,N1+i2),1)
	    i2 = i2 + 1
	    
	 end if
	 
         i=IB2(ib)
	 if(i.ne.0) then
	 
            call dcopy(N,U(1,N+i),1,VB(1,i1),1)
            call dcopy(N,V(1,N+i),1,UB(1,i1),1)
	    i1 = i1 + 1
	    
	 end if
	 
      end do
      
      IQP = 1
      
      do i=1,N
      
         if(i.ne.IB1(IQP)) then
	 
            call dcopy(N,U(1,i),1,UB(1,i1),1)
            call dcopy(N,V(1,i),1,VB(1,i1),1)
            i1 = i1 + 1
	    
	 else
	 
	    IQP=IQP+1
	    
         end if
	 
      end do
      
      IQP = 1
      
      do i=1,N
      
         if(i.ne.IB2(IQP)) then
	 
            call dcopy(N,U(1,N+i),1,UB(1,N1+i2),1)
            call dcopy(N,V(1,N+i),1,VB(1,N1+i2),1)
            i2 = i2 + 1
	    
	 else
	 
	    IQP=IQP+1
	    
         end if
	 
      end do
      
      RETURN
      END
c +-------------------------------------------------------------------+
c +-------------------------------------------------------------------+
      Double precision function RMMP(L,N1,N2,O1_20,O2_20,E_1,E_2)
c
      Implicit Real*8 (A-H,O-Z)
c
      Dimension O1_20(N1,N2),O2_20(N1,N2)
      Dimension E_1(N1),E_2(N2)
      
      sum = 0.0d+00
      
      if (L.eq.0) then 
         do j=1,N2
            do i=1,N1
 
      	       sum = sum + O1_20(i,j)*O2_20(i,j)
 
      	    end do
         end do
      else
         do j=1,N2
            do i=1,N1
	 
	       sum = sum + O1_20(i,j)*O2_20(i,j)/((E_1(i)+E_2(j))**L)
	    
	    end do
         end do
      end if
      
      rmmp = sum
      
      return
      end
c +-------------------------------------------------------------------+
c +-------------------------------------------------------------------+
      Double precision function RMMPE(L,N1,N2,O1_20,O2_20,E_1,E_2,eta)
c
      Implicit Real*8 (A-H,O-Z)
c
      Dimension O1_20(N1,N2),O2_20(N1,N2)
      Dimension E_1(N1),E_2(N2)
      
      sum = 0.0d+00
      
      if (L.eq.0) then 
         do j=1,N2
            do i=1,N1
 
      	       sum = sum + O1_20(i,j)*O2_20(i,j)
 
      	    end do
         end do
      else
         do j=1,N2
            do i=1,N1
	       xx=E_1(i)+E_2(j)
c	       if(xx.lt.0.0) write(6,*) ' ** ',xx
	       sum=sum+O1_20(i,j)*O2_20(i,j)/((xx+eta)**L)
	    
	    end do
         end do
      end if
      
      rmmpe = sum
      
      return
      end
c +-------------------------------------------------------------------+
c +-------------------------------------------------------------------+
      Double precision function RMM(N12,O1_20,O2_20,QUOT)
c
      Implicit Real*8 (A-H,O-Z)
c
      Dimension O1_20(*),O2_20(*)
      Dimension QUOT(*)
            
      sum = 0.0d+00
      
      do i=1,N12
	 
	 sum = sum + O1_20(i)*O2_20(i)*QUOT(i)
	    
      end do
      
      rmm = sum
      
      return
      end
c +-------------------------------------------------------------------+
c +-------------------------------------------------------------------+
      Subroutine QUOT2QPE(N1,N2,E_1,E_2,QUOT)
      Implicit Real*8 (A-H,O-Z)
c
      Dimension QUOT(N1,N2)
      Dimension E_1(N1),E_2(N2)
      
      cutoff = 1.0d+00
      
      do j=1,N2
         do i=1,N1
	 
	    e2qp = max(E_1(i)+E_2(j),cutoff)
	    QUOT(i,j) = 1.0d+00/e2qp
	    
	 end do
      end do
            
      return
      end
c +-------------------------------------------------------------------+
c +-------------------------------------------------------------------+
      Subroutine DBE(N1,N2,Z0,QUOT)
      Implicit Real*8 (A-H,O-Z)
c
      Dimension QUOT(N1,N2),Z0(N1,N2)
      
      do j=1,N2
         do i=1,N1
	 
	    Z0(i,j) = Z0(i,j)*QUOT(i,j)
	    
	 end do
      end do
            
      return
      end
c +-------------------------------------------------------------------+
c +-------------------------------------------------------------------+
      Subroutine ZETA(it,eta,H20,C20,Z0,IX,ETC)
      Implicit real*8 (A-H,O-Z)
      Parameter (NCM=06,NYV=20*NCM,NLAM=NCM-3) ! maximum number of constraints
      Dimension C20(*),Z0(*),H20(*)
      Dimension IX(NCM),ETC(NCM)
c
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC3/ NNC20(4,NCM),NNC11(4,NCM)
      
      isoi = it/3 + 1
      N   = ND (it)
      NN  = ND2(it)
      
      do k=1,NN
         Z0(k) = -eta*H20(k)
      end do
      
      aa = etc(isoi)
      Call DAXPY(NN,aa,C20(NNC20(it,isoi)),1,Z0,1)
      
      do k=3,NCM
         if(ix(k).ne.0) then
	   aa = etc(k)
           Call DAXPY(NN,aa,C20(NNC20(it,k)),1,Z0,1)
	 end if
      end do
      
      return
      end
c +-------------------------------------------------------------------+
c |   Computes                                                        |
c |            U1 D1    U2 D2    V1 D1   V2 D2                        |
c |                                                                   |
c |                                                                   |
c |   U and V are stored as matrices of dimension N x 2N of the form  |
c |                                                                   |
c |      U=(U1,U2)      V=(V1,V2)                                     |
c |                                                                   |
c |   with U1 a matrix N x N1 and U2 a matrix N x N2  (N1+N2=2N)      |
c |                                                                   |
c +-------------------------------------------------------------------+
      Subroutine UVTRD(U,V,N1,N2,N,D1,D2,AUX)
c
      Implicit Real*8 (A-H,O-Z)
c
      Dimension U(N,*),V(N,*)
      Dimension AUX(N,*)
      Dimension D1(N1,N1),D2(N2,N2)
C
      I1 = 1 + N1
      ICOPY = N*(N1+N2)
      
      a  = 1.0d+00 
      b  = 0.0d+00
      call dgemm('n','n',N,N1,N1,a,V      ,N,D1,N1,b,AUX      ,N) ! V1 D1
      call dgemm('n','n',N,N2,N2,a,V(1,I1),N,D2,N2,b,AUX(1,I1),N) ! V2 D2
      
      call dcopy(ICOPY,AUX,1,V,1)
      
      call dgemm('n','n',N,N1,N1,a,U      ,N,D1,N1,b,AUX      ,N) ! U1 D1
      call dgemm('n','n',N,N2,N2,a,U(1,I1),N,D2,N2,b,AUX(1,I1),N) ! U2 D2
      
      call dcopy(ICOPY,AUX,1,U,1)
c
      return
      end
c +-------------------------------------------------------------------+
c +-------------------------------------------------------------------+
      Double precision function trace(A,N)      
c
      Implicit Real*8 (A-H,O-Z)
      Dimension A(N,N)
      s = 0.0d+00
      do i=1,N
         s = s + A(i,i)
      end do
      trace = s
      return
      end	 

c +-------------------------------------------------------------------+
c +-------------------------------------------------------------------+
      Subroutine CPYDIAG(A,N,D)      
      Implicit Real*8 (A-H,O-Z)
      Dimension A(N,N),D(N)
      
      do i=1,N
         D(i) = A(i,i)
      end do
      
      return
      end
