c +-------------------------------------------------------------------+
c |                                                                   |
c |                                                                   |
c |   U and V are stored as matrices of dimension N x 2N of the form  |
c |                                                                   |
c |      U=(U1,U2)      V=(V1,V2)                                     |
c |                                                                   |
c |      RO=(RO1,RO2)                                                 |
c |                                                                   |
c |   with U1 a matrix N x N1 and U2 a matrix N x N2  (N1+N2=2N)      |
c |                                                                   |
c +-------------------------------------------------------------------+
      Subroutine RO(V,N1,N2,N,ROM)
c
      Implicit Real*8 (A-H,O-Z)
c
      Dimension V(N,*)
      Dimension ROM(N,*)
c      
      a = 1.0d+00
      b = 0.0d+00
      call dgemm('n','t',N,N,N2,a,V(1,N1+1),N,V(1,N1+1),N,b,ROM,N)
      call dgemm('n','t',N,N,N1,a,V,N,V,N,b,ROM(1,N+1),N)
c
      RETURN
      END
c +-------------------------------------------------------------------+
c |                                                                   |
c |                                                                   |
c |   U and V are stored as matrices of dimension N x 2N of the form  |
c |                                                                   |
c |      U=(U1,U2)      V=(V1,V2)                                     |
c |                                                                   |
c |   with U1 a matrix N x N1 and U2 a matrix N x N2  (N1+N2=2N)      |
c |                                                                   |
c +-------------------------------------------------------------------+
      Subroutine KAPPA(U,V,N1,N2,N,AKAP)
c
      Implicit Real*8 (A-H,O-Z)
c
      Dimension U(N,*),V(N,*)
      Dimension AKAP(N,N)
c      
      a = 1.0d+00
      b  = 0.0d+00
      call dgemm('n','t',N,N,N2,a,V(1,N1+1),N,U(1,N1+1),N,b,AKAP,N)
c
      RETURN
      END
