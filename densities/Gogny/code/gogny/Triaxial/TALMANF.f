      Real*8 Function DT1(n1,n2,Nmu)
c+---------------------------------------------------------------------+
c|    This function   computes the one-dimensional T1     Coeficient   |
c|                                                                     |
c|                   T  (   N1     N2    Nmu )                         |
c|                    1                                                |
c+---------------------------------------------------------------------+
c
      Implicit real*8 (A-H,O-Z)
      Logical srul
      Parameter (NFAC=150)
      Common /Nfact/dfact(0:Nfac),ddfact(0:Nfac),Nfacm
      if(Nfac.ne.Nfacm) call inigf
c+--------------------------------------------------------------------+
c|    dfact(I)  = ln(   I!   )   |   ln(I!) = dfact (I)               |
c|    ddfact(I) = ln( (2I-1)!! ) |   ln(I!!)= ddfact((I+1)/2)         |
c+--------------------------------------------------------------------+
c
      Nmumin = Iabs(n1-n2)
      Nmumax = n1+n2
      
      srul = (Mod(n1+n2,2).eq.Mod(nmu,2)).and. ! Paridad
     &       (n1 .ge.0).and.                   !
     &       (n2 .ge.0).and.                   !
     &       (nmu.ge.0).and.                   !
     &       (nmu.le.Nmumax).and.              !  nmu =< Nmumax
     &       (nmu.ge.Nmumin)                   !  nmu >= Nmumin
     
      
      if(srul) then
      
         im1= (Nmumin+Nmu)/2
         im2= (Nmumax-Nmu)/2
         im3= (Nmu-Nmumin)/2
         tx = 0.5d+00*(dfact(n2)+dfact(n1)+dfact(Nmu))
     *              -(dfact(im1)+dfact(im2)+dfact(im3))
         DT1 = dexp(tx)
      else
         DT1 = 0.0d+00
      end if
      
      return
      end
      Real*8 Function DJ1BB(n1,n2,Nmu,amu)
c+---------------------------------------------------------------------+
c|    This function   computes the one-dimensional J1 BB  Coeficient   |
c|                                                                     |
c|                   J  (   N1     N2    Nmu )                         |
c|                    1                                                |
c+---------------------------------------------------------------------+
c
      Implicit real*8 (A-H,O-Z)
      Logical srul
      Parameter (NFAC=150,MZTOT=25)
      Common /Nfact /dfact(0:Nfac),ddfact(0:Nfac),Nfacm
      Common /Ngfact/g12(-Nfac:Nfac),fg12(-Nfac:Nfac)
      if(Nfac.ne.Nfacm) call inigf
      
      srul = (Mod(n1+n2,2).eq.Mod(nmu,2)).and.
     &       (n1.ge.0).and.
     &       (n2.ge.0).and.
     &       (nmu.ge.0)
      
      if(srul) then
      
        gam     = 0.5d+00*amu**2
        dlgam   = dlog(gam)
        G       = dlog(1.0d+00+gam)
        pi      = 4.0d+00*atan(1.d+00)
        fa      = amu/(dsqrt(2.0d+00*pi*(1.0d+00+gam))*pi)
c ---------------------------------------------------------------------
        is    = (n1+n2+nmu)/2
        dlfac  = g12(is-n1)+g12(is-n2)
     &         + 0.5d+00*(dfact(n1)+dfact(n2)-dfact(nmu))
        fac  = fa*dexp(dlfac-dfloat(is)*G)*fg12(is-n1)*fg12(is-n2)
	
        sum  = 0.0d+00
        do ir = 0,Min(n1,n2)
           ff	= g12(is-nmu-ir)-(dfact(ir)+dfact(n1-ir)+dfact(n2-ir))
           sum = sum + dexp(ff+dfloat(ir)*dlgam)*fg12(is-nmu-ir)
        end do
	
        DJ1BB =  sum*fac
      else
        DJ1BB = 0.0d+00
      end if
      
      return
      end	
      Real*8 Function DL1(n1,n2,Nmu)
c+---------------------------------------------------------------------+
c|    This function   computes the one-dimensional L1     Coeficient   |
c|                                                                     |
c|                   L  (   N1     N2    Nmu )                         |
c|                    1                                                |
c+---------------------------------------------------------------------+
c
      Implicit real*8 (A-H,O-Z)
      Logical srul
      Parameter (NFAC=150,MZTOT=25)
      Common /Nfact /dfact(0:Nfac),ddfact(0:Nfac),Nfacm
      Common /Ngfact/g12(-Nfac:Nfac),fg12(-Nfac:Nfac)
      if(Nfac.ne.Nfacm) call inigf
      
      srul = (Mod(n1+n2,2).eq.Mod(nmu,2)).and.
     &       (n1.ge.0).and.
     &       (n2.ge.0).and.
     &       (nmu.ge.0)
      
      if(srul) then
      
        pi    = 4.0d+00*atan(1.d+00)
        ff    = 1.0d+00/(dsqrt(2.0d+00)*pi**1.5d+00)
        is    = (n1+n2+nmu)/2
        ss    = g12(is-nmu)+g12(is-n1)+g12(is-n2)
     &        -0.5d+00*(dfact(nmu)+dfact(n1)+dfact(n2))
        phas = fg12(is-nmu)*fg12(is-n1)*fg12(is-n2)
	
        DL1=  ff*dexp(ss)*phas
      else
        DL1 = 0.0d+00
      end if
      
      return
      end	
      Real*8 Function DD1(n1,n2,Nmu)
c+---------------------------------------------------------------------+
c|    This function   computes the one-dimensional D1     Coeficient   |
c|                                                                     |
c|                   D  (   N1  |   N2    Nmu )                        |
c|                    1                                                |
c+---------------------------------------------------------------------+
c
      Implicit real*8 (A-H,O-Z)

      if(n1.ge.0.and.n2.ge.0.and.nmu.ge.0) then
         DD1 = dsqrt(Dfloat(2*  n1)  )*DL1(n1-1,n2,Nmu)-
     &         dsqrt(Dfloat(2*(n1+1)))*DL1(n1+1,n2,Nmu)
      else
         DD1 = 0.0d+00
      end if
      
      return
      end	
      Real*8 Function DK1(n1,n2,Nmu)
c+---------------------------------------------------------------------+
c|    This function   computes the one-dimensional K1     Coeficient   |
c|                                                                     |
c|                   K  (   N1    N2    Nmu )                          |
c|                    1                                                |
c+---------------------------------------------------------------------+
c
      Implicit real*8 (A-H,O-Z)

      if(n1.ge.0.and.n2.ge.0.and.nmu.ge.0) then
         DK1 =(-dsqrt(Dfloat(n1    ))*DL1(n1-1,n2  ,Nmu)
     &         +dsqrt(Dfloat((n1+1)))*DL1(n1+1,n2  ,Nmu)
     &         -dsqrt(Dfloat((n2+1)))*DL1(n1  ,n2+1,Nmu)
     &         +dsqrt(Dfloat(n2    ))*DL1(n1  ,n2-1,Nmu))/dsqrt(2.0d+00)
      else
         DK1 = 0.0d+00
      end if
      
      return
      end	
      Real*8 Function DG1(n1,n2,Nmu)
c+---------------------------------------------------------------------+
c|    This function   computes the one-dimensional G1     Coeficient   |
c|                                                                     |
c|                   G  (   N1    N2    Nmu )                          |
c|                    1                                                |
c+---------------------------------------------------------------------+
c
      Implicit real*8 (A-H,O-Z)

      if(n1.ge.0.and.n2.ge.0.and.nmu.ge.0) then
         DG1 =(-dsqrt(Dfloat(n1    ))*DT1(n1-1,n2  ,Nmu)
     &         +dsqrt(Dfloat((n1+1)))*DT1(n1+1,n2  ,Nmu)
     &         -dsqrt(Dfloat((n2+1)))*DT1(n1  ,n2+1,Nmu)
     &         +dsqrt(Dfloat(n2    ))*DT1(n1  ,n2-1,Nmu))/dsqrt(2.0d+00)
      else
         DG1 = 0.0d+00
      end if
      
      return
      end	
      Real*8 Function EE1(n1,n2,Nmu)
c+---------------------------------------------------------------------+
c|    This function   computes the one-dimensional E1     Coeficient   |
c|                                                                     |
c|       (eta/b)**2  E  (   N1  |   N2    Nmu )                        |
c|                    1                                                |
c+---------------------------------------------------------------------+
c
      Implicit real*8 (A-H,O-Z)

      if(n1.ge.0.and.n2.ge.0.and.nmu.ge.0) then
         EE1 = dsqrt(Dfloat(2*  n1)  )*DL1(n1-1,n2,Nmu)+
     &         dsqrt(Dfloat(2*(n1+1)))*DL1(n1+1,n2,Nmu)
      else
         EE1 = 0.0d+00
      end if
      
      return
      end	
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< INIGF
      Subroutine INIGF()
      Implicit Real*8 (A-H,O-Z)
      Parameter (NFAC=150)
      Common /Nfact   / dfact(0:NFAC),ddfact(0:NFAC),Nfacm
      Common /Ngfact  / g12(-NFAC:NFAC),fg12(-NFAC:NFAC)
c+--------------------------------------------------------------------+
c|    dfact(I)  = ln(   I!   )   |   ln(I!) = dfact (I)               |
c|    ddfact(I) = ln( (2I-1)!! ) |   ln(I!!)= ddfact((I+1)/2)         |
c+--------------------------------------------------------------------+
      Nfacm = Nfac
      dfact(0) = 0.0d+00
      dfact(1) = 0.0d+00
      ddfact(0)= 0.0d+00
      ddfact(1)= 0.0d+00
      do ifac =2,nfacm
           dfact (ifac) = dlog(dfloat(ifac))  + dfact (ifac-1)
           ddfact(ifac) = dlog(dfloat(2*ifac-1))+ ddfact(ifac-1)
      end do
c+--------------------------------------------------------------------+
c|    g12(I) = ln(Abs(Gamma(I+1/2)))   fg12(i) = Sign of Gamma(I+1/2) |
c+--------------------------------------------------------------------+
      pi      = 4.0d+00*atan(1.d+00)
      g12(0)  = 0.5d+00*dlog(pi)
      fg12(0) = 1.0d+00
      do ifac=1,nfacm
         g12( ifac) = dlog(dfloat(ifac)-0.5d+00)+g12(ifac-1)
         g12(-ifac) = dlog(pi) - g12(ifac)
         fg12( ifac) = 1.0d+00
         fg12(-ifac) = Dfloat(1-2*mod(ifac,2))
      end do
      return
      end
