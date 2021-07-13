c +---------------------------------------------------------------------+
c |  Subroutine MOMANG: Computes the matrix elements of                 |
c |                                                                     |
c |     J   iJ    J                                                     |
c |      x    y    z                                                    |
c |                                                                     |
c |  in the triaxial basis.                                             |
c |                                                                     |
c |  Input: bx, by bz ........ Oscillator lengths   ! Common            |
c |         N ................ Dimension                                |   
c |                                                                     |
c |  OUTPUT: Ji(N,N)           Matrix element of the choosen operator   |
c |                                                                     |
c +---------------------------------------------------------------------+
      Subroutine MOMANG(Tipo,AJE,AJO)
      Implicit real*8 (A-H,O-Z)
      Include 'MAXDIM'
      Character*1 Tipo
      Dimension SQ(NSQ)    ! Internal
      Dimension AJE(NP,NP),AJO(NM,NM)
c      
      Save SQ,ICALL
C --------------------------------------------------------- Commons
      Common /CQN/    IQMAX,IQNE(3,NPMAX),IQNO(3,NMMAX),IFILL
      Common /DIMS/   NP,NM,NROP,NROM
      Common /OSCLEN/ bx,by,bz
c
      if(IFILL.eq.0) call setQN()
      
      if(icall.ne.61060) then
         do  i=1,NSQ
            sq(i) = dsqrt(dfloat(i-1))
         end do 
         icall = 61060
      end if
c
      if (Tipo.eq."x".or.Tipo.eq."X") then
         Itipo = 1
      else if (Tipo.eq."y".or.Tipo.eq."Y") then
         Itipo = 2
      else if (Tipo.eq."z".or.Tipo.eq."Z") then
         Itipo = 3
      else
         write(6,*) ' Undefined type in MOMANG ',Tipo
	 stop
      end if
c      
      bxp = by/bz+bz/by
      bxm = by/bz-bz/by
      byp = bz/bx+bx/bz
      bym = bz/bx-bx/bz
      bzp = bx/by+by/bx
      bzm = bx/by-by/bx
c+---------------------------------------------------------------------+
c                  Even parity
c
c  --- Note that the quantum number given by iqne are the real one + 1
c+---------------------------------------------------------------------+
      do ia=1,NP
         nxa = iqne(1,ia)
         jnxa= 1-2*Mod(nxa,2)
         nya = iqne(2,ia)
         jnya= 1-2*Mod(nya,2)
         nza = iqne(3,ia)
         do ib=ia,NP
            aje(ia,ib) = 0.0d+00
            aje(ib,ia) = 0.0d+00
c
            nxb = iqne(1,ib)
            nyb = iqne(2,ib)
            nzb = iqne(3,ib)
	    
            ixs= (nxb-nxa)
            iys= (nyb-nya)
            izs= (nzb-nza)
	    
            ix = Iabs(ixs)
            iy = Iabs(iys)
            iz = Iabs(izs)
c------------------------------------------------------- jx, jy, jz
            if((ix.le.1).and.(iy.le.1).and.(iz.le.1)) then
c--------------------------------------------------------------- jx
            if(Itipo.eq.1.and.ixs.eq.0) then
            zes = 0.0d+00
            if((iys.eq. 1).and.(izs.eq.1)) zes= bxm*sq(nya+1)*sq(nza+1)
            if((iys.eq.-1).and.(izs.eq.1)) zes=-bxp*sq(nya)  *sq(nza+1)
            if((iys.eq.1).and.(izs.eq.-1)) zes=-bxp*sq(nya+1)*sq(nza)
            if((iys.eq.-1).and.(izs.eq.-1))zes= bxm*sq(nya)  *sq(nza)
            if((iys.eq.0).and.(izs.eq.0)) zes=jnxa
            aje(ia,ib) = 0.5d+00*zes
            aje(ib,ia) = 0.5d+00*zes
c	    write(6,*) ' ia .. ',ia,ib,zes
            end if
c------------------------------------------------------------- i*jy
            if(Itipo.eq.2.and.iys.eq.0) then
            zes = 0.0d+00
            if((izs.eq. 1).and.(ixs.eq.1)) zes= bym*sq(nza+1)*sq(nxa+1)
            if((izs.eq.-1).and.(ixs.eq.1)) zes= byp*sq(nza)  *sq(nxa+1)
            if((izs.eq.1).and.(ixs.eq.-1)) zes=-byp*sq(nza+1)*sq(nxa)
            if((izs.eq.-1).and.(ixs.eq.-1))zes=-bym*sq(nza)  *sq(nxa)
            if((izs.eq.0).and.(ixs.eq.0)) zes=jnxa
            aje(ia,ib) = 0.5d+00*zes*jnya*jnxa
            aje(ib,ia) = 0.5d+00*zes*jnya*jnxa
            end if
c--------------------------------------------------------------- jz
            if(Itipo.eq.3.and.izs.eq.0) then
            zes = 0.0d+00
            if((ixs.eq. 1).and.(iys.eq.1)) zes=bzm*sq(nxa+1)*sq(nya+1)
            if((ixs.eq.-1).and.(iys.eq.1)) zes=bzp*sq(nxa)  *sq(nya+1)
            if((ixs.eq.1).and.(iys.eq.-1)) zes=bzp*sq(nxa+1)*sq(nya)
            if((ixs.eq.-1).and.(iys.eq.-1))zes=bzm*sq(nxa)  *sq(nya)
            if((ixs.eq.0).and.(iys.eq.0)) zes=1.d+00
            aje(ia,ib) =-0.5d+00*zes*jnya*jnxa
            aje(ib,ia) =-0.5d+00*zes*jnya*jnxa
            end if
c
            end if ! ix = 1 ,etc condition
         end do    ! ib loop
      end do       ! ia loop
c
c+---------------------------------------------------------------------+
c                  Odd parity
c
c  --- Note that the quantum number given by iqne are the real one + 1
c+---------------------------------------------------------------------+
      do ia=1,NM
         nxa = iqno(1,ia)
         jnxa= 1-2*Mod(nxa,2)
         nya = iqno(2,ia)
         jnya= 1-2*Mod(nya,2)
         nza = iqno(3,ia)
         do ib=ia,NM
            ajo(ia,ib) = 0.0d+00
            ajo(ib,ia) = 0.0d+00
c
            nxb = iqno(1,ib)
            nyb = iqno(2,ib)
            nzb = iqno(3,ib)
	    
            ixs= (nxb-nxa)
            iys= (nyb-nya)
            izs= (nzb-nza)
	    
            ix = Iabs(ixs)
            iy = Iabs(iys)
            iz = Iabs(izs)
c------------------------------------------------------- jx, jy, jz
            if((ix.le.1).and.(iy.le.1).and.(iz.le.1)) then
c--------------------------------------------------------------- jx
            if(Itipo.eq.1.and.ixs.eq.0) then
            zes = 0.0d+00
            if((iys.eq. 1).and.(izs.eq.1)) zes= bxm*sq(nya+1)*sq(nza+1)
            if((iys.eq.-1).and.(izs.eq.1)) zes=-bxp*sq(nya)  *sq(nza+1)
            if((iys.eq.1).and.(izs.eq.-1)) zes=-bxp*sq(nya+1)*sq(nza)
            if((iys.eq.-1).and.(izs.eq.-1))zes= bxm*sq(nya)  *sq(nza)
            if((iys.eq.0).and.(izs.eq.0)) zes=jnxa
            ajo(ia,ib) = 0.5d+00*zes
            ajo(ib,ia) = 0.5d+00*zes
c	    write(6,*) ' ia .. ',ia,ib,zes
            end if
c------------------------------------------------------------- i*jy
            if(Itipo.eq.2.and.iys.eq.0) then
            zes = 0.0d+00
            if((izs.eq. 1).and.(ixs.eq.1)) zes= bym*sq(nza+1)*sq(nxa+1)
            if((izs.eq.-1).and.(ixs.eq.1)) zes= byp*sq(nza)  *sq(nxa+1)
            if((izs.eq.1).and.(ixs.eq.-1)) zes=-byp*sq(nza+1)*sq(nxa)
            if((izs.eq.-1).and.(ixs.eq.-1))zes=-bym*sq(nza)  *sq(nxa)
            if((izs.eq.0).and.(ixs.eq.0)) zes=jnxa
            ajo(ia,ib) = 0.5d+00*zes*jnya*jnxa
            ajo(ib,ia) = 0.5d+00*zes*jnya*jnxa
            end if
c--------------------------------------------------------------- jz
            if(Itipo.eq.3.and.izs.eq.0) then
            zes = 0.0d+00
            if((ixs.eq. 1).and.(iys.eq.1)) zes=bzm*sq(nxa+1)*sq(nya+1)
            if((ixs.eq.-1).and.(iys.eq.1)) zes=bzp*sq(nxa)  *sq(nya+1)
            if((ixs.eq.1).and.(iys.eq.-1)) zes=bzp*sq(nxa+1)*sq(nya)
            if((ixs.eq.-1).and.(iys.eq.-1))zes=bzm*sq(nxa)  *sq(nya)
            if((ixs.eq.0).and.(iys.eq.0)) zes=1.d+00
            ajo(ia,ib) =-0.5d+00*zes*jnya*jnxa
            ajo(ib,ia) =-0.5d+00*zes*jnya*jnxa
            end if
c
            end if ! ix = 1 ,etc condition
         end do    ! ib loop
      end do       ! ia loop
      return
      end
c +---------------------------------------------------------------------+
c |  Subroutine ESPIN: Computes the matrix elements of                  |
c |                                                                     |
c |     S   iS    S                                                     |
c |      x    y    z                                                    |
c |                                                                     |
c |  in the triaxial basis.                                             |
c |                                                                     |
c |  Input: bx, by bz ........ Oscillator lengths   ! Common            |
c |         N ................ Dimension                                |   
c |                                                                     |
c |  OUTPUT: Ji(N,N)           Matrix element of the choosen operator   |
c |                                                                     |
c +---------------------------------------------------------------------+
      Subroutine ESPIN(Tipo,AJE,AJO)
      Implicit real*8 (A-H,O-Z)
      Include 'MAXDIM'
      Character*1 Tipo
      Dimension AJE(NP,NP),AJO(NM,NM)
c      
C --------------------------------------------------------- Commons
      Common /CQN/    IQMAX,IQNE(3,NPMAX),IQNO(3,NMMAX),IFILL
      Common /DIMS/   NP,NM,NROP,NROM
      Common /OSCLEN/ bx,by,bz
c
      if(IFILL.eq.0) call setQN()
      
c
      if (Tipo.eq."x".or.Tipo.eq."X") then
         Itipo = 1
      else if (Tipo.eq."y".or.Tipo.eq."Y") then
         Itipo = 2
      else if (Tipo.eq."z".or.Tipo.eq."Z") then
         Itipo = 3
      else
         write(6,*) ' Undefined type in ESPIN ',Tipo
	 stop
      end if
c      
c+---------------------------------------------------------------------+
c                  Even parity
c
c  --- Note that the quantum number given by iqne are the real one + 1
c+---------------------------------------------------------------------+
      do ia=1,NP
         nxa = iqne(1,ia)
         jnxa= 1-2*Mod(nxa,2)
         nya = iqne(2,ia)
         jnya= 1-2*Mod(nya,2)
         do ib=ia,NP
            aje(ia,ib) = 0.0d+00
            aje(ib,ia) = 0.0d+00
	 end do
c
c --------------------------------------------------------------- Sx
         if(Itipo.eq.1) then
            aje(ia,ia) = 0.5d+00*jnxa
         end if
c ------------------------------------------------------------- i*Sy
         if(Itipo.eq.2) then
            aje(ia,ia) = 0.5d+00*jnya
         end if
c --------------------------------------------------------------- Sz
         if(Itipo.eq.3) then
            aje(ia,ia) =-0.5d+00*jnya*jnxa
         end if
c
      end do       ! ia loop
c
c+---------------------------------------------------------------------+
c                  Odd parity
c
c  --- Note that the quantum number given by iqne are the real one + 1
c+---------------------------------------------------------------------+
      do ia=1,NM
         nxa = iqno(1,ia)
         jnxa= 1-2*Mod(nxa,2)
         nya = iqno(2,ia)
         jnya= 1-2*Mod(nya,2)
         do ib=ia,NM
            ajo(ia,ib) = 0.0d+00
            ajo(ib,ia) = 0.0d+00
	 end do
c
c --------------------------------------------------------------- Sx
         if(Itipo.eq.1) then
            ajo(ia,ia) = 0.5d+00*jnxa
         end if
c ------------------------------------------------------------- i*Sy
         if(Itipo.eq.2) then
            ajo(ia,ia) = 0.5d+00*jnya
         end if
c --------------------------------------------------------------- Sz
         if(Itipo.eq.3) then
            ajo(ia,ia) =-0.5d+00*jnya*jnxa
         end if
c
      end do       ! ia loop
      return
      end
c+---------------------------------------------------------------------+
c|  Subroutine KIN : Computes the matrix elements of                   |
c|                                                                     |
c|                    T                                                |
c|                                                                     |
c|                                                                     |
c|  in the triaxial basis.                                             |
c|                                                                     |
c|  Input: bx, by bz ........ Oscillator lengths   ! Common            |
c|         N ................ Dimension of the matrix                  |  
c|                                                                     |
c|  OUTPUT: TE(N,N)                                                    |
c|                                                                     |
c|  Internal DZ2 ............ One dimensional matrix elements of 2nd   |
c|                           derivative.                               |
c+---------------------------------------------------------------------+
      Subroutine Kin(TE,TO)
      Implicit real*8 (A-H,O-Z)
      Include 'MAXDIM'
      
      Dimension DZ2(NSQ,NSQ)   ! Internal
      
      Dimension TE(NP,NP),TO(NM,NM)
c      
      Save DZ2,ICALL
C --------------------------------------------------------- Commons
      Common /CQN/    IQMAX,IQNE(3,NPMAX),IQNO(3,NMMAX),IFILL
      Common /DIMS/   NP,NM,NROP,NROM
      Common /OSCLEN/ bx,by,bz
c
      if(IFILL.eq.0) call setQN()
c ------------------------------------------------- 
      if(icall.ne.61060) then
      do  i=1,NSQ
         do  j=i,NSQ
           dz2(i,j) = 0.0d+00
           if (i.eq.j) then
              dz2(i,j) = 0.5d+00 -dfloat(i)
           else if (j.eq.i+2) then
              dz2(i,j) = 0.5d+00*dsqrt(dfloat(i)*dfloat(i+1))
           end if
           dz2(j,i) = dz2(i,j)
         end do   ! m loop
      end do      ! n loop
      icall = 61060
      end if
c
      q2 = (bx/by)**2
      p2 = (bx/bz)**2
      !ct= -20.734863d+00/(bx*bx)
      ct= -20.73667622931578756154d+00/(bx*bx)
c+---------------------------------------------------------------------+
c                  Even parity
c
c  --- Note that the quantum number given by iqne are the real one + 1
c+---------------------------------------------------------------------+
      do ia=1,NP
         nxa = iqne(1,ia)
         nya = iqne(2,ia)
         nza = iqne(3,ia)
         do ib=ia,NP
             te(ia,ib) = 0.0d+00
             te(ib,ia) = 0.0d+00
c
            nxb = iqne(1,ib)
            nyb = iqne(2,ib)
            nzb = iqne(3,ib)
	    
            ix = Iabs(nxb-nxa)
            iy = Iabs(nyb-nya)
            iz = Iabs(nzb-nza)
            if((ix.le.2).and.(iy.le.2).and.(iz.le.2)) then
c --------------------------------------------------- kinetic energy
            yf = dfloat(1 - iy)
            ts = 0.0d+00
            if((iy.eq.0).and.(iz.eq.0)) ts=           dz2(nxa,nxb)
            if((ix.eq.0).and.(iz.eq.0)) ts=ts + yf*q2*dz2(nya,nyb)
            if((ix.eq.0).and.(iy.eq.0)) ts=ts +    p2*dz2(nza,nzb)
            te(ia,ib) = ct*ts
            te(ib,ia) = ct*ts
            end if ! ix = 2 ,etc condition
         end do    ! ib loop
      end do       ! ia loop
c
c     Call printd(te,  np,np,np,'   T E ')
c
c+---------------------------------------------------------------------+
c                  Odd parity
c
c  --- Note that the quantum number given by iqne are the real one + 1
c+---------------------------------------------------------------------+
      do ia=1,NM
         nxa = iqno(1,ia)
         nya = iqno(2,ia)
         nza = iqno(3,ia)
         do ib=ia,NM
             to(ia,ib) = 0.0d+00
             to(ib,ia) = 0.0d+00
c
            nxb = iqno(1,ib)
            nyb = iqno(2,ib)
            nzb = iqno(3,ib)
	    
            ix = Iabs(nxb-nxa)
            iy = Iabs(nyb-nya)
            iz = Iabs(nzb-nza)
            if((ix.le.2).and.(iy.le.2).and.(iz.le.2)) then
c --------------------------------------------------- kinetic energy
            yf = dfloat(1 - iy)
            ts = 0.0d+00
            if((iy.eq.0).and.(iz.eq.0)) ts=            dz2(nxa,nxb)
            if((ix.eq.0).and.(iz.eq.0)) ts=ts + yf*q2*dz2(nya,nyb)
            if((ix.eq.0).and.(iy.eq.0)) ts=ts +    p2*dz2(nza,nzb)
            to(ia,ib) = ct*ts
            to(ib,ia) = ct*ts
            end if ! ix = 2 ,etc condition
         end do    ! ib loop
      end do       ! ia loop
c
c     Call printd(to,  nm,nm,nm,'   T O ')
      return
      end
c +-------------------------------------------------------------------+
c |                                                                   |
c |  Subroutine QLMME : Computes the matrix elements of               |
c |                                                                   |
c |                  1    /           m          \                    |
c |           Q   =------|  M   + (-1)  r   M     |                   |
c |            lm   fm    \  lm          m   l-m /                    |
c |                                                                   |
c |   in the triaxial basis.                                          |
c |                                     |  1   m ge 0                 |
c |    fm = sqrt(2)   m ne 0       r  =<                              |
c |    fm = 2         m = 0         m   | -1   m lt 0                 |
c |                                                                   |
c |                                                                   |
c |                                                                   |
c | Input: bx, by bz .... Oscillator lengths   ! Common               |
c |        ndime ndimo .. Dimensions of the parity even and odd parts |
c |                                                                   |
c | OUTPUT: QLME  Matrix elements for even states                     |
c |         QLMO  Matrix elements for odd  states                     |
c |         IS    Signature                                           |
c |                                                                   |
c |    If QLM is a negative parity operator then the corresponding    |
c |    matrix element is in QLME                                      |
c |                                                                   |
c |                                                               n   |
c | Internal: ZN  ............ One dimensional matrix elements of z   |
c |                                                                   |
c +-------------------------------------------------------------------+
c
      Subroutine QLMME(l,m,QLME,QLMO,IS)
      Implicit real*8 (A-H,O-Z)
      Include 'MAXDIM'
      Dimension QLME(NP,*), QLMO(NM,*)
c
      Dimension coef(NCOEFQ),ipow(3,NCOEFQ)
C --------------------------------------------------------- Commons
      Dimension ZN(KMAX1,KMAX1,MAXP)
      
      Save ZN,notcall
      
      Common /CQN/    IQMAX,IQNE(3,NPMAX),IQNO(3,NMMAX),IFILL
      Common /DIMS/   NP,NM,NROP,NROM
      Common /OSCLEN/ bx,by,bz
c
      Common /notc/ icall
c
      if(IFILL.eq.0) call setQN()
      
      if(l.gt.(maxp-1)) then 
         write(6,*) ' Increase MAXP in include MAXDIM (QLMME) to ',L
	 stop
      end if
      
      if(IFILL.eq.0) call setQN()
      
      if(IQMAX.gt.KMAX1) then
         write(6,*) ' Increase KMAX in include MAXDIM (QLMME) to ',IQMAX
	 stop
      end if
      
c
c
c     Call ZNS the first time the routine is invoqued
c
      if(notcall.ne.61060) then  ! quite a bad luck if you find this value
        Nz   = maxp              ! the first time ( it's my birthday !)
        Call ZNS (ZN,Nz)
        notcall = 61060
      end if
c
      ideb = 1
      call cqlm3(l,m,coef,ipow,NCOEFQ,Idim,ideb)
c
      if(m.ge.0) then
         is   = (-1)**m
         iad  = 0
      else
         is   =-(-1)**m
         iad  = 1
      end if
      ipar = 1-2*mod(l,2)
c
      if(ipar.eq.1) then       ! positive parity
c
c
c                  Even parity
c
c  --- Note that the quantum number given by iqne are the real one + 1
c
c
      do ia=1,NP
         nxa = iqne(1,ia)
         nya = iqne(2,ia)
         nza = iqne(3,ia)
         do ib=ia,NP
            qlme(ia,ib) = 0.0d+00
c
            nxb = iqne(1,ib)
            nyb = iqne(2,ib)
            nzb = iqne(3,ib)
c
            if(is.eq.1) then
               nybam = (nyb-nya+iad)/2
               fas   = dfloat((-1)**nybam)
            else
               nybap = (nyb+nya-iad)/2
               fas   = dfloat((-1)**(nybap+nxb))
            end if
c
            do id=1,idim
               iz = ipow(3,id) + 1
               iy = ipow(2,id) + 1
               ix = ipow(1,id) + 1
               qlme(ia,ib) = qlme(ia,ib) + coef(id)*fas*
     *         zn(nxa,nxb,ix)*zn(nya,nyb,iy)*zn(nza,nzb,iz)
            end do
         qlme(ib,ia) = qlme(ia,ib)
         end do    ! ib loop
      end do       ! ia loop
c
c                  Odd  parity
c
      do ia=1,NM
         nxa = iqno(1,ia)
         nya = iqno(2,ia)
         nza = iqno(3,ia)
         do ib=ia,NM
            qlmo(ia,ib) = 0.0d+00
c
            nxb = iqno(1,ib)
            nyb = iqno(2,ib)
            nzb = iqno(3,ib)
c
            if(is.eq.1) then
               nybam = (nyb-nya+iad)/2
               fas   = dfloat((-1)**nybam)
            else
               nybap = (nyb+nya-iad)/2
               fas   = dfloat((-1)**(nybap+nxb))
            end if
c
            do id=1,idim
               iz = ipow(3,id) + 1
               iy = ipow(2,id) + 1
               ix = ipow(1,id) + 1
               qlmo(ia,ib) = qlmo(ia,ib) + coef(id)*fas*
     *         zn(nxa,nxb,ix)*zn(nya,nyb,iy)*zn(nza,nzb,iz)
            end do
         qlmo(ib,ia) = qlmo(ia,ib)
         end do    ! ib loop
      end do       ! ia loop
c
      else if(ipar.eq.(-1)) then          ! negative parity operator
c
      do ia=1,NP
         nxa = iqne(1,ia)
         nya = iqne(2,ia)
         nza = iqne(3,ia)
         do ib=1,NM
            qlme(ia,ib) = 0.0d+00
c
            nxb = iqno(1,ib)
            nyb = iqno(2,ib)
            nzb = iqno(3,ib)
c
            if(is.eq.1) then
               nybam = (nyb-nya+iad)/2
               fas   = dfloat((-1)**nybam)
            else
               nybap = (nyb+nya-iad)/2
               fas   = dfloat((-1)**(nybap+nxb))
            end if
c
            do id=1,idim
               iz = ipow(3,id) + 1
               iy = ipow(2,id) + 1
               ix = ipow(1,id) + 1
               qlme(ia,ib) = qlme(ia,ib) + coef(id)*fas*
     *         zn(nxa,nxb,ix)*zn(nya,nyb,iy)*zn(nza,nzb,iz)
            end do
         end do    ! ib loop
      end do       ! ia loop
c
      end if    ! parity
      return
      end
c +-------------------------------------------------------------------+
C
C  To compute the coeficients of z y and x in Q      any m
C                                              lm
C
C
C      Ndim         Number of coeficients (output)
C      bx,by,bz     Oscillator lenghts
C
c +-------------------------------------------------------------------+
C    Last revision: 31 May 1992
c    Tested       : 10 june 1992
c +-------------------------------------------------------------------+
c
      Subroutine cqlm3(l,m,coef,ipow,ndim,idim,ideb)
      Implicit real*8 (A-H,O-Z)
      Include 'MAXDIM'
      Dimension coef(ndim),ipow(3,ndim)
      Common /fact/   dfact(nfac),ddfact(nfac),nfacm
      Common /OSCLEN/ bx,by,bz
c
      if(nfac.ne.nfacm) call setfact()
      
      if(ideb.ne.0) write(6,100) l,m
c
      if(m.ge.0) then         ! m>0
           imu = (-1)**m
           ms  = 1
      else                    ! m<0
           imu = 1
           ms  = 2
      end if
      
      im  = iabs(m)
      twol = dfloat(l)*dlog(2.0d+00)
      fact =imu*dexp(0.5d+00*(dfact(l-im+1)-dfact(l+im+1))-twol)
      
      if(m.eq.0) then
        fact = fact*0.5d+00
      else
        fact = fact*0.707106781d+00  ! 1/sqrt(2)
      end if
      
      nzin = Mod(l+im,2)+1
      nzfi = l+1
      idim = 0
      
      do nz=nzin,nzfi,2
         nyin = ms
         nyfi = l-nz + 2
         if(nyfi.ge.nyin) then
         do ny=nyin,nyfi,2
           nx = l-nz-ny+3
           idim  = idim + 1
           sum   = 0.0d+00
           nroin = (l-im-nz+1)/2 + 1
           if(nroin.le.0) nroin = 1
           nrofi = (l-im)/2 + 1
           is = (-1)**nroin
           do nro=nroin,nrofi
              is = -is
              ik = (-1)**ms
              kfin = Min(ny,(im+1))
              do k=ms,kfin,2
                 ik = - ik
                 if((nx+k-im).ge.2) then
        slog = dfact(2*l-2*nro+3)+dfact(im+1)
     *  -(dfact((nx+k-im)/2)+dfact(k)+dfact(im-k+2)+dfact(nro-nroin+1)
     *  +dfact((ny-k)/2+1)+dfact(l-nro+2)+dfact(l-im-2*nro+3))
                 sum = sum + dfloat(is*ik)*dexp(slog)
		 end if
              end do
           end do
	   
           if(dabs(sum).ge.1.d-10) then
	   
              coef(idim)=fact*sum*2.0d+00
              IZ = nz-1
              IY = ny-1
              IX = nx-1
              ipow(1,idim) = ix
              ipow(2,idim) = iy
              ipow(3,idim) = iz
c>deb
              if(ideb.ne.0) then
                  write(6,101) iz,iy,ix,coef(idim)
              end if
c>edeb
              coef(idim) = coef(idim)*bz**iz * by**iy * bx**ix
           else   ! of sum if
             idim = idim - 1
           end if
         end do
	 end if
      end do
      if(idim.gt.ndim) then
         write(6,*)' SUBROUTINE CQLM3 IDIM gt NDIM ',idim,ndim
         stop
      end if
  100 format ( ' Coeficients of z  y  x  in  M( ',i2,',',i2,')')
  101 format ( 26x,i2,4x,i2,4x,i2,/
     *  ,5x,d16.10,' * Z ',2x,' Y  ',2x,'  X ')
      return
      end
c+---------------------------------------------------------------------+
c|                                                                     |
c|                                                                     |
c|                                                   n                 |
c|  Subroutine ZNS : To compute matrix elements of  z                  |
c|                                                                     |
c+---------------------------------------------------------------------+
c|  Nzmax ...............  Max value of nz quantum number              |
c|  Nz ..................  Max value of the z power+1 (n=0 is also     |
c|                         included !)                                 |
c|                                                                   n |
c|  ZN (Nzmax,Nzmax,Nz)... Matrix containing the matrix elements of z  |
c+---------------------------------------------------------------------+
c|                                                                     |
c|  Dependencies:                                                      |
c|                                                                     |
c|      COMMON/CONST/         Subroutine SETUP                         |
c|      COMMON/FACT/          Subroutine SETUP                         |
c+---------------------------------------------------------------------+
c|                                                  Date:  28 May 92   |
c+---------------------------------------------------------------------+
      Subroutine ZNS (ZN,Nz)
      Implicit Real*8 (A-H,O-Z)
      Include 'MAXDIM'
      Dimension ZN(KMAX1,KMAX1,MAXP)
      Common /fact / dfact(Nfac),ddfact(Nfac),Nfacm
c
      if(nfac.ne.nfacm) call setfact()
c
c                                 0
c                               Z
c         
      do nza=1,KMAX1
         do nzb=1,KMAX1
            zn(nza,nzb,1) = 0.0d+00
         end do
         zn(nza,nza,1) = 1.0d+00
      end do
c                               n
c ---------------------------- z
      dln2 = dlog(2.0d+00)
      do inn=2,Nz
      
         in = inn - 1
         fa = dfact(in+1)-dfloat(in)*dln2
	 
         do n=1,kmax1
	 
            do m=n,kmax1
	    
                a= 0.0d+00
                imod = mod(n+m+in,2)
c ........ symmetry test
                if ((imod.eq.0).and.(m.le.(n+in))) then
                   imumin = m-n+1
                   imuma0 = m+n-1
                   imumax = min(imuma0,in+1)
                   if (imumax.ge.imumin) then
                       sum = 0.0d+00
                       do imu=imumin,imumax,2
                         i1= (imumin+imu)/2
                         i2= (imuma0-imu)/2 + 1
                         i3= (imu-imumin)/2 + 1
                         i4= (in-imu+1)/2   + 1
                         suml=fa+0.5d+00*(dfact(n)+dfact(m)+
     *                        (imu-1)*dln2)
     *                       -(dfact(i1)+dfact(i2)+dfact(i3)+dfact(i4))
                         sum = sum + dexp(suml)
                       end do  
                       a = sum
                   end if
                end if
                zn(n,m,inn) = a
                zn(m,n,inn) = a
            end do  
         end do  
      end do  
      return
      end
c
c    Matrix elements of x**l+y**l+z**l (L even)
c
      Subroutine RLME(l,RLE,RLO)
      Implicit real*8 (A-H,O-Z)
      Include 'MAXDIM'
      Dimension RLE(NP,*), RLO(NM,*)
c
C --------------------------------------------------------- Commons
      Dimension ZN(KMAX1,KMAX1,MAXP)
      
      Save ZN,notcall
      
      Common /CQN/    IQMAX,IQNE(3,NPMAX),IQNO(3,NMMAX),IFILL
      Common /DIMS/   NP,NM,NROP,NROM
      Common /OSCLEN/ bx,by,bz
c
      if(IFILL.eq.0) call setQN()
      
      if(mod(l,2).ne.0) then
         write(6,*) ' RLME called with odd L ',L
	 stop
      end if
      
      if(l.gt.(maxp-1)) then 
         write(6,*) ' Increase MAXP in include MAXDIM (RLME) to ',L
	 stop
      end if
      
      if(IFILL.eq.0) call setQN()
      
      if(IQMAX.gt.KMAX1) then
         write(6,*) ' Increase KMAX in include MAXDIM (RLME) to ',IQMAX
	 stop
      end if
      
c
c
c     Call ZNS the first time the routine is invoqued
c
      if(notcall.ne.61060) then  ! quite a bad luck if you find this value
        Nz   = maxp              ! the first time ( it's my birthday !)
        Call ZNS (ZN,Nz)
        notcall = 61060
      end if
c
c
c
c
c                  Even parity
c
c  --- Note that the quantum number given by iqne are the real one + 1
c
c
      do ia=1,NP
         nxa = iqne(1,ia)
         nya = iqne(2,ia)
         nza = iqne(3,ia)
         do ib=ia,NP
            rle(ia,ib) = 0.0d+00
c
            nxb = iqne(1,ib)
            nyb = iqne(2,ib)
            nzb = iqne(3,ib)
c
            fas   = dfloat(1-2*Mod(((nyb-nya)/2),2))
c
            rle(ia,ib) =  fas*(
     *      zn(nxa,nxb,l+1)*zn(nya,nyb,  1)*zn(nza,nzb,  1)*bx**l +
     *      zn(nxa,nxb,  1)*zn(nya,nyb,l+1)*zn(nza,nzb,  1)*by**l +
     *      zn(nxa,nxb,  1)*zn(nya,nyb,  1)*zn(nza,nzb,l+1)*bz**l )
            rle(ib,ia) = rle(ia,ib)
         end do    ! ib loop
      end do       ! ia loop
c
c                  Odd  parity
c
      do ia=1,NM
         nxa = iqno(1,ia)
         nya = iqno(2,ia)
         nza = iqno(3,ia)
         do ib=ia,NM
            rlo(ia,ib) = 0.0d+00
c
            nxb = iqno(1,ib)
            nyb = iqno(2,ib)
            nzb = iqno(3,ib)
c
            fas   = dfloat(1-2*Mod(((nyb-nya)/2),2))
c
            rlo(ia,ib) = fas*(
     *      zn(nxa,nxb,l+1)*zn(nya,nyb,  1)*zn(nza,nzb,  1)*bx**l +
     *      zn(nxa,nxb,  1)*zn(nya,nyb,l+1)*zn(nza,nzb,  1)*by**l +
     *      zn(nxa,nxb,  1)*zn(nya,nyb,  1)*zn(nza,nzb,l+1)*bz**l )
            rlo(ib,ia) = rlo(ia,ib)
         end do    ! ib loop
      end do       ! ia loop
c
      return
      end
c +-------------------------------------------------------------------+
c |                                                                   |
c |  Subroutine NECK : Computes the matrix elements of the neck       |
c |  operator                                                         |
c |                        - (z/a)**2                                 |
c |                      e                                            |
c |                                                                   |
c |                                                                   |
c |   in the triaxial basis.                                          |
c |                                                                   |
c | Input: bx, by bz .... Oscillator lengths   ! Common               |
c |        ndime ndimo .. Dimensions of the parity even and odd parts |
c |                                                                   |
c | Output: ANKE  Matrix elements for even states                     |
c |         ANKO  Matrix elements for odd  states                     |
c |                                                                   |
c |                                                         -(z/a)**2 |
c | Internal: ANECKZ .. One dimensional matrix elements of e          |
c |                                                                   |
c +-------------------------------------------------------------------+
c
      Subroutine NECKME(ANKE,ANKO)
      Implicit real*8 (A-H,O-Z)
      Include 'MAXDIM'
      Dimension ANKE(NP,*), ANKO(NM,*)
      Dimension Aneckz(kmax1,kmax1)         ! Internal
      
      Save Aneckz,notcall
c
      Common /CQN/    IQMAX,IQNE(3,NPMAX),IQNO(3,NMMAX),IFILL
      Common /DIMS/   NP,NM,NROP,NROM
      Common /OSCLEN/ bx,by,bz
      Common /neckco/ aneck,rneckp,rneckn
      
      if(IFILL.eq.0) call setQN()
      if(IQMAX.gt.KMAX1) then
         write(6,*) ' Increase KMAX in include MAXDIM (NECKME) to',IQMAX
	 stop
      end if
c
      if(notcall.ne.61060) then  ! quite a bad luck if you find this value
        nzmax = Kmax1
        bbz   = bz
        aa    = aneck    ! necking parameter
c
        Call neckz(nzmax,Aneckz,aa,bbz)
        notcall = 61060
      end if
c
c +-------------------------------------------------------------------+
c
c                  Even parity
c
c  --- Note that the quantum number given by iqne are the real one + 1
c
c +-------------------------------------------------------------------+
c
      do ia=1,NP
         nxa = iqne(1,ia)
         nya = iqne(2,ia)
         nza = iqne(3,ia)
         do ib=ia,NP
            nxb = iqne(1,ib)
            nyb = iqne(2,ib)
            nzb = iqne(3,ib)
c
            anke(ia,ib) = 0.0d+00
      if((nxa.eq.nxb).and.(nya.eq.nyb).and.(Mod(nza+nzb,2).eq.0)) then
            anke(ia,ib) = Aneckz(nza,nzb)
      end if
c
         anke(ib,ia) = anke(ia,ib)
         end do    ! ib loop
      end do       ! ia loop
c
c +-------------------------------------------------------------------+
c                  Odd  parity
c +-------------------------------------------------------------------+
c
      do ia=1,NM
         nxa = iqno(1,ia)
         nya = iqno(2,ia)
         nza = iqno(3,ia)
         do ib=ia,NM
            nxb = iqno(1,ib)
            nyb = iqno(2,ib)
            nzb = iqno(3,ib)
c
            anko(ia,ib) = 0.0d+00
      if((nxa.eq.nxb).and.(nya.eq.nyb).and.(Mod((nza+nzb),2).eq.0)) then
            anko(ia,ib) = Aneckz(nza,nzb)
      end if
c
         anko(ib,ia) = anko(ia,ib)
         end do    ! ib loop
      end do       ! ia loop
c
      return
      end
c+---------------------------------------------------------------------+
c|   Computes the matrix element of the necking operator               |
c|                                                                     |
c|           -(z/a)**2                                                 |
c|    <n|  e            |m>                                            |
c|                                                                     |
c|  where |n> is a harmonic oscillator wave function in the "z"        |
c|  direction                                                          |
c|                                                                     |
c+---------------------------------------------------------------------+
      Subroutine NeckZ(nzmax,Aneckz,aa,bz)
      Implicit real*8 (A-H,O-Z)
c ---------------------------------------------------<< Start    include
      Include 'MAXDIM'
c ---------------------------------------------------<< End      include
      Dimension Aneckz(nzmax,nzmax)
c
      Common /fact/dfact(Nfac),ddfact(Nfac),Nfacm
c
      if(nfac.ne.nfacm) call setfact()
c
      qq  = (bz/aa)**2
      qq1 = dlog(qq/(1.0d+00+qq))
      qqf = 1.0d+00/Dsqrt(1.0d+00+qq)
      dln2= dlog(2.0d+00)
c
      do in =1,nzmax
         do im=in,nzmax
            Aneckz(in,im) = 0.0d+00
            If(Mod(in+im,2).eq.0) then     ! selection rule
              imumin = im-in+1
              imumax = im+in-1
              sum = 0.0d+00
              do imu =imumin,imumax,2
                 im1 = (imumin+imu)/2
                 im2 = (imumax-imu)/2 + 1
                 im3 = (imu-imumin)/2 + 1
                 imu2= (imu-1)/2
       tx =
     * 0.5d+00*(dfact(in)+dfact(im)+Dfloat(imu-1)*(qq1-dln2))
     * + dfact(imu)
     * -(dfact(im1)+dfact(im2)+dfact(im3)+dfact(imu2+1))
                 rx = dexp(tx)*Dfloat(1-2*Mod(imu2,2))
                 sum = sum + rx
              end do
              Aneckz(in,im) = sum * qqf
            end if
            Aneckz(im,in) = Aneckz(in,im)
         end do
      end do
      return
      end
