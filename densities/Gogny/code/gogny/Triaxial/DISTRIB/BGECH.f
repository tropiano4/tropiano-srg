c+---------------------------------------------------------------------+
c|    Exchange and Pairing  contribution coming from the Brink-Boecker |
c|    interaction.                                                     |
c|                                                                     |
c+---------------------------------------------------------------------+
      Subroutine GECH(icero,icerop,GROP,GROM,xtjx,xtjxp,xtjy,xtjz,
     * ixtj,ixtj2,iytj,iytj2,iztj,iztj2,
     * ZOPTI,YZOPTI,gp1p,gp2p,gn1p,gn2p,dP1P,dp2p,dn1p,dn2p,
     *gp1m,gp2m,gn1m,gn2m,dp1m,dp2m,dn1m,dn2m)
      Implicit real*8 (A-H,O-Z)
      Logical LMZE,LMZO
      Logical lx12,lxy12
      Include 'COMDIM'
c                                                     Pos parity
c
c    GROP: rosp, roap, rospw, roapw, akasp, akaap, akaspw, akaapw
c
c
      Dimension GP1P(ngp),gn1p(ngp),gp2p(ngp),gn2p(ngp)
      Dimension dP1P(ngp),dn1p(ngp),dp2p(ngp),dn2p(ngp)
      Dimension GROP(ngp8)
c                                                     Neg parity fields
c
c    GROP: rosm, roam, rosmw, roamw, akasm, akaam, akasmw, akaamw
c
      Dimension GP1m(ngm),gn1m(ngm),gp2m(ngm),gn2m(ngm)
      Dimension dP1m(ngm),dn1m(ngm),dp2m(ngm),dn2m(ngm)
      Dimension GROM(ngm8)
c
c  xtj terms
c
      Dimension xtjx (2,maxtjx)
      Dimension xtjxp(2,maxtjx)
      Dimension xtjy(2,maxtjy)
      Dimension xtjz(2,maxtjz)
c
      Dimension ixtj(nxmax,nxmax),ixtj2(nxmax,nxmax)
      Dimension iytj(nymax,nymax),iytj2(nymax,nymax)
      Dimension iztj(nzmax,nzmax),iztj2(nzmax,nzmax)
c+---------------------------------------------------------------------+
c
c  ZOPTI                  __
C                        \     BB                __
C      ZA(nz1,nz2)   =    >   I                  RO
C                        /__   nz1 nzq nzq' nz2    q q'
C                         __
C                        \     BB                __
C      ZB(nz1,nz2)   =    >   I                  RO
C                        /__   nz1 nzq'nzq  nz2    q q'
C
c      ZB(nz1,nz2) = ZA(nz2,nz1)
c
C    ZA[t][is][i]
c     where                   __         1   2         1   2
c           t type of density RO   1 = ro +ro    2 = ro -ro
c                                       (s)           (a)
c                                  3 = K         4 = K
c
c           is isospin P = protons   N = neutrons
c
c           i=1,2 BB interaction length
c
c   All this quantities are stored in ZOPTI(mazopti)
c
c           mazopti = 16*nzmax*nzmax
c
c                             2      2
c           mayopti = 16*nymax *nzmax
c+---------------------------------------------------------------------+
      Dimension ZOPTI(mazopti)
      Dimension YZOPTI(mayopti)
c
      Dimension rozo(8)
C --------------------------------------------------------- Common Blocks
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
      common /pulga/ideb
c+---------------------------------------------------------------------+
c        Starts calculation of the exchange and pairing
c
c---->    1 =K1 2=K2  K1 >= K2
c
c---->    e sufix stands for even parity
c---->    o sufix stands for odd  parity
c+---------------------------------------------------------------------+
c
c     Sufix 1 stands for
c
c         _ 1                _  1
c        |            and   /_\         k1 >= k2
c          k1,k2                k1,k2
c
c     Sufix 2 stands for
c
c         _ 2                _  1
c        |            and   /_\         k1 >= k2
c          k1,k2                k2,k1
c+---------------------------------------------------------------------+
      if(icero.eq.0) then     ! do I set fields to zero ?
      do igp = 1,ngp
         gp1p(igp) = 0.0d+00
         gp2p(igp) = 0.0d+00
         gn1p(igp) = 0.0d+00
         gn2p(igp) = 0.0d+00
      end do
c
      do igm = 1,ngm
         gp1m(igm) = 0.0d+00
         gp2m(igm) = 0.0d+00
         gn1m(igm) = 0.0d+00
         gn2m(igm) = 0.0d+00
      end do
      end if
      if(icerop.eq.0) then     ! do I set fields to zero ?
      do igp = 1,ngp
         dp1p(igp) = 0.0d+00
         dp2p(igp) = 0.0d+00
         dn1p(igp) = 0.0d+00
         dn2p(igp) = 0.0d+00
      end do
c
      do igm = 1,ngm
         dp1m(igm) = 0.0d+00
         dp2m(igm) = 0.0d+00
         dn1m(igm) = 0.0d+00
         dn2m(igm) = 0.0d+00
      end do
      end if
c
      irop = 0
      irom = 0
      do 1 nqx1=1,nxmax
        nqym1 = MY(nqx1)
        do 11 nqx2=1,nqx1
          icx  = Mod((nqx1+nqx2),2)                     ! type of field
          nqym2 = MY(nqx2)                               ! maximum of nqy
          lx12 = nqx1.eq.nqx2
c ------------------------------------------   setting yzopti to zero
          do iyzo =1,mayopti
             yzopti(iyzo) = 0.0d+00
          end do
c
          do 2 nqy1 = 1,nqym1
            nqxy1 = nqx1+nqy1
            nqzi1e= nzie(nqx1,nqy1)
            nqzi1o= nzio(nqx1,nqy1)
            nqzm1 = MZ(nqx1,nqy1)                      ! maximum of nqz1
c
            if(lx12) nqym2 = nqy1
            do 12 nqy2=1,nqym2
              icy   = Mod((nqy1+nqy2),2)                 ! type of field
              nqxy  = nqx2+nqy2
              nqzi2e= nzie(nqx2,nqy2)
              nqzi2o= nzio(nqx2,nqy2)
              nqzm2 = MZ(nqx2,nqy2)
              lxy12 = lx12.and.(nqy1.eq.nqy2)
c
c---->    Positive parity
c
c
c be careful not to move the if below, the iz1e=iy1e need to be done
c
              do izo =1,mazopti
                  zopti(izo) = 0.0d+00
              end do
c
              if(LMZE(nqx1,nqy1).and.LMZE(nqx2,nqy2)) then
c
              do 30 nqz1 =nqzi1e,nqzm1,2
                if(lxy12) nqzm2=nqz1
                do 31 nqz2 = nqzi2e,nqzm2,2
                  irop  = irop + 1
c
c                 Starts the calculation of Z (nz1,nz2) for each nyq,nxq,etc
c    NEUTRONS
c
c    Ro = Ro              nqz' = nqz1
c           nqz ,nqz'
c
                  rozo(1)= grop (irop)      ! ro (1) + ro (2) pos. parity
                  rozo(2)= grop (irop+1  )  ! ro (1) - ro (2= pos. parity
                  rozo(3)= grop (irop+4   ) ! kappa symmetric   "      "
                  rozo(4)= grop (irop+5   ) ! kappa antysymmetric
c    PROTONS
                  rozo(5)= grop (irop+2   ) ! ro (1) + ro (2) pos. parity
                  rozo(6)= grop (irop+3   ) ! ro (1) - ro (2= pos. parity
                  rozo(7)= grop (irop+6   ) ! kappa symmetric   "      "
                  rozo(8)= grop (irop+7   ) ! kappa antysymmetric
                  irop = irop + 7
c
                 if(lxy12.and.(nqz2.eq.nqz1)) then
                 do kk=1,8
                 rozo(kk) = rozo(kk)*0.5d+00
                 end do
                 end if
c ----------------------------------------------------- GECHZO in line
      call timeit(0,15,'GECH: GECHZO    ')
      imzq = nqz1+nqz2
      nnz1 = 0
      do nz1=1,nzmax
         jxtj = iztj(nz1,nqz1)  ! starting pos of xtj
         nz2i = 1+Mod((imzq+nz1-1),2)
	 nnz2 = (nz2i-1)*16
         do nz2 =nz2i,nzmax,2
            jx = jxtj + iztj2(nz2,nqz2)
            iz12 = nnz1 + nnz2 
            xxtj1 = xtjz(1,jx)
            xxtj2 = xtjz(2,jx)
	    do izo=1,8
	       zopti(iz12+izo  ) = zopti(iz12+izo  ) + xxtj1*rozo(izo)
	       zopti(iz12+izo+8) = zopti(iz12+izo+8) + xxtj2*rozo(izo)
	    end do 
c            call daxpy(8,xxtj1,rozo,1,zopti(iz12+1),1)
c	    call daxpy(8,xxtj2,rozo,1,zopti(iz12+9),1)
c	    call dger(8,2,1.0d+00,rozo,1,xtjz(1,jx),1,zopti(iz12+1),8)
            nnz2 = nnz2 + 32 
         end do
	 nnz1 = nnz1 + 16*nzmax
      end do
      call timeit(1,15,'GECH: GECHZO    ')
c                 call gechzo(iztj,iztj2,nqz1,nqz2,rozo,xtjz,zopti)
c ----------------------------------------------------- GECHZO in line

   31           continue
   30         continue
              end if
c
c---->    Negative parity
c
c
c be careful not to move the if below, the iz1o=iy1o need to be done
c
              if(LMZO(nqx1,nqy1).and.LMZO(nqx2,nqy2)) then
c
              do 40 nqz1 =nqzi1o,nqzm1,2
                if(lxy12) nqzm2=nqz1
                do 41 nqz2 = nqzi2o,nqzm2,2
                  irom  = irom  + 1
c
c                 Starts the calculation of Z (nz1,nz2) for each nyq,nxq,etc
c    NEUTRONS
c
c    Ro = Ro              nqz' = nqz1
c           nqz ,nqz'
c
                  rozo(1)= grom (irom)      ! ro (1) + ro (2) pos. parity
                  rozo(2)= grom (irom+1  )  ! ro (1) - ro (2= pos. parity
                  rozo(3)= grom (irom+4   ) ! kappa symmetric   "      "
                  rozo(4)= grom (irom+5   ) ! kappa antysymmetric
c    PROTONS
                  rozo(5)= grom (irom+2   ) ! ro (1) + ro (2) pos. parity
                  rozo(6)= grom (irom+3   ) ! ro (1) - ro (2= pos. parity
                  rozo(7)= grom (irom+6   ) ! kappa symmetric   "      "
                  rozo(8)= grom (irom+7   ) ! kappa antysymmetric
                  irom = irom + 7
c
c
                 if(lxy12.and.(nqz2.eq.nqz1)) then
                 do kk=1,8
                 rozo(kk) = rozo(kk)*0.5d+00
                 end do
                 end if
c ----------------------------------------------------- GECHZO in line
      call timeit(0,15,'GECH: GECHZO    ')
      imzq = nqz1+nqz2
      nnz1 = 0
      do nz1=1,nzmax
         jxtj = iztj(nz1,nqz1)  ! starting pos of xtj
         nz2i = 1+Mod((imzq+nz1-1),2)
	 nnz2 = (nz2i-1)*16
         do nz2 =nz2i,nzmax,2
            jx = jxtj + iztj2(nz2,nqz2)
            iz12 = nnz1 + nnz2 
            xxtj1 = xtjz(1,jx)
            xxtj2 = xtjz(2,jx)
	    do izo=1,8
	       zopti(iz12+izo  ) = zopti(iz12+izo  ) + xxtj1*rozo(izo)
	       zopti(iz12+izo+8) = zopti(iz12+izo+8) + xxtj2*rozo(izo)
	    end do 
c            call daxpy(8,xxtj1,rozo,1,zopti(iz12+1),1)
c	    call daxpy(8,xxtj2,rozo,1,zopti(iz12+9),1)
c	    call dger(8,2,1.0d+00,rozo,1,xtjz(1,jx),1,zopti(iz12+1),8)
            nnz2 = nnz2 + 32 
         end do
	 nnz1 = nnz1 + 16*nzmax
      end do
      call timeit(1,15,'GECH: GECHZO    ')
c ----------------------------------------------------- GECHZO in line
c                 call gechzo(iztj,iztj2,nqz1,nqz2,rozo,xtjz,zopti)
   41           continue
   40         continue
              end if
	      
              call gechyzo(iytj,iytj2,nqy1,nqy2,xtjy,zopti,yzopti)
   12       continue                                     ! end nqy  loop
    2     continue                                       ! end nqy1 loop

          Call gechxyz(ixtj,ixtj2,nqx1,nqx2,xtjx,
     -  xtjxp,yzopti,gp1p,gn1p,gp2p,gn2p,dp1p,dn1p,dp2p,
     -  dn2p,gp1m,gn1m,gp2m,gn2m,dp1m,dn1m,dp2m,dn2m)
c
   11   continue
    1 continue
      return
      end
c
c     Do the final calculation for the exchange and pairing fields
c
      Subroutine gechxyz(ixtj,ixtj2,nqx1,nqx,xtj
     -  ,xtjp,yzopti,gp1p,gn1p,gp2p,gn2p,dp1p,dn1p,dp2p,
     -  dn2p,gp1m,gn1m,gp2m,gn2m,dp1m,dn1m,dp2m,dn2m)
      Implicit Real*8 (A-H,O-Z)
      Logical lmze,lmzo
      Logical lx12,lxy12,lxyz12,lcx,lcy
      Include 'COMDIM'
      Dimension ixtj(nxmax,nxmax),ixtj2(nxmax,nxmax)
      Dimension yzopti(mayopti)
      Dimension xtj (2,maxtjx)
      Dimension xtjp(2,maxtjx)
c                                                     Pos parity fields
c     ngp size of pos parity symmetric matrices (one side)
      Dimension GP1P(ngp),gn1p(ngp),gp2p(ngp),gn2p(ngp)
      Dimension dP1P(ngp),dn1p(ngp),dp2p(ngp),dn2p(ngp)
c                                                     Neg parity fields
c     ngm size of neg parity symmetric matrices (one side)
      Dimension GP1m(ngm),gn1m(ngm),gp2m(ngm),gn2m(ngm)
      Dimension dP1m(ngm),dn1m(ngm),dp2m(ngm),dn2m(ngm)
c---------------------------------------------------- Internal vectors
      Dimension stf11(8),stf21(8),stf31(8),stf41(8)
      Dimension stf12(8),stf22(8),stf32(8),stf42(8)
c---------------------------------------------------------------------
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
c
      Common /ILOCECH/IXGP(ixmax,ixmax),IXGM(ixmax,ixmax)
      Common/GOGINT/amu(2),xw(2),xh(2),xb(2),xm(2),WLS,t3,alpha,x0,e2
c
c Handy combinations of the BB parameters.
c
      Common/GOGECH/c1(2),c2(2),p1(2),p2(2)
c
      call timeit(0,17,'GECH: GECHXYZ   ')
      c11 = c1(1)
      c12 = c1(2)
      c21 = c2(1)
      c22 = c2(2)
      p11 = p1(1)
      p12 = p1(2)
      p21 = p2(1)
      p22 = p2(2)
      h1 = xh(1)
      h2 = xh(2)
      w1 = xw(1)
      w2 = xw(2)
c ------------------------------------------------------------------
c---->                 1 >= 2
c---->    e sufix stands for even parity
c---->    o sufix stands for odd  parity
c ------------------------------------------------------------------
      nzsh = 16*nzmax*nzmax
      imxq = Mod((nqx1+nqx),2)
      
      do 1 nx1=1,nxmax
        nym1 = MY(nx1)
        jxtj = ixtj(nx1,nqx1)
        kxtj = ixtj(nx1,nqx )
        do 11 nx2=1,nx1
          imx12= Mod((nx1+nx2),2)
          lcx  = imx12.eq.0
          if(imxq.ne.imx12) go to 11       ! selection rule
c
          igp  = ixgp(nx1,nx2)
          igm  = ixgm(nx1,nx2)
          jx   = jxtj+ixtj2(nx2,nqx )
          kx   = kxtj+ixtj2(nx2,nqx1)
          xxtj1 = xtj (1,jx)
          xxtj2 = xtj (2,jx)
          xxtj1p= xtjp(1,jx)
          xxtj2p= xtjp(2,jx)
          xxtk1 = xtj (1,kx)
          xxtk2 = xtj (2,kx)
          xxtk1p= xtjp(1,kx)
          xxtk2p= xtjp(2,kx)
c
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          do 2 ny1 = 1,nym1
            nny1 = (ny1-1)*nymax
            nxy1 = nx1+ny1
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do 12 ny2=1,nym2
              nny2  = (ny2-1)*nymax
              imy12 = Mod((ny1+ny2),2)
              lcy   = imy12.eq.0
              iyz120= nzsh*(nny1+ny2-1)               ! standar
              jyz120= nzsh*(nny2+ny1-1)               ! q q' reversed
              nxy2 = nx2+ny2
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              lxy12 = lx12.and.(ny1.eq.ny2)
c
c---->    Positive parity
c
c
c be careful not to move the if below, the iz1e=iy1e need to be done
c
              if(lmze(nx1,ny1).and.lmze(nx2,ny2)) then
c
              do 30 nz1 =nzi1e,nzm1,2
                nnz1 = (nz1-1)*nzmax
                if(lxy12) nzm2=nz1
                do 31 nz2 = nzi2e,nzm2,2
                  nnz2   = (nz2-1)*nzmax
                  igp    = igp + 1
                  lxyz12 = lxy12.and.(nz1.eq.nz2)
                  feq    = 1.0d+00
c                  if(lxyz12) feq = 0.5d+00
                  iyz12  = iyz120 + 16*(nnz1+nz2-1)             ! pos in yzopti
                  jyz12  = jyz120 + 16*(nnz2+nz1-1)   ! q q' reversed
c
                  do 300 iro=1,8                               ! range 1
                  iyz12      = iyz12 + 1
                  jyz12      = jyz12 + 1
                  stf11(iro) = xxtj1 *yzopti(iyz12)   ! standar
                  stf21(iro) = xxtj1p*yzopti(iyz12)   ! tilde
                  stf31(iro) = xxtk1 *yzopti(jyz12)   ! q q' rev.
                  stf41(iro) = xxtk1p*yzopti(jyz12)   !
  300             continue
c
                  do 301 iro=1,8                              ! range 2
                  iyz12      = iyz12 + 1
                  jyz12      = jyz12 + 1
                  stf12(iro) = xxtj2 *yzopti(iyz12)   ! standar
                  stf22(iro) = xxtj2p*yzopti(iyz12)   ! tilde
                  stf32(iro) = xxtk2 *yzopti(jyz12)   ! q q' rev.
                  stf42(iro) = xxtk2p*yzopti(jyz12)   !
  301             continue
c -----------------------------------------------------------------------
c  stf indexes: 1: ro1+ro2  2: ro1-ro2  3: k sim  4:  k anti   NEUTRONS
c               5: ro1+ro2  6: ro1-ro2  7: k sim  8:  k anti   PROTONS
c
c   c1 = M + 0.5H   c2 = -(B+0.5W)     p1 = W+M -H-B   p2 = W+B-H-M
c   h1 = H          w1 = W
c ------------------------------------------------------------------------
c
c                            +------------------------------+
c                            |        F I E L D   ( 0 )     |
c                            +------------------------------+
c
                  if((lcy).and.(lcx)) then                       ! Field (0)
                  omn1= stf11(1)+stf31(1)
                  omp1= stf11(5)+stf31(5)
                  t1n1=0.5d+00*(stf21(2)+stf41(2))
                  t1p1=0.5d+00*(stf21(6)+stf41(6))
                  omn2= stf12(1)+stf32(1)
                  omp2= stf12(5)+stf32(5)
                  t1n2=0.5d+00*(stf22(2)+stf42(2))
                  t1p2=0.5d+00*(stf22(6)+stf42(6))
c
c ------------------------------------------------------ Exchange
c
                  g1t = c11*(omn1+omp1)+h1*(t1n1+t1p1)+ ! the same for p and n
     *                  c12*(omn2+omp2)+h2*(t1n2+t1p2)
c
                  gp1p(igp) = gp1p(igp) + ( g1t +           ! protons
     *            c21*omp1-w1*t1p1 + c22*omp2-w2*t1p2 )*feq
                  gn1p(igp) = gn1p(igp) + ( g1t +           ! neutrons
     *            c21*omn1-w1*t1n1 + c22*omn2-w2*t1n2 )*feq
c
                  g2t = c11*(omn1+omp1)-h1*(t1n1+t1p1)+ ! the same for p and n
     *                  c12*(omn2+omp2)-h2*(t1n2+t1p2)
                  gp2p(igp) = gp2p(igp) + ( g2t +           ! protons
     *            c21*omp1+w1*t1p1 + c22*omp2+w2*t1p2 )*feq
                  gn2p(igp) = gn2p(igp) + ( g2t +           ! neutrons
     *            c21*omn1+w1*t1n1 + c22*omn2+w2*t1n2 )*feq
c
c ------------------------------------------------------ Pairing
c
                  ddps = 0.5d+00*(p11*(stf11(7)+stf31(7))+
     *                            p12*(stf12(7)+stf32(7)))
                  ddns = 0.5d+00*(p11*(stf11(3)+stf31(3))+
     *                            p12*(stf12(3)+stf32(3)))
                  ddpa =-0.5d+00*(p21*(stf21(8)-stf41(8))+
     *                            p22*(stf22(8)-stf42(8)))
                  ddna =-0.5d+00*(p21*(stf21(4)-stf41(4))+
     *                            p22*(stf22(4)-stf42(4)))
                  dp1p(igp) = dp1p(igp) + feq*(ddps+ddpa)   ! prot   k1,k2
                  dn1p(igp) = dn1p(igp) + feq*(ddns+ddna)   ! neut
                  dp2p(igp) = dp2p(igp) + feq*(ddps-ddpa)   !        k2,k1
                  dn2p(igp) = dn2p(igp) + feq*(ddns-ddna)
c
c                            +------------------------------+
c                            |        F I E L D   ( 1 )     |
c                            +------------------------------+
c
                  else if(lcy.and.(.not.lcx)) then      ! Field (1)
                  t1n1 = 0.5d+00*(stf21(1)+stf41(1))
                  t1p1 = 0.5d+00*(stf21(5)+stf41(5))
                  t2n1 = 0.5d+00*(stf11(2)+stf31(2))
                  t2p1 = 0.5d+00*(stf11(6)+stf31(6))
                  t1n2 = 0.5d+00*(stf22(1)+stf42(1))
                  t1p2 = 0.5d+00*(stf22(5)+stf42(5))
                  t2n2 = 0.5d+00*(stf12(2)+stf32(2))
                  t2p2 = 0.5d+00*(stf12(6)+stf32(6))
                  g1t  = h1*(t1n1+t1p1+t2n1+t2p1)
     *                 + h2*(t1n2+t1p2+t2n2+t2p2)
                  g2t  = h1*(t1n1+t1p1-t2n1-t2p1)
     *                 + h2*(t1n2+t1p2-t2n2-t2p2)
c
c ------------------------------------------------------ Exchange
c
                  gp1p(igp) = gp1p(igp) + feq*(g1t - w1*(t1p1+t2p1)
     *                      -w2*(t1p2+t2p2))
                  gn1p(igp) = gn1p(igp) + feq*(g1t - w1*(t1n1+t2n1)
     *                      -w2*(t1n2+t2n2))
c
                  gp2p(igp) = gp2p(igp) + feq*(g2t - w1*(t1p1-t2p1)
     *                      -w2*(t1p2-t2p2))
                  gn2p(igp) = gn2p(igp) + feq*(g2t - w1*(t1n1-t2n1)
     *                      -w2*(t1n2-t2n2))
c
c ------------------------------------------------------ Pairing
c
                  ddps = 0.5d+00*(p21*(stf21(7)+stf41(7))+
     *                            p22*(stf22(7)+stf42(7)))
                  ddns = 0.5d+00*(p21*(stf21(3)+stf41(3))+
     *                            p22*(stf22(3)+stf42(3)))
                  ddpa =-0.5d+00*(p21*(stf11(8)-stf31(8))+
     *                            p22*(stf12(8)-stf32(8)))
                  ddna =-0.5d+00*(p21*(stf11(4)-stf31(4))+
     *                            p22*(stf12(4)-stf32(4)))
                  dp1p(igp) = dp1p(igp) + feq*(ddps+ddpa)   ! prot   k1,k2
                  dn1p(igp) = dn1p(igp) + feq*(ddns+ddna)   ! neut
                  dp2p(igp) = dp2p(igp) + feq*(ddps-ddpa)   !        k2,k1
                  dn2p(igp) = dn2p(igp) + feq*(ddns-ddna)
c
c                            +------------------------------+
c                            |        F I E L D   ( 2 )     |
c                            +------------------------------+
c
                  else if((.not.lcy).and.lcx) then      ! Field (2)
c
c  YOUR ATENTION PLEASE: Due to the definition of YZOPTI there is
c  a change of sign in the q q' reversed part for Fields (2) and (3)
c

                  omn1= stf11(2)+stf31(2)
                  omp1= stf11(6)+stf31(6)
                  t1n1=0.5d+00*(stf21(1)+stf41(1))
                  t1p1=0.5d+00*(stf21(5)+stf41(5))
                  omn2= stf12(2)+stf32(2)
                  omp2= stf12(6)+stf32(6)
                  t1n2=0.5d+00*(stf22(1)+stf42(1))
                  t1p2=0.5d+00*(stf22(5)+stf42(5))
c
c ------------------------------------------------------ Exchange
c
                  g1t = c11*(omn1+omp1)+h1*(t1n1+t1p1)  ! the same for p and n
     *                 +c12*(omn2+omp2)+h2*(t1n2+t1p2)
c
                  gp1p(igp) = gp1p(igp) + ( g1t +           ! protons
     *            c21*omp1-w1*t1p1 + c22*omp2-w2*t1p2 )*feq
                  gn1p(igp) = gn1p(igp) + ( g1t +           ! neutrons
     *            c21*omn1-w1*t1n1 + c22*omn2-w2*t1n2 )*feq
c
                  g2t =-c11*(omn1+omp1)+h1*(t1n1+t1p1)  ! the same for p and n
     *                 -c12*(omn2+omp2)+h2*(t1n2+t1p2)
                  gp2p(igp) = gp2p(igp) + ( g2t -           ! protons
     *            c21*omp1-w1*t1p1 - c22*omp2-w2*t1p2 )*feq
                  gn2p(igp) = gn2p(igp) + ( g2t -           ! neutrons
     *            c21*omn1-w1*t1n1 - c22*omn2-w2*t1n2 )*feq
c
c ------------------------------------------------------ Pairing
c
                  ddps = 0.5d+00*(p21*(stf21(7)+stf41(7))+
     *                            p22*(stf22(7)+stf42(7)))
                  ddns = 0.5d+00*(p21*(stf21(3)+stf41(3))+
     *                            p22*(stf22(3)+stf42(3)))
                  ddpa =-0.5d+00*(p11*(stf11(8)-stf31(8))+
     *                            p12*(stf12(8)-stf32(8)))
                  ddna =-0.5d+00*(p11*(stf11(4)-stf31(4))+
     *                            p12*(stf12(4)-stf32(4)))
                  dp1p(igp) = dp1p(igp) + feq*(ddps+ddpa)   ! prot   k1,k2
                  dn1p(igp) = dn1p(igp) + feq*(ddns+ddna)   ! neut
                  dp2p(igp) = dp2p(igp) + feq*(ddps-ddpa)   !        k2,k1
                  dn2p(igp) = dn2p(igp) + feq*(ddns-ddna)
c
c                            +------------------------------+
c                            |        F I E L D   ( 3 )     |
c                            +------------------------------+
c
                  else                                  ! Field (3)
                  t1n1 = 0.5d+00*(stf21(2)+stf41(2))
                  t1p1 = 0.5d+00*(stf21(6)+stf41(6))
                  t2n1 = 0.5d+00*(stf11(1)+stf31(1))
                  t2p1 = 0.5d+00*(stf11(5)+stf31(5))
                  t1n2 = 0.5d+00*(stf22(2)+stf42(2))
                  t1p2 = 0.5d+00*(stf22(6)+stf42(6))
                  t2n2 = 0.5d+00*(stf12(1)+stf32(1))
                  t2p2 = 0.5d+00*(stf12(5)+stf32(5))
                  g1t  = h1*(t1n1+t1p1+t2n1+t2p1)
     *                 + h2*(t1n2+t1p2+t2n2+t2p2)
                  g2t  = h1*(t2n1+t2p1-t1n1-t1p1)
     *                 + h2*(t2n2+t2p2-t1n2-t1p2)
c
c ------------------------------------------------------ Exchange
c
                  gp1p(igp) = gp1p(igp) + feq*(g1t - w1*(t1p1+t2p1)
     *                      -w2*(t1p2+t2p2))
                  gn1p(igp) = gn1p(igp) + feq*(g1t - w1*(t1n1+t2n1)
     *                      -w2*(t1n2+t2n2))
c
                  gp2p(igp) = gp2p(igp) + feq*(g2t - w1*(t2p1-t1p1)
     *                      -w2*(t2p2-t1p2))
                  gn2p(igp) = gn2p(igp) + feq*(g2t - w1*(t2n1-t1n1)
     *                      -w2*(t2n2-t1n2))
c
c ------------------------------------------------------ Pairing
c
                  ddps = 0.5d+00*(p21*(stf11(7)+stf31(7))+
     *                            p22*(stf12(7)+stf32(7)))
                  ddns = 0.5d+00*(p21*(stf11(3)+stf31(3))+
     *                            p22*(stf12(3)+stf32(3)))
                  ddpa =-0.5d+00*(p21*(stf21(8)-stf41(8))+
     *                            p22*(stf22(8)-stf42(8)))
                  ddna =-0.5d+00*(p21*(stf21(4)-stf41(4))+
     *                            p22*(stf22(4)-stf42(4)))
                  dp1p(igp) = dp1p(igp) + feq*(ddps+ddpa)   ! prot   k1,k2
                  dn1p(igp) = dn1p(igp) + feq*(ddns+ddna)   ! neut
                  dp2p(igp) = dp2p(igp) + feq*(ddps-ddpa)   !        k2,k1
                  dn2p(igp) = dn2p(igp) + feq*(ddns-ddna)

                  end if                                ! end if FIELDS
   31           continue
   30         continue
c
c---->    Negative parity
c
              end if       ! lmze condition
c
c be careful not to move the if below, the iz1o=iy1o need to be done
c
              if(lmzo(nx1,ny1).and.lmzo(nx2,ny2)) then
c
              do 40 nz1 =nzi1o,nzm1,2
                nnz1 = (nz1-1)*nzmax
                if(lxy12) nzm2=nz1
                do 41 nz2 = nzi2o,nzm2,2
                  nnz2   = (nz2-1)*nzmax
                  igm    = igm + 1
                  lxyz12 = lxy12.and.(nz1.eq.nz2)
                  feq    = 1.0d+00
c                  if(lxyz12) feq = 0.5d+00
                  iyz12  = iyz120 + 16*(nnz1+nz2-1)             ! pos in yzopti
                  jyz12  = jyz120 + 16*(nnz2+nz1-1)             ! q q' reversed
c
                  do 400 iro=1,8                               ! range 1
                  iyz12      = iyz12 + 1
                  jyz12      = jyz12 + 1
                  stf11(iro) = xxtj1 *yzopti(iyz12)   ! standar
                  stf21(iro) = xxtj1p*yzopti(iyz12)   ! tilde
                  stf31(iro) = xxtk1 *yzopti(jyz12)   ! q q' rev.
                  stf41(iro) = xxtk1p*yzopti(jyz12)   !
  400             continue
c
                  do 401 iro=1,8                              ! range 2
                  iyz12      = iyz12 + 1
                  jyz12      = jyz12 + 1
                  stf12(iro) = xxtj2 *yzopti(iyz12)   ! standar
                  stf22(iro) = xxtj2p*yzopti(iyz12)   ! tilde
                  stf32(iro) = xxtk2 *yzopti(jyz12)   ! q q' rev.
                  stf42(iro) = xxtk2p*yzopti(jyz12)   !
  401             continue
c
c
c                            +------------------------------+
c                            |        F I E L D   ( 0 )     |
c                            +------------------------------+
c
                  if((lcy).and.(lcx)) then                       ! Field (0)
                  omn1= stf11(1)+stf31(1)
                  omp1= stf11(5)+stf31(5)
                  t1n1=0.5d+00*(stf21(2)+stf41(2))
                  t1p1=0.5d+00*(stf21(6)+stf41(6))
                  omn2= stf12(1)+stf32(1)
                  omp2= stf12(5)+stf32(5)
                  t1n2=0.5d+00*(stf22(2)+stf42(2))
                  t1p2=0.5d+00*(stf22(6)+stf42(6))
c
c ------------------------------------------------------ Exchange
c
                  g1t = c11*(omn1+omp1)+h1*(t1n1+t1p1)+ ! the same for p and n
     *                  c12*(omn2+omp2)+h2*(t1n2+t1p2)
c
                  gp1m(igm) = gp1m(igm) + ( g1t +           ! protons
     *            c21*omp1-w1*t1p1 + c22*omp2-w2*t1p2 )*feq
                  gn1m(igm) = gn1m(igm) + ( g1t +           ! neutrons
     *            c21*omn1-w1*t1n1 + c22*omn2-w2*t1n2 )*feq
c
                  g2t = c11*(omn1+omp1)-h1*(t1n1+t1p1)+ ! the same for p and n
     *                  c12*(omn2+omp2)-h2*(t1n2+t1p2)
                  gp2m(igm) = gp2m(igm) + ( g2t +           ! protons
     *            c21*omp1+w1*t1p1 + c22*omp2+w2*t1p2 )*feq
                  gn2m(igm) = gn2m(igm) + ( g2t +           ! neutrons
     *            c21*omn1+w1*t1n1 + c22*omn2+w2*t1n2 )*feq
c
c ------------------------------------------------------ Pairing
c
                  ddps = 0.5d+00*(p11*(stf11(7)+stf31(7))+
     *                            p12*(stf12(7)+stf32(7)))
                  ddns = 0.5d+00*(p11*(stf11(3)+stf31(3))+
     *                            p12*(stf12(3)+stf32(3)))
                  ddpa =-0.5d+00*(p21*(stf21(8)-stf41(8))+
     *                            p22*(stf22(8)-stf42(8)))
                  ddna =-0.5d+00*(p21*(stf21(4)-stf41(4))+
     *                            p22*(stf22(4)-stf42(4)))
                  dp1m(igm) = dp1m(igm) + feq*(ddps+ddpa)   ! prot   k1,k2
                  dn1m(igm) = dn1m(igm) + feq*(ddns+ddna)   ! neut
                  dp2m(igm) = dp2m(igm) + feq*(ddps-ddpa)   !        k2,k1
                  dn2m(igm) = dn2m(igm) + feq*(ddns-ddna)
c
c                            +------------------------------+
c                            |        F I E L D   ( 1 )     |
c                            +------------------------------+
c
                  else if(lcy.and.(.not.lcx)) then      ! Field (1)
                  t1n1 = 0.5d+00*(stf21(1)+stf41(1))
                  t1p1 = 0.5d+00*(stf21(5)+stf41(5))
                  t2n1 = 0.5d+00*(stf11(2)+stf31(2))
                  t2p1 = 0.5d+00*(stf11(6)+stf31(6))
                  t1n2 = 0.5d+00*(stf22(1)+stf42(1))
                  t1p2 = 0.5d+00*(stf22(5)+stf42(5))
                  t2n2 = 0.5d+00*(stf12(2)+stf32(2))
                  t2p2 = 0.5d+00*(stf12(6)+stf32(6))
                  g1t  = h1*(t1n1+t1p1+t2n1+t2p1)
     *                 + h2*(t1n2+t1p2+t2n2+t2p2)
                  g2t  = h1*(t1n1+t1p1-t2n1-t2p1)
     *                 + h2*(t1n2+t1p2-t2n2-t2p2)
c
c ------------------------------------------------------ Exchange
c
                  gp1m(igm) = gp1m(igm) + feq*(g1t - w1*(t1p1+t2p1)
     *                      -w2*(t1p2+t2p2))
                  gn1m(igm) = gn1m(igm) + feq*(g1t - w1*(t1n1+t2n1)
     *                      -w2*(t1n2+t2n2))
c
                  gp2m(igm) = gp2m(igm) + feq*(g2t - w1*(t1p1-t2p1)
     *                      -w2*(t1p2-t2p2))
                  gn2m(igm) = gn2m(igm) + feq*(g2t - w1*(t1n1-t2n1)
     *                      -w2*(t1n2-t2n2))
c
c ------------------------------------------------------ Pairing
c
                  ddps = 0.5d+00*(p21*(stf21(7)+stf41(7))+
     *                            p22*(stf22(7)+stf42(7)))
                  ddns = 0.5d+00*(p21*(stf21(3)+stf41(3))+
     *                            p22*(stf22(3)+stf42(3)))
                  ddpa =-0.5d+00*(p21*(stf11(8)-stf31(8))+
     *                            p22*(stf12(8)-stf32(8)))
                  ddna =-0.5d+00*(p21*(stf11(4)-stf31(4))+
     *                            p22*(stf12(4)-stf32(4)))
                  dp1m(igm) = dp1m(igm) + feq*(ddps+ddpa)   ! prot   k1,k2
                  dn1m(igm) = dn1m(igm) + feq*(ddns+ddna)   ! neut
                  dp2m(igm) = dp2m(igm) + feq*(ddps-ddpa)   !        k2,k1
                  dn2m(igm) = dn2m(igm) + feq*(ddns-ddna)
c
c                            +------------------------------+
c                            |        F I E L D   ( 2 )     |
c                            +------------------------------+
c
                  else if((.not.lcy).and.lcx) then      ! Field (2)
c
c  YOUR ATENTION PLEASE: Due to the definition of YZOPTI there is
c  a change of sign in the q q' reversed part for Fields (2) and (3)
c

                  omn1= stf11(2)+stf31(2)
                  omp1= stf11(6)+stf31(6)
                  t1n1=0.5d+00*(stf21(1)+stf41(1))
                  t1p1=0.5d+00*(stf21(5)+stf41(5))
                  omn2= stf12(2)+stf32(2)
                  omp2= stf12(6)+stf32(6)
                  t1n2=0.5d+00*(stf22(1)+stf42(1))
                  t1p2=0.5d+00*(stf22(5)+stf42(5))
c
c ------------------------------------------------------ Exchange
c
                  g1t = c11*(omn1+omp1)+h1*(t1n1+t1p1)  ! the same for p and n
     *                 +c12*(omn2+omp2)+h2*(t1n2+t1p2)
c
                  gp1m(igm) = gp1m(igm) + ( g1t +           ! protons
     *            c21*omp1-w1*t1p1 + c22*omp2-w2*t1p2 )*feq
                  gn1m(igm) = gn1m(igm) + ( g1t +           ! neutrons
     *            c21*omn1-w1*t1n1 + c22*omn2-w2*t1n2 )*feq
c
                  g2t =-c11*(omn1+omp1)+h1*(t1n1+t1p1)  ! the same for p and n
     *                 -c12*(omn2+omp2)+h2*(t1n2+t1p2)
                  gp2m(igm) = gp2m(igm) + ( g2t -           ! protons
     *            c21*omp1-w1*t1p1 - c22*omp2-w2*t1p2 )*feq
                  gn2m(igm) = gn2m(igm) + ( g2t -           ! neutrons
     *            c21*omn1-w1*t1n1 - c22*omn2-w2*t1n2 )*feq
c
c ------------------------------------------------------ Pairing
c
                  ddps = 0.5d+00*(p21*(stf21(7)+stf41(7))+
     *                            p22*(stf22(7)+stf42(7)))
                  ddns = 0.5d+00*(p21*(stf21(3)+stf41(3))+
     *                            p22*(stf22(3)+stf42(3)))
                  ddpa =-0.5d+00*(p11*(stf11(8)-stf31(8))+
     *                            p12*(stf12(8)-stf32(8)))
                  ddna =-0.5d+00*(p11*(stf11(4)-stf31(4))+
     *                            p12*(stf12(4)-stf32(4)))
                  dp1m(igm) = dp1m(igm) + feq*(ddps+ddpa)   ! prot   k1,k2
                  dn1m(igm) = dn1m(igm) + feq*(ddns+ddna)   ! neut
                  dp2m(igm) = dp2m(igm) + feq*(ddps-ddpa)   !        k2,k1
                  dn2m(igm) = dn2m(igm) + feq*(ddns-ddna)
c
c                            +------------------------------+
c                            |        F I E L D   ( 3 )     |
c                            +------------------------------+
c
                  else                                  ! Field (3)
                  t1n1 = 0.5d+00*(stf21(2)+stf41(2))
                  t1p1 = 0.5d+00*(stf21(6)+stf41(6))
                  t2n1 = 0.5d+00*(stf11(1)+stf31(1))
                  t2p1 = 0.5d+00*(stf11(5)+stf31(5))
                  t1n2 = 0.5d+00*(stf22(2)+stf42(2))
                  t1p2 = 0.5d+00*(stf22(6)+stf42(6))
                  t2n2 = 0.5d+00*(stf12(1)+stf32(1))
                  t2p2 = 0.5d+00*(stf12(5)+stf32(5))
                  g1t  = h1*(t1n1+t1p1+t2n1+t2p1)
     *                 + h2*(t1n2+t1p2+t2n2+t2p2)
                  g2t  = h1*(t2n1+t2p1-t1n1-t1p1)
     *                 + h2*(t2n2+t2p2-t1n2-t1p2)
c
c ------------------------------------------------------ Exchange
c
                  gp1m(igm) = gp1m(igm) + feq*(g1t - w1*(t1p1+t2p1)
     *                      -w2*(t1p2+t2p2))
                  gn1m(igm) = gn1m(igm) + feq*(g1t - w1*(t1n1+t2n1)
     *                      -w2*(t1n2+t2n2))
c
                  gp2m(igm) = gp2m(igm) + feq*(g2t - w1*(t2p1-t1p1)
     *                      -w2*(t2p2-t1p2))
                  gn2m(igm) = gn2m(igm) + feq*(g2t - w1*(t2n1-t1n1)
     *                      -w2*(t2n2-t1n2))
c
c ------------------------------------------------------ Pairing
c
                  ddps = 0.5d+00*(p21*(stf11(7)+stf31(7))+
     *                            p22*(stf12(7)+stf32(7)))
                  ddns = 0.5d+00*(p21*(stf11(3)+stf31(3))+
     *                            p22*(stf12(3)+stf32(3)))
                  ddpa =-0.5d+00*(p21*(stf21(8)-stf41(8))+
     *                            p22*(stf22(8)-stf42(8)))
                  ddna =-0.5d+00*(p21*(stf21(4)-stf41(4))+
     *                            p22*(stf22(4)-stf42(4)))
                  dp1m(igm) = dp1m(igm) + feq*(ddps+ddpa)   ! prot   k1,k2
                  dn1m(igm) = dn1m(igm) + feq*(ddns+ddna)   ! neut
                  dp2m(igm) = dp2m(igm) + feq*(ddps-ddpa)   !        k2,k1
                  dn2m(igm) = dn2m(igm) + feq*(ddns-ddna)
c
                  end if                                ! end if FIELDS

   41           continue
   40         continue
            end if        ! lmzo
   12       continue
    2     continue
   11   continue
    1 continue
      call timeit(1,17,'GECH: GECHXYZ   ')
      return
      end
c+---------------------------------------------------------------------+
c|   Computes Y(ny1,nz1;ny2,nz2) for optimization of the exchange term |
c|          ny2-ny1         nqy1-nqy                                   |
c|   The (i)         and (i)         are included                      |
c|                                                                     |
c|                           YZOPTI                                    |
c|  Y(ny1,nz1;ny2,nz2)=                                                |
c|                                                                     |
c|    ___                                                              |
c|    \     BB,y                              ny2+nqy1-ny1-nqy         |
c|     >   I                  Z  (nz1,nz2)  (i)                        |
c|    /__   ny1 nyq nyq' ny2   nqx,nqy;nqx',nqy'                       |
c|  nyq',nyq                                                           |
c|                                                                     |
c+---------------------------------------------------------------------+
      Subroutine Gechyzo(iytj,iytj2,nqy1,nqy,xtj,zopt,yzopt)
      Implicit Real*8 (A-H,O-Z)
      Logical LMZE,LMZO
      Include 'COMDIM'
      Dimension iytj(nymax,nymax),iytj2(nymax,nymax)
      Dimension zopt(mazopti),yzopt(mayopti)
      Dimension xtj(2,maxtjy)
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
c
      call timeit(0,16,'GECH: GECHYZO   ')
      ny2max= 16*nzmax*nzmax
      imyq = Mod((nqy1+nqy),2)
      do 200 ny1=1,nymax
         nny1 = (ny1-1)*nymax
         jxtj = iytj(ny1,nqy1)  ! starting pos of xtj
         nzm1 = MZ(1,ny1)
	 nyi2 = 1 + Mod(imyq+ny1-1,2)
c         do 201 ny2 =1,nymax
         do 201 ny2 =nyi2,nymax,2
c            imy12= Mod((ny1+ny2),2)
c            if(imyq.ne.imy12) goto 201
            nzm2 = MZ(1,ny2)
            phas= Dfloat(1-2*Mod(Iabs(ny2+nqy1-ny1-nqy)/2,2))
            jx = jxtj + iytj2(ny2,nqy )
c            xxtj1 = xtj(1,jx)
c            xxtj2 = xtj(2,jx)
            xxtj1 = xtj(1,jx)*phas
            xxtj2 = xtj(2,jx)*phas
            jyz12 = ny2max*(nny1+ny2-1)    ! ny2max elements each ny1 ny2
c
c It does not compute for all posible nz, just for the maximum allowed
c by the ny+nz .le. nzmax condition ( most unfavorable condition nx=0)
c
            nnz1 = 0
            do 300 nz1=1,nzm1
c               nnz1 = (nz1-1)*nzmax
	       iz12 = nnz1
	       iyz12 = jyz12 + nnz1
               do 301 nz2=1,nzm2
c                  iz12 = (nnz1+nz2-1)*16
c                  iyz12= jyz12+iz12
c! zopt(8,2,nzmax,nzmax)
c		  call daxpy(8,xxtj1,zopt(iz12+1),1,yzopt(iyz12+1),1)
c		  call daxpy(8,xxtj2,zopt(iz12+9),1,yzopt(iyz12+9),1)
c si pones daxpy hay que poner tambien iz12 = iz12 + 16, etc the abajo
                  do iro=1,8              ! alcance 1
                    iz12  = iz12  + 1
                    iyz12 = iyz12 + 1
                    yzopt(iyz12)=yzopt(iyz12)+xxtj1*zopt(iz12)
		  end do
                  do iro=1,8              ! alcance 2
                    iz12  = iz12  + 1
                    iyz12 = iyz12 + 1
                    yzopt(iyz12)=yzopt(iyz12)+xxtj2*zopt(iz12)
		  end do
c               iz12  = iz12  + 16
c	       iyz12 = iyz12 + 16
  301          continue
            nnz1 = nnz1 + 16*nzmax
  300       continue
  201    continue
  200    continue
      call timeit(1,16,'GECH: GECHYZO   ')
      return
      end
c+---------------------------------------------------------------------+
c|   Computes Z(nz1,nz2) for optimization of the exchange term         |
c|                                                                     |
c|                                                                     |
c| ZOPTI                  __                                           |
c|                       \     BB                __                    |
c|     ZA(nz1,nz2)   =    >   I                  RO                    |
c|                       /__   nz1 nzq nzq' nz2    q q'                |
c+---------------------------------------------------------------------+
      Subroutine gechzo(iztj,iztj2,nzq1,nzq2,rozo,xtj,zopt)
      Implicit Real*8 (A-H,O-Z)
      Logical LMZE,LMZO
      Include 'COMDIM'
      Dimension iztj(nzmax,nzmax),iztj2(nzmax,nzmax)
      Dimension xtj(2,maxtjz)
      Dimension zopt(mazopti),rozo(8)
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
c
      call timeit(0,15,'GECH: GECHZO    ')
      imzq = nzq1+nzq2
      nnz1 = 0
      do nz1=1,nzmax
         jxtj = iztj(nz1,nzq1)  ! starting pos of xtj
c         nnz1 = (nz1-1)*nzmax
         nz2i = 1+Mod((imzq+nz1-1),2)
	 nnz2 = (nz2i-1)*16
         do nz2 =nz2i,nzmax,2
            jx = jxtj + iztj2(nz2,nzq2)
            xxtj1 = xtj(1,jx)
            xxtj2 = xtj(2,jx)
c            iz12  = 16*(nnz1+nz2-1)
            iz12 = nnz1 + nnz2 
c	    call dger(8,2,1.0d+00,rozo,1,xtj(1,jx),1,zopt(iz12+1),8)
            call daxpy(8,xxtj1,rozo,1,zopt(iz12+1),1)
	    call daxpy(8,xxtj2,rozo,1,zopt(iz12+9),1)
c            do iro=1,8
c               iz12 = iz12 + 1
c               zopt(iz12) = zopt(iz12) + xxtj1*rozo(iro)
c            end do
c            do iro=1,8
c               iz12 = iz12 + 1
c               zopt(iz12) = zopt(iz12) + xxtj2*rozo(iro)
c            end do
            nnz2 = nnz2 + 32 
         end do
	 nnz1 = nnz1 + 16*nzmax
      end do
      call timeit(1,15,'GECH: GECHZO    ')
      return
      end
