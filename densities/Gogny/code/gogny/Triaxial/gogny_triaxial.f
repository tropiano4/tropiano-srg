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
c+---------------------------------------------------------------------+
c|  First version May   1999                                           |
c+=====================================================================+
c|                                                                     |
c|               S P I N   O R B I T   P A I R I N G                   |
c|                                                                     |
c+---------------------------------------------------------------------+
      Subroutine Prepaiso(isym,Soxyz,GROP,GROM)
c
      Implicit real*8(a-h,o-z)
      Implicit logical (l)
      Include 'COMDIM'
c-------------------------------------------------------------- internal
      Parameter (KAXLSO=((imax*(imax+1))/2)*(imax+3))
      Parameter (KAXKSO= imax*(imax+1)*(imax+3))
      Dimension ALSOMU(kaxlso),AKSOMU(kaxkso)
      Dimension ILSOMU(iMAX,iMAX),IKSOMU(iMAX,iMAX)
      Parameter (NSOZ=12,NSOYZ=18,NSOXYZ=24)
      Dimension Soz (izsrc,NSOZ)
      Dimension Soyz(iyzsrc,NSOYZ)
c--------------------------------------------------------- end internal
      Dimension Soxyz(ndmu,NSOXYZ)
c +---------------------------------------------------------------------+
c |   In  Soxyz the meaning of the second index is                      |
c |=====================================================================|
c |                         nxmu     nymu     nzmu                      |
c |      1,2 .... T yz (2)  even     even      odd                      |
c |      3,4 .... T zy (2)  even      odd     even                      |
c |      5,6 .... T xz (2)   odd      odd      odd                      |
c |      7,8 .... T xy (2)   odd      odd      odd                      |
c +---------------------------------------------------------------------+
c |      9,10.... T xz (1)  even     even      odd                      |
c |      11,12... T xy (3)  even      odd     even                      |
c |      13,14... T yx (3)   odd     even     even                      |
c |      15,16... T zx (1)   odd     even     even                      |
c |      17,18... T yz (1)   odd      odd      odd                      |
c |      19,20... T yx (1)   odd      odd      odd                      |
c |      21,22... T zy (3)   odd      odd      odd                      |
c |      23,24... T zx (3)   odd      odd      odd                      |
c +---------------------------------------------------------------------+
c |      Odd indices are for neutrons while even ones are for protons   |
c +---------------------------------------------------------------------+
      Dimension GROP(nrop8),GROM(nrom8)
c
      Save alsomu,aksomu,ilsomu,iksomu,icall
c      
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      Common/DIMEN/kmax,kmax1,kxmax,kymax,kzmax,nwf2,
     * maxtz1,maxjz1,maxjx1,maxjy1,maxlso,maxdso,
     * ndmu,ndthet,nacou,nrop,nrom,nrop8,nrom8,nherm,nherm38,nlegn
     
      Common/OSCLEN/bx,by,bz
      
c      write(6,*) ' PREPAISO ISYM = ',isym
c-------------------------------------------------------------------
      if(icall.ne.061060) then
      call inigf
      Maxmu = 2*nmax+2
      indls = 1
      do i=1,nmax
         do j=i,nmax
	    ilsomu(i,j) = indls
	    ilsomu(j,i) = indls
            Minmu  = Mod(i+j,2)+1
	    do imu=Minmu,Maxmu,2
	       alsomu(indls) = DL1(i-1,j-1,imu-1)
	       indls = indls + 1
	    end do
	 end do
      end do
      if((indls-1).gt.KAXLSO) then
         write(6,*) ' indls gt KAXLSO in PREPAISO '
	 stop
      end if
      ff = 1.0d+00/dsqrt(2.0d+00)
      Maxmu = 2*nmax+2
      indks = 1
      do i=1,nmax
         do j=1,nmax
	    iksomu(i,j) = indks
            Minmu  = 2-Mod(i+j,2)
	    do imu=Minmu,Maxmu,2
	       aa1 =  dsqrt(dfloat(j-1))*DL1(i-1,j-2,imu-1)
	       aa2 = -dsqrt(dfloat(j  ))*DL1(i-1,j  ,imu-1)
	       aa3 = -dsqrt(dfloat(i-1))*DL1(i-2,j-1,imu-1)
	       aa4 =  dsqrt(dfloat(i  ))*DL1(i  ,j-1,imu-1)
	       aksomu(indks) = ff*(aa1+aa2+aa3+aa4)
	       indks = indks + 1
	    end do
	 end do
      end do
      if((indks-1).gt.KAXKSO) then
         write(6,*) ' indks gt KAXKSO in PREPAISO '
	 stop
      end if
	    
      icall = 061060
      	    
      end if
      
      imuze = 1 + isym
      imuzo = 2 - isym	    
      
      imuye = 1 + isym
      imuyo = 2 - isym	    
      
      mumax = 2*nxmax
      mumay = 2*nymax
      mumaz = 2*nzmax
      
      irop = 0
      irom = 0

      
      do ind=1,NSOXYZ
         do imuxyz=1,ndmu
            soxyz(imuxyz,ind) = 0.0d+0
	 end do
      end do
      
c+------------------------------------------------------------------+
c      isym=0        symmetric kappa 
c           1        skew-symmetric kappa
c      
c+------------------------------------------------------------------+
c             1 -> q'     2 -> q      2<=1 (lower triangle)
c
c     Kappa
c           q q'
c+------------------------------------------------------------------+
      do nx1=1,nxmax
         nym1 = MY(nx1)
         do nx2=1,nx1
            nym2 = MY(nx2)
            lx12 = nx1.eq.nx2
            lpx  = Mod(nx1+nx2,2).eq.0
	  
            do indy=1,NSOYZ
               do imuyz=1,iyzsrc
                  soyz(imuyz,indy) = 0.0d+0
	       end do
            end do
c
          do ny1 = 1,nym1
            nzi1e = nzie(nx1,ny1)
            nzi1o = nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
            if(lx12) nym2 = ny1
            do ny2=1,nym2
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              lxy12 = lx12.and.(ny1.eq.ny2)
              lpy  = Mod(ny1+ny2+isym,2).eq.0
c
              do indz=1,NSOZ
         	 do imuz=1,izsrc
         	    soz(imuz,indz) = 0.0d+0
         	 end do
              end do
c+---------------------------------------------------------------------+
c|                                                 Positive parity     |
c+---------------------------------------------------------------------+
c
c be careful not to move the if below, the iz1e=iy1e need to be done
c
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
              do nz1 =nzi1e,nzm1,2
                if(lxy12) nzm2=nz1
                do nz2 = nzi2e,nzm2,2
c
                 irop  = irop + 1
                 akan = GROP(irop+4+isym) ! kappa 
                 akap = GROP(irop+6+isym) ! kappa 
                 irop   = irop + 7
                 lxyz = lxy12.and.(nz2.eq.nz1)
                 if(lxyz) then
                    akan = 0.5d+00*akan
                    akap = 0.5d+00*akap
                 end if
		 
                 ils = ilsomu(nz2,nz1)
                 iks = iksomu(nz2,nz1)
                 indls= ils
                 indks= iks
c 		      +--------------------------------------------+
c 		      |   No  Field (0)   contribution  	   |
c 		      | 					   |
c 		      +--------------------------------------------+

c 		      +--------------------------------------------+
c 		      | 	Field (1)			   |
c 		      | 					   |
c 		      |    Tyz, Txz, Tyx	    Muz odd	   |
c 		      |    Tzx  		    Muz even	   |
c 		      +--------------------------------------------+
                 if(lpy.and.(.not.lpx)) then
                   imu = 0
                   do muz=imuzo,mumaz,2
                     imu = imu + 1
		     soz(imu,1) = soz(imu,1) + akan*alsomu(indls)
		     soz(imu,2) = soz(imu,2) + akap*alsomu(indls)
                     indls = indls + 1
                   end do  
                   imu = 0
                   do muz=imuze,mumaz,2
                     imu = imu + 1
                     soz (imu,3) =soz (imu,3)+ akan*aksomu(indks )/bz
                     soz (imu,4) =soz (imu,4)+ akap*aksomu(indks )/bz
                     indks = indks + 1
                   end do
c 			+--------------------------------------------+
c 			|	  Field (2)			     |
c 			|					     |
c 			|    Tyz, Txz, Txy	      Muz odd	     |
c 			|    Tzy		      Muz even       |
c 			+--------------------------------------------+
                   else if(.not.lpy.and.lpx) then
                   imu = 0
                   do muz=imuzo,mumaz,2
                     imu = imu + 1
		     soz(imu,5) = soz(imu,5) + akan*alsomu(indls)
		     soz(imu,6) = soz(imu,6) + akap*alsomu(indls)
                     indls = indls + 1
                   end do  
                   imu = 0
                   do muz=imuze,mumaz,2
                     imu = imu + 1
                     soz (imu,7) =soz (imu,7)+ akan*aksomu(indks )/bz
                     soz (imu,8) =soz (imu,8)+ akap*aksomu(indks )/bz
                     indks = indks + 1
                   end do
c 			+--------------------------------------------+
c 			|	  Field (3)			     |
c 			|					     |
c 			|    Tzy, Tzx		      Muz odd	     |
c 			|    Txy, Tyx		      Muz even       |
c 			+--------------------------------------------+
                   else if((.not.lpy).and.(.not.lpx)) then
                   imu = 0
                   do muz=imuzo,mumaz,2
                     imu = imu + 1
		     soz(imu, 9) = soz(imu, 9) + akan*aksomu(indks)/bz
		     soz(imu,10) = soz(imu,10) + akap*aksomu(indks)/bz
                     indks = indks + 1
                   end do  
                   imu = 0
                   do muz=imuze,mumaz,2
                     imu = imu + 1
                     soz (imu,11) =soz (imu,11)+ akan*alsomu(indls )
                     soz (imu,12) =soz (imu,12)+ akap*alsomu(indls )
                     indls = indls + 1
                   end do
c
                   end if ! Fields
c
                end do    ! nz2
              end do      ! nz1
          end if          ! lmze
c+--------------------------------------------------+
c|                                                  |
c|         Paridad negativa                         |
c|                                                  |
c+--------------------------------------------------+
c be careful not to move the if below, the iz1e=iy1e need to be done
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
	      
	      
              do nz1 = nzi1o,nzm1,2
                if(lxy12) nzm2=nz1
                do nz2 = nzi2o,nzm2,2
c
                 irom  = irom + 1
                 akan = GROM(irom+4+isym) ! kappa 
                 akap = GROM(irom+6+isym) ! kappa 
                 irom   = irom + 7
                 lxyz = lxy12.and.(nz2.eq.nz1)
                 if(lxyz) then
                    akan = 0.5d+00*akan
                    akap = 0.5d+00*akap
                 end if
                 ils = ilsomu(nz2,nz1)
                 iks = iksomu(nz2,nz1)
                 indls= ils
                 indks= iks
c      	              +--------------------------------------------+
c 	              |   No  Field (0)   contribution  	   |
c 	              | 					   |
c 	              +--------------------------------------------+

c 		      +--------------------------------------------+
c 		      |	       Field (1)			   |
c 		      | 					   |
c 		      |    Tyz, Txz, Tyx	    Muz odd	   |
c 		      |    Tzx  		    Muz even	   |
c 		      +--------------------------------------------+
c
                 if(lpy.and.(.not.lpx)) then
                   imu = 0
                   do muz=imuzo,mumaz,2
                     imu = imu + 1
		     soz(imu,1) = soz(imu,1) + akan*alsomu(indls)
		     soz(imu,2) = soz(imu,2) + akap*alsomu(indls)
                     indls = indls + 1
                   end do  
                   imu = 0
                   do muz=imuze,mumaz,2
                     imu = imu + 1
                     soz (imu,3) =soz (imu,3)+ akan*aksomu(indks )/bz
                     soz (imu,4) =soz (imu,4)+ akap*aksomu(indks )/bz
                     indks = indks + 1
                   end do
c 			+--------------------------------------------+
c 			|	  Field (2)			     |
c 			|					     |
c 			|    Tyz, Txz, Txy	      Muz odd	     |
c 			|    Tzy		      Muz even       |
c 			+--------------------------------------------+
                   else if(.not.lpy.and.lpx) then
                   imu = 0
                   do muz=imuzo,mumaz,2
                     imu = imu + 1
		     soz(imu,5) = soz(imu,5) + akan*alsomu(indls)
		     soz(imu,6) = soz(imu,6) + akap*alsomu(indls)
                     indls = indls + 1
                   end do  
                   imu = 0
                   do muz=imuze,mumaz,2
                     imu = imu + 1
                     soz (imu,7) =soz (imu,7)+ akan*aksomu(indks )/bz
                     soz (imu,8) =soz (imu,8)+ akap*aksomu(indks )/bz
                     indks = indks + 1
                   end do
c 			+--------------------------------------------+
c 			|	  Field (3)			     |
c 			|					     |
c 			|    Tzy, Tzx		      Muz odd	     |
c 			|    Txy, Tyx		      Muz even       |
c 			+--------------------------------------------+
                   else if((.not.lpy).and.(.not.lpx)) then
                   imu = 0
                   do muz=imuzo,mumaz,2
                     imu = imu + 1
		     soz(imu, 9) = soz(imu, 9) + akan*aksomu(indks)/bz
		     soz(imu,10) = soz(imu,10) + akap*aksomu(indks)/bz
                     indks = indks + 1
                   end do  
                   imu = 0
                   do muz=imuze,mumaz,2
                     imu = imu + 1
                     soz (imu,11) =soz (imu,11)+ akan*alsomu(indls )
                     soz (imu,12) =soz (imu,12)+ akap*alsomu(indls )
                     indls = indls + 1
                   end do
c
                   end if ! Fields
                end do     ! nz2
              end do       ! nz1
            end if ! lmzo
c
c 	    		    	   +------------------------+
c 	    		    	   |	  Y	L O O P     |
c 	    		    	   +------------------------+
c
            ils  = ilsomu(ny2,ny1)
            indls= ils
            iks  = iksomu(ny2,ny1)
            indks= iks
c
c
c 				+--------------------------+
c 				| N O	F I E L D  ( 0 )   |
c 				+--------------------------+

c 				   +---------------------+
c 				   |  F I E L D  ( 1 )   |
c 				   +---------------------+
c
            if(lpy.and.(.not.lpx)) then
	    	    
              fasy = Dfloat(1-2*Mod(Iabs(ny2-ny1+isym)/2,2)) ! (-)**(ny2-ny1)/2
c             
c                                        muy even 
c                                      
              iyzeo = 0
              iyzee = 0
c --------------------------------------  muz odd  Txz (1)
              do muy=imuye,mumay,2
                 sls= alsomu(indls)
                 indls = indls + 1
                 imuz = 0
                 do muz=imuzo,mumaz,2
                    imuz = imuz + 1
                    iyzeo = iyzeo + 1
		    soyz(iyzeo,1)=soyz(iyzeo,1)+fasy*soz(imuz,1)*sls
		    soyz(iyzeo,2)=soyz(iyzeo,2)+fasy*soz(imuz,2)*sls
                 end do  ! muz
c --------------------------------------  muz even Tzx (1)
                 imuz = 0
                 do muz=imuze,mumaz,2
                    imuz = imuz + 1
                    iyzee = iyzee + 1
		    soyz(iyzee,3)=soyz(iyzee,3)+fasy*soz(imuz,3)*sls
		    soyz(iyzee,4)=soyz(iyzee,4)+fasy*soz(imuz,4)*sls
                 end do  ! muz
              end do     ! muy
c             
c                                        muy odd 
c                                      
              iyzoo = 0
c --------------------------------------  muz odd  Tyz (1) Tyx (1)
              do muy=imuyo,mumay,2
                 sks= aksomu(indks)/by
                 indks = indks + 1
                 imuz = 0
                 do muz=imuzo,mumaz,2
                    imuz = imuz + 1
                    iyzoo = iyzoo + 1
		    soyz(iyzoo,5)=soyz(iyzoo,5)+fasy*soz(imuz,1)*sks
		    soyz(iyzoo,6)=soyz(iyzoo,6)+fasy*soz(imuz,2)*sks
                 end do  ! muz
	      end do     ! muy
c
c 				     +------------------+
c 				     | F I E L D  ( 2 ) |
c 				     +------------------+
c
            else if(.not.lpy.and.lpx) then
	                                                       
              fasy= Dfloat(1-2*Mod(Iabs(ny2-ny1+1-isym)/2,2)) !(-)**(ny2-ny1+1)/2
c	      
c                                        muy even 
c                                      
              iyzeo = 0
c --------------------------------------  muz odd  Tyz (2)
              do muy=imuye,mumay,2
                 sks= aksomu(indks)/by
                 indks = indks + 1
                 imuz = 0
                 do muz=imuzo,mumaz,2
                    imuz = imuz + 1
                    iyzeo = iyzeo + 1
		    soyz(iyzeo,7)=soyz(iyzeo,7)+fasy*soz(imuz,5)*sks
		    soyz(iyzeo,8)=soyz(iyzeo,8)+fasy*soz(imuz,6)*sks
                 end do  ! muz
              end do  	 ! muy
c             
c                                        muy odd 
c                                      
              iyzoo = 0
              iyzoe = 0
c --------------------------------------  muz odd  Txz (2) Txy (2)
              do muy=imuyo,mumay,2
                 sls= alsomu(indls)
                 indls = indls + 1
                 imuz = 0
                 do muz=imuzo,mumaz,2
                    imuz = imuz + 1
                    iyzoo = iyzoo + 1
		    soyz(iyzoo, 9)=soyz(iyzoo, 9)+fasy*soz(imuz,5)*sls
		    soyz(iyzoo,10)=soyz(iyzoo,10)+fasy*soz(imuz,6)*sls
                 end do  ! muz
c --------------------------------------  muz even Tzy (2)
                 imuz = 0
                 do muz=imuze,mumaz,2
                    imuz = imuz + 1
                    iyzoe = iyzoe + 1
		    soyz(iyzoe,11)=soyz(iyzoe,11)+fasy*soz(imuz,7)*sls
		    soyz(iyzoe,12)=soyz(iyzoe,12)+fasy*soz(imuz,8)*sls
                 end do  ! muz
              end do  	 ! muy
c
c 				   +--------------------------+
c 				   |	F I E L D  ( 3 )      |
c 				   +--------------------------+
c
            else if((.not.lpy).and.(.not.lpx)) then
	    
              fasy= Dfloat(1-2*Mod(Iabs(ny2-ny1+1-isym)/2,2)) !(-)**(ny2-ny1+1)/2
c	      
c                                        muy even 
c                                      
              iyzee = 0
              do muy=imuye,mumay,2
                 sks= aksomu(indks)/by
                 indks = indks + 1
c --------------------------------------  muz even Tyx (3)
                 imuz = 0
                 do muz=imuze,mumaz,2
                    imuz = imuz + 1
                    iyzee = iyzee + 1
		    soyz(iyzee,13)=soyz(iyzee,13)+fasy*soz(imuz,11)*sks
		    soyz(iyzee,14)=soyz(iyzee,14)+fasy*soz(imuz,12)*sks
                 end do  ! muz
              end do  	 ! muy
c             
c                                        muy odd 
c                                      
              iyzoo = 0
              iyzoe = 0
c --------------------------------------  muz odd  Tzy Tzx (3)
              do muy=imuyo,mumay,2
                 sls= alsomu(indls)
                 indls = indls + 1
                 imuz = 0
                 do muz=imuzo,mumaz,2
                    imuz = imuz + 1
                    iyzoo = iyzoo + 1
		    soyz(iyzoo,15)=soyz(iyzoo,15)+fasy*soz(imuz, 9)*sls
		    soyz(iyzoo,16)=soyz(iyzoo,16)+fasy*soz(imuz,10)*sls
                 end do  ! muz
c --------------------------------------  muz even Txy (3)
                 imuz = 0
                 do muz=imuze,mumaz,2
                    imuz = imuz + 1
                    iyzoe = iyzoe + 1
		    soyz(iyzoe,17)=soyz(iyzoe,17)+fasy*soz(imuz,11)*sls
		    soyz(iyzoe,18)=soyz(iyzoe,18)+fasy*soz(imuz,12)*sls
                 end do  ! muz
              end do  	 ! muy
	      
c ------------------------------------------------------------------
            end if    !    F I E L D S   Y
c ------------------------------------------------------------------
            end do    ! ny2
          end do      ! ny1
c
c 	               +------------------------+
c 	               |      X     L O O P	|
c 	               +------------------------+
c
            ils  = ilsomu(nx2,nx1)
            indls= ils
            iks  = iksomu(nx2,nx1)
            indks= iks
c
            fasx  = Dfloat(1-2*Mod(nx1+1,2))  !   (-)**nxq'
            fasx1 = Dfloat(1-2*Mod(nx1,2))    ! - (-)**nxq'
c  					 +------------+
c  --------------------------------------| nx+nx' par |
c  					 +------------+
            if(lpx) then
c                                                     mux even
              ixyzeeo = 0
              ixyzeoe = 0
              do mux =1,mumax,2
                 sls  = alsomu(indls)
                 indls = indls + 1
                 iyzeo = 0
c                                                     muy even
                 do muy=imuye,mumay,2
c                                                     muz odd        T yz (2)
                    do muz=imuzo,mumaz,2
                       iyzeo = iyzeo + 1
                       ixyzeeo = ixyzeeo + 1
	 soxyz(ixyzeeo,1)=soxyz(ixyzeeo,1)+soyz(iyzeo,7)*sls*fasx
	 soxyz(ixyzeeo,2)=soxyz(ixyzeeo,2)+soyz(iyzeo,8)*sls*fasx
		    end do  ! muz
		 end do     ! muy 
		 
                 iyzoe = 0
c                                                     muy odd
                 do muy=imuyo,mumay,2
c                                                     muz even      T zy (2)
                    do muz=imuze,mumaz,2
                       iyzoe = iyzoe + 1
                       ixyzeoe = ixyzeoe + 1
	 soxyz(ixyzeoe,3)=soxyz(ixyzeoe,3)+soyz(iyzoe,11)*sls*fasx
	 soxyz(ixyzeoe,4)=soxyz(ixyzeoe,4)+soyz(iyzoe,12)*sls*fasx
		    end do  ! muz
		 end do     ! muy
	      end do        ! mux
c
c                                                     mux odd
              ixyzooo = 0
              do mux =2,mumax,2
                 sks  = aksomu(indks)/bx
                 indks = indks + 1
                 iyzoo = 0
c                                                     muy odd
                 do muy=imuyo,mumay,2
c                                               muz odd    Txz  Txy
                    do muz=imuzo,mumaz,2
                       iyzoo = iyzoo + 1
                       ixyzooo = ixyzooo + 1
c     Txz (2)		       
	      soxyz(ixyzooo,5)=soxyz(ixyzooo,5)+soyz(iyzoo, 9)*sks*fasx
	      soxyz(ixyzooo,6)=soxyz(ixyzooo,6)+soyz(iyzoo,10)*sks*fasx
c     Txy (2)         
	      soxyz(ixyzooo,7)=soxyz(ixyzooo,7)+soyz(iyzoo, 9)*sks*fasx
	      soxyz(ixyzooo,8)=soxyz(ixyzooo,8)+soyz(iyzoo,10)*sks*fasx
		    end do  ! muz
		 end do     ! muy
	      end do 	    ! mux

c
c
c 				          +----------------+
c   --------------------------------------| nx+nx' impar   |
c 	                       	          +----------------+
            else if(.not.lpx) then
c                                                         
c                                                     mux even
              ixyzeeo = 0
              ixyzeoe = 0
              do mux =1,mumax,2
                 sks  = aksomu(indks)/bx
                 indks = indks + 1
                 iyzeo = 0
c                                                     muy even
                 do muy=imuye,mumay,2
c                                                     muz odd  T XZ (1)       
                    do muz=imuzo,mumaz,2
                       iyzeo = iyzeo + 1
                       ixyzeeo = ixyzeeo + 1
	soxyz(ixyzeeo, 9)=soxyz(ixyzeeo, 9)+soyz(iyzeo,1)*sks*fasx
	soxyz(ixyzeeo,10)=soxyz(ixyzeeo,10)+soyz(iyzeo,2)*sks*fasx
		    end do  ! muz
		 end do     ! muy
		 
                 iyzoe = 0
c                                                     muy odd
                 do muy=imuyo,mumay,2
c                                                     muz even T xy (3)     
                    do muz=imuze,mumaz,2
                       iyzoe = iyzoe + 1
                       ixyzeoe = ixyzeoe + 1
	  soxyz(ixyzeoe,11)=soxyz(ixyzeoe,11)+soyz(iyzoe,17)*sks
	  soxyz(ixyzeoe,12)=soxyz(ixyzeoe,12)+soyz(iyzoe,18)*sks
		    end do  ! muz
		 end do     ! muy
	      end do 	    ! mux
c
c                                                     mux odd
              ixyzooo = 0
              ixyzoee = 0
              do mux =2,mumax,2
                 sls  = alsomu(indls)
                 indls = indls + 1
                 iyzee = 0
c                                                     muy even
                 do muy=imuye,mumay,2
c                                       muz even T zx (1) T yx (3)   
                    do muz=imuze,mumaz,2
                       iyzee = iyzee + 1
                       ixyzoee = ixyzoee + 1
c   Tyx (3)
	    soxyz(ixyzoee,13)=soxyz(ixyzoee,13)+soyz(iyzee,13)*sls
	    soxyz(ixyzoee,14)=soxyz(ixyzoee,14)+soyz(iyzee,14)*sls
c   Tzx (1)
	    soxyz(ixyzoee,15)=soxyz(ixyzoee,15)+soyz(iyzee,3)*sls*fasx1
	    soxyz(ixyzoee,16)=soxyz(ixyzoee,16)+soyz(iyzee,4)*sls*fasx1
		    end do  ! muz
		 end do     ! muy
                 iyzoo = 0
c                                                     muy odd
                 do muy=imuyo,mumay,2
c                                                     muz odd    
                    do muz=imuzo,mumaz,2
                       iyzoo = iyzoo + 1
                       ixyzooo = ixyzooo + 1
c   T yz (1)
	     soxyz(ixyzooo,17)=soxyz(ixyzooo,17)+soyz(iyzoo,5)*sls*fasx
	     soxyz(ixyzooo,18)=soxyz(ixyzooo,18)+soyz(iyzoo,6)*sls*fasx
c   T yx (1)
	     soxyz(ixyzooo,19)=soxyz(ixyzooo,19)+soyz(iyzoo,5)*sls*fasx1
	     soxyz(ixyzooo,20)=soxyz(ixyzooo,20)+soyz(iyzoo,6)*sls*fasx1
c   T zy (3)
	     soxyz(ixyzooo,21)=soxyz(ixyzooo,21)+soyz(iyzoo,15)*sls
	     soxyz(ixyzooo,22)=soxyz(ixyzooo,22)+soyz(iyzoo,16)*sls
c   T zx (3)
	     soxyz(ixyzooo,23)=soxyz(ixyzooo,23)+soyz(iyzoo,15)*sls
	     soxyz(ixyzooo,24)=soxyz(ixyzooo,24)+soyz(iyzoo,16)*sls
		    end do  ! muz
		 end do     ! muy
	      end do 	    ! mux
c ------------------------------------------------------------------
            end if   ! lpx
c ------------------------------------------------------------------
        end do  ! nx2
      end do    ! nx1
      
      return
      end
c****************************************************************************
c                         P   A   I   S   O
c****************************************************************************
      Subroutine Paiso(isym,grop,grom,dpp,dnp,dpm,dnm,THE,Soxyz)
c
      Implicit real*8 (A-H,O-Z)
      Logical LMZE,LMZO,vmux
      Logical lx12,lxy12,lxyz,lpx,lpy
c
      Include 'COMDIM'
      Parameter (KAXT1=((imax*(imax+1))/2)*(imax+3))
      Parameter (KAXGSO= imax*imax*(imax+3))
c
      Dimension t1d(KAXT1),gsomu(KAXGSO)
c
      Dimension Soxyz(NDMU,24)
c
      Dimension THE(NDTHET,24)
      
      Dimension it1d(imax,imax),igsomu(imax,imax)
      
      Dimension itheeeo(izmax,izmax)
      Dimension itheeoe(izmax,izmax)
      Dimension itheoee(izmax,izmax)
      Dimension itheooo(izmax,izmax)
c
      Dimension dpp(nrop),dnp(nrop)
      Dimension dpm(nrom),dnm(nrom)
c +---------------------------------------------------------------------+
c |   In  Soxyz the meaning of the second index is                      |
c |=====================================================================|
c |                         nxmu     nymu     nzmu                      |
c |      1,2 .... T yz (2)  even     even      odd                      |
c |      3,4 .... T zy (2)  even      odd     even                      |
c |      5,6 .... T xz (2)   odd      odd      odd                      |
c |      7,8 .... T xy (2)   odd      odd      odd                      |
c +---------------------------------------------------------------------+
c |      9,10.... T xz (1)  even     even      odd                      |
c |      11,12... T xy (3)  even      odd     even                      |
c |      13,14... T yx (3)   odd     even     even                      |
c |      15,16... T zx (1)   odd     even     even                      |
c |      17,18... T yz (1)   odd      odd      odd                      |
c |      19,20... T yx (1)   odd      odd      odd                      |
c |      21,22... T zy (3)   odd      odd      odd                      |
c |      23,24... T zx (3)   odd      odd      odd                      |
c +---------------------------------------------------------------------+
c |      Odd indices are for neutrons while even ones are for protons   |
c +---------------------------------------------------------------------+
      Dimension GROP(nrop8),GROM(nrom8)
c ----------------------------------------------------------------
      save t1d,gsomu,it1d,igsomu,icall
c ----------------------------------------------------------------      
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      Common/DIMPDD/ irortmp,idmus0,jdmus0,idmus1,jdmus1,irorz,irory
c
      Common/DIMEN/kmax,kmax1,kxmax,kymax,kzmax,nwf2,
     * maxtz1,maxjz1,maxjx1,maxjy1,maxlso,maxdso,
     * ndmu,ndthet,nacou,nrop,nrom,nrop8,nrom8,nherm,nherm38,nlegn
c
      Common/GOGINT/amu(2),xw(2),xh(2),xb(2),xm(2),WLS,t3,alpha,x0,e2
      Common/OSCLEN/bx,by,bz
c
c      write(6,*) ' PAISO ISYM = ',isym
c-------------------------------------------------------------------
      if(icall.ne.061060) then
      call inigf
      Maxmu = 2*nmax+2
      indt1 = 1
      do i=1,nmax
         do j=i,nmax
	    it1d(i,j) = indt1
	    it1d(j,i) = indt1
            Minmu  = Mod(i+j,2)+1
	    do imu=Minmu,Maxmu,2
	       t1d(indt1) = DT1(i-1,j-1,imu-1)
	       indt1 = indt1 + 1
	    end do
	 end do
      end do
      if((indt1-1).gt.KAXT1) then
         write(6,*) ' indt1 gt KAXT1 in PAISO '
	 stop
      end if
      ff = 1.0d+00/dsqrt(2.0d+00)
      Maxmu = 2*nmax+2
      indgs = 1
      do i=1,nmax
         do j=1,nmax
	    igsomu(i,j) = indgs
            Minmu  = 2-Mod(i+j,2)
	    do imu=Minmu,Maxmu,2
	       aa1 =  dsqrt(dfloat(j-1))*DT1(i-1,j-2,imu-1)
	       aa2 = -dsqrt(dfloat(j  ))*DT1(i-1,j  ,imu-1)
	       aa3 = -dsqrt(dfloat(i-1))*DT1(i-2,j-1,imu-1)
	       aa4 =  dsqrt(dfloat(i  ))*DT1(i  ,j-1,imu-1)
	       gsomu(indgs) = ff*(aa1+aa2+aa3+aa4)
	       indgs = indgs + 1
	    end do
	 end do
      end do
      if((indgs-1).gt.KAXGSO) then
         write(6,*) ' indgs gt KAXGSO in PAISO '
	 stop
      end if
      	    
      icall = 061060
      	    
      end if	    
c
c
      call timeit(0,11,'  PREDIRSO      ')
      call Prepaiso(isym,Soxyz,GROP,GROM)
      call timeit(1,11,'  PREDIRSO      ')
c
      call timeit(0,14,'REST OF GDIR    ')
c
      pi32  = (4.0d+00* datan(1.0d+00))**(1.50d+00)
      fso  = WLS/(4.0d+00*pi32*bx*by*bz)
c
c     Optimization  (THETA)
c
      do nz1=1,nzmax
         do nz2=1,nzmax

            itheeeo(nz1,nz2) = 0
            itheeoe(nz1,nz2) = 0 
            itheoee(nz1,nz2) = 0
            itheooo(nz1,nz2) = 0
	    
	  end do
       end do


      imumax = 2*nxmax
      imumay = 2*nymax
      imumaz = 2*nzmax
      
      imumixe = 1
      imumixo = 2
      imumiye = 1+isym
      imumiyo = 2-isym
      imumize = 1+isym
      imumizo = 2-isym
c
c                                      nz1 + nz2   impar 
c
      iteteeo = 0
      iteteoe = 0
      itetoee = 0
      itetooo = 0
      
      do nz1=1,nzmax
         izmin = 1+MOD((nz1+isym),2) 
         do nz2=izmin,nzmax,2

            itheeeo(nz1,nz2) = iteteeo 
            itheeoe(nz1,nz2) = iteteoe 
            itheoee(nz1,nz2) = itetoee 
            itheooo(nz1,nz2) = itetooo 

	 
            imu = 1
	    indt1 = it1d(nz1,nz2)
            do imux=imumixe,imumax,2
               do imuy=imumiye,imumay,2
	       
                  sumn1     = 0.0d+00
                  sump1     = 0.0d+00
                  sumn2     = 0.0d+00
                  sump2     = 0.0d+00
                  ii = indt1
                  do imuz=imumizo,imumaz,2
                        sumn1 = sumn1 + t1d(ii)*soxyz(imu,1) 
                        sump1 = sump1 + t1d(ii)*soxyz(imu,2)
			
                        sumn2 = sumn2 + t1d(ii)*soxyz(imu,9)
                        sump2 = sump2 + t1d(ii)*soxyz(imu,10)
			imu = imu + 1 
			ii  = ii  + 1 
                  end do                               ! imuz
c                                       even even odd      T yz    T xz
                  iteteeo = iteteeo + 1
                  the (iteteeo,1) = sumn1 ! Tyz (2)
                  the (iteteeo,2) = sump1
                  the (iteteeo,3) = sumn2 ! Txz (1)
                  the (iteteeo,4) = sump2
               end do                                  ! imuy
            end do	                               ! imux
	    
            imu = 1
	    indt1 = it1d(nz1,nz2)
            do imux=imumixo,imumax,2
               do imuy=imumiyo,imumay,2
	       
                  sumn1     = 0.0d+00
                  sump1     = 0.0d+00
                  sumn2     = 0.0d+00
                  sump2     = 0.0d+00
                  sumn3     = 0.0d+00
                  sump3     = 0.0d+00
                  sumn4     = 0.0d+00
                  sump4     = 0.0d+00
                  ii = indt1
                  do imuz=imumizo,imumaz,2
                        sumn1 = sumn1 + t1d(ii)*soxyz(imu,5) ! Txz (2)
                        sump1 = sump1 + t1d(ii)*soxyz(imu,6)
			
                        sumn2 = sumn2 + t1d(ii)*soxyz(imu,21) ! T zy (3)
                        sump2 = sump2 + t1d(ii)*soxyz(imu,22)
			
                        sumn3 = sumn3 + t1d(ii)*soxyz(imu,17) ! Tyz (1)
                        sump3 = sump3 + t1d(ii)*soxyz(imu,18)
			
                        sumn4 = sumn4 + t1d(ii)*soxyz(imu,23) ! T zx (3)
                        sump4 = sump4 + t1d(ii)*soxyz(imu,24)
			
			imu = imu + 1 
			ii  = ii  + 1 
                  end do                               ! imuz
c                           odd odd odd  Txz (2) T zy(3) Tyz(1) Tzx (3)
                  itetooo = itetooo + 1
                  the (itetooo, 5) = sumn1  ! Txz (2)
                  the (itetooo, 6) = sump1
                  the (itetooo, 7) = sumn2  ! T zy (3)
                  the (itetooo, 8) = sump2
                  the (itetooo, 9) = sumn3  ! Tyz (1)
                  the (itetooo,10) = sump3
                  the (itetooo,11) = sumn4  ! T zx (3)
                  the (itetooo,12) = sump4
               end do                                  ! imuy
            end do                                     ! imux
	                                        
            imu = 1
	    indgs = igsomu(nz1,nz2)
            do imux=imumixe,imumax,2
               do imuy=imumiyo,imumay,2
	       
                  sumn1     = 0.0d+00
                  sump1     = 0.0d+00
		  
                  ii = indgs
                  do imuz=imumize,imumaz,2
                        sumn1 = sumn1 + gsomu(ii)*soxyz(imu,11)/bz ! Txy (3)
                        sump1 = sump1 + gsomu(ii)*soxyz(imu,12)/bz
						
			imu = imu + 1 
			ii  = ii  + 1 
                  end do                               ! imuz
c                           even odd even  Txy (3) Tyx (1)
                  iteteoe = iteteoe + 1
                  the (iteteoe,13) = sumn1 ! Txy (3)
                  the (iteteoe,14) = sump1
		  
               end do                                  ! imuy
            end do                                     ! imux
	                                        
            imu = 1
	    indgs = igsomu(nz1,nz2)
            do imux=imumixo,imumax,2
               do imuy=imumiye,imumay,2
	       
                  sumn1     = 0.0d+00
                  sump1     = 0.0d+00
                  ii = indgs
                  do imuz=imumize,imumaz,2
                        sumn1 = sumn1 + gsomu(ii)*soxyz(imu,13)/bz ! Tyx (3)
                        sump1 = sump1 + gsomu(ii)*soxyz(imu,14)/bz
						
			imu = imu + 1 
			ii  = ii  + 1 
                  end do                               ! imuz
c                            odd even even  Tyx (3)
                  itetoee = itetoee + 1
                  the (itetoee,15) = sumn1 ! Tyx (3)
                  the (itetoee,16) = sump1
               end do                                  ! imuy
            end do                                     ! imux
         end do                                        ! nz2
      end do                                           ! nz1

c
c                                      nz1 + nz2   par 
c
      
                                           
      iteteoe = 0
      itetoee = 0
      itetooo = 0
      
      do nz1=1,nzmax
         izmin = 1+MOD((nz1+1-isym),2) 
         do nz2=izmin,nzmax,2
	 

            itheeoe(nz1,nz2) = iteteoe
            itheoee(nz1,nz2) = itetoee
            itheooo(nz1,nz2) = itetooo

            imu = 1
	    indt1 = it1d(nz1,nz2)
            do imux=imumixe,imumax,2
               do imuy=imumiyo,imumay,2
	       
                  sumn1     = 0.0d+00
                  sump1     = 0.0d+00
                  ii = indt1
                  do imuz=imumize,imumaz,2
                        sumn1 = sumn1 + t1d(ii)*soxyz(imu,3) ! Tzy (2)
                        sump1 = sump1 + t1d(ii)*soxyz(imu,4)
			
			imu = imu + 1 
			ii  = ii  + 1 
                  end do                               ! imuz
c                                       even odd even      T zy
                  iteteoe = iteteoe + 1
                  the (iteteoe,17) = sumn1 ! Tzy (2)
                  the (iteteoe,18) = sump1
               end do                                  ! imuy
            end do                                     ! imux
	                                        
            imu = 1
	    indt1 = it1d(nz1,nz2)
            do imux=imumixo,imumax,2
               do imuy=imumiye,imumay,2
	       
                  sumn1     = 0.0d+00
                  sump1     = 0.0d+00
                  ii = indt1
                  do imuz=imumize,imumaz,2
                        sumn1 = sumn1 + t1d(ii)*soxyz(imu,15) ! Tzx (1)
                        sump1 = sump1 + t1d(ii)*soxyz(imu,16)
			
			imu = imu + 1 
			ii  = ii  + 1 
                  end do                               ! imuz
c                           odd even even  Tzx (1) 
                  itetoee = itetoee + 1
                  the (itetoee,19) = sumn1 ! Tzx (1)
                  the (itetoee,20) = sump1
               end do                                  ! imuy
            end do                                     ! imux
	                                        
            imu = 1
	    indgs = igsomu(nz1,nz2)
            do imux=imumixo,imumax,2
               do imuy=imumiyo,imumay,2
	       
                  sumn1     = 0.0d+00
                  sump1     = 0.0d+00
                  sumn2     = 0.0d+00
                  sump2     = 0.0d+00
		  
                  ii = indgs
                  do imuz=imumizo,imumaz,2
                        sumn1 = sumn1 + gsomu(ii)*soxyz(imu,7)/bz ! Txy (2)
                        sump1 = sump1 + gsomu(ii)*soxyz(imu,8)/bz
						
                        sumn2 = sumn2 + gsomu(ii)*soxyz(imu,19)/bz ! Tyx (1)
                        sump2 = sump2 + gsomu(ii)*soxyz(imu,20)/bz
			
			imu = imu + 1 
			ii  = ii  + 1 
                  end do                               ! imuz
c                           odd odd odd  Txy (2) 
                  itetooo = itetooo + 1
                  the (itetooo,21) = sumn1 ! Txy (2)
                  the (itetooo,22) = sump1
		  
                  the (itetooo,23) = sumn2 ! Tyx (1)
                  the (itetooo,24) = sump2
               end do                                  ! imuy
            end do                                     ! imux
	                                        
         end do                                        ! nz2
      end do                                           ! nz1
      
c+-----------------------------------------------------------------+
c|                        T H E                                    |
c+-----------------------------------------------------------------+
c|                     nz1+nz2 odd  (Fields (1) and (2))           |
c+-----------------------------------------------------------------+
c|       1, 2 ......... Tyz(2)          mux even muy even muz odd  |
c|       3, 4 ......... Txz(1)                                     |
c|       5, 6 ......... Txz(2)          mux odd  muy odd  muz odd  |
c|       7, 8 ......... Tzy(3)                                     |
c|       9,10 ......... Tyz(1)                                     |
c|      11,12 ......... Tzx(3)                                     |
c|      13,14 ......... Txy(3)          mux even muy odd  muz even |
c|      15,16 ......... Tyx(3)          mux odd  muy even muz even |
c+-----------------------------------------------------------------+
c|                     nz1+nz2 even (Fields (0) and (3))           |
c+-----------------------------------------------------------------+
c|      17,18 ......... Tzy(2)          mux even muy odd  muz even |
c|      19,20 ......... Tzx(1)          mux odd  muy even muz even |
c|      21,22 ......... Txy(2)          mux odd  muy odd  muz odd  |
c|      23,24 ......... Tyx(1)                                     |
c+-----------------------------------------------------------------+
c
c
c
c
c
c      
c
c    +---------------------------------------------------------+
c    |   Starts calculation of the pairing field               |
c    +---------------------------------------------------------+
c

      irop = 0
      irom = 0
      jrop = 0
      jrom = 0
      
      deltap = 0.0d+00
      deltan = 0.0d+00
      
      do nx1=1,nxmax
        nym1 = MY(nx1)
        do nx2=1,nx1
          fx2  = 1-2*Mod((nx2+1),2)                    ! (-)**x2
          lpx  = Mod((nx1+nx2),2).eq.0                 ! type of field
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          igx  = igsomu(nx1,nx2)                        ! G index
          itx  = it1d  (nx1,nx2)                        ! T index
          do ny1 = 1,nym1
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do ny2=1,nym2
	      lpy  = Mod((ny1+ny2+isym),2).eq.0
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              igy = igsomu(ny1,ny2)                      ! G index
              ity = it1d  (ny1,ny2)                      ! T index
              fy1 = 1-2*Mod(Iabs(ny2-ny1-isym)/2,2)   ! i**(ny2-ny1  -f(-))
              fy2 = 1-2*Mod(Iabs(ny2-ny1+isym)/2,2)   ! i**(ny2-ny1+1-f(+))
              fy3 = 1-2*Mod(Iabs(ny2-ny1+1-isym)/2,2) ! i**(ny2-ny1+1-f(-))
              fy4 = 1-2*Mod(Iabs(ny2-ny1-1+isym)/2,2) ! i**(ny2-ny1  -f(+))
              lxy12 = lx12.and.(ny1.eq.ny2)
c
c---->    Positive parity
c
c
c be careful not to move the if below, the iz1e=iy1e need to be done
c
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
c
              do nz1 =nzi1e,nzm1,2
                if(lxy12) nzm2=nz1
                do nz2 = nzi2e,nzm2,2
c
c                                             F I E L D  (1)
c
            sumnx = 0.0d+00
            sumpx = 0.0d+00
            sumny = 0.0d+00
            sumpy = 0.0d+00
            sumnz = 0.0d+00
            sumpz = 0.0d+00


            if(.not.lpx.and.lpy) then       ! FIELD (1)
	    
            iteteeo = itheeeo(nz1,nz2)
            itetooo = itheooo(nz1,nz2)
            itetoee = itheoee(nz1,nz2)
	    	    
	    inddsx = igx
            do imux=imumixe,imumax,2
               indlsy = ity
               do imuy=imumiye,imumay,2
                  iteteeo = iteteeo + 1
		  aa = gsomu(inddsx)*t1d(indlsy)/bx
c  Tyz(2)		  
                  sumnx = sumnx + aa*the(iteteeo,1)		  
                  sumpx = sumpx + aa*the(iteteeo,2)		  
                  indlsy = indlsy + 1
               end do                                  ! imuy
               inddsx = inddsx + 1
            end do                                     ! imux
	    
	    sumnx = 2.0d+00*sumnx * fx2 * fy2
	    sumpx = 2.0d+00*sumpx * fx2 * fy2
	    
	    indlsx = itx
            do imux=imumixo,imumax,2
               inddsy = igy
               do imuy=imumiyo,imumay,2
                  itetooo = itetooo + 1
		  aa = t1d(indlsx)*gsomu(inddsy)/by
            sumny = sumny + aa*(the(itetooo,11)-the(itetooo,5))
            sumpy = sumpy + aa*(the(itetooo,12)-the(itetooo,6))
                  inddsy = inddsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	    
	    sumny = 2.0d+00*sumny * fx2 * fy2
	    sumpy = 2.0d+00*sumpy * fx2 * fy2
	    
	    indlsx = itx
            do imux=imumixo,imumax,2
               indlsy = ity
               do imuy=imumiye,imumay,2
                  itetoee = itetoee + 1
		  aa = t1d(indlsx)*t1d(indlsy)
                  sumnz = sumnz - aa*the(itetoee,15)		  
                  sumpz = sumpz - aa*the(itetoee,16)		  
                  indlsy = indlsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	                                        
	    sumnz = 2.0d+00 * sumnz * fx2 * fy2
	    sumpz = 2.0d+00 * sumpz * fx2 * fy2
c
c                                             F I E L D  (2)
c
            else if(lpx.and..not.lpy) then       ! FIELD (2)
	    
            iteteeo = itheeeo(nz1,nz2)
            itetooo = itheooo(nz1,nz2)
            iteteoe = itheeoe(nz1,nz2)
	    
	    inddsx = igx
            do imux=imumixo,imumax,2
               indlsy = ity
               do imuy=imumiyo,imumay,2
                  itetooo = itetooo + 1
		  aa = gsomu(inddsx)*t1d(indlsy)/bx
                  sumnx = sumnx + aa*(the(itetooo, 9)*fy3
     &	                 -fy4*the(itetooo,7))
                  sumpx = sumpx + aa*(the(itetooo,10)*fy3
     &		         -fy4*the(itetooo,8))
                  indlsy = indlsy + 1
               end do                                  ! imuy
               inddsx = inddsx + 1
            end do                                     ! imux
	                                         
            sumnx = 2.0d+00*sumnx * fx2
            sumpx = 2.0d+00*sumpx * fx2
	    
	    indlsx = itx
            do imux=imumixe,imumax,2
               inddsy = igy
               do imuy=imumiye,imumay,2
                  iteteeo = iteteeo + 1
		  aa = t1d(indlsx)*gsomu(inddsy)/by
                  sumny = sumny - aa*the(iteteeo,3)
                  sumpy = sumpy - aa*the(iteteeo,4)
                  inddsy = inddsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	                                        
            sumny = 2.0d+00*sumny * fx2 * fy3
            sumpy = 2.0d+00*sumpy * fx2 * fy3
	    
	    indlsx = itx
            do imux=imumixe,imumax,2
               indlsy = ity
               do imuy=imumiyo,imumay,2
                  iteteoe = iteteoe + 1
		  aa = t1d(indlsx)*t1d(indlsy)
                  sumnz = sumnz + aa*the(iteteoe,13)		  
                  sumpz = sumpz + aa*the(iteteoe,14)		  
                  indlsy = indlsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	                                        
            sumnz = sumnz * fx2 * fy4 * 2.0d+00
            sumpz = sumpz * fx2 * fy4 * 2.0d+00 
	    
c
c                                             F I E L D  (3)
c
            else if(.not.lpx.and..not.lpy) then       ! FIELD (3)
	    
            itetoee = itheoee(nz1,nz2)
            itetooo = itheooo(nz1,nz2)
            iteteoe = itheeoe(nz1,nz2)
	    
	    inddsx = igx
            do imux=imumixe,imumax,2
               indlsy = ity
               do imuy=imumiyo,imumay,2
                  iteteoe = iteteoe + 1
		  aa = gsomu(inddsx)*t1d(indlsy)/bx
                  sumnx = sumnx - aa*the(iteteoe,17)
                  sumpx = sumpx - aa*the(iteteoe,18)
                  indlsy = indlsy + 1
               end do                                  ! imuy
               inddsx = inddsx + 1
            end do                                     ! imux
	    
            sumnx = -sumnx  * fy4 * 2.0d+00
            sumpx = -sumpx  * fy4 * 2.0d+00
	    
	    indlsx = itx
            do imux=imumixo,imumax,2
               inddsy = igy
               do imuy=imumiye,imumay,2
                  itetoee = itetoee + 1
		  aa = t1d(indlsx)*gsomu(inddsy)/by
                  sumny = sumny + aa*the(itetoee,19)
                  sumpy = sumpy + aa*the(itetoee,20)
                  inddsy = inddsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	    
            sumny = -sumny  * fy3 * 2.0d+00
            sumpy = -sumpy  * fy3 * 2.0d+00
	    
	    indlsx = itx
            do imux=imumixo,imumax,2
               indlsy = ity
               do imuy=imumiyo,imumay,2
                  itetooo = itetooo + 1
		  aa = t1d(indlsx)*t1d(indlsy)
          sumnz = sumnz + aa*(the(itetooo,21)*fy4-the(itetooo,23)*fy3)
          sumpz = sumpz + aa*(the(itetooo,22)*fy4-the(itetooo,24)*fy3)
                  indlsy = indlsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux 
	                                        
            sumnz = -sumnz  * 2.0d+00
            sumpz = -sumpz  * 2.0d+00
c
           end if  ! fields
c
                irop = irop + 1
	        dpp(irop) = dpp(irop) + fso*(sumpx + sumpy + sumpz ) 
	        dnp(irop) = dnp(irop) + fso*(sumnx + sumny + sumnz )
		
                jrop  = jrop + 1
                akan = GROP(jrop+4+isym) ! kappa 
                akap = GROP(jrop+6+isym) ! kappa 
                jrop   = jrop + 7
                lxyz = lxy12.and.(nz2.eq.nz1)
                if(lxyz) then
                    akan = 0.5d+00*akan
                    akap = 0.5d+00*akap
                end if
		
		deltap = deltap + fso*(sumpx + sumpy + sumpz )*akap
		deltan = deltan + fso*(sumnx + sumny + sumnz )*akan
		 
                end do   ! nz2
              end do     ! nz1
              end if    ! lmze
c
c                                           ---->    Negative parity
c
c
c be careful not to move the if below, the iz1o=iy1o need to be done
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
c
              do nz1 =nzi1o,nzm1,2
                if(lxy12) nzm2=nz1
                do nz2 = nzi2o,nzm2,2
		
c
            sumnx = 0.0d+00
            sumpx = 0.0d+00
            sumny = 0.0d+00
            sumpy = 0.0d+00
            sumnz = 0.0d+00
            sumpz = 0.0d+00


            if(.not.lpx.and.lpy) then       ! FIELD (1)
	    
            iteteeo = itheeeo(nz1,nz2)
            itetooo = itheooo(nz1,nz2)
            itetoee = itheoee(nz1,nz2)
	    	    
	    inddsx = igx
            do imux=imumixe,imumax,2
               indlsy = ity
               do imuy=imumiye,imumay,2
                  iteteeo = iteteeo + 1
		  aa = gsomu(inddsx)*t1d(indlsy)/bx
c  Tyz(2)		  
                  sumnx = sumnx + aa*the(iteteeo,1)		  
                  sumpx = sumpx + aa*the(iteteeo,2)		  
                  indlsy = indlsy + 1
               end do                                  ! imuy
               inddsx = inddsx + 1
            end do                                     ! imux
	    
	    sumnx = 2.0d+00*sumnx * fx2 * fy2
	    sumpx = 2.0d+00*sumpx * fx2 * fy2
	    
	    indlsx = itx
            do imux=imumixo,imumax,2
               inddsy = igy
               do imuy=imumiyo,imumay,2
                  itetooo = itetooo + 1
		  aa = t1d(indlsx)*gsomu(inddsy)/by
            sumny = sumny + aa*(the(itetooo,11)-the(itetooo,5))
            sumpy = sumpy + aa*(the(itetooo,12)-the(itetooo,6))
                  inddsy = inddsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	    
	    sumny = 2.0d+00*sumny * fx2 * fy2
	    sumpy = 2.0d+00*sumpy * fx2 * fy2
	    
	    indlsx = itx
            do imux=imumixo,imumax,2
               indlsy = ity
               do imuy=imumiye,imumay,2
                  itetoee = itetoee + 1
		  aa = t1d(indlsx)*t1d(indlsy)
                  sumnz = sumnz - aa*the(itetoee,15)		  
                  sumpz = sumpz - aa*the(itetoee,16)		  
                  indlsy = indlsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	                                        
	    sumnz = 2.0d+00 * sumnz * fx2 * fy2
	    sumpz = 2.0d+00 * sumpz * fx2 * fy2
c
c                                             F I E L D  (2)
c
            else if(lpx.and..not.lpy) then       ! FIELD (2)
	    
            iteteeo = itheeeo(nz1,nz2)
            itetooo = itheooo(nz1,nz2)
            iteteoe = itheeoe(nz1,nz2)
	    
	    inddsx = igx
            do imux=imumixo,imumax,2
               indlsy = ity
               do imuy=imumiyo,imumay,2
                  itetooo = itetooo + 1
		  aa = gsomu(inddsx)*t1d(indlsy)/bx
                  sumnx = sumnx + aa*(the(itetooo, 9)*fy3
     &	                 -fy4*the(itetooo,7))
                  sumpx = sumpx + aa*(the(itetooo,10)*fy3
     &		         -fy4*the(itetooo,8))
                  indlsy = indlsy + 1
               end do                                  ! imuy
               inddsx = inddsx + 1
            end do                                     ! imux
	                                         
            sumnx = 2.0d+00*sumnx * fx2
            sumpx = 2.0d+00*sumpx * fx2
	    
	    indlsx = itx
            do imux=imumixe,imumax,2
               inddsy = igy
               do imuy=imumiye,imumay,2
                  iteteeo = iteteeo + 1
		  aa = t1d(indlsx)*gsomu(inddsy)/by
                  sumny = sumny - aa*the(iteteeo,3)
                  sumpy = sumpy - aa*the(iteteeo,4)
                  inddsy = inddsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	                                        
            sumny = 2.0d+00*sumny * fx2 * fy3
            sumpy = 2.0d+00*sumpy * fx2 * fy3
	    
	    indlsx = itx
            do imux=imumixe,imumax,2
               indlsy = ity
               do imuy=imumiyo,imumay,2
                  iteteoe = iteteoe + 1
		  aa = t1d(indlsx)*t1d(indlsy)
                  sumnz = sumnz + aa*the(iteteoe,13)		  
                  sumpz = sumpz + aa*the(iteteoe,14)		  
                  indlsy = indlsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	                                        
            sumnz = sumnz * fx2 * fy4 * 2.0d+00
            sumpz = sumpz * fx2 * fy4 * 2.0d+00 
	    
c
c                                             F I E L D  (3)
c
            else if(.not.lpx.and..not.lpy) then       ! FIELD (3)
	    
            itetoee = itheoee(nz1,nz2)
            itetooo = itheooo(nz1,nz2)
            iteteoe = itheeoe(nz1,nz2)
	    
	    inddsx = igx
            do imux=imumixe,imumax,2
               indlsy = ity
               do imuy=imumiyo,imumay,2
                  iteteoe = iteteoe + 1
		  aa = gsomu(inddsx)*t1d(indlsy)/bx
                  sumnx = sumnx - aa*the(iteteoe,17)
                  sumpx = sumpx - aa*the(iteteoe,18)
                  indlsy = indlsy + 1
               end do                                  ! imuy
               inddsx = inddsx + 1
            end do                                     ! imux
	    
            sumnx = -sumnx  * fy4 * 2.0d+00
            sumpx = -sumpx  * fy4 * 2.0d+00
	    
	    indlsx = itx
            do imux=imumixo,imumax,2
               inddsy = igy
               do imuy=imumiye,imumay,2
                  itetoee = itetoee + 1
		  aa = t1d(indlsx)*gsomu(inddsy)/by
                  sumny = sumny + aa*the(itetoee,19)
                  sumpy = sumpy + aa*the(itetoee,20)
                  inddsy = inddsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	    
            sumny = -sumny  * fy3 * 2.0d+00
            sumpy = -sumpy  * fy3 * 2.0d+00
	    
	    indlsx = itx
            do imux=imumixo,imumax,2
               indlsy = ity
               do imuy=imumiyo,imumay,2
                  itetooo = itetooo + 1
		  aa = t1d(indlsx)*t1d(indlsy)
          sumnz = sumnz + aa*(the(itetooo,21)*fy4-the(itetooo,23)*fy3)
          sumpz = sumpz + aa*(the(itetooo,22)*fy4-the(itetooo,24)*fy3)
                  indlsy = indlsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux 
	                                        
            sumnz = -sumnz  * 2.0d+00
            sumpz = -sumpz  * 2.0d+00
c
           end if  ! fields
c
                irom = irom + 1
	        dpm(irom) = dpm(irom) + fso*(sumpx + sumpy + sumpz ) 
	        dnm(irom) = dnm(irom) + fso*(sumnx + sumny + sumnz )
	   
                jrom  = jrom + 1
                akan = GROM(jrom+4+isym) ! kappa 
                akap = GROM(jrom+6+isym) ! kappa 
                jrom   = jrom + 7
                lxyz = lxy12.and.(nz2.eq.nz1)
                if(lxyz) then
                    akan = 0.5d+00*akan
                    akap = 0.5d+00*akap
                end if
		
		deltap = deltap + fso*(sumpx + sumpy + sumpz )*akap
		deltan = deltan + fso*(sumnx + sumny + sumnz )*akan
		 

                end do  ! nz2
              end do    ! nz1
              end if    ! lmzo
            end do  ! ny2
          end do    ! ny1
        end do      ! nx2
      end do        ! nx1
      
      fsym = dfloat( (-1)**isym )
      deltap = fsym*deltap
      deltan = fsym*deltan
      write(6,*) ' Spin-orbit pairing ISYM ',isym,deltap,deltan
      
      call timeit(1,14,'REST OF GDIR    ')
      return
      end
      Subroutine CONSOP (COP)
      Implicit real*8 (A-H,O-Z)
      Parameter (NCM=06,NLAM=NCM-3)          ! maximum number of constraints
      Dimension COP(*)                       ! (NCM-2)*(NP2+NM2)
c COP:   Jx (NP*NP)  Jx (NM*NM)
c        Qlm(NP*NP)  Qlm(NM*NM)     l-> ILamb  m->Imm
c
      Common /CCONS/ cval(ncm),prec(ncm),IC(NCM),
     & ILAMB(NLAM),IMM(NLAM),ISOS(NCM,2)
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
     
      NP2 = ND2(1)
      NM2 = ND2(2)
      
      if(IC(3).ne.0) call momang('X',COP(1),COP(1+NP2))
      
      do k=1,NLAM
         if(IC(k+3).ne.0) then
	    ICOP = k*(NP2+NM2) + 1
	    ICOM = ICOP + NP2
	    if(ILAMB(k).eq.10) then
               call NECKME(COP(ICOP),COP(ICOM))
	    else
               call QLMME(ILAMB(k),IMM(k),COP(ICOP),COP(ICOM),ISI)
c	       write(6,*) ' Operador Qlm ',Ilamb(k),Imm(k),isi
	    end if
	 end if
      end do
      return
      end
c ---------------------------------------------------------------------------      
      Subroutine CONSMV (RO,COP,rval)
      Implicit real*8 (A-H,O-Z)
      Parameter (NCM=06,NLAM=NCM-3)          ! maximum number of constraints
      Dimension COP(*)
      Dimension RO(*)
      Dimension rval(NCM)
      
      Common /CCONS/ cval(ncm),prec(ncm),IC(NCM),
     & ILAMB(NLAM),IMM(NLAM),ISOS(NCM,2)
     
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC1/ NNU(4),NNV(4),NNRO(4),NNKA(4)
      
      do k=1,NCM
         rval(k) = 0.0d+00
      end do
      
      NP2 = ND2(1)
      NM2 = ND2(2)
      
      do it=1,4
      
	 IRO1 = NNRO(it)
	 IRO2 = NNRO(it) + ND2(it)
	 
	 rval(it/3+1) = rval(it/3+1) + trace(RO(IRO1),ND(it))+
     &                                 trace(RO(IRO2),ND(it))
c
c                         < Jx > (time odd -> the minus sign) 
c	 
         if(IC(3).ne.0) then
         ICOP = (1-Mod(it,2))*NP2+1
         rval(3) = rval(3) + ddot(ND2(it),RO(IRO1),1,COP(ICOP),1)-
     &                       ddot(ND2(it),RO(IRO2),1,COP(ICOP),1)
         end if
	 
c
c                         < Qlm > (time even -> the plus sign) 
c	 
         do k=1,NLAM
	 
         if(IC(k+3).ne.0) then
         ICOP = (1-Mod(it,2))*NP2+k*(NP2+NM2)+1
         rval(3+k) = rval(3+k) + ddot(ND2(it),RO(IRO1),1,COP(ICOP),1)+
     &                           ddot(ND2(it),RO(IRO2),1,COP(ICOP),1)
         end if
	 
	 end do ! k
	 
      end do ! it
      return
      end
c---------------------------------------------------------------------------      
      Subroutine CONS20 (UV,COP,C20,AUX)
      Implicit real*8 (A-H,O-Z)
      Parameter (NCM=06,NLAM=NCM-3)          ! maximum number of constraints
      Dimension COP (*)
      Dimension C20 (*)
      Dimension UV  (*)
      Dimension AUX (*)
      
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC1/ NNU(4),NNV(4),NNRO(4),NNKA(4)
      Common /ITLOC3/ NNC20(4,NCM),NNC11(4,NCM)
      
      Common /CCONS/ cval(ncm),prec(ncm),IC(NCM),
     & ILAMB(NLAM),IMM(NLAM),ISOS(NCM,2)
     
      NP2 = ND2(1)
      NM2 = ND2(2)
     
      do it=1,4
      
         IU = NNU(it)
         IV = NNV(it)
	 N  = ND(it)
	 N1 = ND(it)-NBLOCK(it)
	 N2 = ND(it)+NBLOCK(it)
	 ISO= it/3+1 ! 1 Protones 2 Neutrones

c         write(6,*) ' CONS20 N ',it,IU,IV,N,N1,N2,ISO,NNC20(it,ISO)	 
         call SN20 (UV(IU),UV(IV),N1,N2,N,C20(NNC20(it,ISO)))
      
         if(IC(3).ne.0) then 
           ICOP = (1-Mod(it,2))*NP2+1
c         write(6,*) ' CONS20 Jx ',it,IU,IV,N,N1,N2,ISO,NNC20(it,3)	 
           call SO20P(UV(IU),UV(IV),N1,N2,N,
     &                COP(ICOP),AUX,C20(NNC20(it,3)),-1,1)
         end if
	 
         do k=1,NLAM
	 
         if(IC(k+3).ne.0) then
           ICOP = (1-Mod(it,2))*NP2+k*(NP2+NM2)+1
c         write(6,*) ' CONS20 Qlm ',it,IU,IV,N,N1,N2,ISO,NNC20(it,k+3)	 
           call SO20P(UV(IU),UV(IV),N1,N2,N,
     &                COP(ICOP),AUX,C20(NNC20(it,k+3)),1,1)
         end if
	 
	 end do ! k

      end do ! it 
      
      return
      end 
c---------------------------------------------------------------------------      
      Subroutine CONS11 (UV,COP,C11,AUX)
      Implicit real*8 (A-H,O-Z)
      Parameter (NCM=06,NLAM=NCM-3)          ! maximum number of constraints
      Dimension COP (*)
      Dimension C11 (*)
      Dimension UV  (*)
      Dimension AUX (*)
      
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC1/ NNU(4),NNV(4),NNRO(4),NNKA(4)
      Common /ITLOC3/ NNC20(4,NCM),NNC11(4,NCM)
      
      Common /CCONS/ cval(ncm),prec(ncm),IC(NCM),
     & ILAMB(NLAM),IMM(NLAM),ISOS(NCM,2)
     
      NP2 = ND2(1)
      NM2 = ND2(2)
           
      do it=1,4
      
         IU = NNU(it)
         IV = NNV(it)
	 N  = ND(it)
	 N1 = ND(it)-NBLOCK(it)
	 N2 = ND(it)+NBLOCK(it)
	 ISO= it/3+1 ! 1 Protones 2 Neutrones
	 
	 I11_1= NNC11(it,ISO)
	 I11_2= I11_1 + N1*N1
	 
         call SN11 (UV(IU),UV(IV),N1,N2,N,C11(I11_1),C11(I11_2))
      
         if(IC(3).ne.0) then 
	   I11_1= NNC11(it,3)
	   I11_2= I11_1 + N1*N1
           ICOP = (1-Mod(it,2))*NP2+1
           call SO11P(UV(IU),UV(IV),N1,N2,N,
     &                COP(ICOP),AUX,C11(I11_1),C11(I11_2),-1,1)
         end if
	 
         do k=1,NLAM
	 
         if(IC(k+3).ne.0) then
           ICOP = (1-Mod(it,2))*NP2+k*(NP2+NM2)+1
	   I11_1= NNC11(it,k+3)
	   I11_2= I11_1 + N1*N1
           call SO11P(UV(IU),UV(IV),N1,N2,N,
     &                COP(ICOP),AUX,C11(I11_1),C11(I11_2),1,1)
         end if
	 
	 end do ! k

      end do ! it 
      
      return
      end 
c---------------------------------------------------------------------------      
      Subroutine LAMBDA0 (C20,H20,dd,ix,clam,cdlam,shhc)
      Implicit real*8 (A-H,O-Z)
      Parameter (NCM=06,NYV=20*NCM,NLAM=NCM-3) ! maximum number of constraints
      Dimension C20(*),H20(*)
      Dimension DD(NCM)
      Dimension IX(NCM),CLAM(NCM),CDLAM(NCM)
c
      Dimension IPIV(ncm)
      Dimension scn(ncm,ncm),scna(ncm,ncm),hcn(ncm),hcna(ncm)
      Dimension work(NYV)
c
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC3/ NNC20(4,NCM),NNC11(4,NCM)
      Common /ITLOC4/ NNH20(4),NNH11(4),NNEQP(4)
c
      Common /CCONS/ cval(ncm),prec(ncm),IC(NCM),
     & ILAMB(NLAM),IMM(NLAM),ISOS(NCM,2)
c
      ioutlev = 0
      
      do i=1,ncm
         do j=i,ncm
            scn (i,j) = 0.0d+00
            scna(i,j) = 0.0d+00
         end do
         scn(i,i) = 1.0d+00
         hcn(i)   = 0.0d+00
         hcna(i)  = 0.0d+00
         ix (i)   = ic(i)
      end do
c
         do i=1,ncm
            do j=i,ncm
	    
               if(IX(i)*IX(j).eq.1) then
	       
                  FFF = 0.0d+00
                  do it=1,4
      
                    IUVI = NNC20(it,i)
                    IUVJ = NNC20(it,j)
		    NN   = ND2(it)-NBLOCK(it)**2  ! N1*N2
                    FFF = FFF + DDOT(NN,C20(IUVI),1,C20(IUVJ),1)
		    
                  end do
		  
               if(i.eq.j) then
                  if(FFF.gt.1.d-9) then
                     SCN(i,i) = FFF
                     FFH = 0.0d+00
                     do it=1,4
      
                        IUVI = NNC20(it,i)
                        IH20 = NNH20(it)
		        NN   = ND2(it)-NBLOCK(it)**2  ! N1*N2
                        FFH=FFH+DDOT(NN,C20(IUVI),1,H20(IH20),1)
		    
                     end do
                     HCN(i)   = FFH
                  else
                     IX(i)    = 0
            write(6,*) ' Norm zero for the gradient of constraint ',I
                  end if
               else
                  SCN(i,j) = FFF
               end if
            end if                       ! IC(J)*IC(I)
         end do                          !    J
      end do                             !    I

c======================================================================
      do i=1,ncm
         do j=i,ncm
            scna(i,j) = scn(i,j)
         end do
         hcna(i) = hcn(i)
         if(ioutlev.ge.4) write(6,600) i,hcn(i),(scn(k,i),k=1,i)
      end do
  600 format ( ' HCN  ',i2,1x,d12.6,' SCN ',8d13.6)
c+=====================================================================+
c||            C H E M I C A L      P O T E N T I A L S               ||
c+=====================================================================+
c
      call DSYSV('U',ncm,1,SCNA,ncm,IPIV,hcn,ncm,work,nyv,info)
      call dcopy(ncm,hcn,1,clam,1)
C========================================================= GRADIENT NORM
      shhc =  - DDot(ncm,clam,1,hcna,1)
c+=====================================================================+
c|| ADJUSTING CHEMICAL POTENTIALS TO YIELD THE DESIRED VALUE OF CONST.||
c+=====================================================================+
c
      if(ioutlev.ge.4) write(6,601) (DD(i),i=1,ncm)
  601 format ( ' DD   ',8d13.6)
c
c
      do i=1,ncm
         do j=i,ncm
            scna(i,j) = scn(i,j)
         end do
         hcna(i) = dd (i)
      end do
c
      call DSYSV('U',ncm,1,SCNA,ncm,IPIV,hcna,ncm,work,nyv,info)
      call dcopy(ncm,hcna,1,cdlam,1)
c
      return 
      end
c---------------------------------------------------------------------------      
      Subroutine LAMBDA1 (C20,H20,QUOT,dd,ix,clam,cdlam,shhc)
      Implicit real*8 (A-H,O-Z)
      Parameter (NCM=06,NYV=20*NCM,NLAM=NCM-3) ! maximum number of constraints
      Dimension C20(*),H20(*),QUOT(*)
      Dimension DD(NCM)
      Dimension IX(NCM),CLAM(NCM),CDLAM(NCM)
c
      Dimension IPIV(ncm)
      Dimension scn(ncm,ncm),scna(ncm,ncm),hcn(ncm),hcna(ncm)
      Dimension work(NYV)
c
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC3/ NNC20(4,NCM),NNC11(4,NCM)
      Common /ITLOC4/ NNH20(4),NNH11(4),NNEQP(4)
c
      Common /CCONS/ cval(ncm),prec(ncm),IC(NCM),
     & ILAMB(NLAM),IMM(NLAM),ISOS(NCM,2)
c
      ioutlev = 0
      
      do i=1,ncm
         do j=i,ncm
            scn (i,j) = 0.0d+00
            scna(i,j) = 0.0d+00
         end do
         scn(i,i) = 1.0d+00
         hcn(i)   = 0.0d+00
         hcna(i)  = 0.0d+00
         ix (i)   = ic(i)
      end do
c
         do i=1,ncm
            do j=i,ncm
	    
               if(IX(i)*IX(j).eq.1) then
	       
                  FFF = 0.0d+00
                  do it=1,4
      
                    IUVI = NNC20(it,i)
                    IUVJ = NNC20(it,j)
		    IQ   = NNH20(it)
		    NN   = ND2(it)-NBLOCK(it)**2  ! N1*N2
                    FFF = FFF + RMM(NN,C20(IUVI),C20(IUVJ),QUOT(IQ))
		    
                  end do
		  
               if(i.eq.j) then
                  if(FFF.gt.1.d-9) then
                     SCN(i,i) = FFF
                     FFH = 0.0d+00
                     do it=1,4
      
                        IUVI = NNC20(it,i)
                        IH20 = NNH20(it)
		        NN   = ND2(it)-NBLOCK(it)**2  ! N1*N2
                        FFH=FFH+RMM(NN,C20(IUVI),H20(IH20),QUOT(IH20))
		    
                     end do
                     HCN(i)   = FFH
                  else
                     IX(i)    = 0
            write(6,*) ' Norm zero for the gradient of constraint ',I
                  end if
               else
                  SCN(i,j) = FFF
               end if
            end if                       ! IC(J)*IC(I)
         end do                          !    J
      end do                             !    I

c======================================================================
      do i=1,ncm
         do j=i,ncm
            scna(i,j) = scn(i,j)
         end do
         hcna(i) = hcn(i)
         if(ioutlev.ge.4) write(6,600) i,hcn(i),(scn(k,i),k=1,i)
      end do
  600 format ( ' HCN  ',i2,1x,d12.6,' SCN ',8d13.6)
c+=====================================================================+
c||            C H E M I C A L      P O T E N T I A L S               ||
c+=====================================================================+
c
      call DSYSV('U',ncm,1,SCNA,ncm,IPIV,hcn,ncm,work,nyv,info)
      call dcopy(ncm,hcn,1,clam,1)
C========================================================= GRADIENT NORM
      shhc =  - DDot(ncm,clam,1,hcna,1)
c+=====================================================================+
c|| ADJUSTING CHEMICAL POTENTIALS TO YIELD THE DESIRED VALUE OF CONST.||
c+=====================================================================+
c
      if(ioutlev.ge.4) write(6,601) (DD(i),i=1,ncm)
  601 format ( ' DD   ',8d13.6)
c
c
      do i=1,ncm
         do j=i,ncm
            scna(i,j) = scn(i,j)
         end do
         hcna(i) = dd (i)
      end do
c
      call DSYSV('U',ncm,1,SCNA,ncm,IPIV,hcna,ncm,work,nyv,info)
      call dcopy(ncm,hcna,1,cdlam,1)
      call dscal(ncm,0.65d+00,cdlam,1)
c
      return 
      end
c---------------------------------------------------------------------------      
      Subroutine DLAMBDA (C20,dd,ix,cdlam)
      Implicit real*8 (A-H,O-Z)
      Parameter (NCM=06,NYV=20*NCM,NLAM=NCM-3) ! maximum number of constraints
      Dimension C20(*)
      Dimension DD(NCM),IX(NCM),CDLAM(NCM)
c
      Dimension IPIV(ncm)
      Dimension scn(ncm,ncm),scna(ncm,ncm),hcna(ncm)
      Dimension work(NYV)
c
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC3/ NNC20(4,NCM),NNC11(4,NCM)
      Common /ITLOC4/ NNH20(4),NNH11(4),NNEQP(4)
c
      Common /CCONS/ cval(ncm),prec(ncm),IC(NCM),
     & ILAMB(NLAM),IMM(NLAM),ISOS(NCM,2)
c
      do i=1,ncm
         do j=i,ncm
            scn (i,j) = 0.0d+00
            scna(i,j) = 0.0d+00
         end do
         scn(i,i) = 1.0d+00
         hcna(i)  = 0.0d+00
         ix (i)   = ic(i)
      end do
c
         do i=1,ncm
            do j=i,ncm
	    
               if(IX(i)*IX(j).eq.1) then
	       
                  FFF = 0.0d+00
                  do it=1,4
      
                    IUVI = NNC20(it,i)
                    IUVJ = NNC20(it,j)
		    NN   = ND2(it)
                    FFF = FFF + DDOT(NN,C20(IUVI),1,C20(IUVJ),1)
		    
                  end do
		  
               if(i.eq.j) then
                  if(FFF.gt.1.d-9) then
                     SCN(i,i) = FFF
                  else
                     IX(i)    = 0
            write(6,*) ' Norm zero for the gradient of constraint ',I
                  end if
               else
                  SCN(i,j) = FFF
               end if
            end if                       ! IC(J)*IC(I)
         end do                          !    J
      end do                             !    I

c
c+=====================================================================+
c|| ADJUSTING CHEMICAL POTENTIALS TO YIELD THE DESIRED VALUE OF CONST.||
c+=====================================================================+
c
c
      ioutlev = 0
      do i=1,ncm
         do j=i,ncm
            scna(i,j) = scn(i,j)
         end do
         if(ioutlev.ge.4) write(6,600) i,dd(i),(scn(k,i),k=1,i)
  600 format ( ' DD  ',i2,1x,d12.6,' SCN ',8d13.6)
         hcna(i) = dd (i)
      end do
c
      call DSYSV('U',ncm,1,SCNA,ncm,IPIV,hcna,ncm,work,nyv,info)
      call dcopy(ncm,hcna,1,cdlam,1)
      call dscal(ncm,0.5d+00,cdlam,1)
c
      return 
      end
      Subroutine HPRIME(it,Gamma,COP,IX,CLAM)
      Implicit real*8 (A-H,O-Z)
      Parameter (NCM=06,NYV=20*NCM,NLAM=NCM-3) ! maximum number of constraints
      Dimension Gamma(*)
      Dimension COP(*)
      Dimension IX(NCM),CLAM(NCM)
c
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      
      iso = it/3
      N   = ND (it)
      N2  = ND2(it)
      if(iso.eq.0) then
         NP2 = ND2(1)
	 NM2 = ND2(2)
      else
         NP2 = ND2(3)
	 NM2 = ND2(4)
      end if
      
      if((ix(1).ne.0).and.(iso.eq.0)) then ! protons
      
         do i=1,N
	    I1 = i + (i-1)*N
	    I2 = i + (i-1)*N + N2
	    Gamma(I1) = Gamma(I1) - clam(1)
	    Gamma(I2) = Gamma(I2) - clam(1)
	 end do 
	 
      end if
      
      if((ix(2).ne.0).and.(iso.eq.1)) then ! neutrons
      
         do i=1,N
	    I1 = i + (i-1)*N
	    I2 = i + (i-1)*N + N2
	    Gamma(I1) = Gamma(I1) - clam(2)
	    Gamma(I2) = Gamma(I2) - clam(2)
	 end do 
	 
      end if
      
      if(ix(3).ne.0) then ! Jx
         aa1 = -clam(3)
	 aa2 =  clam(3)
	 
	 if(Mod(it,2).eq.1) then
	    ICOP = 1
	 else
	    ICOP = 1 + NP2
	 end if
	 
         Call DAXPY(N2,aa1,COP(ICOP),1,GAMMA,1)
         Call DAXPY(N2,aa2,COP(ICOP),1,GAMMA(1+N2),1)
	 
      end if 

      do k=1,NLAM
      
     	 if(ix(k+3).ne.0) then
     	    aa1 = -clam(k+3)
  	    aa2 = -clam(k+3)
 
  	    if(Mod(it,2).eq.1) then
  	       ICOP = k*(NP2+NM2) + 1
  	    else
  	       ICOP = k*(NP2+NM2) + 1 + NP2
  	    end if
 
     	    Call DAXPY(N2,aa1,COP(ICOP),1,GAMMA,1)
     	    Call DAXPY(N2,aa2,COP(ICOP),1,GAMMA(1+N2),1)
 
     	 end if
      end do
 
      return
      end
      Subroutine ZETA0(it,C20,Z0,IX,ETC)
      Implicit real*8 (A-H,O-Z)
      Parameter (NCM=06,NYV=20*NCM,NLAM=NCM-3) ! maximum number of constraints
      Dimension C20(*),Z0(*)
      Dimension IX(NCM),ETC(NCM)
c
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC3/ NNC20(4,NCM),NNC11(4,NCM)
      
      isoi = it/3 + 1
      N   = ND (it)
      NN  = ND2(it)
      
      do k=1,NN
         Z0(k) = 0.0d+00
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
      Subroutine DOCONS (COP,ROM,C20,UV,IX,dd,Z0,AUX)
      Implicit real*8 (A-H,O-Z)
      Parameter (NCM=06,NYV=20*NCM,NLAM=NCM-3) ! maximum number of constraints
      Dimension C20(*),UV(*)
      Dimension ROM(*),COP(*)
      Dimension Z0(*),AUX(*)
      
      Dimension IX(NCM),DD(NCM),cdlam(NCM),rval(NCM)
      
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC1/ NNU(4),NNV(4),NNRO(4),NNKA(4)
      Common /ITLOC4/ NNH20(4),NNH11(4),NNEQP(4)
      Common /CCONS/ cval(ncm),prec(ncm),IC(NCM),
     & ILAMB(NLAM),IMM(NLAM),ISOS(NCM,2)
     
      do iter=1,4

         do it=1,4
            N  = ND (it)
            call RO(UV(NNV(it)),N,N,N,ROM(NNRO(it)))
         end do

         call CONSMV (ROM,COP,rval)

         do i=1,ncm
             cdlam(i) = 0.0d+00
            if(ic(i).ne.0) then 
               dd(i) = cval(i)-rval(i)
            else
	       dd(i) = 0.0d+00
	    end if
         end do
c      write(6,'(" RVAL   ",6d13.6)') (RVAL(ii),ii=1,ncm)


         write(6,601) iter,(DD(i),i=1,ncm)

         call CONS20 (UV,COP,C20,AUX)

         call DLAMBDA (C20,dd,ix,cdlam)

c      write(6,'(" CDLAM  ",6d13.6)') (CDLAM(ii),ii=1,ncm)
         do it=1,4
c
           N  = ND (it)
	   N1 = N - NBLOCK(it)
	   N2 = N + NBLOCK(it)
	   IU = NNU(it)
	   IV = NNV(it)
	   IZ = NNH20(it)
c+---------------------------------------------------------------------+
c|            T H O U L E S S     Z     M A T R I X                    |
c+---------------------------------------------------------------------+
	   itt =it
           call ZETA0(itt,C20,Z0(IZ),IX,cdlam)
c
c+---------------------------------------------------------------------+
c|            N E W    W A V E    F U N C T I O N S                    |
c+---------------------------------------------------------------------+
c
           call NEWUV(UV(IU),UV(IV),N1,N2,N,Z0(IZ),AUX)
	  
	 end do ! it

      end do ! iter
c
c     Compute again after finishing the calculation
c
      do it=1,4
         N  = ND (it)
         call RO(UV(NNV(it)),N,N,N,ROM(NNRO(it)))
      end do

      call CONSMV (ROM,COP,rval)


      do i=1,ncm
         if(ic(i).ne.0) then 
            dd(i) = cval(i)-rval(i)
         end if
      end do

      call CONS20 (UV,COP,C20,AUX)

      write(6,601) iter,(DD(i),i=1,ncm)
  601 format ( ' DD   ',i3,8d13.6)
      return
      end
c---------------------------------------------------------------------------      
      Subroutine COLLMas (C20,EQP,ix,AMMN)
      Implicit real*8 (A-H,O-Z)
      Parameter (NCM=06,NYV=20*NCM,NLAM=NCM-3) ! maximum number of constraints
      Dimension C20(*),EQP(*)
      Dimension IX(NCM)
      Dimension AMMN(NCM,NCM,0:3)
c
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC3/ NNC20(4,NCM),NNC11(4,NCM)
      Common /ITLOC4/ NNH20(4),NNH11(4),NNEQP(4)
c
      do i=1,ncm
         do j=i,ncm
	 
	    do L=0,3
	       AMMN(i,j,L) = 0.0d+00
	       AMMN(j,i,L) = 0.0d+00
	    end do
         
            if(IX(i)*IX(j).eq.1) then
      	       do L=0,3
            
                  FF = 0.0d+00
                  do it=1,4
 
                    IUVI = NNC20(it,i)
                    IUVJ = NNC20(it,j)
                    N    = ND (it)
                    N1   = N+NBLOCK(it)
                    N2   = N-NBLOCK(it)
                    IE1  = NNEQP(it)
                    IE2  = IE1 + N1
                    FF= FF+
     &		    RMMP(L,N1,N2,C20(IUVI),C20(IUVJ),EQP(IE1),EQP(IE2))
                  end do		      ! IT 
	          AMMN(i,j,L) = 2.0d+00*FF 
	          AMMN(j,i,L) = 2.0d+00*FF 
	       end do		      ! L
	       
            end if 		      ! IC(J)*IC(I)
         end do                          !    J
      end do                             !    I

c
      return 
      end
      
      Subroutine COLLMas0 (C20,EQP,AMMN,eta)
      Implicit real*8 (A-H,O-Z)
      Dimension C20(*),EQP(*)
      Dimension AMMN(0:3)
c
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC4/ NNH20(4),NNH11(4),NNEQP(4)
c
	 
      do L=0,3
         AMMN(L) = 0.0d+00
      end do
      
      do L=0,3
      
         FF = 0.0d+00
         do it=1,4
 
           IUV  = NNH20(it)
           N	= ND (it)
           N1	= N+NBLOCK(it)
           N2	= N-NBLOCK(it)
           IE1  = NNEQP(it)
           IE2  = IE1 + N1
           X=RMMPE(L,N1,N2,C20(IUV),C20(IUV),EQP(IE1),EQP(IE2),eta)
           FF= FF+X
c	   write(6,*) ' it l ',it,l,ff,EQP(IE1),EQP(IE2)
         end do 		     ! IT 
         AMMN(L) = 2.0d+00*FF 
      end do		     ! L	      
c
      return 
      end
c Set of subroutines calling the functions defined in TALMANF.f
      Subroutine TZ1DIM(nzmax,maxtz1,tz1d,itz1d)
      Implicit real*8 (A-H,O-Z)
c ---------------------------------------------------<< Start    include
      Include 'COMDIM'
c ---------------------------------------------------<< End      include
      Dimension tz1d(maxtz1)
      Dimension itz1d(nzmax,nzmax)
c
      Common /fact/dfact(Nfac),ddfact(Nfac),Nfacm
c
      if(nfac.ne.nfacm) then
        write(6,*) ' ERROR IN TZ1DIM:  NFAC,NFACM',NFAC,NFACM
        stop
      end if
c
      indx = 0
c
      do 1 in =1,nzmax
           do 2 im=in,nzmax
                imumin = im-in+1
                imumax = im+in-1
                itz1d(in,im) = indx
                itz1d(im,in) = indx
                do 3 imu =imumin,imumax,2
                     indx=indx+1
                     tz1d(indx)=DT1(in-1,im-1,imu-1)
    3           continue
    2      continue
    1 continue
      if((indx).gt.maxtz1) call errout('ITZ > MAXTZ1    ','TZ1DIM  ',4)
      return
      end
c+---------------------------------------------------------------------+
c|    Subroutine J1BB: Computes the J(n,m,mu) coeficients for the      |
c|    Brink-Boecker interaction.                                       |
c|                                                                     |
c|     - m >= n  and mu = 0,...,2*(Nzmax-1)                            |
c|     - mu must have the same parity of n+m                           |
c+---------------------------------------------------------------------+
c|    Parameters:                                                      |
c|                                                                     |
c|           - Nzmax ..... Maximum value of nz + 1                     |
c|           - Maxjz1 .... Number of elements in J                     |
c|           - Maxtz1 .... Number of elements in T                     |
c|           - Tz1d ...... Vector containing T                         |
c|           - Itz1d ..... Index of T                                  |
c|           - Work ...... Scratch vector                              |
c|           - Nwork ..... Dimension of Work                           |
c|                                                                     |
c|           - Amu(2) .... Interaction lenghts mu1 and  mu2            |
c|           - b ......... Oscillator Parameter                        |
c|                                                                     |
c|           - Ajz1d ..... Output vector containing J                  |
c|           - ijz1d ..... Index of J                                  |
c+---------------------------------------------------------------------+
c|                                             |                       |
c|    Dependencies:  Common /Fact/             |     SETUPN            |
c|                                             |                       |
c+---------------------------------------------------------------------+
      Subroutine J1BB (nzmax,maxjz1,maxtz1,tz1d,itz1d,nmax,work,
     * amu,b,ajz1d,ijz1d)
      Implicit Real*8 (A-H,O-Z)
c ---------------------------------------------------<< Start    include
      Include 'COMDIM'
c ---------------------------------------------------<< End      include
      Dimension ajz1d(2,maxjz1),tz1d(maxtz1)
      Dimension ijz1d(nzmax,nzmax),itz1d(nmax,nmax)
      Dimension amu(2)
c ----------------------------------------------------------------------
      Common /fact/dfact(Nfac),ddfact(Nfac),Nfacm
c ----------------------------------------------------------------------
      Maxmu = 2*nzmax - 1
      Minmu = 1
c ----------------------------------------------------------------------
      indx = 0
c
      do 1 in =1,nzmax
           do 2 im=in,nzmax
                ijz1d(in,im) = indx
                ijz1d(im,in) = indx
                inumin = im-in+1
                inumax = im+in-1
                Minmu  = Mod(in+im,2)+1
                do 3 imu = Minmu,Maxmu,2
                indx = indx + 1
                ajz1d(1,indx) = DJ1BB(in-1,im-1,imu-1,amu(1)/b)
                ajz1d(2,indx) = DJ1BB(in-1,im-1,imu-1,amu(2)/b)
    3           continue
    2      continue
    1 continue
      if((indx).gt.maxjz1) call errout('ITZ > MAXJZ1    ',' J1BB   ',4)
      return
      end
c+---------------------------------------------------------------------+
c³    Subroutine LSOMU: Computes the L(n,m,mu) coeficients for the     ³
c³    Spin-Orbin matrix elements.                                      ³
c³                                                                     ³
c³     - m =0,...,NZMAX                                                ³
c³     - m >= n  and mu = 0,...,2* Nzmax                               ³
c³     - mu must have the same parity of n+m                           ³
c+---------------------------------------------------------------------+
c³    Parameters:                                                      ³
c³                                                                     ³
c³           - Nzmax1..... Maximum value of nz + 2                     ³
c³           - Maxl1 ..... Number of elements in L                     ³
c³           - Maxtz2 .... Number of elements in T                     ³
c³           - Tz2d ...... Vector containing T                         ³
c³           - Itz2d ..... Index of T                                  ³
c³           - Work ...... Scratch vector                              ³
c³           - Nwork ..... Dimension of Work                           ³
c³                                                                     ³
c³                                                                     ³
c³           - Al1d ...... Output vector containing J                  ³
c³           - il1d ...... Index of J                                  ³
c+---------------------------------------------------------------------+
c³                                             ³                       ³
c³    Dependencies:  Common /Fact/             ³     SETUPN            ³
c³                                             ³                       ³
c+---------------------------------------------------------------------+
      Subroutine LSOMU
     *(nzmax1,maxl1,maxtz2,tz2d,itz2d,work,al1d,il1d)
c
      Implicit Real*8 (A-H,O-Z)
c ---------------------------------------------------<< Start    include
      Include 'COMDIM'
c ---------------------------------------------------<< End      include
      Dimension al1d(maxl1),tz2d(maxtz2)
      Dimension il1d(nzmax1,nzmax1),itz2d(nzmax1,nzmax1)
c ----------------------------------------------------------------------
      Common /fact/dfact(Nfac),ddfact(Nfac),Nfacm
c ----------------------------------------------------------------------
      Maxmu = 2*nzmax1 - 1
      Minmu = 1
c ----------------------------------------------------------------------
      indx = 0
c
      do 1 in =1,nzmax1
           do 2 im=in,nzmax1
                il1d(in,im) = indx + 1
                il1d(im,in) = indx + 1
                Minmu  = Mod(in+im,2)+1
                do 3 imu = Minmu,Maxmu,2
                indx = indx + 1
                al1d(indx) = DL1(in-1,im-1,imu-1)
    3           continue
    2      continue
    1 continue
      if((indx).gt.maxl1) call errout(' IL > MAXL1     ',' LSOMU  ',4)
      return
      end
c+---------------------------------------------------------------------+
c³    Subroutine DSOMU: Computes the D(n|m,mu) coeficients for the     ³
c³    Spin-Orbin matrix elements.                                      ³
c³                                                                     ³
c³     - m =0,...,NZMAX-1                                              ³
c³     - m >= n  and mu = 0,...,2* NZMAX                               ³
c³     - mu must have opposite parity to n+m                           ³
c+---------------------------------------------------------------------+
c³    Parameters:                                                      ³
c³                                                                     ³
c³           - Nzmax ..... Maximum value of nz + 1                     ³
c³           - Nzmax1..... Maximum value of nz + 2                     ³
c³           - Maxl1 ..... Number of elements in L                     ³
c³           - MaxD1 ..... Number of elements in D                     ³
c³           - Al1d ...... Vector containing L                         ³
c³           - Il1d ...... Index of L                                  ³
c³                                                                     ³
c³                                                                     ³
c³           - DSOmu...... Output vector containing D                  ³
c³           - idsomu .... Index of D                                  ³
c+---------------------------------------------------------------------+
c³                                             ³                       ³
c³    Dependencies:  Subroutine LSOMU          ³                       ³
c³                                             ³                       ³
c+---------------------------------------------------------------------+
c not modified
      Subroutine DSOMU
     *(nzmax,nzmax1,maxl1,maxd1,al1d,il1d,DSOmuv,iDSOmu)
c
      Implicit Real*8 (A-H,O-Z)
      Logical v1
      Dimension al1d(maxl1),DSOmuv(maxD1)
      Dimension il1d(nzmax1,nzmax1),idsomu(nzmax,nzmax)
c
      indx   = 0
      Maxmu = 2*NZMAX+1
c
      do 1 in =1,nzmax
           v1 = in.ne.1
           fac1 = dsqrt(2.0d+00*dfloat(in-1))
           fac2 = dsqrt(2.0d+00*dfloat(in))
           do 2 im=1,nzmax
                idsomu(in,im) = indx + 1
                if(v1) indt1  = il1d(im,in-1)
                indt2         = il1d(im,in+1)
                Minmu  = 2-Mod(in+im,2)
                do 3 imu = Minmu,Maxmu,2
                   sum = -fac2*al1d(indt2)
                   if(v1) sum = sum + fac1*al1d(indt1)
                   indt1 = indt1+1
                   indt2 = indt2+1
                   indx = indx + 1
                   DSOmuv(indx) = sum
    3           continue
    2      continue
    1 continue
      if((indx).gt.maxD1) call errout(' IL > MAXD1     ',' DSOMU  ',4)
      return
      end
c+---------------------------------------------------------------------+
c
c   February 20, 1992
c
c
c+---------------------------------------------------------------------+
c   Subroutine SXTJ               __
c                                \
c     Computes the contraction    >   T(k1,iq1,imu) J(iq2,k2,imu)
c                                /__
c                                imu
c
c     For k1 >= iq1 and k2 >= iq2
c
c     The vector IXTJ(k1,iq1) gives the position of the above contraction
c     for a given k1,iq1 in the XTJ vector.
c
c     The dimension of XTJ is
c                               2
c     ((nzmax2+1)(nzmax-nzmax2))  +
c                                                   2
c     ( nzmax(nzmax+1)/2 -(nzmax2+1)(nzmax-nzmax2) )
c
c     where nzmax2 = integer_part_of( nzmax/2 )
c
c   The T(k1,iq1,mu) is assumed to be computed with a maximum number
c   of k1,etc equal to nmax
c
c
c+---------------------------------------------------------------------+
      Subroutine Sxtj(NZMAX,nmax,MAXTZ1,TZ1D,ITZ1D,MAXJZ1D,AJZ1D,IJZ1D,
     x                MAXXTJ,XTJ)
      Implicit real*8 (A-H,O-Z)
      Dimension tz1d(maxtz1),ajz1d(2,maxjz1d),xtj(2,maxxtj)
      Dimension itz1d(nmax,nmax),ijz1d(nzmax,nzmax)
c
      ic = 0
      do 10 k1=1,nzmax
         do 11 iq1 = 1,k1
         im1 = Mod((k1+iq1),2)
         it1 = itz1d (k1,iq1)
         imumax = k1+iq1-1
         imumin = k1-iq1+1
         minmu  = im1 +1             ! the starting imu for J
         ish    = (imumin-minmu)/2   ! the shift in do 30
            do 20 k2=1,nzmax
               do 21 iq2=1,k2
               im2 = Mod((k2+iq2),2)
               if (im1.ne.im2) goto 21  ! This can be done in do 21
               ij2 = ijz1d(k2,iq2)+ish
               it2 = it1
               sum1 = 0.0d+00
               sum2 = 0.0d+00
                  do 30 imu = imumin,imumax,2
                     it2 = it2 + 1
                     ij2 = ij2 + 1
                     sum1= sum1+ tz1d(it2)*ajz1d(1,ij2)
                     sum2= sum2+ tz1d(it2)*ajz1d(2,ij2)
   30             continue
               ic  = ic + 1
               xtj(1,ic) = sum1
               xtj(2,ic) = sum2
   21          continue
   20       continue
   11    continue
   10 continue
      if(ic.gt.maxxtj) call errout('IC > MAXXTJ     ','SXTJ    ',4)
      return
      end
c
c   February 17, 1992
c
c
c   Subroutine SXTJP              __
c                       k1+iq1   \
c     Computes      (-1)          >   T(k1,iq1,imu) J(iq2,k2,imu)
c                                /__
c                                imu
c
c     For k1 >= iq1 and k2 >= iq2
c
c     This subroutine is needed for the exchange and pairing part
c
c     The vector IXTJ(k1,iq1) gives the position of the above contraction
c     for a given k1,iq1 in the XTJ vector.
c
c     The dimension of XTJ is
c                               2
c     ((nzmax2+1)(nzmax-nzmax2))  +
c                                                   2
c     ( nzmax(nzmax+1)/2 -(nzmax2+1)(nzmax-nzmax2) )
c
c     where nzmax2 = integer_part_of( nzmax/2 )
c
c
c   The T(k1,iq1,mu) is assumed to be computed with a maximum number
c   of k1,etc equal to nmax
c
c
      Subroutine Sxtjp(NZMAX,nmax,MAXTZ1,TZ1D,ITZ1D,MAXJZ1D,AJZ1D,IJZ1D,
     x                MAXXTJ,XTJ)
      Implicit real*8 (A-H,O-Z)
      Dimension tz1d(maxtz1),ajz1d(2,maxjz1d),xtj(2,maxxtj)
      Dimension itz1d(nmax,nmax),ijz1d(nzmax,nzmax)
c
      ic = 0
      do 10 k1=1,nzmax
         do 11 iq1 = 1,k1
         im1 = Mod((k1+iq1),2)
         it1 = itz1d (k1,iq1)
         phas= dfloat(1-2*im1)
         imumax = k1+iq1-1
         imumin = k1-iq1+1
         minmu  = im1 +1             ! the starting imu for J
         ish    = (imumin-minmu)/2   ! the shift in do 30
            do 20 k2=1,nzmax
               do 21 iq2=1,k2
               im2 = Mod((k2+iq2),2)
               if (im1.ne.im2) goto 21  ! This can be done in do 21
               ij2 = ijz1d(k2,iq2)+ish
               it2 = it1
               sum1= 0.0d+00
               sum2= 0.0d+00
                  do 30 imu = imumin,imumax,2
                     it2 = it2 + 1
                     ij2 = ij2 + 1
                     sum1= sum1+ tz1d(it2)*ajz1d(1,ij2)
                     sum2= sum2+ tz1d(it2)*ajz1d(2,ij2)
   30             continue
               ic  = ic + 1
               xtj(1,ic) = sum1* phas
               xtj(2,ic) = sum2* phas
   21          continue
   20       continue
   11    continue
   10 continue
      if(ic.gt.maxxtj) call errout('IC > MAXXTJ     ','SXTJ    ',4)
      return
      end
c
c   February 14, 1992
c
c   Subroutine SIXTJ:
c
c   Computes two integer matrices needed to locate the elements of XTJ
c   easily.
c
c     The matrix IXTJ(k1,iq1) gives the position in the vector XTJ
c     where the k2,iq2 values start for a given k1,iq1.
c
c     The matrix IXTJ2(k2,iq2) gives the relative position (with
c     respect to the starting point ixtj(k1,iq1) of the element in
c     XTJ
c                             BB
c     The position in XTJ of I                is then given by
c                             k1 iq2 iq1 k2
c
c     IXTJ(k1,iq1)+IXTJ2(k2,iq2)
c
c
      Subroutine Sixtj(NZMAX,IXTJ,IXTJ2)
      Implicit real*8 (A-H,O-Z)
      Dimension IXTJ(nzmax,nzmax),IXTJ2(nzmax,nzmax)
c
      ic = 0
      ice= 0
      ico= 0
      do 10 k1=1,nzmax
         do 11 iq1 = 1,k1
         im1 = Mod((k1+iq1),2)
         if(im1.eq.0) then
            ixtj2(k1,iq1) = ice
            ixtj2(iq1,k1) = ice
            ice = ice + 1
         else
            ixtj2(k1,iq1) = ico
            ixtj2(iq1,k1) = ico
            ico = ico + 1
         end if
         ixtj(k1,iq1) = ic + 1
         ixtj(iq1,k1) = ic + 1
            do 20 k2=1,nzmax
               do 21 iq2=1,k2
               im2 = Mod((k2+iq2),2)
               if (im1.ne.im2) goto 21  ! This can be done in do 21
               ic  = ic + 1
   21          continue
   20       continue
   11    continue
   10 continue
      return
      end
c Falta: 
c      1) Transformar a la forma cuadrada Gamma y Delta
c      2) Anadir la energia cinetica a uno y dos cuerpos
c      3) Definir un gamma y delta en forma cuadrada
c
c+---------------------------------------------------------------------+
c|  ROM and AKAP.... Density matrix and pairing tensor in square form  |
c|  GAMMA, DELTA.... HF field and pairing field in square form         |
c+---------------------------------------------------------------------+
      Subroutine HFBFIELD(ROM,AKAP,GAMMA,DELTA,GRO,GFLD,SS)
      Implicit real*8 (A-H,O-Z)
      Character*8 TEE
      Logical vcom,vcom2
      Include 'DIMPERM'
      Include 'DIMTRIAX'
c --------------------------------------------------- PERMANENT VECTORS
      Parameter (NEEE=30) ! dimension of the output matrix
c -------- The dimension parameters are in DIMPERM
      Dimension tz1d(iaxtz1d)
      Dimension ajx1d(2,iaxjx1d),ajy1d(2,iaxjy1d),ajz1d(2,iaxjz1d)
      Dimension itz1d(imax,imax)
      Dimension ijx1d(ixmax,ixmax),ijy1d(iymax,iymax),ijz1d(izmax,izmax)
      Dimension xtjx(2,iaxtjx),xtjxp(2,iaxtjx)
      Dimension xtjy(2,iaxtjy),xtjz(2,iaxtjz)
      Dimension ixtj(ixmax,ixmax),ixtj2(ixmax,ixmax)
      Dimension iytj(iymax,iymax),iytj2(iymax,iymax)
      Dimension iztj(izmax,izmax),iztj2(izmax,izmax)
c
      Dimension wfho(imax,iherm),xherm(iherm),weher(iherm)
      Dimension wfhonec(imax,iherm),xhnec(iherm)       ! necking
      Dimension wf(iwf2,iherm),xeherm(iherm)
c
      Dimension xleg(ilegn),wleg(ilegn)
c
      Dimension Tz2D(iaxtz2d),itz2d(imax1,imax1)
      Dimension AL1D(iaxlso),il1d(imax1,imax1)
      Dimension DSOMUv(iaxdso),IDSOMU(imax,imax)
c
      Dimension ACOU(iacou),COUM(imax,imax)
c ----------------------------------------------------------------------
      Save tz1d,ajx1d,ajy1d,ajz1d,xtjx,xtjxp,xtjy,xtjz
      Save wfho,wfhonec,xhnec,xherm,weher,wf,xeherm,xleg,wleg,tz2d,al1d
      Save dsomuv,acou,coum
      Save itz1d,ijx1d,ijy1d,ijz1d
      Save ixtj,ixtj2,iytj,iytj2,iztj,iztj2,itz2d,il1d,idsomu
c ----------------------------------------------------------------------
      Dimension ROM  (*), AKAP (*)      
      Dimension GAMMA(*), DELTA(*)
      
      Dimension GRO(*), GFLD(*)
      Dimension SS(*)              ! Scratch
C
c --- C O M M O N     B L O C K S -----------------------------------
c
      Common/DIMPDD/ irortmp,idmus0,jdmus0,idmus1,jdmus1,irorz,irory
c
      Common/DIMEN/nmax,nmax1,nxmax,nymax,nzmax,nwf2,
     * maxtz1,maxjz1,maxjx1,maxjy1,maxlso,maxdso,
     * ndmu,ndtheta,nacou,nrop,nrom,nrop8,nrom8,nherm,nherm38,nlegn
c
      common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
     
      Common /EEEEEE/ TEE(NEEE),EEE(7,NEEE)
      Common /IIIEEE/ mkin,mehf,mepair,mehfb,medd,merea,mecou,
     &                mn,mjx,msx,mq20,mq22,mq40,mneck,mr2,
     &                mbet2,mgamm,mbet4,maxrat,
     &                mfn,mjx2,mjy2,mjz2
c
      common /mfval/ ehfp,ehfn,epaip,epain,errp,errn,ecech,ojxp,ojxn,
     *        op,on,ojz2p,ojz2n,oq20p,oq20n,oq22p,oq22n,or2p,or2n,
     *        osxp,osxn
c
      Common /HFBOPT/Amass,vcom,vcom2,icouech
      
      Common /NOTC  / icall
      Common /CDEB  / ideb(9)
      Common /IFIELD/ idir,iech
      
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC1/ NNU(4),NNV(4),NNRO(4),NNKA(4)
c--------------------------------------------------------------------
      NP2 = ND2(1)
      NM2 = ND2(2)
c-----------------------------------------------------------------------
      if(icall.ne.61060) then    ! first time calculations
c+---------------------------------------------------------------------+
c     Computes dimensions
c+---------------------------------------------------------------------+
        call dimens(ixmax,iymax,izmax,iherm,ilegn,imy,imz)
c+---------------------------------------------------------------------+
c     Initializes permanent vectors
c+---------------------------------------------------------------------+
        call INIT(tz1d,ajx1d,ajy1d,ajz1d,xtjx,xtjxp,xtjy,xtjz,wfhonec,
     * xhnec,wfho,xherm,xeherm,weher,wf,xleg,wleg,tz2d,al1d,dsomuv,acou,
     * coum,SS,itz1d,itz2d,ijx1d,ijy1d,ijz1d,il1d,idsomu,ixtj,ixtj2,
     * iytj,iytj2,iztj,iztj2)

         icall = 61060
      end if
       
       IG = nrop8+nrom8
       
       do ii=1,ig
         gro(ii) = 0.0d+00
         gfld(ii) = 0.0d+00
       end do
c
c
      IG1P= 1
      IG2P= 1 +   NROP
      IG3P= 1 + 2*NROP
      IG4P= 1 + 3*NROP
      IG5P= 1 + 4*NROP
      IG6P= 1 + 5*NROP
      IG7P= 1 + 6*NROP
      IG8P= 1 + 7*NROP
      
      IGP0 = NROP8
      
      IG1M= 1            + IGP0
      IG2M= 1 +   NROM   + IGP0 
      IG3M= 1 + 2*NROM   + IGP0 
      IG4M= 1 + 3*NROM   + IGP0 
      IG5M= 1 + 4*NROM   + IGP0 
      IG6M= 1 + 5*NROM   + IGP0 
      IG7M= 1 + 6*NROM   + IGP0
      IG8M= 1 + 7*NROM   + IGP0
c
c+---------------------------------------------------------------------+
c Transforms the density matrix and the pairing tensor to packed format
c+---------------------------------------------------------------------+
       call ro2pack(rom,akap,gro(ig1p),gro(ig1m))
c
c+---------------------------------------------------------------------+
c|  Computes several density matrices in coordinate representation     |
c|  for the calculation of the density dependent part                  |
c+---------------------------------------------------------------------+
       call timeit(0,18,'DENR            ')
       it1 = 1   +15*NDMU           ! ror
       it2 = it1 +   Nherm38        ! rorz    These are redefined
       it3 = it2 +   irorz          ! rory    afterwards.
       
       call denr (gro(ig1p),gro(ig1m),SS(it1),wfho,wfhonec,
     *            SS(it2),SS(it3))
     
       call timeit(1,18,'DENR            ')
c+---------------------------------------------------------------------+
c   Calculations previous to the direct field calculation for the
c   density dependent part.
c   Densitity dependent terms ( Mean Field, Rearrangement and Coulomb
c   exchange in the Slater approximation )
c+---------------------------------------------------------------------+
      it1 =   1 +15*NDMU           ! ror
      it2 = it1 +   Nherm38        ! rortmp
      it3 = it2 +   irortmp        ! DMUS0
      it4 = it3 +   jdmus0         ! DMUS1
      it5 = 1                      ! DMU
c
      if(idir.ne.0) then
      call timeit(0,13,'  PREDIRDD      ')
      call predirdd(wf,xherm,xeherm,SS(it1),SS(it2),
     *                              SS(it3),SS(it4),SS(it5))
      call timeit(1,13,'  PREDIRDD      ')
      end if
c
      it5 = 1                      ! DMU
      it6 = it5 +15*NDMU           ! sxyz1p
      it7 = it6 + 2*NDMU           ! sxyz1n
      it8 = it7 + 2*NDMU           ! sxyz2p
      it9 = it8 + 2*NDMU           ! sxyz2n
      it10= it9 + 2*NDMU           ! sxyz3p
      it11= it10+ 2*NDMU           ! sxyz3n
      it12= it11+ 2*NDMU           ! sxyz4p
      it13= it12+ 2*NDMU           ! sxyz4n
      it14= it13+ 2*NDMU           ! sxyzc
      it15= it14+   NDMU           ! scf
      it16= it15+   NDMU           ! soxyzp
      it17= it16+10*NDMU           ! soxyzn
      it18= it17+10*NDMU           ! GSUM0EP
      it19= it18+   NDMU           ! GSUM0EN
      it20= it19+   NDMU           ! GSUM0OP
      it21= it20+   NDMU           ! GSUM0ON
      it22= it21+   NDMU           ! GSUM1OP
      it23= it22+   NDMU           ! GSUM1ON
      it24= it23+   NDMU           ! GSUM3OP
      it25= it24+   NDMU           ! GSUM3ON
      it26= it25+   NDMU           ! THE0EP
      it27= it26+   NDTHETa        ! THE0EN
      it28= it27+   NDTHETa        ! THE0OP
      it29= it28+   NDTHETa        ! THE0ON
      it30= it29+   NDTHETa        ! THE1OP
      it31= it30+   NDTHETa        ! THE1ON
      it32= it31+   NDTHETa        ! THE3OP
      it33= it32+   NDTHETa        ! THE3ON
      it34= it33+   NDTHETa        ! opt spin orbit 1
      it35= it34+   NDTHETa        !                2
      it36= it35+   NDTHETa        !                3
      it37= it36+   NDTHETa        !                4
      it38= it37+   NDTHETa        !                5
      it39= it38+   NDTHETa        !                6
      it40= it39+   NDTHETa        !                7
      it41= it40+   NDTHETa        !                8
      it42= it41+   NDTHETa        !                9
      it43= it42+   NDTHETa        !               10
      it44= it43+   NDTHETa        !               11
      it45= it44+   NDTHETa        !               12
      it46= it45+   NDTHETa        !
      IDIGDIR = IT46-1
c      write(6,*) ' DIMENSION OF TEMP IN GDIR  ',IDIGDIR
      if(idir.ne.0) then
       call timeit(0,7,' DIRECT FIELD   ')
c+---------------------------------------------------------------------+
c|                                             D I R E C T     T E R M |
c+---------------------------------------------------------------------+
      Call GDIR(GRO(IG1P),GRO(IG1M),tz1d,ajx1d,ajy1d,ajz1d,al1d,dsomuv,
     * acou,coum,itz1d,ijx1d,ijy1d,ijz1d,il1d,idsomu,
     *SS(it5 ),SS(it6 ),SS(it7 ),SS(it8 ),SS(it9 ),SS(it10),
     *SS(it11),SS(it12),SS(it13),SS(it14),SS(it15),SS(it16),
     *SS(it17),SS(it18),SS(it19),SS(it20),SS(it21),SS(it22),
     *SS(it23),SS(it24),SS(it25),SS(it26),SS(it27),SS(it28),
     *SS(it29),SS(it30),SS(it31),SS(it32),SS(it33),SS(it34),
     *SS(it35),SS(it36),SS(it37),SS(it38),SS(it39),SS(it40),
     *SS(it41),SS(it42),SS(it43),SS(it44),SS(it45),
     * GFLD(IG1P),GFLD(IG2P),GFLD(IG3P),GFLD(IG4P),
     * GFLD(IG1M),GFLD(IG2M),GFLD(IG3M),GFLD(IG4M))
       call timeit(1,7,' DIRECT FIELD   ')
      end if
c
      if(iech.ne.0) then
c+---------------------------------------------------------------------+
c|                 E X C H A N G E   A N D   P A I R I N G             |
c+---------------------------------------------------------------------+
      call timeit(0,8,'EXCHANGE+PAIRING')
c+---------------------------------------------------------------------+
c|               S P I N - O R B I T     P A I R I N G                 |
c+---------------------------------------------------------------------+
      if(ideb(8).eq.0) then
      isym = 0
      call Paiso(isym,gro(Ig1P),gro(Ig1M),GFLD(IG5P),GFLD(IG7P),
     * GFLD(IG5M),GFLD(IG7M),SS(it12),SS(it5))
c     
      isym = 1
      call Paiso(isym,gro(Ig1P),gro(Ig1M),GFLD(IG6P),GFLD(IG8P),
     * GFLD(IG6M),GFLD(IG8M),SS(it12),SS(it5))

      do i=1,nrop
         d1 = gfld(IG5p+i-1)+gfld(IG6p+i-1)
         d2 = gfld(IG5p+i-1)-gfld(IG6p+i-1)
         gfld(IG5p+i-1)=d1
         gfld(IG6p+i-1)=d2
         d1 = gfld(IG7p+i-1)+gfld(IG8p+i-1)
         d2 = gfld(IG7p+i-1)-gfld(IG8p+i-1)
         gfld(IG7p+i-1)=d1
         gfld(IG8p+i-1)=d2
      end do
      do i=1,nrom
         d1 = gfld(IG5m+i-1)+gfld(IG6m+i-1)
         d2 = gfld(IG5m+i-1)-gfld(IG6m+i-1)
         gfld(IG5m+i-1)=d1
         gfld(IG6m+i-1)=d2
         d1 = gfld(IG7m+i-1)+gfld(IG8m+i-1)
         d2 = gfld(IG7m+i-1)-gfld(IG8m+i-1)
         gfld(IG7m+i-1)=d1
         gfld(IG8m+i-1)=d2
      end do
      end if
c+---------------------------------------------------------------------+
c|       B R I N K - B O E C K E R       P A I R I N G                 |
c+---------------------------------------------------------------------+
      jt1 = 1
      jt2 = jt1 + mazopti
      icero = 1
      icerop= 0
      call GECH(icero,icerop,GRO(Ig1P),GRO(IG1M),xtjx,xtjxp,xtjy,xtjz,
     * ixtj,ixtj2,iytj,iytj2,iztj,iztj2,    
     * SS(jt1),SS(jt2),
     * GFLD(IG1P),GFLD(IG2P),GFLD(IG3P),GFLD(IG4P),
     * GFLD(IG5P),GFLD(IG6P),GFLD(IG7P),GFLD(IG8P),
     * GFLD(IG1M),GFLD(IG2M),GFLD(IG3M),GFLD(IG4M),
     * GFLD(IG5M),GFLD(IG6M),GFLD(IG7M),GFLD(IG8M))
       call timeit(1,8,'EXCHANGE+PAIRING')
      end if
      
c
      call field2sq(gfld,gamma,delta)
c
c+---------------------------------------------------------------------+
c     Computes Kinetic Energy
c+---------------------------------------------------------------------+
c
       call kin(SS(1),SS(1+NP2)) ! TE TO
       com = 1.0d+00
       if(vcom) com = 1.0d+00-1.0d+00/Amass        ! Center of mass corr.
c+---------------------------------------------------------------------+
c   Computes the two body correction to the kinetic energy
c+---------------------------------------------------------------------+
       if(vcom2) then
          ikin2G = 1+NP2+NM2
	  ikin2D = ikin2G + 4*(NP2+NM2)
	  call KIN2BC(ROM,AKAP,SS(IKIN2G),SS(IKIN2D))
       end if
c       
c       Adding up the kinet energy and computing the energy	  
c
c
       ehfp = 0.0d+00
       ehfn = 0.0d+00
       
       do it=1,4
       
	 IRO1 = NNRO(it)
	 IRO2 = NNRO(it) + ND2(it)
         IKIN = (1-Mod(it,2))*NP2+1
	 
	 if(vcom2) then 
	   one = 1.0d+00
	   call daxpy(nd2(it),one,SS(IKIN2G+IRO1-1),1,gamma(IRO1),1)
	   call daxpy(nd2(it),one,SS(IKIN2G+IRO2-1),1,gamma(IRO2),1)
	   ykin1 = ddot(ND2(it),ROM(IRO1),1,SS(IKIN2G+IRO1-1),1)
	   ykin2 = ddot(ND2(it),ROM(IRO2),1,SS(IKIN2G+IRO2-1),1)
	   ykin  = 0.5d+00*(ykin1 + ykin2)
           EEE(it,24) =  ykin
	 end if
	 
	 EGR1 = ddot(ND2(it),Gamma(IRO1),1,ROM(IRO1),1)
	 EGR2 = ddot(ND2(it),Gamma(IRO2),1,ROM(IRO2),1)
	 
	 call daxpy(nd2(it),com,SS(IKIN),1,gamma(IRO1),1)
	 call daxpy(nd2(it),com,SS(IKIN),1,gamma(IRO2),1)
	 	 
	 xkin1 = ddot(ND2(it),ROM(IRO1),1,SS(IKIN),1)
	 xkin2 = ddot(ND2(it),ROM(IRO2),1,SS(IKIN),1)
	 xkin  = (xkin1 + xkin2)*com
	 
	 EGR = EGR1+EGR2
	 
c	 xx  = dfloat(it/3)*errn + dfloat(1-it/3)*(errp+ecech)
	 xx  = dfloat(it/3)*errn + dfloat(1-it/3)*(errp)
         if(icouech.eq.1) xx = xx + dfloat(1-it/3)*ecech
	 ehf = xkin + 0.5d+00*egr - 0.5d+00*xx
	 
	 if(it/3.eq.0) then
	    ehfp = ehfp + ehf
	 else 
	    ehfn = ehfn + ehf
	 end if
	 
         EEE(it,mkin) =  xkin
	 EEE(it,mehf) =  ehf

       end do
       
      Epair = 0.0d+00
      Epairp = 0.0d+00
      Epairn = 0.0d+00
      do it=1,4
      
         ID = NNKA(it)
	 
	 if(vcom2) then
	   one = 1.0d+00
	   call daxpy(nd2(it),one,SS(IKIN2D+ID-1),1,Delta(ID),1)
	   ykin = ddot(ND2(it),AKAP(ID),1,SS(IKIN2D+ID-1),1)
           EEE(it,25) =  ykin
	 end if
	 
	 EPP = ddot(ND2(it),Delta(ID),1,akap(ID),1)
         Epair = Epair + EPP
	 if(it/3.eq.0) then
	   Epairp = Epairp + EPP
	 else 
	   Epairn = Epairn + EPP
	 end if
         EEE(it,mepair) = epp
      end do

      do it=1,4
         EEE(it,mehfb)=EEE(it,mehf) + EEE(it,mepair)
      end do
      
c      etot = ehfp+ehfn+epairp+epairn
      
c      write(6,'( " ENERGIES ",5f15.5)') ehfp,ehfn,epairp,epairn,etot
      
      return
      end
c+------------------------------------------------------------------+
c|                                                                  |
c|  Triaxial HFB with 2nd order gradient (approximate)              |
c|                                                                  |
c+------------------------------------------------------------------+
      Program TEST
      Implicit real*8 (A-H,O-Z)
      Logical vcom,vcom2
      Character*8 TEE
      Include 'DIMTRIAX'
      Parameter (NEEE=30) ! dimension of the output matrix
      Parameter (NMAXB=4) ! maximum number of blockings in each channel
      Parameter (NCM=6,NCOP=(NCM-2)*(NP*NP+NM*NM),NLAM=NCM-3)
      Parameter (NAUX=3*(NP*NP+NM*NM))
      Parameter (NGAM=4*(NP*NP+NM*NM)) 
      Parameter (NC20=2*NCM*(NP*NP+NM*NM))
      Parameter (NG20=2*(NP*NP+NM*NM)) ! carefull with Nh20 in dimtriax
      Parameter (NG11=4*(NP*NP+NM*NM+2*NMAXB)) ! Now handles blocking properly
      Parameter (NEQP=4*(NP+NM))
      
      Dimension AMMN(NCM,NCM,0:3)
      
      Dimension COP(NCOP),rval(NCM)
      Dimension C20(NC20),AUX(NAUX)
      Dimension H20(NG20),Z0(NG20)
      Dimension H11(NG11)
      Dimension EQP(NEQP),QUOT(NG20)
      
      Dimension CANNON(NEQP,5) ! 1 = SPE 2=v**2 3=Delta 4=Q20  5=Q22
      
      Dimension ix(NCM),ICC(NCM),CLAM(NCM),CDLAM(NCM),DD(NCM),etc(NCM)
      Dimension KBLO  (4)
      Dimension NNUOLD(4)
      Dimension UVOLD (NUV)
      Dimension UV    (NUV)
      Dimension ROM   (NGAM), AKA   (NGAM)
      Dimension GAMMA (NGAM), DELTA (NGAM)
      Parameter (NSS=NH20+2*(IGP+IGM))
      Dimension SS    (NSS)
      
      Common /EEEEEE/ TEE(NEEE),EEE(7,NEEE)
      Common /IIIEEE/ mkin,mehf,mepair,mehfb,medd,merea,mecou,
     &                mn,mjx,msx,mq20,mq22,mq40,mneck,mr2,
     &                mbet2,mgamm,mbet4,maxrat,
     &                mfn,mjx2,mjy2,mjz2
     
      common /CITER/ epsg,etamax,etamin,dmax,dmin,tshh,maxiter
c
C   UV => U1P U2P V1P V2P U1N U2N V1N V2N      PARIDAD POSITIVA  8NP*NP
C         U1P U2P V1P V2P U1N U2N V1N V2N      PARIDAD NEGATIVA  8NM*NM

C          UV => U1PP U2PP U1PM U2PM  V1PP V2PP V1PM V2PM     
C                U1NP U2NP U1NM U2NM  V1NP V2NP V1NM V2NM     
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC1/ NNU(4),NNV(4),NNRO(4),NNKA(4)
      Common /ITLOC2/ ND2PACK(4),NNGPACK(4),NNDPACK(4)
      Common /ITLOC3/ NNC20(4,NCM),NNC11(4,NCM)
      Common /ITLOC4/ NNH20(4),NNH11(4),NNEQP(4)
      
      Common /CCONS/ cval(ncm),prec(ncm),IC(NCM),
     & ILAMB(NLAM),IMM(NLAM),ISOS(NCM,2)
     
      Common /HFBOPT/ Amass,vcom,vcom2,icouech
      Common /STARTO/ iwf,iol
      Common /IFIELD/ idir,iech
      Common /OSCLEN/ bx,by,bz
      Common /CDEB  / ideb(9)

      itstg = ideb(9)

      call readdata

      idir = 1
      iech = 1
             
      KBLO(1) = 0
      KBLO(2) = 0
      KBLO(3) = 0
      KBLO(4) = 0

      KP = NP
      KM = NM       

      call ITLOC(KP,KM,KBLO)
       
      NNUOLD(1) = 1
      NNUOLD(2) = 1 + 8*NP*NP
      NNUOLD(3) = 1 + 4*NP*NP
      NNUOLD(4) = 1 + 8*NP*NP+4*NM*NM
      
      NL=10
      if (iwf.eq.0) then 
       read (NL) bx0,by0,bz0,UVOLD
 
       do it=1,4
        NC = 2*ND2(it)
        call dcopy(NC,UVOLD(NNUOLD(it)),1,UV(NNU(it)),1)
        call dcopy(NC,UVOLD(NNUOLD(it)+NC),1,UV(NNV(it)),1)
       end do
        write(6,'( //,"  READING WF WITH OLD FORMAT ",//,40("="))')
      else if(iwf.eq.1) then 
      
          read (NL) bx0,by0,bz0,UV
          write(6,'( //,"  READING WF WITH NEW FORMAT ",//,40("="))')
      else
      
         call HOINIT(UV,cval(1),cval(2))
      
      end if
     
      if(iol.eq.1) then
       bx = bx0
       by = by0 
       bz = bz0
       write(6,'( //,"  READING BX, BY, BZ from disk",3x,3f7.4,//)')
     &                                                   bx,by,bz
      else if(iol.eq.0) then
       write(6,'( //,"  READING BX, BY, BZ from data",3x,3f7.4,//)')
     &                                                   bx,by,bz
      end if
c
c
c +---------------------------------------------------------------------+
c |            I T E R A T I O N    P R O C E D U R E                   |
c +---------------------------------------------------------------------+

      call consop(COP)
      
c
c +---------------------------------------------------------------------+
c |    A D J U S T I N G    T H E    C O N S T R A I N T S              |
c +---------------------------------------------------------------------+
      Call DOCONS (COP,ROM,C20,UV,IX,dd,Z0,AUX)
      
      kiter = 0
    1 continue       ! do while (shh.gt.epsg).and.(kiter.lt.maxiter)
      kiter = kiter + 1
      
c +---------------------------------------------------------------------+
c |    C O M P U T I N G    T H E    D E N S I T Y    M A T R I X       |
c +---------------------------------------------------------------------+

      do it=1,4
         N  = ND (it)
         call RO(UV(NNV(it)),N,N,N,ROM(NNRO(it)))
         call KAPPA(UV(NNU(it)),UV(NNV(it)),N,N,N,AKA(NNKA(it)))
      end do

c
c +---------------------------------------------------------------------+
c |            H F B              F I E L D S                           |
c +---------------------------------------------------------------------+
      IG = IGP+IGM
      
      call HFBFIELD(ROM,AKA,GAMMA,DELTA,SS,SS(IG+1),SS(2*IG+1))
c
c       In the first iteration we do not have lambdas
c
      if(kiter.eq.1) then
        
      SHH0 = 0.0d+00
      do it=1,4
         N  = ND (it)
         NU = NNU(it)
         NV = NNV(it)
         NG = NNRO(it)
         NK = NNKA(it)
         N20= NNH20(it)
         call SH20(UV(NU),UV(NV),N,N,N,GAMMA(NG),DELTA(NK),AUX,H20(N20))
         SHH0 = SHH0 + ddot(ND2(it),H20(N20),1,H20(N20),1)
      end do
      
      call CONS20 (UV,COP,C20,AUX)
c
      do i=1,ncm
         dd(i) = 0.0d+00
      end do
c
      call LAMBDA0 (C20,H20,dd,ix,clam,cdlam,shhc)
      
         ioutlev = 0
         if(ioutlev.ge.4) then
           write(6,'(" DD     ",6d13.6)') (dd   (i),i=1,ncm)
           write(6,'(" CLAM   ",6d13.6)') (clam (i),i=1,ncm)
           write(6,'(" CDLAM  ",6d13.6)') (cdlam(i),i=1,ncm)
         end if
      shh = shh0 + shhc
C
      end if ! first iteration       
c +---------------------------------------------------------------------+
c |                                                             11      |
c |                                                           H'        |
c +---------------------------------------------------------------------+
      SHH0 = 0.0d+00
      do it=1,4
      
         N  = ND (it)
	 N1 = N-NBLOCK(it)
	 N2 = N+NBLOCK(it)
         NN_1 = N1*N1
         NN_2 = N2*N2
	 NU = NNU(it)
	 NV = NNV(it)
	 NG = NNRO(it)
	 NK = NNKA(it)
	 N11_1 = NNH11(it)
	 N11_2 = N11_1+NN_1
	 N20= NNH20(it)
	 IEQP1 = NNEQP(it)
	 IEQP2 = IEQP1 + N1
	 
c +=====================================================================+
c ||      h'     = h     - clam(i)* O   (i)                            ||
c +=====================================================================+
c
         icopy = nn_1+nn_2
	 call dcopy(icopy,gamma(NG),1,SS,1)
	 call hprime(it,SS,COP,IX,CLAM)
c +---------------------------------------------------------------------+
c |                                                             11      |
c |                                                           H'        |
c +---------------------------------------------------------------------+
	 	 
         call SH11(UV(NU),UV(NV),N,N,N,SS,DELTA(NK),AUX,
     &	 H11(N11_1),H11(N11_2))
c     
c +=====================================================================+
c ||                                                       ' 11        ||
c ||      D I A G O N A L I Z A T I O N     O F           H            ||
c +=====================================================================+
         call dsyev('V','U',N1,H11(N11_1),N1,EQP(IEQP1),AUX,NAUX,INFO)
	 
	 if(INFO.ne.0) then
	   write(6,'(" PROBLEMS WITH DSYEV, IT, INFO ",2I3)') it,info
	 end if
	 	 
         call dsyev('V','U',N2,H11(N11_2),N2,EQP(IEQP2),AUX,NAUX,INFO)
	 
	 if(INFO.ne.0) then
	   write(6,'(" PROBLEMS WITH DSYEV, IT, INFO ",2I3)') it,info
	 end if
	 
         call QUOT2QPE(N1,N2,EQP(IEQP1),EQP(IEQP2),QUOT(N20))
c +=====================================================================+
c ||  T R A N S F O R M I N G     U  V    T O   Q. P.   B A S I S      ||
c +=====================================================================+
         call UVTRD(UV(NU),UV(NV),N1,N2,N,H11(N11_1),H11(N11_2),AUX)
c +---------------------------------------------------------------------+
c |                                                            20       |
c |                                                           H         |
c +---------------------------------------------------------------------+
         
       call SH20(UV(NU),UV(NV),N1,N2,N,GAMMA(NG),DELTA(NK),AUX,H20(N20))
c	 SHH0 = SHH0 + ddot(ND2(it),H20(N20),1,H20(N20),1)
	 SHH0 = SHH0 + rmm(N1*N2,H20(N20),H20(N20),QUOT(N20))
	 	 
      end do	 
      
c +---------------------------------------------------------------------+
c |                      C O N S T R A I N T S                          |
c +---------------------------------------------------------------------+
c
      
      call CONSMV (ROM,COP,rval)
      call CONS20 (UV,COP,C20,AUX)
c
      do i=1,ncm
         dd(i) = cval(i)-rval(i)
      end do
c
      call LAMBDA1 (C20,H20,QUOT,dd,ix,clam,cdlam,shhc)
      
         ioutlev = 0
         if(ioutlev.ge.4) then
           write(6,'(" DD     ",6d13.6)') (dd   (i),i=1,ncm)
           write(6,'(" CLAM   ",6d13.6)') (clam (i),i=1,ncm)
           write(6,'(" CDLAM  ",6d13.6)') (cdlam(i),i=1,ncm)
         end if
C --------------------------------------------------------- GRADIENT NORM
      shh = shh0 + shhc
C
c +=====================================================================+
c ||   D E T E R M I N I N G    T H E    G R A D I E N T    S T E P    ||
c +=====================================================================+

      ehfb = 0.0d+00
      do it=1,4
         ehfb = ehfb + EEE(it,mehfb)
      end do
      
      if(kiter.eq.1) then
         eta = 0.010d+00
         energy0 = ehfb
         de_pred  = shh0 + shhc
c >>>         eta = 0.0d+00
         if(itstg.eq.1) then
            eta = 1.d-08
            do i=1,ncm
               ix(i) = 0
            end do
            de_pred  = shh0
         end if
         de_pred = -2.0d+00*eta*de_pred
      else 
         de_real = ehfb-energy0
         energy0 = ehfb
         if(itstg.eq.1) then
            write(6,999) de_pred,de_real
            stop
         end if
c
         dtest   = de_real/de_pred
	 ett     = eta
         eta     = etanew(dtest,shh,shhold,ett)
         shhold  = shh
         de_pred  = shh0 + shhc
         de_pred = -2.0d+00*eta*de_pred
      end if
c
      do i=1,ncm
         if(ix(i).ne.0) etc(i) = eta*clam(i)+cdlam(i)
      end do
      
         do it=1,4
c
           N  = ND (it)
	   N1 = N - NBLOCK(it)
	   N2 = N + NBLOCK(it)
	   IU = NNU(it)
	   IV = NNV(it)
	   IZ = NNH20(it)
c +---------------------------------------------------------------------+
c |            T H O U L E S S     Z     M A T R I X                    |
c +---------------------------------------------------------------------+
           call ZETA(it,eta,H20(IZ),C20,Z0(IZ),IX,etc)
c ------------------------------------------------- DIV BY 2 QP ENERGIES
	   
           call DBE(N1,N2,Z0(IZ),QUOT(IZ))
c
c +---------------------------------------------------------------------+
c |            N E W    W A V E    F U N C T I O N S                    |
c +---------------------------------------------------------------------+
c
           call NEWUV(UV(IU),UV(IV),N1,N2,N,Z0(IZ),AUX)
	  
	 end do ! it
      
      if(Mod(kiter,10).eq.0) then 
         Call DOCONS (COP,ROM,C20,UV,IX,dd,Z0,AUX)
	
         write (11) bx,by,bz,UV
	 rewind 11
       end if 
	 
c +--------------------------------------------------------------------+
c |                  E N D       I  T  E  R  A  T  I  O  N             |
c +--------------------------------------------------------------------+
      write(6,'(i5,3D11.4,f7.4,f17.10)')kiter,shh,dtest,de_real,eta,ehfb
      
      
      if((shh.gt.epsg).and.(kiter.lt.maxiter)) goto 1  ! end do while
      
      if (shh.le.epsg) then
          write(6,103) shh,epsg,kiter,ehfb
      else
          write(6,104) shh,epsg,kiter,ehfb
      end if
      
      call COLLMas(C20,EQP,ix,AMMN)

      do L=0,3
         write(6,'( " M (-",I1,")  " )') L
         do i=1,NCM
	    write(6,'( 6d12.5 )') (AMMN(i,j,l),j=1,NCM)
	 end do
      end do

      call RESMV(ROM,AUX)
      call RESFLUC(UV,COP,AUX)

      write(6,900)
      
      do i=1,NEEE
          if(i.ne.mr2) then
             EEE(5,i) = EEE(1,i)+EEE(2,i)
             EEE(6,i) = EEE(3,i)+EEE(4,i)
             EEE(7,i) = EEE(5,i)+EEE(6,i)
	  end if
      end do
c +--------------------------------------------------------------------+
c |                                            Beta_2 and Gamma        |
c +--------------------------------------------------------------------+
c      
      pi = 4.0d+00*datan(1.0d+00)
      f2=dsqrt((4.d+00*pi)/5.0d+00)
      
      do it=1,7
         EEE(it,mbet2) = f2*dsqrt(EEE(it,mq20)**2+EEE(it,mq22)**2)/
     &	                 (EEE(it,mn)*EEE(it,mr2)**2)
         EEE(it,mgamm) = datan(EEE(it,mq22)/(EEE(it,mq20)+1.d-10))
     &            	 *180.0d+00/pi
         EEE(it,mq20 ) = EEE(it,mq20)*1.d-02 ! (barns)
         EEE(it,mq22 ) = EEE(it,mq22)*1.d-02 ! (barns)
         EEE(it,mq40 ) = EEE(it,mq40)*1.d-04 ! (barns**2)
      end do
      
      do i=1,NEEE
          write(6,901) TEE(i),(EEE(k,i),k=1,7)
      end do
c
      write(38) TEE,EEE,rval,clam,ix,ammn
c      
      write(6,902) (rval(i),i=1,ncm),(clam(i),i=1,ncm)

  103 format ( 10x,12('*******'),/,32x,
     &' P R O P E R L Y     F I N I S H E D  ',//,20x,' SHH ',f10.6,
     &' EPSG ',f10.6,' ITER ',i5,' E HFB ',f20.10,/,10x,12('*******'),/)

  104 format ( 10x,12('*******'),/,26x,
     &'M A X I M U M   O F  I T E R A T I O N S   E X C E E D E D',//
     &,20x,' SHH ',f10.6,
     &' EPSG ',f10.6,' ITER ',i5,' E HFB ',f20.10,/,10x,12('*******'),/)

      
  900 format ( 69x,'PROTON',6x,'NEUTRON',7x,'TOTAL',/)

  901 format (3x,A8,7(1x,f12.6))

  902 format ( //,3x,'CONSTR',6(f15.6,2x),//,3x,'LAMBDA',6(f15.10,2x))
  
  999 format ( 4x,5('******'),/,10x,'  G R A D I E N T    T E S T ',/,
     & 10x,' PREDICTED : ',f17.11,/,10x,' REAL      : ',f17.11)


c




c +---------------------------------------------------------------------+
c |                                                            11       |
c |                                                           H         |
c +---------------------------------------------------------------------+
c
      do it=1,4
      
         N  = ND (it)
         NN = N*N
         N1 = N-NBLOCK(it)
         N2 = N+NBLOCK(it)
         NN_1 = N1*N1
         NU = NNU(it)
         NV = NNV(it)
         NG = NNRO(it)
         NK = NNKA(it)
         N11_1 = NNH11(it)
         N11_2 = N11_1+NN_1
         IEQP1 = NNEQP(it)
         IEQP2 = IEQP1 + N1
         IE1   = NNEQP(it)
         IE2   = IE1 + N
	 
         write(6,*) ' Quasiparticle energies ',it
	 
         call hprime(it,gamma(NG),COP,IX,CLAM)
	 	 
         call SH11(UV(NU),UV(NV),N,N,N,GAMMA(NG),DELTA(NK),AUX,
     &	 H11(N11_1),H11(N11_2))
c     
         call dsyev('V','U',N1,H11(N11_1),N1,EQP(IEQP1),AUX,NAUX,INFO)
	 
         if(INFO.ne.0) then
           write(6,'(" PROBLEMS WITH DSYEV, IT, INFO ",2I3)') it,info
         end if
	 
         write(6,'("+i",10f12.5 )') (EQP(IEQP1+kk),kk=0,9)
 
         call dsyev('V','U',N2,H11(N11_2),N2,EQP(IEQP2),AUX,NAUX,INFO)
	 
         if(INFO.ne.0) then
            write(6,'(" PROBLEMS WITH DSYEV, IT, INFO ",2I3)') it,info
         end if
	 
         write(6,'( "-i",10f12.5 )') (EQP(IEQP2+kk),kk=0,9)
	 
	 	 
c+---------------------------------------------------------------------+
c|                                            Spe energies             |
c+---------------------------------------------------------------------+
c
         call dsyev('V','U',N,GAMMA(NG),N,CANNON(IE1,1),AUX,NAUX,INFO)

         a = 1.0d+00
         b = 0.0d+00
	 
         call dgemm('n','n',N,N,N,a,ROM(NG),N,GAMMA(NG),N,b,AUX,N)
         call dgemm('t','n',N,N,N,a,GAMMA(NG),N,AUX,N,b,AUX(NN+1),N)
         call cpydiag(AUX(NN+1),N,CANNON(IE1,2))
c         
         call dgemm('n','n',N,N,N,a,DELTA(NK),N,GAMMA(NG),N,b,AUX,N)
         call dgemm('t','n',N,N,N,a,GAMMA(NG),N,AUX,N,b,AUX(NN+1),N)
         call cpydiag(AUX(NN+1),N,CANNON(IE1,3))
         
         call QLMME(2,0,COP,COP(1+NN),IS)
         
         call dgemm('n','n',N,N,N,a,COP,N,GAMMA(NG),N,b,AUX,N)
         call dgemm('t','n',N,N,N,a,GAMMA(NG),N,AUX,N,b,AUX(NN+1),N)
         call cpydiag(AUX(NN+1),N,CANNON(IE1,4))
         
         call QLMME(2,2,COP(1+2*NN),COP(1+3*NN),IS)
         
         call dgemm('n','n',N,N,N,a,COP(1+2*NN),N,GAMMA(NG),N,b,AUX,N)
         call dgemm('t','n',N,N,N,a,GAMMA(NG),N,AUX,N,b,AUX(NN+1),N)
         call cpydiag(AUX(NN+1),N,CANNON(IE1,5))
c ---------------------------------------------------------neg signature
         NG = NG + NN
         NK = NK + NN
         call dsyev('V','U',N,GAMMA(NG),N,CANNON(IE2,1),AUX,NAUX,INFO)

         a = 1.0d+00
         b = 0.0d+00
	 
         call dgemm('n','n',N,N,N,a,ROM(NG),N,GAMMA(NG),N,b,AUX,N)
         call dgemm('t','n',N,N,N,a,GAMMA(NG),N,AUX,N,b,AUX(NN+1),N)
         call cpydiag(AUX(NN+1),N,CANNON(IE2,2))
c         
         call dgemm('n','n',N,N,N,a,DELTA(NK),N,GAMMA(NG),N,b,AUX,N)
         call dgemm('t','n',N,N,N,a,GAMMA(NG),N,AUX,N,b,AUX(NN+1),N)
         call cpydiag(AUX(NN+1),N,CANNON(IE2,3))
         
         call dgemm('n','n',N,N,N,a,COP(1+NN),N,GAMMA(NG),N,b,AUX,N)
         call dgemm('t','n',N,N,N,a,GAMMA(NG),N,AUX,N,b,AUX(NN+1),N)
         call cpydiag(AUX(NN+1),N,CANNON(IE2,4))
                  
         call dgemm('n','n',N,N,N,a,COP(1+3*NN),N,GAMMA(NG),N,b,AUX,N)
         call dgemm('t','n',N,N,N,a,GAMMA(NG),N,AUX,N,b,AUX(NN+1),N)
         call cpydiag(AUX(NN+1),N,CANNON(IE2,5))
c+---------------------------------------------------------------------+
         
	
	 
      end do	 
      
      call YOCCOZ(UV,EQP,COP,AUX)      
      write (11) bx,by,bz,UV,EQP
      

      write(69) CANNON
      
      stop
      end      
      
      Double precision function etanew(dtest,shh,shhold,ett)
      Implicit real*8 (A-H,O-Z)
      common /CITER/ epsg,etamax,etamin,dmax,dmin,tshh,maxiter
      etfac=1.00
      if (dtest.gt.dmax.and.shh.lt.tshh)    etfac= 2.00d+00
      if (dtest.lt.dmin)                    etfac= 0.50d+00
      if (shh.gt.shhold.and.dtest.lt.0.9)  etfac= 0.25d+00
      eta=ett*etfac
      if (eta.gt.etamax+1.d-05) eta= etamax
      if (eta.lt.etamin-1.d-05) eta= etamin
      etanew = eta
c
      return
      end
      SUBROUTINE HOINIT(UV,AZ,AN)
      IMPLICIT REAL*8(A-H,O-Z)
      Include 'DIMTRIAX'
      Include 'MAXDIM'
      
      Dimension UV(NUV)
      Dimension ee(NP),eo(NM),v2e(NP,2),v2o(NM,2)
      
      Dimension ie(NP),io(NM)
      
      Common /CQN/    IQMAX,IQNE(3,NPMAX),IQNO(3,NMMAX),IFILL
      Common /OSCLEN/ bx,by,bz
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC1/ NNU(4),NNV(4),NNRO(4),NNKA(4)

c
      if(IFILL.eq.0) call setQN()
      
      hc2=197.3286**2
      manuc=938.9
      omex=hc2/manuc/bx**2
      omey=hc2/manuc/by**2
      omez=hc2/manuc/bz**2

      write(6,*) ' HOINIT ***************** ',omex,omey,omez
       write(6,*) ' HOINIT ***************** ',AZ,AN
      write(6,*) ' HOINIT ***************** ',ND(1),ND(2),ND(3),ND(4)
     
c---------------------------------------------------------------------- +
c     TRIAXIAL HARMONIC OSCILLATOR
c+---------------------------------------------------------------------+
c                  Even parity
c
c  --- Note that the quantum number given by iqne are the real one + 1
c+---------------------------------------------------------------------+
      do ia=1,NP
         nxa = iqne(1,ia)
         nya = iqne(2,ia)
         nza = iqne(3,ia)
         
         ee(ia) = omex*(dfloat(nxa-1)+0.5d+00) +
     &            omey*(dfloat(nya-1)+0.5d+00) +    
     &            omez*(dfloat(nza-1)+0.5d+00)
      end do       ! ia loop
c
c+---------------------------------------------------------------------+
c                  Odd parity
c
c  --- Note that the quantum number given by iqne are the real one + 1
c+---------------------------------------------------------------------+
c
      do ia=1,NM
         nxa = iqno(1,ia)
         nya = iqno(2,ia)
         nza = iqno(3,ia)
         
         eo(ia) = omex*(dfloat(nxa-1)+0.5d+00) +
     &            omey*(dfloat(nya-1)+0.5d+00) +    
     &            omez*(dfloat(nza-1)+0.5d+00)
      end do       ! ia loop
c



c
      call indexx(np,ee,ie)
      call indexx(nm,eo,io)
      
      NZ = AZ + 0.1d+00
      NN = AN + 0.1d+00
c      
      ke=1
      ko=1
      do in=1,NZ/2
        if(ee(ie(ke)).le.eo(io(ko)))then
            v2e(ie(ke),1)=1.0d0
c            write(6,*) ' e ',ke,ee(ie(ke))
            ke=ke+1
         else
            v2o(io(ko),1)=1.0d0
c            write(6,*) 'o ',ko,eo(io(ko))
            ko=ko+1
         end if
       end do
c       
      alambdap=max(ee(ie(ke-1)),eo(io(ko-1)))
      write(6,*) ' ALAMBDA P ',alambdap  
      
      ke=1
      ko=1
      do in=1,NN/2
        if(ee(ie(ke)).le.eo(io(ko)))then
            v2e(ie(ke),2)=1.0d0
            ke=ke+1
         else
            v2o(io(ko),2)=1.0d0
            ko=ko+1
         end if
       end do
       
      alambdan=max(ee(ie(ke-1)),eo(io(ko-1)))
      write(6,*) ' ALAMBDA P ',alambdap  
c
c
c      
      DeltaP = 1.0 d+00
      DeltaN = 1.0 d+00
      
      ZZ = 0.0d+00
      do ia=1,NP
         eep = ee(ia)-alambdap
         een = ee(ia)-alambdan
         v2e(ia,1) = 0.5d+00*(1.0d+00 - eep/dsqrt(eep**2+DeltaP**2))
         v2e(ia,2) = 0.5d+00*(1.0d+00 - een/dsqrt(een**2+DeltaN**2))
         ZZ = ZZ + v2e(ia,1)
c         write(6,*) v2e(ia,1)
      end do
     
      do ia=1,NM
         eep = eo(ia)-alambdap
         een = eo(ia)-alambdan
         v2o(ia,1) = 0.5d+00*(1.0d+00 - eep/dsqrt(eep**2+DeltaP**2))
         v2o(ia,2) = 0.5d+00*(1.0d+00 - een/dsqrt(een**2+DeltaN**2))
         ZZ = ZZ + v2o(ia,1)
c         write(6,*) v2o(ia,1)
      end do
   
      write(6,*) ' HOINIT ***** ',ZZ
c
c     CALCULATION OF U AND V MATRICES
c
c          UV => U1PP U2PP U1PM U2PM  V1PP V2PP V1PM V2PM     
c                U1NP U2NP U1NM U2NM  V1NP V2NP V1NM V2NM 
    
      do i=1,NUV
         UV(i) = 0.0d+00
      end do
            
      do ia=1,NP
         vvp = dsqrt(v2e(ia,1))
         uup = dsqrt(max(1.0d+00-v2e(ia,1),0.0d+00))
         
         uv(NNU(1)       +(ia-1)*np+ia)= uup
         uv(NNU(1)+ND2(1)+(ia-1)*np+ia)= uup
         uv(NNV(1)       +(ia-1)*np+ia)= vvp
         uv(NNV(1)+ND2(1)+(ia-1)*np+ia)=-vvp
         
         vvn = dsqrt(v2e(ia,2))
         uun = dsqrt(max(1.0d+00-v2e(ia,2),0.0d+00))
         
         uv(NNU(3)       +(ia-1)*np+ia)= uun
         uv(NNU(3)+ND2(3)+(ia-1)*np+ia)= uun
         uv(NNV(3)       +(ia-1)*np+ia)= vvn
         uv(NNV(3)+ND2(3)+(ia-1)*np+ia)=-vvn
      end do

       do ia=1,NM
         vvp = dsqrt(v2o(ia,1))
         uup = dsqrt(max(1.0d+00-v2o(ia,1),0.0d+00))
         
         uv(NNU(2)       +(ia-1)*nm+ia)= uup
         uv(NNU(2)+ND2(2)+(ia-1)*nm+ia)= uup
         uv(NNV(2)       +(ia-1)*nm+ia)= vvp
         uv(NNV(2)+ND2(2)+(ia-1)*nm+ia)=-vvp
         
         vvn = dsqrt(v2o(ia,2))
         uun = dsqrt(max(1.0d+00-v2o(ia,2),0.0d+00))
         
         uv(NNU(4)       +(ia-1)*nm+ia)= uun
         uv(NNU(4)+ND2(4)+(ia-1)*nm+ia)= uun
         uv(NNV(4)       +(ia-1)*nm+ia)= vvn
         uv(NNV(4)+ND2(4)+(ia-1)*nm+ia)=-vvn
      end do
      
      return
      end
c+---------------------------------------------------------------------+
c|   Numerical Recipes  subroutine INDEXX                              |
c|   (C) Copr. 1986-92 Numerical Recipes Software.                     |
c+---------------------------------------------------------------------+
      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      DOUBLE PRECISION arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      DOUBLE PRECISION a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
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
      Subroutine kin2bc(ROM,AKAP,GAMMA,DELTA)
c
      Implicit real*8 (A-H,O-Z)
      logical vcom,vcom2
c ...............................................................
      Include 'DIMTRIAX' ! NP and NM
      Include 'MAXDIM'
c ...............................................................
c
      Dimension Sq(NSQ)
      Dimension IQNEO(3,2,NP),IQNOE(3,2,NM)
c
c ----------------------------------------------------------------------
      Dimension ROM  (*), AKAP (*)      
      Dimension GAMMA(*), DELTA(*)
c      
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC1/ NNU(4),NNV(4),NNRO(4),NNKA(4)
c
      Common /CQN/    IQMAX,IQNE(3,NPMAX),IQNO(3,NMMAX),IFILL
      Common /HFBOPT/Amass,vcom,vcom2,icouech
      Common /OSCLEN/ bx,by,bz
c
      call setQNK2B(IQNEO,IQNOE)
c
      if(NSQ.lt.IQMAX+1) write(6,*) ' problems in KIN2BC ',NSQ,IQMAX
      do ii=1,NSQ
         sq(ii) = Dsqrt(Dfloat(ii))
      end do
c      
      fx = 20.734863D+00/(bx*bx)/amass
      fy = 20.734863D+00/(by*by)/amass
      fz = 20.734863D+00/(bz*bz)/amass
      
      do it=1,4
         ipar = 2*mod(it,2)-1 
	 itb=it+ipar ! the opposite parity index
         N  = ND (it )
	 NB = ND (itb)
	 
c      write(6,*) ' KIN2BC ',it,itb,ipar,N,NB
      
      do i=1,N
         if(ipar.eq.+1) then
            nxi  = iqne(1,i)
            nyi  = iqne(2,i)
            nzi  = iqne(3,i)
	    ixp1 = IQNEO(1,1,i)
	    ixm1 = IQNEO(1,2,i)
	    iyp1 = IQNEO(2,1,i)
	    iym1 = IQNEO(2,2,i)
	    izp1 = IQNEO(3,1,i)
	    izm1 = IQNEO(3,2,i)
	 else 
            nxi  = iqno(1,i)
            nyi  = iqno(2,i)
            nzi  = iqno(3,i)
	    ixp1 = IQNOE(1,1,i)
	    ixm1 = IQNOE(1,2,i)
	    iyp1 = IQNOE(2,1,i)
	    iym1 = IQNOE(2,2,i)
	    izp1 = IQNOE(3,1,i)
	    izm1 = IQNOE(3,2,i)
	 end if
         do j=1,N
	    if(ipar.eq.+1) then
               nxj  = iqne(1,j)
               nyj  = iqne(2,j)
               nzj  = iqne(3,j)
	       jxp1 = IQNEO(1,1,j)
	       jxm1 = IQNEO(1,2,j)
	       jyp1 = IQNEO(2,1,j)
	       jym1 = IQNEO(2,2,j)
	       jzp1 = IQNEO(3,1,j)
	       jzm1 = IQNEO(3,2,j)
	    else
               nxj  = iqno(1,j)
               nyj  = iqno(2,j)
               nzj  = iqno(3,j)
	       jxp1 = IQNOE(1,1,j)
	       jxm1 = IQNOE(1,2,j)
	       jyp1 = IQNOE(2,1,j)
	       jym1 = IQNOE(2,2,j)
	       jzp1 = IQNOE(3,1,j)
	       jzm1 = IQNOE(3,2,j)
	    end if
	    
	    kg1 = NNRO(it) + i + (j-1)*N -1
	    kg2 = NNRO(it) + N*N + i + (j-1)*N -1
	    kd  = NNKA(it) + i + (j-1)*N -1
	    
	    KBro1 = NNRO(itb) -1
	    KBro2 = NNRO(itb) + NB*NB -1
	    KBkap = NNKA(itb) -1
	    
            phasx = dfloat(1-2*Mod((nxi+nxj+nyi+nyj+1),2))

            GAMMA(kg1) = 0.0d+00
            GAMMA(kg2) = 0.0d+00
            DELTA(kd ) = 0.0d+00

c
c ---------------------------------------------------------- P                         
c                                                             x
c
             g1 = 0.0d+00
	     g2 = 0.0d+00
	     d  = 0.0d+00
             if((ixp1.ne.0).and.(jxp1.ne.0)) then
               SS = -Sq(nxi)*Sq(nxj)*phasx*fx
               g1 = g1 + SS*ROM (KBro2+ixp1+(jxp1-1)*NB)
               g2 = g2 + SS*ROM (KBro1+ixp1+(jxp1-1)*NB)
               d  = d  + SS*AKAP(KBkap+jxp1+(ixp1-1)*NB)
             end if
             if((ixp1.ne.0).and.(jxm1.ne.0)) then
               SS = Sq(nxi)*Sq(nxj-1)*phasx*fx
               g1 = g1 + SS*ROM (KBro2+ixp1+(jxm1-1)*NB)
               g2 = g2 + SS*ROM (KBro1+ixp1+(jxm1-1)*NB)
               d  = d  + SS*AKAP(KBkap+jxm1+(ixp1-1)*NB)
             end if
             if((ixm1.ne.0).and.(jxp1.ne.0)) then
               SS = Sq(nxi-1)*Sq(nxj)*phasx*fx
               g1 = g1 + SS*ROM (KBro2+ixm1+(jxp1-1)*NB)
               g2 = g2 + SS*ROM (KBro1+ixm1+(jxp1-1)*NB)
               d  = d  + SS*AKAP(KBkap+jxp1+(ixm1-1)*NB)
             end if
             if((ixm1.ne.0).and.(jxm1.ne.0)) then
               SS = -Sq(nxi-1)*Sq(nxj-1)*phasx*fx
               g1 = g1 + SS*ROM (KBro2+ixm1+(jxm1-1)*NB)
               g2 = g2 + SS*ROM (KBro1+ixm1+(jxm1-1)*NB)
               d  = d  + SS*AKAP(KBkap+jxm1+(ixm1-1)*NB)
             end if
c
c ---------------------------------------------------------- P
c                                                             y
             if((iyp1.ne.0).and.(jyp1.ne.0)) then
               SS = Sq(nyi)*Sq(nyj)*fy
               g1 = g1 + SS*ROM (KBro1+iyp1+(jyp1-1)*NB)
               g2 = g2 + SS*ROM (KBro2+iyp1+(jyp1-1)*NB)
               d  = d  + SS*AKAP(KBkap+iyp1+(jyp1-1)*NB)
             end if
             if((iyp1.ne.0).and.(jym1.ne.0)) then
               SS = Sq(nyi)*Sq(nyj-1)*fy
               g1 = g1 + SS*ROM (KBro1+iyp1+(jym1-1)*NB)
               g2 = g2 + SS*ROM (KBro2+iyp1+(jym1-1)*NB)
               d  = d  + SS*AKAP(KBkap+iyp1+(jym1-1)*NB)
             end if
             if((iym1.ne.0).and.(jyp1.ne.0)) then
               SS = Sq(nyi-1)*Sq(nyj)*fy
               g1 = g1 + SS*ROM (KBro1+iym1+(jyp1-1)*NB)
               g2 = g2 + SS*ROM (KBro2+iym1+(jyp1-1)*NB)
               d  = d  + SS*AKAP(KBkap+iym1+(jyp1-1)*NB)
             end if
             if((iym1.ne.0).and.(jym1.ne.0)) then
               SS = Sq(nyi-1)*Sq(nyj-1)*fy
               g1 = g1 + SS*ROM (KBro1+iym1+(jym1-1)*NB)
               g2 = g2 + SS*ROM (KBro2+iym1+(jym1-1)*NB)
               d  = d  + SS*AKAP(KBkap+iym1+(jym1-1)*NB)
             end if
c
c ---------------------------------------------------------- P
c                                                             z
             if((izp1.ne.0).and.(jzp1.ne.0)) then
               SS = Sq(nzi)*Sq(nzj)*fz
               g1 = g1 + SS*ROM (KBro1+izp1+(jzp1-1)*NB)
               g2 = g2 + SS*ROM (KBro2+izp1+(jzp1-1)*NB)
               d  = d  + SS*AKAP(KBkap+izp1+(jzp1-1)*NB)
             end if
             if((izp1.ne.0).and.(jzm1.ne.0)) then
               SS = -Sq(nzi)*Sq(nzj-1)*fz
               g1 = g1 + SS*ROM (KBro1+izp1+(jzm1-1)*NB)
               g2 = g2 + SS*ROM (KBro2+izp1+(jzm1-1)*NB)
               d  = d  + SS*AKAP(KBkap+izp1+(jzm1-1)*NB)
             end if
             if((izm1.ne.0).and.(jzp1.ne.0)) then
               SS = -Sq(nzi-1)*Sq(nzj)*fz
               g1 = g1 + SS*ROM (KBro1+izm1+(jzp1-1)*NB)
               g2 = g2 + SS*ROM (KBro2+izm1+(jzp1-1)*NB)
               d  = d  + SS*AKAP(KBkap+izm1+(jzp1-1)*NB)
             end if
             if((izm1.ne.0).and.(jzm1.ne.0)) then
               SS = Sq(nzi-1)*Sq(nzj-1)*fz
               g1 = g1 + SS*ROM (KBro1+izm1+(jzm1-1)*NB)
               g2 = g2 + SS*ROM (KBro2+izm1+(jzm1-1)*NB)
               d  = d  + SS*AKAP(KBkap+izm1+(jzm1-1)*NB)
             end if
	     gamma(kg1) = g1
	     gamma(kg2) = g2
	     delta(kd ) = d
         end do     ! j
      end do        ! i
      end do        ! it
c
      return
      end
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
      ct= -20.734863d+00/(bx*bx)
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
      Subroutine RESMV (RO,OP)
      Implicit real*8 (A-H,O-Z)
      Character*8 TEE
      Parameter (NEEE=30) ! dimension of the output matrix
      
      Dimension OP(*) ! scratch matrix to hold the operator's matrix elements
      Dimension RO(*)
      
      Common /EEEEEE/ TEE(NEEE),EEE(7,NEEE)
      Common /IIIEEE/ mkin,mehf,mepair,mehfb,medd,merea,mecou,
     &                mn,mjx,msx,mq20,mq22,mq40,mneck,mr2,
     &                mbet2,mgamm,mbet4,maxrat,
     &                mfn,mjx2,mjy2,mjz2
     
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC1/ NNU(4),NNV(4),NNRO(4),NNKA(4)
     
      NP2 = ND2(1)
      NM2 = ND2(2)
      
      
      do it=1,4
      
	 IRO1 = NNRO(it)
	 IRO2 = NNRO(it) + ND2(it)
	 
	 EEE(it,mn) = trace(RO(IRO1),ND(it)) + trace(RO(IRO2),ND(it))
	 
      end do
c
c   < Jx > (time odd -> the minus sign) 
c	 

      call momang('X',OP(1),OP(1+NP2))
      
      do it=1,4
      
	 IRO1 = NNRO(it)
	 IRO2 = NNRO(it) + ND2(it)
	 
         ICOP = (1-Mod(it,2))*NP2+1
         EEE(it,mjx ) = ddot(ND2(it),RO(IRO1),1,OP(ICOP),1)-
     &                  ddot(ND2(it),RO(IRO2),1,OP(ICOP),1)
         
      end do
c
c   < Sx > (time odd -> the minus sign) 
c	 

      call espin('X',OP(1),OP(1+NP2))
      
      do it=1,4
      
	 IRO1 = NNRO(it)
	 IRO2 = NNRO(it) + ND2(it)
	 
         ICOP = (1-Mod(it,2))*NP2+1
         EEE(it,msx ) = ddot(ND2(it),RO(IRO1),1,OP(ICOP),1)-
     &                  ddot(ND2(it),RO(IRO2),1,OP(ICOP),1)
         
      end do
c
c   < Qlm > (time even -> the plus sign) 
c	 

c+-------------------------------------------------------------+
c|                              Q                              |
c|                               20                            |
c+-------------------------------------------------------------+     
      call QLMME(2,0,OP,OP(1+NP2),ISI)
      
      do it=1,4
      
	 IRO1 = NNRO(it)
	 IRO2 = NNRO(it) + ND2(it)
	 
         ICOP = (1-Mod(it,2))*NP2+1
         EEE(it,mq20) = ddot(ND2(it),RO(IRO1),1,OP(ICOP),1)+
     &                  ddot(ND2(it),RO(IRO2),1,OP(ICOP),1)
         
      end do
      
c+-------------------------------------------------------------+
c|                              Q                              |
c|                               22                            |
c+-------------------------------------------------------------+     
      call QLMME(2,2,OP,OP(1+NP2),ISI)
      
      do it=1,4
      
	 IRO1 = NNRO(it)
	 IRO2 = NNRO(it) + ND2(it)
	 
         ICOP = (1-Mod(it,2))*NP2+1
         EEE(it,mq22) = ddot(ND2(it),RO(IRO1),1,OP(ICOP),1)+
     &                  ddot(ND2(it),RO(IRO2),1,OP(ICOP),1)
         
      end do

c+-------------------------------------------------------------+
c|                              Q                              |
c|                               40                            |
c+-------------------------------------------------------------+     
      call QLMME(4,0,OP,OP(1+NP2),ISI)
      
      do it=1,4
      
	 IRO1 = NNRO(it)
	 IRO2 = NNRO(it) + ND2(it)
	 
         ICOP = (1-Mod(it,2))*NP2+1
         EEE(it,mq40) = ddot(ND2(it),RO(IRO1),1,OP(ICOP),1)+
     &                  ddot(ND2(it),RO(IRO2),1,OP(ICOP),1)
         
      end do
      
c+-------------------------------------------------------------+
c|                             2                               |
c|                            r                                |
c+-------------------------------------------------------------+     
      call RLME(2,OP,OP(1+NP2))
      
      do it=1,4
      
	 IRO1 = NNRO(it)
	 IRO2 = NNRO(it) + ND2(it)
	 
         ICOP = (1-Mod(it,2))*NP2+1
         EEE(it,mr2 ) = ddot(ND2(it),RO(IRO1),1,OP(ICOP),1)+
     &                  ddot(ND2(it),RO(IRO2),1,OP(ICOP),1)
         
      end do
c
c     Mean Square radius
c
      
      zp  = EEE(1,mn ) + EEE(2,mn )
      zn  = EEE(3,mn ) + EEE(4,mn )
      r2p = EEE(1,mr2) + EEE(2,mr2)
      r2n = EEE(3,mr2) + EEE(4,mr2)
      r2t = (r2p+r2n)/(zp+zn)
      EEE(5,mr2) = dsqrt(r2p/zp)
      EEE(6,mr2) = dsqrt(r2n/zn)
      EEE(7,mr2) = dsqrt(r2t)
      
      do it=1,4
         EEE(it,mr2) = dsqrt(EEE(it,mr2)/EEE(it,mn))
      end do
      
c+-------------------------------------------------------------+
c|                   N E C K                                   |
c|                                                             |
c+-------------------------------------------------------------+     
      call NECKME(OP,OP(1+NP2))
      
      do it=1,4
      
	 IRO1 = NNRO(it)
	 IRO2 = NNRO(it) + ND2(it)
c	 
         ICOP = (1-Mod(it,2))*NP2+1
         EEE(it,mneck) = ddot(ND2(it),RO(IRO1),1,OP(ICOP),1)+
     &                   ddot(ND2(it),RO(IRO2),1,OP(ICOP),1)
c         
      end do
      
      return
      end
c---------------------------------------------------------------------------      
      Subroutine RESFLUC (UV,OP,AUX)
      Implicit real*8 (A-H,O-Z)
      Dimension OP (*)
      Dimension UV  (*)
      Dimension AUX (*)
      Character*8 TEE
      Parameter (NEEE=30) ! dimension of the output matrix
      
      Common /EEEEEE/ TEE(NEEE),EEE(7,NEEE)
      Common /IIIEEE/ mkin,mehf,mepair,mehfb,medd,merea,mecou,
     &                mn,mjx,msx,mq20,mq22,mq40,mneck,mr2,
     &                mbet2,mgamm,mbet4,maxrat,
     &                mfn,mjx2,mjy2,mjz2
      
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC1/ NNU(4),NNV(4),NNRO(4),NNKA(4)
           
      NP2 = ND2(1)
      NM2 = ND2(2)
     
      do it=1,4
      
         IU = NNU(it)
         IV = NNV(it)
	 N  = ND(it)
	 N1 = ND(it)-NBLOCK(it)
	 N2 = ND(it)+NBLOCK(it)
	 N12= N1*N2

         call SN20 (UV(IU),UV(IV),N1,N2,N,AUX)
	 
	 EEE(it,mfn) = ddot(N12,AUX,1,AUX,1)
	 
      end do
      
      call momang('X',OP(1),OP(1+NP2))
	 
      do it=1,4
      
         IU = NNU(it)
         IV = NNV(it)
	 N  = ND(it)
	 N1 = ND(it)-NBLOCK(it)
	 N2 = ND(it)+NBLOCK(it)
	 N12= N1*N2

      
         ICOP = (1-Mod(it,2))*NP2+1
         call SO20P(UV(IU),UV(IV),N1,N2,N,
     &                OP(ICOP),AUX,AUX(1+N*N2),-1,1)
	 
	 EEE(it,mjx2) = ddot(N12,AUX(1+N*N2),1,AUX(1+N*N2),1)
      end do
      
      call momang('Y',OP(1),OP(1+NP2))
	 
      do it=1,4
      
         IU = NNU(it)
         IV = NNV(it)
	 N  = ND(it)
	 N1 = ND(it)-NBLOCK(it)
	 N2 = ND(it)+NBLOCK(it)
	 N11= N1*N1
	 N22= N2*N2
	 I1 = 1 + 2*N*N         ! Q20_1
	 I2 = 1 + 2*N*N+N1*N1   ! Q20_2

      
         ICOP = (1-Mod(it,2))*NP2+1
         call SO20M(UV(IU),UV(IV),N1,N2,N,
     &                OP(ICOP),AUX,AUX(I1),AUX(I2),1,-1)
	 
	 EEE(it,mjy2) = 0.5d+00*(
     &  	 ddot(N11,AUX(I1),1,AUX(I1),1)+
     &         	 ddot(N22,AUX(I2),1,AUX(I2),1) )
      end do
      
      call momang('Z',OP(1),OP(1+NP2))
	 
      do it=1,4
      
         IU = NNU(it)
         IV = NNV(it)
	 N  = ND(it)
	 N1 = ND(it)-NBLOCK(it)
	 N2 = ND(it)+NBLOCK(it)
	 N11= N1*N1
	 N22= N2*N2
	 I1 = 1 + 2*N*N         ! Q20_1
	 I2 = 1 + 2*N*N+N1*N1   ! Q20_2

      
         ICOP = (1-Mod(it,2))*NP2+1
         call SO20M(UV(IU),UV(IV),N1,N2,N,
     &                OP(ICOP),AUX,AUX(I1),AUX(I2),-1,1)
	 
	 EEE(it,mjz2) = 0.5d+00*(
     &  	 ddot(N11,AUX(I1),1,AUX(I1),1)+
     &         	 ddot(N22,AUX(I2),1,AUX(I2),1) )
      end do
      
      return
      end 
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
      Subroutine setfact()
      Implicit real*8 (A-H,O-Z)
      Include 'MAXDIM'
      Common /fact/ dfact(Nfac),ddfact(Nfac),Nfacm
c+--------------------------------------------------------------------+
c|    dfact(I)  = ln( (I-1)! )   |   ln(I!) = dfact (I+1)             |
c|    ddfact(I) = ln( (2I-3)!! ) |   ln(I!!)= ddfact((I+3)/2)         |
c+--------------------------------------------------------------------+
      Nfacm = Nfac
      dfact(1) = 0.0d+00
      dfact(2) = 0.0d+00
      ddfact(1)= 0.0d+00
      ddfact(2)= 0.0d+00
      do ifac =3,nfacm
           dfact (ifac) = dlog(dfloat(ifac-1))  + dfact (ifac-1)
           ddfact(ifac) = dlog(dfloat(2*ifac-3))+ ddfact(ifac-1)
      end do
      return 
      end
c------------------------------------------------------------------------      
      Subroutine setQN()
      Implicit real*8(a-h,o-z)
      Include 'DIMPERM'
      Include 'MAXDIM'
      Common /CQN/    IQMAX,IQNE(3,NPMAX),IQNO(3,NMMAX),IFILL
c
      IFILL = 1
      
      jz1e = 0
      jz1o = 0
      do nx1=1,ixmax
        nym1 = imy(nx1)
          do ny1 = 1,nym1
            nxy1 = nx1+ny1
            nzi1e= 1+Mod(nxy1,2)
            nzi1o= 2-Mod(nxy1,2)
            nzm1 = imz(nx1,ny1)
            if(nzi1e.le.nzm1) then
              do nz1=nzi1e,nzm1,2
                 jz1e = jz1e + 1
                 iqne(1,jz1e) = nx1
                 iqne(2,jz1e) = ny1
                 iqne(3,jz1e) = nz1
              end do
            end if
            if(nzi1o.le.nzm1) then
              do nz1=nzi1o,nzm1,2
                 jz1o = jz1o + 1
                 iqno(1,jz1o) = nx1
                 iqno(2,jz1o) = ny1
                 iqno(3,jz1o) = nz1
              end do
            end if
          end do ! ny1
      end do     ! nx1
c
      IQMAX = MAX(ixmax,iymax,izmax)
      return
      end
c
c
c
      Subroutine setQNK2B(IQNEO,IQNOE)
      Implicit real*8(a-h,o-z)
      Include 'DIMPERM'
      Include 'MAXDIM'
      Dimension IQNEO(3,2,*),IQNOE(3,2,*)
      Common /CQN/    IQMAX,IQNE(3,NPMAX),IQNO(3,NMMAX),IFILL
      Common /DIMS/   NP,NM,NROP,NROM
c
      if(IFILL.ne.1) call setQN()
c
c  iqneo(*,*,i) i=1,NP index in the odd parity states
c
c      
      do i=1,NP
      
      	nxe= iqne(1,i)
      	nye= iqne(2,i)
      	nze= iqne(3,i)
	
	do k=1,3
	  do l=1,2
	     iqneo(k,l,i) = 0
	  end do
	end do
	
	do j=1,NM
	   nxo = iqno(1,j)
	   nyo = iqno(2,j)
	   nzo = iqno(3,j)
	   
      if((nxo.eq.nxe+1).and.(nyo.eq.nye).and.(nzo.eq.nze)) then
         iqneo(1,1,i)=j
      else if((nxo.eq.nxe-1).and.(nyo.eq.nye).and.(nzo.eq.nze)) then
         iqneo(1,2,i)=j
      else if((nxo.eq.nxe).and.(nyo.eq.nye+1).and.(nzo.eq.nze)) then
         iqneo(2,1,i)=j
      else if((nxo.eq.nxe).and.(nyo.eq.nye-1).and.(nzo.eq.nze)) then
         iqneo(2,2,i)=j
      else if((nxo.eq.nxe).and.(nyo.eq.nye).and.(nzo.eq.nze+1)) then
         iqneo(3,1,i)=j
      else if((nxo.eq.nxe).and.(nyo.eq.nye).and.(nzo.eq.nze-1)) then
         iqneo(3,2,i)=j
      end if
         end do
      end do 
      
      do i=1,NM
      
      	nxo= iqno(1,i)
      	nyo= iqno(2,i)
      	nzo= iqno(3,i)
	
	do k=1,3
	  do l=1,2
	     iqnoe(k,l,i) = 0
	  end do
	end do
	
	do j=1,NP
	   nxe = iqne(1,j)
	   nye = iqne(2,j)
	   nze = iqne(3,j)
	   
      if((nxo+1.eq.nxe).and.(nyo.eq.nye).and.(nzo.eq.nze)) then
         iqnoe(1,1,i)=j
      else if((nxo-1.eq.nxe).and.(nyo.eq.nye).and.(nzo.eq.nze)) then
         iqnoe(1,2,i)=j
      else if((nxo.eq.nxe).and.(nyo+1.eq.nye).and.(nzo.eq.nze)) then
         iqnoe(2,1,i)=j
      else if((nxo.eq.nxe).and.(nyo-1.eq.nye).and.(nzo.eq.nze)) then
         iqnoe(2,2,i)=j
      else if((nxo.eq.nxe).and.(nyo.eq.nye).and.(nzo+1.eq.nze)) then
         iqnoe(3,1,i)=j
      else if((nxo.eq.nxe).and.(nyo.eq.nye).and.(nzo-1.eq.nze)) then
         iqnoe(3,2,i)=j
      end if
         end do
      end do 
      
      return
      end
      
c +--------------------------------------------------------------------+
c    Subroutine TIMEIT: To keep track of the CPU time spent in the
c    calculation.
c +--------------------------------------------------------------------+
      Subroutine TIMEIT(iflag,index,text)
      Implicit real*8 (A-H,O-Z)
      Parameter (Nmax=30) ! Maximum number of items to be accounted for
      Character*16 text
      Character*16 texto
      Common /CTIMEIT/ texto(Nmax),prevtim(Nmax),accutim(Nmax),
     *       Ncalls(Nmax),ntot
c
c    index=0     Set storage to 0
c
      if(index.eq.0) then
         do I=1,Nmax
            prevtim(I) = 0.0d+00
            accutim(I) = 0.0d+00
            Ncalls (I) = 0
            texto  (I) = '****************'
         end do
         ntot = 0
      else if ((index.gt.0).and.(index.le.Nmax)) then
c
c    index .ne. 0
c
         if(iflag.eq.0) then
            prevtim(index) = second()
         else if (iflag.eq.1) then
            ttime = second()
            accutim(index) = accutim(index) + ttime - prevtim(index)
            if(ncalls(index).eq.0) texto(index) = text
            ncalls(index) = ncalls(index) + 1
            if(index.gt.ntot) ntot = index
         else
            write(6,*) ' IFLAG value in TIMEIT unknown '
            write(6,*) ' IFLAG ',IFLAG,' INDEX ',INDEX,TEXT
         end if
c
c    index out of range
c
      else
         write(6,*) ' INDEX value in TIMEIT out of range '
         write(6,*) ' IFLAG ',IFLAG,' INDEX ',INDEX,TEXT
      end if
c
      return
      end
c +--------------------------------------------------------------------+
c    Subroutine TIMEREP: To write a report of the cpu time spent in
c    the CHFBJX program
c +--------------------------------------------------------------------+
      Subroutine TIMEREP
      Implicit real*8 (A-H,O-Z)
      Parameter (Nmax=30) ! Maximum number of items to be accounted for
      Character*16 texto
      Common /CTIMEIT/ texto(Nmax),prevtim(Nmax),accutim(Nmax),
     *       Ncalls(Nmax),ntot
      COMMON /GRENZE/ EPSG,SH,EPSN,SN,EPSI,SI,EPSO,SO,SSS,ZZ,LAN,IG,MAXG
      Common/infob/aq0inf,ap0inf,arn0inf,nxmax,nymax,nzmax
c
      nout = 50
      write(nout,100) aq0inf,ap0inf,arn0inf,ig
  100 format( 10x,'  CPU TIME REPORT OF THE CHFBJX PROGRAM ',/,
     *     5x,10('====='),//,
     *    10x,' Q0 = ',f5.2,' P0 = ',f5.2,' N0 = ',f5.2,//,
     *    20x,' NUMBER OF ITERATIONS',I5,//)
      write(nout,101)
  101 format ( 5x,'Subroutine name ',' Number of calls ',
     * ' Total time (s) ',' Time per call (s) ',' % of total ',/,5x,
     * 20('-----'),/)
c
c It is assumed that i=1 is always the total time
c
      timetot = accutim(1)
      do I=1,ntot
         ttime = accutim(i)
         timepc= ttime/dfloat(ncalls(i))
         perc  = ttime/timetot *1.D+2
         write(nout,102) texto(i),ncalls(i),ttime,timepc,perc
      end do
  102 format( 5x,A16,8x,i8,4x,f12.2,7x,f12.5,5x,f5.1)
c
      return
      end
      subroutine errout(Mesage,subr,icode)
c ÉÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ»
c º                                                                      º
c º                                                                      º
c º                                                                      º
c º                                                                      º
c º                                                                      º
c º                                                                      º
c º                                                                      º
c ÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼
      character*16 Mesage
      character*8 subr
      integer*4 icode
      common/fin/ ifin
c
      write(6,100) subr,mesage,icode
  100 format( ' From subroutine ',a8,' >>>>>>>    ',a16,/,
     *        ' Code:  ',i4,/)
c
      if(icode.gt.3) then
           write(6,101)
  101      format(' Code greater than 3 .... Abend       ')
           stop
      end if
      if(icode.eq.3) then
           write(6,102)
  102      format(' Code equal 3 ............ Normal end ')
           ifin = 9999
      end if
c
      return
      end
c---------------------------------------------------------------------------      
      Subroutine YOCCOZ (UV,EQP,OP,AUX)
      Implicit real*8 (A-H,O-Z)
      Dimension OP (*)
      Dimension UV  (*),EQP(*)
      Dimension AUX (*)
      Dimension DJ2(3),DEJ2(3)
      Character*8 TEE
      Parameter (NEEE=30) ! dimension of the output matrix
      
      Common /EEEEEE/ TEE(NEEE),EEE(7,NEEE)
      Common /IIIEEE/ mkin,mehf,mepair,mehfb,medd,merea,mecou,
     &                mn,mjx,msx,mq20,mq22,mq40,mneck,mr2,
     &                mbet2,mgamm,mbet4,maxrat,
     &                mfn,mjx2,mjy2,mjz2
     
      Common /YOCRES/ YZ,EZ,YN,EN,YJ,EJ
      
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC1/ NNU(4),NNV(4),NNRO(4),NNKA(4)
      Common /ITLOC4/ NNH20(4),NNH11(4),NNEQP(4)
           
      NP2 = ND2(1)
      NM2 = ND2(2)
     
      do it=1,4
      
         IU = NNU(it)
         IV = NNV(it)
	 N  = ND(it)
	 N1 = ND(it)-NBLOCK(it)
	 N2 = ND(it)+NBLOCK(it)
	 N12= N1*N2
         IE1  = NNEQP(it)
         IE2  = IE1 + N1

         call SN20 (UV(IU),UV(IV),N1,N2,N,AUX)
	 
	 EEE(it,mfn) = ddot(N12,AUX,1,AUX,1)
	 
	 sum = 0.0d+00
	 do i=1,N1
	    do j=1,N2
	       sum=sum+AUX(i+(j-1)*N1)**2*(EQP(IE1+i-1)+EQP(IE2+j-1))
            end do
	 end do
         
	 EEE(it,27) = sum
	 	 
      end do
      
      call momang('X',OP(1),OP(1+NP2))
	 
      do it=1,4
      
         IU = NNU(it)
         IV = NNV(it)
	 N  = ND(it)
	 N1 = ND(it)-NBLOCK(it)
	 N2 = ND(it)+NBLOCK(it)
	 N12= N1*N2
         IE1  = NNEQP(it)
         IE2  = IE1 + N1

      
         ICOP = (1-Mod(it,2))*NP2+1
         call SO20P(UV(IU),UV(IV),N1,N2,N,
     &                OP(ICOP),AUX,AUX(1+N*N2),-1,1)
		 
	 EEE(it,mjx2) = ddot(N12,AUX(1+N*N2),1,AUX(1+N*N2),1)
	 
	 sum = 0.0d+00
	 do i=1,N1
	    ii = I + N*N2
	    do j=1,N2
	       sum=sum+AUX(ii+(j-1)*N1)**2*(EQP(IE1+i-1)+EQP(IE2+j-1))
            end do
	 end do
         
	 EEE(it,28) = sum
      end do
      
      call momang('Y',OP(1),OP(1+NP2))
	 
      do it=1,4
      
         IU = NNU(it)
         IV = NNV(it)
	 N  = ND(it)
	 N1 = ND(it)-NBLOCK(it)
	 N2 = ND(it)+NBLOCK(it)
	 N11= N1*N1
	 N22= N2*N2
	 I1 = 1 + 2*N*N         ! Q20_1
	 I2 = 1 + 2*N*N+N1*N1   ! Q20_2
         IE1  = NNEQP(it)
         IE2  = IE1 + N1

      
         ICOP = (1-Mod(it,2))*NP2+1
         call SO20M(UV(IU),UV(IV),N1,N2,N,
     &                OP(ICOP),AUX,AUX(I1),AUX(I2),1,-1)
     
	 EEE(it,mjy2) = 0.5d+00*(
     &  	 ddot(N11,AUX(I1),1,AUX(I1),1)+
     &         	 ddot(N22,AUX(I2),1,AUX(I2),1) )
	 
	 sum1 = 0.0d+00
	 do i=1,N1
	    ii = I + I1-1
	    do j=1,N1
	       sum1=sum1+AUX(ii+(j-1)*N1)**2*(EQP(IE1+i-1)+EQP(IE1+j-1))
            end do
	 end do
	 
	 sum2 = 0.0d+00
	 do i=1,N2
	    ii = I + I2-1
	    do j=1,N2
	       sum2=sum2+AUX(ii+(j-1)*N2)**2*(EQP(IE2+i-1)+EQP(IE2+j-1))
            end do
	 end do
         
	 EEE(it,29) = 0.5d+00*(sum1+sum2)
      end do
      
      call momang('Z',OP(1),OP(1+NP2))
	 
      do it=1,4
      
         IU = NNU(it)
         IV = NNV(it)
	 N  = ND(it)
	 N1 = ND(it)-NBLOCK(it)
	 N2 = ND(it)+NBLOCK(it)
	 N11= N1*N1
	 N22= N2*N2
	 I1 = 1 + 2*N*N         ! Q20_1
	 I2 = 1 + 2*N*N+N1*N1   ! Q20_2
         IE1  = NNEQP(it)
         IE2  = IE1 + N1

      
         ICOP = (1-Mod(it,2))*NP2+1
         call SO20M(UV(IU),UV(IV),N1,N2,N,
     &                OP(ICOP),AUX,AUX(I1),AUX(I2),-1,1)
	 
	 EEE(it,mjz2) = 0.5d+00*(
     &  	 ddot(N11,AUX(I1),1,AUX(I1),1)+
     &         	 ddot(N22,AUX(I2),1,AUX(I2),1) )
	 
	 sum1 = 0.0d+00
	 do i=1,N1
	    ii = I + I1-1
	    do j=1,N1
	       sum1=sum1+AUX(ii+(j-1)*N1)**2*(EQP(IE1+i-1)+EQP(IE1+j-1))
            end do
	 end do
	 
	 sum2 = 0.0d+00
	 do i=1,N2
	    ii = I + I2-1
	    do j=1,N2
	       sum2=sum2+AUX(ii+(j-1)*N2)**2*(EQP(IE2+i-1)+EQP(IE2+j-1))
            end do
	 end do
         
	 EEE(it,30) = 0.5d+00*(sum1+sum2)
	 
      end do

      DZ2  = EEE(1,mfn) + EEE(2,mfn)
      DN2  = EEE(3,mfn) + EEE(4,mfn)
      DEZ2 = EEE(1, 27) + EEE(2, 27)
      DEN2 = EEE(3, 27) + EEE(4, 27)
      
      YZ   = DZ2**2/DEZ2
      YN   = DN2**2/DEN2
      
      EZ   = 0.5d+00*DZ2/YZ
      EN   = 0.5d+00*DN2/YN

      SDJ2  = 0.0d+00
      SDJ4  = 0.0d+00
      SDEJ2 = 0.0d+00
      do i=1,3
         DJ2 (i) = 0.0d+00
	 DEJ2(i) = 0.0d+00
	 do it=1,4
	    DJ2 (i) = DJ2 (i) + EEE(it,mjx2+i-1)
	    DEJ2(i) = DEJ2(i) + EEE(it,  28+i-1)
	 end do
	 SDJ2  = SDJ2  + DJ2 (i)
	 SDJ4  = SDJ4  + DJ2 (i)**2
	 SDEJ2 = SDEJ2 + DEJ2(i)
	 write(6,101) i,DJ2(i),DEJ2(i)
      end do
   
      YJ   = SDJ4/SDEJ2
      EJ   = 0.5d+00*SDJ2/YJ
      
      write(6,100) YZ,EZ,YN,EN,YJ,EJ
 100  format( " YOCCOZ Z,N,J ",6f12.3)
 101  format( " YOCCOZ J(",I1,") ",2f12.3)
      
      return
      end 
      Subroutine GDIR(GROP,GROM,tz1d,ajx1d,ajy1d,ajz1d,alsomu,dsomu,
     * acou,coum,itz1d,ijx1d,ijy1d,ijz1d,ilsomu,idsomu,
     * dmu,sxyz1p,sxyz1n,sxyz2p,sxyz2n,
     * sxyz3p,sxyz3n,sxyz4p,sxyz4n,sxyzc,scf,Soxyzp,Soxyzn,
     * GSUM0EP,GSUM0EN,GSUM0OP,GSUM0ON,GSUM1OP,GSUM1ON,GSUM3OP,GSUM3ON,
     * THE0EP,THE0EN,THE0OP,THE0ON,THE1OP,THE1ON,THE3OP,THE3ON,
     * DSO0XP,DSO0XN,DSO0XRP,DSO0XRN,DSO0YP,DSO0YN,DSO0ZP,
     * DSO0ZN,DSO1ZP,DSO1ZN,DSO3YP,DSO3YN,
     * gp1p,gp2p,gn1p,gn2p,gp1m,gp2m,gn1m,gn2m)
c
      Implicit real*8 (A-H,O-Z)
      Logical LMZE,LMZO,vmux
      Logical lx12,lxy12,vcom,vcom2
c
      Include 'COMDIM'
c
      Dimension GROP(NROP8),GROM(NROM8)
c
      Dimension tz1d(maxtz1),alsomu(maxlso),dsomu(maxdso)
      Dimension ajx1d(2,maxjx1),ajy1d(2,maxjy1),ajz1d(2,maxjz1)
      Dimension Acou(Nacou),coum(nmax,nmax)
c
      Dimension DMU(15,NDMU)
      Dimension sxyz1p(2,NDMU),sxyz1n(2,NDMU),sxyz2p(2,NDMU),
     *          sxyz2n(2,NDMU),sxyz3p(2,NDMU),sxyz3n(2,NDMU),
     *          sxyz4p(2,NDMU),sxyz4n(2,NDMU)
c
      Dimension sxyzc(NDMU),scf(NDMU)
c
      Dimension Soxyzp(10,NDMU),Soxyzn(10,NDMU)
c
      Dimension GSUM0EP(NDMU),GSUM0EN(NDMU),GSUM0OP(NDMU),
     *          GSUM0ON(NDMU),GSUM1OP(NDMU),GSUM1ON(NDMU),
     *          GSUM3OP(NDMU),GSUM3ON(NDMU)
c
      Dimension THE0EP(NDTHET),THE0EN(NDTHET),THE0OP(NDTHET),
     *          THE0ON(NDTHET),THE1OP(NDTHET),THE1ON(NDTHET),
     *          THE3OP(NDTHET),THE3ON(NDTHET)
c
      Dimension DSO0XP(NDTHET),DSO0XN(NDTHET),
     *          DSO0XRP(NDTHET),DSO0XRN(NDTHET),
     *          DSO0YP(NDTHET),DSO0YN(NDTHET),
     *          DSO0ZP(NDTHET),DSO0ZN(NDTHET),
     *          DSO1ZP(NDTHET),DSO1ZN(NDTHET),
     *          DSO3YP(NDTHET),DSO3YN(NDTHET)
c
      Dimension itz1d(nmax,nmax)
      Dimension ijx1d(nxmax,nxmax),ijy1d(nymax,nymax),ijz1d(nzmax,nzmax)
      Dimension ilsomu(nmax1,nmax1),idsomu(nmax,nmax)
      Dimension itet0(izmax,izmax),itet1(izmax,izmax),
     *          itet3(izmax,izmax)
      Dimension idsoxy0(izmax,izmax),idsoz0(izmax,izmax),
     *          idsoz1 (izmax,izmax),idsoy3(izmax,izmax)
c
      Dimension gp1p(nrop),gp2p(nrop),gn1p(nrop),gn2p(nrop)
      Dimension gp1m(nrom),gp2m(nrom),gn1m(nrom),gn2m(nrom)
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      Common/DIMPDD/ irortmp,idmus0,jdmus0,idmus1,jdmus1,irorz,irory
c
      Common/DIMEN/kmax,kmax1,kxmax,kymax,kzmax,nwf2,
     * maxtz1,maxjz1,maxjx1,maxjy1,maxlso,maxdso,
     * ndmu,ndthet,nacou,nrop,nrom,nrop8,nrom8,nherm,nherm38,nlegn
c
      Common/GOGINT/amu(2),xw(2),xh(2),xb(2),xm(2),WLS,t3,alpha,x0,e2
      Common/GOGHAN/w2b(2),h2m(2),x01,x02
      Common/OSCLEN/bx,by,bz
      Common /Const/ pi,sqpi,dln2,sq2
      Common /PULGA/ ideb
      Common /HFBOPT/ Acom,vcom,vcom2,icouech
c
c
c    Previous calculations are performed in predir BB, Spin-Orbit
c    and Coulomb direct.
c
      call timeit(0,11,'  PREDIRSO      ')
      call predirso(tz1d,ajx1d,ajy1d,ajz1d,itz1d,ijx1d,ijy1d,ijz1d,
     * alsomu,ilsomu,dsomu,idsomu,sxyz1p,sxyz1n,sxyz2p,sxyz2n,sxyz3p,
     * sxyz3n,sxyz4p,sxyz4n,sxyzc,Soxyzp,Soxyzn,GROP,GROM)
      call timeit(1,11,'  PREDIRSO      ')
c
c   Some additional folding is performed in predirco for the direct
c   Coulomb term     sxyzc --> scf
c
      call timeit(0,12,'  PREDIRCO      ')
      call predirco(sxyzc,scf,coum,acou)
      call timeit(1,12,'  PREDIRCO      ')
c
c    Summing the different contributions (the ones to be folded with T )
c    The sum goes to
c
c        GSUM0E   GSUM0O           Field (0)       Even(E) and Odd(O)
c                 GSUM1O           Field (1)       under time reversal
c                 GSUM3O           Field (3)
c
      call timeit(0,14,'REST OF GDIR    ')
c
      onethird = 1.0d+00/3.0d+00
      pi34  = pi**(0.75d+00)
      pi32  = pi**(1.50d+00)
      fe2   = e2
      fd1   = 8.0d+00*t3/pi34
      
      if(icouech.ne.1) then
         fec = 0.0d+00
      else
         fec   = 8.0d+00*e2/pi34*(3.0d+00/pi)**onethird   ! Coulomb exchange
      end if
      
      fdr   = 0.5d+00*fd1*alpha
      fsox  = WLS/(8.0d+00*pi32*bx*by**2*bz**2)
      fsoy  = WLS/(8.0d+00*pi32*bx**2*by*bz**2)
      fsoz  = WLS/(8.0d+00*pi32*bx**2*by**2*bz)
c
c                                                     Field (0)
c
      mumax = 2*nxmax-1
      mumay = 2*nymax-1
      mumaz = 2*nzmax-1
      i = 0
      do mux =1,mumax,2
         do muy=1,mumay,2
            do muz=1,mumaz,2
               i = i + 1
               s1t1 = sxyz1p(1,i)+sxyz1n(1,i)
               s1t2 = sxyz1p(2,i)+sxyz1n(2,i)
               s2t1 = sxyz2p(1,i)+sxyz2n(1,i)
               s2t2 = sxyz2p(2,i)+sxyz2n(2,i)
               dmu1t= dmu(1,i)+dmu(2,i)
               dmu3t= dmu(3,i)+dmu(4,i)
               dmu10t = dmu(9,i)+dmu(11,i)
               dmu13t = dmu(12,i)+dmu(14,i)
               sxt = Soxyzp(8,i)+ Soxyzn(8,i)
               syt = Soxyzp(9,i)+ Soxyzn(9,i)
               szt = Soxyzp(10,i)+ Soxyzn(10,i)
               dxt = Soxyzp(5,i)+ Soxyzn(5,i)
c                                                   Proton
               GSUM0EP(i) =
     * 0.5d+00*(w2b(1)*s1t1+w2b(2)*s1t2
     *          -h2m(1)*sxyz1p(1,i)-h2m(2)*sxyz1p(2,i))      ! BB
     * + fe2*scf(i)                                          ! Coulomb D
     * + fd1*(x01*dmu1t-x02*dmu(1,i))                        ! DD
     * - fec*dmu(15,i)                                       ! Coulomb E
     * + fdr*(dmu10t*(x01-x02)+2.0d+00*x01*dmu(10,i) +
     * 0.5d+00*(x0-1)*dmu13t + x0*dmu(13,i) )                ! Rearr
     * + (sxt+Soxyzp(8,i))*fsox+
     * (syt+Soxyzp(9,i))*fsoy+(szt+Soxyzp(10,i))*fsoz       ! Spin-Orbit
c                                                   Neutron
               GSUM0EN(i) =
     * 0.5d+00*(w2b(1)*s1t1+w2b(2)*s1t2
     *          -h2m(1)*sxyz1n(1,i)-h2m(2)*sxyz1n(2,i))      ! BB
     * + fd1*(x01*dmu1t-x02*dmu(2,i))                        ! DD
     * + fdr*(dmu10t*(x01-x02)+2.0d+00*x01*dmu(10,i) +
     * 0.5d+00*(x0-1)*dmu13t + x0*dmu(13,i) )                ! Rearr
     * + (sxt+Soxyzn(8,i))*fsox+
     * (syt+Soxyzn(9,i))*fsoy+(szt+Soxyzn(10,i))*fsoz        ! Spin-Orbit
c                                                   Proton
               GSUM0OP(i) =
     * 0.5d+00*(xb(1)*s2t1+xb(2)*s2t2
     *          -xm(1)*sxyz2p(1,i)-xm(2)*sxyz2p(2,i))        ! BB
     * + fd1*0.5* (x0*dmu3t-dmu(3,i))                        ! DD
     * -(dxt+Soxyzp(5,i))*fsox                               !Spin-orbit
c                                                   Neutron
               GSUM0ON(i) =
     * 0.5d+00*(xb(1)*s2t1+xb(2)*s2t2
     *          -xm(1)*sxyz2n(1,i)-xm(2)*sxyz2n(2,i))        ! BB
     * + fd1*0.5* (x0*dmu3t-dmu(4,i))                        ! DD
     * - (dxt+Soxyzn(5,i))*fsox                 ! Spin-orbit
c
c                                                For the J part of SO
c
               gxt = Soxyzp(1,i) + Soxyzn(1,i)
               rxt = Soxyzp(4,i) + Soxyzn(4,i)
               Soxyzp(1,i) = gxt + Soxyzp(1,i)
               Soxyzn(1,i) = gxt + Soxyzn(1,i)
               Soxyzp(4,i) = rxt + Soxyzp(4,i)
               Soxyzn(4,i) = rxt + Soxyzn(4,i)
            end do                               ! muz
         end do                                  ! muy
      end do                                     ! mux
c
c
c                                                     Field (1)
c
      i = 0
      do mux =2,mumax,2
         do muy=1,mumay,2
            do muz=2,mumaz,2
               i = i + 1
               s3t1 = sxyz3p(1,i)+sxyz3n(1,i)
               s3t2 = sxyz3p(2,i)+sxyz3n(2,i)
               dmu5t= dmu(5,i)+dmu(6,i)
               dzt = Soxyzp(7,i)+ Soxyzn(7,i)
c                                                   Proton
               GSUM1OP(i) =
     * 0.5d+00*(xb(1)*s3t1+xb(2)*s3t2
     *          -xm(1)*sxyz3p(1,i)-xm(2)*sxyz3p(2,i))        ! BB
     * + fd1*0.5* (x0*dmu5t-dmu(5,i))                        ! DD
     * + (dzt+Soxyzp(7,i))*fsoz                 ! Spin-orbit
c                                                   Neutron
               GSUM1ON(i) =
     * 0.5d+00*(xb(1)*s3t1+xb(2)*s3t2
     *          -xm(1)*sxyz3n(1,i)-xm(2)*sxyz3n(2,i))        ! BB
     * + fd1*0.5* (x0*dmu5t-dmu(6,i))                        ! DD
     * + (dzt+Soxyzn(7,i))*fsoz                 ! Spin-orbit
c
c                                                For the J part of SO
c
               gzt = Soxyzp(3,i) + Soxyzn(3,i)
               Soxyzp(3,i) = gzt + Soxyzp(3,i)
               Soxyzn(3,i) = gzt + Soxyzn(3,i)
            end do                               ! muz
         end do                                  ! muy
      end do                                     ! mux
c
c                                                     Field (3)
c
      i = 0
      do mux =2,mumax,2
         do muy=2,mumay,2
            do muz=1,mumaz,2
               i = i + 1
               s4t1 = sxyz4p(1,i)+sxyz4n(1,i)
               s4t2 = sxyz4p(2,i)+sxyz4n(2,i)
               dmu7t= dmu(7,i)+dmu(8,i)
               dyt = Soxyzp(6,i)+ Soxyzn(6,i)
c                                                   Proton
               GSUM3OP(i) =
     * 0.5d+00*(xb(1)*s4t1+xb(2)*s4t2
     *          -xm(1)*sxyz4p(1,i)-xm(2)*sxyz4p(2,i))        ! BB
     * + fd1*0.5* (x0*dmu7t-dmu(7,i))                        ! DD
     * - (dyt+Soxyzp(6,i))*fsoy                 ! Spin-orbit
c                                                   Neutron
               GSUM3ON(i) =
     * 0.5d+00*(xb(1)*s4t1+xb(2)*s4t2
     *          -xm(1)*sxyz4n(1,i)-xm(2)*sxyz4n(2,i))        ! BB
     * + fd1*0.5* (x0*dmu7t-dmu(8,i))                        ! DD
     * - (dyt+Soxyzn(6,i))*fsoy                  ! Spin-orbit
c
c                                                For the J part of SO
c
               gyt = Soxyzp(2,i) + Soxyzn(2,i)
               Soxyzp(2,i) = gyt + Soxyzp(2,i)
               Soxyzn(2,i) = gyt + Soxyzn(2,i)
            end do                               ! muz
         end do                                  ! muy
      end do                                     ! mux
c
c     Optimization  (THETA)
c
c                                    FIELD (0)
         imumix = 1
         imumiy = 1
         imumiz = 1
c
      imumax = 2*nxmax-1
      imumay = 2*nymax-1
      imumaz = 2*nzmax-1
      itet = 0
      do nz1=1,nzmax
         izmin = 1+MOD((nz1+imumiz),2)
         do nz2=izmin,nzmax,2
            itz = itz1d(nz1,nz2)
            muzic= Iabs(nz2-nz1) + 1
            muzfc= nz2+nz1-1
            inds = 0
            itet0(nz1,nz2) = itet
            do imux=imumix,imumax,2
               do imuy=imumiy,imumay,2
                  itet = itet + 1
                  the0ep(itet) = 0.0d+00
                  the0op(itet) = 0.0d+00
                  the0en(itet) = 0.0d+00
                  the0on(itet) = 0.0d+00
                  sumep     = 0.0d+00
                  sumop     = 0.0d+00
                  sumen     = 0.0d+00
                  sumon     = 0.0d+00
                  indc = itz
                  do imuz=imumiz,imumaz,2
                     inds = inds + 1
                     if(imuz.ge.muzic.and.imuz.le.muzfc) then
                        indc = indc + 1
                        sumep = sumep + tz1d(indc)*gsum0ep(inds)
                        sumop = sumop + tz1d(indc)*gsum0op(inds)
                        sumen = sumen + tz1d(indc)*gsum0en(inds)
                        sumon = sumon + tz1d(indc)*gsum0on(inds)
                        end if
                  end do                               ! imuz
                  the0ep (itet) = sumep
                  the0op (itet) = sumop
                  the0en (itet) = sumen
                  the0on (itet) = sumon
               end do                                  ! imuy
            end do                                     ! imux
         end do                                        ! nz2
      end do                                           ! nz1
c                                   FIELD (1)
      if(ideb.gt.0) write(6,*) ' GDIR after opt Field (0) itet ',itet
         imumix = 2
         imumiy = 1
         imumiz = 2
c
      itet = 0
      do nz1=1,nzmax
         izmin = 1+MOD((nz1+imumiz),2)
         do nz2=izmin,nzmax,2
            itz = itz1d(nz1,nz2)
            muzic= Iabs(nz2-nz1) + 1
            muzfc= nz2+nz1-1
            inds = 0
            itet1(nz1,nz2) = itet
            do imux=imumix,imumax,2
               do imuy=imumiy,imumay,2
                  itet = itet + 1
                  the1op(itet) = 0.0d+00
                  the1on(itet) = 0.0d+00
                  sumop     = 0.0d+00
                  sumon     = 0.0d+00
                  indc = itz
                  do imuz=imumiz,imumaz,2
                     inds = inds + 1
                     if(imuz.ge.muzic.and.imuz.le.muzfc) then
                        indc = indc + 1
                        sumop = sumop + tz1d(indc)*gsum1op(inds)
                        sumon = sumon + tz1d(indc)*gsum1on(inds)
                        end if
                  end do                               ! imuz
                  the1op (itet) = sumop
                  the1on (itet) = sumon
               end do                                  ! imuy
            end do                                     ! imux
         end do                                        ! nz2
      end do                                           ! nz1
c                                      FIELD (3)
      if(ideb.gt.0) write(6,*) ' GDIR after opt Field (1) itet ',itet
         imumix = 2
         imumiy = 2
         imumiz = 1
c
      itet = 0
      do nz1=1,nzmax
         izmin = 1+MOD((nz1+imumiz),2)
         do nz2=izmin,nzmax,2
            itz = itz1d(nz1,nz2)
            muzic= Iabs(nz2-nz1) + 1
            muzfc= nz2+nz1-1
            inds = 0
            itet3(nz1,nz2) = itet
            do imux=imumix,imumax,2
               do imuy=imumiy,imumay,2
                  itet = itet + 1
                  the3op(itet) = 0.0d+00
                  the3on(itet) = 0.0d+00
                  sumop     = 0.0d+00
                  sumon     = 0.0d+00
                  indc = itz
                  do imuz=imumiz,imumaz,2
                     inds = inds + 1
                     if(imuz.ge.muzic.and.imuz.le.muzfc) then
                        indc = indc + 1
                        sumop = sumop + tz1d(indc)*gsum3op(inds)
                        sumon = sumon + tz1d(indc)*gsum3on(inds)
                        end if
                  end do                               ! imuz
                  the3op (itet) = sumop
                  the3on (itet) = sumon
               end do                                  ! imuy
            end do                                     ! imux
         end do                                        ! nz2
      end do                                           ! nz1
c
c
c     Optimization  Spin orbit J part
c
c                                    FIELD (0)
c                               gamma_x field (2) and ro field (1)
         imumix = 1
         imumiy = 1
         imumiz = 1
c
      itet = 0
      do nz1=1,nzmax
         izmin = 1+MOD(nz1,2)    ! D(nz1|nz2,mzu)
         do nz2=izmin,nzmax,2
            idso= idsomu(nz1,nz2)
            inds = 0
            idsoxy0(nz1,nz2) = itet
            do imux=imumix,imumax,2
               do imuy=imumiy,imumay,2
                  itet = itet + 1
                  dso0xp (itet) = 0.0d+00
                  dso0xn (itet) = 0.0d+00
                  dso0xRp(itet) = 0.0d+00
                  dso0xRn(itet) = 0.0d+00
                  dso0yp (itet) = 0.0d+00
                  dso0yn (itet) = 0.0d+00
                  sumxp      = 0.0d+00
                  sumxRp     = 0.0d+00
                  sumyp      = 0.0d+00
                  sumxn      = 0.0d+00
                  sumxRn     = 0.0d+00
                  sumyn      = 0.0d+00
                  indc = idso
                  do imuz=imumiz,imumaz,2
                     inds = inds + 1
                        sumxp = sumxp + dsomu(indc)*Soxyzp(1,inds)
                        sumxRp= sumxRp+ dsomu(indc)*Soxyzp(4,inds)
                        sumyp = sumyp + dsomu(indc)*Soxyzp(4,inds)
                        sumxn = sumxn + dsomu(indc)*Soxyzn(1,inds)
                        sumxRn= sumxRn+ dsomu(indc)*Soxyzn(4,inds)
                        sumyn = sumyn + dsomu(indc)*Soxyzn(4,inds)
                        indc = indc + 1
                  end do                               ! imuz
                  dso0xp (itet) = sumxp *fsox
                  dso0xrp(itet) = sumxrp*fsox
                  dso0yp (itet) = sumyp *fsoy
                  dso0xn (itet) = sumxn *fsox
                  dso0xrn(itet) = sumxrn*fsox
                  dso0yn (itet) = sumyn *fsoy
               end do                                  ! imuy
            end do                                     ! imux
         end do                                        ! nz2
      end do                                           ! nz1
c
c                                    FIELD (0)
c                                ro field (3)
         imumix = 1
         imumiy = 1
         imumiz = 1
c
      itet = 0
      do nz1=1,nzmax
         izmin = 1+MOD((nz1+1),2)    ! L(nz1,nz2,mzu)
         do nz2=izmin,nzmax,2
            ilso= ilsomu(nz1,nz2)
            inds = 0
            idsoz0(nz1,nz2) = itet
            do imux=imumix,imumax,2
               do imuy=imumiy,imumay,2
                  itet = itet + 1
                  dso0zp(itet) = 0.0d+00
                  dso0zn(itet) = 0.0d+00
                  sumzp     = 0.0d+00
                  sumzn     = 0.0d+00
                  indc = ilso
                  do imuz=imumiz,imumaz,2
                     inds = inds + 1
                        sumzp = sumzp +alsomu(indc)*Soxyzp(4,inds)
                        sumzn = sumzn +alsomu(indc)*Soxyzn(4,inds)
                        indc = indc + 1
                  end do                               ! imuz
                  dso0zp (itet) = sumzp*fsoz
                  dso0zn (itet) = sumzn*fsoz
               end do                                  ! imuy
            end do                                     ! imux
         end do                                        ! nz2
      end do                                           ! nz1
c
c                                   FIELD (1)
c                                   gamma_z
         imumix = 2
         imumiy = 1
         imumiz = 2
c
      itet = 0
      do nz1=1,nzmax
         izmin = 1+MOD(nz1,2)    ! D(nz1|nz2,mzu)
         do nz2=izmin,nzmax,2
            ilso= ilsomu(nz1,nz2)
            inds = 0
            idsoz1(nz1,nz2) = itet
            do imux=imumix,imumax,2
               do imuy=imumiy,imumay,2
                  itet = itet + 1
                  dso1zp(itet) = 0.0d+00
                  dso1zn(itet) = 0.0d+00
                  sumzp     = 0.0d+00
                  sumzn     = 0.0d+00
                  indc = ilso
                  do imuz=imumiz,imumaz,2
                     inds = inds + 1
                        sumzp = sumzp + alsomu(indc)*Soxyzp(3,inds)
                        sumzn = sumzn + alsomu(indc)*Soxyzn(3,inds)
                        indc = indc + 1
                  end do                               ! imuz
                  dso1zp (itet) = sumzp*fsoz
                  dso1zn (itet) = sumzn*fsoz
               end do                                  ! imuy
            end do                                     ! imux
         end do                                        ! nz2
      end do                                           ! nz1
c
c                                      FIELD (3)
c                                      gamma_y
         imumix = 2
         imumiy = 2
         imumiz = 1
c
      itet = 0
      do nz1=1,nzmax
         izmin = 1+MOD(nz1,2)    ! D(nz1|nz2,mzu)
         do nz2=izmin,nzmax,2
            idso= idsomu(nz1,nz2)
            inds = 0
            idsoy3(nz1,nz2) = itet
            do imux=imumix,imumax,2
               do imuy=imumiy,imumay,2
                  itet = itet + 1
                  dso3yp(itet) = 0.0d+00
                  dso3yn(itet) = 0.0d+00
                  sumyp     = 0.0d+00
                  sumyn     = 0.0d+00
                  indc = idso
                  do imuz=imumiz,imumaz,2
                     inds = inds + 1
                        sumyp = sumyp + dsomu(indc)*Soxyzp(2,inds)
                        sumyn = sumyn + dsomu(indc)*Soxyzn(2,inds)
                        indc = indc + 1
                  end do                               ! imuz
                  dso3yp (itet) = sumyp*fsoy
                  dso3yn (itet) = sumyn*fsoy
               end do                                  ! imuy
            end do                                     ! imux
         end do                                        ! nz2
      end do                                           ! nz1
c    +---------------------------------------------------------+
c    |   Starts calculation of the direct field                |
c    +---------------------------------------------------------+
      irop = 0
      irom = 0
      do 1 nx1=1,nxmax
        nym1 = MY(nx1)
        fasx = 1-2*Mod((nx1+1),2)
        do 11 nx2=1,nx1
          icx  = Mod((nx1+nx2),2)                     ! type of field
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          itx = itz1d (nx1,nx2)                        ! tz1d  index
          idx = idsomu(nx1,nx2)                        ! dsomu index
          idxr= idsomu(nx2,nx1)                        ! dsomu index
          ilx = ilsomu(nx1,nx2)                        ! dsomu index
          muxic= Iabs(nx2-nx1) + 1
          muxfc= nx2+nx1-1
          do 2 ny1 = 1,nym1
            nxy1 = nx1+ny1
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do 12 ny2=1,nym2
              icy  = Mod((ny1+ny2),2)                 ! type of field
              nxy2 = nx2+ny2
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              ity = itz1d(ny1,ny2)                       ! tz1d index
              idy = idsomu(ny1,ny2)                      ! dsomu index
              idyr= idsomu(ny2,ny1)                      ! dsomu index
              ily = ilsomu(ny1,ny2)                      ! dsomu index
              muyic= Iabs(ny2-ny1) + 1
              muyfc= ny2+ny1-1
              if(icy.eq.0) then
                 fasy = 1-2*Mod(Iabs(ny2-ny1)/2,2)
              else
                 fasy = 1-2*Mod(Iabs(ny2-ny1+1)/2,2)
              end if
              fasxy = fasx*fasy
              lxy12 = lx12.and.(ny1.eq.ny2)
c
c---->    Positive parity
c
c
c be careful not to move the if below, the iz1e=iy1e need to be done
c
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
c
              do 30 nz1 =nzi1e,nzm1,2
                if(lxy12) nzm2=nz1
                do 31 nz2 = nzi2e,nzm2,2
                  irop = irop + 1
c
            if((icx.eq.0).and.(icy.eq.0)) then       ! FIELD (0)
            imumix=1
            imumiy=1
            ithet0 = itet0(nz1,nz2)
            indcx = itx
            sumep = 0.0d+00
            sumen = 0.0d+00
            sumop = 0.0d+00
            sumon = 0.0d+00
            do imux=imumix,imumax,2
               vmux =(imux.ge.muxic.and.imux.le.muxfc)
               if(vmux) indcx = indcx + 1
               tx = tz1d(indcx)
               indcy = ity
               do imuy=imumiy,imumay,2
                  ithet0 = ithet0 + 1
               if(vmux.and.(imuy.ge.muyic.and.imuy.le.muyfc)) then
               indcy = indcy + 1
               ty = tz1d(indcy)
               sumep = sumep + the0ep(ithet0)*tx*ty
               sumen = sumen + the0en(ithet0)*tx*ty
               sumop = sumop + the0op(ithet0)*tx*ty
               sumon = sumon + the0on(ithet0)*tx*ty
               end if
               end do                                  ! imuy
            end do                                     ! imux
c
           gp1p(irop) = gp1p(irop) + (sumep + fasx*sumop)*fasy
           gp2p(irop) = gp2p(irop) + (sumep - fasx*sumop)*fasy
           gn1p(irop) = gn1p(irop) + (sumen + fasx*sumon)*fasy
           gn2p(irop) = gn2p(irop) + (sumen - fasx*sumon)*fasy
c
            else if((icx.eq.1).and.(icy.eq.0)) then       ! FIELD (1)
            imumix=2                                      ! standar
            imumiy=1
            ithet1 = itet1(nz1,nz2)
            indcx = itx
            sumep = 0.0d+00
            sumen = 0.0d+00
            sumop = 0.0d+00
            sumon = 0.0d+00
            do imux=imumix,imumax,2
               vmux =(imux.ge.muxic.and.imux.le.muxfc)
               if(vmux) indcx = indcx + 1
               tx = tz1d(indcx)
               indcy = ity
               do imuy=imumiy,imumay,2
                  ithet1 = ithet1 + 1
               if(vmux.and.(imuy.ge.muyic.and.imuy.le.muyfc)) then
               indcy = indcy + 1
               ty = tz1d(indcy)
               sumop = sumop + the1op(ithet1)*tx*ty
               sumon = sumon + the1on(ithet1)*tx*ty
               end if
               end do                                  ! imuy
            end do                                     ! imux
c
            imumix=1                                   ! spin orbit J
            imumiy=1
            ithet = idsoxy0(nz1,nz2)
            ithetr= idsoxy0(nz2,nz1)
            indcx = idx
            indcxr= idxr
            susop = 0.0d+00
            suson = 0.0d+00
            do imux=imumix,imumax,2
               dx  = dsomu(indcx )
               dxr = dsomu(indcxr)
               indcy = ily
               do imuy=imumiy,imumay,2
                  ithet = ithet + 1
                  ithetr= ithetr+ 1
               aly = alsomu(indcy)
               susop = susop + aly*(dso0yp(ithet )*dxr
     *                            - dso0yp(ithetr)*dx  )
               suson = suson + aly*(dso0yn(ithet )*dxr
     *                            - dso0yn(ithetr)*dx  )
               indcy = indcy + 1
               end do                                  ! imuy
               indcx = indcx + 1
               indcxr= indcxr+ 1
            end do                                     ! imux
c
           gp1p(irop) = gp1p(irop) + (susop*fasx + sumop) * fasy
           gp2p(irop) = gp2p(irop) + (susop*fasx - sumop) * fasy
           gn1p(irop) = gn1p(irop) + (suson*fasx + sumon) * fasy
           gn2p(irop) = gn2p(irop) + (suson*fasx - sumon) * fasy
c
c
            else if((icx.eq.0).and.(icy.eq.1)) then       ! FIELD (2)
c
            imumix=1                                   ! spin orbit J
            imumiy=1
            ithet = idsoxy0(nz1,nz2)
            ithetr= idsoxy0(nz2,nz1)
            indcx = ilx
            susope = 0.0d+00
            susone = 0.0d+00
            susopo = 0.0d+00
            susono = 0.0d+00
            do imux=imumix,imumax,2
               alx  = alsomu(indcx )
               indcy = idy
               indcyr= idyr
               do imuy=imumiy,imumay,2
                  ithet = ithet + 1
                  ithetr= ithetr+ 1
               dy = dsomu(indcy )
               dyr= dsomu(indcyr)
               susope = susope + alx*(dso0xrp(ithetr)*dy
     *                            - dso0xrp(ithet )*dyr  )
               susone = susone + alx*(dso0xrn(ithetr)*dy
     *                            - dso0xrn(ithet )*dyr  )
               susopo = susopo + alx*(dso0xp(ithetr)*dy
     *                            - dso0xp(ithet )*dyr  )
               susono = susono + alx*(dso0xn(ithetr)*dy
     *                            - dso0xn(ithet )*dyr  )
               indcy = indcy + 1
               indcyr= indcyr+ 1
               end do                                  ! imuy
               indcx = indcx + 1
            end do                                     ! imux
c
c
            imumix=2                                   ! spin orbit J
            imumiy=1
            ithet = idsoz1(nz1,nz2)
            indcx = idx
            indcxr= idxr
            do imux=imumix,imumax,2
               dx  = dsomu(indcx )
               dxr = dsomu(indcxr)
               indcy = idy
               indcyr= idyr
               do imuy=imumiy,imumay,2
               ithet = ithet + 1
               dy = dsomu(indcy )
               dyr= dsomu(indcyr)
               susopo = susopo + dso1zp(ithet)*(dx*dyr-dxr*dy)
               susono = susono + dso1zn(ithet)*(dx*dyr-dxr*dy)
               indcy = indcy + 1
               indcyr= indcyr+ 1
               end do                                  ! imuy
               indcx  = indcx  + 1
               indcxr = indcxr + 1
            end do                                     ! imux
c
c
            imumix=2                                   ! spin orbit J
            imumiy=2
            ithet = idsoy3(nz1,nz2)
            ithetr= idsoy3(nz2,nz1)
            indcx = idx
            indcxr= idxr
            do imux=imumix,imumax,2
               dx  =  dsomu(indcx )
               dxr =  dsomu(indcxr)
               indcy = ily
               do imuy=imumiy,imumay,2
               ithet = ithet + 1
               ithetr= ithetr+ 1
               aly = alsomu(indcy )
               susopo = susopo + aly*(dso3yp(ithet )*dxr
     *                            - dso3yp(ithetr)*dx   )
               susono = susono + aly*(dso3yn(ithet )*dxr
     *                            - dso3yn(ithetr)*dx   )
               indcy = indcy + 1
               end do                                  ! imuy
               indcx = indcx + 1
               indcxr= indcxr+ 1
            end do                                     ! imux
c
           gp1p(irop) = gp1p(irop) + (-susope*fasx + susopo ) * fasy
           gp2p(irop) = gp2p(irop) + (-susope*fasx - susopo ) * fasy
           gn1p(irop) = gn1p(irop) + (-susone*fasx + susono ) * fasy
           gn2p(irop) = gn2p(irop) + (-susone*fasx - susono ) * fasy
c
            else if((icx.eq.1).and.(icy.eq.1)) then       ! FIELD (3)
            imumix=2                               ! standar
            imumiy=2
            ithet3 = itet3(nz1,nz2)
            indcx = itx
            sumep = 0.0d+00
            sumen = 0.0d+00
            sumop = 0.0d+00
            sumon = 0.0d+00
            do imux=imumix,imumax,2
               vmux =(imux.ge.muxic.and.imux.le.muxfc)
               if(vmux) indcx = indcx + 1
               tx = tz1d(indcx)
               indcy = ity
               do imuy=imumiy,imumay,2
                  ithet3 = ithet3 + 1
               if(vmux.and.(imuy.ge.muyic.and.imuy.le.muyfc)) then
               indcy = indcy + 1
               ty = tz1d(indcy)
               sumop = sumop + the3op(ithet3)*tx*ty
               sumon = sumon + the3on(ithet3)*tx*ty
               end if
               end do                                  ! imuy
            end do                                     ! imux
c
            imumix=1                               ! Spin orbit J
            imumiy=1
            ithet = idsoz0(nz1,nz2)
            indcx = idx
            indcxr= idxr
            susop = 0.0d+00
            suson = 0.0d+00
            do imux=imumix,imumax,2
               dx  = dsomu(indcx )
               dxr = dsomu(indcxr)
               indcy = idy
               indcyr= idyr
               do imuy=imumiy,imumay,2
               ithet = ithet + 1
               dy  = dsomu(indcy )
               dyr = dsomu(indcyr)
               susop = susop + dso0zp(ithet)*(dx*dyr-dxr*dy)
               suson = suson + dso0zn(ithet)*(dx*dyr-dxr*dy)
               indcy = indcy + 1
               indcyr= indcyr+ 1
               end do                                  ! imuy
               indcx = indcx + 1
               indcxr= indcxr+ 1
            end do                                     ! imux
c
           gp1p(irop) = gp1p(irop) + (susop*fasx + sumop) * fasxy
           gp2p(irop) = gp2p(irop) + (susop*fasx - sumop) * fasxy
           gn1p(irop) = gn1p(irop) + (suson*fasx + sumon) * fasxy
           gn2p(irop) = gn2p(irop) + (suson*fasx - sumon) * fasxy
c
           end if
c
   31           continue
   30         continue
c
c                                           ---->    Negative parity
c
              end if
c
c be careful not to move the if below, the iz1o=iy1o need to be done
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
c
              do 40 nz1 =nzi1o,nzm1,2
                if(lxy12) nzm2=nz1
                do 41 nz2 = nzi2o,nzm2,2
                  irom = irom + 1
c
            if((icx.eq.0).and.(icy.eq.0)) then       ! FIELD (0)
            imumix=1
            imumiy=1
            ithet0 = itet0(nz1,nz2)
            indcx = itx
            sumep = 0.0d+00
            sumen = 0.0d+00
            sumop = 0.0d+00
            sumon = 0.0d+00
            do imux=imumix,imumax,2
               vmux =(imux.ge.muxic.and.imux.le.muxfc)
               if(vmux) indcx = indcx + 1
               tx = tz1d(indcx)
               indcy = ity
               do imuy=imumiy,imumay,2
                  ithet0 = ithet0 + 1
               if(vmux.and.(imuy.ge.muyic.and.imuy.le.muyfc)) then
               indcy = indcy + 1
               ty = tz1d(indcy)
               sumep = sumep + the0ep(ithet0)*tx*ty
               sumen = sumen + the0en(ithet0)*tx*ty
               sumop = sumop + the0op(ithet0)*tx*ty
               sumon = sumon + the0on(ithet0)*tx*ty
               end if
               end do                                  ! imuy
            end do                                     ! imux
c
           gp1m(irom) = gp1m(irom) + (sumep + fasx*sumop)*fasy
           gp2m(irom) = gp2m(irom) + (sumep - fasx*sumop)*fasy
           gn1m(irom) = gn1m(irom) + (sumen + fasx*sumon)*fasy
           gn2m(irom) = gn2m(irom) + (sumen - fasx*sumon)*fasy
c
            else if((icx.eq.1).and.(icy.eq.0)) then       ! FIELD (1)
            imumix=2
            imumiy=1
            ithet1 = itet1(nz1,nz2)
            indcx = itx
            sumep = 0.0d+00
            sumen = 0.0d+00
            sumop = 0.0d+00
            sumon = 0.0d+00
            do imux=imumix,imumax,2
               vmux =(imux.ge.muxic.and.imux.le.muxfc)
               if(vmux) indcx = indcx + 1
               tx = tz1d(indcx)
               indcy = ity
               do imuy=imumiy,imumay,2
                  ithet1 = ithet1 + 1
               if(vmux.and.(imuy.ge.muyic.and.imuy.le.muyfc)) then
               indcy = indcy + 1
               ty = tz1d(indcy)
               sumop = sumop + the1op(ithet1)*tx*ty
               sumon = sumon + the1on(ithet1)*tx*ty
               end if
               end do                                  ! imuy
            end do                                     ! imux
c
            imumix=1                                   ! spin orbit J
            imumiy=1
            ithet = idsoxy0(nz1,nz2)
            ithetr= idsoxy0(nz2,nz1)
            indcx = idx
            indcxr= idxr
            susop = 0.0d+00
            suson = 0.0d+00
            do imux=imumix,imumax,2
               dx  = dsomu(indcx )
               dxr = dsomu(indcxr)
               indcy = ily
               do imuy=imumiy,imumay,2
                  ithet = ithet + 1
                  ithetr= ithetr+ 1
               aly = alsomu(indcy)
               susop = susop + aly*(dso0yp(ithet )*dxr
     *                            - dso0yp(ithetr)*dx  )
               suson = suson + aly*(dso0yn(ithet )*dxr
     *                            - dso0yn(ithetr)*dx  )
               indcy = indcy + 1
               end do                                  ! imuy
               indcx = indcx + 1
               indcxr= indcxr+ 1
            end do                                     ! imux
c
           gp1m(irom) = gp1m(irom) + (susop*fasx + sumop) * fasy
           gp2m(irom) = gp2m(irom) + (susop*fasx - sumop) * fasy
           gn1m(irom) = gn1m(irom) + (suson*fasx + sumon) * fasy
           gn2m(irom) = gn2m(irom) + (suson*fasx - sumon) * fasy
c
c
            else if((icx.eq.0).and.(icy.eq.1)) then       ! FIELD (2)
c
            imumix=1                                   ! spin orbit J
            imumiy=1
            ithet = idsoxy0(nz1,nz2)
            ithetr= idsoxy0(nz2,nz1)
            indcx = ilx
            susope = 0.0d+00
            susone = 0.0d+00
            susopo = 0.0d+00
            susono = 0.0d+00
            do imux=imumix,imumax,2
               alx  = alsomu(indcx )
               indcy = idy
               indcyr= idyr
               do imuy=imumiy,imumay,2
                  ithet = ithet + 1
                  ithetr= ithetr+ 1
               dy = dsomu(indcy )
               dyr= dsomu(indcyr)
               susope = susope + alx*(dso0xrp(ithetr)*dy
     *                            - dso0xrp(ithet )*dyr  )
               susone = susone + alx*(dso0xrn(ithetr)*dy
     *                            - dso0xrn(ithet )*dyr  )
               susopo = susopo + alx*(dso0xp(ithetr)*dy
     *                            - dso0xp(ithet )*dyr  )
               susono = susono + alx*(dso0xn(ithetr)*dy
     *                            - dso0xn(ithet )*dyr  )
               indcy = indcy + 1
               indcyr= indcyr+ 1
               end do                                  ! imuy
               indcx = indcx + 1
            end do                                     ! imux
c
c
            imumix=2                                   ! spin orbit J
            imumiy=1
            ithet = idsoz1(nz1,nz2)
            indcx = idx
            indcxr= idxr
            do imux=imumix,imumax,2
               dx  = dsomu(indcx )
               dxr = dsomu(indcxr)
               indcy = idy
               indcyr= idyr
               do imuy=imumiy,imumay,2
               ithet = ithet + 1
               dy = dsomu(indcy )
               dyr= dsomu(indcyr)
               susopo = susopo + dso1zp(ithet)*(dx*dyr-dxr*dy)
               susono = susono + dso1zn(ithet)*(dx*dyr-dxr*dy)
               indcy = indcy + 1
               indcyr= indcyr+ 1
               end do                                  ! imuy
               indcx  = indcx  + 1
               indcxr = indcxr + 1
            end do                                     ! imux
c
c
            imumix=2                                   ! spin orbit J
            imumiy=2
            ithet = idsoy3(nz1,nz2)
            ithetr= idsoy3(nz2,nz1)
            indcx = idx
            indcxr= idxr
            do imux=imumix,imumax,2
               dx  =  dsomu(indcx )
               dxr =  dsomu(indcxr)
               indcy = ily
               do imuy=imumiy,imumay,2
               ithet = ithet + 1
               ithetr= ithetr+ 1
               aly = alsomu(indcy )
               susopo = susopo + aly*(dso3yp(ithet )*dxr
     *                            - dso3yp(ithetr)*dx   )
               susono = susono + aly*(dso3yn(ithet )*dxr
     *                            - dso3yn(ithetr)*dx   )
               indcy = indcy + 1
               end do                                  ! imuy
               indcx = indcx + 1
               indcxr= indcxr+ 1
            end do                                     ! imux
c
           gp1m(irom) = gp1m(irom) + (-susope*fasx + susopo ) * fasy
           gp2m(irom) = gp2m(irom) + (-susope*fasx - susopo ) * fasy
           gn1m(irom) = gn1m(irom) + (-susone*fasx + susono ) * fasy
           gn2m(irom) = gn2m(irom) + (-susone*fasx - susono ) * fasy
c
            else if((icx.eq.1).and.(icy.eq.1)) then       ! FIELD (3)
            imumix=2
            imumiy=2
            ithet3 = itet3(nz1,nz2)
            indcx = itx
            sumep = 0.0d+00
            sumen = 0.0d+00
            sumop = 0.0d+00
            sumon = 0.0d+00
            do imux=imumix,imumax,2
               vmux =(imux.ge.muxic.and.imux.le.muxfc)
               if(vmux) indcx = indcx + 1
               tx = tz1d(indcx)
               indcy = ity
               do imuy=imumiy,imumay,2
                  ithet3 = ithet3 + 1
               if(vmux.and.(imuy.ge.muyic.and.imuy.le.muyfc)) then
               indcy = indcy + 1
               ty = tz1d(indcy)
               sumop = sumop + the3op(ithet3)*tx*ty
               sumon = sumon + the3on(ithet3)*tx*ty
               end if
               end do                                  ! imuy
            end do                                     ! imux
c
            imumix=1                               ! Spin orbit J
            imumiy=1
            ithet = idsoz0(nz1,nz2)
            indcx = idx
            indcxr= idxr
            susop = 0.0d+00
            suson = 0.0d+00
            do imux=imumix,imumax,2
               dx  = dsomu(indcx )
               dxr = dsomu(indcxr)
               indcy = idy
               indcyr= idyr
               do imuy=imumiy,imumay,2
               ithet = ithet + 1
               dy  = dsomu(indcy )
               dyr = dsomu(indcyr)
               susop = susop + dso0zp(ithet)*(dx*dyr-dxr*dy)
               suson = suson + dso0zn(ithet)*(dx*dyr-dxr*dy)
               indcy = indcy + 1
               indcyr= indcyr+ 1
               end do                                  ! imuy
               indcx = indcx + 1
               indcxr= indcxr+ 1
            end do                                     ! imux
c
           gp1m(irom) = gp1m(irom) + (susop*fasx + sumop) * fasxy
           gp2m(irom) = gp2m(irom) + (susop*fasx - sumop) * fasxy
           gn1m(irom) = gn1m(irom) + (suson*fasx + sumon) * fasxy
           gn2m(irom) = gn2m(irom) + (suson*fasx - sumon) * fasxy
c
           end if
   41           continue
   40         continue
            end if
   12       continue
    2     continue
   11   continue
    1 continue
      call timeit(1,14,'REST OF GDIR    ')
      return
      end
c+---------------------------------------------------------------------+
c|  Last change February 12 1993                                       |
c+=====================================================================+
c|  Quantities needed in the calculation of the direct field are       |
c|  computed for the following fields:                                 |
c|                                                                     |
c|   -  Brink-Boecker term                                             |
c|                                                                     |
c|   -  Coulomb interaction (Direct  + Exchange in the Slater approx)  |
c|                                                                     |
c|   -  Spin-Orbit                                                     |
c|                                                                     |
c+---------------------------------------------------------------------+
      Subroutine Predirso(tz1d,ajx1d,ajy1d,ajz1d,itz1d,ijx1d,ijy1d,ijz1d
     *,alsomu,ilsomu,dsomu,idsomu,sxyz1p,sxyz1n,sxyz2p,sxyz2n,sxyz3p,
     * sxyz3n,sxyz4p,sxyz4n,sxyzc,Soxyzp,Soxyzn,GROP,GROM)
c
      Implicit real*8(a-h,o-z)
      Implicit logical (l)
      Include 'COMDIM'
      Dimension TZ1D(maxtz1)
      Dimension AJZ1D(2,maxjz1),AJY1D(2,maxjy1),AJX1D(2,maxjx1)
      Dimension ALSOMU(maxlso),DSOMU(maxdso)
      Dimension ITZ1D(NMAX,NMAX)
      Dimension ijx1d(nxmax,nxmax)
      Dimension ijy1d(nymax,nymax)
      Dimension ijz1d(nzmax,nzmax)
      Dimension ILSOMU(NMAX1,NMAX1),IDSOMU(NMAX,NMAX)
c-------------------------------------------------------------- internal
      Dimension scrz1p(2,izsrc),scrz2p(2,izsrc),scrz3p(2,izsrc)
      Dimension scrz1n(2,izsrc),scrz2n(2,izsrc),scrz3n(2,izsrc)
      Dimension scrz4p(2,izsrc),scrz4n(2,izsrc),scrzc(izsrc)
      Dimension scryz1p(2,iyzsrc),scryz2p(2,iyzsrc),scryz3p(2,iyzsrc)
      Dimension scryz1n(2,iyzsrc),scryz2n(2,iyzsrc),scryz3n(2,iyzsrc)
      Dimension scryz4p(2,iyzsrc),scryz4n(2,iyzsrc),scryzc(iyzsrc)
c
      Dimension Sozp (izsrc,10),Sozep(izsrc,6)
      Dimension Sozn (izsrc,10),Sozen(izsrc,6)
      Dimension Soyzp(iyzsrc,10),Soyzep(iyzsrc,6)
      Dimension Soyzn(iyzsrc,10),Soyzen(iyzsrc,6)
c--------------------------------------------------------- end internal
      Dimension sxyz1p(2,ndmu),sxyz2p(2,ndmu),sxyz3p(2,ndmu)
      Dimension sxyz1n(2,ndmu),sxyz2n(2,ndmu),sxyz3n(2,ndmu)
      Dimension sxyz4p(2,ndmu),sxyz4n(2,ndmu),sxyzc(ndmu)
c
      Dimension Soxyzp(10,ndmu)
      Dimension Soxyzn(10,ndmu)
c +---------------------------------------------------------------------+
c |   In Soz, Soyz y Soxyz the meaning of the second index is           |
c +---------------------------------------------------------------------+
c |                         nxmu     nymu     nzmu                      |
c |      1 ..... gamma_x    even     even     even                      |
c |      2 ..... gamma_y     odd      odd     even    With TZ           |
c |      3 ..... gamma_z     odd     even      odd                      |
c |      4 ..... Ro mu      even     even     even                      |
c +---------------------------------------------------------------------+
c |      5 ..... Delta_x    even     even     even    (2)               |
c |      6 ..... Delta_y     odd      odd     even    (2)               |
c |      7 ..... Delta_z     odd     even      odd    (2)               |
c |                                                                     |
c |    The three Sigma (0) components                                   |
c |                                                                     |
c |      8 ..... sigma_x    even     even     even    (2)               |
c |      9 ..... sigma_y    even     even     even    (1)               |
c |     10 ..... sigma_z    even     even     even    (3)               |
c |                                                                     |
c +---------------------------------------------------------------------+
      Dimension GROP(nrop8),GROM(nrom8)
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      Common/DIMEN/kmax,kmax1,kxmax,kymax,kzmax,nwf2,
     * maxtz1,maxjz1,maxjx1,maxjy1,maxlso,maxdso,
     * ndmu,ndthet,nacou,nrop,nrom,nrop8,nrom8,nherm,nherm38,nlegn
c
      nrop1 = nrop
      nrop2 = nrop * 2
      nrop3 = nrop * 3
      nrom1 = nrom
      nrom2 = nrom * 2
      nrom3 = nrom * 3
c
      mumax = 2*nxmax - 1
      mumay = 2*nymax - 1
      mumaz = 2*nzmax - 1
      irop = 0
      irom = 0
c -----------------------------mux par   muy par  muz par Campo (0)
          imuxyz = 0
          do mux=1,mumax,2
            do muy=1,mumay,2
              do muz=1,mumaz,2
                imuxyz = imuxyz + 1
                sxyz1p(1,imuxyz) = 0.0d+0
                sxyz2p(1,imuxyz) = 0.0d+0
                sxyz1p(2,imuxyz) = 0.0d+0
                sxyz2p(2,imuxyz) = 0.0d+0
                sxyzc (imuxyz)   = 0.0d+0
                soxyzp(1,imuxyz) = 0.0d+0
                soxyzp(4,imuxyz) = 0.0d+0
                soxyzp(5,imuxyz) = 0.0d+0
                soxyzp(8,imuxyz) = 0.0d+0
                soxyzp(9,imuxyz) = 0.0d+0
                soxyzp(10,imuxyz)= 0.0d+0
                sxyz1n(1,imuxyz) = 0.0d+0
                sxyz2n(1,imuxyz) = 0.0d+0
                sxyz1n(2,imuxyz) = 0.0d+0
                sxyz2n(2,imuxyz) = 0.0d+0
                soxyzn(1,imuxyz) = 0.0d+0
                soxyzn(4,imuxyz) = 0.0d+0
                soxyzn(5,imuxyz) = 0.0d+0
                soxyzn(8,imuxyz) = 0.0d+0
                soxyzn(9,imuxyz) = 0.0d+0
                soxyzn(10,imuxyz)= 0.0d+0
              end do
            end do
          end do
c -----------------------------mux impar muy par  muz impar Campo (1)
          imuxyz = 0
          do mux=2,mumax,2
            do  muy=1,mumay,2
              do muz=2,mumaz,2
                imuxyz = imuxyz + 1
                sxyz3p(1,imuxyz) = 0.0d+0
                sxyz3p(2,imuxyz) = 0.0d+0
                soxyzp(3,imuxyz ) = 0.0d+0
                soxyzp(7,imuxyz ) = 0.0d+0
                sxyz3n(1,imuxyz) = 0.0d+0
                sxyz3n(2,imuxyz) = 0.0d+0
                soxyzn(3,imuxyz ) = 0.0d+0
                soxyzn(7,imuxyz ) = 0.0d+0
              end do
            end do
          end do
c ------------------------------mux impar muy impar muz par Campo (3)
          imuxyz = 0
          do  mux=2,mumax,2
            do muy=2,mumay,2
              do muz=1,mumaz,2
                imuxyz = imuxyz + 1
                sxyz4p(1,imuxyz) = 0.0d+0
                sxyz4p(2,imuxyz) = 0.0d+0
                soxyzp(2,imuxyz ) = 0.0d+0
                soxyzp(6,imuxyz ) = 0.0d+0
                sxyz4n(1,imuxyz) = 0.0d+0
                sxyz4n(2,imuxyz) = 0.0d+0
                soxyzn(2,imuxyz ) = 0.0d+0
                soxyzn(6,imuxyz ) = 0.0d+0
              end do
            end do
          end do
c ------------------------------------------------------------------
c                              1 -> q'     2 -> q
      do 1 nx1=1,nxmax
        nym1 = MY(nx1)
        do 11 nx2=1,nx1
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          lpx  = Mod(nx1+nx2,2).eq.0
c -------------------------------------- muy par  muz par Campo (0)
          imuyz = 0
          do muy=1,mumay,2
            do muz=1,mumaz,2
              imuyz = imuyz + 1
              scryz1p(1,imuyz) = 0.0d+00
              scryz2p(1,imuyz) = 0.0d+00
              scryz1p(2,imuyz) = 0.0d+00
              scryz2p(2,imuyz) = 0.0d+00
              scryzc (imuyz  ) = 0.0d+00
              soyzp(imuyz,1)   = 0.0d+00
              soyzp(imuyz,4)   = 0.0d+00
              soyzp(imuyz,5)   = 0.0d+00
              soyzp(imuyz,8)   = 0.0d+00
              soyzp(imuyz,9)   = 0.0d+00
              soyzp(imuyz,10)  = 0.0d+00
              soyzep(imuyz,1)  = 0.0d+00
              soyzep(imuyz,4)  = 0.0d+00
              soyzep(imuyz,5)  = 0.0d+00
              soyzep(imuyz,6)  = 0.0d+00
              scryz1n(1,imuyz) = 0.0d+00
              scryz2n(1,imuyz) = 0.0d+00
              scryz1n(2,imuyz) = 0.0d+00
              scryz2n(2,imuyz) = 0.0d+00
              soyzn(imuyz,1)   = 0.0d+00
              soyzn(imuyz,4)   = 0.0d+00
              soyzn(imuyz,5)   = 0.0d+00
              soyzn(imuyz,8)   = 0.0d+00
              soyzn(imuyz,9)   = 0.0d+00
              soyzn(imuyz,10)  = 0.0d+00
              soyzen(imuyz,1)  = 0.0d+00
              soyzen(imuyz,4)  = 0.0d+00
              soyzen(imuyz,5)  = 0.0d+00
              soyzen(imuyz,6)  = 0.0d+00
            end do
          end do
c -------------------------------------- muy par  muz impar Campo (1)
          imuyz = 0
          do muy=1,mumay,2
            do muz=2,mumaz,2
              imuyz = imuyz + 1
              scryz3p(1,imuyz) = 0.0d+00
              scryz3p(2,imuyz) = 0.0d+00
              soyzp(imuyz,3) = 0.0d+00
              soyzp(imuyz,7) = 0.0d+00
              soyzep(imuyz,3)= 0.0d+00
              scryz3n(1,imuyz) = 0.0d+00
              scryz3n(2,imuyz) = 0.0d+00
              soyzn(imuyz,3) = 0.0d+00
              soyzn(imuyz,7) = 0.0d+00
              soyzen(imuyz,3)= 0.0d+00
            end do
          end do
c -------------------------------------- muy impar muz par Campo (3)
          imuyz = 0
          do muy=2,mumay,2
            do muz=1,mumaz,2
              imuyz = imuyz + 1
              scryz4p(1,imuyz ) = 0.0
              scryz4p(2,imuyz ) = 0.0
              soyzp  (imuyz,2)  = 0.0
              soyzp  (imuyz,6)  = 0.0
              soyzep (imuyz,2)  = 0.0
              scryz4n(1,imuyz ) = 0.0
              scryz4n(2,imuyz ) = 0.0
              soyzn  (imuyz,2)  = 0.0
              soyzn  (imuyz,6)  = 0.0
              soyzen (imuyz,2)  = 0.0
            end do
          end do
c ------------------------------------------------------------------
c
          do 2 ny1 = 1,nym1
            nxy1 = nx1+ny1
            nzi1e = nzie(nx1,ny1)
            nzi1o = nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
            if(lx12) nym2 = ny1
            do 12 ny2=1,nym2
              nxy2 = nx2+ny2
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              lxy12 = lx12.and.(ny1.eq.ny2)
              lpy  = Mod(ny1+ny2,2).eq.0
c ------------------------------------- muz par
              imu = 0
              do muz =1,mumaz,2
                 imu=imu+1
                 scrzc(imu) = 0.0d+00
                 scrz1p(1,imu) = 0.0d+00
                 scrz2p(1,imu) = 0.0d+00
                 scrz4p(1,imu) = 0.0d+00
                 scrz1p(2,imu) = 0.0d+00
                 scrz2p(2,imu) = 0.0d+00
                 scrz4p(2,imu) = 0.0d+00
                 sozp(imu,1)   = 0.0d+00
                 sozp(imu,2)   = 0.0d+00
                 sozp(imu,4)   = 0.0d+00
                 sozp(imu,5)   = 0.0d+00
                 sozp(imu,6)   = 0.0d+00
                 sozp(imu,8)   = 0.0d+00
                 sozp(imu,9)   = 0.0d+00
                 sozp(imu,10)  = 0.0d+00
                 sozep(imu,1)  = 0.0d+00
                 sozep(imu,2)  = 0.0d+00
                 sozep(imu,4)  = 0.0d+00
                 sozep(imu,5)  = 0.0d+00
                 sozep(imu,6)  = 0.0d+00
                 scrz1n(1,imu) = 0.0d+00
                 scrz2n(1,imu) = 0.0d+00
                 scrz4n(1,imu) = 0.0d+00
                 scrz1n(2,imu) = 0.0d+00
                 scrz2n(2,imu) = 0.0d+00
                 scrz4n(2,imu) = 0.0d+00
                 sozn(imu,1)   = 0.0d+00
                 sozn(imu,2)   = 0.0d+00
                 sozn(imu,4)   = 0.0d+00
                 sozn(imu,5)   = 0.0d+00
                 sozn(imu,6)   = 0.0d+00
                 sozn(imu,8)   = 0.0d+00
                 sozn(imu,9)   = 0.0d+00
                 sozn(imu,10)  = 0.0d+00
                 sozen(imu,1)  = 0.0d+00
                 sozen(imu,2)  = 0.0d+00
                 sozen(imu,4)  = 0.0d+00
                 sozen(imu,5)  = 0.0d+00
                 sozen(imu,6)  = 0.0d+00
              end do  ! muz
c ------------------------------------- muz impar
              imu = 0
              do muz =2,mumaz,2
                 imu=imu+1
                 scrz3p(1,imu) = 0.0d+00
                 scrz3p(2,imu) = 0.0d+00
                 sozp(imu,3)   = 0.0d+00
                 sozp(imu,7)   = 0.0d+00
                 sozep(imu,3)  = 0.0d+00   ! soz(imu,7)
                 scrz3n(1,imu) = 0.0d+00
                 scrz3n(2,imu) = 0.0d+00
                 sozn(imu,3)   = 0.0d+00
                 sozn(imu,7)   = 0.0d+00
                 sozen(imu,3)  = 0.0d+00   ! soz(imu,7)
              end do   ! muz
c+---------------------------------------------------------------------+
c|                                                 Positive parity     |
c+---------------------------------------------------------------------+
              if(nzi1e.gt.nzm1) goto 20
              do 3 nz1 =nzi1e,nzm1,2
                if(lxy12) nzm2=nz1
                if(nzi2e.gt.nzm2) goto 3
                do 13 nz2 = nzi2e,nzm2,2
c
                 irop  = irop + 1
                 rospen = GROP(irop      )
                 roapen = GROP(irop+1)
                 rospep = GROP(irop+2)
                 roapep = GROP(irop+3)
                 irop   = irop + 7
                 lxyz = lxy12.and.(nz2.eq.nz1)
                 if(.not.lxyz) then
                    rospen = rospen + rospen
                    roapen = roapen + roapen
                    rospep = rospep + rospep
                    roapep = roapep + roapep
                 end if
                 ijz = ijz1d(nz1,nz2)
                 itz = itz1d(nz1,nz2)
                 indx = ijz
                 indc = itz
                 ils = ilsomu(nz1,nz2)
                 ids = idsomu(nz1,nz2)
                 idse= idsomu(nz2,nz1)
                 indls= ils
                 indds= ids
                 inddse= idse
                 muzic= Iabs(nz2-nz1) + 1
                 muzfc= nz2+nz1-1
c+--------------------------------------------+
c|         Field (0)     muz par              |
c|                                            |
c|    scrz1    Brink-Boecker Delta mu         |
c|    scrz2    Brink-Boecker S1    mu         |
c|    scrzc    Coulomb directo     mu         |
c|    soz 1    gamma x             mu         |
c|    soz 4    Ro                  mu         |
c+--------------------------------------------+
                   if(lpy.and.lpx) then
                   imu = 0
                   do 100 muz=1,mumaz,2
                     imu = imu + 1
                     indx = indx + 1
                     scrz1n(1,imu)=scrz1n(1,imu)+rospen*ajz1d(1,indx)
                     scrz2n(1,imu)=scrz2n(1,imu)+roapen*ajz1d(1,indx)
                     scrz1n(2,imu)=scrz1n(2,imu)+rospen*ajz1d(2,indx)
                     scrz2n(2,imu)=scrz2n(2,imu)+roapen*ajz1d(2,indx)
                     scrz1p(1,imu)=scrz1p(1,imu)+rospep*ajz1d(1,indx)
                     scrz2p(1,imu)=scrz2p(1,imu)+roapep*ajz1d(1,indx)
                     scrz1p(2,imu)=scrz1p(2,imu)+rospep*ajz1d(2,indx)
                     scrz2p(2,imu)=scrz2p(2,imu)+roapep*ajz1d(2,indx)
                     if(muz.lt.muzic.or.muz.gt.muzfc) go to 100
                     indc = indc + 1
                     sozn(imu,1) = sozn(imu,1) + roapen*tz1d(indc)
                     sozn(imu,4) = sozn(imu,4) + rospen*tz1d(indc)
                     sozp(imu,1) = sozp(imu,1) + roapep*tz1d(indc)
                     sozp(imu,4) = sozp(imu,4) + rospep*tz1d(indc)
                     scrzc(imu) = scrzc(imu) + rospep*tz1d(indc)
  100              continue
c+--------------------------------------------+
c|         Field (1)   muz impar              |
c|                                            |
c|    scrz3    Brink-Boecker S2    mu         |
c|    soz 3    gamma z             mu         |
c|    soz 9    Sigma y             mu  muz par|
c+--------------------------------------------+
                   else if(lpy.and.(.not.lpx)) then
                   imu = 0
                   do 101 muz=2,mumaz,2
                     imu = imu + 1
                     indx = indx + 1
                     scrz3n(1,imu)=scrz3n(1,imu) +roapen*ajz1d(1,indx)
                     scrz3n(2,imu)=scrz3n(2,imu) +roapen*ajz1d(2,indx)
                     scrz3p(1,imu)=scrz3p(1,imu) +roapep*ajz1d(1,indx)
                     scrz3p(2,imu)=scrz3p(2,imu) +roapep*ajz1d(2,indx)
                     if(muz.lt.muzic.or.muz.gt.muzfc) go to 101
                     indc = indc + 1
                     sozn(imu,3) = sozn(imu,3) + roapen*tz1d(indc)
                     sozp(imu,3) = sozp(imu,3) + roapep*tz1d(indc)
  101              continue
                   imu = 0
                   do muz=1,mumaz,2
                     imu = imu + 1
                     sozn (imu,9) =sozn (imu,9)+ rospen*dsomu(indds )
                     sozen(imu,5) =sozen(imu,5)- rospen*dsomu(inddse)
                     sozp (imu,9) =sozp (imu,9)+ rospep*dsomu(indds )
                     sozep(imu,5) =sozep(imu,5)- rospep*dsomu(inddse)
                     indds = indds + 1
                     inddse= inddse+ 1
                   end do
c+--------------------------------------------+
c|         Field (2)   muz   par              |
c|  Delta   soz (5,6)                         |
c|  sigmax  soz (8)                           |
c|                     muz impar              |
c|  Delta z soz (7)                           |
c+--------------------------------------------+
                   else if(.not.lpy.and.lpx) then
                   imu = 0
                   do 102 muz=1,mumaz,2
                     imu = imu + 1
        sozn (imu,8) =sozn (imu,8)- rospen*dsomu(indds )
        sozen(imu,4) =sozen(imu,4)+ rospen*dsomu(inddse)
        sozn (imu,5) =sozn (imu,5)- roapen*dsomu(indds ) ! delta x
        sozen(imu,1) =sozen(imu,1)+ roapen*dsomu(inddse)
        sozn (imu,6) =sozn (imu,6)+ roapen*dsomu(indds ) ! delta y
        sozen(imu,2) =sozen(imu,2)- roapen*dsomu(inddse)
        sozp (imu,8) =sozp (imu,8)- rospep*dsomu(indds )
        sozep(imu,4) =sozep(imu,4)+ rospep*dsomu(inddse)
        sozp (imu,5) =sozp (imu,5)- roapep*dsomu(indds ) ! delta x
        sozep(imu,1) =sozep(imu,1)+ roapep*dsomu(inddse)
        sozp (imu,6) =sozp (imu,6)+ roapep*dsomu(indds ) ! delta y
        sozep(imu,2) =sozep(imu,2)- roapep*dsomu(inddse)
                     indds = indds + 1
                     inddse= inddse+ 1
  102                continue
                   imu = 0
                   do muz =2,mumaz,2
                     imu = imu + 1
        sozn (imu,7) =sozn (imu,7)+ roapen*alsomu(indls) ! delta z
        sozp (imu,7) =sozp (imu,7)+ roapep*alsomu(indls) ! delta z
                     indls = indls + 1
                   end do
c+--------------------------------------------+
c|         Field (3)    muz par               |
c|                                            |
c|    scrz4    Brink-Boecker S1    mu         |
c|    soz 2    gamma y             mu         |
c|    soz 10   sigma z             mu         |
c+--------------------------------------------+
                   else if((.not.lpy).and.(.not.lpx)) then
                   imu = 0
                   do 103 muz=1,mumaz,2
                     imu = imu + 1
                     indx = indx + 1
                     scrz4n(1,imu)=scrz4n(1,imu)+roapen*ajz1d(1,indx)
                     scrz4n(2,imu)=scrz4n(2,imu)+roapen*ajz1d(2,indx)
                     sozn(imu,10) =sozn(imu,10)+ rospen*alsomu(indls) ! sigma z
                     scrz4p(1,imu)=scrz4p(1,imu)+roapep*ajz1d(1,indx)
                     scrz4p(2,imu)=scrz4p(2,imu)+roapep*ajz1d(2,indx)
                     sozp(imu,10) =sozp(imu,10)+ rospep*alsomu(indls) ! sigma z
                     indls = indls + 1
                     if(muz.lt.muzic.or.muz.gt.muzfc) go to 103
                     indc = indc + 1
                     sozn(imu,2) = sozn(imu,2) + roapen*tz1d(indc)
                     sozp(imu,2) = sozp(imu,2) + roapep*tz1d(indc)
  103                continue
c
                   end if
c
   13           continue
    3         continue
c+--------------------------------------------------+
c|                                                  |
c|         Paridad negativa                         |
c|                                                  |
c+--------------------------------------------------+
   20 continue
              if(nzi1o.gt.nzm1) goto 21
              do 4 nz1 =nzi1o,nzm1,2
                if(lxy12) nzm2=nz1
                if(nzi2o.gt.nzm2) goto 4
                do 14 nz2 = nzi2o,nzm2,2
c
                 irom  = irom + 1
                 rosmen = GROM(irom      )
                 roamen = GROM(irom+1)
                 rosmep = GROM(irom+2)
                 roamep = GROM(irom+3)
                 irom   = irom + 7
                 lxyz = lxy12.and.(nz2.eq.nz1)
                 if(.not.lxyz) then
                    rosmen = rosmen + rosmen
                    roamen = roamen + roamen
                    rosmep = rosmep + rosmep
                    roamep = roamep + roamep
                 end if
                 ijz = ijz1d(nz1,nz2)
                 itz = itz1d(nz1,nz2)
                 indx = ijz
                 indc = itz
                 ils = ilsomu(nz1,nz2)
                 ids = idsomu(nz1,nz2)
                 idse= idsomu(nz2,nz1)
                 indls= ils
                 indds= ids
                 inddse= idse
                 muzic= Iabs(nz2-nz1) + 1
                 muzfc= nz2+nz1-1
c+--------------------------------------------+
c|         Field (0)     muz par              |
c|                                            |
c|    scrz1    Brink-Boecker Delta mu         |
c|    scrz2    Brink-Boecker S1    mu         |
c|    scrzc    Coulomb directo     mu         |
c+--------------------------------------------+
                   if(lpy.and.lpx) then
                   imu = 0
                   do 200 muz=1,mumaz,2
                     imu = imu + 1
                     indx = indx + 1
                     scrz1n(1,imu)=scrz1n(1,imu)+rosmen*ajz1d(1,indx)
                     scrz2n(1,imu)=scrz2n(1,imu)+roamen*ajz1d(1,indx)
                     scrz1n(2,imu)=scrz1n(2,imu)+rosmen*ajz1d(2,indx)
                     scrz2n(2,imu)=scrz2n(2,imu)+roamen*ajz1d(2,indx)
                     scrz1p(1,imu)=scrz1p(1,imu)+rosmep*ajz1d(1,indx)
                     scrz2p(1,imu)=scrz2p(1,imu)+roamep*ajz1d(1,indx)
                     scrz1p(2,imu)=scrz1p(2,imu)+rosmep*ajz1d(2,indx)
                     scrz2p(2,imu)=scrz2p(2,imu)+roamep*ajz1d(2,indx)
                     if(muz.lt.muzic.or.muz.gt.muzfc) go to 200
                     indc = indc + 1
                     sozn(imu,1) = sozn(imu,1) + roamen*tz1d(indc)
                     sozn(imu,4) = sozn(imu,4) + rosmen*tz1d(indc)
                     sozp(imu,1) = sozp(imu,1) + roamep*tz1d(indc)
                     sozp(imu,4) = sozp(imu,4) + rosmep*tz1d(indc)
                     scrzc(imu) = scrzc(imu) + rosmep*tz1d(indc)
  200              continue
c+--------------------------------------------+
c|         Field (1)   muz impar              |
c|                                            |
c|    scrz3    Brink-Boecker S2    mu         |
c+--------------------------------------------+
                   else if(lpy.and.(.not.lpx)) then
                   imu = 0
                   do 201 muz=2,mumaz,2
                     imu = imu + 1
                     indx = indx + 1
                     scrz3n(1,imu)=scrz3n(1,imu)+roamen*ajz1d(1,indx)
                     scrz3n(2,imu)=scrz3n(2,imu)+roamen*ajz1d(2,indx)
                     scrz3p(1,imu)=scrz3p(1,imu)+roamep*ajz1d(1,indx)
                     scrz3p(2,imu)=scrz3p(2,imu)+roamep*ajz1d(2,indx)
                     if(muz.lt.muzic.or.muz.gt.muzfc) go to 201
                     indc = indc + 1
                     sozn(imu,3) = sozn(imu,3) + roamen*tz1d(indc)
                     sozp(imu,3) = sozp(imu,3) + roamep*tz1d(indc)
  201              continue
                   imu = 0
                   do muz=1,mumaz,2
                     imu = imu + 1
                     sozn (imu,9) =sozn (imu,9)+ rosmen*dsomu(indds )
                     sozen(imu,5) =sozen(imu,5)- rosmen*dsomu(inddse)
                     sozp (imu,9) =sozp (imu,9)+ rosmep*dsomu(indds )
                     sozep(imu,5) =sozep(imu,5)- rosmep*dsomu(inddse)
                     indds = indds + 1
                     inddse= inddse+ 1
                   end do
c+--------------------------------------------+
c|         Field (2)   muz   par              |
c|  Delta   soz (5,6)                         |
c|  sigmax  soz (8)                           |
c|                     muz impar              |
c|  Delta z soz (7)                           |
c+--------------------------------------------+
                   else if(.not.lpy.and.lpx) then
                   imu = 0
                   do 202 muz=1,mumaz,2
                     imu = imu + 1
                     sozn (imu,8) =sozn (imu,8)- rosmen*dsomu(indds )
                     sozen(imu,4) =sozen(imu,4)+ rosmen*dsomu(inddse)
                     sozn (imu,5) =sozn (imu,5)- roamen*dsomu(indds ) ! delta x
                     sozen(imu,1) =sozen(imu,1)+ roamen*dsomu(inddse)
                     sozn (imu,6) =sozn (imu,6)+ roamen*dsomu(indds ) ! delta y
                     sozen(imu,2) =sozen(imu,2)- roamen*dsomu(inddse)
                     sozp (imu,8) =sozp (imu,8)- rosmep*dsomu(indds )
                     sozep(imu,4) =sozep(imu,4)+ rosmep*dsomu(inddse)
                     sozp (imu,5) =sozp (imu,5)- roamep*dsomu(indds ) ! delta x
                     sozep(imu,1) =sozep(imu,1)+ roamep*dsomu(inddse)
                     sozp (imu,6) =sozp (imu,6)+ roamep*dsomu(indds ) ! delta y
                     sozep(imu,2) =sozep(imu,2)- roamep*dsomu(inddse)
                     indds = indds + 1
                     inddse= inddse+ 1
  202                continue
                   imu = 0
                   do muz =2,mumaz,2
                     imu = imu + 1
                     sozn (imu,7) =sozn (imu,7)+ roamen*alsomu(indls) ! delta z
                     sozp (imu,7) =sozp (imu,7)+ roamep*alsomu(indls) ! delta z
                     indls = indls + 1
                   end do
c+--------------------------------------------+
c|         Field (3)    muz par               |
c|                                            |
c|    scrz4    Brink-Boecker S1    mu         |
c+--------------------------------------------+
                   else if((.not.lpy).and.(.not.lpx)) then
                   imu = 0
                   do 203 muz=1,mumaz,2
                     imu = imu + 1
                     indx = indx + 1
                     scrz4n(1,imu)=scrz4n(1,imu)+roamen*ajz1d(1,indx)
                     scrz4n(2,imu)=scrz4n(2,imu)+roamen*ajz1d(2,indx)
                     sozn(imu,10) =sozn(imu,10)+ rosmen*alsomu(indls) ! sigma z
                     scrz4p(1,imu)=scrz4p(1,imu)+roamep*ajz1d(1,indx)
                     scrz4p(2,imu)=scrz4p(2,imu)+roamep*ajz1d(2,indx)
                     sozp(imu,10) =sozp(imu,10)+ rosmep*alsomu(indls) ! sigma z
                     indls = indls + 1
                     if(muz.lt.muzic.or.muz.gt.muzfc) go to 203
                     indc = indc + 1
                     sozn(imu,2) = sozn(imu,2) + roamen*tz1d(indc)
                     sozp(imu,2) = sozp(imu,2) + roamep*tz1d(indc)
  203                continue
c
                   end if
c
   14           continue
    4         continue
   21       continue
c
c            +------------------------+
c            |      Y     L O O P     |
c            +------------------------+
c
            ijy = ijy1d(ny1,ny2)
            indx = ijy
            itz = itz1d(ny1,ny2)
            indt = itz
            ils  = ilsomu(ny1,ny2)
            indls= ils
            ids  = idsomu(ny1,ny2)
            indds= ids
            idse = idsomu(ny2,ny1)
            inddse= idse
c
            ny2m1= iabs(ny2-ny1)
            muiy = ny2m1 + 1
            mufy = ny2+ny1-1
c
c
c    +---------------------+
c    |  C A M P O  ( 0 )   |
c    +---------------------+
c
c
            if(lpy.and.lpx) then
c -------------------------------------- muy par  muz par
              fasy = Dfloat(1-2*Mod(Iabs(ny2-ny1)/2,2))
c                                                   (-)**(ny2-ny1)/2
              imuyz = 0
              do 44 muy=1,mumay,2
                indx = indx + 1
                aa1  = ajy1d(1,indx)
                aa2  = ajy1d(2,indx)
                bb   = 0.0d+00
                if (muy.ge.muiy.and.muy.le.mufy) then
                    indt = indt + 1
                    bb = tz1d(indt)
                end if
                imuz = 0
                do 44 muz=1,mumaz,2
                imuz = imuz + 1
                imuyz = imuyz + 1
      scryz1n(1,imuyz)=scryz1n(1,imuyz)+fasy*scrz1n(1,imuz)*aa1
      scryz2n(1,imuyz)=scryz2n(1,imuyz)+fasy*scrz2n(1,imuz)*aa1
      scryz1n(2,imuyz)=scryz1n(2,imuyz)+fasy*scrz1n(2,imuz)*aa2
      scryz2n(2,imuyz)=scryz2n(2,imuyz)+fasy*scrz2n(2,imuz)*aa2
      soyzn  (imuyz,1) = soyzn(imuyz,1) + fasy*sozn(imuz,1)*bb
      soyzn  (imuyz,4) = soyzn(imuyz,4) + fasy*sozn(imuz,4)*bb
      scryz1p(1,imuyz)=scryz1p(1,imuyz)+fasy*scrz1p(1,imuz)*aa1
      scryz2p(1,imuyz)=scryz2p(1,imuyz)+fasy*scrz2p(1,imuz)*aa1
      scryz1p(2,imuyz)=scryz1p(2,imuyz)+fasy*scrz1p(2,imuz)*aa2
      scryz2p(2,imuyz)=scryz2p(2,imuyz)+fasy*scrz2p(2,imuz)*aa2
      soyzp  (imuyz,1) = soyzp(imuyz,1) + fasy*sozp(imuz,1)*bb
      soyzp  (imuyz,4) = soyzp(imuyz,4) + fasy*sozp(imuz,4)*bb
   44 scryzc(imuyz) = scryzc(imuyz)+ fasy*scrzc(imuz)*bb
c
c    +---------------------+
c    |  C A M P O  ( 1 )   |
c    +---------------------+
c
            else if(lpy.and..not.lpx) then
c -------------------------------------- muy par  muz impar
              fasy = Dfloat(1-2*Mod(Iabs(ny2-ny1)/2,2))
c                                                   (-)**(ny2-ny1)/2
              imuyzo = 0
              imuyze = 0
              do 45 muy=1,mumay,2
               indx = indx+1
               imuz = 0
               aa1= ajy1d(1,indx)
               aa2= ajy1d(2,indx)
               sls= alsomu(indls)
               indls = indls + 1
               bb = 0.0d+00
               if (muy.ge.muiy.and.muy.le.mufy) then
                   indt = indt + 1
                   bb = tz1d(indt)
               end if
               do 46 muz=2,mumaz,2
                imuz = imuz + 1
                imuyzo = imuyzo + 1
      soyzn  (imuyzo,3) = soyzn(imuyzo,3) + fasy*sozn(imuz,3)*bb
      scryz3n(1,imuyzo)=scryz3n(1,imuyzo)+fasy*scrz3n(1,imuz)*aa1
      scryz3n(2,imuyzo)=scryz3n(2,imuyzo)+fasy*scrz3n(2,imuz)*aa2
      soyzp  (imuyzo,3) = soyzp(imuyzo,3) + fasy*sozp(imuz,3)*bb
      scryz3p(1,imuyzo)=scryz3p(1,imuyzo)+fasy*scrz3p(1,imuz)*aa1
      scryz3p(2,imuyzo)=scryz3p(2,imuyzo)+fasy*scrz3p(2,imuz)*aa2
   46          continue
c -------------------------------------- muy par  muz par
               imuz = 0
               do 47 muz=1,mumaz,2
                imuz = imuz + 1
                imuyze = imuyze + 1
       soyzn (imuyze,9)=soyzn(imuyze,9) +fasy*sozn (imuz,9)*sls
       soyzen(imuyze,5)=soyzen(imuyze,5)+fasy*sozen(imuz,5)*sls
       soyzp (imuyze,9)=soyzp(imuyze,9) +fasy*sozp (imuz,9)*sls
       soyzep(imuyze,5)=soyzep(imuyze,5)+fasy*sozep(imuz,5)*sls
   47          continue
   45         continue
c
c     +------------------+
c     | C A M P O  ( 2 ) |
c     +------------------+
c
            else if(.not.lpy.and.lpx) then
c -------------------------------------- muy impar  muz impar
              fasy = Dfloat(1-2*Mod(Iabs(ny2-ny1-1)/2,2))
c                                                   (-)**(ny2-ny1-1)/2
              fasy1= Dfloat(1-2*Mod(Iabs(ny2-ny1+1)/2,2))
c                                                   (-)**(ny2-ny1+1)/2
              imuyzo = 0
              imuyze = 0
              do 48 muy=1,mumay,2
               imuz = 0
               sls= alsomu(indls)
               indls = indls + 1
               sds= dsomu (indds)
               indds = indds + 1
               sdse= dsomu (inddse)
               inddse = inddse + 1
               do 49 muz=2,mumaz,2
                imuz = imuz + 1
                imuyzo = imuyzo + 1
        soyzn (imuyzo,7)=soyzn(imuyzo,7) -fasy1*sozn (imuz,7)*sds
        soyzen(imuyzo,3)=soyzen(imuyzo,3)+fasy1*sozn (imuz,7)*sdse
        soyzp (imuyzo,7)=soyzp(imuyzo,7) -fasy1*sozp (imuz,7)*sds
        soyzep(imuyzo,3)=soyzep(imuyzo,3)+fasy1*sozp (imuz,7)*sdse
   49          continue
c -------------------------------------- muy impar  muz par
               imuz = 0
               do muz=1,mumaz,2
                imuz = imuz + 1
                imuyze = imuyze + 1
       soyzn (imuyze,5)=soyzn (imuyze,5)+fasy1*sozn (imuz,5)*sdse
       soyzen(imuyze,1)=soyzen(imuyze,1)+fasy1*sozen(imuz,1)*sds
       soyzn (imuyze,6)=soyzn (imuyze,6)+fasy1*sozn (imuz,6)*sls
       soyzen(imuyze,2)=soyzen(imuyze,2)+fasy1*sozen(imuz,2)*sls
       soyzn (imuyze,8)=soyzn (imuyze,8)+fasy1*sozn (imuz,8)*sdse
       soyzen(imuyze,4)=soyzen(imuyze,4)+fasy1*sozen(imuz,4)*sds
       soyzp (imuyze,5)=soyzp (imuyze,5)+fasy1*sozp (imuz,5)*sdse
       soyzep(imuyze,1)=soyzep(imuyze,1)+fasy1*sozep(imuz,1)*sds
       soyzp (imuyze,6)=soyzp (imuyze,6)+fasy1*sozp (imuz,6)*sls
       soyzep(imuyze,2)=soyzep(imuyze,2)+fasy1*sozep(imuz,2)*sls
       soyzp (imuyze,8)=soyzp (imuyze,8)+fasy1*sozp (imuz,8)*sdse
       soyzep(imuyze,4)=soyzep(imuyze,4)+fasy1*sozep(imuz,4)*sds
               end do ! muz
   48         continue
c
c
c  +--------------------------+
c  |    C A M P O  ( 3 )      |
c  +--------------------------+
c
            else if((.not.lpy).and.(.not.lpx)) then
c -------------------------------------- muy impar muz par
              fasy = Dfloat(1-2*Mod(Iabs(ny2-ny1-1)/2,2))
c                                                   (-)**(ny2-ny1-1)/2
              fasy1= Dfloat(1-2*Mod(Iabs(ny2-ny1+1)/2,2))
c                                                   (-)**(ny2-ny1+1)/2
              imuyz = 0
              do muy=2,mumay,2
               indx = indx+1
               imuz = 0
               aa1= ajy1d(1,indx)
               aa2= ajy1d(2,indx)
               bb = 0.0d+00
               if (muy.ge.muiy.and.muy.le.mufy) then
                   indt = indt + 1
                   bb = tz1d(indt)
               end if
               do 52 muz=1,mumaz,2
                 imuz = imuz + 1
                 imuyz = imuyz + 1
       soyzn  (imuyz,2) = soyzn(imuyz,2) + fasy1*sozn(imuz,2)*bb
       scryz4n(1,imuyz)=scryz4n(1,imuyz)+fasy*scrz4n(1,imuz)*aa1
       scryz4n(2,imuyz)=scryz4n(2,imuyz)+fasy*scrz4n(2,imuz)*aa2
       soyzp  (imuyz,2) = soyzp(imuyz,2) + fasy1*sozp(imuz,2)*bb
       scryz4p(1,imuyz)=scryz4p(1,imuyz)+fasy*scrz4p(1,imuz)*aa1
       scryz4p(2,imuyz)=scryz4p(2,imuyz)+fasy*scrz4p(2,imuz)*aa2
   52          continue
              end do    ! muy
c -------------------------------------- muy par muz par
              imuyz = 0
              do 53 muy=1,mumay,2
               imuz = 0
               sds  = dsomu(indds)
               indds = indds + 1
               sdse = dsomu(inddse)
               inddse = inddse + 1
               do 54 muz=1,mumaz,2
                 imuz = imuz + 1
                 imuyz = imuyz + 1
       soyzn  (imuyz,10)=soyzn(imuyz,10)-fasy1*sozn(imuz,10)*sds
       soyzen (imuyz,6)=soyzen(imuyz,6) +fasy1*sozn(imuz,10)*sdse
       soyzp  (imuyz,10)=soyzp(imuyz,10)-fasy1*sozp(imuz,10)*sds
       soyzep (imuyz,6)=soyzep(imuyz,6) +fasy1*sozp(imuz,10)*sdse
   54          continue
   53        continue
c ------------------------------------------------------------------
            end if    !    C A M P O S   Y
c ------------------------------------------------------------------
   12       continue
    2     continue
c
c            +------------------------+
c            |      X     L O O P     |
c            +------------------------+
c
            ijx   = ijx1d(nx1,nx2)
            indx  = ijx
            itz   = itz1d(nx1,nx2)
            indt  = itz
            ils  = ilsomu(nx1,nx2)
            indls= ils
            ids  = idsomu(nx1,nx2)
            indds= ids
            idse = idsomu(nx2,nx1)
            inddse=idse
c
            fasx  = Dfloat(1-2*Mod(nx1+1,2))  !   (-)**nxq'
            fasx1 = Dfloat(1-2*Mod(nx1,2))    ! - (-)**nxq'
            muix  = Iabs(nx2-nx1) + 1
            mufx  = nx2+nx1-1
c                                       +------------+
c --------------------------------------| nx+nx' par |
c                                       +------------+
            if(lpx) then
c                                                     Campo (0)
              imuxyz = 0
              do 444 mux =1,mumax,2
                indx = indx + 1
                aa1  = ajx1d(1,indx)
                aa2  = ajx1d(2,indx)
                sls  = alsomu(indls)
                indls = indls + 1
                bb   = 0.0d+00
                if (mux.ge.muix.and.mux.le.mufx) then
                    indt = indt + 1
                    bb = tz1d(indt)
                end if
c
                imuyz = 0
                  do 444 muy=1,mumay,2
                    do 444 muz=1,mumaz,2
                    imuyz  = imuyz + 1
                    imuxyz = imuxyz + 1
        soxyzn (1,imuxyz) = soxyzn(1,imuxyz) + fasx1*soyzn(imuyz,1)*bb
        soxyzn (4,imuxyz) = soxyzn(4,imuxyz) +       soyzn(imuyz,4)*bb
        soxyzn (5,imuxyz) = soxyzn(5,imuxyz) +
     *  (soyzn(imuyz,5)+soyzen(imuyz,1))*sls
        soxyzn (8,imuxyz) = soxyzn(8,imuxyz) +
     *  (soyzn(imuyz,8)+soyzen(imuyz,4))*sls*fasx1
        soxyzp (1,imuxyz) = soxyzp(1,imuxyz) + fasx1*soyzp(imuyz,1)*bb
        soxyzp (4,imuxyz) = soxyzp(4,imuxyz) +       soyzp(imuyz,4)*bb
        soxyzp (5,imuxyz) = soxyzp(5,imuxyz) +
     *  (soyzp(imuyz,5)+soyzep(imuyz,1))*sls
        soxyzp (8,imuxyz) = soxyzp(8,imuxyz) +
     *  (soyzp(imuyz,8)+soyzep(imuyz,4))*sls*fasx1
        sxyz1n(1,imuxyz)=sxyz1n(1,imuxyz)+scryz1n(1,imuyz)*aa1
        sxyz2n(1,imuxyz)=sxyz2n(1,imuxyz)+scryz2n(1,imuyz)*aa1*fasx
        sxyz1n(2,imuxyz)=sxyz1n(2,imuxyz)+scryz1n(2,imuyz)*aa2
        sxyz2n(2,imuxyz)=sxyz2n(2,imuxyz)+scryz2n(2,imuyz)*aa2*fasx
        sxyz1p(1,imuxyz)=sxyz1p(1,imuxyz)+scryz1p(1,imuyz)*aa1
        sxyz2p(1,imuxyz)=sxyz2p(1,imuxyz)+scryz2p(1,imuyz)*aa1*fasx
        sxyz1p(2,imuxyz)=sxyz1p(2,imuxyz)+scryz1p(2,imuyz)*aa2
        sxyz2p(2,imuxyz)=sxyz2p(2,imuxyz)+scryz2p(2,imuyz)*aa2*fasx
  444   sxyzc(imuxyz)   =sxyzc(imuxyz)   + scryzc(imuyz)*bb
c
c
c                                                         campo (1)
c
              imuxyz = 0
              ipuxyz = 0
              do 445 mux =2,mumax,2
                sds   = dsomu(indds)
                indds = indds + 1
                sdse   = dsomu(inddse)
                inddse = inddse + 1
                imuyz = 0
c
                do 446 muy=1,mumay,2
                  do 446 muz=2,mumaz,2
                    imuyz  = imuyz  + 1
                    imuxyz = imuxyz + 1
                    soxyzn  (7,imuxyz) = soxyzn(7,imuxyz) +
     *              (soyzn(imuyz,7)*sdse+soyzen(imuyz,3)*sds)
                    soxyzp  (7,imuxyz) = soxyzp(7,imuxyz) +
     *              (soyzp(imuyz,7)*sdse+soyzep(imuyz,3)*sds)
  446             continue
c
c                                        muy impar muz par Campo (3)
c
                imuyz = 0
                do 447 muy=2,mumay,2
                  do 447 muz=1,mumaz,2
                    imuyz  = imuyz  + 1
                    ipuxyz = ipuxyz + 1
                    soxyzn  (6,ipuxyz) = soxyzn(6,ipuxyz) +
     *              (soyzn(imuyz,6)*sdse+soyzen(imuyz,2)*sds)
                    soxyzp  (6,ipuxyz) = soxyzp(6,ipuxyz) +
     *              (soyzp(imuyz,6)*sdse+soyzep(imuyz,2)*sds)
  447             continue
  445           continue
c                                       +----------------+
c --------------------------------------| nx+nx' impar   |
c                                       +----------------+
            else if(.not.lpx) then
c                                                         Campo (1)
              imuxyz = 0
              ipuxyz = 0
              do 448 mux =2,mumax,2
                indx = indx + 1
                imuyz = 0
                aa1  = ajx1d(1,indx)
                aa2  = ajx1d(2,indx)
                bb   = 0.0d+00
                if (mux.ge.muix.and.mux.le.mufx) then
                    indt = indt + 1
                    bb = tz1d(indt)
                end if
c
                do 449 muy=1,mumay,2
                  do 449 muz=2,mumaz,2
                    imuyz  = imuyz  + 1
                    imuxyz = imuxyz + 1
         soxyzn(3,imuxyz)=soxyzn(3,imuxyz)+soyzn  (imuyz,3)*bb
         sxyz3n(1,imuxyz)=sxyz3n(1,imuxyz)+scryz3n(1,imuyz)*aa1
         sxyz3n(2,imuxyz)=sxyz3n(2,imuxyz)+scryz3n(2,imuyz)*aa2
         soxyzp(3,imuxyz)=soxyzp(3,imuxyz)+soyzp  (imuyz,3)*bb
         sxyz3p(1,imuxyz)=sxyz3p(1,imuxyz)+scryz3p(1,imuyz)*aa1
         sxyz3p(2,imuxyz)=sxyz3p(2,imuxyz)+scryz3p(2,imuyz)*aa2
  449             continue
c
c                                        muy impar muz par Campo (3)
c
                imuyz = 0
                do 450 muy=2,mumay,2
                  do 450 muz=1,mumaz,2
                    imuyz  = imuyz  + 1
                    ipuxyz = ipuxyz + 1
                    soxyzn  (2,ipuxyz) = soxyzn(2,ipuxyz) +
     *              fasx1*soyzn(imuyz,2)*bb
                    sxyz4n(1,ipuxyz)=sxyz4n(1,ipuxyz)+
     *               fasx1*scryz4n(1,imuyz)*aa1
                    sxyz4n(2,ipuxyz)=sxyz4n(2,ipuxyz)+
     *               fasx1*scryz4n(2,imuyz)*aa2
                    soxyzp  (2,ipuxyz) = soxyzp(2,ipuxyz) +
     *              fasx1*soyzp(imuyz,2)*bb
                    sxyz4p(1,ipuxyz)=sxyz4p(1,ipuxyz)+
     *               fasx1*scryz4p(1,imuyz)*aa1
                    sxyz4p(2,ipuxyz)=sxyz4p(2,ipuxyz)+
     *               fasx1*scryz4p(2,imuyz)*aa2
  450             continue
  448           continue
c                                                    Campo (0)
              imuxyz = 0
              do 451 mux =1,mumax,2
                sds   = dsomu(indds)
                indds = indds + 1
                sdse  = dsomu(inddse)
                inddse = inddse + 1
c
                imuyz = 0
                  do 451 muy=1,mumay,2
                    do 451 muz=1,mumaz,2
                    imuyz  = imuyz + 1
                    imuxyz = imuxyz + 1
                    soxyzn  (10,imuxyz) = soxyzn(10,imuxyz) +
     *              (soyzn(imuyz,10)*sdse+soyzen(imuyz,6)*sds)
                    soxyzn  (9,imuxyz) = soxyzn(9,imuxyz) +
     *              (soyzn(imuyz, 9)*sdse+soyzen(imuyz,5)*sds)*fasx
                    soxyzp  (10,imuxyz) = soxyzp(10,imuxyz) +
     *              (soyzp(imuyz,10)*sdse+soyzep(imuyz,6)*sds)
                    soxyzp  (9,imuxyz) = soxyzp(9,imuxyz) +
     *              (soyzp(imuyz, 9)*sdse+soyzep(imuyz,5)*sds)*fasx
  451               continue
c ------------------------------------------------------------------
            end if
c ------------------------------------------------------------------
   11   continue
    1 continue
      return
      end
c+---------------------------------------------------------------------+
c|    Previous to the calculation of the Coulomb part                  |
c+---------------------------------------------------------------------+
      Subroutine predirco(sxyzc,scf,coum,acou)
      Implicit real*8 (A-H,O-Z)
      Dimension sxyzc(ndmu),scf(ndmu)
      Dimension acou(nacou)
      dimension coum(nmax,nmax)
C
      Common/DIMEN/nmax,nmax1,nxmax,nymax,nzmax,nwf2,
     * maxtz1,maxjz1,maxjx1,maxjy1,maxlso,maxdso,
     * ndmu,ndthet,nacou,nrop,nrom,nrop8,nrom8,nherm,nherm38,nlegn
c
      Common/OSCLEN/bx,by,bz
      Common/pulga/ideb
c
      !fac = 1.0d+00/(bx*by*bz*2.0d+00*dsqrt(2.0d+00))
      fac = 0.0d+00
c
      mumax = 2*nxmax-1
      mumay = 2*nymax-1
      mumaz = 2*nzmax-1
c
      mumax2 = mumay*mumaz
      indmu=0
      do imux=1,mumax,2
         jmux=(imux+1)/2
         do imuy=1,mumay,2
            jmuy = (imuy+1)/2
            do imuz=1,mumaz,2
               jmuz = (imuz+1)/2
               indmu = indmu+1
c
               indnu=0
               sum  = 0.0d+00
               do inux=1,mumax,2
                  jnux = (inux+1)/2
                  iacx = (imux+inux)/2
                  jacx = (iacx-1)*mumax2
                  cox = coum(jmux,jnux)
                  do inuy=1,mumay,2
                     jnuy = (inuy+1)/2
                     iacy = (imuy+inuy)/2
                     jacy = jacx + (iacy-1)*mumaz
                     coy  = coum(jmuy,jnuy)
                     do inuz=1,mumaz,2
                        jnuz= (inuz+1)/2
                        iacz= (imuz+inuz)/2
                        jacz= jacy + iacz
                        coz = coum(jmuz,jnuz)
                        indnu=indnu+1
                        sum = sum + cox*coy*coz*acou(jacz)*sxyzc(indnu)
                     end do                                  ! inuz
                  end do                                     ! inuy
               end do                                        ! inux
               sum = sum * fac
               scf (indmu) = sum
            end do                                           ! imuz
         end do                                              ! imuy
      end do                                                 ! imux
      return
      end
c+---------------------------------------------------------------------+
c|  Last change: February 24 1993                                      |
c+=====================================================================+
c|                                                                     |
c|  Some quantities required in the calculation of the direct field    |
c|  associated to the density dependent part of the interaction are    |
c|  computed here.                                                     |
c|                                                                     |
c|  Rortmp, dmus0 and dmus1 are temporary vectors                      |
c|                                                                     |
c|                                                                     |
c|                                                                     |
c+---------------------------------------------------------------------+
      subroutine predirdd(wf,xherm,xeherm,ror,rortmp,dmus0,dmus1,dmu)
c
      Implicit real*8(a-h,o-z)
      Implicit logical (l)
c
      Dimension ror(Nherm38),wf(nwf2,Nherm)
      Dimension xherm(Nherm),xeherm(Nherm)
      Dimension DMU(15,ndmu),sumv(15)
      Dimension rortmp(irortmp)
      Dimension dmus0(15,idmus0)
      Dimension dmus1(15,idmus1)
c
      Common/DIMPDD/ irortmp,idmus0,jdmus0,idmus1,jdmus1,irorz,irory
c
      Common/DIMEN/nmax,nmax1,nxmax,nymax,nzmax,nwf2,
     * maxtz1,maxjz1,maxjx1,maxjy1,maxlso,maxdso,
     * ndmu,ndthet,nacou,nrop,nrom,nrop8,nrom8,nherm,nherm38,nlegn
c
      Common/OSCLEN/ bx,by,bz
      Common/GOGINT/amu(2),xw(2),xh(2),xb(2),xm(2),WLS,t3,alpha,x0,e2
      Common/GOGHAN/w2b(2),h2m(2),x01,x02
      Common /Const/ pi,sqpi,dln2,sq2
c
      Common /mfval/ ehfp,ehfn,epaip,epain,errp,errn,ecech,ojxp,ojxn,
     *        op,on,ojz2p,ojz2n,oq20p,oq20n,oq22p,oq22n,or2p,or2n,
     *        osxp,osxn
c
      third  = 1.0d+00/3.0d+00
      fthird = 4.0d+00/3.0d+00
      nherm3 = nherm**3

c+---------------------------------------------------------------------+
c                alfa           1/3
c  Computing   ro      and   ro
c                              Protons
c+---------------------------------------------------------------------+
      nh1 = nherm3
       ih = 0
              do nhx=1,nherm
                 do nhy=1,nherm
                    do nhz=1,nherm
                       ih = ih + 1
                       ih1 = ih+nh1
                       rop = ror(ih)
                       ron = ror(ih+nh1)
                       if(rop.lt.0.0d+00) rop=1.0d-20
                       rotot = rop + ron
                       if(rotot.lt.0.0d+00) rotot = 1.0d-20
                       rortmp(ih    ) = rotot**alpha
                       rortmp(ih+nh1) = rop**third
                    end do
                 end do
              end do
c
c     Optimization: The z part
c
      nh1 = nherm3
      nh2 = nh1 + nherm3
      nh3 = nh2 + nherm3
      nh4 = nh3 + nherm3
      nh5 = nh4 + nherm3
      nh6 = nh5 + nherm3
      nh7 = nh6 + nherm3
c
      imumaz = 2*nzmax-1
           is0 = 0
           do imuz=1,imumaz,2
              ih  = 0
              do nhx=1,nherm
                 do nhy=1,nherm
                    do id=1,15
                      sumv(id) = 0.0d+00
                    end do    ! id
                    do nhz=1,nherm
                       ih = ih + 1
                       ih1 = ih+nh1
                       ih2 = ih+nh2
                       ih3 = ih+nh3
                       ih4 = ih+nh4
                       ih5 = ih+nh5
                       ih6 = ih+nh6
                       ih7 = ih+nh7
                       w = wf(imuz,nhz)
                       rotot = ror(ih    ) + ror(ih+nh1)
                       if(rotot.le.0.0d+00) rotot = 1.d-20
                       roalf   = rortmp(ih)
                       rothird = rortmp(ih+nh1)
                       roalf1= roalf/rotot
                       sumv(1) = sumv(1) + w*(roalf*ror(ih))
                       sumv(2) = sumv(2) + w*(roalf*ror(ih1))
                       sumv(3) = sumv(3) + w*(roalf*ror(ih4))
                       sumv(4) = sumv(4) + w*(roalf*ror(ih5))
                       sumv(7) = sumv(7) + w*(roalf*ror(ih6))
                       sumv(8) = sumv(8) + w*(roalf*ror(ih7))
                       sumv(9) = sumv(9) + w*(roalf1*ror(ih)*ror(ih))
                       sumv(10)= sumv(10)+ w*(roalf1*ror(ih1)*ror(ih))
                       sumv(11)= sumv(11)+ w*(roalf1*ror(ih1)*ror(ih1))
                       sumv(12)= sumv(12)+ w*roalf1*
     *   (ror(ih2)*ror(ih2)+ror(ih4)*ror(ih4)+ror(ih6)*ror(ih6))
                       sumv(13)= sumv(13)+ w*roalf1*
     *   (ror(ih2)*ror(ih3)+ror(ih4)*ror(ih5)+ror(ih6)*ror(ih7))
                       sumv(14)= sumv(14)+ w*roalf1*
     *   (ror(ih3)*ror(ih3)+ror(ih5)*ror(ih5)+ror(ih7)*ror(ih7))
                       sumv(15)= sumv(15)+ w*rothird
                    end do    ! nhz
                 is0 = is0 + 1
                 do id=1,15
                   dmus0(id,is0) = sumv(id)
                 end do       ! id
                 end do       ! nhy
              end do          ! nhx
           end do             ! imuz
c
           is0 = 0
           do imuz=2,imumaz,2
              ih  = 0
              do nhx=1,nherm
                 do nhy=1,nherm
                      sumv(5) = 0.0d+00
                      sumv(6) = 0.0d+00
                    do nhz=1,nherm
                       ih = ih + 1
                       ih1 = ih+nh1
                       ih2 = ih+nh2
                       ih3 = ih+nh3
                       w = wf(imuz,nhz)
                       roalf   = rortmp(ih)
                       sumv(5) = sumv(5) + w*(roalf*ror(ih2))
                       sumv(6) = sumv(6) + w*(roalf*ror(ih3))
                    end do    ! nhz
                 is0 = is0 + 1
                   dmus0(5,is0)  = sumv(5)
                   dmus0(6,is0)  = sumv(6)
                 end do       ! nhy
              end do          ! nhx
           end do             ! imuz
c
c      Optimization. The y part
c
      imumay = 2*nymax-1
        is1 = 0
        do imuy=1,imumay,2
           is0 = 0
           do imuz=1,imumaz,2
              do nhx=1,nherm
                 do id=1,15
                     sumv(id) = 0.0d+00
                 end do
                 do nhy=1,nherm
                       wy = wf(imuy,nhy)
                       w  = wy
                       is0 = is0 + 1
                       sumv(1) = sumv(1) + dmus0(1 ,is0)*w
                       sumv(2) = sumv(2) + dmus0(2 ,is0)*w
                       sumv(3) = sumv(3) + dmus0(3 ,is0)*w
                       sumv(4) = sumv(4) + dmus0(4 ,is0)*w
                       sumv(9) = sumv(9) + dmus0(9 ,is0)*w
                       sumv(10)= sumv(10)+ dmus0(10,is0)*w
                       sumv(11)= sumv(11)+ dmus0(11,is0)*w
                       sumv(12)= sumv(12)+ dmus0(12,is0)*w
                       sumv(13)= sumv(13)+ dmus0(13,is0)*w
                       sumv(14)= sumv(14)+ dmus0(14,is0)*w
                       sumv(15)= sumv(15)+ dmus0(15,is0)*w
                 end do   ! nhy
                 is1 = is1 + 1
                 do id=1,15
                   dmus1(id,is1) = sumv(id)
                 end do
              end do      ! nhx
c
           end do         ! imuz
        end do            ! imuy
c
      is11 = 0
        do imuy=1,imumay,2  ! (1)
           is01 = 0
           do imuz=2,imumaz,2
              do nhx=1,nherm
                 sumv(5) = 0.0d+00
                 sumv(6) = 0.0d+00
                 do nhy=1,nherm
                 wy = wf(imuy,nhy)
                    w  = wy
                    is01 = is01 + 1
                       sumv(5) = sumv(5) + dmus0(5,is01)*w
                       sumv(6) = sumv(6) + dmus0(6,is01)*w
                 end do    ! nhy
                 is11 = is11 + 1
                 dmus1(5,is11) = sumv(5)
                 dmus1(6,is11) = sumv(6)
              end do       ! nhx
           end do          ! imuz
        end do             ! imuy
c
      is13 = 0
        do imuy=2,imumay,2  ! (3)
           is02 = 0
           do imuz=1,imumaz,2
              do nhx=1,nherm
              sumv(7) = 0.0d+00
              sumv(8) = 0.0d+00
                 do nhy=1,nherm
                 wy = wf(imuy,nhy)
                       w = wy
                       is02 = is02 + 1
                       sumv(7) = sumv(7) + dmus0(7,is02)*w
                       sumv(8) = sumv(8) + dmus0(8,is02)*w
                 end do   ! nhy
                 is13 = is13 + 1
                 dmus1(7,is13) = sumv(7)
                 dmus1(8,is13) = sumv(8)
              end do      ! nhx
           end do         ! imuz
        end do            ! imuy
c
c
c     Final  optimization
c
      imumax = 2*nxmax-1
      imu = 0
      do 501 imux=1,imumax,2    ! (0) and (2)
        is1 = 0
        do 502 imuy=1,imumay,2  ! (0)
           do 503 imuz=1,imumaz,2
              imu= imu+1
              do 700 id=1,15
  700            sumv(id) = 0.0d+00
              do 504 nhx=1,nherm
                 wx = wf(imux,nhx)
                       w  = wx
                       is1 = is1 + 1
                       sumv(1) = sumv(1) + dmus1(1 ,is1)*w
                       sumv(2) = sumv(2) + dmus1(2 ,is1)*w
                       sumv(3) = sumv(3) + dmus1(3 ,is1)*w
                       sumv(4) = sumv(4) + dmus1(4 ,is1)*w
                       sumv(9) = sumv(9) + dmus1(9 ,is1)*w
                       sumv(10)= sumv(10)+ dmus1(10,is1)*w
                       sumv(11)= sumv(11)+ dmus1(11,is1)*w
                       sumv(12)= sumv(12)+ dmus1(12,is1)*w
                       sumv(13)= sumv(13)+ dmus1(13,is1)*w
                       sumv(14)= sumv(14)+ dmus1(14,is1)*w
                       sumv(15)= sumv(15)+ dmus1(15,is1)*w
  504         continue
c
              do id=1,15
                 DMU(id,imu) = sumv(id)
              end do
c
  503      continue
  502   continue
  501 continue
c
      imu1 = 0
      imu3 = 0
      do 601 imux=2,imumax,2    ! (1) and (3)
        is11 = 0
        do 602 imuy=1,imumay,2  ! (1)
           do 603 imuz=2,imumaz,2
              sumv(5) = 0.0d+00
              sumv(6) = 0.0d+00
              imu1= imu1+1
              do 604 nhx=1,nherm
                 wx = wf(imux,nhx)
                    w  = wx
                    is11 = is11 + 1
                       sumv(5) = sumv(5) + dmus1(5,is11)*w
                       sumv(6) = sumv(6) + dmus1(6,is11)*w
  604         continue
           DMU(5,imu1) = sumv(5)
           DMU(6,imu1) = sumv(6)
  603      continue
  602   continue
c
        is12 = 0
        do 607 imuy=2,imumay,2  ! (3)
           do 608 imuz=1,imumaz,2
              sumv(7) = 0.0d+00
              sumv(8) = 0.0d+00
              imu3= imu3+1
              do 609 nhx=1,nherm
                 wx = wf(imux,nhx)
                       w = wx
                       is12 = is12 + 1
                       sumv(7) = sumv(7) + dmus1(7,is12)*w
                       sumv(8) = sumv(8) + dmus1(8,is12)*w
  609         continue
           DMU(7,imu3) = sumv(7)
           DMU(8,imu3) = sumv(8)
  608      continue
  607   continue
c
  601 continue

c
              ecch   = 0.0 d+00
              erreap = 0.0 d+00
              errean = 0.0 d+00
              q20p   = 0.0 d+00
              q20n   = 0.0 d+00
              q22p   = 0.0 d+00
              q22n   = 0.0 d+00
              r2p    = 0.0 d+00
              r2n    = 0.0 d+00
              ih  = 0
              do nhx=1,nherm
                 wx = xeherm(nhx)
                 x2 = (xherm(nhx)*bx)**2
                 do nhy=1,nherm
                    wy = xeherm(nhy)
                    y2 = (xherm(nhy)*by)**2
                    do nhz=1,nherm
                       z2 = (xherm(nhz)*bz)**2
                       ih = ih + 1
                       ih1 = ih+nh1
                       ih2 = ih+nh2
                       ih3 = ih+nh3
                       ih4 = ih+nh4
                       ih5 = ih+nh5
                       ih6 = ih+nh6
                       ih7 = ih+nh7
                       w = wx*wy*xeherm(nhz)
                       rotot = ror(ih    ) + ror(ih+nh1)
                       if(rotot.lt.0.0d+00) rotot = 1.0d-20
                       roalf = rortmp(ih)
                       roalf1= roalf/rotot
                       rop     = ror(ih)
                       ron     = ror(ih1)
                       ropft   = rop*rortmp(ih1)
                       q20p    = q20p + w*rop*(2.*z2-x2-y2)
                       q20n    = q20n + w*ron*(2.*z2-x2-y2)
                       q22p    = q22p + w*rop*(x2-y2)
                       q22n    = q22n + w*ron*(x2-y2)
                       r2p     = r2p  + w*rop*(x2+y2+z2)
                       r2n     = r2n  + w*ron*(x2+y2+z2)
                       xxpp    =
     *   (ror(ih2)*ror(ih2)+ror(ih4)*ror(ih4)+ror(ih6)*ror(ih6))
                       xxpn    =
     *   (ror(ih2)*ror(ih3)+ror(ih4)*ror(ih5)+ror(ih6)*ror(ih7))
                       xxnn    =
     *   (ror(ih3)*ror(ih3)+ror(ih5)*ror(ih5)+ror(ih7)*ror(ih7))
                       ecch    = ecch    + w*ropft
                       erreap  = erreap  + w*rop*roalf1*
     *  (x01*(rop*rop+2.0d+00*rop*ron+ron*ron)-x02*(rop*rop+ron*ron)
     *  +0.5d+00*(x0*(xxpp+2.0*xxpn+xxnn)-(xxpp+xxnn)))
                       errean  = errean  + w*ron*roalf1*
     *  (x01*(rop*rop+2.0d+00*rop*ron+ron*ron)-x02*(rop*rop+ron*ron)
     *  +0.5d+00*(x0*(xxpp+2.0*xxpn+xxnn)-(xxpp+xxnn)))
                    end do
                 end do
              end do
      fint  = 8.0d+00*bx*by*bz
      ecech = ecch*fint*0.25d+00*e2*(3.0d+00/pi)**third
      errp  = erreap*fint*0.25d+00*alpha*t3
      errn  = errean*fint*0.25d+00*alpha*t3
      oq20p = q20p*fint
      oq20n = q20n*fint
      oq22p = q22p*fint
      oq22n = q22n*fint
      or2p  = r2p *fint
      or2n  = r2n *fint
c
      return
      end
c+---------------------------------------------------------------------+
c| Initial calculations (the ones needed at the very begining)         |
c|                                                                     |
c|   TZ1D, JX1D, XTJ, XTJP, IXTJ, etc                                  |
c+---------------------------------------------------------------------+
      Subroutine INIT(tz1d,ajx1d,ajy1d,ajz1d,xtjx,xtjxp,xtjy,xtjz,
     *wfhonec,xhnec,
     *wfho,xherm,xeherm,weher,wf,xleg,wleg,tz2d,al1d,dsomuv,acou,coum,
     *work,itz1d,itz2d,ijx1d,ijy1d,ijz1d,il1d,idsomu,ixtj,ixtj2,
     *iytj,iytj2,iztj,iztj2)
      Implicit Real*8 (A-H,O-Z)
      Logical LMZE,LMZO
c ---------------------------------------------------<< Start includes
      Include 'COMDIM'
c ---------------------------------------------------<< End   includes
      Dimension tz1d(maxtz1)
      Dimension ajx1d(2,maxjx1),ajy1d(2,maxjy1),ajz1d(2,maxjz1)
      Dimension itz1d(nmax,nmax)
      Dimension ijx1d(nxmax,nxmax),ijy1d(nymax,nymax),ijz1d(nzmax,nzmax)
      Dimension xtjx(2,maxtjx),xtjxp(2,maxtjx)
      Dimension xtjy(2,maxtjy),xtjz(2,maxtjz)
      Dimension ixtj(nxmax,nxmax),ixtj2(nxmax,nxmax)
      Dimension iytj(nymax,nymax),iytj2(nymax,nymax)
      Dimension iztj(nzmax,nzmax),iztj2(nzmax,nzmax)
c
      Dimension wfho(nmax,Nherm),xherm(nherm),weher(nherm)
      Dimension wfhonec(nmax,Nherm),xhnec(nherm)
      Dimension xeherm(nherm)
      Dimension wf(nwf2,Nherm)
c
      Dimension Work(1000)
      Dimension xleg(Nlegen),wleg(Nlegen)
c
      Dimension Tz2D(maxtz2),itz2d(nmax1,nmax1)
      Dimension AL1D(maxlso),il1d(nmax1,nmax1)
      Dimension DSOMUv(maxdso),IDSOMU(nmax,nmax)
c
      Dimension ACOU(Nacou),COUM(nmax,nmax)
c
      Common/DIMEN/nmax,nmax1,nxmax,nymax,nzmax,nwf2,
     * maxtz1,maxjz1,maxjx1,maxjy1,maxlso,maxdso,
     * ndmu,ndthet,nacou,nrop,nrom,nrop8,nrom8,nherm,nherm38,nlegn
c
      common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
c
      Common /NLEG/ Nlegen,maxtz2
      Common/GOGINT/amu(2),xw(2),xh(2),xb(2),xm(2),WLS,t3,alpha,x0,e2
      Common/GOGHAN/w2b(2),h2m(2),x01,x02
      Common/OSCLEN/bx,by,bz
      Common/NECKCO/ aneck,rneckp,rneckn
      Common /Const/ pi,sqpi,dln2,sq2
      Common /pulga/ideb
      
      call setfact()
c+--------------------------------------------------------------------+
c|    Gauss-Legendre and Gauss-Hermite points and weigths.            |
c|    Only the nodes greater than zero are stored                     |
c+--------------------------------------------------------------------+
      Nlega= Nlegen
      Nher = Nherm
      al   = 0.0d+00
      be   = 0.0d+00
      kpts = 0
      N    = 2*Nlega
      N1   = 3
      N2   = N1 + N
      N3   = N2 + N
c+--------------------------------------------------------------------+
c|---------------------------->> Gauss-Legendre <<--------------------|
c+--------------------------------------------------------------------+
      Kind = 1
      Call Gaussq(Kind,N,al,be,kpts,Work(1),Work(N1),Work(N2),Work(N3))
      do I=1,Nlega
           J1 = N2 + Nlega+ i -1
           J2 = N3 + Nlega+ i -1
           xleg(i) = Work(j1)
           wleg(i) = Work(j2)
      end do
c+--------------------------------------------------------------------+
c|---------------------------->> Gauss-Hermite <<---------------------|
c+--------------------------------------------------------------------+
      N    = 2*Nher
      N1   = 3
      N2   = N1 + N
      N3   = N2 + N
      Kind = 4
      Call Gaussq(Kind,N,al,be,kpts,Work(1),Work(N1),Work(N2),Work(N3))
      do I=1,Nher
           J1 = N2 + Nher + i -1
           J2 = N3 + Nher + i -1
           xherm(i) = Work(j1)
           weher(i) = Work(j2)
           xeherm(i)= Dexp(xherm(i)*xherm(i))*weher(i)
      end do
c+--------------------------------------------------------------------+
c|          Constants                                                 |
c+--------------------------------------------------------------------+
      pi    = 3.1415926535897932384 d+00
      sqpi  = dsqrt(pi)
      dln2  = dlog(2.0d+00)
      sq2   = dsqrt(2.0d+00)
c+--------------------------------------------------------------------+
c|  Brink Boeker stuff                                                |
c+--------------------------------------------------------------------+
      call TZ1DIM(nmax,maxtz1,tz1d,itz1d)
c
      call J1BB
     * (nxmax,maxjx1,maxtz1,tz1d,itz1d,nmax,work,amu,bx,ajx1d,ijx1d)
      call J1BB
     * (nymax,maxjy1,maxtz1,tz1d,itz1d,nmax,work,amu,by,ajy1d,ijy1d)
      call J1BB
     * (nzmax,maxjz1,maxtz1,tz1d,itz1d,nmax,work,amu,bz,ajz1d,ijz1d)
c+--------------------------------------------------------------------+
c|  Spin orbit stuff                                                  |
c+-----------------------------------------------  nmax1 = nmax + 1 --+
      call TZ1DIM(nmax1,maxtz2,tz2d,itz2d)
      Call LSOMU(nmax1,maxlso,maxtz2,tz2d,itz2d,work,al1d,il1d)
      Call DSOMU(nmax,nmax1,maxlso,maxdso,al1d,il1d,DSOmuv,iDSOmu)
c
      Nxmu = 2*nxmax - 1
      Nymu = 2*nymax - 1
      Nzmu = 2*nzmax - 1
c
      nll0 = 1
      nll1 = nlegen+1
      nll2 = 2*nlegen+1
      nll3 = 3*nlegen+1
      nll4 = 4*nlegen+1
      nll5 = 5*nlegen+1
      macou = nacou
      Call COULIN(macou,Nxmu,Nymu,Nzmu,xleg,wleg,work(nll0),
     *work(nll1),work(nll2),work(nll3),work(nll4),work(nll5),Acou)
c
      Call Coulm(coum,nmax)
c+--------------------------------------------------------------------+
c|                                                        Exchange    |
c+--------------------------------------------------------------------+
      call Sxtj(NXMAX,nmax,MAXTZ1,TZ1D,ITZ1D,MAXJX1,AJX1D,IJX1D,
     x                MAXTJX,XTJX)
c Tilde X term
      call Sxtjp(NXMAX,nmax,MAXTZ1,TZ1D,ITZ1D,MAXJX1,AJX1D,IJX1D,
     x                MAXTJX,XTJXP)
c
      call Sxtj(NYMAX,nmax,MAXTZ1,TZ1D,ITZ1D,MAXJY1,AJY1D,IJY1D,
     x                MAXTJY,XTJY)
c
      call Sxtj(NZMAX,nmax,MAXTZ1,TZ1D,ITZ1D,MAXJZ1,AJZ1D,IJZ1D,
     x                MAXTJZ,XTJZ)
c
      call sixtj(nxmax,ixtj,ixtj2)
      call sixtj(nymax,iytj,iytj2)
      call sixtj(nzmax,iztj,iztj2)
c+--------------------------------------------------------------------+
c|   HO wave functions in coordinate representation (Densities)       |
c+--------------------------------------------------------------------+
      iop = 1
      call howf1(wfho,xherm,nmax,Nherm,iop)
c+--------------------------------------------------------------------+
c|   HO wave functions (no exp term) (PREDIRDD)                       |
c+--------------------------------------------------------------------+
      iop = 0                                   ! no exponential term
      call howf1 (wf,xherm,nwf2,Nherm,iop)
      do ih=1,nherm
         do imu=1,nwf2
            wf(imu,ih) = wf(imu,ih)*weher(ih)
         end do    ! imu
      end do       ! ih
c+--------------------------------------------------------------------+
c|   HO wave functions (no exp term) (necking)                        |
c+--------------------------------------------------------------------+
      ff = 1.0d+00/dsqrt(1.d+00 + (bz/aneck)**2)
      do ih=1,nherm
         xhnec(ih) = ff*xherm(ih)
      end do
      iop = 0                                   ! no exponential term
      call howf1 (wfhonec,xhnec,nmax,Nherm,iop)
      do ih=1,nherm
         do imu=1,nmax
            wfhonec(imu,ih) = wfhonec(imu,ih)*Dsqrt(weher(ih))
         end do    ! imu
      end do       ! ih
c
      return
      end
c+--------------------------------------------------------------------+
c|  D I M E N S :         To compute dimensions and handy vectors.    |
c+--------------------------------------------------------------------+
c| Modified in version 4.1                                            |
c+--------------------------------------------------------------------+
      Subroutine DIMENS(jxmax,jymax,jzmax,iherm,ilegn,imy,imz)
      Implicit real*8 (a-h,o-z)
      Logical lx12,lxy12
      Logical LMZE,LMZO
      
      Include 'COMDIM'
      Include 'DIMTRIAX'
c
      Dimension imy(ixmax),imz(ixmax,iymax)
      Dimension INQE (ixmax ,iymax ,izmax ),INQO (ixmax ,iymax ,izmax )
c
      Common/UK2BC/ IEXP1(NP),IEXM1(NP),IEYP1(NP),IEYM1(NP),
     *IEZP1(NP),IEZM1(NP),IOXP1(NM),IOXM1(NM),IOYP1(NM),IOYM1(NM),
     *IOZP1(NM),IOZM1(NM)
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      Common/DIMPDD/ irortmp,idmus0,jdmus0,idmus1,jdmus1,irorz,irory
c
      Common/DIMEN/kmax,kmax1,kxmax,kymax,kzmax,nwf2,
     * maxtz1,maxjz1,maxjx1,maxjy1,maxlso,maxdso,
     * ndmu,ndtheta,nacou,nrop,nrom,nrop8,nrom8,nherm,nherm38,nlegn
c
      common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
c
      Common /NLEG/ Nlegen,maxtz2
      Common/ILOC/MP,MM,ILOCE(ixmax,iymax),ILOCO(ixmax,iymax)
      Common/ILOCECH/IXGP(ixmax,ixmax),IXGM(ixmax,ixmax)
      Common/basin/IQNE(3,inqne),IQNO(3,inqno)
      Common/PULGA/ideb
      Common/infob/aq0inf,ap0inf,arn0inf,nixmax,niymax,nizmax
c
      nxmax = jxmax
      nymax = jymax
      nzmax = jzmax
      nixmax = jxmax
      niymax = jymax
      nizmax = jzmax
      nmax  = Max(nxmax,nymax,nzmax)
      nmax1 = nmax+1
      kmax  = nmax
      kmax1 = nmax1
      kxmax = nxmax
      kymax = nymax
      kzmax = nzmax
c
      nherm = iherm
      nlegn = ilegn
      nlegen= ilegn
c
      do i=1,nxmax
         my(i) = imy(i)
      end do
      
      do j=1,nymax
         do i=1,nxmax
           mz(i,j)=imz(i,j)
         end do
      end do
c
      aq0inf = q0inf
      ap0inf = p0inf
      arn0inf = rn0inf
c
      write(6,100) q0inf,p0inf,rn0inf,nxmax,nymax,nzmax
      write(6,101) (my(ii),ii=1,nxmax)
      do i=1,nxmax
         write(6,102) (mz(i,jj),jj=1,nymax)
      end do
  100 format( ///,15x,' Oscillator Basis parameters ',/,15x,29('-'),/,
     *  15x,' q= ',f7.4,' p= ',f7.4,' N0= ',f8.4,/,
     *  15x,' Nxmax= ',i3,' Nymax ',i3,' Nzmax ',i3)
  101 format( /,15x,' My(nx)    ',20i3,/)
  102 format( /,15x,' MZ(NX,NY) ',20i3)
c
c
c    Computing  LMZE LMZO NZIE NZIO
c
      do nx=1,nxmax
        nym = my(nx)
          do ny=1,nym
             nxy  = nx +ny
             nzie(nx,ny) = 1 + Mod(nxy,2)
             nzio(nx,ny) = 2 - Mod(nxy,2)
             nzm = mz(nx,ny)
             lmze(nx,ny) = .false.
             lmzo(nx,ny) = .false.
             if(nzm.ge.nzie(nx,ny)) lmze(nx,ny) = .true.
             if(nzm.ge.nzio(nx,ny)) lmzo(nx,ny) = .true.
          end do
      end do
c
c--------------------------------------> Setting inqe,o to zero
      do nx=1,nxmax
         do ny=1,nymax
            do nz=1,nzmax
               inqe(nx,ny,nz) = 0
               inqo(nx,ny,nz) = 0
            end do
         end do
      end do
c
      iz1e = 0
      iz1o = 0
      jz1e = 0
      jz1o = 0
      do nx1=1,nxmax
        nym1 = my(nx1)
          do ny1 = 1,nym1
            ILOCE(nx1,ny1) = iz1e
            ILOCO(nx1,ny1) = iz1o
            nxy1 = nx1+ny1
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = mz(nx1,ny1)
            if(lmze(nx1,ny1)) then
              iz1e=iz1e+(nzm1-nzi1e)/2+1  ! pos parity
              do nz1=nzi1e,nzm1,2
                 jz1e = jz1e + 1
                 iqne(1,jz1e) = nx1
                 iqne(2,jz1e) = ny1
                 iqne(3,jz1e) = nz1
                 inqe(nx1,ny1,nz1) = jz1e
              end do
            end if
            if(lmzo(nx1,ny1)) then
              iz1o=iz1o+(nzm1-nzi1o)/2+1  ! neg parity
              do nz1=nzi1o,nzm1,2
                 jz1o = jz1o + 1
                 iqno(1,jz1o) = nx1
                 iqno(2,jz1o) = ny1
                 iqno(3,jz1o) = nz1
                 inqo(nx1,ny1,nz1) = jz1o
              end do
            end if
          end do ! ny1
      end do     ! nx1
c
      MP = iz1e     ! Positive parity matrices dimension NPxNP
      MM = iz1o     ! Negative parity matrices dimension NMxNM
c
      if((MP.ne.NP).or.(MM.ne.NM)) then
        write(6,*) '****** DIMENS: Np and/or Nm inconsistent'
        write(6,*) '****** Np = ',Mp,' Nm ',MM
        stop
      end if
c ---------------------------------------------
      do i=1,np
         nx= iqne(1,i)
         ny= iqne(2,i)
         nz= iqne(3,i)
c
         if((nx+1).le.nxmax) IEXP1(i) = inqo(nx+1,ny  ,nz  )
         if((nx).ne.1)       IEXM1(i) = inqo(nx-1,ny  ,nz  )
         if((ny+1).le.nymax) IEYP1(i) = inqo(nx  ,ny+1,nz  )
         if((ny).ne.1)       IEYM1(i) = inqo(nx  ,ny-1,nz  )
         if((nz+1).le.nzmax) IEZP1(i) = inqo(nx  ,ny  ,nz+1)
         if((nz).ne.1)       IEZM1(i) = inqo(nx  ,ny  ,nz-1)
      end do
c
c     Negative parity
c
      do i=1,nm
         nx= iqno(1,i)
         ny= iqno(2,i)
         nz= iqno(3,i)
c
         if((nx+1).le.nxmax) IOXP1(i) = inqe(nx+1,ny  ,nz  )
         if((nx).ne.1)       IOXM1(i) = inqe(nx-1,ny  ,nz  )
         if((ny+1).le.nymax) IOYP1(i) = inqe(nx  ,ny+1,nz  )
         if((ny).ne.1)       IOYM1(i) = inqe(nx  ,ny-1,nz  )
         if((nz+1).le.nzmax) IOZP1(i) = inqe(nx  ,ny  ,nz+1)
         if((nz).ne.1)       IOZM1(i) = inqe(nx  ,ny  ,nz-1)
      end do
c
c
c
      nrop = (np*(np+1))/2   !   Ro dimensions    Pos parity
      nrom = (nm*(nm+1))/2   !                    Neg parity
      nrop8 = nrop*8
      nrom8 = nrom*8
      ngp   = nrop
      ngm   = nrom
      ngp8  = nrop8
      ngm8  = nrom8
c
c      Define IXGP, IXGM needed in gechxyz
c
      jgp = 0
      jgm = 0
      do nx1=1,nxmax
        nym1 = my(nx1)
        do nx2=1,nx1
          nym2 = my(nx2)
          lx12 = nx1.eq.nx2
c
          ixgp(nx1,nx2) = jgp
          ixgm(nx1,nx2) = jgm
          ixgp(nx2,nx1) = jgp
          ixgm(nx2,nx1) = jgm
c
          do ny1 = 1,nym1
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = mz(nx1,ny1)
c
c   lgst1e    to control the appearence of ghost states (those not
c             appearing due to parity restrictions )
c
            if(lx12) nym2 = ny1
            do ny2=1,nym2
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = mz(nx2,ny2)
              lxy12 = lx12.and.(ny1.eq.ny2)
c
c---->    Positive parity
c
              if(lmze(nx1,ny1).and.lmze(nx2,ny2)) then
c
              do nz1 =nzi1e,nzm1,2
                if(lxy12) nzm2=nz1
                do nz2 = nzi2e,nzm2,2
                jgp = jgp + 1
                end do      ! nz2
              end do        ! nz1
c
              end if
c
c---->    Negative parity
c
c
              if(lmzo(nx1,ny1).and.lmzo(nx2,ny2)) then
c
              do nz1 =nzi1o,nzm1,2
                if(lxy12) nzm2=nz1
                do nz2 = nzi2o,nzm2,2
                jgm = jgm + 1
                end do     ! nz2
              end do       ! nz1
              end if
            end do         ! ny2
          end do           ! ny1
        end do             ! nx2
      end do               ! nx1
c+--------------------------------------------------------------------+
c|     Permanent storage                                              |
c|                                     Direct                         |
c+--------------------------------------------------------------------+
      Maxtz1 = (Nmax*(Nmax**2 + 3*Nmax + 2))/6
      Maxtz2 = ((Nmax+1)*((Nmax+1)**2 + 3*Nmax + 5))/6
      Maxjx1 = ((Nxmax+1)**2/4)*Nxmax + (Nxmax**2/4)*(Nxmax-1)
      Maxjy1 = ((Nymax+1)**2/4)*Nymax + (Nymax**2/4)*(Nymax-1)
      Maxjz1 = ((Nzmax+1)**2/4)*Nzmax + (Nzmax**2/4)*(Nzmax-1)
      Maxlso = ((Nmax+2)**2/4)*(Nmax+1) + ((Nmax+1)**2/4)*Nmax
      Maxdso = ((Nmax+1)**2/2-Nmax)*(Nmax+1)+(Nmax**2/2)*Nmax
      Nacou  = (2*Nxmax-1)*(2*Nymax-1)*(2*Nzmax-1)
      Ncoulm = Nmax**2
      nwf2   = 2*Nmax-1
      Nwavef = Nmax*Nherm
      Nwavefd= nwf2*Nherm
c+--------------------------------------------------------------------+
c|                                     Exchange (SXTJ)                |
c+--------------------------------------------------------------------+
      nxmax2 = nxmax/2
      nymax2 = nymax/2
      nzmax2 = nzmax/2
      maxtjx = ((nxmax2+1)*(nxmax-nxmax2))**2+
     *         ( (nxmax*(nxmax+1))/2-(nxmax2+1)*(nxmax-nxmax2) )**2
      maxtjy = ((nymax2+1)*(nymax-nymax2))**2+
     *         ( (nymax*(nymax+1))/2-(nymax2+1)*(nymax-nymax2) )**2
      maxtjz = ((nzmax2+1)*(nzmax-nzmax2))**2+
     *         ( (nzmax*(nzmax+1))/2-(nzmax2+1)*(nzmax-nzmax2) )**2
      maxx   = 4*maxtjx+2*maxtjy+2*maxtjz
c
      NPERTOT = Maxtz1 + Maxtz2 + 2*Maxjx1 + 2*Maxjy1 + 2*Maxjz1 +
     *          Maxlso + Maxdso + Nacou    + Ncoulm   + Nwavef   +
     *          Nwavefd+ 2*Nherm+  2*Nlegn + maxx
c+--------------------------------------------------------------------+
c|   Temporary                                                        |
c|                                                  Direct            |
c+--------------------------------------------------------------------+
      NDMU    = nxmax*nymax*nzmax
      NDTHETA = Ndmu*(nzmax/2+1)
      nherm38 = 8*nherm**3
      irorz  = Nherm*4
      irory  = Nherm*Nherm*6
      irortmp= 2*Nherm**3
      idmus0 = nherm**2 * nmax
      idmus1 = nherm    * nmax**2
      jdmus0 = 15*idmus0
      jdmus1 = 15*idmus1
c
      NDIRT   = 61 * Ndmu + 20 * NDTHETA + nherm38 + irortmp +
     *          jdmus0 + jdmus1
c+--------------------------------------------------------------------+
c|                                             Exchange and pairing   |
c+--------------------------------------------------------------------+
      mazopti= 16*nzmax*nzmax
      mayopti= 16*nzmax**2 * nymax**2
      NECHT  = mazopti + mayopti
      NTEMP  = MAX(NDIRT,NECHT)
c
      return
      end
cÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
c³    Computes the integrals needed to compute the Coulomb matrix    ³
c³    elements. The integrals are computed by means of a Gauss-      ³
c³    Legendre quadrature with 2*Nleg points.                        ³
cÃÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ´
c³                                                                   ³
c³    Ndim       Dimension of output vector ACOU                     ³
c³               Ndim = Nxmax*Nymax*Nzmax                            ³
c³    Nxmax      Maximum value of K in the x direction + 1           ³
c³    Nymax      Maximum value of K in the y direction + 1           ³
c³    Nzmax      Maximum value of K in the z direction + 1           ³
c³    bx by bz   Oscillator lengths                                  ³
c³    Acou       Output vector                                       ³
cÃÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ´
cÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
      SUBROUTINE COULIN(Ndim,Nxmax,Nymax,Nzmax,xleg,wleg,
     * aux,auy,auz,bux,buy,buz,Acou)
      Implicit Real*8 (A-H,O-Z)
      Dimension Acou(Ndim)
      Dimension aux(Nlegen),auy(Nlegen),auz(Nlegen)
      Dimension bux(Nlegen),buy(Nlegen),buz(Nlegen)
      Dimension xleg(Nlegen),wleg(Nlegen)
c
      Common/DIMEN/kmax,kmax1,kxmax,kymax,kzmax,nwf2,
     * maxtz1,maxjz1,maxjx1,maxjy1,maxlso,maxdso,
     * ndmu,ndthet,nacou,nrop,nrom,nrop8,nrom8,nherm,nherm38,nlegen
c
      Common /OSCLEN/bx,by,bz
      common/Const/ pi,sqpi,Dln2,sq2
      Common/PULGA/ideb
c
      if(ideb.gt.0) then
         write(6,*) ' COULIN ****  NACOU,NLEG ',Ndim,Nlegen
      end if
c
      bi  = Min(bx,by,bz)
      bex = (bi/bx)**2
      bey = (bi/by)**2
      bez = (bi/bz)**2
c
c
      cx  = (1.d+00-bex)/bex
      cy  = (1.d+00-bey)/bey
      cz  = (1.d+00-bez)/bez
      fac = 4.0d+00*bi*bi/sqpi
c
      do il=1,Nlegen
           x2leg   = xleg(il)**2
           dx      = 1.d+00/(1.0d+00+cx*x2leg)/bex
           dy      = 1.d+00/(1.0d+00+cy*x2leg)/bey
           dz      = 1.d+00/(1.0d+00+cz*x2leg)/bez
           bux(il) = x2leg*(dx)
           buy(il) = x2leg*(dy)
           buz(il) = x2leg*(dz)
           aux(il) = 1.0d+00/(x2leg**3 * dsqrt(dx*dy*dz))
      end do
c
      indx = 0
      do Nxl=1,Nxmax
          do il=1,Nlegen
             aux(il) = aux(il)*bux(il)
             auy(il) = aux(il)
          end do
          do Nyl=1,Nymax
              do il=1,Nlegen
                 auy(il) = auy(il)*buy(il)
                 auz(il) = auy(il)
              end do
              do Nzl=1,Nzmax
                  do il=1,Nlegen
                     auz(il) = auz(il)*buz(il)
                  end do
                  sum = 0.0d+00
                  do il=1,Nlegen
                     sum= sum + wleg(il)*auz(il)
                  end do
                  indx = indx + 1
                  acou(indx) = sum*fac
              end do     ! nzl
          end do         ! nyl
      end do             ! nxl
cdebug
c     write(6,*) ' COULIN  indx, nacou  ',indx,nacou
cdebug
      return
      end
      Subroutine coulm(coum,nmax)
      Implicit real*8 (a-h,o-z)
      Include 'COMDIM'
      Dimension coum(nmax,nmax)
      Common /fact/ dfact(Nfac),ddfact(Nfac),nfacm
      Common /Const/ pi,sqpi,dln2,sq2
      imumax = 2*nmax - 1
      inumax = imumax
      indmu = 0
      do imu=1,imumax,2
         indmu = indmu+1
         indnu = 0
         do inu=1,inumax,2
           fas = dfloat(1-2*Mod(Iabs(imu-inu)/2,2))
           indnu = indnu+1
           fac = ddfact((imu+inu)/2)-(0.5d+00*(dfact(imu)+dfact(inu))+
     *      dln2*dfloat((imu+inu-2)/2))
           coum(indmu,indnu) = fas*dexp(fac)
         end do ! inu
      end do ! imu
c      call printd(coum,nmax,nmax,nmax,'  coum  ')
      return
      end
       SUBROUTINE GAUSSQ(KIND,N,ALPHA,BETA,KPTS,ENDPTS,B,T,W)
       IMPLICIT REAL*8(A-H,O-Z)
       REAL*8 MUZERO
       DIMENSION B(N),T(N),W(N),ENDPTS(2)
       CALL CLASS(KIND,N,ALPHA,BETA,B,T,MUZERO)
       IF(KPTS.EQ.0)GO TO 100
       IF(KPTS.EQ.2)GO TO 50
       T(N)=GBSLVE(ENDPTS(1),N,T,B)*B(N-1)**2+ENDPTS(1)
       GO TO 100
  50  GAM=GBSLVE(ENDPTS(1),N,T,B)
       T1=((ENDPTS(1)-ENDPTS(2))/(GBSLVE(ENDPTS(2),N,T,B)-GAM))
       B(N-1)=DSQRT(T1)
       T(N)=ENDPTS(1)+GAM*T1
 100  W(1)=1.D0
       DO 105 I=2,N
 105  W(I)=0.D0
       CALL GBTQL2(N,T,B,W,IERR)
       DO 110 I=1,N
 110  W(I)=MUZERO*W(I)*W(I)
       RETURN
       END
      DOUBLE PRECISION FUNCTION GBSLVE(SHIFT,N,A,B)
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION A(N),B(N)
       ALPHA=A(1)-SHIFT
       NM1=N-1
       DO 10 I=2,NM1
  10  ALPHA=A(I)-SHIFT-B(I-1)**2/ALPHA
       GBSLVE=1.D0/ALPHA
       RETURN
       END
       SUBROUTINE CLASS(KIND,N,ALPHA,BETA,B,A,MUZERO)
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSIONA(N),B(N)
       REAL*8MUZERO
       DATA PI / 3.141592653589793D0 /
       NM1=N-1
       GO TO (10,20,30,40,50,60),KIND
  10  MUZERO=2.D0
       DO 11 I=1,NM1
       A(I)=0.D0
       ABI=I
  11  B(I)=ABI/DSQRT(4.D0*ABI*ABI-1.D0)
       A(N)=0.D0
       RETURN
  20  MUZERO=PI
       DO 21 I=1,NM1
       A(I)=0.D0
  21  B(I)=0.5D0
       B(1)=DSQRT(0.5D0)
       A(N)=0.D0
       RETURN
  30  MUZERO=PI/2.D0
       DO 31 I=1,NM1
       A(I)=0.D0
  31  B(I)=0.5D0
       A(N)=0.D0
       RETURN
  40  MUZERO=DSQRT(PI)
       DO 41 I=1,NM1
       A(I)=0.D0
       DI20=I/2.D0
  41  B(I)=DSQRT(DI20)
       A(N)=0.D0
       RETURN
  50  AB=ALPHA+BETA
       ABI=2.D0+AB
       MUZERO=2.D0**(AB+1.D0)*DGAMMA(ALPHA+1.D0)*DGAMMA(BETA+1.D0)/DGAM
     1MA(ABI)
       A(1)=(BETA-ALPHA)/ABI
       B(1)=DSQRT(4.D0*(1.D0+ALPHA)*(1.D0+BETA)/((ABI+1.D0)*ABI*ABI))
       A2B2=BETA*BETA-ALPHA*ALPHA
       DO 51 I=2,NM1
       ABI=2.D0*I+AB
       A(I)=A2B2/((ABI-2.D0)*ABI)
       FI=I
  51  B(I)=DSQRT(4.D0*FI*(FI+ALPHA)*(FI+BETA)*(FI+AB)/((ABI*ABI-1.D0)*A
     1BI*ABI))
       ABI=2.D0*N+AB
       A(N)=A2B2/((ABI-2.D0)*ABI)
       RETURN
  60  MUZERO=DGAMMA(ALPHA+1.D0)
       DO 61 I=1,NM1
       FI = I
       A(I)=2.D0*FI-1.D0+ALPHA
  61  B(I)=DSQRT(FI*(FI+ALPHA))
       A(N)=2.D0*N-1.D0+ALPHA
       RETURN
      END
      SUBROUTINE GBTQL2 ( N,D,E,Z, IERR )
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER I,J,K,L,M,N,II,MML,IERR
      DIMENSION D(N),E(N),Z(N)
      REAL*8 MACHEP
        MACHEP = 16.0D0**(-14)
      IERR=0
      IF (N .EQ. 1 ) GO TO 1001
      E(N) = 0.0
      DO 240 L= 1,N
        J=0
 105   DO 110 M=L,N
            IF (M .EQ. N) GO TO 120
            IF ( DABS(E(M)) .LE. MACHEP * (DABS(D(M))  + DABS(D(M+1))))
     *         GO TO 120
 110   CONTINUE
 120   P = D(L)
       IF ( M .EQ. L ) GO TO 240
       IF (J .EQ. 30) GO TO 1000
       J=J+1
       G = (D(L+1) - P) /(2.0 * E(L))
       R =DSQRT( G*G+1.0 )
       G = D(M) - P+ E(L) / (G +DSIGN(R , G))
       S = 1.0
       C = 1.0
       P = 0.0
       MML = M - L
       DO 200 II = 1, MML
             I = M - II
             F = S * E(I)
             B = C * E(I)
             IF (DABS(F) .LT.DABS(G)) GO TO 150
             C = G / F
             R =DSQRT(C*C+1.0 )
             E(I+1) = F*R
             S = 1.0 /R
             C = C * S
             GO TO 160
 150         S = F / G
             R =DSQRT(S*S+1.0 )
             E(I+1) = G * R
             C = 1.0  / R
             S = S * C
 160         G = D(I+1) - P
             R = ( D(I) - G) * S + 2.0 * C * B
             P = S * R
             D(I+1) = G + P
             G = C * R - B
             F = Z(I+1)
             Z(I+1) = S * Z(I) + C * F
             Z(I) = C * Z(I) - S * F
 200  CONTINUE
       D(L) = D(L) - P
       E(L) = G
       E(M) = 0.0
       GO TO 105
 240  CONTINUE
      DO 300 II = 2 , N
       I = II - 1
       K = I
       P = D(I)
       DO 260 J = II,N
             IF ( D(J) .GE. P) GO TO 260
             K = J
             P = D(J)
 260   CONTINUE
       IF (K .EQ. I ) GO TO 300
       D(K) = D(I)
       D(I) = P
       P = Z(I)
       Z(I) = Z(K)
       Z(K) = P
 300  CONTINUE
      GO TO 1001
1000  IERR = L
1001  RETURN
      END
      DOUBLE PRECISION    FUNCTION DGAMMA(Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION  A(18)
        A(1)=1.0
        A(2)=.4227843350984678
        A(3)=.4118403304263672
        A(4)=.0815769192502609
         A(5)=.0742490106800904
        A(6)=-.0002669810333484
        A(7)=.0111540360240344
        A(8)=-.0028525821446197
        A(9)=.0021036287024598
        A(10)=-.0009184843690991
        A(11)=.0004874227944768
        A(12)=-.0002347204018919
        A(13)=.0001115339519666
        A(14)=-.0000478747983834
        A(15)=.0000175102727179
        A(16)=-.0000049203750904
        A(17)=.0000009199156407
        A(18)=-.0000000839940496
        IF(Z .LE. 1.0 ) GO TO 10
        IF ( Z .LE. 2.0 ) GO TO 20
        T = Z - 2.0
        GO TO 30
  10    T = Z
        GO TO 30
 20     T = Z-1.0
 30     P = A(18)
        DO 40 K1 = 1,17
        K = 18 -K1
        P = T*P+ A(K)
  40        CONTINUE
        IF(Z .GT. 2.0 ) GO TO 50
        IF (Z .GT. 1.0 ) GO TO 60
        DGAMMA = P /(Z * (Z+1.0))
        RETURN
 60     DGAMMA = P / Z
        RETURN
 50     DGAMMA = P
        RETURN
        END
      subroutine howf1(hz,axh,nzmax,nnh,iop)
      implicit real*8 (a-h,o-u),logical(v)
      dimension hz(nzmax,nnh),axh(nnh)
c-------------------------------------------------------------------
c>> iop = 0 computes phi(n)
c>> iop = 1 computes phi(n)*exp(-0.5*axh*axh))
c>>                       1/2  n    -1/2
c>> where phi (x) = ( (pi)    2  n!)      h  (x)  is the HOWF for b=1
c>>          n                              n
c>>
c------------------------------------------------------------------
      sq2=1.4142135623730950488d+00
      pi4=7.511255444649424828587d-1
c
      viop= iop.eq.0
c
c----------------------------------------------------
c    oscilador armonico en una dimension
c----------------------------------------------------
      do 5 j=1,nnh
      az=axh(j)
      hz(1,j)=pi4
      hz(2,j)=sq2*az*pi4
      ri=1.0d0
c
      do 6 i=3,nzmax
      ri=ri+1.0d0
      ri1=dsqrt(ri-1.0d0)
      aa=1.d0/dsqrt(ri)
    6 hz(i,j)=sq2*aa*az*hz(i-1,j)-ri1*aa*hz(i-2,j)
c
      if(viop) goto 5
      chf = dexp(-0.5*az*az)
      do 11 i=1,nzmax
   11 hz(i,j)=hz(i,j)*chf
    5 continue
      return
      end
      
c+---------------------------------------------------------------------+
c|  Computes ro(r), etc                                                |
c|                                                                     |
c|   Ro           only k1>=k2 is stored in GROP  Pos parity            |
c|     k2 k1                               GROM  Neg parity            |
c|                                                                     |
c|   GROP: rosp, roap, rospw, roapw, akasp, akaap, akaspw, akaapw      |
c|                                                                     |
c|   GROM: rosm, roam, rosmw, roamw, akasm, akaam, akasmw, akaamw      |
c|                                                                     |
c|   ROR : ror0p ror0n ror1p ror1n ror2p ror2n ror3p ror3n             |
c|                                   3                              3  |
c| each of them with a nherm3 = nherm  dimension   nherm38 = 8 nherm   |
c|                                                                     |
c|                                                                     |
c+---------------------------------------------------------------------+
      Subroutine denr(GROP,GROM,ror,wfho,wfhonec,rorz,roryz)
c
      Implicit real*8 (A-H,O-Z)
      Logical lcx,lcy,lx12,lxy12
      Logical LMZE,LMZO
      Include 'COMDIM'
      Dimension sum(8)
c
      Dimension GROP(nrop8)
      Dimension GROM(nrom8)
c
      Dimension ror(nherm38),wfho(nmax,nherm),wfhonec(nmax,nherm)
c
      Dimension rorz (irorz) ! Temporary vectors
      Dimension roryz(irorY) !
c
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      Common/DIMEN/kmax,kmax1,kxmax,kymax,kzmax,nwf2,
     * maxtz1,maxjz1,maxjx1,maxjy1,maxlso,maxdso,
     * ndmu,ndthet,nacou,nrop,nrom,nrop8,nrom8,nherm,nherm38,nlegn
c
      Common/DIMPDD/ irortmp,idmus0,jdmus0,idmus1,jdmus1,irorz,irory
c
      Common/OSCLEN/bx,by,bz
      Common/NECKCO/ aneck,rneckp,rneckn
      Common/PULGA/ideb
c
      bm1= 1.0d+00/(bx*by*bz)
      fneck = 2.0d+00/Dsqrt(1.d+00 + (bz/aneck)**2)
      nherm3 = nherm**3
      nherm4 = nherm *4
      nherm2 = nherm**2
      nherm26= nherm**2 * 6
c
      do jher=1,nherm38
         ror(jher) = 0.0d+00
      end do
c
      nh1 = nherm3
      nh2 = nh1 + nherm3
      nh3 = nh2 + nherm3
      nh4 = nh3 + nherm3
      nh5 = nh4 + nherm3
      nh6 = nh5 + nherm3
      nh7 = nh6 + nherm3
c
c---->        K1 >= K2
c---->    e sufix stands for even parity
c---->    o sufix stands for odd  parity
c
      irop = 0
      irom = 0
      xneckp = 0.0d+00
      xneckn = 0.0d+00
      do 1 nx1=1,nxmax
        nym1 = MY(nx1)
        do 11 nx2=1,nx1
          phasx= Dfloat(1-2*Mod((nx2+1),2))
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          lcx  = Mod((nx1+nx2),2).eq.0
c
              do jher=1,nherm26
                 roryz(jher) = 0.0d+00
              end do
c
          do 2 ny1 = 1,nym1
            nxy1 = nx1+ny1
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do 12 ny2=1,nym2
              lcy  = Mod((ny1+ny2),2) .eq.0
              phas = Dfloat(1-2*Mod(Iabs(ny1-ny2)/2,2))
              phase= Dfloat(1-2*Mod(Iabs(ny1-ny2+1)/2,2))
              nxy2 = nx2+ny2
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              lxy12 = lx12.and.(ny1.eq.ny2)
c
              do jher=1,nherm4
                 rorz(jher) = 0.0d+00
              end do
c                                                  ----> Positive parity
c
c be careful not to move the if below, the iz1e=iy1e need to be done
c
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
c
              do 30 nz1 =nzi1e,nzm1,2
                if(lxy12) nzm2=nz1
                do 31 nz2 = nzi2e,nzm2,2
                feq = 2.0d+00
                if(lxy12.and.(nz1.eq.nz2)) feq = 1.0d+00
c
c    GROP: rosp, roap, rospw, roapw, akasp, akaap, akaspw, akaapw
c
c        sum(1) =  ro1n+ro2n
c        sum(2) =  ro1n-ro2n
c        sum(3) =  ro1p+ro2p
c        sum(4) =  ro1p-ro2p
c        sum(5) =  ak1n+ak2n
c        sum(6) =  ak1n-ak2n
c        sum(7) =  ak1p+ak2p
c        sum(8) =  ak1p-ak2p
c
         call dcopy(8,grop(irop+1),1,sum,1)
         irop = irop + 8
                if((.not.lcy).and.lcx) goto 31     ! no field 2 density
                ro0n = sum(1)
                ro0p = sum(3)
                ro2n = sum(2)
                ro2p = sum(4)
                iher = 1
                do nzher=1,nherm
                   wfz = wfho(nz1,nzher)*wfho(nz2,nzher)
                   rorz(iher)  =rorz(iher  )+ro0p*feq*wfz
                   rorz(iher+1)=rorz(iher+1)+ro0n*feq*wfz
                   rorz(iher+2)=rorz(iher+2)+ro2p*feq*wfz
                   rorz(iher+3)=rorz(iher+3)+ro2n*feq*wfz
                   iher = iher + 4
                end do
c
                if(lxy12) then
                   do nzher=1,nherm
                     wfz = wfhonec(nz1,nzher)*wfhonec(nz2,nzher)
                     xneckp = xneckp + ro0p*feq*wfz
                     xneckn = xneckn + ro0n*feq*wfz
                   end do
                end if
c
   31           continue
   30         continue
c                                                   ----> Negative parity
              end if
c
c be careful not to move the if below, the iz1o=iy1o need to be done
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
c
              do 40 nz1 =nzi1o,nzm1,2
                if(lxy12) nzm2=nz1
                do 41 nz2 = nzi2o,nzm2,2
                feq = 2.0d+00
                if(lxy12.and.(nz1.eq.nz2)) feq = 1.0d+00
c
c    GROM: rosm, roam, rosmw, roamw, akasm, akaam, akasmw, akaamw
c
         call dcopy(8,grom(irom+1),1,sum,1)
         irom = irom + 8
                if((.not.lcy).and.lcx) goto 41   ! no field 2 density
                ro0n = sum(1)
                ro0p = sum(3)
                ro2n = sum(2)
                ro2p = sum(4)
                iher = 1
                do nzher=1,nherm
                   wfz = wfho(nz1,nzher)*wfho(nz2,nzher)
                   rorz(iher)  =rorz(iher  )+ro0p*feq*wfz
                   rorz(iher+1)=rorz(iher+1)+ro0n*feq*wfz
                   rorz(iher+2)=rorz(iher+2)+ro2p*feq*wfz
                   rorz(iher+3)=rorz(iher+3)+ro2n*feq*wfz
                   iher = iher + 4
                end do
c
                if(lxy12) then
                   do nzher=1,nherm
                     wfz = wfhonec(nz1,nzher)*wfhonec(nz2,nzher)
                     xneckp = xneckp + ro0p*feq*wfz
                     xneckn = xneckn + ro0n*feq*wfz
                   end do
                end if
c
c
   41           continue
   40         continue
            end if
c
            iher = 1
            do nyher=1,nherm
               wfy = wfho(ny1,nyher)*wfho(ny2,nyher)
               iherz = 1
               do nzher=1,nherm
                  if(lcy) then
                  roryz(iher  )=roryz(iher  )+rorz(iherz  )*wfy*phas
                  roryz(iher+1)=roryz(iher+1)+rorz(iherz+1)*wfy*phas
                  roryz(iher+2)=roryz(iher+2)+rorz(iherz+2)*wfy*phas
                  roryz(iher+3)=roryz(iher+3)+rorz(iherz+3)*wfy*phas
                  else
                  roryz(iher+4)=roryz(iher+4)+rorz(iherz+2)*wfy*phase
                  roryz(iher+5)=roryz(iher+5)+rorz(iherz+3)*wfy*phase
                  end if
               iherz = iherz + 4
               iher  = iher  + 6
               end do
            end do
c
   12       continue
    2     continue
          iher = 0
          do nxher=1,nherm
             wfx = wfho(nx1,nxher)*wfho(nx2,nxher)
             ihery = 1
             do nyzher=1,nherm2
                iher = iher + 1
                if(lcx) then
                ror(iher    )=ror(iher    )+roryz(ihery  )*wfx
                ror(iher+nh1)=ror(iher+nh1)+roryz(ihery+1)*wfx
                ror(iher+nh4)=ror(iher+nh4)+roryz(ihery+2)*wfx*phasx
                ror(iher+nh5)=ror(iher+nh5)+roryz(ihery+3)*wfx*phasx
                else
                ror(iher+nh2)=ror(iher+nh2)+roryz(ihery+2)*wfx
                ror(iher+nh3)=ror(iher+nh3)+roryz(ihery+3)*wfx
                ror(iher+nh6)=ror(iher+nh6)+roryz(ihery+4)*wfx*phasx
                ror(iher+nh7)=ror(iher+nh7)+roryz(ihery+5)*wfx*phasx
                end if
                ihery = ihery + 6
             end do
          end do
   11   continue
    1 continue
c
      do nher=1,nherm38
         ror(nher) = ror(nher)*bm1
      end do
c
      rneckp = xneckp*fneck
      rneckn = xneckn*fneck
      return
      end
c+---------------------------------------------------------------------+
c    Version 4.0
c
c    This program is used to compute the dimensions needed in the
c    triaxial code and also to create the include files.
c+---------------------------------------------------------------------+
c    The basis is determined by using the energy condition
c
c       a  n    + a  n   + a  n   <  N
c        x  x      y  y     z  z  =   0
c
c
c          1/3   -1/3            2/3   1/3           -1/3  -2/3
c    a  = q    p           a  = p    q         a  = p     q
c     x                     y                   z
c
c                                         Rz           Ry
c     q and p are the axis ratios     q = ---  and p = ---
c                                         Rx           Rx
c
c
c+---------------------------------------------------------------------+
c     In the particular case of p=1 the condition reads
c
c          1                    N0
c         --- n  + (n +n ) <  ------
c          q   z     x  y  =  q^1/3
c
c      which is the same as in the axial case but with 
c
c          N  (triaxial) = q^1/3 N  (AXIAL)
c           0                     0
c+---------------------------------------------------------------------+
c
      Implicit real*8 (A-H,O-Z)
      Parameter (maxdim=30)
      Logical LMZE(maxdim,maxdim),LMZO(maxdim,maxdim)
      Character*1 sep(maxdim)
      Dimension MY(maxdim),MZ(maxdim,maxdim)
      Dimension NZIE(maxdim,maxdim),NZIO(maxdim,maxdim)
c
      do i=1,maxdim
         sep(I)=','
      end do
c
      write(6,*) ' Enter q p and N0  '
      read (5,*) q,p,rn0
c
      open(unit=7,file='DIMPERM',status='NEW')
      open(unit=8,file='COMDIM',status='NEW')
      open(unit=9,file='DIMTRIAX',status='NEW')
c

      onethird = 1.d+00/3.d+00
      twothird = 2.d+00/3.d+00
      ax= q**onethird /(p**onethird)
      ay= p**twothird *(q**onethird)
      az= 1./(p**onethird * q**twothird)
c
      nxmax = rn0/ax
      nxmax = nxmax + 1
      nymax = rn0/ay
      nymax = nymax + 1
      nzmax = rn0/az
      nzmax = nzmax + 1
c
      sep(nxmax)='/'

      Do nx=1,nxmax
         temp   = rn0-(nx-1)*ax
         itemp  = temp/ay
         MY(nx) = itemp + 1
      end do
c
      Do nx=1,nxmax
         nym= MY(nx)
         Do ny=1,nym
            temp  = rn0 -(nx-1)*ax - (ny-1)*ay
            itemp = temp/az
            MZ(nx,ny) = itemp + 1
         end do
      end do
c
      nmax  = Max(nxmax,nymax,nzmax)
      nmax1 = nmax+1
c
      nherm = 2*nmax1
      nlegn = 20
c
c    Computing  LMZE LMZO NZIE NZIO
c
      do nx=1,nxmax
        nym = my(nx)
          do ny=1,nym
             nxy  = nx +ny
             nzie(nx,ny) = 1 + Mod(nxy,2)
             nzio(nx,ny) = 2 - Mod(nxy,2)
             nzm = mz(nx,ny)
             lmze(nx,ny) = .false.
             lmzo(nx,ny) = .false.
             if(nzm.ge.nzie(nx,ny)) lmze(nx,ny) = .true.
             if(nzm.ge.nzio(nx,ny)) lmzo(nx,ny) = .true.
          end do
      end do
c
      iz1e = 0
      iz1o = 0
      jz1e = 0
      jz1o = 0
      do nx1=1,nxmax
        nym1 = my(nx1)
          do ny1 = 1,nym1
            nxy1 = nx1+ny1
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = mz(nx1,ny1)
            if(lmze(nx1,ny1)) then
              iz1e=iz1e+(nzm1-nzi1e)/2+1  ! pos parity
            end if
            if(lmzo(nx1,ny1)) then
              iz1o=iz1o+(nzm1-nzi1o)/2+1  ! neg parity
            end if
          end do ! ny1
      end do     ! nx1
c
      NP = iz1e     ! Positive parity matrices dimension NPxNP
      NM = iz1o     ! Negative parity matrices dimension NMxNM
      NPM= Max(NP,NM)
c
      nrop = (np*(np+1))/2   !   Ro dimensions    Pos parity
      nrom = (nm*(nm+1))/2   !                    Neg parity
      nrop8 = nrop*8
      nrom8 = nrom*8
      ngp   = nrop
      ngm   = nrom
      ngp8  = nrop8
      ngm8  = nrom8
c+--------------------------------------------------------------------+
c|     Permanent storage                                              |
c|                                     Direct                         |
c+--------------------------------------------------------------------+
      Maxtz1 = (Nmax*(Nmax**2 + 3*Nmax + 2))/6
      Maxtz2 = ((Nmax+1)*((Nmax+1)**2 + 3*Nmax + 5))/6
      Maxjx1 = ((Nxmax+1)**2/4)*Nxmax + (Nxmax**2/4)*(Nxmax-1)
      Maxjy1 = ((Nymax+1)**2/4)*Nymax + (Nymax**2/4)*(Nymax-1)
      Maxjz1 = ((Nzmax+1)**2/4)*Nzmax + (Nzmax**2/4)*(Nzmax-1)
      Maxlso = ((Nmax+2)**2/4)*(Nmax+1) + ((Nmax+1)**2/4)*Nmax
      Maxdso = ((Nmax+1)**2/2-Nmax)*(Nmax+1)+(Nmax**2/2)*Nmax
      Nacou  = (2*Nxmax-1)*(2*Nymax-1)*(2*Nzmax-1)
      Ncoulm = Nmax**2
      nwf2   = 2*Nmax-1
      Nwavef = Nmax*Nherm
      Nwavefd= nwf2*Nherm
c+--------------------------------------------------------------------+
c|                                     Exchange (SXTJ)                |
c+--------------------------------------------------------------------+
      nxmax2 = nxmax/2
      nymax2 = nymax/2
      nzmax2 = nzmax/2
      maxtjx = ((nxmax2+1)*(nxmax-nxmax2))**2+
     *         ( (nxmax*(nxmax+1))/2-(nxmax2+1)*(nxmax-nxmax2) )**2
      maxtjy = ((nymax2+1)*(nymax-nymax2))**2+
     *         ( (nymax*(nymax+1))/2-(nymax2+1)*(nymax-nymax2) )**2
      maxtjz = ((nzmax2+1)*(nzmax-nzmax2))**2+
     *         ( (nzmax*(nzmax+1))/2-(nzmax2+1)*(nzmax-nzmax2) )**2
      maxx   = 4*maxtjx+2*maxtjy+2*maxtjz
c
      NPERTOT = Maxtz1 + Maxtz2 + 2*Maxjx1 + 2*Maxjy1 + 2*Maxjz1 +
     *          Maxlso + Maxdso + Nacou    + Ncoulm   + Nwavef   +
     *          Nwavefd+ 2*Nherm+  2*Nlegn + maxx
c+--------------------------------------------------------------------+
c|   Temporary                                                        |
c|                                                  Direct            |
c+--------------------------------------------------------------------+
      NDMU    = nxmax*nymax*nzmax
      NDTHETA = Ndmu*(nzmax/2+1)
      nherm38 = 8*nherm**3
      irorz  = Nherm*4
      irory  = Nherm*Nherm*6
      irortmp= 2*Nherm**3
      idmus0 = nherm**2 * nmax
      idmus1 = nherm    * nmax**2
      jdmus0 = 15*idmus0
      jdmus1 = 15*idmus1
c
      NDIRT1  = 15*NDMU + nherm38 + irorz + irory
      NDIRT2  = 15*NDMU + nherm38 + irortmp + jdmus0 + jdmus1
      NDIRT3  = 61 * Ndmu + 20 * NDTHETA
      NDIRT   = MAX(NDIRT1,NDIRT2,NDIRT3)
      write(6,*) ' TEMP DIM DIR ',NDIRT1,NDIRT2,NDIRT3,NDIRT
c+--------------------------------------------------------------------+
c|                                             Exchange and pairing   |
c+--------------------------------------------------------------------+
      mazopti= 16*nzmax*nzmax
      mayopti= 16*nzmax**2 * nymax**2
c
      kk = 0
      do ny1=1,nymax
        nz1m = MZ(1,ny1)
        do ny2=1,nymax
           nz2m = MZ(1,ny2)
           kk = kk + nz1m*nz2m
        end do
      end do
      kk = kk * 16
      write(6,*) ' mayopti ',mayopti,kk
c
      NECHT  = mazopti + mayopti
      NTEMP  = MAX(NDIRT,NECHT)
c

      NUVP = np*np
      NUVM = nm*nm
      nropt = 16*nrop
      nromt = 16*nrom
      nuvtp = 16*nuvp
      nuvtm = 16*nuvm
      nuvtt = 8*(nuvp+nuvm)
c
      rMb    = 1.0d+00/dfloat(1024*1024)
      rMbropt = dfloat(nropt*8)*rMb
      rMbromt = dfloat(nromt*8)*rMb
      rMBnuvpt= dfloat(nuvtp*8)*rMb
      rMBnuvmt= dfloat(nuvtm*8)*rMb
      rmaz   = 8.0d+00*Dfloat(mazopti)*rMB
      rmay   = 8.0d+00*Dfloat(mayopti)*rMb
c
      NTEMP  = MAX(NTEMP,NUVTT)
c
      write(6,1001) np,nm
      write(6,1002) nrop,nrom
      write(6,1005) rmbropt,rmbromt
      write(6,1006) mazopti,rmaz,mayopti,rmay
      write(6,1007) NDIRT,NECHT,NUVTT,NTEMP
      write(6,1008) nuvtp,rmbnuvpt,nuvtm,rmbnuvmt
      write(6,1010) NUVTT,nrop8,nrom8,ntemp
c
c     Dimperm
c
      write(7,100) nmax,nmax1,nxmax,nymax,nzmax,
     * nwf2,maxtz1,maxtz2,maxjx1,maxjy1,maxjz1,maxlso,
     * maxdso,ndmu,ndtheta,nacou,nherm,nlegn,maxtjx,maxtjy,maxtjz
  100 format ( '       Parameter (imax=',i2,',imax1=',i2,
     * ',ixmax=',i2,',iymax=',i2,',izmax=',i2,/,
     *5x,'* ,iwf2 =',i3,',iaxtz1d=',i6,',iaxtz2d=',i6,',iaxjx1d=',i6,/,
     *5x,'* ,iaxjy1d=',i6,',iaxjz1d=',i6,',iaxlso=',i6,',iaxdso=',i6,/,
     *5x,'* ,idmu=',i8,',idthet=',i8,',iacou=',i6,',iherm=',i2,/,
     *5x,'* ,ilegn= ',i2,',iaxtjx=',i8,',iaxtjy=',i8,',iaxtjz=',i8,')')
c
      write(7,101)
  101 format( '       Dimension iMY(ixmax),iMZ(ixmax,iymax)')
      write(7,2000) (my(k),sep(k),k=1,nxmax)
      write(7,2001) (mz(k,1),k=1,nxmax)
      do j=2,(nymax-1)
          write(7,2002) (mz(k,j),k=1,nxmax)
      end do
          write(7,2003) (mz(k,nymax),sep(k),k=1,nxmax)
c
c   Comdim
c
      izsrc = nzmax
      iyzsrc= nzmax*nymax
      Nfac  = 100
      iqne  = NP
      iqno  = NM
      maxp  = 11
      Kdim  = 30
c
      write(8,200) nmax,nmax1,nxmax,nymax,nzmax,izsrc,iyzsrc,nfac,
     *   iqne,iqno,maxp,kdim,q,p,rn0
  200 format ( '       Parameter (imax=',i2,',imax1=',i2,
     * ',ixmax=',i2,',iymax=',i2,',izmax=',i2,/,
     * 5x'*  ,izsrc=',i2,',iyzsrc=',i5,',nfac=',i3,',inqne=',i4,/,
     * 5x'*  ,inqno=',i4,',maxp=',i2,',kdim=',i2,')',/,
     *  '       Parameter (q0inf=',f7.4,',p0inf=',f7.4,',rn0inf=',
     *  f8.4,')')
c
c    DIMTRIAX
c
      write(9,300) NP,NM,NPM,NUVTT,NTEMP,NROP8,NROM8
  300 format('      Parameter (NP=',i3,',',/,
     *       '     *           NM=',i3,',',/,
     *       '     *          NPM=',i3,',',/,
     *       '     *          NUV=',i8,',',/,
     *       '     *         NH20=',i8,',',/,
     *       '     *          IGP=',i8,',',/,
     *       '     *          IGM=',i8,')')

c  -----------------+ F O R M A T S +------------------------
 1001 format( ' Positive parity states ',i5,
     *        ' Negative parity states ',i5)
c
 1002 format( ' Ro matrix dimensions m>= n ',
     *        ' Positive parity ',i8,
     *        ' Negative parity ',i8)
c
 1005 format( ' For the 16 matrices (densities and fields) ',/,
     * ' Positive parity:',f15.7,'Mb',
     * ' Negative parity:',f15.7,'Mb')
c
 1006 format ( ' Exchange subroutine ',/,
     *         ' Z opt vector        ',i10,'= ',f15.7,' Mb',/,
     *         ' Y opt vector        ',i10,'= ',f15.7,' Mb')
 1007 format ( /,' Size of the TEMP vectors ',/,
     *      ' DIRECT ',I10,' EXCH ',I10,' H20  ',I10,' MAXIMUM ',I10)
c
 1008 format( ' For the 16 matrices (U,V,H20,JX20,N20,JX,JZ) ',/,
     * ' Positive parity, a R*8 vector of ',i10,' and ',f15.7,'Mb',/,
     * ' Negative parity, a R*8 vector of ',i10,' and ',f15.7,'Mb',/)
c
 1010 format( '   D I M E N S I O N S  ',/,
     * ' UV     ',i10,/,
     * ' IGP    ',i10,/,
     * ' IGM    ',i10,/,
     * ' NH20V  ',i10,/)
c
 2000 format('      DATA iMY/',20(i2,A1))
 2001 format('      DATA iMZ/',20(i2,','))
 2002 format('     *        ',20(i2,','))
 2003 format('     *        ',20(i2,A1))
      stop
      end
      Subroutine Gognyf(ipar)
      Implicit double precision (A-H,O-Z)
      CHARACTER*6 FORCE
      Common/GOGINT/amu(2),xw(2),xh(2),xb(2),xm(2),WLS,t3,alpha,x0,e2
      Common/GOGHAN/w2b(2),h2m(2),x01,x02
      Common/GOGECH/c1(2),c2(2),p1(2),p2(2)
c
c+--------------------------------------------------------------------+
c|     electron charge square                                         |
c+--------------------------------------------------------------------+
c
      e2     = 1.44197028 d+00    ! e**2 in Mev/fm
c     e2     = 0.0d+00
c+--------------------------------------------------------------------+
c Gogny's values
c+--------------------------------------------------------------------+
c Spherical
c     e2     = 1.43997865
c+--------------------------------------------------------------------+
c Triaxial
c     e2     = 1.44d+00
c+--------------------------------------------------------------------+
c
      if (ipar.eq.0) then
c+--------------------------------------------------------------------+
c|                                                  D1S  Parameters   |
c+--------------------------------------------------------------------+
      FORCE  = ' D1S  '

      amu(1) = 0.7d+00              ! fm**2
      amu(2) = 1.2d+00              ! fm**2

      xw (1) = -1720.3   d+00
      xw (2) =   103.639 d+00
      xb (1) =  1300.0   d+00
      xb (2) =  -163.483 d+00
      xh (1) = -1813.53  d+00
      xh (2) =   162.812 d+00
      xm (1) =  1397.6   d+00
      xm (2) =  -223.934 d+00

      WLS    =   130.0   d+00

      t3     =  1390.6   d+00

      alpha  = 1.0d+00/3.0d+00
      x0     = 1.0d+00
c
      else if (ipar.eq.1) then
c+--------------------------------------------------------------------+
c|                                                  D1   Parameters   |
c+--------------------------------------------------------------------+
      FORCE  = ' D1   '

      amu(1) = 0.7d+00
      amu(2) = 1.2d+00

      xw (1) =  -402.4   d+00
      xw (2) =   -21.30  d+00
      xb (1) =  -100.0   d+00
      xb (2) =   -11.77  d+00
      xh (1) =  -496.20  d+00
      xh (2) =    37.27  d+00
      xm (1) =   -23.56  d+00
      xm (2) =   -68.81  d+00

      WLS    =   115.0   d+00

      t3     =  1350.0   d+00

      alpha  = 1.0d+00/3.0d+00
      x0     = 1.0d+00
c
      else if (ipar.eq.2) then
c+--------------------------------------------------------------------+
c|                                                  D1'  Parameters   |
c+--------------------------------------------------------------------+
      FORCE  = ' D1 LS'

      amu(1) = 0.7d+00
      amu(2) = 1.2d+00

      xw (1) =  -402.4   d+00
      xw (2) =   -21.30  d+00
      xb (1) =  -100.0   d+00
      xb (2) =   -11.77  d+00
      xh (1) =  -496.20  d+00
      xh (2) =    37.27  d+00
      xm (1) =   -23.56  d+00
      xm (2) =   -68.81  d+00

      WLS    =   130.0   d+00

      t3     =  1350.0   d+00

      alpha  = 1.0d+00/3.0d+00
      x0     = 1.0d+00
c
      else if (ipar.eq.3) then
c+--------------------------------------------------------------------+
c|                                                  D1N  Parameters   |
c| Michel Girod PLB  2008                                             |
c+--------------------------------------------------------------------+
      FORCE  = ' D1N 8'

      amu(1) = 0.8d+00
      amu(2) = 1.2d+00

      xw (1) = -2047.610 d+00
      xw (2) =   293.020 d+00
      xb (1) =  1700.0   d+00
      xb (2) =  -300.780 d+00
      xh (1) = -2414.930 d+00
      xh (2) =   414.590 d+00
      xm (1) =  1519.35  d+00
      xm (2) =  -316.840 d+00

      WLS    =   115.0   d+00

      t3     =  1609.46  d+00

      alpha  = 1.0d+00/3.0d+00
      x0     = 1.0d+00
      
      else if (ipar.eq.4) then
c+--------------------------------------------------------------------+
c|                                                  D1N  Parameters   |
c| Michel Girod June 2006                                             |
c+--------------------------------------------------------------------+
      FORCE  = ' D1N 6'

      amu(1) = 0.8d+00
      amu(2) = 1.2d+00

      xw (1) = -2136.394 d+00
      xw (2) =   309.606 d+00
      xb (1) =  1800.0   d+00
      xb (2) =  -316.813 d+00
      xh (1) = -2544.472 d+00
      xh (2) =   439.596 d+00
      xm (1) =  1585.86  d+00
      xm (2) =  -326.601 d+00

      WLS    =   115.0   d+00

      t3     =  1631.0   d+00

      alpha  = 1.0d+00/3.0d+00
      x0     = 1.0d+00
      
      else if (ipar.eq.5) then
c+--------------------------------------------------------------------+
c|                                                  D1M  Parameters   |
c| Hilaire & Gorieli PRL 2009                                         |
c+--------------------------------------------------------------------+
      FORCE  = ' D1M 9'

      amu(1) = 0.5d+00
      amu(2) = 1.0d+00

      xw (1) = -12797.57 d+00
      xw (2) =    490.95 d+00
      xb (1) =  14048.85 d+00
      xb (2) =   -752.27 d+00
      xh (1) = -15144.43 d+00
      xh (2) =    675.12 d+00
      xm (1) =  11963.89 d+00
      xm (2) =   -693.57 d+00
      
      WLS    =    115.36 d+00
      
      t3     =   1562.22 d+00
      
      alpha  = 1.0d+00/3.0d+00
      x0     = 1.0d+00
      
      end if
c+--------------------------------------------------------------------+
c|     Set them to zero                                               |      
c+--------------------------------------------------------------------+
c      xw (1) =    0.0   d+00
c      xw (2) =   0.0   d+00
c      xb (1) =   0.0   d+00
c      xb (2) =   0.0   d+00
c      xh (1) =   0.0   d+00
c      xh (2) =   0.0   d+00
c      xm (1) =   0.0   d+00
c      xm (2) =   0.0   d+00
c      wls    = 0.0d+00
c      t3     = 0.0d+00
c+--------------------------------------------------------------------+
c|                                   Usefull combinations             |
c+--------------------------------------------------------------------+
      do i=1,2
         w2b(i) = xw(i)*2.0d+00+xb(i)
         h2m(i) = xh(i)*2.0d+00+xm(i)
         c1 (i) = xm(i) + 0.5d+00*xh(i)
         c2 (i) = -(xb(i) + 0.5d+00*xw(i))
         p1 (i) = xw(i)+xm(i)-(xh(i)+xb(i))
         p2 (i) = xw(i)+xb(i)-(xh(i)+xm(i))
      end do
c
      x01 = 1.0d+00+0.5d+00*x0
      x02 = 0.5d+00+x0
c+--------------------------------------------------------------------+
c|                       W R I T I N G    O U T                       |
c+--------------------------------------------------------------------+
      write(6,1000) FORCE
      do i=1,2
         write(6,1001) amu(i),xw(i),xb(i),xh(i),xm(i)
      end do
      write(6,1002) WLS,t3,alpha,x0,e2
 1000 format(20x,' GOGNY FORCE      ',A6,//,5x,
     *  9x,'mu',8x,'W',11x,'B',11x,'H',11x,'M',/)
 1001 format(5x,5(2x,f10.2))
 1002 format(/,3X,' WLS = ',f7.2,'  t3 ',f8.2,
     * ' ALPHA ',f7.4,' x0 ',F5.2,' e2 ',f10.7,//)
c
      return
      end
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
      Subroutine PACK2SQ(OPEPACK,OPOPACK,OPE,OPO,is)
      Implicit real*8(a-h,o-z)
      Logical LMZE,LMZO
      Logical lx12,lxy12
      Include 'COMDIM'
      Dimension OPE(NP,NP),OPO(NM,NM)
      Dimension OPEPACK(NGP),OPOPACK(NGM)
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      Common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
      Common/ILOC/NP,NM,ILOCE(ixmax,iymax),ILOCO(ixmax,iymax)

c
c---->    q'=1  q=2    q' >= q
c
c---->    e sufix stands for even parity
c---->    o sufix stands for odd  parity
      irop = 0
      irom = 0
      do nx1=1,nxmax
        nym1 = MY(nx1)
        do nx2=1,nx1
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          do ny1 = 1,nym1
            iz1e0 = ILOCE(nx1,ny1)
            iz1o0 = ILOCO(nx1,ny1)
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do ny2=1,nym2
              iz2e0 = ILOCE(nx2,ny2)
              iz2o0 = ILOCO(nx2,ny2)
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              lxy12 = lx12.and.(ny1.eq.ny2)
c
c---->    Positive parity
c
c
c be careful not to move the if below, the iz1e=iy1e need to be done
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
c
              iz1e = iz1e0
              do nz1 =nzi1e,nzm1,2
                iz1e = iz1e + 1
                if(lxy12) nzm2=nz1
                iz2e = iz2e0
                do nz2 = nzi2e,nzm2,2
                   iz2e = iz2e + 1
                   irop = irop + 1
                   OPE(iz1e,iz2e)= OPEPACK(irop)
                   if(is.eq.1) then
                      ope(iz2e,iz1e) = ope(iz1e,iz2e)
                   else if (is.eq.-1) then
                     ope(iz2e,iz1e) = - ope(iz1e,iz2e)
                   end if
                end do
              end do
c
c---->    Negative parity
c
              end if
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
c
              iz1o = iz1o0
              do nz1 =nzi1o,nzm1,2
                iz1o = iz1o+ 1
                if(lxy12) nzm2=nz1
                iz2o = iz2o0
                do nz2 = nzi2o,nzm2,2
                   iz2o = iz2o+1
                   irom = irom + 1
                   OPO(iz1o,iz2o) = OPOPACK(irom)
                   if(is.eq.1) then
                     opo(iz2o,iz1o) = opo(iz1o,iz2o)
                   else if (is.eq.-1) then
                     opo(iz2o,iz1o) = - opo(iz1o,iz2o)
                   end if
                end do
              end do
            end if
            end do    ! ny2
          end do      ! ny1
        end do        ! nx2
      end do          ! nx1
      return
      end
      Subroutine SQ2PACK(OPE,OPO,OPEPACK,OPOPACK)
      Implicit real*8(a-h,o-z)
      Logical LMZE,LMZO
      Logical lx12,lxy12
      Include 'COMDIM'
      Dimension OPE(NP,NP),OPO(NM,NM)
      Dimension OPEPACK(NGP),OPOPACK(NGM)
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      Common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
      Common/ILOC/NP,NM,ILOCE(ixmax,iymax),ILOCO(ixmax,iymax)
c
c---->    q'=1  q=2    q' >= q
c
c---->    e sufix stands for even parity
c---->    o sufix stands for odd  parity
      irop = 0
      irom = 0
      do nx1=1,nxmax
        nym1 = MY(nx1)
        do nx2=1,nx1
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          do ny1 = 1,nym1
            iz1e0 = ILOCE(nx1,ny1)
            iz1o0 = ILOCO(nx1,ny1)
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do ny2=1,nym2
              iz2e0 = ILOCE(nx2,ny2)
              iz2o0 = ILOCO(nx2,ny2)
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              lxy12 = lx12.and.(ny1.eq.ny2)
c
c---->    Positive parity
c
c
c be careful not to move the if below, the iz1e=iy1e need to be done
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
c
              iz1e = iz1e0
              do nz1 =nzi1e,nzm1,2
                iz1e = iz1e + 1
                if(lxy12) nzm2=nz1
                iz2e = iz2e0
                do nz2 = nzi2e,nzm2,2
                   iz2e = iz2e + 1
                   irop = irop + 1
                   OPEPACK(irop) = OPE(iz1e,iz2e)
                end do
              end do
c
c---->    Negative parity
c
              end if
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
c
              iz1o = iz1o0
              do nz1 =nzi1o,nzm1,2
                iz1o = iz1o+ 1
                if(lxy12) nzm2=nz1
                iz2o = iz2o0
                do nz2 = nzi2o,nzm2,2
                   iz2o = iz2o+1
                   irom = irom + 1
                   OPOPACK(irom) = OPO(iz1o,iz2o)
                end do
              end do
            end if
            end do    ! ny2
          end do      ! ny1
        end do        ! nx2
      end do          ! nx1
      return
      end
      subroutine d2sq(d1p,d2p,d1m,d2m,aux1p,aux2p,aux1m,aux2m)
      Implicit real*8(a-h,o-z)
      Logical LMZE,LMZO
      Logical lx12,lxy12
      Include 'COMDIM'
      Dimension aux1p(NP,NP),aux1m(NM,NM)
      Dimension aux2p(NP,NP),aux2m(NM,NM)
      Dimension d1p(ngp),d2p(ngp),d1m(ngm),d2m(ngm)
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
c
      Common/ILOC/NP,NM,ILOCE(ixmax,iymax),ILOCO(ixmax,iymax)
c
c---->    q'=1  q=2    q' >= q
c---->    e sufix stands for even parity
c---->    o sufix stands for odd  parity
      irop = 0
      irom = 0
      do nx1=1,nxmax
        nym1 = MY(nx1)
        do nx2=1,nx1
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          do ny1 = 1,nym1
            iz1e0 = ILOCE(nx1,ny1)
            iz1o0 = ILOCO(nx1,ny1)
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do ny2=1,nym2
              iz2e0 = ILOCE(nx2,ny2)
              iz2o0 = ILOCO(nx2,ny2)
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
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
c
              iz1e = iz1e0
              do nz1 =nzi1e,nzm1,2
                iz1e = iz1e + 1
                if(lxy12) nzm2=nz1
                iz2e = iz2e0
                do nz2 = nzi2e,nzm2,2
                iz2e = iz2e + 1
                irop = irop + 1
c		write(6,*) ' d2sq p ',irop,d1p(irop),d2p(irop)
                AUX1P(iz1e,iz2e) = d1p(irop)
                AUX2P(iz2e,iz1e) =-d1p(irop)
                AUX1P(iz2e,iz1e) = d2p(irop)
                AUX2P(iz1e,iz2e) =-d2p(irop)
                end do
              end do
c
c---->    Negative parity
c
              end if
c
c be careful not to move the if below, the iz1o=iy1o need to be done
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
c
              iz1o = iz1o0
              do nz1 =nzi1o,nzm1,2
                iz1o = iz1o+ 1
                if(lxy12) nzm2=nz1
                iz2o = iz2o0
                do nz2 = nzi2o,nzm2,2
                iz2o = iz2o+1
                irom = irom + 1
c		write(6,*) ' d2sq m ',irom,d1m(irom),d2m(irom)
                AUX1M(iz1o,iz2o) = d1m(irom)
                AUX2M(iz2o,iz1o) =-d1m(irom)
                AUX1M(iz2o,iz1o) = d2m(irom)
                AUX2M(iz1o,iz2o) =-d2m(irom)
                end do
              end do
            end if
            end do     ! ny2
          end do       ! ny1
        end do         ! nx2
      end do           ! nx1
c
      if((irop.gt.ngp).or.(irom.gt.ngm))
     *   call errout('Troubles        ','D2SQ    ',4)
c
      return
      end
      subroutine de2sq(d1p,d2p,d1m,d2m,aux1p,aux1m)
      Implicit real*8(a-h,o-z)
      Logical LMZE,LMZO
      Logical lx12,lxy12
      Include 'COMDIM'
      Dimension aux1p(NP,NP),aux1m(NM,NM)
      Dimension d1p(ngp),d2p(ngp),d1m(ngm),d2m(ngm)
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
c
      Common/ILOC/NP,NM,ILOCE(ixmax,iymax),ILOCO(ixmax,iymax)
c
c---->    q'=1  q=2    q' >= q
c---->    e sufix stands for even parity
c---->    o sufix stands for odd  parity
      irop = 0
      irom = 0
      do nx1=1,nxmax
        nym1 = MY(nx1)
        do nx2=1,nx1
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          do ny1 = 1,nym1
            iz1e0 = ILOCE(nx1,ny1)
            iz1o0 = ILOCO(nx1,ny1)
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do ny2=1,nym2
              iz2e0 = ILOCE(nx2,ny2)
              iz2o0 = ILOCO(nx2,ny2)
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
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
c
              iz1e = iz1e0
              do nz1 =nzi1e,nzm1,2
                iz1e = iz1e + 1
                if(lxy12) nzm2=nz1
                iz2e = iz2e0
                do nz2 = nzi2e,nzm2,2
                iz2e = iz2e + 1
                irop = irop + 1
                AUX1P(iz1e,iz2e) = d1p(irop)
                AUX1P(iz2e,iz1e) = d2p(irop)
                end do
              end do
c
c---->    Negative parity
c
              end if
c
c be careful not to move the if below, the iz1o=iy1o need to be done
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
c
              iz1o = iz1o0
              do nz1 =nzi1o,nzm1,2
                iz1o = iz1o+ 1
                if(lxy12) nzm2=nz1
                iz2o = iz2o0
                do nz2 = nzi2o,nzm2,2
                iz2o = iz2o+1
                irom = irom + 1
c		write(6,*) ' d2sq m ',irom,d1m(irom),d2m(irom)
                AUX1M(iz1o,iz2o) = d1m(irom)
                AUX1M(iz2o,iz1o) = d2m(irom)
                end do
              end do
            end if
            end do     ! ny2
          end do       ! ny1
        end do         ! nx2
      end do           ! nx1
c
      if((irop.gt.ngp).or.(irom.gt.ngm))
     *   call errout('Troubles        ','DE2SQ    ',4)
c
      return
      end
      Subroutine JT2disk(ajxe,ajxo,ajye,ajyo,ajze,ajzo,te,to,GROP,GROM)
      Implicit real*8(a-h,o-z)
      Include 'COMDIM'
      Logical LMZE,LMZO
      Logical lx12,lxy12
      Dimension AJXE(NP,NP),AJXO(NM,NM)
      Dimension ajyE(NP,NP),ajyO(NM,NM)
      Dimension ajzE(NP,NP),ajzO(NM,NM)
      Dimension   TE(NP,NP),  TO(NM,NM)
      Dimension GROP(ngp8),GROM(ngm8)
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      Common/ILOC/NP,NM,ILOCE(ixmax,iymax),ILOCO(ixmax,iymax)
c
      Common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
c
      ngp2 = 2*ngp
      ngp3 = 3*ngp
      ngp4 = 4*ngp
      ngm2 = 2*ngm
      ngm3 = 3*ngm
      ngm4 = 4*ngm
c
c---->    q'=1  q=2    q' >= q                <--------
c---->                                        <--------
c---->    e sufix stands for even parity      <--------
c---->    o sufix stands for odd  parity      <--------
      irop = 0
      irom = 0
      do nx1=1,nxmax
        nym1 = MY(nx1)
        do nx2=1,nx1
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          do ny1 = 1,nym1
            iz1e0 = ILOCE(nx1,ny1)
            iz1o0 = ILOCO(nx1,ny1)
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do ny2=1,nym2
              iz2e0 = ILOCE(nx2,ny2)
              iz2o0 = ILOCO(nx2,ny2)
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              lxy12 = lx12.and.(ny1.eq.ny2)
c
c---->    Positive parity
c
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
c
                  iz1e = iz1e0
                  do nz1 =nzi1e,nzm1,2
                    iz1e = iz1e + 1
                    if(lxy12) nzm2=nz1
                    iz2e = iz2e0
                    do nz2 = nzi2e,nzm2,2
                        iz2e = iz2e + 1
                        irop = irop + 1
                        GROP(irop     ) =   TE(iz1e,iz2e)
                        GROP(irop+ngp ) = AJXE(iz1e,iz2e)
                        GROP(irop+ngp2) = AJZE(iz1e,iz2e)
                        GROP(irop+ngp3) = AJYE(iz1e,iz2e)
                    end do    ! nz1
                  end do      ! nz2
              end if      ! positive parity condition
c
c---->    Negative parity
c
c
c be careful not to move the if below, the iz1o=iy1o need to be done
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
c
                  iz1o = iz1o0
                  do nz1 =nzi1o,nzm1,2
                    iz1o = iz1o+ 1
                    if(lxy12) nzm2=nz1
                    iz2o = iz2o0
                    do nz2 = nzi2o,nzm2,2
                       iz2o = iz2o+1
                       irom = irom + 1
                       GROM(irom     ) =   TO(iz1o,iz2o)
                       GROM(irom+ngm ) = AJXO(iz1o,iz2o)
                       GROM(irom+ngm2) = AJZO(iz1o,iz2o)
                       GROM(irom+ngm3) = AJYO(iz1o,iz2o)
                    end do  ! nz2
                  end do    ! nz1
              end if    ! negative parity condition
            end do      ! ny2
          end do        ! ny1
        end do          ! nx2
      end do            ! nx1
c
      REWIND 98
      REWIND 99
c
      WRITE(98) (grop(j),j=1,ngp4)
      WRITE(99) (grom(j),j=1,ngm4)
c
      REWIND 98
      REWIND 99
c
      return
      end
      Subroutine RO2PACK(ROM,AKA,grop,grom)
      Implicit real*8(a-h,o-z)
      Logical LMZE,LMZO
      Logical lx12,lxy12
      Include 'COMDIM'
      Dimension sum(8)
      Dimension ROM(*),AKA(*)
      Dimension grop(NGP8),grom(NGM8)
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      Common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
      Common/ILOC/NP,NM,ILOCE(ixmax,iymax),ILOCO(ixmax,iymax)
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC1/ NNU(4),NNV(4),NNRO(4),NNKA(4)
c
c---->    q'=1  q=2    q' >= q
c
c---->    e sufix stands for even parity
c---->    o sufix stands for odd  parity
      irop = 0
      irom = 0
      do nx1=1,nxmax
        nym1 = MY(nx1)
        do nx2=1,nx1
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          do ny1 = 1,nym1
            iz1e0 = ILOCE(nx1,ny1)
            iz1o0 = ILOCO(nx1,ny1)
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do ny2=1,nym2
              iz2e0 = ILOCE(nx2,ny2)
              iz2o0 = ILOCO(nx2,ny2)
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              lxy12 = lx12.and.(ny1.eq.ny2)
c
c---->    Positive parity
c
c
c be careful not to move the if below, the iz1e=iy1e need to be done
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
c
              iz1e = iz1e0
              do nz1 =nzi1e,nzm1,2
                iz1e = iz1e + 1
                if(lxy12) nzm2=nz1
                iz2e = iz2e0
                do nz2 = nzi2e,nzm2,2
                   iz2e = iz2e + 1
		   
		   ro1n =  ROM(iz2e-1+(iz1e-1)*ND(3)+NNRO(3))       ! ro1 neut
		   ro2n =  ROM(iz2e-1+(iz1e-1)*ND(3)+NNRO(3)+ND2(3))! ro2 neut
		   ro1p =  ROM(iz2e-1+(iz1e-1)*ND(1)+NNRO(1))       ! ro1 prot
		   ro2p =  ROM(iz2e-1+(iz1e-1)*ND(1)+NNRO(1)+ND2(1))! ro2 prot
		   
		   ak1n =  AKA(iz2e-1+(iz1e-1)*ND(3)+NNKA(3))	  ! k 1 neut
		   ak2n =  AKA(iz1e-1+(iz2e-1)*ND(3)+NNKA(3))	  !-k 2 neut
		   ak1p =  AKA(iz2e-1+(iz1e-1)*ND(1)+NNKA(1))	  ! k 1 prot
		   ak2p =  AKA(iz1e-1+(iz2e-1)*ND(1)+NNKA(1))	  !-k 2 prot
		   
                   sum(1) =  ro1n+ro2n
                   sum(2) =  ro1n-ro2n
                   sum(3) =  ro1p+ro2p
                   sum(4) =  ro1p-ro2p
                   sum(5) =  ak1n+ak2n
                   sum(6) =  ak1n-ak2n
                   sum(7) =  ak1p+ak2p
                   sum(8) =  ak1p-ak2p
		   
                do ir=1,8
                   irop = irop + 1
                   GROP(irop) = sum(ir)
                end do
		
                end do
              end do
c
c---->    Negative parity
c
              end if
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
c
              iz1o = iz1o0
              do nz1 =nzi1o,nzm1,2
                iz1o = iz1o+ 1
                if(lxy12) nzm2=nz1
                iz2o = iz2o0
                do nz2 = nzi2o,nzm2,2
                   iz2o = iz2o+1
		   
		   ro1n =  ROM(iz2o-1+(iz1o-1)*ND(4)+NNRO(4))       ! ro1 neut
		   ro2n =  ROM(iz2o-1+(iz1o-1)*ND(4)+NNRO(4)+ND2(4))! ro2 neut
		   ro1p =  ROM(iz2o-1+(iz1o-1)*ND(2)+NNRO(2))       ! ro1 prot
		   ro2p =  ROM(iz2o-1+(iz1o-1)*ND(2)+NNRO(2)+ND2(2))! ro2 prot
		   
		   ak1n =  AKA(iz2o-1+(iz1o-1)*ND(4)+NNKA(4))	  ! k 1 neut
		   ak2n =  AKA(iz1o-1+(iz2o-1)*ND(4)+NNKA(4))	  !-k 2 neut
		   ak1p =  AKA(iz2o-1+(iz1o-1)*ND(2)+NNKA(2))	  ! k 1 prot
		   ak2p =  AKA(iz1o-1+(iz2o-1)*ND(2)+NNKA(2))	  !-k 2 prot
		   
                   sum(1) =  ro1n+ro2n
                   sum(2) =  ro1n-ro2n
                   sum(3) =  ro1p+ro2p
                   sum(4) =  ro1p-ro2p
                   sum(5) =  ak1n+ak2n
                   sum(6) =  ak1n-ak2n
                   sum(7) =  ak1p+ak2p
                   sum(8) =  ak1p-ak2p
		   
                do ir=1,8
                   irom = irom + 1
                   GROM(irom) = sum(ir)
                end do
		
                end do
              end do
            end if
            end do    ! ny2
          end do      ! ny1
        end do        ! nx2
      end do          ! nx1
      return
      end
      Subroutine field2sq(gfield,gam,delta)
      Implicit real*8 (a-h,o-z)
      Dimension gfield(*)
      Dimension gam(*),delta(*)
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC1/ NNU(4),NNV(4),NNRO(4),NNKA(4)
      Common /ITLOC2/ ND2PACK(4),NNGPACK(4),NNDPACK(4)
      
c
c  Transforming gamma to square form
c

       do it=1,3,2
          
	  IGP = NNGPACK(it)
	  IGM = NNGPACK(it+1)
	  IGGP= NNRO(it)
	  IGGM= NNRO(it+1)
c
          call pack2sq(gfield(IGP),gfield(IGM),Gam(IGGP),Gam(IGGM),1)

       end do
       
       do it=1,3,2
          
	  IGP = NNGPACK(it)   + ND2PACK(it)
	  IGM = NNGPACK(it+1) + ND2PACK(it+1)
	  IGGP= NNRO(it)      + ND2(it)
	  IGGM= NNRO(it+1)    + ND2(it+1)
c
          call pack2sq(gfield(IGP),gfield(IGM),Gam(IGGP),Gam(IGGM),1)

       end do
c
c  Transforming delta to square form
c
       do it=1,3,2
          
	  IG1P = NNDPACK(it)
	  IG1M = NNDPACK(it+1)
	  IGG1P= NNKA(it)
	  IGG1M= NNKA(it+1)
	  
	  IG2P = NNDPACK(it)     + ND2PACK(it)
	  IG2M = NNDPACK(it+1)	 + ND2PACK(it+1)
c
          call de2sq(gfield(IG1P),gfield(IG2P),gfield(IG1M),gfield(IG2M)
     *	    ,Delta(IGG1P),Delta(IGG1M))
     
       end do
       return
       end
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
      Subroutine ReadData
      Implicit real*8 (A-H,O-Z)
      Logical vcom,vcom2
      Character*2 cnuc
      Character*8 TEE
      Parameter (NCM=06,NLAM=NCM-3)  ! maximum number of constraints
      Parameter (NEEE=30) ! dimension of the output matrix
      Common /CCONS/ cval(ncm),prec(ncm),IC(NCM),
     & ILAMB(NLAM),IMM(NLAM),ISOS(NCM,2)
      Common /CITER / epsg,etamax,etamin,dmax,dmin,tshh,maxiter
      Common /HFBOPT/ Acom,vcom,vcom2,icouech
      Common/NECKCO/ aneck,rneckp,rneckn
      Common /STARTO/ iwf,iol
      Common /CDEB  / ideb(9)
      Common /OUTLEV/ ioutlev
      Common /OSCLEN/ bx,by,bz
      Common /CNECK / z00,a00
      Common /EEEEEE/ TEE(NEEE),EEE(7,NEEE)
      Common /IIIEEE/ mkin,mehf,mepair,mehfb,medd,merea,mecou,
     &                mn,mjx,msx,mq20,mq22,mq40,mneck,mr2,
     &                mbet2,mgamm,mbet4,maxrat,
     &                mfn,mjx2,mjy2,mjz2
      aneck = 2.0d+00
c
      do i=1,9
         ideb(i) = 0
      end do
c
      ideb(8) = 1 ! 1 = no spin-orbit pairing
      
      IC(1) = 1
      IC(2) = 1
c
      do i=1,ncm
         ISOS(i,1) = 1
         ISOS(i,2) = 1
      end do

      ISOS(1,2) = 0
      ISOS(2,1) = 0
c      
c      ISOS(7,2) = 0
c      ISOS(8,1) = 0
c
      read(5,1) ia,cnuc,iz,icom,icom2,icouech
    1 format (11x,I3,1x,a2,6x,i3,25x,i1,5x,i1,9x,i1)
c
      write(6,101)ia,cnuc,iz,icom,icom2,icouech
  101 format (2x,10('========'),//,30x,'T R I A X I A L     H F B ',//,
     &        32x,'NUCLEUS ',i4,2x,a2,' Z= ',i4,//,30x,
     &        'COM ',i1,' COM2 ',i1,' COU-ECH ',i1,//,
     &        2x,10('========'),//)
      in = ia-iz
      cval(1) = dfloat(iz)
      cval(2) = dfloat(in)
      Acom    = dfloat(ia)
      vcom    = icom .ne.0
      vcom2   = icom2.ne.0
c
      read(5,2) epsg,maxiter,ioutlev,itstg
    2 format(4x,f9.6,9x,i5,20x,i1,11x,i1)
      ideb(9) = itstg
      write(6,102) epsg,maxiter
  102 format (30x,'EPSG ',f9.6,' MAXITER ',i5)
c
      read(5,3) etamax,etamin,dmax,dmin,tshh
    3 format( 5x,f7.4,6x,f7.4,5x,f5.2,5x,f5.2,5x,f6.1)
c
      read(5,4) ipar
    4 format(19x,i1)
c+--------------------------------------------------------------------+
c|+------------------------------------------------------------------+|
c||   ipar = 0     D1S                                               ||
c||   ipar = 1     D1             G O G N Y    F O R C E             ||
c||   ipar = 2     D1'                                               ||
c|+------------------------------------------------------------------+|
c+--------------------------------------------------------------------+
      write(6,112)
c      
      call Gognyf(ipar)
c
      write(6,112)
c      
      read(5,4) iwf
      read(5,5) iol,bx,by,bz
    5 format(19x,i1,26x,f7.4,3x,f7.4,3x,f7.4)
      read(5,6) ic(3),cval(3),prec(3)
    6 format( /,8x,i1,6x,d16.8,5x,d10.2)
      if(ic(3).ne.0) then 
             cval(3) = dsqrt(cval(3)*(cval(3)+1.0d+00))
             write(6,106) cval(3),prec(3)
      end if
  106 format( //,15x,' CONSTRAINT ON <Jx>   ',d16.8,2x,d10.2)
      read(5,7) ILAMB(1),IMM(1),ic(4),cval(4),prec(4)
      if(ic(4).ne.0) write(6,107) ILAMB(1),IMM(1),cval(4),prec(4)
  107 format( /,15x,' CONSTRAINT ON Q',2i1,'    ',d16.8,2x,d10.2)
    7 format( 1x,2i1,5x,i1,6x,d16.8,5x,d10.2)
      read(5,7) ILAMB(2),IMM(2),ic(5),cval(5),prec(5)       
      if(ic(5).ne.0) write(6,108) ILAMB(2),IMM(2),cval(5),prec(5)
  108 format( /,15x,' CONSTRAINT ON Q',2i1,'    ',d16.8,2x,d10.2)
      read(5,7) ILAMB(3),IMM(3),ic(6),cval(6),prec(6)       
      if(ic(6).ne.0) write(6,109) ILAMB(3),IMM(3),cval(6),prec(6)
  109 format( /,15x,' CONSTRAINT ON Q',2i1,'    ',d16.8,2x,d10.2)
c
      write(6,112)
  112 format(//,2x,10('========'),//)
      mkin    = 1
      TEE(mkin  ) = 'Kinetic '
      mehf    = 2
      TEE(mehf  ) = 'HF Ener '
      mepair  = 3
      TEE(mepair) = 'Pairing '
      mehfb   = 4
      TEE(mehfb ) = 'HFB Ener'
      medd    = 5
      TEE(medd  ) = 'DD Ener '
      merea   = 6
      TEE(merea ) = 'Rearrang'
      mecou   = 7
      TEE(mecou ) = 'Coul Ex.'
      mn      = 8
      TEE(mn    ) = '   N    '
      mjx     = 9
      TEE(mjx   ) = '  Jx    '
      msx     = 10
      TEE(msx   ) = '  Sx    '
      mq20    = 11
      TEE(mq20  ) = ' Q20(b) '
      mq22    = 12
      TEE(mq22  ) = ' Q22(b) '
      mq40    = 13
      TEE(mq40  ) = ' Q40(b2)'
      mneck   = 14
      TEE(mneck ) = '  Neck  '
      mr2     = 15
      TEE(mr2   ) = 'MS Rad  '
      mbet2   = 16
      TEE(mbet2 ) = ' Beta 2 '
      mgamm   = 17
      TEE(mgamm ) = ' Gamma  '
      mbet4   = 18
      TEE(mbet4 ) = ' Beta 4 '
      maxrat  = 19
      TEE(maxrat) = '  z/x   '
      mfn     = 20
      TEE(mfn   ) = ' <N^2>  '
      mjx2    = 21
      TEE(mjx2  ) = '<Jx**2> '
      mjy2    = 22
      TEE(mjy2  ) = '<Jy**2> '
      mjz2    = 23
      TEE(mjz2  ) = '<Jz**2> '
      
      TEE(24    ) = ' KIN2 G '
      TEE(25    ) = ' KIN2 D '
      TEE(26    ) = '        '
      TEE(27    ) = '        '
      TEE(28    ) = '        '
      TEE(29    ) = '        '
      TEE(30    ) = '        '
      return
      end
      Double precision  function second()
      Implicit real*8 (A-H,O-Z)
c      i1= mclock()
c   mclock  ------- RS/6000 xlf fortran compiler
c
c     i1 = 0
c
c not known yet for the NDP compiler
c
c      second = secnds(0.0)
      second = 0.0
      return
      end
