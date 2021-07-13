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
      
