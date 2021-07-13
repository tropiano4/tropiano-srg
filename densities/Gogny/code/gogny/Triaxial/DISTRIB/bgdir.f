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
      fac = 1.0d+00/(bx*by*bz*2.0d+00*dsqrt(2.0d+00))
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
