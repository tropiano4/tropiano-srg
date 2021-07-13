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
