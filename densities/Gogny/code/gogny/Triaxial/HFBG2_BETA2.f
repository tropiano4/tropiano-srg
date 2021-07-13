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
         AKA=0.000001d0
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
