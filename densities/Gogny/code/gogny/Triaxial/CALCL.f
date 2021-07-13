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
