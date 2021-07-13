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
