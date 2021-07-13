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
