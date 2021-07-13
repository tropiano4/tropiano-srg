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
