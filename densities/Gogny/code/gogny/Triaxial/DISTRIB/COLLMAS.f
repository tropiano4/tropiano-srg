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
