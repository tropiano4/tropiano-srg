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
      toto=20.73667622931578756154d+00
      fx = toto/(bx*bx)/amass
      fy = toto/(by*by)/amass
      fz = toto/(bz*bz)/amass

c      fx = 20.734863D+00/(bx*bx)/amass
c      fy = 20.734863D+00/(by*by)/amass
c      fz = 20.734863D+00/(bz*bz)/amass
      
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
