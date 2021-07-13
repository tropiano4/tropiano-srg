c Set of subroutines calling the functions defined in TALMANF.f
      Subroutine TZ1DIM(nzmax,maxtz1,tz1d,itz1d)
      Implicit real*8 (A-H,O-Z)
c ---------------------------------------------------<< Start    include
      Include 'COMDIM'
c ---------------------------------------------------<< End      include
      Dimension tz1d(maxtz1)
      Dimension itz1d(nzmax,nzmax)
c
      Common /fact/dfact(Nfac),ddfact(Nfac),Nfacm
c
      if(nfac.ne.nfacm) then
        write(6,*) ' ERROR IN TZ1DIM:  NFAC,NFACM',NFAC,NFACM
        stop
      end if
c
      indx = 0
c
      do 1 in =1,nzmax
           do 2 im=in,nzmax
                imumin = im-in+1
                imumax = im+in-1
                itz1d(in,im) = indx
                itz1d(im,in) = indx
                do 3 imu =imumin,imumax,2
                     indx=indx+1
                     tz1d(indx)=DT1(in-1,im-1,imu-1)
    3           continue
    2      continue
    1 continue
      if((indx).gt.maxtz1) call errout('ITZ > MAXTZ1    ','TZ1DIM  ',4)
      return
      end
c+---------------------------------------------------------------------+
c|    Subroutine J1BB: Computes the J(n,m,mu) coeficients for the      |
c|    Brink-Boecker interaction.                                       |
c|                                                                     |
c|     - m >= n  and mu = 0,...,2*(Nzmax-1)                            |
c|     - mu must have the same parity of n+m                           |
c+---------------------------------------------------------------------+
c|    Parameters:                                                      |
c|                                                                     |
c|           - Nzmax ..... Maximum value of nz + 1                     |
c|           - Maxjz1 .... Number of elements in J                     |
c|           - Maxtz1 .... Number of elements in T                     |
c|           - Tz1d ...... Vector containing T                         |
c|           - Itz1d ..... Index of T                                  |
c|           - Work ...... Scratch vector                              |
c|           - Nwork ..... Dimension of Work                           |
c|                                                                     |
c|           - Amu(2) .... Interaction lenghts mu1 and  mu2            |
c|           - b ......... Oscillator Parameter                        |
c|                                                                     |
c|           - Ajz1d ..... Output vector containing J                  |
c|           - ijz1d ..... Index of J                                  |
c+---------------------------------------------------------------------+
c|                                             |                       |
c|    Dependencies:  Common /Fact/             |     SETUPN            |
c|                                             |                       |
c+---------------------------------------------------------------------+
      Subroutine J1BB (nzmax,maxjz1,maxtz1,tz1d,itz1d,nmax,work,
     * amu,b,ajz1d,ijz1d)
      Implicit Real*8 (A-H,O-Z)
c ---------------------------------------------------<< Start    include
      Include 'COMDIM'
c ---------------------------------------------------<< End      include
      Dimension ajz1d(2,maxjz1),tz1d(maxtz1)
      Dimension ijz1d(nzmax,nzmax),itz1d(nmax,nmax)
      Dimension amu(2)
c ----------------------------------------------------------------------
      Common /fact/dfact(Nfac),ddfact(Nfac),Nfacm
c ----------------------------------------------------------------------
      Maxmu = 2*nzmax - 1
      Minmu = 1
c ----------------------------------------------------------------------
      indx = 0
c
      do 1 in =1,nzmax
           do 2 im=in,nzmax
                ijz1d(in,im) = indx
                ijz1d(im,in) = indx
                inumin = im-in+1
                inumax = im+in-1
                Minmu  = Mod(in+im,2)+1
                do 3 imu = Minmu,Maxmu,2
                indx = indx + 1
                ajz1d(1,indx) = DJ1BB(in-1,im-1,imu-1,amu(1)/b)
                ajz1d(2,indx) = DJ1BB(in-1,im-1,imu-1,amu(2)/b)
    3           continue
    2      continue
    1 continue
      if((indx).gt.maxjz1) call errout('ITZ > MAXJZ1    ',' J1BB   ',4)
      return
      end
c+---------------------------------------------------------------------+
c³    Subroutine LSOMU: Computes the L(n,m,mu) coeficients for the     ³
c³    Spin-Orbin matrix elements.                                      ³
c³                                                                     ³
c³     - m =0,...,NZMAX                                                ³
c³     - m >= n  and mu = 0,...,2* Nzmax                               ³
c³     - mu must have the same parity of n+m                           ³
c+---------------------------------------------------------------------+
c³    Parameters:                                                      ³
c³                                                                     ³
c³           - Nzmax1..... Maximum value of nz + 2                     ³
c³           - Maxl1 ..... Number of elements in L                     ³
c³           - Maxtz2 .... Number of elements in T                     ³
c³           - Tz2d ...... Vector containing T                         ³
c³           - Itz2d ..... Index of T                                  ³
c³           - Work ...... Scratch vector                              ³
c³           - Nwork ..... Dimension of Work                           ³
c³                                                                     ³
c³                                                                     ³
c³           - Al1d ...... Output vector containing J                  ³
c³           - il1d ...... Index of J                                  ³
c+---------------------------------------------------------------------+
c³                                             ³                       ³
c³    Dependencies:  Common /Fact/             ³     SETUPN            ³
c³                                             ³                       ³
c+---------------------------------------------------------------------+
      Subroutine LSOMU
     *(nzmax1,maxl1,maxtz2,tz2d,itz2d,work,al1d,il1d)
c
      Implicit Real*8 (A-H,O-Z)
c ---------------------------------------------------<< Start    include
      Include 'COMDIM'
c ---------------------------------------------------<< End      include
      Dimension al1d(maxl1),tz2d(maxtz2)
      Dimension il1d(nzmax1,nzmax1),itz2d(nzmax1,nzmax1)
c ----------------------------------------------------------------------
      Common /fact/dfact(Nfac),ddfact(Nfac),Nfacm
c ----------------------------------------------------------------------
      Maxmu = 2*nzmax1 - 1
      Minmu = 1
c ----------------------------------------------------------------------
      indx = 0
c
      do 1 in =1,nzmax1
           do 2 im=in,nzmax1
                il1d(in,im) = indx + 1
                il1d(im,in) = indx + 1
                Minmu  = Mod(in+im,2)+1
                do 3 imu = Minmu,Maxmu,2
                indx = indx + 1
                al1d(indx) = DL1(in-1,im-1,imu-1)
    3           continue
    2      continue
    1 continue
      if((indx).gt.maxl1) call errout(' IL > MAXL1     ',' LSOMU  ',4)
      return
      end
c+---------------------------------------------------------------------+
c³    Subroutine DSOMU: Computes the D(n|m,mu) coeficients for the     ³
c³    Spin-Orbin matrix elements.                                      ³
c³                                                                     ³
c³     - m =0,...,NZMAX-1                                              ³
c³     - m >= n  and mu = 0,...,2* NZMAX                               ³
c³     - mu must have opposite parity to n+m                           ³
c+---------------------------------------------------------------------+
c³    Parameters:                                                      ³
c³                                                                     ³
c³           - Nzmax ..... Maximum value of nz + 1                     ³
c³           - Nzmax1..... Maximum value of nz + 2                     ³
c³           - Maxl1 ..... Number of elements in L                     ³
c³           - MaxD1 ..... Number of elements in D                     ³
c³           - Al1d ...... Vector containing L                         ³
c³           - Il1d ...... Index of L                                  ³
c³                                                                     ³
c³                                                                     ³
c³           - DSOmu...... Output vector containing D                  ³
c³           - idsomu .... Index of D                                  ³
c+---------------------------------------------------------------------+
c³                                             ³                       ³
c³    Dependencies:  Subroutine LSOMU          ³                       ³
c³                                             ³                       ³
c+---------------------------------------------------------------------+
c not modified
      Subroutine DSOMU
     *(nzmax,nzmax1,maxl1,maxd1,al1d,il1d,DSOmuv,iDSOmu)
c
      Implicit Real*8 (A-H,O-Z)
      Logical v1
      Dimension al1d(maxl1),DSOmuv(maxD1)
      Dimension il1d(nzmax1,nzmax1),idsomu(nzmax,nzmax)
c
      indx   = 0
      Maxmu = 2*NZMAX+1
c
      do 1 in =1,nzmax
           v1 = in.ne.1
           fac1 = dsqrt(2.0d+00*dfloat(in-1))
           fac2 = dsqrt(2.0d+00*dfloat(in))
           do 2 im=1,nzmax
                idsomu(in,im) = indx + 1
                if(v1) indt1  = il1d(im,in-1)
                indt2         = il1d(im,in+1)
                Minmu  = 2-Mod(in+im,2)
                do 3 imu = Minmu,Maxmu,2
                   sum = -fac2*al1d(indt2)
                   if(v1) sum = sum + fac1*al1d(indt1)
                   indt1 = indt1+1
                   indt2 = indt2+1
                   indx = indx + 1
                   DSOmuv(indx) = sum
    3           continue
    2      continue
    1 continue
      if((indx).gt.maxD1) call errout(' IL > MAXD1     ',' DSOMU  ',4)
      return
      end
c+---------------------------------------------------------------------+
c
c   February 20, 1992
c
c
c+---------------------------------------------------------------------+
c   Subroutine SXTJ               __
c                                \
c     Computes the contraction    >   T(k1,iq1,imu) J(iq2,k2,imu)
c                                /__
c                                imu
c
c     For k1 >= iq1 and k2 >= iq2
c
c     The vector IXTJ(k1,iq1) gives the position of the above contraction
c     for a given k1,iq1 in the XTJ vector.
c
c     The dimension of XTJ is
c                               2
c     ((nzmax2+1)(nzmax-nzmax2))  +
c                                                   2
c     ( nzmax(nzmax+1)/2 -(nzmax2+1)(nzmax-nzmax2) )
c
c     where nzmax2 = integer_part_of( nzmax/2 )
c
c   The T(k1,iq1,mu) is assumed to be computed with a maximum number
c   of k1,etc equal to nmax
c
c
c+---------------------------------------------------------------------+
      Subroutine Sxtj(NZMAX,nmax,MAXTZ1,TZ1D,ITZ1D,MAXJZ1D,AJZ1D,IJZ1D,
     x                MAXXTJ,XTJ)
      Implicit real*8 (A-H,O-Z)
      Dimension tz1d(maxtz1),ajz1d(2,maxjz1d),xtj(2,maxxtj)
      Dimension itz1d(nmax,nmax),ijz1d(nzmax,nzmax)
c
      ic = 0
      do 10 k1=1,nzmax
         do 11 iq1 = 1,k1
         im1 = Mod((k1+iq1),2)
         it1 = itz1d (k1,iq1)
         imumax = k1+iq1-1
         imumin = k1-iq1+1
         minmu  = im1 +1             ! the starting imu for J
         ish    = (imumin-minmu)/2   ! the shift in do 30
            do 20 k2=1,nzmax
               do 21 iq2=1,k2
               im2 = Mod((k2+iq2),2)
               if (im1.ne.im2) goto 21  ! This can be done in do 21
               ij2 = ijz1d(k2,iq2)+ish
               it2 = it1
               sum1 = 0.0d+00
               sum2 = 0.0d+00
                  do 30 imu = imumin,imumax,2
                     it2 = it2 + 1
                     ij2 = ij2 + 1
                     sum1= sum1+ tz1d(it2)*ajz1d(1,ij2)
                     sum2= sum2+ tz1d(it2)*ajz1d(2,ij2)
   30             continue
               ic  = ic + 1
               xtj(1,ic) = sum1
               xtj(2,ic) = sum2
   21          continue
   20       continue
   11    continue
   10 continue
      if(ic.gt.maxxtj) call errout('IC > MAXXTJ     ','SXTJ    ',4)
      return
      end
c
c   February 17, 1992
c
c
c   Subroutine SXTJP              __
c                       k1+iq1   \
c     Computes      (-1)          >   T(k1,iq1,imu) J(iq2,k2,imu)
c                                /__
c                                imu
c
c     For k1 >= iq1 and k2 >= iq2
c
c     This subroutine is needed for the exchange and pairing part
c
c     The vector IXTJ(k1,iq1) gives the position of the above contraction
c     for a given k1,iq1 in the XTJ vector.
c
c     The dimension of XTJ is
c                               2
c     ((nzmax2+1)(nzmax-nzmax2))  +
c                                                   2
c     ( nzmax(nzmax+1)/2 -(nzmax2+1)(nzmax-nzmax2) )
c
c     where nzmax2 = integer_part_of( nzmax/2 )
c
c
c   The T(k1,iq1,mu) is assumed to be computed with a maximum number
c   of k1,etc equal to nmax
c
c
      Subroutine Sxtjp(NZMAX,nmax,MAXTZ1,TZ1D,ITZ1D,MAXJZ1D,AJZ1D,IJZ1D,
     x                MAXXTJ,XTJ)
      Implicit real*8 (A-H,O-Z)
      Dimension tz1d(maxtz1),ajz1d(2,maxjz1d),xtj(2,maxxtj)
      Dimension itz1d(nmax,nmax),ijz1d(nzmax,nzmax)
c
      ic = 0
      do 10 k1=1,nzmax
         do 11 iq1 = 1,k1
         im1 = Mod((k1+iq1),2)
         it1 = itz1d (k1,iq1)
         phas= dfloat(1-2*im1)
         imumax = k1+iq1-1
         imumin = k1-iq1+1
         minmu  = im1 +1             ! the starting imu for J
         ish    = (imumin-minmu)/2   ! the shift in do 30
            do 20 k2=1,nzmax
               do 21 iq2=1,k2
               im2 = Mod((k2+iq2),2)
               if (im1.ne.im2) goto 21  ! This can be done in do 21
               ij2 = ijz1d(k2,iq2)+ish
               it2 = it1
               sum1= 0.0d+00
               sum2= 0.0d+00
                  do 30 imu = imumin,imumax,2
                     it2 = it2 + 1
                     ij2 = ij2 + 1
                     sum1= sum1+ tz1d(it2)*ajz1d(1,ij2)
                     sum2= sum2+ tz1d(it2)*ajz1d(2,ij2)
   30             continue
               ic  = ic + 1
               xtj(1,ic) = sum1* phas
               xtj(2,ic) = sum2* phas
   21          continue
   20       continue
   11    continue
   10 continue
      if(ic.gt.maxxtj) call errout('IC > MAXXTJ     ','SXTJ    ',4)
      return
      end
c
c   February 14, 1992
c
c   Subroutine SIXTJ:
c
c   Computes two integer matrices needed to locate the elements of XTJ
c   easily.
c
c     The matrix IXTJ(k1,iq1) gives the position in the vector XTJ
c     where the k2,iq2 values start for a given k1,iq1.
c
c     The matrix IXTJ2(k2,iq2) gives the relative position (with
c     respect to the starting point ixtj(k1,iq1) of the element in
c     XTJ
c                             BB
c     The position in XTJ of I                is then given by
c                             k1 iq2 iq1 k2
c
c     IXTJ(k1,iq1)+IXTJ2(k2,iq2)
c
c
      Subroutine Sixtj(NZMAX,IXTJ,IXTJ2)
      Implicit real*8 (A-H,O-Z)
      Dimension IXTJ(nzmax,nzmax),IXTJ2(nzmax,nzmax)
c
      ic = 0
      ice= 0
      ico= 0
      do 10 k1=1,nzmax
         do 11 iq1 = 1,k1
         im1 = Mod((k1+iq1),2)
         if(im1.eq.0) then
            ixtj2(k1,iq1) = ice
            ixtj2(iq1,k1) = ice
            ice = ice + 1
         else
            ixtj2(k1,iq1) = ico
            ixtj2(iq1,k1) = ico
            ico = ico + 1
         end if
         ixtj(k1,iq1) = ic + 1
         ixtj(iq1,k1) = ic + 1
            do 20 k2=1,nzmax
               do 21 iq2=1,k2
               im2 = Mod((k2+iq2),2)
               if (im1.ne.im2) goto 21  ! This can be done in do 21
               ic  = ic + 1
   21          continue
   20       continue
   11    continue
   10 continue
      return
      end
