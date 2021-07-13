c+---------------------------------------------------------------------+
c    Version 4.0
c
c    This program is used to compute the dimensions needed in the
c    triaxial code and also to create the include files.
c+---------------------------------------------------------------------+
c    The basis is determined by using the energy condition
c
c       a  n    + a  n   + a  n   <  N
c        x  x      y  y     z  z  =   0
c
c
c          1/3   -1/3            2/3   1/3           -1/3  -2/3
c    a  = q    p           a  = p    q         a  = p     q
c     x                     y                   z
c
c                                         Rz           Ry
c     q and p are the axis ratios     q = ---  and p = ---
c                                         Rx           Rx
c
c
c+---------------------------------------------------------------------+
c     In the particular case of p=1 the condition reads
c
c          1                    N0
c         --- n  + (n +n ) <  ------
c          q   z     x  y  =  q^1/3
c
c      which is the same as in the axial case but with 
c
c          N  (triaxial) = q^1/3 N  (AXIAL)
c           0                     0
c+---------------------------------------------------------------------+
c
      Implicit real*8 (A-H,O-Z)
      Parameter (maxdim=30)
      Logical LMZE(maxdim,maxdim),LMZO(maxdim,maxdim)
      Character*1 sep(maxdim)
      Dimension MY(maxdim),MZ(maxdim,maxdim)
      Dimension NZIE(maxdim,maxdim),NZIO(maxdim,maxdim)
c
      do i=1,maxdim
         sep(I)=','
      end do
c
      write(6,*) ' Enter q p and N0  '
      read (5,*) q,p,rn0
c
      open(unit=7,file='DIMPERM',status='NEW')
      open(unit=8,file='COMDIM',status='NEW')
      open(unit=9,file='DIMTRIAX',status='NEW')
c

      onethird = 1.d+00/3.d+00
      twothird = 2.d+00/3.d+00
      ax= q**onethird /(p**onethird)
      ay= p**twothird *(q**onethird)
      az= 1./(p**onethird * q**twothird)
c
      nxmax = rn0/ax
      nxmax = nxmax + 1
      nymax = rn0/ay
      nymax = nymax + 1
      nzmax = rn0/az
      nzmax = nzmax + 1
c
      sep(nxmax)='/'

      Do nx=1,nxmax
         temp   = rn0-(nx-1)*ax
         itemp  = temp/ay
         MY(nx) = itemp + 1
      end do
c
      Do nx=1,nxmax
         nym= MY(nx)
         Do ny=1,nym
            temp  = rn0 -(nx-1)*ax - (ny-1)*ay
            itemp = temp/az
            MZ(nx,ny) = itemp + 1
         end do
      end do
c
      nmax  = Max(nxmax,nymax,nzmax)
      nmax1 = nmax+1
c
      nherm = 2*nmax1
      nlegn = 20
c
c    Computing  LMZE LMZO NZIE NZIO
c
      do nx=1,nxmax
        nym = my(nx)
          do ny=1,nym
             nxy  = nx +ny
             nzie(nx,ny) = 1 + Mod(nxy,2)
             nzio(nx,ny) = 2 - Mod(nxy,2)
             nzm = mz(nx,ny)
             lmze(nx,ny) = .false.
             lmzo(nx,ny) = .false.
             if(nzm.ge.nzie(nx,ny)) lmze(nx,ny) = .true.
             if(nzm.ge.nzio(nx,ny)) lmzo(nx,ny) = .true.
          end do
      end do
c
      iz1e = 0
      iz1o = 0
      jz1e = 0
      jz1o = 0
      do nx1=1,nxmax
        nym1 = my(nx1)
          do ny1 = 1,nym1
            nxy1 = nx1+ny1
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = mz(nx1,ny1)
            if(lmze(nx1,ny1)) then
              iz1e=iz1e+(nzm1-nzi1e)/2+1  ! pos parity
            end if
            if(lmzo(nx1,ny1)) then
              iz1o=iz1o+(nzm1-nzi1o)/2+1  ! neg parity
            end if
          end do ! ny1
      end do     ! nx1
c
      NP = iz1e     ! Positive parity matrices dimension NPxNP
      NM = iz1o     ! Negative parity matrices dimension NMxNM
      NPM= Max(NP,NM)
c
      nrop = (np*(np+1))/2   !   Ro dimensions    Pos parity
      nrom = (nm*(nm+1))/2   !                    Neg parity
      nrop8 = nrop*8
      nrom8 = nrom*8
      ngp   = nrop
      ngm   = nrom
      ngp8  = nrop8
      ngm8  = nrom8
c+--------------------------------------------------------------------+
c|     Permanent storage                                              |
c|                                     Direct                         |
c+--------------------------------------------------------------------+
      Maxtz1 = (Nmax*(Nmax**2 + 3*Nmax + 2))/6
      Maxtz2 = ((Nmax+1)*((Nmax+1)**2 + 3*Nmax + 5))/6
      Maxjx1 = ((Nxmax+1)**2/4)*Nxmax + (Nxmax**2/4)*(Nxmax-1)
      Maxjy1 = ((Nymax+1)**2/4)*Nymax + (Nymax**2/4)*(Nymax-1)
      Maxjz1 = ((Nzmax+1)**2/4)*Nzmax + (Nzmax**2/4)*(Nzmax-1)
      Maxlso = ((Nmax+2)**2/4)*(Nmax+1) + ((Nmax+1)**2/4)*Nmax
      Maxdso = ((Nmax+1)**2/2-Nmax)*(Nmax+1)+(Nmax**2/2)*Nmax
      Nacou  = (2*Nxmax-1)*(2*Nymax-1)*(2*Nzmax-1)
      Ncoulm = Nmax**2
      nwf2   = 2*Nmax-1
      Nwavef = Nmax*Nherm
      Nwavefd= nwf2*Nherm
c+--------------------------------------------------------------------+
c|                                     Exchange (SXTJ)                |
c+--------------------------------------------------------------------+
      nxmax2 = nxmax/2
      nymax2 = nymax/2
      nzmax2 = nzmax/2
      maxtjx = ((nxmax2+1)*(nxmax-nxmax2))**2+
     *         ( (nxmax*(nxmax+1))/2-(nxmax2+1)*(nxmax-nxmax2) )**2
      maxtjy = ((nymax2+1)*(nymax-nymax2))**2+
     *         ( (nymax*(nymax+1))/2-(nymax2+1)*(nymax-nymax2) )**2
      maxtjz = ((nzmax2+1)*(nzmax-nzmax2))**2+
     *         ( (nzmax*(nzmax+1))/2-(nzmax2+1)*(nzmax-nzmax2) )**2
      maxx   = 4*maxtjx+2*maxtjy+2*maxtjz
c
      NPERTOT = Maxtz1 + Maxtz2 + 2*Maxjx1 + 2*Maxjy1 + 2*Maxjz1 +
     *          Maxlso + Maxdso + Nacou    + Ncoulm   + Nwavef   +
     *          Nwavefd+ 2*Nherm+  2*Nlegn + maxx
c+--------------------------------------------------------------------+
c|   Temporary                                                        |
c|                                                  Direct            |
c+--------------------------------------------------------------------+
      NDMU    = nxmax*nymax*nzmax
      NDTHETA = Ndmu*(nzmax/2+1)
      nherm38 = 8*nherm**3
      irorz  = Nherm*4
      irory  = Nherm*Nherm*6
      irortmp= 2*Nherm**3
      idmus0 = nherm**2 * nmax
      idmus1 = nherm    * nmax**2
      jdmus0 = 15*idmus0
      jdmus1 = 15*idmus1
c
      NDIRT1  = 15*NDMU + nherm38 + irorz + irory
      NDIRT2  = 15*NDMU + nherm38 + irortmp + jdmus0 + jdmus1
      NDIRT3  = 61 * Ndmu + 20 * NDTHETA
      NDIRT   = MAX(NDIRT1,NDIRT2,NDIRT3)
      write(6,*) ' TEMP DIM DIR ',NDIRT1,NDIRT2,NDIRT3,NDIRT
c+--------------------------------------------------------------------+
c|                                             Exchange and pairing   |
c+--------------------------------------------------------------------+
      mazopti= 16*nzmax*nzmax
      mayopti= 16*nzmax**2 * nymax**2
c
      kk = 0
      do ny1=1,nymax
        nz1m = MZ(1,ny1)
        do ny2=1,nymax
           nz2m = MZ(1,ny2)
           kk = kk + nz1m*nz2m
        end do
      end do
      kk = kk * 16
      write(6,*) ' mayopti ',mayopti,kk
c
      NECHT  = mazopti + mayopti
      NTEMP  = MAX(NDIRT,NECHT)
c

      NUVP = np*np
      NUVM = nm*nm
      nropt = 16*nrop
      nromt = 16*nrom
      nuvtp = 16*nuvp
      nuvtm = 16*nuvm
      nuvtt = 8*(nuvp+nuvm)
c
      rMb    = 1.0d+00/dfloat(1024*1024)
      rMbropt = dfloat(nropt*8)*rMb
      rMbromt = dfloat(nromt*8)*rMb
      rMBnuvpt= dfloat(nuvtp*8)*rMb
      rMBnuvmt= dfloat(nuvtm*8)*rMb
      rmaz   = 8.0d+00*Dfloat(mazopti)*rMB
      rmay   = 8.0d+00*Dfloat(mayopti)*rMb
c
      NTEMP  = MAX(NTEMP,NUVTT)
c
      write(6,1001) np,nm
      write(6,1002) nrop,nrom
      write(6,1005) rmbropt,rmbromt
      write(6,1006) mazopti,rmaz,mayopti,rmay
      write(6,1007) NDIRT,NECHT,NUVTT,NTEMP
      write(6,1008) nuvtp,rmbnuvpt,nuvtm,rmbnuvmt
      write(6,1010) NUVTT,nrop8,nrom8,ntemp
c
c     Dimperm
c
      write(7,100) nmax,nmax1,nxmax,nymax,nzmax,
     * nwf2,maxtz1,maxtz2,maxjx1,maxjy1,maxjz1,maxlso,
     * maxdso,ndmu,ndtheta,nacou,nherm,nlegn,maxtjx,maxtjy,maxtjz
  100 format ( '       Parameter (imax=',i2,',imax1=',i2,
     * ',ixmax=',i2,',iymax=',i2,',izmax=',i2,/,
     *5x,'* ,iwf2 =',i3,',iaxtz1d=',i6,',iaxtz2d=',i6,',iaxjx1d=',i6,/,
     *5x,'* ,iaxjy1d=',i6,',iaxjz1d=',i6,',iaxlso=',i6,',iaxdso=',i6,/,
     *5x,'* ,idmu=',i8,',idthet=',i8,',iacou=',i6,',iherm=',i2,/,
     *5x,'* ,ilegn= ',i2,',iaxtjx=',i8,',iaxtjy=',i8,',iaxtjz=',i8,')')
c
      write(7,101)
  101 format( '       Dimension iMY(ixmax),iMZ(ixmax,iymax)')
      write(7,2000) (my(k),sep(k),k=1,nxmax)
      write(7,2001) (mz(k,1),k=1,nxmax)
      do j=2,(nymax-1)
          write(7,2002) (mz(k,j),k=1,nxmax)
      end do
          write(7,2003) (mz(k,nymax),sep(k),k=1,nxmax)
c
c   Comdim
c
      izsrc = nzmax
      iyzsrc= nzmax*nymax
      Nfac  = 100
      iqne  = NP
      iqno  = NM
      maxp  = 11
      Kdim  = 30
c
      write(8,200) nmax,nmax1,nxmax,nymax,nzmax,izsrc,iyzsrc,nfac,
     *   iqne,iqno,maxp,kdim,q,p,rn0
  200 format ( '       Parameter (imax=',i2,',imax1=',i2,
     * ',ixmax=',i2,',iymax=',i2,',izmax=',i2,/,
     * 5x'*  ,izsrc=',i2,',iyzsrc=',i5,',nfac=',i3,',inqne=',i4,/,
     * 5x'*  ,inqno=',i4,',maxp=',i2,',kdim=',i2,')',/,
     *  '       Parameter (q0inf=',f7.4,',p0inf=',f7.4,',rn0inf=',
     *  f8.4,')')
c
c    DIMTRIAX
c
      write(9,300) NP,NM,NPM,NUVTT,NTEMP,NROP8,NROM8
  300 format('      Parameter (NP=',i3,',',/,
     *       '     *           NM=',i3,',',/,
     *       '     *          NPM=',i3,',',/,
     *       '     *          NUV=',i8,',',/,
     *       '     *         NH20=',i8,',',/,
     *       '     *          IGP=',i8,',',/,
     *       '     *          IGM=',i8,')')

c  -----------------+ F O R M A T S +------------------------
 1001 format( ' Positive parity states ',i5,
     *        ' Negative parity states ',i5)
c
 1002 format( ' Ro matrix dimensions m>= n ',
     *        ' Positive parity ',i8,
     *        ' Negative parity ',i8)
c
 1005 format( ' For the 16 matrices (densities and fields) ',/,
     * ' Positive parity:',f15.7,'Mb',
     * ' Negative parity:',f15.7,'Mb')
c
 1006 format ( ' Exchange subroutine ',/,
     *         ' Z opt vector        ',i10,'= ',f15.7,' Mb',/,
     *         ' Y opt vector        ',i10,'= ',f15.7,' Mb')
 1007 format ( /,' Size of the TEMP vectors ',/,
     *      ' DIRECT ',I10,' EXCH ',I10,' H20  ',I10,' MAXIMUM ',I10)
c
 1008 format( ' For the 16 matrices (U,V,H20,JX20,N20,JX,JZ) ',/,
     * ' Positive parity, a R*8 vector of ',i10,' and ',f15.7,'Mb',/,
     * ' Negative parity, a R*8 vector of ',i10,' and ',f15.7,'Mb',/)
c
 1010 format( '   D I M E N S I O N S  ',/,
     * ' UV     ',i10,/,
     * ' IGP    ',i10,/,
     * ' IGM    ',i10,/,
     * ' NH20V  ',i10,/)
c
 2000 format('      DATA iMY/',20(i2,A1))
 2001 format('      DATA iMZ/',20(i2,','))
 2002 format('     *        ',20(i2,','))
 2003 format('     *        ',20(i2,A1))
      stop
      end
