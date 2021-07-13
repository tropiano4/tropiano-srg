      Subroutine setfact()
      Implicit real*8 (A-H,O-Z)
      Include 'MAXDIM'
      Common /fact/ dfact(Nfac),ddfact(Nfac),Nfacm
c+--------------------------------------------------------------------+
c|    dfact(I)  = ln( (I-1)! )   |   ln(I!) = dfact (I+1)             |
c|    ddfact(I) = ln( (2I-3)!! ) |   ln(I!!)= ddfact((I+3)/2)         |
c+--------------------------------------------------------------------+
      Nfacm = Nfac
      dfact(1) = 0.0d+00
      dfact(2) = 0.0d+00
      ddfact(1)= 0.0d+00
      ddfact(2)= 0.0d+00
      do ifac =3,nfacm
           dfact (ifac) = dlog(dfloat(ifac-1))  + dfact (ifac-1)
           ddfact(ifac) = dlog(dfloat(2*ifac-3))+ ddfact(ifac-1)
      end do
      return 
      end
c------------------------------------------------------------------------      
      Subroutine setQN()
      Implicit real*8(a-h,o-z)
      Include 'DIMPERM'
      Include 'MAXDIM'
      Common /CQN/    IQMAX,IQNE(3,NPMAX),IQNO(3,NMMAX),IFILL
c
      IFILL = 1
      
      jz1e = 0
      jz1o = 0
      do nx1=1,ixmax
        nym1 = imy(nx1)
          do ny1 = 1,nym1
            nxy1 = nx1+ny1
            nzi1e= 1+Mod(nxy1,2)
            nzi1o= 2-Mod(nxy1,2)
            nzm1 = imz(nx1,ny1)
            if(nzi1e.le.nzm1) then
              do nz1=nzi1e,nzm1,2
                 jz1e = jz1e + 1
                 iqne(1,jz1e) = nx1
                 iqne(2,jz1e) = ny1
                 iqne(3,jz1e) = nz1
              end do
            end if
            if(nzi1o.le.nzm1) then
              do nz1=nzi1o,nzm1,2
                 jz1o = jz1o + 1
                 iqno(1,jz1o) = nx1
                 iqno(2,jz1o) = ny1
                 iqno(3,jz1o) = nz1
              end do
            end if
          end do ! ny1
      end do     ! nx1
c
      IQMAX = MAX(ixmax,iymax,izmax)
      return
      end
c
c
c
      Subroutine setQNK2B(IQNEO,IQNOE)
      Implicit real*8(a-h,o-z)
      Include 'DIMPERM'
      Include 'MAXDIM'
      Dimension IQNEO(3,2,*),IQNOE(3,2,*)
      Common /CQN/    IQMAX,IQNE(3,NPMAX),IQNO(3,NMMAX),IFILL
      Common /DIMS/   NP,NM,NROP,NROM
c
      if(IFILL.ne.1) call setQN()
c
c  iqneo(*,*,i) i=1,NP index in the odd parity states
c
c      
      do i=1,NP
      
      	nxe= iqne(1,i)
      	nye= iqne(2,i)
      	nze= iqne(3,i)
	
	do k=1,3
	  do l=1,2
	     iqneo(k,l,i) = 0
	  end do
	end do
	
	do j=1,NM
	   nxo = iqno(1,j)
	   nyo = iqno(2,j)
	   nzo = iqno(3,j)
	   
      if((nxo.eq.nxe+1).and.(nyo.eq.nye).and.(nzo.eq.nze)) then
         iqneo(1,1,i)=j
      else if((nxo.eq.nxe-1).and.(nyo.eq.nye).and.(nzo.eq.nze)) then
         iqneo(1,2,i)=j
      else if((nxo.eq.nxe).and.(nyo.eq.nye+1).and.(nzo.eq.nze)) then
         iqneo(2,1,i)=j
      else if((nxo.eq.nxe).and.(nyo.eq.nye-1).and.(nzo.eq.nze)) then
         iqneo(2,2,i)=j
      else if((nxo.eq.nxe).and.(nyo.eq.nye).and.(nzo.eq.nze+1)) then
         iqneo(3,1,i)=j
      else if((nxo.eq.nxe).and.(nyo.eq.nye).and.(nzo.eq.nze-1)) then
         iqneo(3,2,i)=j
      end if
         end do
      end do 
      
      do i=1,NM
      
      	nxo= iqno(1,i)
      	nyo= iqno(2,i)
      	nzo= iqno(3,i)
	
	do k=1,3
	  do l=1,2
	     iqnoe(k,l,i) = 0
	  end do
	end do
	
	do j=1,NP
	   nxe = iqne(1,j)
	   nye = iqne(2,j)
	   nze = iqne(3,j)
	   
      if((nxo+1.eq.nxe).and.(nyo.eq.nye).and.(nzo.eq.nze)) then
         iqnoe(1,1,i)=j
      else if((nxo-1.eq.nxe).and.(nyo.eq.nye).and.(nzo.eq.nze)) then
         iqnoe(1,2,i)=j
      else if((nxo.eq.nxe).and.(nyo+1.eq.nye).and.(nzo.eq.nze)) then
         iqnoe(2,1,i)=j
      else if((nxo.eq.nxe).and.(nyo-1.eq.nye).and.(nzo.eq.nze)) then
         iqnoe(2,2,i)=j
      else if((nxo.eq.nxe).and.(nyo.eq.nye).and.(nzo+1.eq.nze)) then
         iqnoe(3,1,i)=j
      else if((nxo.eq.nxe).and.(nyo.eq.nye).and.(nzo-1.eq.nze)) then
         iqnoe(3,2,i)=j
      end if
         end do
      end do 
      
      return
      end
      
c +--------------------------------------------------------------------+
c    Subroutine TIMEIT: To keep track of the CPU time spent in the
c    calculation.
c +--------------------------------------------------------------------+
      Subroutine TIMEIT(iflag,index,text)
      Implicit real*8 (A-H,O-Z)
      Parameter (Nmax=30) ! Maximum number of items to be accounted for
      Character*16 text
      Character*16 texto
      Common /CTIMEIT/ texto(Nmax),prevtim(Nmax),accutim(Nmax),
     *       Ncalls(Nmax),ntot
c
c    index=0     Set storage to 0
c
      if(index.eq.0) then
         do I=1,Nmax
            prevtim(I) = 0.0d+00
            accutim(I) = 0.0d+00
            Ncalls (I) = 0
            texto  (I) = '****************'
         end do
         ntot = 0
      else if ((index.gt.0).and.(index.le.Nmax)) then
c
c    index .ne. 0
c
         if(iflag.eq.0) then
            prevtim(index) = second()
         else if (iflag.eq.1) then
            ttime = second()
            accutim(index) = accutim(index) + ttime - prevtim(index)
            if(ncalls(index).eq.0) texto(index) = text
            ncalls(index) = ncalls(index) + 1
            if(index.gt.ntot) ntot = index
         else
            write(6,*) ' IFLAG value in TIMEIT unknown '
            write(6,*) ' IFLAG ',IFLAG,' INDEX ',INDEX,TEXT
         end if
c
c    index out of range
c
      else
         write(6,*) ' INDEX value in TIMEIT out of range '
         write(6,*) ' IFLAG ',IFLAG,' INDEX ',INDEX,TEXT
      end if
c
      return
      end
c +--------------------------------------------------------------------+
c    Subroutine TIMEREP: To write a report of the cpu time spent in
c    the CHFBJX program
c +--------------------------------------------------------------------+
      Subroutine TIMEREP
      Implicit real*8 (A-H,O-Z)
      Parameter (Nmax=30) ! Maximum number of items to be accounted for
      Character*16 texto
      Common /CTIMEIT/ texto(Nmax),prevtim(Nmax),accutim(Nmax),
     *       Ncalls(Nmax),ntot
      COMMON /GRENZE/ EPSG,SH,EPSN,SN,EPSI,SI,EPSO,SO,SSS,ZZ,LAN,IG,MAXG
      Common/infob/aq0inf,ap0inf,arn0inf,nxmax,nymax,nzmax
c
      nout = 50
      write(nout,100) aq0inf,ap0inf,arn0inf,ig
  100 format( 10x,'  CPU TIME REPORT OF THE CHFBJX PROGRAM ',/,
     *     5x,10('====='),//,
     *    10x,' Q0 = ',f5.2,' P0 = ',f5.2,' N0 = ',f5.2,//,
     *    20x,' NUMBER OF ITERATIONS',I5,//)
      write(nout,101)
  101 format ( 5x,'Subroutine name ',' Number of calls ',
     * ' Total time (s) ',' Time per call (s) ',' % of total ',/,5x,
     * 20('-----'),/)
c
c It is assumed that i=1 is always the total time
c
      timetot = accutim(1)
      do I=1,ntot
         ttime = accutim(i)
         timepc= ttime/dfloat(ncalls(i))
         perc  = ttime/timetot *1.D+2
         write(nout,102) texto(i),ncalls(i),ttime,timepc,perc
      end do
  102 format( 5x,A16,8x,i8,4x,f12.2,7x,f12.5,5x,f5.1)
c
      return
      end
      subroutine errout(Mesage,subr,icode)
c ษออออออออออออออออออออออออออออออออออออออออออออออออออออออออออออออออออออออป
c บ                                                                      บ
c บ                                                                      บ
c บ                                                                      บ
c บ                                                                      บ
c บ                                                                      บ
c บ                                                                      บ
c บ                                                                      บ
c ศออออออออออออออออออออออออออออออออออออออออออออออออออออออออออออออออออออออผ
      character*16 Mesage
      character*8 subr
      integer*4 icode
      common/fin/ ifin
c
      write(6,100) subr,mesage,icode
  100 format( ' From subroutine ',a8,' >>>>>>>    ',a16,/,
     *        ' Code:  ',i4,/)
c
      if(icode.gt.3) then
           write(6,101)
  101      format(' Code greater than 3 .... Abend       ')
           stop
      end if
      if(icode.eq.3) then
           write(6,102)
  102      format(' Code equal 3 ............ Normal end ')
           ifin = 9999
      end if
c
      return
      end
