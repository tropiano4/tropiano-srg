c+---------------------------------------------------------------------+
c|  Computes ro(r), etc                                                |
c|                                                                     |
c|   Ro           only k1>=k2 is stored in GROP  Pos parity            |
c|     k2 k1                               GROM  Neg parity            |
c|                                                                     |
c|   GROP: rosp, roap, rospw, roapw, akasp, akaap, akaspw, akaapw      |
c|                                                                     |
c|   GROM: rosm, roam, rosmw, roamw, akasm, akaam, akasmw, akaamw      |
c|                                                                     |
c|   ROR : ror0p ror0n ror1p ror1n ror2p ror2n ror3p ror3n             |
c|                                   3                              3  |
c| each of them with a nherm3 = nherm  dimension   nherm38 = 8 nherm   |
c|                                                                     |
c|                                                                     |
c+---------------------------------------------------------------------+
      Subroutine denr(GROP,GROM,ror,wfho,wfhonec,rorz,roryz)
c
      Implicit real*8 (A-H,O-Z)
      Logical lcx,lcy,lx12,lxy12
      Logical LMZE,LMZO
      Include 'COMDIM'
      Dimension sum(8)
c
      Dimension GROP(nrop8)
      Dimension GROM(nrom8)
c
      Dimension ror(nherm38),wfho(nmax,nherm),wfhonec(nmax,nherm)
c
      Dimension rorz (irorz) ! Temporary vectors
      Dimension roryz(irorY) !
c
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      Common/DIMEN/kmax,kmax1,kxmax,kymax,kzmax,nwf2,
     * maxtz1,maxjz1,maxjx1,maxjy1,maxlso,maxdso,
     * ndmu,ndthet,nacou,nrop,nrom,nrop8,nrom8,nherm,nherm38,nlegn
c
      Common/DIMPDD/ irortmp,idmus0,jdmus0,idmus1,jdmus1,irorz,irory
c
      Common/OSCLEN/bx,by,bz
      Common/NECKCO/ aneck,rneckp,rneckn
      Common/PULGA/ideb
c
      bm1= 1.0d+00/(bx*by*bz)
      fneck = 2.0d+00/Dsqrt(1.d+00 + (bz/aneck)**2)
      nherm3 = nherm**3
      nherm4 = nherm *4
      nherm2 = nherm**2
      nherm26= nherm**2 * 6
c
      do jher=1,nherm38
         ror(jher) = 0.0d+00
      end do
c
      nh1 = nherm3
      nh2 = nh1 + nherm3
      nh3 = nh2 + nherm3
      nh4 = nh3 + nherm3
      nh5 = nh4 + nherm3
      nh6 = nh5 + nherm3
      nh7 = nh6 + nherm3
c
c---->        K1 >= K2
c---->    e sufix stands for even parity
c---->    o sufix stands for odd  parity
c
      irop = 0
      irom = 0
      xneckp = 0.0d+00
      xneckn = 0.0d+00
      do 1 nx1=1,nxmax
        nym1 = MY(nx1)
        do 11 nx2=1,nx1
          phasx= Dfloat(1-2*Mod((nx2+1),2))
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          lcx  = Mod((nx1+nx2),2).eq.0
c
              do jher=1,nherm26
                 roryz(jher) = 0.0d+00
              end do
c
          do 2 ny1 = 1,nym1
            nxy1 = nx1+ny1
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do 12 ny2=1,nym2
              lcy  = Mod((ny1+ny2),2) .eq.0
              phas = Dfloat(1-2*Mod(Iabs(ny1-ny2)/2,2))
              phase= Dfloat(1-2*Mod(Iabs(ny1-ny2+1)/2,2))
              nxy2 = nx2+ny2
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              lxy12 = lx12.and.(ny1.eq.ny2)
c
              do jher=1,nherm4
                 rorz(jher) = 0.0d+00
              end do
c                                                  ----> Positive parity
c
c be careful not to move the if below, the iz1e=iy1e need to be done
c
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
c
              do 30 nz1 =nzi1e,nzm1,2
                if(lxy12) nzm2=nz1
                do 31 nz2 = nzi2e,nzm2,2
                feq = 2.0d+00
                if(lxy12.and.(nz1.eq.nz2)) feq = 1.0d+00
c
c    GROP: rosp, roap, rospw, roapw, akasp, akaap, akaspw, akaapw
c
c        sum(1) =  ro1n+ro2n
c        sum(2) =  ro1n-ro2n
c        sum(3) =  ro1p+ro2p
c        sum(4) =  ro1p-ro2p
c        sum(5) =  ak1n+ak2n
c        sum(6) =  ak1n-ak2n
c        sum(7) =  ak1p+ak2p
c        sum(8) =  ak1p-ak2p
c
         call dcopy(8,grop(irop+1),1,sum,1)
         irop = irop + 8
                if((.not.lcy).and.lcx) goto 31     ! no field 2 density
                ro0n = sum(1)
                ro0p = sum(3)
                ro2n = sum(2)
                ro2p = sum(4)
                iher = 1
                do nzher=1,nherm
                   wfz = wfho(nz1,nzher)*wfho(nz2,nzher)
                   rorz(iher)  =rorz(iher  )+ro0p*feq*wfz
                   rorz(iher+1)=rorz(iher+1)+ro0n*feq*wfz
                   rorz(iher+2)=rorz(iher+2)+ro2p*feq*wfz
                   rorz(iher+3)=rorz(iher+3)+ro2n*feq*wfz
                   iher = iher + 4
                end do
c
                if(lxy12) then
                   do nzher=1,nherm
                     wfz = wfhonec(nz1,nzher)*wfhonec(nz2,nzher)
                     xneckp = xneckp + ro0p*feq*wfz
                     xneckn = xneckn + ro0n*feq*wfz
                   end do
                end if
c
   31           continue
   30         continue
c                                                   ----> Negative parity
              end if
c
c be careful not to move the if below, the iz1o=iy1o need to be done
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
c
              do 40 nz1 =nzi1o,nzm1,2
                if(lxy12) nzm2=nz1
                do 41 nz2 = nzi2o,nzm2,2
                feq = 2.0d+00
                if(lxy12.and.(nz1.eq.nz2)) feq = 1.0d+00
c
c    GROM: rosm, roam, rosmw, roamw, akasm, akaam, akasmw, akaamw
c
         call dcopy(8,grom(irom+1),1,sum,1)
         irom = irom + 8
                if((.not.lcy).and.lcx) goto 41   ! no field 2 density
                ro0n = sum(1)
                ro0p = sum(3)
                ro2n = sum(2)
                ro2p = sum(4)
                iher = 1
                do nzher=1,nherm
                   wfz = wfho(nz1,nzher)*wfho(nz2,nzher)
                   rorz(iher)  =rorz(iher  )+ro0p*feq*wfz
                   rorz(iher+1)=rorz(iher+1)+ro0n*feq*wfz
                   rorz(iher+2)=rorz(iher+2)+ro2p*feq*wfz
                   rorz(iher+3)=rorz(iher+3)+ro2n*feq*wfz
                   iher = iher + 4
                end do
c
                if(lxy12) then
                   do nzher=1,nherm
                     wfz = wfhonec(nz1,nzher)*wfhonec(nz2,nzher)
                     xneckp = xneckp + ro0p*feq*wfz
                     xneckn = xneckn + ro0n*feq*wfz
                   end do
                end if
c
c
   41           continue
   40         continue
            end if
c
            iher = 1
            do nyher=1,nherm
               wfy = wfho(ny1,nyher)*wfho(ny2,nyher)
               iherz = 1
               do nzher=1,nherm
                  if(lcy) then
                  roryz(iher  )=roryz(iher  )+rorz(iherz  )*wfy*phas
                  roryz(iher+1)=roryz(iher+1)+rorz(iherz+1)*wfy*phas
                  roryz(iher+2)=roryz(iher+2)+rorz(iherz+2)*wfy*phas
                  roryz(iher+3)=roryz(iher+3)+rorz(iherz+3)*wfy*phas
                  else
                  roryz(iher+4)=roryz(iher+4)+rorz(iherz+2)*wfy*phase
                  roryz(iher+5)=roryz(iher+5)+rorz(iherz+3)*wfy*phase
                  end if
               iherz = iherz + 4
               iher  = iher  + 6
               end do
            end do
c
   12       continue
    2     continue
          iher = 0
          do nxher=1,nherm
             wfx = wfho(nx1,nxher)*wfho(nx2,nxher)
             ihery = 1
             do nyzher=1,nherm2
                iher = iher + 1
                if(lcx) then
                ror(iher    )=ror(iher    )+roryz(ihery  )*wfx
                ror(iher+nh1)=ror(iher+nh1)+roryz(ihery+1)*wfx
                ror(iher+nh4)=ror(iher+nh4)+roryz(ihery+2)*wfx*phasx
                ror(iher+nh5)=ror(iher+nh5)+roryz(ihery+3)*wfx*phasx
                else
                ror(iher+nh2)=ror(iher+nh2)+roryz(ihery+2)*wfx
                ror(iher+nh3)=ror(iher+nh3)+roryz(ihery+3)*wfx
                ror(iher+nh6)=ror(iher+nh6)+roryz(ihery+4)*wfx*phasx
                ror(iher+nh7)=ror(iher+nh7)+roryz(ihery+5)*wfx*phasx
                end if
                ihery = ihery + 6
             end do
          end do
   11   continue
    1 continue
c
      do nher=1,nherm38
         ror(nher) = ror(nher)*bm1
      end do
c
      rneckp = xneckp*fneck
      rneckn = xneckn*fneck
      return
      end
