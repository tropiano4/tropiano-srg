      Subroutine PACK2SQ(OPEPACK,OPOPACK,OPE,OPO,is)
      Implicit real*8(a-h,o-z)
      Logical LMZE,LMZO
      Logical lx12,lxy12
      Include 'COMDIM'
      Dimension OPE(NP,NP),OPO(NM,NM)
      Dimension OPEPACK(NGP),OPOPACK(NGM)
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      Common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
      Common/ILOC/NP,NM,ILOCE(ixmax,iymax),ILOCO(ixmax,iymax)

c
c---->    q'=1  q=2    q' >= q
c
c---->    e sufix stands for even parity
c---->    o sufix stands for odd  parity
      irop = 0
      irom = 0
      do nx1=1,nxmax
        nym1 = MY(nx1)
        do nx2=1,nx1
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          do ny1 = 1,nym1
            iz1e0 = ILOCE(nx1,ny1)
            iz1o0 = ILOCO(nx1,ny1)
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do ny2=1,nym2
              iz2e0 = ILOCE(nx2,ny2)
              iz2o0 = ILOCO(nx2,ny2)
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              lxy12 = lx12.and.(ny1.eq.ny2)
c
c---->    Positive parity
c
c
c be careful not to move the if below, the iz1e=iy1e need to be done
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
c
              iz1e = iz1e0
              do nz1 =nzi1e,nzm1,2
                iz1e = iz1e + 1
                if(lxy12) nzm2=nz1
                iz2e = iz2e0
                do nz2 = nzi2e,nzm2,2
                   iz2e = iz2e + 1
                   irop = irop + 1
                   OPE(iz1e,iz2e)= OPEPACK(irop)
                   if(is.eq.1) then
                      ope(iz2e,iz1e) = ope(iz1e,iz2e)
                   else if (is.eq.-1) then
                     ope(iz2e,iz1e) = - ope(iz1e,iz2e)
                   end if
                end do
              end do
c
c---->    Negative parity
c
              end if
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
c
              iz1o = iz1o0
              do nz1 =nzi1o,nzm1,2
                iz1o = iz1o+ 1
                if(lxy12) nzm2=nz1
                iz2o = iz2o0
                do nz2 = nzi2o,nzm2,2
                   iz2o = iz2o+1
                   irom = irom + 1
                   OPO(iz1o,iz2o) = OPOPACK(irom)
                   if(is.eq.1) then
                     opo(iz2o,iz1o) = opo(iz1o,iz2o)
                   else if (is.eq.-1) then
                     opo(iz2o,iz1o) = - opo(iz1o,iz2o)
                   end if
                end do
              end do
            end if
            end do    ! ny2
          end do      ! ny1
        end do        ! nx2
      end do          ! nx1
      return
      end
      Subroutine SQ2PACK(OPE,OPO,OPEPACK,OPOPACK)
      Implicit real*8(a-h,o-z)
      Logical LMZE,LMZO
      Logical lx12,lxy12
      Include 'COMDIM'
      Dimension OPE(NP,NP),OPO(NM,NM)
      Dimension OPEPACK(NGP),OPOPACK(NGM)
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      Common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
      Common/ILOC/NP,NM,ILOCE(ixmax,iymax),ILOCO(ixmax,iymax)
c
c---->    q'=1  q=2    q' >= q
c
c---->    e sufix stands for even parity
c---->    o sufix stands for odd  parity
      irop = 0
      irom = 0
      do nx1=1,nxmax
        nym1 = MY(nx1)
        do nx2=1,nx1
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          do ny1 = 1,nym1
            iz1e0 = ILOCE(nx1,ny1)
            iz1o0 = ILOCO(nx1,ny1)
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do ny2=1,nym2
              iz2e0 = ILOCE(nx2,ny2)
              iz2o0 = ILOCO(nx2,ny2)
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              lxy12 = lx12.and.(ny1.eq.ny2)
c
c---->    Positive parity
c
c
c be careful not to move the if below, the iz1e=iy1e need to be done
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
c
              iz1e = iz1e0
              do nz1 =nzi1e,nzm1,2
                iz1e = iz1e + 1
                if(lxy12) nzm2=nz1
                iz2e = iz2e0
                do nz2 = nzi2e,nzm2,2
                   iz2e = iz2e + 1
                   irop = irop + 1
                   OPEPACK(irop) = OPE(iz1e,iz2e)
                end do
              end do
c
c---->    Negative parity
c
              end if
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
c
              iz1o = iz1o0
              do nz1 =nzi1o,nzm1,2
                iz1o = iz1o+ 1
                if(lxy12) nzm2=nz1
                iz2o = iz2o0
                do nz2 = nzi2o,nzm2,2
                   iz2o = iz2o+1
                   irom = irom + 1
                   OPOPACK(irom) = OPO(iz1o,iz2o)
                end do
              end do
            end if
            end do    ! ny2
          end do      ! ny1
        end do        ! nx2
      end do          ! nx1
      return
      end
      subroutine d2sq(d1p,d2p,d1m,d2m,aux1p,aux2p,aux1m,aux2m)
      Implicit real*8(a-h,o-z)
      Logical LMZE,LMZO
      Logical lx12,lxy12
      Include 'COMDIM'
      Dimension aux1p(NP,NP),aux1m(NM,NM)
      Dimension aux2p(NP,NP),aux2m(NM,NM)
      Dimension d1p(ngp),d2p(ngp),d1m(ngm),d2m(ngm)
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
c
      Common/ILOC/NP,NM,ILOCE(ixmax,iymax),ILOCO(ixmax,iymax)
c
c---->    q'=1  q=2    q' >= q
c---->    e sufix stands for even parity
c---->    o sufix stands for odd  parity
      irop = 0
      irom = 0
      do nx1=1,nxmax
        nym1 = MY(nx1)
        do nx2=1,nx1
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          do ny1 = 1,nym1
            iz1e0 = ILOCE(nx1,ny1)
            iz1o0 = ILOCO(nx1,ny1)
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do ny2=1,nym2
              iz2e0 = ILOCE(nx2,ny2)
              iz2o0 = ILOCO(nx2,ny2)
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              lxy12 = lx12.and.(ny1.eq.ny2)
c
c---->    Positive parity
c
c
c be careful not to move the if below, the iz1e=iy1e need to be done
c
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
c
              iz1e = iz1e0
              do nz1 =nzi1e,nzm1,2
                iz1e = iz1e + 1
                if(lxy12) nzm2=nz1
                iz2e = iz2e0
                do nz2 = nzi2e,nzm2,2
                iz2e = iz2e + 1
                irop = irop + 1
c		write(6,*) ' d2sq p ',irop,d1p(irop),d2p(irop)
                AUX1P(iz1e,iz2e) = d1p(irop)
                AUX2P(iz2e,iz1e) =-d1p(irop)
                AUX1P(iz2e,iz1e) = d2p(irop)
                AUX2P(iz1e,iz2e) =-d2p(irop)
                end do
              end do
c
c---->    Negative parity
c
              end if
c
c be careful not to move the if below, the iz1o=iy1o need to be done
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
c
              iz1o = iz1o0
              do nz1 =nzi1o,nzm1,2
                iz1o = iz1o+ 1
                if(lxy12) nzm2=nz1
                iz2o = iz2o0
                do nz2 = nzi2o,nzm2,2
                iz2o = iz2o+1
                irom = irom + 1
c		write(6,*) ' d2sq m ',irom,d1m(irom),d2m(irom)
                AUX1M(iz1o,iz2o) = d1m(irom)
                AUX2M(iz2o,iz1o) =-d1m(irom)
                AUX1M(iz2o,iz1o) = d2m(irom)
                AUX2M(iz1o,iz2o) =-d2m(irom)
                end do
              end do
            end if
            end do     ! ny2
          end do       ! ny1
        end do         ! nx2
      end do           ! nx1
c
      if((irop.gt.ngp).or.(irom.gt.ngm))
     *   call errout('Troubles        ','D2SQ    ',4)
c
      return
      end
      subroutine de2sq(d1p,d2p,d1m,d2m,aux1p,aux1m)
      Implicit real*8(a-h,o-z)
      Logical LMZE,LMZO
      Logical lx12,lxy12
      Include 'COMDIM'
      Dimension aux1p(NP,NP),aux1m(NM,NM)
      Dimension d1p(ngp),d2p(ngp),d1m(ngm),d2m(ngm)
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
c
      Common/ILOC/NP,NM,ILOCE(ixmax,iymax),ILOCO(ixmax,iymax)
c
c---->    q'=1  q=2    q' >= q
c---->    e sufix stands for even parity
c---->    o sufix stands for odd  parity
      irop = 0
      irom = 0
      do nx1=1,nxmax
        nym1 = MY(nx1)
        do nx2=1,nx1
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          do ny1 = 1,nym1
            iz1e0 = ILOCE(nx1,ny1)
            iz1o0 = ILOCO(nx1,ny1)
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do ny2=1,nym2
              iz2e0 = ILOCE(nx2,ny2)
              iz2o0 = ILOCO(nx2,ny2)
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              lxy12 = lx12.and.(ny1.eq.ny2)
c
c---->    Positive parity
c
c
c be careful not to move the if below, the iz1e=iy1e need to be done
c
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
c
              iz1e = iz1e0
              do nz1 =nzi1e,nzm1,2
                iz1e = iz1e + 1
                if(lxy12) nzm2=nz1
                iz2e = iz2e0
                do nz2 = nzi2e,nzm2,2
                iz2e = iz2e + 1
                irop = irop + 1
                AUX1P(iz1e,iz2e) = d1p(irop)
                AUX1P(iz2e,iz1e) = d2p(irop)
                end do
              end do
c
c---->    Negative parity
c
              end if
c
c be careful not to move the if below, the iz1o=iy1o need to be done
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
c
              iz1o = iz1o0
              do nz1 =nzi1o,nzm1,2
                iz1o = iz1o+ 1
                if(lxy12) nzm2=nz1
                iz2o = iz2o0
                do nz2 = nzi2o,nzm2,2
                iz2o = iz2o+1
                irom = irom + 1
c		write(6,*) ' d2sq m ',irom,d1m(irom),d2m(irom)
                AUX1M(iz1o,iz2o) = d1m(irom)
                AUX1M(iz2o,iz1o) = d2m(irom)
                end do
              end do
            end if
            end do     ! ny2
          end do       ! ny1
        end do         ! nx2
      end do           ! nx1
c
      if((irop.gt.ngp).or.(irom.gt.ngm))
     *   call errout('Troubles        ','DE2SQ    ',4)
c
      return
      end
      Subroutine JT2disk(ajxe,ajxo,ajye,ajyo,ajze,ajzo,te,to,GROP,GROM)
      Implicit real*8(a-h,o-z)
      Include 'COMDIM'
      Logical LMZE,LMZO
      Logical lx12,lxy12
      Dimension AJXE(NP,NP),AJXO(NM,NM)
      Dimension ajyE(NP,NP),ajyO(NM,NM)
      Dimension ajzE(NP,NP),ajzO(NM,NM)
      Dimension   TE(NP,NP),  TO(NM,NM)
      Dimension GROP(ngp8),GROM(ngm8)
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      Common/ILOC/NP,NM,ILOCE(ixmax,iymax),ILOCO(ixmax,iymax)
c
      Common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
c
      ngp2 = 2*ngp
      ngp3 = 3*ngp
      ngp4 = 4*ngp
      ngm2 = 2*ngm
      ngm3 = 3*ngm
      ngm4 = 4*ngm
c
c---->    q'=1  q=2    q' >= q                <--------
c---->                                        <--------
c---->    e sufix stands for even parity      <--------
c---->    o sufix stands for odd  parity      <--------
      irop = 0
      irom = 0
      do nx1=1,nxmax
        nym1 = MY(nx1)
        do nx2=1,nx1
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          do ny1 = 1,nym1
            iz1e0 = ILOCE(nx1,ny1)
            iz1o0 = ILOCO(nx1,ny1)
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do ny2=1,nym2
              iz2e0 = ILOCE(nx2,ny2)
              iz2o0 = ILOCO(nx2,ny2)
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              lxy12 = lx12.and.(ny1.eq.ny2)
c
c---->    Positive parity
c
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
c
                  iz1e = iz1e0
                  do nz1 =nzi1e,nzm1,2
                    iz1e = iz1e + 1
                    if(lxy12) nzm2=nz1
                    iz2e = iz2e0
                    do nz2 = nzi2e,nzm2,2
                        iz2e = iz2e + 1
                        irop = irop + 1
                        GROP(irop     ) =   TE(iz1e,iz2e)
                        GROP(irop+ngp ) = AJXE(iz1e,iz2e)
                        GROP(irop+ngp2) = AJZE(iz1e,iz2e)
                        GROP(irop+ngp3) = AJYE(iz1e,iz2e)
                    end do    ! nz1
                  end do      ! nz2
              end if      ! positive parity condition
c
c---->    Negative parity
c
c
c be careful not to move the if below, the iz1o=iy1o need to be done
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
c
                  iz1o = iz1o0
                  do nz1 =nzi1o,nzm1,2
                    iz1o = iz1o+ 1
                    if(lxy12) nzm2=nz1
                    iz2o = iz2o0
                    do nz2 = nzi2o,nzm2,2
                       iz2o = iz2o+1
                       irom = irom + 1
                       GROM(irom     ) =   TO(iz1o,iz2o)
                       GROM(irom+ngm ) = AJXO(iz1o,iz2o)
                       GROM(irom+ngm2) = AJZO(iz1o,iz2o)
                       GROM(irom+ngm3) = AJYO(iz1o,iz2o)
                    end do  ! nz2
                  end do    ! nz1
              end if    ! negative parity condition
            end do      ! ny2
          end do        ! ny1
        end do          ! nx2
      end do            ! nx1
c
      REWIND 98
      REWIND 99
c
      WRITE(98) (grop(j),j=1,ngp4)
      WRITE(99) (grom(j),j=1,ngm4)
c
      REWIND 98
      REWIND 99
c
      return
      end
      Subroutine RO2PACK(ROM,AKA,grop,grom)
      Implicit real*8(a-h,o-z)
      Logical LMZE,LMZO
      Logical lx12,lxy12
      Include 'COMDIM'
      Dimension sum(8)
      Dimension ROM(*),AKA(*)
      Dimension grop(NGP8),grom(NGM8)
c
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      Common /DIMECH/maxtjx,maxtjy,maxtjz,mazopti,mayopti,ngp,ngp8,
     *               ngm,ngm8
      Common/ILOC/NP,NM,ILOCE(ixmax,iymax),ILOCO(ixmax,iymax)
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC1/ NNU(4),NNV(4),NNRO(4),NNKA(4)
c
c---->    q'=1  q=2    q' >= q
c
c---->    e sufix stands for even parity
c---->    o sufix stands for odd  parity
      irop = 0
      irom = 0
      do nx1=1,nxmax
        nym1 = MY(nx1)
        do nx2=1,nx1
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          do ny1 = 1,nym1
            iz1e0 = ILOCE(nx1,ny1)
            iz1o0 = ILOCO(nx1,ny1)
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do ny2=1,nym2
              iz2e0 = ILOCE(nx2,ny2)
              iz2o0 = ILOCO(nx2,ny2)
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              lxy12 = lx12.and.(ny1.eq.ny2)
c
c---->    Positive parity
c
c
c be careful not to move the if below, the iz1e=iy1e need to be done
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
c
              iz1e = iz1e0
              do nz1 =nzi1e,nzm1,2
                iz1e = iz1e + 1
                if(lxy12) nzm2=nz1
                iz2e = iz2e0
                do nz2 = nzi2e,nzm2,2
                   iz2e = iz2e + 1
		   
		   ro1n =  ROM(iz2e-1+(iz1e-1)*ND(3)+NNRO(3))       ! ro1 neut
		   ro2n =  ROM(iz2e-1+(iz1e-1)*ND(3)+NNRO(3)+ND2(3))! ro2 neut
		   ro1p =  ROM(iz2e-1+(iz1e-1)*ND(1)+NNRO(1))       ! ro1 prot
		   ro2p =  ROM(iz2e-1+(iz1e-1)*ND(1)+NNRO(1)+ND2(1))! ro2 prot
		   
		   ak1n =  AKA(iz2e-1+(iz1e-1)*ND(3)+NNKA(3))	  ! k 1 neut
		   ak2n =  AKA(iz1e-1+(iz2e-1)*ND(3)+NNKA(3))	  !-k 2 neut
		   ak1p =  AKA(iz2e-1+(iz1e-1)*ND(1)+NNKA(1))	  ! k 1 prot
		   ak2p =  AKA(iz1e-1+(iz2e-1)*ND(1)+NNKA(1))	  !-k 2 prot
		   
                   sum(1) =  ro1n+ro2n
                   sum(2) =  ro1n-ro2n
                   sum(3) =  ro1p+ro2p
                   sum(4) =  ro1p-ro2p
                   sum(5) =  ak1n+ak2n
                   sum(6) =  ak1n-ak2n
                   sum(7) =  ak1p+ak2p
                   sum(8) =  ak1p-ak2p
		   
                do ir=1,8
                   irop = irop + 1
                   GROP(irop) = sum(ir)
                end do
		
                end do
              end do
c
c---->    Negative parity
c
              end if
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
c
              iz1o = iz1o0
              do nz1 =nzi1o,nzm1,2
                iz1o = iz1o+ 1
                if(lxy12) nzm2=nz1
                iz2o = iz2o0
                do nz2 = nzi2o,nzm2,2
                   iz2o = iz2o+1
		   
		   ro1n =  ROM(iz2o-1+(iz1o-1)*ND(4)+NNRO(4))       ! ro1 neut
		   ro2n =  ROM(iz2o-1+(iz1o-1)*ND(4)+NNRO(4)+ND2(4))! ro2 neut
		   ro1p =  ROM(iz2o-1+(iz1o-1)*ND(2)+NNRO(2))       ! ro1 prot
		   ro2p =  ROM(iz2o-1+(iz1o-1)*ND(2)+NNRO(2)+ND2(2))! ro2 prot
		   
		   ak1n =  AKA(iz2o-1+(iz1o-1)*ND(4)+NNKA(4))	  ! k 1 neut
		   ak2n =  AKA(iz1o-1+(iz2o-1)*ND(4)+NNKA(4))	  !-k 2 neut
		   ak1p =  AKA(iz2o-1+(iz1o-1)*ND(2)+NNKA(2))	  ! k 1 prot
		   ak2p =  AKA(iz1o-1+(iz2o-1)*ND(2)+NNKA(2))	  !-k 2 prot
		   
                   sum(1) =  ro1n+ro2n
                   sum(2) =  ro1n-ro2n
                   sum(3) =  ro1p+ro2p
                   sum(4) =  ro1p-ro2p
                   sum(5) =  ak1n+ak2n
                   sum(6) =  ak1n-ak2n
                   sum(7) =  ak1p+ak2p
                   sum(8) =  ak1p-ak2p
		   
                do ir=1,8
                   irom = irom + 1
                   GROM(irom) = sum(ir)
                end do
		
                end do
              end do
            end if
            end do    ! ny2
          end do      ! ny1
        end do        ! nx2
      end do          ! nx1
      return
      end
      Subroutine field2sq(gfield,gam,delta)
      Implicit real*8 (a-h,o-z)
      Dimension gfield(*)
      Dimension gam(*),delta(*)
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC1/ NNU(4),NNV(4),NNRO(4),NNKA(4)
      Common /ITLOC2/ ND2PACK(4),NNGPACK(4),NNDPACK(4)
      
c
c  Transforming gamma to square form
c

       do it=1,3,2
          
	  IGP = NNGPACK(it)
	  IGM = NNGPACK(it+1)
	  IGGP= NNRO(it)
	  IGGM= NNRO(it+1)
c
          call pack2sq(gfield(IGP),gfield(IGM),Gam(IGGP),Gam(IGGM),1)

       end do
       
       do it=1,3,2
          
	  IGP = NNGPACK(it)   + ND2PACK(it)
	  IGM = NNGPACK(it+1) + ND2PACK(it+1)
	  IGGP= NNRO(it)      + ND2(it)
	  IGGM= NNRO(it+1)    + ND2(it+1)
c
          call pack2sq(gfield(IGP),gfield(IGM),Gam(IGGP),Gam(IGGM),1)

       end do
c
c  Transforming delta to square form
c
       do it=1,3,2
          
	  IG1P = NNDPACK(it)
	  IG1M = NNDPACK(it+1)
	  IGG1P= NNKA(it)
	  IGG1M= NNKA(it+1)
	  
	  IG2P = NNDPACK(it)     + ND2PACK(it)
	  IG2M = NNDPACK(it+1)	 + ND2PACK(it+1)
c
          call de2sq(gfield(IG1P),gfield(IG2P),gfield(IG1M),gfield(IG2M)
     *	    ,Delta(IGG1P),Delta(IGG1M))
     
       end do
       return
       end
