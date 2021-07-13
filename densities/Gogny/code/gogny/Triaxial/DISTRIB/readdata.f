      Subroutine ReadData
      Implicit real*8 (A-H,O-Z)
      Logical vcom,vcom2
      Character*2 cnuc
      Character*8 TEE
      Parameter (NCM=06,NLAM=NCM-3)  ! maximum number of constraints
      Parameter (NEEE=30) ! dimension of the output matrix
      Common /CCONS/ cval(ncm),prec(ncm),IC(NCM),
     & ILAMB(NLAM),IMM(NLAM),ISOS(NCM,2)
      Common /CITER / epsg,etamax,etamin,dmax,dmin,tshh,maxiter
      Common /HFBOPT/ Acom,vcom,vcom2,icouech
      Common/NECKCO/ aneck,rneckp,rneckn
      Common /STARTO/ iwf,iol
      Common /CDEB  / ideb(9)
      Common /OUTLEV/ ioutlev
      Common /OSCLEN/ bx,by,bz
      Common /CNECK / z00,a00
      Common /EEEEEE/ TEE(NEEE),EEE(7,NEEE)
      Common /IIIEEE/ mkin,mehf,mepair,mehfb,medd,merea,mecou,
     &                mn,mjx,msx,mq20,mq22,mq40,mneck,mr2,
     &                mbet2,mgamm,mbet4,maxrat,
     &                mfn,mjx2,mjy2,mjz2
      aneck = 2.0d+00
c
      do i=1,9
         ideb(i) = 0
      end do
c
      ideb(8) = 1 ! 1 = no spin-orbit pairing
      
      IC(1) = 1
      IC(2) = 1
c
      do i=1,ncm
         ISOS(i,1) = 1
         ISOS(i,2) = 1
      end do

      ISOS(1,2) = 0
      ISOS(2,1) = 0
c      
c      ISOS(7,2) = 0
c      ISOS(8,1) = 0
c
      read(5,1) ia,cnuc,iz,icom,icom2,icouech
    1 format (11x,I3,1x,a2,6x,i3,25x,i1,5x,i1,9x,i1)
c
      write(6,101)ia,cnuc,iz,icom,icom2,icouech
  101 format (2x,10('========'),//,30x,'T R I A X I A L     H F B ',//,
     &        32x,'NUCLEUS ',i4,2x,a2,' Z= ',i4,//,30x,
     &        'COM ',i1,' COM2 ',i1,' COU-ECH ',i1,//,
     &        2x,10('========'),//)
      in = ia-iz
      cval(1) = dfloat(iz)
      cval(2) = dfloat(in)
      Acom    = dfloat(ia)
      vcom    = icom .ne.0
      vcom2   = icom2.ne.0
c
      read(5,2) epsg,maxiter,ioutlev,itstg
    2 format(4x,f9.6,9x,i5,20x,i1,11x,i1)
      ideb(9) = itstg
      write(6,102) epsg,maxiter
  102 format (30x,'EPSG ',f9.6,' MAXITER ',i5)
c
      read(5,3) etamax,etamin,dmax,dmin,tshh
    3 format( 5x,f7.4,6x,f7.4,5x,f5.2,5x,f5.2,5x,f6.1)
c
      read(5,4) ipar
    4 format(19x,i1)
c+--------------------------------------------------------------------+
c|+------------------------------------------------------------------+|
c||   ipar = 0     D1S                                               ||
c||   ipar = 1     D1             G O G N Y    F O R C E             ||
c||   ipar = 2     D1'                                               ||
c|+------------------------------------------------------------------+|
c+--------------------------------------------------------------------+
      write(6,112)
c      
      call Gognyf(ipar)
c
      write(6,112)
c      
      read(5,4) iwf
      read(5,5) iol,bx,by,bz
    5 format(19x,i1,26x,f7.4,3x,f7.4,3x,f7.4)
      read(5,6) ic(3),cval(3),prec(3)
    6 format( /,8x,i1,6x,d16.8,5x,d10.2)
      if(ic(3).ne.0) then 
             cval(3) = dsqrt(cval(3)*(cval(3)+1.0d+00))
             write(6,106) cval(3),prec(3)
      end if
  106 format( //,15x,' CONSTRAINT ON <Jx>   ',d16.8,2x,d10.2)
      read(5,7) ILAMB(1),IMM(1),ic(4),cval(4),prec(4)
      if(ic(4).ne.0) write(6,107) ILAMB(1),IMM(1),cval(4),prec(4)
  107 format( /,15x,' CONSTRAINT ON Q',2i1,'    ',d16.8,2x,d10.2)
    7 format( 1x,2i1,5x,i1,6x,d16.8,5x,d10.2)
      read(5,7) ILAMB(2),IMM(2),ic(5),cval(5),prec(5)       
      if(ic(5).ne.0) write(6,108) ILAMB(2),IMM(2),cval(5),prec(5)
  108 format( /,15x,' CONSTRAINT ON Q',2i1,'    ',d16.8,2x,d10.2)
      read(5,7) ILAMB(3),IMM(3),ic(6),cval(6),prec(6)       
      if(ic(6).ne.0) write(6,109) ILAMB(3),IMM(3),cval(6),prec(6)
  109 format( /,15x,' CONSTRAINT ON Q',2i1,'    ',d16.8,2x,d10.2)
c
      write(6,112)
  112 format(//,2x,10('========'),//)
      mkin    = 1
      TEE(mkin  ) = 'Kinetic '
      mehf    = 2
      TEE(mehf  ) = 'HF Ener '
      mepair  = 3
      TEE(mepair) = 'Pairing '
      mehfb   = 4
      TEE(mehfb ) = 'HFB Ener'
      medd    = 5
      TEE(medd  ) = 'DD Ener '
      merea   = 6
      TEE(merea ) = 'Rearrang'
      mecou   = 7
      TEE(mecou ) = 'Coul Ex.'
      mn      = 8
      TEE(mn    ) = '   N    '
      mjx     = 9
      TEE(mjx   ) = '  Jx    '
      msx     = 10
      TEE(msx   ) = '  Sx    '
      mq20    = 11
      TEE(mq20  ) = ' Q20(b) '
      mq22    = 12
      TEE(mq22  ) = ' Q22(b) '
      mq40    = 13
      TEE(mq40  ) = ' Q40(b2)'
      mneck   = 14
      TEE(mneck ) = '  Neck  '
      mr2     = 15
      TEE(mr2   ) = 'MS Rad  '
      mbet2   = 16
      TEE(mbet2 ) = ' Beta 2 '
      mgamm   = 17
      TEE(mgamm ) = ' Gamma  '
      mbet4   = 18
      TEE(mbet4 ) = ' Beta 4 '
      maxrat  = 19
      TEE(maxrat) = '  z/x   '
      mfn     = 20
      TEE(mfn   ) = ' <N^2>  '
      mjx2    = 21
      TEE(mjx2  ) = '<Jx**2> '
      mjy2    = 22
      TEE(mjy2  ) = '<Jy**2> '
      mjz2    = 23
      TEE(mjz2  ) = '<Jz**2> '
      
      TEE(24    ) = ' KIN2 G '
      TEE(25    ) = ' KIN2 D '
      TEE(26    ) = '        '
      TEE(27    ) = '        '
      TEE(28    ) = '        '
      TEE(29    ) = '        '
      TEE(30    ) = '        '
      return
      end
