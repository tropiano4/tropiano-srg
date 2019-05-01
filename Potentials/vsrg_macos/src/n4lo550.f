         MODULE newEM_n4lo550
         implicit real*8(a-h,o-z)
         PRIVATE
         PUBLIC :: n4lo550
         CONTAINS
         subroutine n4lo550
c
c******************************************************************
c
c        Final version, December 2016.
c
c******************************************************************
c
c        This code computes the
c
c        Charge-Dependent Chiral NN Potential at Order Five
c        --------------------------------------------------
c        applying a cutoff of 550 MeV: N4LO(550)
c        ---------------------------------------
c
c        this package is self-contained and includes
c        all subroutines needed.
c        only `n4lo550' needs to be called by the user.
c        all codes are consistently in double precision.
c        when working on an UNIX/LINUX system, it is recommended
c        to compile this code with the  "static"  option 
c        ( " -s "  on some compilers).
c        more information on the code is given below.
c
c*******************************************************************
c
c        authors:     D. R. Entem, R. Machleidt, and Y. Nosyk
c                     department of physics
c                     university of idaho
c                     moscow, idaho 83844-0903
c                     u. s. a.
c                     e-mail:  machleid@uidaho.edu
c
c        References:
c        D. R. Entem, N. Kaiser, R. Machleidt, and Y. Nosyk,
c        Phys. Rev. C 91, 014002 (2015).
c                     and
c        D. R. Entem, R. Machleidt, and Y. Nosyk,
c        arXiv:1703:05454 [nucl-th].
c       
c*******************************************************************
c*******************************************************************
c
c
c
c
      implicit real*8 (a-h,o-z)
c
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
c        arguments and values of this subroutine
c
      common /cpot/   v(6),xmev,ymev
      common /cstate/ j,heform,sing,trip,coup,endep,label
      common /cnn/ inn
c
c
c        this has been the end of the common-blocks containing
c        the arguments and values of this subroutine.
c
c        specifications for these common blocks
c
      logical heform,sing,trip,coup,endep
      character*4 label
c
c
c*********************************************************************
c        THE ABOVE FOUR COMMON BLOCKS IS ALL THE USER NEEDS
c        TO BE FAMILIAR WITH.
c*********************************************************************
c
c        here are now some explanations of what those common blocks contain:
c        -------------------------------------------------------------------
c
c        xmev and ymev are the final and initial relative momenta,
c        respectively, in units of mev/c.
c        v is the potential in units of mev**(-2).
c        concerning units, factors of pi, etc.,
c        cf. with the partial-wave Lippmann-Schwinger equation, Eq. (A25),
c        and with the phase shift relation, Eq. (A33), given in Appendix A 
c        of the article: R. Machleidt, PRC 63, 024001 (2001).
c
c        the partial-wave Lippmann-Schwinger equation for the
c        K-matrix reads:
c
c        K(q',q) = V(q',q) + M P \int dk k^2 V(q',k) K(k,q)/(q^2-k^2)
c
c        with M the nucleon mass in MeV and P denoting the principal value;
c        V(q',q) as provided by this code in common block /cpot/;
c        all momenta in MeV.
c
c        the phase-shift relation is:
c
c        tan \delta_L = -(pi/2) M q K_L(q,q)
c
c        with M and q in units of MeV, K_L in MeV**(-2) like V.
c
c
c        if heform=.true., v contains the 6 matrix elements
c        associated with one j in the helicity formalism
c        in the following order:
c        0v, 1v, 12v, 34v, 55v, 66v 
c        (for notation see Appendix A of above article).
c
c        if heform=.false., v contains the 6 matrix elements
c        associated with one j in the lsj formalism
c        in the following order:
c        0v(singlet), 1v(uncoupled triplet), v++, v--, v+-, v-+ (coupled)
c        (for notation, see explanations given in the above article 
c        below Eq. (A31)).
c
c        j is the total angular momentum. there is essentially no upper
c        limit for j.
c        sing, trip, and coup should in general be .true..
c        endep and label can be ignored.
c        it is customary, to set kread=5 and kwrite=6;
c        ignore kpunch and kda(9).
c
c        the meaning of the parameter inn in the common block
c
c                  common /cnn/ inn 
c        is
c                  inn=1  means pp potential,
c                  inn=2  means np potential, and
c                  inn=3  means nn potential.
c
c        the user needs to include this common block in his/her code,
c        and specify which potential he/she wants to use. 
c
c
c        THIS IS ESSENTIALLY ALL THE USER NEEDS TO KNOW.
c
c        if you have further questions, do not hesitate to contact 
c        machleid@uidaho.edu
c
c**********************************************************************
c
c
c        common block for all chi-subroutines
c
      common /cchi/ vj(32,270),c(20,270),fff,ff,f(52),aa(96),ai(19,30),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                gaa(3),fpia(3),ezz1(3),ezz2(3),ct(96),wt(96),
     3                ic(20,270),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(40,3),mggo(40,3),ima(30,40,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,270),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
      common /comlsj/ clsj(15,50),cutlsj(15,50),indlsj
      logical indlsj
c
      common /crrr/ rrr
c
c
c        further specifications
c
      dimension vl(4),adminv(4,4),ldminv(4),mdminv(4)
      dimension vv0(6),vv2(6),vv4(6)
      character*4 nucnuc(3)
      character*4 mesong(40)
      logical index
      logical indmg(40)
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     1            '1-  ','1-t ','1-tt','1-st','1-ss',
     2            'c   ','ss  ','ls  ','sq  ','sk  ',
     3            'sl  ',
     4         24*'    '/
      data index/.false./
      data indmg/40*.false./
      data jj/-1/
      data pi/3.141592653589793d0/
      data innn/-1/
      data nucnuc/'537p','537n','537m'/
c******************    2016.12.09    *********
      save
c
c
c
c
      if (index) go to 10
      index=.true.
c
c
c        call subroutine chipar once and only once
c
c
      call chipar
c     -----------
c
c
c        if you want the potential to be zero for very large momenta,
c        choose rrr=1000.
c        if you want no technical problems in the calculation of the deuteron
c        wave functions, choose rrr=80.
c
      rrr=80.
c
   10 continue
c
c
c
c
      if (inn.lt.1.or.inn.gt.3) then
c        choose the np potential as the default:
      inn=2
      endif
      if (j.lt.0) then
      write (kwrite,19002)
19002 format (////' error in n4lo: total angular momentum j',
     1' is negative.'/' execution terminated.'////)
      stop
      endif
c
c
c
c
c        set the inn dependent parameters
c
      if (inn.eq.innn) go to 30
      innn=inn
      inter=inn
      label=nucnuc(inter)
c
      go to (21,22,23), inter
   21 write (kwrite,10001) 
10001 format (' The pp potential is used.')
      go to 24
   22 write (kwrite,10002) 
10002 format (' The np potential is used.')
      go to 24
   23 write (kwrite,10003) 
10003 format (' The nn potential is used.')
   24 write (kwrite,10004) 
10004 format (' -------------------------'//)
c
c
      iftgo=ift(inter)+1
      dwn=1.d0/wnn(inter)
c
c
c        prepare constant over-all factor
c
      fac=pi/(2.d0*pi)**3*dwn*dwn
c     ---------------------------
c
c
c
      iman=imaa(inter)
      imen=imea(inter)
c
      imanm1=iman-1
c
      iman1=imanm1+1
      iman2=imanm1+2
      iman3=imanm1+3
      iman4=imanm1+4
      iman5=imanm1+5
      iman6=imanm1+6
      iman7=imanm1+7
      iman8=imanm1+8
      iman9=imanm1+9
      imen24=imen-24
      imen23=imen-23
      imen22=imen-22
      imen21=imen-21
      imen15=imen-15
      imen14=imen-14
c
c
c
c
   30 if (j.eq.jj) go to 50
      jj=j
      if (j.eq.0) go to 50
      aj=dble(j)
      aj1=dble(j+1)
      a2j1=dble(2*j+1)
      aaj6=dsqrt(aj*aj1)
c
c        coefficient matrix for the translations into lsj formalism
c
      adminv(1,1)=aj1
      adminv(1,2)=aj
      adminv(1,3)=-aaj6
      adminv(1,4)=-aaj6
      adminv(2,1)=aj
      adminv(2,2)=aj1
      adminv(2,3)=aaj6
      adminv(2,4)=aaj6
      adminv(3,1)=aaj6
      adminv(3,2)=-aaj6
      adminv(3,3)=aj1
      adminv(3,4)=-aj
      adminv(4,1)=aaj6
      adminv(4,2)=-aaj6
      adminv(4,3)=-aj
      adminv(4,4)=aj1
c
c       inversion
c
      call dminv (adminv,4,deter,ldminv,mdminv)
c
c
c
c
c        prepare expressions depending on x and y
c        ----------------------------------------
c        ----------------------------------------
c
c
c
c
   50 xa=xmev*dwn
      ya=ymev*dwn
      indxy=.false.
      x=xa
      xx=x*x
      y=ya
      yy=y*y
      xy2=x*y*2.d0
      xxpyy=xx+yy
      ex=dsqrt(1.d0+xx)
      ey=dsqrt(1.d0+yy)
      eem12=(ex*ey-1.d0)*2.d0
c
c
      xy=xy2*0.5d0
      ee=ex*ey
      ree=dsqrt(ee)
      eem1=ee-1.d0
      eme=ex-ey
      emeh=eme*0.5d0
      emehq=emeh*emeh
      eep1=ee+1.d0
       epe=ex+ey
      xxyy=xx*yy
c
c
      xxpyyh=xxpyy*0.5d0
      xy3=xy*3.d0
      xy4=xy*4.d0
c
c
c
c
      do 63 iv=1,6
      vv0(iv)=0.d0
      vv2(iv)=0.d0
      vv4(iv)=0.d0
   63 v(iv)=0.d0
      do 65 il=iman,imen
      do 65 iv=1,32
   65 vj(iv,il)=0.d0
c
c
c
c
c        prepare over-all factor
c
c
      go to (70,71,72),iftgo
c
c        no additional factor
c
   70 fff=fac
      go to 80
c
c        minimal relativity
c
   71 fff=fac/ree
      go to 80
c
c        factor m/e*m/e
c
   72 fff=fac/ee
c
c
   80 continue
c
c
c
c
c        contributions
c        -------------
c        -------------
c
c
c
c
      do 5995 img=1,mge
      mg=mggo(img,inter)
      if (mg.gt.16) go to 9000
      if (mg.eq.0) go to 8000
      me=mgg(mg,inter)
      go to (9000,9000,9000,9000,9000,9000,9000,9000,9000,9000,
     1       1100,1200,1300,1400,1500,1600),mg
c
c
c
c
c        c   , central force
c        -------------------
c
c
c
c
 1100 mc=1
c
      ff=1.d0
      f(1)=2.d0
      f(2)=0.d0
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-f(1)
      f(8)=f(7)
c
      call chistr(1,1,me)
      go to 5995
c
c
c
c
c        ss  , spin-spin force
c        ---------------------
c
c
c
c
 1200 mc=1
c
      ff=1.d0
      f(1)=-6.d0
      f(2)=0.d0
      f(3)=2.d0
      f(4)=0.d0
      f(5)=0.d0
      f(6)=f(3)
      f(7)=-f(3)
      f(8)=f(7)
c
      call chistr(1,1,me)
      go to 5995
c
c
c
c
c        ls  , spin-orbit force
c        ----------------------
c
c
c
c
 1300 mc=1
c
      ff=1.d0
      f(1)=0.d0
      f(2)=0.d0
      f(3)=0.d0
      f(4)=-xy2
      f(5)=-xy2
      f(6)=0.d0
      f(7)=0.d0
      f(8)=0.d0
      f(9)=0.d0
      f(10)=+xy2
      f(11)=-xy2
c
      call chistr(2,1,me)
      go to 5995
c
c
c
c
c        sq  , sq tensor force (where q denotes the momentum transfer)
c        ---------------------
c
c
c
c
 1400 mc=1
c
      ff=1.d0
      f(1)=-xxpyy*2.0d0
      f(2)=xy*4.d0
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=(xx-yy)*2.0d0
      f(8)=-f(7)
c
      call chistr(1,1,me)
      go to 5995
c
c
c
c
c        sk  , sk tensor force (where k denotes the average momentum)
c        ---------------------
c
c
c
c
 1500 mc=1
c
      ff=0.25d0
      f(1)=-xxpyy*2.0d0
      f(2)=-xy*4.d0
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=(xx-yy)*2.0d0
      f(8)=-f(7)
c
      call chistr(1,1,me)
      go to 5995
c
c
c
c
c        sl  , "quadratic spin-orbit force"
c               or sigma-l operator
c        ----------------------------------
c
c
c
c
 1600 mc=1
c
      ff=1.d0
      f(1)=-xxyy*2.d0
      f(2)=0.d0
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=-f(1)
      f(7)=f(1)
      f(8)=f(7)
      f(9)=f(6)*2.d0
c
      call chistr(4,1,me)
      go to 5995
c
c
c
c
c
c        this has been the end of the contributions of mesons
c        ----------------------------------------------------
c
c
c
c
c        errors and warnings
c        -------------------
c
c
c
c
 9000 if (indmg(mg)) go to 5995
c**** write (kwrite,19000) mesong(mg)
19000 format(1h ////' warning in chinn: contribution ',a4,' does not exi
     1st in this program.'/' contribution ignored. execution continued.'
     2////)
      indmg(mg)=.true.
c
c
c
c
 5995 continue
c
c
c
c
c        add up contributions
c        --------------------
c
c
c
c
 8000 continue
c
c
c        charge-dependent OPE contribution
c        ---------------------------------
c
      if (mod(j,2).eq.1) go to 8020
c
c        j even
c
      v(1)=-vj(1,iman1)+2.d0*vj(1,iman5)
      v(1)=v(1)-vj(1,iman2)+2.d0*vj(1,iman6)
      v(1)=v(1)-vj(1,iman3)+2.d0*vj(1,iman7)
      v(1)=v(1)-vj(1,iman4)+2.d0*vj(1,iman8)
c
      v(2)=-vj(2,iman1)-2.d0*vj(2,iman5)
      v(2)=v(2)-vj(2,iman2)-2.d0*vj(2,iman6)
      v(2)=v(2)-vj(2,iman3)-2.d0*vj(2,iman7)
      v(2)=v(2)-vj(2,iman4)-2.d0*vj(2,iman8)
c
      do 8015 iv=3,6
      v(iv)=-vj(iv,iman1)+2.d0*vj(iv,iman5)
      v(iv)=v(iv)-vj(iv,iman2)+2.d0*vj(iv,iman6)
      v(iv)=v(iv)-vj(iv,iman3)+2.d0*vj(iv,iman7)
      v(iv)=v(iv)-vj(iv,iman4)+2.d0*vj(iv,iman8)
 8015 continue
      go to 8030
c
c        j odd
c
 8020 continue
      v(1)=-vj(1,iman1)-2.d0*vj(1,iman5)
      v(1)=v(1)-vj(1,iman2)-2.d0*vj(1,iman6)
      v(1)=v(1)-vj(1,iman3)-2.d0*vj(1,iman7)
      v(1)=v(1)-vj(1,iman4)-2.d0*vj(1,iman8)
c
      v(2)=-vj(2,iman1)+2.d0*vj(2,iman5)
      v(2)=v(2)-vj(2,iman2)+2.d0*vj(2,iman6)
      v(2)=v(2)-vj(2,iman3)+2.d0*vj(2,iman7)
      v(2)=v(2)-vj(2,iman4)+2.d0*vj(2,iman8)
c
      do 8025 iv=3,6
      v(iv)=-vj(iv,iman1)-2.d0*vj(iv,iman5)
      v(iv)=v(iv)-vj(iv,iman2)-2.d0*vj(iv,iman6)
      v(iv)=v(iv)-vj(iv,iman3)-2.d0*vj(iv,iman7)
      v(iv)=v(iv)-vj(iv,iman4)-2.d0*vj(iv,iman8)
 8025 continue
c
c
 8030 continue
c
c
      if (iman9.gt.imen) go to 8500
c
c
      if (.not.indlsj) then
      do 8105 il=iman9,imen
      do 8105 iv=1,6
 8105 v(iv)=v(iv)+vj(iv,il)
      else
c
c
c        there are contact terms
c        -----------------------
c
      if (iman9.gt.imen24) go to 8200
c
c        the non-contact terms
c
      do 8155 il=iman9,imen24
      do 8155 iv=1,6
 8155 v(iv)=v(iv)+vj(iv,il)
c
c        contact contributions
c        ---------------------
c
 8200 continue
c
c        Q^0 contacts
      do 8205 il=imen23,imen22
      do 8205 iv=1,6
 8205 vv0(iv)=vv0(iv)+vj(iv,il)
c
c        Q^2 contacts
      do 8215 il=imen21,imen15
      do 8215 iv=1,6
 8215 vv2(iv)=vv2(iv)+vj(iv,il)
c
c        Q^4 contacts
      do 8225 il=imen14,imen
      do 8225 iv=1,6
 8225 vv4(iv)=vv4(iv)+vj(iv,il)
c
c
c        ------------------------------------------------------
c        NOTE: partial-wave potentials that add-up to zero need
c        to be cutoff, because they diverge for large momenta.
c        ------------------------------------------------------
c
c        use 3d3 cutoff as default for all j.gt.5 partial waves
c
      if (j.gt.5) then
      if (cutlsj(1,15).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,15))**(2.d0*cutlsj(1,15))
     1    +(ymev/cutlsj(2,15))**(2.d0*cutlsj(1,15))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
c
      do 8275 iv=1,6
      vv0(iv)=vv0(iv)*expexp
      vv2(iv)=vv2(iv)*expexp
 8275 vv4(iv)=vv4(iv)*expexp
      go to 8400
      end if

c
c
c
c
c        look into individual partial waves and
c        multiply with partial-wave dependent cutoffs
c        --------------------------------------------
c
      j1=j+1
      go to (8310,8320,8330,8340,8350,8360),j1
c
c
c        j=0
c        ---
c        ---
c
 8310 continue
c
c        1s0
c        ---
c        Q^0 term
c
      if (cutlsj(1,1).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,1))**(2.d0*cutlsj(1,1))
     1    +(ymev/cutlsj(2,1))**(2.d0*cutlsj(1,1))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv0(1)=vv0(1)*expexp
c
c        Q^2 terms
c
      if (cutlsj(3,1).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(4,1))**(2.d0*cutlsj(3,1))
     1    +(ymev/cutlsj(4,1))**(2.d0*cutlsj(3,1))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv2(1)=vv2(1)*expexp
c
c        Q^4 terms
c
      if (cutlsj(5,1).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(6,1))**(2.d0*cutlsj(5,1))
     1    +(ymev/cutlsj(6,1))**(2.d0*cutlsj(5,1))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv4(1)=vv4(1)*expexp
c
c        3p0
c        ---
c        Q^2 term
c
      if (cutlsj(1,2).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,2))**(2.d0*cutlsj(1,2))
     1    +(ymev/cutlsj(2,2))**(2.d0*cutlsj(1,2))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv2(3)=vv2(3)*expexp
      vv0(3)=vv0(3)*expexp
c
c        Q^4 term
c
      if (cutlsj(3,2).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(4,2))**(2.d0*cutlsj(3,2))
     1    +(ymev/cutlsj(4,2))**(2.d0*cutlsj(3,2))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv4(3)=vv4(3)*expexp
c
      go to 8400
c
c
c        j=1
c        ---
c        ---
c
 8320 continue
c
c        1p1
c        ---
c        Q^2 term
c
      if (cutlsj(1,3).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,3))**(2.d0*cutlsj(1,3))
     1    +(ymev/cutlsj(2,3))**(2.d0*cutlsj(1,3))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv2(1)=vv2(1)*expexp
      vv0(1)=vv0(1)*expexp
c
c        Q^4 term
c
      if (cutlsj(3,3).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(4,3))**(2.d0*cutlsj(3,3))
     1    +(ymev/cutlsj(4,3))**(2.d0*cutlsj(3,3))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv4(1)=vv4(1)*expexp
c
c        3p1
c        ---
c        Q^2 term
c
      if (cutlsj(1,4).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,4))**(2.d0*cutlsj(1,4))
     1    +(ymev/cutlsj(2,4))**(2.d0*cutlsj(1,4))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv2(2)=vv2(2)*expexp
      vv0(2)=vv0(2)*expexp
c
c        Q^4 term
c
      if (cutlsj(3,4).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(4,4))**(2.d0*cutlsj(3,4))
     1    +(ymev/cutlsj(4,4))**(2.d0*cutlsj(3,4))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv4(2)=vv4(2)*expexp
c
c        3s1
c        ---
c        Q^0 term
c
      if (cutlsj(1,5).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,5))**(2.d0*cutlsj(1,5))
     1    +(ymev/cutlsj(2,5))**(2.d0*cutlsj(1,5))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv0(4)=vv0(4)*expexp
c
c        Q^2 terms
c
      if (cutlsj(3,5).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(4,5))**(2.d0*cutlsj(3,5))
     1    +(ymev/cutlsj(4,5))**(2.d0*cutlsj(3,5))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv2(4)=vv2(4)*expexp
c
c        Q^4 terms
c
      if (cutlsj(5,5).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(6,5))**(2.d0*cutlsj(5,5))
     1    +(ymev/cutlsj(6,5))**(2.d0*cutlsj(5,5))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv4(4)=vv4(4)*expexp
c
c        3d1
c        ---
c        Q^4 term
c
      if (cutlsj(1,6).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,6))**(2.d0*cutlsj(1,6))
     1    +(ymev/cutlsj(2,6))**(2.d0*cutlsj(1,6))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv4(3)=vv4(3)*expexp
      vv2(3)=vv2(3)*expexp
      vv0(3)=vv0(3)*expexp
c
c        3s/d1
c        -----
c        Q^2 term
c
      if (cutlsj(1,7).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,7))**(2.d0*cutlsj(1,7))
     1    +(ymev/cutlsj(2,7))**(2.d0*cutlsj(1,7))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv2(5)=vv2(5)*expexp
      vv2(6)=vv2(6)*expexp
      vv0(5)=vv0(5)*expexp
      vv0(6)=vv0(6)*expexp
c
c        Q^4 term
c
      if (cutlsj(3,7).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(4,7))**(2.d0*cutlsj(3,7))
     1    +(ymev/cutlsj(4,7))**(2.d0*cutlsj(3,7))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv4(5)=vv4(5)*expexp
      vv4(6)=vv4(6)*expexp
c
      go to 8400
c
c
c        j=2
c        ---
c        ---
c
 8330 continue
c
c        1d2
c        ---
c        Q^4 term
c
      if (cutlsj(1,8).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,8))**(2.d0*cutlsj(1,8))
     1    +(ymev/cutlsj(2,8))**(2.d0*cutlsj(1,8))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv4(1)=vv4(1)*expexp
      vv2(1)=vv2(1)*expexp
      vv0(1)=vv0(1)*expexp
c
c        3d2
c        ---
c        Q^4 term
c
      if (cutlsj(1,9).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,9))**(2.d0*cutlsj(1,9))
     1    +(ymev/cutlsj(2,9))**(2.d0*cutlsj(1,9))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv4(2)=vv4(2)*expexp
      vv2(2)=vv2(2)*expexp
      vv0(2)=vv0(2)*expexp
c
c        3p2
c        ---
c
c        Q^2 term
c
      if (cutlsj(1,10).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,10))**(2.d0*cutlsj(1,10))
     1    +(ymev/cutlsj(2,10))**(2.d0*cutlsj(1,10))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv2(4)=vv2(4)*expexp
      vv0(4)=vv0(4)*expexp
c
c        Q^4 terms
c
      if (cutlsj(3,10).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(4,10))**(2.d0*cutlsj(3,10))
     1    +(ymev/cutlsj(4,10))**(2.d0*cutlsj(3,10))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv4(4)=vv4(4)*expexp
c
c        3p/f2
c        -----
c        Q^4 term
c
      if (cutlsj(1,12).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,12))**(2.d0*cutlsj(1,12))
     1    +(ymev/cutlsj(2,12))**(2.d0*cutlsj(1,12))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
      vv4(5)=vv4(5)*expexp
      vv4(6)=vv4(6)*expexp
      vv2(5)=vv2(5)*expexp
      vv2(6)=vv2(6)*expexp
      vv0(5)=vv0(5)*expexp
      vv0(6)=vv0(6)*expexp
c
c        3f2
c        ---
c        Q^6 term
c*****************************
c        3f2 p.w. pot calculated directly
      vv0(3)=0.d0
      vv2(3)=0.d0
      vv4(3)=xx*yy*xy*clsj(1,11)*fff/fac
     1 /(dwn**6*1.d20*(2.d0*pi)**3)
c*****************************
c
      if (cutlsj(1,11).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,11))**(2.d0*cutlsj(1,11))
     1    +(ymev/cutlsj(2,11))**(2.d0*cutlsj(1,11))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
c
      vv0(3)=vv0(3)*expexp
      vv2(3)=vv2(3)*expexp
      vv4(3)=vv4(3)*expexp
c
      go to 8400
c
c
c        j=3
c        ---
c        ---
c
 8340 continue
c
c        1f3
c        ---
c        Q^6 term
c*****************************
c        1f3 p.w. pot calculated directly
      vv0(1)=0.d0
      vv2(1)=0.d0
      vv4(1)=xx*yy*xy*clsj(1,13)*fff/fac
     1 /(dwn**6*1.d20*(2.d0*pi)**3)
c*****************************
c
      if (cutlsj(1,13).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,13))**(2.d0*cutlsj(1,13))
     1    +(ymev/cutlsj(2,13))**(2.d0*cutlsj(1,13))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
c
      vv0(1)=vv0(1)*expexp
      vv2(1)=vv2(1)*expexp
      vv4(1)=vv4(1)*expexp
c
c
c
c        3f3
c        ---
c        Q^6 term
c*****************************
c        3f3 p.w. pot calculated directly
      vv0(2)=0.d0
      vv2(2)=0.d0
      vv4(2)=xx*yy*xy*clsj(1,14)*fff/fac
     1 /(dwn**6*1.d20*(2.d0*pi)**3)
c*****************************
c
      if (cutlsj(1,14).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,14))**(2.d0*cutlsj(1,14))
     1    +(ymev/cutlsj(2,14))**(2.d0*cutlsj(1,14))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
c
      vv0(2)=vv0(2)*expexp
      vv2(2)=vv2(2)*expexp
      vv4(2)=vv4(2)*expexp
c
c
c        3d3
c        ---
c        Q^4 term
c
      if (cutlsj(1,15).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,15))**(2.d0*cutlsj(1,15))
     1    +(ymev/cutlsj(2,15))**(2.d0*cutlsj(1,15))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
c
      vv0(4)=vv0(4)*expexp
      vv2(4)=vv2(4)*expexp
      vv4(4)=vv4(4)*expexp
c
c        3g3
c        ---
c        Q^8 term
c*****************************
c        3g3 p.w. pot calculated directly
      vv0(3)=0.d0
      vv2(3)=0.d0
      vv4(3)=xx*xx*yy*yy*clsj(1,16)*fff/fac
     1 /(dwn**8*1.d26*(2.d0*pi)**3)
c*****************************
c
      if (cutlsj(1,16).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,16))**(2.d0*cutlsj(1,16))
     1    +(ymev/cutlsj(2,16))**(2.d0*cutlsj(1,16))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
c
      vv0(3)=vv0(3)*expexp
      vv2(3)=vv2(3)*expexp
      vv4(3)=vv4(3)*expexp
c
c        3g-d3
c        -----
c        Q^6 term
c*****************************
c        3g-d3 p.w. pot calculated directly
      vv0(5)=0.d0
      vv2(5)=0.d0
      vv4(5)=xx*xx*yy*clsj(1,17)*fff/fac
     1 /(dwn**6*1.d20*(2.d0*pi)**3)
c*****************************
c
      if (cutlsj(1,17).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,17))**(2.d0*cutlsj(1,17))
     1    +(ymev/cutlsj(2,17))**(2.d0*cutlsj(1,17))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
c
      vv0(5)=vv0(5)*expexp
      vv2(5)=vv2(5)*expexp
      vv4(5)=vv4(5)*expexp
c
c        3d-g3
c        -----
c        Q^6 term
c*****************************
c        3d-g3 p.w. pot calculated directly
      vv0(6)=0.d0
      vv2(6)=0.d0
      vv4(6)=xx*yy*yy*clsj(1,17)*fff/fac
     1 /(dwn**6*1.d20*(2.d0*pi)**3)
c*****************************
c
      if (cutlsj(1,17).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,17))**(2.d0*cutlsj(1,17))
     1    +(ymev/cutlsj(2,17))**(2.d0*cutlsj(1,17))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
c
      vv0(6)=vv0(6)*expexp
      vv2(6)=vv2(6)*expexp
      vv4(6)=vv4(6)*expexp
c
c
      go to 8400
c
c
c        j=4
c        ---
c        ---
c
 8350 continue
c
c        1g4
c        ---
c        Q^8 term
c*****************************
c        1g4 p.w. pot calculated directly
      vv0(1)=0.d0
      vv2(1)=0.d0
      vv4(1)=xx*xx*yy*yy*clsj(1,18)*fff/fac
     1 /(dwn**8*1.d26*(2.d0*pi)**3)
c*****************************
c
      if (cutlsj(1,18).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,18))**(2.d0*cutlsj(1,18))
     1    +(ymev/cutlsj(2,18))**(2.d0*cutlsj(1,18))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
c
      vv0(1)=vv0(1)*expexp
      vv2(1)=vv2(1)*expexp
      vv4(1)=vv4(1)*expexp
c
c        3g4
c        ---
c        Q^8 term
c*****************************
c        3g4 p.w. pot calculated directly
      vv0(2)=0.d0
      vv2(2)=0.d0
      vv4(2)=xx*xx*yy*yy*clsj(1,19)*fff/fac
     1 /(dwn**8*1.d26*(2.d0*pi)**3)
c*****************************
c
      if (cutlsj(1,19).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,19))**(2.d0*cutlsj(1,19))
     1    +(ymev/cutlsj(2,19))**(2.d0*cutlsj(1,19))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
c
      vv0(2)=vv0(2)*expexp
      vv2(2)=vv2(2)*expexp
      vv4(2)=vv4(2)*expexp
c
c        3f4
c        ---
c        Q^6 term
c*****************************
c        3f4 p.w. pot calculated directly
      vv0(4)=0.d0
      vv2(4)=0.d0
      vv4(4)=xx*yy*xy*clsj(1,20)*fff/fac
     1 /(dwn**6*1.d20*(2.d0*pi)**3)
c*****************************
c
      if (cutlsj(1,20).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,20))**(2.d0*cutlsj(1,20))
     1    +(ymev/cutlsj(2,20))**(2.d0*cutlsj(1,20))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
c
      vv0(4)=vv0(4)*expexp
      vv2(4)=vv2(4)*expexp
      vv4(4)=vv4(4)*expexp
c
c        3h4
c        ---
c
      vv0(3)=0.d0
      vv2(3)=0.d0
      vv4(3)=0.d0
c
c        3h-f4
c        -----
c        Q^8 term
c*****************************
c        3h-f4 p.w. pot calculated directly
      vv0(5)=0.d0
      vv2(5)=0.d0
      vv4(5)=xx*xx*xy*yy*clsj(1,22)*fff/fac
     1 /(dwn**8*1.d26*(2.d0*pi)**3)
c*****************************
c
      if (cutlsj(1,22).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,22))**(2.d0*cutlsj(1,22))
     1    +(ymev/cutlsj(2,22))**(2.d0*cutlsj(1,22))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
c
      vv0(5)=vv0(5)*expexp
      vv2(5)=vv2(5)*expexp
      vv4(5)=vv4(5)*expexp
c
c        3f-h4
c        -----
c        Q^8 term
c*****************************
c        3f-h4 p.w. pot calculated directly
      vv0(6)=0.d0
      vv2(6)=0.d0
      vv4(6)=xx*xy*yy*yy*clsj(1,22)*fff/fac
     1 /(dwn**8*1.d26*(2.d0*pi)**3)
c*****************************
c
      if (cutlsj(1,22).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,22))**(2.d0*cutlsj(1,22))
     1    +(ymev/cutlsj(2,22))**(2.d0*cutlsj(1,22))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
c
      vv0(6)=vv0(6)*expexp
      vv2(6)=vv2(6)*expexp
      vv4(6)=vv4(6)*expexp
c
      go to 8400
c
c
c        j=5
c        ---
c        ---
c
 8360 continue
c
c        3g5
c        ---
c        Q^8 term
c*****************************
c        3g5 p.w. pot calculated directly
      vv0(4)=0.d0
      vv2(4)=0.d0
      vv4(4)=xx*xx*yy*yy*clsj(1,25)*fff/fac
     1 /(dwn**8*1.d26*(2.d0*pi)**3)
c*****************************
c
      if (cutlsj(1,25).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,25))**(2.d0*cutlsj(1,25))
     1    +(ymev/cutlsj(2,25))**(2.d0*cutlsj(1,25))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
c
      vv0(4)=vv0(4)*expexp
      vv2(4)=vv2(4)*expexp
      vv4(4)=vv4(4)*expexp
c
c        3i-g5
c        -----
c        Q^10 term
c*****************************
c        3i-g5 p.w. pot calculated directly
      vv0(5)=0.d0
      vv2(5)=0.d0
      vv4(5)=xx*xx*xx*yy*yy*clsj(1,27)*fff/fac
     1 /(dwn**10*1.d32*(2.d0*pi)**3)
c*****************************
c
      if (cutlsj(1,27).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,27))**(2.d0*cutlsj(1,27))
     1    +(ymev/cutlsj(2,27))**(2.d0*cutlsj(1,27))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
c
      vv0(5)=vv0(5)*expexp
      vv2(5)=vv2(5)*expexp
      vv4(5)=vv4(5)*expexp
c
c        3g-i5
c        -----
c        Q^10 term
c*****************************
c        3g-i5 p.w. pot calculated directly
      vv0(6)=0.d0
      vv2(6)=0.d0
      vv4(6)=xx*xx*yy*yy*yy*clsj(1,27)*fff/fac
     1 /(dwn**10*1.d32*(2.d0*pi)**3)
c*****************************
c
      if (cutlsj(1,27).eq.0.d0) then
      expexp=1.d0
      else
      expo=(xmev/cutlsj(2,27))**(2.d0*cutlsj(1,27))
     1    +(ymev/cutlsj(2,27))**(2.d0*cutlsj(1,27))
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
      end if
c
      vv0(6)=vv0(6)*expexp
      vv2(6)=vv2(6)*expexp
      vv4(6)=vv4(6)*expexp
c
c        set all other j=5 potentials to zero
c
      do 8365 iv=1,3
      vv0(iv)=0.d0
      vv2(iv)=0.d0
      vv4(iv)=0.d0
 8365 continue
c
c
c
c
c        final add up
c        ------------
c
 8400 do 8405 iv=1,6
 8405 v(iv)=v(iv)+vv0(iv)+vv2(iv)+vv4(iv)
c
      end if
c
c
c
c
 8500 if (j.eq.0.or..not.heform) go to 8900
c
c
c         translation into (combinations of) helicity states
c
c
      do 8505 i=1,4
 8505 vl(i)=v(i+2)
c
      do 8520 ii=1,4
      iii=ii+2
      v(iii)=0.d0
c
      do 8515 i=1,4
 8515 v(iii)=v(iii)+adminv(ii,i)*vl(i)
 8520 v(iii)=v(iii)*a2j1
c
c
c
c
 8900 return
      end
      subroutine chipar
c
c        chipar reads, writes, and stores the parameter for all 
c        chi-subroutines.
c
c
      implicit real*8 (a-h,o-z)
c
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
c
c
c        common block for all chi-subroutines
c
      common /cchi/ vj(32,270),c(20,270),fff,ff,f(52),aa(96),ai(19,30),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                gaa(3),fpia(3),ezz1(3),ezz2(3),ct(96),wt(96),
     3                ic(20,270),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(40,3),mggo(40,3),ima(30,40,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,270),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
      common /compar/ cb1a(3),cb2a(3),cb3a(3),cb4a(3),
     1                cd12a(3),cd3a(3),cd5a(3),cd145a(3),
     2                ce14a(3),ce15a(3),ce16a(3),ce17a(3)
c
      common /comlsj/ clsj(15,50),cutlsj(15,50),indlsj
      logical indlsj
c
c
c        further specifications
c
      dimension cc(5),cca(5)
      dimension clec(15,50)
      dimension a(1024),b(32)
      dimension ttab(5,175),tab(5,175)
      dimension topepp(5,2),topenp(5,2),topenn(5,2)
      dimension t1s0pp(5),t1s0np(5),t1s0nn(5)
      character*4 ntab(3,175)
      real*4 eps
      character*4 label
      character*4 name(3)
      integer imga(3)
      character*4 cut,cuta,fun,lsj,lec,end
      character*4 mesong(40)
      logical index
      logical indlec
      logical indca,indlca
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     1            '1-  ','1-t ','1-tt','1-st','1-ss',
     2            'c   ','ss  ','ls  ','sq  ','sk  ',
     3            'sl  ',
     4         24*'    '/
      data index/.false./
      data pi/3.141592653589793d0/
      data eps/1.e-15/
      data cut/'cut '/,cuta/'cuta'/
      data fun/'fun '/,lsj/'lsj '/,lec/'lec '/,end/'end '/
c
c
c
c
c        parameter tables
c        ----------------
c        ----------------
c
c
c        identification table
c        --------------------
c
c
      data ((ntab(i,j),i=1,3),j=1,99)/
     1 'cuta','ll  ','  ',
     2 'sq  ',' ope','p ',
     3 '    ','    ','  ',
     4 '    ','    ','  ',
     5 '    ','    ','  ',
     6 'sq  ',' ope','p ',
     7 'sq  ',' pi-','g ',
     8 'fun ','    ','  ',
     9 '    ','    ','  ',
     * '    ','    ','  ',
     1 'cuta','ll  ','  ',
     2 'c   ',' tpn','1 ',
     3 'fun ','    ','  ',
     4 'sq  ',' tpn','1 ',
     5 'fun ','    ','  ',
     6 'ss  ',' tpn','1 ',
     7 'fun ','    ','  ',
     8 'fun ','    ','  ',
     9 'c   ',' tpn','2 ',
     * 'fun ','    ','  ',
     1 'sq  ',' tpn','2 ',
     2 'fun ','    ','  ',
     3 'ss  ',' tpn','2 ',
     4 'fun ','    ','  ',
     5 'fun ','    ','  ',
     6 'c   ',' tpn','2m',
     7 'fun ','    ','  ',
     8 'c   ',' tpn','2m',
     9 'fun ','    ','  ',
     * 'sq  ',' tpn','2m',
     1 'fun ','    ','  ',
     2 'ss  ',' tpn','2m',
     3 'fun ','    ','  ',
     4 'fun ','    ','  ',
     5 'sq  ',' tpn','2m',
     6 'fun ','    ','  ',
     7 'ss  ',' tpn','2m',
     8 'fun ','    ','  ',
     9 'fun ','    ','  ',
     * 'ls  ',' tpn','2m',
     1 'fun ','    ','  ',
     2 'ls  ',' tpn','2m',
     3 'fun ','    ','  ',
     4 'c   ',' tpn','3 ',
     5 'fun ','    ','  ',
     6 'sq  ',' tpn','3 ',
     7 'fun ','    ','  ',
     8 'ss  ',' tpn','3 ',
     9 'fun ','    ','  ',
     * 'fun ','    ','  ',
     1 'c   ',' tpn','32',
     2 'fun ','    ','  ',
     3 'sq  ',' tpn','32',
     4 'fun ','    ','  ',
     5 'ss  ',' tpn','32',
     6 'fun ','    ','  ',
     7 'fun ','    ','  ',
     8 'sq  ',' tpn','32',
     9 'fun ','    ','  ',
     * 'ss  ',' tpn','32',
     1 'fun ','    ','  ',
     2 'fun ','    ','  ',
     3 'c   ',' tpn','32',
     4 'fun ','    ','  ',
     5 'sq  ',' tpn','32',
     6 'fun ','    ','  ',
     7 'ss  ',' tpn','32',
     8 'fun ','    ','  ',
     9 'fun ','    ','  ',
     * 'c   ',' tpn','32',
     1 'fun ','    ','  ',
     2 'ls  ',' tpc','m ',
     3 'fun ','    ','  ',
     4 'c   ',' tpc','m ',
     5 'fun ','    ','  ',
     6 'fun ','    ','  ',
     7 'sq  ',' tpc','m ',
     8 'fun ','    ','  ',
     9 'ss  ',' tpc','m ',
     * 'fun ','    ','  ',
     1 'fun ','    ','  ',
     2 'ls  ',' tpc','m ',
     3 'fun ','    ','  ',
     4 'c   ',' tpc','m ',
     5 'fun ','    ','  ',
     6 'c   ',' n4l','1 ',
     7 'fun ','    ','  ',
     8 'sq  ',' n4l','1 ',
     9 'fun ','    ','  ',
     * 'ss  ',' n4l','1 ',
     1 'fun ','    ','  ',
     2 'fun ','    ','  ',
     3 'sq  ',' n4l','2 ',
     4 'fun ','    ','  ',
     5 'ss  ',' n4l','2 ',
     6 'fun ','    ','  ',
     7 'fun ','    ','  ',
     8 'sq  ',' n4l','2 ',
     9 'fun ','    ','  '/
c
      data ((ntab(i,j),i=1,3),j=100,175)/
     * 'ss  ',' n4l','2 ',
     1 'fun ','    ','  ',
     2 'fun ','    ','  ',
     3 'c   ',' n4l','2 ',
     4 'fun ','    ','  ',
     5 'c   ',' n4l','2 ',
     6 'fun ','    ','  ',
     7 'c   ',' n43','p1',
     8 'fun ','    ','  ',
     9 'ss  ',' n43','p2',
     * 'fun ','    ','  ',
     1 'fun ','    ','  ',
     2 'sq  ',' n43','p3',
     3 'fun ','    ','  ',
     4 'sq  ',' n43','p3',
     5 'fun ','    ','  ',
     6 'ss  ',' n43','p4',
     7 'fun ','    ','  ',
     8 'fun ','    ','  ',
     9 'sq  ',' n43','p5',
     * 'fun ','    ','  ',
     1 'sq  ',' n43','p5',
     2 'fun ','    ','  ',
     3 'ss  ',' n43','p6',
     4 'fun ','    ','  ',
     5 'ss  ',' n43','p7',
     6 'fun ','    ','  ',
     7 'fun ','    ','  ',
     8 'sq  ',' n43','p8',
     9 'fun ','    ','  ',
     * 'sq  ',' n43','p8',
     1 'fun ','    ','  ',
     2 'ss  ',' n43','p9',
     3 'fun ','    ','  ',
     4 'fun ','    ','  ',
     5 'sq  ',' 43p','10',
     6 'fun ','    ','  ',
     7 'sq  ',' 43p','10',
     8 'fun ','    ','  ',
     9 'cuta','ll  ','  ',
     * 'lsj ',' 1S0','  ',
     1 'lsj ',' 1S0','  ',
     2 'lsj ',' 1S0','  ',
     3 'lsj ',' 1S0','  ',
     4 'lsj ',' 3P0','  ',
     5 'lsj ',' 3P0','  ',
     6 'lsj ',' 1P1','  ',
     7 'lsj ',' 1P1','  ',
     8 'lsj ',' 3P1','  ',
     9 'lsj ',' 3P1','  ',
     * 'lsj ',' 3S1','  ',
     1 'lsj ',' 3S1','  ',
     2 'lsj ',' 3S1','  ',
     3 'lsj ',' 3S1','  ',
     4 'lsj ',' 3D1','  ',
     5 'lsj ',' 3SD','1 ',
     6 'lsj ',' 3SD','1 ',
     7 'lsj ',' 3SD','1 ',
     8 'lsj ',' 1D2','  ',
     9 'lsj ',' 3D2','  ',
     * 'lsj ',' 3P2','  ',
     1 'lsj ',' 3P2','  ',
     2 'lsj ',' 3PF','2 ',
     3 'lsj ',' 3D3','  ',
     4 'lsj ',' 3F2','  ',
     5 'lsj ',' 1F3','  ',
     6 'lsj ',' 3F3','  ',
     7 'lsj ',' 3G3','  ',
     8 'lsj ',' 3D-','G3',
     9 'lsj ',' 1G4','  ',
     * 'lsj ',' 3G4','  ',
     1 'lsj ',' 3F4','  ',
     2 'lsj ',' 3F-','H4',
     3 'lsj ',' 3G5','  ',
     4 'lsj ',' 3G-','I5',
     5 'end ','para','m.'/
c
c
c        parameters
c        ----------
c
c
      data ((tab(i,j),i=1,5),j=1,99)/
     1    6.000000000d0,   0.0d0,    4.0000d0,  550.0d0,   0.0d0,
     2   -1.290000000d0,  92.4d0,  134.9766d0,    0.0d0,   0.0d0,
     3    0.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     4    0.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     5    0.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     6   -1.290000000d0,  92.4d0,  139.5702d0,    0.0d0,   0.0d0,
     7   -0.062170000d0,  92.4d0,  139.5702d0,    0.0d0,   0.0d0,
     8   10.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     9    0.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     *    0.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     1    6.000000000d0,   0.0d0,    2.0000d0,  550.0d0,   0.0d0,
     2    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     3   13.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     4    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     5   14.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     6   -1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     7   11.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     8   14.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     9    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     *   57.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     1    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     2   59.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     3   -1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     4   11.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     5   59.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     6    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     7   81.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     8    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     9   82.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     *    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     1   83.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     2   -1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     3   11.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     4   83.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     5    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     6   84.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     7   -1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     8   11.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     9   84.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     *    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     1   19.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     2    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     3   20.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     4    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     5   30.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     6    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     7   31.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     8   -1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     9   11.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     *   31.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     1    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     2   40.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     3   -1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     4   41.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     5    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     6   41.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     7   11.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     8   -1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     9   42.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     *    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     1   42.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     2   11.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     3    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     4   43.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     5   -1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     6   44.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     7    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     8   44.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     9   11.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     *    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     1   45.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     2    4.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     3   36.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     4    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     5   36.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     6   11.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     7    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     8   37.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     9   -1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     *   37.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     1   11.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     2    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     3   38.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     4    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     5   39.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     6    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     7   46.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     8   -1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     9   47.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     *    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     1   47.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     2   11.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     3    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     4   77.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     5   -1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     6   77.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     7   11.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     8    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     9   78.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0/
c
      data ((tab(i,j),i=1,5),j=100,175)/
     *   -1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     1   78.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     2   11.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     3    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     4   79.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     5    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     6   80.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     7    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     8   67.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     9   -1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     *   68.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     1   11.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     2    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     3   68.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     4    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     5   69.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     6   -1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     7   70.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     8   11.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     9    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     *   70.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     1    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     2   71.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     3    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     4   72.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     5   -1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     6   73.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     7   11.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     8    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     9   73.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     *    1.000000000d0,   0.0d0,  138.0390d0,    0.0d0,  -1.0d0,
     1   74.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     2   -1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     3   75.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     4   11.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     5    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     6   75.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     7    1.000000000d0,   0.0d0,  138.0390d0,    1.0d0,  -1.0d0,
     8   76.000000000d0,  24.0d0,    0.0000d0,    0.0d0,   0.0d0,
     9    0.000000000d0,   0.0d0,    2.0000d0,  550.0d0,   0.0d0,
     *   -0.144393750d0,   3.0d0,  550.0000d0,    0.0d0,   0.0d0,
     1    2.485068852d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     2   -2.642035856d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     3  -11.361531513d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     4    1.228900585d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     5    4.775644065d0,   3.0d0,  550.0000d0,    0.0d0,   0.0d0,
     6    0.237739558d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     7   10.812551448d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     8   -0.831549992d0,   3.0d0,  550.0000d0,    0.0d0,   0.0d0,
     9    4.033160272d0,   3.0d0,  550.0000d0,    0.0d0,   0.0d0,
     *   -0.1397643393d0,  3.0d0,  550.0000d0,    0.0d0,   0.0d0,
     1    0.846182897d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     2   -3.623910172d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     3  -19.563644316d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     4   -5.848244237d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     5    0.547696233d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     6    1.393811021d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     7    1.687849390d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     8   -1.510648724d0,   3.0d0,  550.0000d0,    0.0d0,   0.0d0,
     9   -3.335094881d0,   3.0d0,  550.0000d0,    0.0d0,   0.0d0,
     *   -0.768480318d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     1    5.489917402d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     2    0.090320081d0,   4.0d0,  550.0000d0,    0.0d0,   0.0d0,
     3   -0.565929446d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     4    1.231521499d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     5    0.583772354d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     6    0.144613481d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     7    0.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     8    0.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     9    0.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     *    0.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     1   -0.140917150d0,   2.0d0,  550.0000d0,    0.0d0,   0.0d0,
     2    0.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     3    0.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     4    0.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0,
     5    0.000000000d0,   0.0d0,    0.0000d0,    0.0d0,   0.0d0/
c
c
      data topepp/
     6   -1.290000d0,  92.4d0,  134.9766d0,    0.0d0,   0.0d0,
     7    0.000000d0,   0.0d0,  139.5702d0,    0.0d0,   0.0d0/
c
c
      data topenp/
     6   -1.290000d0,  92.4d0,  139.5702d0,    0.0d0,   0.0d0,
     7   -0.062170d0,  92.4d0,  139.5702d0,    0.0d0,   0.0d0/
c
c
      data topenn/
     6   -1.290000d0,  92.4d0,  134.9766d0,    0.0d0,   0.0d0,
     7    0.000000d0,   0.0d0,  139.5702d0,    0.0d0,   0.0d0/
c
c
      data t1s0pp/
     7   -0.143030188d0,  3.0d0,  550.0000d0,    0.0d0,   0.0d0/
c
c
      data t1s0np/
     7   -0.14439375d0,  3.0d0,  550.0000d0,    0.0d0,   0.0d0/
c
c
      data t1s0nn/
     7   -0.14367715d0,  3.0d0,  550.0000d0,    0.0d0,   0.0d0/
c
c
c        this has been the end of tables
c
c
      save
c
c
c
c
10000 format (2a4,a2,15a4)
10001 format (1h )
10002 format (1h /' contribution  par_1     par_2     mass     iso-spin
     1    iprop'/9x,'fun-type  f u n c t i o n   p a r a m e t e r s')
10003 format (2a4,a2,f15.13,4f10.6)
10004 format (1h ,2a4,a2,f12.6,f10.6,1x,f10.5,2(f7.1,3x))
10005 format (1h ,2a4,a2,f4.1,1x,2f10.5,f13.5,f10.5)
10006 format (2a4,a2,3i3)
10007 format (1h ,2a4,a2,3i3)
10008 format (1h ,47(1h-))
10009 format (2a4,a2,i3,2f10.2)
10010 format (1h ,2a4,a2,i3,2f10.2)
10011 format (//' n4lo550: Charge-Dependent Chiral NN Potential',
     1/' at Order Five using Lambda = 550 MeV, N4LO(550)')
10015 format (' input-parameter-set:'/1h ,20(1h-))
10016 format (1h ,2a4,a2,15a4)
10020 format (1h ,2a4,a2,f28.13,4f9.2)
10021 format (1h ,2a4,a2,5f10.6)
c
c
c
c
      if (index) go to 50
      index=.true.
c
      x=-1.d0
      y=-1.d0
c
c
c
c
c        maxima of certain indices related to the dimension as follows:
c        dimension c(mme,imee),ic(mice,imee),indc(mindce,imee),
c                  mgg(mge,3),mggo(mge,3),mesong(mge),vj(32,imee),
c                  ima(mee,mge,3)
c
      mge=40
      mee=30
      mme=20
      mice=20
      mindce=2
      imb=1
      ime=0
      imee=270
c        mme always ge mice, mindce
c
c        set all parameters and indices to zero or .false.
c
      do 1 int=1,3
      imga(int)=0
      indpar(int)=.false.
      do 1 mgx=1,mge
      mgg(mgx,int)=0
    1 mggo(mgx,int)=0
c
c
      do 2 il=1,imee
      do 2 mm=1,mme
      if (mm.le.mindce) indc(mm,il)=.false.
      if (mm.le.mice) ic(mm,il)=0
    2 c(mm,il)=0.d0
      endep=.false.
c
c
      pi2=pi*pi
      pi4=pi2*pi2
c
c
c
c
c        start
c        -----
c        -----
c
c
c
c        write title
c
   50 continue
      write (kwrite,10011)
      write (kwrite,10008)
      write (kwrite,10008)
c
c
c        store systematically the parameter sets
c        for pp, np, nn.
c
c
      do 8999 inter=1,3
c
c
      indca=.false.
      indlca=.false.
      indlsj=.false.
      indlec=.false.
      ilsj=0
      ilec=0
      do 55 ii=1,50
      do 55 i=1,15
      clsj(i,ii)=0.d0
      cutlsj(i,ii)=0.d0
   55 clec(i,ii)=0.d0
c
c
c        fix index-parameter concerning over-all factor 
c        and the sfr cutoff 
c
      ift(inter)=1
      ezz1(inter)=0.d0
      ezz2(inter)=700.d0
c
c**** write (kwrite,10010) name,ift(inter),ezz1(inter),ezz2(inter)
      iftyp=ift(inter)
c
c
c        set parameters for numerical integration
c
      mint(inter)=8
      maxt(inter)=96
c
c**** write (kwrite,10007) name,mint(inter),maxt(inter)
c
c        nucleon mass
c
      go to (51,52,53), inter
c        mass used for pp
   51 wn=938.2720d0
      go to 54
c        mass used for np
   52 wn=938.9183d0
      go to 54
c        mass used for nn
   53 wn=939.5654d0
   54 continue
c**** write (kwrite,10004) name,wn
      wnq=wn*wn
      dwn=1.d0/wn
      dwnq=dwn*dwn
      wnn(inter)=wn
c
c
c        ga and fpi
c
      ga=1.29d0
      fpi=92.4d0
c
c**** write (kwrite,10004) name,ga,fpi
      ga2=ga*ga
      ga4=ga2*ga2
      ga6=ga4*ga2
      fpi=fpi*dwn
      fpi2=fpi*fpi
      fpi3=fpi2*fpi
      fpi4=fpi2*fpi2
      fpi6=fpi4*fpi2
      gaa(inter)=ga2
      fpia(inter)=fpi2
c
c        fix the LECs of the pi-N Lagrangian
c
c        the c_i LECs
      cc(1)=-1.10d0
      cc(2)=3.57d0
      cc(3)=-5.54d0
      cc(4)=4.17d0
c**** write (kwrite,10021) name,cc
      cb1a(inter)=cc(1)*wn*1.d-3
      cb2a(inter)=cc(2)*wn*1.d-3
      cb3a(inter)=cc(3)*wn*1.d-3
      cb4a(inter)=cc(4)*wn*1.d-3
c
c        the d_i LECs
      cc(1)=6.18d0
      cc(2)=-8.91d0
      cc(3)=0.86d0
      cc(4)=-12.18d0
c**** write (kwrite,10021) name,cc
      cd12a(inter)=cc(1)*wnq*1.d-6
      cd3a(inter)=cc(2)*wnq*1.d-6
      cd5a(inter)=cc(3)*wnq*1.d-6
      cd145a(inter)=cc(4)*wnq*1.d-6
c
c        the e_i LECs
      cc(1)=1.18d0
      cc(2)=0.0d0
      cc(3)=0.0d0
      cc(4)=-0.18d0
c**** write (kwrite,10021) name,cc
      ce14a(inter)=cc(1)*wnq*wn*1.d-9
      ce15a(inter)=cc(2)*wnq*wn*1.d-9
      ce16a(inter)=cc(3)*wnq*wn*1.d-9
      ce17a(inter)=cc(4)*wnq*wn*1.d-9
c
      cb1=cb1a(inter)
      cb2=cb2a(inter)
      cb3=cb3a(inter)
      cb4=cb4a(inter)
      cd12=cd12a(inter)
      cd3=cd3a(inter)
      cd5=cd5a(inter)
      cd145=cd145a(inter)
      ce14=ce14a(inter)
      ce17=ce17a(inter)
c
c
c
c
c        prepare table
c
      do 56 ll=1,175
      do 56 i=1,5
   56 ttab(i,ll)=tab(i,ll)
c
c
c        charge-dependent modifications for pp
c
      if (inter.eq.1) then
      do 57 i=1,5
      ttab(i,6)=topepp(i,1)
      ttab(i,7)=topepp(i,2)
   57 ttab(i,140)=t1s0pp(i)
      end if
c
c
c        charge-dependent modifications for np
c
      if (inter.eq.2) then
      do 58 i=1,5
      ttab(i,6)=topenp(i,1)
      ttab(i,7)=topenp(i,2)
   58 ttab(i,140)=t1s0np(i)
      end if
c
c
c        charge-dependent modifications for nn
c
      if (inter.eq.3) then
      do 59 i=1,5
      ttab(i,6)=topenn(i,1)
      ttab(i,7)=topenn(i,2)
   59 ttab(i,140)=t1s0nn(i)
      end if
c
c
c
c
c        get parameters from tables, line by line
c        ----------------------------------------
c        ----------------------------------------
c
c
c
      line=0
c
   61 line=line+1
      do i=1,5
      if (i.le.3) then
      name(i)=ntab(i,line)
      end if
      cc(i)=ttab(i,line)
      end do
c
c        check if end of input
c
      if (name(1).eq.end) go to 7000
c
c        check if lsj or lec
c
      if (name(1).eq.lsj) go to 6000
      if (name(1).eq.lec) go to 6500
c
c        check if data-card just read contains cut-off or
c        function parameters
c
      if (name(1).eq.cut.or.name(1).eq.fun) go to 70
c
      if (name(1).eq.cuta) then
c**** write (kwrite,10005) name,cc
      indca=.true.
      do i=1,5
      cca(i)=cc(i)
      end do
      go to 61
      end if

c
c
c
c
c        write parameters which are no cut-off or function parameters
c        ------------------------------------------------------------
c
c
c
c
c**** write (kwrite,10004) name,cc
c
c        find out number of contribution mg
c
      do 63 mg=1,mge
      if (name(1).eq.mesong(mg)) go to 64
   63 continue
      go to 9000
c
c
c
c
c        store parameters which are no cut-off or function parameters
c        ------------------------------------------------------------
c
c
c
c
   64 ime=ime+1
      if (ime.gt.imee) go to 9011
      mgg(mg,inter)=mgg(mg,inter)+1
      m=mgg(mg,inter)
      if (m.gt.mee) go to 9001
      ima(m,mg,inter)=ime
      if (m.ne.1) go to 65
      imga(inter)=imga(inter)+1
      mggo(imga(inter),inter)=mg
   65 continue
c
c
      c(1,ime)=cc(1)
c
c
      if (mg.le.10) then
      c(1,ime)=c(1,ime)*4.d0*pi
      end if

      if (mg.le.3.and.cc(2).ne.0.d0) then
      c(1,ime)=(cc(1)/cc(2)*wn)**2
      if (cc(1).lt.0.d0) c(1,ime)=-c(1,ime)
      end if
c
c
      if (mg.ge.6.and.mg.le.10) then
c        store coupling constant f*g
      c(3,ime)=cc(2)*c(1,ime)
c        store coupling constant f**2
      c(2,ime)=cc(2)*c(3,ime)
      if (mg.eq.10)
     1  c(1,ime)=c(1,ime)+c(3,ime)*2.d0+c(2,ime)
      end if
c
c
      if (mg.ge.11.and.cc(2).ne.0.d0) then
      c(1,ime)=(cc(1)/(2.d0*cc(2))*wn)**2
      if (cc(1).lt.0.d0) c(1,ime)=-c(1,ime)
      end if
c
c
c        store meson mass square in units of nucleon mass square
      c(4,ime)=cc(3)*cc(3)*dwnq
c
c        test iso-spin
      icc=cc(4)
      if (icc.ne.0.and.icc.ne.1) go to 9004
c         store isospin as logical constant
      if (icc.eq.1) indc(1,ime)=.true.
c        store and test iprsp
      icc=cc(5)
      ic(1,ime)=icc
      if (iabs(ic(1,ime)).gt.1) go to 9005
c
c        index values for further storing
      mi=4
      mm=5
c
c
c        check if there is a `cutall' cutoff
c
      if (indca) then
      name(1)=cut
      do i=1,5
      cc(i)=cca(i)
      end do
      go to 72
      else
      go to 61
      end if
c
c
c
c
c        write cut-off or function parameters
c        ------------------------------------
c
c
c
c
   70 continue
c**** write (kwrite,10005) name,cc
c
   72 continue
c
c
c
c
c        store parameters
c        ----------------
c
c
c
      ityp=cc(1)
c
      if (ityp.eq.0) go to 5995
      if (ityp.lt.1.or.ityp.gt.84) go to 9002
c
      im=ime
c
c        store typ of cut-off or function
      ic(mi,im)=ityp
c
      if (ityp.le.10) then
c        store and test typ of propagator of cut-off
      ic(mi+1,im)=cc(2)
      if (ic(mi+1,im).lt.0.or.ic(mi+1,im).gt.1) go to 9006
      end if
c
      go to (100,100,300,9002,500,600,9002,9002,9002,1000,
     1 1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,
     2 2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,
     3 3100,3200,3300,3400,3500,3600,3700,3800,3900,4000,
     4 4100,4200,4300,4400,4500,4600,4700,4800,4900,5000,
     5 5100,5200,5300,5400,5500,5600,5700,5800,25900,26000,
     6 26100,26200,26300,26400,26500,26600,26700,26800,26900,27000,
     7 27100,27200,27300,27400,27500,27600,27700,27800,27900,28000,
     8 28100,28200,28300,28400),ityp
c
c
c
c
c        cut-off of dipole type
c        **********************
c
c
c        store and test exponent of cut-off
  100 ic(mi+2,im)=cc(3)
      if (ic(mi+2,im).lt.0) go to 9009
      if (ic(mi+2,im).gt.0) go to 101
c        exponent is zero, omit cut-off
      ic(mi,im)=0
      ic(mi+1,im)=0
      go to 5995
c        store cut-off mass for denominator
  101 c(mm+1,im)=cc(4)*cc(4)*dwnq
c        store numerator of cut-off
      c(mm,im)=c(mm+1,im)
      if (ityp.eq.2)     c(mm,im)=c(mm,im)-c(4,im)
      mi=mi+3
      mm=mm+2
      go to 5995
c
c
c
c
c        exponential form factor of momentum transfer
c        ********************************************
c
c
c        check exponent
  300 if (cc(3).lt.0.d0) go to 9009
      if (cc(3).gt.0.d0) go to 301
c        exponent is zero, omit cutoff
      ic (mi,im)=0
      ic (mi+1,im)=0
      go to 5995
c        store exponent
  301 c(mm+1,im)=cc(3)
c        compute constant factor for argument of exponential function
      c(mm,im)=wnq/(cc(4)*cc(4))
      mi=mi+2
      mm=mm+2
      go to 5995
c
c
c
c
c        sharp cutoff in x and y
c        ***********************
c
c
  500 c(mm,im)=cc(4)*dwn
      mi=mi+2
      mm=mm+1
      go to 5995
c
c
c
c
c        exponential form factor of xx and yy
c        ************************************
c
c
c        check exponent
  600 if (cc(3).lt.0.d0) go to 9009
      if (cc(3).gt.0.d0) go to 601
c        exponent is zero, omit cutoff
      ic (mi,im)=0
      ic (mi+1,im)=0
      go to 5995
c        store exponent
  601 c(mm+1,im)=cc(3)
c        compute constant factor for argument of exponential function
      c(mm,im)=wnq/(cc(4)*cc(4))
      mi=mi+2
      mm=mm+2
      go to 5995
c
c
c
c
c        pi-gamma potential
c        ******************
c
c
 1000 c(mm,im)=cc(3)
      mi=mi+2
      mm=mm+1
      go to 5995
c
c
c
c
c        function q^2 (momentum-transfer squared)
c        ************
c
c
 1100 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function k^2 (average-momentum squared)
c        ************
c
c
 1200 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function 1 for tpn1 (=NLO)
c        **************************
c
c
 1300 c(mm,im)=-1.d0/(384.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function 2 for tpn1
c        *******************
c
c
 1400 c(mm,im)=-3.d0*ga4/(64.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn2 (=N^2LO), function 1
c        *************************
c
c
 1500 c(mm,im)=-3.d0*ga2/(16.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn2, function 2
c        ****************
c
c
 1600 c(mm,im)=-ga2/(128.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn2, function 3
c        ****************
c
c
 1700 c(mm,im)=9.d0*ga4/(512.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn2, function 4
c        ****************
c
c
 1800 c(mm,im)=-ga2/(32.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c
c        tpn2, function 5
c        ****************
c
c
 1900 c(mm,im)=6.d0*ga4/(64.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c        tpn2, function 6
c        ****************
c
c
 2000 c(mm,im)=2.d0*ga2*(1.d0-ga2)/(64.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        function q^4 (momentum-transfer to the power of 4)
c        ************
c
 2100 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c        function k^4 (average-momentum to the power of 4)
c        ************
c
 2200 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c        function +q^2*k^2
c        *****************
c
 2300 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c        function  (\vec q x \vec k)^2
c        *****************************
c
 2400 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c        function xy
c        ***********
c
 2500 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c        function xx+yy
c        **************
c
 2600 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c        function xx*xx+yy*yy
c        ********************
c
 2700 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c        function xx
c        ***********
c
 2800 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c        function yy
c        ***********
c
 2900 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn3 (= N^3LO with one loop), function 1
c        ****************************************
c
c
 3000 c(mm,im)=3.d0/(16.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn3, function 2
c        ****************
c
c
 3100 continue
      c(mm,im)=cb4*cb4/(96.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c        function 1.d0
c        *************
c
 3200 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c        function 1-q^2/8-k^2/2
c        *************************
c
 3300 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c        function 1-q^2/8
c        ****************
c
 3400 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c        function 1+k^2/2
c        ****************
c
 3500 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn3, function 3
c        ****************
c
c
 3600 c(mm,im)=-cb4/(192.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn3, function 4
c        ****************
c
c
 3700 c(mm,im)=cb4/(192.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn3, function 5
c        ****************
c
c
 3800 continue
      c(mm,im)=cb2*ga2/(8.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn3, function 6
c        ****************
c
c
 3900 c(mm,im)=-ga2/(32.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn32 (=N^3LO with 2 loops), function 1
c        ***************************************
c
c
 4000 c(mm,im)= - (-2.d0/pi) * 3.d0*ga4/(4096.d0*pi*fpi6)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn32, function 2
c        *****************
c
c
 4100 c(mm,im)=  (2.d0/pi) * ga4/(4096.d0*pi*fpi6)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn32, function 3
c        *****************
c
c
 4200 c(mm,im)=  (2.d0/pi) * 2.d0*ga6/((8.d0*pi*fpi2)**3) * 
     1                  48.d0*pi2*fpi2/ga4*cd145
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn32, function 4
c        *****************
c
c
 4300 c(mm,im)= - (-2.d0/pi) * 2.d0/3.d0/((8.d0*pi*fpi2)**3)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn32, function 5
c        *****************
c
c
 4400 c(mm,im)=2.d0*ga6/(pi*(8.d0*pi*fpi2)**3)
      ic(mi+1,im)=cc(2)
      mi=mi+2
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn32, function 6
c        *****************
c
c
 4500 c(mm,im)=2.d0/(3.d0*pi*(8.d0*pi*fpi2)**3)
      ic(mi+1,im)=cc(2)
      mi=mi+2
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn4, function 1
c        ****************
c
 4600 c(mm,im) = -(-2.d0/pi) * ( 1.d0 / (((4.d0*fpi)**6)*pi*pi) )
      ic(mi+1,im)=cc(2)
      mi=mi+2
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn4, function 2
c        ****************
c
 4700 c(mm,im) = (2.d0/pi) * ( cb4 * ga2 / (((4.d0*fpi)**6)*pi*pi) )
      ic(mi+1,im)=cc(2)
      mi=mi+2
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn3m (= N^3LO, 1/M^2 terms), function 1
c        ****************************************
c
c
 4800 c(mm,im)=-ga4/(32.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn3m, function 2
c        *****************
c
c
 4900 c(mm,im)=-1.d0/(768.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn3m, function 3
c        *****************
c
c
 5000 c(mm,im)=ga4/(32.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn3m, function 4
c        *****************
c
c
 5100 c(mm,im)=1.d0/(1536.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn3m, function 5
c        *****************
c
c
 5200 c(mm,im)=1.d0/(256.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn3m, function 6
c        *****************
c
c
 5300 c(mm,im)=ga4/(4.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn3m, function 7
c        *****************
c
c
 5400 c(mm,im)=ga4/(32.d0*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn2c (correction for our it 2pi), function 1
c        *********************************************
c
c
 5500 c(mm,im)=ga4/(128.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn2c (correction for our it 2pi), function 2
c        *********************************************
c
c
 5600 c(mm,im)=-ga4/(256.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn2 (=N^2LO), function 1, the non-1/M part of 1500
c        ***************************************************
c
c
 5700 c(mm,im)=-3.d0*ga2/(16.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn2m (=N^2LO), function 1, the 1/M part of 1500
c        ***********************************************
c
c
 5800 c(mm,im)=-3.d0*ga2/(16.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn2, function 4, the non-1/M part of 1800 (tpn2)
c        ************************************************
c
c
25900 c(mm,im)=-ga2/(32.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c        tpn2m, function 4, the 1/M part of 1800 (tpn2m)
c        **********************************************
c
c
26000 c(mm,im)=-ga2/(32.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c
c
c
c
c           <<-- Vc4 BEGIN
c           **************
c
26100 c(mm,im) = 3.d0*ga4/(512.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c           Vc4 END -->>
c           ************
c
c
c
c           <<-- Wc4 BEGIN
c           **************
c
26200 c(mm,im) = ga2/(128.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c           Wc4 END -->>
c           ************
c
c
c
c
c           <<-- Vt4 BEGIN
c           **************
c
26300 c(mm,im) = 9.d0*ga2/(512.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c           Vt4 END -->>
c           ************
c
c
c
c           <<-- Wt4 BEGIN
c           **************
c
26400 c(mm,im) = -ga2/(256.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c           Wt4 END -->>
c           ************
c
c
c
c           <<-- V_LS4 BEGIN
c           ****************
c
26500 c(mm,im) = (-2.d0)*(-3.d0*ga4)/(64.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c           V_LS4 END -->>
c           **************
c
c
c
c           <<-- W_LS4 BEGIN
c           ****************
c
26600 c(mm,im) = (-2.d0)*(-ga2)*(1.d0 - ga2)/(64.d0*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
c
c           W_LS4 END -->>
c           **************
c
c
c
c
c           <<-- Wc13 BEGIN
c           ***************
c     
26700 c(mm,im) = -(-2.d0/pi) * ( ga4*cb4 / (((4.d0*fpi)**6)*pi*pi) )
c           
      ic(mi+1,im)=cc(2)
      mi=mi+2
      mm=mm+1
      go to 5995
c
c           Wc13 END -->>
c           *************
c
c
c
c           <<-- Vs13 BEGIN
c           ***************
c     
26800 c(mm,im) = -(2.d0/pi) * (ga4*cb4 / (((4.d0*fpi)**6)*pi*pi) )
c     
      ic(mi+1,im)=cc(2)
      mi=mi+2
      mm=mm+1
      go to 5995
c
c           Vs13 END -->>
c           *************
c
c
c
c
c           <<-- Vt13 BEGIN
c           ***************
c     
26900 c(mm,im) = -(2.d0/pi) * (ga4*cb4 / (((4.d0*fpi)**6)*pi*pi) )
c     
      ic(mi+1,im)=cc(2)
      mi=mi+2
      mm=mm+1
      go to 5995
c
c           Vt13 END -->>
c           *************
c
c
c
c           <<-- Ws13 BEGIN
c           ***************
c     
27000 c(mm,im) = -(2.d0/pi) * (ga4 / (((4.d0*fpi)**6)*pi*pi) )
c     
      ic(mi+1,im)=cc(2)
      mi=mi+2
      mm=mm+1
      go to 5995
c
c           Ws13 END -->>
c           *************
c
c
c
c           <<-- Wt13 BEGIN
c           ***************
c
c     
27100 c(mm,im) = 1.d0
      ic(mi+1,im)=cc(2)
      mi=mi+2
      mm=mm+1
      go to 5995
c
c           Wt13 END -->>
c           *************
c
c
c
c
c           <<-- Vs14 BEGIN
c           ***************
c     
27200 c(mm,im) = -(-2.d0/pi) * (ga4*cb4 / (((4.d0*fpi)**6)*pi*pi) )
c     
      ic(mi+1,im)=cc(2)
      mi=mi+2
      mm=mm+1
      go to 5995
c
c           Vs14 END -->>
c           *************
c
c
c
c
c           <<-- Vs12 BEGIN
c           ***************
c     
27300 c(mm,im) = -(2.d0/pi) * ( ga2*cb4 / (((4.d0*fpi)**6)*pi*pi) )
c     
      ic(mi+1,im)=cc(2)
      mi=mi+2
      mm=mm+1
      go to 5995
c
c           Vs12 END -->>
c           *************
c
c
c
c           <<-- Vt12 BEGIN
c           ***************
c     
27400 c(mm,im) = -(2.d0/pi) * ( ga2*cb4 / (((4.d0*fpi)**6)*pi*pi) )
c     
      ic(mi+1,im)=cc(2)
      mi=mi+2
      mm=mm+1
      go to 5995
c
c           Vt12 END -->>
c           *************
c
c
c
c           <<-- Ws12 BEGIN
c           ***************
c     
27500 c(mm,im) = -(2.d0/pi) * (ga2 / (((4.d0*fpi)**6)*pi*pi) )
c     
      ic(mi+1,im)=cc(2)
      mi=mi+2
      mm=mm+1
      go to 5995
c
c           Ws12 END -->>
c           *************
c
c
c
c           <<-- Wt12 BEGIN
c           ***************
c     
27600 c(mm,im) = -(2.d0/pi) * (ga2 / (((4.d0*fpi)**6)*pi*pi) )
c     
      ic(mi+1,im)=cc(2)
      mi=mi+2
      mm=mm+1
      go to 5995
c
c           Wt12 END -->>
c           *************
c
c
c
c           <<-- Vt_class2 BEGIN
c           ********************
c     
27700 c(mm,im) = -(2.d0/pi) * (ga4*(cb4-cb3) / (((4.d0*fpi)**6)*pi*pi) )
c     
      mi=mi+1
      mm=mm+1
      go to 5995
c
c           Vt_class2 END -->>
c           ******************
c
c
c
c           <<-- Wt_class2 BEGIN
c           ********************
c     
27800 c(mm,im) = -(2.d0/pi) * (ga2 / (((4.d0*fpi)**6)*pi*pi) )
c     
      mi=mi+1
      mm=mm+1
      go to 5995
c
c           Wt_class2 END -->>
c           ******************
c
c
c
c           <<-- Vc_class2 BEGIN
c           ********************
c     
27900 c(mm,im) = 1.d0
c     
      mi=mi+1
      mm=mm+1
      go to 5995
c
c           Vc_class2 END -->>
c           ******************
c
c
c
c           <<-- Wc_class2 BEGIN
c           ********************
c     
28000 c(mm,im) = 1.d0
c     
      mi=mi+1
      mm=mm+1
      go to 5995
c
c           Wc_class2 END -->>
c           ******************
c
c
c           <<-- Vc leading relativistic corrections BEGIN
c           **********************************************
c     
28100 c(mm,im) = 3.d0 * ga4 / (128.d0*pi*(fpi**4) )
c     
      mi=mi+1
      mm=mm+1
      go to 5995
c
c           Vc leading relativistic corrections END -->>
c           ********************************************
c
c
c
c           <<-- Wc leading relativistic corrections BEGIN
c           **********************************************
c     
28200 c(mm,im) = ga2 / (64.d0*pi*(fpi**4) )
c     
      mi=mi+1
      mm=mm+1
      go to 5995
c
c           Wc leading relativistic corrections END -->>
c           ********************************************
c
c
c
c           <<-- Vt leading relativistic corrections BEGIN
c           **********************************************
c     
28300 c(mm,im) = 3.d0 * ga4 / (256.d0*pi*(fpi**4) )
c     
      mi=mi+1
      mm=mm+1
      go to 5995
c
c           Vt leading relativistic corrections END -->>
c           ********************************************
c
c
c
c           <<-- Wt leading relativistic corrections BEGIN
c           **********************************************
c     
28400 c(mm,im) = ga2 / (128.d0*pi*(fpi**4) )
c     
      mi=mi+1
      mm=mm+1
      go to 5995
c
c           Wt leading relativistic corrections END -->>
c           ********************************************
c
c
c
c
c        end cut-offs and functions
c        **************************
c
c        test dimensions
 5995 if (mi.gt.mice.or.mm.gt.mme) go to 9010
c
      if (indlca) go to 7800
c
      go to 61
c
c
c
c
c        partial wave LEC's
c        ------------------
c
c
 6000 continue
c**** write (kwrite,10020) name,cc
      indlsj=.true.
      ilsj=ilsj+1
      if (ilsj.le.4) iilsj=1
      if (ilsj.ge.5.and.ilsj.le.6) iilsj=2
      if (ilsj.ge.7.and.ilsj.le.8) iilsj=3
      if (ilsj.ge.9.and.ilsj.le.10) iilsj=4
      if (ilsj.ge.11.and.ilsj.le.14) iilsj=5
      if (ilsj.ge.15.and.ilsj.le.15) iilsj=6
      if (ilsj.ge.16.and.ilsj.le.18) iilsj=7
      if (ilsj.ge.19.and.ilsj.le.19) iilsj=8
      if (ilsj.ge.20.and.ilsj.le.20) iilsj=9
      if (ilsj.ge.21.and.ilsj.le.22) iilsj=10
      if (ilsj.ge.23.and.ilsj.le.23) iilsj=12
      if (ilsj.ge.24.and.ilsj.le.24) iilsj=15
      if (ilsj.ge.25.and.ilsj.le.25) iilsj=11
      if (ilsj.ge.26.and.ilsj.le.26) iilsj=13
      if (ilsj.ge.27.and.ilsj.le.27) iilsj=14
      if (ilsj.ge.28.and.ilsj.le.28) iilsj=16
      if (ilsj.ge.29.and.ilsj.le.29) iilsj=17
      if (ilsj.ge.30.and.ilsj.le.30) iilsj=18
      if (ilsj.ge.31.and.ilsj.le.31) iilsj=19
      if (ilsj.ge.32.and.ilsj.le.32) iilsj=20
      if (ilsj.ge.33.and.ilsj.le.33) iilsj=22
      if (ilsj.ge.34.and.ilsj.le.34) iilsj=25
      if (ilsj.ge.35.and.ilsj.le.35) iilsj=27
      if (ilsj.eq.1) iord=0
      if (ilsj.eq.5) iord=0
      if (ilsj.eq.7) iord=0
      if (ilsj.eq.9) iord=0
      if (ilsj.eq.11) iord=0
      if (ilsj.eq.15) iord=0
      if (ilsj.eq.16) iord=0
      if (ilsj.eq.19) iord=0
      if (ilsj.eq.20) iord=0
      if (ilsj.eq.21) iord=0
      if (ilsj.eq.23) iord=0
      if (ilsj.ge.24) iord=0
      iord=iord+1
      clsj(iord,iilsj)=cc(1)
      cutlsj(2*iord-1,iilsj)=cc(2)
      cutlsj(2*iord,iilsj)=cc(3)
      go to 61
c
c
c
c
c        lec LEC's
c        ---------
c
c
 6500 continue
c**** write (kwrite,10020) name,cc
      indlec=.true.
      ilec=ilec+1
      go to (6510,6510,6522,6540,6542,6544),ilec
 6510 do 6515 i=1,5
 6515 clec(ilec,i)=cc(i)
      go to 61
 6522 do 6523 i=1,2
 6523 clec(2,i+5)=cc(i)
      go to 61
 6540 do 6541 i=1,5
 6541 clec(3,i)=cc(i)
      go to 61
 6542 do 6543 i=1,5
 6543 clec(3,i+5)=cc(i)
      go to 61
 6544 do 6545 i=1,5
 6545 clec(3,i+10)=cc(i)
      go to 61
c
c
c
c
c        conclusions
c        -----------
c        -----------
c
c
c        write end
 7000 continue
c**** write (kwrite,10004) name
c**** write (kwrite,10008)
c**** write (kwrite,10008)
c
      if (indlsj) go to 7100
      if (indlec) go to 7500
      go to 8995
c
c
c        determine the low-energy constants (clec)
c        from the partial wave constants (clsj)
c
c
c        LEC's for Q^0 (LO)
c        ------------------
c
 7100 clec(1,1)=(clsj(1,1)+3.d0*clsj(1,5))*0.25d0/(4.d0*pi)
      clec(1,2)=(clsj(1,5)  -   clsj(1,1))*0.25d0/(4.d0*pi)
c
c
c        LEC's for Q^2 (NLO)
c        -------------------
c
c        vector b
c
c        1s0
      b(1)=clsj(2,1)
c        3p0
      b(2)=clsj(1,2)
c        1p1
      b(3)=clsj(1,3)
c        3p1
      b(4)=clsj(1,4)
c        3s1
      b(5)=clsj(2,5)
c        3s-d1
      b(6)=clsj(1,7)
c        3p2
      b(7)=clsj(1,10)
c
c
      do 7205 i=1,7
 7205 b(i)=b(i)/(4.d0*pi)
c
c
c        matrix a for the C parameters
c
c        1. column
      a(1)=1.d0
      a(2)=-2.d0/3.d0
      a(3)=a(2)
      a(4)=a(2)
      a(5)=1.d0
      a(6)=0.d0
      a(7)=a(2)
c
c        2. column
      a(8)=0.25d0
      a(9)=1.d0/6.d0
      a(10)=a(9)
      a(11)=a(9)
      a(12)=0.25d0
      a(13)=0.d0
      a(14)=a(9)
c
c        3. column
      a(15)=-3.d0
      a(16)=-2.d0/3.d0
      a(17)=2.d0
      a(18)=a(16)
      a(19)=1.d0
      a(20)=0.d0
      a(21)=a(16)
c
c        4. column
      a(22)=-0.75d0
      a(23)=1.d0/6.d0
      a(24)=-0.5d0
      a(25)=a(23)
      a(26)=0.25d0
      a(27)=0.d0
      a(28)=a(23)
c
c        5. column
      a(29)=0.d0
      a(30)=-2.d0/3.d0
      a(31)=0.d0
      a(32)=-1.d0/3.d0
      a(33)=0.d0
      a(34)=0.d0
      a(35)=1.d0/3.d0
c
c        6. column
      a(36)=-1.d0
      a(37)=2.d0
      a(38)=2.d0/3.d0
      a(39)=-4.d0/3.d0
      a(40)=1.d0/3.d0
      a(41)=-2.d0*dsqrt(2.d0)/3.d0
      a(42)=0.d0
c
c        7. column
      a(43)=-0.25d0
      a(44)=-0.5d0
      a(45)=-1.d0/6.d0
      a(46)=1.d0/3.d0
      a(47)=1.d0/12.d0
      a(48)=-dsqrt(2.d0)/6.d0
      a(49)=0.d0
c
c
c
c
      call dgelg (b,a,7,1,eps,ier)
c
      if (ier.ne.0) write (kwrite,19500) ier
19500 format (///' warning in chipar. the error index of dgelg is',
     1 ' ier =',i5/' for the calculation of the C parameters.'///)
c

      do 7255 i=1,7
 7255 clec(2,i)=b(i)
c
c
c        LEC's for Q^4 (N^3LO)
c        -------------------
c
c        vector b
c
c        1s0
      b(1)=clsj(3,1)
c        1s0
      b(2)=clsj(4,1)
c        3p0
      b(3)=clsj(2,2)
c        1p1
      b(4)=clsj(2,3)
c        3p1
      b(5)=clsj(2,4)
c        3s1
      b(6)=clsj(3,5)
c        3s1
      b(7)=clsj(4,5)
c        3d1
      b(8)=clsj(1,6)
c        3s-d1
      b(9)=clsj(2,7)
c        3s-d1
      b(10)=clsj(3,7)
c        1d2
      b(11)=clsj(1,8)
c        3d2
      b(12)=clsj(1,9)
c        3p2
      b(13)=clsj(2,10)
c        3p-f2
      b(14)=clsj(1,12)
c        3d3
      b(15)=clsj(1,15)
c
c
      do 7305 i=1,15
 7305 b(i)=b(i)/(4.d0*pi)
c
c
c        matrix a for the D parameters
c
c        1. column
      a(1)=1.d0
      a(2)=10.d0/3.d0
      a(3)=-4.d0/3.d0
      a(4)=a(3)
      a(5)=a(3)
      a(6)=1.d0
      a(7)=a(2)
      a(8)=8.d0/15.d0
      a(9)=0.d0
      a(10)=0.d0
      a(11)=a(8)
      a(12)=a(8)
      a(13)=a(3)
      a(14)=0.d0
      a(15)=a(8)
c
c        2. column
      a(16)=1.d0/16.d0
      a(17)=5.d0/24.d0
      a(18)=1.d0/12.d0
      a(19)=a(18)
      a(20)=a(18)
      a(21)=a(16)
      a(22)=a(17)
      a(23)=1.d0/30.d0
      a(24)=0.d0
      a(25)=0.d0
      a(26)=a(23)
      a(27)=a(23)
      a(28)=a(18)
      a(29)=0.d0
      a(30)=a(23)
c
c        3. column
      a(31)=1.d0/4.d0
      a(32)=1.d0/6.d0
      a(33)=0.d0
      a(34)=0.d0
      a(35)=0.d0
      a(36)=a(31)
      a(37)=a(32)
      a(38)=-2.d0/15.d0
      a(39)=0.d0
      a(40)=0.d0
      a(41)=a(38)
      a(42)=a(38)
      a(43)=0.d0
      a(44)=0.d0
      a(45)=a(38)
c
c        4. column
      a(46)=0.d0
      a(47)=2.d0/3.d0
      a(48)=0.d0
      a(49)=0.d0
      a(50)=0.d0
      a(51)=0.d0
      a(52)=a(47)
      a(53)=-2.d0/15.d0
      a(54)=0.d0
      a(55)=0.d0
      a(56)=a(53)
      a(57)=a(53)
      a(58)=0.d0
      a(59)=0.d0
      a(60)=a(53)
c
c        5. column
      a(61)=-3.d0
      a(62)=-10.d0
      a(63)=-4.d0/3.d0
      a(64)=4.d0
      a(65)=a(63)
      a(66)=1.d0
      a(67)=10.d0/3.d0
      a(68)=8.d0/15.d0
      a(69)=0.d0
      a(70)=0.d0
      a(71)=-8.d0/5.d0
      a(72)=a(68)
      a(73)=a(63)
      a(74)=0.d0
      a(75)=a(68)
c
c        6. column
      a(76)=-3.d0/16.d0
      a(77)=-5.d0/8.d0
      a(78)=1.d0/12.d0
      a(79)=-1.d0/4.d0
      a(80)=a(78)
      a(81)=1.d0/16.d0
      a(82)=5.d0/24.d0
      a(83)=1.d0/30.d0
      a(84)=0.d0
      a(85)=0.d0
      a(86)=-1.d0/10.d0
      a(87)=a(83)
      a(88)=a(78)
      a(89)=0.d0
      a(90)=a(83)
c
c        7. column
      a(91)=-3.d0/4.d0
      a(92)=-1.d0/2.d0
      a(93)=0.d0
      a(94)=0.d0
      a(95)=0.d0
      a(96)=1.d0/4.d0
      a(97)=1.d0/6.d0
      a(98)=-2.d0/15.d0
      a(99)=0.d0
      a(100)=0.d0
      a(101)=2.d0/5.d0
      a(102)=a(98)
      a(103)=0.d0
      a(104)=0.d0
      a(105)=a(98)
c
c        8. column
      a(106)=0.d0
      a(107)=-2.d0
      a(108)=0.d0
      a(109)=0.d0
      a(110)=0.d0
      a(111)=0.d0
      a(112)=2.d0/3.d0
      a(113)=-2.d0/15.d0
      a(114)=0.d0
      a(115)=0.d0
      a(116)=2.d0/5.d0
      a(117)=a(113)
      a(118)=0.d0
      a(119)=0.d0
      a(120)=a(113)
c
c        9. column
      a(121)=0.d0
      a(122)=0.d0
      a(123)=-2.d0/3.d0
      a(124)=0.d0
      a(125)=-1.d0/3.d0
      a(126)=0.d0
      a(127)=0.d0
      a(128)=2.d0/5.d0
      a(129)=0.d0
      a(130)=0.d0
      a(131)=0.d0
      a(132)=2.d0/15.d0
      a(133)=1.d0/3.d0
      a(134)=0.d0
      a(135)=-4.d0/15.d0
c
c        10. column
      a(136)=0.d0
      a(137)=0.d0
      a(138)=-1.d0/6.d0
      a(139)=0.d0
      a(140)=-1.d0/12.d0
      a(141)=0.d0
      a(142)=0.d0
      a(143)=-1.d0/10.d0
      a(144)=0.d0
      a(145)=0.d0
      a(146)=0.d0
      a(147)=-1.d0/30.d0
      a(148)=1.d0/12.d0
      a(149)=0.d0
      a(150)=1.d0/15.d0
c
c        11. column
      a(151)=-1.d0
      a(152)=-10.d0/3.d0
      a(153)=8.d0/3.d0
      a(154)=4.d0/3.d0
      a(155)=-2.d0
      a(156)=1.d0/3.d0
      a(157)=10.d0/9.d0
      a(158)=-4.d0/9.d0
      a(159)=-2.d0*dsqrt(2.d0)/3.d0
      a(160)=-14.d0*dsqrt(2.d0)/9.d0
      a(161)=-8.d0/15.d0
      a(162)=4.d0/5.d0
      a(163)=-2.d0/15.d0
      a(164)=4.d0*dsqrt(6.d0)/15.d0
      a(165)=0.d0
c
c        12. column
      a(166)=-1.d0/4.d0
      a(167)=-1.d0/6.d0
      a(168)=1.d0/3.d0
      a(169)=0.d0
      a(170)=-1.d0/6.d0
      a(171)=1.d0/12.d0
      a(172)=1.d0/18.d0
      a(173)=1.d0/9.d0
      a(174)=-dsqrt(2.d0)/6.d0
      a(175)=dsqrt(2.d0)/18.d0
      a(176)=2.d0/15.d0
      a(177)=-1.d0/5.d0
      a(178)=1.d0/30.d0
      a(179)=-dsqrt(6.d0)/15.d0
      a(180)=0.d0
c
c        13. column
      a(181)=-1.d0/4.d0
      a(182)=-1.d0/6.d0
      a(183)=-1.d0/3.d0
      a(184)=0.d0
      a(185)=1.d0/6.d0
      a(186)=1.d0/12.d0
      a(187)=1.d0/18.d0
      a(188)=1.d0/9.d0
      a(189)=-dsqrt(2.d0)/6.d0
      a(190)=dsqrt(2.d0)/18.d0
      a(191)=2.d0/15.d0
      a(192)=-1.d0/5.d0
      a(193)=-1.d0/30.d0
      a(194)=dsqrt(6.d0)/15.d0
      a(195)=0.d0
c
c        14. column
      a(196)=-1.d0/16.d0
      a(197)=-5.d0/24.d0
      a(198)=-1.d0/6.d0
      a(199)=-1.d0/12.d0
      a(200)=1.d0/8.d0
      a(201)=1.d0/48.d0
      a(202)=5.d0/72.d0
      a(203)=-1.d0/36.d0
      a(204)=-dsqrt(2.d0)/24.d0
      a(205)=-7.d0*dsqrt(2.d0)/72.d0
      a(206)=-1.d0/30.d0
      a(207)=1.d0/20.d0
      a(208)=1.d0/120.d0
      a(209)=-dsqrt(6.d0)/60.d0
      a(210)=0.d0
c
c        15. column
      a(211)=0.d0
      a(212)=-2.d0/3.d0
      a(213)=0.d0
      a(214)=0.d0
      a(215)=0.d0
      a(216)=0.d0
      a(217)=2.d0/9.d0
      a(218)=-16.d0/45.d0
      a(219)=0.d0
      a(220)=2.d0*dsqrt(2.d0)/9.d0
      a(221)=2.d0/15.d0
      a(222)=4.d0/15.d0
      a(223)=0.d0
      a(224)=0.d0
      a(225)=-2.d0/15.d0
c
c
c
c
      call dgelg (b,a,15,1,eps,ier)
c
      if (ier.ne.0) write (kwrite,19501) ier
19501 format (///' warning in chipar. the error index of dgelg is',
     1 ' ier =',i5/' for the calculation of the D parameters.'///)
c
c
      do 7355 i=1,15
 7355 clec(3,i)=b(i)
c
c
c
c        write LEC's
c        -----------
c
c
 7500 continue
c**** write (kwrite,10100)
10100 format (//' Low energy parameters (LEC):'/
     1          ' ----------------------------'/)
c
c        Q^0 (LO)
c**** write (kwrite,10101) (clec(1,i),i=1,2)
10101 format ('lec  CS,CT',2f10.6)
c
c        Q^2 (NLO)
c**** write (kwrite,10102) (clec(2,i),i=1,7)
10102 format ('lec  C_i  ',5f10.6)
c
c        Q^4 (N^3LO)
c**** write (kwrite,10103) (clec(3,i),i=1,15)
10103 format ('lec  D_i  ',5f10.6)
c
c
c
c
c        store LEC's appropriately
c        -------------------------
c
c
      iorder=0
 7600 iorder=iorder+1
c
c
      mg=10
      iterm=0
 7700 iterm=iterm+1
c
c
      if (iorder.eq.1.and.iterm.gt.2) go to 7600
      if (iorder.eq.2.and.iterm.gt.7) go to 7600
c
c
      mg=mg+1
c
      if (iorder.eq.2) then
      if (iterm.eq.2) mg=mg-1
      if (iterm.eq.4) mg=mg-1
      end if
c
      if (iorder.eq.3) then
      if (iterm.eq.2) mg=mg-1
      if (iterm.eq.3) mg=mg-1
      if (iterm.eq.4) mg=mg-1
      if (iterm.eq.6) mg=mg-1
      if (iterm.eq.7) mg=mg-1
      if (iterm.eq.8) mg=mg-1
      if (iterm.eq.10) mg=mg-1
      if (iterm.eq.12) mg=mg-1
      if (iterm.eq.14) mg=mg-1
      end if
c
c
      ime=ime+1
      if (ime.gt.imee) go to 9011
      mgg(mg,inter)=mgg(mg,inter)+1
      m=mgg(mg,inter)
      if (m.gt.mee) go to 9001
      ima(m,mg,inter)=ime
      if (m.eq.1) then
      imga(inter)=imga(inter)+1
      mggo(imga(inter),inter)=mg
      end if
c
c
      c(1,ime)=clec(iorder,iterm)*wnq*1.d-2
      ic(1,ime)=-1
c
c
      mi=4
      mm=5
c
c
      if (indca) then
      indlca=.true.
      name(1)=cut
      do i=1,5
      cc(i)=cca(i)
      end do
      go to 72
      end if
c
c
 7800 indlca=.false.
c
c
      if (iorder.eq.2) then
      c(1,ime)=c(1,ime)*wnq*1.d-6
      if (iterm.le.4) then
      imod=mod(iterm,2)
      if (imod.eq.0) imod=2
      ic(mi,ime)=10+imod
      end if
      end if
c
c
      if (iorder.eq.3) then
      c(1,ime)=c(1,ime)*(wnq*1.d-6)**2
      if (iterm.le.8) then
      imod=mod(iterm,4)
      if (imod.eq.0) imod=4
      ic(mi,ime)=20+imod
      end if
      if (iterm.ge.9.and.iterm.le.14) then
      imod=mod(iterm,2)
      if (imod.eq.0) imod=2
      ic(mi,ime)=10+imod
      end if
      end if
c
c
 7900 if (iterm.lt.15) go to 7700
      if (iorder.lt.3) go to 7600
c
c
 8995 imaa(inter)=imb
      imea(inter)=ime
      imb=ime+1
c
 8999 continue
c        this has been the end of the inter loop
c
      return
c
c
c
c        errors
c        ------
c        ------
c
c
c
c
 9000 write (kwrite,19000) name(1)
19000 format (1h ////' error in chipar:  contribution  ',a4,'   does not
     1 exist in this program.'/' execution terminated.'////)
      go to 9999
c
c
 9001 write (kwrite,19001)
19001 format (1h ////' error in chipar:too many contributions within a g
     1roup with respect to'/' the given dimensions. execution terminated
     2.'////)
      go to 9999
c
c
 9002 write (kwrite,19002) cc(1)
19002 format (1h ////' error in chipar: cut/fun typ',f10.4,'  does not e
     1xist in this program.'/' execution terminated.'////)
      go to 9999
c
c
 9003 write (kwrite,19003) iftyp
19003 format (1h ////' error in chipar: factor typ has the non-permissib
     1le value',i4,' .'/' execution terminated.'////)
      go to 9999
c
c
 9004 write (kwrite,19004) cc(4)
19004 format (1h ////' error in chipar: isospin has the non-permissible
     1value',f10.4,'  .'/' execution terminated.'////)
      go to 9999
c
c
 9005 write (kwrite,19005) cc(5)
19005 format (1h ////' error in chipar: iprop/spe has the non-permissibl
     1e value',f10.4,'  .'/' execution terminated.'////)
      go to 9999
c
c
 9006 write (kwrite,19006) cc(2)
19006 format (1h ////' error in chipar: the index for the propagator of
     1the cut-off has the'/' non-permissible value',f10.4,'  . execution
     2 terminated.'////)
      go to 9999
c
c
 9009 write (kwrite,19009)
19009 format (1h ////' error in chipar: the exponent of the cut-off is l
     1ess than zero.'/' execution terminated.'////)
      go to 9999
c
c
 9010 write (kwrite,19010)
19010 format (1h ////' error in chipar: too many cut/fun parameters with
     1 respect to the given'/' dimensions. execution terminated.'////)
      go to 9999
c
c
 9011 write (kwrite,19011)
19011 format (1h ////' error in chipar:  too many contr. with respect to
     1 the dimensions given'/' to this program. execution terminated.'
     2////)
      go to 9999
c
c
 9999 stop
      end
      subroutine chistr (icase,max,mex)
c
c        chistr computes the structure of one-boson-exchanges
c
c
      implicit real*8 (a-h,o-z)
c
c
c        common blocks
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
c
c
c        common block for all chi-subroutines
c
      common /cchi/ vj(32,270),c(20,270),fff,ff,f(52),aa(96),ai(19,30),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                gaa(3),fpia(3),ezz1(3),ezz2(3),ct(96),wt(96),
     3                ic(20,270),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(40,3),mggo(40,3),ima(30,40,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,270),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c     further specifications
c
      dimension vv(32)
      dimension tt(2,3)
      logical index
      logical indiso
      data jj/-1/
      data index/.false./
      save
c
c
c
c
      if (index) go to 50
      index=.true.
c
c
      do 1 ii=1,3
      tt(1,ii)=1.d0
    1 tt(2,ii)=-3.d0
c
c
c
c
c
   50 do 1095 m=max,mex
      im=ima(m,mg,inter)
c
c
      if (mc.ne.1) go to 60
c
c
c
c
c        call integrals
c        --------------
c
c
c
c
      call chiai
c
c
c
c
   60 if (mc.lt.1) mc=1
c
      if (c(mc,im).eq.0.d0) go to 1095
c
c
c
c
c        nn-nn helicity amplitudes /combinations/
c        ----------------------------------------
c
c
c
c
c        basic structure (a factor of 2 is included in v5 and v6)
c
c
      ive=6
c
      vv(1)=f(1)*ai(1,m)+f(2)*ai(2,m)
      vv(2)=f(3)*ai(1,m)+f(4)*ai(3,m)
      vv(3)=f(5)*ai(1,m)+f(6)*ai(2,m)
      vv(4)=f(4)*ai(1,m)+f(3)*ai(3,m)
      vv(5)=f(7)*ai(4,m)
      vv(6)=f(8)*ai(4,m)
c
c
      go to (1000,120,130,140),icase
c
c
c        additional terms required for the tensor coupling
c        of the rho-meson or for certain operators,
c        like, the spin-orbit operator (`ls  ')
c
c
  120 vv(1)=vv(1)+f(9)*ai(5,m)
      vv(2)=vv(2)+f(10)*ai(2,m)+f(9)*ai(6,m)
      vv(3)=vv(3)+f(10)*ai(5,m)
      vv(4)=vv(4)+f(9)*ai(2,m)+f(10)*ai(6,m)
         e1=f(11)*ai(7,m)
      vv(5)=vv(5)+e1
      vv(6)=vv(6)+e1
      go to 1000
c
c
c        additional terms in case of 2+ mesons
c        not needed here
c
c
  130 continue
      go to 1000
c
c
c        additional terms needed for the sigma-l operator (`sl  ')
c
c
  140 vv(1)=vv(1)+f(6)*ai(5,m)
      vv(2)=vv(2)+f(1)*ai(5,m)+f(9)*ai(6,m)
      vv(3)=vv(3)+f(1)*ai(11,m)
      vv(4)=vv(4)+f(9)*ai(2,m)+f(1)*ai(12,m)
      vv(5)=vv(5)+f(6)*ai(13,m)
      vv(6)=vv(6)+f(6)*ai(13,m)
c
c
c
c
 1000 continue
c
c
c
c
c        set certain cases to zero
c
      if (j.ne.0) go to 1021
      vv(2)=0.d0
      vv(4)=0.d0
      vv(5)=0.d0
      vv(6)=0.d0
c

 1021 mmod=mod(j,2)
      if (.not.sing.or.(mmod.eq.1.and.inter.ne.2)) vv(1)=0.d0
      if (.not.trip.or.(mmod.eq.0.and.inter.ne.2)) vv(2)=0.d0
      if (coup.and.(mmod.eq.0.or.inter.eq.2)) go to 1030
      do 1025 iv=3,6
 1025 vv(iv)=0.d0
c
 1030 continue
c
c
c        transformation into lsj-formalism
c 
      if (j.eq.jj) go to 1035
      jj=j
      aj=dfloat(j)
      aj1=dfloat(j+1)
      d2j1=1.d0/dfloat(2*j+1)
      arjj1=dsqrt(aj*aj1)
c
 1035 v3=vv(3)
      v4=vv(4)
      v5=vv(5)
      v6=vv(6)
      v34=-arjj1*(v3-v4)
      v56=arjj1*(v5+v6)
      vv(3)=d2j1*(aj1*v3+aj*v4-v56)
      vv(4)=d2j1*(aj*v3+aj1*v4+v56)
      vv(5)=d2j1*(v34-aj1*v5+aj*v6)
      vv(6)=d2j1*(v34+aj*v5-aj1*v6)
c
c
c        possible different sign depending on the convention used
      vv(5)=-vv(5)
      vv(6)=-vv(6)
c
c
c
c
c        multiply with factors
c        ---------------------
c
c
c
c
 1040 is=mod(j,2)+1
      it=mod(is,2)+1
      indiso=indc(1,im)
      cmc=c(mc,im)
      fc=fff*ff*cmc
      do 1045 iv=1,ive
c
c        multiply with coupling-constant and factors fff and ff
c
      vv(iv)=vv(iv)*fc
c
c        multiply with isospin factor
c
      if (.not.indiso) go to 1045
      if (iv.eq.2) go to 1043
      vv(iv)=vv(iv)*tt(is,inter)
      go to 1045
 1043 vv(iv)=vv(iv)*tt(it,inter)
c
c     add up in case of several couplings for one meson-exchange
c     and store
 1045 vj(iv,im)=vj(iv,im)+vv(iv)
c
c
 1095 continue
c
c
      return
      end
      subroutine chiai
c
c        chiai integrates over theta
c
c
      implicit real*8 (a-h,o-z)
c
      common /cpot/   v(6),xmev,ymev
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
c
c
c        common block for all chi-subroutines
c
      common /cchi/ vj(32,270),c(20,270),fff,ff,f(52),aa(96),ai(19,30),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                gaa(3),fpia(3),ezz1(3),ezz2(3),ct(96),wt(96),
     3                ic(20,270),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(40,3),mggo(40,3),ima(30,40,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,270),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c
c        further specifications
      dimension gi(7)
c
      dimension pj(7,96)
      real*4 axy2,aomq,am
      logical indj
      data ige/7/
      data nnt/-1/,iinter/-1/,jj/-1/
      save
c
c
c
c
      if (inter.eq.iinter) go to 60
      iinter=inter
      min=mint(inter)
      max=maxt(inter)
c
      igeint=7
c
      wn=wnn(inter)
      dwn=1.d0/wn
      wnq=wn*wn
c
c
c
c
   60 if (j.eq.jj) go to 70
      jj=j
      indj=.false.
c
c
      aj=dfloat(j)
      aj1=dfloat(j+1)
      dj1=1.d0/aj1
      ajdj1=aj*dj1
      aaj=dsqrt(ajdj1)
c
c
      aj2=dfloat(j+2)
      ajm1=dfloat(j-1)
c
c
      ajj1=aj*aj1
      ajj2=ajm1*aj2
      ajjb=aj*ajm1
c
      aajj=0.d0
      if (j.gt.1)
     1aajj=aj/dsqrt(ajj1*ajj2)
c
      aaj1=aajj*ajm1
      aaj2=aajj*aj1
      aaj3=aajj*2.d0
c
      if (j.gt.1) go to 62
      aajj=0.d0
      go to 63
   62 aajj=1.d0/(aj1*dsqrt(ajj2))
c
   63 aaj4=aajj*ajjb
      aaj5=aajj*aj1*2.d0
      aaj6=aajj*(ajj1+2.d0)
      aaj7=aajj*ajj2
c
c
c
c
c        find out appropriate number of gauss-points
c        -------------------------------------------
c
c
   70 c4=c(4,im)
      if (c4.eq.0.d0) then
      c4=(138.*dwn)**2
      end if
      iprsp=ic(1,im)
c
c
c        compute am
c
      axy2=xy2
      if (iprsp.ne.1) go to 91
      aomq=eem12+c4
      go to 92
   91 aomq=xxpyy+c4
c
   92 am=axy2/aomq
c
c
c        compute number of gausspoints (nt)
c
c
      if (am.gt.0.999) go to 94
c
c
      if (am.gt.0.85) am=am**(-alog(1.-am)-0.9)
c
c
      nt=float(min)/(1.-am)+0.9
c
c
      if (nt.gt.max) nt=max
      go to 95
c
c
   94 nt=max
c
c
   95 nt=nt+j
c
c        compute nt, which is suitable for gset
c
      if (nt.le.16) go to 98
      if (nt.gt.24) go to 96
      nt=4*(nt/4)
      go to 98
   96 if (nt.gt.48) go to 97
      nt=8*(nt/8)
      go to 98
   97 nt=16*(nt/16)
      if (nt.gt.96) nt=96
c
   98 if (nt.eq.nnt.and.indj) go to 100
c
c
c
c
c        call gauss-points
c        -----------------
c
c
c
c
      call gset (-1.d0,1.d0,nt,ct,wt)
      nnt=nt
c
c
c
c
c        call legendre-polynoms if necessary
c        -----------------------------------
c
c
c
c
      indxy=.false.
      indj=.true.
      do 99 i=1,nt
      t=ct(i)
      call legp (pj(1,i),pj(3,i),t,j)
      pj(2,i)=pj(1,i)*t
      pj(4,i)=pj(2,i)*t
      pj(6,i)=pj(4,i)*t
      pj(5,i)=pj(3,i)*t
   99 pj(7,i)=pj(5,i)*t
c
c
c
c
c        call integrand
c        --------------
c
c
c
c
  100 call chiaa
c
c
c
c
c        prepare for integration
c
c
c
c
      do 2001 ig=1,igeint
 2001 gi(ig)=0.d0
c
c
c
c
c        integration-loop of theta
c        -------------------------
c
c
c
c
      do 2005 i=1,nt
      do 2005 ig=1,igeint
 2005 gi(ig)=gi(ig)+pj(ig,i)*aa(i)
c
c
c
      if (j.ne.0) go to 2010
      gi(3)=0.d0
      gi(5)=0.d0
      gi(7)=0.d0
c
c
c
c
c        combinations of integrals
c        -------------------------
c
c
c
c
 2010 ai(1,m)=gi(1)
c
      ai(2,m)=gi(2)
      ai(3,m)= ajdj1*gi(2)+dj1*gi(3)
      gi23m  =gi(2)-gi(3)
      ai(4,m)=aaj*gi23m
c
c
      ai(5,m)=gi(4)
      ai(6,m)= ajdj1*gi(4)+dj1*gi(5)
      gi45m  =gi(4)-gi(5)
      ai(7,m)=aaj*gi45m
c
c
      ai( 8,m)= aaj1*gi(4)-aaj2*gi(1)+aaj3*gi(5)
      aai1    = aaj4*gi(4)+aaj5*gi(1)-aaj6*gi(5)
      aai2    = aaj7*gi23m
      ai( 9,m)= aai2+aai1
      ai(10,m)= aai2-aai1
c
c
      ai(11,m)=gi(6)
      ai(12,m)=ajdj1*gi(6)+dj1*gi(7)
      ai(13,m)=aaj*(gi(6)-gi(7))
c
c
      return
      end
      subroutine chiaa
c
c        chiaa computes propagators, cutoffs, and functions
c
c
      implicit real*8 (a-h,o-z)
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
c
c
c        common block for all chi-subroutines
c
      common /cchi/ vj(32,270),c(20,270),fff,ff,f(52),aa(96),ai(19,30),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                gaa(3),fpia(3),ezz1(3),ezz2(3),ct(96),wt(96),
     3                ic(20,270),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(40,3),mggo(40,3),ima(30,40,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,270),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
      common /compar/ cb1a(3),cb2a(3),cb3a(3),cb4a(3),
     1                cd12a(3),cd3a(3),cd5a(3),cd145a(3),
     2                ce14a(3),ce15a(3),ce16a(3),ce17a(3)
c
c
      common /crrr/ rrr
c
c
c
c        further specifications
      dimension deltaq(96,7)
      dimension ell(96),cpa(96),cpaa(96)
      dimension xi(96),wxi(96)
      dimension amu(96),wmu(96)
      dimension amu2(96),amu4(96),amu6(96)
      dimension a3mu(96),w3mu(96)
      dimension a3mu2(96),a3mu4(96),a3mu6(96)
      dimension v40(96), v41(96), v42(96), v43(96)
      dimension v44a(96), v44b(96), v45(96)
      dimension v46(96), v47(96), v67(96), v68(96), v69(96), v70(96)
      dimension v71(96), v72(96), v73(96), v74(96), v75(96), v76(96)
      dimension v77(96), v78(96), v79(96), v80(96)
      logical iv40, iv41, iv42, iv43
      logical iv44, iv45
      logical iv46, iv47, iv67, iv68, iv69, iv70
      logical iv71, iv72, iv73, iv74, iv75, iv76
      logical iv77, iv78, iv79, iv80
      logical indla,indmu, i3ndmu
      data iinter/-1/
      data cc4/-1.d0/,cc445/-1.d0/,cc4mu/-1.d0/,c3c4mu/-1.d0/
      data nnxi/-1/
      data pi/3.141592653589793d0/
      data iv40/.false./,iv41/.false./,iv42/.false./,iv43/.false./
      data iv44/.false./,iv45/.false./
      data iv46/.false./,iv47/.false./,iv67/.false./,iv68/.false./
      data iv69/.false./,iv70/.false./,iv71/.false./,iv72/.false./
      data iv73/.false./,iv74/.false./,iv75/.false./,iv76/.false./
      data iv77/.false./,iv78/.false./,iv79/.false./,iv80/.false./
      save
c
c
c
c
      if (inter.eq.iinter) go to 10
      iinter=inter
      ga2=gaa(inter)
      ga4=ga2*ga2
      fpi2=fpia(inter)
      pih=pi/2.d0
      pi2=pi*pi
      wn=wnn(inter)
      dwn=1.d0/wn
      wnq=wn*wn
      ez2=ezz2(inter)*dwn
c
      cb1=cb1a(inter)
      cb2=cb2a(inter)
      cb3=cb3a(inter)
      cb4=cb4a(inter)
      cd12=cd12a(inter)
      cd3=cd3a(inter)
      cd5=cd5a(inter)
      cd145=cd145a(inter)
      ce14=ce14a(inter)
      ce17=ce17a(inter)
c
      indxy=.false.
      indla=.false.
      indmu=.false.
      i3ndmu=.false.
      iv40=.false.
      iv41=.false.
      iv42=.false.
      iv43=.false.
      iv44=.false.
      iv45=.false.
      iv46=.false.
      iv47=.false.
      iv77=.false.
      iv78=.false.
      iv79=.false.
      iv80=.false.
      iv67=.false.
      iv68=.false.
      iv69=.false.
      iv70=.false.
      iv71=.false.
      iv72=.false.
      iv73=.false.
      iv74=.false.
      iv75=.false.
      iv76=.false.
   10 continue
c
c
c
c
c        delta square
c        ------------
c
c
c
c
      if (indxy) go to 50
      indxy=.true.
      indla=.false.
      do 15 i=1,nt
      xy2t=xy2*ct(i)
c
c
c        function  -q^2 (- momentum-transfer-squared)
c        --------------
c
c        retardation ignored
c
      deltaq(i,1)=xy2t-xxpyy
c
c        retardation incorporated
c
      deltaq(i,2)=xy2t-eem12
c
c
c        function  +k^2 (average-momentum squared)
c        --------------
c
      deltaq(i,3)=(xy2t+xxpyy)*0.25d0
c
c        function  q^4 (momentum-transfer to the power of 4)
c        -------------
c
      deltaq(i,4)=deltaq(i,1)*deltaq(i,1)
c
c        function  k^4 (average-momentum to the power of 4)
c        -------------
c
      deltaq(i,5)=deltaq(i,3)*deltaq(i,3)
c
c        function  +q^2*k^2
c        -----------------
c
      deltaq(i,6)=-deltaq(i,1)*deltaq(i,3)
c
c        function  (\vec q x \vec k)^2
c        -----------------------------
c
      deltaq(i,7)=xx*yy*(1.d0-ct(i)*ct(i))
c
   15 continue
      go to 50
c
c
c
c     calculate ell, cpa, and cpaa
c
   20 indla=.true.
      cc4=c4
c
      wpi=dsqrt(c4)
      wpi2=2.d0*wpi
      wpi2q=wpi2*wpi2
      ez2q=ez2*ez2
      ez3=ez2-wpi2
      ez4=ez2*wpi2
      ez5q=ez2q-wpi2q
      ez5=dsqrt(ez5q)
c
      do 25 i=1,nt
      akk=-deltaq(i,1)
      ak=dsqrt(akk)
      radi=wpi2q+akk
      root=dsqrt(radi)
c
      if (ez2.eq.0.d0) then
      ell(i)=root*dlog((root+ak)/wpi2)/ak
      cpa(i)=datan(ak/wpi2)/(2.d0*ak)
      else
      ell(i)=root*dlog((radi*ez2q+akk*ez5q+2.d0*ez2*ak*root*ez5)
     1       /(wpi2q*(ez2q+akk)))/(2.d0*ak)
      cpa(i)=datan(ak*ez3/(akk+ez4))/(2.d0*ak)
      end if
c
      cpaa(i)=(2.d0*c4+akk)*cpa(i)
   25 continue
      go to 6000
c
c
c        prepare mu integration. 2-pion case
c
   30 indmu=.true.
      iv40=.false.
      iv41=.false.
      iv42=.false.
      iv43=.false.
      iv44=.false.
      iv45=.false.
      iv46=.false.
      iv47=.false.
      iv77=.false.
      iv78=.false.
      iv79=.false.
      iv80=.false.
c
      cc4mu=c4
c
c     pion mass
      wpi=dsqrt(c4)
c
c     set lower integration boundary      
      amulow = 2.d0*wpi
c
      nmu=96
c
      cmu=1.0d0
c
      if (ez2.eq.0.d0) then
c           [2m_pi, +infinity] case
            acut=1.d0
            call gset (0.d0,acut,nmu,amu,wmu)
            do 32 ii=1, nmu
                  xmu=pih*amu(ii)
c
c                 transformed gauss point
                  amu(ii)=dtan(xmu)*cmu + amulow
c
c                 transformed gauss weight
                  dcmu=1.d0/dcos(xmu)
                  wmu(ii)=pih*cmu*dcmu*dcmu*wmu(ii)
   32       continue      
      else
c           [2m_pi, cutoff] case, where cutoff = ez2
            acut=ez2
            call gset (amulow,acut,nmu,amu,wmu)
      end if
c
c     precalculate mu^2, mu^4 and mu^6 to speed up further calculations
      do 35 ii=1,nmu
            amu2(ii)=amu(ii)*amu(ii)
            amu4(ii)=amu2(ii)*amu2(ii)
            amu6(ii)=amu4(ii)*amu2(ii)
   35 continue
      go to 6000
c
c
c        prepare mu integration. 3-pion case.
c
   40 i3ndmu=.true.
      iv67=.false.
      iv68=.false.
      iv69=.false.
      iv70=.false.
      iv71=.false.
      iv72=.false.
      iv73=.false.
      iv74=.false.
      iv75=.false.
      iv76=.false.
c
      c3c4mu=c4
c
c     pion mass
      wpi=dsqrt(c4)
c
c     set lower integration boundary      
      amulow = 3.d0*wpi
c
      n3mu=24
c      
      cmu=1.0d0
c
      if (ez2.eq.0.d0) then
c           [3m_pi, +infinity] case
            acut=1.d0
            call gset (0.d0,acut,n3mu,a3mu,w3mu)
            do 42 ii=1, n3mu
                  xmu=pih*a3mu(ii)
c
c                 transformed gauss point
                  a3mu(ii)=dtan(xmu)*cmu + amulow
c
c                 transformed gauss weight
                  dcmu=1.d0/dcos(xmu)
                  w3mu(ii)=pih*cmu*dcmu*dcmu*w3mu(ii)
   42       continue      
      else
c           [3m_pi, cutoff] case, where cutoff = ez2
            acut=ez2
            call gset (amulow,acut,n3mu,a3mu,w3mu)
      end if
c
c     precalculate mu^2, mu^4 and mu^6 to speed up further calculations
      do 45 ii=1,n3mu
            a3mu2(ii)=a3mu(ii)*a3mu(ii)
            a3mu4(ii)=a3mu2(ii)*a3mu2(ii)
            a3mu6(ii)=a3mu4(ii)*a3mu2(ii)
   45 continue
      go to 6000
c
c
c
c        propagator
c        ----------
c        ----------
c
c
c
c
   50 c4=c(4,im)
      iprsp=ic(1,im)
      if (iprsp.lt.0) go to 60
      iret=iprsp+1
c
c         propagator for the nn case
      do 55 i=1,nt
   55 aa(i)=wt(i)/(c4-deltaq(i,iret))
      go to 80
c
c
c        "no propagator"
c
   60 do 65 i=1,nt
   65 aa(i)=wt(i)
c
c
   80 continue
c
c
c
c
c
c        cut-offs and functions
c        ----------------------
c        ----------------------
c
c
c
c
      mi=4
      mm=5
c
c
 5999 ityp=ic(mi,im)
      if (ityp.eq.0) go to 8000
      if (ityp.le.10) then
      iprspc=ic(mi+1,im)
      iret=iprspc+1
      end if
 6000 go to (100,100,300,9002,500,600,9002,9002,9002,1000,
     1 1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,
     2 2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,
     3 3100,3200,3300,3400,3500,3600,3700,3800,3900,4000,
     4 4100,4200,4300,4400,4500,4600,4700,4800,4900,5000,
     5 5100,5200,5300,5400,5500,5600,5700,5800,25900,26000,
     6 26100,26200,26300,26400,26500,26600,26700,26800,26900,27000,
     7 27100,27200,27300,27400,27500,27600,27700,27800,27900,28000,
     8 28100,28200,28300,28400),ityp
c
c
c
c
c
c
c        cut-off of dipole type
c        **********************
c
c
  100 c5=c(mm,im)
      c6=c(mm+1,im)
      nexp=ic(mi+2,im)
c
      do 105 i=1,nt
c
      aaa=c5/(c6-deltaq(i,iret))
c     -------------------------
c
      do 105 ii=1,nexp
  105 aa(i)=aa(i)*aaa
c
c
      mi=mi+3
      mm=mm+2
      go to 5999
c
c
c
c
c        exponential form factor of momentum transfer
c        ********************************************
c
c
  300 c5=c(mm,im)
      c6=c(mm+1,im)
      do 305 i=1,nt
c
      expo=(c5*dabs(deltaq(i,iret)))**c6
c     ----------------------------
c
      if (expo.gt.rrr) expo=rrr
c
      aa(i)=aa(i)*dexp(-expo)
c     ----------------------
c
  305 continue
      mi=mi+2
      mm=mm+2
      go to 5999
c
c
c
c
c        sharp cutoff of x and y
c        ***********************
c
c
  500 c5=c(mm,im)
c
      if (x.gt.c5.or.y.gt.c5) then
c     ----------------------------
      do 505 i=1,nt
  505 aa(i)=0.d0
      end if
c
      mi=mi+2
      mm=mm+1
      go to 5999
c
c
c
c
c        exponential form factor of xx and yy
c        ************************************
c
c
  600 c5=c(mm,im)
      c6=c(mm+1,im)
c
      expo=(c5*xx)**c6+(c5*yy)**c6
c     ----------------------------
      if (expo.gt.rrr) expo=rrr
      expexp=dexp(-expo)
c     ------------------
c
      do 605 i=1,nt
  605 aa(i)=aa(i)*expexp
      mi=mi+2
      mm=mm+2
      go to 5999
c
c
c
c
c
c        pi-gamma potential
c        ******************
c
c
 1000 c5=c(mm,im)
      do 1055 i=1,nt
      betaq=-deltaq(i,1)/c4
      betaq1=betaq+1.d0
      aaa=-(1.d0-betaq)**2/(2.d0*betaq*betaq)*dlog(betaq1)
     1    +betaq1/(2.d0*betaq)
     2    -2.d0*c5
 1055 aa(i)=aa(i)*aaa
      mi=mi+2
      mm=mm+1
      go to 5999
c
c
c
c
c        function +q^2 (momentum-transfer squared)
c        *************
c
c
 1100 c5=c(mm,im)
      if (c5.eq.0.d0) c5=1.d0
      do 1105 i=1,nt
 1105 aa(i)=-aa(i)*deltaq(i,1)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function k^2 (average-momentum squared)
c        ************
c
c
 1200 c5=c(mm,im)
      if (c5.eq.0.d0) c5=1.d0
      do 1205 i=1,nt
 1205 aa(i)=aa(i)*deltaq(i,3)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function 1 for tpn1
c        *******************
c
c
 1300 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 1305 i=1,nt
      ga4=ga2*ga2
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      brak=4.d0*c4*(5.d0*ga4-4.d0*ga2 -1.d0)
     1    +akk*(23.d0*ga4-10.d0*ga2-1.d0)
     2    +48.d0*ga4*c4*c4/radi
 1305 aa(i)=aa(i)*c5*ell(i)*brak
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function 2 for tpn1
c        *******************
c
c
 1400 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 1405 i=1,nt
 1405 aa(i)=aa(i)*c5*ell(i)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn2, function 1
c        ****************
c
c
 1500 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 1505 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      term1=-ga2*c4**(2.5d0)/(16.d0*radi)
      term2=(2.d0*c4*(2.d0*cb1-cb3)-akk*(cb3+3.d0/16.d0*ga2))
     1     *cpaa(i)
 1505 aa(i)=aa(i)*c5*(term1+term2)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn2, function 2
c        ****************
c
c
 1600 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 1605 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      term1=-3.d0*ga2*c4**(2.5d0)/radi
      term2=(4.d0*c4+2.d0*akk-ga2*(4.d0*c4+3.d0*akk))
     1     *cpaa(i)
 1605 aa(i)=aa(i)*c5*(term1+term2)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn2, function 3
c        ****************
c
c
 1700 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 1705 i=1,nt
 1705 aa(i)=aa(i)*c5*cpaa(i)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn2, function 4
c        ****************
c
c
 1800 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 1805 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      term1=(cb4+0.25d0)*radi
      term2=-ga2/8.d0*(10.d0*c4+3.d0*akk)
 1805 aa(i)=aa(i)*c5*(term1+term2)*cpa(i)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn2, function 5
c        ****************
c
c
 1900 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 1905 i=1,nt
 1905 aa(i)=aa(i)*c5*cpaa(i)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn2, function 6
c        ****************
c
c
 2000 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 2005 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
 2005 aa(i)=aa(i)*c5*radi*cpa(i)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function q^4 (momentum-transfer to the power of 4)
c        ************
c
c
 2100 c5=c(mm,im)
      if (c5.eq.0.d0) c5=1.d0
      do 2105 i=1,nt
 2105 aa(i)=aa(i)*deltaq(i,4)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function k^4 (average-momentum to the power of 4)
c        ************
c
c
 2200 c5=c(mm,im)
      if (c5.eq.0.d0) c5=1.d0
      do 2205 i=1,nt
 2205 aa(i)=aa(i)*deltaq(i,5)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function +q^2*k^2
c        *****************
c
c
 2300 c5=c(mm,im)
      if (c5.eq.0.d0) c5=1.d0
      do 2305 i=1,nt
 2305 aa(i)=aa(i)*deltaq(i,6)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        function  (\vec q x \vec k)^2
c        *****************************
c
c
 2400 c5=c(mm,im)
      if (c5.eq.0.d0) c5=1.d0
      do 2405 i=1,nt
 2405 aa(i)=aa(i)*deltaq(i,7)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c        function xy
c        ***********
c
 2500 c5=c(mm,im)
      if (c5.eq.0.d0) c5=1.d0
      aaxy=xy2*0.5d0*c5
      do 2505 i=1,nt
 2505 aa(i)=aa(i)*aaxy
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c        function xx+yy
c        **************
c
 2600 c5=c(mm,im)
      if (c5.eq.0.d0) c5=1.d0
      do 2605 i=1,nt
 2605 aa(i)=aa(i)*xxpyy*c5
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c        function xx*xx+yy*yy
c        ********************
c
 2700 c5=c(mm,im)
      if (c5.eq.0.d0) c5=1.d0
      aaxy=(xx*xx+yy*yy)*c5
      do 2705 i=1,nt
 2705 aa(i)=aa(i)*aaxy
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c        function xx
c        ***********
c
 2800 c5=c(mm,im)
      if (c5.eq.0.d0) c5=1.d0
      do 2805 i=1,nt
 2805 aa(i)=aa(i)*xx*c5
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c        function yy
c        ***********
c
 2900 c5=c(mm,im)
      if (c5.eq.0.d0) c5=1.d0
      do 2905 i=1,nt
 2905 aa(i)=aa(i)*yy*c5
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn3, function 1
c        ****************
c
c
 3000 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 3005 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      brak=(cb2*radi/6.d0+cb3*(2.d0*c4+akk)-4.d0*cb1*c4)**2
     1    +(cb2*radi)**2/45.d0
 3005 aa(i)=aa(i)*c5*ell(i)*brak
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn3, function 2
c        ****************
c
c
 3100 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 3105 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
 3105 aa(i)=aa(i)*c5*ell(i)*radi
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c        function 1.d0
c        *************
c
 3200 c5=c(mm,im)
      if (c5.eq.0.d0) c5=1.d0
      do 3205 i=1,nt
 3205 aa(i)=aa(i)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c        function 1-q^2/8-k^/2
c        *********************
c
 3300 c5=c(mm,im)
      if (c5.eq.0.d0) c5=1.d0
      do 3305 i=1,nt
 3305 aa(i)=aa(i)*(1.d0+deltaq(i,1)/8.d0-deltaq(i,3)/2.d0)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c        function 1-q^2/8
c        ****************
c
 3400 c5=c(mm,im)
      if (c5.eq.0.d0) c5=1.d0
      do 3405 i=1,nt
 3405 aa(i)=aa(i)*(1.d0+deltaq(i,1)/8.d0)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c        function 1+k^/2
c        ***************
c
 3500 c5=c(mm,im)
      if (c5.eq.0.d0) c5=1.d0
      do 3505 i=1,nt
 3505 aa(i)=aa(i)*(1.d0+deltaq(i,3)/2.d0)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn3, function 3
c        ****************
c
c
 3600 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 3605 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      brak=radi+ga2*(8.d0*c4+5.d0*akk)
 3605 aa(i)=aa(i)*c5*ell(i)*brak
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn3, function 4
c        ****************
c
c
 3700 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 3705 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      brak=radi-ga2*(16.d0*c4+7.d0*akk)
 3705 aa(i)=aa(i)*c5*ell(i)*brak
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn3, function 5
c        ****************
c
c
 3800 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 3805 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
 3805 aa(i)=aa(i)*c5*ell(i)*radi
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn3, function 6
c        ****************
c
c
 3900 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 3905 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      brak=(cb2-6.d0*cb3)*akk*akk
     1     +4.d0*(6.d0*cb1+cb2-3.d0*cb3)*akk*c4
     2     +6.d0*(cb2-2.d0*cb3)*c4*c4
     3     +24.d0*(2.d0*cb1+cb3)*c4*c4*c4/radi
 3905 aa(i)=aa(i)*c5*ell(i)*brak
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn32, function 1
c        *****************
c
c
 4000 continue
c
c     2-pion case      
      if (.not.indmu.or.cc4mu.ne.c4) go to 30
c
      wpi=dsqrt(c4)
c
c     precalculate ImV(mu)
      if (.not.iv40) then
            iv40 = .true.
            do 4035 ii = 1, nmu
                  ANS1 = 2.d0*c4 - amu2(ii)
                  v40(ii) = -ANS1/amu(ii)*
     1         (  (c4 - 2.d0*amu2(ii))*(2.d0*wpi + ANS1/2.d0/amu(ii)*
     2            dlog((amu(ii) + 2.d0*wpi)/(amu(ii) - 2.d0*wpi)))
     3            + 4.d0*ga2*wpi*ANS1    )
 4035       continue
      end if
c
c     iterate over output points aa(i), store results in them
      do 4005 i=1,nt
c           set aakk = q^2 and c5 = overall_const_factor * q^6, see [EM, D.16]
            aakk=-deltaq(i,1)
            c5=c(mm,im)*(aakk**3)
c
c           integrate over mu in [2m_pi, +infinity], and store result in aintmu
            aintmu=0.d0
            do 4015 ii=1,nmu
                  aintmu = aintmu + wmu(ii)*v40(ii) / 
     1                  (amu4(ii)*amu(ii)*(amu2(ii)+aakk))
 4015       continue
c
            aa(i)=aa(i)*c5*aintmu
 4005 continue
c
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn32, function 2
c        *****************
c
c
 4100 continue
c
c     2-pion case      
      if (.not.indmu.or.cc4mu.ne.c4) go to 30
c
      wpi=dsqrt(c4)
c
c     precalculate ImV(mu)
      if (.not.iv41) then
            iv41 = .true.
            do 4135 ii = 1, nmu
                  v41(ii) = (amu2(ii) - 4.d0*c4)*((c4 - amu2(ii)/4.d0)
     1            *dlog((amu(ii) + 2.d0*wpi)/(amu(ii) - 2.d0*wpi))
     2            + (1.d0 + 2.d0*ga2)*amu(ii)*wpi)
 4135       continue
      end if
c
c     iterate over output points aa(i), store results in them
      do 4105 i=1,nt
c           set aakk = q^2 and c5 = overall_const_factor * q^4, see [EM, D.17]
            aakk=-deltaq(i,1)
            c5=c(mm,im)*(aakk**2)
c
c           integrate over mu in [2m_pi, +infinity], and store result in aintmu
            aintmu=0.d0
            do 4115 ii=1,nmu
                  aintmu = aintmu + wmu(ii)*v41(ii) / 
     1                  (amu4(ii)*amu(ii)*(amu2(ii)+aakk))
 4115       continue
c
            aa(i)=aa(i)*c5*aintmu
 4105 continue
c
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn32, function 3
c        *****************
c
c
 4200 continue
c
c     2-pion case      
      if (.not.indmu.or.cc4mu.ne.c4) go to 30
c
      wpi=dsqrt(c4)
c
c     precalculate ImV(mu)
      if (.not.iv42) then
            iv42 = .true.
            do 4235 ii = 1, nmu
                  rootk = dsqrt(amu2(ii)/4.d0 - c4)
                  v42(ii) = 2.d0/3.d0 * amu(ii)*(rootk**3)
 4235       continue
      end if
c
c     iterate over output points aa(i), store results in them
      do 4205 i=1,nt
c           set aakk = q^2 and c5 = overall_const_factor * q^4, see [EM, D.17]
            aakk=-deltaq(i,1)
            c5=c(mm,im)*(aakk**2)
c
c           integrate over mu in [2m_pi, +infinity], and store result in aintmu
            aintmu=0.d0
            do 4215 ii=1,nmu
                  aintmu = aintmu + wmu(ii)*v42(ii) / 
     1                  (amu4(ii)*amu(ii)*(amu2(ii)+aakk))
 4215       continue
c
            aa(i)=aa(i)*c5*aintmu
 4205 continue
c
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn32, function 4
c        *****************
c
c
 4300 continue
c
c     2-pion case      
      if (.not.indmu.or.cc4mu.ne.c4) go to 30
c
      wpi=dsqrt(c4)
c
c     precalculate ImV(mu)
      if (.not.iv43) then
            iv43 = .true.
            do 4335 ii = 1, nmu
                  rootk = dsqrt(amu2(ii)/4.d0 - c4)
                  v43(ii) = rootk/amu(ii)*
     1  (      (  (4.d0*c4*(1.d0 + 2.d0*ga2) 
     2                  - amu2(ii)*(1.d0 + 5.d0*ga2)) * rootk/amu(ii)*
     3            dlog((amu(ii) + 2.d0*rootk)/2.d0/wpi)
     4            + amu2(ii)*(5.d0 + 13.d0*ga2)/12.d0
     5            - 2.d0*c4*(1.d0 + 2.d0*ga2)
     6            + 96.d0*pi2*fpi2 * ( 
     7                  (2.d0*c4 - amu2(ii))*cd12 
     9                  + 4.d0*c4*cd5)  )
     1     *(ga2*(2.d0*c4-amu2(ii))+2.d0/3.d0*(ga2-1.d0)*rootk*rootk)
     2  -192.d0*pi2*fpi2*cd3*rootk*rootk*(ga2/3.d0*(2.d0*c4-amu2(ii))
     3                    +2.d0/5.d0*(ga2-1.d0)*rootk*rootk)        )
c
 4335       continue
      end if
c
c     iterate over output points aa(i), store results in them
      do 4305 i=1,nt
c           set aakk = q^2 and c5 = overall_const_factor * q^6, see [EM, D.16]
            aakk=-deltaq(i,1)
            c5=c(mm,im)*(aakk**3)
c
c           integrate over mu in [2m_pi, +infinity], and store result in aintmu
            aintmu=0.d0
            do 4315 ii=1,nmu
                  aintmu = aintmu + wmu(ii)*v43(ii) / 
     1                  (amu4(ii)*amu(ii)*(amu2(ii)+aakk))
 4315       continue
c
            aa(i)=aa(i)*c5*aintmu
 4305 continue
c 
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn32, function 5
c        *****************
c
c
 4400 continue
c
c     2-pion case      
      if (.not.indmu.or.cc4mu.ne.c4) go to 30
c
      wpi=dsqrt(c4)
      wpi2=2.d0*wpi
c
c     precalculate ImV(mu)
      if (.not.iv44) then
            iv44 = .true.
c
            nxi=ic(mi+1,im)
            if (nxi.ne.nnxi) then
                  nnxi=nxi
                  call gset(0.d0,1.d0,nxi,xi,wxi)
            end if
c
            do 4435 ii = 1, nmu
                  akk=-amu2(ii)
                  radi=4.d0*c4+akk
                  root15=(-radi)**(1.5d0)
c
                  wbarq=1.d0+akk/(4.d0*c4)
                  wbar=dsqrt(-wbarq)
c
                  aintxia=0.d0
                  aintxib=0.d0
                  do 4425 iii=1,nxi
                        xixi=xi(iii)*xi(iii)
                        term1=wbarq*xixi
                        term2=1.d0/term1
                        cterm1=(1.d0-term2)**(1.5d0)
                        cterm2=dsqrt(1.d0-term1)
                        cterm3=dlog(wbar*xi(iii)+cterm2)
                        term4=cterm1*cterm3
                        term2=-term2
                        term4=-term4
                        term24=term2+term4
c
                        term0=-1.d0/6.d0
                        step44a=wxi(iii)*(1.d0-xixi)
     1                        *(term0       )
                        step44b=wxi(iii)*(1.d0-xixi)
     1                        *(      term24)
c
                        aintxia=aintxia+step44a
                        aintxib=aintxib+step44b
 4425             continue
c
                  v44a(ii) = aintxia*root15
                  v44b(ii) = aintxib*root15
 4435       continue
      end if
c
      do 4405 i=1,nt
            aakk=-deltaq(i,1)
            c5=c(mm,im)*(aakk**2)/4.d0
c
            aintmua=0.d0
            aintmub=0.d0
            do 4415 ii=1,nmu
                  aintmua=aintmua
     1                  +v44a(ii)*wmu(ii)
     2                  /(amu4(ii)*(amu2(ii)+aakk))
                  aintmub=aintmub
     1                  +v44b(ii)*wmu(ii)
     2                  /(amu4(ii)*(amu2(ii)+aakk))
 4415       continue
c
            aintmu=aintmua+aintmub
            aa(i)=aa(i)*c5*aintmu
 4405 continue
c
      mi=mi+2
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn32, function 6
c        *****************
c
c
 4500 continue
c
c     2-pion case
      if (.not.indmu.or.cc4mu.ne.c4) go to 30
      wpi=dsqrt(c4)
      wpi2=2.d0*wpi
c
c     precalculate ImV(mu)
      if (.not.iv45) then
            iv45 = .true.
c
            nxi=ic(mi+1,im)
            if (nxi.ne.nnxi) then
                  nnxi=nxi
                  call gset(0.d0,1.d0,nxi,xi,wxi)
            end if
c
            do 4535 ii = 1, nmu
                  akk=-amu2(ii)
                  radi=4.d0*c4+akk
                  root=dsqrt(-radi)
c
                  wbarq=1.d0+akk/(4.d0*c4)
                  wbar=dsqrt(-wbarq)
                  brak1=ga2*(2.d0*c4+akk)
                  brak2=2.d0*(ga2-1.d0)*c4*wbarq
                  brak3=3.d0*c4*wbarq
                  brak4=6.d0*c4*wbar
c
c
                  aintxi=0.d0
                  do 4525 iii=1,nxi
                        xixi=xi(iii)*xi(iii)
                        term1=wbarq*xixi
                        term2=1.d0/term1
                        cterm1=(1.d0-term2)**(1.5d0)
                        cterm2=dsqrt(1.d0-term1)
                        cterm3=dlog(wbar*xi(iii)+cterm2)
                        term4=cterm1*cterm3
                        term24=term2+term4
c
                        term5=cterm2*cterm3*brak4*xi(iii)
                        brak5=ga2*ga2*
     1                        (akk-2.d0*c4*wbarq*xixi+2.d0*c4)
c
c
                        term6=5.d0/6.d0-term24
                        step45=wxi(iii)*(brak1-brak2*xixi)
     1                        *(brak3*xixi+term5-brak5*term6)
c
                        aintxi=aintxi+step45
 4525             continue
                  v45(ii) = aintxi*root
 4535       continue
      end if
c
      do 4505 i=1,nt
            aakk=-deltaq(i,1)
            c5=c(mm,im)*(aakk**3)
c
            aintmu=0.d0
            do 4515 ii=1,nmu
                  aintmu=aintmu
     1                  +v45(ii)*wmu(ii)
     2                  /(amu6(ii)*(amu2(ii)+aakk))
 4515       continue
c
            aa(i)=aa(i)*c5*aintmu
 4505 continue
c
      mi=mi+2
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn4, function 1
c        ****************
c           with analytical formula for ImVc(mu)      
c
 4600 continue
c     prepare mu integration. 2-pion case 
      if (.not.indmu.or.cc4mu.ne.c4) go to 30
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi5 = wpi*c4*c4
c
c     precalculate ImV(mu)
      if (.not.iv46) then
            iv46 = .true.
            do 4635 ii = 1, nmu
c                 unitless u = mu/mpi
                  u = amu(ii)/wpi
                  u2 = u*u
                  u3 = u2*u
                  u4 = u2*u2
                  uroot = dsqrt(u2 - 4.d0)
                  alogfn = dlog((u + uroot)/2.d0)
c
                  ANS1 = ga2*uroot*(5.d0 - 2.d0*u2 - 2.d0/u2)*
     1            (24.d0*cb1 + cb2*(u2 - 4.d0) + 6.d0*cb3*(u2 - 2.d0))
     2            *dlog((u + 2.d0)/(u - 2.d0))
     3            + 8.d0/u*
     4         (  3.d0*(4.d0*cb1 + cb3*(u2 - 2.d0))*
     5            (1.d0 - 10.d0*ga4 + 4.d0*ga4*u2) 
     6            + cb2*(6.d0*ga4*u2 - 10.d0*ga4 - 3.d0)  )*alogfn
c
                  ANS2 = 3.d0*(2.d0 - u2)*(4.d0*cb1 + cb3*(u2 - 2.d0))
     1            + cb2*(7.d0*u2 - 6.d0 - u4)
     2            + 4.d0*ga2/u*(2.d0*u2 - 1.d0)*
     3            (4.d0*(6.d0*cb1 - cb2 - 3.d0*cb3) 
     4                  + (cb2 + 6.d0*cb3)*u2)
c                  
                  ANS3 = ga4*
     1         (  128.d0/(2.d0+u)*(2.d0*cb1 + cb3) + 256.d0/u*
     2            (2.d0*cb1 + cb2/3.d0 - cb3) 
     3                  + 4.d0*(14.d0*cb3 - 5.d0*cb2 - 92.d0*cb1)
     4            + 32.d0*u*(6.d0*cb3 - 5.d0*cb2/3.d0)
     5            + 2.d0*u2/3.d0*(36.d0*cb1 + 13.d0*cb2 - 156.d0*cb3)
     6            + 2.d0*u4/3.d0*(2.d0*cb2 + 9.d0*cb3)               )
c
                  v46(ii) = ANS1 + uroot*(ANS2+ANS3)
 4635       continue
      end if
c
c     iterate over output points aa(i), store results in them
      do 4605 i=1,nt
            aakk=-deltaq(i,1)
            c5=c(mm,im)*wpi5*(aakk**3)
c
c           integrate over mu in [2m_pi, +infinity], and store result in aintmu
            aintmu=0.d0
            do 4615 ii=1,nmu
                  aintmu = aintmu + wmu(ii)*v46(ii) / 
     1                  (amu4(ii)*amu(ii)*(amu2(ii)+aakk))
 4615       continue
c
c
            aa(i)=aa(i)*c5*aintmu
 4605 continue
      mi=mi+2
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn4, function 1
c        ****************
c           with analytical formula for ImVt(mu). Calculates Vt. For Vs multilpy by -q^2 through input data file.
c
 4700 continue
c     prepare mu integration. 2-pion case 
      if (.not.indmu.or.cc4mu.ne.c4) go to 30
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi5 = wpi*c4*c4      
c
c     precalculate ImV(mu)
      if (.not.iv47) then
            iv47 = .true.
            do 4735 ii = 1, nmu
c                 unitless u = mu/mpi
                  u = amu(ii)/wpi
                  u2 = u*u
                  u2m4 = u2 - 4.d0
                  u3 = u2*u
                  u4 = u2*u2
                  uroot = dsqrt(u2m4)
                  alogfn = dlog((u + uroot)/2.d0)
c
                  v47(ii) = 8.d0*ga2*u*(u2 - 5.d0)*alogfn 
     1            - 1.d0/3.d0*u2m4*u2m4*uroot*
     2                        dlog((u + 2.d0)/(u - 2.d0))
     3            + u/3.d0*uroot*(ga2*(u3 - 30.d0*u + 64.d0)
     4                  + 4.d0*u2 - 16.d0)
 4735       continue
      end if
c
c     iterate over output points aa(i), store results in them
      do 4705 i=1,nt
            aakk=-deltaq(i,1)
            c5=c(mm,im)*wpi5*(aakk*aakk)
c
c           integrate over mu in [2m_pi, +infinity], and store result in aintmu
            aintmu = 0.d0
            do 4715 ii = 1, nmu
                  aintmu = aintmu + wmu(ii)*v47(ii) / 
     1                  (amu4(ii)*amu(ii)*(amu2(ii)+aakk))
 4715       continue
c
            aa(i)=aa(i)*c5*aintmu
 4705 continue
      mi=mi+2
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn3m (= N^3LO, 1/M^2 terms), function 1
c        ****************************************
c
c
 4800 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      c44=c4*c4
      c46=c44*c4
      c48=c46*c4
      do 4805 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      radiq=radi*radi
      brak1=c46/(2.d0*radi)
      brak2=(2.d0*c48/radiq+8.d0*c46/radi-akk*akk-2.d0*c44)*ell(i)
 4805 aa(i)=aa(i)*c5*(brak1+brak2)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn3m, function 2
c        *****************
c
c
 4900 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      c44=c4*c4
      c46=c44*c4
      c48=c46*c4
      yyoff=0.5d0*(xx+yy)
      do 4905 i=1,nt
      akk=-deltaq(i,1)
      akkq=akk*akk
      radi=4.d0*c4+akk
      radiq=radi*radi
      brak1=(radi*(akk-4.d0*yyoff)
     1      +8.d0*ga2*(11.d0/4.d0*akkq+5.d0*c4*akk+3.d0*c44
     2       -6.d0*c46/radi-yyoff*(8.d0*c4+5.d0*akk))
     3      +4.d0*ga4*(yyoff*(20.d0*c4+7.d0*akk-16.d0*c44/radi)
     5       +16.d0*c48/radiq+12.d0*c46/radi-27.d0/4.d0*akkq
     6       -11.d0*c4*akk-6.d0*c44))*ell(i)
      brak2=16.d0*ga4*c46/radi
 4905 aa(i)=aa(i)*c5*(brak1+brak2)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn3m, function 3
c        *****************
c
c
 5000 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      c44=c4*c4
      yyoff=0.5d0*(xx+yy)
      do 5005 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      brak=yyoff+3.d0/8.d0*akk+c44/radi
 5005 aa(i)=aa(i)*c5*brak*ell(i)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn3m, function 4
c        *****************
c
c
 5100 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      c44=c4*c4
      do 5105 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      brak=radi-32.d0*ga2*(c4+7.d0/16.d0*akk)
     1     +4.d0*ga4*(7.d0*c4+17.d0/4.d0*akk+4.d0*c44/radi)
 5105 aa(i)=aa(i)*c5*brak*ell(i)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn3m, function 5
c        *****************
c
c
 5200 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      c44=c4*c4
      do 5205 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      brak=-radi+16.d0*ga2*(c4+3.d0/8.d0*akk)
     1     +4.d0/3.d0*ga4*(-9.d0*c4-11.d0/4.d0*akk+4.d0*c44/radi)
 5205 aa(i)=aa(i)*c5*brak*ell(i)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn3m, function 6
c        *****************
c
c
 5300 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      c44=c4*c4
      do 5305 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      brak=11.d0/32.d0*akk+c44/radi
 5305 aa(i)=aa(i)*c5*brak*ell(i)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn3m, function 7
c        *****************
c
c
 5400 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 5405 i=1,nt
 5405 aa(i)=aa(i)*c5*ell(i)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn2c (correction for our it 2pi), function 1
c        *********************************************
c
c
 5500 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 5505 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      term1=dsqrt(c4)*radi
      term2=(2.d0*c4+akk)*cpaa(i)
 5505 aa(i)=aa(i)*c5*(term1+term2)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn2c (correction for our it 2pi), function 2
c        *********************************************
c
c
 5600 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 5605 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      term1=dsqrt(c4)
      term2=radi*cpa(i)
 5605 aa(i)=aa(i)*c5*(term1+term2)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn2, function 1, the non-1/M part of 1500
c        ******************************************
c
c
 5700 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 5705 i=1,nt
      akk=-deltaq(i,1)
      term2=(2.d0*c4*(2.d0*cb1-cb3)-akk*(cb3))
     1     *cpaa(i)
 5705 aa(i)=aa(i)*c5*(term2)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn2m, function 1, the 1/M part of 1500
c        ***************************************
c
c
 5800 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 5805 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      term1=-ga2*c4**(2.5d0)/(16.d0*radi)
      term2=(-akk*(3.d0/16.d0*ga2))
     1     *cpaa(i)
 5805 aa(i)=aa(i)*c5*(term1+term2)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn2, function 4, the non-1/M part of 1800
c        ******************************************
c
c
25900 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 25905 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      term1=(cb4)*radi
25905 aa(i)=aa(i)*c5*(term1)*cpa(i)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c        tpn2m, function 4, the 1/M part of 1800
c        ***************************************
c
c
26000 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 26005 i=1,nt
      akk=-deltaq(i,1)
      radi=4.d0*c4+akk
      term1=(0.25d0)*radi
      term2=-ga2/8.d0*(10.d0*c4+3.d0*akk)
26005 aa(i)=aa(i)*c5*(term1+term2)*cpa(i)
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c
c
c
c           <<-- Vc4 BEGIN
c           **************
c
c
26100 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
c      
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi4 = wpi*wpi
c      
      do 26105 i=1,nt
c            q^2
            akk=-deltaq(i,1)
c            
c           omega^2            
            radi=4.d0*c4+akk
c
c           ANS = {curly bracket}
            ANS = 2.d0*wpi4*wpi/radi 
     1      - 3.d0*(4.d0*wpi4 - akk*akk)*cpa(i)
c
            aa(i)=aa(i)*c5*ANS
26105 continue
      mi=mi+1
      mm=mm+1
      go to 5999
c
c           Vc4 END -->>
c           ************
c
c
c
c
c           <<-- Wc4 BEGIN
c           **************
c
26200 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
c      
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi4 = wpi*wpi
c      
      do 26205 i=1,nt
c            q^2
            akk=-deltaq(i,1)
c            
c           omega^2            
            radi=4.d0*c4+akk
c
c           ANS = {curly bracket}
            ANS = 3.d0*ga2*wpi4*wpi/radi
     1      - (4.d0*c4 + 2.d0*akk - ga2*(7.d0*c4 + 4.5d0*akk))*cpaa(i)
c
            aa(i)=aa(i)*c5*ANS
26205 continue
      mi=mi+1
      mm=mm+1
      go to 5999
c
c           Wc4 END -->>
c           ************
c
c
c
c           <<-- Vt4 BEGIN
c           **************
c
c
26300 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
c      
      do 26305 i=1,nt
c            q^2
            akk=-deltaq(i,1)
c
c           here, pion mass is: c4 = m_pi^2
            ANS = (4.d0*c4 + 1.5d0*akk)*cpa(i)
c                  
            aa(i)=aa(i)*c5*ANS
26305 continue
      mi=mi+1
      mm=mm+1
      go to 5999
c
c           Vt4 END -->>
c           ************
c
c
c
c           <<-- Wt4 BEGIN
c           **************
c
26400 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
c      
      do 26405 i=1,nt
c            q^2
            akk=-deltaq(i,1)
c
c           here, pion mass is: c4 = m_pi^2
            ANS = (8.d0*(1.d0 - ga2)*c4 + (2.d0 - 2.5d0*ga2)*akk)
     1            * cpa(i)
c                  
            aa(i)=aa(i)*c5*ANS
26405 continue
      mi=mi+1
      mm=mm+1
      go to 5999
c
c           Wt4 END -->>
c           ************
c
c
c
c
c           <<-- V_LS4 BEGIN
c           ****************
c
26500 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
c      
      do 26505 i=1,nt
            aa(i)=aa(i)*c5* cpaa(i)
26505 continue
      mi=mi+1
      mm=mm+1
      go to 5999
c
c           V_LS4 END -->>
c           **************
c
c
c
c           <<-- W_LS4 BEGIN
c           ****************
c
26600 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
c      
      do 26605 i=1,nt
c            q^2
            akk=-deltaq(i,1)
c
c           here, pion mass is: c4 = m_pi^2
            ANS = (4.d0*c4 + akk)*cpa(i)
c                  
            aa(i)=aa(i)*c5*ANS
26605 continue
      mi=mi+1
      mm=mm+1
      go to 5999
c
c           W_LS4 END -->>
c           **************
c
c
c
c
c           <<-- Wc13 BEGIN
c           ***************
c
c
26700 continue
c     prepare mu integration. 3-pion case 
      if (.not.i3ndmu.or.c3c4mu.ne.c4) go to 40
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi5 = wpi*c4*c4      
c
c     precalculate ImV(mu)
      if (.not.iv67) then
            iv67 = .true.
            do 26735 ii = 1, n3mu
c                 unitless u = mu/mpi
                  u = a3mu(ii)/wpi
                  u2 = u*u
                  u3 = u2*u
                  u4 = u2*u2
c                 uroot = wurz; uln = logneu                  
                  uroot = dsqrt(u2 - 2.d0*u - 3.d0)
                  uln = dlog((u - 1.d0 + uroot)/2.d0)
c
                  v67(ii) = 8.d0/3.d0*uroot*(u - 1.d0)*
     1            (u - 4.d0 - 2.d0*u2 - u3)
     2            + 32.d0*uln*(u3 - 4.d0*u + 1.d0/u)
26735       continue
      end if
c
c     iterate over output points aa(i), store results in them
      do 26705 i=1,nt
            aakk=-deltaq(i,1)
            c5=c(mm,im)*wpi5*(aakk**3)
c
c
c           integrate over mu in [3m_pi, +infinity], and store result in aintmu
            aintmu = 0.d0
            do 26715 ii = 1, n3mu
                  aintmu = aintmu + w3mu(ii)*v67(ii) / 
     1                  (a3mu4(ii)*a3mu(ii)*(a3mu2(ii)+aakk))
26715       continue
c
            aa(i)=aa(i)*c5*aintmu
26705 continue
c
      mi=mi+2
      mm=mm+1
      go to 5999
c
c           Wc13 END -->>
c           *************
c
c
c
c           <<-- Vs13 BEGIN
c           ***************
c
26800 continue
c     prepare mu integration. 3-pion case 
      if (.not.i3ndmu.or.c3c4mu.ne.c4) go to 40
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi5 = wpi*c4*c4      
c
c     precalculate ImV(mu)
      if (.not.iv68) then
            iv68 = .true.
            do 26835 ii = 1, n3mu
c                 unitless u = mu/mpi
                  u = a3mu(ii)/wpi
                  u2 = u*u
                  u3 = u2*u
                  u4 = u2*u2
                  u5 = u3*u2
                  u6 = u3*u3
c                 uroot = wurz; uln = logneu                  
                  uroot = dsqrt(u2 - 2.d0*u - 3.d0)
                  uln = dlog((u - 1.d0 + uroot)/2.d0)
c
                  ANS = uroot*(u - 1.d0)/24.d0*
     1            (37.d0*u6 + 74.d0*u5 - 251.d0*u4 - 268.d0*u3 
     2                  + 349.d0*u2 - 58.d0*u - 135.d0) 
     3            + 2.d0*uln*(39.d0*u4 - 2.d0 - 52.d0*u2 - 6.d0*u6)
c
                  v68(ii) = ANS/u3
26835       continue
      end if
c      
c     iterate over output points aa(i), store results in them
      do 26805 i=1,nt
            aakk=-deltaq(i,1)
            c5=c(mm,im)*wpi5*aakk*aakk
c
c           integrate over mu in [3m_pi, +infinity], and store result in aintmu
            aintmu = 0.d0
            do 26815 ii = 1, n3mu
                  aintmu = aintmu + w3mu(ii)*v68(ii) / 
     1                  (a3mu4(ii)*a3mu(ii)*(a3mu2(ii)+aakk))
26815       continue
c
            aa(i)=aa(i)*c5*aintmu
26805 continue
c
      mi=mi+2
      mm=mm+1
      go to 5999
c
c           Vs13 END -->>
c           *************
c
c
c
c           <<-- Vt13 BEGIN
c           ***************
c
26900 continue
c     prepare mu integration. 3-pion case 
      if (.not.i3ndmu.or.c3c4mu.ne.c4) go to 40
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi3 = wpi*c4
c
c     precalculate ImV(mu)
      if (.not.iv69) then
            iv69 = .true.
            do 26935 ii = 1, n3mu
c                 unitless u = mu/mpi
                  u = a3mu(ii)/wpi
                  u2 = u*u
                  u3 = u2*u
                  u4 = u2*u2
                  u5 = u3*u2
                  u6 = u3*u3
c                 uroot = wurz; uln = logneu                  
                  uroot = dsqrt(u2 - 2.d0*u - 3.d0)
                  uln = dlog((u - 1.d0 + uroot)/2.d0)
c
                  ANS = uroot*(u - 1.d0)/12.d0*
     1            (5.d0*u6 + 10.d0*u5 - 3.d0*u4 - 252.d0*u3 - 443.d0*u2
     2                  - 58.d0*u - 135.d0)
     3            + 4.d0*uln*(3.d0*u4 + 22.d0*u2 - 2.d0)
c
                  v69(ii) = ANS / u5
26935       continue
      end if
c      
c     iterate over output points aa(i), store results in them
      do 26905 i=1,nt
            aakk=-deltaq(i,1)
            c5=c(mm,im)*wpi3*aakk*aakk
c
c           integrate over mu in [3m_pi, +infinity], and store result in aintmu
            aintmu = 0.d0
            do 26915 ii = 1, n3mu
                  aintmu = aintmu + w3mu(ii)*v69(ii) / 
     1                  (a3mu2(ii)*a3mu(ii)*(a3mu2(ii)+aakk))
26915       continue
c
            aa(i)=aa(i)*c5*aintmu
26905 continue
c
      mi=mi+2
      mm=mm+1
      go to 5999
c
c
c           Vt13 END -->>
c           *************
c
c
c
c           <<-- Ws13 BEGIN
c           ***************
c
27000 continue
c     prepare mu integration. 3-pion case 
      if (.not.i3ndmu.or.c3c4mu.ne.c4) go to 40
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi5 = wpi*c4*c4      
c
c     precalculate ImV(mu)
      if (.not.iv70) then
            iv70 = .true.
            do 27035 ii = 1, n3mu
c                 unitless u = mu/mpi
                  u = a3mu(ii)/wpi
                  u2 = u*u
                  u3 = u2*u
                  u4 = u2*u2
                  u5 = u3*u2
                  u6 = u3*u3
c                 uroot = wurz; uln = logneu                  
                  uroot = dsqrt(u2 - 2.d0*u - 3.d0)
                  uln = dlog((u - 1.d0 + uroot)/2.d0)
c
                  ANS1 = uroot*(u - 1.d0)*
     1         (  2.d0*cb1*u*(5.d0*u3 + 10.d0*u2 - 5.d0*u - 4.d0)
     2            + cb2/48.d0*(135.d0 + 58.d0*u - 277.d0*u2 - 36.d0*u3
     3                  + 147.d0*u4 - 10.d0*u5 - 5.d0*u6)
     4            + cb3/8.d0*(7.d0*u6 + 14.d0*u5 -145.d0*u4 - 20.d0*u3
     5                  + 111.d0*u2 + 18.d0*u + 27.d0)
     6            + cb4/6.d0*(44.d0*u3 + 37.d0*u4 - 14.d0*u5 - 7.d0*u6
     7                  - 3.d0*u2 - 18.d0*u - 27.d0)                  )
c     
                  ANS2 = uln*(24.d0*cb1*(1.d0 + 4.d0*u2 - 3.d0*u4)
     1            + cb2*(2.d0 + 2.d0*u2 - 3.d0*u4)
     2            + 6.d0*cb3*u2*(3.d0*u2 - 2.d0)
     3            + 8.d0*cb4*u2*(u4 - 5.d0*u2 + 5.d0) )
c
                  v70(ii) = (ANS1 + ANS2)/u3
27035       continue
      end if
c      
c     iterate over output points aa(i), store results in them
      do 27005 i=1,nt
            aakk=-deltaq(i,1)
            c5=c(mm,im)*wpi5*aakk*aakk
c
c
c           integrate over mu in [3m_pi, +infinity], and store result in aintmu
            aintmu = 0.d0
            do 27015 ii = 1, n3mu
                  aintmu = aintmu + w3mu(ii)*v70(ii) / 
     1                  (a3mu4(ii)*a3mu(ii)*(a3mu2(ii)+aakk))
27015       continue
c
            aa(i)=aa(i)*c5*aintmu
27005 continue
c
      mi=mi+2
      mm=mm+1
      go to 5999
c
c           Ws13 END -->>
c           *************
c
c
c
c           <<-- Wt13 BEGIN
c           ***************
c
c*************** new 27100 3/14/2015
c
27100 continue
c
c     prepare mu integration. 3-pion case 
      if (.not.i3ndmu.or.c3c4mu.ne.c4) go to 40
c
      wpi=dsqrt(c4)
      wpi3=c4*wpi
      fpi=dsqrt(fpi2)
c      
c
      if (.not.iv71) then
      iv71 = .true.
c
      c5=-2.d0*ga4*wpi3/((4.d0*fpi)**6*pi2*pi)
c
      do 27135 ii = 1, n3mu
c
      u=a3mu(ii)/wpi
      u2=u*u
      u3=u2*u
      u4=u2*u2
      u5=u4*u
      u6=u5*u
      uroot = dsqrt(u2 - 2.d0*u - 3.d0)
      uln = dlog((u - 1.d0 + uroot)/2.d0)
c
      term1=4.d0*cb1*u*(5.d0*u3+10.d0*u2+7.d0*u-4.d0)
      term2=cb2/24.d0*(135.d0+58.d0*u+227.d0*u2+204.d0*u3
     1      +27.d0*u4-10.d0*u5-5.d0*u6)
      term3=cb3/4.d0*(27.d0+18.d0*u-9.d0*u2-68.d0*u3
     1      -121.d0*u4+14.d0*u5+7.d0*u6)
      term4=cb4*(4.d0*u3+19.d0*u4-2.d0*u5-u6-9.d0*u2-6.d0*u-9.d0)
c
      brak1=uroot*(u-1.d0)*(term1+term2+term3+term4)
c
      brak2=2.d0*uln*(24.d0*cb1*(1.d0-3.d0*u4)
     1      +cb2*(2.d0-10.d0*u2-3.d0*u4)
     2      +6.d0*cb3*u2*(3.d0*u2+2.d0)-8.d0*cb4*u4)
c
      v71(ii) = c5/u5*(brak1+brak2)
27135 continue
      end if
c
c     iterate over output points aa(i), store results in them
      do 27105 i=1,nt
c           aakk = q^2 
            aakk=-deltaq(i,1)
c
c           integrate over mu in [2m_pi, +infinity], and store result in aintmu
            aintmu = 0.d0
            do 27115 ii = 1, n3mu
                  aintmu = aintmu + w3mu(ii)*v71(ii) / 
     1                  (a3mu2(ii)*a3mu(ii)*(a3mu2(ii)+aakk))
27115       continue
c
            aa(i)=aa(i)*(aakk**2)*aintmu
27105 continue
c
      mi=mi+2
      mm=mm+1
      go to 5999
c
c           Wt13 END -->>
c           *************
c
c
c
c
c           <<-- Vs14 BEGIN
c           ***************
c
27200 continue
c     prepare mu integration. 3-pion case 
      if (.not.i3ndmu.or.c3c4mu.ne.c4) go to 40
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi5 = wpi*c4*c4      
c
c     precalculate ImV(mu)
      if (.not.iv72) then
            iv72 = .true.
            do 27235 ii = 1, n3mu
c                 unitless u = mu/mpi
                  u = a3mu(ii)/wpi
                  u2 = u*u
                  u3 = u2*u
                  u4 = u2*u2
                  u5 = u3*u2
                  u6 = u3*u3
c                 uroot = wurz; uln = logneu                  
                  uroot = dsqrt(u2 - 2.d0*u - 3.d0)
                  uln = dlog((u - 1.d0 + uroot)/2.d0)
c
                  ANS = uroot*(u - 1.d0)/24.d0*(637.d0*u2 - 58.d0*u
     1            - 135.d0 + 116.d0*u3 - 491.d0*u4 - 22.d0*u5
     2                                                - 11.d0*u6)
     3            + 2.d0*uln*(6.d0*u6 - 9.d0*u4 + 8.d0*u2 - 2.d0)
c
                  v72(ii) = ANS / u3
27235       continue
      end if
c
c     iterate over output points aa(i), store results in them
      do 27205 i=1,nt
            aakk=-deltaq(i,1)
            c5=c(mm,im)*wpi5*(aakk**3)
c
c           integrate over mu in [3m_pi, +infinity], and store result in aintmu
            aintmu = 0.d0
            do 27215 ii = 1, n3mu
                  aintmu = aintmu + w3mu(ii)*v72(ii) / 
     1                  (a3mu4(ii)*a3mu(ii)*(a3mu2(ii)+aakk))
27215       continue
c
            aa(i)=aa(i)*c5*aintmu
27205 continue
c
      mi=mi+2
      mm=mm+1
      go to 5999
c
c           Vs14 END -->>
c           *************
c
c
c
c
c           <<-- Vs12 BEGIN
c           ***************
c
27300 continue
c     prepare mu integration. 3-pion case 
      if (.not.i3ndmu.or.c3c4mu.ne.c4) go to 40
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi5 = wpi*c4*c4      
c
c     precalculate ImV(mu)
      if (.not.iv73) then
            iv73 = .true.
            do 27335 ii = 1, n3mu
c                 unitless u = mu/mpi
                  u = a3mu(ii)/wpi
                  u2 = u*u
                  u3 = u2*u
                  u4 = u2*u2
                  u5 = u3*u2
                  u6 = u3*u3
c                 uroot = wurz; uln = logneu                  
                  uroot = dsqrt(u2 - 2.d0*u - 3.d0)
                  uln = dlog((u - 1.d0 + uroot)/2.d0)
c
                  ANS = uroot/12.d0*(27.d0 + 23.d0*u + 101.d0*u2
     1            - 251.d0*u3 - 85.d0*u4 + 199.d0*u5 - 7.d0*u6 
     2                                                - 7.d0*u6*u)
     3            + 4.d0*uln*(2.d0 + 10.d0*u2 - 9.d0*u4)
c
                  v73(ii) = ANS / u3
27335       continue
      end if
c
c     iterate over output points aa(i), store results in them
      do 27305 i=1,nt
            aakk=-deltaq(i,1)
            c5=c(mm,im)*wpi5*aakk*aakk
c
c           integrate over mu in [3m_pi, +infinity], and store result in aintmu
            aintmu = 0.d0
            do 27315 ii = 1, n3mu
                  aintmu = aintmu + w3mu(ii)*v73(ii) / 
     1                  (a3mu4(ii)*a3mu(ii)*(a3mu2(ii)+aakk))
27315       continue
c
            aa(i)=aa(i)*c5*aintmu
27305 continue
c
      mi=mi+2
      mm=mm+1
      go to 5999
c
c           Vs12 END -->>
c           *************
c
c
c
c
c           <<-- Vt12 BEGIN
c           ***************
c
27400 continue
c     prepare mu integration. 3-pion case 
      if (.not.i3ndmu.or.c3c4mu.ne.c4) go to 40
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi3 = wpi*c4
c
c     precalculate ImV(mu)
      if (.not.iv74) then
            iv74 = .true.
            do 27435 ii = 1, n3mu
c                 unitless u = mu/mpi
                  u = a3mu(ii)/wpi
                  u2 = u*u
                  u3 = u2*u
                  u4 = u2*u2
                  u5 = u3*u2
                  u6 = u3*u3
c                 uroot = wurz; uln = logneu                  
                  uroot = dsqrt(u2 - 2.d0*u - 3.d0)
                  uln = dlog((u - 1.d0 + uroot)/2.d0)
c
                  ANS = uroot*(u - 1.d0)/6.d0*(u6 + 2.d0*u5 - 39.d0*u4
     1            - 12.d0*u3 + 65.d0*u2 - 50.d0*u - 27.d0)
     2            + 8.d0*uln*(3.d0*u4 - 10.d0*u2 + 2.d0)
c
                  v74(ii) = ANS / u5
27435       continue
      end if
c
c     iterate over output points aa(i), store results in them
      do 27405 i=1,nt
            aakk=-deltaq(i,1)
            c5=c(mm,im)*wpi3*aakk*aakk
c
c           integrate over mu in [3m_pi, +infinity], and store result in aintmu
            aintmu = 0.d0
            do 27415 ii = 1, n3mu
                  aintmu = aintmu + w3mu(ii)*v74(ii) / 
     1                  (a3mu2(ii)*a3mu(ii)*(a3mu2(ii)+aakk))
27415       continue
c
            aa(i)=aa(i)*c5*aintmu
27405 continue
c
      mi=mi+2
      mm=mm+1
      go to 5999
c
c           Vt12 END -->>
c           *************
c
c
c
c
c           <<-- Ws12 BEGIN
c           ***************
c
27500 continue
c     prepare mu integration. 3-pion case 
      if (.not.i3ndmu.or.c3c4mu.ne.c4) go to 40
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi5 = wpi*c4*c4      
c
c     precalculate ImV(mu)
      if (.not.iv75) then
            iv75 = .true.
            do 27535 ii = 1, n3mu
c                 unitless u = mu/mpi
                  u = a3mu(ii)/wpi
                  u2 = u*u
                  u3 = u2*u
                  u4 = u2*u2
                  u5 = u3*u2
                  u6 = u3*u3
c                 uroot = wurz; uln = logneu                  
                  uroot = dsqrt(u2 - 2.d0*u - 3.d0)
                  uln = dlog((u - 1.d0 + uroot)/2.d0)
c
                  ANS1 = uroot*(u - 1.d0)*
     1         (  4.d0/3.d0*cb1*u*(u3 + 2.d0*u2 - u + 4.d0)
     2            + cb2/72.d0*(u6 + 2.d0*u5 - 39.d0*u4 - 12.d0*u3
     3                  + 65.d0*u2 - 50.d0*u - 27.d0)
     4            + cb3/12.d0*(u6 + 2.d0*u5 - 31.d0*u4 + 4.d0*u3
     5                  + 57.d0*u2 - 18.d0*u - 27.d0)
     6            + cb4/72.d0*(7.d0*u6 + 14.d0*u5 - 185.d0*u4
     7                  - 100.d0*u3 + 151.d0*u2 + 50.d0*u + 27.d0)  )
c     
                  ANS2 = uln*(16.d0*cb1*(4.d0*u2 - 1.d0 - u4)
     1            + 2.d0*cb2/3.d0*(2.d0 - 10.d0*u2 + 3.d0*u4)
     2            + 4.d0*cb3*u2*(u2 - 2.d0)
     3            + 2.d0*cb4/3.d0*(9.d0*u4 - 10.d0*u2 - 2.d0) )
c
                  v75(ii) = (ANS1 + ANS2)/u3
27535       continue
      end if
c
c     iterate over output points aa(i), store results in them
      do 27505 i=1,nt
            aakk=-deltaq(i,1)
            c5=c(mm,im)*wpi5*aakk*aakk
c
c           integrate over mu in [3m_pi, +infinity], and store result in aintmu
            aintmu = 0.d0
            do 27515 ii = 1, n3mu
                  aintmu = aintmu + w3mu(ii)*v75(ii) / 
     1                  (a3mu4(ii)*a3mu(ii)*(a3mu2(ii)+aakk))
27515       continue
c
            aa(i)=aa(i)*c5*aintmu
27505 continue
c
      mi=mi+2
      mm=mm+1
      go to 5999
c
c           Ws12 END -->>
c           *************
c
c
c
c
c           <<-- Wt12 BEGIN
c           ***************
c
27600 continue
c     prepare mu integration. 3-pion case 
      if (.not.i3ndmu.or.c3c4mu.ne.c4) go to 40
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi3 = wpi*c4
c
c     precalculate ImV(mu)
      if (.not.iv76) then
            iv76 = .true.
            do 27635 ii = 1, n3mu
c                 unitless u = mu/mpi
                  u = a3mu(ii)/wpi
                  u2 = u*u
                  u3 = u2*u
                  u4 = u2*u2
                  u5 = u3*u2
                  u6 = u3*u3
c                 uroot = wurz; uln = logneu                  
                  uroot = dsqrt(u2 - 2.d0*u - 3.d0)
                  uln = dlog((u - 1.d0 + uroot)/2.d0)
c
                  ANS1 = uroot*(u - 1.d0)*
     1         (  16.d0/3.d0*cb1*u*(2.d0 + u - 2.d0*u2 - u3)
     2            + cb2/36.d0*(73.d0*u4 - 6.d0*u5 - 3.d0*u6 + 44.d0*u3
     3                  - 43.d0*u2 - 50.d0*u - 27.d0)
     4            + cb3/2.d0*(19.d0*u4 - 2.d0*u5 - u6 + 4.d0*u3
     5                  - 9.d0*u2 - 6.d0*u - 9.d0)
     6            + cb4/36.d0*(39.d0*u4 - 2.d0*u5 - u6 + 12.d0*u3
     7                  - 65.d0*u2 + 50.d0*u + 27.d0)                )
c
                  ANS2 = 4.d0*uln*(8.d0*cb1*(u4 - 1.d0) 
     1            + cb2*(2.d0/3.d0 - u4) - 2.d0*cb3*u4
     2            + cb4/3.d0*(10.d0*u2 - 2.d0 - 3.d0*u4)  )
c
                  v76(ii) = (ANS1 + ANS2) / u5
27635       continue
      end if
c
c     iterate over output points aa(i), store results in them
      do 27605 i=1,nt
            aakk=-deltaq(i,1)
            c5=c(mm,im)*wpi3*aakk*aakk
c
c           integrate over mu in [3m_pi, +infinity], and store result in aintmu
            aintmu = 0.d0
            do 27615 ii = 1, n3mu
                  aintmu = aintmu + w3mu(ii)*v76(ii) / 
     1                  (a3mu2(ii)*a3mu(ii)*(a3mu2(ii)+aakk))
27615       continue
c
            aa(i)=aa(i)*c5*aintmu
27605 continue
c
      mi=mi+2
      mm=mm+1
      go to 5999
c
c           Wt12 END -->>
c           *************
c
c
c
c
c
c           <<-- Vt_class2 BEGIN
c           ********************

c
27700 continue
c     prepare mu integration. 2-pion case 
      if (.not.indmu.or.cc4mu.ne.c4) go to 30
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi5 = wpi*c4*c4
c
c     precalculate ImV(mu)
      if (.not.iv77) then
            iv77 = .true.
            do 27735 ii = 1, nmu
c                 unitless u = mu/mpi
                  u = amu(ii)/wpi
                  u2 = u*u
                  uroot = dsqrt(u2-4.d0)
c                  
c                 u-expression of eq (1)                  
                  v77(ii) = u*(  uroot*(u*(u2 - 30.d0) + 64.d0) 
     1            + 24.d0*(u2 - 5)*dlog((u + uroot)/2.d0)  )
27735       continue
      end if
c
c     iterate over output points aa(i), store results in them
      do 27705 i=1,nt
            aakk=-deltaq(i,1)
            c5=c(mm,im)*wpi5*aakk*aakk
c
c           integrate over mu in [2m_pi, +infinity], and store result in aintmu
            aintmu = 0.d0
            do 27715 ii = 1, nmu
                  aintmu = aintmu + wmu(ii)*v77(ii) / 
     1                  (amu4(ii)*amu(ii)*(amu2(ii)+aakk))
27715       continue
c
            aa(i)=aa(i)*c5*aintmu
27705 continue
c
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c           Vt_class2 END -->>
c           ******************
c
c
c
c           <<-- Wt_class2 BEGIN
c           ********************
c
27800 continue
c     prepare mu integration. 2-pion case 
      if (.not.indmu.or.cc4mu.ne.c4) go to 30
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi5 = wpi*c4*c4
c
c     precalculate ImV(mu)
      if (.not.iv78) then
            iv78 = .true.
            do 27835 ii = 1, nmu
c                 unitless u = mu/mpi
                  u = amu(ii)/wpi
                  u2 = u*u
                  uroot = dsqrt(u2-4.d0)
c                  
c                 u-expression of eq (2)                  
                  v78(ii) = (u2 - 4.d0)* 
     1         (  cb4/3.d0*
     2            (uroot*(2.d0*u2 - 8.d0)*dlog((u + uroot)/2.d0)
     3                        + u*(8.d0 + 36.d0*ga2 -5.d0/3.d0*u2))
     4            + 128.d0*pi2*fpi2*ce17*u*(u2 - 2.d0)              )
27835       continue
      end if
c
c     iterate over output points aa(i), store results in them
      do 27805 i=1,nt
            aakk=-deltaq(i,1)
            c5=c(mm,im)*wpi5*aakk*aakk
c
c           integrate over mu in [2m_pi, +infinity], and store result in aintmu
            aintmu = 0.d0
            do 27815 ii = 1, nmu
                  aintmu = aintmu + wmu(ii)*v78(ii) / 
     1                  (amu4(ii)*amu(ii)*(amu2(ii)+aakk))
27815       continue
c
            aa(i)=aa(i)*c5*aintmu
27805 continue
c
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c           Wt_class2 END -->>
c           ******************
c
c
c
c           <<-- Vc_class2 BEGIN
c           ********************
c
27900 continue
c     prepare mu integration. 2-pion case 
      if (.not.indmu.or.cc4mu.ne.c4) go to 30
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi5 = wpi*c4*c4
c
c     precalculate ImV(mu)
      if (.not.iv79) then
            iv79 = .true.
c
c           f_pi
            fpi = dsqrt(fpi2)
c      
            c51 = (2.d0/pi) * ga2*wpi5 / ((4.d0*fpi)**6) / pi2
            c52 = (2.d0/pi) * 3.d0*ga2*wpi5 / ((2.d0*fpi)**4)
c
            do 27935 ii = 1, nmu
c                 unitless u = mu/mpi
                  u = amu(ii)/wpi
                  u2 = u*u
                  u2m2 = u2 - 2.d0
                  uroot = dsqrt(u2-4.d0)
c                  
c                 ANS = eq (3) (times const_factor of [EM D.16])
                  v79(ii) = c51*u2m2*(2.d0 - 1.d0/u2)*
     1         (  2.d0*uroot*(24.d0*cb1 + cb2*(u2 - 4.d0) 
     2                  + 6.d0*cb3*u2m2)*dlog((u + uroot)/2)
     3            + u*(cb2*(8.d0 - 5.d0/3.d0*u2) - 6.d0*cb3*u2m2 
     4                                               - 24.0*cb1)  )
     5            + c52*ce14*(u2m2**3)/u
27935       continue
      end if
c
c     iterate over output points aa(i), store results in them
      do 27905 i=1,nt
c           aakk = q^2 
            aakk=-deltaq(i,1)
c
c           integrate over mu in [2m_pi, +infinity], and store result in aintmu
            aintmu = 0.d0
            do 27915 ii = 1, nmu
                  aintmu = aintmu + wmu(ii)*v79(ii) / 
     1                  (amu4(ii)*amu(ii)*(amu2(ii)+aakk))
27915       continue
c
            aa(i)=aa(i)*(aakk**3)*aintmu
27905 continue
c
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c           Vc_class2 END -->>
c           ******************
c
c
c
c           <<-- Wc_class2 BEGIN
c           ********************
c
28000 continue
c     prepare mu integration. 2-pion case 
      if (.not.indmu.or.cc4mu.ne.c4) go to 30
c
      wpi=dsqrt(c4)
      wpi5=c4*c4*wpi
      fpi=dsqrt(fpi2)
c      
c
      if (.not.iv80) then
      iv80 = .true.
c
      c51=2.d0*cb1*wpi5/((2.d0*fpi)**6*pi2*pi)
      c52=2.d0*cb2*wpi5/((2.d0*fpi)**6*pi2*pi)
      c53=2.d0*cb3*wpi5/((4.d0*fpi)**6*pi2*pi)
      c54=2.d0*cb4*ga2*wpi5/((4.d0*fpi)**6*pi2*pi)
c
      do 28035 ii = 1, nmu
c
      u=amu(ii)/wpi
      u2=u*u
      u4=u2*u2
      radi=u2-4.d0
      root=dsqrt(radi)
      alogfun=dlog((u+root)/2.d0)
c
      brak1=(3.d0*ga2+1.d0)/8.d0*root*(2.d0-u2)
     1     +((3.d0*ga2+1.d0)/u-2.d0*ga2*u)*alogfun
      brak2=root/96.d0*(7.d0*u2-6.d0-u4+ga2*(5.d0*u2-6.d0-2.d0*u4))
     1     +(ga2*u2-1.d0-ga2)/(4.d0*u)*alogfun
      brak3=2.d0/9.d0*root*(3.d0*(7.d0*u2-6.d0-u4)
     1     +4.d0*ga2*(32.d0/u-12.d0-20.d0*u+7.d0*u2-u4)
     2     +ga4*(114.d0-512.d0/u+368.d0*u
     3     -169.d0*u2+7.d0*u4+192.d0/(u+2.d0)))
     4     +16.d0/(3.d0*u)*(ga4*(6.d0*u4-30.d0*u2+35.d0)
     5     +ga2*(6.d0*u2-8.d0)-3.d0)*alogfun
      brak4=2.d0/9.d0*root*(30.d0-128.d0/u+80.d0*u-13.d0*u2-2.d0*u4
     1     +ga2*(512.d0/u-114.d0-368.d0*u
     2     +169.d0*u2-7.d0*u4-192.d0/(u+2.d0)))
     3     +16.d0/(3.d0*u)*(5.d0-3.d0*u2
     4     +ga2*(30.d0*u2-35.d0-6.d0*u4))*alogfun
c
      term1=c51*brak1
      term2=c52*brak2
      term3=c53*brak3
      term4=c54*brak4
      term=term1+term2+term3+term4
c
      v80(ii) = term
28035 continue
      end if
c
c     iterate over output points aa(i), store results in them
      do 28005 i=1,nt
c           aakk = q^2 
            aakk=-deltaq(i,1)
c
c           integrate over mu in [2m_pi, +infinity], and store result in aintmu
            aintmu = 0.d0
            do 28015 ii = 1, nmu
                  aintmu = aintmu + wmu(ii)*v80(ii) / 
     1                  (amu4(ii)*amu(ii)*(amu2(ii)+aakk))
28015       continue
c
            aa(i)=aa(i)*(aakk**3)*aintmu
28005 continue
c
      mi=mi+1
      mm=mm+1
      go to 5999
c
c
c           Wc_class2 END -->>
c           ******************
c
c
c
c           <<-- Vc leading relativistic corrections BEGIN
c           **********************************************
c
28100 continue
      if (.not.indla.or.cc4.ne.c4) go to 20
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi5 = wpi*c4*c4      
c
      c5=c(mm,im)
c
c     iterate over output points aa(i), store results in them
      do 28105 i=1,nt
c           aakk = q^2 
            aakk=-deltaq(i,1)
            omega2 = 4.d0*c4 + aakk
c
            aa(i)=aa(i) * c5 *
     1   (  wpi5/2.d0/omega2 + (2.d0*c4 + aakk)*(aakk - c4)*cpa(i) )
28105 continue
c
      mi=mi+1
      mm=mm+1
      go to 5999
c
c           Vc leading relativistic corrections END -->>
c           ********************************************
c
c
c
c           <<-- Wc leading relativistic corrections BEGIN
c           **********************************************
c
28200 continue
      if (.not.indla.or.cc4.ne.c4) go to 20
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
      wpi5 = wpi*c4*c4      
c
      c5=c(mm,im)
c
c     iterate over output points aa(i), store results in them
      do 28205 i=1,nt
c           aakk = q^2 
            aakk=-deltaq(i,1)
            omega2 = 4.d0*c4 + aakk
c
            aa(i)=aa(i) * c5 *
     1   (  3.d0*ga2*wpi5/2.d0/omega2 + 
     2            (ga2*(3.d0*c4 + 2.d0*aakk) - 2.d0*c4 - aakk)*
     3            (2.d0*c4 + aakk)*cpa(i)      )
28205 continue
c
      mi=mi+1
      mm=mm+1
      go to 5999
c
c           Wc leading relativistic corrections END -->>
c           ********************************************
c
c
c
c           <<-- Vt leading relativistic corrections BEGIN
c           **********************************************
c
28300 continue
      if (.not.indla.or.cc4.ne.c4) go to 20
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
c
      c5=c(mm,im)
c
c     iterate over output points aa(i), store results in them
      do 28305 i=1,nt
c           aakk = q^2 
            aakk=-deltaq(i,1)
c
            aa(i)=aa(i) * c5 * (5.d0*c4 + 2.d0*aakk)*cpa(i)
28305 continue
c
      mi=mi+1
      mm=mm+1
      go to 5999
c
c           Vt leading relativistic corrections END -->>
c           ********************************************
c
c
c
c           <<-- Wt leading relativistic corrections BEGIN
c           **********************************************
c
28400 continue
      if (.not.indla.or.cc4.ne.c4) go to 20
c
c     pion mass, c4 = m_pi^2
      wpi=dsqrt(c4)
c
      c5=c(mm,im)
c
c     iterate over output points aa(i), store results in them
      do 28405 i=1,nt
c           aakk = q^2 
            aakk=-deltaq(i,1)
            omega2 = 4.d0*c4 + aakk
c
            aa(i)=aa(i) * c5 * 
     1            (ga2*(3.d0*c4 + aakk) - omega2)*cpa(i)
28405 continue
c
      mi=mi+1
      mm=mm+1
      go to 5999
c
c           Wt leading relativistic corrections END -->>
c           ********************************************
c
c
c
c     <<-- Error section BEGIN
c     ************************
c
c
 9002 write (kwrite,19002) ityp
19002 format (1h ////' error in chiaa: invalid value of ityp =',I5,  
     1 /' execution terminated.'////)
      stop
c      
c     Error section END -->>
c     **********************
c
c
 8000 return
      end
      subroutine legp (pj,pjm1,x,j)
c
c
c        subroutine legp   computes the legendre polynominals
c
      real*8 pj,pjm1,x,a,b
c
c
c
c        compute legendre polynom for j equals zero
c
c
      if (j.gt.0) go to 1
      pj=1.d0
      pjm1=0.d0
      if (j.lt.0) pj=0.d0
      return
c
c
c
c        compute legendre polynoms for j equals one
c
c
c
    1 pj=x
      pjm1=1.d0
      if (j.eq.1) return
c
c
c
c        compute legendre polynom for j greater or equal two
c
c
c
      do 2 i=2,j
      a=x*pj
      b=a-pjm1
      pjm1=pj
    2 pj=-b/dfloat(i)+b+a
c
c
      return
      end
      subroutine gset(ax,bx,n,z,w)
c
c
c        this code has been obtained from the CERN computer library
c        in the year of the lord 1972.
c
c
      implicit real*8 (a-h,o-z)
c
c     n-point gauss zeros and weights for the interval (ax,bx) are
c           stored in  arrays z and w respectively.
c
      dimension     a(273),x(273),ktab(96)
      dimension z(2),w(2)
c
c-----table of initial subscripts for n=2(1)16(4)96
      data ktab(2)/1/
      data ktab(3)/2/
      data ktab(4)/4/
      data ktab(5)/6/
      data ktab(6)/9/
      data ktab(7)/12/
      data ktab(8)/16/
      data ktab(9)/20/
      data ktab(10)/25/
      data ktab(11)/30/
      data ktab(12)/36/
      data ktab(13)/42/
      data ktab(14)/49/
      data ktab(15)/56/
      data ktab(16)/64/
      data ktab(20)/72/
      data ktab(24)/82/
      data ktab(28)/82/
      data ktab(32)/94/
      data ktab(36)/94/
      data ktab(40)/110/
      data ktab(44)/110/
      data ktab(48)/130/
      data ktab(52)/130/
      data ktab(56)/130/
      data ktab(60)/130/
      data ktab(64)/154/
      data ktab(68)/154/
      data ktab(72)/154/
      data ktab(76)/154/
      data ktab(80)/186/
      data ktab(84)/186/
      data ktab(88)/186/
      data ktab(92)/186/
      data ktab(96)/226/
c
c-----table of abscissae (x) and weights (a) for interval (-1,+1).
c
c**** n=2
      data x(1)/0.577350269189626  d0/, a(1)/1.000000000000000  d0/
c**** n=3
      data x(2)/0.774596669241483  d0/, a(2)/0.555555555555556  d0/
      data x(3)/0.000000000000000  d0/, a(3)/0.888888888888889  d0/
c**** n=4
      data x(4)/0.861136311594053  d0/, a(4)/0.347854845137454  d0/
      data x(5)/0.339981043584856  d0/, a(5)/0.652145154862546  d0/
c**** n=5
      data x(6)/0.906179845938664  d0/, a(6)/0.236926885056189  d0/
      data x(7)/0.538469310105683  d0/, a(7)/0.478628670499366  d0/
      data x(8)/0.000000000000000  d0/, a(8)/0.568888888888889  d0/
c**** n=6
      data x(9)/0.932469514203152  d0/, a(9)/0.171324492379170  d0/
      data x(10)/0.661209386466265 d0/, a(10)/0.360761573048139 d0/
      data x(11)/0.238619186083197 d0/, a(11)/0.467913934572691 d0/
c**** n=7
      data x(12)/0.949107912342759 d0/, a(12)/0.129484966168870 d0/
      data x(13)/0.741531185599394 d0/, a(13)/0.279705391489277 d0/
      data x(14)/0.405845151377397 d0/, a(14)/0.381830050505119 d0/
      data x(15)/0.000000000000000 d0/, a(15)/0.417959183673469 d0/
c**** n=8
      data x(16)/0.960289856497536 d0/, a(16)/0.101228536290376 d0/
      data x(17)/0.796666477413627 d0/, a(17)/0.222381034453374 d0/
      data x(18)/0.525532409916329 d0/, a(18)/0.313706645877887 d0/
      data x(19)/0.183434642495650 d0/, a(19)/0.362683783378362 d0/
c**** n=9
      data x(20)/0.968160239507626 d0/, a(20)/0.081274388361574 d0/
      data x(21)/0.836031107326636 d0/, a(21)/0.180648160694857 d0/
      data x(22)/0.613371432700590 d0/, a(22)/0.260610696402935 d0/
      data x(23)/0.324253423403809 d0/, a(23)/0.312347077040003 d0/
      data x(24)/0.000000000000000 d0/, a(24)/0.330239355001260 d0/
c**** n=10
      data x(25)/0.973906528517172 d0/, a(25)/0.066671344308688 d0/
      data x(26)/0.865063366688985 d0/, a(26)/0.149451349150581 d0/
      data x(27)/0.679409568299024 d0/, a(27)/0.219086362515982 d0/
      data x(28)/0.433395394129247 d0/, a(28)/0.269266719309996 d0/
      data x(29)/0.148874338981631 d0/, a(29)/0.295524224714753 d0/
c**** n=11
      data x(30)/0.978228658146057 d0/, a(30)/0.055668567116174 d0/
      data x(31)/0.887062599768095 d0/, a(31)/0.125580369464905 d0/
      data x(32)/0.730152005574049 d0/, a(32)/0.186290210927734 d0/
      data x(33)/0.519096129206812 d0/, a(33)/0.233193764591990 d0/
      data x(34)/0.269543155952345 d0/, a(34)/0.262804544510247 d0/
      data x(35)/0.000000000000000 d0/, a(35)/0.272925086777901 d0/
c**** n=12
      data x(36)/0.981560634246719 d0/, a(36)/0.047175336386512 d0/
      data x(37)/0.904117256370475 d0/, a(37)/0.106939325995318 d0/
      data x(38)/0.769902674194305 d0/, a(38)/0.160078328543346 d0/
      data x(39)/0.587317954286617 d0/, a(39)/0.203167426723066 d0/
      data x(40)/0.367831498998180 d0/, a(40)/0.233492536538355 d0/
      data x(41)/0.125233408511469 d0/, a(41)/0.249147045813403 d0/
c**** n=13
      data x(42)/0.984183054718588 d0/, a(42)/0.040484004765316 d0/
      data x(43)/0.917598399222978 d0/, a(43)/0.092121499837728 d0/
      data x(44)/0.801578090733310 d0/, a(44)/0.138873510219787 d0/
      data x(45)/0.642349339440340 d0/, a(45)/0.178145980761946 d0/
      data x(46)/0.448492751036447 d0/, a(46)/0.207816047536889 d0/
      data x(47)/0.230458315955135 d0/, a(47)/0.226283180262897 d0/
      data x(48)/0.000000000000000 d0/, a(48)/0.232551553230874 d0/
c**** n=14
      data x(49)/0.986283808696812 d0/, a(49)/0.035119460331752 d0/
      data x(50)/0.928434883663574 d0/, a(50)/0.080158087159760 d0/
      data x(51)/0.827201315069765 d0/, a(51)/0.121518570687903 d0/
      data x(52)/0.687292904811685 d0/, a(52)/0.157203167158194 d0/
      data x(53)/0.515248636358154 d0/, a(53)/0.185538397477938 d0/
      data x(54)/0.319112368927890 d0/, a(54)/0.205198463721296 d0/
      data x(55)/0.108054948707344 d0/, a(55)/0.215263853463158 d0/
c**** n=15
      data x(56)/0.987992518020485 d0/, a(56)/0.030753241996117 d0/
      data x(57)/0.937273392400706 d0/, a(57)/0.070366047488108 d0/
      data x(58)/0.848206583410427 d0/, a(58)/0.107159220467172 d0/
      data x(59)/0.724417731360170 d0/, a(59)/0.139570677926154 d0/
      data x(60)/0.570972172608539 d0/, a(60)/0.166269205816994 d0/
      data x(61)/0.394151347077563 d0/, a(61)/0.186161000015562 d0/
      data x(62)/0.201194093997435 d0/, a(62)/0.198431485327111 d0/
      data x(63)/0.000000000000000 d0/, a(63)/0.202578241925561 d0/
c**** n=16
      data x(64)/0.989400934991650 d0/, a(64)/0.027152459411754 d0/
      data x(65)/0.944575023073233 d0/, a(65)/0.062253523938648 d0/
      data x(66)/0.865631202387832 d0/, a(66)/0.095158511682493 d0/
      data x(67)/0.755404408355003 d0/, a(67)/0.124628971255534 d0/
      data x(68)/0.617876244402644 d0/, a(68)/0.149595988816577 d0/
      data x(69)/0.458016777657227 d0/, a(69)/0.169156519395003 d0/
      data x(70)/0.281603550779259 d0/, a(70)/0.182603415044924 d0/
      data x(71)/0.095012509837637 d0/, a(71)/0.189450610455069 d0/
c**** n=20
      data x(72)/0.993128599185094 d0/, a(72)/0.017614007139152 d0/
      data x(73)/0.963971927277913 d0/, a(73)/0.040601429800386 d0/
      data x(74)/0.912234428251325 d0/, a(74)/0.062672048334109 d0/
      data x(75)/0.839116971822218 d0/, a(75)/0.083276741576704 d0/
      data x(76)/0.746331906460150 d0/, a(76)/0.101930119817240 d0/
      data x(77)/0.636053680726515 d0/, a(77)/0.118194531961518 d0/
      data x(78)/0.510867001950827 d0/, a(78)/0.131688638449176 d0/
      data x(79)/0.373706088715419 d0/, a(79)/0.142096109318382 d0/
      data x(80)/0.227785851141645 d0/, a(80)/0.149172986472603 d0/
      data x(81)/0.076526521133497 d0/, a(81)/0.152753387130725 d0/
c**** n=24
      data x(82)/0.995187219997021 d0/, a(82)/0.012341229799987 d0/
      data x(83)/0.974728555971309 d0/, a(83)/0.028531388628933 d0/
      data x(84)/0.938274552002732 d0/, a(84)/0.044277438817419 d0/
      data x(85)/0.886415527004401 d0/, a(85)/0.059298584915436 d0/
      data x(86)/0.820001985973902 d0/, a(86)/0.073346481411080 d0/
      data x(87)/0.740124191578554 d0/, a(87)/0.086190161531953 d0/
      data x(88)/0.648093651936975 d0/, a(88)/0.097618652104113 d0/
      data x(89)/0.545421471388839 d0/, a(89)/0.107444270115965 d0/
      data x(90)/0.433793507626045 d0/, a(90)/0.115505668053725 d0/
      data x(91)/0.315042679696163 d0/, a(91)/0.121670472927803 d0/
      data x(92)/0.191118867473616 d0/, a(92)/0.125837456346828 d0/
      data x(93)/0.064056892862605 d0/, a(93)/0.127938195346752 d0/
c**** n=32
      data x(94)/0.997263861849481 d0/, a(94)/0.007018610009470 d0/
      data x(95)/0.985611511545268 d0/, a(95)/0.016274394730905 d0/
      data x(96)/0.964762255587506 d0/, a(96)/0.025392065309262 d0/
      data x(97)/0.934906075937739 d0/, a(97)/0.034273862913021 d0/
      data x(98)/0.896321155766052 d0/, a(98)/0.042835898022226 d0/
      data x(99)/0.849367613732569 d0/, a(99)/0.050998059262376 d0/
      data x(100)/0.794483795967942d0/, a(100)/0.058684093478535d0/
      data x(101)/0.732182118740289d0/, a(101)/0.065822222776361d0/
      data x(102)/0.663044266930215d0/, a(102)/0.072345794108848d0/
      data x(103)/0.587715757240762d0/, a(103)/0.078193895787070d0/
      data x(104)/0.506899908932229d0/, a(104)/0.083311924226946d0/
      data x(105)/0.421351276130635d0/, a(105)/0.087652093004403d0/
      data x(106)/0.331868602282127d0/, a(106)/0.091173878695763d0/
      data x(107)/0.239287362252137d0/, a(107)/0.093844399080804d0/
      data x(108)/0.144471961582796d0/, a(108)/0.095638720079274d0/
      data x(109)/0.048307665687738d0/, a(109)/0.096540088514727d0/
c**** n=40
      data x(110)/0.998237709710559d0/, a(110)/0.004521277098533d0/
      data x(111)/0.990726238699457d0/, a(111)/0.010498284531152d0/
      data x(112)/0.977259949983774d0/, a(112)/0.016421058381907d0/
      data x(113)/0.957916819213791d0/, a(113)/0.022245849194166d0/
      data x(114)/0.932812808278676d0/, a(114)/0.027937006980023d0/
      data x(115)/0.902098806968874d0/, a(115)/0.033460195282547d0/
      data x(116)/0.865959503212259d0/, a(116)/0.038782167974472d0/
      data x(117)/0.824612230833311d0/, a(117)/0.043870908185673d0/
      data x(118)/0.778305651426519d0/, a(118)/0.048695807635072d0/
      data x(119)/0.727318255189927d0/, a(119)/0.053227846983936d0/
      data x(120)/0.671956684614179d0/, a(120)/0.057439769099391d0/
      data x(121)/0.612553889667980d0/, a(121)/0.061306242492928d0/
      data x(122)/0.549467125095128d0/, a(122)/0.064804013456601d0/
      data x(123)/0.483075801686178d0/, a(123)/0.067912045815233d0/
      data x(124)/0.413779204371605d0/, a(124)/0.070611647391286d0/
      data x(125)/0.341994090825758d0/, a(125)/0.072886582395804d0/
      data x(126)/0.268152185007253d0/, a(126)/0.074723169057968d0/
      data x(127)/0.192697580701371d0/, a(127)/0.076110361900626d0/
      data x(128)/0.116084070675255d0/, a(128)/0.077039818164247d0/
      data x(129)/0.038772417506050d0/, a(129)/0.077505947978424d0/
c**** n=48
      data x(130)/0.998771007252426d0/, a(130)/0.003153346052305d0/
      data x(131)/0.993530172266350d0/, a(131)/0.007327553901276d0/
      data x(132)/0.984124583722826d0/, a(132)/0.011477234579234d0/
      data x(133)/0.970591592546247d0/, a(133)/0.015579315722943d0/
      data x(134)/0.952987703160430d0/, a(134)/0.019616160457355d0/
      data x(135)/0.931386690706554d0/, a(135)/0.023570760839324d0/
      data x(136)/0.905879136715569d0/, a(136)/0.027426509708356d0/
      data x(137)/0.876572020274247d0/, a(137)/0.031167227832798d0/
      data x(138)/0.843588261624393d0/, a(138)/0.034777222564770d0/
      data x(139)/0.807066204029442d0/, a(139)/0.038241351065830d0/
      data x(140)/0.767159032515740d0/, a(140)/0.041545082943464d0/
      data x(141)/0.724034130923814d0/, a(141)/0.044674560856694d0/
      data x(142)/0.677872379632663d0/, a(142)/0.047616658492490d0/
      data x(143)/0.628867396776513d0/, a(143)/0.050359035553854d0/
      data x(144)/0.577224726083972d0/, a(144)/0.052890189485193d0/
      data x(145)/0.523160974722233d0/, a(145)/0.055199503699984d0/
      data x(146)/0.466902904750958d0/, a(146)/0.057277292100403d0/
      data x(147)/0.408686481990716d0/, a(147)/0.059114839698395d0/
      data x(148)/0.348755886292160d0/, a(148)/0.060704439165893d0/
      data x(149)/0.287362487355455d0/, a(149)/0.062039423159892d0/
      data x(150)/0.224763790394689d0/, a(150)/0.063114192286254d0/
      data x(151)/0.161222356068891d0/, a(151)/0.063924238584648d0/
      data x(152)/0.097004699209462d0/, a(152)/0.064466164435950d0/
      data x(153)/0.032380170962869d0/, a(153)/0.064737696812683d0/
c**** n=64
      data x(154)/0.999305041735772d0/, a(154)/0.001783280721696d0/
      data x(155)/0.996340116771955d0/, a(155)/0.004147033260562d0/
      data x(156)/0.991013371476744d0/, a(156)/0.006504457968978d0/
      data x(157)/0.983336253884625d0/, a(157)/0.008846759826363d0/
      data x(158)/0.973326827789910d0/, a(158)/0.011168139460131d0/
      data x(159)/0.961008799652053d0/, a(159)/0.013463047896718d0/
      data x(160)/0.946411374858402d0/, a(160)/0.015726030476024d0/
      data x(161)/0.929569172131939d0/, a(161)/0.017951715775697d0/
      data x(162)/0.910522137078502d0/, a(162)/0.020134823153530d0/
      data x(163)/0.889315445995114d0/, a(163)/0.022270173808383d0/
      data x(164)/0.865999398154092d0/, a(164)/0.024352702568710d0/
      data x(165)/0.840629296252580d0/, a(165)/0.026377469715054d0/
      data x(166)/0.813265315122797d0/, a(166)/0.028339672614259d0/
      data x(167)/0.783972358943341d0/, a(167)/0.030234657072402d0/
      data x(168)/0.752819907260531d0/, a(168)/0.032057928354851d0/
      data x(169)/0.719881850171610d0/, a(169)/0.033805161837141d0/
      data x(170)/0.685236313054233d0/, a(170)/0.035472213256882d0/
      data x(171)/0.648965471254657d0/, a(171)/0.037055128540240d0/
      data x(172)/0.611155355172393d0/, a(172)/0.038550153178615d0/
      data x(173)/0.571895646202634d0/, a(173)/0.039953741132720d0/
      data x(174)/0.531279464019894d0/, a(174)/0.041262563242623d0/
      data x(175)/0.489403145707052d0/, a(175)/0.042473515123653d0/
      data x(176)/0.446366017253464d0/, a(176)/0.043583724529323d0/
      data x(177)/0.402270157963991d0/, a(177)/0.044590558163756d0/
      data x(178)/0.357220158337668d0/, a(178)/0.045491627927418d0/
      data x(179)/0.311322871990210d0/, a(179)/0.046284796581314d0/
      data x(180)/0.264687162208767d0/, a(180)/0.046968182816210d0/
      data x(181)/0.217423643740007d0/, a(181)/0.047540165714830d0/
      data x(182)/0.169644420423992d0/, a(182)/0.047999388596458d0/
      data x(183)/0.121462819296120d0/, a(183)/0.048344762234802d0/
      data x(184)/0.072993121787799d0/, a(184)/0.048575467441503d0/
      data x(185)/0.024350292663424d0/, a(185)/0.048690957009139d0/
c**** n=80
      data x(186)/0.999553822651630d0/, a(186)/0.001144950003186d0/
      data x(187)/0.997649864398237d0/, a(187)/0.002663533589512d0/
      data x(188)/0.994227540965688d0/, a(188)/0.004180313124694d0/
      data x(189)/0.989291302499755d0/, a(189)/0.005690922451403d0/
      data x(190)/0.982848572738629d0/, a(190)/0.007192904768117d0/
      data x(191)/0.974909140585727d0/, a(191)/0.008683945269260d0/
      data x(192)/0.965485089043799d0/, a(192)/0.010161766041103d0/
      data x(193)/0.954590766343634d0/, a(193)/0.011624114120797d0/
      data x(194)/0.942242761309872d0/, a(194)/0.013068761592401d0/
      data x(195)/0.928459877172445d0/, a(195)/0.014493508040509d0/
      data x(196)/0.913263102571757d0/, a(196)/0.015896183583725d0/
      data x(197)/0.896675579438770d0/, a(197)/0.017274652056269d0/
      data x(198)/0.878722567678213d0/, a(198)/0.018626814208299d0/
      data x(199)/0.859431406663111d0/, a(199)/0.019950610878141d0/
      data x(200)/0.838831473580255d0/, a(200)/0.021244026115782d0/
      data x(201)/0.816954138681463d0/, a(201)/0.022505090246332d0/
      data x(202)/0.793832717504605d0/, a(202)/0.023731882865930d0/
      data x(203)/0.769502420135041d0/, a(203)/0.024922535764115d0/
      data x(204)/0.744000297583597d0/, a(204)/0.026075235767565d0/
      data x(205)/0.717365185362099d0/, a(205)/0.027188227500486d0/
      data x(206)/0.689637644342027d0/, a(206)/0.028259816057276d0/
      data x(207)/0.660859898986119d0/, a(207)/0.029288369583267d0/
      data x(208)/0.631075773046871d0/, a(208)/0.030272321759557d0/
      data x(209)/0.600330622829751d0/, a(209)/0.031210174188114d0/
      data x(210)/0.568671268122709d0/, a(210)/0.032100498673487d0/
      data x(211)/0.536145920897131d0/, a(211)/0.032941939397645d0/
      data x(212)/0.502804111888784d0/, a(212)/0.033733214984611d0/
      data x(213)/0.468696615170544d0/, a(213)/0.034473120451753d0/
      data x(214)/0.433875370831756d0/, a(214)/0.035160529044747d0/
      data x(215)/0.398393405881969d0/, a(215)/0.035794393953416d0/
      data x(216)/0.362304753499487d0/, a(216)/0.036373749905835d0/
      data x(217)/0.325664370747701d0/, a(217)/0.036897714638276d0/
      data x(218)/0.288528054884511d0/, a(218)/0.037365490238730d0/
      data x(219)/0.250952358392272d0/, a(219)/0.037776364362001d0/
      data x(220)/0.212994502857666d0/, a(220)/0.038129711314477d0/
      data x(221)/0.174712291832646d0/, a(221)/0.038424993006959d0/
      data x(222)/0.136164022809143d0/, a(222)/0.038661759774076d0/
      data x(223)/0.097408398441584d0/, a(223)/0.038839651059051d0/
      data x(224)/0.058504437152420d0/, a(224)/0.038958395962769d0/
      data x(225)/0.019511383256793d0/, a(225)/0.039017813656306d0/
c**** n=96
      data x(226)/0.999689503883230d0/, a(226)/0.000796792065552d0/
      data x(227)/0.998364375863181d0/, a(227)/0.001853960788946d0/
      data x(228)/0.995981842987209d0/, a(228)/0.002910731817934d0/
      data x(229)/0.992543900323762d0/, a(229)/0.003964554338444d0/
      data x(230)/0.988054126329623d0/, a(230)/0.005014202742927d0/
      data x(231)/0.982517263563014d0/, a(231)/0.006058545504235d0/
      data x(232)/0.975939174585136d0/, a(232)/0.007096470791153d0/
      data x(233)/0.968326828463264d0/, a(233)/0.008126876925698d0/
      data x(234)/0.959688291448742d0/, a(234)/0.009148671230783d0/
      data x(235)/0.950032717784437d0/, a(235)/0.010160770535008d0/
      data x(236)/0.939370339752755d0/, a(236)/0.011162102099838d0/
      data x(237)/0.927712456722308d0/, a(237)/0.012151604671088d0/
      data x(238)/0.915071423120898d0/, a(238)/0.013128229566961d0/
      data x(239)/0.901460635315852d0/, a(239)/0.014090941772314d0/
      data x(240)/0.886894517402420d0/, a(240)/0.015038721026994d0/
      data x(241)/0.871388505909296d0/, a(241)/0.015970562902562d0/
      data x(242)/0.854959033434601d0/, a(242)/0.016885479864245d0/
      data x(243)/0.837623511228187d0/, a(243)/0.017782502316045d0/
      data x(244)/0.819400310737931d0/, a(244)/0.018660679627411d0/
      data x(245)/0.800308744139140d0/, a(245)/0.019519081140145d0/
      data x(246)/0.780369043867433d0/, a(246)/0.020356797154333d0/
      data x(247)/0.759602341176647d0/, a(247)/0.021172939892191d0/
      data x(248)/0.738030643744400d0/, a(248)/0.021966644438744d0/
      data x(249)/0.715676812348967d0/, a(249)/0.022737069658329d0/
      data x(250)/0.692564536642171d0/, a(250)/0.023483399085926d0/
      data x(251)/0.668718310043916d0/, a(251)/0.024204841792364d0/
      data x(252)/0.644163403784967d0/, a(252)/0.024900633222483d0/
      data x(253)/0.618925840125468d0/, a(253)/0.025570036005349d0/
      data x(254)/0.593032364777572d0/, a(254)/0.026212340735672d0/
      data x(255)/0.566510418561397d0/, a(255)/0.026826866725591d0/
      data x(256)/0.539388108324357d0/, a(256)/0.027412962726029d0/
      data x(257)/0.511694177154667d0/, a(257)/0.027970007616848d0/
      data x(258)/0.483457973920596d0/, a(258)/0.028497411065085d0/
      data x(259)/0.454709422167743d0/, a(259)/0.028994614150555d0/
      data x(260)/0.425478988407300d0/, a(260)/0.029461089958167d0/
      data x(261)/0.395797649828908d0/, a(261)/0.029896344136328d0/
      data x(262)/0.365696861472313d0/, a(262)/0.030299915420827d0/
      data x(263)/0.335208522892625d0/, a(263)/0.030671376123669d0/
      data x(264)/0.304364944354496d0/, a(264)/0.031010332586313d0/
      data x(265)/0.273198812591049d0/, a(265)/0.031316425596861d0/
      data x(266)/0.241743156163840d0/, a(266)/0.031589330770727d0/
      data x(267)/0.210031310460567d0/, a(267)/0.031828758894411d0/
      data x(268)/0.178096882367618d0/, a(268)/0.032034456231992d0/
      data x(269)/0.145973714654896d0/, a(269)/0.032206204794030d0/
      data x(270)/0.113695850110665d0/, a(270)/0.032343822568575d0/
      data x(271)/0.081297495464425d0/, a(271)/0.032447163714064d0/
      data x(272)/0.048812985136049d0/, a(272)/0.032516118713868d0/
      data x(273)/0.016276744849602d0/, a(273)/0.032550614492363d0/
c
c
c-----test n
      alpha=0.5d0*(ax+bx)
      beta=0.5d0*(bx-ax)
      if( n.lt.1 .or. n.gt.96 ) go to 100
      if(n.ne.1) go to 1
      z(1)=alpha
      w(1)=bx-ax
      return
c
    1 if (n.le.16) go to 3
      if (n.gt.24) go to 4
      n=4*(n/4)
      go to 3
    4 if (n.gt.48) go to 5
      n=8*(n/8)
      go to 3
    5 n=16*(n/16)
c
c----- set k equal to initial subscript and store results
    3 k=ktab(n)
      m=n/2
      do 2 j=1,m
      jtab=k-1+j
      wtemp=beta*a(jtab)
      delta=beta*x(jtab)
      z(j)=alpha-delta
      w(j)=wtemp
      jp=n+1-j
      z(jp)=alpha+delta
      w(jp)=wtemp
    2 continue
      if((n-m-m).eq.0) return
      z(m+1)=alpha
      jmid=k+m
      w(m+1)=beta*a(jmid)
      return
c
  100 zn=n
      write(6,200) zn
  200 format(/////' error in gset. n has the non-permissible value',
     1e11.3/' execution terminated.')
      stop
      end
c name:    dgelg
c        programmbibliothek rhrz bonn        02/02/81       dgelg
c                                            fortran iv     ibm 370/168
c
c purpose:
c
c to solve a general system of simultaneous linear equations.
c
c usage:   call dgelg(r,a,m,n,eps,ier)
c
c parameters:
c
c r:       double precision m by n right hand side matrix
c          (destroyed). on return r contains the solutions
c          of the equations.
c
c a:       double precision m by m coefficient matrix
c          (destroyed).
c
c m:       the number of equations in the system.
c
c n:       the number of right hand side vectors.
c
c eps:     single precision input constant which is used as
c          relative tolerance for test on loss of
c          significance.
c
c ier:     resulting error parameter coded as follows
c           ier=0  - no error,
c           ier=-1 - no result because of m less than 1 or
c                   pivot element at any elimination step
c                   equal to 0,
c           ier=k  - warning due to possible loss of signifi-
c                   cance indicated at elimination step k+1,
c                   where pivot element was less than or
c                   equal to the internal tolerance eps times
c                   absolutely greatest element of matrix a.
c
c remarks: (1) input matrices r and a are assumed to be stored
c              columnwise in m*n resp. m*m successive storage
c              locations. on return solution matrix r is stored
c              columnwise too.
c          (2) the procedure gives results if the number of equations m
c              is greater than 0 and pivot elements at all elimination
c              steps are different from 0. however warning ier=k - if
c              given indicates possible loss of significance. in case
c              of a well scaled matrix a and appropriate tolerance eps,
c              ier=k may be interpreted that matrix a has the rank k.
c              no warning is given in case m=1.
c
c method:
c
c solution is done by means of gauss-elimination with
c complete pivoting.
c
c programs required:
c          none
c
c access:
c
c load module:    sys3.fortlib(dgelg)
c source module:  sys3.symlib.fortran(dgelg)
c description:    sys3.infolib(dgelg)
c
c author:         ibm, ssp iii
c installation:   ibm 370/168, mvs-jes2, fortran iv (h ext. enh.)
c
c**********************************************************************
      subroutine dgelg(r,a,m,n,eps,ier)
c
c
      implicit real*8 (a-h,o-z)
      dimension a(1),r(1)
      real*4 eps
c
c
c
c
      if(m)23,23,1
c
c     search for greatest element in matrix a
    1 ier=0
      piv=0.d0
      mm=m*m
      nm=n*m
      do 3 l=1,mm
      tb=dabs(a(l))
      if(tb-piv)3,3,2
    2 piv=tb
      i=l
    3 continue
      tol=eps*piv
c     a(i) is pivot element. piv contains the absolute value of a(i).
c
c
c     start elimination loop
      lst=1
      do 17 k=1,m
c
c     test on singularity
      if(piv)23,23,4
    4 if(ier)7,5,7
    5 if(piv-tol)6,6,7
    6 ier=k-1
    7 pivi=1.d0/a(i)
      j=(i-1)/m
      i=i-j*m-k
      j=j+1-k
c     i+k is row-index, j+k column-index of pivot element
c
c     pivot row reduction and row interchange in right hand side r
      do 8 l=k,nm,m
      ll=l+i
      tb=pivi*r(ll)
      r(ll)=r(l)
    8 r(l)=tb
c
c     is elimination terminated
      if(k-m)9,18,18
c
c     column interchange in matrix a
    9 lend=lst+m-k
      if(j)12,12,10
   10 ii=j*m
      do 11 l=lst,lend
      tb=a(l)
      ll=l+ii
      a(l)=a(ll)
   11 a(ll)=tb
c
c     row interchange and pivot row reduction in matrix a
   12 do 13 l=lst,mm,m
      ll=l+i
      tb=pivi*a(ll)
      a(ll)=a(l)
   13 a(l)=tb
c
c     save column interchange information
      a(lst)=j
c
c     element reduction and next pivot search
      piv=0.d0
      lst=lst+1
      j=0
      do 16 ii=lst,lend
      pivi=-a(ii)
      ist=ii+m
      j=j+1
      do 15 l=ist,mm,m
      ll=l-j
      a(l)=a(l)+pivi*a(ll)
      tb=dabs(a(l))
      if(tb-piv)15,15,14
   14 piv=tb
      i=l
   15 continue
      do 16 l=k,nm,m
      ll=l+j
   16 r(ll)=r(ll)+pivi*r(l)
   17 lst=lst+m
c     end of elimination loop
c
c
c     back substitution and back interchange
   18 if(m-1)23,22,19
   19 ist=mm+m
      lst=m+1
      do 21 i=2,m
      ii=lst-i
      ist=ist-lst
      l=ist-m
      l=a(l)+.5d0
      do 21 j=ii,nm,m
      tb=r(j)
      ll=j
      do 20 k=ist,mm,m
      ll=ll+1
   20 tb=tb-a(k)*r(ll)
      k=j+l
      r(j)=r(k)
   21 r(k)=tb
   22 return
c
c
c     error return
   23 ier=-1
      return
      end
c********************************************************************
c name:    dminv
c        programmbibliothek rhrz bonn        28/11/78       dminv
c                                            fortran iv     ibm 370/168
c
c purpose:
c
c invert a matrix
c
c usage:   call dminv (a,n,d,l,m)
c
c parameters:
c
c a:       input matrix, destroyed in computation and replaced by
c          resultant inverse.
c          double precision required.
c
c n:       order of matrix a
c
c d:       resultant determinant
c          double precision required.
c
c l:       work vector of length n
c
c m:       work vector of length n
c
c remarks: matrix a must be a general matrix
c
c method:
c
c the standard gauss-jordan method is used. the determinant
c is also calculated. a determinant of zero indicates that
c the matrix is singular.
c
c programs required:
c          none
c
c author:  ibm, ssp iii
c
c**********************************************************************
      subroutine dminv (a,n,d,l,m)
      implicit real*8 (a-h,o-z)
      dimension a(1),l(1),m(1)

c
c
c        search for largest element
c
      d=1.d0
      nk=-n
      do 80 k=1,n
      nk=nk+n
      l(k)=k
      m(k)=k
      kk=nk+k
      biga=a(kk)
      do 20 j=k,n
      iz=n*(j-1)
      do 20 i=k,n
      ij=iz+i
   10 if (dabs(biga)-dabs(a(ij)))  15,20,20
   15 biga=a(ij)
      l(k)=i
      m(k)=j
   20 continue
c
c        interchange rows
c
      j=l(k)
      if(j-k) 35,35,25
   25 ki=k-n
      do 30 i=1,n
      ki=ki+n
      hold=-a(ki)
      ji=ki-k+j
      a(ki)=a(ji)
   30 a(ji) =hold
c
c        interchange columns
c
   35 i=m(k)
      if(i-k) 45,45,38
   38 jp=n*(i-1)
      do 40 j=1,n
      jk=nk+j
      ji=jp+j
      hold=-a(jk)
      a(jk)=a(ji)
   40 a(ji) =hold
c
c        divide column by minus pivot (value of pivot element is
c        contained in biga)
c
   45 if(biga) 48,46,48
   46 d=0.d0
      return
   48 do 55 i=1,n
      if(i-k) 50,55,50
   50 ik=nk+i
      a(ik)=a(ik)/(-biga)
   55 continue
c
c        reduce matrix
c
      do 65 i=1,n
      ik=nk+i
      hold=a(ik)
      ij=i-n
      do 65 j=1,n
      ij=ij+n
      if(i-k) 60,65,60
   60 if(j-k) 62,65,62
   62 kj=ij-i+k
      a(ij)=hold*a(kj)+a(ij)
   65 continue
c
c        divide row by pivot
c
      kj=k-n
      do 75 j=1,n
      kj=kj+n
      if(j-k) 70,75,70
   70 a(kj)=a(kj)/biga
   75 continue
c
c        product of pivots
c
      d=d*biga
c
c        replace pivot by reciprocal
c
      a(kk)=1.d0/biga
   80 continue
c
c        final row and column interchange
c
      k=n
  100 k=(k-1)
      if(k) 150,150,105
  105 i=l(k)
      if(i-k) 120,120,108
  108 jq=n*(k-1)
      jr=n*(i-1)
      do 110 j=1,n
      jk=jq+j
      hold=a(jk)
      ji=jr+j
      a(jk)=-a(ji)
  110 a(ji) =hold
  120 j=m(k)
      if(j-k) 100,100,125
  125 ki=k-n
      do 130 i=1,n
      ki=ki+n
      hold=a(ki)
      ji=ki-k+j
      a(ki)=-a(ji)
  130 a(ji) =hold
      go to 100
  150 return
      end
c**************** this is the end of the program n4lo550 ********************
      END MODULE
