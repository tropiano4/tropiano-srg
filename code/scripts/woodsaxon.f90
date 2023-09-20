!-------------------------------------------------------------------------------
! Use Numerov algorithm with bisection on energies to find the single particle 
! orbitals and single particle energies for a Wood-Saxon well
!
! References: arXiv:0709.3525 [nucl-th]
!             https://nucracker.volya.net
!
! Version history
!  1 - 2016-04-13 - Lovato    - original code
!  2 - 2017-02-09 - Lonardoni - module, updates and bug fix
!  3 - 2017-02-14 - Lonardoni - fix ls potential and parameters
!  4 - 2017-02-15 - Lonardoni - vws and v0 sign change
!  5 - 2017-02-17 - Lonardoni - rdiv added and tables fixed
!  6 - 2017-02-27 - Lonardoni - minor changes, prntorb added
!  7 - 2017-07-19 - Lonardoni - wine-bottle added
!  8 - 2017-10-11 - Lonardoni - wine-bottle fixed
!  9 - 2017-10-23 - Lonardoni - reduced output pot table
! 10 - 2023-09-20 - Tropiano  - increased nmax and lmax for Pb208
!-------------------------------------------------------------------------------
! The employed potential is given by:
! - Wood-Saxon + Wine-Bottle:
!     vc = vws * ( 1/(1+exp((r-rws)/aws)) + kwb*exp(-(r/awb)**2) )
! - Spin-Orbit (derivative of a Wood-Saxon like potential):
!     vso = vls / (2*m**2*r*als)*exp((r-rls)/als)/(1+exp((r-rls)/als))**2
! - Coulomb:
!     r > rem  vem = (Z-1)*hc*alpha/r
!     r < rem  vem = (Z-1)*0.5*hc*alpha*(3-r**2/rem**2)/rem
!
!   where:
!     rws = r0ws*A^1/3
!     rls = r0ls*A^1/3
!     rem = r0em*A^1/3
!     vws = -vws
!     vls = -kls*v0
!
! The subroutine returns the orbital tables and corresponding orbital names for
! a given set of tz,n,l,j
!     ntau = type of orbitals: 1 - nucleons with no Coulomb
!                              2 - distinguish protons and neutrons
!        A = number of nucleons
!        Z = number of protons
!     ntab = number of points for orbital tables
!     rmax = maximum r for orbital tables
!      prm = Wood-Saxon parameters prm(2,9): prm(1,:) - protons
!                                            prm(2,:) - neutrons
!     rdiv = divide orbital by r             true: get R(r); false: get u(r)=r R(r)
!      orb = orbital tables
!  orbname = orbital name tables
!     prnt = print summary, potentials and densities
!  prntorb = print orbital tables
!
! Note: element prm(2,9) not used
!       column  prm(2,:) ignored if ntau = 1
!
! If prnt = .true. the following output files are produced:
! - ws_log: summary with the list of bound states ordered with increasing energy
! - ws_pot: potential tables
! - ws_rho: density tables
! - p(n).nXX.lXX.jXX.orb: orbital tables
!
! Parameters (according to the references above):
!  v0 = 49.6                 als = 0.7
!  aws = 0.7                 r0ls p = 1.32
!  r0ws p = 1.275            r0ws n = 1.31
!  r0ws n = 1.347            kls p = 36
!  r0em   = 1.275            kls n = 35
!  vws (p & n) =
!     H3:  96.3168  51.8011
!    He3:  51.8011  96.3168
!    He4:  76.8412  76.8412
!    He6:  79.6234  46.2366
!    Li6:  68.4945  68.4945
!    O16:  58.0611  58.0611
!   Ca40:  54.3051  54.3051
!   Ca48:  59.4522  46.9322
! kwb = 0.0
! awb = 1.0
!-------------------------------------------------------------------------------

subroutine ws(ntau,A,Z,ntab,rmax,nrad,orbws,norb,lorb,jorb,prm,rdiv,prnt,prntorb,dens)
   implicit none
   integer, parameter :: i4=selected_int_kind(9)
   integer, parameter :: r8=selected_real_kind(15,9)
!-------------------------------------------------------------------------------
! parameters according to Argonne potential
!-------------------------------------------------------------------------------
   real(kind=r8), parameter :: mp=938.27231_r8     ! MeV
   real(kind=r8), parameter :: mn=939.56563_r8     ! MeV
   real(kind=r8), parameter :: mnuc=0.5_r8*(mp+mn) ! MeV
   real(kind=r8), parameter :: amu=931.494095_r8   ! MeV
   real(kind=r8), parameter :: mpi=138.0363_r8     ! MeV
   real(kind=r8), parameter :: hc=197.327053_r8    ! MeV * fm
   real(kind=r8), parameter :: alpha=1.0_r8/137.03599_r8
   real(kind=r8), parameter :: pi=acos(-1.0_r8)
!-------------------------------------------------------------------------------
!                              
! parameters Numerov
!-------------------------------------------------------------------------------
   integer(kind=i4), parameter :: nstep=12000
   integer(kind=i4), parameter :: nmatch=nstep/1000
   integer(kind=i4), parameter :: enstep=400
   integer(kind=i4), parameter :: lmax=6  ! good up to Pb208
   integer(kind=i4), parameter :: nmax=3  ! good up to Pb208
   real(kind=r8),    parameter :: emax=0.0_r8  ! Emax=0 if bound states only
   real(kind=r8),    parameter :: emin=-100.0_r8
   real(kind=r8),    parameter :: elim=0.00001_r8
   real(kind=r8),    parameter :: small=1.0e-10
!-------------------------------------------------------------------------------
   integer(kind=i4) :: ntau,nrad,A,Z,ntab,itau,l,ei,i,k,ir,isz,nn,n,tfn
   integer(kind=i4) :: red,bound(2),npart(2),deg
   integer(kind=i4) :: norb(nrad),lorb(nrad),jorb(nrad)
   integer(kind=i4), allocatable :: nsrt(:,:),lsrt(:,:),jsrt(:,:)
   real(kind=r8) :: rmax,prm(2,9)
   real(kind=r8) :: eup,edn,etrial
   real(kind=r8) :: hr,hrr,sz,j,ls,norm,match,nprot,nneut,rprot,rneut
   real(kind=r8) :: vempn(2),mass(2),mcore(2)
   real(kind=r8), allocatable :: v0(:),vws(:),rws(:),aws(:),kwb(:),awb(:)
   real(kind=r8), allocatable :: kls(:),rls(:),als(:),rem(:)
   real(kind=r8), allocatable :: vc(:,:),vso(:,:),vem(:,:)
   real(kind=r8), allocatable :: r(:),ri(:),kr(:),u(:),ul(:),ur(:)
   real(kind=r8), allocatable :: rho(:,:),elm1jms(:),elm1jps(:)
   real(kind=r8), allocatable :: el(:,:,:,:),eunsrt(:,:,:,:),esrt(:,:)
   real(kind=r8), allocatable :: orb(:,:,:,:,:)
   real(kind=r8) :: orbws(2,nrad,ntab), rr(ntab)  ! here 
   character(len=17), allocatable :: orbname(:,:,:,:)
   character(len=17) :: line
   logical :: prnt,prntorb,rdiv, dens

!-------------------------------------------------------------------------------
! inizialize
!-------------------------------------------------------------------------------
   if (nstep .le. ntab .or. mod(nstep,ntab) .ne. 0) then
      write(6,*) "Error in Wood Saxon, ntab must be smaller than ntab and nstep/ntab must be an integer"
      write(6,*) "ntab, nstep=", ntab, nstep
      stop
   endif

   do itau=1,ntau
      do k=1,9
         if (prm(itau,k).lt.0.0_r8.and.k.ne.4) prm(itau,k)=abs(prm(itau,k)) ! all parameters but #4 must be positive
      enddo
   enddo
   npart(1)=Z              ! proton number
   npart(2)=A-Z            ! neutron number 
   hr=rmax/real(nstep,r8)  ! step for internal grid
   hrr=rmax/real(ntab,r8)  ! step for output grid
   do i=1,ntab
      rr(i) = i * hrr
   enddo      
   red=int(hrr/hr, i4)                      !nstep/ntab          ! reduced grid factor  
 
   mcore(:)=(A-1)*amu      ! A-1 core masse
   if (ntau.eq.1) then     ! proton .eq. neutron
      mass(:)=1.0_r8/(1.0_r8/mnuc+1.0_r8/mcore(:))  ! reduced masses
   else                    ! proton .ne. neutron
      mass(1)=1.0_r8/(1.0_r8/mp+1.0_r8/mcore(1))    ! reduced masses
      mass(2)=1.0_r8/(1.0_r8/mn+1.0_r8/mcore(2))
   endif
   allocate(v0(ntau),vws(ntau),rws(ntau),aws(ntau),kwb(ntau),awb(ntau))
   allocate(kls(ntau),rls(ntau),als(ntau),rem(ntau))
   allocate(vc(ntau,nstep),vso(ntau,nstep),vem(ntau,nstep))
   allocate(r(nstep),ri(nstep),kr(nstep),u(nstep),elm1jms(0:lmax),elm1jps(0:lmax))
   allocate(el(ntau,0:nmax,0:lmax,2*lmax+1),eunsrt(ntau,0:nmax,0:lmax,2*lmax+1))
   allocate(ul(nstep),ur(nstep),rho(2,ntab))
   if (allocated(orb)) deallocate(orb)
   if (allocated(orbname)) deallocate(orbname)
   allocate(orb(2,0:nmax,0:lmax,2*lmax+1,ntab))
   allocate(orbname(2,0:nmax,0:lmax,2*lmax+1))
   v0=49.6_r8
   do itau=1,ntau
      vws(itau)=prm(itau,1)  ! Wood-Saxon: vws 
      rws(itau)=prm(itau,2)  ! Wood-Saxon: r0ws
      aws(itau)=prm(itau,3)  ! Wood-Saxon: aws
      kwb(itau)=prm(itau,4)  ! Wine-Bottle: kwb
      awb(itau)=prm(itau,5)  ! Wine-Bottle: awb
      kls(itau)=prm(itau,6)  ! Spin-Orbit: kls 
      rls(itau)=prm(itau,7)  ! Spin-Orbit: r0ls
      als(itau)=prm(itau,8)  ! Spin-Orbit: als
      rem(itau)=prm(itau,9)  ! Coulomb: r0em
      aws(itau)=max(aws(itau),small)
      awb(itau)=max(awb(itau),small)
      als(itau)=max(als(itau),small)
   enddo
   if (prnt) then
      open(unit=20,file='ws_log',status='unknown',access='append')
      rewind 20
      write(20,'(''Wood-Saxon orbitals: summary'',/)')
      if (ntau.eq.1) then
         if (kls(1).eq.0.0_r8) write(20,'(''Wood-Saxon'')')
         if (kls(1).ne.0.0_r8) write(20,'(''Wood-Saxon + Spin-Orbit'')')
         if (kwb(1).eq.0.0_r8) write(20,'(''no Wine-Bottle'')')
         write(20,'(/,8x,''prot & neut'')')
      else
         if ((kls(1)+kls(2)).eq.0.0_r8.and.rem(1).eq.0.0_r8) write(20,'(''Wood-Saxon'')')
         if ((kls(1)+kls(2)).eq.0.0_r8.and.rem(1).ne.0.0_r8) write(20,'(''Wood-Saxon + Coulomb'')')
         if ((kls(1)+kls(2)).ne.0.0_r8.and.rem(1).eq.0.0_r8) write(20,'(''Wood-Saxon + Spin-Orbit'')')
         if ((kls(1)+kls(2)).ne.0.0_r8.and.rem(1).ne.0.0_r8) write(20,'(''Wood-Saxon + Spin-Orbit + Coulomb'')')
         if ((kwb(1)+kwb(2)).eq.0.0_r8) write(20,'(''no Wine-Bottle'')')
         write(20,'(/,12x,''protons'',7x,''neutrons'')')
      endif
      write(20,'(''v0   ='',t10,f10.6,5x,f10.6)') (v0(i), i=1,ntau) 
      write(20,'(''vws  ='',t10,f10.6,5x,f10.6)') (vws(i),i=1,ntau) 
      write(20,'(''r0ws ='',t10,f10.6,5x,f10.6)') (rws(i),i=1,ntau)
      write(20,'(''rws  ='',t10,f10.6,5x,f10.6)') (rws(i)*A**(1.0_r8/3.0_r8),i=1,ntau)
      write(20,'(''aws  ='',t10,f10.6,5x,f10.6)') (aws(i),i=1,ntau)
      write(20,'(''kwb  ='',t10,f10.6,5x,f10.6)') (kwb(i),i=1,ntau)
      write(20,'(''awb  ='',t10,f10.6,5x,f10.6)') (awb(i),i=1,ntau)
      write(20,'(''kls  ='',t10,f10.6,5x,f10.6)') (kls(i),i=1,ntau)
      write(20,'(''r0ls ='',t10,f10.6,5x,f10.6)') (rls(i),i=1,ntau)
      write(20,'(''rls  ='',t10,f10.6,5x,f10.6)') (rls(i)*A**(1.0_r8/3.0_r8),i=1,ntau)
      write(20,'(''als  ='',t10,f10.6,5x,f10.6)') (als(i),i=1,ntau)
      if (ntau.eq.2) then
         write(20,'(''r0em ='',t10,f10.6,5x,f10.6)') (rem(i),i=1,ntau)
         write(20,'(''rem  ='',t10,f10.6,5x,f10.6)') (rem(i)*A**(1.0_r8/3.0_r8),i=1,ntau)
      endif
      write(20,'(/,''rmax orbitals ='',f8.3)') rmax
      write(20,'(''ntab orbitals ='',i8)') ntab
      write(20,'(''divide by r =  '',l8/)') rdiv
      write(20,'(''protons  ='',i4)')   npart(1)
      write(20,'(''neutrons ='',i4,/)') npart(2)
      write(20,'(''proton  reduced mass ='',f12.6)')   mass(1)
      write(20,'(''neutron reduced mass ='',f12.6,/)') mass(2)
      open(unit=30,file='ws_pot',status='unknown',access='append')
      open(unit=40,file='ws_rho',status='unknown',access='append')
      write(30,'(''Wood-Saxon orbitals: potentials'')')
      if (ntau.eq.1) then
         write(30,'(4(6x,a5,6x))') '  r  ','  vc ',' vls ',' vem '
      else
         write(30,'(6(6x,a5,6x))') '  r  ','vc_p ','vls_p','vc_n ','vls_n',' vem '
      endif
      write(40,'(''Wood-Saxon orbitals: densities'')')
      write(40,'(3(6x,a5,6x))') '  r  ','rho_p','rho_n'
   endif
   rws=rws*A**(1.0_r8/3.0_r8)
   rls=rls*A**(1.0_r8/3.0_r8)
   rem=rem*A**(1.0_r8/3.0_r8)
   vws=-vws
   kls=-kls*v0

!-------------------------------------------------------------------------------
! define radial step, radial coordinate and potentials
!-------------------------------------------------------------------------------
   vc=0.0_r8
   vso=0.0_r8
   vem=0.0_r8
   vempn(1)=1.0_r8
   vempn(2)=0.0_r8
   ri=1.0_r8
   do i=1,nstep
      r(i)=i*hr
      if (rdiv) ri(i)=1.0_r8/r(i)
      do itau=1,ntau
         vc(itau,i)=vws(itau)*(1.0_r8/(1.0_r8+exp((r(i)-rws(itau))/aws(itau)))+kwb(itau)*exp(-(r(i)/awb(itau))**2))
         if (kls(itau).ne.0.0_r8) vso(itau,i)=kls(itau)/(2.0_r8*r(i)*als(itau)*(mass(itau)/hc)**2) &
            *exp((r(i)-rls(itau))/als(itau))/(1.0_r8+exp((r(i)-rls(itau))/als(itau)))**2
         if (rem(itau).ne.0.0_r8.and.ntau.eq.2) then
            if (r(i).gt.rem(itau)) then
               vem(itau,i)=vempn(itau)*real(Z-1,r8)*hc*alpha/r(i)
            else
               vem(itau,i)=vempn(itau)*real(Z-1,r8)*0.5_r8*hc*alpha*(3.0_r8-r(i)**2/rem(itau)**2)/rem(itau)
            endif
         endif
      enddo
      if (prnt.and.(mod(i,50).eq.0)) write(30,'(6(e15.8,2x))') r(i),(vc(k,i),vso(k,i),k=1,ntau),vem(1,i)
   enddo
   deallocate(v0,vws,rws,aws,kwb,awb,kls,rls,als,rem)

!-------------------------------------------------------------------------------
! Numerov algorithm
!-------------------------------------------------------------------------------
   bound=0
   el=0.0_r8
   orb=0.0_r8
   orbname='notbound'
   do itau=1,ntau  ! loop on isospin: 1=protons, 2=neutrons
      do l=0,lmax  ! loop on angular momentum l
         do isz=0,1  ! loop on spin projection sz
            sz=real(isz,r8)-0.5_r8
            j=real(l,r8)+sz 
            ls=0.5_r8*(j*(j+1.0_r8)-real(l,r8)*real(l+1,r8)-0.75_r8)
            if (j.ge.0.0_r8) then  ! check j > 0
               if (l.eq.0) then
                  edn=emin
               else
                  edn=min(elm1jms(l-1),elm1jps(l-1))  ! El >= El-1 (faster code)
               endif
               eup=emax
               tfn=0
               do n=0,nmax  ! loop on nodes n: the larger n, the larger the energy
                  do ei=1,enstep
                     etrial=(edn+eup)/2.0_r8  ! bisection on energies
                     do i=1,nstep
                        kr(i)=2.0_r8*mass(itau)/hc**2*(etrial-vc(itau,i)-vso(itau,i)*ls-vem(itau,i)) &
                             -real(l,r8)*(real(l,r8)+1.0_r8)/r(i)**2
                     enddo
                     nn=0
                     ul(1)=r(1)**(l+1)  ! left boundary conditions
                     ul(2)=r(2)**(l+1)
                     do i=2,nstep-1
                        ul(i+1)=(2.0_r8*(1.0_r8-5.0_r8/12.0_r8*hr**2*kr(i))*ul(i)-(1.0_r8+hr**2/12.0_r8*kr(i-1))*ul(i-1)) &
                               /(1.0_r8+hr**2/12.0_r8*kr(i+1))  ! left wave functions
                        if (ul(i-1)*ul(i).lt.0.0_r8) then  ! count nodes number
                           nn=nn+1
                        endif
                     enddo
                     u=ul
                     if (abs(eup-edn).le.elim) then  ! convergence: store energies for 0 nodes
                        if (n.eq.0.and.isz.eq.0) then
                           elm1jms(l)=etrial  ! j=l+s
                        elseif (n.eq.0.and.isz.eq.1) then
                           elm1jps(l)=etrial  ! j=l-s
                        endif
                        exit
                     endif
                     if (nn.gt.n) then  ! check on nodes number
                        eup=etrial      ! nn > n : reduce etrial
                     else
                        edn=etrial      ! nn < n : increase etrial
                     endif
                     if (abs(abs(etrial)-abs(emax)).lt.elim) then  ! etrial = emax: no solutions for given n
                        if (prnt) write(20,'(a20,1x,i2,2x,a15,2x,i2,2x,a9,2x,i2)') &
                           ' impossibile to have',n,'nodes with  l =',l,'and  2j =',nint(2.0_r8*j)
                        tfn=1  ! too few nodes: exit loop on n
                        exit      
                     endif
                  enddo  ! end loop ei
                  if (tfn.eq.1) exit
                  ur(nstep)=exp(-sqrt(-2.0_r8*mass(itau)/hc**2*etrial)*r(nstep))  ! right boundary conditions
                  ur(nstep-1)=exp(-sqrt(-2.0_r8*mass(itau)/hc**2*etrial)*r(nstep-1))
                  do i=nstep-1,2,-1
                     ur(i-1)=(2.0_r8*(1.0_r8-5.0_r8/12.0_r8*hr**2*kr(i))*ur(i)-(1.0_r8+hr**2/12.0_r8*kr(i+1))*ur(i+1)) &
                            /(1.0_r8+hr**2/12.0_r8*kr(i-1))  ! right wave functions
                  enddo
                  match=ul(nstep/nmatch)/ur(nstep/nmatch)  ! match right with left wave functions
                  ur=ur*match
                  do i=nstep/nmatch,nstep
                     u(i)=ur(i)
                  enddo
                  norm=0.0_r8  ! normalize wave function
                  do i=1,nstep
                     norm=norm+r(i)**2*(u(i)/r(i))**2
                  enddo
                  norm=hr*norm
                  u=u/sqrt(norm)
                  if (etrial.lt.0.0_r8) then  ! store bound states' energies and wave functions
                     if (ntau.eq.1)               write(line,'(  ''n'',i2,''.l'',i2,''.j'',i2,''.orb'')') n,l,nint(2.0_r8*j)
                     if (ntau.eq.2.and.itau.eq.1) write(line,'(''p.n'',i2,''.l'',i2,''.j'',i2,''.orb'')') n,l,nint(2.0_r8*j)
                     if (ntau.eq.2.and.itau.eq.2) write(line,'(''n.n'',i2,''.l'',i2,''.j'',i2,''.orb'')') n,l,nint(2.0_r8*j)
                     call stripspaces(line)
                     orbname(itau,n,l,nint(2.0_r8*j))=trim(line)
                     bound(itau)=bound(itau)+1
                     el(itau,n,l,nint(2.0_r8*j))=etrial
                     if (prntorb) open(unit=50,file=trim(line),status='unknown')
                     !write(1996,*), "HERE ", red
                     do i=red,nstep,red
                        ir=i/red
                        orb(itau,n,l,nint(2.0_r8*j),ir)=u(i)*ri(i)     
                        if (prntorb) write(50,'(1p,e15.8,4x,e15.8)') r(i),orb(itau,n,l,nint(2.0_r8*j),ir)
                     enddo
                     if (prntorb) close(50)
                  endif
                  if (prnt) then
                     write(20,'(60(''-''))')
                     write(20,'('' itau ='',i2,4x,''n = '',i2,4x,''l = '',i2,4x,''2j = '',i2)') itau,n,l,nint(2.0_r8*j)
                     write(20,'('' E ='',f15.8)') el(itau,n,l,nint(2.0_r8*j))
                  endif
                  eup=emax
               enddo  ! end loop n
            endif  ! end check j>0
         enddo  ! end loop sz
      enddo  ! end loop l
   enddo  ! end loop itau
   deallocate(r,vc,vso,vem,kr,u,ul,ur,elm1jms,elm1jps)
!-------------------------------------------------------------------------------
! sort ground state proton and neutron energies
!-------------------------------------------------------------------------------
   n=maxval(bound)
   allocate(nsrt(2,n),lsrt(2,n),jsrt(2,n),esrt(2,n))
   eunsrt=el
   esrt=0.0_r8
   orbws=0.0_r8
   do itau=1,ntau
      if (prnt) then
         if (ntau.eq.1) write(20,'(/,''prot & neut'')')
         if (ntau.eq.2.and.itau.eq.1) write(20,'(/,''protons'')')
         if (ntau.eq.2.and.itau.eq.2) write(20,'(/,''neutrons'')')
         write(20,'(4(a6,1x),a12,3x,a6)')'  itau',' n','l',' 2j','energy','deg'  
      endif
      i=1
      do while (i.le.bound(itau))
         do n=0,nmax
            do l=0,lmax
               do isz=0,1
                  sz=real(isz,r8)-0.5_r8
                  j=real(l,r8)+sz 
                  if (j.gt.0.0_r8) then
                     if (eunsrt(itau,n,l,nint(2.0_r8*j)).lt.esrt(itau,i)) then
                        esrt(itau,i)=eunsrt(itau,n,l,nint(2.0_r8*j))
                        nsrt(itau,i)=n
                        lsrt(itau,i)=l
                        jsrt(itau,i)=nint(2.0_r8*j)
                     endif
                  endif
               enddo
            enddo
         enddo
         eunsrt(itau,nsrt(itau,i),lsrt(itau,i),jsrt(itau,i))=0.0_r8
         if (i.le.nrad) then
            if (prnt) then
               write(6,*)'wood saxon ordering'
               write(6,*)'n=',nsrt(itau,i), norb(i)
               write(6,*)'l=',lsrt(itau,i), lorb(i)
               write(6,*)'j=',jsrt(itau,i), jorb(i)
            endif
!            orbws(itau,i,:) = orb(itau,nsrt(itau,i),lsrt(itau,i),jsrt(itau,i), :) / rr(:)**lsrt(itau,i)
         endif

         if (esrt(itau,i).eq.0.0_r8) exit
         if (prnt) write(20,'(4(i6,1x),1x,f12.6,1x,i6)') &
            itau,nsrt(itau,i),lsrt(itau,i),jsrt(itau,i),esrt(itau,i),jsrt(itau,i)+1
         i=i+1
      enddo
   enddo
   do i=1,nrad
      orbws(1,i,:) = orb(1,norb(i)-1,lorb(i),jorb(i), :) / rr(:)**lorb(i)
   enddo

   if (ntau.eq.1) then
      bound(2)=bound(1)
      orb(2,:,:,:,:)=orb(1,:,:,:,:)
      orbws(2,:,:)=orbws(1,:,:)
      orbname(2,:,:,:)=orbname(1,:,:,:)
      esrt(2,:)=esrt(1,:)
      nsrt(2,:)=nsrt(1,:)
      lsrt(2,:)=lsrt(1,:)
      jsrt(2,:)=jsrt(1,:)
   endif
   if (prnt) then 
      write(20,'(/,''proton  bound states ='',i4)') bound(1)
      write(20,'(  ''neutron bound states ='',i4)') bound(2)
   endif
   deallocate(el,eunsrt,esrt)
   if (.not. dens) then
      deallocate(ri,nsrt,lsrt,jsrt,rho)
      return
   endif

!-------------------------------------------------------------------------------
! calculate proton and neutron densities, radii and particle number
!-------------------------------------------------------------------------------
   nprot=0.0_r8
   nneut=0.0_r8
   rprot=0.0_r8 
   rneut=0.0_r8
   rho=0.0_r8
   ri=1.0_r8
   do ir=1,ntab
      if (.not.rdiv) ri(ir)=1.0_r8/(ir*hrr)
      do itau=1,2
         n=0
         i=1
         do while (n.lt.npart(itau))
            if (i.gt.bound(itau)) then
               if (prnt) write(20,'(/,''Warning: not enough bound states ('',i3,'' ) to fit '',i3,'' particles'')') &
                      bound(itau),npart(itau)
               goto 10
            endif
            deg=jsrt(itau,i)+1
            if (n+deg.lt.npart(itau)) then
               rho(itau,ir)=rho(itau,ir)+real(deg,r8)*0.25_r8/pi*(ri(ir)*orb(itau,nsrt(itau,i),lsrt(itau,i),jsrt(itau,i),ir))**2
            else
               rho(itau,ir)=rho(itau,ir)+real(npart(itau)-n,r8)*0.25_r8/pi*(ri(ir) &
                   *orb(itau,nsrt(itau,i),lsrt(itau,i),jsrt(itau,i),ir))**2
            endif
            n=n+deg
            i=i+1
            if (n.ge.npart(itau)) exit
         enddo
      enddo
      if (prnt) write(40,'(3(e15.8,2x))') rr(ir),rho(1,ir),rho(2,ir)
      nprot=nprot+rr(ir)**2*rho(1,ir)
      nneut=nneut+rr(ir)**2*rho(2,ir)
      rprot=rprot+rr(ir)**4*rho(1,ir)
      rneut=rneut+rr(ir)**4*rho(2,ir)
   enddo
   nprot=nprot*4.0_r8*pi*hrr
   nneut=nneut*4.0_r8*pi*hrr
   rprot=rprot*4.0_r8*pi*hrr/real(nprot,r8)
   rneut=rneut*4.0_r8*pi*hrr/real(nneut,r8)
   if (prnt) then 
      write(20,'(/,''norm proton ='',f10.6,5x,''norm neutron ='',f10.6)') nprot,nneut
      write(20,'(  ''rms  proton ='',f10.6,5x,''rms  neutron ='',f10.6)') sqrt(rprot),sqrt(rneut)
   endif
10 if (prnt) then
      write(20,*)
      write(30,*)
      write(40,*)
      close(20)
      close(30)
      close(40)
   endif
   deallocate(ri,nsrt,lsrt,jsrt,rho)

!   do ir=1,ntab
!      write(2008,*) rr(ir), orbws(1,1,ir)
!   enddo
!   flush(2008)
   return


end subroutine ws

!-------------------------------------------------------------------------------
! remove empty spaces from string 
!-------------------------------------------------------------------------------
subroutine stripspaces(str)
   implicit none
   integer, parameter :: i4=selected_int_kind(9)
   integer, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4) :: strlen,last,actual
   character(len=*) :: str
   strlen=len(str)
   last=1
   actual=1
   do while (actual.lt.strlen)
      if (str(last:last).eq.' ') then
         actual=actual+1
         str(last:last)=str(actual:actual)
         str(actual:actual)=' '
      else
         last=last+1
         if (actual.lt.last) actual=last
      endif
   end do
   return
end subroutine stripspaces
