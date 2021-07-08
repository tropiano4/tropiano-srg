!
!=========================================================================
!
!
! To compile the code HFBRAD under linux using the Intel
! fortran compiler version 8.x, you must use the option "-pc80".
!
!
!=========================================================================
!
module cste
  !
  !! This module contains the constantes which are used throughout
  !! the code by all subroutines and functions.
  !
  implicit none
  public
  !
  ! pr = kind of real, you can switch it to extented precision
  !      by requiring more significant digits
  !
  integer, parameter :: pr = selected_real_kind( p = 12 )
  !
  !! The place where you want to store the output, potential
  !! and force files (this latter is unused in the present version
  !! of the code)
  !
  !character ( len = * ), parameter :: out = "out/", potdir = "pot/", &
  !     forces = "forces/"
  character ( len = * ), parameter :: out = "", potdir = "", forces = ""
  !
  !! Turn debug information on/off and memory allocation informations
  !
  logical :: debug = .false., memlog = .false.
  !
  !! Switch to turn the Coulomb interaction off
  !
  logical, parameter :: coulomb = .true.
  !
  !! Unit numbers for the output files are set here
  !
  integer, parameter :: ulog   = 10  ! log file
  integer, parameter :: uwave  = 12  ! wave functions, densities
  integer, parameter :: uspe   = 13  ! spectra
  integer, parameter :: uinput = 14  ! main input file
  integer, parameter :: usumup = 16  ! summary
  !
  !! Something small
  !
  real ( kind = pr ) :: small
  !
  !! Cubic root of the smallest real number of kind = pr
  !
  real ( kind = pr ) :: crtiny
  !
  !! qp = 4 pi
  !! pi = pi
  !! coef, echarg: the elementary charge used for the Coulomb
  !!               direct (echarg) and exchange (coef) fields.
  !!               Use the values of zero to switch the Coulomb
  !!               force off.
  !! maxa = log( biggest real of kind pr )
  !! eps = smallest real with  1 + eps /= 1
  !! t13 = 1/3
  !! t43 = 4/3
  !
  real ( kind = pr ) :: qp, coef, echarg, maxa, eps, pi, t13, t43
  !
  !! id2, id3 = 2x2 and 3x3 identity matrices
  !
  real ( kind = pr ), dimension(2,2), parameter :: &
       id2 = reshape( (/ 1.0_pr, 0.0_pr, 0.0_pr, 1.0_pr /), (/ 2, 2/) )
  real ( kind = pr ), dimension(3,3), parameter :: &
       id3 = reshape( (/ 1.0_pr, 0.0_pr, 0.0_pr, &
       &                 0.0_pr, 1.0_pr, 0.0_pr, &
       &                 0.0_pr, 0.0_pr, 1.0_pr /), (/ 3, 3/) )
  !
  !! Seven points vector used to compute derivatives
  !
  real ( kind = pr ), dimension(:), allocatable :: deriv
  !
  !! time
  !
  real ( kind = pr ) :: time0, time_previous
  !
  !! For diagonalization
  !
  real ( kind = pr ) :: diaeps
  real ( kind = pr ) :: diatol

  !......................................................................

contains

  !
  !!
  !

  subroutine set_cste()
    implicit none
    !
    !! Set constantes
    !
    t13 = 1.0_pr / 3.0_pr
    t43 = 4.0_pr / 3.0_pr
    crtiny = 100.0_pr * tiny(1.0_pr)**t13
    qp = 16.0_pr * atan(1.0_pr)
    pi =  4.0_pr * atan(1.0_pr)
    if (coulomb) then
       !echarg = 1.44_pr    ! physical constant
       !echarg = 1.43986_pr ! fluctuations
       echarg = 1.439978_pr
       coef = - 0.75_pr * ( 3 / pi )**t13 * echarg
    else
       coef = 0.0_pr
       echarg = 0.0_pr
       print '(16x,47("*"),/,16x,a,/,16x,47("*"))', &
            '************** NO COULOMB FIELDS **************'
    end if
    maxa = log( huge( 1.0_pr ) )
    eps = spacing( 1.0_pr )
    small = 10 * 10.0_pr**int( log ( 10000.0_pr * eps ) / log(10.0_pr) )
    !
    !! For diagonalization
    !
    diaeps = 40000 * spacing(1.0_pr)
    diatol = spacing(1.0_pr)**2.3_pr
    !
  end subroutine set_cste

  !
  !!
  !

end module cste
!
!=========================================================================
!
!
!=========================================================================
!
module param
  !
  !! Various parameters and global variables used throughout the code
  !
  use cste
  implicit none
  public
  save
  !
  !! Name of the input file
  !
  character ( len = 128 ) :: hfb_input
  !
  !! Skyrme force
  !
  character ( len = 4 ) :: force
  !
  !! Geometry
  !
  real ( kind = pr ) :: h     ! Step of integration.
  real ( kind = pr ) :: h12   ! Definition for Numerov
  real ( kind = pr ) :: h_12  ! For derivations
  real ( kind = pr ) :: hh_12 !      ''
  real ( kind = pr ) :: h_60  !      ''
  real ( kind = pr ) :: h_120 !      ''
  real ( kind = pr ) :: rbox  ! size of the sperical box
  !
  !! Number of particles.
  !
  integer :: n_nuclei
  integer :: proton, neutron
  real ( kind = pr ) :: rnumpart(2)
  !
  !! cut off method
  !
  integer :: i_cut
  !
  !! Dripline is a flag used, when non 0, to compute the dripline
  !! position. It is raised when one particle number read from input
  !! is negative.
  !
  integer :: dripline
  !
  !! Extension for the output files
  !
  character ( len = 8 ) :: extn
  !
  !! r, r**2 and 1/(4*pi*r**2) on the mesh points
  !
  real ( kind = pr ), dimension(:), allocatable :: r, r2, qpr2
  !
  !! Energy scale
  !
  real ( kind = pr ) :: dmshb, hb
  !
  !! Gaps, Fermi energies, mass...
  !
  real ( kind = pr ) :: amb(2), del(2), mass
  !
  !! Global energies
  !
  real ( kind = pr )               :: field_energy
  real ( kind = pr )               :: coulomb_energy
  real ( kind = pr )               :: ex_coulomb_energy
  real ( kind = pr )               :: spin_orbit_energy
  real ( kind = pr )               :: rearangment_energy
  real ( kind = pr )               :: total_energy
  real ( kind = pr )               :: energy_per_nucleon
  real ( kind = pr ), dimension(3) :: kinetic_energy
  real ( kind = pr ), dimension(3) :: pairing_energy
  real ( kind = pr ), dimension(3) :: pairing_kinetic_energy
  !
  !! Switches for the blocked states !.....................UNUSED
  !
  integer :: iing(2), ning(2), ling(2), jing(2)
  !
  !! Occupation probability of the blocked states !........UNUSED
  !
  real ( kind = pr ) :: vkblo(2)
  !
  !! Steering parameters
  !
  real ( kind = pr ) :: cut_diffuseness
  integer :: it_max
  integer :: pairing_force, boundary_condition, npr(2)
  logical :: bogolyubov(2)
  logical :: fixed_ph
  logical :: fixed_pp = .false.
  logical :: regularization
  integer :: j_max(2)
  !
  !! The pairing field is neglected fo r > r_cut
  !
  real ( kind = pr ) :: r_cut
  !
  !! Output control
  !
  logical :: quasiparticles, meanfields, densities, canonical_states
  logical :: cano
  !
  real ( kind = pr ) :: eps_energy ! accuracy for the converged total energy
  real ( kind = pr ) :: max_delta  !    "      "   "      "     mean gaps
  !
  ! .......................... Potentials
  !
  character ( len = 40 ) :: pot
  !
  ! .......................... ph mean field
  !
  real ( kind = pr ), dimension(:,:), allocatable :: v1f
  !
  ! .......................... ph kinetic
  !
  real ( kind = pr ), dimension(:,:), allocatable :: v2f
  !
  ! .......................... ph spin-orbit
  !
  real ( kind = pr ), dimension(:,:), allocatable :: v3f
  !
  ! .......................... pp mean field
  !
  real ( kind = pr ), dimension(:,:), allocatable :: v1p
  !
  ! .......................... pp kinetic
  !
  real ( kind = pr ), dimension(:,:), allocatable :: v2p
  !
  ! .......................... pp spin-orbit
  !
  real ( kind = pr ), dimension(:,:), allocatable :: v3p
  !
  ! .......................... Coulomb field
  !
  real ( kind = pr ), dimension(:),   allocatable ::  vc
  !
  !! Effective mass
  !
  real ( kind = pr ), dimension(:,:), allocatable :: v12f, v22f, v12p, v22p
  !
  !! Density
  !
  real ( kind = pr ), dimension(:,:), allocatable :: density
  !
  !! Number of states found in a given (it,l,j) block,
  !! dimension is 2*(j_max+1)
  !
  integer, dimension(:), allocatable :: nfobl
  !
  !! Logical flag which is lowered when then "accelereted" method
  !! can be used to find the solution in a given block
  !
  logical, dimension(:), allocatable :: brutal
  !
  !! Wave functions ff(:,:,:) and there derivatives dff(:,:,:)
  !!   indice 1 = up / down
  !!   indice 2 = function index (in a given l,j,it block)
  !!   indice 3 = point
  !! Since it is not possible to determine how many functions will
  !! be computed, a number NALLOC of functions will be allocated,
  !! this number will be automatically enlarged if needed (by 250)
  !! as many times as necessary.
  !! This method may allocate unused memory, but it is easier and
  !! faster to handle than a linked list.
  !
  real ( kind = pr ), dimension(:,:,:), allocatable :: ff, dff
  !
  !! Canonical states
  !
  real ( kind = pr ), dimension(:,:,:), allocatable :: canwf(:,:)
  !
  integer :: nalloc = 500     ! Initial guess
  integer :: nalloc_add = 250 ! increment for nalloc
  !
  !! Arraies containing the quasi-particle properties, these arraies
  !! are handled the same way as ff
  !!
  !
  real ( kind = pr ), dimension(:), allocatable :: ehfb, qpv2, qprad, ehf, &
       dhf, vhf
  real ( kind = pr ), dimension(:,:), allocatable :: e_up, e_down, e_mid
  integer, dimension(:), allocatable :: numnodes
  !
  !! Array containing informations about the canonical states: occupation,
  !! energy (diagonal m.e.), etc...
  !
  real ( kind = pr ), dimension(:), allocatable :: v2can, ecan, deltacan
  real ( kind = pr ), dimension(:), allocatable :: mecan
  !
  !! Particles number fluctuations
  !
  real ( kind = pr ) :: dispersion(2)
  !
  !! Densities
  !
  real ( kind = pr ), dimension(:,:), allocatable :: rho, tau, cur, &
       rho_p, tau_p, cur_p, geff, rega, regb
  !
  !! cut_off = energy cutoff,
  !! estep0 = initial step in energy for the search of the solutions
  !!          recommended value: 0.7 for the box size r = h * npt = 20 fm.
  !! xmu = part of the old potential taken to determine the new
  !!       potential in each iteration,
  !!       recommended value: 0.7 for the forces with the mass effective 1.
  !! xpmu = idem for the pairing part of the field (I have separated these
  !!        two, because you can be more brutal with the pairing field).
  !
  real ( kind = pr ) :: cut_off, estep0, estep0min, xmu, xmu0, xpmu, xpmu0
  !......................................................................

  interface allocatefct
     module procedure allocatefct
  end interface

  interface reallocatefct
     module procedure reallocatefct
  end interface

  !......................................................................

contains

  !
  !! Allocate memory for the solutions
  !

  subroutine allocatefct( npt, dimnfobl )
    implicit none
    integer, intent(in) :: npt, dimnfobl
    integer :: ifail
    if ( .not. allocated(vc) ) then
       if ( memlog ) &
            write( ulog, '(" ALLOCATEFCT()")' )
       allocate( vc(npt), &
            v1f( npt, 2 ),    v2f( npt, 2 ),  v3f( npt, 2 ), &
            v1p( npt, 2 ),    v2p( npt, 2 ),  v3p( npt, 2 ), &
            v12f( npt, 2 ),   v22f( npt, 2 ), &
            v12p( npt, 2 ),   v22p( npt, 2 ), &
            density(npt,npt), &
            stat = ifail )
       if ( ifail /= 0 ) stop 'Not enough memory... (param/1)'
       allocate( ff(2,nalloc,1:npt), dff(2,nalloc,1:npt), &
            canwf(nalloc,1:npt),                          &
            ehfb(nalloc), qpv2(nalloc), qprad(nalloc),    &
            ehf(nalloc), dhf(nalloc), vhf(nalloc),        &
            v2can(nalloc), mecan(nalloc), ecan(nalloc), deltacan(nalloc), &
            e_up(nalloc,2), e_down(nalloc,2), e_mid(nalloc,2),     &
            rho(npt,2), tau(npt,2), cur(npt,2),                    &
            rho_p(npt,2), tau_p(npt,2), cur_p(npt,2),              &
            geff(npt,2), rega(npt,2), regb(npt,2),                 &
            numnodes(nalloc),                                      &
            stat = ifail )
       if ( ifail /= 0 ) stop 'Not enough memory...(param/2)'
    end if
    !
    !! The maximum value of j_max can be different for different systems
    !! so the following arraies have to be reallocated
    !
    if ( allocated(nfobl) ) deallocate( nfobl, brutal )
    if ( memlog ) write( ulog, '(" ALLOCATEFCT() for nfobl and brutal")' )
    allocate( nfobl(dimnfobl), brutal(dimnfobl), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...(param/3)'
    !
  end subroutine allocatefct

  !
  !! Deallocate memory for wave functions and fields
  !

  subroutine deallocatefct()
    implicit none
    !
    deallocate( vc, v1f, v2f, v3f, v1p, v2p, v3p, v12f, v22f, &
         v12p, v22p, density )
    deallocate( ff, dff, canwf, ehfb, qpv2, qprad, ehf, dhf, vhf,  &
         v2can, mecan, ecan, deltacan, e_up, e_down, e_mid,        &
         rho, tau, cur, rho_p, tau_p, cur_p, geff, rega, regb,     &
         numnodes )
    deallocate( nfobl, brutal )
    !
  end subroutine deallocatefct

  !
  !! In case of necessity, reallocate more memory
  !

  subroutine reallocatefct(npt)
    implicit none
    integer, intent(in) :: npt
    integer :: ifail, n
    real ( kind = pr ), dimension(:,:,:), allocatable :: tmp
    real ( kind = pr ), dimension(:), allocatable :: tmpe
    real ( kind = pr ), dimension(:,:), allocatable :: tmp2
    integer, dimension(:), allocatable :: itmp
    !
    if ( memlog ) write( ulog, '(" REALLOCATEFCT()")' )
    if ( memlog ) then
       write( ulog, &
            '(" ** Reallocating memory for the wave functions...")' )
       write( ulog, &
            '(" ** ",i5,"    --->   ",i5)' ) nalloc, nalloc + nalloc_add
    end if
    n = nalloc
    nalloc = nalloc + nalloc_add
    allocate( tmp(2,n,1:npt), tmpe(n), tmp2(n,2), itmp(n), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    !
    itmp = numnodes
    deallocate(numnodes)
    allocate( numnodes(nalloc), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    numnodes(1:n) = itmp
    numnodes(n+1:nalloc) = 0
    !
    tmp = ff
    deallocate(ff)
    allocate( ff(2,nalloc,1:npt), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    ff(:,1:n,:) = tmp(:,:,:)
    ff(:,n+1:nalloc,:) = 0.0_pr
    !
    tmp = dff
    deallocate(dff)
    allocate( dff(2,nalloc,1:npt), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    dff(:,1:n,:) = tmp(:,:,:)
    dff(:,n+1:nalloc,:) = 0.0_pr
    !
    tmp(1,:,:) = canwf(:,:)
    deallocate(canwf)
    allocate( canwf(nalloc,1:npt), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    canwf(1:n,:) = tmp(1,:,:)
    canwf(n+1:nalloc,:) = 0.0_pr
    !
    tmpe = ehfb
    deallocate(ehfb)
    allocate( ehfb(nalloc), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    ehfb = 0.0_pr
    ehfb(1:n) = tmpe(:)
    ehfb(n+1:nalloc) = 0.0_pr
    !
    tmp2 = e_up
    deallocate(e_up)
    allocate( e_up(nalloc,2), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    e_up(1:n,:) = tmp2(:,:)
    e_up(n+1:nalloc,:) = 0.0_pr
    !
    tmp2 = e_down
    deallocate(e_down)
    allocate( e_down(nalloc,2), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    e_down(1:n,:) = tmp2(:,:)
    e_down(n+1:nalloc,:) = 0.0_pr
    !
    tmp2 = e_mid
    deallocate(e_mid)
    allocate( e_mid(nalloc,2), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    e_mid(1:n,:) = tmp2(:,:)
    e_mid(n+1:nalloc,:) = 0.0_pr
    !
    tmpe = qpv2
    deallocate(qpv2)
    allocate( qpv2(nalloc), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    qpv2(1:n) = tmpe(:)
    qpv2(n+1:nalloc) = 0.0_pr
    !
    tmpe = qprad
    deallocate(qprad)
    allocate( qprad(nalloc), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    qprad(1:n) = tmpe(:)
    qprad(n+1:nalloc) = 0.0_pr
    !
    tmpe = ehf
    deallocate(ehf)
    allocate( ehf(nalloc), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    ehf(1:n) = tmpe(:)
    ehf(n+1:nalloc) = 0.0_pr
    !
    tmpe = dhf
    deallocate(dhf)
    allocate( dhf(nalloc), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    dhf(1:n) = tmpe(:)
    dhf(n+1:nalloc) = 0.0_pr
    !
    tmpe = vhf
    deallocate(vhf)
    allocate( vhf(nalloc), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    vhf(1:n) = tmpe(:)
    vhf(n+1:nalloc) = 0.0_pr
    !
    tmpe = v2can
    deallocate(v2can)
    allocate( v2can(nalloc), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    v2can(1:n) = tmpe(:)
    v2can(n+1:nalloc) = 0.0_pr
    !
    tmpe = ecan
    deallocate(ecan)
    allocate( ecan(nalloc), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    ecan(1:n) = tmpe(:)
    ecan(n+1:nalloc) = 0.0_pr
    !
    tmpe = mecan
    deallocate(mecan)
    allocate( mecan(nalloc), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    mecan(1:n) = tmpe(:)
    mecan(n+1:nalloc) = 0.0_pr
    !
    tmpe = deltacan
    deallocate(deltacan)
    allocate( deltacan(nalloc), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory...'
    deltacan(1:n) = tmpe(:)
    deltacan(n+1:nalloc) = 0.0_pr
    !
    deallocate( tmp, tmpe, tmp2, itmp )
    if ( memlog ) write( ulog, '("    Done.")' )
    !
  end subroutine reallocatefct

  !
  !! It's probably useless to flush the allocated array
  !! but it does not hurt.
  !

  subroutine flush()
    implicit none
    !
    if ( memlog ) write( ulog, '(" FLUSH()")' )
    v1f      = 0.0_pr
    v2f      = 0.0_pr
    v3f      = 0.0_pr
    v1p      = 0.0_pr
    v2p      = 0.0_pr
    v3p      = 0.0_pr
    vc       = 0.0_pr
    v12f     = 0.0_pr
    v22f     = 0.0_pr
    v12p     = 0.0_pr
    v22p     = 0.0_pr
    density  = 0.0_pr
    ehfb     = 0.0_pr
    qpv2     = 0.0_pr
    qprad    = 0.0_pr
    ehf      = 0.0_pr
    dhf      = 0.0_pr
    vhf      = 0.0_pr
    e_up     = 0.0_pr
    e_down   = 0.0_pr
    e_mid    = 0.0_pr
    v2can    = 0.0_pr
    ecan     = 0.0_pr
    deltacan = 0.0_pr
    mecan    = 0.0_pr
    rho      = 0.0_pr
    tau      = 0.0_pr
    cur      = 0.0_pr
    rho_p    = 0.0_pr
    tau_p    = 0.0_pr
    cur_p    = 0.0_pr
    geff     = 0.0_pr
    rega     = 0.0_pr
    regb     = 0.0_pr
    ff       = 0.0_pr
    dff      = 0.0_pr
    numnodes = 0
    nfobl    = 0
    brutal   = .true.
    !
  end subroutine flush

  !
  !!
  !

end module param
!
!=========================================================================
!
!
!=========================================================================
!
module linalg
  !
  !! Several linear algebra routines
  !
  use cste
  implicit none

  private
  public :: matinv22, matinv33, matinv44, diagon

  !......................................................................

  interface matinv22
     module procedure matinv22
  end interface

  interface matinv33
     module procedure matinv33
  end interface

  interface matinv44
     module procedure matinv44
  end interface

  interface diagon
     module procedure diagon
  end interface

  !......................................................................

contains

  subroutine matinv22( a, ai )
    !
    !! Inverse of a 2x2 matrix:
    !!   a = matrix, ai = inverse (if not singular)
    !
    implicit none
    real ( kind = pr ), intent(in) :: a(2,2)
    real ( kind = pr ), intent(inout) :: ai(2,2)
    real ( kind = pr ) :: det
    det = a(1,1) * a(2,2) - a(1,2) * a(2,1)
    if ( det == 0.0_pr ) stop "Singular 2x2 matrix..."
    ai(1,1) = a(2,2)
    ai(2,2) = a(1,1)
    ai(1,2) = - a(1,2)
    ai(2,1) = - a(2,1)
    ai = ai / det
    !
  end subroutine matinv22

  !
  !!
  !

  subroutine matinv33( a, ai )
    !
    !! Inverse of a 3x3 matrix:
    !!   a = matrix, ai = inverse (if not singular).
    !!   If the determinant is zero, the subroutine sets
    !!   the inverse to something big.
    !
    implicit none
    real ( kind = pr ), intent(in) :: a(3,3)
    real ( kind = pr ), intent(inout) :: ai(3,3)
    real ( kind = pr ) :: det
    !
    det =  a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
         - a(1,2) * ( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) &
         + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )
    !det =  a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
    !     + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) &
    !     - a(1,2) * ( a(2,1) * a(3,3) - a(2,3) * a(3,1) )
    if ( det == 0.0_pr ) then
       !
       !! Sometimes you can have a singular matrix accidentaly,
       !! don't worry too much if it happens only once in a
       !! while during the iterations.
       !
       print '("Singular 3x3 matrix...")'
!print *, a
!print *
!print *, a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) )
!print *, - a(1,2) * ( a(2,1) * a(3,3) - a(2,3) * a(3,1) )
!print *, a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )
!print *
!print *, a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
!- a(1,2) * ( a(2,1) * a(3,3) - a(2,3) * a(3,1) )
!print *, a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )
!print *
!print *, det
!print *, ( a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
!+ a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) ) &
!- a(1,2) * ( a(2,1) * a(3,3) - a(2,3) * a(3,1) )
       ai = 1.e10_pr
    else
       ai(1,1) = a(2,2) * a(3,3) - a(2,3) * a(3,2)
       ai(2,1) = - ( a(2,1) * a(3,3) - a(2,3) * a(3,1) )
       ai(3,1) = a(2,1) * a(3,2) - a(2,2) * a(3,1)
       ai(1,2) = - ( a(1,2) * a(3,3) - a(1,3) * a(3,2) )
       ai(2,2) = a(1,1) * a(3,3) - a(1,3) * a(3,1)
       ai(3,2) = - ( a(1,1) * a(3,2) - a(1,2) * a(3,1) )
       ai(1,3) = a(1,2) * a(2,3) - a(1,3) * a(2,2)
       ai(2,3) = - ( a(1,1) * a(2,3) - a(1,3) * a(2,1) )
       ai(3,3) = a(1,1) * a(2,2) - a(1,2) * a(2,1)
       ai = ai / det
    end if
    !
  end subroutine matinv33

  !
  !! Warning: this subroutine does not return the inverse !
  !

  subroutine matinv44( a, ai, det )
    !
    !! Inverse of a 4x4 matrix
    !!   a = matrix, ai = inverse * determinant (if not singular)
    !!                              -----------
    !!   det = determinant
    !
    implicit none
    real ( kind = pr ), intent(in) :: a(4,4)
    real ( kind = pr ), intent(inout) :: ai(4,4)
    real ( kind = pr ), intent(out) :: det
    integer :: s, j, k, j1, j2, j3, k1, k2, k3
    !
    det = 0.0_pr
    s = 1
    j1 = 2
    j2 = 3
    j3 = 4
    do j = 1, 4
       if ( j == 2 ) j1 = 1
       if ( j == 3 ) j2 = 2
       if ( j == 4 ) j3 = 3
       k1 = 2
       k2 = 3
       k3 = 4
       do k = 1, 4
          if ( k == 2 ) k1 = 1
          if ( k == 3 ) k2 = 2
          if ( k == 4 ) k3 = 3
          ai(k,j)= a(j1,k1) &
               * ( a(j2,k2) * a(j3,k3) - a(j2,k3) * a(j3,k2) ) &
               - a(j1,k2) &
               * ( a(j2,k1) * a(j3,k3) - a(j2,k3) * a(j3,k1) ) &
               + a(j1,k3) &
               * ( a(j2,k1) * a(j3,k2) - a(j2,k2) * a(j3,k1) )
          ai(k,j) = ai(k,j) * s
          s = - s
       end do
       s = - s
       det = det + a(j,1) * ai(1,j)
    end do
    !
  end subroutine matinv44

  !-------------------------------------------------------------------------
  !
  !! Diagonalization
  !
  !-------------------------------------------------------------------------

  !
  !! Subroutines to swap scalars or vectors
  !

  subroutine swaps( a, b )
    implicit none
    real ( kind = pr ), intent(inout) :: a, b
    real ( kind = pr ) :: w
    !
    w = a
    a = b
    b = w
    !
  end subroutine swaps

  subroutine swapv( a, b, n )
    implicit none
    integer, intent(in) :: n
    real ( kind = pr ), intent(inout) :: a(n), b(n)
    real ( kind = pr ) :: w(n)
    !
    w = a
    a = b
    b = w
    !
  end subroutine swapv

  !
  !! Matrix digonalization...
  !

  subroutine diagon( n, a, d, z )
    implicit none
    integer, intent(in) :: n
    real ( kind = pr ), dimension(n,n), intent(in) :: a
    real ( kind = pr ), dimension(n,n), intent(inout) :: z
    real ( kind = pr ), dimension(n), intent(inout) :: d
    real ( kind = pr ), dimension(n) :: e
    integer :: i, l, j, k, j1, m
    real ( kind = pr ) :: h, scale, f, g, hh, b, p
    real ( kind = pr ) :: q, c, s
    !
    z = a
    if ( n <= 1 ) return
    do i = n, 2, -1
       l = i - 1
       h = 0.0_pr
       scale = 0.0_pr
       if ( l /= 1 ) scale = sum( abs( z(i,1:l) ) )
       if ( scale <= diatol ) then
          e(i) = z(i,l)
          d(i) = h
          cycle
       end if
       z(i,1:l) = z(i,1:l) / scale
       h = dot_product( z(i,1:l), z(i,1:l) )
       f = z(i,l)
       g = - sign( sqrt(h), f )
       e(i) = g * scale
       h = h - f * g
       z(i,l) = f - g
       f = 0.0_pr
       do j = 1, l
          z(j,i) = z(i,j) / ( h * scale )
          g = dot_product( z(j,1:j), z(i,1:j) )
          j1 = j + 1
          if ( j1  <= l ) g = g + dot_product( z(j1:l,j), z(i,j1:l) )
          e(j) = g / h
          f = f + e(j) * z(i,j)
       end do
       hh = f / ( h + h )
       do j = 1, l
          f = z(i,j)
          g = e(j) - hh * f
          e(j) = g
          z(j,1:j) = z(j,1:j) - f * e(1:j) - g * z(I,1:j)
       end do
       z(i,1:l) = scale * z(i,1:l)
       d(i) = h
    end do
    d(1) = 0.0_pr
    e(1) = 0.0_pr
    do i = 1, n
       l = i - 1
       if ( abs(d(i)) >= 1.e-16_pr .and. l /= 0 ) then
          do j = 1, l
             g = dot_product( z(i,1:l), z(1:l,j) )
             z(1:l,j) = z(1:l,j) - g * z(1:l,i)
          end do
       end if
       d(i) = z(i,i)
       z(i,i) = 1.0_pr
       if ( l == 0 ) cycle
       z(i,1:l) = 0.0_pr
       z(1:l,i) = 0.0_pr
    end do
    !
    e = eoshift( e, 1 )
    b = 0.0_pr
    f = 0.0_pr
    do l = 1, n
       j = 0
       h = diaeps * ( abs(d(l)) + abs(e(l)) )
       if ( b < h ) b = h
       do m = l, n
          if ( abs(e(m)) - b <= 0.0_pr ) exit
       end do
       if ( m /= l ) then
          do
             if ( j == 30 ) stop '1111'
             j = j + 1
             p = ( d(l+1) - d(l) ) / ( 2 * e(l) )
             q = sqrt( p * p + 1.0_pr )
             h = d(l) - e(l) / ( p + sign( q, p ) )
             d(l:n) = d(l:n) - h
             f = f + h
             p = d(m)
             c = 1.0_pr
             s = 0.0_pr
             do i = m - 1, l, -1
                g = c * e(i)
                h = c * p
                if ( abs(p) - abs(e(i)) >= 0.0_pr ) then
                   c = e(i) / p
                   q = sqrt( c * c + 1.0_pr )
                   e(i+1) = s * p * q
                   s = c / q
                   c = 1.0_pr / q
                else
                   c = p / e(i)
                   q = sqrt( c * c + 1.0_pr )
                   e(i+1) = s * e(i) * q
                   s = 1.0_pr / q
                   c = c / q
                end if
                p = c * d(i) - s * g
                d(i+1) = h + s * ( c * g + s * d(i) )
                do k = 1, n
                   h = z(k,i+1)
                   z(k,i+1) = s * z(k,i) + c * h
                   z(k,i) = c * z(k,i) - s * h
                end do
             end do
             e(l) = s * p
             d(l) = c * p
             if ( abs(e(l)) - b <= 0.0_pr ) exit
          end do
       end if
       d(l) = d(l) + f
    end do
    do i = 1, n - 1
       k = i
       p = d(i)
       do j = i + 1, n
          if ( d(j) - p >= 0.0_pr ) cycle
          k = j
          p = d(j)
       end do
       if ( k == i ) cycle
       call swaps( d(k), d(i) )
       call swapv( z(1:n,i), z(1:n,k), n )
    end do

  end subroutine diagon

  !
  !!
  !

end module linalg
!
!=========================================================================
!

!
!=========================================================================
!
module eqdifstatic

  use cste
  real ( kind = pr ), dimension(:,:), allocatable :: yfu, yfd
  real ( kind = pr ), dimension(:,:), allocatable :: ybu, ybd
  real ( kind = pr ), dimension(:), allocatable :: deta
  real ( kind = pr ), dimension(:), allocatable :: d2f1, d2f2, dif1, dif2

  interface allocateeqdif
     module procedure allocateeqdif
  end interface

contains

  subroutine allocateeqdif(n)
    implicit none
    integer, intent(in) :: n
    integer :: ifail
    allocate( yfu(2,n), yfd(2,n), ybu(2,n), ybd(2,n), &
         deta(n), d2f1(n), d2f2(n), dif1(n), dif2(n), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory... (eqdif)'
  end subroutine allocateeqdif

  !
  !!
  !

  subroutine deallocateeqdif()
    implicit none
    deallocate( yfu, yfd, ybu, ybd, deta, d2f1, d2f2, dif1, dif2 )
  end subroutine deallocateeqdif

  !
  !!
  !

end module eqdifstatic

!
!=========================================================================
!

module eqdif

  use cste
  implicit none

  !......................................................................

  interface integrate
     module procedure integrate
  end interface

  interface derivative
     module procedure derivative_sca
     module procedure derivative_vec
  end interface

  interface d1_and_d2
     module procedure d1_and_d2
  end interface

  interface accuracy_check
     module procedure accuracy_check
  end interface

  interface count_nodes
     module procedure count_nodes
  end interface

  !......................................................................

contains

  !
  !!
  !

  subroutine d1_and_d2( f, df, d2f, n )
    !
    use param, only : h_12, hh_12
    implicit none
    !
    !! This subroutine computes the first and second derivative of
    !! of function evaluated on the meshpoints 1,...,npt.
    !! The input is the function f with extrapolated values in -1, 0.
    !
    integer, intent(in) :: n
    real ( kind = pr ), intent(in) :: f(-1:n)
    real ( kind = pr ), intent(inout) :: df(1:n)
    real ( kind = pr ), intent(inout), optional :: d2f(1:n)
    integer :: i
    !
    do i = 1, n - 2
       df(i) = ( 8 * ( f(i+1) - f(i-1) ) - f(i+2) + f(i-2) ) / h_12
       if ( present(d2f) ) &
            d2f(i) = ( - 30 * f(i) + 16 * ( f(i+1) + f(i-1) ) &
            - f(i+2) - f(i-2) ) / hh_12
    end do
    !
    df(n-1) = ( - f(n-4) + 6 * f(n-3) - 18 * f(n-2) &
         + 10 * f(n-1) + 3 * f(n) ) / h_12
    df(n) = ( 3 * f(n-4) - 16 * f(n-3) + 36 * f(n-2) &
         - 48 * f(n-1) + 25 * f(n) ) / h_12
    if ( present(d2f) ) then
       d2f(n-1) = ( - f(n-4) + 4 * f(n-3) + 6 * f(n-2) &
            - 20 * f(n-1) + 11 * f(n) ) / hh_12
       d2f(n) = ( 11 * f(n-4) - 56 * f(n-3) + 114 * f(n-2) &
            - 104 * f(n-1) + 35 * f(n) ) / hh_12
    end if
    !
  end subroutine d1_and_d2

  !
  !!
  !

  subroutine integrate( n0, l0, n, eh12d, h12v, h12w, det1, edet2,  &
       fall1, fall2, nmatch, l, det, f1, f2 )
    use param, only : deriv
    use eqdifstatic
    use linalg
    implicit none
    !
    !! This subroutine integrates the coupled differential equations:
    !!
    !!    ( -d2/dr2 + v ) * f1 + w * f2 = E/d * f1
    !!    w * f1 - ( -d2/dr2 + v ) * f2 = E/d * f2
    !!
    !! The result is det = determinant of the mathing matrix in nmatch.
    !! If the optional arguments f1 and f2 are present they will
    !! contain the solutions at the end.
    !
    integer, intent(in) :: n, nmatch, l, n0
    real ( kind = pr ), intent(inout), optional :: f1(n), f2(n)
    real ( kind = pr ), intent(out) :: det
    real ( kind = pr ), intent(in) :: fall1, fall2, l0
    real ( kind = pr ), dimension(n), intent(in) :: eh12d, h12v, h12w, &
         det1, edet2
    integer :: nm, idm, i, iii
    real ( kind = pr ) :: mat(4,4), matinv(4,4)
    real ( kind = pr ) :: g0(2), g1(2), a, b, c, d, mc
    real ( kind = pr ) :: cond(3), d3, d0
    real ( kind = pr ), dimension(3,3) :: m3, mi31, mi32, mi33, mi34
    real ( kind = pr ) :: detb, f11, f12, f21, f22, t11, t12, t21, t22, &
         t, g11, g21, x
    !
    !! Integration:
    !!  we integration four times with different conditions at
    !!  the origin and at the box radius.
    !
    deta = det1 - edet2
    !
    !! Forward integration
    !
    if ( l <= 1 .or. n0 == 1 ) then
       g0 = (/ l0, 0.0_pr /)
       g1 = (/ h12v(1) + eh12d(1), h12w(1) /)
       yfu(:,1) = (/ 1.0_pr, 0.0_pr /)
    else
       do i = 1, n0
          yfu(:,i) = (/ ( real( i, pr ) / ( n0 - 1 ) )**(l+1), 0.0_pr /)
       end do
       g0 = 0.0_pr
       g1 = (/ h12v(n0) + eh12d(n0), h12w(n0) /) * yfu(1,n0)
    end if
    nm = nmatch + 3
    call numerov( yfu, g0, g1, eh12d, h12v, h12w, deta, n0+1, nm, n )
    !
    if ( l <= 1 .or. n0 == 1 ) then
       g0 = (/ 0.0_pr, l0 /)
       g1 = (/ -h12w(1), h12v(1) - eh12d(1) /)
       yfd(:,1) = (/ 0.0_pr, 1.0_pr /)
    else
       do i = 1, n0
          yfd(:,i) = (/ 0.0_pr, ( real( i, pr ) / ( n0 - 1 ) )**(l+1) /)
       end do
       g0 = 0.0_pr
       g1 = (/ - h12w(n0), h12v(n0) - eh12d(n0) /) * yfd(2,n0)
    end if
    call numerov( yfd, g0, g1, eh12d, h12v, h12w, deta, n0+1, nm, n )
    !
    !! Backward integration
    !
    nm = nmatch - 3
    ybu(:,n) = (/ 1.0_pr, 0.0_pr /)
    g1 = (/ h12v(n) + eh12d(n), h12w(n) /)
    g0 = g1 * fall1
    call numerov( ybu, g0, g1, eh12d, h12v, h12w, deta, n-1, nm, n )
    !
    ybd(:,n) = (/ 0.0_pr, 1.0_pr /)
    g1 = (/ -h12w(n), h12v(n) - eh12d(n) /)
    g0 = g1 * fall2
    call numerov( ybd, g0, g1, eh12d, h12v, h12w, deta, n-1, nm, n )
    !
    !! Matching matrix
    !
    nm = nmatch
    mat(1:2,1) = yfu(1:2,nm)
    mat(1:2,2) = yfd(1:2,nm)
    mat(1:2,3) = - ybu(1:2,nm)
    mat(1:2,4) = - ybd(1:2,nm)
    !
    mat(3,1) = dot_product( deriv(-3:3), yfu(1,nm-3:nm+3) )
    mat(3,2) = dot_product( deriv(-3:3), yfd(1,nm-3:nm+3) )
    mat(3,3) = - dot_product( deriv(-3:3), ybu(1,nm-3:nm+3) )
    mat(3,4) = - dot_product( deriv(-3:3), ybd(1,nm-3:nm+3) )
    !
    mat(4,1) = dot_product( deriv(-3:3), yfu(2,nm-3:nm+3) )
    mat(4,2) = dot_product( deriv(-3:3), yfd(2,nm-3:nm+3) )
    mat(4,3) = - dot_product( deriv(-3:3), ybu(2,nm-3:nm+3) )
    mat(4,4) = - dot_product( deriv(-3:3), ybd(2,nm-3:nm+3) )
    !
    !! Construction of the inverse
    !
    nm = nmatch
    call matinv44( mat, matinv, det )
    !
    !! Build the solutions, if needed
    !
    if ( present(f1) .and. present(f2) ) then
       !
       !! Test of accuracy at the matching point
       !
       detb = mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1)
       f11 = ( mat(1,2) * mat(2,3) - mat(2,2) * mat(1,3) ) / detb
       f12 = ( mat(1,2) * mat(2,4) - mat(2,2) * mat(1,4) ) / detb
       f21 = ( mat(2,1) * mat(1,3) - mat(1,1) * mat(2,3) ) / detb
       f22 = ( mat(2,1) * mat(1,4) - mat(1,1) * mat(2,4) ) / detb
       t11 = - mat(3,3) - mat(3,1) * f11 - mat(3,2) * f21
       t12 = - mat(3,4) - mat(3,1) * f12 - mat(3,2) * f22
       t21 = - mat(4,3) - mat(4,1) * f11 - mat(4,2) * f21
       t22 = - mat(4,4) - mat(4,1) * f12 - mat(4,2) * f22
       t = t11 * t22 - t21 * t12
       if ( abs(t11) + abs(t12) < abs(t22) + abs(t21) ) then
          g11 =   t22
          g21 = - t21
          x = mat(1,3) * g11 + mat(1,4) * g21
       else
          g11 = - t12
          g21 =   t11
          x = mat(2,3) * g11 + mat(2,4) * g21
       end if
       if ( x == 0.0_pr ) print '(" WRONG MATCHING:  l = ",i2)', l
       !
       idm = 1
       m3 = mat(1:3,1:3)
       call matinv33( m3, mi31 )
       d3 = maxval( abs( id3 - matmul( m3, mi31 ) ) )
       d0 = d3
       !
       m3(1,1) = mat(1,1)
       m3(2:3,1) = mat(3:4,1)
       m3(1,2:3) = mat(1,3:4)
       m3(2:3,2:3) = mat(3:4,3:4)
       call matinv33( m3, mi32 )
       d3 = maxval( abs( id3 - matmul( m3, mi32 ) ) )
       if ( d3 < d0 ) then
          idm = 2
          d0 = d3
       end if
       !
       m3(1:2,1:2) = mat(1:2,1:2)
       m3(1:2,3) = mat(1:2,4)
       m3(3,1:2) = mat(4,1:2)
       m3(3,3) = mat(4,4)
       call matinv33( m3, mi33 )
       d3 = maxval( abs( id3 - matmul( m3, mi33 ) ) )
       if ( d3 < d0 ) then
          idm = 3
          d0 = d3
       end if
       !
       m3 = mat(2:4,2:4)
       call matinv33( m3, mi34 )
       d3 = maxval( abs( id3 - matmul( m3, mi34 ) ) )
       if ( d3 < d0 ) then
          idm = 4
          d0 = d3
       end if
       !
       !! Warning if the matching condition is not acurately fulfilled
       !
       if ( d0 > 1.e9 * eps ) then
          print &
               '(" Warning: matching not accurate  Err =",e12.4,3x,a,i2)', &
               d0, 'for  l = ', l
       end if
       !
       iii = idm
       d0 = 1.e9_pr
       !
       select case (idm)
          !
       case (1)
          cond = - matmul( mi31, (/ mat(1,4), mat(2,4), mat(3,4) /) )
          a = cond(1)
          b = cond(2)
          c = cond(3)
          d = 1.0_pr
       case (2)
          cond = - matmul( mi32, (/ mat(1,2), mat(3,2), mat(4,2) /) )
          a = cond(1)
          b = 1.0_pr
          c = cond(2)
          d = cond(3)
       case (3)
          cond = - matmul( mi33, (/ mat(1,3), mat(2,3), mat(4,3) /) )
          a = cond(1)
          b = cond(2)
          c = 1.0_pr
          d = cond(3)
       case (4)
          cond = - matmul( mi34, (/ mat(2,1), mat(3,1), mat(4,1) /) )
          a = 1.0_pr
          b = cond(1)
          c = cond(2)
          d = cond(3)
          !
       end select
       !
       mc = max( abs(a), abs(b), abs(c), abs(d) )
       a = a / mc
       b = b / mc
       c = c / mc
       d = d / mc
       !
       f1(1:nm-2) = a * yfu(1,1:nm-2) + b * yfd(1,1:nm-2)
       f1(nm+2:n) = c * ybu(1,nm+2:n) + d * ybd(1,nm+2:n)
       f1(nm-1:nm+1) = ( a * yfu(1,nm-1:nm+1) + b * yfd(1,nm-1:nm+1)     &
            &          + c * ybu(1,nm-1:nm+1) + d * ybd(1,nm-1:nm+1) ) / 2
       !
       f2(1:nm-2) = a * yfu(2,1:nm-2) + b * yfd(2,1:nm-2)
       f2(nm+2:n) = c * ybu(2,nm+2:n) + d * ybd(2,nm+2:n)
       f2(nm-1:nm+1) = ( a * yfu(2,nm-1:nm+1) + b * yfd(2,nm-1:nm+1)     &
            &          + c * ybu(2,nm-1:nm+1) + d * ybd(2,nm-1:nm+1) ) / 2
       !
    end if
    !
  end subroutine integrate

  !
  !!
  !

  subroutine numerov( y, g0i, g1i, eh12d, h12v, h12w, deta, &
       nstart, nstop, n )
    use linalg
    implicit none
    !----------- Integration by the Numerov method -----------
    !
    !! Solves:  y''(x) + k^2(x) * y(x) = 0, between x0 and x1.
    !! The function y is a two components vector and k^2 is a 2x2 matrix
    !! defined by:
    !!         ...
    !! The Numerov algorithm leads to:
    !!
    !! (1+h^2/12*k^2_{n+1}) * y_{n+1} - 2 * (1-5h^2/12*k^2_n) * y_n
    !!                                + (1+h^2/12*k^2_{n-1}) * y_{n-1} = 0
    !!
    !! where h^2 is the square of the integration step.
    !!
    !! Since this subroutine is called many times, it is very time
    !! consuming. It can be optimized by remarking that the iterative
    !! relation is of the form:
    !!
    !!   A_{n+1} * y_{n+1} = B_n * y_n - A_{n-1} * y_{n-1}
    !!
    !! The matrix A can be simply inverted:
    !!
    !!
    !!             1       /  h12v - e * h12d         h12w        \
    !! Inv[A] = -------- * |                                      |
    !!           det(A)    \      -h12w          h12v + e * h12d  /
    !!
    !! with:
    !!
    !!         det(A) = (h12v)^2 + (h12w)^2 - e^2 * (h12d)^2
    !!                = det1 - e^2 * det2 
    !!
    !!          det1 = (h12v)^2 + (h12w)^2 
    !!
    !!        det2 = (h12d)^2    ( and   edet2 = e^2 * det2  )
    !!
    !!        h12v = 1 - h^2/12 * v 
    !!
    !!        h12d = h^2/12 / d    ( and   eh12d = e * h12d  )
    !!
    !!        h12w = h^2/12 * w 
    !!
    !!  The equation can be rewritten:
    !!
    !!              G_{n+1} = 12 * y_n - 10 * G_n - G_{n-1},
    !!              y_{n+1} = A^{-1}_{n+1} G_{n+1}
    !!
    !! with  G_n = ( 1 - h^2/12 * k^2_n ) y_n.
    !!
    !!  So, starting with G_0, G_1 and y_1 we can compute G_2 and y_2
    !!    and iterate.
    !
    integer, intent(in) :: n, nstart, nstop
    real ( kind = pr ), intent(inout) :: y(2,n)
    real ( kind = pr ), intent(in), dimension(n) :: eh12d, h12v, h12w, deta
    real ( kind = pr ), intent(in) :: g0i(2), g1i(2)
    real ( kind = pr ) :: g0(2), g1(2), g2(2)
    integer :: i
    !
    g0 = g0i
    g1 = g1i
    if ( nstart > nstop ) then
       do i = nstart, nstop, -1
          g2 = 12 * y(:,i+1) - 10 * g1 - g0
          y(1,i) = &
               ( ( h12v(i) - eh12d(i) ) * g2(1) + h12w(i) * g2(2) ) / deta(i)
          y(2,i) = &
               ( ( h12v(i) + eh12d(i) ) * g2(2) - h12w(i) * g2(1) ) / deta(i)
          g0 = g1
          g1 = g2
       end do
    else
       do i = nstart, nstop
          g2 = 12 * y(:,i-1) - 10 * g1 - g0
          y(1,i) = &
               ( ( h12v(i) - eh12d(i) ) * g2(1) + h12w(i) * g2(2) ) / deta(i)
          y(2,i) = &
               ( ( h12v(i) + eh12d(i) ) * g2(2) - h12w(i) * g2(1) ) / deta(i)
          g0 = g1
          g1 = g2
       end do
    end if
    !
  end subroutine numerov

  !
  !! Derivative of a...
  !

  ! ...salar function

  subroutine derivative_sca( n, f, df, l )
    use param, only : h_12, h_60, boundary_condition
    implicit none
    integer, intent(in) :: n
    real ( kind = pr ), intent(in) :: f(n)
    real ( kind = pr ), intent(inout) :: df(n)
    integer, intent(in) :: l
    real ( kind = pr ) :: sig
    integer :: i
    !
    sig = ( modulo(l,2) - 0.5_pr ) * 2
    df(1) = ( 8.0_pr * f(2) - f(3) + sig * f(1) ) / h_12
    df(2) = ( 45.0_pr * ( f(3) - f(1) ) - 9.0_pr * f(4) &
         + f(5) - sig * f(1) ) / h_60
    df(3) = ( 45.0_pr * ( f(4) - f(2) ) - 9.0_pr * ( f(5) - f(1) ) &
         + f(6) ) / h_60
    !
    if ( boundary_condition == 0 .or.  &
         boundary_condition == 2 .and. l == 2 * ( l / 2 ) .or. &
         boundary_condition == 3 .and. l /= 2 * ( l / 2 ) ) then
       df(n) = ( -8.0_pr * f(n-1) + f(n) + f(n-2) ) / h_12
       df(n-1) = ( 45.0_pr * ( f(n) - f(n-2) ) + 9.0_pr * f(n-3) &
            - f(n) - f(n-4) ) / h_60
       df(n-2) = ( 45.0_pr * ( f(n-1) - f(n-3) ) &
            - 9.0_pr * ( f(n) - f(n-4) ) - f(n-5) ) / h_60
    end if
    if ( boundary_condition == 1 .or.  &
         boundary_condition == 2 .and. l /= 2 * ( l / 2 ) .or. &
         boundary_condition == 3 .and. l == 2 * ( l / 2 ) ) then
       df(n) = ( - 54 * f(n-1) + 45 * f(n) + 10 * f(n-2) - f(n-3) ) / h_60
       df(n-1) = ( 36 * f(n) + f(n-1) - 45 * f(n-2) &
            + 9 * f(n-3) - f(n-4) ) / h_60
       df(n-2) = ( 45 * ( f(n-1) - f(n-3) ) - 8 * f(n)  &
            + 9 * f(n-4) - f(n-5) ) / h_60
    end if
    !
    do i = 4, n - 3
       df(i) = ( 45.0_pr * ( f(i+1) - f(i-1) )  &
            - 9.0_pr * ( f(i+2) - f(i-2) ) &
            + f(i+3) - f(i-3) ) / h_60
    end do
    !
  end subroutine derivative_sca

  ! ...vector function

  subroutine derivative_vec( n, f, df, l )
    use param, only : h_12, h_60, boundary_condition
    implicit none
    integer, intent(in) :: n
    real ( kind = pr ), intent(in) :: f(2,n)
    real ( kind = pr ), intent(inout) :: df(2,n)
    integer, intent(in) :: l
    real ( kind = pr ) :: sig
    integer :: i
    !
    sig = ( modulo(l,2) - 0.5_pr ) * 2
    df(:,1) = ( 8.0_pr * f(:,2) - f(:,3) + sig * f(:,1) ) / h_12
    df(:,2) = ( 45.0_pr * ( f(:,3) - f(:,1) ) - 9.0_pr * f(:,4) &
         + f(:,5) - sig * f(:,1) ) / h_60
    df(:,3) = ( 45.0_pr * ( f(:,4) - f(:,2) )  &
         - 9.0_pr * ( f(:,5) - f(:,1) )        &
         + f(:,6) ) / h_60
    !
    if ( boundary_condition == 0 .or.  &
         boundary_condition == 2 .and. l == 2 * ( l / 2 ) .or. &
         boundary_condition == 3 .and. l /= 2 * ( l / 2 ) ) then
       df(:,n) = ( -8.0_pr * f(:,n-1) + f(:,n) + f(:,n-2) ) / h_12
       df(:,n-1) = ( 45.0_pr * ( f(:,n) - f(:,n-2) ) + 9.0_pr * f(:,n-3) &
            - f(:,n) - f(:,n-4) ) / h_60
       df(:,n-2) = ( 45.0_pr * ( f(:,n-1) - f(:,n-3) ) &
            - 9.0_pr * ( f(:,n) - f(:,n-4) ) - f(:,n-5) ) / h_60
    end if
    if ( boundary_condition == 1 .or.  &
         boundary_condition == 2 .and. l /= 2 * ( l / 2 ) .or. &
         boundary_condition == 3 .and. l == 2 * ( l / 2 ) ) then
       df(:,n) = ( - 54 * f(:,n-1) + 45 * f(:,n)  &
            + 10 * f(:,n-2) - f(:,n-3) ) / h_60
       df(:,n-1) = ( 36 * f(:,n) + f(:,n-1) - 45 * f(:,n-2) &
            + 9 * f(:,n-3) - f(:,n-4) ) / h_60
       df(:,n-2) = ( 45 * ( f(:,n-1) - f(:,n-3) ) - 8 * f(:,n)  &
            + 9 * f(:,n-4) - f(:,n-5) ) / h_60
    end if

    !
    do i = 4, n - 3
       df(1,i) = ( 45.0_pr * ( f(1,i+1) - f(1,i-1) ) &
            - 9.0_pr * ( f(1,i+2) - f(1,i-2) ) + f(1,i+3) - f(1,i-3) ) / h_60
       df(2,i) = ( 45.0_pr * ( f(2,i+1) - f(2,i-1) ) &
            - 9.0_pr * ( f(2,i+2) - f(2,i-2) ) + f(2,i+3) - f(2,i-3) ) / h_60
    end do
    !
  end subroutine derivative_vec

  !
  !! Check the accuracy of the solutions
  !!  *** UNUSED IN THE PRESENT VERSION OF THE CODE ***
  !

  subroutine accuracy_check( v, w, d, e, f1, f2, it, l, j, node, n, &
       f1ext, f2ext )
    use param, only : hh_12, bogolyubov, boundary_condition
    use eqdifstatic
    implicit none
    integer, intent(in) :: n, l, it, j, node
    real ( kind = pr ), dimension(1:n), intent(in) :: v, w, d, f1, f2
    real ( kind = pr ), intent(in) :: e, f1ext, f2ext
    real ( kind = pr ) :: sig, sn1, sn2, norm1, norm2, nx1, nx2
    real ( kind = pr ) :: norm, er1, er2
    real ( kind = pr ), parameter :: ermax = 0.1_pr
    integer :: i, iier
    !
    !! WARNING:
    !! This test is not relevant in the HF case when the energy
    !! is very close to 0.
    !! Better skip it in this case...
    !
    if ( .not. bogolyubov(it) .and. e < 0.001_pr ) return
    !
    sig = ( modulo(l,2) - 0.5_pr ) * 2
    !
    !! If l is even, we derivate f1 and f2 like odd functions
    !! If l is odd, we derivate f1 and f2 like even functions
    !
    d2f1(1) = - ( 30 + sig ) * f1(1) + 16 * f1(2) - f1(3)
    d2f2(1) = - ( 30 + sig ) * f2(1) + 16 * f2(2) - f2(3)
    !
    d2f1(2) = 16 * ( f1(1) + f1(3) ) - 30 * f1(2) - f1(4)
    d2f2(2) = 16 * ( f2(1) + f2(3) ) - 30 * f2(2) - f2(4)
    !
    if ( boundary_condition == 0 ) then
       d2f1(n) = - 29 * f1(n) + 16 * f1(n-1) - f1(n-2)
       d2f2(n) = - 29 * f2(n) + 16 * f2(n-1) - f2(n-2)
       d2f1(n-1) = 16 * ( f1(n) + f1(n-2) ) - 30 * f1(n-1) - f1(n-3)
       d2f2(n-1) = 16 * ( f2(n) + f2(n-2) ) - 30 * f2(n-1) - f2(n-3)
    else
       d2f1(n) = 16 * ( f1(n-2) + f1(n) ) - 30 * f1(n-1) - f1(n-3) - f1ext
       d2f2(n) = 16 * ( f2(n-2) + f2(n) ) - 30 * f2(n-1) - f2(n-3) - f2ext
       d2f1(n-1) = - 30 * f1(n-1) + 16 * ( f1(n-2) + f1(n) ) - f1(n-3) - f1ext
       d2f2(n-1) = - 30 * f2(n-1) + 16 * ( f2(n-2) + f2(n) ) - f2(n-3) - f2ext
    end if
    !
    do i = 3, n - 2
       d2f1(i) = - 30 * f1(i) + 16 * ( f1(i-1) + f1(i+1) ) - f1(i-2) - f1(i+2)
       d2f2(i) = - 30 * f2(i) + 16 * ( f2(i-1) + f2(i+1) ) - f2(i-2) - f2(i+2)
    end do
    d2f1 = d2f1 / hh_12
    d2f2 = d2f2 / hh_12
    !
    dif1 = ( - d2f1 + v * f1 + w * f2 ) * d
    dif2 = ( d2f2 - v * f2 + w * f1 ) * d
    !
    sn1 = dot_product( dif1, dif1 )
    sn2 = dot_product( dif2, dif2 )
    norm1 = dot_product( f1, f1 )
    norm2 = dot_product( f2, f2 )
    nx1 = dot_product( dif1, f1 )
    nx2 = dot_product( dif2, f2 )
    er1 = 0.0_pr
    er2 = 0.0_pr
    norm = norm1 + norm2
    if ( norm1 / norm > 1.e-7_pr ) &
         er1 = 1.0_pr - ( nx1 / sn1 ) * ( nx1 / norm1 )
    if ( norm2 / norm > 1.e-7_pr ) &
         er2 = 1.0_pr - ( nx2 / sn2 ) * ( nx2 / norm2 )
    print '(1x,a,/,1x,"IT = ",i1,"  N =",i2,"  l =",i2,    &
         &"  j = ",i2,"/2  E = ",f10.6,"  er",i1," =",e16.8,2x, &
         &"norms =",2e11.3)', &
         'ACCURACY OF THE WAVE FUNCTION:', &
         it, node + 1, l, j, e, iier, max( er1, er2 ), norm1, norm2

    if ( er1 > ermax .or. er2 > ermax ) then
       iier = 1
       if ( er2 > er1 ) iier = 2
       !
       !! In the HF case, there is no need to worry if the state
       !! is unoccuppied !
       !
       print '(1x,a,/,1x,"IT = ",i1,"  N =",i2,"  l =",i2,    &
            &"  j = ",i2,"/2  E = ",f10.6,"  er",i1," =",e16.8,2x, &
            &"norms =",2e11.3)', &
            'INSUFFICIENT ACCURACY OF THE WAVE FUNCTION:', &
            it, node + 1, l, j, e, iier, max( er1, er2 ), norm1, norm2
    end if
    if (debug) &
         print '(" it = ",i1," l = ",i2," 2j = ",i2,"  at E = ",&
         &f10.6,"   err = ",2e12.4)', it, l, j, e, er1, er2
    !
  end subroutine accuracy_check

  !
  !!
  !

  function count_nodes( f1, f2, npt ) result(cnode)
    implicit none
    integer, intent(in) :: npt
    real ( kind = pr ), intent(in) :: f1(npt), f2(npt)
    integer :: node1, node2, i
    integer :: cnode
    !
    node1 = 0
    node2 = 0
    do i = 2, npt
       !
       if ( f1(i) * f1(i-1) < 0.0_pr )  node1 = node1 + 2
       if ( f1(i) * f1(i-1) == 0.0_pr ) node1 = node1 + 1
       if ( f2(i) * f2(i-1) < 0.0_pr )  node2 = node2 + 2
       if ( f2(i) * f2(i-1) == 0.0_pr ) node2 = node2 + 1
       !
    end do
    cnode = node1 / 2
    if ( dot_product( f2, f2 ) > dot_product( f1, f1 ) ) &
         cnode = node2 / 2
    !
  end function count_nodes

  !
  !!
  !

end module eqdif
!
!=========================================================================
!
!
!=========================================================================
!
!! Module FORCES
!! this module contains the parameters of the Skyrme interactions:
!! SIII, SLy4, SLy5, Skm*  and SkP.
!
!=========================================================================
!

module forcesstatic

  use cste
  real ( kind = pr ), dimension(:), allocatable :: rho_tot, tau_tot, cur_tot, &
       rho_gamma, rho_gammap, rho_gamma2
  real ( kind = pr ), dimension(:), allocatable :: tmp
  real ( kind = pr ), dimension(:,:), allocatable :: drho, d2rho, &
       drho_p, d2rho_p, dcur
  real ( kind = pr ), dimension(:), allocatable :: drho_tot, d2rho_tot, &
       dcur_tot
  real ( kind = pr ), dimension(:,:), allocatable :: v1fx, v2fx, v3fx, v12fx, &
       v22fx, v1px, v12px, v2px, v3px, v22px

  interface allocateforces
     module procedure allocateforces
  end interface

contains

  subroutine allocateforces(n)
    implicit none
    integer, intent(in) :: n
    integer :: ifail
    allocate( rho_tot(n), tau_tot(n), cur_tot(n),               &
         rho_gamma(n), rho_gammap(n), rho_gamma2(n),            &
         tmp(-1:n), drho(1:n,2), d2rho(1:n,2),                  &
         drho_p(1:n,2), d2rho_p(1:n,2), dcur(1:n,2),            &
         drho_tot(n), d2rho_tot(n), dcur_tot(n),                &
         v1fx(1:n,2), v2fx(1:n,2), v3fx(1:n,2), v12fx(1:n,2),   &
         v22fx(1:n,2), v1px(1:n,2), v12px(1:n,2), v2px(1:n,2),  &
         v3px(1:n,2), v22px(1:n,2), &
         stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory... (forces)'
  end subroutine allocateforces

  !
  !!
  !

  subroutine deallocateforces()
    implicit none
    deallocate( rho_tot, tau_tot, cur_tot,            &
         rho_gamma, rho_gammap, rho_gamma2,           &
         tmp, drho, d2rho,                            &
         drho_p, d2rho_p, dcur,                       &
         drho_tot, d2rho_tot, dcur_tot,               &
         v1fx, v2fx, v3fx, v12fx,                     &
         v22fx, v1px, v12px, v2px,                    &
         v3px, v22px )
  end subroutine deallocateforces

  !
  !!
  !

end module forcesstatic

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

module skforces

  use cste
  use param
  implicit none

  !
  !! Skyrme force
  !
  !character ( len = 4 ) :: force
  !
  ! Parameters of the Skyrme force:
  !
  ! V(1,2) = t0 ( 1 + x0 P_sigma ) delta(1,2)
  !        + 1/2 t1 ( 1 + x1 P_sigma ) ( k'^2 delta(1,2) + delta(1,2) k^2 )
  !        + t2 ( 1 + x2 P_sigma ) k'.delta(1,2) k
  !        + 1/6 t3 ( 1 + x3 P_sigma ) rho^gamma  delta(1,2)
  !!       + 1/6 t4 ( 1 + x4 P_sigma ) rho^gamma2 delta(1,2)  !... FOR TESTS
  !        + i W0 ( sigma1 + sigma2 ).( k' x delta(1,2) k )
  !
  real ( kind = pr ), private, save :: t0, x0, t1, x1, t2, x2
  real ( kind = pr ), private, save :: t3, x3, gamma
  real ( kind = pr ), private, save :: t4, x4, gamma2
  real ( kind = pr ), private, save :: w
  !
  !! The coef. can be different in the pp channel
  !
  real ( kind = pr ), private, save :: t0p, x0p, t1p, x1p, t2p, x2p
  real ( kind = pr ), private, save :: t3p, x3p, gammap
  real ( kind = pr ), private, save :: wp
  !
  !! Coef. for the extended (spin-orbit) functionnal
  !
  real ( kind = pr ), private, save :: b_4, b_4p
  !
  !! Flags related with the Skyrme force
  !
  integer, private, save :: ixtls, ifecm, j2terms
  !
  !! (Force dependant) energy scale
  !
  real ( kind = pr ), private, save :: dmshb0
  !
  !! Auxiliary coefficients
  !
  real ( kind = pr ), private, save :: hqp,                     &
       af1, af2, af3, af4, af5, af6, af7, af8, af9, af10, af11, &
       ap1, ap2, ap3, ap4, ap5,                                 &
       bf1, bf2, bf3, bf4, bf5, bf6, bf7, bf8, bf9, bf10, bf11, &
       bf13, bf14,                                              &
       bp1, bp2, bp3, bp4, bp5, af3b, af4b, ap2b, bf3b, bf14b,  &
       bp2b, bf4b

  real ( kind = pr ), private, save :: a0, b0, a1, b1, a2, b2, &
       a3, b3, a4, b4

  interface set_force
     module procedure set_force
  end interface

  interface energy
     module procedure energy
  end interface

  interface updatefields
     module procedure updatefields
  end interface

  interface modify_force
     module procedure modify_force
  end interface

  interface set_init_param
     module procedure set_init_param
  end interface

contains

  !
  !! Set hb = hbar**2/2m and calculate auxilliary constantes
  !

  subroutine subparam()
    !
    !! This subroutine sets the values of different usefull parameters
    !! used to compute the different parts of the energy.
    !
    implicit none
    character ( len = 3 ) :: name1, name2
    integer :: it

    !
    !! Initialize several variables
    !
    vkblo = 0.0_pr
    hqp = h * qp
    !
    !! Switch between different schemes for c.m. correction
    !
    if ( ifecm == 0 ) then
       dmshb = dmshb0 / ( 1.0_pr - 1.0_pr / mass )
    else
       dmshb = dmshb0
    end if
    hb = 1.0_pr / dmshb
    !
    !! Fields energy coefficients
    !
    a0 = t0 * ( 2 + x0 ) / 4
    b0 = - t0 * ( 2 * x0 + 1 ) / 4
    a1 = ( t1 * ( 2 + x1 ) + t2 * ( 2 + x2 ) ) / 8
    b1 = - ( t1 * ( 2 * x1 + 1 ) - t2 * ( 2 * x2 + 1 ) ) / 8
    a2 = ( 3 * t1 * ( 2 + x1 ) - t2 * ( 2 + x2 ) ) / 32
    b2 = - ( 3 * t1 * ( 2 * x1 + 1 ) + t2 * ( 2 * x2 + 1 ) ) / 32
    a3 = t3 * ( 2 + x3 ) / 24
    b3 = - t3 * ( 2 * x3 + 1 ) / 24
    a4 = t4 * ( 2 + x4 ) / 24
    b4 = - t4 * ( 2 * x4 + 1 ) / 24
    !
    af1 = t0 * ( 2 + x0 ) / 4.0_pr
    af2 = - t0 * ( 2 * x0 + 1 ) / 4.0_pr
    af3 = t3 * ( 2 + x3 ) / 24.0_pr
    af3b = t4 * ( 2 + x4 ) / 24.0_pr
    af4 = - t3 * ( 2 * x3 + 1 ) / 24.0_pr
    af4b = - t4 * ( 2 * x4 + 1 ) / 24.0_pr
    af5 = ( t1 * ( 2 + x1 ) + t2 * ( 2 + x2 ) ) / 8.0_pr
    af6 = - ( 3 * t1 * ( 2 + x1 ) - t2 * ( 2 + x2 ) ) / 32.0_pr
    af7 = - ( t1 * ( 2 * x1 + 1 ) - t2 * ( 2 * x2 + 1 ) ) / 8.0_pr
    af8 = ( 3 * t1 * ( 2 * x1 + 1 ) + t2 * ( 2 * x2 + 1 ) ) / 32.0_pr
    !
    if ( j2terms == 1 ) then
       af9 = - ( t1 * x1 + t2 * x2 ) / 16.0_pr
       af10 = ( t1 - t2 ) / 16.0_pr
    else
       af9 = 0.0_pr
       af10 = 0.0_pr
    end if
    !
    !!
    !
    af11 = - w / 2.0_pr
    !
    !! Pairing energy coefficients  (FAUX)
    !
    ap1 = t0 * ( 1 - x0 ) / 4.0_pr
    ap2 = t3 * ( 1 - x3 ) / 24.0_pr
    ap2b = t4 * ( 1 - x4 ) / 24.0_pr
    ap3 = t1 * ( 1 - x1 ) / 4.0_pr
    ap4 = - ap3 / 4.0_pr
    ap5 = ( t2 * ( 1 + x2 ) + 2 * w ) / 8.0_pr
    !
    !!
    !
    if ( pairing_force /= 0 ) then
       ap1 = t0p / 4.0_pr
       ap2 = t3p / 24.0_pr
       ap2b = t3p / 24.0_pr
       ap3 = 0.0_pr
       ap4 = 0.0_pr
       ap5 = 0.0_pr
    end if
    !
    !!
    !
    bf1 = 2 * af1
    bf2 = 2 * af2
    bf3 = ( 2 + gamma ) * af3
    bf3b = ( 2 + gamma2 ) * af3b
    bf4 = 2 * af4
    bf4b = 2 * af4b
    bf5 = af5
    bf6 = 2 * af6
    bf7 = af7
    bf8 = 2 * af8
    bf9 = 2 * af9
    bf10 = 2 * af10
    bf11 = - af11
    bf13 = gammap * ap2
    bf14 = gamma * af4
    bf14b = gamma2 * af4b
    !
    !! Pairing potential coefficients
    !
    bp1 = 2 * ap1
    bp2 = 2 * ap2
    bp2b = 2 * ap2b
    bp3 = ap3
    bp4 = 2 * ap4
    bp5 = 2 * ap5
    write( name1, '(i3)' ) npr(1)
    write( name2, '(i3)' ) npr(2)
    extn = "_"//trim(adjustl(name1))//"_"//trim(adjustl(name2))
    !
    del = 0.0_pr
    do it = 1, 2
       if ( bogolyubov(it) ) then
          del(it) = 0.5_pr
          !
          !! If the nucleus is magic, there is no need to start
          !! with a too strong pairing gap
          !
          if ( npr(it) == 8   .or. &
               npr(it) == 20  .or. &
               npr(it) == 28  .or. &
               npr(it) == 50  .or. &
               npr(it) == 82  .or. &
               npr(it) == 126 .or. &
               npr(it) == 184 ) then ! and after 184 ?...
             del(it) = del(it) / 4.0_pr
          end if
          if ( npr(it) == 40 ) then
             del(it) = del(it) / 2.75_pr ! ...
          end if
          if ( npr(it) == 16 ) then
             del(it) = del(it) / 1.5_pr ! ...
          end if
          !
          if ( dripline == it ) del(it) = 0.5_pr
          !
       end if
    end do
    !
    !! Fermi levels
    !
    amb(1) = -8.0_pr
    amb(2) = -8.0_pr
    !
  end subroutine subparam

  !
  !!
  !

  subroutine set_force()
    !
    !!
    !
    implicit none
    !
    if ( memlog ) write( ulog, '(" SET_FORCE()")' )
    !
    gamma = 1.0_pr
    dmshb0 = 0.04823_pr
    w = 0.0_pr
    b_4 = 0.0_pr
    b_4p = 0.0_pr
    !
    !! These three parameters are for testing purposes only...
    !
    t4 = 0.0_pr
    x4 = 0.0_pr
    gamma2 = 2.0_pr / 3.0_pr
    !
    select case (trim(force))
       !
    case ("SIII")
       !
       !! SIII SET OF PARAMETERS
       !
       j2terms = 0
       !
       t0 = -1128.750000_pr
       x0 =     0.450000_pr
       t1 =   395.000000_pr
       x1 =     0.000000_pr
       t2 =   -95.000000_pr
       x2 =     0.000000_pr
       t3 = 14000.000000_pr
       x3 =     1.000000_pr
       t4 =     0.000000_pr
       x4 =     0.000000_pr
       w  =   120.000000_pr
       gamma =   1.000000_pr
       !
       ixtls = 0
       ifecm = 0
       !
       x0p =    0.000000_pr
       t1p =    0.000000_pr
       x1p =    0.000000_pr
       t2p =    0.000000_pr
       x2p =    0.000000_pr
       x3p =    0.000000_pr
       wp  =    0.000000_pr
       gammap = 1.000000_pr
       !
       select case (pairing_force)
          !
       case (0)
          stop ' Non sense !'
       case (1) ! Volume
          t0p =  -234.0_pr
          if ( i_cut == 1 ) t0p = -218.5_pr
          if ( .not. regularization ) t0p = -159.6_pr
          if ( i_cut == 2 .and. regularization ) t0p = -197.6_pr
          t3p = 0.0_pr
       case (2) ! Surface
          t0p = -483.2_pr
          if ( i_cut == 2 .and. regularization ) t0p = -807.0_pr
          t3p = - 6 * t0p / 0.16_pr
       case (3) ! Mixed
          if ( .not. regularization ) t0p = -248.5_pr
          if ( i_cut == 2 .and. regularization ) t0p = -316.9_pr
          t3p = - 6 * t0p / 0.32_pr
       case default
          stop ' Non sense !'
          !
       end select
       !
    case ("SLY4")
       !
       !! SLY4 SET OF PARAMETERS
       !
       j2terms = 0
       !
       t0 = -2488.913000_pr
       x0 =     0.834000_pr
       t1 =   486.818000_pr
       x1 =    -0.344000_pr
       t2 =  -546.395000_pr
       x2 =    -1.000000_pr
       t3 = 13777.000000_pr
       x3 =     1.354000_pr
       t4 =     0.000000_pr
       x4 =     0.000000_pr
       w  =   123.000000_pr
       gamma =   1.0_pr / 6.0_pr
       !
       dmshb0 = 1.0_pr / 20.7355300000_pr
       !
       ixtls = 0
       ifecm = 0
       x0p =    0.000000_pr
       t1p =    0.000000_pr
       x1p =    0.000000_pr
       t2p =    0.000000_pr
       x2p =    0.000000_pr
       x3p =    0.000000_pr
       wp  =    0.000000_pr
       gammap = 1.000000_pr
       !
       select case (pairing_force)
          !
       case (0)
          stop ' Non sense !'
       case (1) ! Volume
          t0p =  -234.0_pr
          if ( i_cut == 1 ) t0p = -218.5_pr
          if ( .not. regularization ) t0p = -186.5_pr
          if ( i_cut == 2 .and. regularization ) t0p = -233.0_pr
          t3p = 0.0_pr
       case (2) ! Surface
          t0p = -509.6_pr
          if ( i_cut == 2 .and. regularization ) t0p = -914.2_pr
          t3p = - 6 * t0p / 0.16_pr
       case (3) ! Mixed
          if ( .not. regularization ) t0p = -283.33_pr
          if ( i_cut == 2 .and. regularization ) t0p = -370.2_pr
          t3p = - 6 * t0p / 0.32_pr
       case default
          stop ' Non sense !'
          !
       end select
       !
    case ("SLY5")
       !
       !! SLY5 SET OF PARAMETERS
       !
       j2terms = 1
       !
       t0 = -2483.450000_pr
       x0 =     0.776000_pr
       t1 =   484.230000_pr
       x1 =    -0.317000_pr
       t2 =  -556.690000_pr
       x2 =    -1.000000_pr
       t3 = 13757.000000_pr
       x3 =     1.263000_pr
       t4 =     0.000000_pr
       x4 =     0.000000_pr
       w  =   125.000000_pr
       !
       gamma =   1.0_pr / 6.0_pr
       !
       dmshb0 = 1.0_pr / 20.7355300000_pr
       !
       ixtls = 0
       ifecm = 0
       x0p =    0.000000_pr
       t1p =    0.000000_pr
       x1p =    0.000000_pr
       t2p =    0.000000_pr
       x2p =    0.000000_pr
       x3p =    0.000000_pr
       wp  =    0.000000_pr
       gammap = 1.000000_pr
       !
       select case (pairing_force)
          !
       case (0)
          stop ' Non sense !'
       case (1) ! Volume
          t0p =  -234.0_pr
          if ( .not. regularization ) t0p = -179.9_pr
          if ( i_cut == 2 .and. regularization ) t0p = -222.7_pr
          t3p =    0.0_pr
       case (2) ! Surface
          t0p = -504.9_pr
          if ( i_cut == 2 .and. regularization ) t0p = -901.9_pr
          t3p = - 6 * t0p / 0.16_pr
       case (3) ! Mixed
          if ( .not. regularization ) t0p = -275.8_pr
          if ( i_cut == 2 .and. regularization ) t0p = -356.6_pr
          t3p = - 6 * t0p / 0.32_pr
       case default
          stop ' Non sense !'
          !
       end select
       !
    case ("SKP")
       !
       !! SKP SET OF PARAMETERS
       !
       j2terms = 1
       !
       t0 = -2931.6960000_pr
       x0 =     0.2921515_pr
       t1 =   320.6182000_pr
       x1 =     0.6531765_pr
       t2 =  -337.4091000_pr
       x2 =    -0.5373230_pr
       t3 = 18708.9600000_pr
       x3 =     0.1810269_pr
       t4 =     0.0000000_pr
       x4 =     0.0000000_pr
       w  =   100.0000000_pr
       b_4  = 50.0_pr
       b_4p = 50.0_pr
       gamma = 0.1666667_pr
       !
       dmshb0 = 0.04823_pr
       !
       ixtls = 1
       ifecm = 0
       !
       x0p =    0.000000_pr
       t1p =    0.000000_pr
       x1p =    0.000000_pr
       t2p =    0.000000_pr
       x2p =    0.000000_pr
       x3p =    0.000000_pr
       wp  =    0.000000_pr
       gammap = 1.000000_pr
       !
       t0p =  -160.0_pr
       t3p = 24175.0_pr
       !
       if ( pairing_force == 1 ) t3p = 0.0_pr
       if ( pairing_force == 2 ) t0p = t0
       select case (pairing_force)
          !
       case (0)
          t0p = t0
          x0p = x0
          t1p = t1
          x1p = x1
          t2p = t2
          x2p = x2
          t3p = t3
          x3p = x3
          gammap = gamma
          wp = w
       case (1) ! Volume
          t0p = -131.6_pr
          if ( i_cut == 2 .and. regularization ) t0p = -186.7_pr
          t3p = 0.0_pr
       case (2) ! Surface
          t0p = -417.5_pr
          if ( i_cut == 2 .and. regularization ) t0p = -848.0_pr
          t3p = - 6 * t0p / 0.16_pr
       case (3) ! Mixed
          if ( .not. regularization ) t0p = -211.3_pr
          if ( i_cut == 2 .and. regularization ) t0p = -308.6_pr
          t3p = - 6 * t0p / 0.32_pr
       end select
       !
    case ("SKM*")
       !
       j2terms = 0
       !
       t0 = -2645.00_pr
       x0 =     0.09_pr
       t1 =   410.00_pr
       x1 =     0.00_pr
       t2 =  -135.00_pr
       x2 =     0.00_pr
       t3 = 15595.00_pr
       x3 =     0.00_pr
       t4 =     0.00_pr
       x4 =     0.00_pr
       w  =   130.00_pr
       gamma = 1.0_pr / 6.0_pr
       !
       ixtls = 0
       ifecm = 0
       x0p =    0.000000_pr
       t1p =    0.000000_pr
       x1p =    0.000000_pr
       t2p =    0.000000_pr
       x2p =    0.000000_pr
       x3p =    0.000000_pr
       wp  =    0.000000_pr
       gammap = 1.000000_pr
       !
       t0p =  -173.0_pr
       t3p = 21591.0_pr
       !
       if ( pairing_force == 1 ) t3p = 0.0_pr
       if ( pairing_force == 2 ) t0p = t0
       if ( pairing_force == 3 ) then
          t0p = -283.7_pr
          t3p = - 6 * t0p / 0.16_pr
       end if
       x0p =    0.000000_pr
       t1p =    0.000000_pr
       x1p =    0.000000_pr
       t2p =    0.000000_pr
       x2p =    0.000000_pr
       x3p =    0.000000_pr
       wp  =    0.000000_pr
       gammap = 1.000000_pr
       !
       select case (pairing_force)
          !
       case (0)
          stop ' Non sense !'
       case (1) ! Volume
          t0p =  -234.0_pr
          if ( .not. regularization ) t0p = -148.6_pr
          if ( i_cut == 2 .and. regularization ) t0p = -184.7_pr
          t3p =    0.0_pr
       case (2) ! Surface
          t0p = -452.6_pr
          if ( i_cut == 2 .and. regularization ) t0p = -798.4_pr
          t3p = - 6 * t0p / 0.16_pr
       case (3) ! Mixed
          if ( .not. regularization ) t0p = -233.9_pr
          if ( i_cut == 2 .and. regularization ) t0p = -300.7_pr
          t3p = - 6 * t0p / 0.32_pr
       case default
          stop ' Non sense !'
          !
       end select
       !
    case default
       !
       print *, "Unknown force..."
       stop
       !
    end select
    !
  end subroutine set_force

  !
  !!
  !

  subroutine print_out_force_parameters()
    implicit none
    !
    if ( ixtls == 0 ) then
       print &
            '(1x,"     t0 =",f13.6,"   x0 =",f11.6,   &
            &    "    t1 =",f13.6,"   x1 =",f11.6,/,  &
            & 1x,"     t2 =",f13.6,"   x2 =",f11.6,   &
            &    "    t3 =",f13.6,"   x3 =",f11.6,/,  &
            & 2x," gamma =",f13.6,"    W =",f13.6,/,  &
            & 2x,"   t0'' =",f13.6,"  t3'' =",f13.6,  &
            &"   gamma'' =",f13.6)',  &
            t0, x0, t1, x1, t2, x2, &
            t3, x3, gamma, w, t0p, t3p, gammap
    else
       print '(1x,"     t0 =",f13.6,"   x0 =",f11.6,  &
            &     "    t1 =",f13.6,"   x1 =",f11.6,   &
            &/,1x,"     t2 =",f13.6,"   x2 =",f11.6,  &
            &     "    t3 =",f13.6,"   x3 =",f11.6,   &
            &/,2x," gamma =",f13.6,"   b_4= ",f10.6,  &
            &     "    b_4'' =",f11.6,                &
            &/,2x,"   t0'' =",f13.6,"   t3''=",f13.6, &
            &"   gamma''=",f13.6)',    &
            t0, x0, t1, x1, t2, x2, &
            t3, x3, gamma, b_4, b_4p, t0p, t3p, gammap
    end if
    if ( t4 /= 0.0_pr ) then
       print '(1x,"   t4 =",f13.6,"    x4 =",f11.6, &
            &/,2x," gamma2=",f13.6)',    &
            t4, x4, gamma2
    end if
    if ( ifecm == 0 ) then
       print '(2x,"hb2/2m = ",f10.7,3x,a)', 1.0_pr / dmshb0, &
            "(with 1-body c.m. correction)"
    else
       print '(2x,"hb2/2m = ",f10.7,3x,a)', 1.0_pr / dmshb0, &
            "(no c.m. correction)"
    end if
    !
  end subroutine print_out_force_parameters

  !
  !! Modify the parameters of the Skyrme force
  !

  subroutine set_init_param( skt0p, skt3p )
    implicit none
    real ( kind = pr ), intent(out) :: skt0p, skt3p
    skt0p = t0p
    skt3p = t3p
  end subroutine set_init_param

  !
  !!
  !

  subroutine modify_force( skt0p, skt3p )
    implicit none
    real ( kind = pr ), intent(in), optional :: skt0p, skt3p
    !
    if ( present(skt0p) ) t0p = skt0p
    if ( present(skt0p) ) t3p = skt3p
    !
  end subroutine modify_force

  !----------------------------------------------------------------------
  !
  !! Summations of energies and fields
  !
  !----------------------------------------------------------------------

  subroutine energy(npt)
    !
    !!
    !
    use forcesstatic
    use eqdif
    implicit none
    integer, intent(in) :: npt
    real ( kind = pr ) :: sum1, sum2
    integer :: i, isospin
    !
    !! various components of the energy density
    !
    real ( kind = pr ), dimension(npt) :: h_rho, h_rhotau, h_drho,  &
         h_gamma, h_gamma2, h_so
    !
    rho_tot = rho(:,1) + rho(:,2)
    rho_gamma = rho_tot**gamma
    rho_gamma2 = rho_tot**gamma2
    rho_gammap = rho_tot**gammap
    tau_tot = tau(:,1) + tau(:,2)
    cur_tot = cur(:,1) + cur(:,2)
    if ( regularization ) then
       do isospin = 1, 2
          geff(:,isospin) = ( t0p + t3p / 6 * rho_gammap )   &
               / ( 1 - ( rega(:,isospin) + regb(:,isospin) ) &
               * ( t0p + t3p / 6 * rho_gammap ) )
       end do
    end if
    !
    !! The densities rho and \tilde\rho are extrapolated to 0 using
    !! a forth order polynomial P(x) with the assumptions:
    !!   P'(0) = 0
    !!   P(-h) = P(h)
    !! Outside the box, they are set to 0.
    !
    !! For J and \tilde J, the extrapolation is:
    !!   P(0) = 0
    !!   P(-h) = -P(h)
    !
    !! The derivatives of the densities nabla(rho), nabla(\tilde\rho),
    !! Delta(rho) and Delta(\tilde\rho)
    !! are computed using a 5 points formula.
    !
    do isospin = 1, 2
       !
       tmp(1:npt) = rho(:,isospin)
       tmp(-1) = tmp(1)
       tmp(0) = ( 15 * tmp(1) - 6 * tmp(2) + tmp(3) ) / 10.0_pr
       call d1_and_d2( tmp, drho(:,isospin), d2rho(:,isospin), npt )
       !
       tmp(1:npt) = rho_p(:,isospin)
       tmp(-1) = tmp(1)
       tmp(0) = ( 15 * tmp(1) - 6 * tmp(2) + tmp(3) ) / 10.0_pr
       call d1_and_d2( tmp, drho_p(:,isospin), d2rho_p(:,isospin), npt )
       !
       tmp(1:npt) = cur(:,isospin)
       tmp(-1) = - tmp(1)
       tmp(0) = 0.0_pr
       call d1_and_d2( tmp, dcur(:,isospin), n = npt )
       !
    end do
    !
    tmp(1:npt) = rho_tot(:)
    tmp(-1) = tmp(1)
    tmp(0) = ( 15 * tmp(1) - 6 * tmp(2) + tmp(3) ) / 10.0_pr
    call d1_and_d2( tmp, drho_tot(:), d2rho_tot(:), npt )
    !
    tmp(1:npt) = cur_tot(:)
    tmp(-1) = - tmp(1)
    tmp(0) = 0.0_pr
    call d1_and_d2( tmp, dcur_tot(:), n = npt )
    !
    !! Coulomb potential
    !
    sum1 = 0.0_pr
    sum2 = 0.0_pr
    do i = 1, npt
       sum1 = sum1 + qp * r2(i) * rho(i,2)
       sum2 = sum2 + qp * r(i) * rho(i,2)
       vc(i) = sum1 / r(i) - sum2
    end do
    vc = echarg * h * ( vc + sum2 )
    !
    !!............................................. Energy
    !
    if ( minval(rho_tot) < -1.e-10_pr ) then
       print '(" Negative density in subroutine ENERGY in points:")'
       do i = 1, npt
          if ( rho_tot(i) < 1.e-10_pr ) print *, i
       end do
    end if
    !
    !!............................................. Direct Coulomb energy
    !
    coulomb_energy = sum( vc * rho(:,2) * r2 ) / 2.0_pr
    !
    !!............................................. Exchange Coulomb energy
    !
    ex_coulomb_energy = sum( rho(:,2)**t43 * r2 ) * coef
    !
    !!............................................. Kinetic energy
    !
    do isospin = 1, 2
       kinetic_energy(isospin) = hb * sum( tau(:,isospin) * r2 )
       !
       if ( regularization ) then
          pairing_energy(isospin) = &
               sum( r2 * ( rho_p(:,isospin)**2 * geff(:,isospin) ) ) / 4
       else
          pairing_energy(isospin) = sum( r2 * ( rho_p(:,isospin)**2  &
               * ( ap1 + ap2 * rho_gammap )                          &
               + ap3 * rho_p(:,isospin) * tau_p(:,isospin)           &
               - ap4 * drho_p(:,isospin)**2                          &
               + ap5 * cur_p(:,isospin)**2 ) )
       end if
       !
       pairing_kinetic_energy(isospin) = sum( r2 * (              &
            + ap3 * rho_p(:,isospin) * tau_p(:,isospin)           &
            - ap4 * drho_p(:,isospin)**2                          &
            + ap5 * cur_p(:,isospin)**2 ) )
    end do
    kinetic_energy(3) = kinetic_energy(1) + kinetic_energy(2)
    pairing_energy(3) = pairing_energy(1) + pairing_energy(2)
    pairing_kinetic_energy(3) =                                   &
         pairing_kinetic_energy(1) + pairing_kinetic_energy(2)
    !
    !!............................................. Spin orbit energy
    !
    if ( ixtls == 0 ) then
       spin_orbit_energy = sum( r2 * af11 * ( cur_tot * drho_tot &
            + cur(:,1) * drho(:,1) + cur(:,2) * drho(:,2) ) )
    else
       spin_orbit_energy = - sum( r2 * ( b_4 * cur_tot * drho_tot &
            + b_4p * ( cur(:,1) * drho(:,1) + cur(:,2) * drho(:,2) ) ) )
    end if
    !
    !!............................................. Field energy
    !
    h_rho = a0 * rho_tot**2 + b0 * ( rho(:,1)**2 + rho(:,2)**2 )
    !
    h_rhotau = a1 * rho_tot * tau_tot  &
         + b1 * ( rho(:,1) * tau(:,1) + rho(:,2) * tau(:,2 ) )
    !
    h_drho = a2 * drho_tot**2 + b2 * ( drho(:,1)**2 + drho(:,2)**2 )
    !
    h_gamma = a3 * rho_gamma * rho_tot**2  &
         + b3 * rho_gamma * ( rho(:,1)**2 + rho(:,2)**2 )
    !
    h_gamma2 = a4 * rho_gamma2 * rho_tot**2  &
         + b4 * rho_gamma2 * ( rho(:,1)**2 + rho(:,2)**2 )
    !
    h_so = af9 * cur_tot**2                                          &
         + af10 * ( cur(:,1)**2 + cur(:,2)**2 )
    !
    field_energy = sum( r2 * ( h_rho + h_rhotau + h_drho   &
         + h_gamma + h_gamma2 + h_so ) )
    !
    !!............................................. Rearrangement energy
    !
    rearangment_energy = sum( r2 * rho_gamma                            &
         * ( af3 * rho_tot**2 + af4 * ( rho(:,1)**2 + rho(:,2)**2 ) ) ) &
         * gamma / 2.0_pr
    rearangment_energy = rearangment_energy + sum( r2 * rho_gamma2        &
         * ( af3b * rho_tot**2 + af4b * ( rho(:,1)**2 + rho(:,2)**2 ) ) ) &
         * gamma2 / 2.0_pr
    rearangment_energy = rearangment_energy + sum( r2 * rho_gammap      &
         * ap2 * ( rho_p(:,1)**2 + rho_p(:,2)**2 ) )                    &
         * gammap / 2.0_pr
    rearangment_energy = rearangment_energy - ex_coulomb_energy / 3.0_pr
    !
    !!
    !
    kinetic_energy         = kinetic_energy         * h * qp ! 3
    field_energy           = field_energy           * h * qp ! 1
    spin_orbit_energy      = spin_orbit_energy      * h * qp ! 1
    coulomb_energy         = coulomb_energy         * h * qp ! 1
    ex_coulomb_energy      = ex_coulomb_energy      * h * qp ! 1
    pairing_energy         = pairing_energy         * h * qp ! 3
    pairing_kinetic_energy = pairing_kinetic_energy * h * qp ! 3
    rearangment_energy     = rearangment_energy     * h * qp ! 1
    !
    !!............................................. Total energy and E/A
    !
    total_energy = kinetic_energy(3) + field_energy &
         + spin_orbit_energy + coulomb_energy &
         + ex_coulomb_energy + pairing_energy(3)
    energy_per_nucleon = total_energy / mass
    !
  end subroutine energy

  !
  !!
  !

  subroutine updatefields()
    !
    !!
    !
    use forcesstatic
    use eqdif
    implicit none
    real ( kind = pr ) :: ymu, ypmu, q
    integer :: isospin
    !
    !! Calculation of the new field and pairing potentials and
    !! effective masses
    !
    ymu  = 1.0_pr - xmu
    ypmu = 1.0_pr - xpmu
    do isospin = 1, 2
       q = real( isospin, kind = pr ) - 1.0_pr
       v1fx(:,isospin) = rho_tot &
            * ( bf1 + bf3 * rho_gamma + bf3b * rho_gamma2 )                &
            + rho(:,isospin)                                               &
            * ( bf2 + bf4 * rho_gamma + bf4b * rho_gamma2 )                &
            + rho_gamma * ( bf14 * ( rho(:,1)**2 + rho(:,2)**2 )           &
            ) / ( rho_tot + eps )                                          &
            + rho_gamma2 * ( bf14b * ( rho(:,1)**2 + rho(:,2)**2 )         &
            ) / ( rho_tot + eps )                                          &
            + rho_gammap * ( bf13 * ( rho_p(:,1)**2                        &!
            + rho_p(:,2)**2 ) ) / ( rho_tot + eps )                        &!
            + bf5 * tau_tot                                                &
            + bf6 * ( d2rho_tot + 2 * drho_tot / r )                       &
            + bf7 * tau(:,isospin)                                         &
            + bf8 * ( d2rho(:,isospin) + 2 * drho(:,isospin) / r )         &
            + q * ( vc + t43 * coef * rho(:,2)**t13 )
       !
       !!................................................. Effective mass
       !
       v2fx(:,isospin) = bf5 * rho_tot + bf7 * rho(:,isospin) + hb
       !
       !!................................................. Spin orbit
       !
       if ( ixtls == 0 ) then
          v1fx(:,isospin) = v1fx(:,isospin)          &
               + bf11 * ( dcur_tot + dcur(:,isospin) &
               + 2 * ( cur_tot + cur(:,isospin) ) / r )
          v3fx(:,isospin) = bf9 * cur_tot + bf10 * cur(:,isospin) &
               - bf11 * ( drho_tot + drho(:,isospin) )
       else
          v1fx(:,isospin) = v1fx(:,isospin) &
               + b_4  * dcur_tot             &
               + b_4p * dcur(:,isospin)      &
               + 2 * ( b_4 * cur_tot + b_4p * cur(:,isospin) ) / r
          v3fx(:,isospin) = bf9 * cur_tot + bf10 * cur(:,isospin) &
               - b_4 * drho_tot  - b_4p * drho(:,isospin)
       end if
       !
       !!................................................. Pairing field
       !
       if ( regularization .and. bogolyubov(isospin) ) then
          v1px(:,isospin) = geff(:,isospin) * rho_p(:,isospin) / 2
       else
          v1px(:,isospin) = rho_p(:,isospin) * ( bp1 + bp2 * rho_gammap )  &
               + bp3 * tau_p(:,isospin)                                    &
               + bp4 * ( d2rho_p(:,isospin) + 2 * drho_p(:,isospin) / r )
       end if
       !
       v2px(:,isospin)  = bp3 * rho_p(:,isospin)
       v3px(:,isospin)  = bp5 * cur_p(:,isospin)
       v12fx(:,isospin) = bf5 * drho_tot  + bf7 * drho(:,isospin)
       v22fx(:,isospin) = bf5 * d2rho_tot + bf7 * d2rho(:,isospin)
       v12px(:,isospin) = bp3 * drho_p(:,isospin)
       v22px(:,isospin) = bp3 * d2rho_p(:,isospin)
       v1fx(:,isospin)  = v1fx(:,isospin) + v12fx(:,isospin) / r
       v1px(:,isospin)  = v1px(:,isospin) + v12px(:,isospin) / r
       !
    end do
    !
    if ( .not.fixed_ph ) then
       v1f  = xmu * v1f  + ymu * v1fx
       v2f  = xmu * v2f  + ymu * v2fx
       v3f  = xmu * v3f  + ymu * v3fx
       v12f = xmu * v12f + ymu * v12fx
       v22f = xmu * v22f + ymu * v22fx
    end if
    !
    if ( .not.fixed_pp ) then
       v1p  = xpmu * v1p  + ypmu * v1px
       v2p  = xpmu * v2p  + ypmu * v2px
       v3p  = xpmu * v3p  + ypmu * v3px
       v12p = xpmu * v12p + ypmu * v12px
       v22p = xpmu * v22p + ypmu * v22px
       !
       !! cut the pairing fields for r > r_cut
       !
       where ( r > r_cut )
          v1p(:,1)  = 0.0_pr
          v1p(:,2)  = 0.0_pr
          v2p(:,1)  = 0.0_pr
          v2p(:,2)  = 0.0_pr
          v3p(:,1)  = 0.0_pr
          v3p(:,2)  = 0.0_pr
          v12p(:,1) = 0.0_pr
          v12p(:,2) = 0.0_pr
          v22p(:,1) = 0.0_pr
          v22p(:,2) = 0.0_pr
       end where
       !
    end if
    !
    !! Calculate particle numbers, radii and meanfield pairing gaps
    !
    do isospin = 1, 2
       del(isospin) = sum( ( v2p(:,isospin) * tau(:,isospin) &
            + v3p(:,isospin) * cur(:,isospin)              &
            + v1p(:,isospin) * rho(:,isospin) ) * r2 )
       rnumpart(isospin) = sum( rho(:,isospin) * r2 )
    end do
    del = del / rnumpart
    rnumpart = rnumpart * qp * h
    !
    !! The (1-body) center of mass correction must be updated
    !! if the particle number is not fixed (dripline)
    !
    if ( dripline /= 0 .and. ifecm == 0 ) then
       mass = sum(rnumpart)
       dmshb = dmshb0 / ( 1.0_pr - 1.0_pr / mass )
    end if
    !
  end subroutine updatefields

  !
  !!
  !

end module skforces
!
!=========================================================================
!

!
!==================================================================
!
!! Module IO
!!
!! Input and output subroutines.
!
!==================================================================
!

module io

  use cste
  use param
  use skforces
  implicit none

  interface read_input
     module procedure read_input
  end interface

  interface print_out_wf
     module procedure print_out_wf
  end interface

  interface print_out_can
     module procedure print_out_can
  end interface

contains

  subroutine open_files()
    implicit none
    integer :: v(8)
    character ( len = 4 ) :: year
    character ( len = 2 ) :: day, month, hour, min, sec
    character ( len = 17 ) :: ename
    !
    call date_and_time( values = v )
    write( year, '(i4)' ) v(1)
    if ( v(2) <= 9 ) then
       write( month, '(2i1)' ) 0, v(2)
    else
       write( month, '(i2)' ) v(2)
    end if
    if ( v(3) <= 9 ) then
       write( day, '(2i1)' ) 0, v(3)
    else
       write( day, '(i2)' ) v(3)
    end if
    if ( v(5) <= 9 ) then
       write( hour, '(2i1)' ) 0, v(5)
    else
       write( hour, '(i2)' ) v(5)
    end if
    if ( v(6) <= 9 ) then
       write( min, '(2i1)' ) 0, v(6)
    else
       write( min, '(i2)' ) v(6)
    end if
    if ( v(7) <= 9 ) then
       write( sec, '(2i1)' ) 0, v(7)
    else
       write( sec, '(i2)' ) v(7)
    end if
    ename = day//'.'//month//'.'//year//'-'//hour//min//sec
    if ( memlog ) &
         open( unit = ulog,   file = out//'hfb.log',     status = 'unknown' )
    open( unit = uinput, file = hfb_input,          status = 'unknown' )
    open( unit = usumup, file = out//'hfb.summary', status = 'unknown' )
    !
    write( usumup, '("#  Format: 2i4,2i3,f6.1,2f9.4,2f7.4,2f11.5,2f10.5,",&
         &"f11.5,f10.5,f10.5,f15.7,4f8.5")' )
    write( usumup, '("# N   Z  J_max E_cut   Fermi energies",&
         &"   |Mean gaps|     Kinetic energies    Pairing energies",&
         &"    SO energy  Coulomb energies       Total        r_n     r_p",&
         &"    r_tot   r_ch")' )
    write( usumup, '("#",7("-"),2(1x,5("-")),1x,17("-"),1x,13("-"),&
         &1x,21("-"),1x,19("-"),1x,10("-"),1x,19("-"),1x,14("-"),&
         &4(1x,7("-")))' )
    !
  end subroutine open_files

  !
  !!
  !

  subroutine init_clock()
    implicit none
    integer :: v(8)
    character ( len = 3 ), parameter :: month(12) = &
         (/ "jan", "feb", "mar", "apr", "may", "jun", &
         "jul", "aug", "sep", "oct", "nov", "dec" /)
    character ( len = 2 ) :: min, sec
    !
    call date_and_time( values = v )
    if ( v(6) < 10 ) then
       write(min,'("0",i1)') v(6)
    else
       write(min,'(i2)') v(6)
    end if
    if ( v(7) < 10 ) then
       write(sec,'("0",i1)') v(7)
    else
       write(sec,'(i2)') v(7)
    end if
    if ( memlog ) &
         write( ulog, &
         '("Calculation start on ",a,"-",i2,", ",i4," at ",i2,"h", &
         &a,":",a)') &
         month(v(2)), v(3), v(1), v(5), min, sec
    time0 = v(5) * 60 * 60 + v(6) * 60 + v(7) + v(8) / 1000.0_pr
    time_previous = time0
    !
  end subroutine init_clock

  !
  !!
  !

  subroutine how_long()
    implicit none
    integer :: v(8)
    real ( kind = pr ) :: time1
    !
    call date_and_time( values = v )
    time1 = v(5) * 60 * 60 + v(6) * 60 + v(7) + v(8) / 1000.0_pr
    if ( memlog ) &
         write( ulog, &
         '(" Total time =",f10.3," sec.   Time since last call =",&
         & f10.3," sec.")') time1 - time0, time1 - time_previous
    time_previous = time1
    !
  end subroutine how_long

  !
  !!
  !

  subroutine read_input( npt, bye )
    implicit none
    integer, intent(inout) :: npt
    logical, intent(out), optional :: bye
    !
    !! parameters that are valid for all runs are saved
    !
    real ( kind = pr ), save :: integ_step
    integer, save :: mesh_points, nucl
    !
    !! Others are not
    !
    character ( len = 40 ), save :: read_pot
    !    integer :: iormax
    integer :: nn, nl, nj, pn, pl, pj       ! Blocking.....UNUSED
    character ( len = 3 ) :: ntype, ptype   ! .............UNUSED
    real ( kind = pr ) :: e_step
    !
    !! Local variables
    !
    integer :: ifail, i
    !
    !! Local parameters of the Skyrme force
    !
    real ( kind = pr ), save :: skt0p, skt3p
    !
    !! Files
    !
    logical, save :: files_are_opened
    character ( len = 72 ) :: str ! dummy string
    namelist /input/ force, mesh_points, it_max, bogolyubov,      &
         eps_energy, integ_step, pairing_force, fixed_ph,         &
         boundary_condition, regularization, xmu, max_delta
    namelist /nucleus/ proton, neutron, j_max, read_pot, cut_off, &
         skt0p, skt3p, r_cut, e_step, densities, quasiparticles,  &
         meanfields, canonical_states, cut_diffuseness
    namelist /blocking/ ntype, nn, nl, nj, ptype, pn, pl, pj
    !
    boundary_condition = 0
    inquire( unit = uinput, opened = files_are_opened )
    if ( .not. files_are_opened ) then
       call open_files()
       !
       !! Scan the input file to see how many nuclei I have to deal with
       !
       n_nuclei = 0
       do
          read( uinput, '(a)', end = 10 ) str
          str = adjustl(str)
          if ( str(1:8) == "&nucleus" ) n_nuclei = n_nuclei + 1
       end do
10     rewind(uinput)
       if ( n_nuclei == 0 ) then
          print '(/,10x,i4," NOTHING TO DO...",/)'
          stop
       end if
       if ( n_nuclei == 1 ) then
          print '(/,28x,"1 nucleus to compute",/)'
       else
          print '(/,25x,i4," nuclei to compute",/)', n_nuclei
       end if
       !
       !! Read the parameters for the run
       !
       !! default values
       !
       force = "SLY4"
       fixed_ph = .false.
       xmu = 0.0_pr
       max_delta = 5.e-7_pr
       !
       read( unit = uinput, nml = input )
       !
       if ( max_delta < eps_energy ) max_delta = max( eps_energy, 5.e-7_pr )
       if ( regularization ) then
          i_cut = 2
          cut_diffuseness = 0.0_pr
       else
          i_cut = 3
          cut_diffuseness = 1.0_pr
       end if
       !
       call set_xmu( 0.8_pr )
       call set_force()
       !
       if ( eps_energy < 1.e2_pr * eps ) then
          print '(/,18x,41("*"),/,18x, &
               & "WARNING: You ask for a too high precision")'
          eps_energy = 2.e2_pr * eps
          print &
               '(18x,"         Reset eps_energy to ",e7.1,/,18x,41("*"),/)', &
               eps_energy
       end if
       !
       if ( force == "" ) force = "SLY4"
       print '(29x,"Skyrme Force = ",a)', force
       nucl = 1 ! first nucleus to read
       cut_off = 0.0_pr
       read_pot = ""
       j_max = (/ 21, 21 /)
       r_cut = 30.0_pr
       e_step = 0.0_pr
       quasiparticles   = .false.
       canonical_states = .false.
       meanfields       = .false.
       densities        = .false.
       !
       call set_init_param( skt0p, skt3p )
       read( uinput, nml = nucleus )
       cano = canonical_states
       if ( i_cut == 3 .and. cut_off == 0.0_pr ) cut_off = 60.0_pr
       call modify_force( skt0p, skt3p )
       !
       !! Dripline ?
       !
       dripline = 0
       if ( proton <= 0 .and. neutron <= 0 ) stop 'Non sense !...'
       if ( neutron <= 0 ) dripline = 1
       if ( proton <= 0 ) dripline = 2
       proton = iabs(proton)
       neutron = iabs(neutron)
    else
       nucl = nucl + 1
       if ( nucl > n_nuclei ) then
          if ( present(bye) ) bye = .true.
          return
       else
          e_step = 0.0_pr
          quasiparticles   = .false.
          canonical_states = .false.
          meanfields       = .false.
          densities        = .false.
          !
          read( uinput, nml = nucleus )
          if ( i_cut == 3 .and. cut_off == 0.0_pr ) cut_off = 60.0_pr
          cano = canonical_states
          call modify_force( skt0p, skt3p )
          !
          !! Dripline ?
          !
          dripline = 0
          if ( proton  <= 0 .and. neutron <= 0 ) stop 'Non sense !...'
          if ( neutron <= 0 ) dripline = 1
          if ( proton  <= 0 ) dripline = 2
          proton = iabs(proton)
          neutron = iabs(neutron)
          if ( present(bye) ) bye = .false.
       end if
    end if
    npr(1) = neutron
    npr(2) = proton
    print '(/,21x," --  N = ",i3,"    --    Z =",i3,"    --",/)', npr
    !
    !******************************************************************
    !
    !! The following lines handle the blocking,
    !! this possibility is not described in the CPC article
    !! and has not been fully tested yet...
    !! So, don't use it ! And, anyway, the time odd terms have not
    !! been implemented in the functionnal yet.
    !
    !******************************************************************
    !
    !! Blocking
    !
    iing = 0
    ning = 0
    ling = 0
    jing = 0
    if ( neutron /= 2 * ( neutron / 2 ) .or. &
         proton /= 2 * ( proton / 2 ) ) then
       read( uinput, nml = blocking )
       print '(6x,46("="))'
       if ( neutron == 2 * ( neutron / 2 ) .or. &
            proton == 2 * ( proton / 2 ) ) then
          print '(6x,"|",14x,"Odd-Even Nucleus",14x,"|")'
       else
          print '(6x,"|",14x,"Odd-Odd Nucleus",15x,"|")'
       end if
       if ( ntype == "1qp" .or. ntype == "1QP" ) then
          iing(1) = 3
          ning(1) = nn
          ling(1) = nl
          jing(1) = nj
          print '(6x,"|  Neutrons: ",a,", n =",i3,",  l =",  &
               & i3,"  j =",i3,"/2  |")', &
               ntype, nn, nl, nj
       end if
       if ( ptype == "1qp" .or. ptype == "1QP" ) then
          iing(2) = 3
          ning(2) = pn
          ling(2) = pl
          jing(2) = pj
          print '(6x,"|   Protons: ",a,", n =",i3,",  l =",  &
               & i3,"  j =",i3,"/2  |")', &
               ptype, pn, pl, pj
       end if
       print '(6x,46("="),/)'
    end if
    !
    pot = read_pot
    !
    print '(1x,32("-"),"  INPUT DATA  ",33("-"))'
    print '(2x,"it_max =",i4,           &
         &7x,"eps energ. =",1pe10.3,7x,"max delta =",1pe10.3)', &
         it_max, eps_energy, max_delta
    print '(2x,"Boundary_condition =",i2)', boundary_condition
    write( *, '(2x,"Pairing force: " )', advance = "no" )
    select case(pairing_force)
    case (0)
       write( *, '("Full Skyrme interaction   ")', advance = "no" )
    case (1)
       write( *, '("Volume pairing            ")', advance = "no" )
    case (2)
       write( *, '("Surface pairing           ")', advance = "no" )
    case (3)
       write( *, '("Volume + Surface pairing  ")', advance = "no" )
    case default
       print *
       stop "Don't known what you want..."
    end select
    !
    if ( .not. bogolyubov(1) .and. .not. bogolyubov(2) ) print &
         '("No Bogolyubov transformation")'
    if ( bogolyubov(1) .and. .not. bogolyubov(2) ) print &
         '("Bogolyubov transf. for neutrons only")'
    if ( bogolyubov(2) .and. .not. bogolyubov(1) ) print &
         '("Bogolyubov transf. for protons only")'
    if ( bogolyubov(1) .and. bogolyubov(2) ) print &
         '("Bogolyubov transf. for both isospins")'
    !
    npt = mesh_points
    rbox = npt * integ_step
    estep0 = 11110.0_pr / rbox**3 - 550.0_pr / rbox**2 + 12.0_pr / rbox
    if ( rbox >= 45.0_pr ) estep0 = estep0 / 1.25_pr
    estep0 = min( 0.5_pr, estep0 * 90.0_pr / sum(npr) )
    if ( e_step /= 0.0_pr ) estep0 = e_step
    estep0min = estep0 / 4
    !
    h = integ_step
    h12 = h * h / 12.0_pr ! For Numerov
    h_12 = 12 * h         ! \
    hh_12 = 12 * h * h    !  | For derivations
    h_60 = 60 * h         ! /
    h_120 = 120 * h       !/
    !
    if ( .not. allocated(r) ) then
       allocate( r(npt), r2(npt), qpr2(npt), deriv(-3:3), stat = ifail )
       if ( ifail /= 0 ) stop 'Memory problem...'
    end if
    r = (/ ( i * h, i = 1, npt ) /) 
    r2 = r * r
    qpr2 = 1.0_pr / qp / r2
    deriv(-3:3) = (/ -1.0_pr, 9.0_pr, -45.0_pr, 0.0_pr, &
         &           45.0_pr, -9.0_pr, 1.0_pr /) / 60.0_pr
    deriv = deriv / h
    !
    print &
         '(2x,"npt =",i4,5x,"integ. step =",f6.3," fm",&
         &5x,"Rbox =",f6.2," fm",5x,"r_cut = ",f6.2," fm")', &
         npt, h, rbox, r_cut
    print '(1x,79("-"))'
    print '(" Skyrme force parameters:")'
    call print_out_force_parameters()
    print '(1x,79("-"))'
    print '(2x,"xmu =",f5.2,7x,"2Jmax =",i3," (neutrons),",2x,&
         &i3," (protons)")', xmu, j_max
    print '(2x,"cut_off =",f8.3," MeV",2x,"(with diff. ",f6.3," MeV)", &
         & 6x,"energ. step = ",f6.3," MeV")', cut_off, cut_diffuseness, &
         estep0
    if ( regularization ) then
       print '(2x,"Use the regularization method of pairing")'
    end if
    print '(1x,79("-"))'
    !
  end subroutine read_input

  !
  !!
  !

  subroutine deallocate_r()
    implicit none
    deallocate( r, r2, qpr2, deriv ) 
  end subroutine deallocate_r

  !
  !!
  !

  subroutine set_xmu(x)
    implicit none
    real ( kind = pr ), intent(in) :: x
    !
    if ( xmu == 0.0_pr ) xmu = x
    xmu0 = xmu
    xpmu = 0.9_pr * xmu ! the change in the pairing field can be a
    xpmu0 = xpmu        ! little bit faster than in the normal field.
    !
  end subroutine set_xmu

  !
  !!
  !

  subroutine print_out_wf( n0, n1, npt, it, l, j )
    implicit none
    integer, intent(in) :: n0, n1, npt
    integer, intent(in) :: l
    integer, intent(in), optional :: it, j
    integer :: i, k
    character ( len = 4 ) :: name
    !
    do i = n0, n0 + n1 - 1
       if ( i < 10 ) write(name,'("000",i1)') i
       if ( i > 9 .and. i < 100 ) write(name,'("00",i2)') i
       if ( i > 99 .and. i < 1000 ) write(name,'("0",i3)') i
       if ( i > 999 ) write(name,'(i4)') i
       open( unit = uwave, file = out//'qp'//name//trim(extn)//'.gfx', &
            status = 'unknown' )
       if ( present(it) .and. present(j) ) &
            write( uwave, &
            '("# it = ",i1,",  l = ",i2,",  j = ",i2,"/2")' ) it, l, j
       write( uwave, '("# E = ",f12.8)' ) ehfb(i)
       write( uwave, '("#",/,"#  r (fm)       psi_1          psi_2",&
            &"          dpsi_1/dr       dpsi_2/dr")' )
       write( uwave, '("# -------",4(" ---------------"))' )
       write( uwave, '(f8.3,1x,4e16.8)' ) 0.0_pr, 0.0_pr, 0.0_pr, &
            ( 600 * ff(:,i,1) - 600 * ff(:,i,2) &
            + 400 * ff(:,i,3) - 150 * ff(:,i,4) + 24 * ff(:,i,5) ) / h_120
       do k = 1, npt
          write( uwave, '(f8.3,1x,4e16.8)' ) r(k), ff(1,i,k), ff(2,i,k), &
               dff(1,i,k), dff(2,i,k)
       end do
       if ( boundary_condition == 0 .or.  &
            boundary_condition == 2 .and. l == 2 * ( l / 2 ) .or. &
            boundary_condition == 3 .and. l /= 2 * ( l / 2 ) ) then
          write( uwave, '(f8.3,1x,4e16.8)' ) ( npt + 1 ) * h, 0.0_pr, 0.0_pr, &
               ( - 600 * ff(:,i,npt) + 600 * ff(:,i,npt-1) &
               - 400 * ff(:,i,npt-2) + 150 * ff(:,i,npt-3) &
               - 24 * ff(:,i,npt-4) ) / h_120
       end if
       if ( boundary_condition == 1 .or.  &
            boundary_condition == 2 .and. l /= 2 * ( l / 2 ) .or. &
            boundary_condition == 3 .and. l == 2 * ( l / 2 ) ) then
          write( uwave, '(f8.3,1x,4e16.8)' ) (npt+1)*h, ff(1,i,npt), &
               ff(2,i,npt), - dff(1,i,npt), - dff(2,i,npt)
       end if
       close(uwave)
    end do
    !
  end subroutine print_out_wf

  !
  !!
  !

  subroutine print_out_can( n0, n1, npt, it, l, j )
    implicit none
    integer, intent(in) :: n0, n1, npt
    integer, intent(in), optional :: it, l, j
    integer :: i, k
    character ( len = 3 ) :: name
    !
    do i = n0, n0 + n1 - 1
       if ( mecan(i) < 60.0_pr ) then
          if ( i < 10 ) write(name,'("00",i1)') i
          if ( i > 9 .and. i < 100 ) write(name,'("0",i2)') i
          if ( i > 99 ) write(name,'(i3)') i
          open( unit = uwave, file = out//'can'//name//trim(extn)//'.gfx', &
               status = 'unknown' )
          if ( present(it) .and. present(l) .and. present(j) ) &
               write( uwave, &
               '("# it = ",i1,",  l = ",i2,",  j = ",i2,"/2")' ) it, l, j
          write( uwave, '("# E = ",f12.8)' ) ecan(i)
          write( uwave, '("#",/,"#  r (fm)        w.f.")' )
          write( uwave, '("# -------    -------------")' )
          write( uwave, '(f8.3,e18.7)' ) 0.0_pr, 0.0_pr
          do k = 1, npt
             write( uwave, '(f8.3,e18.7)' ) k * h, canwf(i,k)
          end do
          close(uwave)
       end if
    end do
    !
  end subroutine print_out_can

  !
  !!
  !

end module io
!
!==================================================================
!
!
!=========================================================================
!
module fields
  use cste
  use param
  implicit none

  interface init_sw
     module procedure init_sw
  end interface

  interface reduce_pairing
     module procedure reduce_pairing
  end interface

contains

  !
  !!
  !

  subroutine init_sw( ma, jz, n )
    !-----
    !  Determination of the Woods Saxon potentials for starting
    !  point in the first iteration
    !-----
    ! ma = total number of particle
    ! jz = proton number
    ! n = number of mesh points
    !-----
    implicit none
    integer, intent(in) :: ma, jz, n
    real ( kind = pr ) :: r0, v1, v2, aa, ab, ac, ad, ae, af, xx
    real ( kind = pr ) :: rad, rr, e, f, z, rc, sym, hr
    real ( kind = pr ), dimension(:), allocatable :: dp
    integer :: i, nr, ifail
    !
    if ( memlog ) &
         write( ulog, '(" INIT_SW()")' )
    allocate( dp(n), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory... (fields)'
    !
    select case (pot)
       !
    case ( "" )
       !
       r0 = 1.25_pr
       v1 = 65.0_pr
       v2 = 31.0_pr
       aa = 0.85_pr
       ab = 0.75_pr
       ac = 0.6_pr
       ae = 0.0_pr
       af = 3.0_pr
       xx = 2.0_pr**(1.0_pr/ac)
       ad = 2.0_pr * ab * ac * ( xx - 1.0_pr ) / xx
       rad = r0 * real( ma, kind = pr )**(1.0_pr/3.0_pr)
       rr = rad - ad * log( 2.0_pr**(1.0_pr/ac) - 1.0_pr )
       rc = rad * 1.09_pr / r0
       sym = real( ma - 2 * jz, kind = pr ) / ma
       !
       !! Coulomb potential
       !
       dp = 0.0_pr
       where ( abs( ( r - rc ) / 0.55_pr ) < maxa )
          dp = 1.0_pr / ( 1.0_pr + exp( ( r - rc ) / 0.55_pr ) )
       end where
       e = 0.0_pr
       f = 0.0_pr
       z = sum( r2 * dp )
       dp = dp * real( jz, kind = pr ) / ( qp * z * h )
       do i = 1, n
          e = e + qp * r2(i) * dp(i)
          f = f + qp * r(i) * dp(i)
          vc(i) = e / r(i) - f
       end do
       vc = echarg * h * ( vc + f )
       !
       !!
       !
       v1f = 0.0_pr
       v3f = 0.0_pr
       where ( abs( ( r - rr ) / ad ) < maxa )
          v1f(:,1) = &
               - 1.0_pr / ( 1.0_pr + exp( ( r - rr ) / ad ) )**ac &
               * ( v1 - v2 * sym )
          v1f(:,2) = &
               - 1.0_pr / ( 1.0_pr + exp( ( r - rr ) / ad ) )**ac &
               * ( v1 + v2 * sym )
          v3f(:,1) = hb * aa * ac * exp( ( r - rr ) / ad ) &
               / ( ad * ( 1.0_pr + exp( ( r - rr ) / ad ) )**(1.0_pr+ac) )
          v3f(:,2) = v3f(:,1)
       end where
       !
       v1f(:,2) = v1f(:,2) + vc(:)
       v2f(:,:) = hb
       !
       !! pairing fields
       !
       v1p = 0.0_pr
       v2p = 0.0_pr
       v3p = 0.0_pr
       where ( abs( ( r - rr ) / ad ) < maxa )
          v1p(:,1) = - del(1) / ( 1.0_pr + exp( ( r - rr ) / ad ) )**ac
          v1p(:,2) = - del(2) / ( 1.0_pr + exp( ( r - rr ) / ad ) )**ac
       end where
       !
       !! Effective mass
       !
       v12f = 0.0_pr
       v22f = 0.0_pr
       v12p = 0.0_pr
       v22p = 0.0_pr
       !
    case default
       !
       print '(/,19x,"*** READ POTENTIALS IN ",a," ***",/)', trim(pot)
       open( unit = uwave, file = potdir//pot, status = 'old', &
            form = 'unformatted' )
       read(uwave) hr
       if ( h /= hr ) then
          print '(" *")'
          print '(" *  The integration step h is not the same in the input")'
          print '(" *  file and in ",a)', trim(pot)
          print '(" *",/)'
          stop
       end if
       read(uwave) nr, del, amb
       read(uwave)   &
            ( v1f(i,1), v2f(i,1), v3f(i,1), &
            v1p(i,1), v2p(i,1), v3p(i,1),   &
            v12f(i,1), v22f(i,1), v12p(i,1), v22p(i,1), i = 1, min( nr, n ) )
       read(uwave) &
            ( v1f(i,2), v2f(i,2), v3f(i,2), &
            v1p(i,2), v2p(i,2), v3p(i,2), vc(i), &
            v12f(i,2), v22f(i,2), v12p(i,2), v22p(i,2), i = 1, min( nr, n ) )
       close(uwave)
       if ( n > nr ) then
          do i = nr + 1, n
             v1f(i,:)  = v1f(nr,:)
             v2f(i,:)  = v2f(nr,:)
             v3f(i,:)  = v3f(nr,:)
             v1p(i,:)  = v1p(nr,:)
             v2p(i,:)  = v2p(nr,:)
             v3p(i,:)  = v3p(nr,:)
             v12f(i,:) = v12f(nr,:)
             v22f(i,:) = v22f(nr,:)
             v12p(i,:) = v12p(nr,:)
             v22p(i,:) = v22p(nr,:)
             vc(i) = vc(nr)
          end do
       end if
       !
       if ( del(1) == 0.0_pr .and. bogolyubov(1) ) then
          del(1) = 0.5_pr
          v1p(:,1) = 0.03_pr * v1f(:,1)
       end if
       if ( del(2) == 0.0_pr .and. bogolyubov(2) ) then
          del(2) = 0.5_pr
          v1p(:,2) = 0.03_pr * v1f(:,2)
       end if
       !
    end select
    !
    deallocate(dp)
    !
    where ( abs(v1p)  < spacing(10.0_pr) ) v1p  = 0.0_pr
    where ( abs(v2p)  < spacing(10.0_pr) ) v2p  = 0.0_pr
    where ( abs(v3p)  < spacing(10.0_pr) ) v3p  = 0.0_pr
    where ( abs(v12p) < spacing(10.0_pr) ) v12p = 0.0_pr
    where ( abs(v22p) < spacing(10.0_pr) ) v22p = 0.0_pr
    where ( r > 18.0_pr )
       v1p(:,1) = 0.0_pr
       v2p(:,1) = 0.0_pr
       v3p(:,1) = 0.0_pr
       v12p(:,1) = 0.0_pr
       v22p(:,1) = 0.0_pr
       v1p(:,2) = 0.0_pr
       v2p(:,2) = 0.0_pr
       v3p(:,2) = 0.0_pr
       v12p(:,2) = 0.0_pr
       v22p(:,2) = 0.0_pr
    end where
    !
  end subroutine init_sw

  !
  !!
  !

  subroutine reduce_pairing( it, a )
    implicit none
    integer, intent(in) :: it, a
    !
    v1p(:,it) = v1p(:,it) / ( 0.2_pr * a + 12.0_pr )
    v3p(:,it) = v3p(:,it) / ( 0.2_pr * a + 12.0_pr )
    !
  end subroutine reduce_pairing

  !
  !!
  !

end module fields
!
!=========================================================================
!
!
!=========================================================================
!
module bogostatic
  !
  !! It seems better to define static (shared) arraies
  !! rather than local (automatic) ones. But this can be different
  !! on other architectures and/or with other compilers.
  !
  use cste
  real ( kind = pr ), dimension(:), allocatable :: vf, vp, am, ab, am1, &
       ab1, am2, ab2, d, v, w, x, y, v1, v2, cosa, sina,                &
       cosa2d, sina2d, f1, f2, h12d, h12v, h12w, det1, det2,            &
       eh12d, edet2
  real ( kind = pr ), dimension(:), allocatable :: valp, rho_d
  real ( kind = pr ), dimension(:,:), allocatable :: vecp
  integer, dimension(:), allocatable :: sib

  interface allocatebogo
     module procedure allocatebogo
  end interface

contains

  subroutine allocatebogo(n)
    implicit none
    integer, intent(in) :: n
    integer :: ifail
    !
    if ( memlog ) &
         write( ulog, '(" Allocate memory for the subroutine bogo().")' )
    allocate( vf(n), vp(n), am(n), ab(n), am1(n), ab1(n), am2(n), &
         ab2(n), d(n), v(n), w(n), x(n), y(n), v1(n), v2(n),      &
         cosa(n), sina(n), cosa2d(n), sina2d(n), f1(n), f2(n),    &
         h12d(n), h12v(n), h12w(n), det1(n), det2(n), eh12d(n),   &
         edet2(n), valp(n), rho_d(n), sib(n), vecp(n,n),          &
         stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory... (bogo)'
    !
  end subroutine allocatebogo

  !
  !!
  !

  subroutine deallocatebogo()
    implicit none
    !
    if ( memlog ) &
         write( ulog, &
         '(" Deallocate memory used in the subroutine bogo().")' )
    deallocate( vf, vp, am, ab, am1, ab1, am2,  &
         ab2, d, v, w, x, y, v1, v2,            &
         cosa, sina, cosa2d, sina2d, f1, f2,    &
         h12d, h12v, h12w, det1, det2, eh12d,   &
         edet2, valp, rho_d, sib, vecp )
    !
  end subroutine deallocatebogo

  !
  !!
  !

end module bogostatic

!
!! Module bogo
!

module bogo

  use cste
  use param
  use skforces
  use fields
  implicit none

  real ( kind = pr ), save :: ehfbmax(2)

  interface skyhfb
     module procedure skyhfb
  end interface

  interface quabcs
     module procedure quabcs
  end interface

  interface sfeden
     module procedure sfeden
  end interface

  interface print_out_fields
     module procedure print_out_fields
  end interface

contains

  subroutine skyhfb( it, l, j, npt, alamb, index, nfcttot, nfound, i0, &
       fin, ninteg, ifail )
    use eqdif
    use linalg
    use bogostatic
    implicit none
    integer, intent(in) :: it, l, j ! quantum numbers
    integer, intent(in) :: npt      ! size of the automatic arrays
    integer, intent(in) :: nfcttot  ! tot number of wave-fct found upto now
    real ( kind = pr ), intent(inout) :: alamb ! Fermi energy
    integer, intent(in) :: index    ! index of the (it,l,j) block
    !
    !! e_mid is an array containing the energies obtained during
    !! the previous iteration, deltae is a real number. If the brutal
    !! method is not used, the search for the new solutions will
    !! be made by starting at energies e_mid() - deltae instead
    !! of a systematic and time consuming search.
    !
    !! number wave functions found in this block
    !! ninteg is the number of times the subroutine "integrate" is called
    integer, intent(out) :: nfound
    integer, intent(inout) :: ninteg
    integer, intent(in) :: i0
    logical, intent(in) :: fin
    !
    !! Flag raised if iterations on the energy don't progress
    !! ifail = 1: node number disagreament
    !
    integer, intent(inout) :: ifail
    !
    !!----------------------------------------------------------------
    !
    ! This subroutine solves the system of two linear differential
    ! HFB equations:
    !   (-d/dr M d/dr + VV ) G = E G
    ! where M and VV are 2x2 matrices and G is a spinor wave function:
    !       / am  ab \         / vf-lambda    vp      \        / g1 \
    !   M = |        |    VV = |                      |    G = |    |
    !       \ ab -am /         \   vp      -vf+lambda /        \ g2 /
    ! By means of the linear transformation G = A F, the system above
    ! id transformed into a system without derivative in the coupling
    ! terms:
    !    ( -d2/dr2 + v ) * f1 + w * f2 = E/d * f1
    !    w * f1 - ( -d2/dr2 + v ) * f2 = E/d * f2
    ! and the latter is solved by using the Numerov and the systematic
    ! search for the eigenvalues.
    !
    !!----------------------------------------------------------------
    !
    integer :: ll1
    real ( kind = pr ) :: fspin, e, e0, fall1, fall2, vdepth, vdepth0, &
         epsene
    real ( kind = pr ) :: ecut
    integer :: nmatch
    ! Dummy variables
    integer :: i, k, m, ncur, tot_nodes, last_fct_nodes
    integer :: iter, n0, nrac, nr0(1), nrm(1)
    real ( kind = pr ) :: estep, det, det0, e00, norm, err,         &
         ecano, dcano, laste, e1, det01, det00, emax, eeq,          &
         xe0, ye0, xdet0, ydet0, l0
    !
    !! Various flags
    !
    logical :: flag, ipart, smaller_step, forcedic, jfo
    !
    !
    !! Brutal flag:
    !
    !! Flags used to choose the method for searching the solutions
    !! if .true., then a systematic search from E=0 is made with
    !! and energy step given by estep,
    !! if .false., the the code starts from the energies found
    !! during the previous iteration minus "something".
    !! Each time a "Node number disagreement" accident occurs,
    !! the flag is switched to .true., it is .flase. otherwise.
    !
    !...................................................................
    !
    !! Starting point for integration (to avoid the sigularity at
    !! the origin)
    !
    ll1 = l * ( l + 1 )
    n0 = int( sqrt( ll1 / 30.0_pr / h ) ) + 1
    !
    !! Different initial condition for l = 1
    !
    if ( l /= 1 ) then
       l0 = 0.0_pr
    else
       l0 = - 1.0_pr / 6.0_pr
    end if
    !
    nfound = 0 ! number of functions found in this (l,j,it) block
    estep = estep0
    fspin = ( ll1 - j * ( j + 2 ) / 4.0_pr + 0.75_pr )
    !
    !! Matching points
    !
    if ( it == 1 ) then
       !
       nrm = maxloc( abs( v1f(:,it) ) )
       nr0 = minloc( abs( v1f(nrm(1):,it)  &
            + abs( v1f(nrm(1),it) / 2.1_pr ) ) )
       nrac = nr0(1) + nrm(1) + 1
       nmatch = int( 1.2_pr * l / h )
    else
       !
       nrm = maxloc( abs( v1f(:,it) ) )
       nr0 = minloc( abs( v1f(nrm(1):,it)  &
            + abs( v1f(nrm(1),it) / 2.0_pr ) ) )
       nrac = nr0(1) + nrm(1) - 1
       nmatch = int( 1.3_pr * l / h )
    end if
    if ( nmatch < nrac     ) nmatch = nrac
    if ( nmatch > 2 * nrac ) nmatch = 2 * nrac
    if ( nmatch > npt - 4  ) nmatch = npt - 4
    !
    !! Calculation of the effectives potentials and inertia parameters
    !
    call set_veff( it, ll1, fspin, alamb )
    !
    epsene = 2 * eps
    forcedic = .false.
    vdepth0 = - min( minval(v1f(:,it)), 0.0_pr )
    select case (i_cut)
       !
    case (1)
       !
       !! Depth of the effective potential for this
       !! isopsin and partial wave
       !
       vdepth = - min( minval(vf), 0.0_pr )
       ecut = vdepth + cut_off
       !
    case (2)
       !
       !! Depth of the central potential for this isospin
       !
       vdepth = vdepth0
       ecut = max( vdepth0, cut_off + max( 1.2_pr * abs(alamb), 1.0_pr ) )
       !
    case (3)
       !
       !! ecut > cut_off + |lambda| since the cut off will be
       !! made in the equivalent energy spectrum
       !
       vdepth = vdepth0
       ecut = max( vdepth0, cut_off + max( 1.2_pr * abs(alamb), 1.0_pr ) )
       !
    case default
       !
       stop 'not implemented yet...'
       !
    end select
    !
    e = minval(vf)
    if ( e < 0.0_pr ) e = 0.0_pr
    e0 = 0.0_pr
    det00 = 0.0_pr
    det0 = 0.0_pr
    ipart = .false.
    flag = .false.
    smaller_step = .false.
    laste = -100.0_pr
    emax = vdepth
    tot_nodes = 0
    jfo = .false.
    !
    mainloop: do
       !
       if ( smaller_step .and. e < vdepth0 ) estep = estep0 / 150.0_pr
       !
       if ( .not.brutal(index) ) then
          if ( debug .or. ( l == 100 .and. it == 2 .and. j == 1 ) ) &
               print *, 'Search in intervalle previously determined...'
          det00 = 0.0_pr
          e0 = e_down( i0 + nfound, 2 )
          call set_boundary_values( l, e0, fall1, fall2 )
          call integrate( n0, l0, npt, eh12d, h12v, h12w, det1, edet2, &
               fall1, fall2, nmatch, l, det0 )
          ninteg = ninteg + 1
          e = e_up( i0 + nfound, 2 )
          call set_boundary_values( l, e, fall1, fall2 )
          call integrate( n0, l0, npt, eh12d, h12v, h12w, det1, edet2, &
               fall1, fall2, nmatch, l, det )
          ninteg = ninteg + 1
          det01 = det
          !
          if ( det0 * det01 > 0.0_pr ) then
             brutal(index) = .true.
             e = 0.0_pr
             nfound = 0
             det0 = 0.0_pr
          else
             xe0 = max( e_mid( i0 + nfound, 2 ) - estep / 100, &
                  e_down( i0 + nfound, 2 ) )
             call set_boundary_values( l, xe0, fall1, fall2 )
             call integrate( n0, l0, npt, eh12d, h12v, h12w, det1, edet2, &
                  fall1, fall2, nmatch, l, xdet0 )
             ninteg = ninteg + 1
             if ( debug .or. ( l == 100 .and. it == 2 .and. j == 1 ) ) &
                  print &
                  '(" E =",f14.9,3x,"Det =",e11.4,4x,"estep = ",&
                  &1f8.5,2x,l1)', &
                  e, det, estep, smaller_step
             ye0 = min( e_mid( i0 + nfound, 2 ) + estep / 100, &
                  e_up( i0 + nfound, 2 ) )
             call set_boundary_values( l, ye0, fall1, fall2 )
             call integrate( n0, l0, npt, eh12d, h12v, h12w, det1, edet2, &
                  fall1, fall2, nmatch, l, ydet0 )
             ninteg = ninteg + 1
             if ( debug .or. ( l == 100 .and. it == 2 .and. j == 1 ) ) &
                  print &
                  '(" E =",f14.9,3x,"Det =",e11.4,4x,"estep = ",&
                  &1f8.5,2x,l1)', &
                  e, det, estep, smaller_step
             if ( xdet0 * ydet0 < 0.0_pr ) then
                e0 = xe0
                det0 = xdet0
                e = ye0
                det = ydet0
             end if
             !
          end if
          !
       end if
       !
       if ( brutal(index) ) then
          if ( debug ) &
               print *, 'Brutal search !'
          !
          !! Initialization of the boundary values
          !
          call set_boundary_values( l, e, fall1, fall2 )
          !
          !! Numerical integration
          !
          call integrate( n0, l0, npt, eh12d, h12v, h12w, det1, edet2, &
               fall1, fall2, nmatch, l, det )
          ninteg = ninteg + 1
          if ( debug ) &
               print &
               '(" E =",f14.9,3x,"Det =",e11.4,4x,"estep = ",1f8.5,2x,2l1)', &
               e, det, estep, smaller_step, jfo
          if ( ( det - det0 ) * ( det0 - det00 ) < 0.0_pr .and. &
               .not.smaller_step .and. .not.jfo ) then
             if (debug) &
                  print *, ' changement de pente...', e
             if ( ( det0 > 0.0_pr .and. det > det0 ) .or. &
                  ( det0 < 0.0_pr .and. det < det0 ) ) then
                if (debug) &
                     print *, 'Possibilite de racine...'
                det = det00
                emax = e
                e = e - 3 * estep
                smaller_step = .true.
             end if
          end if
       end if
       !
       jfo = .false.
       !
       if ( det * det0 < 0.0_pr ) then
          iter = 0
          e00 = e
          e1 = e
          det01 = det
          do
             iter = iter + 1
             if ( iter == 150 ) then
                !
                !! If convergence is not progressing well, it seems better
                !! to switch back to dichotomy
                !
                epsene = 5 * epsene
                forcedic = .true.
                if ( debug ) &
                     print '(" Force dichotomy...")'
             end if
             if ( iter == 1000 ) then
                print *, 'Too many iterations...'
                print '(" ** it = ",i1,",  l = ",i2,",  j = ",i2,"/2")', &
                     it, l, j
                print '(" ** E=",f15.10,3x,"Det=",e12.4,&
                     &3x,"E0=",e12.4,3x,"E1=",e12.4,3x,"D(E)=",e10.3)',&
                     e, det, e0, e1, abs( ( e1 - e0 ) / e1 )
                ifail = 2
                if ( debug ) stop
                return
             end if
             !
             if ( abs( e1 - e0 ) > estep .or. forcedic ) then
                e = ( e0 + e1 ) / 2
             else
                e = ( det01 * e0 - det0 * e1 ) / ( det01 - det0 )
             end if
             !
             if ( abs( ( e00 - e ) / e ) < epsene .or. det01 == det0 ) exit
             call set_boundary_values( l, e, fall1, fall2 )
             call integrate( n0, l0, npt, eh12d, h12v, h12w, det1, edet2, &
                  fall1, fall2, nmatch, l, det )
             ninteg = ninteg + 1
             if ( det * det0 < 0.0_pr ) then
                e1 = e
                det01 = det
             else
                e0 = e
                det0 = det
             end if
             e00 = e
             if ( debug .or. ( l == 100 .and. it == 2 .and. j == 1 ) ) &
                  print '(" ** E =",f19.14,3x,"Det =",e18.10,&
                  &3x,"E0 =",e18.8,3x,"E1 =",e18.8,3x,"D(E)",e12.4)',&
                  e, det, e0, e1, abs( ( e1 - e0 ) / e1 )
          end do
          !
          !! A solution has been found, ask to build it...
          !
          if ( debug ) &
               print *, "Solution found..."
          epsene = 5 * eps
          forcedic = .false.
          call set_boundary_values( l, e, fall1, fall2 )
          call integrate( n0, l0, npt, eh12d, h12v, h12w, det1, edet2,  &
               fall1, fall2, nmatch, l, det, f1, f2 )
          ninteg = ninteg + 1
          !
          !! Raise the exit flag
          !
          ipart = .true.
          !
          !! Check its accuracy (not used)
          !
          !call accuracy_check( v, w, d, e, f1, f2, it, l, j, nfound, &
          !     npt, alamb, f1ext0, f2ext0 )
          !
          !! store it...
          !
          nfound = nfound + 1
          if ( nfound >= nfobl(index) ) brutal(index) = .true.
          ncur = nfcttot + nfound
          if ( ncur > nalloc ) call reallocatefct(npt)
          ff(1,ncur,:) = cosa2d(:) * f1(:) - sina2d(:) * f2(:)
          ff(2,ncur,:) = sina2d(:) * f1(:) + cosa2d(:) * f2(:)
          !
          !! Norm it...
          !
          if ( debug ) print *, 'Norm the solution...'
          norm = sqrt( ( dot_product( ff(1,ncur,:), ff(1,ncur,:) ) &
               + dot_product( ff(2,ncur,:), ff(2,ncur,:) ) ) * h )
          ff(:,ncur,:) = ff(:,ncur,:) / norm
          !
          !! Check its node number...
          !
          if ( debug ) print *, 'Chech nodes number...'
          last_fct_nodes = count_nodes( ff(1,ncur,:), ff(2,ncur,:), npt )
          tot_nodes = tot_nodes + last_fct_nodes
          numnodes(ncur) = last_fct_nodes + 1
          !
          !! Calculate its derivative...
          !
          if ( debug ) print *, 'Compute derivative...'
          call derivative( npt, ff(:,ncur,:), dff(:,ncur,:), l )
          ehfb(ncur) = e
          !
          !! Occupation of the lower component
          !
          if ( debug ) print *, 'Compute occupation factor...'
          qpv2(ncur) = sum( ff(2,ncur,:) * ff(2,ncur,:) ) * h
          !
          !! Radius
          !
          if ( debug ) print *, '      ... radius...'
          qprad(ncur) = 0.0_pr
          if ( qpv2(ncur) /= 0.0_pr ) qprad(ncur) = &
               sqrt( sum( ff(2,ncur,:)**2 * r2 * h ) / qpv2(ncur) )
          !
          !! Set next energy and reset determinant
          !
          if ( debug ) print *, ' And search for the next solution...'
          estep = estep0
          if ( e < vdepth0 ) then
             e = e + estep / 1000.0_pr
          else
             e = e + estep / 5.0_pr
          end if
          eh12d = e * h12d
          edet2 = e * e * det2
          if ( debug .or. ( l == 100 .and. it == 2 .and. j == 1 ) ) &
               print '("Integration just above the previous solution:")'
          call set_boundary_values( l, e, fall1, fall2 )
          call integrate( n0, l0, npt, eh12d, h12v, h12w, det1, edet2, &
               fall1, fall2, nmatch, l, det )
          if ( debug .or. ( l == 100 .and. it == 2 .and. j == 1 ) ) &
               print '("     E =",f14.6,4x,"Det = ",e18.6)', e, det
          ninteg = ninteg + 1
          jfo = .true.
       end if
       !
       !! Restore the original value for estep when we are out of
       !! the region where the slope changes
       !
       if ( e >= emax ) then
          smaller_step = .false.
          estep = estep0
       end if
       e0 = e
       if (brutal(index)) det00 = det0
       det0 = det
       e = e + estep
       !
       if ( e >= ecut .and. ipart ) then
          if ( debug ) print *, ' Enough functions have been found !...'
          !
          !! Check the number of nodes of the last function
          !
          !   if ( tot_nodes /= sum( (/ ( i, i = 0, nfound - 1 ) /) ) ) then
          !      ! Say nothing, but worries a little...
          !      if ( nfound - 1 /= last_fct_nodes ) then
          !         print '(" Node number disagreement !",&
          !              &1x,"it = ",i1,"  l = ",i2,"  j = ",i2,"/2")', &
          !              it, l, j
          !         print '(5x,"nrac = ",i3,4x,"nmatch = ",i3)', nrac, nmatch
          !         ifail = 1
          !         estep0 = estep0 / 2
          !         if ( estep0 < estep0min / 4 ) estep0 = estep0min
          !         brutal(index) = .true.
          !      end if
          !   end if
          !
          !! Remarque: As it is discussed in the CPC article, the number of
          !! nodes is sensitive to the "modified asymptotic behaviour"
          !! of the HFB spinors. So this test is meaningless for weakly
          !! bound nuclei.
          !
          if ( debug ) print *, ' Exit main loop...    *****'
          exit
       end if
       if ( debug ) print *, ' End of main loop...  *****'
    end do mainloop
    !
    e_down(nfcttot+1,1) = 0.0_pr
    e_up(ncur,1) = e + 5 * estep
    if ( ncur > nfcttot + 1 ) then
       e_down(nfcttot+2:ncur,1) = &
            ( ehfb(nfcttot+1:ncur-1) + ehfb(nfcttot+2:ncur) ) / 2
       e_up(nfcttot+1:ncur-1,1) = &
            ( ehfb(nfcttot+1:ncur-1) + ehfb(nfcttot+2:ncur) ) / 2
    end if
    e_mid(nfcttot+1:ncur,1) = ehfb(nfcttot+1:ncur)
    brutal(index) = ( nfobl(index) /= nfound )
    !
    !! One particle canonical state blocking comes here ...
    !
    ! ... but it is not implemented...
    !
    if ( fin .and. cano ) then
       !
       !! Build density
       !
       if ( debug ) print *, ' Build density...'
       do k = 1, npt
          do i = 1, npt
             select case (i_cut)
             case (1:2)
                density(i,k) = dot_product( ff(2,nfcttot+1:ncur,i), &
                     ff(2,nfcttot+1:ncur,k) )
             case (3)
                density(i,k) = ff(2,nfcttot+1,i) * ff(2,nfcttot+1,k)
                do m = nfcttot + 2, ncur
                   eeq = alamb + ( 1 - 2 * qpv2(m) ) * ehfb(m)
                   if ( eeq > cut_off ) exit
                   density(i,k) = density(i,k) + ff(2,m,i) * ff(2,m,k)
                end do
             end select
          end do
       end do
       density = density * h
       !
       !! Diagonalize it
       !
       if ( debug ) print *, 'Last iteration, diagonalize density...'
       call diagon( npt, density, valp, vecp )
       !
       !! Blocking    ?
       !
       ! .............. to do.
       !
       !! Norm the canonical states and test the accuracy
       !
       do i = 1, npt
          vecp(:,i) = vecp(:,i) &
               / sqrt( dot_product( vecp(:,i), vecp(:,i) ) * h )
          err = maxval( matmul( density, vecp(:,i) ) - valp(i) * vecp(:,i) )
          if ( err > 20000 * eps ) then ! There is no need in asking
             !                          ! too high precision...
             print '("Bad accuracy in canonical states, ERR =",e12.4)', err
          end if
       end do
       !
       !! Diagonal matrix element of the particle and pairing fields
       !! in the canonical basis
       !
       do i = 1, nfound
          canwf(nfcttot+i,:) = vecp(:,npt-i+1) ! save canonical wave functions
          call derivative( npt, vecp(:,npt-i+1), rho_d(:), l )
          v2can(nfcttot+i) = valp(npt-i+1)
          ! avoid Floating underflow :
          where ( abs(vecp(:,npt-i+1)) < crtiny ) vecp(:,npt-i+1) = 0.0_pr
          !
          ecano = sum( rho_d(:)**2 * am(:) + vecp(:,npt-i+1)**2 * vf(:) ) * h
          dcano = sum( rho_d(:)**2 * ab(:) + vecp(:,npt-i+1)**2 * vp(:) ) * h
          if ( l == 0 ) then
             ! Correction of integrals du to a nonvanishing
             ! contributions of kinetic terms at the origin
             ecano = ecano + ( 4 * rho_d(1) - rho_d(2) )**2 &
                  * ( 4 * am(1) - am(2) ) / 54.0_pr * h
             dcano = dcano + ( 4 * rho_d(1) - rho_d(2) )**2 &
                  * ( 4 * ab(1) - ab(2) ) / 54.0_pr * h
          end if
          mecan(nfcttot+i) = ecano
          ecano = ecano - alamb
          dcano = dcano
          ecan(nfcttot+i) = sqrt( ecano * ecano + dcano * dcano )
       end do
    end if
    !
    if ( debug ) print *, 'End of SKYHFB...'
    !
  end subroutine skyhfb

  !
  !!
  !

  subroutine set_veff( it, ll1, fspin, alamb )
    !
    !! The effective potential for a given partial wave is set here.
    !
    use bogostatic
    implicit none
    integer, intent(in) :: it, ll1
    real ( kind = pr ), intent(in) :: fspin, alamb
    !
    vf(:)  = v1f(:,it) + v2f(:,it) * ll1 / r2 + v3f(:,it) * fspin / r
    am(:)  = v2f(:,it)
    am1(:) = v12f(:,it)
    am2(:) = v22f(:,it)
    vp(:)  = v1p(:,it) + v2p(:,it) * ll1 / r2 + v3p(:,it) * fspin / r
    ab(:)  = v2p(:,it)
    ab1(:) = v12p(:,it)
    ab2(:) = v22p(:,it)
    d = sqrt( am * am + ab * ab )
    x = ( am * ab1 - ab * am1 ) / 4.0_pr  / d / d
    y = ( am * am1 + ab * ab1 ) / 4.0_pr  / d / d
    v1 = vf - alamb + am2 / 2 + ab1 * x - am1 * y
    v2 = vp         + ab2 / 2 - am1 * x - ab1 * y
    cosa = am / d
    sina = ab / d
    v = ( v1 * cosa + v2 * sina ) / d
    w = ( v2 * cosa - v1 * sina ) / d
    sib = 1
    where ( ab < 0.0_pr ) sib = -1
    where ( ab == 0.0_pr )
       cosa2d = 1.0_pr / sqrt(d)
       sina2d = 0.0_pr
       elsewhere
       cosa2d = sqrt( ( 1.0_pr + am / d ) / 2.0_pr / d )
       sina2d = 0.0_pr
       where ( 1.0_pr - am / d > 0.0_pr ) &
            sina2d = sqrt( ( 1.0_pr - am / d ) / 2.0_pr / d ) * sib
    end where
    h12d = h12 / d
    h12v = 1.0_pr - h12 * v
    h12w = h12 * w
    det1 = h12v**2 + h12w**2
    det2 = h12d**2
    !
  end subroutine set_veff

  !
  !!
  !

  subroutine set_boundary_values( l, e, fall1, fall2 )
    !
    use bogostatic, only : eh12d, edet2, det2, h12d
    implicit none
    integer, intent(in) :: l
    real ( kind = pr ), intent(in) :: e
    real ( kind = pr ), intent(out) :: fall1, fall2
    !
    eh12d = e * h12d
    edet2 = e * e * det2
    !
    if ( boundary_condition == 0 .or.  &
         boundary_condition == 2 .and. l == 2 * ( l / 2 ) .or. &
         boundary_condition == 3 .and. l /= 2 * ( l / 2 ) ) then
       fall1 = 0.0_pr
       fall2 = 0.0_pr
    end if
    if ( boundary_condition == 1 .or.  &
         boundary_condition == 2 .and. l /= 2 * ( l / 2 ) .or. &
         boundary_condition == 3 .and. l == 2 * ( l / 2 ) ) then
       fall1 = 1.0_pr
       fall2 = 1.0_pr
    end if
    !
  end subroutine set_boundary_values

  !
  !!
  !

  subroutine quabcs( fin )
    !
    !! Calculation of the "equivalent" single-particle energies
    !! and state-dependent pairing gaps
    !
    implicit none
    logical, intent(in) :: fin
    integer :: it, l, j ! quantum numbers
    real ( kind = pr ) :: alamb
    integer :: i, n, ll, i0, i1, n0, n1, np, nn, k
    integer, parameter :: itdichomax = 150
    integer :: itdicho
    real ( kind = pr ) :: mqv2e, ehfbmin, x, y, a, b, vh
    real ( kind = pr ) :: xinf, xsup, einf, esup
    real ( kind = pr ) :: degeneracy, e, f, epa
    real ( kind = pr ), parameter :: epseqp = 1.e-4_pr, &
         epslam = 1.e-7_pr, epspar = 1.e-9_pr
    logical :: icze
    !
    if ( debug ) print *, 'Step in QUABCS...'
    !
    ehfbmax = 0.0_pr
    !
    nn = sum(nfobl(1:j_max(1)+1))                   ! number of neutron states
    np = sum(nfobl(j_max(1)+2:j_max(1)+j_max(2)+2)) ! idem for protons
    !
    !! Loop over isospin
    !
    do it = 1, 2
       alamb = amb(it)
       if ( it == 1 ) then         ! Boundary values of the loops
          i0 = 1                   !
          i1 = nn                  !
          n0 = 1                   !
          n1 = j_max(it) + 1       !
       else                        !
          i0 = nn + 1              !
          i1 = nn + np             !
          n0 = j_max(1) + 2        !
          n1 = sum(j_max) + 2      !
       end if                      !
       do i = i0, i1
          mqv2e = ( 1 - 2 * qpv2(i) ) * ehfb(i)
          ehf(i) = alamb + mqv2e
          dhf(i) = sqrt( abs( ehfb(i) * ehfb(i) - mqv2e * mqv2e ) )
       end do
       !
       !! Compute the Fermi energy for the next iteration
       !! (see appendix B in the DFT paper)
       !
       if ( debug ) print *, &
            '     Compute the Fermi energy for the next iteration...'
       ehfbmin = minval(ehfb(i0:i1))
       x = alamb
       ehfbmax(it) = maxval( ehfb(i0:i1) )
       icze = .false.
       xinf = -100.0_pr
       xsup = 80.0_pr
       einf = 0.0_pr
       esup = 1.0_pr
       itdicho = 0
       epa = 0.0_pr

       do
          e = 0.0_pr
          f = 0.0_pr
          i = i0
          j = 1
          ll = 1
          l = 0
          st1: do n = n0, n1
             degeneracy = j + 1
             do i = 1, nfobl(n)
                vh = 0.0_pr
                k = sum(nfobl(1:n)) - nfobl(n) + i
                y = ehf(k) - x
                a = y * y + dhf(k) * dhf(k)
                b = sqrt(a)
                if ( b > 0.0_pr ) vh = 0.5_pr * ( 1.0_pr - y / b )
                !
                !! This statement enforces the partial occupation of the
                !! HF orbit when there is no pairing
                !
                if ( icze .and. (.not.bogolyubov(it)) .and. b < epseqp ) &
                     vh = - einf / ( esup - einf )
                if ( vh < 0.0_pr ) vh = 0.0_pr
                if ( vh > 1.0_pr ) vh = 1.0_pr
                if ( cut_diffuseness /= 0.0_pr ) then
                   !
                   !! Smoothed cut off
                   !
                   if ( i_cut == 3 ) &
                        vh = vh / ( 1 + exp( ( ehf(k) - cut_off )  &
                        / cut_diffuseness ) )
                else
                   !
                   !! or sharp
                   !
                   if ( ehf(k) > cut_off .and. i_cut == 2 ) vh = 0.0_pr
                   !
                end if
                !
                !! Here is the effect of one quasiparticle blocking
                !
                if ( it + 2 == iing(it)  & ! 1 qp is blocked in this fluid
                     .and. i == ning(it) & ! that is this qp
                     .and. l == ling(it) & ! with l and 2*j
                     .and. j == jing(it) ) then
                   e = e + vh * ( degeneracy - 1 ) + 1 - vh
                   if ( b > 0.0_pr ) f = f + dhf(k) * dhf(k) / a / b &
                        * ( degeneracy - 2 )
                   !     !
                else     !! if there is no blocking
                   !     !
                   e = e + vh * degeneracy
                   if ( b > 0.0_pr ) f = f + dhf(k) * dhf(k) / a / b &
                        * degeneracy
                end if
                vhf(k) = vh
             end do
             ll = ll + 2
             if ( ll == 5 ) then
                j = j + 2
                ll = 1
             end if
             l = ( j - 2 + ll ) / 2
          end do st1
          f = f / 2
          e = e - npr(it)
          if ( abs( e / npr(it) ) - epspar <= 0.0_pr .and.  &
               itdicho > 0 ) exit
          !
          !! The occupation probability of the blocked state
          !! is taken into account
          !
          if ( e == 0.0_pr ) exit
          if ( e < 0.0_pr ) then
             xinf = max( xinf, x )
             einf = e
          else
             xsup = min( xsup, x )
             esup = e
          end if
          x = x - e / ( f + 1.e-20_pr )
          if ( xsup - xinf <= epslam ) icze = .true.
          if ( x < xinf .or. x > xsup ) x = ( xinf + xsup ) / 2
          !
          if ( it == dripline ) exit
          !
          itdicho = itdicho + 1
          if ( itdicho == itdichomax ) exit ! sometimes it failes to converge !
          !
          if ( debug ) print *, ' x =', x
          !
       end do
       amb(it) = x
       !
       if ( it == dripline ) amb(it) = 0.0_pr
       !
    end do ! it
    if (fin) call print_out_qp()
    !
  end subroutine quabcs

  !
  !!
  !

  subroutine print_out_qp()
    !
    !! Print out the summary of quasiparticle states properties
    !
    implicit none
    integer :: it, l, j ! quantum numbers
    integer :: i, n, ll, i0, i1, n0, n1, np, nn, k, degeneracy
    !
    nn = sum(nfobl(1:j_max(1)+1))            ! number of neutron states
    np = sum(nfobl(j_max(1)+2:sum(j_max)+2)) ! idem for protons
    !
    !! Loop over isospin
    !
    do it = 1, 2
       !
       if ( it == 1 ) then         ! Boundary values of the loops
          i0 = 1                   !
          i1 = nn                  !
          n0 = 1                   !
          n1 = j_max(1) + 1        !
       else                        !
          i0 = nn + 1              !
          i1 = nn + np             !
          n0 = j_max(1) + 2        !
          n1 = sum(j_max) + 2      !
       end if                      !
       !
       if ( it == 1 ) then
          write( uspe, '("  Fermi energies (MeV):",/,2x,21("-"))' )
          write( uspe, &
               '(5x,"lambda_n = ",f10.6,/,5x,"lambda_p = ",f10.6,/)' ) amb
       end if
       write( uspe, '(/,32x,27("-"))' )
       if ( it == 1 ) then
          write( uspe, '(32x,"|  Spectrum for neutrons  |")' )
       else
          write( uspe, '(32x,"|  Spectrum for protons   |")' )
       end if
       write( uspe, '(32x,27("-"),/)' )
       write( uspe, '(32x,"( lambda = ",f10.6," MeV )")' ) amb(it)
       write( uspe, '(//,34x,"Quasiparticles properties",/,34x,25("-"),/)' )
       j = 1
       ll = 1
       l = 0
       write( uspe, '(" it 2j   l   n    N",6x,"E_qp",6x,"N_qp", &
            &7x,"E_eq      D_eq      N_eq",7x,"r.m.s   nodes")' )
       write( uspe, '(90("-"))' )
       loopsummaryqp: do n = n0, n1
          degeneracy = j + 1
          do i = 1, nfobl(n)
             k = sum(nfobl(1:n)) - nfobl(n) + i
             !
             write( uspe, '(1x,i1,2(2x,i2),1x,i3,1x,i4,2x,f10.6,&
                  & 1x,f9.6,1x,f10.6, &
                  & 2(1x,f9.6),1x,f10.6,2x,i3)' ) &
                  it, j, l, i, k, &
                  ehfb(k), qpv2(k), ehf(k), dhf(k), vhf(k), qprad(k), &
                  numnodes(k)
             !
          end do
          ll = ll + 2
          if ( ll == 5 ) then
             j = j + 2
             if ( j > j_max(it) ) exit loopsummaryqp
             ll = 1
          end if
          l = ( j - 2 + ll ) / 2
       end do loopsummaryqp
       !
       !! fluctuation of particles numbers
       !
       dispersion(it) = 0.0_pr
       !
       !! Canonical states properties
       !
       if ( cano ) then
          write( uspe, &
               '(//,26x,"Canonical states properties",/,26x,27("-"),/)' )
          j = 1
          ll = 1
          l = 0
          write( uspe, '(" it 2j   l   n   N",6x,"E_can",6x,"v2", &
               &8x,"E_eq")' )
          write( uspe, '(51("-"))' )
          loopsummarycan: do n = n0, n1
             do i = 1, nfobl(n)
                k = sum(nfobl(1:n)) - nfobl(n) + i
                !
                if ( mecan(k) < cut_off ) then
                   write( uspe, &
                        '(1x,i1,2x,i2,2x,i2,2(1x,i3),2x,f10.6,1x,&
                        &f9.6,1x,f10.6)')&
                        it, j, l, i, k, ecan(k), v2can(k), mecan(k)
                end if
                !
                dispersion(it) = dispersion(it) &
                     + ( j + 1 ) * v2can(k) * ( 1.0_pr - v2can(k) )
             end do
             ll = ll + 2
             if ( ll == 5 ) then
                j = j + 2
                if ( j > j_max(it) ) exit loopsummarycan
                ll = 1
             end if
             l = ( j - 2 + ll ) / 2
          end do loopsummarycan
          !
          dispersion(it) = 2 * abs(dispersion(it))
          !
       end if
       !
    end do ! it
    !
  end subroutine print_out_qp

  !
  !!
  !

  subroutine sfeden( fin, npt )
    !
    implicit none
    logical, intent(in) :: fin
    integer, intent(in) :: npt
    !
    !! This subroutine sums up the densities
    !
    integer :: it, ll1, n0, n1, i0, i1, nn, np, i, n, k, ll, j, l
    integer :: degeneracy, p, ifail
    integer, dimension(:), allocatable :: compx, compy
    real ( kind = pr ), dimension(:), allocatable :: gv2, gu2, guv
    real ( kind = pr ) :: alamb, fspin, fu2, fv2, fuv
    real ( kind = pr ) :: x, y, z, t, xx, yy
    real ( kind = pr ), parameter :: epsano = 1.e-8_pr
    !
    !! For regularization
    !
    real ( kind = pr ) :: vf, kc, lc, ec, a, b
    real ( kind = pr ) :: kf
    !
    nn = sum(nfobl(1:j_max(1)+1))            ! number of neutron states
    np = sum(nfobl(j_max(1)+2:sum(j_max)+2)) ! idem for protons
    !
    allocate( compx(nn+np), compy(nn+np), gu2(nn+np), gv2(nn+np), &
         guv(nn+np), stat = ifail )
    if ( ifail /= 0 ) stop 'Not enough memory... (sfeden in bogo)'
    !
    !! Flush densities
    !
    rho   = 0.0_pr
    tau   = 0.0_pr
    cur   = 0.0_pr
    rho_p = 0.0_pr
    tau_p = 0.0_pr
    cur_p = 0.0_pr
    !
    !! Loop over isospin
    !
    do it = 1, 2
       alamb = amb(it)
       if ( it == 1 ) then         ! Boundary values of the loops
          i0 = 1                   !
          i1 = nn                  !
          n0 = 1                   !
          n1 = j_max(1) + 1        !
       else                        !
          i0 = nn + 1              !
          i1 = nn + np             !
          n0 = j_max(1) + 2        !
          n1 = sum(j_max) + 2      !
       end if                      !
       !
       j = 1
       ll = 1
       l = 0
       do n = n0, n1
          degeneracy = j + 1
          ll1 = l * ( l + 1 )
          fspin = ( ll1 - j * ( j + 2 ) / 4.0_pr + 0.75_pr )
          do i = 1, nfobl(n)
             k = sum(nfobl(1:n)) - nfobl(n) + i
             compx(k) = 2
             compy(k) = 1
             !
             !! ...
             !
             if ( qpv2(k) <= 1.e-7_pr .and. (.not.bogolyubov(it)) ) then
                compx(k) = 1
                fv2 = vhf(k)
                fu2 = 1.d0 - vhf(k)
             else
                if ( 1.0_pr - qpv2(k) <= 1.e-7_pr &
                     .and. (.not.bogolyubov(it)) ) then
                   compy(k) = 2
                   fv2 = vhf(k)
                   fu2 = 1.0_pr - vhf(k)
                else
                   !
                   !! ...
                   !
                   if ( qpv2(k) > epsano ) then
                      fv2 = vhf(k) / qpv2(k)
                   else
                      fv2 = 1.0_pr
                   end if
                   if ( 1.0_pr - qpv2(k) > epsano ) then
                      fu2 = ( 1.0_pr - vhf(k) ) / ( 1.0_pr - qpv2(k) )
                   else
                      fu2 = 1.0_pr
                   end if
                end if
             end if
             fuv = sqrt( fv2 * fu2 )
             if ( compx(k) == compy(k) ) fuv = - fuv
             gv2(k) = fv2
             gu2(k) = fu2
             guv(k) = fuv
             !
          end do ! i (k)
          !
          !! loop over mesh points
          !
          do p = 1, npt
             !
             !! loop over states in a given block
             !
             x = 0.0_pr
             y = 0.0_pr
             z = 0.0_pr
             t = 0.0_pr
             !
             do i = 1, nfobl(n)
                k = sum(nfobl(1:n)) - nfobl(n) + i
                if ( gv2(k) <= 0.0_pr ) cycle

                xx = dff( compx(k), k, p ) - ff( compx(k), k, p ) / r(p)
                yy = dff( compy(k), k, p ) - ff( compy(k), k, p ) / r(p)
                x = x + ff( compx(k), k, p ) * ff( compx(k), k, p ) &
                     * gv2(k) * degeneracy
                y = y + xx * xx * gv2(k) * degeneracy
                z = z - ff( compx(k), k, p ) * ff( compy(k), k, p ) &
                     * guv(k) * degeneracy
                t = t - xx * yy * guv(k) * degeneracy
                !
                !! One quasiparticle blocking
                !
                if ( it + 2  == iing(it) .and. i == ning(it) .and. &
                     l == ling(it) .and. j .eq. jing(it) ) then
                   x = x + ff( compy(k), k, p ) * ff( compy(k), k, p ) &
                        * gu2(k)                                       &
                        - ff( compx(k), k, p ) * ff( compx(k), k, p )  &
                        * gv2(k)
                   y = y + yy * yy * gu2(k) - xx * xx * gv2(k)
                   stop 'qp blocking + regularization'
                   z = z + ff( compx(k), k, p ) * ff( compy(k), k, p ) &
                        * guv(k) * 2
                   t = t + xx * yy * guv(k) * 2
                end if
             end do ! i (k)
             rho(p,it) = rho(p,it) + x * qpr2(p)
             tau(p,it) = tau(p,it) + ( y + ll1 * x / r2(p) ) * qpr2(p)
             cur(p,it) = cur(p,it) + x * qpr2(p) * fspin / r(p)
             if ( bogolyubov(it) ) then
                rho_p(p,it) = rho_p(p,it) + z * qpr2(p)
                tau_p(p,it) = tau_p(p,it) &
                     + ( t + ll1 * z / r2(p) ) * qpr2(p)
                cur_p(p,it) = cur_p(p,it) + z * qpr2(p) * fspin / r(p)
             end if
          end do ! p
          ll = ll + 2
          if ( ll == 5 ) then
             j = j + 2
             if ( j > j_max(it) ) exit
             ll = 1
          end if
          l = ( j - 2 + ll ) / 2
       end do
       !
       if ( regularization ) then
          !
          !!
          !! Regularization of the abnormal density:
          !!
          !
          ec = cut_off - amb(it)
          do p = 1, npt
             !
             vf = v1f(p,it)
             if ( vf + ec - amb(it) < 0.0_pr ) then
                !
                kf = sqrt( ( amb(it) - vf ) / v2f(p,it) )
                kc = sqrt( ( amb(it) + ec - vf ) / v2f(p,it) )
                lc = sqrt( ( amb(it) - ec - vf ) / v2f(p,it) )
                a = kc / 4 / pi**2 / v2f(p,it) &
                     * ( 1 - kf / 2 / kc * log( ( kc + kf ) / ( kc - kf ) ) )
                b = lc / 4 / pi**2 / v2f(p,it) &
                     * ( 1 - kf / 2 / lc * log( ( lc + kf ) / ( kf - lc ) ) )
                rega(p,it) = a
                regb(p,it) = b
                !........................ geff(p,it) = 1 / ( 1 / g - a - b )
             else
                if ( vf - amb(it) < 0.0_pr ) then
                   !
                   kf = sqrt( ( amb(it) - vf ) / v2f(p,it) )
                   kc = sqrt( ( amb(it) + ec - vf ) / v2f(p,it) )
                   a = kc / 4 / pi**2 / v2f(p,it) * &
                        ( 1 - kf / 2 / kc * log( ( kc + kf ) / ( kc - kf ) ) )
                   rega(p,it) = a
                   regb(p,it) = 0.0_pr
                   !..................... geff(p,it) = 1 / ( 1 / g - a )
                else
                   if ( vf - ec - amb(it) < 0.0_pr ) then
                      !
                      kf = sqrt( - ( amb(it) - vf ) / v2f(p,it) )
                      kc = sqrt( ( amb(it) + ec - vf ) / v2f(p,it) )
                      a = kc / 4 / pi**2 / v2f(p,it) &
                           * ( 1 + kf / kc * atan( kf / kc ) )
                      rega(p,it) = a
                      regb(p,it) = 0.0_pr
                      !.................. geff(p,it) = 1 / ( 1 / g - a )
                   else
                      !
                      rega(p,it) = 0.0_pr
                      regb(p,it) = 0.0_pr
                      !.................. geff(p,it) = g
                   end if
                end if
             end if
             !
          end do
       end if
       !
    end do ! it
    deallocate( compx, compy, gu2, gv2, guv )
    if ( fin ) call print_out_fields( npt )
    !
  end subroutine sfeden

  !
  !!
  !

  subroutine print_out_fields(npt)
    implicit none
    integer, intent(in) :: npt
    integer :: k
    !
    !! Print out mean-fields.
    !
    write( uspe, '(//,35x,"Neutrons mean fields")')
    write( uspe, '(35x,"--------------------")')
    write( uspe, '(16x,"Particle-hole channel",&
         &15x,"Particle-particle channel")')
    write( uspe, '(12x,35("-"),1x,35("-"))')
    write( uspe, '(" it    r      Average     Kinetic    Spin-orbit",&
         &"   Average     Kinetic    Spin-orbit")')
    write( uspe, '(" -- ",7("-"),1x,71("-"))' )
    do k = 1, npt
       write( uspe, '(2x,i1,1x,f6.2,1x,6(1pe12.4))' ) 1,  &
            r(k), v1f(k,1), v2f(k,1), v3f(k,1), v1p(k,1), v2p(k,1), v3p(k,1)
    end do
    if (meanfields) then
       open( unit = uwave, file = out//'neutron'//trim(extn)//'.mf', &
            status = 'unknown' )
       write( uwave, '("#                  Particle-hole channel       ",&
            &"            Particle-particle channel")')
       write( uwave, '("#        -----------------------------------------",&
            &"  ----------------------------------------")')
       write( uwave, '("#   r       Average       Kinetic      Spin-orbit",&
            &"     Average       Kinetic      Spin-orbit")')
       write( uwave, '("#-----------------------------------------",&
            &"--------------------------------------------------")')
       do k = 1, npt
          write( uwave, '(1x,f6.2,1x,6e14.6)' )    &
               r(k), v1f(k,1), v2f(k,1), v3f(k,1), &
               v1p(k,1), v2p(k,1), v3p(k,1)
       end do
       close(uwave)
    end if
    write( uspe, '(///)' )
    !
    !! ... and densities.
    !
    write( uspe, '(37x,"Neutrons densities")')
    write( uspe, '(37x,"------------------")')
    write( uspe, '(18x,"Particle-hole channel", &
         &15x,"Particle-particle channel")')
    if ( regularization ) then
       write( uspe, '(12x,34("-"),2x,46("-"))')
       write( uspe, '(" it   r      Average     Kinetic    Spin-orbit", &
            &"   Average     Kinetic    Spin-orbit    Reg.")')
       write( uspe, '(1x,93("-"))')
       do k = 1, npt
          write( uspe, '(1x,i1,1x,f6.2,1x,7(1pe12.4))' ) &
               1, r(k), rho(k,1), tau(k,1), cur(k,1),    &
               rho_p(k,1), tau_p(k,1), cur_p(k,1),       &
               rho_p(k,1) + 2 * v1p(k,1) * rega(k,1)
       end do
    else
       write( uspe, '(12x,34("-"),2x,34("-"))')
       write( uspe, '(" it   r      Average     Kinetic    Spin-orbit", &
            &"   Average     Kinetic    Spin-orbit")')
       write( uspe, '(1x,81("-"))')
       do k = 1, npt
          write( uspe, '(1x,i1,1x,f6.2,1x,6(1pe12.4))' ) &
               1, r(k), rho(k,1), tau(k,1), cur(k,1),    &
               rho_p(k,1), tau_p(k,1), cur_p(k,1)
       end do
    end if
    if (densities) then
       open( unit = uwave, file = out//'neutron'//trim(extn)//'.dens', &
            status = 'unknown' )
       if ( regularization ) then
          write( uwave, '("#   r   ",&
               &2("     Average       Kinetic      Spin-orbit"),&
               &"   Regularized    Regulator")' )
          write( uwave, '("# ------ ",8(" -------------"))' )
          do k = 1, npt
             write( uwave, '(1x,f7.3,1x,8e14.6)' )    &
                  r(k), rho(k,1),   tau(k,1),   cur(k,1),       &
                  &     rho_p(k,1), tau_p(k,1), cur_p(k,1),    &
                  &     rho_p(k,1) + 2 * v1p(k,1) * rega(k,1), &
                  &     rega(k,1)
          end do
       else
          write( uwave, '("#   r   ",&
               &2("     Average       Kinetic      Spin-orbit"))' )
          write( uwave, '("# ------ ",6(" -------------"))' )
          do k = 1, npt
             write( uwave, '(1x,f7.3,1x,6e14.6)' )           &
                  r(k), rho(k,1),   tau(k,1),   cur(k,1),    &
                  &     rho_p(k,1), tau_p(k,1), cur_p(k,1)
          end do
       end if
       close(uwave)
    end if
    write( uspe, '(///)' )
    !
    !!
    !
    write( uspe, '(//,36x,"Protons mean fields")')
    write( uspe, '(36x,"-------------------")')
    write( uspe, '(17x,"Particle-hole channel",15x,&
         &"Particle-particle channel")')
    write( uspe, '(12x,35("-"),2x,34("-"))')
    write( uspe, '(" it    r      Average     Kinetic    Spin-orbit",&
         &"   Average     Kinetic    Spin-orbit   Coulomb")')
    write( uspe, '(1x,"-- ",7("-"),1x,71("-"),2x,10("-"))')
    do k = 1, npt
       write( uspe, '(2x,i1,1x,f6.2,1x,7(1pe12.4))' ) &
            2, r(k), v1f(k,2), v2f(k,2), v3f(k,2), &
            v1p(k,2), v2p(k,2), v3p(k,2), vc(k)
    end do
    if (meanfields) then
       open( unit = uwave, file = out//'proton'//trim(extn)//'.mf', &
            status = 'unknown' )
       write( uwave, '("#                   Particle-hole channel       ",&
            &"            Particle-particle channel")')
       write( uwave, '("#         -----------------------------------------",&
            &"  ----------------------------------------")')
       write( uwave, '("#    r       Average       Kinetic      Spin-orbit",&
            &"     Average       Kinetic      Spin-orbit     Coulomb")')
       write( uwave, '("#-------------------------------------------------",&
            &"---------------------------------------------- ----------")')
       do k = 1, npt
          write( uwave, '(1x,f7.3,1x,7e14.6)' )    &
               r(k), v1f(k,2), v2f(k,2), v3f(k,2), &
               v1p(k,2), v2p(k,2), v3p(k,2), vc(k)
       end do
       close(uwave)
    end if
    write( uspe, '(///)' )
    !
    !!
    !
    if (densities) &
         open( unit = uwave, file = out//'proton'//trim(extn)//'.dens', &
         status = 'unknown' )
    write( uspe, '(36x,"Protons densities")')
    write( uspe, '(36x,"-----------------")')
    write( uspe, '(17x,"Particle-hole channel", &
         &15x,"Particle-particle channel")')
    if ( regularization ) then
       write( uspe, '(12x,34("-"),2x,46("-"))')
       write( uspe, '(" it   r      Average     Kinetic    Spin-orbit", &
            &"   Average     Kinetic    Spin-orbit    Reg.")')
       write( uspe, '(1x,93("-"))')
       do k = 1, npt
          write( uspe, '(1x,i1,1x,f6.2,1x,7(1pe12.4))' )  &
               2, r(k), rho(k,2), tau(k,2), cur(k,2),     &
               rho_p(k,2), tau_p(k,2), cur_p(k,2),        &
               rho_p(k,2) + 2 * v1p(k,2) * rega(k,2)
       end do
    else
       write( uspe, '(12x,34("-"),2x,34("-"))')
       write( uspe, '(" it   r      Average     Kinetic    Spin-orbit", &
            &"   Average     Kinetic    Spin-orbit")')
       write( uspe, '(1x,81("-"))')
       do k = 1, npt
          write( uspe, '(1x,i1,1x,f6.2,1x,6(1pe12.4))' ) &
               2, r(k), rho(k,2), tau(k,2), cur(k,2),    &
               rho_p(k,2), tau_p(k,2), cur_p(k,2)
       end do
    end if
    if (densities) then
       do k = 1, npt
          write( uwave, '(1x,f7.3,1x,6e14.6)' )    &
               r(k), rho(k,2), tau(k,2), cur(k,2), &
               rho_p(k,2), tau_p(k,2), cur_p(k,2)
       end do
       close(uwave)
    end if
    write( uspe, '(///)' )
    !
  end subroutine print_out_fields

  !
  !!
  !

  subroutine savepot(npt)
    !
    implicit none
    integer, intent(in) :: npt
    integer :: k
    !
    !! Save the potentials to use it in an other run
    !
    open( unit = uwave, file = potdir//'pot'//trim(extn)//'.mf', &
         status = 'unknown', form = 'unformatted' )
    write(uwave) h
    write(uwave) npt, del, amb
    write(uwave)  &
         ( v1f(k,1), v2f(k,1), v3f(k,1),  &
         v1p(k,1), v2p(k,1), v3p(k,1),    &
         v12f(k,1), v22f(k,1), v12p(k,1), v22p(k,1), k = 1, npt )
    write(uwave)      &
         ( v1f(k,2), v2f(k,2), v3f(k,2),      &
         v1p(k,2), v2p(k,2), v3p(k,2), vc(k), &
         v12f(k,2), v22f(k,2), v12p(k,2), v22p(k,2), k = 1, npt )
    close(uwave)

  end subroutine savepot

  !
  !!
  !

end module bogo
!
!=========================================================================
!
!
!=========================================================================
!     *    "Consider a spherical cow"                            *
!     *    (John Harte, University of California, Berkeley.)     *
!=========================================================================
!
!=========================================================================
!
!!   **** Sperical HFB code in coordinate representation ****
!! Mostly based on the code by Dobaczewski, Flocard and Treiner,
!! Nucl. Phys. A 422 (1994) 103.
!
!=========================================================================
!
!! Main program
!!
!! driver for the different subroutines
!
!=========================================================================
!

program hfb
  use cste ;  use io ;  use fields ;  use skforces ;  use bogo
  use bogostatic ;  use eqdifstatic ;  use forcesstatic
  implicit none
  integer :: npt       ! Number of points in the box..
  integer :: nfcttot   ! total number of wave functions
  integer :: it, j, ll, l, iter, nfound
  real ( kind = pr ) :: alamb ! Fermi energy for the current iteration
  ! dummy
  integer :: i, ninteg, i0, ifail
  logical :: bye, done, failed
  !
  !! Energy and energy relative evolution
  !
  real ( kind = pr ) :: previous_total_energy, delta_tot_e, &
       previous_tot_del, tot_del, delta_tot_del
  !
  hfb_input = 'hfb.input' ! you can change the name of the input file.
  !
  print '(34("*"),1x,a,i2,1x,35("*"),/)', "PR = ", pr
  call set_cste()
  call read_input(npt)
  call init_clock()
  call allocatefct( npt, dimnfobl = sum(j_max) + 2 )
  !
  do                   ! loop over nuclei to compute
     !
     call allocatebogo(npt)
     call allocateeqdif(npt)
     call allocateforces(npt)
     !
     call flush()      ! (re)set everything to 0
     mass = real( sum(npr), kind = pr )
     call subparam()   ! set system dependant parameters
     open( unit = uspe, file = out//'hfb'//trim(extn)//'.spe', &
          status = 'unknown' )
     !
     call writeparticlenumbers()
     !
     call init_sw( ma = int(mass), jz = npr(2), n = npt )
     previous_total_energy = 0.0_pr
     previous_tot_del = 0.0_pr
     delta_tot_e = 1.0_pr
     delta_tot_del = 1.0_pr
     brutal = .true.
     nfobl = 0.0_pr
     xmu = xmu0
     done = .false.
     !
     !! Beginning of the iteration loop
     !
     call firstlines()
     iter = 0
     do
        ifail = 0
        iter = iter + 1
        if ( iter < 4 ) brutal = .true.
        !
        !! Loop over the isospin, it = 1 = neutron, it = 2 = proton
        !
        failed = .false.
        ninteg = 0
        nfcttot = 0
        i0 = 1 ! index of the first function in the block in the
        !      ! previous iteration
        i = 1
        do it = 1, 2
           !
           !! Loops over angular momentum j and orbital momentum l
           !
           alamb = amb(it)
           do j = 1, j_max(it), 2
              do ll = 1, 3, 2
                 l = ( j - 2 + ll ) / 2
                 !
                 !! Solution of the HFB equations in a given block
                 !
                 ifail = 0
                 call skyhfb( it, l, j, npt, alamb, i, nfcttot, nfound, i0, &
                      done, ninteg, ifail )
                 if ( failed .or. ifail == 1 ) then
                    failed = .true.
                    !done = .false.
                 end if
                 i0 = i0 + nfobl(i)
                 nfobl(i) = nfound
                 if ( done .or. iter == it_max ) call writeoutwf() ! then
                 nfcttot = nfcttot + nfound
                 i = i + 1 ! next (it,l,j) block
              end do ! ll
           end do ! j
        end do ! it
        !
        if ( ifail > 1 ) exit ! skyhfb was not able to converge
        !
        !! Store energies
        !
        e_down(:,2) = e_down(:,1)
        e_mid(:,2) = e_mid(:,1)
        e_up(:,2) = e_up(:,1)
        !
        call quabcs( done )
        !
        !! Update densities
        !
        call sfeden( done, npt )
        !
        !! Compute energy and mean-fields
        !
        call energy(npt)
        delta_tot_e = &
             abs( ( total_energy - previous_total_energy ) / total_energy )
        if ( iter < it_max ) then
           if ( total_energy < 0.0_pr .and. .not.failed ) then
              xmu  = xmu0
              xpmu = xpmu0
              call updatefields()
           else
              xmu  = 0.95_pr
              xpmu = 0.95_pr
              call updatefields()
              if ( total_energy > 0.0_pr ) print *, ' Wrong matching ?...'
           end if
        else
           !
           call updatefields()
           !
        end if
        if ( done ) call savepot(npt)
        !
        !! The total energy is not the best quantity to determine if the
        !! equations have converged, it's interesting to take into account
        !! the paring as well (i.e. mean deltas)
        !
        tot_del = sum(abs(del))
        if ( tot_del > 1.e-4_pr ) then
           delta_tot_del = abs( ( tot_del - previous_tot_del ) / tot_del )
        else
           delta_tot_del = 0.0_pr
        end if
        if ( done .and. iter > 15 ) call firstlines()
        print '(1x,i3,2x,i6,f11.6,1x,f11.6,2f11.6,1x,f14.7,e12.4)', &
             iter, ninteg, amb, - del, total_energy, delta_tot_e
        !
        !! After one third of the iteration, check if convergence is not
        !! too slow (it can happend when pairing is very small).
        !
        if ( iter == it_max / 3 .or. iter == 2 * it_max / 3 &
             .or. iter == it_max / 2 ) then
           if ( dripline == 0 ) call doesitconverge()
        end if
        previous_tot_del = tot_del
        previous_total_energy = total_energy
        if ( iter >= it_max ) then
           ifail = 2
           exit
        end if
        if ( done ) then
           call happyend()
           exit
        end if
        if ( delta_tot_e < eps_energy  .and. &
             delta_tot_del < max_delta ) done = .true.
     end do ! iter
     if ( ifail > 1 ) call failtoconverge()
     close(uspe)
     call how_long()
     call read_input( npt, bye )
        call deallocatebogo()
        call deallocateeqdif()
        call deallocateforces()
     if (bye) then
        close(uinput)
        exit
     end if
     call allocatefct( npt, dimnfobl = sum(j_max) + 2 )
     !
  end do
  !
  call deallocatefct()
  call deallocate_r()
  !
  !=========================================================================
  !
contains

  !
  !! Several local routines written outside of the main program for
  !! sake of clarity.
  !

  subroutine doesitconverge()
    !
    implicit none
    integer :: ido
    logical :: once
    !
    once = .false.
    if ( delta_tot_e < 200 * eps_energy  &
         .and. tot_del < previous_tot_del ) then
       do ido = 1, 2
          if ( abs(del(ido)) <= 0.075_pr .and. bogolyubov(ido) ) then
             if ( .not. once ) then
                print '(3x,"warning: Iterations are not progressing well",&
                     &" and pairing seems to be small...")'
                once = .true.
             end if
             call reduce_pairing( ido, npr(ido) )
             print '(12x,"Trying to force convergence for it = ",i1)', ido
          end if
       end do
    end if
  end subroutine doesitconverge

  !
  !!
  !

  subroutine firstlines()
    implicit none
    print '(/,1x,"iter  Nint       Fermi Energies    ", &
         &"       Mean Delta       Total Energy  | D(E)/E |")'
    print   '(1x,"---- ------ ---- N --------- Z ----", &
         &" ---- N ------- Z ---- -------------- ----------")'
  end subroutine firstlines

  !
  !!
  !

  subroutine failtoconverge()
    implicit none
    !
    print '(/,8x,"**** Unable to converge for this system... ****")'
    print '(/,/,2x,82("-"),/)'
    call print_out_qp()
    call print_out_fields( npt )
    call happyend()
    !
  end subroutine failtoconverge

  !
  !!
  !

  subroutine writeparticlenumbers()
    implicit none
    character ( len = 20 ) :: tail
    !
    if ( memlog ) &
         write( ulog, '(/,10x,"N = ",i3,10x,"Z =",i3,/)' ) npr
    write( uspe, '(/,2x,"Nucleus:",/,2x,"--------")' )
    write( uspe, '(5x,"N = ",i3,10x,"Z =",i3,/)' ) npr
    write( uspe, '(/,2x,"Force: ",a,/,2x,"------")' ) force
    if ( bogolyubov(1) .or. bogolyubov(2) ) then
       tail = ""
       if ( regularization ) tail = " with regularization"
       select case (pairing_force)
       case (1)
          write( uspe, '(/,2x,"Pairing force: ",a,/,2x,"--------------",/)' ) &
               "Volume pairing"//tail
       case (2)
          write( uspe, '(/,2x,"Pairing force: ",a,/,2x,"--------------",/)' ) &
               "Surface pairing"//tail
       case (3)
          write( uspe, '(/,2x,"Pairing force: ",a,/,2x,14("-"),/)' ) &
               "Volume + Surface pairing"//tail
       end select
    end if
    !
  end subroutine writeparticlenumbers

  !
  !!
  !

  subroutine writeoutwf()
    implicit none
    if ( quasiparticles ) &
         call print_out_wf( nfcttot + 1, nfound, npt, it, l, j )
    if ( canonical_states ) &
         call print_out_can( nfcttot + 1, nfound, npt, it, l, j )
  end subroutine writeoutwf

  !
  !!
  !

  subroutine happyend()
    !
    implicit none
    real ( kind = pr ) :: part_number(2), radius(2)
    !
    print '(/,21x,"Neutrons        Protons          Total")'
    print '(1x,62("."))'
    print '(1x,"Fermi Ener.  =",2f16.8)', amb
    print '(1x,"Mean Gaps    =",2f16.8)', del
    print '(1x,"Kinetic En.  =",3f16.8)', kinetic_energy
    print '(1x,"Pairing En.  =",3f16.8)', pairing_energy
    print '(1x,"Pair.Kin.En. =",3f16.8)', pairing_kinetic_energy
    print '(/,1X,"Energies (in MeV):",/,2x,"E/A =",f11.6,4x,&
         &"Variat. =",f11.6,6x,"=====>  Etot =",f13.6)', &
         energy_per_nucleon, abs( delta_tot_e * total_energy ), &
         total_energy
    print '(1x,"Contributions:")'
    print '(2x,"Field =",f13.6,2x,"Spin-Or. =",f12.6,   &
         &  5x,"Coul. =",f12.6)', &
         field_energy, spin_orbit_energy , coulomb_energy
    print '(1x,"(Rear. =",f13.6,")",25x,"Coul.Ex. =",f12.6)', &
         rearangment_energy, ex_coulomb_energy
    write( usumup, '(2i4,2i3,f6.1,2f9.4,2f7.4,2f11.5,2f10.5,&
         &f11.5,f10.5,f10.5,f15.7)', advance = "no" ) &
         npr, j_max, cut_off, amb, abs(del),&
         kinetic_energy(1:2), pairing_energy(1:2), &
         spin_orbit_energy , coulomb_energy, ex_coulomb_energy, &
         total_energy
    part_number(1) = sum( r2 * rho(:,1) ) * qp * h
    part_number(2) = sum( r2 * rho(:,2) ) * qp * h
    radius(1) = sqrt( sum( r2 * r2 * rho(:,1) ) * h * qp / part_number(1) )
    radius(2) = sqrt( sum( r2 * r2 * rho(:,2) ) * h * qp / part_number(2) )
    !
    print '(/,20x,"Neutrons       Protons          Total         Charge")'
    print '(" Part. numbers: ",f12.6,3x,f12.6)', part_number
    if ( cano ) &
         print '(" Fluctuations:  ",f12.6,3x,f12.6)', dispersion
    print '(" Radii:         ",f12.6,3(3x,f12.6))', radius, &
         sqrt( sum( radius**2 * part_number ) / sum( part_number ) ), &
         sqrt( radius(2)**2 + 0.64_pr )
    write( usumup, '(4f8.5)' ) radius, &
         sqrt( sum( radius**2 * part_number ) / sum( part_number ) ), &
         sqrt( radius(2)**2 + 0.64_pr )
    print '(/,2x,82("-"),/)'

  end subroutine happyend

  !
  !!
  !

end program hfb

!
!=========================================================================
!
