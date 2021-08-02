PROGRAM generate_vsrg_vlowk

!***************************************************************************
! Stand-alone main program for the new vsrg/vlowk codes to generate smooth 
!  (and sharp) potentials with various regulators.
!
! This code can be used to generate files with vlowk potentials in specified
!  partial waves, or it can be a prototype for directly calling the vlowk
!  subroutines.
!
! Main Programmer for vlowk routines:
!   Scott Bogner  bogner@mps.ohio-state.edu
!
! Sous Chef (this main program)
!   Dick Furnstahl  furnstahl.1@osu.edu
!
! Revision history:
!   06-Aug-2006 --- original version, based on previous codes
!   04-Oct-2006 --- upgraded to include Similarity RG option.
!   18-Nov-2006 --- upgraded to output bare potential
!   29-Nov-2006 --- output vlowk with zeros beyond nmod
!   02-Feb-2007 --- new vlowk/vsrg codes from Scott
!   19-Dec-2014 --- added Weinberg eigenvalue routines
!
!***************************************************************************

!***************************************************************************
! Include relevant modules
!
USE nrtype                      ! constants and KIND definitions 
USE vnn_module, ONLY: vnnmompw   ! calculate starting NN potential
  ! 
USE vlowkmodule, ONLY: VlowkType,ConstructVlowkType,GetVlowk,GetVlowkSize, &
                               GetVlowkMethod    
  ! 
USE vlowkperipherals, ONLY: MeshType,ConstructMeshType,GetMeshWts, &
                             GetNmod,GetNtot,GetKmax,GetLambdaeff,GetLambda, &
                             ChannelType,ConstructChannelType,GetCoupled, &
			     GetL,GetS,GetJ,GetT,GetHb2m, &
                             MethodType,ConstructMethodType
  ! input parameters for vlowk type, mesh, regulator
USE vlowk_input, ONLY: Input_parameters, kvnn, Lambda, fac, smooth, ireg, &
                        nsmooth, rsmooth, imethod, iherm, &
                        ntot, nmod, kmax, &
			Input_channel, l, s, jt, itz, it  

USE convergence, ONLY: weinberg       
!
!***************************************************************************

!***************************************************************************
! Type definitions
IMPLICIT NONE

! All components of the following derived types are PRIVATE. 
!  You have to use the PUBLIC "constructor" and "accessor" functions/subroutines 
!  declared in the relevant MODULES to initialize them and to gain access 
!  to the individual components.
TYPE(MeshType) :: kmesh        ! derived type for mesh/wts and relevant 
                               !   parameters (see vlowkperipherals)
TYPE(ChannelType) :: channel   ! derived type for nn-channel quantum #'s, etc. 
                               !   (see vlowkperipherals)
TYPE(MethodType) :: methods    ! derived type specifying the regulators, 
                               !   hermitization, etc (see vlowkperipherals)
TYPE(VlowkType) :: vlowk       ! derived type for vlowk (herm, non-herm, 
                               !   matrix size, etc)  (see vlowkmodule)
 
! Local variables 
!!INTEGER(I4B) :: l,s,jt,itz,it  ! channel variables
INTEGER(I4B) :: i,j            ! dummy index variables
INTEGER(I4B) :: nout           ! actual number of mesh points
REAL(DP) :: convert, convert_vnn            ! factor to convert vlowk units
REAL(DP),ALLOCATABLE :: xk(:), wk(:),vnn(:,:),hnn(:,:) ,vlk(:,:)
REAL(DP) :: small, zero

REAL(DP) :: Ecm, Pcm, kfermi
COMPLEX(16),ALLOCATABLE :: eta(:), psi(:,:)
REAL(DP),ALLOCATABLE :: vnn_scaled(:,:), vlk_scaled(:,:)
	      
!
!***************************************************************************
!
! Get the vlowk parameters from an input file (given 1st on the command line)
CALL Input_parameters

! Get the channel parameters from an input file (given 2nd on command line)
CALL Input_channel

IF (imethod.EQ.3) THEN
   nout = ntot
 ELSE
   nout = nmod
ENDIF 

! put all equal to nout  
nout = ntot
	      
!
!***************************************************************************
! Relevant NN channel data for the ChannelType structure "channel"
!
!   itz --- specify np or pp or nn [integer]
!     l --- orbital angular momentum [integer]
!     s --- spin [integer]
!    jt --- total j [integer]
!
! 1S0 np example 
!   itz = 0
!   l = 0         
!   s = 0 
!   jt = 0
!
! 3S1-3D1 np example 
!  itz = 0
!  l = 0         
!  s = 1 
!  jt = 1
!
! 3P2-3F2 np example
!  itz = 0
!  l = 1
!  s = 1
!  jt = 2
!
! 3D2 np example
!  itz = 0
!  l = 2
! s = 1
!  jt = 2
!
!***************************************************************************

! Here is the calling sequence to build a potential in a given channel 

! 1. Allocate space for arrays (allow for coupled channel)
ALLOCATE(xk(ntot),wk(ntot),vnn(2*ntot,2*ntot),hnn(2*ntot,2*ntot))

! 2. Construct the MeshType structure kmesh
CALL ConstructMeshType(kmesh,Lambda,fac*Lambda,kmax,nmod,ntot)   

! 3. Retrieve the xk(:) and wk(:) PRIVATE components of kmesh
CALL GetMeshWts(kmesh,xk,wk)  

! 4. Initialize the MethodsType structure imethod
methods=ConstructMethodType(imethod,iherm,ireg,rsmooth,nsmooth,lambda)  

! 5. Initialize the ChannelType structure 
!  (you can also use an overloaded version w/ nwave,itz as arguments)
channel = ConstructChannelType(l,s,jt,itz)    
it = GetT(channel)   ! get the corresponding isospin

! 6. Construct the bare vnn potential and Hamiltonian (hnn)
CALL vnnmompw(GetCoupled(channel),kvnn,l,s,jt,it,1,-1,xk,wk,ntot,hnn,vnn) 

! 7. Construct the VlowkType structure vlowk
CALL ConstructVlowkType(kmesh,channel,methods,vnn,vlowk) 

! 8. Retrieve the PRIVATE matrix elements of the vlowk structure 
!  (returned in vlk(:,:), which is first allocated and initialized to 0)
ALLOCATE(vlk(1:GetVlowkSize(vlowk),1:GetVlowkSize(vlowk)));  vlk=0_dp
CALL GetVlowk(vlowk,vlk)  

! Here we choose to simply print the results to files for later use.

! Print out the mesh to an output file
WRITE(6,*) "Writing mesh to generate_vsrg_mesh.dat . . ."
OPEN (UNIT=17, FILE='generate_vsrg_vlowk_mesh.dat', STATUS='REPLACE')
   
!WRITE(17,'(''# kvnn='',i2,''  Lambda='', f5.2, ''  factor='', &
!       f4.2, '' regulator='', i2, '' integer smooth='', i2)') &
!       kvnn, lambda, fac, ireg, nsmooth
!WRITE(17,'(''#  real smooth='', f8.4, '' method='', i2, &
!       '' hermitization='', i2)') rsmooth, imethod, iherm
DO i = 1,nout
  WRITE(17,'(1pe22.13,3x,1pe22.13)') xk(i), wk(i)
ENDDO       

CLOSE(17)


! Print out the potentials to output files
!  generate_vsrg_vlowk_bare.out --- initial V_nn potential
!  generate_vsrg_vlowk.out --- V_lowk or V_srg
!  generate_vsrg_vlowk_diag.out --- diagonal V_lowk or V_srg elements
!
WRITE(6,*) "Writing output to generate_vsrg.out . . ."
OPEN (UNIT=12, FILE='generate_vsrg_vlowk.out', STATUS='REPLACE')
OPEN (UNIT=13, FILE='generate_vsrg_vlowk_diag.out', STATUS='REPLACE')
OPEN (UNIT=77, FILE='generate_vsrg_vlowk_bare.out', STATUS='REPLACE')
OPEN (UNIT=78, FILE='generate_vsrg_vlowk_bare_diag.out', STATUS='REPLACE')
! EDIT BY AT: Need file with Gaussian weights: open and save a file
!OPEN (UNIT=79, FILE='gw.dat', STATUS='REPLACE')
   
WRITE(12,'(''# vlowk itz = '', i2, '' l = '', i2, '' s = '', i2, &
               '' jt = '', i2, '' kvnn =  '', i2, ''  Lambda = '', &
	       f6.2 )') itz, l, s, jt, kvnn, lambda 
WRITE(12,'(''#  fac='', &
       f4.2, '' reg.='', i2, '' int. smooth='', i2,  &
       '' real smooth='', f8.4, '' meth='', i2, '' herm='', i2)') &
       fac, ireg, nsmooth, rsmooth, imethod, iherm

WRITE(13,'(''# vlowk itz = '', i2, '' l = '', i2, '' s = '', i2, &
               '' jt = '', i2, '' kvnn =  '', i2, ''  Lambda = '', &
	       f6.2 )') itz, l, s, jt, kvnn, lambda 
WRITE(13,'(''#  fac='', &
       f4.2, '' reg.='', i2, '' int. smooth='', i2,  &
       '' real smooth='', f8.4, '' meth='', i2, '' herm='', i2)') &
       fac, ireg, nsmooth, rsmooth, imethod, iherm
   
WRITE(77,'(''# vlowk itz = '', i2, '' l = '', i2, '' s = '', i2, &
               '' jt = '', i2, '' kvnn =  '', i2, ''  Lambda = '', &
	       f6.2 )') itz, l, s, jt, kvnn, lambda 
WRITE(77,'(''#  fac='', &
       f4.2, '' reg.='', i2, '' int. smooth='', i2,  &
       '' real smooth='', f8.4, '' meth='', i2, '' herm='', i2)') &
       fac, ireg, nsmooth, rsmooth, imethod, iherm
   
WRITE(78,'(''# vlowk itz = '', i2, '' l = '', i2, '' s = '', i2, &
               '' jt = '', i2, '' kvnn =  '', i2, ''  Lambda = '', &
         f6.2 )') itz, l, s, jt, kvnn, lambda 
WRITE(78,'(''#  fac='', &
       f4.2, '' reg.='', i2, '' int. smooth='', i2,  &
       '' real smooth='', f8.4, '' meth='', i2, '' herm='', i2)') &
       fac, ireg, nsmooth, rsmooth, imethod, iherm

convert = (pi/2._dp)/GetHb2m(channel)   ! \pi/(2*hbar^2/M)
IF (GetCoupled(channel)) THEN   ! coupled channel
   WRITE(12,'(''#     k            kp           V11             V12      '', &
              ''        V21             V22 '')')
   WRITE(77,'(''#     k            kp           V11             V12      '', &
              ''        V21             V22 '')')
   DO i = 1,nout
      ! EDIT BY AT: Write Gaussian weights to file here
      !WRITE(79,'(f11.6)') wk(i)
      DO j = 1,nout
         zero = 0.0_dp;
         IF (imethod.EQ.3) THEN
           WRITE(12,'(2f11.6,1x,4(1pe21.12))') xk(i), xk(j), & 
	         small(convert*vlk(i,j)), &
                 small(convert*vlk(i,j+nout)), small(convert*vlk(i+nout,j)), &
	         small(convert*vlk(i+nout,j+nout))
         ELSE IF ((i.LE.nmod).AND.(j.LE.nmod)) THEN
           WRITE(12,'(2f11.6,1x,4(1pe21.12))') xk(i), xk(j), & 
	         small(convert*vlk(i,j)), &
                 small(convert*vlk(i,j+nmod)), small(convert*vlk(i+nmod,j)), &
	         small(convert*vlk(i+nmod,j+nmod))
         ELSE
           WRITE(12,'(2f11.6,1x,4(1pe21.12))') xk(i), xk(j), & 
	         zero, zero, zero, zero
         ENDIF
          
         convert_vnn = convert/xk(i)/xk(j)/sqrt(wk(i)*wk(j))
         WRITE(77,'(2f11.6,1x,4(1pe21.12))') xk(i), xk(j), & 
	       small(convert_vnn*vnn(i,j)), &
               small(convert_vnn*vnn(i,j+nout)), &
	       small(convert_vnn*vnn(i+nout,j)), &
	       small(convert_vnn*vnn(i+nout,j+nout))
      ENDDO
   ENDDO
 ELSE                           ! uncoupled channel
   WRITE(12,'(''#     k            kp           V '')')
   WRITE(77,'(''#     k            kp           V '')')
   DO i = 1,nout
      DO j = 1,nout
         convert_vnn = convert/xk(i)/xk(j)/sqrt(wk(i)*wk(j))
         WRITE(77,'(2f11.6,1x,1(1pe21.12))') xk(i), xk(j), &
	    small(convert_vnn*vnn(i,j))
	    
         zero = 0.0_dp;
         IF (imethod.EQ.3) THEN
           WRITE(12,'(2f11.6,1x,1(1pe21.12))') xk(i), xk(j), &
	      small(convert*vlk(i,j))
         ELSE IF ((i.LE.nmod).AND.(j.LE.nmod)) THEN
           WRITE(12,'(2f11.6,1x,1(1pe21.12))') xk(i), xk(j), &
	      small(convert*vlk(i,j))
         ELSE
           WRITE(12,'(2f11.6,1x,1(1pe21.12))') xk(i), xk(j), & 
	         zero
         ENDIF
      ENDDO
   ENDDO
ENDIF

IF (GetCoupled(channel)) THEN   ! coupled channel
   WRITE(13,'(''#     k               V11             V12      '', &
              ''        V21             V22 '')')
   WRITE(78,'(''#     k               V11             V12      '', &
              ''        V21             V22 '')')
   DO i = 1,nout
      WRITE(13,'(1f11.6,1x,4(1pe15.6))') xk(i), convert*vlk(i,i), &
                     convert*vlk(i,i+nout), convert*vlk(i+nout,i), &
		     convert*vlk(i+nout,i+nout)
      convert_vnn = convert/xk(i)/xk(i)/sqrt(wk(i)*wk(i))
      WRITE(78,'(1f11.6,1x,4(1pe15.6))') xk(i), convert_vnn*vnn(i,i), &
                convert_vnn*vnn(i,i+nout), convert_vnn*vnn(i+nout,i), &
                convert_vnn*vnn(i+nout,i+nout)
   ENDDO
 ELSE                           ! uncoupled channel
   WRITE(13,'(''#     k               V '')')
   WRITE(78,'(''#     k               V '')')
   DO i = 1,nout
      WRITE(13,'(1f11.6,1x,1(f15.8))') xk(i), small(convert*vlk(i,i))
      convert_vnn = convert/xk(i)/xk(i)/sqrt(wk(i)*wk(i))
      WRITE(78,'(1f11.6,1x,1(f15.8))') xk(i), small(convert_vnn*vnn(i,i))
   ENDDO
ENDIF

CLOSE(12)
CLOSE(13)
CLOSE(14)
CLOSE(77)
CLOSE(78)
! EDIT BY AT: Close gw.dat
!CLOSE(79)

! Now we try out the Weinberg eigenvalue calculatation
!  The calling sequence is
!     CALL weinberg(coupled, P, kf, lambda, Ecm, xk, wk, vnn, np, eta, psi)
!  where
!    coupled = true/false
!    P = COM mom. (I've only really ever looked at P=0)
!    kf = fermi mom. (set to 0.0 for free space, 
!                      otherwise it uses pauli blocked propagators)
!    lambda = UV cutoff (either the vlowk one, or whatever the maximum mesh 
!                        is for smooth cutoffs)
!    xk = momentum mesh
!    wk = momentum wts
!    vnn = potential m.e.'s in fm  (np x np for uncoupled, 
!            2np x 2np for coupled channel)
!    Ecm = COM energy in MeV
!    eta = complex vector of length np (uncoupled) 2np (coupled) 
!            with the eigenvalues
!    psi = complex array with the weinberg states
    
ALLOCATE(vlk_scaled(1:GetVlowkSize(vlowk),1:GetVlowkSize(vlowk)))

vlk_scaled = convert * vlk

IF (GetCoupled(channel)) THEN   ! coupled channel

   ALLOCATE(psi(2*ntot,1:2*ntot),eta(2*ntot))
   ALLOCATE(vnn_scaled(2*ntot,2*ntot))

   DO i = 1,nout
      DO j = 1,nout
         convert_vnn = convert/xk(i)/xk(j)/sqrt(wk(i)*wk(j))
         vnn_scaled(i,j) = convert_vnn * vnn(i,j)
         vnn_scaled(i,j+nout) = convert_vnn * vnn(i,j+nout)
         vnn_scaled(i+nout,j) = convert_vnn * vnn(i+nout,j)
         vnn_scaled(i+nout,j+nout) = convert_vnn * vnn(i+nout,j+nout)
      ENDDO
   ENDDO
 ELSE                           ! uncoupled channel

   ALLOCATE(psi(ntot,ntot),eta(ntot))
   ALLOCATE(vnn_scaled(ntot,ntot))

   DO i = 1,nout
      DO j = 1,nout
         convert_vnn = convert/xk(i)/xk(j)/sqrt(wk(i)*wk(j))
         vnn_scaled(i,j) = convert_vnn * vnn(i,j)
      ENDDO
   ENDDO
ENDIF


Pcm = 0.0d0
kfermi = 0.0d0
kmax = GetKmax(kmesh)

OPEN (UNIT=37, FILE='generate_vsrg_vlowk_weinberg.out', STATUS='REPLACE')
OPEN (UNIT=39, FILE='generate_vsrg_vlowk_weinberg2.out', STATUS='REPLACE')

Ecm = 0.0d0
CALL weinberg(GetCoupled(channel),Pcm,kfermi,kmax,Ecm,xk,wk,vnn_scaled,nout,eta,psi)
WRITE(37,'(f7.2,2x,5(f8.4,f8.4,1x,f8.4,2x))') Ecm, ( eta(i), abs(eta(i)), i=1,5 )

Ecm = 25.0d0
CALL weinberg(GetCoupled(channel),Pcm,kfermi,kmax,Ecm,xk,wk,vnn_scaled,nout,eta,psi)
WRITE(37,'(f7.2,2x,5(f8.4,f8.4,1x,f8.4,2x))') Ecm, ( eta(i), abs(eta(i)), i=1,5 )

Ecm = 66.0d0
CALL weinberg(GetCoupled(channel),Pcm,kfermi,kmax,Ecm,xk,wk,vnn_scaled,nout,eta,psi)
WRITE(37,'(f7.2,2x,5(f8.4,f8.4,1x,f8.4,2x))') Ecm, ( eta(i), abs(eta(i)), i=1,5 )

Ecm = 100.0d0
CALL weinberg(GetCoupled(channel),Pcm,kfermi,kmax,Ecm,xk,wk,vnn_scaled,nout,eta,psi)
WRITE(37,'(f7.2,2x,5(f8.4,f8.4,1x,f8.4,2x))') Ecm, ( eta(i), abs(eta(i)), i=1,5 )

Ecm = 150.0d0
CALL weinberg(GetCoupled(channel),Pcm,kfermi,kmax,Ecm,xk,wk,vnn_scaled,nout,eta,psi)
WRITE(37,'(f7.2,2x,5(f8.4,f8.4,1x,f8.4,2x))') Ecm, ( eta(i), abs(eta(i)), i=1,5 )

Ecm = 200.0d0
CALL weinberg(GetCoupled(channel),Pcm,kfermi,kmax,Ecm,xk,wk,vnn_scaled,nout,eta,psi)
WRITE(37,'(f7.2,2x,5(f8.4,f8.4,1x,f8.4,2x))') Ecm, ( eta(i), abs(eta(i)), i=1,5 )

Ecm = 250.0d0
CALL weinberg(GetCoupled(channel),Pcm,kfermi,kmax,Ecm,xk,wk,vnn_scaled,nout,eta,psi)
WRITE(37,'(f7.2,2x,5(f8.4,f8.4,1x,f8.4,2x))') Ecm, ( eta(i), abs(eta(i)), i=1,5 )

Ecm = 300.0d0
CALL weinberg(GetCoupled(channel),Pcm,kfermi,kmax,Ecm,xk,wk,vnn_scaled,nout,eta,psi)
WRITE(37,'(f7.2,2x,5(f8.4,f8.4,1x,f8.4,2x))') Ecm, ( eta(i), abs(eta(i)), i=1,5 )

DO j = 1,76
  Ecm = DBLE(j-1)*2.d0
  CALL weinberg(GetCoupled(channel),Pcm,kfermi,kmax,Ecm,xk,wk,vnn_scaled,nout,eta,psi)
  WRITE(39,'(f7.2,2x,5(f8.4,f8.4,1x,f8.4,2x))') Ecm, ( eta(i), abs(eta(i)), i=1,5 )
ENDDO  

CLOSE(37)
CLOSE(39)

!Ecm = -2.224
!CALL weinberg(GetCoupled(channel),Pcm,kfermi,kmax,Ecm,xk,wk,vnn_scaled,nout,eta,psi)
!DO i = 1,5  ! nout
!   WRITE(38,'(2(f15.8),3x)') eta(i)
!ENDDO




! Recover the array memory
DEALLOCATE(xk,wk,vnn,hnn,vlk)  

! That's all folks!

STOP
END

!***************************************************************************
!***************************************************************************

   REAL*8 FUNCTION small(x)
! 
! Simple function to eliminate output problem with Intel ifort compiler
!  when printing formatted output with number less than 1.e-99 in magnitude.
!  (The "E" gets dropped when there is a third digit in the exponent.)
!
   USE nrtype                      ! constants and KIND definitions 
   REAL*8 x
   IF (abs(x) .LT. 1.e-99_dp) THEN 
     small = 0.
   ELSE IF (abs(x) .GT. 1.e20_dp) THEN
     small = 0.  
   ELSE
     small = x
   ENDIF
   
   RETURN
   END

!***************************************************************************
!***************************************************************************

