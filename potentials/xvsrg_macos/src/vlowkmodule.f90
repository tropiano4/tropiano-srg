! CAUTION: Documentation of this code is incomplete. 

! Coded by Scott Bogner 
! This code is the front-end wrapper module to calculate 
! low-momentum interactions ("V_lowk") and the related 
! similarity renormalization group ("V_SRG") interactions. 
!
! 
! ############################### CHANGE LOG ##############################################            
!
!  1/25/07         Fixed a bug in the wrapper routines. The wrong value of isospin
!                  was being passed to the bare vnn subroutine (it was always setting T = 0).

!  1/10/07         Added some new items to the VlowkType structure (which are accessible by
!                  the PUBLIC "accessor" functions listed below). For example, if you create
!                  an object of TYPE(VlowkType), you can retrieve the vlowk matrix elements
!                  with the gaussian weights and measure factors either attached or 
!                  unattached.  Also, I've defined a new derived type structure called VlowkInputsType 
!                  to streamline passing of the many input parameters (see below) the user must specify.
!
!  1/4/07          Split out all the computational routines associated with sharp and smooth cutoff 
!                  V_lowk and V_SRG calculations into separate modules (vlowk_sharp_cutoff.f90, 
!                  vlowk_smooth_cutoff.f90, v_srg.f90, hermitization.f90). Now vlowkmodule.f90 is 
!                  just a front-end wrapper.
!   
!  1/4/07         Added two wrapper subroutines (VlowkWrapperAndreas and VlowkWrapper) 
!                 to access the codes, in addition to the existing calling sequence. All of 
!                 these calling sequences will be illustrated below. 
! ############################### END CHANGE LOG ##############################################

! ############################ CALLING SEQUENCE INSTRUCTIONS ##################################
!
!                                     WAY # 1 
!
!    NOTE-- This assumes that the user is either using our bare vnn code (allvnnmodels.f), or at 
!           least has his or her vnn codes set up to have a wrapper interface like our vnnmompw 
!           subroutine, i.e.
!   
!        INTERFACE         
!           SUBROUTINE vnnmompw(icoup,kvnn,ll,is,jt,it,itz1,itz2,xkp,xwkp,nkp,hnn,vnn)
!              implicit real*8 (a-h,o-z)
!              implicit integer*4 (i-n)
!              LOGICAL::icoup
!              dimension xkp(:),xwkp(:)
!              dimension vnn(:,:),hnn(:,:)
!           END SUBROUTINE
!        END INTERFACE
!   
!   Given a set of k-space mesh points and weights and a specified partial wave
!
!            icoup      - logical (.true. for coupled channel, .false. for uncoupled channel)
!
!  a) In the calling subroutine place the following USE
!     statements and declare the derived-type structures
!     VlowkInputs and nnChannel:
!
!            USE vlowkmodule, ONLY:VlowkWrapperAndreas, VlowkInputsType, ConstructVlowkInputsType
!            USE vlowkperipherals, ONLY:ChannelType, ConstructChannelType
!            USE vnn_module, ONLY:vnnmompw                 
!            TYPE(VlowkInputsType) :: VlowkInputs
!            TYPE(ChannelType) :: nnChannel
!
!  b) Initialize the 2 data structures that are passed to the vlowk codes:
!
!            VlowkInputs = ConstructVlowkInputsType(imethod, iherm, ireg, nsmooth, rsmooth, ntot, &
!                              nmod, kmax, lambda, lambdaeff,kvnn)      
!
!            nnChannel = ConstructChannelType(l12, s12, j12, Tz)

!      where the above input parameters are (list and explain the various inputs below...)
!      
!  c) And finally... 
!
!           CALL VlowkWrapperAndreas( VlowkInputs, nnChannel, vnnmompw, xk, wk, vlowk )   
!           
!  where     
!             xk(1:ntot) and wk(1:ntot) are the k-space mesh points/weights that vlowk is calculated on,     
!             and the vlowk matrix elements (in units of MeV-fm**3 ) are returned in vlowk(:,:).
! 
!             For the vlowk calculations (sharp or smooth), vlowk is defined on the lowest NMOD
!             mesh points. I.e., for an uncoupled channel (e.g., 1S0), the code returns calculated
!             matrix elements in vlowk(1:nmod,1:nmod), while for a coupled channel (e.g., 3S1-3D1) 
!             the vlowk elements are returned in vlowk(1:2*nmod,1:2*nmod)  [i.e., vlowk contains the 
!             sub-blocks of 3S1-3S1, 3S1-3D1, 3D1-3S1, 3D1-3D1, where each block is defined on
!             NMOD points]. You can encompass both cases by dimensioning as vlowk(1:2*nmod, 1:2*nmod).
!
!             For the SRG calculations, the V_srg is defined on the full momentum space grid (NTOT points). 
!             Therefore, the uncoupled V_srg elements are returned in vlowk(1:ntot,1:ntot), while a coupled
!             channel would return elements in vlowk(1:2*ntot, 1:2*ntot).     
!************************************************************************************************************!

MODULE VlowkModule
 USE nrtype
 USE VlowkPeripherals
 USE vlowk_sharp_cutoff, ONLY:leesuzukivlowk
 USE vlowk_smooth_cutoff, ONLY:smoothvlowk
 USE v_srg, ONLY:CalcVsrg
 USE hermitization, ONLY:hermitize
 IMPLICIT NONE
 PRIVATE

     INTEGER(I4B), PARAMETER :: ndim = 400

     
   TYPE VlowkType
      PRIVATE
      INTEGER(I4B) :: theVlowkSize
      TYPE(MeshType) :: theVlowkKmeshType
      TYPE(ChannelType) :: theVlowkChannel
      TYPE(MethodType) :: theVlowkMethod
      REAL(DP) :: theNonHermVlowk(1:ndim,1:ndim) 
      REAL(DP) :: theVlowk(1:ndim,1:ndim) 
      REAL(DP) :: theVlowk_w_kkww(1:ndim,1:ndim)
      REAL(DP) :: theVlowkMesh(1:ndim)
      REAL(DP) :: theVlowkWts(1:ndim)
   END TYPE VlowkType

   TYPE VlowkInputsType
      PRIVATE
      INTEGER(I4B) :: theKvnn
      INTEGER(I4B) :: theMethod
      INTEGER(I4B) :: theHerm
      INTEGER(I4B) :: theReg
      INTEGER(I4B) :: theIntSmooth
      REAL(DP)     :: theRealSmooth
      INTEGER(I4B) :: theNtot
      INTEGER(I4B) :: theNmod
      REAL(DP)     :: theKmax
      REAL(DP)     :: theLambda
      REAL(DP)     :: theLambdaeff
   END Type VlowkInputsType
   

   PUBLIC :: VlowkType,GetVlowkSize,ConstructVlowkType,GetVlowk,GetNonHermVlowk, &
             GetVlowkKmeshType, GetVlowkChannel, GetVlowkMethod  &
             ,VlowkInputsType, ConstructVlowkInputsType,  &
              GetVlowk_w_kkww
           



   PUBLIC :: VlowkWrapperAndreas, VlowkWrapper

CONTAINS

  
  SUBROUTINE VlowkWrapperAndreas( Inputs, nnChannel, vnnmompw, xk, wk, vlowk)
     TYPE(VlowkInputsType) :: Inputs
     TYPE(ChannelType) :: nnChannel
     TYPE(MethodType)  :: vlkmethod 
     TYPE(MeshType)    :: kmesh
     TYPE(VlowkType)   :: vlowkstruct
     LOGICAL(LGT) :: coupled
     INTEGER(I4B) :: kvnn, l, s, jt, it, tz1, tz2, ntot, tz 
     REAL(DP) :: xk(:), wk(:), vlowk(:,:)
     REAL(DP), ALLOCATABLE :: vnn(:,:), hnn(:,:) 
     INTENT(IN) :: Inputs, nnChannel
     INTENT(OUT) :: xk, wk, vlowk

!  MAKE SURE TO MODIFY THE INTERFACE BELOW TO CONFORM WITH YOUR VERSION OF vnnmompw
INTERFACE
     SUBROUTINE vnnmompw(icoup,kvnn,ll,is,jt,it,itz1,itz2,xkp,xwkp,nkp,hnn,vnn)
!     SUBROUTINE vnnmompw(icoup,ll,is,jt,it,itz1,itz2,xkp,xwkp,nkp,hnn,vnn) ! uncomment this for andreas's vnnmompw routine
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      LOGICAL::icoup
      dimension xkp(:),xwkp(:)
      dimension vnn(:,:),hnn(:,:)
     END SUBROUTINE
END INTERFACE


! build MeshType structure kmesh 
     CALL ConstructMeshType( kmesh, Inputs%theLambda, Inputs%theLambdaeff, Inputs%theKmax, &
                            Inputs%theNmod, Inputs%theNtot )  
     CALL GetMeshWts(kmesh, xk, wk)
! build MethodType structure vlkmethod (holds various inputs for the type of Vlowk/SRG calc. to be done)
     vlkmethod = ConstructMethodType( Inputs%theMethod, Inputs%theHerm, Inputs%theReg, Inputs%theRealSmooth, &
                                     Inputs%theIntSmooth, Inputs%theLambda)


! allocate space for the "bare" vnn model matrix and initialize channel needed to call vnnmompw(....)
     ALLOCATE(vnn(2*Inputs%theNtot, 2*Inputs%theNtot)); vnn = 0_dp     
     ALLOCATE(hnn(2*Inputs%theNtot, 2*Inputs%theNtot)); hnn = 0_dp     
          
     l  = GetL(nnChannel) ; s  = GetS(nnChannel) ; it  = GetT(nnChannel) 
     tz = GetTz(nnChannel); jt = GetJ(nnChannel)
     
     IF(tz.eq.1)THEN
        tz1 = 1; tz2 = 1
     ELSE IF(tz.eq.0)THEN
        tz1 = 1; tz2 = -1
     ELSE IF(tz.eq.-1)THEN
        tz1 = -1; tz2 = -1
     ELSE
        WRITE(*,*)'invalid tz value must have tz = 0, 1 , or -1'
        CALL ABORT
     ENDIF

     kvnn = Inputs%theKvnn     !which vnn model to use
     coupled = GetCoupled(nnchannel)
     ntot = Inputs%theNtot
 
     CALL vnnmompw( coupled, kvnn, l, s, jt, it, tz1, tz2, xk, wk, ntot, hnn, vnn) !my bare vnn subroutine
!    CALL vnnmompw( coupled, l, s, jt, it, tz1, tz2, xk, wk, ntot, hnn, vnn) ! uncomment this for andreas's vnnmompw routine



     CALL ConstructVlowkType(kmesh, nnChannel, vlkmethod, vnn, vlowkstruct) 
     vlowk = 0_dp 
     CALL GetVlowk(vlowkstruct, vlowk(1:vlowkstruct%theVlowkSize, 1:vlowkstruct%theVlowkSize))

  
     DEALLOCATE(vnn,hnn)


  END SUBROUTINE VlowkWrapperAndreas
 
  SUBROUTINE VlowkWrapper( Inputs, nnChannel, Getvnn, vlowk_kmesh, vlowk_kwts, vlowk ) 
     TYPE(VlowkInputsType) :: Inputs
     TYPE(ChannelType) :: nnChannel
     TYPE(MethodType)  :: vlkmethod 
     TYPE(MeshType)    :: kmesh
     TYPE(VlowkType)   :: vlowkstruct
     REAL(DP) :: vlowk_kmesh(:), vlowk_kwts(:), vlowk(:,:), trace_eff, trace_bare
     REAL(DP), ALLOCATABLE :: vnn(:,:) , xxk(:), wwk(:)    
     INTEGER(I4B) :: i
     INTENT(IN) :: Inputs, nnChannel
     INTENT(OUT) :: vlowk_kwts, vlowk_kmesh, vlowk

INTERFACE 
  SUBROUTINE getvnn( channel, ntot, xk, wk, vnn )
   USE vnn_module
   USE vlowkperipherals
   IMPLICIT NONE
   TYPE(ChannelType) :: channel
   INTEGER :: iT, l, s, jt, itz, ntot
   DOUBLE PRECISION :: xk(:), wk(:), vnn(:,:), temp(400,400)
   LOGICAL :: coupled

  END SUBROUTINE getvnn
END INTERFACE

     ALLOCATE(xxk(1:Inputs%theNtot), wwk(1:Inputs%theNtot))
     CALL ConstructMeshType( kmesh, Inputs%theLambda, Inputs%theLambdaeff, Inputs%theKmax, &
                            Inputs%theNmod, Inputs%theNtot )  
     CALL GetMeshWts(kmesh, xxk, wwk)
 
     vlkmethod = ConstructMethodType( Inputs%theMethod, Inputs%theHerm, Inputs%theReg, Inputs%theRealSmooth, &
                                     Inputs%theIntSmooth, Inputs%theLambda)


     ALLOCATE(vnn(2*Inputs%theNtot, 2*Inputs%theNtot)); vnn = 0_dp     
     CALL GetVnn( nnchannel, Inputs%theNtot, xxk, wwk, vnn)  


     CALL ConstructVlowkType(kmesh, nnChannel, vlkmethod, vnn, vlowkstruct) 
     vlowk = 0_dp 
     CALL GetVlowk(vlowkstruct, vlowk(1:vlowkstruct%theVlowkSize, 1:vlowkstruct%theVlowkSize))

     IF( GetCoupled(nnChannel) )THEN
           vlowk_kmesh(1:vlowkstruct%theVlowkSize/2) = xxk(1:vlowkstruct%theVlowkSize/2)     
           vlowk_kwts(1:vlowkstruct%theVlowkSize/2) = wwk(1:vlowkstruct%theVlowkSize/2)  
     ELSE
           vlowk_kmesh(1:vlowkstruct%theVlowkSize) = xxk(1:vlowkstruct%theVlowkSize)     
           vlowk_kwts(1:vlowkstruct%theVlowkSize) = wwk(1:vlowkstruct%theVlowkSize)  
     ENDIF
  
     DEALLOCATE(xxk, wwk, vnn)
  END SUBROUTINE VlowkWrapper
 
  SUBROUTINE ConstructVlowkType(kmesh,channel,vlkmethod,vnn,vlowk)
     TYPE(VlowkType) :: vlowk
     TYPE(MeshType) :: kmesh
     TYPE(ChannelType) :: channel
     TYPE(MethodType) :: vlkmethod
     REAL(DP) :: vnn(:,:), fkkww
     INTEGER(I4B) :: nsizeV, nsizeMesh, i, j, ii, jj
     REAL(DP), ALLOCATABLE :: vherm(:,:), vnonherm(:,:)
     INTENT(IN)::kmesh,channel,vnn,vlkmethod
     INTENT(OUT)::vlowk


     nsizeV = GetNmod(kmesh) 
     nsizeMesh = GetNmod(kmesh)

     IF(GetMethod(vlkmethod).eq.3) THEN  !V_SRG defined on full momentum space
          nsizeV = GetNtot(kmesh) 
          nsizeMesh = GetNtot(kmesh)
     ENDIF

     IF(GetCoupled(channel))nsizeV=2*nsizeV
     Vlowk%theVlowkSize = nsizeV

     ALLOCATE(vherm(1:nsizeV,1:nsizeV),vnonherm(1:nsizeV,1:nsizeV)) 
     vherm=0_dp;vnonherm=0_dp

     Vlowk%theVlowk = 0_dp ; Vlowk%theNonHermVlowk = 0_dp
     CALL GetMeshWts(kmesh, Vlowk%theVlowkMesh(1:nsizeMesh), Vlowk%theVlowkWts(1:nsizeMesh))
 
     CALL calcvlowk(kmesh,channel,vlkmethod,vnn,vnonherm,vherm)
          Vlowk%theVlowkKmeshType = kmesh
          Vlowk%theVlowkMethod = vlkmethod 
          Vlowk%theVlowkChannel = channel
          Vlowk%theNonHermVlowk(1:nsizeV,1:nsizeV) = vnonherm(1:nsizeV,1:nsizeV)
          Vlowk%theVlowk(1:nsizeV,1:nsizeV) = vherm(1:nsizeV,1:nsizeV)

!    vlowk with weights and integration measure factors attached [units = MeV]
          DO i = 1, nsizeV
              DO j = 1, nsizeV
                   ii = i
                   jj = j
                   IF(i .gt. nsizeMesh) ii = i - nsizeMesh
                   IF(j .gt. nsizeMesh) jj = j - nsizeMesh
                   fkkww = sqrt(Vlowk%theVlowkWts(ii)*Vlowk%theVlowkWts(jj))* &
                            Vlowk%theVlowkMesh(ii)*Vlowk%theVlowkMesh(jj)      
                   Vlowk%theVlowk_w_kkww(i,j) = vherm(i,j)*fkkww
              ENDDO
         ENDDO
       
     DEALLOCATE(vherm,vnonherm)
  END SUBROUTINE ConstructVlowkType
   

   SUBROUTINE calcvlowk(kmesh,channel,vlkmethod,vnn,vnonherm,vherm)
     TYPE(MeshType) :: kmesh
     TYPE(ChannelType) :: channel
     TYPE(MethodType) :: vlkmethod
     REAL(DP) :: vnn(:,:), vnonherm(:,:),vherm(:,:)
     REAL(DP) :: srg_lam
       
  IF(GetMethod(vlkmethod).eq.1)THEN
         CALL leesuzukivlowk( kmesh, channel, vnn, vnonherm, vherm )
         IF(GetHerm(vlkmethod).ne.2)THEN  ! only for GS or Cholesky (Lee-Suzuki does Okubo internally)
            CALL hermitize( vlkmethod, kmesh, channel, vnonherm, vherm )
         ENDIF
  ELSE IF(GetMethod(vlkmethod).eq.2)THEN
         CALL smoothvlowk( vlkmethod, kmesh, channel, vnn, vnonherm )
         CALL hermitize( vlkmethod, kmesh, channel, vnonherm, vherm )
  ELSE IF(GetMethod(vlkmethod).eq.3)THEN

         srg_lam = GetLambda(kmesh) 
         CALL CalcVsrg(kmesh, srg_lam, vnn, channel, vherm)
         vherm = GetHb2m(channel)*vherm*2.0_dp/pi_d
      ENDIF


   END SUBROUTINE calcvlowk
  FUNCTION ConstructVlowkInputsType(meth, herm, ireg, nsmooth, rsmooth, n_tot, n_mod, k_max, lam, lam_eff, kvnn)
     TYPE(VlowkInputsType) :: ConstructVlowkInputsType
     INTEGER(I4B) :: meth, herm, ireg, nsmooth, n_tot, n_mod
     INTEGER(I4B), OPTIONAL :: kvnn
     REAL(DP) :: rsmooth, k_max, lam, lam_eff

     IF(PRESENT(kvnn))THEN 
     ConstructVlowkInputsType = VlowkInputsType(kvnn, meth, herm, ireg, nsmooth, rsmooth, n_tot, n_mod, k_max, &
                                               lam, lam_eff)
     ELSE
     ConstructVlowkInputsType = VlowkInputsType(0, meth, herm, ireg, nsmooth, rsmooth, n_tot, n_mod, k_max, &
                                               lam, lam_eff)
     END IF

  END FUNCTION ConstructVlowkInputsType

  FUNCTION GetVlowkKmeshType(vlowk)
     TYPE(MeshType) :: GetVlowkKmeshType
     TYPE(VlowkType) :: vlowk
     INTENT(IN) :: vlowk
   
     GetVlowkKmeshType = vlowk%theVlowkKmeshType     

     RETURN
  END FUNCTION GetVlowkKmeshType


  FUNCTION GetVlowkChannel(vlowk)
     TYPE(ChannelType) :: GetVlowkChannel
     TYPE(VlowkType) :: vlowk
     INTENT(IN) :: vlowk
   
     GetVlowkChannel = vlowk%theVlowkChannel     

     RETURN
  END FUNCTION GetVlowkChannel

  FUNCTION GetVlowkMethod(vlowk)
     TYPE(MethodType) :: GetVlowkMethod
     TYPE(VlowkType) :: vlowk
     INTENT(IN) :: vlowk
   
     GetVlowkMethod = vlowk%theVlowkMethod    

     RETURN
  END FUNCTION GetVlowkMethod

 
  SUBROUTINE GetVlowk(vlowk,vkk)
     TYPE(VlowkType) :: vlowk
     REAL(DP) :: vkk(:,:)
     INTEGER(I4B) :: nsize
     INTENT(IN) :: vlowk
     INTENT(OUT) :: vkk
     vkk = 0_dp 
     nsize = vlowk%theVlowkSize
     vkk(1:nsize,1:nsize) = vlowk%theVlowk(1:nsize,1:nsize)
  END SUBROUTINE GetVlowk

  SUBROUTINE GetVlowk_w_kkww(vlowk,vkk)
     TYPE(VlowkType) :: vlowk
     REAL(DP) :: vkk(:,:)
     INTEGER(I4B) :: nsize
     INTENT(IN) :: vlowk
     INTENT(OUT) :: vkk
     vkk = 0_dp 
     nsize = vlowk%theVlowkSize
     vkk(1:nsize,1:nsize) = vlowk%theVlowk_w_kkww(1:nsize,1:nsize)
  END SUBROUTINE GetVlowk_w_kkww

  SUBROUTINE GetNonHermVlowk(vlowk,vkk)
     TYPE(VlowkType) :: vlowk
     REAL(DP) :: vkk(:,:)
     INTEGER(I4B) :: nsize
     INTENT(IN) :: vlowk
     INTENT(OUT) :: vkk
     vkk = 0_dp 
     nsize = vlowk%theVlowkSize
     vkk(1:nsize,1:nsize) = vlowk%theNonHermVlowk(1:nsize,1:nsize)
  END SUBROUTINE GetNonHermVlowk
     

  INTEGER(I4B) FUNCTION GetVlowkSize(vlowk)
       TYPE(VlowkType) :: vlowk
       GetVlowkSize = vlowk%theVlowkSize
  END FUNCTION GetVlowkSize

   

END MODULE VlowkModule 
