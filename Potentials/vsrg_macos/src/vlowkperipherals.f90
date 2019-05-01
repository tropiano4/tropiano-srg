! 1.26.07 modified the ConstructMeshType subroutine to break the high k 
!         region into smaller subintervals.

! 6/05/06 Programmed by Scott K. Bogner
! This module is part of my first pass at re-writing the vlowk codes
! to implement some of the (primitive) object-oriented 
! features that Fortran 90/95 offers (e.g., derived-type data structures, 
! localization of data and data-hiding via modules, 
! procedure and operator overloading, etc.) 
! The current module defines 3 (this will grow as I think of more) 
! convenient data-structures that are used in the calculation of vlowk, and that
! are also repeatedly encountered in many applications of vlowk. 
! All individual components of the derived-types are declared PRIVATE, which is to 
! say one has to use the PUBLIC "constructor" and "accessor" functions and
! subroutines listed below to access the individual components in a given external
! procedure.

MODULE VlowkPeripherals
USE nrtype
IMPLICIT NONE
    PRIVATE
    INTEGER(I4B), PARAMETER :: nkdim = 200

    TYPE MeshType    ! derived type with all the (private) data specifying the momentum space mesh/wts
         PRIVATE
         INTEGER(I4B) :: theNkF    ! # of points from 0-thekF
         INTEGER(I4B) :: theNmod   ! # of points from 0-theLambdaeff
         INTEGER(I4B) :: theNtot   ! # of points from 0-theKmax
         REAL(DP) :: thekF         ! # fermi momentum (or any momentum) with thekF < theLambda 
         REAL(DP) :: theLambda     ! cutoff in vlowk
         REAL(DP) :: theLambdaeff  ! effective cutoff (e.g., for smooth cutoff, typically chosen as ~1.25*lambda)
         REAL(DP) :: theKmax       ! maximum momentum
         REAL(DP) :: theKmesh(1:nkdim) ! mesh vector
         REAL(DP) :: theWts(1:nkdim)   ! weight vector

    END TYPE MeshType

    TYPE ChannelType !derived type for the different nn-channels
        PRIVATE
        INTEGER(I4B) :: theNwave
        INTEGER(I4B) :: theL
        INTEGER(I4B) :: theS
        INTEGER(I4B) :: theJ
        INTEGER(I4B) :: theT
        INTEGER(I4B) :: theTz
        REAL(DP) :: theHb2m   
        LOGICAL(LGT) :: theCoupled
        LOGICAL(LGT) :: theDeut     
    END TYPE ChannelType

    TYPE MethodType!derived type for choice of calculational method, hermitization method, and regulator function
        PRIVATE
        INTEGER(I4B) :: theVlowkMethod    ! 1 = lee-suzuki (sharp cutoff only), 2 = 3-step  (smooth or sharp cutoff)       
        INTEGER(I4B) :: theHerm      ! 1 = Gram-Schmidt, 2 = Okubo, 3 = Cholesky, 4 = Kato
        INTEGER(I4B) :: theRegulator ! 1 = theta function, 2 = exponential, 3 = woods-saxon, 4 = tanh , 5 = powerlaw   
        REAL(DP) :: theRealWidth     ! real smoothness parameter (used by strutinsky, woods-saxon,...) 
        INTEGER(I4B) :: theIntWidth  ! int smoothness parameter (used by exponential, power-law, ...)
        REAL(DP) :: theLambda        ! cutoff
    END TYPE MethodType 
 

!  INTERFACES for overloaded procedures
    INTERFACE ConstructMeshType
         MODULE PROCEDURE ConstructMeshType, ConstructMeshType_medium, CosThetaMesh
    END INTERFACE ConstructMeshType 
         
    INTERFACE ConstructChannelType
        MODULE PROCEDURE ConstructChannel_nwave, ConstructChannel_lsjt
    END INTERFACE 

    INTERFACE GetMeshWts
        MODULE PROCEDURE GetMeshWts, GetNpointsMeshWts
    END INTERFACE GetMeshWts
!  List of PUBLIC subroutines and functions to construct the derived-type 
!  data structures and to access their individual components:

    PUBLIC :: MeshType, ConstructMeshType, GetMeshWts, GetNkF, &
              GetNmod, GetNtot, GetkF, GetKmax, GetLambdaeff, GetLambda

    PUBLIC :: ChannelType, ConstructChannelType, GetCoupled, &
              GetL, GetS, GetJ, GetT, GetTz, GetNwave, GetHb2m, GetDeut

    PUBLIC :: MethodType, ConstructMethodType, GetRealWidth,GetIntWidth,&
              GetHerm, GetMethod, smoothP, smoothQ, gauleg

CONTAINS

! Constructor function for MethodType derived-type structures 
  FUNCTION ConstructMethodType(imethod,iherm,ireg,rsmoothparam,nsmoothparam,lambda)
     TYPE(MethodType) :: ConstructMethodType
     INTEGER(I4B) :: imethod
     INTEGER(I4B) :: iherm, ireg, nsmoothparam
     REAL(DP) :: rsmoothparam,lambda
     INTENT(IN) :: imethod, iherm, rsmoothparam, nsmoothparam
     
     IF(imethod.eq.1.and.ireg.ne.1)ireg = 1 ! Lee-Suzuki only defined for sharp cutoffs
     ConstructMethodType = MethodType(imethod,iherm,ireg,rsmoothparam,nsmoothparam,lambda)

  END FUNCTION ConstructMethodType


!!!!!!!!!  BEGIN overloaded constructor subroutines to create MeshType structure KMESH !!!!!!!!!! 
! The various constructor subroutines are ALL invoked by CALL ConstructMeshType(...) -- i.e., use
! the generic procedure name. The compiler sorts out which implementation to use based on the 
! arguments in the given CALL statement.


   SUBROUTINE ConstructMeshType(kmesh,lambda,lambdaeff,kmax,nmod,ntot)
       TYPE(MeshType) :: kmesh
       REAL(DP), ALLOCATABLE :: xk(:), wk(:)
       REAL(DP) :: lambda,lambdaeff,kmax
       INTEGER(I4B) :: nmod,ntot
       INTENT(IN) :: lambda,lambdaeff,kmax,nmod,ntot
       INTENT(OUT) :: kmesh

        
       ALLOCATE(xk(ntot),wk(ntot))

       CALL gauleg(0._dp,lambdaeff,xk(1:nmod),wk(1:nmod),nmod) 
       CALL gauleg(lambdaeff,kmax,xk(nmod+1:ntot),wk(nmod+1:ntot),ntot-nmod) 

       kmesh%theWts = 0_dp ; kmesh%theKmesh = 0_dp

       kMesh%theNmod = nmod
       kMesh%theNtot = ntot
       kMesh%theLambdaeff = lambdaeff
       kMesh%theLambda = lambda
       kMesh%theKmax = kmax
       kMesh%theKmesh(1:ntot) = xk(1:ntot)
       kMesh%theWts(1:ntot) = wk(1:ntot)     
       DEALLOCATE(xk,wk)
       RETURN
    END SUBROUTINE ConstructMeshType   
   SUBROUTINE CosThetaMesh(kmesh,lambdaminus,lambda,lambdaeff,kmax,nmod,ntot)
       TYPE(MeshType) :: kmesh
       REAL(DP), ALLOCATABLE :: xk(:), wk(:)
       REAL(DP) :: lambda,lambdaeff,kmax, lambdaminus
       INTEGER(I4B) :: nmod,ntot
       INTENT(IN) :: lambda,lambdaeff,kmax,nmod,ntot,lambdaminus
       INTENT(OUT) :: kmesh

        
       ALLOCATE(xk(ntot),wk(ntot))

       CALL gauleg(lambdaminus,lambdaeff,xk(1:nmod),wk(1:nmod),nmod) 
       CALL gauleg(lambdaeff,kmax,xk(nmod+1:ntot),wk(nmod+1:ntot),ntot-nmod) 

       kmesh%theWts = 0_dp ; kmesh%theKmesh = 0_dp

       kMesh%theNmod = nmod
       kMesh%theNtot = ntot
       kMesh%theLambdaeff = lambdaeff
       kMesh%theLambda = lambda
       kMesh%theKmax = kmax
       kMesh%theKmesh(1:ntot) = xk(1:ntot)
       kMesh%theWts(1:ntot) = wk(1:ntot)     
       DEALLOCATE(xk,wk)
       RETURN
    END SUBROUTINE CosThetaMesh   
!  Slightly Embellished version (convenient for calcs. where you want a fixed number of
!  mesh points from 0-xkf, where xkf < lambda. For example, nuclear matter calculations...) 
   SUBROUTINE ConstructMeshType_medium(kmesh,xkf,lambda,lambdaeff,kmax,nkf,nmod,ntot)
       TYPE(MeshType) :: kmesh
       REAL(DP), ALLOCATABLE :: xk(:), wk(:)
       REAL(DP) :: lambda,lambdaeff,kmax,xkf
       INTEGER(I4B) :: nmod,ntot,nkf
       INTENT(IN) :: lambda,lambdaeff,kmax,nmod,ntot,xkf,nkf
       INTENT(OUT) :: kmesh

        
       ALLOCATE(xk(ntot),wk(ntot))
       
       CALL gauleg(0._dp,xkf,xk(1:nkf),wk(1:nkf),nkf) 
       CALL gauleg(xkf,lambdaeff,xk(nkf+1:nmod),wk(nkf+1:nmod),nmod-nkf) 
       CALL gauleg(lambdaeff,kmax,xk(nmod+1:ntot),wk(nmod+1:ntot),ntot-nmod) 

       kmesh%theWts = 0_dp ; kmesh%theKmesh = 0_dp

       kMesh%theNkf = nkf
       kMesh%theNmod = nmod
       kMesh%theNtot = ntot
       kMesh%thekF = xkf
       kMesh%theLambdaeff = lambdaeff
       kMesh%theLambda = lambda
       kMesh%theKmax = kmax
       kMesh%theKmesh(1:ntot) = xk(1:ntot)
       kMesh%theWts(1:ntot) = wk(1:ntot)     

       DEALLOCATE(xk,wk)
       RETURN
    END SUBROUTINE ConstructMeshType_medium   

!!!!!!!!!!! END of overloaded ConstructMeshType subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!!!!!!!!!! BEGIN PUBLIC accessor functions for PRIVATE components of TYPE(Mesh) structures !!!!!!!!!!!!!!!!!!!!
  SUBROUTINE GetMeshWts(kmesh,xk,wk)
      TYPE(MeshType) :: kmesh
      REAL(DP) :: xk(:),wk(:)
      INTEGER(I4B) :: n
      INTENT(IN) :: kmesh
      INTENT(OUT) :: xk,wk

      n = kmesh%theNtot
      xk(1:n) = kmesh%theKmesh(1:n)
      wk(1:n) = kmesh%theWts(1:n)
  END SUBROUTINE GetMeshWts

  SUBROUTINE GetNpointsMeshWts(kmesh,n,xk,wk)
      TYPE(MeshType) :: kmesh
      REAL(DP) :: xk(:),wk(:)
      INTEGER(I4B) :: n
      INTENT(IN) :: kmesh,n
      INTENT(OUT) :: xk,wk

!      n = kmesh%theNtot
      xk(1:n) = kmesh%theKmesh(1:n)
      wk(1:n) = kmesh%theWts(1:n)
  END SUBROUTINE GetNpointsMeshWts

  INTEGER(I4B) FUNCTION GetNkF(kmesh)
      TYPE(MeshType) :: kmesh
      INTENT(IN) :: kmesh
      GetNkF = kmesh%theNkF
  END FUNCTION GetNkF
 
  INTEGER(I4B) FUNCTION GetNmod(kmesh)
      TYPE(MeshType) :: kmesh
      INTENT(IN) :: kmesh
      GetNmod = kmesh%theNmod
  END FUNCTION GetNmod 

  INTEGER(I4B) FUNCTION GetNtot(kmesh)
      TYPE(MeshType) :: kmesh
      INTENT(IN) :: kmesh
      GetNtot = kmesh%theNtot
  END FUNCTION GetNtot 

  REAL(DP) FUNCTION GetkF(kmesh)
      TYPE(MeshType) :: kmesh
      INTENT(IN) :: kmesh
      GetkF = kmesh%thekF
  END FUNCTION GetkF

  REAL(DP) FUNCTION GetLambdaeff(kmesh)
      TYPE(MeshType) :: kmesh
      INTENT(IN) :: kmesh
      GetLambdaeff = kmesh%theLambdaeff
  END FUNCTION GetLambdaeff

  REAL(DP) FUNCTION GetLambda(kmesh)
      TYPE(MeshType) :: kmesh
      INTENT(IN) :: kmesh
      GetLambda = kmesh%theLambda
  END FUNCTION GetLambda

  REAL(DP) FUNCTION GetKmax(kmesh)
      TYPE(MeshType) :: kmesh
      INTENT(IN) :: kmesh
      GetKmax = kmesh%theKmax
  END FUNCTION GetKmax
!!!!!!!!!!! END public accessor functions for PRIVATE components of TYPE(Mesh) derived-type !!!!!!!!!!!!!!!!!!!!

  SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER(I4B) :: n,i,j,m
      REAL(DP) :: x1,x2,x(n),w(n),EPS
      PARAMETER (EPS=3.e-14_dp)
      REAL(DP) :: p1,p2,p3,pp,xl,xm,z,z1
      INTENT(IN) :: x1,x2,n
      INTENT(OUT) :: x,w                                                                                                                                      
      m=(n+1)/2
      xm=0.5_dp*(x2+x1)
      xl=0.5_dp*(x2-x1)
      do 12 i=1,m
        z=cos(pi_d*(i-.25_dp)/(n+.5_dp))
1       continue
          p1=1.0_dp
          p2=0.0_dp
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.0_dp)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.0_dp*xl/((1.0_dp-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      END SUBROUTINE gauleg


    FUNCTION ConstructChannel_nwave(nwave,itz)
        TYPE(ChannelType) :: ConstructChannel_nwave
        INTEGER(I4B) :: nwave,itz,jt,l,s,it,la(29),ja(29),isa(29) 
        LOGICAL(LGT) :: coupled,deut
        REAL(DP) :: hb2m 
        INTENT(IN) :: nwave,itz
        DATA la /0,0,1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9/
        DATA ja /0,1,1,0,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9/
        DATA isa/0,1,0,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1/
             
             if(itz.eq.0)hb2m = 41.47105
             if(itz.eq.1)hb2m = 41.49960
             if(itz.eq.-1)hb2m = 41.4425
             
             l=la(nwave);jt=ja(nwave);s=isa(nwave)
             it=1
             if((-1)**(it+l+s).ge.0)it=0

             coupled=.false.
             if(s.eq.1.and.jt.gt.l)coupled=.true. 
            
             deut=.false.
             if(nwave.eq.2.and.itz.eq.0)deut=.true. 
            
             ConstructChannel_nwave = ChannelType(nwave,l,s,jt,it,itz,hb2m,coupled,deut)              
        RETURN
    END FUNCTION ConstructChannel_nwave

    FUNCTION ConstructChannel_lsjt(l,s,jt,itz)
        TYPE(ChannelType) :: ConstructChannel_lsjt
        INTEGER(I4B) :: l,s,jt,itz,it,nwave 
        LOGICAL(LGT) :: coupled,deut
        REAL(DP) :: hb2m 
!        DATA la /0,0,1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9/
!        DATA ja /0,1,1,0,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9/
!        DATA isa/0,1,0,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1/

        IF(l.eq.0)THEN
           IF(jt.eq.0)nwave = 1
           IF(jt.eq.1)nwave = 2
        ELSE IF (l.eq.1)THEN
           IF(s.eq.0)nwave = 3
           IF(s.eq.1)THEN
                 IF(jt.eq.0)nwave = 4
                 IF(jt.eq.1)nwave = 5
                 IF(jt.eq.2)nwave = 6
           ENDIF
       ELSE IF (l.eq.2)THEN
           IF(s.eq.0)nwave = 7
           IF(s.eq.1)THEN
                 IF(jt.eq.2)nwave = 8
                 IF(jt.eq.3)nwave = 9
           ENDIF
       ELSE IF (l.eq.3)THEN
           IF(s.eq.0)nwave = 10
           IF(s.eq.1)THEN
               IF(jt.eq.3)nwave = 11
               IF(jt.eq.4)nwave = 12
           ENDIF
       ELSE IF (l.eq.4)THEN
           IF(s.eq.0)nwave = 13
           IF(s.eq.1)THEN
               IF(jt.eq.4)nwave = 14
               IF(jt.eq.5)nwave = 15
           ENDIF
       ELSE IF(l.eq.5)THEN
           IF(s.eq.0)nwave = 16
           IF(s.eq.1)THEN
               IF(jt.eq.5)nwave = 17
               IF(jt.eq.6)nwave = 18
           ENDIF
       ELSE IF(l.eq.6)THEN
           IF(s.eq.0)nwave = 19
           IF(s.eq.1)THEN
              IF(jt.eq.6)nwave = 20
              IF(jt.eq.7)nwave = 21
           ENDIF
       ELSE IF(l.eq.7)THEN
           IF(s.eq.0)nwave = 22
           IF(s.eq.1)THEN
              IF(jt.eq.7)nwave = 23
              IF(jt.eq.8)nwave = 24
           ENDIF
       ELSE IF(l.eq.8)THEN
           IF(s.eq.0)nwave = 25
           IF(s.eq.1)THEN
              IF(jt.eq.8)nwave = 26
              IF(jt.eq.9)nwave = 27
           ENDIF
       ELSE IF(l.eq.9)THEN
           IF(s.eq.0)nwave = 28
           IF(s.eq.1)THEN
              IF(jt.eq.9)nwave = 29
              IF(jt.eq.10)nwave = 30
           ENDIF
       ELSE IF(l.eq.10)THEN
           IF(s.eq.0)nwave = 31
           IF(s.eq.1)THEN
              IF(jt.eq.10)nwave = 32
              IF(jt.eq.11)nwave = 33
           ENDIF
       ELSE IF(l.eq.11)THEN
           IF(s.eq.0)nwave = 34
           IF(s.eq.1)THEN
              IF(jt.eq.11)nwave = 35
              IF(jt.eq.12)nwave = 36
           ENDIF
       ELSE IF(l.eq.12)THEN
           IF(s.eq.0)nwave = 37
           IF(s.eq.1)THEN
              IF(jt.eq.12)nwave = 38
              IF(jt.eq.13)nwave = 39
           ENDIF
       ENDIF

        if(itz.eq.0)hb2m = 41.47105
        if(itz.eq.1)hb2m = 41.49960
        if(itz.eq.-1)hb2m = 41.4425

        it=1
        if((-1)**(it+l+s).ge.0)it=0

        coupled=.false.
        if(s.eq.1.and.jt.gt.l)coupled=.true. 

        deut=.false.
        if(s.eq.1.and.l.eq.0.and.jt.eq.1.and.itz.eq.0)deut=.true. 

        ConstructChannel_lsjt = ChannelType(nwave,l,s,jt,it,itz,hb2m,coupled,deut)
        RETURN
    END FUNCTION ConstructChannel_lsjt

    INTEGER(I4B) FUNCTION GetNwave(channel)
        TYPE(ChannelType) :: channel
        GetNwave = channel%theNwave
    END FUNCTION GetNwave

    INTEGER(I4B) FUNCTION GetTz(channel)
        TYPE(ChannelType) :: channel
        GetTz = channel%theTz
    END FUNCTION GetTz   

    LOGICAL(LGT) FUNCTION GetCoupled(channel)
        TYPE(ChannelType) :: channel
        GetCoupled = channel%theCoupled
        RETURN
    END FUNCTION GetCoupled 
  
    LOGICAL(LGT) FUNCTION GetDeut(channel)
        TYPE(ChannelType) :: channel
        GetDeut = channel%theDeut
        RETURN
    END FUNCTION GetDeut
 
    REAL(DP) FUNCTION GetHb2m(channel)
        TYPE(ChannelType) :: channel
        GetHb2m = channel%theHb2m
        RETURN
    END FUNCTION GetHb2m   
    
    INTEGER(I4B) FUNCTION GetL(channel)
        TYPE(ChannelType) :: channel
        GetL = channel%theL
        RETURN
    END FUNCTION GetL

    INTEGER(I4B) FUNCTION GetS(channel)
        TYPE(ChannelType) :: channel
        GetS = channel%theS
        RETURN
    END FUNCTION GetS

    INTEGER(I4B) FUNCTION GetJ(channel)
        TYPE(ChannelType) :: channel
        GetJ = channel%theJ
        RETURN
    END FUNCTION GetJ

    INTEGER(I4B) FUNCTION GetT(channel)
        TYPE(ChannelType) :: channel
        GetT = channel%theT
        RETURN
    END FUNCTION GetT



 REAL(DP) FUNCTION GetRealWidth(vlkmethod)
     TYPE(MethodType) :: vlkmethod
     GetRealWidth = vlkmethod%theRealWidth
     RETURN
 END FUNCTION GetRealWidth

 INTEGER(I4B) FUNCTION GetIntWidth(vlkmethod)
     TYPE(MethodType) :: vlkmethod
     GetIntWidth = vlkmethod%theIntWidth
     RETURN
 END FUNCTION GetIntWidth

 INTEGER(I4B) FUNCTION GetHerm(vlkmethod)
     TYPE(MethodType) :: vlkmethod
     GetHerm = vlkmethod%theHerm
 END FUNCTION GetHerm

 INTEGER(I4B) FUNCTION GetMethod(vlkmethod)
     TYPE(MethodType) :: vlkmethod
     GetMethod = vlkmethod%theVlowkMethod
 END FUNCTION GetMethod


! Smooth regulator function
 REAL(DP) FUNCTION SmoothP(x,vlkmethod)
        TYPE(MethodType) :: vlkmethod
        REAL(DP) :: x,lambda,eps
        INTEGER(I4B) :: n,ireg
        INTENT(IN) :: vlkmethod, x
 
        ireg = vlkmethod%theRegulator 
        n = vlkmethod%theIntWidth
        lambda = vlkmethod%theLambda
        eps = vlkmethod%theRealWidth
        
        IF(ireg.eq.1)THEN     !theta
             SmoothP = 0_dp
             IF(x.le.lambda)SmoothP = 1._dp
             RETURN
        ELSE IF(ireg.eq.2)THEN  !exponential
             SmoothP = dexp(-(x/lambda)**(2.0_dp*n)) 
             RETURN
        ELSE IF(ireg.eq.3)THEN  ! woods-saxon
             SmoothP=1.0_dp/(1.0_dp+dexp((x**2-lambda**2)/eps**2)) 
             RETURN
        ELSE IF(ireg.eq.4)THEN  ! hyperbolic tangent
             SmoothP=0.5_dp*(1.0_dp +tanh((lambda**2-x**2)/(lambda*x*eps)))
             RETURN
        ELSE IF(ireg.eq.5)THEN  ! power-law 
             SmoothP=1.0_dp/(1.0_dp +(x/lambda)**n) 
             RETURN
        ELSE
             write(11,*)'problem in SMoothP ireg =', ireg
             call abort
        ENDIF 
        RETURN
 END FUNCTION smoothp
 
! "complement" 1-SmoothP**2 regulator function
 REAL(DP) FUNCTION SmoothQ(x,vlkmethod)
        TYPE(MethodType) :: vlkmethod
        REAL(DP) :: x
        INTENT(IN) :: x, vlkmethod

        SmoothQ = 1_dp - SmoothP(x,vlkmethod)**2
 
 END FUNCTION SmoothQ


END MODULE VlowkPeripherals



