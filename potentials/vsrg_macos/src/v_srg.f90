MODULE v_srg
 USE nrtype
 USE VlowkPeripherals

 IMPLICIT NONE
 PRIVATE

     INTEGER(I4B), PARAMETER :: iode_solver = 1 !(1 = pred/corr, 2 = odepack (not active), 3-5 different rk methods)
     REAL(DP), PARAMETER :: rtol = 1.0e-8_dp, atol = 1.0e-8_dp

     REAL(DP) :: finaldeut, baredeut, srgdeut, hb2m
     REAL(DP), ALLOCATABLE :: xk(:), wk(:)
     INTEGER(I4B) :: n_ode_eqns, nsizeV, nkpoint
     
     PUBLIC :: CalcVsrg

CONTAINS

SUBROUTINE CalcVsrg(kmesh, srg_lam, vnn, channel, vsrg)
   USE ode_shampine_and_gordon
   USE rksuite
   TYPE(MeshType) :: kmesh
   TYPE(ChannelType) :: channel
   REAL(DP) :: vnn(:,:), vsrg(:,:), kkww, srg_lam, h_0, hstart, tstopped
   REAL(DP), ALLOCATABLE :: v_vec(:), ham(:,:), eval(:), work(:), RWORK(:),vp_vec(:), vmax_vec(:),thres(:)
   INTEGER(I4B) ::  i, j, ii, jj, ij, nok, nbad, info, lwork, indx(400)
   INTEGER(I4B) :: ITOL, ITASK, ISTATE, IOPT, LRW, LIW, MF
   INTEGER, ALLOCATABLE :: IWORK(:)
   LOGICAL(LGT) :: MESSAGE, ERRASS

!!!!!!!!!!!! set module work variables !!!!!!!!!!
   hb2m = GetHb2m(channel)
   nkpoint = GetNtot(kmesh)       
   nsizeV  = nkpoint              
   IF(GetCoupled(Channel))nsizeV = 2*nkpoint
   n_ode_eqns = nsizeV*(nsizeV+1)/2

   ALLOCATE(xk(nkpoint), wk(nkpoint))
   CALL GetMeshWts(kmesh, xk, wk)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ALLOCATE(v_vec(n_ode_eqns)); v_vec = 0_dp
  

! put vnn(k,k') into a 1-d vector
     ij = 0
     DO j = 1, nsizeV 
         DO i = j, nsizeV
             ii = i ; jj = j
             IF(i .gt. nkpoint) ii = i - nkpoint
             IF(j .gt. nkpoint) jj = j - nkpoint 
             kkww = sqrt(wk(ii)*wk(jj))*xk(ii)*xk(jj)
             ij = ij+1
             v_vec(ij) = vnn(i,j)/kkww/hb2m*pi_d/2.0_dp
         ENDDO
     ENDDO


!      call ODE stepper

SELECT CASE(iode_solver)  
   CASE(1)     !shampine_gordon code using predictor/corrector method
        lrw  = 100 + 21*n_ode_eqns    ! length of rwork(:) array 
        liw   = 5                     !length of iwork(:) array
        istate = 1  !normal operation
        ALLOCATE(rwork(lrw), iwork(liw))

        CALL ode(derivs, n_ode_eqns, v_vec, GetKmax(kmesh), srg_lam, rtol, atol, istate, rwork, iwork)
        WRITE(*,*)'istate',istate
        IF(istate.ne.2)write(11,*)'problem in ode, istate=',istate
        IF(istate.ne.2)write(*,*)'problem in ode, istate=',istate
        DEALLOCATE(RWORK, IWORK)

   CASE(2)  !ODEPACK code using  predictor/corrector method with scalar error tolerances 

        itol =  1   ! specifies scalar absolute error tolerances
        itask = 1 
        lrw  = 20 + 16*n_ode_eqns    ! length of rwork(:) array 
        liw   = 20                   !length of iwork(:) array
        ALLOCATE(rwork(lrw), iwork(liw))
        iwork(6) = 500000  !increase max # steps
        iopt = 1    ! 0 => no optional inputs used, 1 => optional inputs used
        istate = 1  !1st call 
        mf = 10   !non-stiff predictor/corrector
!        CALL DLSODE (xderivs, n_ode_eqns, v_vec, GetKmax(kmesh), srg_lam, ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, &
!                 RWORK, LRW, IWORK, LIW, JACOBIAN, MF)
        DEALLOCATE(RWORK, IWORK)
        IF(istate.ne.2)WRITE(*,*)'problem in dlsode istate=',istate 

    CASE(3:5) !rksuite Runge-Kutta Methods 
         SELECT CASE (iode_solver)
              CASE(3)        !RK (2,3) method
                  MF = 1
              CASE(4)        !RK (4,5) method
                  MF = 2
              CASE(5)        !RK (7,8) method
                  MF = 3
         END SELECT          

         ALLOCATE(THRES(n_ode_eqns), vp_vec(n_ode_eqns), vmax_vec(n_ode_eqns))
         THRES(:) = 1.0e-4     ! threshold for ith element of v_vec
         ERRASS = .false.      ! estimate true error (costly)  
         HSTART = 0.0_dp       ! 1st step
         LRW = 32*n_ode_eqns   ! size of work array
         MESSAGE = .true.      ! output diagnostic messages to std output
         ALLOCATE(rwork(lrw))
         tstopped = GetKmax(kmesh)
!         CALL SETUP(n_ode_eqns,GetKmax(kmesh),v_vec,srg_lam,RTOL,THRES,MF,'U',ERRASS,HSTART,RWORK,LRW,MESSAGE)
!         CALL UT(rkderivs,srg_lam,tstopped,v_vec,vp_vec,vmax_vec,RWORK,istate)
         WRITE(*,*)'istate',istate
         WRITE(11,*)'ran rksuite', srg_lam, tstopped
         DEALLOCATE (THRES, rwork, vp_vec, vmax_vec)
END SELECT





! put evolved v_vec(:) back into matrix form vsrg(:,:)
     ij = 0
     DO j = 1, nsizeV 
         DO i = j, nsizeV 
              ij = ij+1
               vsrg(i,j) = v_vec(ij) 
               vsrg(j,i) = v_vec(ij) 
         ENDDO
     ENDDO
            
     IF(GetCoupled(channel).and.(GetL(channel).eq.0))THEN   ! check deuteron B.E.
         ALLOCATE(ham(nsizeV, nsizeV), eval(nsizeV))
         lwork = 6*nsizeV-1; info = 100; ALLOCATE(work(lwork))

!        set up SRG hamiltonian to be diagonalized
         DO i = 1, nsizeV 
            DO j = 1, nsizeV 
                 ii = i; jj = j
                 IF(i .gt. nkpoint) ii = i - nkpoint 
                 IF(j .gt. nkpoint) jj = j - nkpoint 
                 h_0 = 0.0_dp
                 IF(i.eq.j) h_0 = hb2m*xk(ii)**2
                 kkww = hb2m*xk(ii)*xk(jj)*sqrt(wk(ii)*wk(jj))
                 ham(i, j) = h_0 + kkww*vsrg(i,j)*2.0_dp/pi_d
            ENDDO
         ENDDO 

         CALL dsyev('V','U',nsizeV,ham,nsizeV,eval,work,lwork,info)

         IF(info.ne.0)THEN
               write(*,*)'PROBLEM IN DSYEV INSIDE SRGVLOWK'
               CALL abort
         ENDIF
         
!        CALL indexx(nkpoint,eval,indx)
        srgdeut = eval(1)
        finaldeut = srgdeut
!        WRITE(14,'(2f15.8)')(xk(i), ham(i,indx(1)), i = 1, ntot)
!        WRITE(14,*) ' '  
!       repeat for the "bare" hamiltonian

         ham = 0_dp; eval = 0_dp; work = 0_dp
         DO i = 1, nsizeV 
            DO j = 1, nsizeV
                 ii = i; jj = j
                 IF(i .gt. nkpoint) ii = i - nkpoint 
                 IF(j .gt. nkpoint) jj = j - nkpoint 
                 h_0 = 0_dp
                 IF(i.eq.j)h_0 = hb2m*xk(ii)**2
                 ham(i, j) = h_0 + vnn(i,j)
            ENDDO
         ENDDO 

         CALL dsyev('V','U',nsizeV,ham,nsizeV,eval,work,lwork,info)

         IF(info.ne.0)THEN
               write(*,*)'PROBLEM IN DSYEV INSIDE SRGVLOWK'
               CALL abort
         ENDIF

 !       CALL indexx(nkpoint,eval,indx)
       baredeut = eval(1)
        write(11,*)'bare srg deut', baredeut, srgdeut
! EDITED BY AT: Include 5 lowest eigenvalues in case spurious bound state is lower than deuteron BE
        write(11,*)'lowest eigenvalues', eval(1), eval(2), eval(3), eval(4), eval(5)
        DEALLOCATE(work, ham, eval)
     ENDIF

     DEALLOCATE( v_vec , xk, wk)
END SUBROUTINE CalcVsrg


  SUBROUTINE xderivs(neqn, lam, v, dv)
!  SUBROUTINE derivs(lam, v, dv)
      USE nrtype
      INTEGER(I4B) ::  neqn, ij, i, j, ii, jj, kk, k
      REAL(DP) :: ans, lam, v(neqn), dv(neqn), k1, k2, p
      REAL(DP), ALLOCATABLE :: vij(:,:)

      ALLOCATE(vij(nsizeV,nsizeV)) 

      ij = 0
      DO j = 1, nsizeV 
         DO i = j, nsizeV 
             ij = ij + 1
             vij(i,j) = v(ij)
             vij(j,i) = v(ij)
         ENDDO
      ENDDO

      ij = 0
      DO j = 1, nsizeV 
         DO i = j, nsizeV 
            ii = i; jj = j    
            IF(i .gt. nkpoint) ii = i - nkpoint 
            IF(j .gt. nkpoint) jj = j - nkpoint 
            ij = ij + 1
            ans = 0_dp
                 DO k = 1, nsizeV 
                      kk = k
                      IF(k .gt. nkpoint) kk = k - nkpoint 
                      k1 = xk(ii)**2 
                      k2 = xk(jj)**2 
                      p  = xk(kk)**2 
                      ans = ans +  2.0_dp/pi_d*(k1+k2-2.0_dp*p)*xk(kk)**2*wk(kk)* &
                            vij(i,k)*vij(k,j) 
                 ENDDO
            dv(ij) =-4.0_dp/lam**5*( ans - (xk(ii)**2-xk(jj)**2)*(k1-k2)*vij(i,j))
         ENDDO
      ENDDO
     
      DEALLOCATE(vij)
  END SUBROUTINE xderivs     

  SUBROUTINE rkderivs(lam, v, dv)
      USE nrtype
      INTEGER(I4B) ::  neqn, ij, i, j, ii, jj, kk, k
      REAL(DP) :: ans, lam, v(*), dv(*), k1, k2, p
      REAL(DP), ALLOCATABLE :: vij(:,:)

      ALLOCATE(vij(nsizeV,nsizeV)) 

      ij = 0
      DO j = 1, nsizeV 
         DO i = j, nsizeV 
             ij = ij + 1
             vij(i,j) = v(ij)
             vij(j,i) = v(ij)
         ENDDO
      ENDDO

      ij = 0
      DO j = 1, nsizeV 
         DO i = j, nsizeV 
            ii = i; jj = j    
            IF(i .gt. nkpoint) ii = i - nkpoint 
            IF(j .gt. nkpoint) jj = j - nkpoint 
            ij = ij + 1
            ans = 0_dp
                 DO k = 1, nsizeV 
                      kk = k
                      IF(k .gt. nkpoint) kk = k - nkpoint 
                      k1 = xk(ii)**2 
                      k2 = xk(jj)**2 
                      p  = xk(kk)**2 
                      ans = ans +  2.0_dp/pi_d*(k1+k2-2.0_dp*p)*xk(kk)**2*wk(kk)* &
                            vij(i,k)*vij(k,j) 
                 ENDDO
            dv(ij) =-4.0_dp/lam**5*( ans - (xk(ii)**2-xk(jj)**2)*(k1-k2)*vij(i,j))
         ENDDO
      ENDDO
     
      DEALLOCATE(vij)
  END SUBROUTINE rkderivs     

!  SUBROUTINE derivs(neqn, lam, v, dv)
  SUBROUTINE derivs(lam, v, dv)
      USE nrtype
      INTEGER(I4B) ::  neqn, ij, i, j, ii, jj, kk, k
      REAL(DP) :: ans, lam, v(:), dv(:), k1, k2, p
      REAL(DP), ALLOCATABLE :: vij(:,:)

      ALLOCATE(vij(nsizeV,nsizeV)) 

      ij = 0
      DO j = 1, nsizeV 
         DO i = j, nsizeV 
             ij = ij + 1
             vij(i,j) = v(ij)
             vij(j,i) = v(ij)
         ENDDO
      ENDDO

      ij = 0
      DO j = 1, nsizeV 
         DO i = j, nsizeV 
            ii = i; jj = j    
            IF(i .gt. nkpoint) ii = i - nkpoint 
            IF(j .gt. nkpoint) jj = j - nkpoint 
            ij = ij + 1
            ans = 0_dp
                 DO k = 1, nsizeV 
                      kk = k
                      IF(k .gt. nkpoint) kk = k - nkpoint 
                      k1 = xk(ii)**2 
                      k2 = xk(jj)**2 
                      p  = xk(kk)**2 
                      ans = ans +  2.0_dp/pi_d*(k1+k2-2.0_dp*p)*xk(kk)**2*wk(kk)* &
                            vij(i,k)*vij(k,j) 
                 ENDDO
            dv(ij) =-4.0_dp/lam**5*( ans - (xk(ii)**2-xk(jj)**2)*(k1-k2)*vij(i,j))
         ENDDO
      ENDDO
     
      DEALLOCATE(vij)
  END SUBROUTINE derivs     


  SUBROUTINE JACOBIAN (NEQ, T, Y, ML, MU, PD, NROWPD)
      USE nrtype
      INTEGER(I4B) :: NEQ, ML, MU, NROWPD
      REAL(DP) :: T, Y(:), PD(:)

  END SUBROUTINE jacobian 
END MODULE v_srg 
