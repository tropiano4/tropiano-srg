MODULE vlowk_smooth_cutoff
 USE nrtype
 USE VlowkPeripherals
 USE interpolation, ONLY:interp2d,vkku399lag

 IMPLICIT NONE
 PRIVATE

     INTEGER(I4B), PARAMETER :: ndim = 400, nit = 10 
     REAL(DP), PARAMETER :: svd_eps = 1.0e-8_dp 
     TYPE(MethodType) :: reg 
     REAL(DP) :: p_on, edeut, lambda, lambdaeff, kmax, finaldeut, hb2m 
     REAL(DP),ALLOCATABLE :: xk(:), wk(:), p_in(:), pwt(:), biorth(:,:)
     INTEGER(I4B) :: nkpt, nmod, np_on, nktot 
     LOGICAL(LGT) :: coupled, deut

     
  PUBLIC :: smoothvlowk
CONTAINS

SUBROUTINE smoothvlowk(vlkmethod, kmesh, channel, xvnn, vnh)
     TYPE(MeshType) :: kmesh, pmesh
     TYPE(ChannelType) :: channel
     TYPE(MethodType) :: vlkmethod
     REAL(DP) :: xvnn(:,:), vnh(:,:),  edeut_sc, fkkw
     REAL(DP),ALLOCATABLE :: veff(:,:), p_sc(:), vec_sc(:,:), veffE(:,:,:), vnn(:,:)
     INTEGER(I4B)::i, j, k, ip, ii, jj, nn
     INTENT(IN):: xvnn, kmesh, vlkmethod, channel
     INTENT(OUT):: vnh

!!!!!!! set module working variables and allocate module arrays !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      reg = vlkmethod
      coupled = GetCoupled(channel)
      hb2m = GetHb2m(channel)
      deut = GetDeut(channel)

      nkpt = GetNtot(kmesh) 
      nmod = GetNmod(kmesh)
      lambda = GetLambda(kmesh)
      lambdaeff = GetLambdaeff(kmesh)
      kmax = GetKmax(kmesh)

      np_on = nmod; nktot = nkpt
      if(coupled)np_on = 2*nmod
      if(coupled)nktot = 2*nkpt

      IF(ALLOCATED(xk))DEALLOCATE(xk)
      IF(ALLOCATED(wk))DEALLOCATE(wk)
      ALLOCATE(xk(1:nkpt),wk(1:nkpt)); xk=0_dp ; wk=0_dp
      CALL GetMeshWts(kmesh,xk,wk)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE(vnn(size(xvnn(1:,1)),size(xvnn(1:,1))))
      vnn=xvnn

!     generate the onshell momentum mesh p_in(1:np_on) for the bloch-horowitz eqn.
      Call ConstructMeshType(pmesh,lambda,lambdaeff,lambdaeff,np_on,np_on)
        nn = GetNtot(pmesh)
        allocate(p_in(1:nn),pwt(1:nn))
        call GetMeshWts(pmesh,p_in,pwt)
        p_in(1:nn)=p_in(1:nn)+.75_dp*xk(1)                                

!       remove factors of sqrt(wt_i*wt_j)*k_i*k_j from vnn and change units to fm by dividing out hb2m
            do i=1,nktot
               do j=1,nktot
                 ii=i;jj=j
                 if(i.gt.nkpt)ii=i-nkpt
                 if(j.gt.nkpt)jj=j-nkpt 
                 fkkw=hb2m*dsqrt(wk(ii)*wk(jj))*xk(ii)*xk(jj)*2.0_dp/pi_d
                 vnn(i,j)=vnn(i,j)/fkkw
               enddo
            enddo

        allocate(veff(np_on,np_on),veffE(np_on,np_on,np_on),vec_sc(np_on,np_on),p_sc(np_on));veffE=0_dp;vec_sc=0_dp;p_sc=0_dp
         
!       calculate energy-dependent veff and self-consistent eigenvectors/energies 
        do ip=1,np_on    !note that np_on=nmod for uncoupled and 2*nmod for coupled channel
              veff=0.0_dp
              p_on=p_in(ip)
              call blochhorwtz(lambda,vnn,veff)
              call eigenvectors(lambda,veff,vnn,p_sc(ip),vec_sc(1:np_on,ip))
              veffE(:,:,ip)=veff
        enddo
!       treat the E<0 case (i.e., the deuteron) separately 
        if(deut)then
              edeut=-2.224/hb2m
              veff=0.0_dp
              call blochhorwtz_deut(lambda,vnn,veff)
              call eigenvectors_deut(lambda,veff,vnn,edeut_sc,vec_sc(1:np_on,1))
              veffE(:,:,1)=veff
        endif

!           DO i = 1,np_on
!             if(p_sc(i).le.2.0)then
!             write(88,*)'psc',p_sc(i) 
!             write(88,'(2f14.5)')(xk(j),veffE(j,j,i),j=1,np_on/2)
!            write(88,*) ' '
!            endif
!         ENDDO
!       attach ki*kj*sqrt(wi*wj) factors to the matrix elements
           do i=1,np_on
              do j=1,np_on
                 ii=i;jj=j
                 if(i.gt.nmod)ii=i-nmod
                 if(j.gt.nmod)jj=j-nmod
                 fkkw=dsqrt(wk(ii)*wk(jj))*xk(ii)*xk(jj)
                 veffE(i,j,:)=veffE(i,j,:)*fkkw
              enddo
           enddo

!        write(11,'(a,3f14.6)')('k   p_in   p_sc', xk(i), p_in(i), p_sc(i),i=1,np_on)


!       calculate energy-independent vlowk
        call vlowknh(p_sc,vec_sc,vnh,veffE)  !non-hermitian

!        write(11,*)'after vlowknh'
!        IF(GetHerm(reg).eq.1)THEN
!           call gramschmidtherm(vec_sc,vnh,vherm)        !hermitian   
!        ELSE
     
!        ENDIF
!         call hermitize(vec_sc,vnh,vherm)
!        call recalcvlowknh(vec_sc,vnh,veffE) 
!         call hermitize(vec_sc,vnh,vherm)
!      vherm = vherm*2.0_dp/pi_d*hb2m
!      vnh = vnh*2.0_dp/pi_d*hb2m

     DEALLOCATE(veff, veffE, vec_sc, p_sc, vnn)
     DEALLOCATE(xk, wk, p_in, pwt, biorth)
END SUBROUTINE smoothvlowk
    

   SUBROUTINE vlowknh(p_sc,vec_sc,vnh,veff)
        USE nrtype
        IMPLICIT NONE
        REAL(DP),ALLOCATABLE::u(:,:),w(:),v(:,:),work(:) !,ui(:,:),wi(:),vi(:,:)
        INTEGER(I4B)::info,i,j,ii,jj,lwork
        REAL(DP)::p_sc(:),vec_sc(:,:),vnh(:,:),veff(:,:,:),fkkw
        INTENT(IN)::p_sc,veff,vec_sc
        INTENT(OUT)::vnh

        if(allocated(biorth))deallocate(biorth)
        allocate(biorth(np_on,np_on),u(np_on,np_on),v(np_on,np_on),w(np_on));biorth=0.0_dp;u=0.0_dp;v=0.0_dp;w=0.0_dp

!   use SV decomposition to build bi-orthogonal complement vectors
               u=vec_sc(:,:);lwork=40*np_on;allocate(work(lwork))
               u=transpose(vec_sc) !temp
               call dgesvd('O','A',np_on,np_on,u,np_on,w,u,np_on,v,np_on,work,lwork,info )
               v=transpose(v); deallocate(work) 
               if(info.ne.0)then
                    write(11,*)'problem in svd inside vlowknh info=',info
                    call abort
               endif

!  set 1/w(i)=0 if w(i)=0 ("solution of minimum norm" trick) and set recalcvecs=.true. if needed
               do i=1,np_on 
!                   write(11,'(a,f14.8)')'<k|psi_p> sing value',w(i)
                   if(w(i).gt.svd_eps)then
                          v(:,i)=v(:,i)/w(i)
                   else
                          write(11,*)'small singular value in vlowknh' 
                          v(:,i)=v(:,i)*0.0_dp
                   endif                        
               enddo
               biorth(1:np_on,1:np_on)=matmul(v,transpose(u))
               biorth=transpose(biorth) !temp
! test biorth
!               u=0_dp
!               u=matmul(vec_sc,biorth)
!               write(88,*)'******* inside vlowknh *******'
!               do i=1,np_on
!                   do j=1,np_on
!                     if(dabs(u(i,j)).gt.0.0001_dp) write(88,'(2i4,f14.6)')i,j,u(i,j)
!                   enddo
!               enddo

!  evaluate u(k,p)= Sum_k"{<k|veff(p^2)|k"><k"|Psi_p>}
               u=0.0_dp 
               do i=1,np_on                          
                   do j=1,np_on 
                       u(i,j)=dot_product(veff(i,:,j),vec_sc(:,j))
                   enddo
               enddo


!   evaluate non-herm vlowk(k,k')=Sum_p{u(k,p)<biorth_p|k'>}
               vnh=matmul(u,biorth)  

!  remove the weights and factors of k_i*k_j
                do i=1,np_on
                 ii=i
                   if(i.gt.nmod)ii=i-nmod  
                     do j=1,np_on
                        jj=j
                           if(j.gt.nmod)jj=j-nmod
                           fkkw=dsqrt(wk(ii)*wk(jj))*xk(ii)*xk(jj)
                        vnh(i,j)=vnh(i,j)/fkkw        
                     enddo
                enddo
        vnh = vnh*hb2m*2_dp/pi_d 
        deallocate(u,v,w)   
        RETURN
   END SUBROUTINE vlowknh

   SUBROUTINE eigenvectors(lambda,veff,vnn,p_sc,vec_sc)
        USE nrtype
        IMPLICIT NONE
        REAL(DP)::veff(:,:),vnn(:,:),vec_sc(:),dia,fkkw,lambda,p_sc
        REAL(DP),ALLOCATABLE::heff(:,:),eval(:),work(:),temp(:),tempvec(:,:)
        INTEGER(I4B)::lwork,info,i,j,n,ii,jj,iter,indx(100000)        
        INTENT(INOUT)::veff
        INTENT(IN)::vnn,lambda
        INTENT(OUT)::p_sc,vec_sc
!        nit=6
        do iter=1,nit !iteration loop for self-consistent solutions of H_eff(E)         
              lwork=40*np_on
              if(allocated(heff))deallocate(heff)
              if(allocated(eval))deallocate(eval)
              if(allocated(work))deallocate(work)
              allocate(heff(np_on,np_on),eval(np_on),work(lwork))

!  set up heff to be diagonalized (note that heff matrix is replaced by eigenvectors on return from dsyev)
                    do i=1,np_on
                       do j=1,np_on
                            ii=i;jj=j
                            if(i.gt.nmod)ii=i-nmod
                            if(j.gt.nmod)jj=j-nmod
                            fkkw=xk(ii)*xk(jj)*dsqrt(wk(ii)*wk(jj))
                            dia=0.0_dp
                            if(i.eq.j)dia=xk(ii)**2  
                            heff(i,j)=dia+2.0_dp/pi_d*fkkw*veff(i,j)
                       enddo
                    enddo

! diagonalize using LAPACK routine dsyev
                    eval=0.0_dp;work=0.0_dp
                    call dsyev('V','U',np_on,heff,np_on,eval,work,lwork,info)
                          if(info.ne.0)then
                              write(11,'(a,i5)')'problem in dysev in eigenvectors, info=',info
                              call abort
                          endif

!  sort eval in ascending order and store pointer in array indx(:)
                    indx=0
                    call indexx(np_on,eval,indx(1:np_on))
                    allocate(temp(np_on),tempvec(np_on,np_on));temp=eval;tempvec=heff
                         do i=1,np_on
                            eval(i)=temp(indx(i))
                            heff(:,i)=tempvec(:,indx(i))
                         enddo
                    deallocate(temp,tempvec)

! locate where input energy p_on**2 falls in eval(:) values to make the next iteration of p_on**2
	            call hunt(eval,p_on**2,j)
                         n=j+1
                         if(abs(eval(j)-p_on**2).le.abs(eval(j+1)-p_on**2))n=j
                         if(j.eq.0)n=1;if(j.eq.np_on)n=np_on
                         p_on=0.5_dp*(p_on+dsqrt(eval(n)))  !new guess for onshell p_on
                         call blochhorwtz(lambda,vnn,veff) ! recalculate veff at new p_on 
        enddo ! self-consistency loop 

! return the self-consistent p_on and corresponding self-consistent eigenvector
        vec_sc(1:np_on)=heff(1:np_on,n)
        p_sc=p_on

        deallocate(work,eval,heff) 
        RETURN
   END SUBROUTINE eigenvectors


   SUBROUTINE eigenvectors_deut(lambda,veff,vnn,edeut_sc,vec_sc)
! temporary fix to handle deuteron; a more sophisticated fix is to use overloading...
        USE nrtype
        IMPLICIT NONE
        REAL(DP)::veff(:,:),vnn(:,:),vec_sc(:),dia,fkkw,lambda,edeut_sc
        REAL(DP),ALLOCATABLE::heff(:,:),eval(:),work(:)
        INTEGER(I4B)::lwork,info,i,j,n,ii        
        INTENT(INOUT)::veff
        INTENT(IN)::vnn,lambda
        INTENT(OUT)::edeut_sc,vec_sc
!        nit=6
               do ii=1,nit         
                    lwork=40*np_on
                    if(allocated(heff))deallocate(heff)
                    if(allocated(eval))deallocate(eval)
                    if(allocated(work))deallocate(work)
                    allocate(heff(2*nmod,2*nmod),eval(2*nmod),work(lwork))
!  set up heff to be diagonalized (note that heff matrix is replaced by eigenvectors on return from dsyev)
                    do i=1,nmod
                       do j=1,nmod
                            fkkw=xk(i)*xk(j)*dsqrt(wk(i)*wk(j))
                            dia=0.0_dp
                            if(i.eq.j)dia=xk(i)**2  
                            heff(i,j)=dia+2.0_dp/pi_d*fkkw*veff(i,j)
                            heff(i+nmod,j+nmod)=dia+2.0_dp/pi_d*fkkw*veff(i+nmod,j+nmod)
                            heff(i,j+nmod)=2.0_dp/pi_d*fkkw*veff(i,j+nmod)
                            heff(i+nmod,j)=2.0_dp/pi_d*fkkw*veff(i+nmod,j)
                       enddo
                    enddo

                     eval=0.0_dp;work=0.0_dp
                     call dsyev('V','U',2*nmod,heff,2*nmod,eval,work,lwork,info)
                          if(info.ne.0)then
                              write(11,'(a,i5)')'problem in dysev in eigenvectors_deut, info=',info
                              call abort
                          endif
!                    write(11,'(a,i5,f14.5)')('nit',ii,hb2m*eval(i),i=1,2*nmod)
!                    write(11,*)' '
	            call hunt(eval,edeut,j)
                         n=j+1
                         if(abs(eval(j)-edeut).le.abs(eval(j+1)-edeut))n=j
                         if(ii.eq.1)write(11,'(2f15.6)')hb2m*edeut,hb2m*eval(n)
                         if(ii.eq.nit)write(11,'(2f15.6)')hb2m*edeut,hb2m*eval(n)

                         edeut=0.5_dp*(edeut+(eval(n)))  !new guess for onshell p_on
                         veff=0_dp
                         call blochhorwtz_deut(lambda,vnn,veff) ! recalculate veff at new p_on 
               enddo

! return the self-consistent p_on and corresponding self-consistent eigenvector
               vec_sc(1:2*nmod)=heff(1:2*nmod,n)
               edeut_sc=edeut

        deallocate(work,eval,heff) 
        RETURN
   END SUBROUTINE eigenvectors_deut
 

   SUBROUTINE blochhorwtz_deut(lambda,vnn,veffwformfactor)
! a stupid (temporary) fix to handle the boundstate problem; eventually i should use overloading  
        IMPLICIT NONE
        REAL(DP)::vnn(:,:),veffwformfactor(:,:),ans,zero,lambda
        REAL(DP),ALLOCATABLE :: kern(:,:),veff(:,:)
        INTEGER(I4B)::i,j,info,ipiv(600)
        INTENT(IN)::vnn,lambda
        INTENT(OUT)::veffwformfactor

        allocate(kern(2*nkpt,2*nkpt),veff(2*nkpt,2*nkpt));kern=0.0_dp;veff=0_dp
! set up kernel=1-gV to solve veff=(1-gV)^-1*V  (no principal value subtraction here)          
               do i=1,nkpt
                   do j=1,nkpt
                       kern(i,j)=-2.0_dp/pi_d*vnn(i,j)*xk(j)**2*wk(j)/(edeut-xk(j)**2)*smoothq(xk(j),reg)  !ss
                       if(i.eq.j)kern(i,j)=kern(i,j)+1.0_dp

                       kern(i+nkpt,j+nkpt)=-2.0_dp/pi_d*vnn(i+nkpt,j+nkpt)*xk(j)**2*wk(j)/(edeut-xk(j)**2)* &  !dd
                       smoothq(xk(j),reg)
                       if(i.eq.j)kern(i+nkpt,j+nkpt)=kern(i+nkpt,j+nkpt)+1.0_dp

                       kern(i,j+nkpt)=-2.0_dp/pi_d*vnn(i,j+nkpt)*xk(j)**2*wk(j)/(edeut-xk(j)**2)* &  !sd
                       smoothq(xk(j),reg)

                       kern(i+nkpt,j)=-2.0_dp/pi_d*vnn(i+nkpt,j)*xk(j)**2*wk(j)/(edeut-xk(j)**2)* &  !ds
                       smoothq(xk(j),reg)
                   enddo
               enddo

! solve for veff using LAPACK routine dgesv
               veff(1:2*nkpt,1:2*nkpt)=vnn(1:2*nkpt,1:2*nkpt)
               call dgesv(2*nkpt,2*nkpt,kern,2*nkpt,ipiv(1:2*nkpt),veff,2*nkpt,info)
               if(info.ne.0)then
                  write(6,'(a,i5)')'problem in dgesv called in blochhorwtz, info=',info
                  call abort
               endif

! attach the external smooth form factors
               do i=1,nkpt
                   do j=1,nkpt
                       veff(i,j)=veff(i,j)*smoothp(xk(i),reg)*smoothp(xk(j),reg)
                       veff(i+nkpt,j+nkpt)=veff(i+nkpt,j+nkpt)*smoothp(xk(i),reg)*smoothp(xk(j),reg)
                       veff(i,j+nkpt)=veff(i,j+nkpt)*smoothp(xk(i),reg)*smoothp(xk(j),reg)
                       veff(i+nkpt,j)=veff(i+nkpt,j)*smoothp(xk(i),reg)*smoothp(xk(j),reg)
                   enddo
               enddo

               veffwformfactor(1:nmod,1:nmod)=veff(1:nmod,1:nmod)
               veffwformfactor(nmod+1:2*nmod,nmod+1:2*nmod)=veff(nkpt+1:nkpt+nmod,nkpt+1:nkpt+nmod)
               veffwformfactor(1:nmod,nmod+1:2*nmod)=veff(1:nmod,nkpt+1:nkpt+nmod)
               veffwformfactor(nmod+1:2*nmod,1:nmod)=veff(nkpt+1:nkpt+nmod,1:nmod)
               deallocate(kern,veff)
               RETURN
    END SUBROUTINE blochhorwtz_deut

    SUBROUTINE blochhorwtz(lambda,vnn,veffwformfactor)
        IMPLICIT NONE
        REAL(DP)::vnn(:,:),veffwformfactor(:,:),ans,zero,lambda
        REAL(DP),ALLOCATABLE :: kern(:,:),vbig(:,:),veff(:,:),xxk(:)
        INTEGER(I4B)::i,j,info,ipiv(600)
        INTENT(IN)::vnn,lambda
        INTENT(OUT)::veffwformfactor
       
! evaluate "zero" to be subtracted from the principal value to smooth the singularity   
        zero=0.0_dp
        do i=1,nkpt
            zero=zero+wk(i)/(p_on**2-xk(i)**2)
        enddo
            zero=zero-.5_dp/p_on*dlog((kmax+p_on)/(kmax-p_on))
        veffwformfactor=0.0_dp
        if(.not.coupled)then !uncoupled channel case
              allocate(kern(nkpt+1,nkpt+1),vbig(nkpt+1,nkpt+1),veff(nkpt+1,nkpt+1),xxk(nkpt+1));vbig=0.0_dp;kern=0.0_dp;veff=0.0_dp

! interpolate vnn onto a bigger mesh that includes the pole for the principal value subtraction  
               xxk(1:nkpt)=xk(1:nkpt);xxk(nkpt+1)=p_on
	       call vkku399lag(nkpt,vnn,nkpt+1,xxk,vbig)

! set up kernel=1-gV to solve veff=(1-gV)^-1*V            
               do i=1,nkpt+1
                   do j=1,nkpt
                       kern(i,j)=-2/pi_d*vbig(i,j)*xk(j)**2*wk(j)/(p_on**2-xk(j)**2)*smoothq(xk(j),reg)
                       if(i.eq.j)kern(i,j)=kern(i,j)+1.0_dp
                   enddo
                   kern(i,nkpt+1)=2/pi_d*vbig(i,nkpt+1)*p_on**2*zero*smoothq(p_on,reg)
                   if(i.eq.nkpt+1)kern(nkpt+1,nkpt+1)=kern(nkpt+1,nkpt+1)+1.0_dp
               enddo

! solve for veff using LAPACK routine dgesv
               call dgesv(nkpt+1,nkpt+1,kern,nkpt+1,ipiv(1:nkpt+1),vbig,nkpt+1,info)
               if(info.ne.0)then
                  write(6,'(a,i5)')'problem in dgesv called in blochhorwtz, info=',info
                  call abort
               endif
               veff(1:nkpt,1:nkpt)=vbig(1:nkpt,1:nkpt)

! attach the external smooth form factors
               do i=1,nkpt
                   do j=1,nkpt
                       veff(i,j)=veff(i,j)*smoothp(xk(i),reg)*smoothp(xk(j),reg)
                   enddo
               enddo
               veffwformfactor(1:nmod,1:nmod)=veff(1:nmod,1:nmod)
               deallocate(xxk,vbig,kern,veff)      
        else  !coupled channel case
               allocate(kern(2*(nkpt+1),2*(nkpt+1)),vbig(2*(nkpt+1),2*(nkpt+1)),veff(2*nkpt+2,2*nkpt+2),xxk(1:nkpt+1))
               vbig=0.0_dp;kern=0.0_dp;veff=0.0_dp 
               xxk(1:nkpt)=xk(1:nkpt);xxk(nkpt+1)=p_on

! interpolate vnn (for each ll' sub-block) onto bigger mesh that includes the pole for the principal value subtraction  

	       call vkku399lag(nkpt,vnn(1:nkpt,1:nkpt),nkpt+1,xxk,vbig(1:nkpt+1,1:nkpt+1))
	       call vkku399lag(nkpt,vnn(nkpt+1:2*nkpt,nkpt+1:2*nkpt),nkpt+1,xxk,vbig(nkpt+2:2*nkpt+2,nkpt+2:2*nkpt+2))
	       call vkku399lag(nkpt,vnn(1:nkpt,nkpt+1:2*nkpt),nkpt+1,xxk,vbig(1:nkpt+1,nkpt+2:2*nkpt+2))
	       call vkku399lag(nkpt,vnn(nkpt+1:2*nkpt,1:nkpt),nkpt+1,xxk,vbig(nkpt+2:2*nkpt+2,1:nkpt+1))


! set up kernel=1-gV to solve veff=(1-gV)^-1*V            
               do i=1,nkpt+1
                   do j=1,nkpt
                       kern(i,j)=-2.0_dp/pi_d*vbig(i,j)*xk(j)**2*wk(j)/(p_on**2-xk(j)**2)*smoothq(xk(j),reg)  !ss
                       if(i.eq.j)kern(i,j)=kern(i,j)+1.0_dp

                       kern(i+nkpt+1,j+nkpt+1)=-2.0_dp/pi_d*vbig(i+nkpt+1,j+nkpt+1)*xk(j)**2*wk(j)/(p_on**2-xk(j)**2)* &  !dd
                       smoothq(xk(j),reg)
                       if(i.eq.j)kern(i+nkpt+1,j+nkpt+1)=kern(i+nkpt+1,j+nkpt+1)+1.0_dp

                       kern(i,j+nkpt+1)=-2.0_dp/pi_d*vbig(i,j+nkpt+1)*xk(j)**2*wk(j)/(p_on**2-xk(j)**2)* &  !sd
                       smoothq(xk(j),reg)

                       kern(i+nkpt+1,j)=-2.0_dp/pi_d*vbig(i+nkpt+1,j)*xk(j)**2*wk(j)/(p_on**2-xk(j)**2)* &  !ds
                       smoothq(xk(j),reg)
                   enddo
                   kern(i,nkpt+1)=2.0_dp/pi_d*vbig(i,nkpt+1)*p_on**2*zero*smoothq(p_on,reg)
                   if(i.eq.nkpt+1)kern(nkpt+1,nkpt+1)=kern(nkpt+1,nkpt+1)+1.0_dp

                   kern(i+nkpt+1,2*nkpt+2)=2.0_dp/pi_d*vbig(i+nkpt+1,2*nkpt+2)*p_on**2*zero*smoothq(p_on,reg)
                   if(i.eq.nkpt+1)kern(2*nkpt+2,2*nkpt+2)=kern(2*nkpt+2,2*nkpt+2)+1.0_dp

                   kern(i,2*nkpt+2)=2.0_dp/pi_d*vbig(i,2*nkpt+2)*p_on**2*zero*smoothq(p_on,reg)
                   kern(i+nkpt+1,nkpt+1)=2.0_dp/pi_d*vbig(i+nkpt+1,nkpt+1)*p_on**2*zero*smoothq(p_on,reg)
               enddo

! solve for veff using LAPACK routine dgesv

               call dgesv(2*(nkpt+1),2*(nkpt+1),kern,2*(nkpt+1),ipiv(1:2*nkpt+2),vbig,2*(nkpt+1),info)
               if(info.ne.0)then
                  write(6,'(a,i5)')'problem in dgesv called in blochhorwtz, info=',info
                  call abort
               endif
               veff(1:nkpt,1:nkpt)=vbig(1:nkpt,1:nkpt)
               veff(nkpt+1:2*nkpt,nkpt+1:2*nkpt)=vbig(nkpt+2:2*nkpt+1,nkpt+2:2*nkpt+1)
               veff(1:nkpt,nkpt+1:2*nkpt)=vbig(1:nkpt,nkpt+2:2*nkpt+1)
               veff(nkpt+1:2*nkpt,1:nkpt)=vbig(nkpt+2:2*nkpt+1,1:nkpt)

! attach the external smooth form factors
               do i=1,nkpt
                   do j=1,nkpt
                       veff(i,j)=veff(i,j)*smoothp(xk(i),reg)*smoothp(xk(j),reg)
                       veff(i+nkpt,j+nkpt)=veff(i+nkpt,j+nkpt)*smoothp(xk(i),reg)*smoothp(xk(j),reg)
                       veff(i,j+nkpt)=veff(i,j+nkpt)*smoothp(xk(i),reg)*smoothp(xk(j),reg)
                       veff(i+nkpt,j)=veff(i+nkpt,j)*smoothp(xk(i),reg)*smoothp(xk(j),reg)
                   enddo
               enddo
               veffwformfactor(1:nmod,1:nmod)=veff(1:nmod,1:nmod)
               veffwformfactor(nmod+1:2*nmod,nmod+1:2*nmod)=veff(nkpt+1:nkpt+nmod,nkpt+1:nkpt+nmod)
               veffwformfactor(1:nmod,nmod+1:2*nmod)=veff(1:nmod,nkpt+1:nkpt+nmod)
               veffwformfactor(nmod+1:2*nmod,1:nmod)=veff(nkpt+1:nkpt+nmod,1:nmod)

        deallocate(xxk,vbig,kern,veff)      
        endif
   END SUBROUTINE blochhorwtz

	SUBROUTINE hunt(xx,x,jlo)
	USE nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(INOUT) :: jlo
	REAL(DP), INTENT(IN) :: x
	REAL(DP), DIMENSION(:), INTENT(IN) :: xx
	INTEGER(I4B) :: n,inc,jhi,jm
	LOGICAL :: ascnd
	n=size(xx)
	ascnd = (xx(n) >= xx(1))
	if (jlo <= 0 .or. jlo > n) then
		jlo=0
		jhi=n+1
	else
		inc=1
		if (x >= xx(jlo) .eqv. ascnd) then
			do
				jhi=jlo+inc
				if (jhi > n) then
					jhi=n+1
					exit
				else
					if (x < xx(jhi) .eqv. ascnd) exit
					jlo=jhi
					inc=inc+inc
				end if
			end do
		else
			jhi=jlo
			do
				jlo=jhi-inc
				if (jlo < 1) then
					jlo=0
					exit
				else
					if (x >= xx(jlo) .eqv. ascnd) exit
					jhi=jlo
					inc=inc+inc
				end if
			end do
		end if
	end if
	do
		if (jhi-jlo <= 1) then
			if (x == xx(n)) jlo=n-1
			if (x == xx(1)) jlo=1
			exit
		else
			jm=(jhi+jlo)/2
			if (x >= xx(jm) .eqv. ascnd) then
				jlo=jm
			else
				jhi=jm
			end if
		end if
	end do
	END SUBROUTINE hunt
         

      SUBROUTINE indexx(n,arr,indx)
      USE nrtype
      INTEGER(I4B) n,indx(n),M,NSTACK
      REAL(DP):: arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER(I4B) i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL(DP):: a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END SUBROUTINE indexx
END MODULE vlowk_smooth_cutoff 
