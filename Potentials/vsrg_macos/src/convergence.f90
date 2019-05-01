 module convergence
      use myinterpolation, only: vkku399lag
      implicit none
      private
      logical :: inmedium
      real(8), parameter :: hb2m = 41.47
      complex(16), allocatable :: psivec(:,:), etas(:) 
      public::weinberg !,diagonalizeK,weinberg_u,weinberg_c

 contains


 subroutine weinberg(coupled,P,xkf,lambda,Ecm,xk,wk,vnn,np,eta,psi)
      logical:: coupled
      integer:: np
      real(8):: P, xkf, lambda, Ecm, xk(:), wk(:), vnn(:,:)
      complex(16) :: eta(size(vnn(1,:))), psi(size(vnn(1,:)),size(vnn(1,:)))
      intent(in) :: coupled, P, xkf, lambda, Ecm, xk, wk, vnn, np
      intent(out) :: eta, psi

      inmedium=.FALSE.
      IF(xkf .gt. 0.0001d0) inmedium = .TRUE.
       
      IF(coupled)ALLOCATE(psivec(1:2*np,1:2*np),etas(1:2*np))
      IF(.not.coupled)ALLOCATE(psivec(1:np,1:np),etas(1:np))

      IF(Ecm .gt. 0)THEN
           
          IF(coupled)THEN
              
              CALL weinberg_c(P,xkf,lambda, Ecm, xk, wk, vnn,np)
          ELSE 
              CALL weinberg_u(P,xkf,lambda, Ecm, xk, wk, vnn(1:np,1:np),np)
          ENDIF
  
      ELSE

        IF(coupled)THEN  
           CALL diagonalizeK(xkf,P,lambda,coupled,Ecm,1,xk,wk,vnn(1:2*np,1:2*np),np)
        ELSE
          CALL diagonalizeK(xkf,P,lambda,coupled,Ecm,1,xk,wk,vnn(1:np,1:np),np)
        ENDIF
      ENDIF

      eta = etas
      psi = psivec

      DEALLOCATE(psivec,etas)

 end subroutine weinberg


 subroutine weinberg_u(P,xkf,lambda,Elab,xk,wk,vnn,np)
!      use vnn_module, ONLY:hb2m
      real(8)::zero,denom,Elab,pole,vnn(:,:),xk(:),wk(:),pi,kkww,lambda,xE,P,xkf
      real(8)::zerokf
      real(8),allocatable::rwork(:),vbig(:,:),xxk(:),wwk(:)
      complex*16,allocatable::kernel(:,:),rvec(:,:),lvec(:,:),evalr(:),work(:),repulse(:),attract(:)
      integer::izero,nkf,ir,ia,npoint,info,lwork,i,j,ie,np,ii,jj
      character*1 jobvl,jobvr
     
      nkf=0
      do i = 1,np
        if(xk(i).le. xkf)nkf=nkf+1
      enddo
!      write(8,*)'nkf in weinberg',nkf

      pi=2.*asin(1.0)
      npoint=np+1

      lwork=10*npoint
!      rwork=0. 
      pole=sqrt(Elab/hb2m) !Elab is w/com subtracted out
      write(6,*)'pole',pole
           allocate(kernel(npoint,npoint),rvec(npoint,npoint),lvec(npoint,npoint),&
           rwork(2*npoint),evalr(npoint),work(lwork),xxk(npoint),wwk(npoint),vbig(npoint,npoint),&
           attract(npoint),repulse(npoint))
           kernel=0.;rvec=0.;lvec=0.;evalr=0.;work=0.;attract=0.;repulse=0.
 
      xxk(1:np)=xk(1:np)
      wwk(1:np)=wk(1:np)
      xxk(npoint)=pole
!     set up principal value subtraction
      zero=0.
      do j=1,np
         denom=(pole**2-xxk(j)**2)
          zero=zero+wwk(j)/denom
      enddo

      zerokf=0.
      do j=1,nkf
         denom=(pole**2-xxk(j)**2)
          zerokf=zerokf+wwk(j)/denom
      enddo

      if(xkf.gt.0.01.and.pole.gt.xkf)then      
      zero=zero-0.5*dlog((lambda+pole)/(lambda-pole))/pole & 
              -0.5*dlog(abs(pole-xkf)/(pole+xkf))/pole*2. 
      zerokf = zerokf -.5*dlog((xkf+pole)/(xkf-pole))/pole
      else 
      zero=zero-0.5*dlog((lambda+pole)/(lambda-pole))/pole  
      endif
      call vkku399lag(np,vnn(1:np,1:np),npoint,xxk,vbig(1:npoint,1:npoint))

!      write(8,'(3f14.6)')(xxk(i),vnn(i,i),vbig(i,i),i=1,npoint)   
      jobvl='N'
      jobvr='V'
           do i=1,npoint
              do j=1,np !npoint
!               kernel(i,j)=cmplx(2./pi*pauli(xkf,P,xxk(j))*vbig(i,j)*xxk(j)**2*wwk(j)/(pole**2-xxk(j)**2),0.d0)
               kernel(i,j)=cmplx(2./pi*pauli(xkf,P,xxk(j))*vbig(i,j)*xxk(j)**2*wwk(j)/(pole**2-xxk(j)**2),0.d0)
              enddo
           enddo
       
          do i=1,npoint
          kernel(i,npoint)=cmplx(-2./pi*pauli(xkf,P,pole)*vbig(i,npoint)* &
                        pole**2*zero,-pole*vbig(i,npoint))

          enddo
!      do ie=1,nen         
      info=1000
      call zgeev(jobvl,jobvr,npoint,kernel(:,:),npoint,evalr(:) &
               ,rvec(:,:),npoint,lvec(:,:),npoint,work,lwork,rwork,info)
      if(info.ne.0)then
           write(6,*)'problem in dgeev info=',info
           call abort
      endif
!      call separablecomplex(vnn,xk,wk,npoint,rvec,lvec,evalr)
          ia=1;ir=1   
      do i=1,npoint
          if(dble(evalr(i)).gt.0.001)then
             attract(ia)=evalr(i)
             ia=ia+1
          else if(dble(evalr(i)).lt.-0.001)then 
             repulse(ir)=evalr(i)
             ir=ir+1
          endif
      enddo
      ia=ia-1
      ir=ir-1     

      izero=1
      if(xkf.gt.0.01)then
      do i=1,npoint
         if(xnorm(evalr(i)).ge.0.001)goto 10
         izero=izero+1
      enddo
10    continue   
      endif

      psivec = rvec
      etas = evalr

!       write(8,'(f9.4,1f14.5)')(elab,xnorm(evalr(i)),i=izero,izero+1)

!      write(8,'(f8.3,8f12.6)')elab,evalr(izero),evalr(izero+1)
!     1 ,evalr(izero+2),evalr(izero+3) !,xnorm(evalr(izero+2))
!        write(8,'(f8.3,4f12.6)')elab,xnorm(evalr(izero)),xnorm(evalr(izero+1)),&
!                               xnorm(evalr(izero+2)),xnorm(evalr(izero+3))
!      write(8,'(2f14.4)')(evalr(i),i=1,4)
!      if(dble(evalr(izero)).lt.0)write(8,'(f8.3,7f11.4)')elab,evalr(izero),evalr(izero+1)
!     1 ,xnorm(evalr(izero)),xnorm(evalr(izero+1)),xnorm(evalr(izero+2))
!      if(dble(evalr(izero)).gt.0)write(8,'(f8.3,7f11.4)')elab,evalr(izero+1),evalr(izero)
!     1 ,xnorm(evalr(izero+1)),xnorm(evalr(izero)),xnorm(evalr(izero+2))
!      write(8,'(2f10.3,2f14.6)')lambda,Elab,
!     1 xnorm(evalr(1)),xnorm(evalr(2))
!      write(8,*)' ' 
!      write(8,'(2f10.3,3f14.6)')(lambda,Elab,attract(i),
!     1 xnorm(attract(i)), i=1,min(6,ia))
!      write(8,*)' '    
!      write(8,'(4f10.3,3f14.6)')(lambda,xkf,P,Elab,evalr(i),
!     1 xnorm(evalr(i)), i=1,10)
!      write(8,'(2f10.3,6f11.4)')lambda,Elab,xnorm(evalr(1)),
!     1 xnorm(evalr(2)),xnorm(evalr(3)),xnorm(evalr(4)),xnorm(evalr(5)),xnorm(evalr(6))
      deallocate(rwork,attract,repulse,kernel,evalr,rvec,lvec,work,xxk,wwk,vbig)
      return
     end subroutine 



      subroutine weinberg_c(P,xkf,lambda,Elab,xk,wk,vnn,np)
!      use vnn_module, ONLY:hb2m
      real(8)::zero,denom,Elab,pole,vnn(:,:),xk(:),wk(:),pi,kkww,lambda,xE,P,xkf
      real(8),allocatable::rwork(:),vbig(:,:),xxk(:),wwk(:)
      complex*16,allocatable::kernel(:,:),rvec(:,:),lvec(:,:),evalr(:),work(:),repulse(:),attract(:)
      integer::npp1,npp2,ir,ia,npoint,info,lwork,i,j,ie,np,ii,jj,izero
      character*1 jobvl,jobvr

      pi=2.*asin(1.0)
      npoint=2*(np+1)
      npp1=np+1
      npp2=np+2
 
      lwork=10*npoint
!      rwork=0. 
      pole=sqrt(Elab/hb2m) !Elab is w/com subtracted out
      write(6,*)'pole',pole
           allocate(kernel(npoint,npoint),rvec(npoint,npoint),lvec(npoint,npoint), &
           rwork(2*npoint),evalr(npoint),work(lwork),xxk(npoint),wwk(npoint),vbig(npoint,npoint),&
           attract(npoint),repulse(npoint))
           kernel=0.;rvec=0.;lvec=0.;evalr=0.;work=0.;attract=0.;repulse=0.
 
      xxk(1:np)=xk(1:np)
      wwk(1:np)=wk(1:np)
      xxk(npp1)=pole
!     set up principal value subtraction
      zero=0.
      do j=1,np
         denom=(pole**2-xxk(j)**2)
          zero=zero+wwk(j)/denom
      enddo

      if(xkf.gt.0.and.pole.gt.xkf)then      
      zero=zero-0.5*dlog((lambda+pole)/(lambda-pole))/pole  
!     1         +0.5*dlog((pole+xkf)/(pole-xkf))/pole !note this is redundant since "zero" above has it
      else 
      zero=zero-0.5*dlog((lambda+pole)/(lambda-pole))/pole  
      endif

      call vkku399lag(np,vnn(1:np,1:np),npp1,xxk,vbig(1:npp1,1:npp1))  !ss block
      call vkku399lag(np,vnn(npp1:2*np,npp1:2*np),npp1,xxk,vbig(npp2:npoint,npp2:npoint)) !dd block
      call vkku399lag(np,vnn(npp1:2*np,1:np),npp1,xxk,vbig(npp2:npoint,1:npp1))  !ds block
      call vkku399lag(np,vnn(1:np,npp1:2*np),npp1,xxk,vbig(1:npp1,npp2:npoint))  !sd block


!      write(8,'(3f14.6)')(xxk(i),vnn(i,i),vbig(i,i),i=1,npoint)   

! ss block      
           do i=1,npp1 !npoint
              do j=1,np !npoint
               kernel(i,j)=cmplx(2./pi*pauli(xkf,P,xxk(j))*vbig(i,j)*xxk(j)**2*wwk(j)/(pole**2-xxk(j)**2),0.d0)
              enddo
           enddo
       
          do i=1,npp1 !npoint
          kernel(i,npp1)=cmplx(-2./pi*pauli(xkf,P,pole)*vbig(i,npp1)* &
                        pole**2*zero,-pole*pauli(xkf,P,pole)*vbig(i,npp1))
          enddo
!dd block
           do i=1,npp1 !npoint
              do j=1,np !npoint
               kernel(i+npp1,j+npp1)=cmplx(2./pi*pauli(xkf,P,xxk(j))* &
                   vbig(i+npp1,j+npp1)*xxk(j)**2*wwk(j)/(pole**2-xxk(j)**2),0.d0)
             enddo
           enddo
       
          do i=1,npp1 !npoint
          kernel(i+npp1,npoint)=cmplx(-2./pi*pauli(xkf,P,pole)*vbig(i+npp1,npoint)* &
                        pole**2*zero,-pole*pauli(xkf,P,pole)*vbig(i+npp1,npoint))
          enddo

!ds block
           do i=1,npp1 !npoint
              do j=1,np !npoint
               kernel(i+npp1,j)=cmplx(2./pi*pauli(xkf,P,xxk(j))* &
                   vbig(i+npp1,j)*xxk(j)**2*wwk(j)/(pole**2-xxk(j)**2),0.d0)
              enddo
           enddo
       
          do i=1,npp1 !npoint
          kernel(i+npp1,npp1)=cmplx(-2./pi*pauli(xkf,P,pole)*vbig(i+npp1,npp1)*&
                        pole**2*zero,-pole*pauli(xkf,P,pole)*vbig(i+npp1,npp1))
          enddo
!sd block
           do i=1,npp1 !npoint
              do j=1,np !npoint
               kernel(i,j+npp1)=cmplx(2./pi*pauli(xkf,P,xxk(j))* &
                   vbig(i,j+npp1)*xxk(j)**2*wwk(j)/(pole**2-xxk(j)**2),0.d0)
              enddo
           enddo
       
          do i=1,npp1 !npoint
          kernel(i,npoint)=cmplx(-2./pi*pauli(xkf,P,pole)*vbig(i,npoint)* &
                        pole**2*zero,-pole*pauli(xkf,P,pole)*vbig(i,npoint))
          enddo



      jobvl='N'
      jobvr='V'
!      do ie=1,nen         
      info=1000
!      WRITE(*,*)'before zgeev'
      call zgeev(jobvl,jobvr,npoint,kernel(:,:),npoint,evalr(:) &
               ,lvec(:,:),npoint,rvec(:,:),npoint,work,lwork,rwork,info)
!      WRITE(*,*)'after zgeev'
      if(info.ne.0)then
           write(6,*)'problem in dgeev info=',info
           call abort
      endif
      
          ia=1;ir=1   
      do i=1,npoint
          if(dble(evalr(i)).gt.0.001)then
             attract(ia)=evalr(i)
             ia=ia+1
          else if(dble(evalr(i)).lt.-0.001)then 
             repulse(ir)=evalr(i)
             ir=ir+1
          endif
      enddo
      ia=ia-1
      ir=ir-1     

!      if(dble(evalr(1)).le.0)write(8,'(2f10.3,2f14.6)')lambda,Elab,
!     1 xnorm(evalr(1)),xnorm(evalr(2))
!      if(dble(evalr(1)).gt.0)write(8,'(2f10.3,2f14.6)')lambda,Elab,
!     1 xnorm(evalr(2)),xnorm(evalr(1))
!      write(8,*)' ' 
!      write(8,'(2f10.3,3f14.6)')(lambda,Elab,attract(i),
!     1 xnorm(attract(i)), i=1,min(6,ia))
!      write(8,*)' '    
      izero=1
      if(xkf.gt.0.01)then
      do i=1,npoint
         if(xnorm(evalr(i)).ge.0.001)goto 10
         izero=izero+1
      enddo
10    continue   
      endif

      psivec = rvec
      etas = evalr
!        write(8,'(f8.3,4f12.6)')elab,xnorm(evalr(izero)),xnorm(evalr(izero+1)),
!     1                          xnorm(evalr(izero+2)),xnorm(evalr(izero+3))
!      write(8,'(2f19.9)')(evalr(i),i=1,npoint)
!      write(8,'(f9.4,6f14.5)')elab,evalr(izero),evalr(izero+1),xnorm(evalr(izero)), &
!                                  xnorm(evalr(izero+1))
!       write(8,'(f9.4,1f14.5)')(elab,xnorm(evalr(i)),i=izero,izero+1)
!      write(8,'(f8.3,7f11.4)')elab,evalr(izero),evalr(izero+1)
!     1 ,xnorm(evalr(izero)),xnorm(evalr(izero+1)),xnorm(evalr(izero+2))
!      if(dble(evalr(izero)).lt.0)write(8,'(f8.3,i3,7f11.4)')elab,izero,evalr(izero),evalr(izero+1)
!     1 ,xnorm(evalr(izero)),xnorm(evalr(izero+1)),xnorm(evalr(izero+2))
!      if(dble(evalr(izero)).gt.0)write(8,'(f8.3,i3,7f11.4)')elab,izero,evalr(izero+1),evalr(izero)
!     1 ,xnorm(evalr(izero+1)),xnorm(evalr(izero)),xnorm(evalr(izero+2))
!      write(8,'(2f10.3,2f14.6)')lambda,Elab,
!      write(8,'(f8.3,6f11.4)')lambda,evalr(izero),evalr(izero+1)
!     1 ,xnorm(evalr(izero)),xnorm(evalr(izero+1))
!      write(8,*)' '
!      write(8,'(4f10.3,3f14.6)')(lambda,xkf,P,Elab,evalr(i),
!     1 xnorm(evalr(i)), i=1,10)
      deallocate(attract,repulse,kernel,evalr,rvec,lvec,work,rwork,xxk,wwk,vbig)
      return
      end subroutine 

      real(8) function xnorm(eval)
      complex*16::eval
      xnorm=dsqrt(dble(eval)**2+aimag(eval)**2)
      return
      end function   

      subroutine diagonalizeK(xkf,P,lambda,coupled,E,nen,xk,wk,vnn,np)
!      use vnn_module, ONLY:hb2m
      logical::coupled
      real(8)::vnn(:,:),xk(:),wk(:),pi,kkww,lambda,xE,xkf,E,P
!      complex*16::E
      real(8),allocatable::rwork(:)
      complex*16,allocatable::kernel(:,:),rvec(:,:),lvec(:,:),evalr(:),work(:)
      integer::izero,nen,npoint,info,lwork,i,j,ie,np,ii,jj
      character*1 jobvl,jobvr
!      write(*,*)'inside diagK'
      pi=2.*asin(1.0)
      npoint=np
      if(coupled)npoint=2*np
      lwork=10*npoint
!      rwork=0. 
      if(.not.coupled)then
           allocate(kernel(npoint,npoint),rvec(npoint,npoint),lvec(npoint,npoint),&
           rwork(2*npoint),evalr(npoint),work(lwork))
           kernel=0.;rvec=0.;lvec=0.;evalr=0.;work=0.

           do i=1,npoint
              do j=1,npoint
                   kkww=sqrt(wk(i)*wk(j))*xk(i)*xk(j)
                   kernel(i,j)=2./pi*vnn(i,j)*kkww*pauli(xkf,P,xk(i))/(E/hb2m-xk(i)**2)
              enddo
           enddo

      jobvl='n'
      jobvr='v'
      info=1000
      call zgeev(jobvl,jobvr,npoint,kernel(:,:),npoint,evalr(:) &
               ,lvec(:,:),npoint,rvec(:,:),npoint,work,lwork,rwork,info)
      if(info.ne.0)then
           write(6,*)'problem in dgeev info=',info
           call abort
      endif
!      call getseparable(vnn,xk,wk,npoint,rvec,evalr)
      else if(coupled)then
      
           allocate(kernel(npoint,npoint),rvec(npoint,npoint),lvec(npoint,npoint),&
           rwork(2*npoint),evalr(npoint),work(lwork))
           kernel=0.;rvec=0.;lvec=0.;evalr=0.;work=0.

           do i=1,npoint
              ii=i
              if(ii.gt.np)ii=i-np
              do j=1,npoint
              jj=j
              if(jj.gt.np)jj=j-np
                   kkww=sqrt(wk(ii)*wk(jj))*xk(ii)*xk(jj)
                   kernel(i,j)=2./pi*vnn(i,j)*pauli(xkf,P,xk(ii))*kkww/(E/hb2m-xk(ii)**2)
              enddo
           enddo

      jobvl='N'
      jobvr='V'
!      do ie=1,nen         
      info=1000
      call zgeev(jobvl,jobvr,npoint,kernel(:,:),npoint,evalr(:) &
               ,lvec(:,:),npoint,rvec(:,:),npoint,work,lwork,rwork,info)
      if(info.ne.0)then
           write(6,*)'problem in dgeev info=',info
           call abort
      endif

      endif
      izero=1
      do i=1,npoint
         if(xnorm(evalr(i)).ge.0.001)goto 10
         izero=izero+1
      enddo
10    continue   
      psivec = rvec
      etas = evalr
!      write(8,'(2f18.9)')(evalr(i),i=1,npoint)
!      write(8,'(1f10.5,8f13.6)')lambda,evalr(1),evalr(2),evalr(3),evalr(4)
!       write(8,'(f9.4,2f14.5)')(lambda,xnorm(evalr(i)),i=izero,izero+1)
!      write(8,'(f8.3,9f11.4)')lambda,evalr(izero),evalr(izero+1),evalr(izero+2)
!     1 ,xnorm(evalr(izero)),xnorm(evalr(izero+1)),xnorm(evalr(izero+2))
!       write(8,'(f9.4,2f14.5)')(lambda,xnorm(evalr(i)),i=izero,izero+1)
!      write(8,'(f8.3,6f12.7)')lambda,real(evalr(izero)),real(evalr(izero+1)),real(evalr(izero+2)), real(evalr(izero+3)),
!     $                               real(evalr(izero+4)), real(evalr(izero+5))
!      if(dble(evalr(izero)).lt.0)write(8,'(f8.3,6f11.4)')lambda,evalr(izero),evalr(izero+1)
!     1 ,xnorm(evalr(izero)),xnorm(evalr(izero+1))
!      if(dble(evalr(izero)).gt.0)write(8,'(f8.3,6f11.4)')lambda,evalr(izero),evalr(izero+1)
!     1 ,xnorm(evalr(izero+1)),xnorm(evalr(izero))
!      write(8,'(f8.3,4f11.4)')lambda,evalr(izero),evalr(izero+1)
!     1 ,xnorm(evalr(izero)),xnorm(evalr(izero+1))
!      if(dble(evalr(1)).lt.0)write(8,'(2f8.3,f10.4,2f14.5)')lambda,xkf,E,dble(evalr(1)),dble(evalr(2))
!      if(dble(evalr(1)).gt.0)write(8,'(2f8.3,f10.4,2f14.5)')lambda,xkf,E,dble(evalr(2)),dble(evalr(1))
!      write(8,'(2f8.3,2f11.4)')lambda,E,xnorm(evalr(1)),xnorm(evalr(2))
!      write(8,'(4f10.3,3f14.6)')(lambda,xkf,P,E,evalr(i),
!     1 xnorm(evalr(i)), i=1,2)
!      write(8,'(4f10.3,3f14.6)')(lambda,xkf,P,E,evalr(i),
!     1 xnorm(evalr(i)), i=1,10)
      deallocate(kernel,rvec,lvec,evalr,work,rwork)
!      write(*,*)'leaving diagK'
      return
      end subroutine 


      real(8) function pauli(xkf,P_cm,xk)
      real(8):: P_cm,xk,xkf, xk1,xk2
      pauli=1.d0
!      return
!Scott      IF(inmedium==.TRUE.)THEN
      IF(inmedium)THEN
      pauli=0.
      if(abs(P_cm) .ge.0.00001)then 
      		if(xk**2.lt.(xkf**2-.25*P_cm**2))return
      		if(xk.ge.(xkf+.5*P_cm))pauli=1.0
      	if((xk**2.gt.(xkf**2-.25*P_cm**2)).and.(xk**2.le.(xkf+.5*P_cm)**2))pauli=1.0/(P_cm*xk)*(.25*P_cm**2+xk**2-xkf**2) 
      else 
      		if(xk.ge.xkf)pauli=1.
      endif
      ENDIF
      return
      end function
!      real function pauli(xkf,xk,P)
end module
