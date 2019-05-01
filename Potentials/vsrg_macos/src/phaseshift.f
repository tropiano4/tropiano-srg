! *** EDITED BY RJF on 08/10/06 ***
! Changes include:
!  * switching to DO/ENDDO constructions
!  * uppercase for fortran keywords, lowercase (in general) for variables
!  * explicit declarations of variables (not complete yet)
!  * explicit variables (with more precision) for pi and hb2m
!  * split lines with multiple statements into one line per statement
!  * consistent indentation for loops and if statments
! Still to do:
!  * check that the algorithm still works for smooth vlowk, including regulated
!  * set up global constants for pi
!  * add comments to document variables and subroutine logic
!  * explicit declarations for all variables; eliminate "implicit"
!  * switch to more informative variable/routine names
!  * convert to fortran 90

! *** Input potential V is in MeV-fm^3
!***************************************************************************
!***************************************************************************

      MODULE phaseshift    

!	USE myCONSTANTS	
      !USE MATRIX_ROUTINES 
      USE vlowkPeripherals,  ONLY:MethodType, SmoothP
      PRIVATE
      PUBLIC:: phase_uncoup,phase_coup,phase_uncoup_reg,phase_coup_reg
      CONTAINS

!***************************************************************************

      SUBROUTINE phase_uncoup(nfull,npont,xk,wk,vkkless,xelab,xkmod,thos,delta)

      USE vnn_module, ONLY:hb2m

      implicit real*8(a-h,o-z)
      INTEGER:: nfull,npont,ndim,nplus   
      REAL:: xk(:),wk(:),vkkless(:,:),thos(:),delta
      
      REAL,ALLOCATABLE::vkk2(:,:),fker(:,:),xxk(:),wwk(:),dum(:,:)
      REAL:: pi, topi  ! ADDED RJF 08/10/06
      
      INTENT(IN)::xk,wk,nfull,npont,vkkless,xkmod,xelab
      INTENT(OUT)::thos,delta

      nppt = npont+1 
      ALLOCATE(vkk2(nppt,nppt),fker(nppt,nppt),xxk(nppt),wwk(nppt),dum(nppt,nppt))
      xxk(:npont) = xk(:npont)
      wwk(:npont) = wk(:npont)
      xxk(nppt) = 0.
      wwk(nppt) = 0.

      pi = 3.14159265358979323846d0  
      topi = 2.d0/pi

      xkz = sqrt(xelab/(2.*hb2m))  
      xkpole = xkz
      vkk2 = 0.
	
      ! Set up principal value subtraction       
     	 
      czz = 0.d0
       	 
      DO i = 1,npont 
        cdeno = (xxk(i)**2-xkz**2)
        czz = czz+wwk(i)/cdeno
      ENDDO       
        
      IF (npont.LT.nfull) THEN
        czz = czz + 0.5*log((xkmod+xkpole)/(xkmod-xkpole))/xkpole
      ENDIF  

      nplus = npont+1     ! onshell point k = sqrt(Elab/83)
      xxk(nplus) = xkz     

!     Interpolate vkkless onto the extra onshell meshpt-->vkk2       
      CALL vkku399lag(npont,vkkless,nplus,XXK,vkk2)
       
!      write(*,*)'xkz,vkk2(nplus)',xkz,vkk2(nplus,nplus)

!     Divide vkk2 by topi to conform with T = V+2/Pi VGT, 
!        where -Tan(delta)/k =  T      
      vkk2(:,:) = vkk2(:,:)/topi/hb2m 

!     set up the linear equation fkern*T = v and solve for T

      DO i = 1,nplus
        DO j = 1,npont
          dij = 0.d0
          IF (i.EQ.j) dij = 1.d0
          deno = xxk(j)**2 - xkz**2    
          fker(i,j) = dij+topi*vkk2(i,j)*wwk(j)*xxk(j)**2/deno
        ENDDO
      ENDDO        

    
      DO i = 1,nplus
        diag = 0.d0
        IF (i.EQ.nplus) diag = 1.d0
        fker(i,nplus) = diag-topi*czz*vkk2(i,nplus)*xkz**2
      ENDDO

      
!      CALL SLZBMINV(fker,nplus,nplus)
      dum(:,:) = fker(:,:)
      CALL  INVMAT(dum,nplus,nplus,fker)
!      CALL dLINRG(nplus,fker,nplus,fker,nplus)
      vkk2(:,:) = MATMUL(fker,vkk2)
              
      cunit = 180.d0/acos(-1.0)
!      write(*,*)'xkz,vkk2(nplus)',xkz,vkk2(nplus,nplus)

      delta = cunit*atan(-xkz*vkk2(nplus,nplus))

!      write(6,'(a,2e12.4)')' uncoup elab delta',xelab,delta
!      write(14,'(a,2e12.4)')' uncoup elab delta',xelab,delta
!      write(6,*)'--------------------------------------------------'

      DO ix = 1,npont
!        xkv = xk(ix)

!        CALL interp2d(vkk2,xk,xk,npont,ndim,4,xkv,xkz,tkkk)
        thos(ix) = vkk2(ix,nplus) !tkkk  !7.23.99 t half on shell
      ENDDO
               
      DEALLOCATE(FKER,VKK2,XXK,WWK,dum)  
      RETURN
      END SUBROUTINE phase_uncoup

!***************************************************************************
!  coupled channel phase shifts (bar conventions)

      SUBROUTINE phase_coup(nfull,npont,rk,wrk,vkkless,xelab,xkmod,dela,delb,epsj)
                
      USE vnn_module, ONLY:hb2m        
    	
      implicit real(a-h,o-z)
	 	
      REAL::rk(:),wrk(:),vkkless(:,:)
      REAL:: pi, topi  !, hb2m  
      REAL,ALLOCATABLE::fker(:,:),vkkbig(:,:),xk(:),wk(:),dum(:,:)
      
      INTENT(IN)::nfull,npont,rk,wrk,vkkless,xelab,xkmod
      INTENT(OUT)::dela,delb,epsj

      pi = 3.14159265358979323846d0  
      topi = 2.d0/pi
      !hb2m = 41.47105   !  np only!  

      npont2 = 2*npont
      nplus = npont+1
      nplus2 = 2*nplus
    
      ALLOCATE(fker(nplus2,nplus2),vkkbig(nplus2,nplus2),xk(nplus),wk(nplus),dum(nplus2,nplus2))
      xk(1:npont) = rk(1:npont)
      xk(nplus) = 0.
      wk(1:npont) = wrk(1:npont)
      wk(nplus) = 0.

      xkz = sqrt(xelab/(2.*hb2m))  
      xkpole = xkz
      xk(nplus) = xkz  
    
!  set up principal value subtraction
      czz = 0.
      
      DO i = 1,npont
        czz = czz+wk(i)/(xk(i)**2-xkz**2)
      ENDDO
     
      IF (npont.LT.nfull) THEN
        czz = czz+0.5*log((xkmod+xkpole)/(xkmod-xkpole))/xkpole
      ENDIF
                    
      CALL vkku399lag(npont,vkkless(1:npont,1:npont),nplus,xk,vkkbig(1:nplus,1:nplus))
      CALL vkku399lag(npont,vkkless(1:npont,nplus:npont2),nplus,xk,vkkbig(1:nplus,nplus+1:nplus2))
      CALL vkku399lag(npont,vkkless(nplus:npont2,1:npont),nplus,xk,vkkbig(nplus+1:nplus2,1:nplus))
      CALL vkku399lag(npont,vkkless(nplus:npont2,nplus:npont2),nplus,xk,vkkbig(nplus+1:nplus2,nplus+1:nplus2))

      vkkbig(:,:) = vkkbig(:,:)/topi/hb2m

!   set up kr = v and solve for r
      nplusm = nplus-1
      nplusp = nplus+1
      nplus2m = nplus2-1

      DO i = 1,nplus2  !i = row
   	ii = i
        IF (i.GT.nplus) ii = i-nplus

        DO j = 1,nplusm
	  jj = j;IF (j.GT.nplus)jj = j-nplus
	  dij = 0.;IF (i.EQ.j) dij = 1.0
	  deno = xk(jj)**2-xkz**2    
          fker(i,j) = dij+topi*vkkbig(i,j)*wk(jj)*xk(jj)**2/deno
        ENDDO  

        DO j = nplusp,nplus2m
          jj = j
          IF (j.GT.nplus) jj = j-nplus

          dij = 0.
          IF (i.EQ.j) dij = 1.0

          deno = xk(jj)**2 - xkz**2
          fker(i,j) = dij+topi*vkkbig(i,j)*wk(jj)*xk(jj)**2/deno    
        ENDDO      
      ENDDO   

      DO i = 1,nplus2
        diag1 = 0.0
        diag2 = 0.0
        IF (i.EQ.nplus)diag1 = 1.0
        IF (i.EQ.nplus2)diag2 = 1.0

        fker(i,nplus) = diag1-topi*czz*vkkbig(i,nplus)*xkz**2
        fker(i,nplus2) = diag2-topi*czz*vkkbig(i,nplus2)*xkz**2
      ENDDO

      dum(:,:) = fker(:,:)

!      CALL SLZBMINV(fker,nplus2,nplus2)
      CALL  INVMAT(dum,Nplus2,nplus2,fker)
!      CALL dlinrg(nplus2,fker,nplus2,fker,nplus2)
      vkkbig(:,:) = MATMUL(fker,vkkbig)
                          
      cunit = 180.0/pi 
      r11 = vkkbig(nplus,nplus)  
      r22 = vkkbig(nplus2,nplus2)
      r12 = vkkbig(nplus,nplus2)

      eps = 0.5*atan(2.0*r12/(r11-r22))
      reps = (r11-r22)/cos(2.0*eps)

      dela = -atan(0.5*xkz*(r11+r22+reps))
      delb = -atan(0.5*xkz*(r11+r22-reps))
      epsbar = .5*asin(sin(2.*eps)*sin(dela-delb))
      delbar1 = .5*(dela+delb+asin(tan(2*epsbar)/tan(2*eps)))
      delbar2 = .5*(dela+delb-asin(tan(2*epsbar)/tan(2*eps)))
      dela = delbar1*cunit
      delb = delbar2*cunit
      epsj = epsbar*cunit 
                                     
      DEALLOCATE(fker,xk,wk,vkkbig,dum)

      RETURN
      END SUBROUTINE 

!***************************************************************************
!***************************************************************************

      SUBROUTINE phase_uncoup_reg(cut,method,nfull,npont,xk,wk,vkkless,xelab,xkmod,thos,delta)

      use vnn_module, ONLY:hb2m        
 !        use smoothvlowk, ONLY:smoothp

      implicit real*8(a-h,o-z)

      TYPE(MethodType) :: method
      INTEGER:: nfull,npont,ndim,nplus,ipiv(1000)   
      REAL:: xk(:),wk(:),vkkless(:,:),thos(:)

      REAL,ALLOCATABLE::vkk2(:,:),fker(:,:),xxk(:),wwk(:),dum(:,:)

      INTENT(IN)::xk,wk,nfull,npont,vkkless,xkmod,xelab,method
      INTENT(OUT)::thos,delta
 !        external smoothp	

      nppt=npont+1 
      ALLOCATE(vkk2(nppt,nppt),fker(nppt,nppt),xxk(nppt),wwk(nppt),dum(nppt,nppt))
      xxk(:npont)=xk(:npont)
      wwk(:npont)=wk(:npont)
      xxk(nppt)=0.
      wwk(nppt)=0.

      topi=2.d0/acos(-1.d0)   !pi

      xkz=sqrt(xelab/(2.d0*hb2m))
      xkpole=xkz
      vkk2=0.

      !  SET UP PRINCIPAL VALUE subtraction       

      czz=0.d0
       	 
      DO i = 1,npont 
        cdeno = (xxk(i)**2-xkz**2)
        czz = czz+wwk(i)/cdeno
      ENDDO       
        
      IF (npont.LT.nfull) THEN
        czz = czz + 0.5*log((xkmod+xkpole)/(xkmod-xkpole))/xkpole
      ENDIF  

      nplus = npont+1     ! onshell point k = sqrt(Elab/83)
      xxk(nplus) = xkz     

!     Interpolate vkkless onto the extra onshell meshpt-->vkk2       
      CALL vkku399lag(npont,vkkless,nplus,XXK,vkk2)
       
!      write(*,*)'xkz,vkk2(nplus)',xkz,vkk2(nplus,nplus)

!     Divide vkk2 by topi to conform with T = V+2/Pi VGT, 
!        where -Tan(delta)/k =  T      
      vkk2(:,:) = vkk2(:,:)/topi/hb2m 

!     set up the linear equation fkern*T = v and solve for T

      DO i = 1,nplus
        DO j = 1,npont
          dij = 0.d0
          IF (i.EQ.j) dij = 1.d0
          deno = xxk(j)**2 - xkz**2    
          fker(i,j) = dij+topi*vkk2(i,j)*wwk(j)*xxk(j)**2/deno
        ENDDO
      ENDDO        

    
      DO i = 1,nplus
        diag = 0.d0
        IF (i.EQ.nplus) diag = 1.d0
        fker(i,nplus) = diag-topi*czz*vkk2(i,nplus)*xkz**2
      ENDDO

      
!      CALL SLZBMINV(fker,nplus,nplus)
      dum(:,:) = fker(:,:)
      CALL  INVMAT(dum,nplus,nplus,fker)
!      CALL dLINRG(nplus,fker,nplus,fker,nplus)
      vkk2(:,:) = MATMUL(fker,vkk2)

!!s      CALL DGESV(NPLUS,NPLUS,DUM,NPLUS,IPIV,VKK2,NPLUS,INFO)
!!s             if(info.ne.0)then
!!s                 write(6,*)'problem in dgesv in phaseshift',info 
!!s                 call abort
!!s             endif

      cunit = 180.d0/acos(-1.d0)
!      write(*,*)'xkz,vkk2(nplus)',xkz,vkk2(nplus,nplus)
      delta = cunit*atan(-xkz*vkk2(nplus,nplus)*SmoothP(xkz,method)**2)  !!smoothp(xkz,cut)**2)
!        write(6,'(a,2e12.4)')' uncoup elab delta',xelab,delta
!        write(14,'(a,2e12.4)')' uncoup elab delta',xelab,delta
!        write(6,*)'--------------------------------------------------'

      DO ix = 1,npont
!        xkv = xk(ix)

!        CALL interp2d(vkk2,xk,xk,npont,ndim,4,xkv,xkz,tkkk)
        thos(ix) = vkk2(ix,nplus) !tkkk  !7.23.99 t half on shell
      ENDDO
               
      DEALLOCATE(FKER,VKK2,XXK,WWK,dum)  
      RETURN
      END SUBROUTINE phase_uncoup_reg


!***************************************************************************
!***************************************************************************

      SUBROUTINE phase_coup_reg(cut,method,nfull,npont,rk,wrk,vkkless,xelab,xkmod,dela,delb,epsj)

      use vnn_module, ONLY:hb2m

 !        use smoothvlowk, ONLY:smoothp
      implicit real(a-h,o-z)

      TYPE(MEthodType) :: method
      REAL::rk(:),wrk(:),vkkless(:,:)

      REAL,ALLOCATABLE::fker(:,:),vkkbig(:,:),xk(:),wk(:),dum(:,:)

      INTENT(IN)::nfull,npont,rk,wrk,vkkless,xelab,xkmod,method
      INTENT(OUT)::dela,delb,epsj
 !       external smoothp     	 

      topi = 2.0/acos(-1.0) !3.1415 !PI
      h2m = hb2m 
      
      npont2 = 2*npont
      nplus = npont+1
      nplus2 = 2*nplus
    
      ALLOCATE(fker(nplus2,nplus2),vkkbig(nplus2,nplus2),xk(nplus),wk(nplus),dum(nplus2,nplus2))
      xk(1:npont) = rk(1:npont)
      xk(nplus) = 0.
      wk(1:npont) = wrk(1:npont)
      wk(nplus) = 0.

      xkz = sqrt(xelab/(2.*hb2m))  
      xkpole = xkz
      xk(nplus) = xkz  
    
!  set up principal value subtraction
      czz = 0.
      
      DO i = 1,npont
        czz = czz+wk(i)/(xk(i)**2-xkz**2)
      ENDDO
     
      IF (npont.LT.nfull) THEN
        czz = czz+0.5*log((xkmod+xkpole)/(xkmod-xkpole))/xkpole
      ENDIF
                    
      CALL vkku399lag(npont,vkkless(1:npont,1:npont),nplus,xk,vkkbig(1:nplus,1:nplus))
      CALL vkku399lag(npont,vkkless(1:npont,nplus:npont2),nplus,xk,vkkbig(1:nplus,nplus+1:nplus2))
      CALL vkku399lag(npont,vkkless(nplus:npont2,1:npont),nplus,xk,vkkbig(nplus+1:nplus2,1:nplus))
      CALL vkku399lag(npont,vkkless(nplus:npont2,nplus:npont2),nplus,xk,vkkbig(nplus+1:nplus2,nplus+1:nplus2))

      vkkbig(:,:)=vkkbig(:,:)/topi/hb2m


!   set up kr=v and solve for r
      nplusm = nplus-1
      nplusp = nplus+1
      nplus2m = nplus2-1

      DO i = 1,nplus2  !i = row
   	ii = i
        IF (i.GT.nplus) ii = i-nplus

        DO j = 1,nplusm
	  jj = j;IF (j.GT.nplus)jj = j-nplus
	  dij = 0.;IF (i.EQ.j) dij = 1.0
	  deno = xk(jj)**2-xkz**2    
          fker(i,j) = dij+topi*vkkbig(i,j)*wk(jj)*xk(jj)**2/deno
        ENDDO  

        DO j = nplusp,nplus2m
          jj = j
          IF (j.GT.nplus) jj = j-nplus

          dij = 0.
          IF (i.EQ.j) dij = 1.0

          deno = xk(jj)**2 - xkz**2
          fker(i,j) = dij+topi*vkkbig(i,j)*wk(jj)*xk(jj)**2/deno    
        ENDDO      
      ENDDO   

      DO i = 1,nplus2
        diag1 = 0.0
        diag2 = 0.0
        IF (i.EQ.nplus)diag1 = 1.0
        IF (i.EQ.nplus2)diag2 = 1.0

        fker(i,nplus) = diag1-topi*czz*vkkbig(i,nplus)*xkz**2
        fker(i,nplus2) = diag2-topi*czz*vkkbig(i,nplus2)*xkz**2
      ENDDO

      dum(:,:) = fker(:,:)

!     CALL SLZBMINV(fker,nplus2,nplus2)
      CALL  INVMAT(dum,Nplus2,nplus2,fker)
!     	CALL dlinrg(nplus2,fker,nplus2,fker,nplus2)
      vkkbig(:,:)=MATMUL(fker,vkkbig)
                          
      cunit = 180.0/acos(-1.0) !3.1415 !PI
      r11 = vkkbig(nplus,nplus)*SmoothP(xkz,method)**2 !smoothp(xkz,cut)**2  
      r22 = vkkbig(nplus2,nplus2)*SmoothP(xkz,method)**2 !smoothp(xkz,cut)**2  
      r12 = vkkbig(nplus,nplus2)*SmoothP(xkz,method)**2 !smoothp(xkz,cut)**2  

      eps = 0.5*atan(2.0*r12/(r11-r22))
      reps = (r11-r22)/cos(2.0*eps)

      dela = -atan(0.5*xkz*(r11+r22+reps))
      delb = -atan(0.5*xkz*(r11+r22-reps))
      epsbar = .5*asin(sin(2.*eps)*sin(dela-delb))
      delbar1 = .5*(dela+delb+asin(tan(2*epsbar)/tan(2*eps)))
      delbar2 = .5*(dela+delb-asin(tan(2*epsbar)/tan(2*eps)))
      dela = delbar1*cunit
      delb = delbar2*cunit
      epsj = epsbar*cunit 
                                     
      DEALLOCATE(fker,xk,wk,vkkbig,dum)
      
      RETURN
      END SUBROUTINE 


!***************************************************************************
!***************************************************************************

      SUBROUTINE interp2d(z,x,y,n,ip,x0,y0,ans)
!      z = z(x,y) = <x\z\y>
!      ip = no of points used for interpolation

      implicit real(a-h,o-z)
      REAL:: z(:,:),y(:),x(:)
      REAL,ALLOCATABLE::zrow(:),col(:)

      ALLOCATE(zrow(SIZE(Z)),col(SIZE(Z)));zrow = 0.;col = 0.

      DO jcol = 1,n
        DO i = 1,n
          col(i) = z(i,jcol)
        ENDDO

        CALL interp1d(col,x,n,ip,x0,zrow(jcol))
      ENDDO

      CALL interp1d(zrow,y,n,ip,y0,ans)

      DEALLOCATE(ZROW,COL)      

      RETURN
      END SUBROUTINE


!***************************************************************************
!***************************************************************************

      SUBROUTINE interp1d(y,x,n,ip,x0,y0)

      implicit real (a-h,o-z)
!      y = y(x)
!      ip = no of points used for interpolation
      REAL:: y(:),x(:)
      logical offpt

      offpt = .true.
      y00 = 0.d0
      y0 = 0.d0
      i0 = 0

      DO i = 1,n
         xi = x(i)
         xi3 = x(i0+ip-1)
         IF (xi.EQ.x0) THEN
           y0 = y(i)
           offpt = .false.
         ELSE IF (i0.LT.n-ip) THEN
           IF (x0.GT.xi3) THEN
             i0 = i0 + 1
           ENDIF
         ENDIF
      ENDDO

      IF (offpt) THEN
        DO  k = 1,ip
          f = 1.d0
          DO l = 1,k-1
            f = f*(x0-x(i0+l))/(x(i0+k)-x(i0+l))
          ENDDO 

          IF (k.NE.ip) THEN
            DO l = k+1,ip
              f = f*(x0-x(i0+l))/(x(i0+k)-x(i0+l))
            ENDDO 
          ENDIF

          y00 = y00 + f*y(i0+k)
        ENDDO    
        y0 = y00
      ENDIF

      RETURN 
      END SUBROUTINE


!***************************************************************************
!***************************************************************************

      SUBROUTINE vkku399lag(nless,vkkless,nr,ra,vkk2)

      implicit real(a-h,o-z)
      REAL::vkkless(:,:),vkk2(:,:),ra(:)
      INTENT(IN)::nless,vkkless,ra,nr
      INTENT(OUT)::vkk2     

!     3/20/99 last point by interpolation
      
      DO j = 1,nr
        DO i = 1,nr

          IF ( (i.NE.nr).and.(j.NE.nr) ) THEN
            vkk2(i,j) = vkkless(i,j)
          ENDIF

          IF ( (i.EQ.nr).OR.(j.EQ.nr) ) THEN
            xi = ra(i)
            xj = ra(j)
            CALL interp2d (vkk2,ra,ra,nless,4,xi,xj,ans)
            vkk2(i,j) = ans
          ENDIF
          
        ENDDO
      ENDDO      

      RETURN
      END SUBROUTINE

!***************************************************************************
!***************************************************************************

      SUBROUTINE INVMAT(A,N,lda,Ainv)
      
      DIMENSION A(lda,lda),Ainv(lda,lda),ipiv(N)
      
      Ainv = 0
      DO i = 1,n
        Ainv(i,i) = 1.
      ENDDO

      CALL DGESV(N,N,A,lda,ipiv,Ainv,lda,info)

      IF (info.NE.0) PAUSE 'error in invmat'

      RETURN
      END SUBROUTINE

!***************************************************************************
!***************************************************************************
      
      END module 
