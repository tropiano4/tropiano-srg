MODULE hermitization 
 USE nrtype
 USE VlowkPeripherals

 IMPLICIT NONE
 PRIVATE
     REAL(DP), PARAMETER :: svd_eps = 1.0e-8_dp 
     TYPE(MethodType) :: reg 
     REAL(DP) ::  finaldeut, hb2m 
     REAL(DP),ALLOCATABLE :: xk(:), wk(:), biorth(:,:) !, vec_sc(:,:)
     INTEGER(I4B) :: nkpt, nmod, np_on, nktot 
     LOGICAL(LGT) :: coupled, deut


  PUBLIC :: hermitize

CONTAINS
      
SUBROUTINE hermitize(vlkmethod, pmesh, channel, vnh,vherm)
        IMPLICIT NONE
        TYPE(MethodType) :: vlkmethod
        TYPE(MeshType) :: pmesh
        TYPE(ChannelType) :: channel
        REAL(DP)::vnh(:,:),vherm(:,:),fkkw,dia
        REAL(DP),ALLOCATABLE::zinv(:,:),z(:,:),hnh(:,:),temp(:,:),evalr(:),evali(:),vecr(:,:),work(:)
        REAL(DP),ALLOCATABLE::u(:,:), vec_sc(:,:)
        INTEGER(I4B)::i,j,n,info,lwork,ii,jj,indx(100000)
        INTENT(IN)::vnh, vlkmethod, pmesh, channel
        INTENT(OUT)::vherm
        
!!!!!!! set module working variables and allocate module arrays !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      reg = vlkmethod
      coupled = GetCoupled(channel)
      hb2m = GetHb2m(channel)
      deut = GetDeut(channel)

      nkpt = GetNtot(pmesh) 
      nmod = GetNmod(pmesh)

      np_on = nmod; nktot = nkpt
      if(coupled)np_on = 2*nmod
      if(coupled)nktot = 2*nkpt

      IF(ALLOCATED(xk))DEALLOCATE(xk)
      IF(ALLOCATED(wk))DEALLOCATE(wk)
      ALLOCATE(xk(1:nkpt),wk(1:nkpt)); xk=0_dp ; wk=0_dp
      CALL GetMeshWts(pmesh,xk,wk)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        allocate(zinv(np_on,np_on),z(np_on,np_on),hnh(np_on,np_on), vec_sc(np_on,np_on)) 
        zinv=0.0_dp;z=0.0_dp;hnh=0.0_dp; vec_sc = 0.0_dp
        call recalcvectors(vnh,vec_sc)  ! calculates right vectors of hnh -> vec_sc(:)  and biorth complements -> biorth(:)

! attach kkww factors to build hnh
        do i=1,np_on
             do j=1,np_on
                   ii=i;jj=j
                   if(i.gt.nmod)ii=i-nmod
                   if(j.gt.nmod)jj=j-nmod  
                   dia=0_dp
                   if(i.eq.j)dia=hb2m*xk(ii)**2      
                   fkkw=dsqrt(wk(ii)*wk(jj))*xk(ii)*xk(jj)
                   hnh(i,j)=dia+vnh(i,j)*fkkw
             enddo
        enddo                 

! calculate the orthogonalizaton z and zinv transformation then calculate vherm 
        if(GetHerm(reg).eq.1)then
               call ztransgramschmidt(np_on,vec_sc,z,zinv)
        else
               call ztransP_plus_wdw(np_on,vec_sc,z,zinv)
        endif
        

        vherm=matmul(z,matmul(hnh,zinv)) !actually h_herm

!!!!!!!!!!!!!!!!!!!!!!!!! check the deuteron be is preserved!!!!!!!!!!!!!!!!!1
           if(deut)then
               lwork=40*np_on;info=100
               allocate(temp(np_on,np_on),evalr(np_on),evali(np_on),vecr(np_on,np_on),work(lwork),u(np_on,np_on))
               temp=hnh  !non-hermitian 
               call dgeev('N','V',np_on,temp,np_on,evalr,evali,u,np_on,vecr,np_on,work,lwork,info )
               if(info.ne.0)call abort 
               call indexx(np_on,evalr,indx(1:np_on))
               write(11,*)'non hermitian deuteron',evalr(indx(1))
! repeat w/hermitian hamiltonian
               temp=vherm;evalr=0_dp;indx=0
               call dgeev('N','V',np_on,temp,np_on,evalr,evali,u,np_on,vecr,np_on,work,lwork,info )
               if(info.ne.0)call abort 
               call indexx(np_on,evalr,indx(1:np_on))
               write(11,*)'hermitian deuteron',evalr(indx(1))
!               write(22,*)'eigenvector'
!               write(22,'(2f14.6)')(xk(i),vecr(i,indx(1))/xk(i)/dsqrt(wk(i)),i=1,nmod) 
               finaldeut=evalr(indx(1))    
               deallocate(temp,evalr,evali,vecr,work,u)
           endif
!!!!!!!!!!!!!! end of deuteron be check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! subtract off kinetic energy        
        do i=1,np_on
           ii=i
           if(i.gt.nmod)ii=i-nmod   
           vherm(i,i)=vherm(i,i)-hb2m*xk(ii)**2  !subtract k.e.
        enddo               

! remove kkww factors for final vherm
        do i=1,np_on
             do j=1,np_on
                  ii=i;jj=j
                     if(i.gt.nmod)ii=i-nmod
                     if(j.gt.nmod)jj=j-nmod  
                     fkkw=dsqrt(wk(ii)*wk(jj))*xk(ii)*xk(jj)
                     vherm(i,j)=vherm(i,j)/fkkw
             enddo
        enddo

        deallocate(hnh,z,zinv)
        deallocate(vec_sc, biorth, xk, wk)
   END SUBROUTINE hermitize 
    
   SUBROUTINE recalcvectors(vnh,vec_sc)
        USE nrtype
        IMPLICIT NONE
        REAL(DP)::vnh(:,:),vec_sc(:,:), fkkw,dia
        REAL(DP),ALLOCATABLE::hnh(:,:),evalr(:),evali(:),vecr(:,:),work(:),u(:,:),v(:,:),w(:)
        INTEGER(I4B)::i,j,ii,jj,info,lwork,indx(100000)
        INTENT(IN)::vnh
        INTENT(OUT) :: vec_sc

        lwork=40*np_on;info=100
        ALLOCATE(hnh(np_on,np_on),evalr(np_on),evali(np_on),vecr(np_on,np_on),work(lwork))
        ALLOCATE(u(np_on,np_on),v(np_on,np_on),w(np_on));w=0_dp;u=0_dp;v=0_dp

!   set up non-hermitian hamiltonian 
        do i=1,np_on
            do j=1,np_on
               ii=i;jj=j
               if(i.gt.nmod)ii=i-nmod
               if(j.gt.nmod)jj=j-nmod
               dia=0.0_dp
               if(i.eq.j)dia=hb2m*xk(ii)**2
               fkkw=dsqrt(wk(ii)*wk(jj))*xk(ii)*xk(jj) !*2.0_dp/pi_d
               hnh(i,j)=dia+vnh(i,j)*fkkw
!                   write(33,*)i,j,vnh(i,j)
            enddo
        enddo

!   solve for right eigenvectors of the non-hermitian eigenvalue problem           
        call dgeev('N','V',np_on,hnh,np_on,evalr,evali,u,np_on,vecr,np_on,work,lwork,info )
        if(info.ne.0)then
             write(11,*)'problem in recalcvectors info',info  
             call abort 
        endif  

!   sort the eigenvalues/vectors in ascending order
        call indexx(np_on,evalr,indx(1:np_on))
!        write(11,*)'lowest eigenvalue',evalr(indx(1))*hb2m
        do i=1,np_on
              vec_sc(:,i)=vecr(:,indx(i))
        enddo

!   use SV decomposition to build bi-orthogonal complement vectors
        if(allocated(biorth))deallocate(biorth)
        allocate(biorth(1:np_on,1:np_on));biorth=0_dp
        biorth = 0_dp
        u=vec_sc(:,:)
        u=transpose(vec_sc(:,:)) !temp
        call dgesvd('O','A',np_on,np_on,u,np_on,w,u,np_on,v,np_on,work,lwork,info )
        v=transpose(v) 

               if(info.ne.0)then
                    write(11,*)'problem in svd inside recalcvectors info=',info
                    call abort
               endif
                    do i=1,np_on 
!                         write(11,'(a,f14.8)')'recalcvector sing value',w(i)
                         if(w(i).le.svd_eps)write(11,*)'small singular value in nonhermeigvec'
                         if(w(i).gt.svd_eps)v(:,i)=v(:,i)/w(i)
                         if(w(i).le.svd_eps)v(:,i)=v(:,i)*0.0_dp
                    enddo
                    biorth(1:np_on,1:np_on)=matmul(v,transpose(u))
                    biorth=transpose(biorth) !temp
! test biorth
!               u=0_dp
!               u=matmul(vec_sc,biorth)
!               write(88,*)'******* inside recalcvectors *******'
!               do i=1,np_on
!                   do j=1,np_on
!                     if(i.eq.j) write(88,'(2i4,f14.6)')i,j,u(i,j)
!                   enddo
!               enddo

!        do i = 1,np_on
!           write(11,'(a,f14.6)')'eigvects', evalr(indx(i))
!           write(11,'(3f14.6)')(xk(j),vecr(j,indx(i)),vecr(j+nmod,indx(i)),j=1,nmod)
!           write(11,*) ' '
!       enddo
              
        deallocate(w,v,u,hnh,evalr,evali,vecr,work)
   END SUBROUTINE recalcvectors



   SUBROUTINE ztransP_plus_wdw(n,vec_sc,z,zinv)
       REAL(DP) :: vec_sc(:,:), z(:,:), zinv(:,:),ans1,ans2
       INTEGER(I4B) :: n,lwork,info,i,j,k,ipiv(400)
       REAL(DP), ALLOCATABLE :: P_plus_wdw(:,:), eval(:), evec(:,:),work(:),temp(:,:)      
       INTENT(IN) :: n, vec_sc
       INTENT(OUT) :: z, zinv

       ALLOCATE(P_plus_wdw(n,n),eval(n),evec(n,n),temp(n,n))


       P_plus_wdw = Matmul(Transpose(biorth),biorth)

       IF(GetHerm(reg).eq.2)THEN !okubo
           evec = P_plus_wdw ! input matrix destroyed; returns eigenvectors
           lwork = 40*n; ALLOCATE(work(lwork)) 
           call dsyev('V','U',n,evec,n,eval,work,lwork,info)
              if(info.ne.0)then
                write(11,'(a,i5)')'problem in dysev in ztransP_plus_wdw, info=',info
                call abort
              endif

             do i = 1, n
               do j = 1,n
                    ans1 =0_dp ; ans2 = 0_dp
                    do k =1, n
                       ans1 = ans1 + evec(i,k)/dsqrt(eval(k))*evec(j,k)
                       ans2 = ans2 + evec(i,k)*dsqrt(eval(k))*evec(j,k)
                    enddo
                    zinv(i,j) = ans1
                    z(i,j) = ans2 
               enddo
            enddo
            deallocate(work)
        ELSE IF (GetHerm(reg).eq.3)THEN !cholesky
            
             call dpotrf('L', n, p_plus_wdw, n, info)
                  if(info.ne.0)then
                     write(6,*)'trouble in dpotrf in ztrans',info
                     call abort
                  endif
             z = 0_dp; zinv=0_dp
             do i = 1, n
                do j = 1, i
                    z(i,j) = p_plus_wdw(i,j)
                enddo
             enddo
             z = transpose(z)
             do i = 1, n
                zinv(i,i) = 1.0_dp
             enddo
          
             temp = z
             call dgesv(n,n,temp,n,ipiv(1:n),zinv,n,info)
                 if(info.ne.0)then
                     write(6,*)'problem in dgesv in ztransP_plus_wdw', info
                     call abort
                 endif


        ELSE IF (GetHerm(reg).eq.4)THEN ! kato
             zinv = 0_dp; z =0_dp
             do i = 1, n
                zinv(i,i) = 1.0_dp
                z(i,i) = 1.0_dp
             enddo
             temp = p_plus_wdw 
             call dgesv(n,n,temp,n,ipiv(1:n),zinv,n,info)
                 if(info.ne.0)then
                     write(6,*)'problem in dgesv in ztransP_plus_wdw', info
                     call abort
                 endif
        ENDIF

       DEALLOCATE(P_plus_wdw,eval,evec,temp)
   END SUBROUTINE ztransP_plus_wdw 

   SUBROUTINE ztransgramschmidt(n,vec_sc,z,zinv)
        USE nrtype
        IMPLICIT NONE
        REAL(DP)::vec_sc(:,:),z(:,:),zinv(:,:)
        REAL(DP),ALLOCATABLE::temp(:,:),vmat(:,:),v(:,:),w(:),work(:)
        INTEGER(I4B)::i,j,n,info,ipiv(100000),lwork
        INTENT(IN)::n,vec_sc
        INTENT(OUT)::z,zinv

        allocate(temp(n,n),w(n),vmat(n,n),v(n,n));temp=0.0_dp;vmat=0.0_dp;v=0.0_dp;w=0.0_dp
         
! run modified gram-schmidt to get orthogonal set v from vec_sc
        call mgsreorth(vec_sc,n,v)

! construct m.e.'s of the z-transformation that connects vec_sc to v
        z=matmul(v,biorth)

! m.e.'s of the inverse z-transformation using svd
        temp=z;info=200;ipiv=0;zinv=0.0_dp
             lwork=40*n;allocate(work(lwork))
             call dgesvd('O','A',n,n,temp,n,w,temp,n,vmat,n,work,lwork,info )
             deallocate(work)
             vmat=transpose(vmat) 
               if(info.ne.0)then
                    write(11,*)'problem in svd inside ztrans info=',info
                    call abort
               endif
! take care of any small w(i)'s by setting 1/0=0
             do i=1,n
                 write(11,'(a,f18.10)')'singular value in ztrans',w(i)
                 if(w(i).lt.svd_eps)then
                    vmat(:,i)=vmat(:,i)*0.0_dp
                    write(11,'(a,f18.10)')'ABORT: small singular value in ztrans',w(i)
                    call abort
                 else
                     vmat(:,i)=vmat(:,i)/w(i)
                 endif
             enddo
        zinv=matmul(vmat,transpose(temp))  
        deallocate(temp,vmat,v,w)                
        RETURN
   END SUBROUTINE ztransgramschmidt

   SUBROUTINE mgs(a,n,q)
        USE nrtype
        IMPLICIT NONE
        INTEGER(I4B)::n,j,k
        REAL(DP)::a(:,:),q(:,:),xnorm
        REAL(DP),ALLOCATABLE::a1(:,:)
        INTENT(IN)::a,n
        INTENT(OUT)::q        
        allocate(a1(n,n));a1=0.0_dp

        do j=1,n
           a1(:,j)=a(:,j)
           do k=1,j-1
              a1(:,j)=a1(:,j)-Dot_Product(a1(:,j),q(:,k))*q(:,k)
           enddo
           xnorm=dsqrt(Dot_Product(a1(:,j),a1(:,j)))       
!           write(31,*)j,xnorm
!           if(dabs(xnorm).le.0.00001_dp)write(31,*)'small norm'
           q(:,j)=a1(:,j)/xnorm
        enddo  

        deallocate(a1)
        RETURN
   END SUBROUTINE mgs


   SUBROUTINE cgs(a,n,q)
        USE nrtype
        IMPLICIT NONE
        INTEGER(I4B)::n,j,k
        REAL(DP)::a(:,:),q(:,:)
        REAL(DP),ALLOCATABLE::a1(:,:)
        INTENT(IN)::a,n
        INTENT(OUT)::q        
        allocate(a1(n,n));a1=0.0_dp

        do j=1,n
           a1(:,j)=a(:,j)
           do k=1,j-1
              a1(:,j)=a1(:,j)-Dot_Product(a(:,j),q(:,k))*q(:,k)
           enddo
           q(:,j)=a1(:,j)/dsqrt(Dot_Product(a1(:,j),a1(:,j)))       
        enddo  

        deallocate(a1)
        RETURN
   END SUBROUTINE cgs

   SUBROUTINE mgsreorth(a,n,q)
        USE nrtype
        IMPLICIT NONE
        INTEGER(I4B)::n,j,k,i
        REAL(DP)::a(:,:),q(:,:),xnorm
        REAL(DP),ALLOCATABLE::a1(:,:,:)
        INTENT(IN)::a,n
        INTENT(OUT)::q        
        allocate(a1(n,n,0:2));a1=0.0_dp

        a1(:,:,0)=a(:,:)

        do j=1,n
           do i=1,2
              a1(:,j,i)=a1(:,j,i-1)
              do k=1,j-1
                 a1(:,j,i)=a1(:,j,i)-Dot_Product(a1(:,j,i),q(:,k))*q(:,k)
              enddo
           enddo
           xnorm=dsqrt(Dot_Product(a1(:,j,2),a1(:,j,2)))      
!           write(31,*)xnorm 
!           if(dabs(xnorm).le.0.0001_dp)write(31,*)'small norm'
           q(:,j)=a1(:,j,2)/xnorm
        enddo  

        deallocate(a1)
        RETURN
   END SUBROUTINE mgsreorth

   SUBROUTINE cgsreorth(a,n,q)
        USE nrtype
        IMPLICIT NONE
        INTEGER(I4B)::n,j,k,i
        REAL(DP)::a(:,:),q(:,:),xnorm
        REAL(DP),ALLOCATABLE::a1(:,:,:)
        INTENT(IN)::a,n
        INTENT(OUT)::q        
        allocate(a1(n,n,0:2));a1=0.0_dp

        a1(:,:,0)=a(:,:)

        do j=1,n
           do i=1,2
              a1(:,j,i)=a1(:,j,i-1)
              do k=1,j-1
                 a1(:,j,i)=a1(:,j,i)-Dot_Product(a1(:,j,i-1),q(:,k))*q(:,k)
              enddo
           enddo
           xnorm=dsqrt(Dot_Product(a1(:,j,2),a1(:,j,2)))      
!           write(31,*)xnorm 
!           if(dabs(xnorm).le.0.0001_dp)write(31,*)'small norm'
           q(:,j)=a1(:,j,2)/xnorm
        enddo  

        deallocate(a1)
        RETURN
   END SUBROUTINE cgsreorth

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
 
END MODULE hermitization 
