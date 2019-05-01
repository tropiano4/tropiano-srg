MODULE vlowk_sharp_cutoff
 USE nrtype
 USE vlowkperipherals

 IMPLICIT NONE
 PRIVATE

     REAL(DP),ALLOCATABLE :: xk(:), wk(:), vnntemp(:,:), h0(:)
     INTEGER(I4B) :: nkpt, nmod, np_on, nktot, nqsize
     LOGICAL(LGT) :: coupled, deut

   
     INTERFACE leesuzukivlowk
        MODULE PROCEDURE leesuzuki_long_argument_list, leesuzuki_short_list
     END INTERFACE

     PUBLIC :: leesuzukivlowk
         
CONTAINS


SUBROUTINE leesuzuki_long_argument_list(ntotpts, nmodpts, nsizeV, nsizeVlk, &
                                    xxk, wwk, hb2m, vnn, vlk_nonherm, vlk_herm)
! Use this form to CALL leesuzuki(....) if you want to generate sharp cutoff vlowk 
! without using the derived type data structures (e.g., MeshType, ChannelType, etc...) 
 
      REAL(DP) :: diag, kkww, vlk_herm(:,:), vlk_nonherm(:,:), vnn(:,:), xxk(:), wwk(:), hb2m
      INTEGER(I4B) :: i, j, ii, jj, ntotpts, nmodpts, nsizeV, nsizeVlk
      REAL(DP), ALLOCATABLE :: h_full(:,:)
      INTENT(IN) :: vnn, hb2m, xxk, wwk, ntotpts, nmodpts, nsizeV, nsizeVlk
      INTENT(OUT) :: vlk_herm, vlk_nonherm 


!!!!!!!!!!!! set module variables and allocate module arrays!!!!!!!!!!
     nkpt  = ntotpts
     nmod  = nmodpts
     nktot = nsizeV
     np_on = nsizeVlk
     nqsize = nktot - np_on

     coupled = .false.
     IF(nktot.eq.2*nkpt)coupled = .true.

     IF(ALLOCATED(xk))DEALLOCATE(xk)
     IF(ALLOCATED(wk))DEALLOCATE(wk)
     IF(ALLOCATED(vnntemp))DEALLOCATE(vnntemp)
     IF(ALLOCATED(h0))DEALLOCATE(h0)

     ALLOCATE(vnntemp(nktot,nktot),h0(nktot), xk(nkpt), wk(nkpt))

     xk(1:nkpt) = xxk(1:nkpt)
     wk(1:nkpt) = wwk(1:nkpt)

     CALL prepare_vnn(hb2m, vnn) !generates module arrays vnntemp and h0 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     IF(ALLOCATED(h_full))DEALLOCATE(h_full)
     ALLOCATE(h_full(nktot,nktot))

          DO i=1,nktot
             DO j=1,nktot
                 diag=0.d0
                 IF(i.EQ.j)diag=h0(i)
                 h_full(i,j)=diag+vnntemp(i,j)
             ENDDO
          ENDDO

     CALL noniterls(h_full,vlk_nonherm,vlk_herm)  
    

     DEALLOCATE(vnntemp, h0, h_full, xk, wk)

END SUBROUTINE leesuzuki_long_argument_list


SUBROUTINE leesuzuki_short_list(Kmesh, Channel, vnn, vlk_nonherm, vlk_herm)
! Use this form to CALL leesuzuki(....) if you want to generate sharp cutoff vlowk 
! using the derived type data structures (e.g., MeshType, ChannelType, etc...) 

      TYPE(MeshType) :: Kmesh
      TYPE(ChannelType) :: Channel 
      REAL(DP) :: diag, kkww, vlk_herm(:,:), vlk_nonherm(:,:), vnn(:,:) 
      INTEGER(I4B) :: i, j, ii, jj 
      REAL(DP), ALLOCATABLE :: h_full(:,:)
      INTENT(IN) :: vnn, Kmesh, Channel
      INTENT(OUT) :: vlk_herm, vlk_nonherm 


!!!!!!!!!!!! set module variables and allocate module arrays!!!!!!!!!!
     nkpt  = GetNtot(Kmesh)
     nmod  = GetNmod(Kmesh) 

     nktot = nkpt
     np_on = nmod
     coupled = GetCoupled(Channel)
     IF(coupled)THEN
         nktot = 2*nkpt 
         np_on = 2*nmod
     ENDIF
     
     nqsize = nktot - np_on

     IF(ALLOCATED(xk))DEALLOCATE(xk)
     IF(ALLOCATED(wk))DEALLOCATE(wk)
     IF(ALLOCATED(vnntemp))DEALLOCATE(vnntemp)
     IF(ALLOCATED(h0))DEALLOCATE(h0)

     ALLOCATE(vnntemp(nktot,nktot),h0(nktot), xk(nkpt), wk(nkpt))

     CALL GetMeshWts(Kmesh,xk,wk)
     
     CALL prepare_vnn(GetHb2m(Channel), vnn) !generates module arrays vnntemp and h0 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     IF(ALLOCATED(h_full))DEALLOCATE(h_full)
     ALLOCATE(h_full(nktot,nktot))

          DO i=1,nktot
             DO j=1,nktot
                 diag=0.d0
                 IF(i.EQ.j)diag=h0(i)
                 h_full(i,j)=diag+vnntemp(i,j)
             ENDDO
          ENDDO

     CALL noniterls(h_full,vlk_nonherm,vlk_herm)  
    

     DEALLOCATE(vnntemp, h0, h_full, xk, wk)
END SUBROUTINE leesuzuki_short_list

  SUBROUTINE prepare_vnn(hb2m, xvnn)
      INTEGER(I4B) :: nq,i,j,ii
      REAL(DP) :: xvnn(:,:), hb2m
      INTENT(in)::xvnn

      nq=nkpt-nmod

      IF(.not.coupled)THEN
           vnntemp(1:nktot,1:nktot) = xvnn(1:nktot,1:nktot)
           DO i = 1, nktot
               h0(i)=hb2m*xk(i)**2
           ENDDO

      ELSE 
!   shuffle the coupled channel vnn into PVP, PVQ, QVQ, QVQ blocks
!   PVP block
           DO i=1,nmod
              DO j=1,nmod
                  vnntemp(i,j)=xvnn(i,j)  !SS
                  vnntemp(i,j+nmod)=xvnn(i,j+nkpt) !SD
                  vnntemp(i+nmod,j)=xvnn(i+nkpt,j) !DS
                  vnntemp(i+nmod,j+nmod)=xvnn(i+nkpt,j+nkpt) !DD
              ENDDO
           ENDDO
!   QVQ block
           DO i=1,nq
              DO j=1,nq
                  vnntemp(i+2*nmod,j+2*nmod)=xvnn(i+nmod,j+nmod) !SS
                  vnntemp(i+2*nmod,j+2*nmod+nq)=xvnn(i+nmod,j+nkpt+nmod) !SD
                  vnntemp(i+2*nmod+nq,j+2*nmod)=xvnn(nkpt+nmod+i,nmod+j) !DS
                  vnntemp(i+2*nmod+nq,j+2*nmod+nq)=xvnn(nkpt+nmod+i,nkpt+nmod+j)!DD
              ENDDO
           ENDDO

!   QVP block
           DO i=1,nq
              DO j=1,nmod
                  vnntemp(i+2*nmod,j)=xvnn(i+nmod,j) !SS
                  vnntemp(i+2*nmod,j+nmod)=xvnn(i+nmod,j+nkpt)!SD
                  vnntemp(i+2*nmod+nq,j)=xvnn(i+nkpt+nmod,j) !DS
                  vnntemp(i+2*nmod+nq,j+nmod)=xvnn(i+nkpt+nmod,j+nkpt)!DD
              ENDDO
           ENDDO

!   PVQ block
           DO i=1,nmod
              DO j=1,nq
                  vnntemp(i,j+2*nmod)=xvnn(i,j+nmod) !SS
                  vnntemp(i,j+2*nmod+nq)=xvnn(i,j+nkpt+nmod) !SD
                  vnntemp(i+nmod,j+2*nmod)=xvnn(i+nkpt,j+nmod) !DS
                  vnntemp(i+nmod,j+2*nmod+nq)=xvnn(i+nkpt,j+nkpt+nmod)!DD
              ENDDO
           ENDDO
!   kinetic energy 
           DO i=1,2*nmod
                ii=i
                IF(ii.GT.nmod)ii=i-nmod
                h0(i)=hb2m*xk(ii)**2  !note PH0P has Ph0P doubled (SS,DD)
           ENDDO
           
           DO i=1,2*nkpt-2*nmod
                ii=i
                IF(ii.GT.nq)ii=i-nq
                h0(i+2*nmod)=hb2m*xk(ii+nmod)**2 !note that QH0Q has SS and DD 
           ENDDO
      ENDIF 
      RETURN
  END SUBROUTINE prepare_vnn



  SUBROUTINE noniterls(hnn,vnonherm,vherm)
     INTEGER(I4B) :: lwork,info,nsize,np,nq,i,j,jj,ii
     REAL(DP) :: hnn(:,:),vnonherm(:,:),vherm(:,:),kkww 
     REAL(DP), ALLOCATABLE :: p_psi_n(:,:), q_psi_n(:,:),omega(:,:),work(:),evec(:,:),eval(:),one(:,:),heff(:,:)
     INTEGER(I4B), ALLOCATABLE :: indx(:),ipiv(:)
     CHARACTER*1 jobz,uplo
     INTENT(IN) :: hnn
     INTENT(OUT) :: vnonherm,vherm

        nq=nqsize
        np=np_on
        nsize=nktot

        ALLOCATE(p_psi_n(np,np),q_psi_n(nq,np),omega(nq,np),evec(nsize,nsize),eval(nsize),indx(nsize),heff(np,np))
 
!   calculate eigenvalues/vectors of the bare 2body hamiltonian 
        lwork=3*nsize-1;info=100;ALLOCATE(work(lwork));jobz='V';uplo='U'
            evec(1:nsize,1:nsize)=hnn(1:nsize,1:nsize)  ! input matrix destroyed returning eigenvectors
                CALL dsyev(jobz,uplo,nsize,evec,nsize,eval,work,lwork,info)
                   DEALLOCATE(work)

                   if(info.ne.0)then
                        write(*,*)'PROBLEM IN DSYEV INSIDE NONITERLS'
                        call abort
                    endif

!  sort eval in ascending order and store pointer in array indx(:)
        call indexx(nsize,eval,indx)
        write(11,'(f14.6)')(eval(indx(i)), i = 1,np_on)
!  set up p_psi_n and q_psi_n matrices
        do i=1,np
           p_psi_n(1:np,i)=evec(1:np,indx(i))
              q_psi_n(1:nq,i)=evec(np+1:nsize,indx(i))
                 enddo


!   write(97,'(2f14.8)')(xk(i),evec(i,indx(1))/xk(i)/sqrt(wk(i)),i=1,nkpt)
!  unit matrix for RHS; written over by inverse of <p|Psi_n> after dgesv 
        allocate(one(np,np),ipiv(np))
        one(:,:)=0.0
        do i=1,np
           one(i,i)=1.0
        enddo
! solve for inverse of <p|Psi_n> matrix, returned in one(:,:)  
             info=100
             CALL dgesv(np,np,p_psi_n,np,ipiv,one,np,info) 
             if(info.ne.0)then
                 write(*,*)'problem in DGESV inside NONITERLS',info
                 call abort
             endif

! construct wave operator w(q,p) = Sum_n <q|Psi_n><p|Psi_n>**-1
        omega(1:nq,1:np)=MatMul(q_psi_n(1:nq,1:np),one(1:np,1:np))

! construct (non-hermitian) lowk hamiltonian
        heff(1:np,1:np)=hnn(1:np,1:np)+Matmul(hnn(1:np,np+1:nsize),omega(1:nq,1:np))

! construct non-hermitian vlowk by subtracting h0(:)
        vnonherm = heff
        do i=1,np
            vnonherm(i,i)=vnonherm(i,i)-h0(i)
        enddo

! remove kkww factors from non-hermitian vlowk
        do i=1,np
             ii=i;IF(ii.GT.nmod)ii=i-nmod
             do j=1,np
                  jj=j;IF(jj.GT.nmod)jj=j-nmod
                  kkww=DSQRT(wk(ii)*wk(jj))*xk(ii)*xk(jj)
                  vnonherm(i,j)=vnonherm(i,j)/kkww
             enddo
        enddo


! hermitize the lowk hamiltonian
        CALL veff_okubo(np,nq,omega,heff)     !okubo  
!        IF(herm)CALL veff_herm(np,nq,omega,heff)     !andreozzi based on cholesky

! construct hermitian vlowk by subtracting h0(:)
        vherm = heff
        do i=1,np
            vherm(i,i)=heff(i,i)-h0(i)
        enddo

! remove kkww factors from hermitian vlowk
        do i=1,np
             ii=i;IF(ii.GT.nmod)ii=i-nmod
             do j=1,np
                  jj=j;IF(jj.GT.nmod)jj=j-nmod
                  kkww=DSQRT(wk(ii)*wk(jj))*xk(ii)*xk(jj)
                  vherm(i,j)=vherm(i,j)/kkww
             enddo
        enddo


        deallocate(heff,p_psi_n,q_psi_n,omega,evec,eval,indx,ipiv,one)
        return
        end subroutine

        SUBROUTINE veff_okubo(np,nq,omega,pheff)
!       Hermitize via the Okamoto/Suzuki method ( which is eq. to Okubo )
        INTEGER(I4B) :: np,nq,i,j,lwork,info
        REAL(DP) :: omega(:,:),pheff(:,:),diag,z
        REAL(DP),ALLOCATABLE :: wdw(:,:),vec(:,:),eval(:)
        REAL(DP),ALLOCATABLE :: pheff_eig(:,:), D_L(:,:),D_R(:,:),work(:)
        CHARACTER*1 jobz,uplo
        INTENT(IN)::np,nq,omega
        INTENT(INOUT)::pheff

        ALLOCATE(wdw(np,np),vec(np,np),pheff_eig(np,np),eval(np))
        ALLOCATE( D_L(np,np),D_R(np,np))
        wdw=0;vec=0;pheff_eig=0;eval=0
        wdw(1:np,1:np)=MATMUL(TRANSPOSE(omega(:nq,:np)),omega(1:nq,1:np))
        
       lwork=3*np-1
        info=100
        ALLOCATE(work(lwork))
        jobz='V'
        uplo='U'
        CALL dsyev(jobz,uplo,np,wdw,np,eval,work,lwork,info)
        DEALLOCATE(work)
!        CALL devcsf(np,wdw,np,eval,vec,np)
        vec(:np,:np)=wdw(:np,:np)
        pheff_eig(:np,:np)=MATMUL(TRANSPOSE(vec(:np,:np)),MATMUL(pheff(:np,:np),vec(:np,:np)))

        DO i=1,np
	DO j=1,np
		D_L(i,j)=SQRT(1.+eval(j))*vec(i,j)
        D_R(i,j)=1./SQRT(1.+eval(j))*vec(i,j)
	ENDDO
        ENDDO
        pheff(1:np,1:np)=MATMUL(D_L(1:np,1:np),MATMUL(pheff_eig(1:np,1:np),TRANSPOSE(D_R(1:np,1:np))))
        DEALLOCATE(wdw,vec,pheff_eig,eval,D_L,D_R)
        RETURN
   END SUBROUTINE veff_okubo

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	


END MODULE vlowk_sharp_cutoff 
