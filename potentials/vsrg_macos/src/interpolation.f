      module interpolation
      use nrtype
      implicit none
      
      contains
	  subroutine vkku399lag(nless,vkkless,nr,ra,vkk2)
!      implicit real(a-h,o-z)
      REAL(DP) :: VKKLESS(:,:),VKK2(:,:),RA(:), xi, xj, ans
      INTEGER(I4B):: i,j, nr, nless 
      INTENT(IN)::NLESS,VKKLESS,RA,NR
	  INTENT(OUT)::VKK2     
!     3/20/99 last point by interpolation
      do 40 j=1,nr
      do 50 i=1,nr
      if(i.ne.nr.and.j.ne.nr)then
      vkk2(i,j)=vkkless(i,j)
      endif

      if(i.eq.nr.or.j.eq.nr)then
      xi=ra(i)
      xj=ra(j)
      call interp2d(vkk2,ra,ra,nless,4,xi,xj,ans)
      vkk2(i,j)=ans
      endif
50    continue
40    continue      

      return
      end subroutine

!  	subroutine nonsquareinterp2d(z,x,y,nx,ny,ip,x0,y0,ans)
!      z=z(x,y)=<x\z\y>
!      ip=no of points used for interpolation
!      implicit real(a-h,o-z)

!      REAL:: z(:,:),y(:),x(:)
!	  REAL,ALLOCATABLE::zrow(:),col(:)
!      ALLOCATE(ZROW(SIZE(Z)),COL(SIZE(Z)));ZROW=0.;COL=0.

!	  do jcol=1,ny

!      do i=1,nx
!      col(i)=z(i,jcol)
!      enddo

!      call myinterp1d(col,x,nx,ip,x0,zrow(jcol))
!      enddo

!      call myinterp1d(zrow,y,ny,ip,y0,ans)
!      DEALLOCATE(ZROW,COL)      
!      return
!      end subroutine

  	subroutine interp2d(z,x,y,n,ip,x0,y0,ans)
!      z=z(x,y)=<x\z\y>
!      ip=no of points used for interpolation
!      implicit real(a-h,o-z)
      REAL(DP) :: z(:,:),y(:),x(:), ans, x0, y0
      INTEGER(I4B) :: n, ip, i, jcol
      REAL(DP),ALLOCATABLE::zrow(:),col(:)
      ALLOCATE(ZROW(SIZE(Z)),COL(SIZE(Z)));ZROW=0.;COL=0.

	  do jcol=1,n

      do i=1,n
      col(i)=z(i,jcol)
      enddo

      call myinterp1d(col,x,n,ip,x0,zrow(jcol))
      enddo

      call myinterp1d(zrow,y,n,ip,y0,ans)
      DEALLOCATE(ZROW,COL)      
      return
      end subroutine

      subroutine myinterp1d(yy,xx,npoint,iorder,x0,y0)
      REAL(DP):: yy(:),xx(:),x0,y0,dy
      INTEGER(I4B) :: npoint,iorder,k,k2,j

!      call locate(xx,npoint,x0,j)
      call hunt(xx,npoint,x0,j)
      k=min(max(j-(iorder-1)/2,1),npoint+1-iorder)

      k2=k+iorder-1
      call polint(xx(k:k2),yy(k:k2),iorder,x0,y0,dy)
!      call ratint(xx(k:k2),yy(k:k2),iorder,x0,y0,dy)

      return
      end subroutine

      SUBROUTINE locate(xx,n,x,j)
      INTEGER(I4B):: j,n
      REAL(DP):: x,xx(n)
      INTEGER(I4B):: jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      if(x.eq.xx(1))then
        j=1
      else if(x.eq.xx(n))then
        j=n-1
      else
        j=jl
      endif
      return
      END subroutine 
      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER(I4B) :: jlo,n
      REAL(DP):: x,xx(n)
      INTEGER(I4B):: inc,jhi,jm
      LOGICAL(LGT):: ascnd
      ascnd=xx(n).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)then
        if(x.eq.xx(n))jlo=n-1
        if(x.eq.xx(1))jlo=1
        return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END subroutine
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER(I4B):: n,NMAX
      REAL(DP) :: dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER(I4B)::i,m,ns
      REAL(DP):: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END subroutine 
      SUBROUTINE ratint(xa,ya,n,x,y,dy)
      INTEGER (I4B) :: n,NMAX
      REAL (DP) :: dy,x,y,xa(n),ya(n),TINY
      PARAMETER (NMAX=10,TINY=1.e-25)
      INTEGER (I4B)::i,m,ns
      REAL(DP):: dd,h,hh,t,w,c(NMAX),d(NMAX)
      ns=1
      hh=abs(x-xa(1))
      do 11 i=1,n
        h=abs(x-xa(i))
        if (h.eq.0.)then
          y=ya(i)
          dy=0.0
          return
        else if (h.lt.hh) then
          ns=i
          hh=h
        endif
        c(i)=ya(i)
        d(i)=ya(i)+TINY
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          w=c(i+1)-d(i)
          h=xa(i+m)-x
          t=(xa(i)-x)*d(i)/h
          dd=t-c(i+1)
          if(dd.eq.0.)pause 'failure in ratint'
          dd=w/dd
          d(i)=c(i+1)*dd
          c(i)=t*dd
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END subroutine
      end module
