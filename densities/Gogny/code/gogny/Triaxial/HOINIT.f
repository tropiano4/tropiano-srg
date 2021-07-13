      SUBROUTINE HOINIT(UV,AZ,AN)
      IMPLICIT REAL*8(A-H,O-Z)
      Include 'DIMTRIAX'
      Include 'MAXDIM'
      
      Dimension UV(NUV)
      Dimension ee(NP),eo(NM),v2e(NP,2),v2o(NM,2)
      
      Dimension ie(NP),io(NM)
      
      Common /CQN/    IQMAX,IQNE(3,NPMAX),IQNO(3,NMMAX),IFILL
      Common /OSCLEN/ bx,by,bz
      Common /ITLOC0/ ND(4),ND2(4),NBLOCK(4)
      Common /ITLOC1/ NNU(4),NNV(4),NNRO(4),NNKA(4)

c
      if(IFILL.eq.0) call setQN()
      
      hc2=197.32891d0**2
      manuc=938.9059d0
      omex=8.0d0 !hc2/manuc/bx**2
      omey=8.0d0 !hc2/manuc/by**2
      omez=8.0d0 !hc2/manuc/bz**2

      write(6,*) ' HOINIT omex, omey, omez  ',omex,omey,omez
      write(6,*) ' HOINIT ***************** ',AZ,AN
      write(6,*) ' HOINIT ***************** ',ND(1),ND(2),ND(3),ND(4)
      write(6,*) ' HOINIT ***************** ',NP,NM,NUV  
c---------------------------------------------------------------------- +
c     TRIAXIAL HARMONIC OSCILLATOR
c+---------------------------------------------------------------------+
c                  Even parity
c
c  --- Note that the quantum number given by iqne are the real one + 1
c+---------------------------------------------------------------------+
      do ia=1,NP
         nxa = iqne(1,ia)
         nya = iqne(2,ia)
         nza = iqne(3,ia)
         
         ee(ia) = omex*(dfloat(nxa-1)+0.5d+00) +
     &            omey*(dfloat(nya-1)+0.5d+00) +    
     &            omez*(dfloat(nza-1)+0.5d+00)
      end do       ! ia loop
c
c+---------------------------------------------------------------------+
c                  Odd parity
c
c  --- Note that the quantum number given by iqne are the real one + 1
c+---------------------------------------------------------------------+
c
      do ia=1,NM
         nxa = iqno(1,ia)
         nya = iqno(2,ia)
         nza = iqno(3,ia)
         
         eo(ia) = omex*(dfloat(nxa-1)+0.5d+00) +
     &            omey*(dfloat(nya-1)+0.5d+00) +    
     &            omez*(dfloat(nza-1)+0.5d+00)
      end do       ! ia loop
c



c
      call indexx(np,ee,ie)
      call indexx(nm,eo,io)
      
      NZ = AZ + 0.1d+00
      NN = AN + 0.1d+00
c      
      ke=1
      ko=1
      do in=1,NZ/2
        if(ee(ie(ke)).le.eo(io(ko)))then
            v2e(ie(ke),1)=1.0d0
c            write(6,*) ' e ',ke,ee(ie(ke))
            ke=ke+1
         else
            v2o(io(ko),1)=1.0d0
c            write(6,*) 'o ',ko,eo(io(ko))
            ko=ko+1
         end if
       end do
c       
      alambdap=max(ee(ie(ke-1)),eo(io(ko-1)))
      write(6,*) ' ALAMBDA P ',alambdap  
      
      ke=1
      ko=1
      do in=1,NN/2
        if(ee(ie(ke)).le.eo(io(ko)))then
            v2e(ie(ke),2)=1.0d0
            ke=ke+1
         else
            v2o(io(ko),2)=1.0d0
            ko=ko+1
         end if
       end do
       
      alambdan=max(ee(ie(ke-1)),eo(io(ko-1)))
      write(6,*) ' ALAMBDA P ',alambdap  
c
c
c      
      DeltaP = 1.0 d+00
      DeltaN = 1.0 d+00
      
      ZZ = 0.0d+00
      do ia=1,NP
         eep = ee(ia)-alambdap
         een = ee(ia)-alambdan
         v2e(ia,1) = 0.5d+00*(1.0d+00 - eep/dsqrt(eep**2+DeltaP**2))
         v2e(ia,2) = 0.5d+00*(1.0d+00 - een/dsqrt(een**2+DeltaN**2))
         ZZ = ZZ + v2e(ia,1)
c         write(6,*) v2e(ia,1)
      end do
     
      do ia=1,NM
         eep = eo(ia)-alambdap
         een = eo(ia)-alambdan
         v2o(ia,1) = 0.5d+00*(1.0d+00 - eep/dsqrt(eep**2+DeltaP**2))
         v2o(ia,2) = 0.5d+00*(1.0d+00 - een/dsqrt(een**2+DeltaN**2))
         ZZ = ZZ + v2o(ia,1)
c         write(6,*) v2o(ia,1)
      end do
   
      write(6,*) ' HOINIT ***** ',ZZ
c
c     CALCULATION OF U AND V MATRICES
c
c          UV => U1PP U2PP U1PM U2PM  V1PP V2PP V1PM V2PM     
c                U1NP U2NP U1NM U2NM  V1NP V2NP V1NM V2NM 
    
      do i=1,NUV
         UV(i) = 0.0d+00
      end do
            
      do ia=1,NP
         ja = ia - 1
         vvp = dsqrt(v2e(ia,1))
         uup = dsqrt(max(1.0d+00-v2e(ia,1),0.0d+00))
         
         uv(NNU(1)        + ja*np + ja)= uup
         uv(NNU(1)+ND2(1) + ja*np + ja)= uup
         uv(NNV(1)        + ja*np + ja)= vvp
         uv(NNV(1)+ND2(1) + ja*np + ja)=-vvp
         
         vvn = dsqrt(v2e(ia,2))
         uun = dsqrt(max(1.0d+00-v2e(ia,2),0.0d+00))
         
         uv(NNU(3)        + ja*np + ja)= uun
         uv(NNU(3)+ND2(3) + ja*np + ja)= uun
         uv(NNV(3)        + ja*np + ja)= vvn
         uv(NNV(3)+ND2(3) + ja*np + ja)=-vvn
      end do

       do ia=1,NM
         ja = ia -1 
         vvp = dsqrt(v2o(ia,1))
         uup = dsqrt(max(1.0d+00-v2o(ia,1),0.0d+00))
         
         uv(NNU(2)        + ja*nm + ja)= uup
         uv(NNU(2)+ND2(2) + ja*nm + ja)= uup
         uv(NNV(2)        + ja*nm + ja)= vvp
         uv(NNV(2)+ND2(2) + ja*nm + ja)=-vvp
         
         vvn = dsqrt(v2o(ia,2))
         uun = dsqrt(max(1.0d+00-v2o(ia,2),0.0d+00))
         
         uv(NNU(4)        + ja*nm + ja)= uun
         uv(NNU(4)+ND2(4) + ja*nm + ja)= uun
         uv(NNV(4)        + ja*nm + ja)= vvn
         uv(NNV(4)+ND2(4) + ja*nm + ja)=-vvn
      end do
      
      return
      end
c+---------------------------------------------------------------------+
c|   Numerical Recipes  subroutine INDEXX                              |
c|   (C) Copr. 1986-92 Numerical Recipes Software.                     |
c+---------------------------------------------------------------------+
      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      DOUBLE PRECISION arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      DOUBLE PRECISION a
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
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
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
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
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
5       indx(l)=indx(j)
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
      END
