May, 2017

Modifying allvnnmodels_r8.f for the new Entem-Machleidt code.

* Is empotkk supposed to be commented out or not????
* Does this depend on the code?

* Already added modules to all of the individual sets of routines provided for each 
   order and each cutoff (15 in total).

* Module vnn_module [4666-24769] is everything except top cdbonn module
!     vnnmompw returns matrix vnn(i,j)=<i|v|j>*sqrt(xwkp(i)*xwkp(j))*xkp(i)*xkp(j)
!     where the units are such that vnn(i,j) is in MeV. The important variables are
!     icoup= .true./.false. for coupled/uncoupled channels
!     kvnn= integer controls which vnn model is used (see below for assignments)
!     ll,is,jt,it,itz1,itz2= partial wave quantum numbers 
!     xkp(i)= momentum-space grid points (in 1/fm)
!     xwkp(i)= corresponding momentum-space quadrature weights (in 1/fm) 

! for both uncoupled and coupled channels
        if (kvnn.eq.10) call n3loidaho(icoup,it,itz,ll,is,jt,nkp,xkp,xwkp) 
        if (kvnn.eq.70) call wrapper_lo450(icoup,it,itz,ll,is,jt,nkp,xkp,xwkp)  

* change this to
    if ((kvnn.ge.70).AND.(kvnn.le.84)) call wrapper_newEM(icoup,it,itz,ll,is,jt,nkp,xkp,xwkp,kvnn)     
  and create wrapper_newEM just like n3loidaho except passes kvnn to interior subroutine POT_newEM_LSJ

      if((kvnn.ne.7).and.(kvnn.ne.6).and.(kvnn.ne.3).and.(kvnn.ne.9).and.(kvnn.ne.10).and.
     1   (kvnn.ne.11).and.(kvnn.ne.12).and.(kvnn.ne.13).and.(kvnn.ne.15)
     1   .and.((kvnn.lt.70).or.(kvnn.gt.84))                 
     1   )
     1   call svdivide1(v,nkp,xkp,xwkp,nkp,vnn_kkww,nd1)!v(i)to vbare(i,j)


! vnnmompw vs. vnnmompwn3lo
!   we only need vnnmompw


