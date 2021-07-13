c+---------------------------------------------------------------------+
c|  First version May   1999                                           |
c+=====================================================================+
c|                                                                     |
c|               S P I N   O R B I T   P A I R I N G                   |
c|                                                                     |
c+---------------------------------------------------------------------+
      Subroutine Prepaiso(isym,Soxyz,GROP,GROM)
c
      Implicit real*8(a-h,o-z)
      Implicit logical (l)
      Include 'COMDIM'
c-------------------------------------------------------------- internal
      Parameter (KAXLSO=((imax*(imax+1))/2)*(imax+3))
      Parameter (KAXKSO= imax*(imax+1)*(imax+3))
      Dimension ALSOMU(kaxlso),AKSOMU(kaxkso)
      Dimension ILSOMU(iMAX,iMAX),IKSOMU(iMAX,iMAX)
      Parameter (NSOZ=12,NSOYZ=18,NSOXYZ=24)
      Dimension Soz (izsrc,NSOZ)
      Dimension Soyz(iyzsrc,NSOYZ)
c--------------------------------------------------------- end internal
      Dimension Soxyz(ndmu,NSOXYZ)
c +---------------------------------------------------------------------+
c |   In  Soxyz the meaning of the second index is                      |
c |=====================================================================|
c |                         nxmu     nymu     nzmu                      |
c |      1,2 .... T yz (2)  even     even      odd                      |
c |      3,4 .... T zy (2)  even      odd     even                      |
c |      5,6 .... T xz (2)   odd      odd      odd                      |
c |      7,8 .... T xy (2)   odd      odd      odd                      |
c +---------------------------------------------------------------------+
c |      9,10.... T xz (1)  even     even      odd                      |
c |      11,12... T xy (3)  even      odd     even                      |
c |      13,14... T yx (3)   odd     even     even                      |
c |      15,16... T zx (1)   odd     even     even                      |
c |      17,18... T yz (1)   odd      odd      odd                      |
c |      19,20... T yx (1)   odd      odd      odd                      |
c |      21,22... T zy (3)   odd      odd      odd                      |
c |      23,24... T zx (3)   odd      odd      odd                      |
c +---------------------------------------------------------------------+
c |      Odd indices are for neutrons while even ones are for protons   |
c +---------------------------------------------------------------------+
      Dimension GROP(nrop8),GROM(nrom8)
c
      Save alsomu,aksomu,ilsomu,iksomu,icall
c      
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      Common/DIMEN/kmax,kmax1,kxmax,kymax,kzmax,nwf2,
     * maxtz1,maxjz1,maxjx1,maxjy1,maxlso,maxdso,
     * ndmu,ndthet,nacou,nrop,nrom,nrop8,nrom8,nherm,nherm38,nlegn
     
      Common/OSCLEN/bx,by,bz
      
c      write(6,*) ' PREPAISO ISYM = ',isym
c-------------------------------------------------------------------
      if(icall.ne.061060) then
      call inigf
      Maxmu = 2*nmax+2
      indls = 1
      do i=1,nmax
         do j=i,nmax
	    ilsomu(i,j) = indls
	    ilsomu(j,i) = indls
            Minmu  = Mod(i+j,2)+1
	    do imu=Minmu,Maxmu,2
	       alsomu(indls) = DL1(i-1,j-1,imu-1)
	       indls = indls + 1
	    end do
	 end do
      end do
      if((indls-1).gt.KAXLSO) then
         write(6,*) ' indls gt KAXLSO in PREPAISO '
	 stop
      end if
      ff = 1.0d+00/dsqrt(2.0d+00)
      Maxmu = 2*nmax+2
      indks = 1
      do i=1,nmax
         do j=1,nmax
	    iksomu(i,j) = indks
            Minmu  = 2-Mod(i+j,2)
	    do imu=Minmu,Maxmu,2
	       aa1 =  dsqrt(dfloat(j-1))*DL1(i-1,j-2,imu-1)
	       aa2 = -dsqrt(dfloat(j  ))*DL1(i-1,j  ,imu-1)
	       aa3 = -dsqrt(dfloat(i-1))*DL1(i-2,j-1,imu-1)
	       aa4 =  dsqrt(dfloat(i  ))*DL1(i  ,j-1,imu-1)
	       aksomu(indks) = ff*(aa1+aa2+aa3+aa4)
	       indks = indks + 1
	    end do
	 end do
      end do
      if((indks-1).gt.KAXKSO) then
         write(6,*) ' indks gt KAXKSO in PREPAISO '
	 stop
      end if
	    
      icall = 061060
      	    
      end if
      
      imuze = 1 + isym
      imuzo = 2 - isym	    
      
      imuye = 1 + isym
      imuyo = 2 - isym	    
      
      mumax = 2*nxmax
      mumay = 2*nymax
      mumaz = 2*nzmax
      
      irop = 0
      irom = 0

      
      do ind=1,NSOXYZ
         do imuxyz=1,ndmu
            soxyz(imuxyz,ind) = 0.0d+0
	 end do
      end do
      
c+------------------------------------------------------------------+
c      isym=0        symmetric kappa 
c           1        skew-symmetric kappa
c      
c+------------------------------------------------------------------+
c             1 -> q'     2 -> q      2<=1 (lower triangle)
c
c     Kappa
c           q q'
c+------------------------------------------------------------------+
      do nx1=1,nxmax
         nym1 = MY(nx1)
         do nx2=1,nx1
            nym2 = MY(nx2)
            lx12 = nx1.eq.nx2
            lpx  = Mod(nx1+nx2,2).eq.0
	  
            do indy=1,NSOYZ
               do imuyz=1,iyzsrc
                  soyz(imuyz,indy) = 0.0d+0
	       end do
            end do
c
          do ny1 = 1,nym1
            nzi1e = nzie(nx1,ny1)
            nzi1o = nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
            if(lx12) nym2 = ny1
            do ny2=1,nym2
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              lxy12 = lx12.and.(ny1.eq.ny2)
              lpy  = Mod(ny1+ny2+isym,2).eq.0
c
              do indz=1,NSOZ
         	 do imuz=1,izsrc
         	    soz(imuz,indz) = 0.0d+0
         	 end do
              end do
c+---------------------------------------------------------------------+
c|                                                 Positive parity     |
c+---------------------------------------------------------------------+
c
c be careful not to move the if below, the iz1e=iy1e need to be done
c
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
              do nz1 =nzi1e,nzm1,2
                if(lxy12) nzm2=nz1
                do nz2 = nzi2e,nzm2,2
c
                 irop  = irop + 1
                 akan = GROP(irop+4+isym) ! kappa 
                 akap = GROP(irop+6+isym) ! kappa 
                 irop   = irop + 7
                 lxyz = lxy12.and.(nz2.eq.nz1)
                 if(lxyz) then
                    akan = 0.5d+00*akan
                    akap = 0.5d+00*akap
                 end if
		 
                 ils = ilsomu(nz2,nz1)
                 iks = iksomu(nz2,nz1)
                 indls= ils
                 indks= iks
c 		      +--------------------------------------------+
c 		      |   No  Field (0)   contribution  	   |
c 		      | 					   |
c 		      +--------------------------------------------+

c 		      +--------------------------------------------+
c 		      | 	Field (1)			   |
c 		      | 					   |
c 		      |    Tyz, Txz, Tyx	    Muz odd	   |
c 		      |    Tzx  		    Muz even	   |
c 		      +--------------------------------------------+
                 if(lpy.and.(.not.lpx)) then
                   imu = 0
                   do muz=imuzo,mumaz,2
                     imu = imu + 1
		     soz(imu,1) = soz(imu,1) + akan*alsomu(indls)
		     soz(imu,2) = soz(imu,2) + akap*alsomu(indls)
                     indls = indls + 1
                   end do  
                   imu = 0
                   do muz=imuze,mumaz,2
                     imu = imu + 1
                     soz (imu,3) =soz (imu,3)+ akan*aksomu(indks )/bz
                     soz (imu,4) =soz (imu,4)+ akap*aksomu(indks )/bz
                     indks = indks + 1
                   end do
c 			+--------------------------------------------+
c 			|	  Field (2)			     |
c 			|					     |
c 			|    Tyz, Txz, Txy	      Muz odd	     |
c 			|    Tzy		      Muz even       |
c 			+--------------------------------------------+
                   else if(.not.lpy.and.lpx) then
                   imu = 0
                   do muz=imuzo,mumaz,2
                     imu = imu + 1
		     soz(imu,5) = soz(imu,5) + akan*alsomu(indls)
		     soz(imu,6) = soz(imu,6) + akap*alsomu(indls)
                     indls = indls + 1
                   end do  
                   imu = 0
                   do muz=imuze,mumaz,2
                     imu = imu + 1
                     soz (imu,7) =soz (imu,7)+ akan*aksomu(indks )/bz
                     soz (imu,8) =soz (imu,8)+ akap*aksomu(indks )/bz
                     indks = indks + 1
                   end do
c 			+--------------------------------------------+
c 			|	  Field (3)			     |
c 			|					     |
c 			|    Tzy, Tzx		      Muz odd	     |
c 			|    Txy, Tyx		      Muz even       |
c 			+--------------------------------------------+
                   else if((.not.lpy).and.(.not.lpx)) then
                   imu = 0
                   do muz=imuzo,mumaz,2
                     imu = imu + 1
		     soz(imu, 9) = soz(imu, 9) + akan*aksomu(indks)/bz
		     soz(imu,10) = soz(imu,10) + akap*aksomu(indks)/bz
                     indks = indks + 1
                   end do  
                   imu = 0
                   do muz=imuze,mumaz,2
                     imu = imu + 1
                     soz (imu,11) =soz (imu,11)+ akan*alsomu(indls )
                     soz (imu,12) =soz (imu,12)+ akap*alsomu(indls )
                     indls = indls + 1
                   end do
c
                   end if ! Fields
c
                end do    ! nz2
              end do      ! nz1
          end if          ! lmze
c+--------------------------------------------------+
c|                                                  |
c|         Paridad negativa                         |
c|                                                  |
c+--------------------------------------------------+
c be careful not to move the if below, the iz1e=iy1e need to be done
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
	      
	      
              do nz1 = nzi1o,nzm1,2
                if(lxy12) nzm2=nz1
                do nz2 = nzi2o,nzm2,2
c
                 irom  = irom + 1
                 akan = GROM(irom+4+isym) ! kappa 
                 akap = GROM(irom+6+isym) ! kappa 
                 irom   = irom + 7
                 lxyz = lxy12.and.(nz2.eq.nz1)
                 if(lxyz) then
                    akan = 0.5d+00*akan
                    akap = 0.5d+00*akap
                 end if
                 ils = ilsomu(nz2,nz1)
                 iks = iksomu(nz2,nz1)
                 indls= ils
                 indks= iks
c      	              +--------------------------------------------+
c 	              |   No  Field (0)   contribution  	   |
c 	              | 					   |
c 	              +--------------------------------------------+

c 		      +--------------------------------------------+
c 		      |	       Field (1)			   |
c 		      | 					   |
c 		      |    Tyz, Txz, Tyx	    Muz odd	   |
c 		      |    Tzx  		    Muz even	   |
c 		      +--------------------------------------------+
c
                 if(lpy.and.(.not.lpx)) then
                   imu = 0
                   do muz=imuzo,mumaz,2
                     imu = imu + 1
		     soz(imu,1) = soz(imu,1) + akan*alsomu(indls)
		     soz(imu,2) = soz(imu,2) + akap*alsomu(indls)
                     indls = indls + 1
                   end do  
                   imu = 0
                   do muz=imuze,mumaz,2
                     imu = imu + 1
                     soz (imu,3) =soz (imu,3)+ akan*aksomu(indks )/bz
                     soz (imu,4) =soz (imu,4)+ akap*aksomu(indks )/bz
                     indks = indks + 1
                   end do
c 			+--------------------------------------------+
c 			|	  Field (2)			     |
c 			|					     |
c 			|    Tyz, Txz, Txy	      Muz odd	     |
c 			|    Tzy		      Muz even       |
c 			+--------------------------------------------+
                   else if(.not.lpy.and.lpx) then
                   imu = 0
                   do muz=imuzo,mumaz,2
                     imu = imu + 1
		     soz(imu,5) = soz(imu,5) + akan*alsomu(indls)
		     soz(imu,6) = soz(imu,6) + akap*alsomu(indls)
                     indls = indls + 1
                   end do  
                   imu = 0
                   do muz=imuze,mumaz,2
                     imu = imu + 1
                     soz (imu,7) =soz (imu,7)+ akan*aksomu(indks )/bz
                     soz (imu,8) =soz (imu,8)+ akap*aksomu(indks )/bz
                     indks = indks + 1
                   end do
c 			+--------------------------------------------+
c 			|	  Field (3)			     |
c 			|					     |
c 			|    Tzy, Tzx		      Muz odd	     |
c 			|    Txy, Tyx		      Muz even       |
c 			+--------------------------------------------+
                   else if((.not.lpy).and.(.not.lpx)) then
                   imu = 0
                   do muz=imuzo,mumaz,2
                     imu = imu + 1
		     soz(imu, 9) = soz(imu, 9) + akan*aksomu(indks)/bz
		     soz(imu,10) = soz(imu,10) + akap*aksomu(indks)/bz
                     indks = indks + 1
                   end do  
                   imu = 0
                   do muz=imuze,mumaz,2
                     imu = imu + 1
                     soz (imu,11) =soz (imu,11)+ akan*alsomu(indls )
                     soz (imu,12) =soz (imu,12)+ akap*alsomu(indls )
                     indls = indls + 1
                   end do
c
                   end if ! Fields
                end do     ! nz2
              end do       ! nz1
            end if ! lmzo
c
c 	    		    	   +------------------------+
c 	    		    	   |	  Y	L O O P     |
c 	    		    	   +------------------------+
c
            ils  = ilsomu(ny2,ny1)
            indls= ils
            iks  = iksomu(ny2,ny1)
            indks= iks
c
c
c 				+--------------------------+
c 				| N O	F I E L D  ( 0 )   |
c 				+--------------------------+

c 				   +---------------------+
c 				   |  F I E L D  ( 1 )   |
c 				   +---------------------+
c
            if(lpy.and.(.not.lpx)) then
	    	    
              fasy = Dfloat(1-2*Mod(Iabs(ny2-ny1+isym)/2,2)) ! (-)**(ny2-ny1)/2
c             
c                                        muy even 
c                                      
              iyzeo = 0
              iyzee = 0
c --------------------------------------  muz odd  Txz (1)
              do muy=imuye,mumay,2
                 sls= alsomu(indls)
                 indls = indls + 1
                 imuz = 0
                 do muz=imuzo,mumaz,2
                    imuz = imuz + 1
                    iyzeo = iyzeo + 1
		    soyz(iyzeo,1)=soyz(iyzeo,1)+fasy*soz(imuz,1)*sls
		    soyz(iyzeo,2)=soyz(iyzeo,2)+fasy*soz(imuz,2)*sls
                 end do  ! muz
c --------------------------------------  muz even Tzx (1)
                 imuz = 0
                 do muz=imuze,mumaz,2
                    imuz = imuz + 1
                    iyzee = iyzee + 1
		    soyz(iyzee,3)=soyz(iyzee,3)+fasy*soz(imuz,3)*sls
		    soyz(iyzee,4)=soyz(iyzee,4)+fasy*soz(imuz,4)*sls
                 end do  ! muz
              end do     ! muy
c             
c                                        muy odd 
c                                      
              iyzoo = 0
c --------------------------------------  muz odd  Tyz (1) Tyx (1)
              do muy=imuyo,mumay,2
                 sks= aksomu(indks)/by
                 indks = indks + 1
                 imuz = 0
                 do muz=imuzo,mumaz,2
                    imuz = imuz + 1
                    iyzoo = iyzoo + 1
		    soyz(iyzoo,5)=soyz(iyzoo,5)+fasy*soz(imuz,1)*sks
		    soyz(iyzoo,6)=soyz(iyzoo,6)+fasy*soz(imuz,2)*sks
                 end do  ! muz
	      end do     ! muy
c
c 				     +------------------+
c 				     | F I E L D  ( 2 ) |
c 				     +------------------+
c
            else if(.not.lpy.and.lpx) then
	                                                       
              fasy= Dfloat(1-2*Mod(Iabs(ny2-ny1+1-isym)/2,2)) !(-)**(ny2-ny1+1)/2
c	      
c                                        muy even 
c                                      
              iyzeo = 0
c --------------------------------------  muz odd  Tyz (2)
              do muy=imuye,mumay,2
                 sks= aksomu(indks)/by
                 indks = indks + 1
                 imuz = 0
                 do muz=imuzo,mumaz,2
                    imuz = imuz + 1
                    iyzeo = iyzeo + 1
		    soyz(iyzeo,7)=soyz(iyzeo,7)+fasy*soz(imuz,5)*sks
		    soyz(iyzeo,8)=soyz(iyzeo,8)+fasy*soz(imuz,6)*sks
                 end do  ! muz
              end do  	 ! muy
c             
c                                        muy odd 
c                                      
              iyzoo = 0
              iyzoe = 0
c --------------------------------------  muz odd  Txz (2) Txy (2)
              do muy=imuyo,mumay,2
                 sls= alsomu(indls)
                 indls = indls + 1
                 imuz = 0
                 do muz=imuzo,mumaz,2
                    imuz = imuz + 1
                    iyzoo = iyzoo + 1
		    soyz(iyzoo, 9)=soyz(iyzoo, 9)+fasy*soz(imuz,5)*sls
		    soyz(iyzoo,10)=soyz(iyzoo,10)+fasy*soz(imuz,6)*sls
                 end do  ! muz
c --------------------------------------  muz even Tzy (2)
                 imuz = 0
                 do muz=imuze,mumaz,2
                    imuz = imuz + 1
                    iyzoe = iyzoe + 1
		    soyz(iyzoe,11)=soyz(iyzoe,11)+fasy*soz(imuz,7)*sls
		    soyz(iyzoe,12)=soyz(iyzoe,12)+fasy*soz(imuz,8)*sls
                 end do  ! muz
              end do  	 ! muy
c
c 				   +--------------------------+
c 				   |	F I E L D  ( 3 )      |
c 				   +--------------------------+
c
            else if((.not.lpy).and.(.not.lpx)) then
	    
              fasy= Dfloat(1-2*Mod(Iabs(ny2-ny1+1-isym)/2,2)) !(-)**(ny2-ny1+1)/2
c	      
c                                        muy even 
c                                      
              iyzee = 0
              do muy=imuye,mumay,2
                 sks= aksomu(indks)/by
                 indks = indks + 1
c --------------------------------------  muz even Tyx (3)
                 imuz = 0
                 do muz=imuze,mumaz,2
                    imuz = imuz + 1
                    iyzee = iyzee + 1
		    soyz(iyzee,13)=soyz(iyzee,13)+fasy*soz(imuz,11)*sks
		    soyz(iyzee,14)=soyz(iyzee,14)+fasy*soz(imuz,12)*sks
                 end do  ! muz
              end do  	 ! muy
c             
c                                        muy odd 
c                                      
              iyzoo = 0
              iyzoe = 0
c --------------------------------------  muz odd  Tzy Tzx (3)
              do muy=imuyo,mumay,2
                 sls= alsomu(indls)
                 indls = indls + 1
                 imuz = 0
                 do muz=imuzo,mumaz,2
                    imuz = imuz + 1
                    iyzoo = iyzoo + 1
		    soyz(iyzoo,15)=soyz(iyzoo,15)+fasy*soz(imuz, 9)*sls
		    soyz(iyzoo,16)=soyz(iyzoo,16)+fasy*soz(imuz,10)*sls
                 end do  ! muz
c --------------------------------------  muz even Txy (3)
                 imuz = 0
                 do muz=imuze,mumaz,2
                    imuz = imuz + 1
                    iyzoe = iyzoe + 1
		    soyz(iyzoe,17)=soyz(iyzoe,17)+fasy*soz(imuz,11)*sls
		    soyz(iyzoe,18)=soyz(iyzoe,18)+fasy*soz(imuz,12)*sls
                 end do  ! muz
              end do  	 ! muy
	      
c ------------------------------------------------------------------
            end if    !    F I E L D S   Y
c ------------------------------------------------------------------
            end do    ! ny2
          end do      ! ny1
c
c 	               +------------------------+
c 	               |      X     L O O P	|
c 	               +------------------------+
c
            ils  = ilsomu(nx2,nx1)
            indls= ils
            iks  = iksomu(nx2,nx1)
            indks= iks
c
            fasx  = Dfloat(1-2*Mod(nx1+1,2))  !   (-)**nxq'
            fasx1 = Dfloat(1-2*Mod(nx1,2))    ! - (-)**nxq'
c  					 +------------+
c  --------------------------------------| nx+nx' par |
c  					 +------------+
            if(lpx) then
c                                                     mux even
              ixyzeeo = 0
              ixyzeoe = 0
              do mux =1,mumax,2
                 sls  = alsomu(indls)
                 indls = indls + 1
                 iyzeo = 0
c                                                     muy even
                 do muy=imuye,mumay,2
c                                                     muz odd        T yz (2)
                    do muz=imuzo,mumaz,2
                       iyzeo = iyzeo + 1
                       ixyzeeo = ixyzeeo + 1
	 soxyz(ixyzeeo,1)=soxyz(ixyzeeo,1)+soyz(iyzeo,7)*sls*fasx
	 soxyz(ixyzeeo,2)=soxyz(ixyzeeo,2)+soyz(iyzeo,8)*sls*fasx
		    end do  ! muz
		 end do     ! muy 
		 
                 iyzoe = 0
c                                                     muy odd
                 do muy=imuyo,mumay,2
c                                                     muz even      T zy (2)
                    do muz=imuze,mumaz,2
                       iyzoe = iyzoe + 1
                       ixyzeoe = ixyzeoe + 1
	 soxyz(ixyzeoe,3)=soxyz(ixyzeoe,3)+soyz(iyzoe,11)*sls*fasx
	 soxyz(ixyzeoe,4)=soxyz(ixyzeoe,4)+soyz(iyzoe,12)*sls*fasx
		    end do  ! muz
		 end do     ! muy
	      end do        ! mux
c
c                                                     mux odd
              ixyzooo = 0
              do mux =2,mumax,2
                 sks  = aksomu(indks)/bx
                 indks = indks + 1
                 iyzoo = 0
c                                                     muy odd
                 do muy=imuyo,mumay,2
c                                               muz odd    Txz  Txy
                    do muz=imuzo,mumaz,2
                       iyzoo = iyzoo + 1
                       ixyzooo = ixyzooo + 1
c     Txz (2)		       
	      soxyz(ixyzooo,5)=soxyz(ixyzooo,5)+soyz(iyzoo, 9)*sks*fasx
	      soxyz(ixyzooo,6)=soxyz(ixyzooo,6)+soyz(iyzoo,10)*sks*fasx
c     Txy (2)         
	      soxyz(ixyzooo,7)=soxyz(ixyzooo,7)+soyz(iyzoo, 9)*sks*fasx
	      soxyz(ixyzooo,8)=soxyz(ixyzooo,8)+soyz(iyzoo,10)*sks*fasx
		    end do  ! muz
		 end do     ! muy
	      end do 	    ! mux

c
c
c 				          +----------------+
c   --------------------------------------| nx+nx' impar   |
c 	                       	          +----------------+
            else if(.not.lpx) then
c                                                         
c                                                     mux even
              ixyzeeo = 0
              ixyzeoe = 0
              do mux =1,mumax,2
                 sks  = aksomu(indks)/bx
                 indks = indks + 1
                 iyzeo = 0
c                                                     muy even
                 do muy=imuye,mumay,2
c                                                     muz odd  T XZ (1)       
                    do muz=imuzo,mumaz,2
                       iyzeo = iyzeo + 1
                       ixyzeeo = ixyzeeo + 1
	soxyz(ixyzeeo, 9)=soxyz(ixyzeeo, 9)+soyz(iyzeo,1)*sks*fasx
	soxyz(ixyzeeo,10)=soxyz(ixyzeeo,10)+soyz(iyzeo,2)*sks*fasx
		    end do  ! muz
		 end do     ! muy
		 
                 iyzoe = 0
c                                                     muy odd
                 do muy=imuyo,mumay,2
c                                                     muz even T xy (3)     
                    do muz=imuze,mumaz,2
                       iyzoe = iyzoe + 1
                       ixyzeoe = ixyzeoe + 1
	  soxyz(ixyzeoe,11)=soxyz(ixyzeoe,11)+soyz(iyzoe,17)*sks
	  soxyz(ixyzeoe,12)=soxyz(ixyzeoe,12)+soyz(iyzoe,18)*sks
		    end do  ! muz
		 end do     ! muy
	      end do 	    ! mux
c
c                                                     mux odd
              ixyzooo = 0
              ixyzoee = 0
              do mux =2,mumax,2
                 sls  = alsomu(indls)
                 indls = indls + 1
                 iyzee = 0
c                                                     muy even
                 do muy=imuye,mumay,2
c                                       muz even T zx (1) T yx (3)   
                    do muz=imuze,mumaz,2
                       iyzee = iyzee + 1
                       ixyzoee = ixyzoee + 1
c   Tyx (3)
	    soxyz(ixyzoee,13)=soxyz(ixyzoee,13)+soyz(iyzee,13)*sls
	    soxyz(ixyzoee,14)=soxyz(ixyzoee,14)+soyz(iyzee,14)*sls
c   Tzx (1)
	    soxyz(ixyzoee,15)=soxyz(ixyzoee,15)+soyz(iyzee,3)*sls*fasx1
	    soxyz(ixyzoee,16)=soxyz(ixyzoee,16)+soyz(iyzee,4)*sls*fasx1
		    end do  ! muz
		 end do     ! muy
                 iyzoo = 0
c                                                     muy odd
                 do muy=imuyo,mumay,2
c                                                     muz odd    
                    do muz=imuzo,mumaz,2
                       iyzoo = iyzoo + 1
                       ixyzooo = ixyzooo + 1
c   T yz (1)
	     soxyz(ixyzooo,17)=soxyz(ixyzooo,17)+soyz(iyzoo,5)*sls*fasx
	     soxyz(ixyzooo,18)=soxyz(ixyzooo,18)+soyz(iyzoo,6)*sls*fasx
c   T yx (1)
	     soxyz(ixyzooo,19)=soxyz(ixyzooo,19)+soyz(iyzoo,5)*sls*fasx1
	     soxyz(ixyzooo,20)=soxyz(ixyzooo,20)+soyz(iyzoo,6)*sls*fasx1
c   T zy (3)
	     soxyz(ixyzooo,21)=soxyz(ixyzooo,21)+soyz(iyzoo,15)*sls
	     soxyz(ixyzooo,22)=soxyz(ixyzooo,22)+soyz(iyzoo,16)*sls
c   T zx (3)
	     soxyz(ixyzooo,23)=soxyz(ixyzooo,23)+soyz(iyzoo,15)*sls
	     soxyz(ixyzooo,24)=soxyz(ixyzooo,24)+soyz(iyzoo,16)*sls
		    end do  ! muz
		 end do     ! muy
	      end do 	    ! mux
c ------------------------------------------------------------------
            end if   ! lpx
c ------------------------------------------------------------------
        end do  ! nx2
      end do    ! nx1
      
      return
      end
c****************************************************************************
c                         P   A   I   S   O
c****************************************************************************
      Subroutine Paiso(isym,grop,grom,dpp,dnp,dpm,dnm,THE,Soxyz)
c
      Implicit real*8 (A-H,O-Z)
      Logical LMZE,LMZO,vmux
      Logical lx12,lxy12,lxyz,lpx,lpy
c
      Include 'COMDIM'
      Parameter (KAXT1=((imax*(imax+1))/2)*(imax+3))
      Parameter (KAXGSO= imax*imax*(imax+3))
c
      Dimension t1d(KAXT1),gsomu(KAXGSO)
c
      Dimension Soxyz(NDMU,24)
c
      Dimension THE(NDTHET,24)
      
      Dimension it1d(imax,imax),igsomu(imax,imax)
      
      Dimension itheeeo(izmax,izmax)
      Dimension itheeoe(izmax,izmax)
      Dimension itheoee(izmax,izmax)
      Dimension itheooo(izmax,izmax)
c
      Dimension dpp(nrop),dnp(nrop)
      Dimension dpm(nrom),dnm(nrom)
c +---------------------------------------------------------------------+
c |   In  Soxyz the meaning of the second index is                      |
c |=====================================================================|
c |                         nxmu     nymu     nzmu                      |
c |      1,2 .... T yz (2)  even     even      odd                      |
c |      3,4 .... T zy (2)  even      odd     even                      |
c |      5,6 .... T xz (2)   odd      odd      odd                      |
c |      7,8 .... T xy (2)   odd      odd      odd                      |
c +---------------------------------------------------------------------+
c |      9,10.... T xz (1)  even     even      odd                      |
c |      11,12... T xy (3)  even      odd     even                      |
c |      13,14... T yx (3)   odd     even     even                      |
c |      15,16... T zx (1)   odd     even     even                      |
c |      17,18... T yz (1)   odd      odd      odd                      |
c |      19,20... T yx (1)   odd      odd      odd                      |
c |      21,22... T zy (3)   odd      odd      odd                      |
c |      23,24... T zx (3)   odd      odd      odd                      |
c +---------------------------------------------------------------------+
c |      Odd indices are for neutrons while even ones are for protons   |
c +---------------------------------------------------------------------+
      Dimension GROP(nrop8),GROM(nrom8)
c ----------------------------------------------------------------
      save t1d,gsomu,it1d,igsomu,icall
c ----------------------------------------------------------------      
      Common /FLOCAR/nmax,nmax1,nxmax,nymax,nzmax,
     *  my(ixmax),mz(ixmax,iymax),nzie(ixmax,iymax),nzio(ixmax,iymax),
     *  lmze(ixmax,iymax),lmzo(ixmax,iymax)
c
      Common/DIMPDD/ irortmp,idmus0,jdmus0,idmus1,jdmus1,irorz,irory
c
      Common/DIMEN/kmax,kmax1,kxmax,kymax,kzmax,nwf2,
     * maxtz1,maxjz1,maxjx1,maxjy1,maxlso,maxdso,
     * ndmu,ndthet,nacou,nrop,nrom,nrop8,nrom8,nherm,nherm38,nlegn
c
      Common/GOGINT/amu(2),xw(2),xh(2),xb(2),xm(2),WLS,t3,alpha,x0,e2
      Common/OSCLEN/bx,by,bz
c
c      write(6,*) ' PAISO ISYM = ',isym
c-------------------------------------------------------------------
      if(icall.ne.061060) then
      call inigf
      Maxmu = 2*nmax+2
      indt1 = 1
      do i=1,nmax
         do j=i,nmax
	    it1d(i,j) = indt1
	    it1d(j,i) = indt1
            Minmu  = Mod(i+j,2)+1
	    do imu=Minmu,Maxmu,2
	       t1d(indt1) = DT1(i-1,j-1,imu-1)
	       indt1 = indt1 + 1
	    end do
	 end do
      end do
      if((indt1-1).gt.KAXT1) then
         write(6,*) ' indt1 gt KAXT1 in PAISO '
	 stop
      end if
      ff = 1.0d+00/dsqrt(2.0d+00)
      Maxmu = 2*nmax+2
      indgs = 1
      do i=1,nmax
         do j=1,nmax
	    igsomu(i,j) = indgs
            Minmu  = 2-Mod(i+j,2)
	    do imu=Minmu,Maxmu,2
	       aa1 =  dsqrt(dfloat(j-1))*DT1(i-1,j-2,imu-1)
	       aa2 = -dsqrt(dfloat(j  ))*DT1(i-1,j  ,imu-1)
	       aa3 = -dsqrt(dfloat(i-1))*DT1(i-2,j-1,imu-1)
	       aa4 =  dsqrt(dfloat(i  ))*DT1(i  ,j-1,imu-1)
	       gsomu(indgs) = ff*(aa1+aa2+aa3+aa4)
	       indgs = indgs + 1
	    end do
	 end do
      end do
      if((indgs-1).gt.KAXGSO) then
         write(6,*) ' indgs gt KAXGSO in PAISO '
	 stop
      end if
      	    
      icall = 061060
      	    
      end if	    
c
c
      call timeit(0,11,'  PREDIRSO      ')
      call Prepaiso(isym,Soxyz,GROP,GROM)
      call timeit(1,11,'  PREDIRSO      ')
c
      call timeit(0,14,'REST OF GDIR    ')
c
      pi32  = (4.0d+00* datan(1.0d+00))**(1.50d+00)
      fso  = WLS/(4.0d+00*pi32*bx*by*bz)
c
c     Optimization  (THETA)
c
      do nz1=1,nzmax
         do nz2=1,nzmax

            itheeeo(nz1,nz2) = 0
            itheeoe(nz1,nz2) = 0 
            itheoee(nz1,nz2) = 0
            itheooo(nz1,nz2) = 0
	    
	  end do
       end do


      imumax = 2*nxmax
      imumay = 2*nymax
      imumaz = 2*nzmax
      
      imumixe = 1
      imumixo = 2
      imumiye = 1+isym
      imumiyo = 2-isym
      imumize = 1+isym
      imumizo = 2-isym
c
c                                      nz1 + nz2   impar 
c
      iteteeo = 0
      iteteoe = 0
      itetoee = 0
      itetooo = 0
      
      do nz1=1,nzmax
         izmin = 1+MOD((nz1+isym),2) 
         do nz2=izmin,nzmax,2

            itheeeo(nz1,nz2) = iteteeo 
            itheeoe(nz1,nz2) = iteteoe 
            itheoee(nz1,nz2) = itetoee 
            itheooo(nz1,nz2) = itetooo 

	 
            imu = 1
	    indt1 = it1d(nz1,nz2)
            do imux=imumixe,imumax,2
               do imuy=imumiye,imumay,2
	       
                  sumn1     = 0.0d+00
                  sump1     = 0.0d+00
                  sumn2     = 0.0d+00
                  sump2     = 0.0d+00
                  ii = indt1
                  do imuz=imumizo,imumaz,2
                        sumn1 = sumn1 + t1d(ii)*soxyz(imu,1) 
                        sump1 = sump1 + t1d(ii)*soxyz(imu,2)
			
                        sumn2 = sumn2 + t1d(ii)*soxyz(imu,9)
                        sump2 = sump2 + t1d(ii)*soxyz(imu,10)
			imu = imu + 1 
			ii  = ii  + 1 
                  end do                               ! imuz
c                                       even even odd      T yz    T xz
                  iteteeo = iteteeo + 1
                  the (iteteeo,1) = sumn1 ! Tyz (2)
                  the (iteteeo,2) = sump1
                  the (iteteeo,3) = sumn2 ! Txz (1)
                  the (iteteeo,4) = sump2
               end do                                  ! imuy
            end do	                               ! imux
	    
            imu = 1
	    indt1 = it1d(nz1,nz2)
            do imux=imumixo,imumax,2
               do imuy=imumiyo,imumay,2
	       
                  sumn1     = 0.0d+00
                  sump1     = 0.0d+00
                  sumn2     = 0.0d+00
                  sump2     = 0.0d+00
                  sumn3     = 0.0d+00
                  sump3     = 0.0d+00
                  sumn4     = 0.0d+00
                  sump4     = 0.0d+00
                  ii = indt1
                  do imuz=imumizo,imumaz,2
                        sumn1 = sumn1 + t1d(ii)*soxyz(imu,5) ! Txz (2)
                        sump1 = sump1 + t1d(ii)*soxyz(imu,6)
			
                        sumn2 = sumn2 + t1d(ii)*soxyz(imu,21) ! T zy (3)
                        sump2 = sump2 + t1d(ii)*soxyz(imu,22)
			
                        sumn3 = sumn3 + t1d(ii)*soxyz(imu,17) ! Tyz (1)
                        sump3 = sump3 + t1d(ii)*soxyz(imu,18)
			
                        sumn4 = sumn4 + t1d(ii)*soxyz(imu,23) ! T zx (3)
                        sump4 = sump4 + t1d(ii)*soxyz(imu,24)
			
			imu = imu + 1 
			ii  = ii  + 1 
                  end do                               ! imuz
c                           odd odd odd  Txz (2) T zy(3) Tyz(1) Tzx (3)
                  itetooo = itetooo + 1
                  the (itetooo, 5) = sumn1  ! Txz (2)
                  the (itetooo, 6) = sump1
                  the (itetooo, 7) = sumn2  ! T zy (3)
                  the (itetooo, 8) = sump2
                  the (itetooo, 9) = sumn3  ! Tyz (1)
                  the (itetooo,10) = sump3
                  the (itetooo,11) = sumn4  ! T zx (3)
                  the (itetooo,12) = sump4
               end do                                  ! imuy
            end do                                     ! imux
	                                        
            imu = 1
	    indgs = igsomu(nz1,nz2)
            do imux=imumixe,imumax,2
               do imuy=imumiyo,imumay,2
	       
                  sumn1     = 0.0d+00
                  sump1     = 0.0d+00
		  
                  ii = indgs
                  do imuz=imumize,imumaz,2
                        sumn1 = sumn1 + gsomu(ii)*soxyz(imu,11)/bz ! Txy (3)
                        sump1 = sump1 + gsomu(ii)*soxyz(imu,12)/bz
						
			imu = imu + 1 
			ii  = ii  + 1 
                  end do                               ! imuz
c                           even odd even  Txy (3) Tyx (1)
                  iteteoe = iteteoe + 1
                  the (iteteoe,13) = sumn1 ! Txy (3)
                  the (iteteoe,14) = sump1
		  
               end do                                  ! imuy
            end do                                     ! imux
	                                        
            imu = 1
	    indgs = igsomu(nz1,nz2)
            do imux=imumixo,imumax,2
               do imuy=imumiye,imumay,2
	       
                  sumn1     = 0.0d+00
                  sump1     = 0.0d+00
                  ii = indgs
                  do imuz=imumize,imumaz,2
                        sumn1 = sumn1 + gsomu(ii)*soxyz(imu,13)/bz ! Tyx (3)
                        sump1 = sump1 + gsomu(ii)*soxyz(imu,14)/bz
						
			imu = imu + 1 
			ii  = ii  + 1 
                  end do                               ! imuz
c                            odd even even  Tyx (3)
                  itetoee = itetoee + 1
                  the (itetoee,15) = sumn1 ! Tyx (3)
                  the (itetoee,16) = sump1
               end do                                  ! imuy
            end do                                     ! imux
         end do                                        ! nz2
      end do                                           ! nz1

c
c                                      nz1 + nz2   par 
c
      
                                           
      iteteoe = 0
      itetoee = 0
      itetooo = 0
      
      do nz1=1,nzmax
         izmin = 1+MOD((nz1+1-isym),2) 
         do nz2=izmin,nzmax,2
	 

            itheeoe(nz1,nz2) = iteteoe
            itheoee(nz1,nz2) = itetoee
            itheooo(nz1,nz2) = itetooo

            imu = 1
	    indt1 = it1d(nz1,nz2)
            do imux=imumixe,imumax,2
               do imuy=imumiyo,imumay,2
	       
                  sumn1     = 0.0d+00
                  sump1     = 0.0d+00
                  ii = indt1
                  do imuz=imumize,imumaz,2
                        sumn1 = sumn1 + t1d(ii)*soxyz(imu,3) ! Tzy (2)
                        sump1 = sump1 + t1d(ii)*soxyz(imu,4)
			
			imu = imu + 1 
			ii  = ii  + 1 
                  end do                               ! imuz
c                                       even odd even      T zy
                  iteteoe = iteteoe + 1
                  the (iteteoe,17) = sumn1 ! Tzy (2)
                  the (iteteoe,18) = sump1
               end do                                  ! imuy
            end do                                     ! imux
	                                        
            imu = 1
	    indt1 = it1d(nz1,nz2)
            do imux=imumixo,imumax,2
               do imuy=imumiye,imumay,2
	       
                  sumn1     = 0.0d+00
                  sump1     = 0.0d+00
                  ii = indt1
                  do imuz=imumize,imumaz,2
                        sumn1 = sumn1 + t1d(ii)*soxyz(imu,15) ! Tzx (1)
                        sump1 = sump1 + t1d(ii)*soxyz(imu,16)
			
			imu = imu + 1 
			ii  = ii  + 1 
                  end do                               ! imuz
c                           odd even even  Tzx (1) 
                  itetoee = itetoee + 1
                  the (itetoee,19) = sumn1 ! Tzx (1)
                  the (itetoee,20) = sump1
               end do                                  ! imuy
            end do                                     ! imux
	                                        
            imu = 1
	    indgs = igsomu(nz1,nz2)
            do imux=imumixo,imumax,2
               do imuy=imumiyo,imumay,2
	       
                  sumn1     = 0.0d+00
                  sump1     = 0.0d+00
                  sumn2     = 0.0d+00
                  sump2     = 0.0d+00
		  
                  ii = indgs
                  do imuz=imumizo,imumaz,2
                        sumn1 = sumn1 + gsomu(ii)*soxyz(imu,7)/bz ! Txy (2)
                        sump1 = sump1 + gsomu(ii)*soxyz(imu,8)/bz
						
                        sumn2 = sumn2 + gsomu(ii)*soxyz(imu,19)/bz ! Tyx (1)
                        sump2 = sump2 + gsomu(ii)*soxyz(imu,20)/bz
			
			imu = imu + 1 
			ii  = ii  + 1 
                  end do                               ! imuz
c                           odd odd odd  Txy (2) 
                  itetooo = itetooo + 1
                  the (itetooo,21) = sumn1 ! Txy (2)
                  the (itetooo,22) = sump1
		  
                  the (itetooo,23) = sumn2 ! Tyx (1)
                  the (itetooo,24) = sump2
               end do                                  ! imuy
            end do                                     ! imux
	                                        
         end do                                        ! nz2
      end do                                           ! nz1
      
c+-----------------------------------------------------------------+
c|                        T H E                                    |
c+-----------------------------------------------------------------+
c|                     nz1+nz2 odd  (Fields (1) and (2))           |
c+-----------------------------------------------------------------+
c|       1, 2 ......... Tyz(2)          mux even muy even muz odd  |
c|       3, 4 ......... Txz(1)                                     |
c|       5, 6 ......... Txz(2)          mux odd  muy odd  muz odd  |
c|       7, 8 ......... Tzy(3)                                     |
c|       9,10 ......... Tyz(1)                                     |
c|      11,12 ......... Tzx(3)                                     |
c|      13,14 ......... Txy(3)          mux even muy odd  muz even |
c|      15,16 ......... Tyx(3)          mux odd  muy even muz even |
c+-----------------------------------------------------------------+
c|                     nz1+nz2 even (Fields (0) and (3))           |
c+-----------------------------------------------------------------+
c|      17,18 ......... Tzy(2)          mux even muy odd  muz even |
c|      19,20 ......... Tzx(1)          mux odd  muy even muz even |
c|      21,22 ......... Txy(2)          mux odd  muy odd  muz odd  |
c|      23,24 ......... Tyx(1)                                     |
c+-----------------------------------------------------------------+
c
c
c
c
c
c      
c
c    +---------------------------------------------------------+
c    |   Starts calculation of the pairing field               |
c    +---------------------------------------------------------+
c

      irop = 0
      irom = 0
      jrop = 0
      jrom = 0
      
      deltap = 0.0d+00
      deltan = 0.0d+00
      
      do nx1=1,nxmax
        nym1 = MY(nx1)
        do nx2=1,nx1
          fx2  = 1-2*Mod((nx2+1),2)                    ! (-)**x2
          lpx  = Mod((nx1+nx2),2).eq.0                 ! type of field
          nym2 = MY(nx2)
          lx12 = nx1.eq.nx2
          igx  = igsomu(nx1,nx2)                        ! G index
          itx  = it1d  (nx1,nx2)                        ! T index
          do ny1 = 1,nym1
            nzi1e= nzie(nx1,ny1)
            nzi1o= nzio(nx1,ny1)
            nzm1 = MZ(nx1,ny1)
c
            if(lx12) nym2 = ny1
            do ny2=1,nym2
	      lpy  = Mod((ny1+ny2+isym),2).eq.0
              nzi2e= nzie(nx2,ny2)
              nzi2o= nzio(nx2,ny2)
              nzm2 = MZ(nx2,ny2)
              igy = igsomu(ny1,ny2)                      ! G index
              ity = it1d  (ny1,ny2)                      ! T index
              fy1 = 1-2*Mod(Iabs(ny2-ny1-isym)/2,2)   ! i**(ny2-ny1  -f(-))
              fy2 = 1-2*Mod(Iabs(ny2-ny1+isym)/2,2)   ! i**(ny2-ny1+1-f(+))
              fy3 = 1-2*Mod(Iabs(ny2-ny1+1-isym)/2,2) ! i**(ny2-ny1+1-f(-))
              fy4 = 1-2*Mod(Iabs(ny2-ny1-1+isym)/2,2) ! i**(ny2-ny1  -f(+))
              lxy12 = lx12.and.(ny1.eq.ny2)
c
c---->    Positive parity
c
c
c be careful not to move the if below, the iz1e=iy1e need to be done
c
              if(LMZE(nx1,ny1).and.LMZE(nx2,ny2)) then
c
              do nz1 =nzi1e,nzm1,2
                if(lxy12) nzm2=nz1
                do nz2 = nzi2e,nzm2,2
c
c                                             F I E L D  (1)
c
            sumnx = 0.0d+00
            sumpx = 0.0d+00
            sumny = 0.0d+00
            sumpy = 0.0d+00
            sumnz = 0.0d+00
            sumpz = 0.0d+00


            if(.not.lpx.and.lpy) then       ! FIELD (1)
	    
            iteteeo = itheeeo(nz1,nz2)
            itetooo = itheooo(nz1,nz2)
            itetoee = itheoee(nz1,nz2)
	    	    
	    inddsx = igx
            do imux=imumixe,imumax,2
               indlsy = ity
               do imuy=imumiye,imumay,2
                  iteteeo = iteteeo + 1
		  aa = gsomu(inddsx)*t1d(indlsy)/bx
c  Tyz(2)		  
                  sumnx = sumnx + aa*the(iteteeo,1)		  
                  sumpx = sumpx + aa*the(iteteeo,2)		  
                  indlsy = indlsy + 1
               end do                                  ! imuy
               inddsx = inddsx + 1
            end do                                     ! imux
	    
	    sumnx = 2.0d+00*sumnx * fx2 * fy2
	    sumpx = 2.0d+00*sumpx * fx2 * fy2
	    
	    indlsx = itx
            do imux=imumixo,imumax,2
               inddsy = igy
               do imuy=imumiyo,imumay,2
                  itetooo = itetooo + 1
		  aa = t1d(indlsx)*gsomu(inddsy)/by
            sumny = sumny + aa*(the(itetooo,11)-the(itetooo,5))
            sumpy = sumpy + aa*(the(itetooo,12)-the(itetooo,6))
                  inddsy = inddsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	    
	    sumny = 2.0d+00*sumny * fx2 * fy2
	    sumpy = 2.0d+00*sumpy * fx2 * fy2
	    
	    indlsx = itx
            do imux=imumixo,imumax,2
               indlsy = ity
               do imuy=imumiye,imumay,2
                  itetoee = itetoee + 1
		  aa = t1d(indlsx)*t1d(indlsy)
                  sumnz = sumnz - aa*the(itetoee,15)		  
                  sumpz = sumpz - aa*the(itetoee,16)		  
                  indlsy = indlsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	                                        
	    sumnz = 2.0d+00 * sumnz * fx2 * fy2
	    sumpz = 2.0d+00 * sumpz * fx2 * fy2
c
c                                             F I E L D  (2)
c
            else if(lpx.and..not.lpy) then       ! FIELD (2)
	    
            iteteeo = itheeeo(nz1,nz2)
            itetooo = itheooo(nz1,nz2)
            iteteoe = itheeoe(nz1,nz2)
	    
	    inddsx = igx
            do imux=imumixo,imumax,2
               indlsy = ity
               do imuy=imumiyo,imumay,2
                  itetooo = itetooo + 1
		  aa = gsomu(inddsx)*t1d(indlsy)/bx
                  sumnx = sumnx + aa*(the(itetooo, 9)*fy3
     &	                 -fy4*the(itetooo,7))
                  sumpx = sumpx + aa*(the(itetooo,10)*fy3
     &		         -fy4*the(itetooo,8))
                  indlsy = indlsy + 1
               end do                                  ! imuy
               inddsx = inddsx + 1
            end do                                     ! imux
	                                         
            sumnx = 2.0d+00*sumnx * fx2
            sumpx = 2.0d+00*sumpx * fx2
	    
	    indlsx = itx
            do imux=imumixe,imumax,2
               inddsy = igy
               do imuy=imumiye,imumay,2
                  iteteeo = iteteeo + 1
		  aa = t1d(indlsx)*gsomu(inddsy)/by
                  sumny = sumny - aa*the(iteteeo,3)
                  sumpy = sumpy - aa*the(iteteeo,4)
                  inddsy = inddsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	                                        
            sumny = 2.0d+00*sumny * fx2 * fy3
            sumpy = 2.0d+00*sumpy * fx2 * fy3
	    
	    indlsx = itx
            do imux=imumixe,imumax,2
               indlsy = ity
               do imuy=imumiyo,imumay,2
                  iteteoe = iteteoe + 1
		  aa = t1d(indlsx)*t1d(indlsy)
                  sumnz = sumnz + aa*the(iteteoe,13)		  
                  sumpz = sumpz + aa*the(iteteoe,14)		  
                  indlsy = indlsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	                                        
            sumnz = sumnz * fx2 * fy4 * 2.0d+00
            sumpz = sumpz * fx2 * fy4 * 2.0d+00 
	    
c
c                                             F I E L D  (3)
c
            else if(.not.lpx.and..not.lpy) then       ! FIELD (3)
	    
            itetoee = itheoee(nz1,nz2)
            itetooo = itheooo(nz1,nz2)
            iteteoe = itheeoe(nz1,nz2)
	    
	    inddsx = igx
            do imux=imumixe,imumax,2
               indlsy = ity
               do imuy=imumiyo,imumay,2
                  iteteoe = iteteoe + 1
		  aa = gsomu(inddsx)*t1d(indlsy)/bx
                  sumnx = sumnx - aa*the(iteteoe,17)
                  sumpx = sumpx - aa*the(iteteoe,18)
                  indlsy = indlsy + 1
               end do                                  ! imuy
               inddsx = inddsx + 1
            end do                                     ! imux
	    
            sumnx = -sumnx  * fy4 * 2.0d+00
            sumpx = -sumpx  * fy4 * 2.0d+00
	    
	    indlsx = itx
            do imux=imumixo,imumax,2
               inddsy = igy
               do imuy=imumiye,imumay,2
                  itetoee = itetoee + 1
		  aa = t1d(indlsx)*gsomu(inddsy)/by
                  sumny = sumny + aa*the(itetoee,19)
                  sumpy = sumpy + aa*the(itetoee,20)
                  inddsy = inddsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	    
            sumny = -sumny  * fy3 * 2.0d+00
            sumpy = -sumpy  * fy3 * 2.0d+00
	    
	    indlsx = itx
            do imux=imumixo,imumax,2
               indlsy = ity
               do imuy=imumiyo,imumay,2
                  itetooo = itetooo + 1
		  aa = t1d(indlsx)*t1d(indlsy)
          sumnz = sumnz + aa*(the(itetooo,21)*fy4-the(itetooo,23)*fy3)
          sumpz = sumpz + aa*(the(itetooo,22)*fy4-the(itetooo,24)*fy3)
                  indlsy = indlsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux 
	                                        
            sumnz = -sumnz  * 2.0d+00
            sumpz = -sumpz  * 2.0d+00
c
           end if  ! fields
c
                irop = irop + 1
	        dpp(irop) = dpp(irop) + fso*(sumpx + sumpy + sumpz ) 
	        dnp(irop) = dnp(irop) + fso*(sumnx + sumny + sumnz )
		
                jrop  = jrop + 1
                akan = GROP(jrop+4+isym) ! kappa 
                akap = GROP(jrop+6+isym) ! kappa 
                jrop   = jrop + 7
                lxyz = lxy12.and.(nz2.eq.nz1)
                if(lxyz) then
                    akan = 0.5d+00*akan
                    akap = 0.5d+00*akap
                end if
		
		deltap = deltap + fso*(sumpx + sumpy + sumpz )*akap
		deltan = deltan + fso*(sumnx + sumny + sumnz )*akan
		 
                end do   ! nz2
              end do     ! nz1
              end if    ! lmze
c
c                                           ---->    Negative parity
c
c
c be careful not to move the if below, the iz1o=iy1o need to be done
c
              if(LMZO(nx1,ny1).and.LMZO(nx2,ny2)) then
c
              do nz1 =nzi1o,nzm1,2
                if(lxy12) nzm2=nz1
                do nz2 = nzi2o,nzm2,2
		
c
            sumnx = 0.0d+00
            sumpx = 0.0d+00
            sumny = 0.0d+00
            sumpy = 0.0d+00
            sumnz = 0.0d+00
            sumpz = 0.0d+00


            if(.not.lpx.and.lpy) then       ! FIELD (1)
	    
            iteteeo = itheeeo(nz1,nz2)
            itetooo = itheooo(nz1,nz2)
            itetoee = itheoee(nz1,nz2)
	    	    
	    inddsx = igx
            do imux=imumixe,imumax,2
               indlsy = ity
               do imuy=imumiye,imumay,2
                  iteteeo = iteteeo + 1
		  aa = gsomu(inddsx)*t1d(indlsy)/bx
c  Tyz(2)		  
                  sumnx = sumnx + aa*the(iteteeo,1)		  
                  sumpx = sumpx + aa*the(iteteeo,2)		  
                  indlsy = indlsy + 1
               end do                                  ! imuy
               inddsx = inddsx + 1
            end do                                     ! imux
	    
	    sumnx = 2.0d+00*sumnx * fx2 * fy2
	    sumpx = 2.0d+00*sumpx * fx2 * fy2
	    
	    indlsx = itx
            do imux=imumixo,imumax,2
               inddsy = igy
               do imuy=imumiyo,imumay,2
                  itetooo = itetooo + 1
		  aa = t1d(indlsx)*gsomu(inddsy)/by
            sumny = sumny + aa*(the(itetooo,11)-the(itetooo,5))
            sumpy = sumpy + aa*(the(itetooo,12)-the(itetooo,6))
                  inddsy = inddsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	    
	    sumny = 2.0d+00*sumny * fx2 * fy2
	    sumpy = 2.0d+00*sumpy * fx2 * fy2
	    
	    indlsx = itx
            do imux=imumixo,imumax,2
               indlsy = ity
               do imuy=imumiye,imumay,2
                  itetoee = itetoee + 1
		  aa = t1d(indlsx)*t1d(indlsy)
                  sumnz = sumnz - aa*the(itetoee,15)		  
                  sumpz = sumpz - aa*the(itetoee,16)		  
                  indlsy = indlsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	                                        
	    sumnz = 2.0d+00 * sumnz * fx2 * fy2
	    sumpz = 2.0d+00 * sumpz * fx2 * fy2
c
c                                             F I E L D  (2)
c
            else if(lpx.and..not.lpy) then       ! FIELD (2)
	    
            iteteeo = itheeeo(nz1,nz2)
            itetooo = itheooo(nz1,nz2)
            iteteoe = itheeoe(nz1,nz2)
	    
	    inddsx = igx
            do imux=imumixo,imumax,2
               indlsy = ity
               do imuy=imumiyo,imumay,2
                  itetooo = itetooo + 1
		  aa = gsomu(inddsx)*t1d(indlsy)/bx
                  sumnx = sumnx + aa*(the(itetooo, 9)*fy3
     &	                 -fy4*the(itetooo,7))
                  sumpx = sumpx + aa*(the(itetooo,10)*fy3
     &		         -fy4*the(itetooo,8))
                  indlsy = indlsy + 1
               end do                                  ! imuy
               inddsx = inddsx + 1
            end do                                     ! imux
	                                         
            sumnx = 2.0d+00*sumnx * fx2
            sumpx = 2.0d+00*sumpx * fx2
	    
	    indlsx = itx
            do imux=imumixe,imumax,2
               inddsy = igy
               do imuy=imumiye,imumay,2
                  iteteeo = iteteeo + 1
		  aa = t1d(indlsx)*gsomu(inddsy)/by
                  sumny = sumny - aa*the(iteteeo,3)
                  sumpy = sumpy - aa*the(iteteeo,4)
                  inddsy = inddsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	                                        
            sumny = 2.0d+00*sumny * fx2 * fy3
            sumpy = 2.0d+00*sumpy * fx2 * fy3
	    
	    indlsx = itx
            do imux=imumixe,imumax,2
               indlsy = ity
               do imuy=imumiyo,imumay,2
                  iteteoe = iteteoe + 1
		  aa = t1d(indlsx)*t1d(indlsy)
                  sumnz = sumnz + aa*the(iteteoe,13)		  
                  sumpz = sumpz + aa*the(iteteoe,14)		  
                  indlsy = indlsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	                                        
            sumnz = sumnz * fx2 * fy4 * 2.0d+00
            sumpz = sumpz * fx2 * fy4 * 2.0d+00 
	    
c
c                                             F I E L D  (3)
c
            else if(.not.lpx.and..not.lpy) then       ! FIELD (3)
	    
            itetoee = itheoee(nz1,nz2)
            itetooo = itheooo(nz1,nz2)
            iteteoe = itheeoe(nz1,nz2)
	    
	    inddsx = igx
            do imux=imumixe,imumax,2
               indlsy = ity
               do imuy=imumiyo,imumay,2
                  iteteoe = iteteoe + 1
		  aa = gsomu(inddsx)*t1d(indlsy)/bx
                  sumnx = sumnx - aa*the(iteteoe,17)
                  sumpx = sumpx - aa*the(iteteoe,18)
                  indlsy = indlsy + 1
               end do                                  ! imuy
               inddsx = inddsx + 1
            end do                                     ! imux
	    
            sumnx = -sumnx  * fy4 * 2.0d+00
            sumpx = -sumpx  * fy4 * 2.0d+00
	    
	    indlsx = itx
            do imux=imumixo,imumax,2
               inddsy = igy
               do imuy=imumiye,imumay,2
                  itetoee = itetoee + 1
		  aa = t1d(indlsx)*gsomu(inddsy)/by
                  sumny = sumny + aa*the(itetoee,19)
                  sumpy = sumpy + aa*the(itetoee,20)
                  inddsy = inddsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux
	    
            sumny = -sumny  * fy3 * 2.0d+00
            sumpy = -sumpy  * fy3 * 2.0d+00
	    
	    indlsx = itx
            do imux=imumixo,imumax,2
               indlsy = ity
               do imuy=imumiyo,imumay,2
                  itetooo = itetooo + 1
		  aa = t1d(indlsx)*t1d(indlsy)
          sumnz = sumnz + aa*(the(itetooo,21)*fy4-the(itetooo,23)*fy3)
          sumpz = sumpz + aa*(the(itetooo,22)*fy4-the(itetooo,24)*fy3)
                  indlsy = indlsy + 1
               end do                                  ! imuy
               indlsx = indlsx + 1
            end do                                     ! imux 
	                                        
            sumnz = -sumnz  * 2.0d+00
            sumpz = -sumpz  * 2.0d+00
c
           end if  ! fields
c
                irom = irom + 1
	        dpm(irom) = dpm(irom) + fso*(sumpx + sumpy + sumpz ) 
	        dnm(irom) = dnm(irom) + fso*(sumnx + sumny + sumnz )
	   
                jrom  = jrom + 1
                akan = GROM(jrom+4+isym) ! kappa 
                akap = GROM(jrom+6+isym) ! kappa 
                jrom   = jrom + 7
                lxyz = lxy12.and.(nz2.eq.nz1)
                if(lxyz) then
                    akan = 0.5d+00*akan
                    akap = 0.5d+00*akap
                end if
		
		deltap = deltap + fso*(sumpx + sumpy + sumpz )*akap
		deltan = deltan + fso*(sumnx + sumny + sumnz )*akan
		 

                end do  ! nz2
              end do    ! nz1
              end if    ! lmzo
            end do  ! ny2
          end do    ! ny1
        end do      ! nx2
      end do        ! nx1
      
      fsym = dfloat( (-1)**isym )
      deltap = fsym*deltap
      deltan = fsym*deltan
      write(6,*) ' Spin-orbit pairing ISYM ',isym,deltap,deltan
      
      call timeit(1,14,'REST OF GDIR    ')
      return
      end
