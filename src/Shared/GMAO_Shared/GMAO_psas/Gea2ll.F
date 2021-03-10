	subroutine Gea2ll(
     &	    nG,          nkvG,  ktmapG,  kpmapG,
     &	    nvecs,
     &	         mxG,   mxkvG,    Gvec,
     &	         idm,     jdm,     kdm,
     &		         HHai,    uuai,    vvai,    qqai,
     &	                slpai,    usai,    vsai )

	implicit none

c   Input parameters
c  ==================

c     ..Define a G-vec type structure
	integer     nG	! the first dimension of Gvec
	  integer   nkvG	! the second dimension of Gvec
	  integer ktmapG(nkvG)	! variable type table of 2-d fields
	  integer kpmapG(nkvG)	! vertical level table of 2-d fields

	integer nvecs	! number of G-vector in Gvec
	  integer   mxG	! the size of the first dimenstion of Gvec
	  integer mxkvG	! the size of the second dimension of Gvec
	  real     Gvec(mxG,mxkvG,nvecs)	! the vectors

c     ..Define a 3-d arrye type structure for all analysis variables
	integer    idm	! number of longitude grid points of ??ai
	integer    jdm	! number of latitude grid points
	integer	   kdm	! number of vertical levels for upper-air

c   Output parameters, still part of the output data structure
c  ===================
	  real	  HHai(idm,jdm,kdm,nvecs)	! height
	  real	  uuai(idm,jdm,kdm,nvecs)	! u-wind, upper-air
	  real	  vvai(idm,jdm,kdm,nvecs)	! v-wind, upper-air
	  real	  qqai(idm,jdm,kdm,nvecs)	! mixing ratio
	  real	 slpai(idm,jdm,nvecs)	! sea-level pressure
	  real	  usai(idm,jdm,nvecs)	! sea-level u-wind
	  real	  vsai(idm,jdm,nvecs)	! sea-level v-wind

c   Local variables
c  =================
	integer ivec	! index of vectors
	integer kv	! index of 2-d fields
	integer kp	! level index in a 3-d grid

c   Pparameters
c  =============
	include "ktmax.h"	! define ktHH, ktuu, etc.

c	==================================
c	..Loop over all level in the G-vec

	do ivec=1,nvecs
	  do kv=1,nkvG
	    kp=kpmapG(kv)	! where is this level in a 3-d grid

c	  ..All upper-air variables

	    if(ktmapG(kv).eq.ktHH) then
	      call ea2ll( HHai(1,1,kp,ivec),idm,jdm,Gvec(1,kv,ivec),nG)
	    elseif(ktmapG(kv).eq.ktuu) then
	      call ea2ll( uuai(1,1,kp,ivec),idm,jdm,Gvec(1,kv,ivec),nG)
	    elseif(ktmapG(kv).eq.ktvv) then
	      call ea2ll( vvai(1,1,kp,ivec),idm,jdm,Gvec(1,kv,ivec),nG)
	    elseif(ktmapG(kv).eq.ktqq) then
	      call ea2ll( qqai(1,1,kp,ivec),idm,jdm,Gvec(1,kv,ivec),nG)

c	  ..All surface variables

	    elseif(ktmapG(kv).eq.ktslp) then
	      call ea2ll(slpai(1,1,ivec),idm,jdm,Gvec(1,kv,ivec),nG)
	    elseif(ktmapG(kv).eq.ktus) then
	      call ea2ll( usai(1,1,ivec),idm,jdm,Gvec(1,kv,ivec),nG)
	    elseif(ktmapG(kv).eq.ktvs) then
	      call ea2ll( vsai(1,1,ivec),idm,jdm,Gvec(1,kv,ivec),nG)
	    endif


	  end do
	end do
	end
c.
