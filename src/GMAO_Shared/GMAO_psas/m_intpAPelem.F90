!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_intpAPelem - Interpolations from a GEOS A-pressure grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_intpAPelem
      implicit none
      private	! except

      public :: intpAP

      interface intpAP; module procedure	&
	dintpvec_,	&	! scatterd data attributes for vectors
	dintpscl_,	&	! scatterd data attributes for scalars
	gintpvec_,	&	! gridded data output for vectors
	gintpscl_		! gridded data output for scalars
      end interface

! !REVISION HISTORY:
! 	14Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_intpAP'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: dintpscl_ - interpolate a A-grid scalar field
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine dintpscl_(ndat,alat,alon,alev,adat,		&
	nlon,nlat,nlev,pres,grd, stat		)

      use m_die,      only : die,perr
      use m_mall,     only : mall_ison,mall_mci,mall_mco
      implicit none

		! About the insitu scalar data for output

      integer,intent(in) :: ndat	! number of insitu data
      real,dimension(:), intent(in)  :: alat	! latitudes
      real,dimension(:), intent(in)  :: alon	! longitudes
      real,dimension(:), intent(in)  :: alev	! pres. levels
      real,dimension(:), intent(out) :: adat	! insitu data

		! About the grid point scalar data as the input

      integer,intent(in) :: nlat	! number of A-grid lat.
      integer,intent(in) :: nlon	! number of longitudes
      integer,intent(in) :: nlev	! number of levels
      real,dimension(:),    intent(in) :: pres	! pres. levels.
      real,dimension(:,:,:),intent(in) :: grd	! a grided scalar field

		! Optional arguments

      integer, optional, intent(out):: stat	! Status code

! !REVISION HISTORY:
! 	19Mar97 - Jing Guo <guo@dao> - revised from intp_ap_vect_()
!EOP ___________________________________________________________________
!#define _DEBUG	!

	! Parameters

  character(len=*),parameter :: myname_=myname//'::dintpscl_'
  logical,parameter	     :: nearest=.true.

  real, parameter :: s_lat = -90., n_lat= +90.
  real, parameter :: w_lon =-180., e_lon=+180.

	! Local variables

      integer :: ierr
      integer :: i,i1
      integer :: j,j1
      integer :: k,k1
      integer :: l,ix

      real :: x
      real :: dlon,dlat
      real :: wtlon,wtlat,wtlev

      real :: am,ap,akp

	! Temporary working space

      real,              allocatable :: wght(:)
      integer,           allocatable :: glev(:)

	! Basic sanity checks

  if(present(stat)) stat=0

  if(nlon < 2 .or. nlat < 2 .or. nlev < 1) then
    if(nlon < 2) call perr(myname_,'invalid nlon',nlon)
    if(nlat < 2) call perr(myname_,'invalid nlat',nlat)
    if(nlev < 1) call perr(myname_,'invalid nlev',nlev)
    if(.not.present(stat)) call die(myname_)
    stat=1
    return
  endif

	! For a GEOS-DAS A grid

		! Again, the input values for alat(:) are expected to
		! be within [-90.,+90.]; for alon(:) are expected to be
		! within [-180.,360.)

        dlon = (e_lon-w_lon)/nlon
        dlat = (n_lat-s_lat)/(nlat-1)

        allocate(glev(ndat),wght(ndat),stat=ierr)
        if(ierr /= 0) then
	  call perr(myname_,'allocate()',ierr)
	  if(.not.present(stat)) call die(myname_)
	  stat=ierr
	  return
        endif
	if(mall_ison()) then
	  call mall_mci(glev,myname)
	  call mall_mci(wght,myname)
	endif

		! Indexing to the vertical level grid values

        call slogtab(.not.nearest,nlev,pres,	&
              ndat,alev,glev,wght		)
       
        do l = 1, ndat
     
		! Longitudinal direction indexing

          x = (alon(l) - w_lon)/dlon
	  ix = int(x)
	  if(ix > x) ix=ix-1	! a faster ix = floor(x)
          wtlon =  x - ix

	 	! ix = modulo(ix,nlon), works only upto 360.+180.
		! ix is in [0,nlon-1].

          if(ix >= nlon) ix=ix-nlon

		! the indices at j, i and i1

	  i =ix+1
	  i1=i +1
	  if(i1 > nlon) i1=1	! Circular indexing
             
		! Latitudinal direction indexing

          x = (alat(l) - s_lat)/dlat
	  ix = int(x)
	  if(ix > x) ix=ix-1	! a faster ix = floor(x)
          wtlat = x - ix

	  j =ix+1
	  if(j >= nlat) then
	    j=nlat-1
	    wtlat=1.
	  endif
	  j1=j +1

        k = glev(l)
	wtlev=wght(l)
	k1=k+1

		! Bilinear interpolation of a scalar field at the
		! lower level.

	am =grd(i ,j ,k )+wtlon*(grd(i1,j ,k )-grd(i ,j ,k ))
	ap =grd(i ,j1,k )+wtlon*(grd(i1,j1,k )-grd(i ,j1,k ))
	adat(l)=am+wtlat*(ap-am)

	if(k1 <= nlev) then

		! Bilinear interpolation of a scalar field at the
		! upper level.

	  am =grd(i ,j ,k1)+wtlon*(grd(i1,j ,k1)-grd(i ,j ,k1))
	  ap =grd(i ,j1,k1)+wtlon*(grd(i1,j1,k1)-grd(i ,j1,k1))
	  akp=am+wtlat*(ap-am)

	  adat(l)=adat(l)+wtlev*(akp-adat(l))
	endif

      enddo   ! ndat

	if(mall_ison()) then
	  call mall_mco(glev,myname)
	  call mall_mco(wght,myname)
	endif
      deallocate(glev,wght,stat=ierr)
      if(ierr /= 0) then
	call perr(myname_,'deallocate()',ierr)
	if(.not.present(stat)) call die(myname_)
	stat=ierr
	return
      endif

      if(present(stat)) stat=0
end subroutine dintpscl_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: gintpscl_ - 2-d grid interpolate a 3-d A-grid scalar field
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gintpscl_(nlon_ot,nlat_ot,plev_ot,glev,	&
	nlon,nlat,nlev,pres,grd, stat		)

      use m_die,      only : die,perr
      implicit none

		! About the insitu scalar data for output

      integer,intent(in) :: nlon_ot
      integer,intent(in) :: nlat_ot
      real,   intent(in) :: plev_ot
      real,dimension(:,:),intent(out) :: glev

		! About the grid point scalar data as the input

      integer,intent(in) :: nlon	! number of longitudes
      integer,intent(in) :: nlat	! number of A-grid lat.
      integer,intent(in) :: nlev	! number of levels
      real,dimension(:),    intent(in) :: pres	! pres. levels.
      real,dimension(:,:,:),intent(in) :: grd	! a grided scalar field

		! Optional arguments

      integer, optional, intent(out):: stat	! Status code

! !REVISION HISTORY:
! 	14Jan99 - Jing Guo <guo@dao> -
!		. revised from intp_ap_.F90
!		. wrote a new prolog
!EOP ___________________________________________________________________

!#define _DEBUG	!

	! Parameters

  character(len=*),parameter :: myname_=myname//'::gintpscl_'
  logical,parameter	     :: NEAREST=.true.

  real, parameter :: s_lat_ot = -90., n_lat_ot= +90.
  real, parameter :: w_lon_ot =-180., e_lon_ot=+180.

  real, parameter :: s_lat = -90., n_lat= +90.
  real, parameter :: w_lon =-180., e_lon=+180.

	! Local variables

      integer :: ierr
      integer :: i,i1
      integer :: j,j1
      integer :: k,k1
      integer :: l,ix
      integer :: i_ot,j_ot

      real :: x
      real :: dlon,dlat
      real :: dlon_ot,dlat_ot
      real :: alon_ot,alat_ot
      real :: wtlon,wtlat,wtlev

      real :: am,ap,akm,akp

	! Basic sanity checks

  if(present(stat)) stat=0

  if(nlon_ot<=0 .or. nlat_ot<=0) then
    if(nlon_ot<=0) call perr(myname_,'invalid nlon_ot',nlon_ot)
    if(nlat_ot<=1) call perr(myname_,'invalid nlat_ot',nlat_ot)
    if(.not.present(stat)) call die(myname_)
    stat=1
    return
  endif

  if(nlon < 2 .or. nlat < 2 .or. nlev < 1) then
    if(nlon < 2) call perr(myname_,'invalid nlon',nlon)
    if(nlat < 2) call perr(myname_,'invalid nlat',nlat)
    if(nlev < 1) call perr(myname_,'invalid nlev',nlev)
    if(.not.present(stat)) call die(myname_)
    stat=1
    return
  endif

	! For a GEOS-DAS A grid for the input 3-d field

		! Again, the input values for alat(:) are expected to
		! be within [-90.,+90.]; for alon(:) are expected to be
		! within [-180.,360.)

      dlon = (e_lon-w_lon)/nlon	! long. is a period grid
      dlat = (n_lat-s_lat)/(nlat-1)

!----------------------------------------
      dlon_ot = (e_lon_ot-w_lon_ot)/nlon_ot
      dlat_ot = (n_lat_ot-s_lat_ot)/(nlat_ot-1)

		! Indexing to the vertical level grid values

      call slogtab(.not.NEAREST,nlev,pres,	&
        1,plev_ot,k,wtlev			)
      k1=k+1

!----------------------------------------
      do j_ot = 1, nlat_ot
	alat_ot=(j_ot-1)*dlat_ot + s_lat_ot

		! Latitudinal direction indexing

		! Compute the deg. offset in [0,180.], then
		! mearue it by unit dlat

        x = (alat_ot - s_lat)/dlat

		! Index it as ix=floor(x), which is in [0,nlat]

	ix = int(x)
	if(ix > x) ix=ix-1	! a faster ix = floor(x)

		! Compute the linear weight

        wtlat = x - ix

	j =ix+1
	if(j >= nlat) then
	  wtlat=wtlat+j-nlat+1.
	  j=nlat-1
	endif
	j1=j +1

	do i_ot = 1, nlon_ot
	  alon_ot=(i_ot-1)*dlon_ot + w_lon_ot
     
		! Longitudinal direction indexing, assuming a periodic
		! grid.

		! Compute the deg. offset in [0,360), then
		! measure it by unit dlon.

	  x = modulo(alon_ot - w_lon, 360.)/dlon

		! Index it as ix=floor(x), which is in [0,nlon-1]

	  ix = int(x)
	  if(ix > x) ix=ix-1	! a faster ix = floor(x)

		! Compute the linear weight

          wtlon =  x - ix

		! the indices at j, i and i1

	  i =ix+1
	  i1=i +1
	  if(i1 > nlon) i1=1	! Circular indexing for periodic long.
             
		! Bilinear interpolation of a scalar field at the
		! lower level.

	am =grd(i ,j ,k )+wtlon*(grd(i1,j ,k )-grd(i ,j ,k ))
	ap =grd(i ,j1,k )+wtlon*(grd(i1,j1,k )-grd(i ,j1,k ))
	glev(i_ot,j_ot)=am+wtlat*(ap-am)

	if(k1 <= nlev) then

		! Bilinear interpolation of a scalar field at the
		! upper level.

	  am =grd(i ,j ,k1)+wtlon*(grd(i1,j ,k1)-grd(i ,j ,k1))
	  ap =grd(i ,j1,k1)+wtlon*(grd(i1,j1,k1)-grd(i ,j1,k1))
	  akp=am+wtlat*(ap-am)

	  glev(i_ot,j_ot)=glev(i_ot,j_ot)+wtlev*(akp-glev(i_ot,j_ot))
	endif

      enddo	! i_ot
    end do	! j_ot

      if(present(stat)) stat=0

end subroutine gintpscl_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: dintpvec_ - interpolate a A-grid vector field
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine dintpvec_(ndat,alat,alon,alev,udat,vdat,		&
	nlon,nlat,nlev,pres,ugrd,vgrd, stat	)

      use m_die,      only : die,perr
      use m_mall,     only : mall_ison,mall_mci,mall_mco
      implicit none

		! About the insitu vector data for output

      integer,intent(in) :: ndat	! number of insitu data
      real,dimension(:), intent(in)  :: alat	! latitudes
      real,dimension(:), intent(in)  :: alon	! longitudes
      real,dimension(:), intent(in)  :: alev	! pres. levels
      real,dimension(:), intent(out) :: udat	! insitu u-vect
      real,dimension(:), intent(out) :: vdat	! insitu v-vect

		! About the grid point vector data as the input

      integer,intent(in) :: nlat	! number of A-grid lat.
      integer,intent(in) :: nlon	! number of longitudes
      integer,intent(in) :: nlev	! number of levels
      real,dimension(:),    intent(in) :: pres	! pres. levels.
      real,dimension(:,:,:),intent(in) :: ugrd	! a grided u-vect
      real,dimension(:,:,:),intent(in) :: vgrd	! a grided v-vect

		! Optional arguments

      integer, optional, intent(out):: stat	! Status code

! !REVISION HISTORY:
!	24Apr00	- Jing Guo
!		. Created based on dintp_() of this module.
! 	19Mar97 - Jing Guo <guo@dao> - revised from intp_ap_vect_()
!EOP ___________________________________________________________________

	! Parameters

  character(len=*),parameter :: myname_=myname//'::dintpvec_'
  logical,parameter	     :: nearest=.true.

  real, parameter :: s_lat = -90., n_lat= +90.
  real, parameter :: w_lon =-180., e_lon=+180.
  real :: deg

	! Local variables

      integer :: ierr
      integer :: i,i1
      integer :: j,j1
      integer :: k,k1
      integer :: l,ix

      real :: x
      real :: dlon,dlat
      real :: wtlon,wtlat,wtlev
      real :: rlon,rlon1
      real :: rlat,rlat1
      real :: wxm,wym,wzm
      real :: wxp,wyp,wzp

      real :: ukm,ukp
      real :: vkm,vkp

	! Temporory working space

      real,              allocatable :: wght(:)
      integer,           allocatable :: glev(:)

	! Basic sanity checks

  if(present(stat)) stat=0

  if(nlon < 2 .or. nlat < 2 .or. nlev < 1) then
    if(nlon < 2) call perr(myname_,'invalid nlon',nlon)
    if(nlat < 2) call perr(myname_,'invalid nlat',nlat)
    if(nlev < 1) call perr(myname_,'invalid nlev',nlev)
    if(.not.present(stat)) call die(myname_)
    stat=1
    return
  endif

      deg=4.*atan(1.)/180.

	! For a GEOS-DAS A grid

		! The input values for alat(:) are expected to
		! be within [-90.,+90.]; for alon(:) are expected to be
		! within [-180.,360.)

      dlon = (e_lon-w_lon)/nlon
      dlat = (n_lat-s_lat)/(nlat-1)

      allocate(glev(ndat),wght(ndat),stat=ierr)
      if(ierr /= 0) then
	  call perr(myname_,'allocate()',ierr)
	  if(.not.present(stat)) call die(myname_)
	  stat=ierr
	  return
      endif
      if(mall_ison()) then
	call mall_mci(glev,myname)
	call mall_mci(wght,myname)
      endif

		! Indexing to the vertical level grid values

      call slogtab(.not.nearest,nlev,pres,	&
              ndat,alev,glev,wght		)
       
      do l = 1, ndat
     
		! Longitudinal direction indexing

        x = (alon(l) - w_lon)/dlon
	ix = int(x)
	if(ix > x) ix=ix-1	! a faster ix = floor(x)
        wtlon =  x - ix

	 	! ix = modulo(ix,nlon), works only upto 360.+180.
		! ix is in [0,nlon-1].

        if(ix >= nlon) ix=ix-nlon

		! the indices at j, i and i1

	i =ix+1
	i1=i +1
	if(i1 > nlon) i1=1	! Circular indexing
             
		! Latitudinal direction indexing

        x = (alat(l) - s_lat)/dlat
	ix = int(x)
	if(ix > x) ix=ix-1	! a faster ix = floor(x)
        wtlat = x - ix

	j =ix+1
	if(j >= nlat) then
	  j=nlat-1
	  wtlat=1.
	endif
	j1=j +1

	rlon =(i -1)*dlon+w_lon
	rlon1=(i1-1)*dlon+w_lon

	rlat =(j -1)*dlat+s_lat
	rlat1=(j1-1)*dlat+s_lat

        k = glev(l)
	wtlev=wght(l)
	k1=k+1

		! Bilinear interpolation of a vector field at the
		! lower level.

        call Linear2Ang_(wxp,wyp,wzp, wtlon  ,rlat1*deg,	&
		ugrd(i ,j1,k ),vgrd(i ,j1,k ),rlon *deg,	&
		ugrd(i1,j1,k ),vgrd(i1,j1,k ),rlon1*deg		)

        call Linear2Ang_(wxm,wym,wzm, wtlon  ,rlat *deg,	&
		ugrd(i ,j ,k ),vgrd(i ,j ,k ),rlon *deg,	&
		ugrd(i1,j ,k ),vgrd(i1,j ,k ),rlon1*deg		)

        call Ang2Linear_(udat(l),vdat(l),	&
		alon(l)*deg,alat(l)*deg,	&
		wtlat,wxm,wym,wzm,wxp,wyp,wzp)

	if(k1 <= nlev) then

		! Bilinear interpolation of a vector field at the
		! upper level.

          call Linear2Ang_(wxp,wyp,wzp, wtlon,rlat1*deg,	&
		ugrd(i ,j1,k1),vgrd(i ,j1,k1),rlon *deg,	&
		ugrd(i1,j1,k1),vgrd(i1,j1,k1),rlon1*deg		)

          call Linear2Ang_(wxm,wym,wzm, wtlon,rlat *deg,	&
		ugrd(i ,j ,k1),vgrd(i ,j ,k1),rlon *deg,	&
		ugrd(i1,j ,k1),vgrd(i1,j ,k1),rlon1*deg		)

          call Ang2Linear_(ukp,vkp,alon(l)*deg,alat(l)*deg,	&
		wtlat,wxm,wym,wzm,wxp,wyp,wzp)

	  udat(l)=udat(l)+wtlev*(ukp-udat(l))
	  vdat(l)=vdat(l)+wtlev*(vkp-vdat(l))

	endif

      enddo   ! ndat

	if(mall_ison()) then
	  call mall_mco(glev,myname)
	  call mall_mco(wght,myname)
	endif
      deallocate(glev,wght,stat=ierr)
      if(ierr /= 0) then
	call perr(myname_,'deallocate()',ierr)
	if(.not.present(stat)) call die(myname_)
	stat=ierr
	return
      endif

      if(present(stat)) stat=0

contains
	!--------------------------------------------------------
  subroutine Ang2Linear_(u,v,rlon,rlat, wt,wxm,wym,wzm,wxp,wyp,wzp)
    implicit none
    real, intent(out) :: u,v	! the resultant vector
    real, intent(in)  :: rlon,rlat	! Coordinates of (u,v)
    real, intent(in)  :: wt	! weight of the increment
    real, intent(in)  :: wxm,wym,wzm	! ang-vect. at the lower corner
    real, intent(in)  :: wxp,wyp,wzp	! ang-vect. at the upper corner

    real :: wx,wy,wz
    real :: coslon,sinlon
    real :: coslat,sinlat
    real :: emx,emy,emz
    real :: elx,ely

    wx=wxm+wt*(wxp-wxm)
    wy=wym+wt*(wyp-wym)
    wz=wzm+wt*(wzp-wzm)

    coslon=cos(rlon)
    sinlon=sin(rlon)
    coslat=cos(rlat)
    sinlat=sin(rlat)

    emx=-coslon*sinlat
    emy=-sinlon*sinlat
    emz=        coslat
    elx=-sinlon
    ely= coslon

    u=  wx*emx+wy*emy+wz*emz
    v=-(wx*elx+wy*ely)

  end subroutine Ang2Linear_
	!--------------------------------------------------------
  subroutine Linear2Ang_(wx,wy,wz, wt,rlat, um,vm,rlonm,up,vp,rlonp)
    implicit none
    real,intent(out) :: wx,wy,wz	! interpolated ang-vect.
    real,intent(in)  :: wt	! latitudinal weight
    real,intent(in)  :: rlat	! interpolation along the same lat.
    real,intent(in)  :: um,vm,rlonm	! vector at the lower corner
    real,intent(in)  :: up,vp,rlonp	! vector at the upper corner

    real :: wxp,wyp,wzp
    real :: coslon,sinlon
    real :: coslat,sinlat
    real :: emx,emy,emz
    real :: elx,ely

    coslon=cos(rlonm)
    sinlon=sin(rlonm)
    sinlat=sin(rlat)
    coslat=cos(rlat)

    emx=-coslon*sinlat
    emy=-sinlon*sinlat
    emz=        coslat
    elx=-sinlon
    ely= coslon

    wx=um*emx-vm*elx
    wy=um*emy-vm*ely
    wz=um*emz

    coslon=cos(rlonp)
    sinlon=sin(rlonp)
    sinlat=sin(rlat)
    coslat=cos(rlat)

    emx=-coslon*sinlat
    emy=-sinlon*sinlat
    emz=        coslat
    elx=-sinlon
    ely= coslon

    wxp=up*emx-vp*elx
    wyp=up*emy-vp*ely
    wzp=up*emz

    wx=wx+wt*(wxp-wx)
    wy=wy+wt*(wyp-wy)
    wz=wz+wt*(wzp-wz)
  end subroutine Linear2Ang_

end subroutine dintpvec_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: gintpvec_ - 2-d grid interpolate a 3-d A-grid vector field
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gintpvec_(nlon_ot,nlat_ot,plev_ot,ulev,vlev,	&
	nlon,nlat,nlev,pres,ugrd,vgrd, stat		)

      use m_die,      only : die,perr
      implicit none

		! About the insitu vector data for output

      integer,intent(in) :: nlon_ot
      integer,intent(in) :: nlat_ot
      real,   intent(in) :: plev_ot
      real,dimension(:,:),intent(out) :: ulev
      real,dimension(:,:),intent(out) :: vlev

		! About the grid point vector data as the input

      integer,intent(in) :: nlon	! number of longitudes
      integer,intent(in) :: nlat	! number of A-grid lat.
      integer,intent(in) :: nlev	! number of levels
      real,dimension(:),    intent(in) :: pres	! pres. levels.
      real,dimension(:,:,:),intent(in) :: ugrd	! a grided vector field
      real,dimension(:,:,:),intent(in) :: vgrd

		! Optional arguments

      integer, optional, intent(out):: stat	! Status code

! !REVISION HISTORY:
! 	14Jan99 - Jing Guo <guo@dao> -
!		. revised from intp_ap_.F90
!		. wrote a new prolog
!EOP ___________________________________________________________________

!#define _DEBUG	!

	! Parameters

  character(len=*),parameter :: myname_=myname//'::gintpvec_'
  logical,parameter	     :: NEAREST=.true.

  real, parameter :: s_lat_ot = -90., n_lat_ot= +90.
  real, parameter :: w_lon_ot =-180., e_lon_ot=+180.

  real, parameter :: s_lat = -90., n_lat= +90.
  real, parameter :: w_lon =-180., e_lon=+180.

	! Local variables

      integer :: ierr
      integer :: i,i1
      integer :: j,j1
      integer :: k,k1
      integer :: l,ix
      integer :: i_ot,j_ot

      real :: x
      real :: dlon,dlat
      real :: dlon_ot,dlat_ot
      real :: alon_ot,alat_ot
      real :: wtlon,wtlat,wtlev

      real :: wxm,wym,wzm
      real :: wxp,wyp,wzp
      real :: ukp,vkp
      real :: rlon,rlon1
      real :: rlat,rlat1
      real :: deg

	! Basic sanity checks

  if(present(stat)) stat=0

  if(nlon_ot<=0 .or. nlat_ot<=0) then
    if(nlon_ot<=0) call perr(myname_,'invalid nlon_ot',nlon_ot)
    if(nlat_ot<=1) call perr(myname_,'invalid nlat_ot',nlat_ot)
    if(.not.present(stat)) call die(myname_)
    stat=1
    return
  endif

  if(nlon < 2 .or. nlat < 2 .or. nlev < 1) then
    if(nlon < 2) call perr(myname_,'invalid nlon',nlon)
    if(nlat < 2) call perr(myname_,'invalid nlat',nlat)
    if(nlev < 1) call perr(myname_,'invalid nlev',nlev)
    if(.not.present(stat)) call die(myname_)
    stat=1
    return
  endif

	! For a GEOS-DAS A grid for the input 3-d field

		! Again, the input values for alat(:) are expected to
		! be within [-90.,+90.]; for alon(:) are expected to be
		! within [-180.,360.)

      dlon = (e_lon-w_lon)/nlon	! long. is a period grid
      dlat = (n_lat-s_lat)/(nlat-1)

      deg = 4.*atan(1.)/180.

!----------------------------------------
      dlon_ot = (e_lon_ot-w_lon_ot)/nlon_ot
      dlat_ot = (n_lat_ot-s_lat_ot)/(nlat_ot-1)

		! Indexing to the vertical level grid values

      call slogtab(.not.NEAREST,nlev,pres,	&
        1,plev_ot,k,wtlev			)
      k1=k+1

!----------------------------------------
      do j_ot = 1, nlat_ot
	alat_ot=(j_ot-1)*dlat_ot + s_lat_ot

		! Latitudinal direction indexing

		! Compute the deg. offset in [0,180.], then
		! mearue it by unit dlat

        x = (alat_ot - s_lat)/dlat

		! Index it as ix=floor(x), which is in [0,nlat]

	ix = int(x)
	if(ix > x) ix=ix-1	! a faster ix = floor(x)

		! Compute the linear weight

        wtlat = x - ix

	j =ix+1
	if(j >= nlat) then
	  wtlat=wtlat+j-nlat+1.
	  j=nlat-1
	endif
	j1=j +1

	rlat =(j -1)*dlat+s_lat
	rlat1=(j1-1)*dlat+s_lat

	do i_ot = 1, nlon_ot
	  alon_ot=(i_ot-1)*dlon_ot + w_lon_ot
     
		! Longitudinal direction indexing, assuming a periodic
		! grid.

		! Compute the deg. offset in [0,360), then
		! measure it by unit dlon.

	  x = modulo(alon_ot - w_lon, 360.)/dlon

		! Index it as ix=floor(x), which is in [0,nlon-1]

	  ix = int(x)
	  if(ix > x) ix=ix-1	! a faster ix = floor(x)

		! Compute the linear weight

          wtlon =  x - ix

		! the indices at j, i and i1

	  i =ix+1
	  i1=i +1
	  if(i1 > nlon) i1=1	! Circular indexing for periodic long.

	  rlon =(i -1)*dlon+w_lon
	  rlon1=(i1-1)*dlon+w_lon
             
		! Bilinear interpolation of a vector field at the
		! lower level.

        call Linear2Ang_(wxp,wyp,wzp, wtlon  ,rlat1*deg,	&
		ugrd(i ,j1,k ),vgrd(i ,j1,k ),rlon *deg,	&
		ugrd(i1,j1,k ),vgrd(i1,j1,k ),rlon1*deg		)

        call Linear2Ang_(wxm,wym,wzm, wtlon  ,rlat *deg,	&
		ugrd(i ,j ,k ),vgrd(i ,j ,k ),rlon *deg,	&
		ugrd(i1,j ,k ),vgrd(i1,j ,k ),rlon1*deg		)

        call Ang2Linear_(ulev(i_ot,j_ot),vlev(i_ot,j_ot),	&
		alon_ot*deg,alat_ot*deg,			&
		wtlat,wxm,wym,wzm,wxp,wyp,wzp)

	if(k1 <= nlev) then

		! Bilinear interpolation of a vector field at the
		! lower level.

          call Linear2Ang_(wxp,wyp,wzp, wtlon  ,rlat1*deg,	&
		ugrd(i ,j1,k ),vgrd(i ,j1,k ),rlon *deg,	&
		ugrd(i1,j1,k ),vgrd(i1,j1,k ),rlon1*deg		)

          call Linear2Ang_(wxm,wym,wzm, wtlon  ,rlat *deg,	&
		ugrd(i ,j ,k ),vgrd(i ,j ,k ),rlon *deg,	&
		ugrd(i1,j ,k ),vgrd(i1,j ,k ),rlon1*deg		)

          call Ang2Linear_(ukp,vkp,		&
		alon_ot*deg,alat_ot*deg,	&
		wtlat,wxm,wym,wzm,wxp,wyp,wzp)

	  ulev(i_ot,j_ot)=ulev(i_ot,j_ot)+wtlev*(ukp-ulev(i_ot,j_ot))
	  vlev(i_ot,j_ot)=vlev(i_ot,j_ot)+wtlev*(vkp-vlev(i_ot,j_ot))
	endif

      enddo	! i_ot
    end do	! j_ot

      if(present(stat)) stat=0

contains
	!--------------------------------------------------------
  subroutine Ang2Linear_(u,v,rlon,rlat, wt,wxm,wym,wzm,wxp,wyp,wzp)
    implicit none
    real, intent(out) :: u,v	! the resultant vector
    real, intent(in)  :: rlon,rlat	! Coordinates of (u,v)
    real, intent(in)  :: wt	! weight of the increment
    real, intent(in)  :: wxm,wym,wzm	! ang-vect. at the lower corner
    real, intent(in)  :: wxp,wyp,wzp	! ang-vect. at the upper corner

    real :: wx,wy,wz
    real :: coslon,sinlon
    real :: coslat,sinlat
    real :: emx,emy,emz
    real :: elx,ely

    wx=wxm+wt*(wxp-wxm)
    wy=wym+wt*(wyp-wym)
    wz=wzm+wt*(wzp-wzm)

    coslon=cos(rlon)
    sinlon=sin(rlon)
    coslat=cos(rlat)
    sinlat=sin(rlat)

    emx=-coslon*sinlat
    emy=-sinlon*sinlat
    emz=        coslat
    elx=-sinlon
    ely= coslon

    u=  wx*emx+wy*emy+wz*emz
    v=-(wx*elx+wy*ely)

  end subroutine Ang2Linear_
	!--------------------------------------------------------
  subroutine Linear2Ang_(wx,wy,wz, wt,rlat, um,vm,rlonm,up,vp,rlonp)
    implicit none
    real,intent(out) :: wx,wy,wz	! interpolated ang-vect.
    real,intent(in)  :: wt	! latitudinal weight
    real,intent(in)  :: rlat	! interpolation along the same lat.
    real,intent(in)  :: um,vm,rlonm	! vector at the lower corner
    real,intent(in)  :: up,vp,rlonp	! vector at the upper corner

    real :: wxp,wyp,wzp
    real :: coslon,sinlon
    real :: coslat,sinlat
    real :: emx,emy,emz
    real :: elx,ely

    coslon=cos(rlonm)
    sinlon=sin(rlonm)
    sinlat=sin(rlat)
    coslat=cos(rlat)

    emx=-coslon*sinlat
    emy=-sinlon*sinlat
    emz=        coslat
    elx=-sinlon
    ely= coslon

    wx=um*emx-vm*elx
    wy=um*emy-vm*ely
    wz=um*emz

    coslon=cos(rlonp)
    sinlon=sin(rlonp)
    sinlat=sin(rlat)
    coslat=cos(rlat)

    emx=-coslon*sinlat
    emy=-sinlon*sinlat
    emz=        coslat
    elx=-sinlon
    ely= coslon

    wxp=up*emx-vp*elx
    wyp=up*emy-vp*ely
    wzp=up*emz

    wx=wx+wt*(wxp-wx)
    wy=wy+wt*(wyp-wy)
    wz=wz+wt*(wzp-wz)
  end subroutine Linear2Ang_

end subroutine gintpvec_
end module m_intpAPelem
