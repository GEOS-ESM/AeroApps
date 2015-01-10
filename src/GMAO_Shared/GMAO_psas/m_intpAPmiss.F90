!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_intpAPmiss - Interpolations from a GEOS A-pressure grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_intpAPmiss
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

  character(len=*),parameter :: myname='m_intpAPmiss'

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
	nlon,nlat,nlev,pres,grd, missing, stat		)

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

      real,   intent(in) :: missing	! the value of missing-scalar

		! Optional arguments

      integer, optional, intent(out):: stat	! Status code

! !REVISION HISTORY:
! 	19Mar97 - Jing Guo <guo@dao> - revised from intp_ap_vect_()
!EOP ___________________________________________________________________

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

      real :: am,ap,akm,akp

	! Temperory working space

      real,              allocatable :: wght(:)
      integer,           allocatable :: glev(:)

	! Function

      real    :: threshold	! The threshold value of missing

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

	! Set the threadhold, it is the missing value less a possible
	! truncation error, to prevent the situation where the system
	! created data has a different machine real representation from
	! the system compiled the program.

      threshold=abs(missing)*(1.-epsilon(1.))

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

	  i =ix+1	! i is in [1,nlon]
	  i1=i +1	! i1 is in {[2,nlon],1}
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

		! The basic approach is that for any 1-d interpolation
		! (3d=1d+1d+1d), if one of grid point value is missing,
		! the other will be taken; if both are missing, a 
		! missing value will be taken; if none is missing, a
		! regular linearly interpolated value is taken.
		!
		! If the all of 8 points in the 3d interpolation are
		! missing, then the final value is missing.

	akp=missing	! let it be missing, in case of k1 > nlev.
	if(k1 <= nlev) then

		! Bilinear interpolation of a scalar field with
		! possiblly missing values, at the upper level.

	  call BiLinear_(akp, wtlon,wtlat,	&
	    grd(i ,j ,k1), grd(i1,j ,k1),	&
	    grd(i ,j1,k1), grd(i1,j1,k1)	)

	endif

		! Bilinear interpolation of a scalar field with
		! possiblly missing values, at the lower level.

	call BiLinear_(akm, wtlon,wtlat,	&
	    grd(i ,j ,k ), grd(i1,j ,k ),	&
	    grd(i ,j1,k ), grd(i1,j1,k )	)

		! Final interpolation between two levels.

	call Linear_(adat(l), wtlev, akm,akp)

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
  function notMissing_(a)
    implicit none
    real, intent(in) :: a
    logical :: notMissing_	! the result

    notMissing_ = abs(a)< threshold
  end function notMissing_
	!--------------------------------------------------------
  subroutine Linear_(a, wt,am,ap)
    implicit none
    real, intent(out) :: a	! the resultant scalar
    real, intent(in)  :: wt	! weight of the increment
    real, intent(in)  :: am	! the scalar at the lower corner
    real, intent(in)  :: ap	! the scalar at the upper corner

    a=ap
    if(   notMissing_(am) ) then
      a=am
      if( notMissing_(ap) ) a=a+wt*(ap-a)
    endif
  end subroutine Linear_
	!--------------------------------------------------------
  subroutine BiLinear_(a, wx_,w_x, amm,apm, amp,app)
    implicit none
    real, intent(out) :: a	! the resultant scalar
    real, intent(in) :: wx_,w_x	! weights of the increments
    real, intent(in) :: amm	! the scalar at the  left-lower corner
    real, intent(in) :: apm	! the scalar at the right-lower corner
    real, intent(in) :: amp	! the scalar at the  left-upper corner
    real, intent(in) :: app	! the scalar at the right-upper corner

    real :: am,ap

    call Linear_(ap, wx_,amp,app)
    call Linear_(am, wx_,amm,apm)
    call Linear_(a,  w_x,am, ap )

  end subroutine BiLinear_
	!--------------------------------------------------------
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
	nlon,nlat,nlev,pres,grd, missing, stat		)

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

      real,   intent(in) :: missing	! the value of missing-scalar

		! Optional arguments

      integer, optional, intent(out):: stat	! Status code

! !REVISION HISTORY:
! 	14Jan99 - Jing Guo <guo@dao> -
!		. revised from intp_ap_.F90
!		. wrote a new prolog
!EOP ___________________________________________________________________

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

	! Function

      real    :: threshold	! The threshold value of missing

	! Basic sanity checks

  if(present(stat)) stat=0

  if(nlon_ot<=0 .or. nlat_ot<=0) then
    if(nlon_ot<=0) call perr(myname_,'invalid nlon_ot',nlon_ot)
    if(nlat_ot<=0) call perr(myname_,'invalid nlat_ot',nlat_ot)
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

	! Set the threadhold, it is the missing value less a possible
	! truncation error, to prevent the situation where the system
	! created data has a different machine real representation from
	! the system compiled the program.

      threshold=abs(missing)*(1.-epsilon(1.))

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
             
		! The basic approach is that for any 1-d interpolation
		! (3d=1d+1d+1d), if one of grid point value is missing,
		! the other will be taken; if both are missing, a 
		! missing value will be taken; if none is missing, a
		! regular linearly interpolated value is taken.
		!
		! If the all of 8 points in the 3d interpolation are
		! missing, then the final value is missing.

	akp=missing	! let it be missing, in case of k1 > nlev.
	if(k1 <= nlev) then

		! Bilinear interpolation of a scalar field with
		! possiblly missing values, at the upper level.

	  call BiLinear_(akp, wtlon,wtlat,	&
	    grd(i ,j ,k1), grd(i1,j ,k1),	&
	    grd(i ,j1,k1), grd(i1,j1,k1)	)

	endif

		! Bilinear interpolation of a scalar field with
		! possiblly missing values, at the lower level.

	call BiLinear_(akm, wtlon,wtlat,	&
	    grd(i ,j ,k ), grd(i1,j ,k ),	&
	    grd(i ,j1,k ), grd(i1,j1,k )	)

		! Final interpolation between two levels.

	call Linear_(glev(i_ot,j_ot), wtlev, akm,akp)

      enddo	! i_ot
    end do	! j_ot

      if(present(stat)) stat=0

contains
	!--------------------------------------------------------
  function notMissing_(a)
    implicit none
    real, intent(in) :: a
    logical :: notMissing_	! the result

    notMissing_ = abs(a)< threshold
  end function notMissing_
	!--------------------------------------------------------
  subroutine Linear_(a, wt,am,ap)
    implicit none
    real, intent(out) :: a	! the resultant scalar
    real, intent(in)  :: wt	! weight of the increment
    real, intent(in)  :: am	! the scalar at the lower corner
    real, intent(in)  :: ap	! the scalar at the upper corner

    a=ap
    if(   notMissing_(am) ) then
      a=am
      if( notMissing_(ap) ) a=a+wt*(ap-a)
    endif
  end subroutine Linear_
	!--------------------------------------------------------
  subroutine BiLinear_(a, wx_,w_x, amm,apm, amp,app)
    implicit none
    real, intent(out) :: a	! the resultant scalar
    real, intent(in) :: wx_,w_x	! weights of the increments
    real, intent(in) :: amm	! the scalar at the  left-lower corner
    real, intent(in) :: apm	! the scalar at the right-lower corner
    real, intent(in) :: amp	! the scalar at the  left-upper corner
    real, intent(in) :: app	! the scalar at the right-upper corner

    real :: am,ap

    call Linear_(ap, wx_,amp,app)
    call Linear_(am, wx_,amm,apm)
    call Linear_(a,  w_x,am, ap )

  end subroutine BiLinear_
	!--------------------------------------------------------
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
	nlon,nlat,nlev,pres,ugrd,vgrd, umiss,vmiss, stat	)

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

      real,   intent(in) :: umiss	! the value of missing-u
      real,   intent(in) :: vmiss	! the value of missing-v

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

	! Temperory working space

      real,              allocatable :: wght(:)
      integer,           allocatable :: glev(:)

	! Threshold values of either u or v being missing

      real    :: uepsilon,vepsilon

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

	! Set the threadhold values of being missing.

      uepsilon=abs(umiss)*epsilon(1.)
      vepsilon=abs(vmiss)*epsilon(1.)

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
	  call mall_mco(glev,myname)
	  call mall_mco(wght,myname)
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

		! The basic approach is that for any 1-d interpolation
		! (3d=1d+1d+1d), if one of grid point value is missing,
		! the other will be taken; if both are missing, a 
		! missing value will be taken; if none is missing, a
		! regular linearly interpolated value is taken.
		!
		! If the all of 8 points in the 3d interpolation are
		! missing, then the final value is missing.

	ukp=umiss	! let it be missing, in case of k1 > nlev.
	vkp=vmiss
	if(k1 <= nlev) then

		! Bilinear interpolation of a vector field with
		! possiblly missing values, at the upper level.

          call Linear2Ang_(wxp,wyp,wzp, wtlon,rlat1*deg,	&
		ugrd(i ,j1,k1),vgrd(i ,j1,k1),rlon *deg,	&
		ugrd(i1,j1,k1),vgrd(i1,j1,k1),rlon1*deg		)

          call Linear2Ang_(wxm,wym,wzm, wtlon,rlat *deg,	&
		ugrd(i ,j ,k1),vgrd(i ,j ,k1),rlon *deg,	&
		ugrd(i1,j ,k1),vgrd(i1,j ,k1),rlon1*deg		)

          call Ang2Linear_(ukp,vkp,alon(l)*deg,alat(l)*deg,	&
		wtlat,wxm,wym,wzm,wxp,wyp,wzp)

	endif

		! Bilinear interpolation of a vector field with
		! possiblly missing values, at the lower level.

        call Linear2Ang_(wxp,wyp,wzp, wtlon  ,rlat1*DEG,	&
		ugrd(i ,j1,k ),vgrd(i ,j1,k ),rlon *DEG,	&
		ugrd(i1,j1,k ),vgrd(i1,j1,k ),rlon1*DEG		)

        call Linear2Ang_(wxm,wym,wzm, wtlon  ,rlat *DEG,	&
		ugrd(i ,j ,k ),vgrd(i ,j ,k ),rlon *DEG,	&
		ugrd(i1,j ,k ),vgrd(i1,j ,k ),rlon1*DEG		)

        call Ang2Linear_(ukm,vkm,alon(l)*DEG,alat(l)*DEG,	&
		wtlat,wxm,wym,wzm,wxp,wyp,wzp)

		! Final interpolation between two levels.

	call Linear_(udat(l),vdat(l), wtlev, ukm,ukp,vkm,vkp)

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
  function notMissing_(u,v)
    implicit none
    real, intent(in) :: u,v
    logical :: notMissing_	! the result

		! == epsilon is not used in case that epsilon
		! values are zeroes (either umiss==0 or vmiss==0).

    notMissing_ = abs(u-umiss)>uepsilon .and.	&
		  abs(v-vmiss)>vepsilon

		! Note that both components must be not-missing for a
		! vector to be not-missing.

  end function notMissing_
	!--------------------------------------------------------
  subroutine Linear_(u,v, wt,um,up,vm,vp)
    implicit none
    real, intent(out) :: u,v	! the resultant vector
    real, intent(in)  :: wt	! weight of the increment
    real, intent(in)  :: um,vm	! the vector at the lower corner
    real, intent(in)  :: up,vp	! the vector at the upper corner

    u=up
    v=vp
    if(   notMissing_(um,vm) ) then
      u=um
      v=vm
      if( notMissing_(up,vp) ) then
	u=u+wt*(up-u)
	v=v+wt*(vp-v)
      endif
    endif
  end subroutine Linear_
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

    wx=wxp
    wy=wyp
    wz=wzp
    if(notMissing_(wxm,wym)) then
      wx=wxm
      wy=wym
      wz=wzm
      if(notMissing_(wxp,wyp)) then
	wx=wx+wt*(wxp-wx)
	wy=wy+wt*(wyp-wy)
	wz=wz+wt*(wzp-wz)
      endif
    endif

    u=umiss
    v=vmiss
    if(notMissing_(wx,wy)) then

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

    endif
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

    wxp=umiss
    wyp=vmiss
    wzp=0.
    if(notMissing_(up,vp)) then

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
    endif

    wx=wxp
    wy=wyp
    wz=wzp
    if(notMissing_(um,vm)) then

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

      if(notMissing_(up,vp)) then
	wx=wx+wt*(wxp-wx)
	wy=wy+wt*(wyp-wy)
	wz=wz+wt*(wzp-wz)
      endif
    endif
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
	nlon,nlat,nlev,pres,ugrd,vgrd, umiss,vmiss, stat	)

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

      real,   intent(in) :: umiss	! the value of missing-vector
      real,   intent(in) :: vmiss

		! Optional arguments

      integer, optional, intent(out):: stat	! Status code

! !REVISION HISTORY:
! 	14Jan99 - Jing Guo <guo@dao> -
!		. revised from intp_ap_.F90
!		. wrote a new prolog
!EOP ___________________________________________________________________

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
      real :: ukm,vkm
      real :: ukp,vkp
      real :: rlon,rlon1
      real :: rlat,rlat1
      real :: deg

	! The threshold value of missing

      real :: uepsilon,vepsilon

	! Basic sanity checks

  if(present(stat)) stat=0

  if(nlon_ot<=0 .or. nlat_ot<=0) then
    if(nlon_ot<=0) call perr(myname_,'invalid nlon_ot',nlon_ot)
    if(nlat_ot<=0) call perr(myname_,'invalid nlat_ot',nlat_ot)
    if(.not.present(stat)) call die(myname_)
    stat=1
    return
  endif

  if(nlon < 2 .or. nlat < 2 .or. nlev < 1) then
    if(nlon < 2) call perr(myname_,'invalid nlon =',nlon)
    if(nlat < 2) call perr(myname_,'invalid nlat =',nlat)
    if(nlev < 1) call perr(myname_,'invalid nlev =',nlev)
    if(.not.present(stat)) call die(myname_)
    stat=1
    return
  endif

      deg=4.*atan(1.)/180.

	! Set the threadhold values of being missing.

      uepsilon=abs(umiss)*epsilon(1.)
      vepsilon=abs(vmiss)*epsilon(1.)

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

	rlat =(j -1)*dlat + s_lat
	rlat1=(j1-1)*dlat + s_lat

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

	  rlon =(i -1)*dlon + w_lon
	  rlon1=(i1-1)*dlon + w_lon
             
		! The basic approach is that for any 1-d interpolation
		! (3d=1d+1d+1d), if one of grid point value is missing,
		! the other will be taken; if both are missing, a 
		! missing value will be taken; if none is missing, a
		! regular linearly interpolated value is taken.
		!
		! If the all of 8 points in the 3d interpolation are
		! missing, then the final value is missing.

	ukp=umiss	! let it be missing, in case of k1 > nlev.
	vkp=vmiss
	if(k1 <= nlev) then

		! Bilinear interpolation of a vector field with
		! possiblly missing values, at the upper level.

          call Linear2Ang_(wxp,wyp,wzp, wtlon,rlat1*deg,	&
		ugrd(i ,j1,k1),vgrd(i ,j1,k1),rlon *deg,	&
		ugrd(i1,j1,k1),vgrd(i1,j1,k1),rlon1*deg		)

          call Linear2Ang_(wxm,wym,wzm, wtlon,rlat *deg,	&
		ugrd(i ,j ,k1),vgrd(i ,j ,k1),rlon *deg,	&
		ugrd(i1,j ,k1),vgrd(i1,j ,k1),rlon1*deg		)

          call Ang2Linear_(ukp,vkp,alon_ot*deg,alat_ot*deg,	&
		wtlat,wxm,wym,wzm,wxp,wyp,wzp)

	endif

		! Bilinear interpolation of a vector field with
		! possiblly missing values, at the lower level.

        call Linear2Ang_(wxp,wyp,wzp, wtlon  ,rlat1*DEG,	&
		ugrd(i ,j1,k ),vgrd(i ,j1,k ),rlon *DEG,	&
		ugrd(i1,j1,k ),vgrd(i1,j1,k ),rlon1*DEG		)

        call Linear2Ang_(wxm,wym,wzm, wtlon  ,rlat *DEG,	&
		ugrd(i ,j ,k ),vgrd(i ,j ,k ),rlon *DEG,	&
		ugrd(i1,j ,k ),vgrd(i1,j ,k ),rlon1*DEG		)

        call Ang2Linear_(ukm,vkm,alon_ot*DEG,alat_ot*DEG,	&
		wtlat,wxm,wym,wzm,wxp,wyp,wzp)

		! Final interpolation between two levels.

	call Linear_(ulev(i_ot,j_ot),vlev(i_ot,j_ot),	&
		wtlev, ukm,ukp,vkm,vkp)

      enddo	! i_ot
    end do	! j_ot

      if(present(stat)) stat=0

contains
	!--------------------------------------------------------
  function notMissing_(u,v)
    implicit none
    real, intent(in) :: u,v
    logical :: notMissing_	! the result

		! == epsilon is not used in case that epsilon
		! values are zeroes (either umiss==0 or vmiss==0).

    notMissing_ = abs(u-umiss)>uepsilon .and.	&
		  abs(v-vmiss)>vepsilon

		! Note that both components must be not-missing for a
		! vector to be not-missing.

  end function notMissing_
	!--------------------------------------------------------
  subroutine Linear_(u,v, wt,um,up,vm,vp)
    implicit none
    real, intent(out) :: u,v	! the resultant vector
    real, intent(in)  :: wt	! weight of the increment
    real, intent(in)  :: um,vm	! the vector at the lower corner
    real, intent(in)  :: up,vp	! the vector at the upper corner

    u=up
    v=vp
    if(   notMissing_(um,vm) ) then
      u=um
      v=vm
      if( notMissing_(up,vp) ) then
	u=u+wt*(up-u)
	v=v+wt*(vp-v)
      endif
    endif
  end subroutine Linear_
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

    wx=wxp
    wy=wyp
    wz=wzp
    if(notMissing_(wxm,wym)) then
      wx=wxm
      wy=wym
      wz=wzm
      if(notMissing_(wxp,wyp)) then
	wx=wx+wt*(wxp-wx)
	wy=wy+wt*(wyp-wy)
	wz=wz+wt*(wzp-wz)
      endif
    endif

    u=umiss
    v=vmiss
    if(notMissing_(wx,wy)) then

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

    endif
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

    wxp=umiss
    wyp=vmiss
    wzp=0.
    if(notMissing_(up,vp)) then

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
    endif

    wx=wxp
    wy=wyp
    wz=wzp
    if(notMissing_(um,vm)) then

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

      if(notMissing_(up,vp)) then
	wx=wx+wt*(wxp-wx)
	wy=wy+wt*(wyp-wy)
	wz=wz+wt*(wzp-wz)
      endif
    endif
  end subroutine Linear2Ang_
	!--------------------------------------------------------
end subroutine gintpvec_
end module m_intpAPmiss
