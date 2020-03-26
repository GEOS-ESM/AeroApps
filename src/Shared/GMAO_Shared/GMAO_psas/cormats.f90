module cormats
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: cormats - interface to correlation matrix subroutines
!
! !INTERFACE:
!	module cormats
!	  public fDDcor1,fDDcorx,fHDcorx,fHDcorx,fQQcor1,fQQcorx,
!	  public oHHcor1,oHHcorx,uHHcor1
!	end module cormats
!
! !DESCRIPTION:
!	Module `cormats' defines the interface to correlation matrix
!	subroutines.  Currently, these subroutines are all written in
!	Fortram 77 to ensure their performance on some system.
!
!	Subroutines defined in this interface module include:
!
!	  fDDcor1()	for a diagonal matrix off either of two height
!			gradiant components.  Forecast errors.
!	  fDDcorx()	for a off-diagonal matrix of two height gradiant
!			components, either the same or different
!			components.  Forecast errors.
!	  fHDcorx()	for a off-diagonal matrix of height and either
!			of two its gradiant components.  Forecast
!			errors.
!	  fHHcor1()	for a diagonal matrix of height.  Forecast
!			error.
!	  fHHcorx()	for a off-diagonal matrix of height.  Forecast
!			error.
!	  fQQcor1()	for a diagonal matrix of Q, mixing ratio.
!			Forecast error.
!	  fQQcorx()	for a off-diagonal matrix of Q, mixing ratio.
!			Forecast error.
!	  oHHcor1()	for a diagonal matrix of height.  Correlated
!			(the same kx, data source, only) observation
!			error.
!	  oHHcorx()	for a off-diagonal matrix of height.  Correlated
!			(the same kx data source only) observation
!			error.
!	  uHHcor1()	for a diagonal matrix of height.  Uncorrelated
!			(vertical or the same profile only) observation
!			error.
!
! !REVISION HISTORY:
! 	27Feb96 - J. Guo	- subroutines have been changed with
!		different arguments.  Changes include,
!
!		  1) removed linear interpolation weights (wtab, vtab)
!		  2) removed Hint array for the original mass-balanced
!		     wind correlation functions
!		  3) removed the sign from some of the subrotines.
!
! 	25Jan96 - J. Guo	- initial version for PSAS version 1
!		correlation models.  With explicitly specified Hint[]
!		array for mass-balanced mass-wind correlations, and a
!		sign (+-) for on some of the subroutines.
!_______________________________________________________________________
implicit none
private

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		! Only callable names
public :: fDDcor1, fDDcorx, fHDcorx, fHHcor1, fHHcorx
public :: oHHcor1, oHHcorx, uHHcor1
!_______________________________________________________________________

interface
!=======================================================================
  subroutine fDDcor1(kind_cov, len, irstr, iistr, coords, iattr, idx, corrF	)
    implicit none

  ! Input arguments:

    integer, intent(in)	:: kind_cov	! type of correlation matrix
    integer, intent(in) :: len
    integer, intent(in) :: irstr, iistr
    real,    intent(in) :: coords(irstr,len)
    integer, intent(in) :: iattr(iistr,len)
    integer, intent(in) :: idx ! control - should be set to either iqm_x or iql_x

  ! Output argument:
    real,   intent(out)	:: corrF(len*(len+1)/2)  ! <D,D> forecast corr.
!!$    real corrF(size(coords,2)*(size(coords,2)+1)/2)  ! <D,D> forecast correlation

  end subroutine fDDcor1
!=======================================================================
  subroutine fDDcorx(kind_cov,					&
		leni, irstr_i, iistr_i, coords_i, iattr_i, idx_i, &
		lenj, irstr_j, iistr_j, coords_j, iattr_j, idx_j, &
		corrF						)
    implicit none

  ! Input arguments:

    integer,intent(in)	:: kind_cov	! type of correlation matrix

    integer,intent(in)	:: leni		! row dimension of the matrix
    integer, intent(in) :: irstr_i, iistr_i
    real, intent(in) :: coords_i(irstr_i,leni)   ! attribute data
    integer, intent(in) :: iattr_i(iistr_i,leni)
    integer, intent(in) :: idx_i

    integer,intent(in)	:: lenj		! column dimension of the matrix
    integer, intent(in) :: irstr_j, iistr_j
    real, intent(in) :: coords_j(irstr_j,lenj)   ! attribute data
    integer, intent(in) :: iattr_j(iistr_j,lenj)
    integer, intent(in) :: idx_j

  ! Output argument:
    real,   intent(out)	:: corrF(leni,lenj)  ! <L,M> forecast corr.

  end subroutine fDDcorx
!=======================================================================
  subroutine fHDcorx(kind_cov, &
       leni, irstr_i, iistr_i, coords_i, iattr_i,  &
       lenj, irstr_j, iistr_j, coords_j, iattr_j, idx_j, &
       corrF					    )
    implicit none

  ! Input arguments:

    integer	kind_cov	! type of correlation matrix

    integer	leni		! row dimension of the matrix
    integer	iistr_i, irstr_i
    integer	iistr_j, irstr_j
    real  coords_i(irstr_i,leni)
    integer iattr_i(iistr_i,leni)

    integer	lenj		! column dimension of the matrix
    real coords_j(irstr_j,lenj)
    integer iattr_j(iistr_j,lenj)
    integer idx_j

  ! Output argument:
    real	corrF(leni,lenj)  ! <H,D> forecast correlation

  end subroutine fHDcorx
!=======================================================================
  subroutine fHHcor1(kind_cov, coords, iattr, corrF)
    implicit none

!  Input arguments:
! =================
    integer	kind_cov	! type of correlation matrix
    real,    intent(in) :: coords(:,:)
    integer, intent(in) :: iattr(:,:)

!  Output arguments:
! ==================
!!$    real		corrF(len*(len+1)/2)	! upper-triangle matrix
    real corrF(size(coords,2)*(size(coords,2)+1)/2)  ! <D,D> forecast correlation

  end subroutine fHHcor1
!=======================================================================
!!$  subroutine fHHcorx(kind_cov,leni,stride_i,coords_i, ktabi, &
!!$		     lenj,stride_j,coords_j, ktabj, corrF)
  subroutine fHHcorx(kind_cov, &
       leni, irstr_i, iistr_i, coords_i, iattr_i, &
       lenj, irstr_j, iistr_j, coords_j, iattr_j, corrF)
    implicit none

! Input arguments:
!=================

    integer	kind_cov	! type of correlation matrix
    
    integer	leni		! row dimension of the matrix
    integer, intent(in) :: irstr_i, iistr_i
    integer, intent(in) :: irstr_j, iistr_j
    real    :: coords_i(irstr_i,leni)
    integer :: iattr_i(iistr_i,leni)
    integer	lenj		! column dimension of the matrix
    real    :: coords_j(irstr_i,lenj)
    integer :: iattr_j(irstr_j,lenj)

! Output argument:
!=================
    real	corrF(leni,lenj)  ! <H,H> forecast correlation
!!$	real    :: corrF(size(coords_i,2),size(coords_j,2))

  end subroutine fHHcorx
!=======================================================================
!!$  subroutine fQQcor1(len,qr_x,qr_y,qr_z,ktab, corrF)
  subroutine fQQcor1(coords, iattr, corrF)
    implicit none

!  Input arguments
! =================
    real,    intent(in) :: coords(:,:)
    integer, intent(in) :: iattr(:,:)

!  Output argument
! =================
!!$    real		corrF(len*(len+1)/2)	! packed matrix
    real corrF(size(coords,2)*(size(coords,2)+1)/2)  ! <D,D> forecast correlation
  end subroutine fQQcor1
!=======================================================================
  subroutine fQQcorx(leni, irstr_i, iistr_i, coords_i, iattr_i, &
		     lenj, irstr_j, iistr_j, coords_j, iattr_j, corrF)
    implicit none

! Input arguments:
!=================
    integer, intent(in) :: leni
    integer, intent(in) :: irstr_i, iistr_i
    integer, intent(in) :: irstr_j, iistr_j
    real    :: coords_i(irstr_i,leni)
    integer :: iattr_i(iistr_i,leni)
    integer	lenj		! column dimension of the matrix
    real    :: coords_j(irstr_i,lenj)
    integer :: iattr_j(irstr_j,lenj)

! Output argument:
!=================
    real	corrF(leni,lenj)  ! <Q,Q> forecast correlation
  end subroutine fQQcorx
!=======================================================================
!!$  subroutine oHHcor1(len,kx,qr_x,qr_y,qr_z,ktab, Mtyp,corrO)
  subroutine oHHcor1(coords, iattr, Mtyp,corrO)
    implicit none

! Inputs:
!========
!!$    integer		len		! dimension of the matrix
!!$    integer		kx(len)		! data source id
!!$    real		qr_x(len)	! clat*clon
!!$    real		qr_y(len)	! clat*slon
!!$    real		qr_z(len)	! slat
    real :: coords(:,:)
!!$    integer		ktab(len)	! reference to voecH table
    integer :: iattr(:,:)

! Outputs:
!==========
    character*1	Mtyp		! "I" or "U", or ..
!!$    real		corrO(len*(len+1)/2)	! packed matrix
    real corrO(size(coords,2)*(size(coords,2)+1)/2)  ! <D,D> forecast correlation

  end subroutine oHHcor1
!=======================================================================
!!$  subroutine oHHcorx(leni,kxi,qri_x,qri_y,qri_z,ktabi,	&
!!$		     lenj,kxj,qrj_x,qrj_y,qrj_z,ktabj,	&
!!$		     Mtyp,corrO				)
  subroutine oHHcorx(coords_i, iattr_i,  &
		     coords_j, iattr_j,	&
		     Mtyp,corrO				)
    implicit none

! Input arguments:
!=================

!!$    integer	leni		! row dimension of the matrix
!!$    integer kxi(leni)	! data source IDs
!!$    real	qri_x(leni)	! x of qri (location vector of rows)
!!$    real	qri_y(leni)	! y of qri
!!$    real	qri_z(leni)	! z of qri
    real, intent(in) :: coords_i(:,:)
!!$    integer	ktabi(leni)	! indices of rows to vfecH table
    integer, intent(in) :: iattr_i(:,:)

!!$    integer	lenj		! column dimension of the matrix
!!$    integer kxj(lenj)	! data source IDs
!!$    real	qrj_x(lenj)	! x of qri (location vector of columns)
!!$    real	qrj_y(lenj)	! y of qri
!!$    real	qrj_z(lenj)	! z of qri
    real, intent(in) :: coords_j(:,:)
!!$    integer	ktabj(lenj)	! indices of columns to vfecH table
    integer, intent(in) :: iattr_j(:,:)

! Output argument:
!=================
    character*1 Mtyp	! "N"ormal or "Z"ero
!!$    real	corrO(leni,lenj)  ! <H,H> forecast correlation
	real    :: corrO(size(coords_i,2),size(coords_j,2))

  end subroutine oHHcorx

!!$  subroutine uHHcor1(len,kx,ks,ktab, Mtyp,corrU)
  subroutine uHHcor1(iattr, Mtyp,corrU)
    implicit none

!  Inputs:
! ========
!!$    integer		len		! dimension of the matrix
!!$    integer		kx(len)		! data source IDs
!!$    integer		ks(len)		! data sounding indices
!!$    integer		ktab(len)	! reference to voecH table
    integer, intent(in) :: iattr(:,:)
!!$    real		wtab(len)	! weight to voecH table

!  Outputs:
! =========
    character*1	Mtyp		! "I" or "U"
!!$    real		corrU(len*(len+1)/2)	! packed matrix
    real corrU(size(iattr,2)*(size(iattr,2)+1)/2)  ! <D,D> forecast correlation

  end subroutine uHHcor1
end interface
end module cormats
!.
