module cordriv
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: cordriv - a F90 module of PSAS' correlation function drivers
!
! !INTERFACE:
!	use cordriv
!
! !DESCRIPTION:
!	cordriv() is a F90 module of PSAS' correlation function driver
!	subroutines.  It defines the interfaces of the following
!	correlation function driver subroutines,
!
!	  diagcorF: forming a diagonal block of the forecast error
!		    correlation matrix;
!	  diagcorO: forming a diagonal block of the observation error
!		    correlation matrix with correlated horizontal
!		    correlation functions;
!	  diagcorU: forming a diagonal block of the observation error
!		    correlation matrix with uncorrelated horizontal
!		    correlation functions;
!	  offdcorF: forming a off-diagonal block of the forecast error
!		    correlation matrix;
!	  offdcorO: forming a off-diagonal block of the observation
!		    error correlation matrix with correlation horizontal
!		    correlation functions.
!
! !BUGS:
!	Several changes are already planned for the interface of the
!	four existent subroutines, as while as a group of new entries
!	to be included in to this group.
!	
!	Although it is convenient to use the module to insure the
!	correctness of the calling interface, it does not help much
!	if the actual subroutines/functions are coded different from
!	this interface.  So it is important to make sure the actual
!	interfaces are accurately recorded in the this module.  An
!	automatic process is desireble.
!
! !SEE ALSO:
!	diagcorF(), diagcorO(), diagcorU(), offdcorF(), and offdcorO().
!
! !REVISION HISTORY:
!	25Jan95 - J. Guo	- added diagcorU() interface.  Some
!		corrections also made according to the changes in other
!		include subroutines.
!
! 	15Dec95 - J. Guo	- programed and added the prolog
!_______________________________________________________________________
interface
!=======================================================================
!!$  subroutine diagcorF(	kind_cov,				&
!!$			kt,len,	qr_x,qr_y,qr_z,			&
!!$				qm_x,qm_y,qm_z,ql_x,ql_y,	&
!!$				ktab,				&
!!$			Mtyp,corF,istat				)
  subroutine diagcorF(	kind_cov,				&
			kt, iattr, rattr, &
			Mtyp,corF,istat				)
!-----------------------------------------------------------------------
    integer, intent(in)	:: kind_cov
    integer, intent(in)	:: kt
    real,    intent(in) :: rattr(:,:)
    integer, intent(in) :: iattr(:,:)
    character*1, intent(out)	:: Mtyp
!!$    real,    intent(out)	:: corF(len*(len+1)/2)
    real,    intent(out)	:: corF(size(iattr,2)*(size(iattr,2)+1)/2)
    integer, intent(out)	:: istat
  end subroutine diagcorF
!=======================================================================
  subroutine diagcorO(	kt, iattr, rattr,  &
			Mtyp,corO, istat	)
!-----------------------------------------------------------------------
    integer, intent(in)	:: kt
    real, intent(in)    :: rattr(:,:)
    integer, intent(in) :: iattr(:,:)
    character*1, intent(out)	:: Mtyp
!!$    real,    intent(out)	:: corO(len*(len+1)/2)
    real, intent(out) :: corO(size(iattr,2)*(size(iattr,2)+1)/2)
    integer, intent(out)	:: istat
  end subroutine diagcorO
!=======================================================================
!!$  subroutine diagcorU(	kt,len,	kx,ks,ktab,Mtyp,corU,istat)
  subroutine diagcorU(	kt, iattr, Mtyp, corU, istat)
!-----------------------------------------------------------------------
    integer, intent(in)	:: kt		! data type
    integer, intent(in) :: iattr(:,:)
    character*1, intent(out)	:: Mtyp
!!$    real,    intent(out)	:: corU(len*(len+1)/2)
    real, intent(out) :: corU(size(iattr,2)*(size(iattr,2)+1)/2)
    integer, intent(out)	:: istat
  end subroutine diagcorU
!=======================================================================
!!$  subroutine offdcorF(kind_cov,					&
!!$		      kti,leni,	ktabi, stride_i, coords_i,                &
!!$		      ktj,lenj,	ktabj, stride_j, coords_j,                &
!!$		      Mtyp,corF,istat				)
  subroutine offdcorF(kind_cov,					&
		      kti,rattr_i, iattr_i,                &
		      ktj,rattr_j, iattr_j,                &
		      Mtyp,corF,istat				)
!-----------------------------------------------------------------------
    integer, intent(in)	:: kind_cov

    integer, intent(in)	:: kti
    integer, intent(in) :: iattr_i(:,:)
    real,    intent(in) :: rattr_i(:,:)

    integer, intent(in)	:: ktj
    integer, intent(in) :: iattr_j(:,:)
    real,    intent(in) :: rattr_j(:,:)

    character*1, intent(out)	:: Mtyp
!!$    real,    intent(out)	:: corF(leni*lenj)
    real, intent(out) :: corf(:)
    integer, intent(out)	:: istat
  end subroutine offdcorF
!=======================================================================
!!$  subroutine offdcorO(kti,leni,kxi,qri_x,qri_y,qri_z,	&
!!$				ktabi,				&
!!$		      ktj,lenj,kxj,qrj_x,qrj_y,qrj_z,	&
!!$				ktabj,				&
!!$		      Mtyp,corO,istat				)
  subroutine offdcorO(kti,rattr_i, iattr_i, &
                      ktj,rattr_j, iattr_j, &
		      Mtyp,corO,istat				)
!-----------------------------------------------------------------------
    integer, intent(in)	:: kti
!!$    integer, intent(in)	:: leni
!!$    integer, intent(in)	::  kxi(leni)
!!$    real,    intent(in)	:: qri_x(leni)
!!$    real,    intent(in)	:: qri_y(leni)
!!$    real,    intent(in)	:: qri_z(leni)
!!$    integer, intent(in)	:: ktabi(leni)
    integer, intent(in) :: iattr_i(:,:)
    real, intent(in) :: rattr_i(:,:)

    integer, intent(in)	:: ktj
!!$    integer, intent(in)	:: lenj
!!$    integer, intent(in)	::  kxj(lenj)
!!$    real,    intent(in)	:: qrj_x(lenj)
!!$    real,    intent(in)	:: qrj_y(lenj)
!!$    real,    intent(in)	:: qrj_z(lenj)
!!$    integer, intent(in)	:: ktabj(lenj)
    integer, intent(in) :: iattr_j(:,:)
    real, intent(in) :: rattr_j(:,:)

    character*1, intent(out)	:: Mtyp
!!$    real,    intent(out)	:: corO(leni*lenj)
    real, intent(out) :: corO(size(iattr_i,2)*size(iattr_j,2))
    integer, intent(out)	:: istat
  end subroutine offdcorO
!=======================================================================
end interface
end module cordriv
