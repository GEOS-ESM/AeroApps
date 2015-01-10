module FEalpha_imat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: FEalpha_imat - indexed matrix/transpose of ALPHA operator
!
! !INTERFACE:
!	use FEalpha_imat[, only : Aum_imat,Avm_imat,Aul_imat,Avl_imat]
!
! !DESCRIPTION:
!	Indexed matrix/transpose of ALPHA operator, for wind to height
!	gradiant fields transformation.  All imats are referenced with
!	respect to a lookup table of pveclev(nveclev) for the first
!	dimension and to ilat = 1 + (rlat-90.) * (MXlat_imat-1)/180. for
!	the second dimension.
!
! !EXAMPLES: (to do)
! !BUGS: (to do)
!
! !SEE ALSO: mvAH(), mvAX(), levtabl.h
!
! !REVISION HISTORY:
! 	09Feb96 - J. Guo	- (to do)
!	19Sep96 - J. Guo	- add `save' to all public variables
!_______________________________________________________________________

  use config, only : MXveclev, MXveclat
  implicit none
  private
	! except these

  public Aum_imat,Avm_imat,Aul_imat,Avl_imat,FEalpha_imat0

  real, save :: Aum_imat(MXveclev,MXveclat)
  real, save :: Avm_imat(MXveclev,MXveclat)
  real, save :: Aul_imat(MXveclev,MXveclat)
  real, save :: Avl_imat(MXveclev,MXveclat)

Contains

  subroutine FEalpha_imat0
  end subroutine FEalpha_imat0

end module FEalpha_imat
!.
