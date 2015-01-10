!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_cache - cache and optimization related parameters
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_cache
      implicit none
      private	! except

      public :: JBANDS		! see m_Cxpy()
      public :: JMINOR		! see m_Cxpy()
      public :: BSIZE		! see m_Cxpy()
      public :: NBIG_BLOCK_COLS	! see m_Cxpy()

! !REVISION HISTORY:
! 	20Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_cache'

#ifdef sysIRIX64

#ifdef _JBANDS_
  integer,parameter :: JBANDS = _JBANDS_
#else
  integer,parameter :: JBANDS = 16
#endif

#ifdef _JMINOR_
  integer,parameter :: JMINOR = _JMINOR_
#else
  integer,parameter :: JMINOR = 3
#endif

#ifdef	_NBIGBC_
  integer,parameter :: NBIG_BLOCK_COLS	= _NBIGBC_
#else
  integer,parameter :: NBIG_BLOCK_COLS	= 1
#endif

#ifdef _BSIZE_
  integer,parameter :: BSIZE=_BSIZE_
#else
  integer,parameter :: BSIZE=-1
#endif

#else

  integer,parameter :: JBANDS = 1
  integer,parameter :: JMINOR = 1
  integer,parameter :: NBIG_BLOCK_COLS	= 1
  integer,parameter :: BSIZE=-1
#endif

end module m_cache
