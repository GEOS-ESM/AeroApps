!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_psasrc - PSAS default resource handle
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_psasrc
      implicit none
      private	! except

      public :: PSAS_DEFRC	! the default resource filename
      public :: psasrc_name	! the current resource filename
      public :: psasrc_open	! open the current resource
      public :: psasrc_close	! close the current resource

      interface psasrc_open; module procedure	&
	open_,	&
	allopen_
      end interface

      interface psasrc_close; module procedure	&
	close_
      end interface

! !REVISION HISTORY:
! 	20Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- Split the data in psasrc.h/psasrcbd() to m_psasrc
!		  and m_psasverssion.
!		- Removed psasrc.h and psasrcbd.f, an effort to clean
!		  BLOCKDATA structures.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_psasrc'

!	(psasrc.h)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE:	psasrc.h - data module defines standard PSAS resources
!
! !SYNOPSIS:	#include "psasrc.h"
!
! !DESCRIPTION:
!	"psasrc.h" defines following standard PSAS resources,
!
!		character*4 psasname	(='PSAS')
!		character*6 version	(='x.x.x[abd]')
!		character*7 def_psasrc	(='psas.rc')
!
!	as well as a global variable "psasrc" for a shared resource
!	filename.  The variable value may be defined by a user.
!
!	The identities therefore are reserved where the module is
!	included.
!
!	"psasname" and "version" should be used in output to mark the
!	version of PSAS library in used.  "def_psasrc" is used to
!	define a default filename for PSAS resource input.  An
!	alternative may be defined through an environment variable or
!	possiblely by user specification.  See !EXAMPLS.
!
! !EXAMPLES:
!			external psasrcbd	! a blockdata unit
!		#include "psasrc.h"
!		...
!			if(psasrc.eq.' ') call getenv('PSASRC',psasrc)
!			if(psasrc.eq.' ') psasrc=def_psasrc
! !BUGS:
!	When there is a change in the values of the identities, remember
!	to update the length of the string.
!
! !SEE ALSO: psasrcbd.f
!
! !REVISION HISTORY:
!	17Oct96 - J. Guo	- 1.1.4d with inpak90 etc.
!	21Sep95 - Jing G.	- included common/rsrcfile/
! 	14Sep95 - J. Guo	- Added the prolog and modified comments
!	01Feb95 - Jing G.	- Created to store standard resources
!_______________________________________________________________________

!	..Shared label for information output
!	=====================================

	! The default resource filename

  character(len=*),parameter :: PSAS_DEFRC='psas.rc'

	! On-line resource filename

  character(len=128),save    :: psasrc_name=PSAS_DEFRC

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: open_ - determine the current resource filename
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine open_(rsrc,stat)
      use m_die,    only : perr,die
      use m_inpak90,only : i90_loadf
      use m_mpout,  only : mpout_log
      implicit none
      character(len=*),optional,intent(in) :: rsrc
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
!	12Jan01	- Jing Guo
!		. Coded specially for single PE situation.  The earlier
!		  implementation is not correct.
! 	03Jan00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::open_'
  integer :: ier

		! Inquering the resource handle

  if(present(stat)) stat=0

  if(present(rsrc)) then

		! Loading the resource handle

    call I90_LoadF(rsrc,ier)
	if(ier/=0) then
	  call perr(myname_,'i90_Loadf("'//trim(rsrc)//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	call mpout_log(myname_, 'using "'//		&
		trim(rsrc)//'" for the runtime resource input')

  else
    psasrc_name=""
    call getenv('PSASRC',psasrc_name)		! Unix binding
    if(psasrc_name.eq.' ') psasrc_name=PSAS_DEFRC	! default name

		! Loading the resource handle

    call I90_LoadF(psasrc_name,ier)
	if(ier/=0) then
	  call perr(myname_,	&
		'i90_LoadF("'//trim(psasrc_name)//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	call mpout_log(myname_, 'using "'//		&
		trim(psasrc_name)//'" for the runtime resource input')

  endif

end subroutine open_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allopen_ - determine the current resource filename
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allopen_(root,comm,rsrc,stat)
      use m_die,    only : perr,die
      use m_inpak90,only : i90_allLoadf
      use m_mpout,  only : mpout_log
      implicit none
      integer,intent(in) :: root
      integer,intent(in) :: comm
      character(len=*),optional,intent(in) :: rsrc
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	03Jan00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allopen_'
  integer :: ier

		! Inquering the resource handle

  if(present(stat)) stat=0

  if(present(rsrc)) then

		! Loading the resource handle

    call I90_allLoadF(rsrc,root,comm,ier)
	if(ier/=0) then
	  call perr(myname_,'i90_allLoadf("'//trim(rsrc)//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	call mpout_log(myname_, 'using "'//		&
		trim(rsrc)//'" for the runtime resource input')

  else
    psasrc_name=""
    call getenv('PSASRC',psasrc_name)		! Unix binding
    if(psasrc_name.eq.' ') psasrc_name=PSAS_DEFRC	! default name

		! Loading the resource handle

    call I90_allLoadF(psasrc_name,root,comm,ier)
	if(ier/=0) then
	  call perr(myname_,	&
		'i90_allLoadF("'//trim(psasrc_name)//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	call mpout_log(myname_, 'using "'//		&
		trim(psasrc_name)//'" for the runtime resource input')

  endif

end subroutine allopen_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: close_ - close the resource file
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine close_(stat)
      use m_die,    only : perr,die
      use m_inpak90,only : i90_release
      implicit none
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	03Jan00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::close_'
  integer :: ier

  if(present(stat)) stat=0

  call i90_release(ier)
	if(ier/=0) then
	  call perr(myname_,'i90_release()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine close_
end module m_psasrc
!.
