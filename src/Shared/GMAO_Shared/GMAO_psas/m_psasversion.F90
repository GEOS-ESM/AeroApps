!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_psasversion - PSAS version module
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_psasversion
      implicit none
      private	! except

      public :: PSAS_name
      public :: PSAS_vers

! !REVISION HISTORY:
! 	20Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- Split the data in psasrc.h/psasrcbd() to m_psasrc
!		  and m_psasversion.
!		- Removed psasrc.h and psasrcbd.f, an effort to clean
!		  BLOCKDATA structures.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_psasversion'

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


  character(len=*),parameter :: PSAS_name='PSAS'
  character(len=*),parameter :: PSAS_vers ='1.3.04'

end module m_psasversion
!.
