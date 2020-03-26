**************************************************************************
*                                                                        *
*                  NASA/Goddard Space Flight Center                      *
*                     Laboratory for Atmospheres                         *
*                      Data Assimilation Office                          *
*                            Code 910.3                                  *
*                                                                        *
*            PHYSICAL-SPACE STATISTICAL ANALYSIS SYSTEM (PSAS)           *
*                                                                        *
**************************************************************************


*........................... ROUTINE PROLOGUE .............................
* !BOP
*
* !FILE: kxname0.f
*
* !ROUTINE: kxname0
* 
* !DESCRIPTION: 
*
*        Initialize array 'kxname' with the known data source names,
*  and array 'kxrank' for superobing in 'proxelim' package.
*
* !CALLING SEQUENCE:
*
*        call KXNAME0
*
* !INPUT PARAMETERS: none
*
* !OUTPUT PARAMETERS: 
*
*       Modifies array 'kxname' in common block.
*
* !SEE ALSO: kxname.h, kxmax.h
*
* !SYSTEM ROUTINES: none.
*
* !FILES USED: none.
*
* !WRITTEN BY: Jim Pf, 02/05/93 
* 
* !REVISION HISTORY:
*
*    05oct94  Data sources > 87 disabled for now (A. da S.)
*    27sep94  New KX list with longer names and more sources provided
*              by D. Lamich. Added initialization of kx ranking array
*              (A. da Silva) 
*
* !EOP
*.........................................................................

      subroutine kxname0
      use m_die,  only : die
      use m_stdio,only : stdout

c.......................................................................
c.... Table kxname to hold data source names.

	use m_output, only : output_ObsErrCov
	use m_stdio
	implicit none
      include     "ktmax.h"
      include     "kxmax.h"
      include     "kxtabl.h"

c  Locals
c ========
	character*7 myname
	parameter  (myname='kxname0')

	integer ikx,istat

c	..Read the table
c	================
	call rdkxtbl(RC_kx,kxmax,kxclas,kxrank,kxdesc,istat)
	if(istat.ne.0) call die(myname,'rdkxtbl()',istat)

c	..Echo input in its original form
c	=================================
	if(output_ObsErrCov) then

	  write(stdout,'(/a)') RC_kx
	  do ikx=1,kxmax
	    write(stdout,'(i3,2x,a,2x,i7,2(x,a))') ikx,
     &	      kxclas(ikx),kxrank(ikx),kxdesc(ikx),'$'
	  end do
	  write(stdout,'(a)') '::'

	endif

      return
      end
