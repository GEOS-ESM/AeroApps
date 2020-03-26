
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.... Initialize logical array kl with a particular choice of data.

      subroutine initkl( verbose, luverb, nobs, kx, kt, kl )

c	20Jan95 - Jing G.	- Added ktneeded filter to remove any
c		data type not needed by a given analysis run.
c	06Jan95 - Jing G.	- Removed the use of kxkttog variable.
c				- Changed `kxkttog.h' to `kxktmask.h'.
 
c.......................................................................
c.... Argument declarations.

      logical      verbose
      integer      luverb
      integer      nobs
      integer      kx(nobs)
      integer      kt(nobs)
      logical      kl(nobs)

c.......................................................................
c.... Table of toggle switches for data source/type combinations.

      include     "kxmax.h"
      include     "ktmax.h"
      include     "kttabl.h"
      include	  "kxtabl.h"
      include	  "ktwanted.h"

c	..Local vars., flag of certain data type is needed by this
c	analysis, according the output list ktwanted[].
c	==========================================================
	logical ktneeded(ktmax)
	integer ikt,jkt
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c	..Elimiting irrelevent data types (not correlated with any
c	`wanted' analysis data types.
c	==========================================================
	do ikt=1,ktmax
	  ktneeded(ikt)=.false.
	  do jkt=1,ktmax
	    ktneeded(ikt)=ktneeded(ikt).or.
     &	      (ktwanted(jkt).and.ktmvar(ikt,jkt))
	  end do
	end do

c.......................................................................
c.... Set logical array kl to reflect the data choices in kxtmask

      kount = 0
      do 100 n = 1, nobs

         kl(n) = kxtmask(kx(n),kt(n)).ne.0 .and. ktneeded(kt(n))
         if( kl(n) ) kount = kount + 1

  100 continue

      if( verbose ) write(luverb,920) kount, nobs
  920 format('   initkl: kl is on in ',I7,' cases out of ',I7/
     $       '   initkl:')

      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end

