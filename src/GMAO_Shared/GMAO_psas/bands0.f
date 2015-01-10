C***********************************************************************
C                          subroutine bands0
C***********************************************************************
C
C WRITTEN BY:              Jim Pfaendtner
C
C!DESCRIPTION:             Set parameters for conjugate gradients 
C                          descent
C
C CALLED FROM:             startup
C
C SYSTEM ROUTINES USED:    none
C
C SUBROUTINES CALLED:      none
C
C!INPUT PARAMETERS:        none
C
C!OUTPUT PARAMETERS:       none
C
C!INPUT/OUTPUT PARAMETERS: none
C
C FILES USED:              none
C
C!REVISION HISTORY:
C  10oct94 - A. da S. - Introduced freq to output Anal. incs.
C  04oct94 - A. da S. - back from blockdata to subroutine
C  21sep94 - Jing G.  - retyped as a blockdata module
C  23jun94 - Jim Pf.  - Removed parameters from .h file
C  09apr94 - Jim Pf.  - Added prologue, put parameters in .h file
C  02oct93 - Jim Pf.  - Original code
!  25jul96 - Jing Guo - changed minmax to minpass and criter(:) to
!			criter(:,2)
!  15Apr99 - Jing Guo - replaced stdio.h with m_stdio
!		      - added BLOCKDATA to avoid warning messages
C
C***********************************************************************

      blockdata bands_bd

      include     "bands.h"

      DATA seplim /      .00,      .00,      .00,    26.50,    58.25/
      DATA maxpass/       30,       16,        8,        8,        8/
      DATA minpass /       3,        3,        3,        3,        3/
      DATA criter /       .1,       .1,       .1,       .1,       .1,
     &                    .3,       .3,       .3,       .3,       .3/
      DATA cgname /'conjgr1','conjgr2','conjgr3','conjgr4','conjgr5'/
      DATA msmall /       16/

      end blockdata bands_bd

      subroutine bands0

	use m_inpak90, only : lablin,intget,fltget
	use m_stdio,   only : stderr,stdout
	use m_die,only : die
c.......................................................................

*     Standard I/O units
*     ------------------
      include     "bands.h"
      include	  "kind_mats.h"

      character*6 myname
      parameter(myname='bands0')

      real :: deg

	! These are old variables which may be removed sometime.

      method = 1
      nbandcg=kind_5mat
      diagscal=1.

      call LABLIN ( 'conjgr_separation_limits:' )
      do 10 i = 1, NBANDS
         seplim(i) = FLTGET ( seplim(i) )
 10   continue

      call LABLIN ( 'conjgr_maximum_no_iterations:' )
      do 20 i = 1, NBANDS
         maxpass(i) = INTGET ( maxpass(i) )
 20   continue

      call LABLIN ( 'conjgr_minimum_no_iterations:' )
      do 40 i = 1, NBANDS
         minpass(i) = INTGET ( minpass(i) )
 40   continue

      call LABLIN ( 'conjgr_max_tolerances:' )
      do i = 1, NBANDS
         criter(i,1) = FLTGET ( criter(i,1) )
      end do

      call LABLIN ( 'conjgr_min_tolerances:' )
      do i = 1, NBANDS
         criter(i,2) = FLTGET ( criter(i,2) )
      end do

       call LABLIN ( 'conjgr_M_small:' )
       msmall = INTGET ( msmall )

c   ..Input verbose options for each level(1 for the top level, 0 o.w.)
      call lablin('conjgr_verbose:')
      do i=1,nbands
	if(i.eq.nbandcg) then
	  cgverb(i)=intget(1).gt.0
	else
	  cgverb(i)=intget(0).gt.0
	endif
      end do

	deg=4.*atan(1.)/180.
	do i=1,nbands
	  cosseplim(i)=cos(seplim(i)*deg)
	end do


      return
      end
