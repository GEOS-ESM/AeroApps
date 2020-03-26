C***********************************************************************
C              subroutine restrict
C***********************************************************************
C
C WRITTEN BY:  Jim Searle
C
C!DESCRIPTION: Restricts data to be used in analysis
C
C CALLED FROM: analdrv
C
C SYSTEM ROUTINES USED: none
C
C SUBROUTINES CALLED:
C  blklist  - Excludes data contained in a union of 
C             (latitude) by (logitude) by (pressure layer) boxes
C  initkl   - Initializes logical vector kl with toggled data choices
C  llboxes  - Restricts data to that contained in a union of
C             (latitude) by (logitude) by (pressure layer) boxes
C  prtdata  - Prints vectors associated with selected data
C  settog   - Toggles data on/off by kt=data type and kx=data source
C  tofront  - Compresses unselected data out of data vectors
C
C!INPUT PARAMETERS:
C  luverb   - integer unit number for verbose ascii output
C  nobs     - integer number of observations in data vectors on input
C  turnoff  - logical switch for print of final data selection
C  verbose  - logical switch for verbose ascii output
C
C!OUTPUT PARAMETERS:
C  nnobs    - integer length of data vectors after all selections
C             have been done 
C
C!INPUT/OUTPUT PARAMETERS: 
C             NOTE: all of these vectors are of length nobs on input
C  del      - real innovation vector (of length nnobs) with the 
C             (observation - first guess) differences, the units of 
C             the k'th component depends on the data type kt(k)
C               kt(k)    data_type                  units
C                 1        near-surface u             m/s
C                 2        near-surface v             m/s
C                 3        sea level pressure         mb
C                 4        upper-air u                m/s
C                 5        upper-air v                m/s
C                 6        geopotential height        m
C                 7        water vapor mixing ratio   g/kg
C  kt       - integer vector (of length nnobs) with the data type
C             indices for the observations
C  kx       - integer vector (of length nnobs) with the data source
C             indices for the observations
C  qcflag   - real vector (of length nnobs) with data quality control 
C             flag,  qcflag(k) has the following meanings:
C               qcflag = 1.0  ---> 
C  rlats    - real vector (of length nobs/nobs) with the latitudes of the 
C             observations in degrees (-90.0 <= rlats(k) <= 90.0)
C  rlevs    - real vector (of length nnobs) with the pressure level of
C             the observations in mb
C  rlons    - real vector (of length nnobs) with the longitudes of the 
C             observations in degrees (-180.0 <= rlats(k) <= 180.0)
C  sigF     - real vector of estimated background error standard 
C             deviations, units are the same as for the innovation 
C             vector (del)
C  sigO     - real vector of estimated observation error standard 
C             deviations, units are the same as for the innovation 
C             vector (del)
C  tstamp   - real vector (of length nnobs) with time stamps for the 
C             individual observations, the value given is the offset 
C             in minutes from the synoptic time 
C               -180.0 <= tstamp(k) <= 180.0
C
C FILES USED:
C  print *  - standard ascii output
C  luverb   - ascii output to unit luverb
C
C!REVISION HISTORY:
!
!	20Dec99	- Jing Guo
!		. revised to fit module m_boxes.
!  30Apr99 - Todling  - increase max number of boxes (nbmax) from 24 to 100 
!  08Oct98 - G. Brin  - changed the number of boxes (nbmax) from 9 to 24 
!                       to accomodate DAO TOVS boxes.
!  09May96 - J. Guo   - New setbox()/llboxes() with interface changes
!  26Apr96 - J. Guo   - Add mark_nsig() that marks off negative sigU or
!			sigO values.
C  04oct94 - A. da S. - Calls to blklist commented out.
C                       Changed call to SETBOX and dimensions of boxes.
C  26apr94 - Jim Pf   - Added prologue and call to blklist
C  02apr93 - James S. - Code modifications based on analy.f from IRIS
C
C***********************************************************************

      subroutine RESTRICT (verbose, luverb,   nobs,   turnoff,
     1			   nboxes,  boxes,
     1                     rlats,   rlons,   rlevs,   kx,      
     1                     kt,      del,     sigU,    sigO,    sigF,   
     1                     tstamp,  nnobs             )

      use m_die,  only : die
      use m_stdio,only : stdout

      implicit none

*     Control parameter for output.
*     -----------------------------
      integer      luverb
      logical      verbose
      logical      turnoff

	integer,intent(in) :: nboxes
	real,dimension(2,6,nboxes),intent(in) :: boxes

*     Size parameters for database.
*     -----------------------------
      include     "ktmax.h"
      include     "kxmax.h"
      integer      nobs
      integer      nnobs


*     Data items.
*     -----------
      real         rlats(nobs)
      real         rlons(nobs)
      real         rlevs(nobs)
      integer      kx(nobs)
      integer      kt(nobs)
      real         del(nobs)
      real         sigU(nobs)
      real         sigO(nobs)
      real         sigF(nobs)
      real         tstamp(nobs)


*     Local storage.
*     -------------

	logical, allocatable :: kl(:)
	integer istat

	character(len=*), parameter :: myname = 'restrict'

*     with subroutine llboxes which is used below.  For diagnostic
*     purposes, the formats in some subroutines have assumed the number
*     boxes (nboxes) has a single digit, or nboxes < 10.
*     ------------------------------------------------------------------

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

*....................................................................


*     At this point no data are turned on for use in the experiment.
*     Select the data wanted for the current experiment by kx and kt
*     --------------------------------------------------------------
	allocate(kl(nobs), stat=istat)
	if(istat.ne.0) call die(myname,'allocate()',istat)

*     Set definitions of geographic lat/lon/lev boxes for data selection

*     Initialize logical array kl with this choice of data.  It first
*     check kx-kt combination against ktmask(kx) to see if the data
*     type is valid.
*     ---------------------------------------------------------------
      ! call INITKL ( verbose, luverb, nobs, kx, kt, kl )

	! kx-kt mask checking removed.

      kl(1:nobs)=.true.


*     Restrict the choice to data within a union of lat/lon/lev boxes.
*     We want to set  kl(n) = kl(n).and.(box1.or.box2.or....or.boxk).
*     ---------------------------------------------------------------
      call LLBOXES ( verbose, luverb, nboxes, boxes,
     1               nobs,    rlons,  rlats,  rlevs, kx, kt, tstamp,
     &		     kl )

!-----------------------------------------------------------------------
!	! Mark off observations with either of negative sigU/sigO/sigF
!
!	call mark_nsig(nobs,sigU,sigO,kl)
!
!-----------------------------------------------------------------------
	! Makr off regional

!	call mark_reg(nobs,rlats,rlons, kl)

!     Slide the selected data to the front of the arrays.
!     --------------------------------------------------
      call TOFRONT ( nobs,  kx,    kt,    kl,   
     1               rlats, rlons, rlevs,
     1               del,   sigU,   sigO,  sigF,  tstamp, nnobs )

!     --------------------------------------------------
	write(stdout,'(a,2(a,i6))') myname,': from ',nobs,' to ',nnobs


!--- begin(tbrm)
*     Print the selected data.
*     -----------------------
      call PRTDATA ( turnoff, luverb, nnobs,   kx,      kt,     kl,
     1               rlats,   rlons,  rlevs,   del,     sigO,   sigF,
     1               tstamp                                          )
!--- end(tbrm)


	deallocate(kl)
      return
      end
