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
* !FILE: proxel.F
*
* !ROUTINE: proxel0
* 
* !DESCRIPTION: 
*
*     Initializes Proximity Elimination package.
*
* !CALLING SEQUENCE:
*
*     call PROXEL0
*
* !INPUT PARAMETERS: none
*
* !OUTPUT PARAMETERS: none
*
* !SEE ALSO: startup
*
* !SYSTEM ROUTINES: 
*
* !FILES USED: 
*
* !WRITTEN BY:  Arlindo da Silva, 5oct94
* 
* !REVISION HISTORY:
* ! a bug is fixed. use of temporary storage for ireglen added.
* ! got rid of cpp/gpp directives: f90 use of dynamic space.
*    G. Brin Mar 96.
*
* !EOP
*.........................................................................

      subroutine PROXEL0

      use m_inpak90, only : lablin, fltget,chrget
      use m_stdio,   only : stdout

      include     'kxmax.h'
      include     'proxel.h'

      character*1 answer

      write(stdout,*)
      write(stdout,*) '  +------------------------------------+'
      write(stdout,*) '  |  Proximity Elimination Parameters  |'
      write(stdout,*) '  +------------------------------------+'
      write(stdout,*)

!      write(stdout,*)
!      write(stdout,*) 'Enter radius in km for proximity check : '
      call LABLIN ( 'radius_in_km_for_proximity_check:' )
      Rkm = FLTGET ( 80. )
      Rtresh = Rkm / Rerth    ! Rerth is in km
      Rtresh2 = Rtresh**2     ! save in common for later use

!      write(stdout,*) 'Enter vertical distance (in Ln-p) for ',
!     .   'proximity check'
      call LABLIN ( 'del_lnp_for_vertical_proximity_check:' )
      dellnp = abs ( FLTGET ( 0.1 ) )
   
!      write(stdout,*)
!      write(stdout,*) 'Do you want to print a detailed superob list?'
      call LABLIN ( 'do_you_want_detailed_superob_list:' )
      answer = CHRGET ( 'n' )
      debug = .false.
      if ( answer .eq. 'y' .or. answer .eq. 'Y' ) then
           debug = .true.
      end if

      write(stdout,*)
      write(stdout,*) '               --- // ---'
      write(stdout,*)

      return
      end


*.........................................................................
* !BOP
*
* !FILE: proxel.f
*
* !ROUTINE: proxel
* 
* !DESCRIPTION: 
*
*    This routine scans the input innovation vectors elimitating those
*  data points that are "too close". The general algorithm follows:
*
*  \begin{enumerate}
*
*  \item The input data comes in sorted by region, type, source, 
*        lat, lon, level.   For each data point, its region 
*  and nearby regions are scanned for other data of same type (kt) 
*  that are horizontally close (the threshhold distance 
*  Rkm is set  in 'proxel.h') and those data points are saved on a 
*  black-list (local data structure defined in 'proxel.h' as well).
*  Next the black-list is processed:
*
*  \begin{itemize}
*   \item the highest ranking obs. source (kx) is selected. There is an
*         array ('kxrank(kxmax)') which defines this ranking.
*         This array is initialized by routine 'kxname0' and must be updated 
*         every time a new data source is included.
*   \item all data in the black-list with other sources (kx) are marked 
*         for elimination
*   \item the selected data (with highest kx) is then averaged
*         horizontally level by level. If the observations are too
*         close in the vertical, they are vertically averaged as well
*         (this vertical distance is done in log-p, see parameter DELLNP)
*         The following quantities are averaged:
*         \begin{itemize}
*               \item horizontal position (lon, lat)
*               \item level (from log-p average)
*               \item innovations (del)
*               \item forecast error std (sigF)
*               \item obs. uncorrelated error std (sig_Ou)  
*         \end{itemize}
*         Since averaging longitudes is not a well posed problem,
*         the following algorithm is adopted. Spherical coordinates
*         are converted to cartesian coordinates on the unit sphere
*         and the average of the (x,y,z) components is computed.
*         Since the average position is not necessarily on the surface
*         of the sphere, we DEFINE the (lon,lat) of the average position as 
*         the longitude/latitude of the closest point over the sphere
*         to the average (x,y,z) vector. This is accomplished by
*         simply normalizing this vector to have L-2 norm 1.
*  \end{itemize}
*
*  \item After the superobing step is completed the master observation data
*  base is updated, marking the eliminated observations and replacing
*  them with a single superobed profile. Finally, the data is resorted\
*  to remove gaps, the total number of obs adjusted and the region
*  and type/region pointers reset. 
*
*  \end{enumerate}
*
*  One of the limitations of this approach is that data thar are
*  horizontally close but vertically apart would be mistakenly
*  discarded.  (The consideration of horizontal distance only is
*  important to preserve profiles). The fix is to reexamine the
*  observations marked for elimination after the superobing is done and 
*  "release" them if they are vertically distant from the superob
*  (TO DO). However, this should not present a problem for the
*  current observing system. To be kept in mind.
*
*  {\bf NOTE:} Chordal distances are used to defined proximity. 
*
* !CALLING SEQUENCE:
*
*      call PROXEL ( verbose, luverb,  nnobs,   kx,    kt,   kl,
*     .              rlats,   rlons,   rlevs,   del,   sig_Oc, sig_Ou,  sigF,
*     .              tstamp,  maxreg,  iregbeg, ireglen, 
*     .              ktmax, ityplen, nobs                                   )
*
*
* !INPUT PARAMETERS:
*
*  logical      verbose               -> verbose flag
*  integer      luverb                -> output unit number for msg
*  integer      nnobs                 -> number of observations
*  integer      kx(nnobs)             -> observation data source
*  integer      kt(nnobs)             -> observation type
*  real         rlats(nnobs)          -> latitude of obs in deg
*  real         rlons(nnobs)          -> longitude of obs in deg
*  real         rlevs(nnobs)          -> levels in mb
*  real         del(nnobs)            -> innovations
*  real         sig_Oc(nnobs)         -> obs. error std (correlated part)
*  real         sig_Ou(nnobs)         -> obs. error std (uncorrelated part)
*  real         sigF(nnobs)           -> forecast error std.
*  integer      maxreg                -> number of regions
*  integer      iregbeg(maxreg)       -> Pointer to first ob in region
*  integer      ireglen(maxreg)       -> No. obs in region
*  integer      ktmax                 -> No. of obs. types
*  integer      ityplen(ktmax,maxreg) -> No. of obs. of type in region
*
* !OUTPUT PARAMETERS:
*
*  integer      nnobs                 -> number of observations
*  integer      kx(nnobs)             -> observation data source
*  integer      kt(nnobs)             -> observation type
*  logical      kl(nnobs)             -> logical flag
*  real         rlats(nnobs)          -> latitude of obs in deg
*  real         rlons(nnobs)          -> longitude of obs in deg
*  real         rlevs(nnobs)          -> levels in mb
*  real         del(nnobs)            -> innovations
*  real         sig_Oc(nnobs)         -> obs. error std (correlated part)
*  real         sig_Ou(nnobs)         -> obs. error std (uncorrelated part)
*  real         sigF(nnobs)           -> forecast error std.
*  real         tstamp(nnobs)         -> time stamp (not used)
*  real         qcflag(nnobs)         -> qc flag (not used)
*  integer      iregbeg(maxreg)       -> Pointer to first ob in region
*  integer      ireglen(maxreg)       -> No. obs in region
*  integer      ityplen(ktmax,maxreg) -> No. of obs. of type in region
*
* !SEE ALSO: proxel.h
*
* !SYSTEM ROUTINES: atan2
*
* !FILES USED: none
*
* !WRITTEN BY: Arlindo da Silva, 20sep94
* 
* !REVISION HISTORY:
*
C  02Feb95  - Jing G.	- Changed CRAY to _UNICOS for consistency and
C			  to follow the guide lines.
*  08Nov94  A. da Silva  Radiosondes made "transparent".
*  08mar96  G.Brin included sig_Oc and sig_Ou in the args,
*                replaced sigO_lst with sigOc_lst, sigOu_lst
*                ( see this change also in proxel.h)
* !EOP
*.........................................................................

      subroutine PROXEL ( verbose, luverb,  nnobs,   kx,    kt,   kl,
     .                  rlats, rlons, rlevs, del, sig_Oc, sig_Ou, sigF,
     .                  tstamp,  maxreg,  iregbeg, ireglen, 
     .                  ktmax1, ityplen,  nprox )

        use m_Spherical_Partition, only : GetRegion
        use m_spherical_Triangle,  only : 
     .          cossepang => CosSeparation, SEPANG_MIN
        use m_spherical_geometry, only : ll2xyz

      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      use m_stdio,only : stdout
      include 'ktmax.h'

      logical      verbose               ! verbose flag
      integer      luverb                ! output unit number for msg
      integer      nnobs                 ! number of observations
      integer      kx(nnobs)             ! observation data source
      integer      kt(nnobs)             ! observation type
      logical      kl(nnobs)             ! logical flag
      real         rlats(nnobs)          ! latitude of obs in deg
      real         rlons(nnobs)          ! longitude of obs in deg
      real         rlevs(nnobs)          ! levels in mb
      real         del(nnobs)            ! innovations
      real         sig_Oc(nnobs)         ! obs. error std (corr. part)
      real         sig_Ou(nnobs)         ! obs. error std (uncorr. part)
      real         sigF(nnobs)           ! forecast error std.
      real         tstamp(nnobs)         ! time stamp (not used)
c-    real         qcflag(*)             ! qc flag (not used)
      integer      maxreg                ! number of regions
      integer      iregbeg(maxreg)       ! Pointer to first ob in region
      integer      ireglen(maxreg)       ! No. obs in region
      integer      ktmax1                ! No. of obs. types
      integer      ityplen(ktmax,maxreg) ! No. of obs. of type in region

!-----------------------------------------------------------------------
*     Local work space
*     ----------------
      logical,allocatable :: tag(:) ! like kl but not reset after super
      real,   allocatable :: xobs(:),yobs(:),zobs(:)	! x-y-z
      integer,allocatable :: ireglen1(:)	! Local copy of ireglen


*     Local data structure
*     --------------------
      include 'kxmax.h'
      include 'proxel.h'

*     Data source rank (kxrank) defined in here
*     -----------------------------------------
      include 'kxtabl.h'
      include 'kttabl.h'

      character*2 ch
      logical did_it
      integer ierr

	character*6 myname
	parameter(myname='proxel')

*     Statement function: chordal distance squared
*     --------------------------------------------
      dist2(i,j) = (xobs(i)-xobs(j))**2 + 
     .             (yobs(i)-yobs(j))**2 + 
     .             (zobs(i)-zobs(j))**2 

!-----------------------------------------------------------------------

	allocate(tag(nnobs),xobs(nnobs),yobs(nnobs),zobs(nnobs),
     &	  ireglen1(maxreg), stat=ierr)
	if(ierr.ne.0) call die(myname,'allocate()',ierr)

	if(mall_ison()) then
	  call mall_mci(tag ,myname)
	  call mall_mci(xobs,myname)
	  call mall_mci(yobs,myname)
	  call mall_mci(zobs,myname)
	  call mall_mci(ireglen1,myname)
	endif

*     Copy region pointer arrays. Necessary to correctly execute
*     the main region loop.

      ireglen1(:) = ireglen(:)

*     Nothing to do if Rkm = 0. and dellnp = 0.
*     -----------------------------------------
      if ( Rtresh2 .eq. 0. .and. dellnp .eq. 0. ) then
     
         nprox = 0
         
         if ( verbose ) then
            write(luverb, *) ' proxel: '
            write(luverb, *) ' proxel: ', nprox,' too close obs found'
            write(luverb, *) ' proxel: nobs reset from ', nnobs,
     .           ' to ', nnobs
            if ( luverb .ne. stdout ) then
               write(stdout, *) ' proxel: '
               write(stdout, *) 
     .              ' proxel: ', nprox,' too close obs found'
               write(stdout, *) ' proxel: nobs reset from ', nnobs,
     .              ' to ', nnobs
            end if
         end if

		if(mall_ison()) then
		  call mall_mco(tag ,myname)
		  call mall_mco(xobs,myname)
		  call mall_mco(yobs,myname)
		  call mall_mco(zobs,myname)
		  call mall_mco(ireglen1,myname)
		endif
	 deallocate(tag,xobs,yobs,zobs,ireglen1)
         return

      end if


*     Initialize counter of performance statistics
*     --------------------------------------------
      do 5 k = 1, kxmax
         kxelim(k) = 0
 5    continue


*     Set flags to .true. (obs is OK)
*     -------------------------------
      do 10 i = 1, nnobs
         kl(i) = .true.
         tag(i) = .true.
 10   continue


*     Compute cartesian (x,y,z) coord. on unity sphere:
*            x = cos(lat) * cos(lon)
*            y = cos(lat) * sin(lon)
*            z = sin(lat)
*     ------------------------------------------------
      call LL2XYZ ( rlons, rlats, nnobs, xobs, yobs, zobs, ierr )
      if(ierr.ne.0) call die(myname,'LL2XYZ()',ierr)

      nobs = nnobs
      nprox = 0


*     Loop over regions
*     -----------------
      do 20 ireg1 = 1, maxreg

*        Skip if no data in region
*        --------------------------
         if( ireglen(ireg1).eq.0 ) goto 20

*        Loop over obs in reference region 1
*        -----------------------------------
         ibeg1 = iregbeg(ireg1)
         ilen1 = ireglen(ireg1)
         do 30 n1 = ibeg1, ibeg1+ilen1-1

*           Radiosondes are transparent
*           ---------------------------
!           if ( kx(n1) .eq. 7 ) go to 30
            if ( kxrank(kx(n1)) .le. 0 ) go to 30

*           Skip if obs. is black-listed
*           ----------------------------
            if ( .not. tag(n1) ) go to 30

*           Initialize list counter
*           -----------------------
            ilist = 0       

*           Look for "close" obs in this and other regions
*           ----------------------------------------------
            do 35 ireg2 = ireg1, maxreg

*             Skip if no data in region 2
*             ---------------------------
              if( ireglen(ireg2).eq.0 ) goto 35

*              Skip if regions are too far apart
*              (SEPLIM defined in proxel.h)
*              ----------------------------------            
               if( cossepang(GetRegion(ireg1),GetRegion(ireg2),
     .                       SEPANG_MIN,cosSEPLIM )) go to 35

*              Compute obs. range in this region
*              ---------------------------------
               if ( ireg2 .eq. ireg1 ) then
                    n2beg = n1+1
                    n2end = ibeg1+ilen1-1
               else
                    n2beg = iregbeg(ireg2)
                    n2end = n2beg + ireglen(ireg2) - 1 
               end if

*            Loop over obs. in region 2
*            --------------------------
             do 40 n2 = n2beg, n2end

*              Radiosondes are transparent
*              ---------------------------
!              if ( kx(n2) .eq. 7 ) go to 40
               if ( kxrank(kx(n2)) .le. 0 ) go to 40

*              Skip if second obs. is not of same type
*              ---------------------------------------
               if ( kt(n1) .ne. kt(n2) ) go to 40

*              Skip if second obs. has already been black-listed
*              --------------------------------------------------
               if ( .not. tag(n2) ) go to 40

*              Add close enough observation of same type to black-list
*              -------------------------------------------------------
               if ( dist2(n1,n2) .lt. RTRESH2 ) then


*                 Add first obs to list if necessary
*                 ----------------------------------
                  if ( tag(n1) ) then
                     ilist = ilist + 1
		     if (ilist.le.nlist) then	! keep acounting
                        idx_lst(ilist)    = n1
                        kx_lst(ilist)     = kx(n1)
                        kt_lst(ilist)     = kt(n1)
                        ireg_lst(ilist)   = ireg1
                        kl_lst(ilist)     = .true.
                        rlats_lst(ilist)  = rlats(n1)
                        rlons_lst(ilist)  = rlons(n1)
                        rlevs_lst(ilist)  = rlevs(n1)
                        x_lst(ilist)      = xobs(n1)
                        y_lst(ilist)      = yobs(n1)
                        z_lst(ilist)      = zobs(n1)
                        del_lst(ilist)    = del(n1)
                        sigOc_lst(ilist)  = sig_Oc(n1)
                        sigOu_lst(ilist)  = sig_Ou(n1)
                        sigF_lst(ilist)   = sigF(n1)
                        tag(n1)           = .false.  
		     endif
                  end if

*                 Add second obs. to black list
*                 -----------------------------
                  ilist = ilist + 1
		  if (ilist.le.nlist) then	! keep acounting
                     idx_lst(ilist)    = n2
                     kx_lst(ilist)     = kx(n2)
                     kt_lst(ilist)     = kt(n2)
                     ireg_lst(ilist)   = ireg2
                     kl_lst(ilist)     = .true.
                     rlats_lst(ilist)  = rlats(n2)
                     rlons_lst(ilist)  = rlons(n2)
                     rlevs_lst(ilist)  = rlevs(n2)
                     x_lst(ilist)      = xobs(n2)
                     y_lst(ilist)      = yobs(n2)
                     z_lst(ilist)      = zobs(n2)
                     del_lst(ilist)    = del(n2)
                     sigOc_lst(ilist)  = sig_Oc(n1)
                     sigOu_lst(ilist)  = sig_Ou(n1)
                     sigF_lst(ilist)   = sigF(n2)
                     tag(n2)           = .false.  
		  endif

               end if

 40          continue

 35         continue

*           No cluster found, go on...
*           --------------------------
            if ( ilist .lt. 2 ) go to 30
            if ( ilist .gt. nlist ) then
	       write(luverb,'(1x,3a,i5,a)') 'proxel: ',
     &		  'Not enough work space to hold black-list.  ',
     &		  'nlist =',ilist,' is expected.'
	       if (stdout.ne.luverb) then
	          write(stdout,'(1x,3a,i5,a)') 'proxel: ',
     &		     'Not enough work space to hold black-list.  ',
     &		     'nlist =',ilist,' is expected.'
	       endif
	       call die(myname)
            end if


*           Now the black-list has a cluster of close enough 
*            observations.  First find the highest ranking observation
*            source (all obs in list assumed to have the same type kt)
*           -----------------------------------------------------------
            kx_hi = IBIG
            kxrank_hi = IBIG
            do 100 i = 1, ilist
               if ( kxrank ( kx_lst(i) ) .lt. kxrank_hi ) then
                    kx_hi = kx_lst(i)
                    kxrank_hi = kxrank ( kx_lst(i) ) 
               end if
 100        continue
            if(kx_hi.eq.IBIG) call die(myname,'invalid kx_hi',kxhi)

*           Now we "superob" profiles of highest ranking source kx_hi
*           ---------------------------------------------------------
            call PRXSOB ( ierr, kx_hi, did_it, tag, nnobs )
            if(ierr.ne.0) call die(myname,'PRXSOB()',ierr)

*           If no elimination has been done, go on...
*           -----------------------------------------
            if ( .not. did_it ) go to 30


*           Update master list with superobed black-list
*           --------------------------------------------
            if ( debug ) write(luverb,1100)  
            do 150 i = 1, ilist

               n = idx_lst(i)

*              Superobed observation
*              ---------------------
               if ( kl_lst(i) ) then

                    kl(n) = .true.

*                   Write result of superobing
*                   --------------------------
                    if ( debug ) then
                     ch = '  '
                     if ( del_lst(i) .ne. del(n) ) ch = ' o'
                     write(luverb,1100) n, ktname(kt_lst(i)),
     .                kx_lst(i), kxdesc(kx_lst(i)), 
     .                kxrank ( kx_lst(i) ), rlons_lst(i), rlats_lst(i),
     .                rlevs_lst(i), del_lst(i), sigF_lst(i), ch
 1100                 format(' proxel: ',i6,a9,i4,1x,a25,i7,
     .                       2F8.1,F8.0,2F8.3,a2)
                     if ( ch .eq. ' o' )
     .                write(luverb,1100) n, ktname(kt(n)),
     .                 kx(n), kxdesc(kx(n)), 
     .                 kxrank ( kx(n) ), rlons(n), rlats(n),
     .                 rlevs(n), del(n), sigF(n), ' *'
                    end if

*                   update master list
*                   ------------------
                    rlons(n) = rlons_lst(i)
                    rlats(n) = rlats_lst(i)
                    rlevs(n) = rlevs_lst(i)
                    del(n) = del_lst(i)
                    sig_Ou(n) = sigOu_lst(i)
                    sigF(n) = sigF_lst(i)

*              Eliminated observation
*              ----------------------
               else

                    kl(n) = .false.

*                   Write result of black listing
*                   -----------------------------
                    if ( debug ) then
                     write(luverb,1100) n, ktname(kt_lst(i)),
     .               kx_lst(i), kxdesc(kx_lst(i)), 
     .               kxrank ( kx_lst(i) ), rlons_lst(i), rlats_lst(i),
     .               rlevs_lst(i), del_lst(i),
     .               sigF_lst(i), ' *'

                    end if

*                   Decrement region pointers
*                   -------------------------
                    ireg = ireg_lst(i)
                    kti = kt_lst(i)
                    ireglen1(ireg) = ireglen1(ireg) - 1
                    ityplen(kti,ireg) = ityplen(kti,ireg) - 1
                    nprox = nprox + 1
                    kxelim(kx(n)) = kxelim(kx(n)) + 1


               end if 

 150        continue   ! end of updating loop

 30      continue      ! end of reference obs. loop

 20   continue         ! end of region loop



*     Slide the superobed data to the front of the arrays
*     ---------------------------------------------------
      call TOFRONT ( nnobs, kx,    kt,    kl,   
     .                 rlats, rlons, rlevs,
     .                 del, sig_Ou, sig_Oc,  sigF,  tstamp, nobs  )


      if ( verbose ) then
          write(luverb, *) ' proxel: '
          write(luverb, *) ' proxel: ', nprox,' too close obs found'
          write(luverb, *) ' proxel: nobs reset from ', nnobs,
     .                     ' to ', nobs
          if ( luverb .ne. stdout ) then
            write(stdout, *) ' proxel: '
            write(stdout, *) ' proxel: ', nprox,' too close obs found'
            write(stdout, *) ' proxel: nobs reset from ', nnobs,
     .                     ' to ', nobs
          end if
      end if


*     No clusters found, all done
*     ---------------------------
      if ( nprox .eq. 0 ) then
		if(mall_ison()) then
		  call mall_mco(tag ,myname)
		  call mall_mco(xobs,myname)
		  call mall_mco(yobs,myname)
		  call mall_mco(zobs,myname)
		  call mall_mco(ireglen1,myname)
		endif
	deallocate(tag,xobs,yobs,zobs,ireglen1)
	return
      endif


*     Reset number of observations
*     ----------------------------
      nnobs = nobs


*     Save new region lengths in permanent location.
*     Reset region pointers.
*     ----------------------------------------------

      ksofar = 0
      do 600 ireg = 1, maxreg
         ireglen(ireg) = ireglen1(ireg)
         if( ireglen(ireg).gt.0 ) then
            iregbeg(ireg) = ksofar + 1
            ksofar = ksofar + ireglen(ireg)
         endif
  600 continue


*     Report elimination statistics
*     -----------------------------
      if ( verbose .and. nprox .gt. 0 ) then
         write(luverb,901)
         if ( luverb .ne. stdout ) write(stdout,901)
 901     format('  proxel: ',/,'  proxel: ','  KX ',2x,
     .   '  DATA SOURCE (kxdesc)  ',2x,' No.')
         do 900 k = 1, kxmax
            if ( kxelim(k) .ne. 0 ) then
               write(luverb,902) k, kxdesc(k), kxelim(k)
               if ( luverb .ne. stdout ) then
                    write(stdout,902) k, kxdesc(k), kxelim(k)
               end if
            end if
 900     continue
 902     format(1x,' proxel: ',i4,2x,a25,2x,i4)
         write(luverb,*) ' proxel: '
         if (luverb.ne.stdout) write(stdout,*) ' proxel: '
      end if

		if(mall_ison()) then
		  call mall_mco(tag ,myname)
		  call mall_mco(xobs,myname)
		  call mall_mco(yobs,myname)
		  call mall_mco(zobs,myname)
		  call mall_mco(ireglen1,myname)
		endif
      deallocate(tag,xobs,yobs,zobs,ireglen1)
      return
      end


*........................... ROUTINE PROLOGUE .............................
* !BOP
*
* !FILE: proxel.f
*
* !ROUTINE: prxsob
* 
* !DESCRIPTION: 
*
*    Superobs observation with type KX_HI discarding all other 
* observations on black-list.
*
* !CALLING SEQUENCE: call PRXSOB ( ierr, kx_hi,did_it,tag,nnobs )
*
* !INPUT PARAMETERS: 
*
*    kx_hi - observation type to superob.
*
* !OUTPUT PARAMETERS: 
*
*    ierr  error condition. All is fine if ierr=0 on return.
*      *   black-list data structure (see proxel.h)
*
* !SEE ALSO: proxel.h
*
* !SYSTEM ROUTINES: none
*
* !FILES USED: none
*
* !WRITTEN BY: Arlindo da Silva, 26sep94
* 
* !REVISION HISTORY:
*
* !EOP
*.........................................................................

      subroutine PRXSOB ( ierr, kx_hi, did_it, tag, nnobs )

      integer ierr, kx_hi
      logical tag(nnobs)

      logical did_it, vclose

*     Local data structure
*     --------------------
      include 'kxmax.h'
      include 'proxel.h'

      parameter ( pi     = 3.141592654 )
      parameter ( pid2   = pi / 2. )
      parameter ( r2d = 180.0 / pi )


*     Eliminate other obs. types
*     --------------------------
      did_it = .false.
      do 10 i = 1, ilist
         kl_lst(i) = ( kx_lst(i) .eq. kx_hi )  
         if ( .not. kl_lst(i) ) did_it = .true. !(?)
 10   continue


*     Average close enough obs
*     ------------------------
      do 20 i1 = 1, ilist

         if ( .not. kl_lst(i1) ) go to 20

         rlnp1 = alog ( rlevs_lst(i1) )

         xmean = x_lst(i1)
         ymean = y_lst(i1)
         zmean = z_lst(i1)
         pmean = rlnp1
         dmean = del_lst(i1) 
         fmean = sigF_lst(i1)
         n = 1

         do 30 i2 = i1+1, ilist

            if ( .not. kl_lst(i2) ) go to 30

            rlnp2 = alog ( rlevs_lst(i2) )            
            if ( abs(rlnp2-rlnp1) .lt. DELLNP ) then

               xmean = xmean + x_lst(i2)
               ymean = ymean + y_lst(i2)
               zmean = zmean + z_lst(i2)

               pmean = pmean + rlnp2
               dmean = dmean + del_lst(i2) 
               fmean = fmean + sigF_lst(i2)
               kl_lst(i2) = .false.
               n = n + 1

            end if

 30      continue

*        Assign average to i1-th obs. in black-list
*        ------------------------------------------
         if ( n .gt. 1 ) then

              did_it = .true.

*             From mean (x,y,z) cartesian coordinates,
*              compute nearest point on the unit sphere
*              by normalizing the mean vector to 1
*              (division by n folded into rnorm). 
*             ------------------------------------------
              rnorm = sqrt ( xmean*xmean + ymean*ymean 
     .                       + zmean*zmean )
              if ( rnorm .eq. 0. ) then       ! very unlikely
                 ierr = 99
                 return
              end if

              xmean = xmean / rnorm
              ymean = ymean / rnorm
              zmean = zmean / rnorm
            
*             Get lat/lon from (x,y,z). Recall:
*                    x = cos(lat) * cos(lon)
*                    y = cos(lat) * sin(lon)
*                    z = sin(lat)
*             --------------------------------
              rlats_lst(i1) = asin ( zmean )
              rlons_lst(i1) = atan2 ( ymean, xmean )

*             lat/lon in degrees
*             ------------------
              rlats_lst(i1) = r2d * rlats_lst(i1)
              rlons_lst(i1) = r2d * rlons_lst(i1)

              rlevs_lst(i1) = exp ( pmean / n )
              del_lst(i1) = dmean / n
              sigF_lst(i1) = fmean / n

* obs of the same data type have the same unc. error std. ( we add variances,
* n*var/n**2) ===> std =  std/sqrt(n)

              sigOu_lst(i1) = sigOu_lst(i1) / sqrt(float(n))

         end if

 20   continue

      if ( .not. did_it ) return


*     Further check: if obs with kx .ne. kx_hi is not
*        vertically close to superob just let it in
*        (this is to avoid elimination of cloud track
*         winds by aircraft reports at different levels)
*     --------------------------------------------------
      did_it = .false.
      do 150 i1 = 1, ilist

         if ( kx_lst(i1) .ne. kx_hi ) then

            rlnp1 = alog ( rlevs_lst(i1) )
            vclose = .false.
            
            do 160 i2 = 1, ilist
               
               if ( .not. kl_lst(i2) ) go to 160
               
               rlnp2 = alog ( rlevs_lst(i2) )
               
               if ( abs(rlnp2-rlnp1) .le. dellnp )  then
                  vclose = .true.
               end if
               
 160        continue

*           If not vertically close, release observation
*           --------------------------------------------
            if ( .not. vclose ) then
               n = idx_lst(i1)
               kl_lst(i1) = .true.   ! release it
               tag(n) = .true.       ! untag it
            end if

         end if

         if ( .not. kl_lst(i1) ) did_it = .true.

 150  continue

      return
      end
* !FILE: proxel.f
* !subroutine prx_tofront
*
C  02Feb95  - Jing G.	- Changed CRAY to _UNICOS for consistency and
C			  to follow the guide lines.
c  03oct94 - Implemented CRAY specifics with IFDEFs.
c  modification for dynamic storage on CRAY                     05/28/93
c
*  11mar96 - G.Brin made a copy of TOFRONT suitable for proxel .
*            i.e. uses an extra argument - sig_Ou.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.... Slide the selected data to the front of the data arrays.

      subroutine PRX_TOFRONT ( nobs,  kx,    kt,    kl,   
     $                     rlats, rlons, rlevs,
     $                     del, sig_Oc, sig_Ou,  sigF,  tstamp, newnr)

c.......................................................................
c.... Argument declarations.

      integer      nobs
      integer      kx(nobs)
      integer      kt(nobs)
      logical      kl(nobs)
      real         rlats(nobs)
      real         rlons(nobs)
      real         rlevs(nobs)
      real         del(nobs)
      real         sig_Oc(nobs)
      real         sig_Ou(nobs)
      real         sigF(nobs)
      real         tstamp(nobs)
c-    real         qcflag(nobs) 	! no longer used
      integer      newnr

c.......................................................................
c.... Local storage

      integer      iperm(nobs)   ! dynamic allocation
      real         rsort(nobs)

c.......................................................................
c.... Initialize the vector to be sorted

      newnr = 0

      do 100 n = 1, nobs

         if( kl(n) ) then
            newnr    = newnr + 1
            rsort(n) = - float(newnr)
         else
            rsort(n) = float(n)
         endif

  100 continue

      if( newnr.gt.0 ) then

c.......................................................................
c....... Perform the sort.

         flip = - float(newnr)
         do 110 n = 1, nobs
            if( kl(n) ) rsort(n) = flip - rsort(n)
  110    continue

         call INDEXXR ( nobs, rsort, iperm )

c.......................................................................
c....... Apply the sorting permutation to all arrays.

         call PERMUTI ( kx,    iperm, nobs, kx    )
         call PERMUTI ( kt,    iperm, nobs, kt    )
         call PERMUTL ( kl,    iperm, nobs, kl    )
         call PERMUTR ( rlats, iperm, nobs, rlats )
         call PERMUTR ( rlons, iperm, nobs, rlons )
         call PERMUTR ( rlevs, iperm, nobs, rlevs )
         call PERMUTR ( del,   iperm, nobs, del   )
         call PERMUTR ( sig_Oc,iperm, nobs, sig_Oc)
         call PERMUTR ( sig_Ou,iperm, nobs, sig_Ou)
         call PERMUTR ( sigF,  iperm, nobs, sigF  )
         call PERMUTR ( tstamp,iperm, nobs, tstamp)
c-       call PERMUTR ( qcflag,iperm, nobs, qcflag)

      endif

      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end
