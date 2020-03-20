!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOI
!
! !TITLE: An Application for Gridding Innovation\\ Time Averaged Statistics (GrITAS v2.00)
!
! !AUTHORS: Arlindo da Silva and Guang Ping Lou
!
! !AFFILIATION: Data Assimilation Office, NASA/GSFC, Greenbelt, MD 20771
!
! !DATE: September 1, 1997 (Revised June 9, 2004)
!
! !INTRODUCTION: System Overview
!
!  GrITAS is a FORTRAN 77/90 program which reads ODS or GEOS-1 style
!  del files and produces global grids of O-F (observation minus forecast),
!  O-A (observation minus analysis) or  O (observed value) with time mean,
!  standard deviation and number of observations in each grid box.
!  Grids are produced for user specified combination of data types
!  (e.g., u, v-winds, heights) and data sources (e.g., radiosondes,
!  TOVS retrievals, ships, etc.) 
!
!
!EOI
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GrITAS()
! 
! !DESCRIPTION: Driver for the {\em Gridding of Innovation Time Average
!               Statistics} system (GrITAS).
!
! !INTERFACE:
!
       Program GrITAS

! !USES:

       Use m_odsmeta, only : KTMAX
       Use m_odsmeta, only : KXMAX
       Use m_odsmeta, only : X_PRE_BAD  ! GSI  passive data
       Use m_odsmeta, only : X_PASSIVE  ! PSAS passive data

       Use m_gritas_grids, only : im, jm, km
       Use m_gritas_grids, only : GrITAS_Pint
       Use m_gritas_grids, only : GrITAS_Accum
       Use m_gritas_grids, only : GrITAS_Norm

       Use m_ods

!      Use m_obs, only : Reducer

       Use m_stdio, only : stdout
       Use m_die

! !REMARKS:
!
!   1) In case of calculating statistics of (specified) sigo and obias the standard deviation 
!      output is replaced by root-mean-square; that is to say that position t=1 in the grads
!      file will have the mean as it normally has, and t=2 will have the RMS instead of the STDV.
!
!   2) In case of Desroziers et al. diagnostics, the mean and standard deviations output 
!      are replaced by cross-variance and the cross standard deviations, respectively.
!
! !REVISION HISTORY: 
!
!  01Sep97  da Silva   Initial design & code.
!  17Oct97  G.P. Lou   Modifications. (RT: What modifications?)
!  15Mar98  da Silva   Minor touch ups (comments only)
!  11Dec98  da Silva   Included ODS input along with O-A/obs option.
!  15Dec98  da Silva   Added GFIO output.
!  23dec00  da Silva   Made part of fvDAS.
!  09Jun04  Todling    - Using m_ods
!                      - Resolution independent executable
!                      - Changed action when no-obs found in syn time
!  24Jun04  Todling    Fixed time frequency of output (RUC compatible).
!  30Jul04  Todling    Added reducer capability
!  12oct05  da Silva   Clean ODS each time
!  07Nov06  Todling    - Implemented Desroziers et al. diagnostics
!                      - Implemented Jo option.
!  23Jan08  Todling    Added -obsimp
!  12Sep08  Todling    Added -nonorm to facilitate obs impact calculation
!  07Jul09  Redder     Added ob_flag as an argument in the call to the
!                      routine GrITAS_Accum
!  20Jul09  Redder     Added the code to generate the descriptions for
!                      for each entry in the kt list.
!  21Jul09  Redder     Fixed bug in call to the routine, GrITAS_Init.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      Implicit NONE

      character(len=*), parameter :: myname = 'GrITAS'

      include 'gritas.h'


!     User defined vertical levels
!     ----------------------------
      real, allocatable :: p1(:), p2(:)  ! pressure intervals for gridding


!     Workspace for Holding observations for each synoptic time
!     Note: lat/lon/etc cannot be pointers because ods is r8 vs lat_r4
!     ----------------------------------------------------------------
      integer                  nobs      ! actual number of obs for this synoptic time
      integer                  nobs_good ! no. of observations after reducer operation
      type     ( ods_vect ) :: ods
      integer, allocatable ::  kx (:)
      integer, allocatable ::  kt (:)
      real,    allocatable ::  lat(:)
      real,    allocatable ::  lon(:)
      real,    allocatable ::  lev(:)
      real,    allocatable ::  qc (:)          ! not used, but needed
      real,    allocatable ::  del(:,:)
      integer, allocatable ::  ob_flag (:)     ! = 0 if ob was used to accumulate
                                               !   bias/RMS, = nonzero otherwise

!     Workspace for holding the surface gridded fields
!     ------------------------------------------------
      real, allocatable ::  bias_2d(:,:,:,:)   ! first  time mean
      real, allocatable ::  rtms_2d(:,:,:,:)   ! root time mean square
      real, allocatable ::  stdv_2d(:,:,:,:)   ! standard deviation
      real, allocatable ::  xcov_2d(:,:,:,:)   ! x-rms and x-variance
      real, allocatable ::  nobs_2d(:,:,:)     ! number of obs per grid box
      integer list_2d(l2d_max)                 ! corresponding list item

!      equivalence ( rtms_2d(1,1,1), stdv_2d(1,1,1) )


!     Workspace for holding the upper-air gridded fields
!     --------------------------------------------------
      real, allocatable ::  bias_3d(:,:,:,:,:) ! first  time mean
      real, allocatable ::  rtms_3d(:,:,:,:,:) ! root time mean square
      real, allocatable ::  stdv_3d(:,:,:,:,:) ! standard deviation
      real, allocatable ::  xcov_3d(:,:,:,:,:) ! x-standard deviation
      real, allocatable ::  nobs_3d(:,:,:,:)   ! number of obs per grid box
      integer list_3d(l3d_max)                 ! corresponding list item


!     Data structure for user selection of data type/sources
!     ------------------------------------------------------
      integer       listsz               ! size of the list (< lmax)
      integer       kt_list(lmax)        ! each grid has one data type (kt)  
      integer       kx_list(kxmax,lmax)  ! but several data sources (kx)
      character*11  var_list(lmax)       ! variable name for output 
                                         ! (e.g., uraob)
      character*256 var_descr(lmax)      ! ... and its description.
!     Note: Negative kt values will be treated as surface variables

!     Data type/data source table. Gives the grid number for each (kx,kt)
!     This table is derived from the user lists above. Notice that 
!     surface grids are associated with kt<4.
!     -------------------------------------------------------------------
      integer       kxkt_table(kxmax,ktmax)

      integer         l2d                      ! actual number of 2d grids
      integer         l3d                      ! actual number of 3d grids

!     Input File names
!     ----------------
      integer          ndel                      ! actual number of del-files
      character*(255)  del_ifname(NDEL_MAX)      ! O-F (del) data files

!     Output File names
!     -----------------
      character*(255)  grid_ofname               ! gridded O-F (del) file

!     Local variables
!     ---------------
      integer idel, isyn, ksyn, ioptn, iopt, nstat, nr, i, j, ir
      integer nymd, nhms, timinc, nsyn, ier, fid, rc
      integer nymd_min, nymd_max, nymd_out
      character*8 ftype, dtype, otype
      character(len=255) RC_red, lb_tallyfile
      logical DO_red, first
      logical ncf, split_tgroups
      logical iampassive, iwantpassive, nrmlz
 
      nymd_min =  99999999
      nymd_max = -99999999
      nymd_out = -1         !  date of output

      first = .true.
      nrmlz = .true.

!...........................................................................


!     User interface
!     --------------
      call GrITAS_Init ( NDEL_MAX, LMAX, KXMAX,
     &                   ftype, dtype, otype, ioptn, ncf, nr,
     &                   split_tgroups, lb_tallyfile,
     &                   nymd_out, DO_red, RC_red,
     &                   ndel, del_ifname, grid_ofname, 
     &                   kt_list, kx_list, var_list, var_descr,
     &                   listsz, nrmlz )

      iopt = abs(ioptn)

!     Allocate space for working arrays
!     ---------------------------------
      allocate ( bias_2d(im,jm,l2d_max,nr),  rtms_2d(im,jm,l2d_max,nr),
     &           stdv_2d(im,jm,l2d_max,nr),
     &           xcov_2d(im,jm,l2d_max,nr),  nobs_2d(im,jm,l2d_max),
     &           stat=ier )
        if(ier/=0) call die (myname,'Error in Alloc(2d)')

      allocate ( bias_3d(im,jm,km,l3d_max,nr), rtms_3d(im,jm,km,l3d_max,nr),
     &           stdv_3d(im,jm,km,l3d_max,nr),
     &           xcov_3d(im,jm,km,l3d_max,nr), nobs_3d(im,jm,km,l3d_max),
     &           stat=ier )
        if(ier/=0) call die (myname,'Error in Alloc(3d)')

      allocate ( p1(km), p2(km),
     &           stat=ier )
        if(ier/=0) call die (myname,'Error in Alloc(pressure)')


!     Create kx-kt table from user defined list
!     -----------------------------------------
      call GrITAS_KxKt ( KTMAX, KXMAX, L2D_MAX, L3D_MAX,
     &                   listsz, kt_list, kx_list, 
     &                   kxkt_table, list_2d, list_3d, l2d, l3d )


!     Set pressure Intervals for gridding
!     -----------------------------------
      call GrITAS_Pint ( p1, p2 )


!     Initialize grids
!     ----------------
      bias_2d(1:im,1:jm,     1:l2d_max,1:nr) = 0.0
      rtms_2d(1:im,1:jm,     1:l2d_max,1:nr) = 0.0
      stdv_2d(1:im,1:jm,     1:l2d_max,1:nr) = 0.0
      xcov_2d(1:im,1:jm,     1:l2d_max,1:nr) = 0.0
      nobs_2d(1:im,1:jm,     1:l2d_max)      = 0.0
      bias_3d(1:im,1:jm,1:km,1:l3d_max,1:nr) = 0.0
      rtms_3d(1:im,1:jm,1:km,1:l3d_max,1:nr) = 0.0
      stdv_3d(1:im,1:jm,1:km,1:l3d_max,1:nr) = 0.0
      xcov_3d(1:im,1:jm,1:km,1:l3d_max,1:nr) = 0.0
      nobs_3d(1:im,1:jm,1:km,1:l3d_max)      = 0.0


!     Figure out range of synoptic time loop
!     --------------------------------------
      ksyn = 32767
      if(ncf) ksyn = 1

!     Loop over del files
!     -------------------
      do idel = 1, ndel

!         Loop over synoptic time on this del file...
!         -------------------------------------------
          do 10 isyn = 1, ksyn     


!            Read del file for this synoptic time
!            ------------------------------------
             if ( trim(ftype) .eq. 'del' ) then

                allocate ( lat(NOBSMAX), lon(NOBSMAX), lev(NOBSMAX),
     &                     qc (NOBSMAX), del(NOBSMAX,nr),
     &                     kx(NOBSMAX),  kt(NOBSMAX),
     &                     ob_flag(NOBSMAX),
     &                     stat=ier )
                   if(ier/=0) call die (myname,
     &                                ' Obs Arrays Error, alloc()')

                call Read_Del ( trim(del_ifname(idel)), NOBSMAX, ioptn,
     &                          lat, lon, lev, del(1,1), kx, kt, qc, nobs,
     &                          nymd, nhms, ier )
                ob_flag ( : nobs ) = 0
             else

                nymd = -1           ! get data for next synoptic time on file
                nhms =  0
                call ODSNxTime ( trim(del_ifname(idel)), nymd, nhms )
                if ( nymd .eq. -1 ) then
!                    print *, 'End-Of-File'
                     exit     ! skip to next file
                end if
                if ( first ) then
                     first = .false.
                else
                   call ods_clean ( ods, ier )
                end if
                call ODS_Get ( trim(del_ifname(idel)), nymd, nhms, ftype, ods, ier, ncf=ncf )

                if ( ier .gt. 0 ) then
                     print *, 'ODS_Get error: ier = ', ier
                     exit     ! skip to next file
                else
                   if ( ods%data%nobs > 0 ) then
                     write(stdout,'(3a,i8,a,i6.6)') 
     &                    '[x] Read ', trim(del_ifname(idel)), ' on ', nymd, ' at ', nhms
                   end if
                end if


!               Set number of observations found in the file
!               --------------------------------------------
                nobs   =  ods%data%nobs

                if ( nobs .eq. 0 ) then
!!!                     print *, 'No data for this synoptic time'
                     cycle    ! skip to next synoptic time
                end if

!               Reduce input data as requested
!               ------------------------------
!               if ( DO_red ) then
!                    call Reducer ( nymd, nhms, nobs, ods, nobs_good, rcfile=RC_red )
!                    ods%data%nobs = nobs_good
!                    nobs          = nobs_good
!                    print *, 'Remaining obs after Reducer: ', nobs
!               end if

                allocate ( lat(nobs), lon(nobs), lev(nobs), qc (nobs), 
     &                     del(nobs,nr),   kx(nobs), kt (nobs),
     &                     ob_flag (nobs), stat=ier )
                   if(ier/=0) call die (myname,
     &                                ' Obs Arrays Error, alloc()')

                lat     = ods%data%lat   (1:nobs)
                lon     = ods%data%lon   (1:nobs)
                lev     = ods%data%lev   (1:nobs)
                kt      = ods%data%kt    (1:nobs)
                kx      = ods%data%kx    (1:nobs)
                qc      = real(ods%data%qcexcl(1:nobs))
                ob_flag ( : nobs ) = 0

!               Force y2K compliance
!               --------------------
                if ( nymd .lt. 19000000 ) nymd = nymd + 19000000
                nsyn = ods%meta%nsyn

!               Collect requested residual (field) (O-F, O-A, obs, etc)
!               -------------------------------------------------------
                nstat = 3     ! writes out mean, std, and nobs count
                if ( iopt .eq. 0 .or. iopt .eq. 1 ) then
                     del(:,1) = ods%data%omf(1:nobs)
                else if ( iopt .eq. 2 .or. iopt .eq. 3 ) then
                     del(:,1) = ods%data%oma(1:nobs)
                else if ( iopt .eq. 4 .or. iopt .eq. 5 ) then
                     del(:,1) = ods%data%obs(1:nobs)
                else if ( iopt .eq. 6 .or. iopt .eq. 7 ) then    ! prescribed sigo
                     del(:,1) = ods%data%xvec(1:nobs)
                else if ( iopt .eq. 8 .or. iopt .eq. 9 ) then    ! observation bias
                     del(:,1) = ods%data%xm(1:nobs)
                else if ( iopt .eq.10 .or. iopt .eq.11 ) then    ! X-cov <AmF,OmF> ~ HBH'
                     del(:,1) = ods%data%oma(1:nobs) - ods%data%omf(1:nobs)
                     del(:,2) = ods%data%omf(1:nobs)
                else if ( iopt .eq.12 .or. iopt .eq.13 ) then    ! X-cov <AmF,OmA> ~ HPaH'
                     del(:,1) = ods%data%oma(1:nobs) - ods%data%omf(1:nobs)
                     del(:,2) = ods%data%oma(1:nobs)
                else if ( iopt .eq.14 .or. iopt .eq.15 ) then    ! X-cov <OmA,OmF> ~ R
                     del(:,1) = ods%data%oma(1:nobs)
                     del(:,2) = ods%data%omf(1:nobs)
                else if ( iopt .eq.16 .or. iopt .eq.17 ) then    ! Jo(OmF)/p
                     del(:,1) = ods%data%omf(1:nobs)
                     del(:,2) = ods%data%xvec(1:nobs)
                     nstat = 1
                else if ( iopt .eq.18 .or. iopt .eq.19 ) then    ! Jo(OmA)/p
                     del(:,1) = ods%data%oma(1:nobs)
                     del(:,2) = ods%data%xvec(1:nobs)
                     nstat = 1
                else if ( iopt.eq.20 .or. iopt.eq.21 ) then      ! original observation
                     del(:,1) = ods%data%obs(1:nobs) + ods%data%xm(1:nobs)
                else if ( iopt.eq.22 .or. iopt.eq.23 ) then      ! observation impacts
                     del(:,1) = ods%data%xvec(1:nobs)
                else
                     call die (myname, 'iopt not valid')
                endif

!               Now that all atrributes have been read, eliminate bad obs
!                and those that have passed QC flag, if user desires. 
!                Notice that passive data types are kept if ioptn>0.
!               --------------------------------------------------------
                if ( mod(iopt,2) .ne. 0 ) then
                   do i = 1, nobs
                      iampassive   = ods%data%qcexcl(i) .eq. X_PASSIVE 
     &                          .or. ods%data%qcexcl(i) .eq. X_PRE_BAD
                      iwantpassive = ioptn .gt. 0 .and. iampassive
                      if ( ods%data%qcexcl(i) .eq. 0 .or.
     &                     iwantpassive ) then
                         ob_flag ( i ) = 0
                      else 
                         ob_flag ( i ) = 1
                      end if
                   end do
                   if(iopt==17 .or. iopt==19) del(1:nobs,1) = del(1:nobs,1)/del(1:nobs,2)  ! omf/sigo for Joa/Job calculation
                end if

             end if

!            Enforce Y2K compliance, etc.
!            ----------------------------
             if ( nymd .lt. 19000000 ) nymd = nymd + 19000000
             nymd_max = max ( nymd_max, nymd )
             nymd_min = min ( nymd_min, nymd )

!            Accumulate grids
!            ----------------
             call GrITAS_Accum ( lat, lon, lev, kx, kt, del(1:nobs,1:nr), nobs, nr,
     &                           kxkt_table, p1, p2,
     &                           ob_flag,
     &                           l2d, bias_2d, rtms_2d, xcov_2d(:,:,:,1),   nobs_2d,
     &                           l3d, bias_3d, rtms_3d, xcov_3d(:,:,:,:,1), nobs_3d )
             if ( trim(ftype)        .ne. 'del' .and.
     &            trim(lb_tallyfile) .ne.  ''  .and.
     &            nobs         .gt.  0 ) then

                call GrITAS_LBTally ( del_ifname(idel), nymd, nhms,
     &                                ob_flag, ods,  rc )
             end if

             deallocate ( lat, lon, lev, qc, del, kx, kt, ob_flag )

 10       continue
 11       continue


       end do

!      Normalize grids
!      ---------------
!      if ( nrmlz  ) then
            call GrITAS_Norm ( nr, nrmlz,
     &                    l2d, bias_2d, rtms_2d, stdv_2d, xcov_2d, nobs_2d,
     &                    l3d, bias_3d, rtms_3d, stdv_3d, xcov_3d, nobs_3d)
!      endif


!      STDV not always makes sense ...
!      -------------------------------
       if ( iopt==6 .or.iopt==7 .or. iopt==8 .or.iopt==9 ) then
            stdv_2d = rtms_2d
            stdv_3d = rtms_3d
       endif

!      In case of x-covariances, write x-rms in the bias slot, and x-variances in stdv slot
!      ------------------------------------------------------------------------------------
       if ( iopt==10 .or.iopt==11 .or. iopt==12 .or.iopt==13 .or. iopt==14 .or.iopt==15 ) then
            bias_2d(:,:,:,  1) = xcov_2d(:,:,:,  1)
            stdv_2d(:,:,:,  1) = xcov_2d(:,:,:,  2)
            bias_3d(:,:,:,:,1) = xcov_3d(:,:,:,:,1)
            stdv_3d(:,:,:,:,1) = xcov_3d(:,:,:,:,2)
       endif

!      In case of Jo(OmA) and Jo(OmF) store Jo in bias array
!      -----------------------------------------------------
       if ( iopt==16 .or.iopt==17 .or. iopt==18 .or.iopt==19 ) then
            bias_2d(:,:,:,  1) = rtms_2d(:,:,:,  1)
            bias_3d(:,:,:,:,1) = rtms_3d(:,:,:,:,1)
       endif
       

!      Write del grids to file
!      -----------------------
       if ( nymd_out < 0 ) then
          nymd = ( nymd_min + nymd_max ) / 2
       else
          nymd = nymd_out ! user specified
       endif
       nhms = 120000
       timinc = 240000 / nsyn ! frequency of ouput
       if ( trim(otype) .eq. 'hdf' ) then
          if ( split_tgroups ) nstat = -1 * nstat
          call GFIO_Output   ( grid_ofname, 
     &                         var_list, var_descr, listsz, timinc,
     &                         l2d, list_2d, bias_2d(1,1,1,1)  , stdv_2d(1,1,1,1)  , nobs_2d,
     &                         l3d, list_3d, bias_3d(1,1,1,1,1), stdv_3d(1,1,1,1,1), nobs_3d,
     &                         nstat, nymd, nhms, fid ) 

       else

          grid_ofname = trim(grid_ofname)//'.dat' 
          call GrADS_Output (  trim(grid_ofname), 
     &                         var_list, var_descr, listsz,
     &                         l2d, list_2d, bias_2d(1,1,1,1)  , stdv_2d(1,1,1,1)  , nobs_2d,
     &                         l3d, list_3d, bias_3d(1,1,1,1,1), stdv_3d(1,1,1,1,1), nobs_3d, nstat)

       end if

!      Clean up
!      --------
       deallocate ( bias_3d, rtms_3d, stdv_3d, xcov_3d, nobs_3d, stat=ier )
         if(ier/=0) call die (myname,'Error in Alloc(3d)')

       deallocate ( bias_2d, rtms_2d, stdv_2d, xcov_2d, nobs_2d,stat=ier )
         if(ier/=0) call die (myname,'Error in Dealloc(2d)')

       call GrITAS_Clean ()

!      All done.
!      --------
       if ( trim(otype) .eq. 'hdf' ) then
          call GFIO_Close ( fid, rc )   ! close GFIO file
       end if

       stop
       end

!EOC

!...........................................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GrITAS_LBTally --- Write summary of obs not included in statistics
! 
! !INTERFACE:
!

      subroutine GrITAS_LBTally ( ods_file, date, time,
     &                            ob_flag,  ods,  rc )
!
! !USES:
      use m_ods,   only : ods_vect, 
     &                    ODS_Merge, ODS_Select, ODS_Tally, ODS_Clean
      use m_stdio, only : stdout
!
! !INPUT PARAMETERS:
      character (len=*), intent (in)  ::
     &   ods_file
      integer,           intent (in)  ::
     &   date,
     &   time,
     &   ob_flag (*)
      type (ods_vect),   intent (in)  ::
     &   ods
!
! !OUTPUT PARAMETERS:
      integer,           intent (out) ::
     &   rc

! !DESCRIPTION:
!
! !REVISION HISTORY: 
!     10Jul09  Redder   Initial code.
!EOP
!-------------------------------------------------------------------------

      type ( ods_vect ) :: ods_out
      integer :: iob, nobs, nsel, lu
      integer, dimension (:), pointer :: qchist

      lu = stdout
      nobs   =  ods     % data % nobs

!     Transfer obs to local space
!     ---------------------------
      call ODS_Merge ( (/ ods /), 1, ods_out, rc )
      if ( rc .ne. 0 ) then
         call gritas_perror ( 'GrITAS_LBTally: error status '
     .                     // 'returned from ODS_Merge.' )
         return
      end if

!     Put flag values to local ods structure
!     --------------------------------------
      nobs   =  ods     % data % nobs
      qchist => ods_out % data % qchist ( : nobs )
      ods_out % data % qchist ( : nobs ) = ob_flag ( : nobs )

!     Select only the obs not included in the statistical computations
!     ----------------------------------------------------------------
      call ODS_Select ( ods_out, nobs, nsel, rc,
     &                  qchist = 1 )
      if ( rc .ne. 0 ) then
         call gritas_perror ( 'GrITAS_LBTally: error status '
     .                     // 'returned from ODS_Select.' )
         return
      end if

!     Write header and ...
!     --------------------
      write (lu, '(4(a/)a,i8/a,i6)',iostat = rc)
     &    'ods_tally:===========================================',
     &    'ods_tally:',
     &    'ods_tally: Talley of obs not used in the computation '
     &                               //  ' of statistics',
     &    'ods_tally:   ods file = ' // trim ( ods_file ),
     &    'ods_tally:   date     = ', date,
     &    'ods_tally:   time     = ', time
      if ( rc .ne. 0 ) then
         call gritas_perror ( 'GrITAS_LBTally: error status '
     .                     // 'returned from write statement.' )
         return
      end if

!     ... summary of selected obs
!     ---------------------------
      call ODS_Tally ( lu, ods_out, nsel, rc )
      if ( rc .ne. 0 ) then
         call gritas_perror ( 'GrITAS_LBTally: error status '
     .                     // 'returned from ODS_Tally.' )
         return
      end if

!     Clean up
!     --------
      call ODS_Clean ( ods_out, rc )
 
      return
      end 

!...........................................................................


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GrITAS_Init --- Command line and resource user interface
! 
! !INTERFACE:
!
      subroutine GrITAS_Init ( NDEL_MAX, LMAX, KXMAX,
     &                         ftype, dtype, otype, iopt, ncf, nr,
     &                         split_tgroups, lb_tallyfile,
     &                         nymd_out, DO_red, RC_red,
     &                         ndel, del_ifname, grid_ofname,
     &                         kt_list, kx_list, var_list, var_descr,
     &                         listsz, nrmlz )


! !USES:
 
      Use m_gritas_grids, only : grids_init
      Use m_gritas_grids, only : im, jm, km
      Use m_gritas_grids, only : nlevs
      Use m_gritas_grids, only : plevs
      Use m_inpak90
      Use m_stdio, only : stderr
      Use m_die

!
! !INPUT PARAMETERS:
!
      implicit        NONE
      integer         NDEL_MAX             ! max number of del files    
      integer         LMAX                 ! max list size
      integer         KXMAX                ! max number of kx's

!
! !OUTPUT PARAMETERS:
!

      character(len=*) ftype                 ! File type: 'ods' or 'del'
      character(len=*) dtype                 ! Data type: 'omf', 'oma' or 'obs'
      character(len=*) otype                 ! Output type: 'hdf' or 'eee'
      integer          nr                    ! number of distinct residual series
      integer          iopt                  ! =1 for O-F
                                             ! =3 for OMA 
                                             ! =5 for obs
      integer          nymd_out              ! date for output file
      logical          ncf                   ! Non-compliant format (to allow
                                             !   handling of GSI-output files)
      logical         split_tgroups          ! = true to split the tgroups,
                                             ! bias, stdv and nobs into 
                                             ! separate files.
      character(len=*) lb_tallyfile          ! name of file containing the
                                             !   summary of ODS obs not
                                             !   included in the summary.
                                             !   If blank name then no
                                             !   summary is printed.  Other-
                                             !   wise the filenam is set to
                                             !   stdout for standard output
                                             !   for this implementation.
      integer         ndel                   ! actual number of del files
      character*255   del_ifname(NDEL_MAX)   ! del file names
      character*255   grid_ofname            ! output grid file name

      integer         listsz                 ! actual list size (<LMAX)
      integer         kt_list(LMAX)          ! data type for the grid
      integer         kx_list(KXMAX,LMAX)    ! data sources for the grid
      character*11    var_list(LMAX)         ! variable name for the grid
      character*256   var_descr(LMAX)        ! ... and its description
      logical          DO_red                ! Defines whether to apply reducer or not
      character(len=*) RC_red                ! Recucer resource filename

! !INPUT/OUTPUT PARAMETERS:

       logical         nrmlz                 ! Controls normalization of accumalated

! !DESCRIPTION: This routine parses the command line for the name of
!               the input/output file names and read the resource file
!  file 'gritas.rc' for defining the combination of data type/data
!  sources for each of the output grids.
!
! !REVISION HISTORY: 
!
!  03Sep97  da Silva   Initial code.
!  ???????  G. P. Lou  Modifications.  
!  15Dec98  da Silva   Added ods/gfio, o-a, obs support.
!  25mar99  da Silva   Added -rc option
!  16nov99  da Silva   Added -nopassive
!  09Jun04  Todling    Added -res; updated to use m_gritas_grids
!  14Jun04  Todling    Changed reading of RC file kt/kx table to allow ":"
!  30Jul04  Todling    Added reducer option; bug fix: added i90_release
!  28Oct06  Todling    Added -sigo, and -obias options
!  07Nov06  Todling    Added various opts: joa/job/esigo/hbh/hah
!  23Jan08  Todling    Added -obsimp
!  12Sep08  Todling    Added -nonorm
!  28May09  Redder     Added capability to define kt as a surface variable
!                      by added the prefix, sfckt_mkr, to its name in the
!                      input rc file
!  08Jun09  Redder     Added the command line option, -ospl, and
!                      the subroutine output parameter, split_tgroups
!  10Jul09  Redder     Added the command line option, -lb and
!                      the subroutine output parameter, lb_tallyfile
!  20Jul09  Redder     Added the subroutine output parameter, var_descr.
!  21Jul09  Redder     Fixed bug in generating the list, var_descr
!EOP
!-------------------------------------------------------------------------
!BOC
      character(len=*), parameter :: myname = 'GrITAS_Init'

!     Prefix for variable name denoting a surface variable
!     ----------------------------------------------------
      character (len=*), parameter :: sfckt_mkr = '_'

      integer        iarg, argc, iargc
      character*255  argv, res, rc_ifname, token
      integer        i, j, nkx, lt, ii, jj, ilist, isfckt_mkr
      integer        iret, ios
      integer        kx1, kx2, kxnext
      real p
      logical passive, long_name_found, long_name


      passive = .true.  ! include passive data in calc of stats
      ncf   = .false.
      ftype = 'ods'
      dtype = 'omf'
      otype = 'hdf'
      iopt  = 1
      nr    = 1
      im=72; jm=46; km=26  ! default resolution
      DO_red = .false.
      RC_red = 'reducer.rc'
      nymd_out = -1

!     Parse command line
!     ------------------
      rc_ifname = 'gritas.rc'
      grid_ofname = 'gritas'      
      split_tgroups = .false.
      lb_tallyfile  = ' '
      argc =  iargc()
!      print *, "begin to read command-line input"
      if ( argc .lt. 1 ) call gritas_usage()
      ndel = 0
      iarg = 0
      do i = 1, 32767
         iarg = iarg + 1
!      print *, "iarg = ", iarg
!      print *, "argc = ", argc
         if ( iarg .gt. argc ) go to 111
         call GetArg ( iArg, argv )
!      print *, "argv = ", argv
         if      ( trim(argv) .eq. '-del' )  then
            ftype = 'del'
         else if ( trim(argv) .eq. '-nopassive' )  then
            passive = .false.
         else if ( trim(argv) .eq. '-omf' )  then
            dtype = 'omf'
            iopt  = 1
         else if ( trim(argv) .eq. '-oma' )  then
            dtype = 'oma'
            iopt  = 3
         else if ( trim(argv) .eq. '-obs' )  then
            dtype = 'obs'
            iopt  = 5
         else if ( trim(argv) .eq. '-sigo' )  then
            dtype = 'sigo'
            iopt  = 7
         else if ( trim(argv) .eq. '-obias' )  then
            dtype = 'obias'
            iopt  = 9
         else if ( trim(argv) .eq. '-hbh' )  then
            dtype = 'hbh'
            iopt  = 11
            nr    = 2
         else if ( trim(argv) .eq. '-hah' )  then
            dtype = 'hah'
            iopt  = 13
            nr    = 2
         else if ( trim(argv) .eq. '-esigo' )  then
            dtype = 'esigo'
            iopt  = 15
            nr    = 2
         else if ( trim(argv) .eq. '-job' )  then
            dtype = 'job'
            iopt  = 17
            nr    = 2
         else if ( trim(argv) .eq. '-joa' )  then
            dtype = 'joa'
            iopt  = 19
            nr    = 2
         else if ( trim(argv) .eq. '-oobs' )  then
            dtype = 'oobs'
            iopt  = 21
         else if ( trim(argv) .eq. '-obsimp' )  then
            dtype = 'obsimp'
            iopt  = 23
         else if ( trim(argv) .eq. '-reduce' )  then
            DO_red = .true.
         else if ( trim(argv) .eq. '-nonorm' )  then
            nrmlz = .false.
         else if ( trim(argv) .eq. '-ncf' )  then
            ncf = .true.
         else if ( trim(argv) .eq. '-hdf' )  then
            otype = 'hdf'
         else if ( trim(argv) .eq. '-ieee' )  then
            otype = 'eee'
         else if ( trim(argv) .eq. '-o' ) then
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, grid_ofname )
         else if ( trim(argv) .eq. '-ospl' ) then
            split_tgroups = .true.
         else if ( trim(argv) .eq. '-lb' ) then
            lb_tallyfile = 'stdout'
         else if ( trim(argv) .eq. '-rc' ) then
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, rc_ifname )
         else if ( trim(argv) .eq.'-rcred' ) then
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, RC_red )
         else if ( trim(argv) .eq. '-nlevs' ) then
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, token )
            read(token,*,iostat=ios) km
            if ( ios /= 0 ) then
                 print *, 'cannot parse nlevs ...'
                 call exit(1)
            end if
         else if ( trim(argv) .eq. '-date' ) then
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, token )
            read(token,*,iostat=ios) nymd_out
            if ( ios /= 0 ) then
                 print *, 'cannot parse date ...'
                 call exit(1)
            end if
         else if ( trim(argv) .eq. '-res' ) then
            if ( iarg+1 .gt. argc ) call gritas_usage()
            iarg = iarg + 1
            call GetArg ( iArg, res )
            if ( res(1:1) == "a" ) then
                 im=72; jm=46
             else if ( res(1:1) == "b" ) then
                 im=144; jm=91
             else if ( res(1:1) == "c" ) then
                 im=288; jm=181
             else if ( res(1:1) == "d" ) then
                 im=576; jm=361
             else if ( res(1:1) == "e" ) then
                 im=1152; jm=721
             else
                 print *, 'unknown resolution ... taking default'
             endif
         else
           ndel = ndel + 1
           del_ifname(ndel) = trim(argv)
         end if
      end do
 111  continue
!       print *, "end of command line input"

      if ( ndel .lt. 1 ) call gritas_usage()               

      if ( .not. passive ) iopt = - abs(iopt)   ! exclude passive data

!     Consistency check: del files have only O-F
!     ------------------------------------------
      if ( ftype .eq. 'del' .and. iopt .ne. 1 ) then
         print *
         print *, 'Note: DEL files only have O-F, ignoring -oma/-obs'
         dtype = 'omf'
         iopt = 1
      end if


!     Echo the parameters
!     -------------------
      print *
      print *, '      Resolution: ', im, jm, km
      print *, '   Residual type: ', trim(dtype)
      print *, 'Input  File type: ', trim(ftype)
      print *, 'Output File type: ', trim(otype)
      print *, '       Del files: ', ndel
      do i = 1, ndel
         print *, '     o ', trim(del_ifname(i))
      end do
      
      print *
      print *, 'Grid output file: ', trim(grid_ofname)
      print *

      if ( split_tgroups ) then
         print *
         print *, ' Bias, stdv and nobs will be placed in the files '
         print *, ' with the basenames, '
         print *, '   ' // trim (grid_ofname) // '.bias,'
         print *, '   ' // trim (grid_ofname) // '.stdv and '
         print *, '   ' // trim (grid_ofname) // '.nobs'
         print *
      end if

!     Initialize grid
!     ---------------
      call grids_init ( im, jm, km, iret)
        if(iret/=0) call die(myname,'error from grids_init')

!     Read resource file
!     ------------------
      call i90_LoadF ( trim(rc_ifname), iret )
      if ( iret .ne. 0 ) then
         print *, 'GrITAS: cannot open resource file ',
     &            trim(rc_ifname)
         call exit(1)
      end if

!     vertical levels
!     ---------------
      call i90_label ( 'GrITAS*Vertical_Levels:', iret )
      if ( iret .ne. 0 ) call gritas_perror
     &   ('GrITAS: cannot find vertical levels on rc file' )
      nlevs = 0
      do i = 1, 32767
         p = i90_gfloat(iret)
         if ( iret .ne. 0 ) go to 222
         nlevs = nlevs + 1
         if ( nlevs .gt. km ) 
     &      call gritas_perror('GrITAS: too many vertical levels')
         plevs(nlevs) = p
      end do
 222  continue
      if ( nlevs .eq. 0 )
     &      call gritas_perror('GrITAS: no vertical levels')

      if ( plevs(1) < plevs(2) ) plevs(1:nlevs) = plevs(nlevs:1:-1)

      print *
      print *, 'Vertical Levels: ', ( plevs(i), i=1,nlevs )

!     Initialize list
!     ---------------
      listsz = 0
      do i = 1, LMAX
         do j = 1, KXMAX
            kx_list(j,i) = -1
         end do
      end do

!     kx-kt list
!     ----------
      long_name_found = .false.
      kx_list = -1
      call i90_label ( 'GrITAS*List::', iret )
      if ( iret .ne. 0 ) call
     &   gritas_perror('GrITAS: cannot find kx-kt table on rc file')
      do 10 i = 1, 32767
         call i90_gline ( iret )
         if ( iret .eq. -1 ) 
     &        call gritas_perror ('GrITAS: premature end of table')
         if ( iret .eq. 1 ) go to 11
         listsz = listsz + 1
         if(listsz .gt. LMAX )
     &      call gritas_perror('GrITAS: kx-kt list is too large')
         kt_list(listsz) = i90_gint(iret)
         if (iret.ne.0) call gritas_perror('GrITAS: table error: kt')
         call i90_gtoken ( argv, iret )
         if (iret.ne.0) call gritas_perror('GrITAS: table error: var')
         isfckt_mkr = scan ( argv, sfckt_mkr )
         if ( isfckt_mkr .eq. 1 ) then
            kt_list(listsz) = -1 * kt_list(listsz)
            var_list(listsz) = argv(isfckt_mkr+1:)
         else
            var_list(listsz) = argv
         end if
         var_descr(listsz) = ' '

         j = 0
         do jj = 1, KXMAX
            call I90_GToken(token, iret )
!            print *, trim(token)
            if (iret/=0) then      ! read error
                exit
!               call die(myname,'I90_GToken error, iret=',iret)
            end if
            long_name = verify ( token, '-0123456780: ') .gt. 0 .and.
     &                  jj .eq. 1
            if ( long_name ) then
               var_descr ( listsz ) =  token
               long_name_found      = .true.
            else
               ii = index(token,':') ! token is single kx or range of kx's
               lt = len_trim(token)
               if (ii==0) then    ! no colon, therefore single kx
                  read(token,*) kx1
                  kx2 = kx1
               else                  ! colon, therefore kx1:kx2
                  read(token(1:ii-1),*) kx1
                  read(token(ii+1:lt),*) kx2
               end if
               if (kx1>kx2) then      ! range error
                  write(stderr,'(2a,i5)')
     &               myname, ': Invalid range: ', token
                  call die(myname)
               end if
               do kxnext = kx1, kx2
                  if (j==KXMAX) then    ! check space
                     write(stderr,'(2a,i5)') myname,': increase KXMAX'
                     call die(myname)
                  else if (iret==0) then
                     j = j + 1
                     kx_list(j,listsz) = kxnext
                  end if
               end do
            end if
         end do ! while

 21      continue
 10   continue
 11   continue

      if ( .not. long_name_found ) then
         do ilist = 1, listsz
            var_descr ( ilist ) = var_list ( ilist )
         end do
      end if

!     Echo the parameters
!     -------------------
      print *
      print *, 'Kx-Kt table...'
      do i = 1, listsz

         nkx = 0
         do j = 1, KXMAX
            if ( kx_list(j,i) .lt. 0 ) go to 51
            nkx = j
         end do
 51      continue


         print *
         print *, 'Variable ', i, ': ', trim(var_list(i))
         print *, '   o kt = ', abs(kt_list(i))
         print *, '   o kx = ', (kx_list(j,i),j=1,nkx)

      end do


!     Release resources
!     -----------------
      call I90_Release()

!     All done
!     --------
      return
      end
        

!...........................................................................
      subroutine GrITAS_Clean ()
      use m_gritas_grids, only : grids_clean
      use m_die
      implicit none
      integer ier
      call grids_clean(ier)
        if(ier/=0) call die('GrITAS_Clean','error from grids_clean')
      end subroutine GrITAS_Clean
!...........................................................................

      subroutine gritas_usage()
      print *
      print *,'Usage:  gritas.x  [-rc fn] [-del]  [-oma|-obs]  [-ieee] [-reduce]'
      print *,'                  [-rcred fn] [-o fn] [-ncf] [-res RES] obs_file(s)'
      print *
      print *,  'where'
      print *
      print *, '-rc fname   resource file name (default: gritas.rc)'
      print *, '-del        input files are GEOS-1 DEL files'
      print *, '             (by default input files are ODS)'
      print *, '-nopassive  passive data types will not be included'
      print *, '-oma        produces gridded O-A residuals instead'
      print *, '              of O-F residuals'
      print *, '-obs        produces gridded observed values instead'
      print *, '              of O-F residuals'
      print *, '-obias      produces gridded observed bias values instead'
      print *, '              of O-F residuals'
      print *, '-obsimp     produces observation impacts'
      print *, '-sigo       produces observation error values'
      print *, '-esigo      produces estimate of sigo from residuals'
      print *, '-hbh        produces gridded O-A x O-B residual statistics'
      print *, '-hah        produces gridded A-F x O-A residual statistics'
      print *, '-job        produces gridded Jo(O-F)/p'
      print *, '-joa        produces gridded Jo(O-A)/p'
      print *, '-oobs       produces gridded for original obs (not bias corrected)'
      print *, '-ieee       output file is in IEEE binary format'
      print *, '              instead of HDF/GFIO format; both formats'
      print *, '              are GrADS compatible'
      print *, '-reduce     to apply reducer (default RC file: reducer.rc)'
      print *, '-rcred fn   specify alternative name for reducer resource file'
      print *, '-ncf        non-compliant format (apply for GSI-diag-conv files)'
      print *, '-nonorm     skips normalization of accumulated summs (good for Obs Impact)'
      print *, '-res RES    resolution definition (X=a, b, or c) '
      print *, '            (default: a)'
      print *, '-nlevs km   Number of vertical levels'
      print *, '            (default: 26)'
      print *, '-date nymd  Date for output file; default is mean of input'
      print *, '-ospl       Split bias, stdv and nobs into separate files'
      print *, '            (only for sdf files)'
      print *, '-lb         Write tally of obs not included in'
      print *, '            computation of statistics'
      print *, '-o bname    specifies output file base name, optional '
      print *, '            (default: gritas)'
      print *
      print *, 'obs_files   ODS, DEL, or DIAG_CONV file names (no default)'
      print *
      print *
      print *, ' NOTES: '
      print *
      write(6,*) ' 1) In case of calculating statistics of (specified) sigo and obias the standard deviation '
      write(6,*) '    output is replaced by root-mean-square; that is to say that position t=1 in the grads '
      write(6,*) '    file will have the mean as it normally has, and t=2 will have the RMS instead of the STDV. '
      print *
      write(6,*)  ' 2) In case of Desroziers et al. diagnostics, the mean and standard deviations output '
      write(6,*) '    are replaced by cross-variance and the cross standard deviations, respectively. '
      print *
      print *, 'Revised: 12 September 2008 (R. Todling)'
      print *, 'Created:  1 September 1997 (A. da Silva)'
      print *
      call gritas_perror('GrITAS: not enough input parameters')
      end


!EOC




