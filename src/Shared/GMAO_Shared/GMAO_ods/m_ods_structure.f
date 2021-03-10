!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  ODS_Structure --- Defines the ODS Structure
!
! !INTERFACE:
!

      module  m_ods_structure

! !USES:

      use  m_odsmeta

      implicit none

! !PUBLIC TYPES:

      PRIVATE
      PUBLIC  ods_attributes_type
      PUBLIC  ods_description_type
      PUBLIC  ods_structure

! !PUBLIC MEMBER FUNCTIONS:

      PUBLIC  Create_ODS_Structure
      PUBLIC  Destroy_ODS_Structure
      PUBLIC  Create_ODS_Attributes
      PUBLIC  Destroy_ODS_Attributes
      PUBLIC  Create_ODS_Description
      PUBLIC  Destroy_ODS_Description
      PUBLIC  ODS_GEOS4to2
      PUBLIC  ODS_GEOS2to4
      PUBLIC  Tally
      PUBLIC  ODS_Reorder
      PUBLIC  ODS_MoveUp

! !DESCRIPTION: Defines all ODS attributes as a f90 structure.
!
! !REVISION HISTORY:
!
!  16Nov98   Todling   Initial code.
!  15Jan99   Todling   Added description data type.
!  28Dec99   Todling   Added odsmeta and removed  kt_max.
!  13Apr00   Todling   Made Desc public.
!  07Feb01   Todling   Added Xvec to structure.
!  09May02   Dee       Moved ODS_Reorder, ODS_MoveUp to this module
!
!EOP
!-------------------------------------------------------------------------

!     Surface state variables
!     ------------------------
      type    ods_attributes_type

                integer             nobs      ! number of observations

                real,    pointer :: lat(:)    ! latitute     of obs (degrees)
                real,    pointer :: lon(:)    ! longitude    of obs (degrees)
                real,    pointer :: lev(:)    ! level        of obs (hPa)
                integer, pointer :: kx(:)     ! data source index
                integer, pointer :: kt(:)     ! data type   index
                integer, pointer :: ks(:)     ! sounding    index
                real,    pointer :: xm(:)     ! metadata
                integer, pointer :: time(:)   ! time

                real,    pointer :: obs(:)    ! observation value
                                              !   (units depend on kt)
                real,    pointer :: OmF(:)    ! obs minus forecast (O-F) innovations
                real,    pointer :: OmA(:)    ! obs minus analysis (O-A) residuals
                real,    pointer :: Xvec(:)   ! PSAS CG solution vector
                integer, pointer :: qcexcl(:) ! On-line QC exclusion flag
                integer, pointer :: qchist(:) ! On-line QC history flag

                integer, pointer :: kid(:)    ! Obs identification index


      end type  ods_attributes_type

      type    ods_description_type

                        !  kt-related entries
                integer  ktmax
                logical          , pointer ::  ktmvar(:,:)
                character(len=8) , pointer ::  ktname(:)
                character(len=8) , pointer ::  ktunit(:)
                character(len=32), pointer ::  ktdesc(:)

                        !  kx-related entries
                integer  kxmax

                integer  , pointer :: kxrank(:)     ! quality ranks, [kxmax]
                integer  , pointer :: kxtmask(:,:)  ! mask data types of data sources
                integer  , pointer :: i_kxclas(:)   ! indices to module OEclass_tbl
                integer  , pointer :: i_hoecHc(:)   ! indices to module OEhcor_tbl
                integer  , pointer :: i_voecHc(:)   ! indices to module OEvcor_tbl
                integer  , pointer :: i_voecHu(:)   ! indices to module OEvcor_tbl

                character(len=8) , pointer ::  kxclas(:)
                character(len=32), pointer ::  kxdesc(:)

                        !  metadata-related entries
                integer  nkx_meta          ! number of metadata pieces
                character(len=32), pointer ::  kxmeta(:)

                        !  qcflag-related entries
!_SOON          integer  nqcXmax
!_SOON          character(len=32), pointer :: qcXnames(:)

      end type  ods_description_type

      type  ods_structure
                type ( ods_attributes_type  ) attr
                type ( ods_description_type ) desc
      end type  ods_structure


      CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Create_ODS_Structure --- Allocate space for all ODS structure
!
! !INTERFACE:
!
      subroutine Create_ODS_Structure ( ods, nobs_max )

! !USES:

      implicit none
      include 'ods_stdio.h'

! !INPUT PARAMETERS:

      integer nobs_max

! !INPUT/OUTPUT PARAMETERS:

      type ( ods_structure  ) ods

! !DESCRIPTION: Allocates space for all ODS structure.
!
! !REVISION HISTORY:
!
!  16Nov98   Todling   Initial code
!  16Feb00   Todling   rename stdio.h to ods_stdio.h
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname = 'Create_ODS_Structure'

      call Create_ODS_Attributes  ( ods%attr, nobs_max )
      call Create_ODS_Description ( ods%desc )

      return
      end subroutine Create_ODS_Structure

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Create_ODS_Attributes --- Allocate space for ODS attributes
!
! !INTERFACE:
!
      subroutine Create_ODS_Attributes ( attr, nobs_max )

! !USES:

      implicit none
      include 'ods_stdio.h'

! !INPUT PARAMETERS:

      integer nobs_max

! !INPUT/OUTPUT PARAMETERS:

      type ( ods_attributes_type  ) attr

! !DESCRIPTION: Allocates space for ODS attributes.
!
! !REVISION HISTORY:
!
!  16Nov98   Todling   Initial code
!  16Feb00   Todling   rename stdio.h to ods_stdio.h
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname = 'Create_ODS_Attributes'
      integer        ier

!     Initialize number of observations
!     ---------------------------------
      attr%nobs = 0

!     Allocate space for ODS-structure
!     --------------------------------
      allocate ( attr%lat(nobs_max),    attr%lon(nobs_max), attr%lev(nobs_max),
     .            attr%kx(nobs_max),     attr%kt(nobs_max),  attr%ks(nobs_max),
     .            attr%xm(nobs_max),   attr%time(nobs_max), attr%obs(nobs_max),
     .           attr%OmF(nobs_max),    attr%OmA(nobs_max), attr%Xvec(nobs_max),
     .        attr%qcexcl(nobs_max), attr%qchist(nobs_max),
     .           attr%kid(nobs_max),
     .          stat = ier )
      if ( ier .ne. 0 ) then
          write(stderr,'(2a,i5)') myname,': allocate(ods_struct) error,
     .          stat =', ier
          call exit(2)
      endif

!     All done
!     --------
      write(stdout,'(a)') 'Done allocating space for ODS attributes'
      call flush ( stdout )
      return
      end subroutine Create_ODS_Attributes


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Create_ODS_Description --- Allocates and defines ODS desc vars
!
! !INTERFACE:
!
      subroutine Create_ODS_Description ( desc )

! !USES:

      implicit none

      include 'ods_stdio.h'

! !INPUT/OUTPUT PARAMETERS:

      type ( ods_description_type  ) desc

! !DESCRIPTION: Allocates space for ODS attributes.
!
! !REVISION HISTORY:
!
!  16Nov98   Todling   Initial code
!  16Feb00   Todling   rename stdio.h to ods_stdio.h
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname = 'Create_ODS_Description'
      integer        ier, iret, jret
      integer        j
      integer        kx_next
      integer        nkx_active
      integer        kx_active(kxmax)
      character*255  word_in

!     Set description quantities
!     --------------------------
      desc%ktmax    = ktmax
      desc%kxmax    = kxmax


!     Allocate space for ODS-description quantities
!     ---------------------------------------------
      allocate ( desc%ktmvar(ktmax,ktmax),
     .           desc%ktname(ktmax), desc%ktunit(ktmax), desc%ktdesc(ktmax),
     .          stat = ier )
      if ( ier .ne. 0 ) then
          write(stderr,'(2a,i5)') myname,': allocate(desc_kt) error,
     .          stat =', ier
          call exit(2)
      endif

      allocate ( desc%kxclas(kxmax),   desc%kxdesc(kxmax),
     .           desc%kxrank(kxmax),   desc%kxtmask(kxmax,ktmax),
     .           desc%i_kxclas(kxmax),
     .           desc%i_hoecHc(kxmax),
     .           desc%i_voecHc(kxmax), desc%i_voecHu(kxmax),
     .           desc%kxmeta(kxmax),
     .          stat = ier )
      if ( ier .ne. 0 ) then
          write(stderr,'(2a,i5)') myname,': allocate(desc_kx) error,
     .          stat =', ier
          call exit(2)
      endif

      return
      end subroutine Create_ODS_Description


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_GEOS4to2 - Convert GEOS-4 obs vector to GEOS-2 equivalent
!
! !INTERFACE:
!
      subroutine ODS_GEOS4to2 ( y, ods, null )

! !USES:

      use m_ods

      implicit none

! !INPUT PARAMETERS:

      type(obs_vect), intent(in)             ::  y  ! New (GEOS-4) obs vector

      logical, intent(in), optional          :: null ! Set to .t. to nullify

! !OUTPUT PARAMETERS:

      type(ods_attributes_type), intent(out) :: ods ! Old (GEOS-2) obs vector
!
!
! !DESCRIPTION: Converts a new (GEOS-4) observation vector to the old
!               (GEOS-4) observation vector type ({\tt ods\_attributes\_type
!  in GEOS-2 parlance).
!
! !REVISION HISTORY:
!
!  12oct1999  da Silva   Initial code.
!  07Feb2001  Todling    Added Xvec.
!  17Feb2001  Todling    Added optional parameter to nullify pointers
!
!
!EOP
!-------------------------------------------------------------------------

      if (present(null)) then
          if ( null ) then
               nullify(ods%lat)
               nullify(ods%lon)
               nullify(ods%lev)
               nullify(ods%kx)
               nullify(ods%kt)
               nullify(ods%ks)
               nullify(ods%xm)
               nullify(ods%time)
               nullify(ods%obs)
               nullify(ods%omf)
               nullify(ods%oma)
               nullify(ods%Xvec)
               nullify(ods%qcexcl)
               nullify(ods%qchist)
               nullify(ods%kid)
          return
          end if
      end if

!     Number of observations
!     ----------------------
      ods%nobs = y%nobs
      if ( y%nobs .eq. 0 ) return

!     Reassign pointers
!     -----------------
      ods%lat    => y%lat
      ods%lon    => y%lon
      ods%lev    => y%lev
      ods%kx     => y%kx
      ods%kt     => y%kt
      ods%ks     => y%ks
      ods%xm     => y%xm
      ods%time   => y%time
      ods%obs    => y%obs
      ods%omf    => y%omf
      ods%oma    => y%oma
      ods%Xvec   => y%Xvec
      ods%qcexcl => y%qcexcl
      ods%qchist => y%qchist
      ods%kid    => y%kid


      end subroutine ODS_GEOS4to2


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_GEOS2to4 - Convert GEOS-2 obs vector to GEOS-4 equivalent
!
! !INTERFACE:
!
      subroutine ODS_GEOS2to4 ( ods, y )

! !USES:

      use m_ods

      implicit none

! !INPUT PARAMETERS:

      type(ods_attributes_type), intent(in) :: ods ! Old (GEOS-2) obs vector

! !OUTPUT PARAMETERS:

      type(obs_vect), intent(out)           ::  y  ! New (GEOS-4) obs vector
!
!
! !DESCRIPTION: Converts a old (GEOS-2/3) observation vector to the new
!               (GEOS-4) observation vector type.
!
! !REVISION HISTORY:
!
!  15Feb2001  Todling    Created based on GEOS4to2
!  09Mar2001  Todling    Bug fix in xvec pointer.
!
!
!EOP
!-------------------------------------------------------------------------

!     Number of observations
!     ----------------------
      y%nobs = ods%nobs
      if ( ods%nobs .eq. 0 ) return

!     Reassign pointers
!     -----------------
      y%lat    => ods%lat
      y%lon    => ods%lon
      y%lev    => ods%lev
      y%kx     => ods%kx
      y%kt     => ods%kt
      y%ks     => ods%ks
      y%xm     => ods%xm
      y%time   => ods%time
      y%obs    => ods%obs
      y%omf    => ods%omf
      y%oma    => ods%oma
      y%Xvec   => ods%Xvec
      y%qcexcl => ods%qcexcl
      y%qchist => ods%qchist
      y%kid    => ods%kid


      end subroutine ODS_GEOS2to4

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tally --- Summarizes information on observations
!
! !INTERFACE:
!

      subroutine Tally ( ods, nobs, nymd, nhms, label )

! !USES:

      Implicit NONE


      include 'ods_stdio.h'

! !INPUT PARAMETERS:

      type ( ods_structure ) ods

      integer   nobs                 ! total no. of observations

                                     ! Synoptic time:
      integer   nymd                 !  Year-month-day, e.g.,  19971012
      integer   nhms                 !  hour-minute-sec, e.g.,   120000

      character*(*) label            ! label for tally table

!
! !DESCRIPTION: Summarizes information on observations.
!               Mimics OI tally routine.
!
! !REVISION HISTORY:
!
!  15Aug98   Todling    Initial code.
!  14Dec98   Todling    Passing ods structure.
!  16Feb00   Todling   rename stdio.h to ods_stdio.h
!
!EOP
!-------------------------------------------------------------------------


      character(len=*), parameter ::  myname = 'Tally'

      integer  n_obs(kxmax)    ! no of obs per kx
      integer  indx(kxmax)     ! index for actual kx's
      integer  t_obs_good      ! total number of good obs

                               ! no of obs per kx
      integer  n_excl          !   obs excluded
      integer  n_slp           !   surf pressure obs
      integer  n_swnd          !   surf wind obs
      integer  n_uwnd          !   upper-air wind obs
      integer  n_hght          !   upper-air height obs
      integer  n_mixr          !   mixing ratio obs

                               ! no of obs for all kx's
      integer  t_excl          !   obs excluded
      integer  t_slp           !   surf pressure obs
      integer  t_swnd          !   surf wind obs
      integer  t_uwnd          !   upper-air wind obs
      integer  t_hght          !   upper-air height obs
      integer  t_mixr          !   mixing ratio obs

      integer  i, k, n

!     Make these into pointers at first during implementation of ods.struc
!     --------------------------------------------------------------------
      real   , pointer ::   lat(:)      ! latitute     of obs (degrees)
      real   , pointer ::   lon(:)      ! longitude    of obs (degrees)
      real   , pointer ::   lev(:)      ! level        of obs (hPa)
      integer, pointer ::   kx(:)       ! data source index
      integer, pointer ::   kt(:)       ! data type   index
      integer, pointer ::  qcexcl(:)    ! exclusion flags

      lat    => ods%attr%lat(1:nobs)
      lon    => ods%attr%lon(1:nobs)
      lev    => ods%attr%lev(1:nobs)
      kx     => ods%attr%kx(1:nobs)
      kt     => ods%attr%kt(1:nobs)
      qcexcl => ods%attr%qcexcl(1:nobs)


!     Nothing to do
!     -------------
      if ( nobs .eq. 0 ) return

!     Count number of excluded observations per kx
!     --------------------------------------------
      n = 0
      do k = 1, kxmax
         n_obs(k)  = count(kx .eq. k, dim=1)
         if ( n_obs(k) .ne. 0 ) then
            n = n + 1
            indx(n) = k
         end if
      end do

      write(stdout,'(80("*"))')
      write(stdout,'(20x,a,3x,i8.8,2x,i6.6)')
     &             'Observation Count Summary:', nymd, nhms
      write(stdout,'(20x,3a)')  myname, ': ', label
      write(stdout,'(1x,a,1x,a,16x,a,6(1x,a))')
     &             'kx'    , 'Data Type',
     &             'SurfP ', 'SurfW ',
     &             'UpAirW', 'UpAirH', 'UpMixR',
     &             ' Excl ', 'Total '
      t_excl = 0
      t_slp  = 0
      t_swnd = 0
      t_uwnd = 0
      t_hght = 0
      t_mixr = 0
      do i = 1, n
         k = indx(i)
         n_excl = count( kx .eq. k .and. qcexcl .ne. 0, dim=1)
         n_swnd = count((kt .eq. ktUS .or.  kt .eq. ktVS) .and.
     &                   kx .eq. k .and. qcexcl .eq. 0, dim=1)
         n_slp  = count( kt .eq. ktSLP .and.
     &                   kx .eq. k .and. qcexcl .eq. 0, dim=1)
         n_uwnd = count((kt .eq. ktUU .or.  kt .eq. ktVV) .and.
     &                   kx .eq. k .and. qcexcl .eq. 0, dim=1)
         n_hght = count( kt .eq. ktHH .and.
     &                   kx .eq. k .and. qcexcl .eq. 0, dim=1)
         n_mixr = count( kt .eq. ktWW .and.
     &                   kx .eq. k .and. qcexcl .eq. 0, dim=1)
         write(stdout,'(i3,1x,a23,7(1x,i6))')
     &         k, ods%desc%kxdesc(k), n_slp , n_swnd,
     &                                n_uwnd, n_hght, n_mixr,
     &                                n_excl, n_obs(k)

         t_excl = t_excl + n_excl
         t_slp  = t_slp  + n_slp
         t_swnd = t_swnd + n_swnd
         t_uwnd = t_uwnd + n_uwnd
         t_hght = t_hght + n_hght
         t_mixr = t_mixr + n_mixr
      end do
      t_obs_good = count( qcexcl .eq. 0, dim=1)
      write(stdout,'(80("-"))')
      write(stdout,'(2x,a25,7(1x,i6))')
     &             'Totals', t_slp , t_swnd,
     &                       t_uwnd, t_hght, t_mixr,
     &                       t_excl, t_obs_good
      write(stdout,'(80("-"))')
      write(stdout,'(2x,a,/,2x,6("-"))') 'Notes: '
      write(stdout,'(2(8x,a,/),80("*"))')
     &      '1. Total & Exclusion counts for all kts per data-type',
     &      '2. All other counts for good obs only'

      do k = 1, ktmax
         i = count(kt.eq.k,dim=1)
         if(i.ne.0) print*, 'kt = ', k, ' found = ', i
      end do
      return
      end subroutine tally

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Destroy_ODS_Structure - Deallocated all ODS structure
!
! !INTERFACE:
!
      subroutine destroy_ods_structure ( ods )

! !USES:

      implicit none

! !INPUT/OUTPUT PARAMETERS:

      type ( ods_structure ) ods

! !DESCRIPTION: Wipe out ODS structure memory allocation.
!
! !REVISION HISTORY:
!
!  15Jan99   Todling   Initial code
!
!EOP
!-------------------------------------------------------------------------

      call destroy_ods_description (ods%desc)
      call destroy_ods_attributes  (ods%attr)

      return
      end subroutine  destroy_ods_structure

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Destroy_ODS_Attributes - Deallocated ODS attribute arrays
!
! !INTERFACE:
!
      subroutine Destroy_ODS_Attributes ( vars )

! !USES:

      implicit none

! !INPUT/OUTPUT PARAMETERS:

      type ( ods_attributes_type ) vars

! !DESCRIPTION: Wipe out ODS attributes memory allocation.
!
! !REVISION HISTORY:
!
!  16Nov98   Todling   Initial code
!
!EOP
!-------------------------------------------------------------------------

      if(associated( vars%lat     )) deallocate ( vars%lat   )
      if(associated( vars%lon     )) deallocate ( vars%lon   )
      if(associated( vars%lev     )) deallocate ( vars%lev   )
      if(associated( vars%kx      )) deallocate ( vars%kx    )
      if(associated( vars%kt      )) deallocate ( vars%kt    )
      if(associated( vars%ks      )) deallocate ( vars%ks    )
      if(associated( vars%xm      )) deallocate ( vars%xm    )
      if(associated( vars%time    )) deallocate ( vars%time  )
      if(associated( vars%obs     )) deallocate ( vars%obs   )
      if(associated( vars%OmF     )) deallocate ( vars%OmF   )
      if(associated( vars%OmA     )) deallocate ( vars%OmA   )
      if(associated( vars%Xvec    )) deallocate ( vars%Xvec  )
      if(associated( vars%qcexcl  )) deallocate ( vars%qcexcl)
      if(associated( vars%qchist  )) deallocate ( vars%qchist)

      if(associated( vars%kid     )) deallocate ( vars%kid   )

      return
      end subroutine  Destroy_ODS_Attributes


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Destroy_ODS_Description - Deallocated ODS description arrays
!
! !INTERFACE:
!
      subroutine Destroy_ODS_Description ( vars )

! !USES:

      implicit none

! !INPUT/OUTPUT PARAMETERS:

      type ( ods_description_type ) vars

! !DESCRIPTION: Wipe out ODS description memory allocation.
!
! !REVISION HISTORY:
!
!  15Jan99   Todling   Initial code
!
!EOP
!-------------------------------------------------------------------------


      if(associated( vars%ktdesc  )) deallocate ( vars%ktdesc )
      if(associated( vars%ktunit  )) deallocate ( vars%ktunit )
      if(associated( vars%ktname  )) deallocate ( vars%ktname )
      if(associated( vars%ktmvar  )) deallocate ( vars%ktmvar )

      if(associated( vars%kxmeta   )) deallocate ( vars%kxmeta   )
      if(associated( vars%i_voecHu )) deallocate ( vars%i_voecHu )
      if(associated( vars%i_voecHc )) deallocate ( vars%i_voecHc )
      if(associated( vars%i_hoecHc )) deallocate ( vars%i_hoecHc )
      if(associated( vars%i_kxclas )) deallocate ( vars%i_kxclas )
      if(associated( vars%kxtmask  )) deallocate ( vars%kxtmask  )
      if(associated( vars%kxrank   )) deallocate ( vars%kxrank   )
      if(associated( vars%kxdesc   )) deallocate ( vars%kxdesc   )
      if(associated( vars%kxclas   )) deallocate ( vars%kxclas   )


      return
      end subroutine  Destroy_ODS_Description

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_Reorder() --- Reorder data in ODS structure.
!
! !INTERFACE:
!
      subroutine ODS_Reorder ( nobs, ods, indx )

!
! !INPUT PARAMETERS:
!
      integer , intent(in)    ::  nobs          ! number of obs
      integer , intent(in)    ::  indx(nobs)    ! sorting index
!
! !INPUT/OUTPUT PARAMETERS:
!
      type ( ods_attributes_type ) , intent (inout) :: ods
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  04Feb99 (Dee)      Initial code
!  07Feb01 (Todling)  Added Xvec.
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname='ODS_Reorder'

      integer i

!     Pointers to ods attributes:
!     --------------------------
      integer, pointer ::   kid(:)      ! observation index
      integer, pointer ::   kt(:)       ! data type index
      integer, pointer ::   kx(:)       ! data source index
      integer, pointer ::   ks(:)       ! sounding index
      real   , pointer ::   lon(:)      ! longitude of obs (degrees)
      real   , pointer ::   lat(:)      ! latitute of obs (degrees)
      real   , pointer ::   lev(:)      ! pressure level of obs (hPa)
      integer, pointer ::   time(:)     ! time
      real   , pointer ::   obs(:)      ! observation
      real   , pointer ::   OmF(:)      ! obs-minus-fcst (O-F)
      real   , pointer ::   OmA(:)      ! obs-minus-anal (O-A)
      real   , pointer ::   Xvec(:)     ! PSAS CG solution vector
      real   , pointer ::   xm(:)       ! metadata
      integer, pointer ::   qcx(:)      ! exclusion mark
      integer, pointer ::   qch(:)      ! history mark

      kid    => ods%kid
      kt     => ods%kt
      kx     => ods%kx
      ks     => ods%ks
      lon    => ods%lon
      lat    => ods%lat
      lev    => ods%lev
      time   => ods%time
      obs    => ods%obs
      OmF    => ods%OmF
      OmA    => ods%OmA
      Xvec   => ods%Xvec
      xm     => ods%xm
      qcx    => ods%qcexcl
      qch    => ods%qchist

!     Nothing to do
!     -------------
      if ( nobs==0 ) return

!     Apply sorting permutation to data arrays
!     ----------------------------------------
      kid (1:nobs) = kid ( (/ (indx(i), i=1,nobs) /) )
      kt  (1:nobs) = kt  ( (/ (indx(i), i=1,nobs) /) )
      kx  (1:nobs) = kx  ( (/ (indx(i), i=1,nobs) /) )
      ks  (1:nobs) = ks  ( (/ (indx(i), i=1,nobs) /) )
      lon (1:nobs) = lon ( (/ (indx(i), i=1,nobs) /) )
      lat (1:nobs) = lat ( (/ (indx(i), i=1,nobs) /) )
      lev (1:nobs) = lev ( (/ (indx(i), i=1,nobs) /) )
      time(1:nobs) = time( (/ (indx(i), i=1,nobs) /) )
      obs (1:nobs) = obs ( (/ (indx(i), i=1,nobs) /) )
      OmF (1:nobs) = OmF ( (/ (indx(i), i=1,nobs) /) )
      OmA (1:nobs) = OmA ( (/ (indx(i), i=1,nobs) /) )
      Xvec(1:nobs) = Xvec( (/ (indx(i), i=1,nobs) /) )
      xm  (1:nobs) = xm  ( (/ (indx(i), i=1,nobs) /) )
      qcx (1:nobs) = qcx ( (/ (indx(i), i=1,nobs) /) )
      qch (1:nobs) = qch ( (/ (indx(i), i=1,nobs) /) )

      end subroutine ODS_Reorder

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_MoveUp() --- Move selected observations to front.
!
! !INTERFACE:
!
      subroutine ODS_MoveUp ( ods, nobs, nmoved,
     &                             attrname, attrvalues, nvalues )

! !USES:
!

      implicit NONE
      include 'ods_stdio.h'
!
! !INPUT PARAMETERS:
!
      integer , intent(in)      ::  nobs       ! total number of obs
      character*(*), intent(in) ::  attrname   ! name of attribute
      integer , intent(in)      ::  nvalues    ! number of selection values
      integer , intent(in)      ::  attrvalues(nvalues)  ! selection values
!
! !INPUT/OUTPUT PARAMETERS:
!
      type ( ods_attributes_type ) , intent (in out) :: ods
!
! !OUTPUT PARAMETERS:
!
      integer , intent(out)   ::  nmoved         ! number of data moved up
!
! !DESCRIPTION: Moves observations with selected attribute values to the
!     front of the list. Only works for integer attributes.
!
! !TO DO:
!     make this more efficient.
!
! !REVISION HISTORY:
!
!  04Feb99 (Dee)  Initial code
!  08Jan01 (AdS)  Removed trim from inside select to avoid Linux/PGI bug.
!  09May02 (Dee)  Moved here from psas_qc.f; removed m_MergeSort dependency
!
!EOP
!-------------------------------------------------------------------------

!     Local work space
!     ----------------
      integer, allocatable ::  indx(:)
      logical, allocatable ::  selected(:)

      character(len=*), parameter :: myname='ODS_MoveUp'
      integer ios, i

!     Pointers to ods attributes:
!     --------------------------
      integer, pointer ::   kid(:)      ! observation index
      integer, pointer ::   kt(:)       ! data type index
      integer, pointer ::   kx(:)       ! data source index
      integer, pointer ::   ks(:)       ! sounding index
      integer, pointer ::   time(:)     ! time
      integer, pointer ::   qcx(:)      ! exclusion mark
      integer, pointer ::   qch(:)      ! history mark

      kid    => ods%kid
      kt     => ods%kt
      kx     => ods%kx
      ks     => ods%ks
      time   => ods%time
      qcx    => ods%qcexcl
      qch    => ods%qchist

!     Nothing to do
!     -------------
      if ( nobs == 0 ) then
          nmoved = 0
          return
      end if

!     Allocate memory
!     ---------------
      allocate ( indx(nobs), selected(nobs), stat=ios )
      if ( ios/=0 ) then
          write(stderr,'(2a)') myname,': allocate() error'
          call exit(2)
      end if


      selected = .false.
      select case (attrname)
      case ('kid')
           do i = 1, nvalues
              where (kid ==attrvalues(i)) selected = .true.
           end do
      case ('kt')
           do i = 1, nvalues
              where (kt  ==attrvalues(i)) selected = .true.
           end do
      case ('kx')
           do i = 1, nvalues
              where (kx  ==attrvalues(i)) selected = .true.
           end do
      case ('ks')
           do i = 1, nvalues
              where (ks  ==attrvalues(i)) selected = .true.
           end do
      case ('time')
           do i = 1, nvalues
              where (time==attrvalues(i)) selected = .true.
           end do
      case ('qcexcl')
           do i = 1, nvalues
              where (qcx ==attrvalues(i)) selected = .true.
           end do
      case ('qchist')
           do i = 1, nvalues
              where (qch ==attrvalues(i)) selected = .true.
           end do
      case default
          write(stderr,'(5a)') myname, ': ',
     &          'invalid attribute name (', attrname, ')'
          call exit(2)
      end select

      indx = (/ (i,i=1,nobs) /)
      indx = (/ pack(indx,mask=selected), pack(indx,mask=.not. selected) /)
      call ODS_Reorder ( nobs, ods, indx )

!     Figure out how many were moved up
!     ---------------------------------
      nmoved = count(selected)

!     Deallocate memory
!     ---------------
      deallocate ( selected, indx )

      end subroutine ODS_MoveUp


      end module  m_ods_structure
