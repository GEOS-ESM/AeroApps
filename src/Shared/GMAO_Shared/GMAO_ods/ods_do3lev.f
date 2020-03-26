!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 610.1, GEOS/DAS      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!BOP
! !ROUTINE: ods_do3lev
!
! !INTERFACE:
!      
      subroutine ods_do3lev ( var, data, ninfo, nobs, 
     .                       ods,  mobs, kobs, iks, rc )

      use m_odsmeta
      use m_ods
      use m_ods_obsdiags, only : ods_obsdiags
      use m_ods_obsdiags, only : ods_obsdiags_getparam
      
      implicit NONE

! !INPUT PARAMETERS:

      integer,          intent(in)    :: mobs         ! max obs in ODS array
      integer,          intent(in)    :: ninfo, nobs  ! data dims
      real(4),          intent(in)    :: data(ninfo,nobs)
      character(len=*), intent(in)    :: var          ! obs type
             
! !INPUT/OUTPUT PARAMETERS:

      integer,          intent(inout) :: kobs, iks
      type(ods_vect),   intent(inout) :: ods


! !OUTPUT PARAMETERS:

      integer,          intent(out)   :: rc
 	            
! !DESCRIPTION: Converts ozone level data from GSI diag_ files to ODS.
!
! !REVISION HISTORY:
!	4Jan2007 - Sienkiewicz   Initial code
!       9Apr2008 - Sienkiewicz   modified for new GSI-2006
!
!EOP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
            
      character(len=*), parameter :: myname = 'ods_do3lev'

      integer, parameter   :: r_single = selected_real_kind(6)
      real(r_single)  zero_single,tiny_single
      real small_num
      
      integer i, iu, iv, kobs0, ks0

      real, parameter :: undef = 1.e15

      logical lobsdiagsave
      
!     Pointers to ods attributes:
!     --------------------------
      integer, pointer :: kt(:)       ! data type index
      integer, pointer :: kx(:)       ! data source index
      integer, pointer :: ks(:)       ! sounding index
      real   , pointer :: lon(:)      ! longitude of obs (degrees)
      real   , pointer :: lat(:)      ! latitude of obs (degrees)
      real   , pointer :: lev(:)      ! pressure level of obs (hPa)
      integer, pointer :: time(:)     ! time
      real   , pointer :: obs(:)      ! observation
      real   , pointer :: OmF(:)      ! obs-minus-fcst (O-F)
      real   , pointer :: OmA(:)      ! obs-minus-ana  (O-A)
      real   , pointer :: sigo(:)     ! obs error std dev
      real   , pointer :: xm(:)       ! metadata
      integer, pointer :: qcx(:)      ! exclusion mark
      integer, pointer :: qch(:)      ! history mark

!     Assign pointers to ODS attributes
!     ---------------------------------
      kt   => ods%data%kt
      kx   => ods%data%kx
      ks   => ods%data%ks
      lon  => ods%data%lon
      lat  => ods%data%lat
      lev  => ods%data%lev
      time => ods%data%time
      obs  => ods%data%obs
      OmF  => ods%data%OmF
      OmA  => ods%data%OmA
      sigo => ods%data%Xvec
      xm   => ods%data%xm
      qcx  => ods%data%qcexcl
      qch  => ods%data%qchist

      tiny_single = tiny(zero_single)
      small_num = 1.e3 * tiny_single

      rc = 0  ! all is well to start with

      call ods_obsdiags_getparam ( 'lobsdiagsave', lobsdiagsave )


      kobs0 = kobs

      ks0 = data(11,1)

      do i = 1, nobs

         kobs         = kobs + 1
         if(kobs>mobs) then
            rc = 99
            print *, myname, 'dim violation, var = ', var
            exit
         endif
         kt(kobs)     = kto3
         kx(kobs)     = nint(data(1,i)) 
         
         lat(kobs)    = data(2,i)
         lon(kobs)    = data(3,i)
         lev(kobs)    = data(4,i)
         time(kobs)   = data(5,i) * 60 ! convert to minutes from ana time
         obs(kobs)    = data(8,i)             ! obs
         omf(kobs)    = (data(8,i)-data(9,i)) ! omf
         oma(kobs)    = undef   ! no info available
         xm(kobs)     = data(13,i)       ! try quality flag
         
         qch(kobs)    = 0       ! no info available
         qcx(kobs)    = 0


         if (data(15,i) < 0. .or. data(15,i) >= 100. ) then
            qcx(kobs) = X_PASSIVE
         endif


         if (data(7,i)>small_num) then
            sigo(kobs) = 1./data(7,i)                ! sigo 
            if(data(6,i)<1) qcx(kobs) = X_NCEP_NLNQC
         else                                         ! rejected by QC
            sigo(kobs) = undef
            qcx(kobs) = 2   
         endif

         if (data(11,i) /= ks0 ) then    ! check sounding index
            iks = iks + 1
            ks0 = data(11,i)
         endif

         ks(kobs)     = iks

c$$$         weight(kobs) = data(6,i) ! fractional weight

      end do

      
      
!     Nullify pointers
!     -----------------
      if(associated(lat)   ) nullify(lat)
      if(associated(lon)   ) nullify(lon)
      if(associated(lev)   ) nullify(lev)
      if(associated(time)  ) nullify(time)
      if(associated(kt)    ) nullify(kt)
      if(associated(kx)    ) nullify(kx)
      if(associated(ks)    ) nullify(ks)
      if(associated(xm)    ) nullify(xm)
      if(associated(qch)   ) nullify(qch)
      if(associated(qcx)   ) nullify(qcx)
      if(associated(obs)   ) nullify(obs)
      if(associated(omf)   ) nullify(omf)
      if(associated(oma)   ) nullify(oma)
      if(associated(sigo)  ) nullify(sigo)
	  
      end subroutine ods_do3lev
