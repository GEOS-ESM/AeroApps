!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  NASA/GSFC, Global Modeling and Assimilation Office, 610.1 GEOS/DAS   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!BOP
! !ROUTINE: ods_dpcp
!
! !INTERFACE:
!      
      subroutine ods_dpcp ( obstype, isat, iinfo, pinfo,
     .                      iint, ireal,
     .                      ods, mobs, kobs, iks, rc )
      use m_odsmeta
      use m_ods
      use m_ods_obsdiags, only : ods_obsdiags
      use m_ods_obsdiags, only : ods_obsdiags_getparam

      implicit NONE

! !INPUT PARAMETERS:

      character(len=*), intent(in)    :: obstype
      integer,          intent(in)    :: isat         ! offset satid value for KX
      integer,          intent(in)    :: ireal        ! # real vars in pinfo
      integer,          intent(in)    :: iint         ! # int  vars in iinfo
      integer,          intent(in)    :: mobs         ! size of ODS arrays
      integer,          intent(in)    :: iks          ! ks for data type
      integer,          intent(in)    :: iinfo(iint)  ! integer diagnostics
      real(4),          intent(in)    :: pinfo(ireal) ! real diagnostics
             
! !INPUT/OUTPUT PARAMETERS:

      integer,          intent(inout) :: kobs      ! current obs counter
      type(ods_vect),   intent(inout) :: ods

! !OUTPUT PARAMETERS:

      integer,          intent(out)   :: rc

! !DESCRIPTION: Converts precip data from GSI diag_ files to ODS.
!
! !REVISION HISTORY:
!	07Nov2005 - Sienkiewicz - Initial code based on 'ods_drad'
!       01Feb2006 - Sienkiewicz, Zhang - corrected time value
!       01Feb2007 - Sienkiewicz - new usage flag
!       19Dec2007 - Todling     - Add ability to read in observation impacts
!       06Feb2009 - Sienkiewicz/RT - correct obs time (mult by 360)
!       17Feb2009 - Todling     - handle sensitivity as well as impact
!       17Feb2009 - Todling     - update to handle sensitivity from GSI-Dec08
!       10Mar2012 - Todling     - add use_threshold to sync w/ setuppcp
!                               - fix calc of omf; now as setuppcp does it
!
!EOP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
            
      character(len=*), parameter :: myname = 'ods_dpcp'

      integer, parameter   :: r_single = selected_real_kind(6)
      real(r_single)  zero_single,tiny_single
      real small_num
      real nlomx, tlomx, obimp
      logical lobsdiagsave, lobssens, passed
      
      integer kidsat,idia

      real, parameter :: undef = 1.e15
      real, parameter :: use_threshold = 1.e-6

      real, parameter :: pcpconv = 24.  !  convert mm/hr to mm/day
       
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
      small_num = 10. * tiny_single


      rc = 0  ! all is well to start with
      call ods_obsdiags_getparam ( 'lobsdiagsave', lobsdiagsave )
      call ods_obsdiags_getparam ( 'lobssens', lobssens )

      kidsat = iinfo(1)          ! "sat id" is 1st element of array
            

        kobs = kobs + 1
        if (kobs>mobs) then
           rc = 99
           print *, myname, ' dim violation, obstype = ', obstype, kobs, mobs
           return
        endif
        kt(kobs)   = ktpr        ! precip rate
        kx(kobs)   = kidsat + isat
        ks(kobs)   = iks
  
        lat(kobs)  = pinfo(1)
        lon(kobs)  = pinfo(2)
        lev(kobs)  = 1
        time(kobs) = pinfo(3) * 360.                   ! analysis divides by 3hrs
        obs(kobs)  = pinfo(5) * pcpconv                ! obs w/o bias corr.
        omf(kobs)  = log(1.+pinfo(5))-log(1.+pinfo(18))! bias-corrected omf
        oma(kobs)  = undef                             ! no info available
        xm(kobs)   = (pinfo(18)-pinfo(17)) * pcpconv   ! bias correction
  
        qch(kobs)  = 0                                 ! no info available
        qcx(kobs)  = 0
        if (iinfo(3) < 0 ) then
           qcx(kobs) = X_PASSIVE
        endif

        if (pinfo(20)>use_threshold) then
           sigo(kobs) = 1.0/sqrt(pinfo(20)) ! sigo
        else                                ! rejected by GSI
           sigo(kobs) = undef
           qcx(kobs) = 2    
!
! Check pinfo(19) value to assign history flag
!    pinfo(20) ~0 and pinfo(19) > 0 :  eliminated by gross check (qch=40)
!                     pinfo(19) ~0  :  ice, or sensitivity/gradient problem
!
!    iinfo(4) =2 : sfcflg indicates ice  (qch=52)
!    otherwise     eliminated due to sensitivity/gradient problem (qch=31)
!
           if (pinfo(19) > small_num) then
              qch(kobs) = 40
           else if ( iinfo(4) == 2 ) then
              qch(kobs) = 52 
           else
              qch(kobs) = 31
           endif
        endif

        if (lobsdiagsave) then              ! place obs sensitivity in sigo slot
            idia = 22
            obimp =  0.0
            sigo(kobs)=obimp
            call ods_obsdiags ( nlomx, tlomx, obimp, reshape(pinfo,(/ireal,1/)), idia, 1, ireal, 1, undef, passed )
            if (passed) then
                sigo(kobs)=obimp
                qcx (kobs)=0
                if(lobssens) sigo(kobs)=obimp/(1.+obs(kobs))
            endif
        endif

      
!     Nullify pointers
!     -----------------
      if(associated(lat) ) nullify(lat)
      if(associated(lon) ) nullify(lon)
      if(associated(lev) ) nullify(lev)
      if(associated(time)) nullify(time)
      if(associated(kt)  ) nullify(kt)
      if(associated(kx)  ) nullify(kx)
      if(associated(ks)  ) nullify(ks)
      if(associated(xm)  ) nullify(xm)
      if(associated(qch) ) nullify(qch)
      if(associated(qcx) ) nullify(qcx)
      if(associated(obs) ) nullify(obs)
      if(associated(omf) ) nullify(omf)
      if(associated(oma) ) nullify(oma)
      if(associated(sigo)) nullify(sigo)
  
      end subroutine ods_dpcp
