!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!BOP
! !ROUTINE: ods_drad
!
! !INTERFACE:
!      
      subroutine ods_drad ( obstype, isat, iuse, iinfo, pinfo, pdata,
     .                      varchn, nuchan,
     .                      ndiag, iint, ireal, ipchan, npred, nchanl,
     .  		    ods, mobs, kobs, iks, ioff, rc )
      use m_odsmeta
      use m_ods
      use m_ods_obsdiags, only : ods_obsdiags
      use m_ods_obsdiags, only : ods_obsdiags_getparam

      implicit NONE

! !INPUT PARAMETERS:

      character(len=*), intent(in)    :: obstype
      integer,          intent(in)    :: isat
      integer,          intent(in)    :: ireal, ipchan, npred, nchanl, iint
      integer,          intent(in)    :: mobs, iks 
      integer,          intent(in)    :: ndiag
      integer,          intent(in)    :: ioff
      integer,          intent(in)    :: iuse(nchanl)
      integer,          intent(in)    :: iinfo(iint)
      real(4),          intent(in)    :: pinfo(ireal)
      real(4),          intent(in)    :: pdata(ndiag,nchanl)
      integer,          intent(in)    :: nuchan(nchanl)
      real(4),          intent(in)    :: varchn(nchanl)
             
! !INPUT/OUTPUT PARAMETERS:

      integer,          intent(inout) :: kobs      ! current obs counter
      type(ods_vect),   intent(inout) :: ods

! !OUTPUT PARAMETERS:

      integer,          intent(out)   :: rc
 	            
! !DESCRIPTION: Converts radiance data from GSI diag_ files to ODS.
!
! !REVISION HISTORY:
!	20Dec2004 - D. Dee - Initial code.
!       02Mar2005 - Dee - Updated for latest GSI
!       20May2005 - Dee - Handle passive radiance data
!       03Nov2005 - Sienkiewicz - addition for SSM/I radiance 
!                                  (has some QC marks -> history marks)
!       24Jan2006 - Todling - Bug fix; defining xvec w/ undef's
!       31Jan2006 - Todling - Filling in time attribute 
!        9Mar2006 - Sienkiewicz - restore xvec ( = sigO) and modify
!                                 sigO screening criteria
!       24Jul2006 - Sienkiewicz - Corrected bug with time attribute 
!       17Oct2006 - Sienkiewicz - add SSU
!       14Dec2006 - Todling     - update to gsi-2006_09
!       22Dec2006 - Sienkiewicz - update SSM/I history mark assignment
!       18Jan2008 - Todling     - Add ability to read in observation impacts
!       17Feb2012 - Todling     - Update to Sep2011 GSI
!       08Mar2012 - Todling     - Use pdata(5,i)<0 to identify passive data
!       29Mar2012 - Bloom/RT    - bug; sigo's were getting sqrt'd twice!
!       02Apr2014 - Todling     - when desired, write out table sigo to file
!       16Dec2016 - Todling     - add lreduced to indicate reduced diag
!                                 (use real chan number)
!EOP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
            
      character(len=*), parameter :: myname = 'ods_drad'

      integer, parameter   :: r_single = selected_real_kind(6)
      real(r_single)  zero_single,tiny_single
      real small_num,nlomx,tlomx

      logical satknown,lobsdiagsave,passed,ladjsigo,lreduced
      integer i, n, kidsat, kobs0, miter

      real, parameter :: undef = 1.e15
       
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
      
      call ods_obsdiags_getparam ( 'ladjsigo', ladjsigo )
      call ods_obsdiags_getparam ( 'lobsdiagsave', lobsdiagsave )
      if ( lobsdiagsave ) then
           call ods_obsdiags_getparam ( 'miter', miter )
           if (ndiag/=ipchan+npred+2+4*miter+1) then
              rc = 99
              print *, myname, ': dim violation(obs sens), obstype = ', obstype
              return
           endif
           if(ioff.ne.ipchan+npred+2) then
              rc = 99
              print *, myname, ': sensitivity index appears to be incorrect = ', ioff,ipchan+npred+2
              return
           endif
      endif
      call ods_obsdiags_getparam ( 'reduce_diag', lreduced )

      satknown = .false.
      kidsat   = isat
      do n = 1, nsats
         if(trim(sats(n))==trim(obstype)) then
            satknown = .true.
            kidsat   = kidsat+idsats(n)
            exit
         else
           cycle
         endif
      end do

      if (.not.satknown .or. kidsat==0) then
        rc = 98
        print *, myname, ': Cannot identify satellite type = ', obstype, ' kidsat = ', kidsat
        return
      endif
      
      kobs0 = kobs
      do i = 1, nchanl

	kobs	   = kobs + 1
        if (kobs>mobs) then
           rc = 98
           print *, myname, ': dim violation, obstype = ', obstype
           exit
        endif
        kt(kobs)   = ktTb		  ! brightness temperature
        kx(kobs)   = kidsat
        ks(kobs)   = iks
  
        lat(kobs)  = pinfo(1)
        lon(kobs)  = pinfo(2)
        if (lreduced) then
           lev(kobs)  = nuchan(i) ! in case of reduced diag output, need to use 
                                  ! true channel number
        else
           lev(kobs)  = i ! nuchan(i) ---> back to index for now (RTodling)
        endif
        time(kobs) = pinfo(4) * 60	  ! time
        obs(kobs)  = pdata(1,i)           ! observed brightness temperature
        omf(kobs)  = pdata(2,i)   	  ! bias-corrected omf
        oma(kobs)  = undef		  ! no info available
        xm(kobs)   = pdata(3,i)-pdata(2,i)! bias correction
  
        qch(kobs)  = 0  		  ! no info available
        qcx(kobs)  = 0
	if (iuse(i)==-1.or.pdata(5,i)<0) qcx(kobs) = 1    ! passive data
	if (pdata(4,i)>small_num) then
           if (ladjsigo) then
              sigo(kobs) = 1.0/pdata(4,i)                 ! adjusted   sigo
           else
              sigo(kobs) = sqrt(varchn(i))                ! prescribed sigo
           endif
	else                                      ! rejected by GSI
	   sigo(kobs) = undef
	   if (qcx(kobs)==0) qcx(kobs) = 2    
	endif

!       Place obs impacts in sigo slot
!       ------------------------------
        if ( lobsdiagsave ) then
             sigo(kobs) = 0.0
             call ods_obsdiags ( nlomx, tlomx, sigo(kobs), pdata, ioff, i, ndiag, nchanl, undef, passed )
             if(passed) qcx(kobs) = 0
        endif

      end do
      
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
	  
      end subroutine ods_drad
