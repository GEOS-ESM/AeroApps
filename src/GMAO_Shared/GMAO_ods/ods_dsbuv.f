!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!BOP
! !ROUTINE: ods_dsbuv
!
! !INTERFACE:
!      
      subroutine ods_dsbuv ( isat, obstype, pinfo, pdata, 
     .                       pobs, tnoise, iouse,
     .                       ireal, nlevs, ndiag,
     .  		     ods, mobs, kobs, iks, rc )
      use m_odsmeta
      use m_ods
      use m_ods_obsdiags, only : ods_obsdiags
      use m_ods_obsdiags, only : ods_obsdiags_getparam

      implicit NONE

! !INPUT PARAMETERS:

      integer,          intent(in)    :: isat
      integer,          intent(in)    :: ireal, nlevs, ndiag
      integer,          intent(in)    :: mobs, iks 
      integer,          intent(in)    :: iouse(nlevs)
      real(4),          intent(in)    :: pinfo(ireal)
      real(4),          intent(in)    :: pdata(ndiag,nlevs)
      real(4),          intent(in)    :: pobs(nlevs)
      real(4),          intent(in)    :: tnoise(nlevs)
      character(len=*), intent(in)    :: obstype
             
! !INPUT/OUTPUT PARAMETERS:

      integer,          intent(inout) :: kobs      ! current obs counter
      type(ods_vect),   intent(inout) :: ods

! !OUTPUT PARAMETERS:

      integer,          intent(out)   :: rc
 	            
! !DESCRIPTION: Converts sbuv data from GSI diag_ files to ODS.
!
! !REVISION HISTORY:
!	06Jan2005 - D. Dee - Initial code.
!       02May2005 - D. Dee - Fix for sbuv/2
!       24Jan2006 - Todling - Bug fix: xvec now defined w/ undef's
!        9Mar2006 - Sienkiewicz   take out xvec change, modify test
!                                   of sigO value
!       14Dec2006 - Todling - update to gsi-2006_09
!       08Apr2006 - Sienkiewicz - obs_sens test now in if(lobsdiagsave) block
!       26Aug2010 - Todling - update to Mar2010 GSI (obs impact handle)
!
!EOP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
            
      character(len=*), parameter :: myname = 'ods_dsbuv'


      integer, parameter   :: r_single = selected_real_kind(6)
      real(r_single)  zero_single,tiny_single
      real small_num, nlomx, tlomx
      
      logical satknown, lobsdiagsave, passed
      integer i, n, kobs0, kidsat, ioff, miter

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

      call ods_obsdiags_getparam ( 'lobsdiagsave', lobsdiagsave )
      if ( lobsdiagsave ) then
        call ods_obsdiags_getparam ( 'miter', miter )
        ioff=6
        if (ndiag/=ioff+4*miter+1) then
            rc = 99
            print *,'ndiag is ',ndiag
            print *,'ioff+4*miter+1 is ',ioff+4*miter+1
            print *, myname, ': dim violation(obs sens), obstype = ', obstype
            return
        endif
      endif

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
      do i = 1, nlevs
      
	kobs	   = kobs + 1
        if (kobs>mobs) then
           rc = 99
           print *, myname, 'dim violation, obstype = ', 'sbuv2'
           exit
        endif
	select case(obstype)
	case('o3lev','mls','mls20','mls22','mls30','mls55')
         kt(kobs)   = kto3mx		  ! ozone mixing ratio (ppmv)
	 lev(kobs)  = pdata(4,i)
	case default
	 if (pobs(i)>0) then
           kt(kobs)   = kto3		  ! ozone
           lev(kobs)  = pobs(i)           ! bottom edge of layer
	   if (i>1) then
	      xm(kobs) = pobs(i-1)        ! top edge of layer
	   else                           ! NOTE: assumes profile order
	      xm(kobs) = 0.0
	   end if
	 else
	   kt(kobs)   = kttco3		  ! total column ozone
	   lev(kobs)  = undef
	 end if
	end select

        kx(kobs)   = kidsat
        ks(kobs)   = iks
  
        lat(kobs)  = pinfo(1)
        lon(kobs)  = pinfo(2)
        time(kobs) = int(pinfo(3)*60.0)   ! minutes from ana time	 
        obs(kobs)  = pdata(1,i)           ! obs
        omf(kobs)  = pdata(2,i)	          ! omf
        oma(kobs)  = undef		  ! no info available
  
        qch(kobs)  = 0  		  ! no info available
        qcx(kobs)  = 0
	if (pdata(3,i)>small_num .and. iouse(i)>0) then
           sigo(kobs) = 1.0/pdata(3,i)    ! sigo
	else				  ! rejected by QC
	   sigo(kobs) = undef
	   qcx(kobs) = 2   
	endif
		
        if ( lobsdiagsave ) then
           sigo(kobs) = 0.0
           call ods_obsdiags ( nlomx, tlomx, sigo(kobs), pdata, ioff, i, ndiag, nlevs, undef, passed )
           if(passed) qcx(kobs) = 0
        endif

!	write(*,'(i2,i3,i5,i6,i5,3(1x,f8.3),1x,i3,6(1x,e8.3e1))') qcx(kobs),i,kobs,kx(kobs),kt(kobs),
!     &		lat(kobs),lon(kobs),lev(kobs),time(kobs),obs(kobs),omf(kobs),oma(kobs),sigo(kobs)
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
	  
      end subroutine ods_dsbuv
