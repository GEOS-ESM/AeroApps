!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!BOP
! !ROUTINE: ods_dconv
!
! !INTERFACE:
!      
      subroutine ods_dconv ( var, cdiagbuf, data, ninfo, nobs, 
     .                       ods, station, weight, mobs, kobs, iks, ioff, rc )

      use m_odsmeta
      use m_ods
      use m_ods_obsdiags, only : ods_obsdiags
      use m_ods_obsdiags, only : ods_obsdiags_getparam
      use m_fpe, only: isInf
      
      implicit NONE

! !INPUT PARAMETERS:

      integer,          intent(in)    :: mobs
      integer,          intent(in)    :: ninfo, nobs, ioff
      real(4),          intent(in)    :: data(ninfo,nobs)
      character(len=*), intent(in)    :: var
             
! !INPUT/OUTPUT PARAMETERS:

      integer,          intent(inout) :: kobs, iks
      character(len=*), intent(inout) :: cdiagbuf(nobs)
      type(ods_vect),   intent(inout) :: ods
      real(4),          intent(inout) :: weight(mobs)
      character(len=*), intent(inout) :: station(mobs)

! !OUTPUT PARAMETERS:

      integer,          intent(out)   :: rc
 	            
! !DESCRIPTION: Converts conventional data from GSI diag_ files to ODS.
!
! !REVISION HISTORY:
!	19Apr2004 - R. Todling - Initial code.
!       03Dec2004 - D. Dee     - remove kx conversion; various corrections;
!                                pass fractional weight
!       16Mar2005 - Meta       - Added sst
!       25May2005 - D. Dee     - Changed qcx values
!       07Nov2005 - Todling    - Updated to handle PW
!       24Jan2006 - Todling    - Bug fix; defining xvec w/ undef's
!       31Jan2006 - Todling    - Added check on varQC
!       24Feb2006 - Steve B, Meta - Modify 'ps' due to change in diag file.
!        8Mar2006 - Meta       - undo xvec undef, modify check on obs error
!       15Mar2006 - Meta        isNaN check on 'sst' sigO work around beta6's
!                                bug in SST diag file write 
!       30Oct2006 - Meta       - added 'muse' to diag output, use (if
!                                  available) to determine passive obs
!       14Dec2006 - Todling    - update to work w/ GSI-2006_09;
!                                Note: obs slot keeps bias corrected obs
!       19Dec2007 - Todling    - Add ability to read in observation impacts
!       17Feb2009 - Todling    - Add handle to store sensitivies rather than impacts
!       17Mar2009 - Meunier    - Add handle for Lagragian-type observations
!       05Apr2009 - Todling    - Add GPS refractivity and bending angle
!       29Jan2010 - Meta       - modified GPS (o-f) for change in diag file in Jun2009 GSI
!       01Apr2010 - McCarty/RT - add doppler wind lidar 
!       26Jul2010 - Todling    - add tropical cycle knob
!       14Feb2012 - Todling    - idia from 19 to 20 (per GSI update)
!       28Jun2012 - Meta       - screen gpsbnd OmF for invalid values
!        2Aug2012 - Meta       - identify drifting buoy from subtype & stn ID
!        4Sep2012 - Meta       - switch gpsbnd test to use m_fpe isInf
!       17May2013 - Meta       - write temperature bias correction to xm slot
!       25Jun2013 - Todling    - upd idia(t)=20 after addition of phase-of-flight
!                                   to temperature diags
!       21Jan2014 - Todling    - option to choose type of sigo to go into ODS file
!       28Jan2014 - Todling    - use sensitivity slot (idia) from header
!
!EOP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
            
      character(len=*), parameter :: myname = 'ods_dconv'

      integer, parameter   :: r_single = selected_real_kind(6)
      real(r_single)  zero_single,tiny_single
      real small_num, nlomx(2), tlomx(2), obimp(2)
      logical :: lobsdiagsave, passed, lobssens, ladjsigo
      
      integer i, iu, iv, kobs0, idia, ilat, ilon, ios, iwmo, n311, isigo_slot

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
      small_num = 1.e3 * tiny_single

      rc = 0  ! all is well to start with

      call ods_obsdiags_getparam ( 'lobsdiagsave', lobsdiagsave )
      call ods_obsdiags_getparam ( 'lobssens', lobssens )
      call ods_obsdiags_getparam ( 'ladjsigo', ladjsigo )

      if ( ladjsigo ) then
          isigo_slot=16 ! GSI-adjusted ob error
      else
          isigo_slot=14 ! prepbufr (prescribed) ob error
      endif

      kobs0 = kobs
      select case (var)
      
      case (' dw')  ! doppler wind lidar
      
          do i = 1, nobs

	      kobs         = kobs + 1
              if(kobs>mobs) then
                 rc = 99
                 print *, myname, ' dim violation, var = ', var
                 exit
              endif
              kt(kobs)     = ktdw
              kx(kobs)     = nint(data(1,i)) 
  
              lat(kobs)    = data(3,i)
              lon(kobs)    = data(4,i)
              lev(kobs)    = data(6,i)
              time(kobs)   = data(8,i) * 60       ! minutes from ana time
              obs(kobs)    = data(17,i)           ! observation
              omf(kobs)    = data(18,i)           ! omf used in the analysis
              oma(kobs)    = undef
              xm(kobs)     = obs(kobs)-omf(kobs)+data(19,i) ! original obs w/o bias correction
  
	      qch(kobs)    = 0                    ! no info available
              qcx(kobs)    = 0

              if (ninfo >= 12) then
                 if (data(12,i) <  0. ) then
                    qcx(kobs) = X_PASSIVE
                 endif
              endif

	      if (data(16,i)>small_num) then
                 sigo(kobs) = 1.0/data(isigo_slot,i) ! sigo used in analysis
                 if(data(13,i)<1) qcx(kobs) = X_NCEP_NLNQC
	      else                                ! rejected by QC
	         sigo(kobs) = undef
		 qcx(kobs) = 2   
	      endif

              ks(kobs)     = iks
              station(kobs)= cdiagbuf(i)
	      weight(kobs) = data(13,i)           ! fractional weight

              if (lobsdiagsave) then              ! place obs sensitivity in sigo slot
                  idia  = ioff
                  obimp = 0.0
                  sigo(kobs) = 0.0
                  call ods_obsdiags ( nlomx(1), tlomx(1), obimp(1), data, idia, i, ninfo, nobs, undef, passed )  
                  if(passed) then
                     sigo(kobs) = obimp(1)
                     qcx (kobs) = 0
                  endif
              endif

          end do

      case (' ps')  ! surface pressure
                    ! Note: It is really ln(ps) that is analyzed
      
          do i = 1, nobs

	      kobs         = kobs + 1
              if(kobs>mobs) then
                 rc = 99
                 print *, myname, 'dim violation, var = ', var
                 exit
              endif
              kt(kobs)     = ktps2m               ! surface pressure
              kx(kobs)     = nint(data(1,i))

              if (kx(kobs) == 180 .and. nint(data(2,i)) == 562) then
                 read(cdiagbuf(i),*,iostat=ios) iwmo
                 if (ios == 0 .and. iwmo > 0) then
                    if (mod(iwmo,1000) >=500) then
                       kx(kobs) = 199
                    end if
                 end if
              end if
  
              lat(kobs)    = data(3,i)
              lon(kobs)    = data(4,i)
              lev(kobs)    = undef
              time(kobs)   = data(8,i) * 60       ! minutes from ana time
              obs(kobs)    = data(17,i)           ! observation
              omf(kobs)    = data(18,i)           ! omf used in analysis
              oma(kobs)    = undef                ! no info available
              xm(kobs)     = data(7,i)            ! surface height
  
              qch(kobs)    = 0                    ! no info available
              qcx(kobs)    = 0

              if (ninfo >= 12) then
                 if (data(12,i) < 0. ) then
                    qcx(kobs) = X_PASSIVE
                 endif
              endif

	      if (data(16,i)>small_num) then
                 sigo(kobs) = 1.0/data(isigo_slot,i) ! sigo (hPa)
                 if(data(13,i)<1) qcx(kobs) = X_NCEP_NLNQC
	      else                                ! rejected by QC
	         sigo(kobs) = undef
		 qcx(kobs) = 2   
	      endif
	      
              ks(kobs)     = iks
              station(kobs)= cdiagbuf(i)
	      weight(kobs) = data(13,i)           ! fractional weight

              if (lobsdiagsave) then              ! place obs sensitivity in sigo slot
                  idia  = ioff
                  obimp = 0.0
                  sigo(kobs) = 0.0
                  call ods_obsdiags ( nlomx(1), tlomx(1), obimp(1), data, idia, i, ninfo, nobs, undef, passed )  
                  if(passed) then
                     sigo(kobs) = obimp(1)
                     qcx (kobs) = 0
                  endif
              endif

          end do
    
      case ('tcp')  ! tropical cyclone simulated surface pressure
      
          do i = 1, nobs

	      kobs         = kobs + 1
              if(kobs>mobs) then
                 rc = 99
                 print *, myname, 'dim violation, var = ', var
                 exit
              endif
              kt(kobs)     = ktps2m               ! surface pressure
              kx(kobs)     = nint(data(1,i)) 
  
              lat(kobs)    = data(3,i)
              lon(kobs)    = data(4,i)
              lev(kobs)    = undef
              time(kobs)   = data(8,i) * 60       ! minutes from ana time
              obs(kobs)    = data(17,i)           ! observation
              omf(kobs)    = data(18,i)           ! omf used in analysis
              oma(kobs)    = undef                ! no info available
              xm(kobs)     = data(7,i)            ! surface height
  
              qch(kobs)    = 0                    ! no info available
              qcx(kobs)    = 0

              if (ninfo >= 12) then
                 if (data(12,i) < 0. ) then
                    qcx(kobs) = X_PASSIVE
                 endif
              endif

	      if (data(16,i)>small_num) then
                 sigo(kobs) = 1.0/data(isigo_slot,i) ! sigo (hPa)
                 if(data(13,i)<1) qcx(kobs) = X_NCEP_NLNQC
	      else                                ! rejected by QC
	         sigo(kobs) = undef
		 qcx(kobs) = 2   
	      endif
	      
              ks(kobs)     = iks
              station(kobs)= cdiagbuf(i)
	      weight(kobs) = data(13,i)           ! fractional weight

              if (lobsdiagsave) then              ! place obs sensitivity in sigo slot
                  idia  = ioff
                  obimp = 0.0
                  sigo(kobs) = 0.0
                  call ods_obsdiags ( nlomx(1), tlomx(1), obimp(1), data, idia, i, ninfo, nobs, undef, passed )  
                  if(passed) then
                     sigo(kobs) = obimp(1)
                     qcx (kobs) = 0
                  endif
              endif

          end do
    
      case ('  q')  ! specific humidity
                    ! Note: it is really pseudo-rh that is analyzed

          do i = 1, nobs

	      kobs         = kobs + 1
              if(kobs>mobs) then
                 rc = 99
                 print *, myname, 'dim violation, var = ', var
                 exit
              endif
              kt(kobs)     = ktqq                 ! specific humidity
              kx(kobs)     = nint(data(1,i)) 

              if (kx(kobs) == 180 .and. nint(data(2,i)) == 562) then
                 read(cdiagbuf(i),*,iostat=ios) iwmo
                 if (ios == 0 .and. iwmo > 0) then
                    if (mod(iwmo,1000) >=500) then
                       kx(kobs) = 199
                    end if
                 end if
              end if
  
              lat(kobs)    = data(3,i)
              lon(kobs)    = data(4,i)
              lev(kobs)    = data(6,i)
              time(kobs)   = data(8,i) * 60       ! convert to minutes from ana time
              obs(kobs)    = 1e3*data(17,i)       ! observation
              omf(kobs)    = 1e3*data(18,i)       ! omf used in analysis [g/kg]
              oma(kobs)    = undef                ! no info available
              xm(kobs)     = 1e3*data(20,i)       ! guess saturated mixing ratio [g/kg]
  
	      qch(kobs)    = 0                    ! no info available
              qcx(kobs)    = 0

              if (ninfo >= 12) then
                 if (data(12,i) < 0. ) then
                    qcx(kobs) = X_PASSIVE
                 endif
              endif

	      if (data(16,i)>small_num) then
                 sigo(kobs) = 1e3/data(isigo_slot,i) ! actual sigo used in analysis (for sp hum, in g/kg)
                 if(data(13,i)<1) qcx(kobs) = X_NCEP_NLNQC
	      else                                ! rejected by QC
	         sigo(kobs) = undef
		 qcx(kobs) = 2   
	      endif

              ks(kobs)     = iks
              station(kobs)= cdiagbuf(i)
	      weight(kobs) = data(13,i)           ! fractional weight

              if (lobsdiagsave) then              ! place obs sensitivity in sigo slot
                  idia  = ioff
                  obimp = 0.0
                  sigo(kobs) = 0.0
                  call ods_obsdiags ( nlomx(1), tlomx(1), obimp(1), data, idia, i, ninfo, nobs, undef, passed )  
                  if(passed) then
                     sigo(kobs) = obimp(1)
                     qcx (kobs) = 0
                     if(lobssens) sigo(kobs)=sigo(kobs)*1.e-3 ! convert to 1/[g/kg]
                  endif 
              endif

          end do

      case ('spd')  ! wind speed
      
          do i = 1, nobs
              
	      kobs         = kobs + 1
              if(kobs>mobs) then
                 rc = 99
                 print *, myname, ' dim violation, var = ', var
                 exit
              endif
              kt(kobs)     = ktus10               ! define it as 10m speeds
              kx(kobs)     = nint(data(1,i)) 
  
              lat(kobs)    = data(3,i)
              lon(kobs)    = data(4,i)
              lev(kobs)    = undef
              time(kobs)   = data(8,i) * 60       ! minutes from ana time
              obs(kobs)    = data(17,i)           ! observation
              omf(kobs)    = data(18,i)           ! omf used in analysis
              oma(kobs)    = undef                ! no info available
              xm(kobs)     = data(7,i)            ! elevation
  
	      qch(kobs)    = 0                    ! no info available
              qcx(kobs)    = 0

              if (ninfo >= 12) then
                 if (data(12,i) <  0. ) then
                    qcx(kobs) = X_PASSIVE
                 endif
              endif

	      if (data(16,i)>small_num) then
                 sigo(kobs) = 1.0/data(isigo_slot,i) ! sigo used in analysis
                 if(data(13,i)<1) qcx(kobs) = X_NCEP_NLNQC
	      else                                ! rejected by QC
	         sigo(kobs) = undef
		 qcx(kobs) = 2   
	      endif

              ks(kobs)     = iks
              station(kobs)= cdiagbuf(i)
	      weight(kobs) = data(13,i)           ! fractional weight

              if (lobsdiagsave) then              ! place obs sensitivity in sigo slot
                  idia  = ioff
                  obimp = 0.0
                  sigo(kobs) = 0.0
                  call ods_obsdiags ( nlomx(1), tlomx(1), obimp(1), data, idia, i, ninfo, nobs, undef, passed )  
                  if(passed) then
                     sigo(kobs) = obimp(1)
                     qcx (kobs) = 0
                  endif
              endif

          end do
      
      case ('  t')  ! virtual temperature
      
	  n311=0
          do i = 1, nobs

	      kobs         = kobs + 1
              if(kobs>mobs) then
                 rc = 99
                 print *, myname, ' dim violation, var = ', var
                 exit
              endif
              kt(kobs)     = ktTv
              kx(kobs)     = nint(data(1,i)) 
  
              if (kx(kobs) == 180 .and. nint(data(2,i)) == 562) then
                 read(cdiagbuf(i),*,iostat=ios) iwmo
                 if (ios == 0 .and. iwmo > 0) then
                    if (mod(iwmo,1000) >=500) then
                       kx(kobs) = 199
                    end if
                 end if
	      elseif(kx(kobs)==311) then
	         n311=n311+1
	         kx(kobs)=304		! a fix/hack for MLS temperature
              end if

              lat(kobs)    = data(3,i)
              lon(kobs)    = data(4,i)
              lev(kobs)    = data(6,i)
              time(kobs)   = data(8,i) * 60       ! minutes from ana time
              obs(kobs)    = data(17,i)           ! observation
              omf(kobs)    = data(18,i)           ! omf used in the analysis
              oma(kobs)    = undef
              xm(kobs)     = data(19,i) - data(18,i) ! omgnbc - omgbc
  
	      qch(kobs)    = 0                    ! no info available
              qcx(kobs)    = 0

              if (ninfo >= 12) then
                 if (data(12,i) <  0. ) then
                    qcx(kobs) = X_PASSIVE
                 endif
              endif

	      if (data(16,i)>small_num) then
                 sigo(kobs) = 1.0/data(isigo_slot,i) ! sigo used in analysis
                 if(data(13,i)<1) qcx(kobs) = X_NCEP_NLNQC
	      else                                ! rejected by QC
	         sigo(kobs) = undef
		 qcx(kobs) = 2   
	      endif

              ks(kobs)     = iks
              station(kobs)= cdiagbuf(i)
	      weight(kobs) = data(13,i)           ! fractional weight

              if (lobsdiagsave) then              ! place obs sensitivity in sigo slot
                  idia  = ioff
                  obimp = 0.0
                  sigo(kobs) = 0.0
                  call ods_obsdiags ( nlomx(1), tlomx(1), obimp(1), data, idia, i, ninfo, nobs, undef, passed )  
                  if(passed) then
                     sigo(kobs) = obimp(1)
                     qcx (kobs) = 0
                  endif
              endif

!	if(kx(kobs)==304) then
!	  write(*,'(i2,i3,i5,i6,i5,3(1x,f8.3),1x,i3,6(1x,e8.3e1))') qcx(kobs),i,kobs,kx(kobs),kt(kobs),
!     &		lat(kobs),lon(kobs),lev(kobs),time(kobs),obs(kobs),omf(kobs),oma(kobs),sigo(kobs)
!	endif
          end do
	  write(*,*) 'n311 =',n311

      case (' uv')  ! vector wind
      
          do i = 1, nobs

	     kobs = kobs + 1
             if(kobs>mobs) then
                rc = 99
                print *, myname, ' dim violation, var = ', var
                exit
             endif
             iu   = kobs
	     kobs = kobs + 1
             if(kobs>mobs) then
                rc = 99
                print *, myname, ' dim violation, var = ', var
                exit
             endif
             iv   = kobs
		   
             kt(iu)     = ktuu
             kx(iu)     = nint(data(1,i))
             lat(iu)    = data(3,i)
             lon(iu)    = data(4,i)
             lev(iu)    = data(6,i)
             time(iu)   = data(8,i) * 60         ! minutes from ana time
             obs(iu)    = data(17,i)             ! observation
             omf(iu)    = data(18,i)             ! omf used in analysis
             oma(iu)    = undef
             xm(iu)     = data(7,i)              ! elevation
	     
	     qch(iu)    = 0                      ! no info available
             qcx(iu)    = 0

	      
             ks(iu)     = iks
             station(iu)= cdiagbuf(i)
	     weight(iu) = data(13,i)             ! fractional weight

             kt(iv)     = ktvv
             kx(iv)     = nint(data(1,i))
             lat(iv)    = data(3,i)
             lon(iv)    = data(4,i)
             lev(iv)    = data(6,i)
             time(iv)   = data(8,i) * 60         ! minutes from ana time
             obs(iv)    = data(20,i)             ! observation
             omf(iv)    = data(21,i)             ! omf used in analysis
             oma(iv)    = undef
             xm(iv)     = data(7,i)              ! elevation
	     
	     qch(iv)    = 0                      ! no info available
             qcx(iv)    = 0

             ks(iv)     = iks
             station(iv)= cdiagbuf(i)
	     weight(iv) = data(13,i)             ! fractional weight

             if (kx(iu) == 280 .and. nint(data(2,i)) == 562) then
                read(cdiagbuf(i),*,iostat=ios) iwmo
                if (ios == 0 .and. iwmo > 0) then
                   if (mod(iwmo,1000) >=500) then
                      kx(iu) = 299
                      kx(iv) = 299
                   end if
                end if
             end if
             
             if (ninfo >= 12) then
                if (data(12,i) <  0. ) then
                   qcx(iu) = X_PASSIVE
                   qcx(iv) = X_PASSIVE
                endif
             endif
	     if (data(16,i)>small_num) then
                sigo(iu) = 1.0/data(isigo_slot,i) ! sigo
                sigo(iv) = 1.0/data(isigo_slot,i) ! sigo
                 if(data(13,i)<1) then
                    qcx(iu) = X_NCEP_NLNQC
                    qcx(iv) = X_NCEP_NLNQC
                 endif
	     else                                ! rejected by QC
	        sigo(iu) = undef
	        sigo(iv) = undef
	        qcx(iu)  = 2   
	        qcx(iv)  = 2   
	     endif

              if (lobsdiagsave) then
                  idia  = ioff
                  obimp = 0.0
                  sigo(iu) = 0.0
                  sigo(iv) = 0.0
                  call ods_obsdiags ( nlomx, tlomx, obimp, data, idia, i, ninfo, nobs, undef, passed )  
                  if(passed) then
                    sigo(iu) = obimp(1)             ! place obs sensitivity in sigo slot
                    sigo(iv) = obimp(2)             ! place obs sensitivity in sigo slot
                    qcx (iu) = 0                    ! overwrite qc flag
                    qcx (iv) = 0                    ! overwrite qc flag
                  endif
              endif

           end do
	   
      case ('sst')     ! sea surface temperature

          do i = 1, nobs

	      kobs         = kobs + 1
              if(kobs>mobs) then
                 rc = 99
                 print *, myname, ' dim violation, var = ', var
                 exit
              endif
              kt(kobs)     = ktSST
              kx(kobs)     = nint(data(1,i)) 
  
              lat(kobs)    = data(3,i)
              lon(kobs)    = data(4,i)
              lev(kobs)    = undef
              time(kobs)   = data(8,i) * 60       ! minutes from ana time
              obs(kobs)    = data(17,i)           ! observation
              omf(kobs)    = data(18,i)           ! omf used in analysis
              oma(kobs)    = undef
              xm(kobs)     = data(5,i)            ! station elevation
  
	      qch(kobs)    = 0                    ! no info available
              qcx(kobs)    = 0

              if (ninfo >= 12) then
                 if (data(12,i) <  0. ) then
                    qcx(kobs) = X_PASSIVE
                 endif
              endif

	      if (data(16,i)>small_num) then
                 sigo(kobs) = 1.0/data(isigo_slot,i) ! sigo used in analysis
                 if(data(13,i)<1) qcx(kobs) = X_NCEP_NLNQC
	      else                                ! rejected by QC
	         sigo(kobs) = undef
		 qcx(kobs) = 2   
	      endif

              ks(kobs)     = iks
              station(kobs)= cdiagbuf(i)
	      weight(kobs) = data(13,i)           ! fractional weight

              if (lobsdiagsave) then              ! place obs sensitivity in sigo slot
                  idia  = ioff
                  obimp = 0.0
                  sigo(kobs) = 0.0
                  call ods_obsdiags ( nlomx(1), tlomx(1), obimp(1), data, idia, i, ninfo, nobs, undef, passed )  
                  if(passed) then
                     sigo(kobs) = obimp(1)
                     qcx (kobs) = 0
                  endif
              endif

          end do
             
      case (' pw')     ! total column water
                                                                                                                       
          do i = 1, nobs
                                                                                                                       
              kobs         = kobs + 1
              if(kobs>mobs) then
                 rc = 99
                 print *, myname, ' dim violation, var = ', var
                 exit
              endif
              kt(kobs)     = ktTPW
              kx(kobs)     = nint(data(1,i))
                                                                                                                       
              lat(kobs)    = data(3,i)
              lon(kobs)    = data(4,i)
              lev(kobs)    = data(6,i)            ! observation pressure
              time(kobs)   = data(8,i) * 60       ! minutes from ana time
              obs(kobs)    = data(17,i)           ! observation
              omf(kobs)    = data(18,i)           ! omf use in analysis
              oma(kobs)    = undef
              xm(kobs)     = data(9,i)            ! input (prepbufr) qc mark

              qch(kobs)    = 0                    ! no info available
              qcx(kobs)    = 0

              if (ninfo >= 12) then
                 if (data(12,i) <  0. ) then
                    qcx(kobs) = X_PASSIVE
                 endif
              endif

              if (data(16,i)>small_num) then
                 sigo(kobs) = 1.0/data(isigo_slot,i) ! sigo
              else                                ! rejected by QC
                 sigo(kobs) = undef
                 qcx(kobs) = 2
              endif

              ks(kobs)     = iks
              station(kobs)= cdiagbuf(i)
              weight(kobs) = data(13,i)           ! fractional weight

              if (lobsdiagsave) then              ! place obs sensitivity in sigo slot
                  idia  = ioff
                  obimp = 0.0
                  sigo(kobs) = 0.0
                  call ods_obsdiags ( nlomx(1), tlomx(1), obimp(1), data, idia, i, ninfo, nobs, undef, passed )  
                  if(passed) then
                     sigo(kobs) = obimp(1)
                     qcx (kobs) = 0
                  endif
              endif

          end do

      case ('lag')

          do i = 1, nobs

              kobs = kobs + 1  
              if(kobs>mobs) then
                 rc = 99 
                 print *, myname, ' dim violation, var = ', var
                 exit
              endif
              ilat = kobs
              kobs = kobs + 1
              if(kobs>mobs) then
                 rc = 99
                 print *, myname, ' dim violation, var = ', var
                 exit
              endif
              ilon = kobs
          
!             Take lat observation first
!             --------------------------
              kt(ilat)     = ktlat
              kx(ilat)     = nint(data(1,i)) 
              station(ilat)= cdiagbuf(i)
              
              lon(ilat)    = data(4,i)
              lat(ilat)    = data(5,i)
              lev(ilat)    = data(6,i)
              time(ilat)   = data(7,i) * 60 ! convert to minutes from ana time
              obs(ilat)    = data(5,i)      ! obs
              omf(ilat)    = data(17,i)     ! omf
              oma(ilat)    = undef   ! no info available
              xm(ilat)     = undef
              
              qch(ilat)    = 0       ! no info available
              qcx(ilat)    = 0
 
              if (data(8,i) < 0.) then
                 qcx(ilat) = X_PASSIVE
              endif
 
              if ( qcx(ilat) == 0 ) then
                 sigo(ilat) = 1./data(12,i)                ! sigo         
              else                                         ! rejected by QC
                 sigo(ilat) = undef
                 qcx(ilat)  = 2   
              endif
 
              ks(ilat)     = nint(data(3,i))   ! balloon identifier


!             Take lon observation next
!             -------------------------
              kt(ilon)     = ktlon
              kx(ilon)     = nint(data(1,i))
              station(ilon)= cdiagbuf(i)
              
              lon(ilon)    = data(4,i)
              lat(ilon)    = data(5,i)
              lev(ilon)    = data(6,i)
              time(ilon)   = data(7,i) * 60 ! convert to minutes from ana time
              obs(ilon)    = data(4,i)      ! obs
              omf(ilon)    = data(16,i)     ! omf
              oma(ilon)    = undef   ! no info available
              xm(ilon)     = undef
              
              qch(ilon)    = 0       ! no info available
              qcx(ilon)    = 0

              if (data(8,i) < 0.) then
                 qcx(ilon) = X_PASSIVE
              endif


              if ( qcx(ilon) == 0 ) then
                 sigo(ilon) = 1./data(11,i)                ! sigo         
              else                                         ! rejected by QC
                 sigo(ilon) = undef
                 qcx(ilon) = 2
              endif

              ks(ilon)     = nint(data(3,i))   ! balloon identifier

!   OBSERVATION SENSITIVIES NOT HANDLED YET

          end do

      case ('gps')     ! GPS refractivity and bending angle (KX should tell)
                                                                                                                       
          do i = 1, nobs
                                                                                                                       
              kobs         = kobs + 1
              if(kobs>mobs) then
                 rc = 99
                 print *, myname, ' dim violation, var = ', var
                 exit
              endif
              if(nint(data(20,i))==0)then
                 kt(kobs)  = ktGPSr    ! refractivity
                 xm(kobs)  = data(9,i) ! height above model terrain
              else
                 kt(kobs)  = ktGPSb    ! bending angle
                 xm(kobs)  = data(7,i) ! impact height in meters
              endif
              kx(kobs)     = nint(data(1,i))
                                                                                                                       
              lat(kobs)    = data(3,i)
              lon(kobs)    = data(4,i)
              lev(kobs)    = data(6,i)            ! observation pressure
              time(kobs)   = data(8,i) * 60       ! minutes from ana time
              obs(kobs)    = data(17,i)           ! observation
              if (isInf(data(5,i))) then
                 omf(kobs) = undef                ! skip problem GPS OmF
              else
                 omf(kobs) = data(17,i)*data(5,i) ! omf use in analysis
              endif

              oma(kobs)    = undef

              qch(kobs)    = 0                    ! no info available
              qcx(kobs)    = 0

              if (ninfo >= 12) then
                 if (data(12,i) <  0. ) then
                    qcx(kobs) = X_PASSIVE
                 endif
              endif

              if (data(16,i)>small_num) then
                 sigo(kobs) = 1.0/data(16,i)      ! sigo
              else                                ! rejected by QC
                 sigo(kobs) = undef
                 qcx(kobs) = 2
              endif

              ks(kobs)     = iks
              station(kobs)= cdiagbuf(i)
              weight(kobs) = data(13,i)           ! fractional weight

              if (lobsdiagsave) then              ! place obs sensitivity in sigo slot
                  idia  = ioff
                  obimp(1) = 0.0
                  sigo(kobs)=0.0
                  call ods_obsdiags ( nlomx(1), tlomx(1), obimp(1), data, idia, i, ninfo, nobs, undef, passed )  
                  if (passed) then
                    sigo(kobs)=obimp(1)
                    qcx(kobs)=0
                  endif
              endif

          end do

      case default  ! everything else
      
           print *, myname, ' can''t handle var = ', var, nobs, ' observations.'
	   rc = 98
	   
      end select

      
      
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
	  
      end subroutine ods_dconv
