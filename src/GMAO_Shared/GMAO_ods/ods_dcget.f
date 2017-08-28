!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ods_dcget:  get data from diag\_ file and convert to ODS
!
! !INTERFACE:

      subroutine ods_dcget ( fname, nymd, nhms, ods, rc )
 
! !USES

      use m_odsmeta
      use m_ods
      use m_odsxsup, only : getodsmeta 
      use m_odsxsup, only : ncepQCXval 
      use m_Sndx,    only : setSndx
      use m_ods_obsdiags, only : ods_obsdiags_getparam
      
      Implicit NONE

! !INPUT PARAMETERS:
 
      character(len=*), intent(in)   :: fname   ! GSI diag_ file name
      integer, intent(in)            :: nymd    ! year-month-day, e.g., 19990701
      integer, intent(in)            :: nhms    ! hour-min-sec,   e.g., 120000

! !OUTPUT PARAMETERS:
 
      type(ods_vect), intent(inout)  :: ods     ! ODS vector

      integer, intent(out)           :: rc      ! Error return code:

! !DESCRIPTION: get data from GSI diag\_ files and convert to ODS
!
! !REVISION HISTORY:
!       14Apr2004 - R. Todling - Initial code.
!       06Dec2004 - D. Dee     - removed kt/kx/qc mappings
!       14Dec2004 - D. Dee     - handle satellite data
!       07Feb2005 - D. Dee     - set ks to actual station id for raobs
!       02May2005 - D. Dee     - fix for sbuv/2
!       20May2005 - D. Dee     - handle passive radiance data
!       03Nov2005 - Sienkiewicz - pass idiagbuf (incl. sfctype) to ods_drad
!       07Nov2005 - Sienkiewicz - handle rain rate obs
!       07Aug2006 - Sienkiewicz - don't use station ID <=0 for 'ks'
!       14Dec2006 - Todling     - update to GSI-2006_09
!       08Apr2008 - Sienkiewicz - changes for OMI and MLSoz ozone
!       31Jan2010 - Todling     - re-defined how ks is determined for 
!                                 Raob and Lagragian ballon obs
!       14Mar2010 - Todling     - update sbuv ozone read to Mar2010 version
!       30Jul2010 - Sienkiewicz - fix sbuv ozone (Mar2010 version)
!        3Nov2010 - Sienkiewicz - fix o3lev ozone (2010 version)
!       17Feb2012 - Todling     - update to Sep2011 GSI
!       19Jul2012 - Sienkiewicz - add 'tirosn' to getsatid_
!       19Oct2012 - Sienkiewicz - add npp and meteosat to getsatid_
!                                 extend for metop-b,c
!	18Mar2012 - Guo         - changes for OMIeff ozone.
!       02May2013 - K Wargan    - Add TOMS ozone     
!        1Jul2013 - Sienkiewicz - add precip types to getsatid_
!       23Sep2013 - Guo         - replaced the use of getsatid_() with
!                                 pcp_getisat_(), rad_getisat_(), and
!                                 ozone_getisat_().
!       28Jan2014 - Todling     - read sensitivity slot indicator (ioff) from header
!       02Apr2014 - Todling     - gather ob variance error for possible output option
!
!EOP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      character(len=*), parameter :: myname_ = 'ods_dcget'
      logical,          parameter :: verbose = .true.

      integer  ODS_Handle
      integer  ios, ier, lu, i, k, id, iks, ii
      integer  mobs,ioff
      integer  idiag,angord,iversion_radiag,inewpc
      integer  nymdf, nhmsf                  ! date and time found in file
      logical  conv
      logical  lobsdiagsave
      logical  lextra
      integer  nobs_ods, nqcmax
      character(len=10) satype
      character(len=10) pcptype             ! shorter string for precip data

      character(len=3)                         :: var
      character(len=5)                         :: statid
      integer                                  :: statks
      integer                                  :: idate
      integer                                  :: nchar,ninfo,nobs,nlevs,mype
      integer                                  :: isat,nchanl,npred,iint,ireal,ipchan,iextra,jextra,jiter,ndiag
      real(4)                                  :: freq4,pol4,wave4,tlap4
      integer                                  :: miter
      integer,     allocatable, dimension(:)   :: nuchan
      integer,     allocatable, dimension(:)   :: ich
      integer,     allocatable, dimension(:)   :: iuse_rad
      real(4),     allocatable, dimension(:)   :: weight
      character(8),allocatable, dimension(:)   :: cdiagbuf
      character(8),allocatable, dimension(:)   :: station
      integer,     allocatable, dimension(:)   :: idiagbuf
      integer,     allocatable, dimension(:,:) :: idiagbuf2
      real(4),     allocatable, dimension(:)   :: varch4
      real(4),     allocatable, dimension(:)   :: diagbuf
      real(4),     allocatable, dimension(:,:) :: diagbuf2
      real(4),     allocatable, dimension(:,:) :: diagbufex
      real(4),     allocatable, dimension(:,:,:) :: diagbuf3
      integer,     allocatable, dimension(:)   :: iouse
      real(4),     allocatable, dimension(:)   :: pobs,gross,tnoise

      character(len=20)  isis                ! sensor/instrument/satellite id  ex.amsua_n15
      character(len=10)  dplat               ! sattelite platform id

      
      rc = 0 ! all is well to start with
      
!     Inquire whether this is diag file w/ history of omx, osens, ...
!     --------------------------------------------------------------- 
      call ods_obsdiags_getparam ( 'lobsdiagsave', lobsdiagsave )
      call ods_obsdiags_getparam ( 'miter', miter )

!     Scan the file for date, number of obs, and data type(s)
!     -------------------------------------------------------
      call ods_dcscan ( verbose, fname, idate, mobs, conv, satype, ier )
      if (ier/=0) then
         rc = 1
         print *, myname_, ': ods_dcscan() error'
         return
      end if

      nymdf = idate / 100
      nhmsf = 10000 * (idate - 100 * nymd)
      if (nymdf/=nymd .or. nhmsf/=nhms) then
         rc = 1
         print*, 'Date/Time requested not found in file'
         print*, ' nymd , nhms =',nymd ,nhms  
         print*, ' nymdf, nhmsf=',nymdf,nhmsf
         return 
      end if

!     Create ODS structure and define metadata
!     ----------------------------------------
      nqcmax = 17
      call ods_init ( ods, mobs, ier, nqc=nqcmax )
      if (ier/=0) then
         rc = 2
         print *, myname_, ' ods_init() error'
         return
      end if

      call getodsmeta( ods )
      call ncepQCXval( ods )    
     
!     Open the data file:
!     ------------------
      lu = ODS_Handle(fname,ier)
      if (ier/=0) then
         rc = 99
         print*, myname_, ': cannot get file handle'
         return
      end if
      open(lu,file=fname,form='unformatted')

      if (conv) then      
      
!        Conventional data:
!        -----------------
         allocate ( station(mobs), weight(mobs), stat=ier )
         if (ier/=0) then
            rc = 3
            print *, myname_, ': cannot alloc()'
            return
         end if
     
         read(lu)      ! very first record is date

         ios      = 0
         nobs_ods = 0
         iks      = -999
         do while (ios==0)

            read (lu,iostat=ios) var,nchar,ninfo,nobs,mype,ioff   ! header

            if (ios==0) then

               if (verbose) print *, 'will read for ', var, ' nobs = ', nobs, ' ninfo = ', ninfo,
     .                               ' ioff= ', ioff
               allocate ( cdiagbuf(nobs), diagbuf2(ninfo,nobs) )

               read(lu) cdiagbuf, diagbuf2

               call ods_dconv ( var, cdiagbuf, diagbuf2, ninfo, nobs, 
     .                          ods, station, weight, mobs, nobs_ods, iks, ioff, ier )
               if (ier/=0) then
                  rc = 4
                  print*, myname_,': ods_dconv() error'
                  return
               end if

               deallocate ( cdiagbuf, diagbuf2 )
         
            end if ! < ios >
         end do

!        Note: For testing purposes, we can store 'weight' in one of the ods attributes
!        ----

!        Randomly assign ks by station id
!        --------------------------------
         call setsndx (ods%data%ks(1:nobs_ods),ods%data%kx(1:nobs_ods),station(1:nobs_ods))

!        Set ks to actual station id for radiosondes
!        -------------------------------------------
         do i = 1, nobs_ods
            if (ods%data%kx(i)==120 .OR. ods%data%kx(i)==220  .or.  ! radiosondes
     .          ods%data%kx(i)==901 .OR. ods%data%kx(i)==902) then  ! lagragian ballon data

!              Normally the station id is a 5-digit number ...
!              -------------------------------------------
               statid = station(i)(1:5)
               read(statid,'(i5)',err=999) statks
               if (statks > 0 ) then
                  ods%data%ks(i) = statks
               else         !  if invalid numeric id, treat like non-numeric id
                            !  RT: code will never be here!
                  do ii =1,len(station(i))
                     ods%data%ks(i) = ods%data%ks(i)+ii*iachar(station(i)(ii:ii)) !  arbitrary scheme to use a number better than 99999
                  enddo
                   print*, 'ods_dcget: YES CODE PASSES HERE SOMETIMES!' 
               endif
               cycle

!              ... but if not, make sure this ks is distinct from any other station id
!              -----------------------------------------------------------------------
  999          continue

               do ii =1,len(station(i))
                  ods%data%ks(i) = ods%data%ks(i)+ii*iachar(station(i)(ii:ii)) !  arbitrary scheme to use a number better than 99999
               enddo
               if (verbose) print *,'Non-numeric station id ', statid, ' for kx = ', ods%data%kx(i), 'assigned ks =', ods%data%ks(i)  

            end if
         end do

!        Clean up
!        --------
         deallocate ( station, weight )

      else if ( any(satype==pcpt) ) then

!        precip data
!        -------------
         ios      = 0
         nobs_ods = 0
         iks      = 0

         read(lu,iostat=ios) isis,dplat,pcptype,jiter,idate,iint,ireal,iextra ! header

         if (ios==0) then

            call getsatid_(isat)
            !call pcp_getisat_(pcptype,dplat,isis,isat)
            if(verbose) print *, 'will read for ', pcptype, ' sat = ', dplat,' isat = ',isat
            allocate( idiagbuf(iint), diagbuf(ireal), stat=ier)

      	    if (ier/=0) then
               rc = 3
               print *, myname_, ': cannot alloc()'
               return
            end if

            do while (ios==0)    ! read data records until eof
               read(lu,iostat=ios) mype,idiagbuf, diagbuf
               if (ios==0) then
                  iks = iks + 1
                  call ods_dpcp( satype, isat, idiagbuf, diagbuf, iint, ireal,
     .                           ods, mobs, nobs_ods, iks, ier )
                  if (ier/=0) then
                     rc = 4
                     print*, myname_, ': ods_dpcp() error'
                     return
                  end if
               end if
            end do

            deallocate( idiagbuf, diagbuf )
         end if

      else if (trim(satype)=='o3lev_old') then

!        o3 level data
!        -------------
         ios      = 0
         nobs_ods = 0
         iks      = 0

         read(lu,iostat=ios) isis,dplat,satype,jiter,nlevs,idate,ninfo
         if (ios/=0) then
            rc = 4
            print*, myname_, ': error reading (1) o3lev data'
            return
         end if
         ios = 0

         do while (ios == 0) 

            read(lu,iostat=ios) nobs
            if (ios/= 0) exit          !end of file

            if (nobs > 0) then
               allocate( diagbuf2(ninfo,nobs), stat=ier)

               if (ier/=0) then
                  rc = 3
                  print *, myname_, 'cannot alloc(), satype = "',trim(satype),'"'
                  print *, myname_, '                 ninfo = ',ninfo
                  print *, myname_, '                 nlevs = ',nlevs
                  print *, myname_, '                  nobs = ',nobs
                  print *, myname_, '                   ier = ',ier
                  return
               end if

               read(lu,iostat=ios) diagbuf2
	       if(.true.) then
		  print*, myname_,': cannot call ods_do3lev() for now'
		  rc=9
		  return
	       endif
               if (ios==0) then
                  iks = iks + 1
                  call ods_do3lev( satype, diagbuf2, ninfo,nobs,
     .                    ods, mobs, nobs_ods, iks, ier )
                  if (ier/=0) then
                     rc = 4
                     print*, myname_, ' ods_do3lev() error'
                     return
                  end if
               end if

               deallocate( diagbuf2 )
               
            else     ! if no observations
               read(lu,iostat=ios)
            end if


         end do

	 
      else if (trim(satype)/='sbuv2' .and. trim(satype)/='omi' .and.
     &         trim(satype)/='mls'   .and. trim(satype)/='mls20' .and.
     &         trim(satype)/='mls22' .and. trim(satype)/='mls30' .and.
     &         trim(satype)/='mls55' .and. trim(satype)/='tomseff' .and.
     &	       trim(satype)/='omieff'.and. trim(satype)/='o3lev' .and.
     &         trim(satype)/='gome'  ) then
      
!        Radiance data:
!        -------------
         ios      = 0
         nobs_ods = 0
         iks      = 0

         read (lu,iostat=ios) isis,dplat,satype,jiter,nchanl,npred,idate,ireal,ipchan,iextra,jextra,
     &                        idiag,angord,iversion_radiag,inewpc,ioff
	 
         lextra=iextra>0
         if (ios==0) then

	    call getsatid_(isat)
	    !call rad_getisat_(satype,dplat,isis,isat)
            if(verbose) print *, 'will read *'//trim(dplat)//'* for ', satype, ' isat = ', isat, ' nchanl = ', nchanl

            iint = 0
            ndiag = ipchan+npred+2
            if (lobsdiagsave) ndiag=ndiag+4*miter+1
            allocate ( iuse_rad(nchanl), idiagbuf(iint), diagbuf(ireal), diagbuf2(ndiag,nchanl),
     &                 varch4(nchanl), nuchan(nchanl), ich(nchanl), stat=ier )	       
      	    if (ier/=0) then
               rc = 3
               print *, myname_, 'cannot alloc()'
               return
            end if
            if (lextra) allocate(diagbufex(iextra,jextra))
	    
	    do i = 1, nchanl
               read (lu) freq4,pol4,wave4,varch4(i),tlap4,iuse_rad(i),nuchan(i),ich(i)
	    end do

	    do while (ios==0)        ! read data records until eof
	    
               if (lextra) then
	          read (lu,iostat=ios) diagbuf,diagbuf2 
               else
	          read (lu,iostat=ios) diagbuf,diagbuf2,diagbufex
               endif
	       
	       if (ios==0) then
		  iks = iks + 1
                  idiagbuf = 0  ! dummy array in this case

                  call ods_drad ( satype, isat, iuse_rad, 
     .                            idiagbuf, diagbuf, diagbuf2, varch4, nuchan, ndiag, iint,
     .                            ireal, ipchan, npred, nchanl,
     .   	                  ods, mobs, nobs_ods, iks, ioff, ier )
                  if (ier/=0) then
        	     rc = 4
        	     print*, myname_,': ods_drad() error'
        	     return
                  end if
	       end if
	    end do
	    
	    deallocate ( iuse_rad, idiagbuf, diagbuf, diagbuf2, varch4, ich )
            if (lextra) deallocate(diagbufex)
	    
	 end if
	 
      else
      
!        Ozone data (SBUV2,OMI, OMIeff (and whatever survived upto this point):
!        ----------------------
         ios      = 0
         nobs_ods = 0
         iks      = 0

         read (lu,iostat=ios) isis,dplat,satype,jiter,nlevs,idate,iint,ireal,iextra ! header
	 
         if (ios/=0) then
            rc = 4
            print*, myname_, ': error reading (1) sbuv2/omi/omieff/tomseff/o3lev/mlsXX data, satype = "',trim(satype),'"'
            return	    
	 end if

         call ozone_getisat_(satype,dplat,isis,isat)
         if(verbose) print *, 'will read *'//trim(dplat)//'* for ', satype, ' isat = ', isat, ' nchanl = ', nchanl
	 	 
         ios = 0	    
         ndiag = 6
         if (lobsdiagsave) ndiag=ndiag+4*miter+1
         allocate ( pobs(nlevs), gross(nlevs), tnoise(nlevs), iouse(nlevs), stat=ier )
         if (ier/=0) then
            rc = 3
            print *, myname_, ':  cannot alloc(pobs), ier= ', ier
            return
         end if
	 	! set default values for "o3lev" case.
	 pobs  (:)=-1.
	 gross (:)= 0.
	 tnoise(:)= 0.
	 iouse (:)=+1.
         
	 select case(satype)
	 case('o3lev','mls','mls20','mls22','mls30','mls55')
	 case default
           read (lu,iostat=ios) pobs,gross,tnoise,iouse
           if (ios/=0) then
              rc = 4
              print*, myname_, ': error reading (2) sbuv2/omi/omieff/tomseff/o3lev/mlsXX data, satype = "',trim(satype),'"'
              return	    
           end if
	 end select

         ios = 0

	 do             ! read until end-of-file
	    
	    read (lu,iostat=ios) nobs
            if (ios/=0) exit   ! end-of-file

	    allocate ( idiagbuf2(iint,nobs), diagbuf2(ireal,nobs), diagbuf3(ndiag,nlevs,nobs), 
     .   	       stat=ier )
      	    if (ier/=0) then
               rc = 3
               print *, myname_, ':  cannot alloc(diagbug2), ier= ', ier
               return
            end if
	    read (lu,iostat=ios) idiagbuf2,diagbuf2,diagbuf3
            if (ios/=0) then
               rc = 4
               print*, myname_, ': error reading (3) sbuv2/omi/omieff/tomseff/o3lev/mlsXX data, satype = "',trim(satype),'"'
               return	    
	    end if
	
	    do i = 1, nobs
	       iks = iks + 1
	       	
               call ods_dsbuv ( isat, satype, diagbuf2(:,i), diagbuf3(:,:,i), 
     .   			pobs, tnoise, iouse,
     .   			ireal, nlevs, ndiag,
     .   			ods, mobs, nobs_ods, iks, ier )
               if (ier/=0) then
         	   rc = 4
         	   print*, myname_,': ods_dsbuv() error, with satype = "',trim(satype),'"'
         	   return
               end if
	    end do
	    	
	    deallocate ( idiagbuf2, diagbuf2, diagbuf3 )
       	 
	 end do

         deallocate (  pobs, gross, tnoise, iouse )
	 	 	       
      end if
         
      close (lu)
      
      if (nobs_ods==0) then
         rc = 10
         print*, myname_, ' no data'
         return
      end if
		 	 
      ods%data%nobs = nobs_ods      

!     Fix longitudes so they become ODS-compatible
!     -------------------------------------------
      where ( ods%data%lon(1:nobs_ods) > 180 ) 
              ods%data%lon(1:nobs_ods) = ods%data%lon(1:nobs_ods) - 360.
      end where  	  

      print *, 'nobs to write to ODS file: ', nobs_ods
      id = count(ods%data%kt(1:nobs_ods)<=0.or.ods%data%kt(1:nobs_ods)>ktmax)
      if (id>0) then
          rc = 5
          print *, myname_,' kt out of range'
          return
      end if

      return

      contains
         subroutine getsatid_(myisat)
! RT: some heck that needs more work
! RT: this needs serious attention as it is becoming a huge heck now (3/30/09)
            implicit none
            integer, intent(out) :: myisat
            integer ios

            myisat     = 0  ! take fixed sat index as in idsats

!	    select case( trim(ladjust(dplat)) )
!	    case ('aura')
!              myisat = 999
!	      return
!	    end select

! first handle precip types
            i = index('pcp',isis(1:3))
            if (i>0) then
               i = index('trmm',dplat(1:4))
               if (i>0) then
                  if (dplat(6:8) == 'lnd') then
                     myisat = 1
                  else if (dplat(6:8) == 'ocn') then
                     myisat = 2
                  else
                     myisat = 0
                  end if
                  return
               else
                  read(dplat(5:6),'(i2)',iostat=ios)myisat
                  return
               end if
            end if

! Need to distinguish between AMSUA from AQUA and METOP for example 
! (NOAA sats are already distinguished)
            i = index('metop',dplat(1:5))
            if(i>0)then
               myisat = 25 + iachar(dplat(7:7)) - iachar('a')  
               return
            endif
            i = index('tiros',dplat(1:5))
            if(i>0)then
               myisat = 5
               return
            endif
!           i = index('aqua',dplat(1:4))
!           if(i>0)then
!              return
!           endif
! the "word" fgnm stands for the platforms dmsp/goes/noaa/meteosat this needs generalization
! e.g. dmsp -> f15   goes -> g12   noaa -> n18  meteosat -> m09
            i = index('fgnm',dplat(1:1))
            if(i>0)then; read(dplat(2:3),'(i2)',iostat=ios)myisat; endif

         end subroutine getsatid_

      subroutine pcp_getisat_(dtype,dplat,dsis,myisat)
        implicit none
	character(len=*),intent(in):: dtype,dplat,dsis
	integer,intent(out):: myisat

	myisat=0	! for a default isat value.

	if(dsis(1:3)=='pcp'.or.dtype(1:3)=='pcp') then
	   select case(dplat)
	   case('trmm','dmsp')
              myisat = 0
	   case('trmm_lnd')
              myisat = 1
	   case('trmm_ocn')
              myisat = 2
	   case default
	      ! case("trmm[0-9][0-9]")
	      if(dplat(1:4)=='trmm' .and. verify(dplat(5:6),'0123456789')==0) then
                 read(dplat(5:6),'(i2)') myisat
	      endif
	   end select
	endif
      end subroutine pcp_getisat_

      subroutine rad_getisat_(dtype,dplat,dsis,myisat)
        implicit none
	character(len=*),intent(in):: dtype,dplat,dsis
	integer,intent(out):: myisat

	myisat=0	! for a default isat value.

	if(len_trim(dplat)==3 .and.
     .     verify(dplat(1:1),'fgnm')==0 .and.
     .     verify(dplat(2:3),'0123456789')==0) then
! the "word" fgnm stands for the platforms dmsp/goes/noaa/meteosat this needs generalization
! e.g. dmsp -> f15   goes -> g12   noaa -> n18  meteosat -> m09
           read(dplat(2:3),'(i2)') myisat

	else
! Need to distinguish between AMSUA from AQUA and METOP for example 
! (NOAA sats are already distinguished)
           select case(dplat)
	   case('tirosn')
	      myisat = 5
	   case('metop-a')
	      myisat = 25
	   case('metop-b')
	      myisat = 26
	   case default
	      if(dplat(1:6)=='metop-') myisat = 25 + iachar(dplat(7:7)) - iachar('a')  
	   end select
	endif
      end subroutine rad_getisat_

      subroutine ozone_getisat_(dtype,dplat,dsis,myisat)
        implicit none
	character(len=*),intent(in):: dtype,dplat,dsis
	integer,intent(out):: myisat

	myisat=0	! for a default isat value.

        if(len_trim(dplat)==3 .and.
     &     verify(dplat(1:1),'fgnm')==0 .and.
     &     verify(dplat(2:3),'0123456789')==0) then
! the "word" fgnm stands for the platforms dmsp/goes/noaa/meteosat this needs generalization
! e.g. dmsp -> f15   goes -> g12   noaa -> n18  meteosat -> m09
           read(dplat(2:3),'(i2)') myisat

        else
	   select case(dplat)
	   case('aura')
	      myisat = 0
	   case('nim07')
	      myisat = 1
	   case('ep')
	      myisat = 2
	   case('metop-a')
	      myisat = 25
	   case('metop-b')
	      myisat = 26
	   case default
	      if(dplat(1:6)=='metop-') myisat = 25 + iachar(dplat(7:7)) - iachar('a')  
	   end select
        endif
      end subroutine ozone_getisat_
 
      end subroutine ods_dcget
