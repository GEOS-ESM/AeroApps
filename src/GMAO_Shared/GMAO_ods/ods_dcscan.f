!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!BOP
! !ROUTINE: ods_dcscan: echoes information on diag\_ files
!
! !DESCRIPTION:
!
! !INTERFACE:
!
      subroutine ods_dcscan ( verbose, fname, idate, mobs, conv, satype, rc )

      use m_odsmeta

      implicit NONE

! !INPUT PARAMETERS:
 
      logical,          intent(in)  :: verbose
      character(len=*), intent(in)  :: fname

! !OUTPUT PARAMETERS:

      integer, intent(out)           :: mobs
      integer, intent(out)           :: idate
      logical, intent(out)           :: conv   ! .true. for conventional data
      character(len=10), intent(out) :: satype ! sensor id if conv==.false.
      integer, intent(out)           :: rc

! !DESCRIPTION: echoes information on GSI diag\_ files
!
! !REVISION HISTORY:
!	15Apr2004 - R. Todling - Initial code.
!       14Dec2004 - D. Dee     - Handle sat files (changed interface)  
!	24Feb2005 - Todling    - Compilation errors/vars double declare
!       03Mar2005 - Dee        - Updated to support additional sensors
!       29Apr2005 - Dee        - Fix for sbuv2
!       08Nov2005 - Sienkiewicz - add pcp (rain rate) data
!       17Oct2006 - Sienkiewicz - add SSU
!       14Dec2006 - Todling     - update to handle GSI 2006_09
!        9Apr2008   Meta          handle OMI with SBUV, add 'mlsoz'
!       05Apr2009   Todling       add GPS (refractivity and bending angle)
!       14Mar2010   Todling       update sbuv to Mar2010 version of GSI
!       26Jul2010   Todling       add tropical cyclone knob (tcp)
!       30Jul2010 - Sienkiewicz - fix sbuv ozone (Mar2010 version)
!        3Nov2010 - Sienkiewicz - fix 'o3lev' ozone (2010 version)
!       02May2013 - K Wargan    - Add TOMS ozone      
!       28Jan2014 - Todling     - upd hear of conv and rad
!
!EOP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      character(len=*), parameter :: myname = 'ods_dcscan'

      character(len=3)  var
      integer           i, icount, lu, is(1), ios, maxnobs, ier
      integer           nchar,ninfo,nobs,nlevs,mype
      integer           jiter,isat,nchanl,npred,iint,ireal,ipchan,iextra,jextra
      integer           idiag,angord,iversion_radiag,inewpc,ioff
      integer     ODS_Handle
      external    ODS_Handle
      character(len=10)  pcptype             ! shorter string for precip data
      character(len=20)  isis                ! sensor/instrument/satellite id  ex.amsua_n15
      character(len=10)  dplat               ! sattelite platform id

      integer, parameter :: nvars = 12
      character(len=*), parameter :: vars(nvars)=(/
     .                 ' dw', ' ps', ' pw', ' rw', 
     .                 '  q', 'spd', '  t', ' uv', 
     .                 'sst', 'gps', 'lag', 'tcp' /)
      integer ivars(nvars)
      logical lvars(nvars)
      character(len=16) :: junk
      
      integer maxloc

! Initialize parameters
! ---------------------
      rc      = 0
      ios     = 0
      icount  = 0
      maxnobs = 0
      mobs    = 0
      lvars   = .false.
      conv    = .false.
      do i=1, nvars; ivars(i) = i; end do

      lu = ODS_Handle (fname, ier)  ! just get a file unit (avoid depend on mpeu)
      if (ier/=0) then
         rc = 1
         print*, myname, ': Cannot obtain a handle for file ', fname
         return
      end if
      open(lu,file=fname,form='unformatted')
      
!     Read first record
!     -----------------      
      read (lu,iostat=ios) idate
      if (ios/=0) then
         rc = 2
         print*, myname, ': Cannot read first record on file ', fname
	 close(lu)
         return
      end if
	 
!     Try conventional data first
!     ---------------------------
      do while (ios==0)
                                        ! try to read the header record
         read (lu,iostat=ios) var,nchar,ninfo,nobs,mype,ioff  
	                               
	 if (ios==0) then              ! we appear to have conventional data
            read (lu,iostat=ios)       ! try read a data record
	    is = maxloc(ivars,vars==var)
            if (is(1)>0 .and. is(1)<=nvars) then
	       lvars(is(1)) = .true.	 
               maxnobs = max(maxnobs,nobs)
               mobs    = mobs   + nobs
	       if (var==vars(8)) mobs = mobs + nobs   ! wind vectors
	       if (var==vars(11)) mobs = mobs + nobs   ! lagrangean data
               icount  = icount + 1
	    end if
	 end if
      end do
      
      if (icount>0) conv = .true.   ! successfully read conventional data
	 
      if (icount==0) then
      
!     Try radiances next
!     ------------------ 
         rewind (lu) 
	 ios = 0     
	                            ! try to read the first header record
         read (lu,iostat=ios) isis,dplat,satype,jiter,nchanl,npred,idate,ireal,ipchan,iextra,jextra,
     &                        idiag,angord,iversion_radiag,inewpc,ioff
	 
	 if (ios==0) then 	

            print*, myname, ': Trying satype as one of sats(:) ...'
            if (.not. any(satype==sats)) then
               print*, myname, ': Unrecognized sensor on diag_sat file: satype = ', trim(satype)
	       goto 301
	    endif
	    do i = 1, nchanl          ! try to read additional header records
               read (lu,iostat=ios) 
	       if (ios/=0) then
                  rc = 3
                  print*, myname, ': Cannot read header record from diag_sat file'
                  return
               endif	         
	    end do
	    do while (ios==0)        ! read data records until eof
 	       read (lu,iostat=ios)
 	       if (ios==0) icount = icount + 1
	    end do
	    maxnobs = nchanl
	    mobs = icount*nchanl	 
	 end if
      end if
 301  continue

      if (icount==0) then
      
!     Try ozone next
!     --------------                         
         rewind (lu) 
	 ios = 0
! try to read the header record
         read (lu,iostat=ios) isis,dplat,satype,jiter,nlevs,idate,iint,ireal,iextra ! header
	 
         if (ios==0) then     
	                                                         
            print*, myname, ': Trying satype as one of sbuv2/omi/omieff/tomseff ...'
            if ( .not. trim(satype)=='sbuv2'  .and. 
     &           .not. trim(satype)=='omieff' .and.
     &           .not. trim(satype)=='omi'    .and.
     &           .not. trim(satype)=='tomseff' .and.
     &           .not. trim(satype)=='gome' ) then
               print*, myname, ': Unrecognized sensor on diag_ file: satype = ', trim(satype)
	       goto 302
            endif	    
            read (lu,iostat=ios) ! header 2
            if (ios/=0) then
               rc = 3
               print*, myname, ': Cannot read second header record from diag_sbuv file'
               return
            endif
            do while (ios==0)   ! read nobs and data records until eof
               read (lu,iostat=ios) nobs
               if (ios==0) then
                  read(lu,iostat=ios) ! data
                  if (ios/=0) then
                     rc = 3
                     print*, myname,
     &                    ': Problem reading data record from sbuv/omi/omieff/tomseff file, satype ="',trim(satype),'"'
                     return
                  endif
                  icount = icount + 1
                  maxnobs = max(maxnobs,nlevs*nobs)
                  mobs = mobs + nlevs*nobs
               endif
            end do
         end if
      end if
 302  continue

!     Try rain rate (pcp) next
!     ------------------------

      if (icount==0) then
         rewind(lu)
         ios = 0
         do while (ios==0)    ! read until end-of-file
                              ! try to read the first header record
            read(lu,iostat=ios) isis,dplat,pcptype,jiter,idate,iint,ireal,iextra ! header

            if (ios == 0) then
               if ( .not. any(pcptype==pcpt)) exit  ! if no match, leave loop

               satype = pcptype

	       do while (ios==0)        ! read data records until eof
 	          read (lu,iostat=ios)
 	          if (ios==0) icount = icount + 1
	       end do

               maxnobs = 1
               mobs = icount

            end if
         end do
      endif
 303  continue

!     Try ozone level (MLS) data next
!     ------------------------

      if (icount==0) then
         rewind(lu)
         ios = 0
! Try to read the header record
         read(lu,iostat=ios) isis,dplat,satype,jiter,nlevs,idate,ninfo
         if (ios==0) then
            print*, myname, ': Trying satype as one of mlsXX/o3lev types ...'
	    select case(satype)
	    case('o3lev','mls','mls20','mls22','mls30','mls55')
	    case default
               print *, myname, 
     &              ': Unrecognized sensor on diag_ file: satype = ', trim(satype)
	       goto 304
            end select
         end if

         do while (ios==0)      ! read nobs and data records until eof
            read(lu,iostat=ios)  nobs
            if (ios==0) then
               read(lu,iostat=ios) !  data
               if (ios/=0) then
                  rc = 3
                  print *, myname,
     &                 ': Cannot read data record from diag_'//trim(satype)//' file'
                  return
               endif
               icount = icount + 1
               maxnobs = max(maxnobs,nobs)
               mobs = mobs + nobs
            endif
         end do
      endif
 304  continue

      close (lu)
      
      if (icount==0) then
      
!     Give up
!     -------
         rc = 3
         print*, myname, ': All attempt failed'
         print*, myname, ': Unrecognized header on file ', fname
         return
	    
      end if

! Summarize what's found so far
! -----------------------------
      if ( verbose ) then
         print *, 'Date on file: ', idate
         if (conv) then
	    if (count(lvars)>0) then
               print *, ' variables found on file: '
               do i = 1, nvars
                  if(lvars(i)) print *, trim(vars(i))
               end do
            endif
	 else
	    print *, ' sensor found on file: ', satype
	 endif   
         print *, ' total number of data records found: ', icount
         print *, ' max number of obs per data record:  ', maxnobs
         print *, ' total number of observations:       ', mobs
      endif
      
      return

      end subroutine ods_dcscan
