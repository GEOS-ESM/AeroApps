!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  ods_scan.x - This stand-alone program scans an ODS file.

      program ods_scan

!
! !DESCRIPTION:  The default behavior of this program will open an ODS 
!                and print the number of observations for each synoptic 
!                time.  An assumption is made that the ODS file contains 
!                data from a single day.
!
! !INTERFACE:
!
!      scan_ods.x [ -t synoptic_hour ] [ -l ] filename
!
!
! !OPTIONS:
!
!          -t synoptic_hour  Print the number of obs for the 
!                            specified hour
!
!          -l                Print the last synoptic time for
!                            which the number of obs is not 0.
!                            (Valid times: 00, 06, 12, or 18).
!
!          -x                Print the number of observations 
!                            for each KX.
!
!
! !RETURNS:
!       Unix return code:
!            0 - no problems
!           99 - some error occured opening or reading file
!
! !NOTES:
!     This code uses three SGI extensions to Fortran 90:
!         subroutine exit()
!         subroutine getarg()
!         integer function iargc()
!
! !REVISION HISTORY:
!     03/10/2000 Lucchesi - Original version.
!
!EOP
!-------------------------------------------------------------------------

      implicit none

      character*180 filename
      character*180, allocatable :: input_buffer(:)
      character*180, allocatable :: kx_names(:)
      integer, allocatable :: kx(:), kxtable(:)
      integer jbeg, jend, jday, ODS_caldat, ODS_NGet
      integer ymd, nobs, syn_hour(4)
      integer user_hour, nkx
      integer ncid, ier, i, iargc, num_options, j, k
      logical DEFAULT
      logical T_FLAG, L_FLAG, X_FLAG
      data syn_hour /0,6,12,18/
 
      DEFAULT=.true.
      T_FLAG=.false.
      L_FLAG=.false.
      X_FLAG=.false.

!***********************************************************
!  Process input arguments
!***********************************************************
 
      ! Check for impossible number of command line arguments.

      num_options = iargc()
      if (num_options .LT. 1 .OR. num_options .GT. 3) call usage ()
      
 
      ! Allocate memory for arguments
 
      allocate ( input_buffer ( num_options) )

      ! Loop through arguments

      i=1
      do while ( i .LE. num_options )
        call getarg (i, input_buffer (i))

        if ( trim(input_buffer(i)) .EQ. '-t' ) then  
          if ( num_options .LT. 3 ) call usage ()
          if ( L_FLAG ) call usage()

          ! Argument after -t will be the synoptic hour.
          ! Entering a non interger will crash the program.

          call getarg (i+1, input_buffer(i+1))
          read (input_buffer(i+1), 99) user_hour
99        format (I2)
          T_FLAG = .true.
          DEFAULT = .false.
          i = i + 1

        else if ( trim(input_buffer(i)) .EQ. '-l' ) then
          if ( num_options .NE. 2 ) call usage ()
          if ( T_FLAG ) call usage ()
          L_FLAG = .true.
          DEFAULT = .false.

        else if ( trim(input_buffer(i)) .EQ. '-x' ) then
          if ( num_options .NE. 2 ) call usage ()
          if ( X_FLAG ) call usage ()
          X_FLAG = .true.
          DEFAULT = .false.

        else
          filename = input_buffer(i)
        endif
        i = i + 1
      enddo

!***********************************************************
!  Open and read ODS file 
!***********************************************************

      call ODS_Open ( ncid, trim(filename), 'r', ier ) ! open the file
      if ( ier .ne. 0 ) then
         print *, 'could not open ods file ',trim(filename)
         call exit (99)
      endif

      call ODS_IGet ( ncid, 'syn_beg:first_julian_day',  jbeg, ier )
      call ODS_IGet ( ncid, 'syn_beg:latest_julian_day', jend, ier )
      call ODS_IGet ( ncid, 'nkx', nkx, ier )

      if (jbeg .NE. jend) then
        print *, 'file has more than one day ',trim(filename)
        call exit (99)
      endif

      ymd = ODS_CalDat ( jbeg )
 
      ! Depending on arguments specified, take appropriate action.

      if (DEFAULT) then
        print *, ' '
        print *, 'Scanning ',trim(filename)
        print *, ' '
        print *, 'Date=',ymd
        print *, ' '
        print *, 'Synoptic Hour    Number of Obs'
        print *, '-------------    -------------'
        do i=1,4
           nobs = ODS_NGet ( ncid, jbeg, syn_hour(i), ier )
           if ( ier .ne. 0 ) then
              print *, 'could not read ods file ',trim(filename)
              call exit (99)
           end if
           write (6,100) syn_hour(i), nobs
100        format ('       ',I2.2,'         ',I8)
        enddo
      else if (T_FLAG) then
        nobs = ODS_NGet ( ncid, jbeg, user_hour,  ier )
        print *, nobs
      else if (L_FLAG) then
        do i=4,1,-1
          nobs = ODS_NGet ( ncid, jbeg, syn_hour(i), ier )
          if ( nobs .GT. 0 ) then
             write (6,110) syn_hour(i)
110          format (I2.2)
             call ODS_close ( ncid, filename, ier )
             stop
          endif
        enddo
      else if (X_FLAG) then

        allocate ( kxtable ( nkx ) )
        do j=1,nkx
           kxtable(j) = 0
        enddo

        allocate ( kx_names ( nkx ) )
        call ODS_GetList ( ncid, "kx_names", nkx, kx_names, ier)
        if ( ier .ne. 0 ) then
           print *, 'could not get kx names from ',trim(filename)
           call exit (99)
        end if

        print *, ' '
        print *, 'Scanning ',trim(filename)
        print *, ' '
        do i=1,4
          nobs = ODS_NGet ( ncid, jbeg, syn_hour(i), ier )
          if ( ier .ne. 0 ) then
            print *, 'could not read ods file ',trim(filename)
            call exit (99)
          end if
          if ( nobs .LE. 0 ) then
             cycle
          endif
          print *, ' '
          print *, 'Date=',ymd,' Synoptic Hour=', syn_hour(i)
          print *, 'Total Number of OBS=',nobs
          print *, ' '
          allocate ( kx ( nobs ) )
          call ODS_GetI ( ncid, "kx", jbeg,  syn_hour(i), nobs,
     .                    kx, ier )
          if ( ier .ne. 0 ) then
            print *, 'could not kx from ods file ',trim(filename)
            call exit (99)
          end if

          do k=1,nobs
            if ( kx(k) .LE. nkx .AND. kx(k) .GT. 0 ) then
              kxtable( kx(k) ) = kxtable( kx(k) ) + 1
            endif
          enddo

          deallocate ( kx )

          do j=1,nkx
             if ( kxtable(j) .NE. 0 ) then
               write (6,120) j, kx_names(j), kxtable(j)
               kxtable(j) = 0
             endif
120          format (I3.3,' ',A60,I8)
          enddo
        enddo
        deallocate ( kx_names )
        deallocate ( kxtable )
      endif

!***********************************************************
!  Clean up 
!***********************************************************

      if ( ier .ne. 0 ) then
         print *, 'could not read ods file ',trim(filename)
         call exit (99)
      endif
      call ODS_close ( ncid, filename, ier )

      stop
      end

!***********************************************************
!  Usage message
!***********************************************************

      subroutine usage ()

         print *, 'Usage: ods_scan [-t synoptic_hour] [ -l ] [ -x ]',
     .            ' filename'
         print *, 'Options: -t synoptic_hour prints number of obs for'
         print *, '            the requested hour'
         print *, '         -l prints the last hour with non-zero obs'
         print *, '         -x dumps obs counts by KX'
         call exit (99)

      return
      end
