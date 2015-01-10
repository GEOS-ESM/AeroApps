
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Make_Fname --- Creates file name from tokens
! 
! !INTERFACE:

      subroutine Make_Fname ( prefix, nymd, nhms, filename )
!
! !USES:
!
      implicit NONE

! !INPUT PARAMETERS: 
!
      character*(*)  prefix    ! run identifier
      integer        nymd      ! year-month-day
      integer        nhms      ! hour-minute-sec

! !OUTPUT PARAMETERS:
!
      character*255  filename  ! file name generated

! !DESCRIPTION: Creates a finame from toikens and dates.
!
! !REVISION HISTORY: 
!
!  06Mar98   da Silva   Revised prologue.
!
!EOP
!-------------------------------------------------------------------------

      character*32   timecode

      integer  julday
      integer  year, month, day, jday, hour

      year  = int(nymd/10000)
      month = int((nymd-year*10000)/100)
      day   = int(nymd-year*10000-month*100)
      if (year .lt. 100) year = year + 1900
      hour  = int(nhms/10000)
      
      write (timecode,'(".t",i4.4,i2.2,i2.2,".",i2.2,"Z")') 
     &       year, month, day, hour
      filename = trim(prefix)//trim(timecode)

      return
      end

      subroutine gritas_perror(msg)
      character*(*) msg
      print *
      print *, msg
      print *
      call exit(7)
      end

