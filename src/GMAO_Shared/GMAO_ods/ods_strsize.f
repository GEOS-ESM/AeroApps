

*..............................................................


      integer function ODS_StrSize    ( String )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:      ODS_StrSize
!
! !DESCRIPTION: 
!
!     Determines the index number of the last non-blank character
!     in the input string.  This number is defined as the length
!     of string which may contain embedded blanks.
!
! !INTERFACE: StrSize = ODS_StrSize ( String )
!
! !INPUT PARAMETER:
      character  String * (*)       ! The string to be examined
!
! !REVISION HISTORY: 
!     17May1996   Redder   Origional version
!     01Nov1999   Redder   Rewrote algorithm to match that for I90
!
!-------------------------------------------------------------------------

      integer  StrLen  ! length of input string
      integer  StrSize ! temporary storage for output
      integer  iStr    ! index variable for do loop

      StrLen  = len ( String )
      StrSize = 1  ! Default for null string
      do iStr = StrLen, 1, -1
         StrSize = iStr
         if ( String ( iStr:iStr ) .ne. ' ' ) go to 11

      end do
      StrSize = StrSize - 1
 11   continue

      ODS_StrSize = StrSize

      return
      end 
