

*..............................................................


      character * ( * ) function ODS_Case ( string, case )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:      ODS_Case
! 
! !DESCRIPTION: 
!     Sets all alphabetical character strings to upper or lower
!     case as determined by the calling routine.  For the case
!     when the lengths of the input and outstring differ, this
!     routine processes the data according the ANSI standard
!     FORTRAN rules for character assignments.  That is, if the
!     output string length is larger than the input string length,
!     then  modified string is stored in the left portion of the
!     output string and the remaining portion is padded with
!     blanks.  If the output string length is smaller, then only
!     the leftmost characters of the input string are modified
!     and stored.
!
! !INTERFACE: new_string = ODS_Case  ( string, case )
!
! !INPUT PARAMETERS:
      character * (*) string ! The string to be modified
      character * (*) case   ! = 'upper' if upper case is desired
                             ! = 'lower' if lower case is desired
                             ! Any other options produce no changes
                             ! The routine is case sensitive to
                             ! this input parameter.
!
! !REVISION HISTORY: 
!     14Feb1996   Redder   Origional version
!     01Nov1999   Redder   Revised code to prevent subscript errors
!                          in character strings
!     24Nov1999   Redder   Insured that any portion of the output
!                          string not assigned by the modified
!                          input string is padded with blanks.
!
!-------------------------------------------------------------------------

*     storage for old and new element in string
*     -----------------------------------------
      character old_element * ( 1 )
      character new_element * ( 1 )

*     maximum and minimum index values for
*     for upper and lower case letters
*     ------------------------------------
      integer   lower_min
      integer   lower_max
      integer   upper_min
      integer   upper_max

*     Other variables
*     ---------------
      integer   i_element     ! index value for character element
      integer   delta_case    ! difference of between the index
                              !   values of the upper and lower
                              !   case of the same letter
      integer   InLen         ! length of input string
      integer   OutLen        ! length of output string
      integer   StrLen        ! minimum of InLen and OutLen
      integer   iStr          ! index variable for do loop
      character case_ * ( 5 ) ! temporary array for input argument

*     Determine upper and lower range of indicies of
*     upper and lower case alphabetical characters
*     ----------------------------------------------
      lower_max   = max ( ichar ( 'a' ), ichar ( 'z' ) )
      lower_min   = min ( ichar ( 'a' ), ichar ( 'z' ) )
      upper_max   = max ( ichar ( 'A' ), ichar ( 'Z' ) )
      upper_min   = min ( ichar ( 'A' ), ichar ( 'Z' ) )

*     Miscellaneous variables
*     -----------------------
      case_       = case
      delta_case  = ichar ( 'A' ) - ichar ( 'a' )
      InLen       = len ( string   )
      OutLen      = len ( ODS_Case )
      StrLen      = min ( InLen, OutLen )

*     For each character element in string
*     ------------------------------------
      do 10, iStr = 1, StrLen

         old_element = string ( iStr: iStr )
         new_element = old_element
         i_element   = ichar ( old_element )

*        if upper case letters are desired ...
*        -------------------------------------
         if   ( case_ ( 1:5 ) .eq. 'upper' ) then

*        and the element is a lower case letter ...
*        ------------------------------------------
         if ( i_element .ge. lower_min .and.
     .        i_element .le. lower_max ) then

*           then make the element an upper case letter 
*           ------------------------------------------
            new_element = char ( ichar ( old_element ) + delta_case )
         end if

*        if lower case letters are desired ...
*        -------------------------------------
         else if ( case_ ( 1:5 ) .eq. 'lower' ) then

*        and the element is an upper case letter ...
*        -------------------------------------------
         if ( i_element .ge. upper_min .and.
     .        i_element .le. upper_max ) then

*           then make the element an lower case letter 
*           ------------------------------------------
            new_element = char ( ichar ( old_element ) - delta_case )
         end if

         end if

         ODS_Case ( iStr: iStr ) = new_element

 10   continue

*     Pad with blanks any portion of the output string
*     that is not defined by the modified input string.
*     -------------------------------------------------
      if ( OutLen .gt. InLen ) ODS_Case ( InLen + 1 : ) = ' '

      return
      end
