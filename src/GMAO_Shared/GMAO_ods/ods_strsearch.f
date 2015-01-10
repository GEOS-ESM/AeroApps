

*..............................................................


      integer function ODS_StrSearch ( SearchStr,
     .                                 NStrings,
     .                                 Strings,
     .                                 ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:      ODS_StrSearch
!
! !DESCRIPTION: 
!
!     Searches for the string, SearchStr, among a list of strings
!     and returns the smallest list index number that is identical
!     to SearchStr.  If none of the list entries contain the 
!     string, SearchStr, then the function returns a value of 
!     Not_Found ( =0 ).
!
! !INTERFACE: list_entry = ODS_StrSearch ( SearchStr,
!                                          NStrings,
!                                          Strings,
!                                          ierr )
!
! !INPUT PARAMETERS:

      character  SearchStr * (*)    ! The string to be located
      integer    NStrings           ! The number of strings to be
                                    !   searched
      character  Strings  ( NStrings ) * (*)
                                    ! The array of strings to be
                                    !   searched

! !OUTPUT PARAMETER:
      integer         ierr          ! The returned error code
!
! !FILES USED:
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_stdio.h, a header file, for defining standard input/output
!            unit numbers
!
! !REVISION HISTORY: 
!     09Apr1996   Redder   Origional version
!     06Apr1998   Redder   Replaced the algorithm for determining
!                          the string length
!     10Jan2003   Redder   Replaced algorithm for determining maximum
!                          string length to avoid array index error
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_stdio.h'

*     Function referenced
*     -------------------
      integer     ODS_StrSize   ! Determines the string size
                                !   excludeing trailing blanks

*     Other variables
*     ---------------
      integer     iString       ! index for do loop
      integer     len_SearchStr ! Length of the string, SearchStr
      integer     len_String    ! Length of a string in Strings
      integer     len_StringMax ! Maximum size of len_String
      integer     list_entry    ! Index number of the list entry return
                                !   from the search string function,
                                !   ODS_StrSearch.  If no string is
                                !   found then the value is equal to
                                !   the parameter, Not_Found
      integer     Not_Found     ! Symbolic constant for the case when
                                !   the search was not successful
      parameter ( Not_Found = 0 )

      ierr = NCNoErr   ! Default error code

*     Determine the length of the search string
*     -----------------------------------------
      len_SearchStr = ODS_StrSize ( SearchStr )
      len_StringMax = len_SearchStr ! Default maximum value

      list_entry    = Not_Found     ! default index value

*     For each string in Strings ...
*     ------------------------------
      do 10, iString = NStrings, 1, - 1

*        Determine the length of the string being examined and its max
*        -------------------------------------------------------------
         len_String    = ODS_StrSize ( Strings ( iString ) )
         len_StringMax =         len ( Strings ( iString ) )

*        Do not even perform the strings comparisons if the
*        lengths do not even match
*        --------------------------------------------------
         if ( len_String .eq. len_SearchStr ) then

*           Examine the ith string to determine if there is a match
*           -------------------------------------------------------
            if ( Strings ( iString ) ( : len_String ) .eq.
     .           SearchStr           ( : len_SearchStr ) )
     .         list_entry = iString
         end if
 10   continue

*     If the length of the search string is greater than
*     the maximum size of each entry in strings then
*     return with an error message and a value of Not_Found
*     -----------------------------------------------------
      if ( len_SearchStr .gt. len_StringMax ) then
         ierr          = NCSysErr
         write ( stderr, 901 ) SearchStr
         ODS_StrSearch = Not_Found
         return
      end if

      ODS_StrSearch = list_entry

      return
*     ------

 901  format ( /, ' ODS_StrSearch : The size of the search string ',
     .         /, '                 exceeds the maximum of each ',
     .         /, '                 string to be searched. ',
     .         /, '                 The search string = ', a )

      end
