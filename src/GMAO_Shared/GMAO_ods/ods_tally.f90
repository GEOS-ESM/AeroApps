!
! Simple program to create an empty ODS file for a given day.
! This is handy when standard files are missing from the operational
! stream.
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: odslist:  Prints contents of ods files to stdout.
!
! !INTERFACE:

      program odslist

!
! !USAGE: see the routine usage() below
!
! !DESCRIPTION: Lists contents of an ODS file in ASC format.
!
! !TODO:
!
!  option for selecting which attributes to print

! !USES:

    use m_ods

    implicit NONE

!
! !REVISION HISTORY:
!
!     ???????   Lucchsi - created
!     15Jun2004 Todling - removed restriction on synoptic times;
!                         added support for GSI-output files.
!
!EOP
!BOC


    character(len=*), parameter :: myname = 'ods_tally'

    type(ods_vect)     :: ods
    integer            :: rc, nymd, nhms
    character(len=255) :: ifname, ftype
    logical            :: ktonly = .true.

    integer, external  :: iargc, ods_caldat
    integer            :: argc, jbeg, jend, nobs, i, j, isyn, ksyn, ia
    logical	       :: ncf
    character(len=255) :: argv

!....................................................................
    ncf = .FALSE.     ! default: handling ODS file
!....................................................................

!   Parse command line
!   ------------------
    argc = iargc()
    if ( argc < 1 ) call usage_()
    call GetArg(1, ifname)
    do ia = 1, argc
       call GetArg(ia, argv)
       if ( trim(argv) .eq. 'long' .or. trim(argv) .eq. 'LONG' ) &
            ktonly = .false.
       if ( trim(argv) .eq. '-ncf' ) ncf = .true.
    end do

!   Figure out range of synoptic time loop
!   --------------------------------------
    ksyn = 32767
    if(ncf) ksyn = 1

!   For each synoptic time
!   ----------------------
    do isyn = 1, ksyn

          nymd = -1           ! get data for next synoptic time on file
          nhms =  0
          call ODSNxTime ( ifname, nymd, nhms )
          if ( nymd .eq. -1 ) then
               print *, 'End-Of-File'
               exit     ! skip to next file
          end if

          print *
          write(6,'(a,i8,a,i6)') & 
         '>>> Date: ', nymd, ', Time: ', nhms
          write(6,'(56x,a)') '  No. OBS    No. OBS'
          write(6,'(56x,a)') ' CLEAR QCX   ALL QCX'
          write(6,'(56x,a)') ' ---------  --------'

!         Read one synoptic time
!         ----------------------
          call ODS_Get ( trim(ifname), nymd, nhms, ftype, ods, rc, ncf=ncf )
          if ( rc .ne. 0 ) then
             print *, 'Error reading ', trim(ifname)
             print *, 'rc = ', rc
             call exit(7)
          else
             nobs = ods%data%nobs
          end if

!         Print tally
!         -----------
          call ods_tally ( 6, ods, nobs, rc, kt_only=ktonly )

    end do

!   All done
!   --------
    call exit(0)

CONTAINS

!...............................................................

    subroutine usage_()

print *
print *, "NAME"
print *, "    ods_tally - Prints out summary by  kx/kt"
print *, ""
print *, "SYNOPSYS"
print *, "    ods_tally.x  filename.ods  [long] [-ncf]"
print *, ""
print *, "DESCRIPTION"
print *, "    This utility reads an existing ODS file (sample.ods)"
print *, "    and prints a summary with number of observations for"
print *, "    for each data type (default) and for each data source (kx)"
print *, "    and data type when 'long' is specified. Summaries are"
print *, "    for all dates and synoptic times on file."
print *, "        "
print *, "  -ncf   used for handling non-ODS files (GSI-output)   "
print *, "        "
print *, "        "
print *, "EXAMPLE"
print *, "    ods_tally.x qscat.ods.flk.t20010511 long"
print *, ""
print *, "AUTHOR"
print *, "    Arlindo da Silva (dasilva@gsfc.nasa.gov)"
print *

  call exit(1)

  end subroutine usage_
    
!.................................................................

      Subroutine Days_ ( fname, jbeg, jend )

      character(len=*) fname
      integer jbeg, jend

!
!     Find beginning and ending Julian days on file.
!
      integer ncid, ier

      call ODS_Open ( ncid, trim(fname), 'r', ier ) ! open the file
      if ( ier .ne. 0 ) then
         print *,  myname//': could not open ods file '// trim(fname) 
         call exit(1)
      end if
      call ODS_IGet ( ncid, 'syn_beg:first_julian_day',  jbeg, ier )
      call ODS_IGet ( ncid, 'syn_beg:latest_julian_day', jend, ier )
      if ( ier .ne. 0 ) then
         print *,  myname//': could not read ods file '// trim(fname) 
         call exit(1)
      end if
      call ODS_close ( ncid, myname, ier )

      end Subroutine Days_

end 



