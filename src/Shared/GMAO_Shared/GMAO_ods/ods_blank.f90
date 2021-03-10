!
! Simple program to create an empty ODS file for a given day.
! This is handy when standard files are missing from the operational
! stream.
! Todling 2010-10-28 - handle synoptic hour individually (no daily file
!                      assumption)
!

    use m_ods

    implicit NONE

    character(len=*), parameter :: myname = 'ods_blank'

    type(ods_vect) ods
    integer rc, nymd, nhms, nymd_, nhms_
    integer ihours,syn_incr,nhms_check
    character(len=256) ifname, ftype, ofname

    integer, external :: iargc, ods_caldat
    integer argc, jbeg, jend,nsyn
    character(len=256) argv

!....................................................................

!   Parse command line
!   ------------------
    argc = iargc()
    if ( argc .ne. 4 ) call usage_()
    call GetArg(1, ifname)
    call GetArg(2, argv);  read(argv,*) nymd
    call GetArg(3, argv);  read(argv,*) nhms
    call GetArg(4, ofname)

!   Find days on file
!   ------------------
    call Days_ ( trim(ifname), jbeg, jend )
    nymd_ = ODS_CalDat ( jbeg )
    nhms_ = 0

!   Read one synoptic time
!   ----------------------
    call ODS_Get ( ifname, nymd_, nhms_, ftype, ods, rc )
    if ( rc .ne. 0 ) then
       print *, 'Error reading ', trim(ifname)
       print *, 'rc = ', rc
       call exit(7)
    end if

    if ( .not. ( trim(ftype) == 'pre_anal' .OR. trim(ftype) == 'post_anal' ) ) then
       ftype = 'pre_anal'
    end if
    ods%data%nobs = 0

!   Create empty file (Note: not RUC friendly)
!   ------------------------------------------
!   do nhms = 0, 060000, 120000, 180000
      nsyn = ods%meta%nsyn
      ihours = 240000
      syn_incr = ihours/nsyn
      nhms_check = mod(nhms,syn_incr)
      if(nhms_check .ne. 0) then
       print *, "============================================================"
       print *, "     ******        ERROR    ******                       "
       print *, " Input blank ods file contains ",nsyn," synoptic times " 
       print *, " and trying to write ",nymd,nhms
       print *, "============================================================"
       call exit(100)
      endif


    print *, 'Writing ods at ', nymd, nhms
       call ODS_Put ( ofname, ftype, nymd, nhms, ods, rc )
       if ( rc .ne. 0 ) then
          print *, 'Error writing ', trim(ofname)
          print *, 'rc = ', rc
          call exit(7)
       end if
!   end do

!   All done
!   --------
    call exit(0)

CONTAINS

!...............................................................

    subroutine usage_()

print *
print *, "NAME"
print *, "    ods_blank - Creates empty ODS file for one day"
print *, ""
print *, "SYNOPSYS"
print *, "    ods_blank.x  sample.ods  nymd  nhms output.ods"
print *, ""
print *, "DESCRIPTION"
print *, "    This utility reads an existing ODS file (sample.ods)"
print *, "    and creates a similar empty ODS file (output.ods)"
print *, "    for a different date. This useful to replace files"
print *, "    missing from the operational stream."
print *, "        "
print *, "EXAMPLE"
print *, "    ods_blank.x qscat.ods.flk.t20010511 20010512 qscat.ods.flk.t20010512 "
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



