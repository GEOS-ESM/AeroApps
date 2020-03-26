!
! Unit test for m_ods and portions of the ODS library
!

    use m_ods

    implicit NONE

    type(ods_vect) ods
    integer rc, nymd, nhms
    character*80 fname, ftype

    print *, 'Starting m_ods unit test....'

    print *, 'Enter nymd, nhms: '
    read *, nymd, nhms

!    nymd = 19980106
!    nhms = 180000

    fname = 'test.ods'

!   Read one synoptic time
!   ----------------------
    call ODS_Get ( fname, nymd, nhms, ftype, ods, rc )
    if ( rc .ne. 0 ) then
       print *, 'Error reading ', fname(1:60)
       print *, 'rc = ', rc
       call exit(7)
    else
       print *, 'ut_ods: nymd, nhms, nobs: ', nymd, nhms, ods%data%nobs
       print *, 'ut_ods: sucessful reading of ', fname(1:60)
    end if

!   Write it
!   --------
    fname = 'output.ods'
    call ODS_Put ( fname, ftype, nymd, nhms, ods, rc )
    if ( rc .ne. 0 ) then
       print *, 'Error writing ', fname(1:60)
       print *, 'rc = ', rc
       call exit(7)
    else
       print *, 'ut_ods: sucessful writing of ', fname(1:60)
    end if


!   All done
!   --------
    end

