!
!  Tests RAOBV obs.
!
 Program ut_mraob

    use m_mraob

    type(raob_oms) oms, oms2
    integer rc

    nsta = 20

    call Init ( oms, nsta, rc, nsyn = 4 )
    if ( rc .ne. 0 ) then
       print *, 'error: init'
    else
       print *, 'init: ok'
    end if

    print *, 'nsta = ', oms%nsta

    oms%xm     = (/ (1.0*i, i=1,nsta) /)
    oms%itype  = (/ (2*i, i=1,nsta) /)
    oms%istatn = (/ (2*i, i=1,nsta) /)
    oms%radcor = (/ (3*i, i=1,nsta) /)
    oms%lon = 180.
    oms%lat = -45.
    
    nymd = 20020107
    nhms = 0
    rc =0

!   First create file for 0z
!   ------------------------
    call Put ( 'test.oms', nymd, nhms, oms, rc, new=.true. )
    if ( rc .ne. 0 ) then
       print *, 'error: put on ', nymd, nhms
    else
       print *, 'Put: ok on ', nymd, nhms
    end if

!   then add other hours
!   --------------------
    do ih = 6, 18, 6
        nhms = ih * 10000
        rc = 0
        call Put ( 'test.oms', nymd, nhms, oms, rc )
        if ( rc .ne. 0 ) then
           print *, 'error: put on ', nymd, nhms
        else
           print *, 'Put: ok on ', nymd, nhms
        end if
    end do


!   Read data back
!   --------------
    do ih = 0, 18, 6
        nhms = ih * 10000
        call clean ( oms, rc ) ! not needed, but just in case
        call Get ( 'test.oms', nymd, nhms, oms, rc )
        if ( rc .ne. 0 ) then
           print *, 'error: get on ', nymd, nhms
        else
           print *, 'Get: ok on ', nymd, nhms, oms%nsta
           print *, 'Get: lon/lat = ', oms%lon(5), oms%lat(5)
           print *, 'Get: xm      = ', oms%xm(1:nsta:5)
           print *, 'Get: itype   = ', oms%itype(1:nsta:5)
           print *, 'Get: radcor  = ', oms%radcor(1:nsta:5)
        end if
    end do

    call Clean ( oms, rc )
    if ( rc .ne. 0 ) then
       print *, 'error: clean'
    else
       print *, 'Clean: ok'
    end if

  end
