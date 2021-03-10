!
! Simple program to created ODS files from ASCII ODT files.
! See usage_() for additioonal details.
!
! 24Feb2003  Gaspari - fixed usage for option -t
!

    use m_ods

    implicit NONE

    character(len=*), parameter :: myname = 'ods_maker'

    type(ods_vect) mods, ods
    integer rc, nymd, nhms, nymd_o, nhms_o, mobs
    character(len=255) ifname, mfname, ftype, mftype, ofname

    integer, external :: iargc, ods_caldat
    integer:: twindow(2), toffset
    integer i, n, iarg, argc, jbeg, jend
    character(len=255) argv

    integer, parameter :: NOBS_MAX = 450 * 1000
    integer :: lu = 10

    logical :: verbose 

!....................................................................

!   Parse command line
!   ------------------
    call Init_()

!   A clean slate
!   -------------
    call ods_init ( ods, mobs, rc )

!   Main loop
!   ---------
    n = 0
    open(lu, file=trim(ifname), status='old' )
    do

       n = n + 1
       if ( n > ods%data%nvct ) then
          print *, myname//': too many obs, ', n, mobs
          print *, 'rc = ', rc
          call exit(7)
       end if

       call nxtobs_ ( lu, ftype, nhms_o, twindow, toffset, &
                      ods%data%kx(n),   ods%data%ks(n), ods%data%kt(n), &
                      ods%data%time(n), ods%data%lev(n), ods%data%lon(n), &
                      ods%data%lat(n),  ods%data%obs(n), ods%data%omf(n), &
                      ods%data%oma(n),  ods%data%qcexcl(n), ods%data%qchist(n), &
                      ods%data%xm (n),  ods%data%xvec(n), rc )
      
!     End of file
!     -----------
      if ( rc < 0 ) then
           ods%data%nobs = n - 1
           exit
      end if

    end do

!   Set ODS metadata
!   ----------------
    call set_meta_()

!   Create empty file (Note: not RUC friendly)
!   ------------------------------------------
    call ODS_Put ( ofname, ftype, nymd_o, nhms_o, ods, rc )
    if ( rc .ne. 0 ) then
       print *, myname//': error writing ', trim(ofname)
       print *, 'rc = ', rc
       call exit(7)
    end if

!   All done
!   --------
    call exit(0)

CONTAINS

!...............................................................

    subroutine init_()

    integer ic,lt
    real    swap
    character*255 SS

!   Defaults
!   --------
    ifname = 'odslist.odt'
    ofname = 'ods_maker.ods'
    ftype  = 'post_analysis'
    mobs   = NOBS_MAX
    verbose = .false.
    twindow(1) = -9999999
    twindow(2) =  9999999
    toffset    =  0

!   Parse command line
!   ------------------
    argc = iargc()
    if ( argc < 2 ) call usage_()

    iarg = 0
    do i = 1, argc
       iarg = iarg + 1
       if ( iarg .gt. argc ) exit
       call GetArg ( iArg, argv )
       if (index(argv,'-i' ) .gt. 0 ) then
          if ( iarg+1 .gt. argc ) call usage_()
          iarg = iarg + 1
          call GetArg ( iArg, ifname )
       else if (index(argv,'-o' ) .gt. 0 ) then
          if ( iarg+1 .gt. argc ) call usage_()
          iarg = iarg + 1
          call GetArg ( iArg, ofname )
       else if (index(argv,'-t ' ) .gt. 0 ) then
          if ( iarg+1 .gt. argc ) call usage_()
          iarg = iarg + 1
          call GetArg ( iArg, ftype )
       else if (index(argv,'-n' ) .gt. 0 ) then
          if ( iarg+1 .gt. argc ) call usage_()
          iarg = iarg + 1
          call GetArg ( iArg, argv ); read(argv,*) mobs
       else if (index(argv,'-v' ) .gt. 0 ) then
          verbose = .true.
       elseif (index(argv,'-toffset') .gt. 0 ) then
          if ( iarg+1 .gt. argc ) call usage_()
          iarg = iarg + 1
          call GetArg ( iArg, SS )
          read(SS,*) toffset
       elseif (index(argv,'-time') .gt. 0 ) then
          if ( iarg+1 .gt. argc ) call usage_()
          iarg = iarg + 1
          call GetArg ( iArg, SS )
          ic = index(SS,':')     ! string is t1 or t1:t2
          lt = len_trim(SS)
          if (ic .eq. 0) then    ! no colon, therefore t1
              read(SS,*) twindow(1)
              twindow(2) = twindow(1)
          else                   ! colon, therefore t1:t2
              read(SS(1:ic-1) ,*) twindow(1)
              read(SS(ic+1:lt),*) twindow(2)
          end if
          if (twindow(2)<twindow(1)) then    ! let's be nice to the user..
              swap = twindow(2)
              twindow(2)   = twindow(1)
              twindow(1)   = swap
          end if
       else
          if ( iArg+2 > argc ) call usage_()
          call GetArg ( iArg,   mfname )
          call GetArg ( iArg+1, argv ); read(argv,*) nymd_o
          call GetArg ( iArg+2, argv ); read(argv,*) nhms_o
          exit
       end if
    end do

!   Echo the parameters
!   -------------------
    if ( verbose ) then
      print *
      print *, "ods_maker - Creates ODS file from ASCII ODT file"
      print *
      print *, '     Input ASCII file: ', trim(ifname)
      print *, 'ODS file for Metadata: ', trim(mfname)
      print *, '    Output   ODS file: ', trim(ofname)
      print *, '    Output  file type: ', trim(ftype)
      print *, '   Date/Time/nobs_max: ', nymd_o, nhms_o, mobs
      print *, '   Within time window: ', twindow
      print *, '     With time offset: ', toffset
      print *
    end if

    end subroutine init_

!...............................................................

    subroutine usage_()

print *
print *, "NAME"
print *, " ods_maker - Creates ODS file from ASCII ODT file"
print *, ""
print *, "SYNOPSYS"
print *, "    ods_maker.x  [...OPTIONS...]  metafn nymd nhms"
print *, ""
print *, "DESCRIPTION"
print *, "    This utility reads a simple rc file containing an ASCII"
print *, "    specification of an ODS 'Observation Vector' (lat, lon, ...)"
print *, "    and creates an ODS file for a given date and time, taking "
print *, "    metadata from a pre-existing ODS file (metafn). This utility is"
print *, "    useful for creating ODS files for 1-observation experiments."
print *, "    This utility is designed to work with (stripped) output from"
print *, "    utility odsselect. The default time is nhms = 00000. "
print *, "        "
print *, "OPTIONS"
print *, "    -i fname     input ASCII ODT file name (default: odslist.odt)"
print *, "    -o fname     output file name (default: ods_maker.ods) "
print *, "    -t ftype     ODS file type: 'pre_anal' or 'post_anal' (default)"
print *, "    -n nobs_max  maximum number of obs possible, default = ", NOBS_MAX
print *, "    -v           verbose: echoes obs to screen    "
print *, "        "
print *, "EXAMPLE"
print *, "    ods_maker.x qscat.ods.flk.t20010511 20010512 120000"
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
         print *,  myname//': could not open meta ods file '// trim(fname) 
         call exit(1)
      end if
      call ODS_IGet ( ncid, 'syn_beg:first_julian_day',  jbeg, ier )
      call ODS_IGet ( ncid, 'syn_beg:latest_julian_day', jend, ier )
      if ( ier .ne. 0 ) then
         print *,  myname//': could not read meta ods file '// trim(fname) 
         call exit(1)
      end if
      call ODS_close ( ncid, myname, ier )

      end Subroutine Days_


!.................................................................

     subroutine nxtobs_ ( lu, ftype, nhms, twindow, toffset, &
                          kx, ks, kt, time, lev, lon, lat, &
                          obs, omf, oma, qcx, qch, xm, xvec, rc )

     integer, intent(in)  :: lu
     character(len=255), intent(in) :: ftype
     integer, intent(in)  :: twindow(2), toffset
     integer, intent(in)  :: nhms
     integer, intent(out) :: kx, ks, kt, time
     real,    intent(out) :: lev, lon, lat
     integer, intent(out) :: qcx, qch
     real,    intent(out) :: obs, omf, oma, xm, xvec
     integer, intent(out) :: rc

     character(len=255) line
     integer              :: kkx, kks, kkt, ttime
     real                 :: llev, llon, llat
     integer              :: qqcx, qqch
     real                 :: oobs, oomf, ooma, xxm, xxvec
     integer ls, ios

     omf  = obs_missing
     oma  = obs_missing
     xvec = obs_missing

     do 

!        Read next line on file
!        --------------
         read(lu,'(a)', iostat=rc) line ! read next line

!        Premature EOF
!        -------------
         if ( rc < 0 ) then
              exit
         else if ( rc > 0 ) then
              print *, myname//': error reading '//trim(ifname)
              call exit(1)
         end if
           
!        Cycle if comment or empty line
!        ------------------------------
         ls = min(254,len_trim(line)) + 1
         line(ls:ls) = '#'
         ls = index(line,'#' ) - 1    ! line length
         if ( ls < 1 ) cycle

!        Read observation
!        ----------------
         if ( ftype(1:4) .eq. 'post' .or. ftype(1:4) .eq. 'POST' ) then
              read(line(1:ls),*,iostat=ios) &
                          kkx, kks, kkt, ttime, llev, llon, llat, & 
                          oobs, oomf, ooma, qqcx, qqch, xxm, xxvec
              if ( ios /= 0 ) then
                 print *, myname//': error decoding post-analysis line:'
                 print *, '|'//trim(line(1:ls))//'|'
                 call exit(1)
              end if
              if ( twindow(1)<=ttime .and. ttime<=twindow(2)) then
                   kx   = kkx
                   ks   = kks
                   kt   = kkt
                   time = ttime - toffset
                   lev  = llev
                   lon  = llon
                   lat  = llat
                   obs  = oobs
                   omf  = oomf
                   oma  = ooma
                   qcx  = qqcx
                   qch  = qqch
                   xm   = xxm
                   xvec = xxvec
              else
                  cycle
              endif
        else if ( ftype(1:3) .eq. 'pre' .or. ftype(1:3) .eq. 'PRE' ) then
              read(line(1:ls),*,iostat=ios) &
                          kkx, kks, kkt, ttime, llev, llon, llat, &
                          oobs, qqcx, qqch, xxm
              if ( ios /= 0 ) then
                 print *, myname//': error decoding pre-analysis line:'
                 print *, '|'//trim(line(1:ls))//'|'
                 call exit(1)
              end if
              if ( twindow(1)<=ttime .and. ttime<=twindow(2)) then
                   kx   = kkx
                   ks   = kks
                   kt   = kkt
                   time = ttime - toffset
                   lev  = llev
                   lon  = llon
                   lat  = llat
                   obs  = oobs
                   qcx  = qqcx
                   qch  = qqch
                   xm   = xxm
              endif
       else
             print *, myname//': unknown file type '//trim(ftype)
             call exit(1)
       end if

       if ( verbose ) then
          print *, kx, ks, kt, time, lev, lon, lat, &
                   obs, omf, oma, qcx, qch, xm, xvec
       end if
                   
!      Convert ECMWF kx/kt to GSI-like
!      -------------------------------
       call eckt2gsikt ( kt )
       call eckx2gsikx ( kx, kt )

       exit

     end do

     end subroutine nxtobs_

     subroutine eckx2gsikx ( kx, kt )
     implicit none
     integer,intent(inout) :: kx
     integer,intent(in)    :: kt
     if ( kx==326  ) then
          if(kt==4 .or. kt==5) then
             kx=229         ! PILOT Land Report (Winds)
          else
             kx=129         ! should never be the case (not NCEP kx)
          endif
          return
     endif
     if ( kx==346  ) then
          if(kt==4 .or. kt==5) then
             kx=223         ! American Wind Profiler
          else
             kx=123         ! should never be the case (not NCEP kx)
          endif
          return
     endif
     if ( kx==355  ) then
          if(kt==4 .or. kt==5) then
             kx=220 !281    ! TEMP Land Report
          else
             kx=120 !181    ! TEMP Land Report
          endif
          return
     endif
     if ( kx==365  ) then
          kx=183         ! TEMP SHIP Report
          return
     endif
     if ( kx==1316 ) then
          if(kt==4 .or. kt==5) then
             kx=228         ! Japanese Wind Profiler
          else
             kx=128         ! should never be the case (not NCEP kx)
          endif
          return
     endif
     if ( kx==1346 ) then
          if(kt==4 .or. kt==5) then
             kx=223         ! European Wind Profiler
          else
             kx=123         ! should never be the case (not NCEP kx)
          endif
          return
     endif
     if ( kx==1355 ) then
          if(kt==4 .or. kt==5) then
            kx=220         ! TEMP Dropsonde Report 
          else
            kx=120         ! TEMP Dropsonde Report 
          endif
          return
     endif
     end subroutine eckx2gsikx

     subroutine eckt2gsikt ( kt )
     use m_odsmeta
     implicit none
     integer,intent(inout) :: kt
     if ( kt==2 ) then
          kt=ktTv 
          return
     endif
     if ( kt==3 ) then
          kt=ktuu
          return
     endif
     if ( kt==4 ) then
          kt=ktvv
          return
     endif
     if ( kt==7 ) then
          kt=ktqq
          return
     endif
     if ( kt==41) then
          kt=ktus10
          return
     endif
     if ( kt==42) then
          kt=ktTs10
          return
     endif
     end subroutine eckt2gsikt

!.........................................................................

    subroutine set_meta_()

!   Find days on meta ODS file
!   --------------------------
    call Days_ ( trim(mfname), jbeg, jend )
    nymd = ODS_CalDat ( jbeg )
    nhms = 0        

!   Read one synoptic time
!   ----------------------
    call ODS_Get ( mfname, nymd, nhms, mftype, mods, rc )
    if ( rc .ne. 0 ) then
       print *, myname//':error reading ', trim(ifname)
       print *, 'rc = ', rc
       call exit(7)
    end if

    ods%meta = mods%meta  ! shallow copy

  end subroutine set_meta_


end 



