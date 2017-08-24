!-------------------------------------------------------------------------
!         NASA/GSFC, Global Modelling Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOI
!
! !TITLE: An Utility to trasnform X synoptic ods file to Y synoptic ods file.
!
! !AUTHORS: Ravi C. Govindaraju                 
!
! !AFFILIATION: Global Modelling Assimilation Office, NASA/GSFC, Greenbelt, MD 20771
!
! !DATE: August 4th 2015 
!
! !INTRODUCTION: System Overview
!     
!         Transform 4 synoptic ods file to 8 synoptic ods file.
!
!
!EOI
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!         NASA/GSFC, Global Modelling Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  ods_synX2Y()
! 
! !DESCRIPTION: Utility to transform the ods files.! !INTERFACE:
!
       Program ods_synX2Y

! !USES:

       Use m_ods
       Use m_StrTemplate
       Use m_stdio, only : stdout
       Use m_die

       Implicit NONE

      character(len=*), parameter :: myname = 'ods_synX2Y'
      character*(255) ifname
      character*(255) ofname,token
      character(len=2048)      argv
      character*8 ftype, dtype, otype
      logical     ncf,verbose

      integer                  nobs,nsyn,nymd,nhms,isyn,ksyn,ios
      integer                  i,inc_out
      integer                  nhms0,hhmmss(24),jsyn,hhmmss0
      integer                  argc,ier,rc,iarg
      integer, external     :: iargc
      integer                  nkount,nkx,nkt
      type     ( ods_vect )    ods
      type     ( ods_vect )    ods_out
      type     ( ods_vect )    ods_temp
      type     ( ods_vect )    mods
!     data hhmmss /000000,030000,060000,090000,120000,150000,180000,210000/
!     data nkount,nkt,nkx/10000,400,600/



       ftype = 'ods'
       otype = 'nc4'
       ncf   = .false.

       call init_()

       nsyn = 24/inc_out   
       hhmmss0 = 0
       do i = 1,nsyn
        hhmmss(i) = hhmmss0 *10000
        hhmmss0 = hhmmss0 + inc_out
       end do
      
       ksyn = 32767
       isyn = 1
       call ODS_Get ( trim(ifname), nymd, nhms, ftype, ods, ier, ncf=ncf )
       ksyn = ods%meta%nsyn
       nkt = ods%meta%nkt
       nkx    = ods%meta%nkx

       nobs   =  ods%data%nobs

       ods_out%meta%nsyn = nsyn
       print *,'nkx: ',nkx,' nkt: ',nkt
       call ods_init(ods_out,nkount,rc,nkx=nkx,nkt=nkt,nsyn=nsyn)
       call ods_init(ods_temp,nkount,rc,nkx=nkx,nkt=nkt,nsyn=nsyn)

       do  isyn = 1, ksyn
        nhms0 = hhmmss(isyn)
        call ODS_Get ( trim(ifname), nymd, nhms, ftype, ods, ier, ncf=ncf )
        nobs   =  ods%data%nobs
        if(nobs == 0 ) then
         ods_temp = ods
        endif
       end do

       do  isyn = 1, nsyn
         nhms0 = hhmmss(isyn)
         call ODS_Get ( trim(ifname), nymd, nhms0, ftype, ods, ier, ncf=ncf )
         nobs   =  ods%data%nobs
         jsyn = ods%meta%nsyn
         print *,' isyn,nymd,nhms0,nhms,jsyn,nobs: ',isyn,nymd,nhms0,nhms,jsyn,nobs

         if(nobs .ne. 0 .and. nhms0 .ne. nhms ) then
          ods_out = ods_temp
         endif

         if(nobs .ne. 0 .and. nhms0 == nhms ) then
          ods_out = ods
         else
          ods_out = ods_temp
         endif

         if ( nobs .eq. 0 ) then
                     print *, 'No data for this synoptic time'
                     call ods_clean ( ods, ier )
!                    cycle    ! skip to next synoptic time

         endif 

          ods_out%meta%nsyn = nsyn
          ods_out%data%nobs = 0

         if ( nhms0 == nhms ) then
          ods_out%meta%nsyn = nsyn
          ods_out%data%nobs = nobs
          call ODS_Put ( trim(ofname), 'post_analysis', nymd, nhms0, ods_out, rc )
          if ( rc .ne. 0 ) then
           call die ( myname, 'could not write ODS file '//trim(ofname) )
          endif
         endif

        end do

   CONTAINS
!.................................................................
         subroutine init_()
                                                                                                                             
!          --------------------------------------
!            Get the user given information.
!          --------------------------------------
                                                                                                                             
!   Defaults
!   --------
         ofname    = 'aod_obs.ods'
         inc_out   = 03
         verbose  = .false.
         argc = iargc()
         print *,' argc ==> ',argc
         if ( argc < 1 ) call usage_()
         iarg = 0
         do i = 1, 32767
            iarg = iarg + 1
            if ( iarg .gt. argc ) exit
            call GetArg ( iArg, argv )
            if(index(argv,'-inc') .gt. 0 ) then
               if ( iarg+1 .gt. argc ) call usage_()
               iarg = iarg + 1
               call GetArg ( iArg, argv ); read(argv,*) inc_out
            else if (index(argv,'-v' ) .gt. 0 ) then
               verbose = .true.
            else 

              if ( iarg .gt. argc ) call usage_()
              call GetArg ( iArg, ifname )
              iarg = iarg + 1
              if ( iarg .gt. argc ) call usage_()
              call GetArg ( iArg, ofname )
              iarg = iarg + 1
              if ( iarg .gt. argc ) call usage_()
              call GetArg ( iArg, argv );read(argv,*) nymd
              iarg = iarg + 1
              if ( iarg .gt. argc ) call usage_()
              call GetArg ( iArg, argv );read(argv,*) nhms
            endif
         end do
 
!   Echo the parameters
!   -------------------
       if ( verbose ) then
         print *
         print *, "ods_synX2Y- rewrites an ODS file with inc_out bined ODS file"
         print *
         print *, '    -inc : ', inc_out
         print *, '     input ODS: ', trim(ifname)
         print *, '     onput ODS: ', trim(ofname)
         print *, '     nymd     : ', nymd
         print *, '     nhms     : ', nhms
 
      end if
 
      print *
    end subroutine init_
!
!.................................................................
subroutine usage_()
 
!
! Prints usage notice.
!
 
print *, "NAME"
print *, "   ods_synX2Y  Transferring input ods file to user given synoptic binned "
print *, "                ods file."           
print *, "   USAGE:  ods_synX2Y.x -inc HH in_ods_file out_ods_file yyyymmdd hhmmss"
print *
print *, "SYNOPSIS"
print *, "   ods_synX2Y  [options] input_ods_fname output_ods_file name nymd nhms"
print *
print *, "OPTIONS"
print *
 
print *, "  -inc                       Synoptic increment HH (default: 08)"
print *, "  -v     verbose:           echoes obs to screen    "
print *
 
    stop
 
    end subroutine usage_
!..............................................................

  end Program ods_synX2Y

