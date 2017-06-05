      PROGRAM trainfo
      
c     ***********************************************************************
c     * Get infos for a trajectory file                                     *
c     * Michael Sprenger / Spring, summer 2010                              *
c     ***********************************************************************

      implicit none
      
c     ----------------------------------------------------------------------
c     Declaration of variables
c     ----------------------------------------------------------------------

c     Input file
      character*80                           inpfile     ! Input filename
      character*80                           mode        ! Mode

c     Trajectories
      integer                                ntra        ! Number of trajectories
      integer                                ntim        ! Number of times
      integer                                ncol        ! Number of columns
      character*80                           vars(100)   ! Variable names
      integer                                refdate(6)  ! Reference date
      real,allocatable, dimension (:,:,:) :: tra         ! Trajectories (ntra,ntim,ncol)

c     Auxiliary variables
      integer                                inpmode
      integer                                stat
      integer                                fid
      integer                                i,j,n
      integer                                old(5)
      integer                                new(5)
      integer                                hour,min
      character*20                           datestr
      character*120                          str
      character*4                            str1
      character*2                            str2,str3,str4,str5,str6
      integer                                isok
      integer                                nleaving

c     ----------------------------------------------------------------------
c     Do the reformating
c     ----------------------------------------------------------------------

c     Read parameters
      open(10,file='trainfo.param')
       read(10,*) inpfile
       read(10,*) mode
      close(10)
      
c     Determine the formats
      call mode_tra(inpmode,inpfile)
      if (inpmode.eq.-1) inpmode=1

c     Get the dimension of the trajectory file
      call info_tra(inpfile,ntra,ntim,ncol,inpmode)

c     Get haeder information
      call ropen_tra(fid,inpfile,ntra,ntim,ncol,refdate,vars,inpmode)
      call close_tra(fid,inpmode)

c     Write dimensions
      if ( (mode.eq.'dim').or.(mode.eq.'all') ) then
         print*,ntra,ntim,ncol
      endif

c     Write single dimensions
      if ( (mode.eq.'ntra').or.(mode.eq.'all') ) then
         print*,ntra
      endif
      if ( (mode.eq.'ntim').or.(mode.eq.'all') ) then
         print*,ntim
      endif
      if ( (mode.eq.'ncol').or.(mode.eq.'all') ) then
         print*,ncol
      endif

c     Write variable names
      if ( (mode.eq.'vars').or.(mode.eq.'all') ) then
         print*,(trim(vars(i))//' ',i=1,ncol)         
      endif

c     Write reference date 
      if ( (mode.eq.'refdate').or.(mode.eq.'all') ) then
         
c        Concatenate date string
         min = refdate(5)         
         call datestring(datestr,
     >             refdate(1),refdate(2),refdate(3),refdate(4) )
         if ( min.eq.0 ) then
            datestr = trim(datestr)//'00'
         elseif (min.lt.10) then
            datestr = trim(datestr)//'0'//char(ichar('0')+min)
         else
            datestr = trim(datestr)//
     >                char(ichar('0')+int(min/10))//
     >                char(ichar('0')+mod(min,10))
         endif
                    
c        Write date string
         print*,trim(datestr)

      endif

c     Load all trajectory times if necessary
      if ( (mode.eq.'all'      ).or.
     >     (mode.eq.'times'    ).or.
     >     (mode.eq.'startdate').or.
     >     (mode.eq.'enddate'  ) )
     >then
         allocate(tra(ntra,ntim,ncol),stat=stat)
         call ropen_tra(fid,inpfile,ntra,ntim,ncol,refdate,vars,inpmode)
         call read_tra (fid,tra,ntra,ntim,ncol,inpmode)
         call close_tra(fid,inpmode)
      endif

c     Write list of all times
      if ( (mode.eq.'times').or.(mode.eq.'all') ) then
         write(*,'(100f8.2)') (tra(1,i,1),i=1,ntim)
      endif

c     Write time range
      if ( (mode.eq.'timerange').or.(mode.eq.'all') ) then
         write(*,'(i6)') refdate(6)
      endif
      
c     Write firstdate (reference date + first time)
      if ( (mode.eq.'startdate').or.(mode.eq.'all') ) then
         
c        Set the time shift of first time relative to reference date
         hour = int(tra(1,1,1))
         min  = mod(nint(100.*tra(1,1,1)),100) + refdate(5)
         if (min.gt.60) then
            min  = min - 60
            hour = hour + 1
         endif
         if (min.lt.0) then
            min  = min + 60
            hour = hour - 1
         endif
         if (min.lt.0) then
            min  = min + 60
            hour = hour - 1
         endif

c        Get new date (hours and minutes)
         old(1) = refdate(1)
         old(2) = refdate(2)
         old(3) = refdate(3)
         old(4) = refdate(4)
         old(5) = 0
         call newdate(old,real(hour),new)

c        Concatenate the date string
         call datestring(datestr,
     >             new(1),new(2),new(3),new(4) )
         if ( min.eq.0 ) then
            datestr = trim(datestr)//'00'
         elseif (min.lt.10) then
            datestr = trim(datestr)//'0'//char(ichar('0')+min)
         else
            datestr = trim(datestr)//
     >                char(ichar('0')+int(min/10))//
     >                char(ichar('0')+mod(min,10))
         endif
                    
c        Write date string
         print*,trim(datestr)

      endif

c     Write enddate (reference date + last time)
      if ( (mode.eq.'enddate').or.(mode.eq.'all') ) then
         
c        Set the time shift of first time relative to reference date
         hour = int(tra(1,ntim,1))
         min  = mod(nint(100.*tra(1,ntim,1)),100) + refdate(5)
         if (min.gt.60) then
            min  = min - 60
            hour = hour + 1
         endif
         if (min.lt.0) then
            min  = min + 60
            hour = hour - 1
         endif
         if (min.lt.0) then
            min  = min + 60
            hour = hour - 1
         endif

c        Get new date (hours and minutes)
         old(1) = refdate(1)
         old(2) = refdate(2)
         old(3) = refdate(3)
         old(4) = refdate(4)
         old(5) = 0
         call newdate(old,real(hour),new)

c        Concatenate the date string
         call datestring(datestr,
     >             new(1),new(2),new(3),new(4) )
         if ( min.eq.0 ) then
            datestr = trim(datestr)//'00'
         elseif (min.lt.10) then
            datestr = trim(datestr)//'0'//char(ichar('0')+min)
         else
            datestr = trim(datestr)//
     >                char(ichar('0')+int(min/10))//
     >                char(ichar('0')+mod(min,10))
         endif
                    
c        Write date string
         print*,trim(datestr)

      endif

c     Write trajectories to screen
      if ( (mode.eq.'list') ) then

c        Read the complete trajectory file
         allocate(tra(ntra,ntim,ncol),stat=stat)
         call ropen_tra(fid,inpfile,ntra,ntim,ncol,refdate,vars,inpmode)
         call read_tra (fid,tra,ntra,ntim,ncol,inpmode)
         call close_tra(fid,inpmode)

c        Get the strings for output
         write(str1,'(i4)') refdate(1)
         write(str2,'(i2)') refdate(2)
         write(str3,'(i2)') refdate(3)
         write(str4,'(i2)') refdate(4)
         write(str5,'(i2)') refdate(5)
         if (refdate(2).eq. 0) str2(1:1)='0'
         if (refdate(3).eq. 0) str3(1:1)='0'
         if (refdate(4).eq. 0) str4(1:1)='0'
         if (refdate(5).eq. 0) str5(1:1)='0'
         if (refdate(2).lt.10) str2(1:1)='0'
         if (refdate(3).lt.10) str3(1:1)='0'
         if (refdate(4).lt.10) str4(1:1)='0'
         if (refdate(5).lt.10) str5(1:1)='0'

c        Write the time specification
         write(*,'(a15,a4,a2,a2,a1,a2,a2,a13,i8,a4)') 
     >          'Reference date ',
     >           str1,str2,str3,'_',str4,str5,
     >          ' / Time range',refdate(6), ' min'
         write(*,*)

c        Write variable names
         str=''
         do i=1,ncol
            str=trim(str)//trim(vars(i))
         enddo
         write(*,'(a6,a9,a8,a6,100a10)') (trim(vars(i)),i=1,ncol)
         write(*,'(a6,a9,a8,a6,100a10)') 
     >              '------','---------','--------','------',
     >              ('----------',i=5,ncol)

         do n=1,ntra
            write(*,*)
            do i=1,ntim
               write(*,'(1f7.2,f9.2,f8.2,i6,100f10.3)') 
     >               (tra(n,i,j),j=1,3),             ! time, lon, lat
     >               nint(tra(n,i,4)),               ! p
     >               (tra(n,i,j),j=5,ncol)           ! fields
            enddo
         enddo

      endif
         
c     Get number of trajectories leaving domain
      if ( (mode.eq.'leaving') ) then

c        Read the complete trajectory file
         allocate(tra(ntra,ntim,ncol),stat=stat)
         call ropen_tra(fid,inpfile,ntra,ntim,ncol,refdate,vars,inpmode)
         call read_tra (fid,tra,ntra,ntim,ncol,inpmode)
         call close_tra(fid,inpmode)

c        Determine the number of trajectories leaving domain
         do i=1,ntra
            isok = 1
            do j=1,ntim
               if ( tra(i,j,4).lt.0. ) isok = 0
            enddo         
            if ( isok.eq.0 ) then
               nleaving = nleaving + 1
            endif
         enddo

c        Write output
	 print*,nleaving     

      endif
         
      end


      

      
