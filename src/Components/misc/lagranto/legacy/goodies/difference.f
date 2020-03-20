
      PROGRAM difference
      
c     ***********************************************************************
c     * Get the difference between tow trajectory files                     *
c     * Michael Sprenger / Spring, summer, autumn 2010                      *
c     ***********************************************************************

      implicit none
      
c     ----------------------------------------------------------------------
c     Declaration of variables
c     ----------------------------------------------------------------------

c     Field name and mode for difference calculation
      character*80                           mode                     ! Difference mode
      character*80                           fieldname                ! Name of differencing

c     Input and output format for trajectories (see iotra.f)
      character*80                           inpfile1                 ! Input filename 1
      character*80                           inpfile2                 ! Input filename 2
      character*80                           outfile                  ! Output filename

c     Input trajectories
      integer                                ntra1      ,ntra2        ! Number of trajectories
      integer                                ntim1      ,ntim2        ! Number of times
      integer                                ncol1      ,ncol2        ! Number of columns
      real,allocatable, dimension (:,:,:) :: trainp1    ,trainp2      ! Trajectories (ntra,ntim,ncol)
      character*80                           vars1(100) ,vars2(100)   ! Variable names
      integer                                refdate1(6),refdate2(6)  ! Reference date
      
c     Output/comparison trajectory
      integer                                ntra                    ! Number of trajectories
      integer                                ntim                    ! Number of times
      integer                                ncol                    ! Number of columns
      real,allocatable, dimension (:,:,:) :: traout                  ! Trajectories (ntra,ntim,ncol)
      character*80                           vars(100)               ! Variable names
      integer                                refdate(6)              ! Reference date

c     Auxiliary variables
      integer                                inpmode1,inpmode2
      integer                                outmode
      integer                                stat
      integer                                fid
      integer                                i,j,k
      integer                                ind1,ind2
      integer                                isok
      character                              ch
      real,allocatable, dimension (:) ::     diff
      integer                                outind

c     Externals 
      real,external                   ::     sdis
      
c     ----------------------------------------------------------------------
c     Preparations
c     ----------------------------------------------------------------------

c     Read parameters
      open(10,file='difference.param')
       read(10,*) inpfile1
       read(10,*) inpfile2
       read(10,*) outfile
       read(10,*) ntra1,ntim1,ncol1
       read(10,*) ntra2,ntim2,ncol2
       read(10,*) mode
       read(10,*) fieldname
      close(10)
      
c     Determine the formats
      call mode_tra(inpmode1,inpfile1)
      if (inpmode1.eq.-1) inpmode1=1

      call mode_tra(inpmode2,inpfile2)
      if (inpmode2.eq.-1) inpmode2=1

      call mode_tra(outmode,outfile)
      if (outmode.eq.-1) outmode=1

c     Set number of trajectories for output
      if ( ntra1.lt.ntra2) then
         ntra = ntra1
      else
         ntra = ntra2
      endif

c     Set number of times for output
      if ( mode.eq.'single' ) then
         ntim = ntim1
      else
         ntim = 1
      endif

c     Set the column names for output 
      if ( fieldname.eq.'LATLON') then         
         ncol    = 1 + 3 + 3     + 1
         vars(1) = 'time'
         vars(2) = 'lon[1]'
         vars(3) = 'lat[1]'
         vars(4) = 'p[1]'
         vars(5) = 'lon[2]'
         vars(6) = 'lat[2]'
         vars(7) = 'p[2]'
         vars(8) = 'SDIS'
      else
         ncol = 1 + 3 + 3 + 2 + 1
         vars( 1) = 'time'
         vars( 2) = 'lon[1]'
         vars( 3) = 'lat[1]'
         vars( 4) = 'p[1]'
         vars( 5) = 'lon[2]'
         vars( 6) = 'lat[2]'
         vars( 7) = 'p[2]'
         vars( 8) = trim(fieldname)//'[1]'
         vars( 9) = trim(fieldname)//'[2]'
         vars(10) = trim(fieldname)//'[1-2]'
      endif

c     Allocate memory
      allocate(trainp1(ntra1,ntim1,ncol1),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array trainp1   ***' 
      allocate(trainp2(ntra2,ntim2,ncol2),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array trainp2   ***' 
      allocate(traout(ntra,ntim,ncol),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array traout    ***' 
      allocate(diff(ntim1),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array diff      ***' 

c     Read inpufiles
      call ropen_tra(fid,inpfile1,ntra1,ntim1,ncol1,
     >                   refdate1,vars1,inpmode1)
      call read_tra (fid,trainp1,ntra1,ntim1,ncol1,inpmode1)
      call close_tra(fid,inpmode1)

      call ropen_tra(fid,inpfile2,ntra2,ntim2,ncol2,
     >                   refdate2,vars2,inpmode2)
      call read_tra (fid,trainp2,ntra2,ntim2,ncol2,inpmode2)
      call close_tra(fid,inpmode2)

c     Check dimensions of the two trajectory files (#tim hard check)
      if (ntim1.ne.ntim2) then
         print*,'Trajectoy files have different time dimensions... Stop'
         stop
      endif

c     Check dimensions of the two trajectory files (#tra soft check)
      if (ntra1.ne.ntra2) then
         print*,'Differing number of trajectories... proceed [y/n]'
         read*,ch
         if (ch.eq.'n') stop
      endif

c     Check whether difference field is available on both files 
      if ( fieldname.ne.'LATLON') then
         ind1 = 0
         ind2 = 0
         do i=1,ncol
            if ( fieldname.eq.vars1(i) ) ind1 = i
            if ( fieldname.eq.vars2(i) ) ind2 = i
         enddo
         if ( (ind1.eq.0).or.(ind2.eq.0) ) then
            print*,'Field ',trim(fieldname),' not available... Stop'
            stop
         endif
      endif

c     Check reference dates (soft check)
      isok = 1
      do i=1,6
         if ( refdate1(i).ne.refdate2(i) ) isok = 0
      enddo
      if ( isok.eq.0 ) then
         print*,'Warning: reference dates differ... proceed [y/n]'
         read*,ch
         if (ch.eq.'n') stop
      endif

c     Check trajectory times (soft check)
      isok = 1
      do i=1,ntim
         if ( trainp1(1,i,1).ne.trainp2(1,i,1) ) isok = 0
      enddo
      if ( isok.eq.0 ) then
         print*,'Warning: trajectory times differ... proceed [y/n]'
         read*,ch
         if (ch.eq.'n') stop
      endif

c     Copy reference date to output
      do i=1,6
         refdate(i) = refdate1(i)
      enddo
      
c     ----------------------------------------------------------------------
c     Calculate the difference (depending on mode)
c     ----------------------------------------------------------------------

c     Loop over all trajectories
      do i=1,ntra

c        Calculate difference for all times
         do j=1,ntim1
            
c           Calculate the difference (distance or absolute value)
            if (fieldname.eq.'LATLON') then
               diff(j) = sdis( trainp1(i,j,2),trainp1(i,j,3),
     >                         trainp2(i,j,2),trainp2(i,j,3) )         
            else
               diff(j) = abs(trainp1(i,j,ind1) - trainp2(i,j,ind2))
            endif

         enddo

c        Save output for each time
         if ( mode.eq.'single' ) then

            do j=1,ntim
            
               if ( fieldname.eq.'LATLON' ) then
                  traout(i,j, 1) = trainp1(i,j,1)    ! time
                  traout(i,j, 2) = trainp1(i,j,2)    ! lon[1]
                  traout(i,j, 3) = trainp1(i,j,3)    ! lat[1]
                  traout(i,j, 4) = trainp1(i,j,4)    ! p[1]
                  traout(i,j, 5) = trainp2(i,j,2)    ! lon[2]
                  traout(i,j, 6) = trainp2(i,j,3)    ! lat[2]
                  traout(i,j, 7) = trainp2(i,j,4)    ! p[2]
                  traout(i,j, 8) = diff(j)           ! SDIS(j)
               else
                  traout(i,j, 1) = trainp1(i,j,1)    ! time
                  traout(i,j, 2) = trainp1(i,j,2)    ! lon[1]
                  traout(i,j, 3) = trainp1(i,j,3)    ! lat[1]
                  traout(i,j, 4) = trainp1(i,j,4)    ! p[1]
                  traout(i,j, 5) = trainp2(i,j,2)    ! lon[2]
                  traout(i,j, 6) = trainp2(i,j,3)    ! lat[2]
                  traout(i,j, 7) = trainp2(i,j,4)    ! p[2]
                  traout(i,j, 8) = trainp1(i,j,ind1) ! field[1]
                  traout(i,j, 9) = trainp2(i,j,ind2) ! field[2]
                  traout(i,j,10) = diff(j)           ! SDIS(j)
               endif

            enddo

c        Save only maximum
         elseif ( mode.eq.'max') then
            
            outind  = 1
            do j=2,ntim1
               if ( diff(j).gt.diff(outind) ) outind = j
            enddo

            if ( fieldname.eq.'LATLON' ) then
               traout(i,1, 1) = trainp1(i,outind,1)     ! time
               traout(i,1, 2) = trainp1(i,outind,2)     ! lon[1]
               traout(i,1, 3) = trainp1(i,outind,3)     ! lat[1]
               traout(i,1, 4) = trainp1(i,outind,4)     ! p[1]
               traout(i,1, 5) = trainp2(i,outind,2)     ! lon[2]
               traout(i,1, 6) = trainp2(i,outind,3)     ! lat[2]
               traout(i,1, 7) = trainp2(i,outind,4)     ! p[2]
               traout(i,1, 8) = diff(outind)            ! SDIS
            else
               traout(i,1, 1) = trainp1(i,outind,1)     ! time
               traout(i,1, 2) = trainp1(i,outind,2)     ! lon[1]
               traout(i,1, 3) = trainp1(i,outind,3)     ! lat[1]
               traout(i,1, 4) = trainp1(i,outind,4)     ! p[1]
               traout(i,1, 5) = trainp2(i,outind,2)     ! lon[2]
               traout(i,1, 6) = trainp2(i,outind,3)     ! lat[2]
               traout(i,1, 7) = trainp2(i,outind,4)     ! p[2]
               traout(i,1, 8) = trainp1(i,outind,ind1)  ! field[1]
               traout(i,1, 9) = trainp2(i,outind,ind2)  ! field[2]
               traout(i,1,10) = diff(outind)            ! SDIS(j)
            endif
            
         endif

      enddo
    
c     ----------------------------------------------------------------------
c     Write output
c     ----------------------------------------------------------------------

c     Write output file
      call wopen_tra(fid,outfile,ntra,ntim,ncol,refdate,vars,outmode)
      call write_tra(fid,traout,ntra,ntim,ncol,outmode)
      call close_tra(fid,outmode)
      
      end




c     ***********************************************************************
c     * Subroutines                                                         *
c     ***********************************************************************

c     ----------------------------------------------------------------------
c     Spherical distance
c     ----------------------------------------------------------------------

      real function sdis(xp,yp,xq,yq)
c
c     calculates spherical distance (in km) between two points given
c     by their spherical coordinates (xp,yp) and (xq,yq), respectively.
c
      real      pi180
      parameter (pi180=3.14159/180.)
      real      re
      parameter (re=6370.)
      real      degkm
      parameter (degkm=111.1775)
      real      xp,yp,xq,yq,arg

      if ( (abs(xp-xq).gt.0.05).and.(abs(yp-yq).gt.0.05) ) then
         arg=sin(pi180*yp)*sin(pi180*yq)+
     >       cos(pi180*yp)*cos(pi180*yq)*cos(pi180*(xp-xq))
         if (arg.lt.-1.) arg=-1.
         if (arg.gt.1.) arg=1.
         sdis=re*acos(arg)
      else
         sdis= (yp-yq)**2 + ( (xp-xq) * cos( pi180*0.5*(yp+yq) ) )**2
         sdis = deg2km * sqrt(sdis)
      endif

      end


      
