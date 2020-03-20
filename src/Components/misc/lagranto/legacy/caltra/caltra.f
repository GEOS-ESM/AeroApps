      PROGRAM caltra

C     ********************************************************************
C     *                                                                  *
C     * Calculates trajectories                                          *
C     *                                                                  *
C     *	Heini Wernli	   first version:	April 1993               *
C     * Michael Sprenger   major upgrade:       2008-2009                *
C     *                                                                  *
C     ********************************************************************

      implicit none      

c     --------------------------------------------------------------------
c     Declaration of parameters
c     --------------------------------------------------------------------

c     Maximum number of levels for input files
      integer   nlevmax
      parameter	(nlevmax=200)

c     Maximum number of input files (dates, length of trajectories)
      integer   ndatmax
      parameter	(ndatmax=500)

c     Numerical epsilon (for float comparison)
      real      eps
      parameter (eps=0.001)

c     Distance in m between 2 lat circles 
      real	deltay
      parameter	(deltay=1.112E5)
      
c     Numerical method for the integration (0=iterative Euler, 1=Runge-Kutta)
      integer   imethod
      parameter (imethod=1)

c     Number of iterations for iterative Euler scheme
      integer   numit
      parameter (numit=3)

c     Input and output format for trajectories (see iotra.f)
      integer   inpmode
      integer   outmode

c     Filename prefix for primary and secondary files (typically 'P' and 'S')
      character*1 prefix
      parameter   (prefix='P')
      character*1 srefix
      parameter   (srefix='S')

c     Physical constants - needed to compute potential temperature
      real      rdcp,tzero
      data      rdcp,tzero /0.286,273.15/

c     --------------------------------------------------------------------
c     Declaration of variables
c     --------------------------------------------------------------------

c     Input parameters
      integer                                fbflag          ! Flag for forward/backward mode
      integer                                numdat          ! Number of input files
      character*11                           dat(ndatmax)    ! Dates of input files
      real                                   timeinc         ! Time increment between input files
      real                                   per             ! Periodicity (=0 if none)
      integer                                ntra            ! Number of trajectories
      character*80                           cdfname         ! Name of output files
      real                                   ts              ! Time step
      real                                   tst,ten         ! Shift of start and end time relative to first data file
      integer                                deltout         ! Output time interval (in minutes)
      integer                                jflag           ! Jump flag (if =1 ground-touching trajectories reenter atmosphere)
      real                                   wfactor         ! Factor for vertical velocity field
      character*80                           strname         ! File with start positions
      character*80                           timecheck       ! Either 'yes' or 'no'
      character*80                           isen            ! Isentropic trajectories ('yes' or 'no')
      integer                                thons           ! Isentropic mode: is TH availanle on S
      character*80                           modlev          ! 2D (model level) trajectories ('yes' or 'no')

c     Trajectories
      integer                                ncol            ! Number of columns for insput trajectories
      real,allocatable, dimension (:,:,:) :: trainp          ! Input start coordinates (ntra,1,ncol)
      real,allocatable, dimension (:,:,:) :: traout          ! Output trajectories (ntra,ntim,4)
      integer                                reftime(6)      ! Reference date
      character*80                           vars(200)       ! Field names
      real,allocatable, dimension (:)     :: xx0,yy0,pp0     ! Position of air parcels
      integer,allocatable, dimension (:)  :: leftflag        ! Flag for domain-leaving
      real                                   xx1,yy1,pp1     ! Updated position of air parcel
      integer                                leftcount       ! Number of domain leaving trajectories
      integer                                ntim            ! Number of output time steps
      real,allocatable, dimension (:)     :: theta           ! Potential temperature for isentropic trajectories
      real,allocatable, dimension (:)     :: zindex          ! Vertical index for model-level (2D) trajectories

c     Meteorological fields
      real,allocatable, dimension (:)     :: spt0,spt1       ! Surface pressure
      real,allocatable, dimension (:)     :: uut0,uut1       ! Zonal wind
      real,allocatable, dimension (:)     :: vvt0,vvt1       ! Meridional wind
      real,allocatable, dimension (:)     :: wwt0,wwt1       ! Vertical wind
      real,allocatable, dimension (:)     :: p3t0,p3t1       ! 3d-pressure 
      real,allocatable, dimension (:)     :: tht0,tht1       ! 3d potential temperature
      real,allocatable, dimension (:)     :: sth0,sth1       ! Surface potential temperature

c     Grid description
      real                                   pollon,pollat   ! Longitude/latitude of pole
      real                                   ak(nlevmax)     ! Vertical layers and levels
      real                                   bk(nlevmax) 
      real                                   xmin,xmax       ! Zonal grid extension
      real                                   ymin,ymax       ! Meridional grid extension
      integer                                nx,ny,nz        ! Grid dimensions
      real                                   dx,dy           ! Horizontal grid resolution
      integer                                hem             ! Flag for hemispheric domain
      real                                   mdv             ! Missing data value

c     Auxiliary variables                 
      real                                   delta,rd
      integer	                             itm,iloop,i,j,k,filo,lalo
      integer                                ierr,stat
      integer                                cdfid,fid
      real	                                 tstart,time0,time1,time
      real                                   reltpos0,reltpos1
      real                                   xind,yind,pind,pp,sp,stagz
      character*80                           filename,varname
      integer                                reftmp(6)
      character                              ch
      real                                   frac,tload
      integer                                itim
      integer                                wstep
      real                                   x1,y1,p1
      real                                   thetamin,thetamax
      real                                   zindexmin,zindexmax

c     Externals
      real                                   int_index4
      external                               int_index4

c     --------------------------------------------------------------------
c     Start of program, Read parameters
c     --------------------------------------------------------------------

c     Write start message
      print*,'========================================================='
      print*,'              *** START OF PROGRAM CALTRA ***'
      print*

c     Open the parameter file
      open(9,file='caltra.param')

c     Read flag for forward/backward mode (fbflag)
      read(9,*) fbflag

c     Read number of input files (numdat)
      read(9,*) numdat
      if (numdat.gt.ndatmax) then
        print*,' ERROR: too many input files ',numdat,ndatmax
        goto 993
      endif

c     Read list of input dates (dat, sort depending on forward/backward mode)
      if (fbflag.eq.1) then
        do itm=1,numdat
          read(9,'(a11)') dat(itm)
        enddo
      else
        do itm=numdat,1,-1
          read(9,'(a11)') dat(itm)
        enddo
      endif

c     Read time increment between input files (timeinc)
      read(9,*) timeinc

C     Read if data domain is periodic and its periodicity
      read(9,*) per

c     Read the number of trajectories and name of position file
      read(9,*) strname
      read(9,*) ntra
      read(9,*) ncol 
      if (ntra.eq.0) goto 991

C     Read the name of the output trajectory file and set the constants file
      read(9,*) cdfname

C     Read the timestep for trajectory calculation (convert from minutes to hours)
      read(9,*) ts
      ts=ts/60.        

C     Read shift of start and end time relative to first data file
      read(9,*) tst
      read(9,*) ten

C     Read output time interval (in minutes)
      read(9,*) deltout

C     Read jumpflag (if =1 ground-touching trajectories reenter the atmosphere)
      read(9,*) jflag

C     Read factor for vertical velocity field
      read(9,*) wfactor

c     Read the reference time and the time range
      read(9,*) reftime(1)                  ! year
      read(9,*) reftime(2)                  ! month
      read(9,*) reftime(3)                  ! day
      read(9,*) reftime(4)                  ! hour
      read(9,*) reftime(5)                  ! min
      read(9,*) reftime(6)                  ! time range (in min)

c     Read flag for 'no time check'
      read(9,*) timecheck

c     Read flag for isentropic trajectories
      read(9,*) isen, thons

c     Read flag for model-level trajectories (2D mode)
      read(9,*) modlev

c     Close the input file
      close(9)

c     Calculate the number of output time steps
      ntim = abs(reftime(6)/deltout) + 1

c     Set the formats of the input and output files
      call mode_tra(inpmode,strname)
      call mode_tra(outmode,cdfname)
      if (outmode.eq.-1) outmode=1

c     Write some status information
      print*,'---- INPUT PARAMETERS -----------------------------------'
      print* 
      print*,'  Forward/Backward       : ',fbflag
      print*,'  #input files           : ',numdat
      print*,'  First/last input file  : ',trim(dat(1)),' ... ',
     >                                     trim(dat(numdat))
      print*,'  time increment         : ',timeinc
      print*,'  Output file            : ',trim(cdfname)
      print*,'  Time step (min)        : ',60.*ts
      write(*,'(a27,f7.2,f7.2)') '   Time shift (start,end) : ',tst,ten
      print*,'  Output time interval   : ',deltout
      print*,'  Jump flag              : ',jflag
      print*,'  Vertical wind (scale)  : ',wfactor
      print*,'  Trajectory pos file    : ',trim(strname)
      print*,'  # of trajectories      : ',ntra
      print*,'  # of output timesteps  : ',ntim
      if ( inpmode.eq.-1) then
         print*,'  Input format           : (lon,lat,p)-list'
      else
         print*,'  Input format           : ',inpmode
      endif
      print*,'  Output format          : ',outmode
      print*,'  Periodicity            : ',per
      print*,'  Time check             : ',trim(timecheck)
      print*,'  Isentropic trajectories: ',trim(isen),thons
      print*,'  Model-level trajs (2D) : ',trim(modlev)
      print*

      if ( (isen.eq.'yes').and.(modlev.eq.'yes') ) then
         print*,
     >    '  WARNING: isentropic and 2D mode chosen -> 2D accepted'
         print*
         isen = 'no'
      endif

c     Init missing data value
      mdv = -999.

      print*,'---- FIXED NUMERICAL PARAMETERS -------------------------'
      print*
      print*,'  Numerical scheme       : ',imethod
      print*,'  Number of iterations   : ',numit
      print*,'  Filename prefix        : ',prefix
      print*,'  Missing data value     : ',mdv
      print*

c     --------------------------------------------------------------------
c     Read grid parameters, checks and allocate memory
c     --------------------------------------------------------------------

c     Read the constant grid parameters (nx,ny,nz,xmin,xmax,ymin,ymax,pollon,pollat)
c     The negative <-fid> of the file identifier is used as a flag for parameter retrieval  
      filename = prefix//dat(1)
      varname  = 'U'
      nx       = 1
      ny       = 1
      nz       = 1
      tload    = -tst
      call input_open (fid,filename)
      call input_grid (-fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >                 tload,pollon,pollat,rd,rd,nz,rd,rd,rd,timecheck)
      call input_close(fid)

C     Check if the number of levels is too large
      if (nz.gt.nlevmax) goto 993

C     Set logical flag for periodic data set (hemispheric or not)
      hem = 0
      if (per.eq.0.) then
         delta=xmax-xmin-360.
         if (abs(delta+dx).lt.eps) then               ! Program aborts: arrays must be closed
            goto 992
        else if (abs(delta).lt.eps) then              ! Periodic and hemispheric
           hem=1
           per=360.
        endif
      else                                            ! Periodic and hemispheric
         hem=1
      endif
      
C     Allocate memory for some meteorological arrays
      allocate(spt0(nx*ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array spt0 ***'   ! Surface pressure
      allocate(spt1(nx*ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array spt1 ***'
      allocate(uut0(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array uut0 ***'   ! Zonal wind
      allocate(uut1(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array uut1 ***'
      allocate(vvt0(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array vvt0 ***'   ! Meridional wind
      allocate(vvt1(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array vvt1 ***'
      allocate(wwt0(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array wwt0 ***'   ! Vertical wind
      allocate(wwt1(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array wwt1 ***'
      allocate(p3t0(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array p3t0 ***'   ! Pressure
      allocate(p3t1(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array p3t1 ***'
      allocate(tht0(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tht0 ***'   ! Potential temperature
      allocate(tht1(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tht1 ***'
      allocate(sth0(nx*ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array spt0 ***'   ! Surface potential temperature
      allocate(sth1(nx*ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array sth1 ***'

C     Get memory for trajectory arrays
      allocate(trainp(ntra,1,ncol),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array trainp   ***' ! Input start coordinates
      allocate(traout(ntra,ntim,4),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array traout   ***' ! Output trajectories
      allocate(xx0(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array xx0      ***' ! X position (longitude)
      allocate(yy0(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array yy0      ***' ! Y position (latitude)
      allocate(pp0(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array pp0      ***' ! Pressure
      allocate(leftflag(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array leftflag ***' ! Leaving-domain flag
      allocate(theta(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array theta ***'    ! Potential temperature for isentropic trajectories
      allocate(zindex(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array kind ***'     ! Vertical index for model-level trajectories

c     Write some status information
      print*,'---- CONSTANT GRID PARAMETERS ---------------------------'
      print*
      print*,'  xmin,xmax     : ',xmin,xmax
      print*,'  ymin,ymax     : ',ymin,ymax
      print*,'  dx,dy         : ',dx,dy
      print*,'  pollon,pollat : ',pollon,pollat
      print*,'  nx,ny,nz      : ',nx,ny,nz
      print*,'  per, hem      : ',per,hem
      print*

c     --------------------------------------------------------------------
c     Initialize the trajectory calculation
c     --------------------------------------------------------------------

c     Read start coordinates from file - Format (lon,lat,lev)
      if (inpmode.eq.-1) then
         open(fid,file=strname)
          do i=1,ntra
             read(fid,*) xx0(i),yy0(i),pp0(i)
          enddo
         close(fid)

c     Read start coordinates from trajectory file - check consistency of ref time
      else
         call ropen_tra(cdfid,strname,ntra,1,ncol,reftmp,vars,inpmode)
         call read_tra (cdfid,trainp,ntra,1,ncol,inpmode)
         do i=1,ntra
            time   = trainp(i,1,1)
            xx0(i) = trainp(i,1,2) 
            yy0(i) = trainp(i,1,3) 
            pp0(i) = trainp(i,1,4) 
         enddo
         call close_tra(cdfid,inpmode)

         if ( ( reftime(1).ne.reftmp(1) ).or.
     >        ( reftime(2).ne.reftmp(2) ).or.
     >        ( reftime(3).ne.reftmp(3) ).or.
     >        ( reftime(4).ne.reftmp(4) ).or.
     >        ( reftime(5).ne.reftmp(5) ) )
     >   then
            print*,' WARNING: Inconsistent reference times'
            write(*,'(5i8)') (reftime(i),i=1,5)
            write(*,'(5i8)') (reftmp (i),i=1,5)
            print*,'Enter a key to proceed...'
            stop
         endif
      endif

c     Set sign of time range
      reftime(6) = fbflag * reftime(6)
         
c     Write some status information
      print*,'---- REFERENCE DATE---------- ---------------------------'
      print*
      print*,' Reference time (year)  :',reftime(1)
      print*,'                (month) :',reftime(2)
      print*,'                (day)   :',reftime(3)
      print*,'                (hour)  :',reftime(4)
      print*,'                (min)   :',reftime(5)
      print*,' Time range             :',reftime(6),' min'
      print*

C     Save starting positions 
      itim = 1
      do i=1,ntra
         traout(i,itim,1) = 0.
         traout(i,itim,2) = xx0(i)
         traout(i,itim,3) = yy0(i)
         traout(i,itim,4) = pp0(i)
      enddo

c     Init the flag and the counter for trajectories leaving the domain
      leftcount=0
      do i=1,ntra
         leftflag(i)=0
      enddo

C     Convert time shifts <tst,ten> from <hh.mm> into fractional time
      call hhmm2frac(tst,frac)
      tst = frac
      call hhmm2frac(ten,frac)
      ten = frac

c     Get 3D and surface pressure from first data file 
      filename = prefix//dat(1)
      call input_open (fid,filename)
      if ( modlev.eq.'no' ) then 
         varname = 'P'
         call input_grid
     >       (fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >        tload,pollon,pollat,p3t1,spt1,nz,ak,bk,stagz,timecheck)
      else  
         varname = 'P.ML'
         call input_grid
     >       (fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >        tload,pollon,pollat,p3t1,spt1,nz,ak,bk,stagz,timecheck)
      endif
      call input_close(fid)
      
c     Check that all starting positions are above topography    
      do i=1,ntra

C       Interpolate surface pressure to actual position (from first input file)
        x1 = xx0(i)
        y1 = yy0(i)
        call get_index4 (xind,yind,pind,x1,y1,1050.,0.,
     >                   p3t1,p3t1,spt1,spt1,3,
     >                   nx,ny,nz,xmin,ymin,dx,dy,mdv)
        sp = int_index4 (spt1,spt1,nx,ny,1,xind,yind,1.,0.,mdv)

c       Decide whether to keep the trajectory
        if ( pp0(i).gt.sp ) then
            write(*,'(a30,3f10.2)')
     >               'WARNING: starting point below topography ',
     >               xx0(i),yy0(i),pp0(i)
            leftflag(i) = 1
        endif

      enddo

c     Special handling for isentropic trajectories - read potential
c     temperature from S file or calculate it based on temperature and
c     pressure from P file; then, calculate for each trajectory its
c     potential temperature - which will stay fixed over time
      if ( isen.eq.'yes' ) then

c         Get potential temperature from S file
          if ( thons.eq.1 ) then
              filename = srefix//dat(1)
              print*,' TH <- ',trim(filename)
              call input_open (fid,filename)
              varname='TH'
              call input_wind
     >            (fid,varname,tht1,tload,stagz,mdv,
     >              xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)
              call input_close(fid)

c         Calculate potential temperature from P file
          else
              filename = prefix//dat(1)
              print*,' TH = T * (1000/P)^RDCP <- ',trim(filename)
              call input_open (fid,filename)
              varname='T'
              call input_wind
     >            (fid,varname,tht1,tload,stagz,mdv,
     >              xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)
              call input_close(fid)
              do i=1,nx*ny*nz
                if (tht1(i).lt.100.) then
                   tht1(i)=(tht1(i)+tzero)*( (1000./p3t1(i))**rdcp )
                else
                   tht1(i)=tht1(i)*( (1000./p3t1(i))**rdcp )
                endif
              enddo
          endif

c         Take surface potential temperature from lowest level
          do i=1,nx*ny
            sth1(i) = tht1(i)
          enddo

c         Calculate now the potential temperature of all trajectories
          do i=1,ntra

             x1 = xx0(i)
             y1 = yy0(i)
             p1 = pp0(i)
             call get_index4 (xind,yind,pind,x1,y1,p1,0.,
     >                 p3t1,p3t1,spt1,spt1,3,
     >                 nx,ny,nz,xmin,ymin,dx,dy,mdv)
             theta(i) =
     >        int_index4(tht1,tht1,nx,ny,nz,xind,yind,pind,0.,mdv)

          enddo

c         Write info about theta range of starting positions
          thetamin = theta(1)
          thetamax = theta(1)
          do i=2,ntra
            if ( theta(i).gt.thetamax ) thetamax = theta(i)
            if ( theta(i).lt.thetamin ) thetamin = theta(i)
          enddo

c         Write some status information
          print*
          print*,
     >     '---- THETA RANGE OF ISENTROPIC TRAJECTORIES -------------'
          print*
          print*,' Theta(min)             :',thetamin
          print*,' Theta(max)             :',thetamax

      endif

c     Special handling for model-level (2D) trajectories - get the
c     vertical index for each trajectory - which will remain fixed
      if ( modlev.eq.'yes' ) then
          do i=1,ntra
             x1 = xx0(i)
             y1 = yy0(i)
             p1 = pp0(i)
             call get_index4 (xind,yind,pind,x1,y1,p1,0.,
     >                 p3t1,p3t1,spt1,spt1,3,
     >                 nx,ny,nz,xmin,ymin,dx,dy,mdv)
             zindex(i) = pind
          enddo

          do i=1,nz
            print*,i,p3t1(189+(141-1)*nx+(i-1)*nx*ny)
          enddo
             print*,x1,y1,p1
             print*,xind,yind,pind

c         Write info about zindex range of starting positions
          zindexmin = zindex(1)
          zindexmax = zindex(1)
          do i=2,ntra
            if ( zindex(i).gt.zindexmax ) zindexmax = zindex(i)
            if ( zindex(i).lt.zindexmin ) zindexmin = zindex(i)
          enddo

c         Write some status information
          print*
          print*,
     >     '---- INDEX RANGE OF MODEL-LEVEL TRAJECTORIES ------------'
          print*
          print*,' Zindex(min)            :',zindexmin
          print*,' Zindex(max)            :',zindexmax


      endif

c     -----------------------------------------------------------------------
c     Loop to calculate trajectories
c     -----------------------------------------------------------------------

c     Write some status information
      print*
      print*,'---- TRAJECTORIES ----------- ---------------------------'
      print*    

C     Set the time for the first data file (depending on forward/backward mode)
      if (fbflag.eq.1) then
        tstart = -tst
      else
        tstart = tst
      endif

c     Set the minute counter for output
      wstep = 0

c     Read wind fields and vertical grid from first file
      filename = prefix//dat(1)

      call frac2hhmm(tstart,tload)

      write(*,'(a16,a20,f9.2)') '  (file,time) : ',
     >                       trim(filename),tload

      call input_open (fid,filename)

      varname='U'                                      ! U
      call input_wind 
     >    (fid,varname,uut1,tload,stagz,mdv,
     >     xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)

      varname='V'                                      ! V
      call input_wind 
     >    (fid,varname,vvt1,tload,stagz,mdv,
     >     xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)

      if ( (modlev.eq.'no').and.(isen.eq.'no') ) then
         varname='OMEGA'                                  ! OMEGA
         call input_wind
     >       (fid,varname,wwt1,tload,stagz,mdv,
     >        xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)
      endif

      if ( modlev.eq.'no' ) then
         call input_grid                                  ! GRID - AK,BK -> P
     >       (fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >        tload,pollon,pollat,p3t1,spt1,nz,ak,bk,stagz,timecheck)
      else
         varname='P.ML'                                   ! GRID - P,PS
         call input_grid                                  !
     >       (fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >        tload,pollon,pollat,p3t1,spt1,nz,ak,bk,stagz,timecheck)
      endif

      call input_close(fid)

c     Special handling for isentropic trajectories - read potential
c     temperature from S file or calculate it based on temperature and
c     pressure from P file
      if ( isen.eq.'yes' ) then

c         Get potential temperature from S file
          if ( thons.eq.1 ) then
              filename = srefix//dat(1)
              print*,' TH <- ',trim(filename)
              call input_open (fid,filename)
              varname='TH'
              call input_wind
     >            (fid,varname,tht1,tload,stagz,mdv,
     >              xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)
              call input_close(fid)

c         Calculate potential temperature from P file
          else
              filename = prefix//dat(1)
              print*,' TH = T * (1000/P)^RDCP <- ',trim(filename)
              call input_open (fid,filename)
              varname='T'
              call input_wind
     >            (fid,varname,tht1,tload,stagz,mdv,
     >              xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)
              call input_close(fid)
              do i=1,nx*ny*nz
                if (tht1(i).lt.100.) then
                   tht1(i)=(tht1(i)+tzero)*( (1000./p3t1(i))**rdcp )
                else
                   tht1(i)=tht1(i)*( (1000./p3t1(i))**rdcp )
                endif
              enddo
          endif

c         Take surface potential temperature from lowest level
          do i=1,nx*ny
            sth1(i) = tht1(i)
          enddo
      endif

c     Loop over all input files (time step is <timeinc>)
      do itm=1,numdat-1

c       Calculate actual and next time
        time0 = tstart+real(itm-1)*timeinc*fbflag
        time1 = time0+timeinc*fbflag

c       Copy old velocities and pressure fields to new ones       
        do i=1,nx*ny*nz
           uut0(i)=uut1(i)
           vvt0(i)=vvt1(i)
           wwt0(i)=wwt1(i)
           p3t0(i)=p3t1(i)
           tht0(i)=tht1(i)
        enddo
        do i=1,nx*ny
           spt0(i)=spt1(i)
           sth0(i)=sth1(i)
        enddo

c       Read wind fields and surface pressure at next time
        filename = prefix//dat(itm+1)

        call frac2hhmm(time1,tload)
        write(*,'(a16,a20,f9.2)') '  (file,time) : ',
     >                          trim(filename),tload

        call input_open (fid,filename)

        varname='U'                                     ! U
        call input_wind 
     >       (fid,varname,uut1,tload,stagz,mdv,
     >        xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)

        varname='V'                                     ! V
        call input_wind 
     >       (fid,varname,vvt1,tload,stagz,mdv,
     >        xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)

        if ( (modlev.eq.'no').and.(isen.eq.'no') ) then
           varname='OMEGA'                              ! OMEGA
           call input_wind
     >          (fid,varname,wwt1,tload,stagz,mdv,
     >           xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)
        endif

        if ( modlev.eq.'no' ) then
           call input_grid                                  ! GRID - AK,NK -> P
     >          (fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >           tload,pollon,pollat,p3t1,spt1,nz,ak,bk,stagz,timecheck)
        else
           varname='P.ML'                                   ! GRID - P,PS
           call input_grid                                  !
     >       (fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >        tload,pollon,pollat,p3t1,spt1,nz,ak,bk,stagz,timecheck)
        endif

        call input_close(fid)

c     Special handling for isentropic trajectories - read potential
c     temperature from S file or calculate it based on temperature and
c     pressure from P file
      if ( isen.eq.'yes' ) then

c         Get TH from S file
          if ( thons.eq.1 ) then
              filename = srefix//dat(itm+1)
              print*,' TH <- ',trim(filename)
              call input_open (fid,filename)
              varname='TH'
              call input_wind
     >            (fid,varname,tht1,tload,stagz,mdv,
     >              xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)
              call input_close(fid)

c         Calculate potential temperature from P file
          else
              filename = prefix//dat(itm+1)
              print*,' TH = T * (1000/P)^RDCP <- ',trim(filename)
              call input_open (fid,filename)
              varname='T'
              call input_wind
     >            (fid,varname,tht1,tload,stagz,mdv,
     >              xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)
              call input_close(fid)
              do i=1,nx*ny*nz
                if (tht1(i).lt.100.) then
                   tht1(i)=(tht1(i)+tzero)*( (1000./p3t1(i))**rdcp )
                else
                   tht1(i)=tht1(i)*( (1000./p3t1(i))**rdcp )
                endif
              enddo
          endif

c         Take surface potential temperature from lowest level
          do i=1,nx*ny
            sth1(i) = tht1(i)
          enddo
      endif
        
C       Determine the first and last loop indices
        if (numdat.eq.2) then
          filo = nint(tst/ts)+1
          lalo = nint((timeinc-ten)/ts) 
        elseif ( itm.eq.1 ) then
          filo = nint(tst/ts)+1
          lalo = nint(timeinc/ts)
        else if (itm.eq.numdat-1) then
          filo = 1
          lalo = nint((timeinc-ten)/ts)
        else
          filo = 1
          lalo = nint(timeinc/ts)
        endif

c       Split the interval <timeinc> into computational time steps <ts>
        do iloop=filo,lalo

C         Calculate relative time position in the interval timeinc (0=beginning, 1=end)
          reltpos0 = ((real(iloop)-1.)*ts)/timeinc
          reltpos1 = real(iloop)*ts/timeinc

c         Timestep for all trajectories
          do i=1,ntra

C           Check if trajectory has already left the data domain
            if (leftflag(i).ne.1) then	

c             3D: Iterative Euler timestep (x0,y0,p0 -> x1,y1,p1)
              if ( (imethod.eq.1   ).and.
     >             (isen   .eq.'no').and.
     >             (modlev .eq.'no') )
     >        then
                 call euler_3d(
     >               xx1,yy1,pp1,leftflag(i),
     >               xx0(i),yy0(i),pp0(i),reltpos0,reltpos1,
     >               ts*3600,numit,jflag,mdv,wfactor,fbflag,
     >               spt0,spt1,p3t0,p3t1,uut0,uut1,vvt0,vvt1,wwt0,wwt1,
     >               xmin,ymin,dx,dy,per,hem,nx,ny,nz)

c             3D: Runge-Kutta timestep (x0,y0,p0 -> x1,y1,p1)
              else if ( (imethod.eq.2   ).and.
     >                  (isen   .eq.'no').and.
     >                  (modlev .eq.'no') )
     >        then
                 call runge(
     >               xx1,yy1,pp1,leftflag(i),
     >               xx0(i),yy0(i),pp0(i),reltpos0,reltpos1,
     >               ts*3600,numit,jflag,mdv,wfactor,fbflag,
     >               spt0,spt1,p3t0,p3t1,uut0,uut1,vvt0,vvt1,wwt0,wwt1,
     >               xmin,ymin,dx,dy,per,hem,nx,ny,nz)

c             ISENTROPIC: Iterative Euler timestep (x0,y0,p0 -> x1,y1,p1)
              else if ( (imethod.eq.1    ).and.
     >                  (isen   .eq.'yes').and.
     >                  (modlev .eq.'no' ) )
     >        then
                 call euler_isen(
     >               xx1,yy1,pp1,leftflag(i),
     >               xx0(i),yy0(i),pp0(i),theta(i),reltpos0,reltpos1,
     >               ts*3600,numit,jflag,mdv,wfactor,fbflag,
     >               spt0,spt1,p3t0,p3t1,uut0,uut1,vvt0,vvt1,
     >               sth0,sth1,tht0,tht1,
     >               xmin,ymin,dx,dy,per,hem,nx,ny,nz)

c             MODEL-LEVEL (2D): Iterative Euler timestep (x0,y0,p0 -> x1,y1,p1)
              else if ( (imethod.eq.1    ).and.
     >                  (isen   .eq.'no' ).and.
     >                  (modlev .eq.'yes') )
     >        then
                 call euler_2d(
     >               xx1,yy1,pp1,leftflag(i),
     >               xx0(i),yy0(i),pp0(i),zindex(i),reltpos0,reltpos1,
     >               ts*3600,numit,jflag,mdv,wfactor,fbflag,
     >               spt0,spt1,p3t0,p3t1,uut0,uut1,vvt0,vvt1,
     >               xmin,ymin,dx,dy,per,hem,nx,ny,nz)

              endif

c             Update trajectory position, or increase number of trajectories leaving domain
              if (leftflag(i).eq.1) then
                leftcount=leftcount+1
                if ( leftcount.lt.10 ) then
                   print*,'     -> Trajectory ',i,' leaves domain'
                elseif ( leftcount.eq.10 ) then
                   print*,'     -> N>=10 trajectories leave domain'
                endif
              else
                xx0(i)=xx1      
                yy0(i)=yy1
                pp0(i)=pp1
              endif

c          Trajectory has already left data domain (mark as <mdv>)
           else
              xx0(i)=mdv
              yy0(i)=mdv
              pp0(i)=mdv
              
           endif

          enddo

C         Save positions only every deltout minutes
          delta = aint(iloop*60*ts/deltout)-iloop*60*ts/deltout
          if (abs(delta).lt.eps) then
c          wstep = wstep + abs(ts)
c          if ( mod(wstep,deltout).eq.0 ) then
            time = time0+reltpos1*timeinc*fbflag
            itim = itim + 1
            if ( itim.le.ntim ) then
              do i=1,ntra
                 call frac2hhmm(time,tload)
                 traout(i,itim,1) = tload
                 traout(i,itim,2) = xx0(i)
                 traout(i,itim,3) = yy0(i)
                 traout(i,itim,4) = pp0(i)
              enddo
            endif
          endif

        enddo

      enddo

c     Write trajectory file
      vars(1)  ='time'
      vars(2)  ='lon'
      vars(3)  ='lat'
      vars(4)  ='p'
      call wopen_tra(cdfid,cdfname,ntra,ntim,4,reftime,vars,outmode)
      call write_tra(cdfid,traout,ntra,ntim,4,outmode)
      call close_tra(cdfid,outmode)   

c     Write some status information, and end of program message
      print*  
      print*,'---- STATUS INFORMATION --------------------------------'
      print*
      print*,'  #leaving domain    ', leftcount
      print*,'  #staying in domain ', ntra-leftcount
      print*
      print*,'              *** END OF PROGRAM CALTRA ***'
      print*,'========================================================='

      stop

c     ------------------------------------------------------------------
c     Exception handling
c     ------------------------------------------------------------------

 991  write(*,*) '*** ERROR: all start points outside the data domain'
      call exit(1)
      
 992  write(*,*) '*** ERROR: close arrays on files (prog. closear)'
      call exit(1)

 993  write(*,*) '*** ERROR: problems with array size'
      call exit(1)

      end 


c     *******************************************************************
c     * Time step : either Euler or Runge-Kutta                         *
c     *******************************************************************

C     Time-step from (x0,y0,p0) to (x1,y1,p1)
C
C     (x0,y0,p0) input	coordinates (long,lat,p) for starting point
C     (x1,y1,p1) output	coordinates (long,lat,p) for end point
C     deltat	 input	timestep in seconds
C     numit	 input	number of iterations
C     jump	 input  flag (=1 trajectories don't enter the ground)
C     left	 output	flag (=1 if trajectory leaves data domain)

c     -------------------------------------------------------------------
c     Iterative Euler time step (KINEMATIC 3D TRAJECTORIES)
c     -------------------------------------------------------------------

      subroutine euler_3d(x1,y1,p1,left,x0,y0,p0,reltpos0,reltpos1,
     >                deltat,numit,jump,mdv,wfactor,fbflag,
     >		          spt0,spt1,p3d0,p3d1,uut0,uut1,vvt0,vvt1,wwt0,wwt1,
     >                xmin,ymin,dx,dy,per,hem,nx,ny,nz)

      implicit none

c     Declaration of subroutine parameters
      integer      nx,ny,nz
      real         x1,y1,p1
      integer      left
      real	   x0,y0,p0
      real         reltpos0,reltpos1
      real   	   deltat
      integer      numit
      integer      jump
      real         wfactor
      integer      fbflag
      real     	   spt0(nx*ny)   ,spt1(nx*ny)
      real         uut0(nx*ny*nz),uut1(nx*ny*nz)
      real 	   vvt0(nx*ny*nz),vvt1(nx*ny*nz)
      real         wwt0(nx*ny*nz),wwt1(nx*ny*nz)
      real         p3d0(nx*ny*nz),p3d1(nx*ny*nz)
      real         xmin,ymin,dx,dy
      real         per
      integer      hem
      real         mdv

c     Numerical and physical constants
      real         deltay
      parameter    (deltay=1.112E5)  ! Distance in m between 2 lat circles
      real         pi                       
      parameter    (pi=3.1415927)    ! Pi

c     Auxiliary variables
      real         xmax,ymax
      real	   xind,yind,pind
      real	   u0,v0,w0,u1,v1,w1,u,v,w,sp
      integer	   icount
      character    ch

c     Externals    
      real         int_index4
      external     int_index4

c     Reset the flag for domain-leaving
      left=0

c     Set the esat-north bounray of the domain
      xmax = xmin+real(nx-1)*dx
      ymax = ymin+real(ny-1)*dy

C     Interpolate wind fields to starting position (x0,y0,p0)
      call get_index4 (xind,yind,pind,x0,y0,p0,reltpos0,
     >                 p3d0,p3d1,spt0,spt1,3,
     >                 nx,ny,nz,xmin,ymin,dx,dy,mdv)
      u0 = int_index4(uut0,uut1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)
      v0 = int_index4(vvt0,vvt1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)
      w0 = int_index4(wwt0,wwt1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)

c     Force the near-surface wind to zero
      if (pind.lt.1.) w0=w0*pind

C     For first iteration take ending position equal to starting position
      x1=x0
      y1=y0
      p1=p0

C     Iterative calculation of new position
      do icount=1,numit

C        Calculate new winds for advection
         call get_index4 (xind,yind,pind,x1,y1,p1,reltpos1,
     >                    p3d0,p3d1,spt0,spt1,3,
     >                    nx,ny,nz,xmin,ymin,dx,dy,mdv)
         u1 = int_index4(uut0,uut1,nx,ny,nz,xind,yind,pind,reltpos1,mdv)
         v1 = int_index4(vvt0,vvt1,nx,ny,nz,xind,yind,pind,reltpos1,mdv)
         w1 = int_index4(wwt0,wwt1,nx,ny,nz,xind,yind,pind,reltpos1,mdv)

c        Force the near-surface wind to zero
         if (pind.lt.1.) w1=w1*pind
 
c        Get the new velocity in between
         u=(u0+u1)/2.
         v=(v0+v1)/2.
         w=(w0+w1)/2.
         
C        Calculate new positions
         x1 = x0 + fbflag*u*deltat/(deltay*cos(y0*pi/180.))
         y1 = y0 + fbflag*v*deltat/deltay
         p1 = p0 + fbflag*wfactor*w*deltat/100.

c       Handle pole problems (crossing and near pole trajectory)
        if ((hem.eq.1).and.(y1.gt.90.)) then
          y1=180.-y1
          x1=x1+per/2.
        endif
        if ((hem.eq.1).and.(y1.lt.-90.)) then
          y1=-180.-y1
          x1=x1+per/2.
        endif
        if (y1.gt.89.99) then
           y1=89.99
        endif

c       Handle crossings of the dateline
        if ((hem.eq.1).and.(x1.gt.xmin+per-dx)) then
           x1=xmin+amod(x1-xmin,per)
        endif
        if ((hem.eq.1).and.(x1.lt.xmin)) then
           x1=xmin+per+amod(x1-xmin,per)
        endif

C       Interpolate surface pressure to actual position
        call get_index4 (xind,yind,pind,x1,y1,1050.,reltpos1,
     >                   p3d0,p3d1,spt0,spt1,3,
     >                   nx,ny,nz,xmin,ymin,dx,dy,mdv)
        sp = int_index4 (spt0,spt1,nx,ny,1,xind,yind,1.,reltpos1,mdv)

c       Handle trajectories which cross the lower boundary (jump flag)
        if ((jump.eq.1).and.(p1.gt.sp)) p1=sp-10.
 
C       Check if trajectory leaves data domain
        if ( ( (hem.eq.0).and.(x1.lt.xmin)    ).or.
     >       ( (hem.eq.0).and.(x1.gt.xmax-dx) ).or.
     >         (y1.lt.ymin).or.(y1.gt.ymax).or.(p1.gt.sp) )
     >  then
          left=1
          goto 100
        endif

      enddo

c     Exit point for subroutine
 100  continue

      return

      end

c     -------------------------------------------------------------------
c     Iterative Euler time step (ISENTROPIC)
c     -------------------------------------------------------------------

      subroutine euler_isen(x1,y1,p1,left,x0,y0,p0,theta,
     >                reltpos0,reltpos1,
     >                deltat,numit,jump,mdv,wfactor,fbflag,
     >                spt0,spt1,p3d0,p3d1,uut0,uut1,vvt0,vvt1,
     >                sth0,sth1,tht0,tht1,
     >                xmin,ymin,dx,dy,per,hem,nx,ny,nz)

      implicit none

c     Declaration of subroutine parameters
      integer      nx,ny,nz
      real         x1,y1,p1
      integer      left
      real         x0,y0,p0
      real         reltpos0,reltpos1
      real         deltat
      integer      numit
      integer      jump
      real         wfactor
      integer      fbflag
      real         spt0(nx*ny)   ,spt1(nx*ny)
      real         sth0(nx*ny)   ,sth1(nx*ny)
      real         uut0(nx*ny*nz),uut1(nx*ny*nz)
      real         vvt0(nx*ny*nz),vvt1(nx*ny*nz)
      real         p3d0(nx*ny*nz),p3d1(nx*ny*nz)
      real         tht0(nx*ny*nz),tht1(nx*ny*nz)
      real         xmin,ymin,dx,dy
      real         per
      integer      hem
      real         mdv
      real         theta

c     Numerical and physical constants
      real         deltay
      parameter    (deltay=1.112E5)  ! Distance in m between 2 lat circles
      real         pi
      parameter    (pi=3.1415927)    ! Pi

c     Auxiliary variables
      real         xmax,ymax
      real         xind,yind,pind
      real         u0,v0,w0,u1,v1,w1,u,v,w,sp
      integer      icount
      character    ch

c     Externals
      real         int_index4
      external     int_index4

c     Reset the flag for domain-leaving
      left=0

c     Set the esat-north bounray of the domain
      xmax = xmin+real(nx-1)*dx
      ymax = ymin+real(ny-1)*dy

C     Interpolate wind fields to starting position (x0,y0,p0)
      call get_index4 (xind,yind,pind,x0,y0,p0,reltpos0,
     >                 p3d0,p3d1,spt0,spt1,3,
     >                 nx,ny,nz,xmin,ymin,dx,dy,mdv)
      u0 = int_index4(uut0,uut1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)
      v0 = int_index4(vvt0,vvt1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)

C     For first iteration take ending position equal to starting position
      x1=x0
      y1=y0
      p1=p0

C     Iterative calculation of new position
      do icount=1,numit

C        Calculate new winds for advection
         call get_index4 (xind,yind,pind,x1,y1,p1,reltpos1,
     >                    p3d0,p3d1,spt0,spt1,3,
     >                    nx,ny,nz,xmin,ymin,dx,dy,mdv)
         u1 = int_index4(uut0,uut1,nx,ny,nz,xind,yind,pind,reltpos1,mdv)
         v1 = int_index4(vvt0,vvt1,nx,ny,nz,xind,yind,pind,reltpos1,mdv)

c        Get the new velocity in between
         u=(u0+u1)/2.
         v=(v0+v1)/2.

C        Calculate new positions
         x1 = x0 + fbflag*u*deltat/(deltay*cos(y0*pi/180.))
         y1 = y0 + fbflag*v*deltat/deltay

c        Get the pressure on the isentropic surface at the new position
         call get_index4 (xind,yind,pind,x1,y1,theta,reltpos1,
     >                    tht0,tht1,sth0,sth1,1,
     >                    nx,ny,nz,xmin,ymin,dx,dy,mdv)
         p1 = int_index4(p3d0,p3d1,nx,ny,nz,xind,yind,pind,reltpos1,mdv)

c       Handle pole problems (crossing and near pole trajectory)
        if ((hem.eq.1).and.(y1.gt.90.)) then
          y1=180.-y1
          x1=x1+per/2.
        endif
        if ((hem.eq.1).and.(y1.lt.-90.)) then
          y1=-180.-y1
          x1=x1+per/2.
        endif
        if (y1.gt.89.99) then
           y1=89.99
        endif

c       Handle crossings of the dateline
        if ((hem.eq.1).and.(x1.gt.xmin+per-dx)) then
           x1=xmin+amod(x1-xmin,per)
        endif
        if ((hem.eq.1).and.(x1.lt.xmin)) then
           x1=xmin+per+amod(x1-xmin,per)
        endif

C       Interpolate surface pressure to actual position
        call get_index4 (xind,yind,pind,x1,y1,1050.,reltpos1,
     >                   p3d0,p3d1,spt0,spt1,3,
     >                   nx,ny,nz,xmin,ymin,dx,dy,mdv)
        sp = int_index4 (spt0,spt1,nx,ny,1,xind,yind,1.,reltpos1,mdv)

c       Handle trajectories which cross the lower boundary (jump flag)
        if ((jump.eq.1).and.(p1.gt.sp)) p1=sp-10.

C       Check if trajectory leaves data domain
        if ( ( (hem.eq.0).and.(x1.lt.xmin)    ).or.
     >       ( (hem.eq.0).and.(x1.gt.xmax-dx) ).or.
     >         (y1.lt.ymin).or.(y1.gt.ymax).or.(p1.gt.sp) )
     >  then
          left=1
          goto 100
        endif

      enddo

c     Exit point for subroutine
 100  continue

      return

      end

c     -------------------------------------------------------------------
c     Iterative Euler time step (MODEL-LEVEL, 2D)
c     -------------------------------------------------------------------

      subroutine euler_2d(x1,y1,p1,left,x0,y0,p0,zindex,
     >                reltpos0,reltpos1,
     >                deltat,numit,jump,mdv,wfactor,fbflag,
     >                spt0,spt1,p3d0,p3d1,uut0,uut1,vvt0,vvt1,
     >                xmin,ymin,dx,dy,per,hem,nx,ny,nz)

      implicit none

c     Declaration of subroutine parameters
      integer      nx,ny,nz
      real         x1,y1,p1
      integer      left
      real         x0,y0,p0
      real         reltpos0,reltpos1
      real         deltat
      integer      numit
      integer      jump
      real         wfactor
      integer      fbflag
      real         spt0(nx*ny)   ,spt1(nx*ny)
      real         uut0(nx*ny*nz),uut1(nx*ny*nz)
      real         vvt0(nx*ny*nz),vvt1(nx*ny*nz)
      real         p3d0(nx*ny*nz),p3d1(nx*ny*nz)
      real         xmin,ymin,dx,dy
      real         per
      integer      hem
      real         mdv
      real         zindex

c     Numerical and physical constants
      real         deltay
      parameter    (deltay=1.112E5)  ! Distance in m between 2 lat circles
      real         pi
      parameter    (pi=3.1415927)    ! Pi
      real         eps
      parameter    (eps=0.001)

c     Auxiliary variables
      real         xmax,ymax
      real         xind,yind,pind
      real         u0,v0,w0,u1,v1,w1,u,v,w,sp
      integer      icount
      character    ch

c     Externals
      real         int_index4
      external     int_index4

c     Reset the flag for domain-leaving
      left=0

c     Set the esat-north bounray of the domain
      xmax = xmin+real(nx-1)*dx
      ymax = ymin+real(ny-1)*dy

C     Interpolate wind fields to starting position (x0,y0,p0)
      call get_index4 (xind,yind,pind,x0,y0,p0,reltpos0,
     >                 p3d0,p3d1,spt0,spt1,3,
     >                 nx,ny,nz,xmin,ymin,dx,dy,mdv)
      u0 = int_index4(uut0,uut1,nx,ny,nz,xind,yind,zindex,reltpos0,mdv)
      v0 = int_index4(vvt0,vvt1,nx,ny,nz,xind,yind,zindex,reltpos0,mdv)

C     For first iteration take ending position equal to starting position
      x1=x0
      y1=y0
      p1=p0

C     Iterative calculation of new position
      do icount=1,numit

C        Calculate new winds for advection
         call get_index4 (xind,yind,pind,x1,y1,p1,reltpos1,
     >                    p3d0,p3d1,spt0,spt1,3,
     >                    nx,ny,nz,xmin,ymin,dx,dy,mdv)
         u1=int_index4(uut0,uut1,nx,ny,nz,xind,yind,zindex,reltpos1,mdv)
         v1=int_index4(vvt0,vvt1,nx,ny,nz,xind,yind,zindex,reltpos1,mdv)
         if ( abs(u1-mdv).lt.eps ) then
             left = 1
             goto 100
         endif
         if ( abs(v1-mdv).lt.eps ) then
             left = 1
             goto 100
         endif

c        Get the new velocity in between
         u=(u0+u1)/2.
         v=(v0+v1)/2.

C        Calculate new positions
         x1 = x0 + fbflag*u*deltat/(deltay*cos(y0*pi/180.))
         y1 = y0 + fbflag*v*deltat/deltay

c        Get the pressure on the model surface at the new position
         xind = (x1 - xmin ) / dx + 1.
         yind = (y1 - ymin ) / dy + 1.
         p1 =
     >     int_index4(p3d0,p3d1,nx,ny,nz,xind,yind,zindex,reltpos1,mdv)
         if ( abs(p1-mdv).lt.eps ) then
             left = 1
             goto 100
         endif

c       Handle pole problems (crossing and near pole trajectory)
        if ((hem.eq.1).and.(y1.gt.90.)) then
          y1=180.-y1
          x1=x1+per/2.
        endif
        if ((hem.eq.1).and.(y1.lt.-90.)) then
          y1=-180.-y1
          x1=x1+per/2.
        endif
        if (y1.gt.89.99) then
           y1=89.99
        endif

c       Handle crossings of the dateline
        if ((hem.eq.1).and.(x1.gt.xmin+per-dx)) then
           x1=xmin+amod(x1-xmin,per)
        endif
        if ((hem.eq.1).and.(x1.lt.xmin)) then
           x1=xmin+per+amod(x1-xmin,per)
        endif

C       Interpolate surface pressure to actual position
        call get_index4 (xind,yind,pind,x1,y1,1050.,reltpos1,
     >                   p3d0,p3d1,spt0,spt1,3,
     >                   nx,ny,nz,xmin,ymin,dx,dy,mdv)
        sp = int_index4 (spt0,spt1,nx,ny,1,xind,yind,1.,reltpos1,mdv)

c       Handle trajectories which cross the lower boundary (jump flag)
        if ((jump.eq.1).and.(p1.gt.sp)) p1=sp-10.

C       Check if trajectory leaves data domain
        if ( ( (hem.eq.0).and.(x1.lt.xmin)    ).or.
     >       ( (hem.eq.0).and.(x1.gt.xmax-dx) ).or.
     >         (y1.lt.ymin).or.(y1.gt.ymax).or.(p1.gt.sp) )
     >  then
          left=1
          goto 100
        endif

      enddo

c     Exit point for subroutine
 100  continue

      return

      end

c     -------------------------------------------------------------------
c     Runge-Kutta (4th order) time-step
c     -------------------------------------------------------------------

      subroutine runge(x1,y1,p1,left,x0,y0,p0,reltpos0,reltpos1,
     >                 deltat,numit,jump,mdv,wfactor,fbflag,
     >		       spt0,spt1,p3d0,p3d1,uut0,uut1,vvt0,vvt1,wwt0,wwt1,
     >                 xmin,ymin,dx,dy,per,hem,nx,ny,nz)

      implicit none

c     Declaration of subroutine parameters
      integer      nx,ny,nz
      real         x1,y1,p1
      integer      left
      real	   x0,y0,p0
      real         reltpos0,reltpos1
      real   	   deltat
      integer      numit
      integer      jump
      real         wfactor
      integer      fbflag
      real     	   spt0(nx*ny)   ,spt1(nx*ny)
      real         uut0(nx*ny*nz),uut1(nx*ny*nz)
      real 	   vvt0(nx*ny*nz),vvt1(nx*ny*nz)
      real         wwt0(nx*ny*nz),wwt1(nx*ny*nz)
      real         p3d0(nx*ny*nz),p3d1(nx*ny*nz)
      real         xmin,ymin,dx,dy
      real         per
      integer      hem
      real         mdv

c     Numerical and physical constants
      real         deltay
      parameter    (deltay=1.112E5)  ! Distance in m between 2 lat circles
      real         pi                       
      parameter    (pi=3.1415927)    ! Pi

c     Auxiliary variables
      real         xmax,ymax
      real	   xind,yind,pind
      real	   u0,v0,w0,u1,v1,w1,u,v,w,sp
      integer	   icount,n
      real         xs,ys,ps,xk(4),yk(4),pk(4)
      real         reltpos

c     Externals    
      real         int_index4
      external     int_index4

c     Reset the flag for domain-leaving
      left=0

c     Set the esat-north bounray of the domain
      xmax = xmin+real(nx-1)*dx
      ymax = ymin+real(ny-1)*dy

c     Apply the Runge Kutta scheme
      do n=1,4
 
c       Get intermediate position and relative time
        if (n.eq.1) then
          xs=0.
          ys=0.
          ps=0.
          reltpos=reltpos0
        else if (n.eq.4) then
          xs=xk(3)
          ys=yk(3)
          ps=pk(3)
          reltpos=reltpos1
        else
          xs=xk(n-1)/2.
          ys=yk(n-1)/2.
          ps=pk(n-1)/2.
          reltpos=(reltpos0+reltpos1)/2.
        endif
        
C       Calculate new winds for advection
        call get_index4 (xind,yind,pind,x0+xs,y0+ys,p0+ps,reltpos,
     >                   p3d0,p3d1,spt0,spt1,3,
     >                   nx,ny,nz,xmin,ymin,dx,dy,mdv)
        u = int_index4 (uut0,uut1,nx,ny,nz,xind,yind,pind,reltpos,mdv)
        v = int_index4 (vvt0,vvt1,nx,ny,nz,xind,yind,pind,reltpos,mdv)
        w = int_index4 (wwt0,wwt1,nx,ny,nz,xind,yind,pind,reltpos,mdv)
         
c       Force the near-surface wind to zero
        if (pind.lt.1.) w1=w1*pind
 
c       Update position and keep them
        xk(n)=fbflag*u*deltat/(deltay*cos(y0*pi/180.))
        yk(n)=fbflag*v*deltat/deltay
        pk(n)=fbflag*w*deltat*wfactor/100.

      enddo
 
C     Calculate new positions
      x1=x0+(1./6.)*(xk(1)+2.*xk(2)+2.*xk(3)+xk(4))
      y1=y0+(1./6.)*(yk(1)+2.*yk(2)+2.*yk(3)+yk(4))
      p1=p0+(1./6.)*(pk(1)+2.*pk(2)+2.*pk(3)+pk(4))

c     Handle pole problems (crossing and near pole trajectory)
      if ((hem.eq.1).and.(y1.gt.90.)) then
         y1=180.-y1
         x1=x1+per/2.
      endif
      if ((hem.eq.1).and.(y1.lt.-90.)) then
         y1=-180.-y1
         x1=x1+per/2.
      endif
      if (y1.gt.89.99) then
         y1=89.99
      endif
      
c     Handle crossings of the dateline
      if ((hem.eq.1).and.(x1.gt.xmin+per-dx)) then
         x1=xmin+amod(x1-xmin,per)
      endif
      if ((hem.eq.1).and.(x1.lt.xmin)) then
         x1=xmin+per+amod(x1-xmin,per)
      endif
      
C     Interpolate surface pressure to actual position
      call get_index4 (xind,yind,pind,x1,y1,1050.,reltpos1,
     >                 p3d0,p3d1,spt0,spt1,3,
     >                 nx,ny,nz,xmin,ymin,dx,dy,mdv)
      sp = int_index4 (spt0,spt1,nx,ny,1,xind,yind,1,reltpos,mdv)

c     Handle trajectories which cross the lower boundary (jump flag)
      if ((jump.eq.1).and.(p1.gt.sp)) p1=sp-10.
      
C     Check if trajectory leaves data domain
      if ( ( (hem.eq.0).and.(x1.lt.xmin)    ).or.
     >     ( (hem.eq.0).and.(x1.gt.xmax-dx) ).or.
     >       (y1.lt.ymin).or.(y1.gt.ymax).or.(p1.gt.sp) )
     >then
         left=1
         goto 100
      endif
      
c     Exit point fdor subroutine
 100  continue

      return
      end
