PROGRAM trace

  ! ********************************************************************
  ! *                                                                  *
  ! * Trace fields along trajectories                                  *
  ! *                                                                  *
  ! * April 1993: First version (Heini Wernli)                         *
  ! * 2008-2009 : Major upgrades (Michael Sprenger)                    *
  ! * Mar 2012: Clustering option (Bojan Skerlak)                      *  
  ! * Nov 2012: Circle options (")                                     *
  ! * Jul 2013: user-defined PV,TH @ clustering mode (")               *
  ! *                                                                  *
  ! ********************************************************************

  implicit none

  ! --------------------------------------------------------------------
  ! Declaration of parameters
  ! --------------------------------------------------------------------

  ! Maximum number of levels for input files
  integer :: nlevmax
  parameter (nlevmax=100)

  ! Maximum number of input files (dates, length of trajectories)
  integer :: ndatmax
  parameter (ndatmax=500)

  ! Numerical epsilon (for float comparison)
  real :: eps
  parameter (eps=0.001)

  ! Conversion factors
  real :: pi180                                   ! deg -> rad
  parameter (pi180=3.14159/180.)
  real :: deg2km                                  ! deg -> km (at equator)
  parameter (deg2km=111.)
  real :: pir
  parameter (pir=255032235.95489)                 ! 2*Pi*R^2 Bojan

  ! Prefix for primary and secondary fields
  character :: charp
  character :: chars
  parameter (charp='P')
  parameter (chars='S')

  ! --------------------------------------------------------------------
  ! Declaration of variables
  ! --------------------------------------------------------------------

  ! Input and output format for trajectories (see iotra.f)
  integer :: inpmode
  integer :: outmode

  ! Input parameters
  character(len=80) :: inpfile         ! Input trajectory file
  character(len=80) :: outfile         ! Output trajectory file
  integer :: ntra                      ! Number of trajectories
  integer :: ncol                      ! Number of columns (including time, lon, lat, p)
  integer :: ntim                      ! Number of times per trajectory
  integer :: ntrace0                   ! Number of trace variables
  character(len=80) :: tvar0(200)      ! Tracing variable (with mode specification)
  character(len=80) :: tvar(200)       ! Tracing variable name (only the variable)
  character(len=80) :: tfil(200)       ! Filename prefix
  real :: fac(200)                     ! Scaling factor
  real :: shift_val(200)               ! Shift in space and time relative to trajectory position
  character(len=80) :: shift_dir(200)  ! Direction of shift
  integer :: compfl(200)               ! Computation flag (1=compute)
  integer :: numdat                    ! Number of input files
  character(len=11) :: dat(ndatmax)    ! Dates of input files
  real :: timeinc                      ! Time increment between input files
  real :: tst                          ! Time shift of start relative to first data file
  real :: ten                          ! Time shift of end relatiev to first data file
  character(len=20) :: startdate       ! First time/date on trajectory
  character(len=20) :: enddate         ! Last time/date on trajectory
  integer :: ntrace1                   ! Count trace and additional variables
  character(len=80) :: timecheck       ! Either 'yes' or 'no'
  character(len=80) :: intmode         ! Interpolation mode ('normal', 'nearest') Bojan ('clustering','circle_avg','circle_max','circle_min')

  ! Trajectories
  real,allocatable, dimension (:,:,:) :: trainp          ! Input trajectories (ntra,ntim,ncol)
  real,allocatable, dimension (:,:,:) :: traint          ! Internal trajectories (ntra,ntim,ncol+ntrace1)
  real,allocatable, dimension (:,:,:) :: traout          ! Output trajectories (ntra,ntim,ncol+ntrace0)
  integer :: reftime(6)                                  ! Reference date
  character(len=80) :: varsinp(200)                      ! Field names for input trajectory
  character(len=80) :: varsint(200)                      ! Field names for internal trajectory
  character(len=80) :: varsout(200)                      ! Field names for output trajectory
  integer :: fid,fod                                     ! File identifier for inp and out trajectories
  real :: x0,y0,p0                                       ! Position of air parcel (physical space)
  real :: reltpos0                                       ! Relative time of air parcel
  real :: xind,yind,pind                                 ! Position of air parcel (grid space)
  integer :: fbflag                                      ! Flag for forward (1) or backward (-1) trajectories
  integer :: fok(200)                                    ! Flag whether field is ready

  ! Meteorological fields
  real,allocatable, dimension (:)     :: spt0,spt1       ! Surface pressure
  real,allocatable, dimension (:)     :: p3t0,p3t1       ! 3d-pressure
  real,allocatable, dimension (:)     :: f3t0,f3t1       ! 3d field for tracing
  character(len=80) :: svars(100)                        ! List of variables on S file
  character(len=80) :: pvars(100)                        ! List of variables on P file
  integer :: n_svars                                     ! Number of variables on S file
  integer :: n_pvars                                     ! Number of variables on P file

  ! Grid description
  real :: pollon,pollat   ! Longitude/latitude of pole
  real :: ak(100)         ! Vertical layers and levels
  real :: bk(100)
  real :: xmin,xmax       ! Zonal grid extension
  real :: ymin,ymax       ! Meridional grid extension
  integer :: nx,ny,nz     ! Grid dimensions
  real :: dx,dy           ! Horizontal grid resolution
  integer :: hem          ! Flag for hemispheric domain
  integer :: per          ! Flag for periodic domain
  real :: stagz           ! Vertical staggering
  real :: mdv             ! Missing data value

  ! Auxiliary variables
  integer :: i,j,k,l,m,n
  real :: rd
  character(len=80) :: filename,varname
  real :: time0,time1,reltpos
  integer :: itime0,itime1
  integer :: stat
  real :: tstart
  integer :: iloaded0,iloaded1
  real :: f0
  real :: frac
  real :: tload,tfrac
  integer :: isok
  character :: ch
  integer :: ind
  integer :: ind1,ind2,ind3,ind4,ind5
  integer :: ind6,ind7,ind8,ind9,ind0
  integer :: noutside
  real :: delta
  integer :: itrace0
  character(len=80) :: string
  integer err_c1,err_c2,err_c3

  ! Bojan
  real    :: dist,circlesum,circlemax,circlemin,circleavg,radius       ! distance (great circle), sum/max/min/avg in circle, radius of circle
  integer :: ist,jst,kst,sp,ml,mr,nd,nu                                ! ijk in stack, sp=stack counter, ml (left), mr (right), nd (down), nu (up)
  integer :: lci,lcj,xindb,xindf                                       ! label count i and j, xind back, xind forward
  integer :: yindb,yindf,pindb,pindf                                   ! yind back, yind forward, pind back, pind forward
  integer :: pvpos,thpos                                               ! position of variables PV and TH in trajectory
  real    :: tropo_pv,tropo_th                                         ! values of PV and TH at dynamical tropopause
  integer, allocatable, dimension (:)   :: stackx,stacky               ! lon/lat of stack
  integer, allocatable, dimension (:)   :: lblcount                    ! counter for label
  integer, allocatable, dimension (:,:) :: connect                     ! array that keeps track of the visited grid points
  real, allocatable, dimension (:)      :: lbl                         ! label
  real, allocatable, dimension (:)      :: circlelon,circlelat         ! value of f, lon and lat in circle
  real, allocatable, dimension (:)      :: circlef,circlearea          ! value of f, lon and lat in circle
  real, allocatable, dimension (:)      :: longrid,latgrid             ! arrays of longitude and latitude of grid points
  ! Bojan

  ! Externals
  real :: int_index4
  external                               int_index4
  real :: sdis ! Bojan: need function sdis (calculates great circle distance)
  external                               sdis ! Bojan: need function sdis

  ! --------------------------------------------------------------------
  ! Start of program, Read parameters, get grid parameters
  ! --------------------------------------------------------------------

  ! Write start message
  print*,'========================================================='
  print*,'              *** START OF PROGRAM TRACE ***'
  print*

  ! Read parameters
  open(10,file='trace.param')
   read(10,*) inpfile
   read(10,*) outfile
   read(10,*) startdate
   read(10,*) enddate
   read(10,*) fbflag
   read(10,*) numdat
   if ( fbflag.eq.1) then
      do i=1,numdat
         read(10,'(a11)') dat(i)
      enddo
   else
      do i=numdat,1,-1
         read(10,'(a11)') dat(i)
      enddo
   endif
   read(10,*) timeinc
   read(10,*) tst
   read(10,*) ten
   read(10,*) ntra
   read(10,*) ntim
   read(10,*) ncol
   read(10,*) ntrace0
   do i=1,ntrace0
      read(10,*) tvar(i), fac(i), compfl(i), tfil(i)
   enddo
   read(10,*) n_pvars
   do i=1,n_pvars
      read(10,*) pvars(i)
   enddo
   read(10,*) n_svars
   do i=1,n_svars
      read(10,*) svars(i)
   enddo
   read(10,*) timecheck
   read(10,*) intmode
   read(10,*) radius
   read(10,*) tropo_pv
   read(10,*) tropo_th
  close(10)

  ! Bojan: error if radius < 0
  if (((intmode .eq. "circle_avg") .or. (intmode .eq. "circle_min") .or. (intmode .eq. "circle_max")) .and. (radius .lt. 0)) then
     print*,'ERROR (circle): radius < 0!'
     stop
  endif

  ! Remove commented tracing fields
  itrace0 = 1
  do while ( itrace0.le.ntrace0)
     string = tvar(itrace0)
     if ( string(1:1).eq.'#' ) then
        do i=itrace0,ntrace0-1
           tvar(i)   = tvar(i+1)
           fac(i)    = fac(i+1)
           compfl(i) = compfl(i+1)
           tfil(i)   = tfil(i+1)
        enddo
        ntrace0 = ntrace0 - 1
     else
        itrace0 = itrace0 + 1
     endif
  enddo

  ! Save the tracing variable (including all mode specifications)
  do i=1,ntrace0
     tvar0(i) = tvar(i)
  enddo

  ! Set the formats of the input and output files
  call mode_tra(inpmode,inpfile)
  if (inpmode.eq.-1) inpmode=1
  call mode_tra(outmode,outfile)
  if (outmode.eq.-1) outmode=1

  ! Convert time shifts <tst,ten> from <hh.mm> into fractional time
  call hhmm2frac(tst,frac)
  tst = frac
  call hhmm2frac(ten,frac)
  ten = frac

  ! Set the time for the first data file (depending on forward/backward mode)
  if (fbflag.eq.1) then
    tstart = -tst
  else
    tstart = tst
  endif

  ! Read the constant grid parameters (nx,ny,nz,xmin,xmax,ymin,ymax,pollon,pollat)
  ! The negative <-fid> of the file identifier is used as a flag for parameter retrieval
  filename = charp//dat(1)
  varname  = 'U'
  call input_open (fid,filename)
  call input_grid (-fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,tstart,pollon,pollat,rd,rd,nz,rd,rd,rd,timecheck)
  call input_close(fid)

  ! Allocate memory for some meteorological arrays
  allocate(spt0(nx*ny),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array spt0 ***'   ! Surface pressure
  allocate(spt1(nx*ny),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array spt1 ***'
  allocate(p3t0(nx*ny*nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array p3t0 ***'   ! Pressure
  allocate(p3t1(nx*ny*nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array p3t1 ***'
  allocate(f3t0(nx*ny*nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array p3t0 ***'   ! Tracing field
  allocate(f3t1(nx*ny*nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array p3t1 ***'

  ! Get memory for trajectory arrays
  allocate(trainp(ntra,ntim,ncol),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array tra      ***'

  ! Bojan
  ! allocate memory for clustering mode 
  if (intmode .eq. 'clustering') then
     allocate(lbl(8),stat=stat)
     if (stat.ne.0) print*,'*** error allocating array lbl      ***'
     allocate(lblcount(5),stat=stat)
     if (stat.ne.0) print*,'*** error allocating array lblcount ***'
  endif
  ! allocate memory for circle mode
  if ( (intmode.eq.'circle_avg') .or. (intmode.eq.'circle_min') .or. (intmode.eq.'circle_max') ) then
     allocate(connect(nx,ny),stat=stat)
     if (stat.ne.0) print*,'*** error allocating connect ***'
     allocate(stackx(nx*ny),stat=stat)
     if (stat.ne.0) print*,'*** error allocating stackx ***'
     allocate(stacky(nx*ny),stat=stat)
     if (stat.ne.0) print*,'*** error allocating stacky ***'
     allocate(circlelon(nx*ny),stat=stat)
     if (stat.ne.0) print*,'*** error allocating circlelon ***'
     allocate(circlelat(nx*ny),stat=stat)
     if (stat.ne.0) print*,'*** error allocating circlelat ***'
     allocate(circlef(nx*ny),stat=stat)
     if (stat.ne.0) print*,'*** error allocating circlef ***'
     allocate(circlearea(nx*ny),stat=stat)
     if (stat.ne.0) print*,'*** error allocating circlearea ***'
     allocate(longrid(nx),stat=stat)
     if (stat.ne.0) print*,'*** error allocating longrid ***'
     allocate(latgrid(ny),stat=stat)
     if (stat.ne.0) print*,'*** error allocating latgrid ***'
     do m=1,nx
        longrid(m)=xmin+dx*(m-1)
     enddo
     do n=1,ny
        latgrid(n)=ymin+dy*(n-1)
     enddo
  endif
  ! Bojan

  ! Set the flags for periodic domains
  if ( abs(xmax-xmin-360.).lt.eps ) then
     per = 1
  elseif ( abs(xmax-xmin-360.+dx).lt.eps ) then
     per = 2
  else
     per = 0
  endif

  ! Set logical flag for periodic data set (hemispheric or not)
  hem = 0
  if (per.eq.0.) then
     delta=xmax-xmin-360.
     if (abs(delta+dx).lt.eps) then               ! Program aborts: arrays must be closed
        print*,' ERROR: arrays must be closed... Stop'
     else if (abs(delta).lt.eps) then ! Periodic and hemispheric
       hem=1
       per=360.
    endif
  else                                            ! Periodic and hemispheric
     hem=1
  endif

  ! Write some status information
  print*,'---- INPUT PARAMETERS -----------------------------------'
  print*
  print*,'  Input trajectory file  : ',trim(inpfile)
  print*,'  Output trajectory file : ',trim(outfile)
  print*,'  Format of input file   : ',inpmode
  print*,'  Format of output file  : ',outmode
  print*,'  Forward/backward       : ',fbflag
  print*,'  #tra                   : ',ntra
  print*,'  #col                   : ',ncol
  print*,'  #tim                   : ',ntim
  print*,'  No time check          : ',trim(timecheck)
  print*,'  Interpolation mode     : ',trim(intmode)
  ! Bojan
  if (trim(intmode) .eq. "clustering") then
     print*,'  Tropopause PV [pvu]    : ',tropo_pv
     print*,'  Tropopause TH [K]      : ',tropo_th
  endif
  do i=1,ntrace0
     if (compfl(i).eq.0) then
        print*,'  Tracing field          : ', trim(tvar(i)), fac(i), ' 0 ', tfil(i)
     else
        print*,'  Tracing field          : ', trim(tvar(i)), fac(i), '  1 ', tfil(i)
     endif
  enddo
  print*
  print*,'---- INPUT DATA FILES -----------------------------------'
  print*
  call frac2hhmm(tstart,tload)
  print*,'  Time of 1st data file  : ',tload
  print*,'  #input files           : ',numdat
  print*,'  time increment         : ',timeinc
  call frac2hhmm(tst,tload)
  print*,'  Shift of start         : ',tload
  call frac2hhmm(ten,tload)
  print*,'  Shift of end           : ',tload
  print*,'  First/last input file  : ',trim(dat(1)), ' ... ', trim(dat(numdat))
  print*,'  Primary variables      : ',trim(pvars(1))
  do i=2,n_pvars
     print*,'                         : ',trim(pvars(i))
  enddo
  if ( n_svars.ge.1 ) then
     print*,'  Secondary variables    : ',trim(svars(1))
     do i=2,n_svars
        print*,'                         : ',trim(svars(i))
     enddo
  endif
  print*
  print*,'---- CONSTANT GRID PARAMETERS ---------------------------'
  print*
  print*,'  xmin,xmax     : ',xmin,xmax
  print*,'  ymin,ymax     : ',ymin,ymax
  print*,'  dx,dy         : ',dx,dy
  print*,'  pollon,pollat : ',pollon,pollat
  print*,'  nx,ny,nz      : ',nx,ny,nz
  print*,'  per, hem      : ',per,hem
  print*

  ! --------------------------------------------------------------------
  ! Load the input trajectories
  ! --------------------------------------------------------------------

  ! Read the input trajectory file
  call ropen_tra(fid,inpfile,ntra,ntim,ncol,reftime,varsinp,inpmode)
  call read_tra (fid,trainp,ntra,ntim,ncol,inpmode)
  call close_tra(fid,inpmode)

  ! Check that first four columns correspond to time,lon,lat,p
  if ( (varsinp(1).ne.'time' ).or. &
       (varsinp(2).ne.'xpos' ).and.(varsinp(2).ne.'lon' ).or. &
       (varsinp(3).ne.'ypos' ).and.(varsinp(3).ne.'lat' ).or. &
       (varsinp(4).ne.'ppos' ).and.(varsinp(4).ne.'p'   ) ) then
     print*,' ERROR: problem with input trajectories ...'
     stop
  endif
  varsinp(1) = 'time'
  varsinp(2) = 'lon'
  varsinp(3) = 'lat'
  varsinp(4) = 'p'

  ! Write some status information of the input trajectories
  print*,'---- INPUT TRAJECTORIES ---------------------------------'
  print*
  print*,' Start date             : ',trim(startdate)
  print*,' End date               : ',trim(enddate)
  print*,' Reference time (year)  : ',reftime(1)
  print*,'                (month) : ',reftime(2)
  print*,'                (day)   : ',reftime(3)
  print*,'                (hour)  : ',reftime(4)
  print*,'                (min)   : ',reftime(5)
  print*,' Time range (min)       : ',reftime(6)
  do i=1,ncol
     print*,' Var                    :',i,trim(varsinp(i))
  enddo
  print*

  ! Check that first time is 0 - otherwise the tracing will produce
  ! wrong results because later in the code only absolute times are
  ! considered: <itime0   = int(abs(tfrac-tstart)/timeinc) + 1>. This
  ! will be changed in a future version.
  if ( abs( trainp(1,1,1) ).gt.eps ) then
     print*,' ERROR: First time of trajectory must be 0, i.e. '
     print*,'     correspond to the reference date. Otherwise'
     print*,'     the tracing will give wrong results... STOP'
     stop
  endif

  ! --------------------------------------------------------------------
  ! Check dependencies for trace fields which must be calculated
  ! --------------------------------------------------------------------

  ! Set the counter for extra fields
  ntrace1 = ntrace0

  ! Loop over all tracing variables
  i = 1
  do while (i.le.ntrace1)

     ! Skip fields which must be available on the input files
     if (i.le.ntrace0) then
        if (compfl(i).eq.0) goto 100
     endif

     ! Get the dependencies for potential temperature (TH)
     if ( tvar(i).eq.'TH' ) then
        varname='P'                                   ! P
        call add2list(varname,tvar,ntrace1)
        varname='T'                                   ! T
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for potential temperature (TH)
     elseif ( tvar(i).eq.'RHO' ) then
        varname='P'                                   ! P
        call add2list(varname,tvar,ntrace1)
        varname='T'                                   ! T
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for relative humidity (RH)
     elseif ( tvar(i).eq.'RH' ) then
        varname='P'                                   ! P
        call add2list(varname,tvar,ntrace1)
        varname='T'                                   ! T
        call add2list(varname,tvar,ntrace1)
        varname='Q'                                   ! Q
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for equivalent potential temperature (THE)
     elseif ( tvar(i).eq.'THE' ) then
        varname='P'                                   ! P
        call add2list(varname,tvar,ntrace1)
        varname='T'                                   ! T
        call add2list(varname,tvar,ntrace1)
        varname='Q'                                   ! Q
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for latent heating rate (LHR)
     elseif ( tvar(i).eq.'LHR' ) then
        varname='P'                                   ! P
        call add2list(varname,tvar,ntrace1)
        varname='T'                                   ! T
        call add2list(varname,tvar,ntrace1)
        varname='Q'                                   ! Q
        call add2list(varname,tvar,ntrace1)
        varname='OMEGA'                               ! OMEGA
        call add2list(varname,tvar,ntrace1)
        varname='RH'                                  ! RH
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for wind speed (VEL)
     elseif ( tvar(i).eq.'VEL' ) then
        varname='U'                                   ! U
        call add2list(varname,tvar,ntrace1)
        varname='V'                                   ! V
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for wind direction (DIR)
     elseif ( tvar(i).eq.'DIR' ) then
        varname='U'                                   ! U
        call add2list(varname,tvar,ntrace1)
        varname='V'                                   ! V
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for du/dx (DUDX)
     elseif ( tvar(i).eq.'DUDX' ) then
        varname='U:+1DLON'                            ! U:+1DLON
        call add2list(varname,tvar,ntrace1)
        varname='U:-1DLON'                            ! U:-1DLON
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for dv(dx (DVDX)
     elseif ( tvar(i).eq.'DVDX' ) then
        varname='V:+1DLON'                            ! V:+1DLON
        call add2list(varname,tvar,ntrace1)
        varname='V:-1DLON'                            ! V:-1DLON
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for du/dy (DUDY)
     elseif ( tvar(i).eq.'DUDY' ) then
        varname='U:+1DLAT'                            ! U:+1DLAT
        call add2list(varname,tvar,ntrace1)
        varname='U:-1DLAT'                            ! U:-1DLAT
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for dv/dy (DVDY)
     elseif ( tvar(i).eq.'DVDY' ) then
        varname='V:+1DLAT'                            ! V:+1DLAT
        call add2list(varname,tvar,ntrace1)
        varname='V:-1DLAT'                            ! V:-1DLAT
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for du/dp (DUDP)
     elseif ( tvar(i).eq.'DUDP' ) then
        varname='U:+1DP'                              ! U:+1DP
        call add2list(varname,tvar,ntrace1)
        varname='U:-1DP'                              ! U:-1DP
        call add2list(varname,tvar,ntrace1)
        varname='P:+1DP'                              ! P:+1DP
        call add2list(varname,tvar,ntrace1)
        varname='P:-1DP'                              ! P:-1DP
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for dv/dp (DVDP)
     elseif ( tvar(i).eq.'DVDP' ) then
        varname='V:+1DP'                              ! V:+1DP
        call add2list(varname,tvar,ntrace1)
        varname='V:-1DP'                              ! V:-1DP
        call add2list(varname,tvar,ntrace1)
        varname='P:+1DP'                              ! P:+1DP
        call add2list(varname,tvar,ntrace1)
        varname='P:-1DP'                              ! P:-1DP
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for dt/dx (DTDX)
     elseif ( tvar(i).eq.'DTDX' ) then
        varname='T:+1DLON'                            ! T:+1DLON
        call add2list(varname,tvar,ntrace1)
        varname='T:-1DLON'                            ! T:-1DLON
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for dth/dy (DTHDY)
     elseif ( tvar(i).eq.'DTHDY' ) then
        varname='T:+1DLAT'                            ! T:+1DLON
        call add2list(varname,tvar,ntrace1)
        varname='T:-1DLAT'                            ! T:-1DLON
        call add2list(varname,tvar,ntrace1)
        varname='P'                                   ! P
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for dth/dx (DTHDX)
     elseif ( tvar(i).eq.'DTHDX' ) then
        varname='T:+1DLON'                            ! T:+1DLON
        call add2list(varname,tvar,ntrace1)
        varname='T:-1DLON'                            ! T:-1DLON
        call add2list(varname,tvar,ntrace1)
        varname='P'                                   ! P
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for dt/dy (DTDY)
     elseif ( tvar(i).eq.'DTDY' ) then
        varname='T:+1DLAT'                            ! T:+1DLON
        call add2list(varname,tvar,ntrace1)
        varname='T:-1DLAT'                            ! T:-1DLON
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for dt/dp (DTDP)
     elseif ( tvar(i).eq.'DTDP' ) then
        varname='T:+1DP'                              ! T:+1DP
        call add2list(varname,tvar,ntrace1)
        varname='T:-1DP'                              ! T:-1DP
        call add2list(varname,tvar,ntrace1)
        varname='P:+1DP'                              ! P:+1DP
        call add2list(varname,tvar,ntrace1)
        varname='P:-1DP'                              ! P:-1DP
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for dth/dp (DTHDP)
     elseif ( tvar(i).eq.'DTHDP' ) then
        varname='T:+1DP'                              ! T:+1DP
        call add2list(varname,tvar,ntrace1)
        varname='T:-1DP'                              ! T:-1DP
        call add2list(varname,tvar,ntrace1)
        varname='T'                                   ! T
        call add2list(varname,tvar,ntrace1)
        varname='P:+1DP'                              ! P:+1DP
        call add2list(varname,tvar,ntrace1)
        varname='P:-1DP'                              ! P:-1DP
        call add2list(varname,tvar,ntrace1)
        varname='P'                                   ! P
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for squared Brunt Vaiäläa frequency (NSQ)
     elseif ( tvar(i).eq.'NSQ' ) then
        varname='DTHDP'                                ! DTHDP
        call add2list(varname,tvar,ntrace1)
        varname='TH'                                   ! TH
        call add2list(varname,tvar,ntrace1)
        varname='RHO'                                  ! RHO
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for relative vorticity (RELVORT)
     elseif ( tvar(i).eq.'RELVORT' ) then
        varname='U'                                    ! U
        call add2list(varname,tvar,ntrace1)
        varname='DUDY'                                 ! DUDY
        call add2list(varname,tvar,ntrace1)
        varname='DVDX'                                 ! DVDX
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for relative vorticity (ABSVORT)
     elseif ( tvar(i).eq.'ABSVORT' ) then
        varname='U'                                    ! U
        call add2list(varname,tvar,ntrace1)
        varname='DUDY'                                 ! DUDY
        call add2list(varname,tvar,ntrace1)
        varname='DVDX'                                 ! DVDX
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for divergence (DIV)
     elseif ( tvar(i).eq.'DIV' ) then
        varname='V'                                    ! U
        call add2list(varname,tvar,ntrace1)
        varname='DUDX'                                 ! DUDX
        call add2list(varname,tvar,ntrace1)
        varname='DVDY'                                 ! DVDY
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for deformation (DEF)
     elseif ( tvar(i).eq.'DEF' ) then
        call add2list(varname,tvar,ntrace1)
        varname='DUDX'                                 ! DUDX
        call add2list(varname,tvar,ntrace1)
        varname='DVDY'                                 ! DVDY
        call add2list(varname,tvar,ntrace1)
        varname='DUDY'                                 ! DUDY
        call add2list(varname,tvar,ntrace1)
        varname='DVDX'                                 ! DVDX
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for potential vorticity (PV)
     elseif ( tvar(i).eq.'PV' ) then
        varname='ABSVORT'                              ! ABSVORT
        call add2list(varname,tvar,ntrace1)
        varname='DTHDP'                                ! DTHDP
        call add2list(varname,tvar,ntrace1)
        varname='DUDP'                                 ! DUDP
        call add2list(varname,tvar,ntrace1)
        varname='DVDP'                                 ! DVDP
        call add2list(varname,tvar,ntrace1)
        varname='DTHDX'                                ! DTHDX
        call add2list(varname,tvar,ntrace1)
        varname='DTHDY'                                ! DTHDY
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for Richardson number (RI)
     elseif ( tvar(i).eq.'RI' ) then
        varname='DUDP'                                 ! DUDP
        call add2list(varname,tvar,ntrace1)
        varname='DVDP'                                 ! DVDP
        call add2list(varname,tvar,ntrace1)
        varname='NSQ'                                  ! NSQ
        call add2list(varname,tvar,ntrace1)
        varname='RHO'                                  ! RHO
        call add2list(varname,tvar,ntrace1)

     ! Get the dependencies for Ellrod&Knapp's turbulence index (TI)
     elseif ( tvar(i).eq.'TI' ) then
        varname='DEF'                                  ! DEF
        call add2list(varname,tvar,ntrace1)
        varname='DUDP'                                 ! DUDP
        call add2list(varname,tvar,ntrace1)
        varname='DVDP'                                 ! DVDP
        call add2list(varname,tvar,ntrace1)
        varname='RHO'                                  ! RHO
        call add2list(varname,tvar,ntrace1)

     endif

     ! Exit point for handling additional fields
 100   continue
     i = i + 1

  enddo

  ! Save the full variable name (including shift specification)
  do i=1,ncol
     varsint(i)      = varsinp(i)
  enddo
  do i=1,ntrace1
     varsint(i+ncol) = tvar(i)
  enddo

  ! Bojan: check that PV and TH are on trajectory
  if (intmode .eq. 'clustering') then
     pvpos=-1
     thpos=-1
     do i=1,ncol+ntrace1
        if (varsint(i) .eq. 'TH') then
        thpos=i
        print*,'Clustering: Found TH at position:',thpos
        endif
        if (varsint(i) .eq. 'PV') then
        pvpos=i
        print*,'Clustering: Found PV at position:',pvpos
        endif
     enddo
     if (thpos .eq. -1) then
        print*,'WARNING (clustering): Did not find TH'
        stop
     endif
     if (pvpos .eq. -1) then
        print*,'WARNING (clustering): Did not find PV'
        stop
     endif
  endif
  ! Bojan

  ! Split the tracing variables
  do i=1,ntrace0
    call splitvar(tvar(i),shift_val(i),shift_dir(i) )
  enddo


  ! Split the variable name and set flags
  do i=ntrace0+1,ntrace1

     ! Set the scaling factor
     fac(i) = 1.

     ! Set the base variable name, the shift and the direction
     call splitvar(tvar(i),shift_val(i),shift_dir(i) )

     ! Set the prefix of the file name for additional fields
     tfil(i)='*'
     do j=1,n_pvars
        if ( tvar(i).eq.pvars(j) ) tfil(i)=charp
     enddo
     do j=1,n_svars
        if ( tvar(i).eq.svars(j) ) tfil(i)=chars
     enddo

     ! Set the computational flag
     if ( (tvar(i).eq.'P'   ).or. &
          (tvar(i).eq.'PLAY').or. &
          (tvar(i).eq.'PLEV') ) then
        compfl(i) = 0
        tfil(i)   = charp
     elseif ( ( tfil(i).eq.charp ).or.( tfil(i).eq.chars ) ) then
        compfl(i) = 0
     else
        compfl(i) = 1
     endif

  enddo

  ! Check whether the shift modes are supported
  do i=1,ntrace1
     if ( ( shift_dir(i).ne.'nil'     ).and. &
          ( shift_dir(i).ne.'DLON'    ).and. &
          ( shift_dir(i).ne.'DLAT'    ).and. &
          ( shift_dir(i).ne.'DP'      ).and. &
          ( shift_dir(i).ne.'HPA'     ).and. &
          ( shift_dir(i).ne.'HPA(ABS)').and. &
          ( shift_dir(i).ne.'KM(LON)' ).and. &
          ( shift_dir(i).ne.'KM(LAT)' ).and. &
          ( shift_dir(i).ne.'H'       ).and. &
          ( shift_dir(i).ne.'MIN'     ).and. &
          ( shift_dir(i).ne.'INDP'    ) ) then
        print*,' ERROR: shift mode ',trim(shift_dir(i)), ' not supported'
        stop
     endif
  enddo

  ! Write status information
  print*
  print*,'---- COMPLETE TABLE FOR TRACING -------------------------'
  print*
  do i=1,ntrace1
     if ( ( shift_dir(i).ne.'nil' ) ) then
        write(*,'(i4,a4,a8,f10.2,a8,3x,a4,i5)') i,' : ',trim(tvar(i)), &
             shift_val(i),trim(shift_dir(i)),tfil(i),compfl(i)
     else
        write(*,'(i4,a4,a8,10x,8x,3x,a4,i5)') &
             i,' : ',trim(tvar(i)),tfil(i),compfl(i)
     endif
  enddo

  ! --------------------------------------------------------------------
  ! Prepare the internal and output trajectories
  ! --------------------------------------------------------------------

  ! Allocate memory for internal trajectories
  allocate(traint(ntra,ntim,ncol+ntrace1),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array traint   ***'

  ! Copy input to output trajectory
  do i=1,ntra
     do j=1,ntim
        do k=1,ncol
           traint(i,j,k)=trainp(i,j,k)
        enddo
     enddo
  enddo

  ! Set the flags for ready fields/colums - at begin only read-in fields are ready
  do i=1,ncol
     fok(i) = 1
  enddo
  do i=ncol+1,ntrace1
     fok(i) = 0
  enddo

  ! --------------------------------------------------------------------
  ! Trace the fields (fields available on input files)
  ! --------------------------------------------------------------------

  print*
  print*,'---- TRACING FROM PRIMARY AND SECONDARY DATA FILES ------'

  ! Loop over all tracing fields
  do i=1,ntrace1

      ! Skip fields which must be computed (compfl=1), will be handled later
      if (compfl(i).ne.0)  goto 110

      ! Write some status information
      print*
      print*,' Now tracing             : ', trim(tvar(i)),shift_val(i),trim(shift_dir(i)),compfl(i),' ',trim(tfil(i))

      ! Set the flag for ready field/column
      fok(ncol+i) = 1

      ! Reset flags for load manager
      iloaded0 = -1
      iloaded1 = -1

      ! Reset the counter for fields outside domain
      noutside = 0
      err_c1   = 0
      err_c2   = 0
      err_c3   = 0

      ! Loop over all times
      do j=1,ntim

         ! Convert trajectory time from hh.mm to fractional time
         call hhmm2frac(trainp(1,j,1),tfrac)

         ! Shift time if requested
         if ( shift_dir(i).eq.'H' ) then
            tfrac = tfrac + shift_val(i)
         elseif ( shift_dir(i).eq.'MIN' ) then
            tfrac = tfrac + shift_val(i)/60.
         endif

         ! Get the times which are needed
         itime0   = int(abs(tfrac-tstart)/timeinc) + 1
         time0    = tstart + fbflag * real(itime0-1) * timeinc
         itime1   = itime0 + 1
         time1    = time0 + fbflag * timeinc
         if ( itime1.gt.numdat ) then
            itime1 = itime0
            time1  = time0
         endif

         ! Load manager: Check whether itime0 can be copied from itime1
         if ( itime0.eq.iloaded1 ) then
            f3t0     = f3t1
            p3t0     = p3t1
            spt0     = spt1
            iloaded0 = itime0
         endif

         ! Load manager: Check whether itime1 can be copied from itime0
         if ( itime1.eq.iloaded0 ) then
            f3t1     = f3t0
            p3t1     = p3t0
            spt1     = spt0
            iloaded1 = itime1
         endif

         ! Load manager:  Load first time (tracing variable and grid)
         if ( itime0.ne.iloaded0 ) then

            filename = trim(tfil(i))//trim(dat(itime0))
            call frac2hhmm(time0,tload)
            varname  = tvar(i)
            write(*,'(a23,a20,a3,a5,f7.2)') '    ->  loading          : ', trim(filename),'  ',trim(varname),tload
            call input_open (fid,filename)
            call input_wind &
                 (fid,varname,f3t0,tload,stagz,mdv, &
                 xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)
            call input_grid &
                 (fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny, &
                 tload,pollon,pollat,p3t0,spt0,nz,ak,bk,stagz, &
                 timecheck)
            call input_close(fid)

            iloaded0 = itime0

         endif

         ! Load manager: Load second time (tracing variable and grid)
         if ( itime1.ne.iloaded1 ) then

            filename = trim(tfil(i))//trim(dat(itime1))
            call frac2hhmm(time1,tload)
            varname  = tvar(i)
            write(*,'(a23,a20,a3,a5,f7.2)') '    ->  loading          : ', trim(filename),'  ',trim(varname),tload
            call input_open (fid,filename)
            call input_wind &
                 (fid,varname,f3t1,tload,stagz,mdv, &
                 xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)
            call input_grid &
                 (fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny, &
                 tload,pollon,pollat,p3t1,spt1,nz,ak,bk,stagz, &
                 timecheck)
            call input_close(fid)

            iloaded1 = itime1

         endif

         ! Loop over all trajectories
         do k=1,ntra

            ! Set the horizontal position where to interpolate to
            x0       = traint(k,j,2)                          ! Longitude
            y0       = traint(k,j,3)                          ! Latitude

            ! Set the vertical position where to interpolate to
            if ( nz.gt.1 ) then
               p0 = traint(k,j,4)                             ! Pressure (3D tracing)
            else
               p0 = 1050.                                     ! Lowest level (2D tracing)
            endif

            ! Set negative pressures to mdv
            if (p0.lt.0.) then
	        f0 = mdv
                goto 109
            endif

            ! Set the relative time
            call hhmm2frac(traint(k,j,1),tfrac)
            reltpos0 = fbflag * (tfrac-time0)/timeinc

            ! Make adjustments depending on the shift flag
            if ( shift_dir(i).eq.'DLON' ) then                         ! DLON
               x0 = x0 + shift_val(i)

            elseif  ( shift_dir(i).eq.'DLAT' ) then                    ! DLAT
               y0 = y0 + shift_val(i)

            elseif ( shift_dir(i).eq.'KM(LON)' ) then                  ! KM(LON)
               x0 = x0 + shift_val(i)/deg2km * 1./cos(y0*pi180)

            elseif ( shift_dir(i).eq.'KM(LAT)' ) then                  ! KM(LAT)
               y0 = y0 + shift_val(i)/deg2km

            elseif ( shift_dir(i).eq.'HPA' ) then                      ! HPA
               p0 = p0 + shift_val(i)

            elseif ( shift_dir(i).eq.'HPA(ABS)' ) then                 ! HPA(ABS)
               p0 = shift_val(i)

            elseif ( shift_dir(i).eq.'DP' ) then                       ! DP
               call get_index4 (xind,yind,pind,x0,y0,p0,reltpos0, &
                    p3t0,p3t1,spt0,spt1,3,nx,ny,nz,xmin,ymin,dx,dy,mdv)
               pind = pind - shift_val(i)
               p0   = int_index4(p3t0,p3t1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)

            elseif ( shift_dir(i).eq.'INDP' ) then
               p0   = int_index4(p3t0,p3t1,nx,ny,nz,xind,yind,shift_val(i),reltpos0,mdv)

            endif

            ! Handle periodic boundaries in zonal direction
            if ( (x0.gt.xmax).and.(per.ne.0) ) x0 = x0 - 360.
            if ( (x0.lt.xmin).and.(per.ne.0) ) x0 = x0 + 360.

            ! Handle pole problems for hemispheric data (taken from caltra.f)
            if ((hem.eq.1).and.(y0.gt.90.)) then
               print*,'WARNING: y0>90 ',y0,' => setting to 180-y0 ',180.-y0
               y0=180.-y0
               x0=x0+per/2.
            endif
            if ((hem.eq.1).and.(y0.lt.-90.)) then
               print*,'WARNING: y0<-90 ',y0,' => setting to -180-y0 ',-180.-y0
               y0=-180.-y0
               x0=x0+per/2.
            endif

            ! Get the index where to interpolate (x0,y0,p0)
            if ( (abs(x0-mdv).gt.eps).and. (abs(y0-mdv).gt.eps) ) then
               call get_index4 (xind,yind,pind,x0,y0,p0,reltpos0, &
                    p3t0,p3t1,spt0,spt1,3,nx,ny,nz,xmin,ymin,dx,dy,mdv)
            else
               xind = mdv
               yind = mdv
               pind = mdv
            endif

           ! Check if point is within grid (keep indices if ok)
            if ( (xind.ge.1.).and.(xind.le.real(nx)).and. &
                 (yind.ge.1.).and.(yind.le.real(ny)).and. &
                 (pind.ge.1.).and.(pind.le.real(nz)) ) then
		 xind = xind
		 yind = yind
		 pind = pind

            ! Check if pressure is outside, but rest okay => adjust to lowest or highest level
            elseif ( (xind.ge.1.).and.(xind.le.real(nx)).and. (yind.ge.1.).and.(yind.le.real(ny)) ) then ! only vertical problem

               if ( pind.gt.nz ) then ! pressure too low, index too high
                 err_c1 = err_c1 + 1
                 if ( err_c1.lt.10 ) then
                    write(*,*) ' WARNING: pressure too low (pind = ',pind,') => adjusted to highest level (pind=nz.)'
                    print*,'(x0,y0,p0)=',x0,y0,p0
                    pind = real(nz)
                 elseif ( err_c1.eq.10 ) then
                    print*,' WARNING: more pressures too low -> adjusted to highest level '
                    pind = real(nz)
                 else
                    pind = real(nz)
                 endif
                 
               elseif (pind.lt.1.) then ! pressure too high, index too low
                 err_c2 = err_c2 + 1
                 if ( err_c2.lt.10 ) then
                    write(*,*) ' WARNING: pressure too high (pind = ',pind,') => adjusted to lowest level, (pind=1.)'
                    print*,'(x0,y0,p0)=',x0,y0,p0
                    pind = 1.
                 elseif ( err_c2.eq.10 ) then
                    print*,' WARNING: more pressures too high -> adjusted to lowest level '
                    pind = 1.
                 else
                    pind = 1.
                 endif

              endif

            ! Grid point is outside!
            else
               
               err_c3 = err_c3 + 1
               if ( err_c3.lt.10 ) then
                  print*,'ERROR: point is outside grid (horizontally)'
                  print*,'   Trajectory # ',k
                  print*,'   Position     ',x0,y0,p0
                  print*,'  (xind,yind):  ',xind,yind
                  xind          = mdv
                  yind          = mdv
                  pind          = mdv
                  traint(k,j,2) = mdv  
                  traint(k,j,3) = mdv  
                  traint(k,j,4) = mdv  
               elseif ( err_c3.eq.10 ) then
                  print*,'ERROR: more points outside grid (horizontally)'
                  xind          = mdv
                  yind          = mdv
                  pind          = mdv
                  traint(k,j,2) = mdv  
                  traint(k,j,3) = mdv  
                  traint(k,j,4) = mdv  
               else
                  xind          = mdv
                  yind          = mdv
                  pind          = mdv
                  traint(k,j,2) = mdv  
                  traint(k,j,3) = mdv  
                  traint(k,j,4) = mdv  
               endif

            endif

            ! ------------------------ NEAREST mode ------------------------------- 
            ! Interpolate to nearest grid point
            if ( intmode.eq.'nearest') then

               xind = real( nint(xind) )
               yind = real( nint(yind) )
               pind = real( nint(pind) )

               if ( xind.lt.1.  ) xind = 1.
               if ( xind.gt.nx  ) xind = real(nx)
               if ( yind.lt.1.  ) yind = 1.
               if ( yind.gt.ny  ) yind = real(ny)

               if ( pind.lt.1.  ) pind = 1.
               if ( pind.gt.nz  ) pind = real(nz)

               if ( abs(reltpos0).ge.eps ) then
                  print*,'ERROR (nearest): reltpos != 0',reltpos0
                  stop
               endif

               ! interpolate
               f0 = int_index4(f3t0,f3t1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)

            ! ------------------------ end NEAREST mode ------------------------------- 

            ! ------------------------ CLUSTERING mode ------------------------------- Bojan
            elseif (intmode.eq.'clustering') then
               if (varname.ne.'LABEL' ) then
                  print*,'ERROR (clustering): varname is not LABEL'
                  stop
               endif

            ! Get indices of box around the point
            xindb=floor(xind)
            xindf=ceiling(xind)
            yindb=floor(yind)
            yindf=ceiling(yind)
            pindb=floor(pind)
            pindf=ceiling(pind)

            ! Make sure all points are within grid
            if ( xindb.lt.1 ) xindb = 1
            if ( xindf.lt.1 ) xindf = 1
            if ( xindb.gt.nx ) xindb = nx
            if ( xindf.gt.nx ) xindf = nx
            if ( yindb.lt.1 ) yindb = 1
            if ( yindf.lt.1 ) yindf = 1
            if ( yindb.gt.ny ) yindb = ny
            if ( yindf.gt.ny ) yindf = ny
            if ( pindb.lt.1 ) pindb = 1
            if ( pindf.lt.1 ) pindf = 1
            if ( pindb.gt.nz ) pindb = nz
            if ( pindf.gt.nz ) pindf = nz

            ! Shift one point if they are equal
            if ( xindb.eq.xindf ) then
               if ( xindf.eq.nx ) then
                  xindb=nx-1
               else
                  xindf=xindb+1
               endif
            endif
            if ( yindb.eq.yindf ) then
               if ( yindf.eq.ny ) then
                  yindb=ny-1
               else
                  yindf=yindb+1
               endif
            endif
            if ( pindb.eq.pindf ) then
               if ( pindf.eq.nz ) then
                  pindb=nz-1
               else
                  pindf=pindb+1
               endif
            endif
            ! Give warnings and stop if problems occur
            if ( xindb.eq.xindf ) then
               print*,'ERROR (clustering): xindb=xindf'
               print*,xind,xindb,xindf
               stop
            endif
            if ( yindb.eq.yindf ) then
               print*,'ERROR (clustering): yindb=yindf'
               print*,yind,yindb,yindf
               stop
            endif
            if ( pindb.eq.pindf ) then
               print*,'ERROR (clustering): pindb=pindf'
               print*,pind,pindb,pindf
               stop
            endif
            if ( ( xindb.lt.1 ).or.( xindf.gt.nx ) ) then
               print*,'ERROR (clustering): xindb/f outside'
               print*,xind,xindb,xindf
               stop
            endif
            if ( ( yindb.lt.1 ).or.( yindf.gt.ny ) ) then
               print*,'ERROR (clustering): yindb/f outside'
               print*,yind,yindb,yindf
               stop
            endif
            if ( ( pindb.lt.1 ).or.( pindf.gt.nz ) ) then
               print*,'ERROR (clustering): pindb/f outside'
               print*,pind,pindb,pindf
               stop
            endif
            if ( abs(reltpos0).ge.eps ) then
               print*,'ERROR (clustering): reltpos != 0',reltpos0
               stop
            endif

            ! Get Value in Box
            lblcount=(/0,0,0,0,0/)

            lbl(1) = f3t0( xindb + nx*(yindb-1) + nx*ny*(pindb-1) )
            lbl(2) = f3t0( xindf + nx*(yindb-1) + nx*ny*(pindb-1) )
            lbl(3) = f3t0( xindb + nx*(yindf-1) + nx*ny*(pindb-1) )
            lbl(4) = f3t0( xindf + nx*(yindf-1) + nx*ny*(pindb-1) )
            lbl(5) = f3t0( xindb + nx*(yindb-1) + nx*ny*(pindf-1) )
            lbl(6) = f3t0( xindf + nx*(yindb-1) + nx*ny*(pindf-1) )
            lbl(7) = f3t0( xindb + nx*(yindf-1) + nx*ny*(pindf-1) )
            lbl(8) = f3t0( xindf + nx*(yindf-1) + nx*ny*(pindf-1) )

            ! Count the number of times every label appears
            do lci=1,5
               do lcj=1,8
                  if ( abs(lbl(lcj)-lci).lt.eps ) then
                     lblcount(lci)=lblcount(lci)+1
                  endif
               enddo
            enddo

            ! Set to -9 to detect if no label was assigned in the end
            f0=-9

            ! Stratosphere (PV)
            if ( abs(traint(k,j,pvpos)) .ge. tropo_pv ) then
               if ( (lblcount(2).ge.lblcount(3)).and. (lblcount(2).ge.lblcount(5)) ) then
                  f0=2
               elseif ( lblcount(3).ge.lblcount(5) ) then
                  f0=3
               elseif ( lblcount(5).gt.lblcount(3) ) then
                  f0=5
               endif
            endif

            ! Troposphere (PV)
            if ( abs(traint(k,j,pvpos)) .lt. tropo_pv ) then
               if ( lblcount(1).ge.lblcount(4) ) then
               f0=1
               elseif ( lblcount(4).gt.lblcount(1) ) then
               f0=4
               endif
            endif

            ! Stratosphere (TH)
            if ( traint(k,j,thpos) .ge. tropo_th ) then
               f0=2
            endif

            if (f0.eq.-9) then
               print*,'ERROR (Clustering): No label assigned!'
               stop
            endif
            ! ------------------------ end CLUSTERING mode -------------------------------

            ! ------------------------ CIRCLE modes ------------------------------- Bojan
            ! elseif (not clustering but one of the possible circle modes)
            elseif ( (intmode.eq.'circle_avg') .or. (intmode.eq.'circle_min') .or. (intmode.eq.'circle_max') ) then

            ! reset arrays for this point
            connect=0
            stackx=0
            stacky=0
            circlelon=0
            circlelat=0
            circlef=0
            circlearea=0

            ! Get indices of one coarse grid point within search radius (nint=round to next integer)
            if ( sdis(x0,y0,longrid(nint(xind)),latgrid(nint(yind))) .gt. radius) then
               print*,'ERROR (circle): Search radius is too small... (1). r =',radius
               print*,'Distance to nint grid point (minimum search radius)=',sdis(x0,y0,longrid(nint(xind)),latgrid(nint(yind)))
               stop
            endif
 
            ! Initialize stack with nint(xind),nint(yind)
            kst=0 ! counts the number of points in circle
            stackx(1)=nint(xind)
            stacky(1)=nint(yind)
            sp=1 ! stack counter
            do while (sp.ne.0)

            ! Get an element from stack
             ist=stackx(sp)
             jst=stacky(sp)
             sp=sp-1

            ! Get distance from reference point
             dist=sdis(x0,y0,longrid(ist),latgrid(jst))

            ! Check whether distance is smaller than search radius: connected
             if (dist.lt.radius) then

            ! Increase total stack index
              kst=kst+1
              circlelon(kst)=longrid(ist)
              circlelat(kst)=latgrid(jst)
  
            ! Interpolate field to position of point (interpolation in time!) 
              circlef(kst) = int_index4(f3t0,f3t1,nx,ny,nz,real(ist),real(jst),pind,reltpos0,mdv)

            ! Calculate area of point (for circle_avg mode only)
              if ( intmode .eq. 'circle_avg' ) then
                 circlearea(kst) = pir/(nx-1)*(sin(pi180*abs(circlelat(kst)))-sin(pi180*(abs(circlelat(kst))-dy)))
              endif

            ! Mark this point as visited
              connect(ist,jst)=1

            ! Get coordinates of neighbouring points and implement periodicity
              mr=ist+1
              if (mr.gt.nx) mr=1
              ml=ist-1
              if (ml.lt.1) ml=nx
              nu=jst+1
              if (nu.gt.ny) nu=ny
              nd=jst-1
              if (nd.lt.1) nd=1

            ! Update stack with neighbouring points
              if (connect(mr,jst).ne. 1) then
                 connect(mr,jst)=1
                 sp=sp+1
                 stackx(sp)=mr
                 stacky(sp)=jst
              endif
              if (connect(ml,jst).ne. 1) then
                 connect(ml,jst)=1
                 sp=sp+1
                 stackx(sp)=ml
                 stacky(sp)=jst
              endif
              if (connect(ist,nd).ne. 1) then
                 connect(ist,nd)=1
                 sp=sp+1
                 stackx(sp)=ist
                 stacky(sp)=nd
              endif
              if (connect(ist,nu).ne. 1) then
                 connect(ist,nu)=1
                 sp=sp+1
                 stackx(sp)=ist
                 stacky(sp)=nu
              endif

             endif ! endif radius is smaller => end of updating stack

            end do ! end working on stack 

            if (kst.ge.1) then
               ! Choose output depending on intmode
               if ( intmode .eq. 'circle_avg' ) then
                  ! calculate area-weighted average of f in circle
                  circlesum=0.
                  do l=1,kst
                     circlesum=circlesum+circlef(l)*circlearea(l)
                  enddo
                  circleavg=circlesum/sum(circlearea(1:kst))
                  !print*,'area-weighted average of f in circle=',circleavg
                  f0=circleavg
               elseif ( intmode .eq. 'circle_min' ) then
                  ! calculate minimum in circle
                  circlemin=circlef(1)
                  do l=1,kst
                     if (circlef(l) .lt. circlemin) then
                        circlemin=circlef(l)
                     endif
                  enddo
                  !print*,'minimum of f in circle=',circlemin       
                  f0=circlemin
               elseif ( intmode .eq. 'circle_max' ) then             
                  ! calculate maximum in circle
                  circlemax=circlef(1)
                  do l=1,kst
                     if (circlef(l) .gt. circlemax) then
                        circlemax=circlef(l)
                     endif
                  enddo
                  !print*,'maximum of f in circle=',circlemax
                  f0=circlemax
               else 
                  print*,'ERROR (circle): intmode not valid!'
                  stop
               endif
            else
               print*,'ERROR (circle): Search radius is too small... (2). r =',radius
               stop
            endif

            ! ------------------------ end CIRCLE modes -------------------------------

            ! ------------------------ NORMAL mode -------------------------------
            else ! not clustering nor circle: NORMAL mode

            ! Check if point is within grid
!            if ( (xind.ge.1.).and.(xind.le.real(nx)).and. &
!                 (yind.ge.1.).and.(yind.le.real(ny)).and. &
!                 (pind.ge.1.).and.(pind.le.real(nz)) ) then
!
            ! Do the interpolation: everthing is ok
               f0 = int_index4(f3t0,f3t1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)

!            ! Check if pressure is outside, but rest okay: adjust to lowest or highest level
!            elseif ( (xind.ge.1.).and.(xind.le.real(nx)).and. (yind.ge.1.).and.(yind.le.real(ny)) ) then ! only vertical problem
!               if ( pind.gt.nz ) then ! pressure too low, index too high
!                 pind = real(nz)
!                 print*,' Warning: pressure too low -> adjusted to highest level, pind=nz.'
!                 print*,'(x0,y0,p0)=',x0,y0,p0
!               elseif (pind.lt.1.) then ! pressure too high, index too low
!                 pind = 1.
!                 print*,' Warning: pressure too high -> adjusted to lowest level, pind=1.'
!                 print*,'(x0,y0,p0)=',x0,y0,p0
!              endif
!              f0 = int_index4(f3t0,f3t1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)

!            ! Less than 10 outside
!            elseif (noutside.lt.10) then
!               print*,' ',trim(tvar(i)),' @ ',x0,y0,p0,'outside'
!               f0       = mdv
!               noutside = noutside + 1
!
!            ! More than 10 outside
!            elseif (noutside.eq.10) then
!               print*,' ...more than 10 outside...'
!               f0       = mdv
!               noutside = noutside + 1

!            ! Else (not everything okay and also not 'tolerated cases') set to missing data
!            else
!               f0       = mdv
!            endif

            ! ------------------------ end NORMAL mode -------------------------------
            endif ! end if nearest case

           ! Exit for loop over all trajectories and times -save interpolated value
 109        continue

            ! Save the new field
            if ( abs(f0-mdv).gt.eps) then
               traint(k,j,ncol+i) = f0
            else
               traint(k,j,ncol+i) = mdv
            endif

         enddo ! end loop over all trajectories

      enddo ! end loop over all times

      ! Exit point for loop over all tracing variables
 110   continue

   enddo ! end loop over all variables

  ! --------------------------------------------------------------------
  ! Calculate additional fields along the trajectories
  ! --------------------------------------------------------------------

  print*
  print*,'---- CALCULATE ADDITIONAL FIELDS FROM TRAJECTORY TABLE --'

  ! Loop over all tracing fields
  do i=ntrace1,1,-1

      ! Skip fields which must not be computed (compfl=0)
      if (compfl(i).eq.0) goto 120

      ! Write some status information
      print*
      write(*,'(a10,f10.2,a5,i3,3x,a2)') &
           trim(tvar(i)),shift_val(i),trim(shift_dir(i)),compfl(i),trim(tfil(i))

      ! Loop over trajectories and times
      do j=1,ntra
      do k=1,ntim

         ! Potential temperature (TH)
         if  ( varsint(i+ncol).eq.'TH' ) then

            varname='T'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='p'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)

            call calc_TH  (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2) )

         ! Density (RHO)
         elseif  ( varsint(i+ncol).eq.'RHO' ) then

            varname='T'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='p'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)

            call calc_RHO (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2) )

         ! Relative humidity (RH)
         elseif  ( varsint(i+ncol).eq.'RH' ) then

            varname='T'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='p'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='Q'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)

            call calc_RH (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3) )

         ! Equivalent potential temperature (THE)
         elseif  ( varsint(i+ncol).eq.'THE' ) then

            varname='T'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='p'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='Q'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)

            call calc_THE (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3) )

         ! Latent heating rate (LHR)
         elseif  ( varsint(i+ncol).eq.'LHR' ) then

            varname='T'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='p'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='Q'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)
            varname='OMEGA'
            call list2ind (ind4,varname,varsint,fok,ncol+ntrace1)
            varname='RH'
            call list2ind (ind5,varname,varsint,fok,ncol+ntrace1)

            call calc_LHR (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3),traint(j,k,ind4),traint(j,k,ind5) )

         ! Wind speed (VEL)
         elseif  ( varsint(i+ncol).eq.'VEL' ) then

            varname='U'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='V'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)

            call calc_VEL (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2) )

         ! Wind direction (DIR)
         elseif  ( varsint(i+ncol).eq.'DIR' ) then

            varname='U'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='V'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)

            call calc_DIR (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2) )

         ! Zonal gradient of U (DUDX)
         elseif  ( varsint(i+ncol).eq.'DUDX' ) then

            varname='U:+1DLON'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='U:-1DLON'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='lat'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)

            call calc_DUDX (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3) )

         ! Zonal gradient of V (DVDX)
         elseif  ( varsint(i+ncol).eq.'DVDX' ) then

            varname='V:+1DLON'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='V:-1DLON'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='lat'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)

            call calc_DVDX (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3) )

         ! Zonal gradient of T (DTDX)
         elseif  ( varsint(i+ncol).eq.'DVDX' ) then

            varname='T:+1DLON'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='T:-1DLON'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='lat'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)

            call calc_DTDX (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3) )

         ! Zonal gradient of TH (DTHDX)
         elseif  ( varsint(i+ncol).eq.'DTHDX' ) then

            varname='T:+1DLON'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='T:-1DLON'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='P'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)
            varname='lat'
            call list2ind (ind4,varname,varsint,fok,ncol+ntrace1)

            call calc_DTHDX (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3),traint(j,k,ind4) )

         ! Meridional gradient of U (DUDY)
         elseif  ( varsint(i+ncol).eq.'DUDY' ) then

            varname='U:+1DLAT'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='U:-1DLAT'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)

            call calc_DUDY (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2) )

         ! Meridional gradient of V (DVDY)
         elseif  ( varsint(i+ncol).eq.'DVDY' ) then

            varname='V:+1DLAT'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='V:-1DLAT'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)

            call calc_DVDY (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2) )

         ! Meridional gradient of T (DTDY)
         elseif  ( varsint(i+ncol).eq.'DTDY' ) then

            varname='T:+1DLAT'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='T:-1DLAT'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)

            call calc_DTDY (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2) )

         ! Meridional gradient of TH (DTHDY)
         elseif  ( varsint(i+ncol).eq.'DTHDY' ) then

            varname='T:+1DLAT'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='T:-1DLAT'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='P'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)

            call calc_DTDY (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3) )


         ! Vertical wind shear DU/DP (DUDP)
         elseif  ( varsint(i+ncol).eq.'DUDP' ) then

            varname='U:+1DP'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='U:-1DP'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='P:+1DP'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)
            varname='P:-1DP'
            call list2ind (ind4,varname,varsint,fok,ncol+ntrace1)

            call calc_DUDP (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3),traint(j,k,ind4) )

         ! Vertical wind shear DV/DP (DVDP)
         elseif  ( varsint(i+ncol).eq.'DVDP' ) then

            varname='V:+1DP'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='V:-1DP'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='P:+1DP'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)
            varname='P:-1DP'
            call list2ind (ind4,varname,varsint,fok,ncol+ntrace1)

            call calc_DVDP (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3),traint(j,k,ind4) )

         ! Vertical derivative of T (DTDP)
         elseif  ( varsint(i+ncol).eq.'DTDP' ) then

            varname='T:+1DP'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='T:-1DP'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='P:+1DP'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)
            varname='P:-1DP'
            call list2ind (ind4,varname,varsint,fok,ncol+ntrace1)

            call calc_DTDP (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3),traint(j,k,ind4) )

         ! Vertical derivative of TH (DTHDP)
         elseif  ( varsint(i+ncol).eq.'DTHDP' ) then

            varname='T:+1DP'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='T:-1DP'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='P:+1DP'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)
            varname='P:-1DP'
            call list2ind (ind4,varname,varsint,fok,ncol+ntrace1)
            varname='P'
            call list2ind (ind5,varname,varsint,fok,ncol+ntrace1)
            varname='T'
            call list2ind (ind6,varname,varsint,fok,ncol+ntrace1)

            call calc_DTHDP (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3),traint(j,k,ind4),traint(j,k,ind5),traint(j,k,ind6) )

         ! Squared Brunt-Vaisäla frequency (NSQ)
         elseif  ( varsint(i+ncol).eq.'NSQ' ) then

            varname='DTHDP'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='TH'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='RHO'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)

            call calc_NSQ (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3))

         ! Relative vorticity (RELVORT)
         elseif  ( varsint(i+ncol).eq.'RELVORT' ) then

            varname='DUDY'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='DVDX'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='U'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)
            varname='lat'
            call list2ind (ind4,varname,varsint,fok,ncol+ntrace1)

            call calc_RELVORT (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3),traint(j,k,ind4))

         ! Absolute vorticity (ABSVORT)
         elseif  ( varsint(i+ncol).eq.'ABSVORT' ) then

            varname='DUDY'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='DVDX'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='U'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)
            varname='lat'
            call list2ind (ind4,varname,varsint,fok,ncol+ntrace1)

            call calc_ABSVORT (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3),traint(j,k,ind4))

         ! Divergence (DIV)
         elseif  ( varsint(i+ncol).eq.'DIV' ) then

            varname='DUDX'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='DVDY'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='V'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)
            varname='lat'
            call list2ind (ind4,varname,varsint,fok,ncol+ntrace1)

            call calc_DIV (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3),traint(j,k,ind4))

         ! Deformation (DEF)
         elseif  ( varsint(i+ncol).eq.'DEF' ) then

            varname='DUDX'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='DVDX'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='DUDY'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)
            varname='DVDY'
            call list2ind (ind4,varname,varsint,fok,ncol+ntrace1)

            call calc_DEF (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3),traint(j,k,ind4))

         ! Potential Vorticity (PV)
         elseif  ( varsint(i+ncol).eq.'PV' ) then

            varname='ABSVORT'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='DTHDP'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='DUDP'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)
            varname='DVDP'
            call list2ind (ind4,varname,varsint,fok,ncol+ntrace1)
            varname='DTHDX'
            call list2ind (ind5,varname,varsint,fok,ncol+ntrace1)
            varname='DTHDY'
            call list2ind (ind6,varname,varsint,fok,ncol+ntrace1)

            call calc_PV (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3),traint(j,k,ind4),traint(j,k,ind5),traint(j,k,ind6) )

         ! Richardson number (RI)
         elseif  ( varsint(i+ncol).eq.'RI' ) then

            varname='DUDP'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='DVDP'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='NSQ'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)
            varname='RHO'
            call list2ind (ind4,varname,varsint,fok,ncol+ntrace1)

            call calc_RI (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3),traint(j,k,ind4) )

         ! Ellrod and Knapp's turbulence idicator (TI)
         elseif  ( varsint(i+ncol).eq.'TI' ) then

            varname='DEF'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='DUDP'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)
            varname='DVDP'
            call list2ind (ind3,varname,varsint,fok,ncol+ntrace1)
            varname='RHO'
            call list2ind (ind4,varname,varsint,fok,ncol+ntrace1)

            call calc_TI (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,k,ind3),traint(j,k,ind4) )

         ! Spherical distance from starting position (DIST0)
         elseif  ( varsint(i+ncol).eq.'DIST0' ) then

            varname='lon'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='lat'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)

            call calc_DIST0 (traint(j,k,ncol+i), traint(j,k,ind1), &
                 traint(j,k,ind2),traint(j,1,ind1),traint(j,1,ind2) )

         ! Spherical distance length of trajectory (DIST)
         elseif  ( varsint(i+ncol).eq.'DIST' ) then

            varname='lon'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='lat'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)

            if ( k.eq.1 ) then
               traint(j,k,ncol+i) = 0.
            else
               call calc_DIST0 (delta, traint(j,k  ,ind1), &
                   traint(j,k  ,ind2),traint(j,k-1,ind1),traint(j,k-1,ind2) )
               traint(j,k,ncol+i) = traint(j,k-1,ncol+i) + delta
            endif

         ! Heading of the trajectory (HEAD)
         elseif  ( varsint(i+ncol).eq.'HEAD' ) then

            varname='lon'
            call list2ind (ind1,varname,varsint,fok,ncol+ntrace1)
            varname='lat'
            call list2ind (ind2,varname,varsint,fok,ncol+ntrace1)

            if (k.eq.ntim) then
               traint(j,k,ncol+i) = mdv
            else
               call calc_HEAD (traint(j,k,ncol+i), &
                    traint(j,k  ,ind1),traint(j,k  ,ind2),traint(j,k+1,ind1),traint(j,k+1,ind2) )
            endif


         ! Invalid tracing variable
         else

            print*,' ERROR: invalid tracing variable ', trim(varsint(i+ncol))
            stop


         endif

      ! End loop over all trajectories and times
      enddo
      enddo

      ! Set the flag for a ready field/column
      fok(ncol+i) = 1


      ! Exit point for loop over all tracing fields
 120   continue

   enddo

  ! --------------------------------------------------------------------
  ! Write output to output trajectory file
  ! --------------------------------------------------------------------

  ! Write status information
  print*
  print*,'---- WRITE OUTPUT TRAJECTORIES --------------------------'
  print*

  ! Allocate memory for internal trajectories
  allocate(traout(ntra,ntim,ncol+ntrace0),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array traout   ***'

  ! Copy input to output trajectory (apply scaling of output)
  do i=1,ntra
     do j=1,ntim
        do k=1,ncol+ntrace0
           if ( k.le.ncol ) then
              traout(i,j,k) = traint(i,j,k)
           elseif ( abs(traint(i,j,k)-mdv).gt.eps ) then
              traout(i,j,k) = fac(k-ncol) * traint(i,j,k)
           else
              traout(i,j,k) = mdv
           endif
        enddo
     enddo
  enddo

  ! Set the variable names for output trajectory
  do i=1,ncol+ntrace0
     varsout(i)      = varsint(i)
  enddo

  ! Write trajectories
  call wopen_tra(fod,outfile,ntra,ntim,ncol+ntrace0, &
       reftime,varsout,outmode)
  call write_tra(fod,traout ,ntra,ntim,ncol+ntrace0,outmode)
  call close_tra(fod,outmode)

  ! Write some status information, and end of program message
  print*
  print*,'---- STATUS INFORMATION --------------------------------'
  print*
  print*,' ok'
  print*
  print*,'              *** END OF PROGRAM TRACE ***'
  print*,'========================================================='


end program trace



! ******************************************************************
! * SUBROUTINE SECTION                                             *
! ******************************************************************

! ------------------------------------------------------------------
! Add a variable to the list if not yet included in this list
! ------------------------------------------------------------------

subroutine add2list (varname,list,nlist)

  implicit none

  ! Declaration of subroutine parameters
  character(len=80) :: varname
  character(len=80) :: list(200)
  integer :: nlist

  ! Auxiliray variables
  integer :: i,j
  integer :: isok

  ! Expand the list, if necessary
  isok = 0
  do i=1,nlist
     if ( list(i).eq.varname ) isok = 1
  enddo
  if ( isok.eq.0 ) then
     nlist       = nlist + 1
     list(nlist) = varname
  endif

  ! Check for too large number of fields
  if ( nlist.ge.200) then
     print*,' ERROR: too many additional fields for tracing ...'
     stop
  endif

end subroutine add2list


! ------------------------------------------------------------------
! Get the index of a variable in the list
! ------------------------------------------------------------------

subroutine list2ind (ind,varname,list,fok,nlist)

  implicit none

  ! Declaration of subroutine parameters
  integer :: ind
  character(len=80) :: varname
  character(len=80) :: list(200)
  integer :: fok(200)
  integer :: nlist

  ! Auxiliray variables
  integer :: i,j
  integer :: isok

  ! Get the index - error message if not found
  ind = 0
  do i=1,nlist
     if ( varname.eq.list(i) ) then
        ind = i
        goto 100
     endif
  enddo

  if ( ind.eq.0) then
     print*
     print*,' ERROR: cannot find ',trim(varname),' in list ...'
     do i=1,nlist
        print*,i,trim(list(i))
     enddo
     print*
     stop
  endif

  ! Exit point
 100   continue

  ! Check whether the field/column is ready
  if ( fok(ind).eq.0 ) then
     print*
     print*,' ERROR: unresolved dependence : ',trim(list(ind))
     print*
     stop
  endif

end subroutine list2ind


! ------------------------------------------------------------------
! Split the variable name into name, shift and direction
! ------------------------------------------------------------------

subroutine splitvar (tvar,shift_val,shift_dir)

  implicit none

  ! Declaration of subroutine parameters
  character(len=80) :: tvar
  real :: shift_val
  character(len=80) :: shift_dir

  ! Auxiliary variables
  integer :: i,j
  integer :: icolon,inumber
  character(len=80) :: name
  character :: ch
  integer      isabsval

  ! Save variable name
  name = tvar

  ! Search for colon
  icolon=0
  do i=1,80
     if ( (name(i:i).eq.':').and.(icolon.ne.0) ) goto 100
     if ( (name(i:i).eq.':').and.(icolon.eq.0) ) icolon=i
  enddo

  ! If there is a colon, split the variable name
  if ( icolon.ne.0 ) then

     tvar = name(1:(icolon-1))

     ! Get the index for number
     do i=icolon+1,80
        ch = name(i:i)
        if ( ( ch.ne.'0' ).and. ( ch.ne.'1' ).and.( ch.ne.'2' ).and. &
             ( ch.ne.'3' ).and. ( ch.ne.'4' ).and.( ch.ne.'5' ).and. &
             ( ch.ne.'6' ).and. ( ch.ne.'7' ).and.( ch.ne.'8' ).and. &
             ( ch.ne.'9' ).and. ( ch.ne.'+' ).and.( ch.ne.'-' ).and. &
             ( ch.ne.'.' ).and. ( ch.ne.' ' )  ) then
           inumber = i
           exit
        endif
     enddo

     ! Get the number
     read(name( (icolon+1):(inumber-1) ),*) shift_val

     ! Decide whether it is a shift relatiev to trajectory or absolute value
     ! If the number starts with + or -, it is relative to the trajectory
     isabsval = 1
     do i=icolon+1,inumber-1
       ch = name(i:i)
       if ( (ch.eq.'+').or.(ch.eq.'-') ) isabsval = 0
     enddo

     ! Get the unit/shift axis
     shift_dir = name(inumber:80)
     if ( isabsval.eq.1 ) then
       shift_dir=trim(shift_dir)//'(ABS)'
     endif

  else

     shift_dir = 'nil'
     shift_val = 0.

  endif

  return

  ! Error handling
 100   continue

  print*,' ERROR: cannot split variable name ',trim(tvar)
  stop

end subroutine splitvar
