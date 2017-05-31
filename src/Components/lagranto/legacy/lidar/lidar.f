      PROGRAM trace

C     ********************************************************************
C     *                                                                  *
C     * Pseudo-lidar plots along trajectories                             *
C     *                                                                  *
C     * Heini Wernli       first version:       April 1993               *
C     * Michael Sprenger   major upgrade:       2008-2009                *
C     *                                                                  *
C     ********************************************************************

      implicit none
      
c     --------------------------------------------------------------------
c     Declaration of parameters
c     --------------------------------------------------------------------

c     Maximum number of levels for input files
      integer   nlevmax
      parameter (nlevmax=100)

c     Maximum number of input files (dates, length of trajectories)
      integer   ndatmax
      parameter (ndatmax=5000)

c     Numerical epsilon (for float comparison)
      real      eps
      parameter (eps=0.001)

c     Conversion factors
      real      pi180                                   ! deg -> rad
      parameter (pi180=3.14159/180.)
      real      deg2km                                  ! deg -> km (at equator)
      parameter (deg2km=111.)

c     Prefix for primary and secondary fields
      character charp
      character chars
      parameter (charp='P')
      parameter (chars='S')

c     --------------------------------------------------------------------
c     Declaration of variables
c     --------------------------------------------------------------------

c     Input and output format for trajectories (see iotra.f)
      integer   inpmode

c     Input parameters
      character*80                           inpfile         ! Input trajectory file
      character*80                           outfile         ! Output netCDF file
      character*80                           outmode         ! Output mode (sum,mean)
      integer                                ntra            ! Number of trajectories
      integer                                ncol            ! Number of columns (including time, lon, lat, p)
      integer                                ntim            ! Number of times per trajectory
      integer                                ntrace0         ! Number of trace variables
      character*80                           tvar(200)       ! Tracing variable name (only the variable)
      character*1                            tfil(200)       ! Filename prefix 
      real                                   fac(200)        ! Scaling factor 
      integer                                compfl(200)     ! Computation flag (1=compute)
      integer                                numdat          ! Number of input files
      character*11                           dat(ndatmax)    ! Dates of input files
      real                                   timeinc         ! Time increment between input files
      real                                   tst             ! Time shift of start relative to first data file
      real                                   ten             ! Time shift of end relatiev to first data file  
      character*20                           startdate       ! First time/date on trajectory
      character*20                           enddate         ! Last time/date on trajectory
      character*80                           timecheck       ! Either 'yes' or 'no'
      character*80                           intmode         ! Interpolation mode ('normal', 'nearest')
	  real                                   pmin,pmax       ! Pressure range for output grid
	  integer                                npre            ! Number of pressure levels in output grid
	  character*80                           centering       ! Centering around trajectory position ('yes','no')
	  character*80                       direction       ! Direction of lidar (vertical,lat,lon,normal)
      character*80                           dumpcoord       ! Dumping coordinates ('yes','no')

c     Trajectories
      real,allocatable, dimension (:,:,:) :: trainp                ! Input trajectories (ntra,ntim,ncol)
      integer                                reftime(6)            ! Reference date
      character*80                           varsinp(100)          ! Field names for input trajectory
      integer                                fid,fod               ! File identifier for inp and out trajectories
      real                                   x0_tra,y0_tra,p0_tra  ! Position of air parcel (physical space)
      real                                   reltpos0              ! Relative time of air parcel
      real                                   xind,yind,pind        ! Position of air parcel (grid space)
      integer                                fbflag                ! Flag for forward (1) or backward (-1) trajectories

c     Meteorological fields from input file
      real,allocatable, dimension (:)     :: spt0,spt1       ! Surface pressure
      real,allocatable, dimension (:)     :: p3t0,p3t1       ! 3d-pressure 
      real,allocatable, dimension (:)     :: f3t0,f3t1       ! 3d field for tracing 
      character*80                           svars(100)      ! List of variables on S file
      character*80                           pvars(100)      ! List of variables on P file
      integer                                n_svars         ! Number of variables on S file
      integer                                n_pvars         ! Number of variables on P file
      
c     Input grid description
      real                                   pollon,pollat   ! Longitude/latitude of pole
      real                                   ak(100)         ! Vertical layers and levels
      real                                   bk(100) 
      real                                   xmin,xmax       ! Zonal grid extension
      real                                   ymin,ymax       ! Meridional grid extension
      integer                                nx,ny,nz        ! Grid dimensions
      real                                   dx,dy           ! Horizontal grid resolution
      integer                                hem             ! Flag for hemispheric domain
      integer                                per             ! Flag for periodic domain
      real                                   stagz           ! Vertical staggering
      real                                   mdv             ! Missing data value

c	  Output grid  and fields
      real                                   levels(1000)    ! Ouput levels
      real                                   times (5000)    ! Output times
	  real,allocatable, dimension (:,:)   :: out_pos         ! Position of trajectories
	  real,allocatable, dimension (:,:)   :: out_val         ! Output lidar field
	  real,allocatable, dimension (:,:)   :: out_cnt         ! # output lidar sum ups


c     Auxiliary variables
      integer                                i,j,k,l,n,m
      real                                   rd
      character*80                           filename
      real                                   time0,time1,reltpos
      integer                                itime0,itime1
      integer                                stat
      real                                   tstart
      integer                                iloaded0,iloaded1
      real                                   f0
      real                                   frac
      real                                   tload,tfrac
      integer                                isok
      character                              ch
      integer                                ind
      integer                                ind1,ind2,ind3,ind4,ind5
      integer                                ind6,ind7,ind8,ind9,ind0
      integer                                noutside
      real                                   delta
      integer                                itrace0
      character*80                           string
      character*80                           cdfname
      character*80                           varname
      real                                   time
      character*80                           longname
      character*80                           unit
      integer                                ind_time
      integer                                ind_pre
      real                                   rlat,rlon
      real                                   x0,y0,p0
      real                                   vx0,vy0,vx1,vy1
      real                                   rotation,lon,lat

c     Externals 
      real                                   int_index4
      external                               int_index4

c     --------------------------------------------------------------------
c     Start of program, Read parameters, get grid parameters
c     --------------------------------------------------------------------

c     Write start message
      print*,'========================================================='
      print*,'              *** START OF PROGRAM LIDAR ***'
      print*

c     Read parameters
      open(10,file='trace.param')
       read(10,*) inpfile
       read(10,*) outfile
       read(10,*) outmode
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
       read(10,*) pmin,pmax,npre
       read(10,*) centering
       read(10,*) direction
       read(10,*) dumpcoord
      close(10)

c     Check that the direction is ok
      if ( ( direction.ne.'vertical' ).and.
     >     ( direction.ne.'lat'      ).and.
     >     ( direction.ne.'lon'      ).and.
     >     ( direction.ne.'normal'   ) ) 
     >then 
         print*,' ERROR: invalid direction ',trim(direction)
         stop
      endif

c     Remove commented tracing fields
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

c     Set the formats of the input  files
      call mode_tra(inpmode,inpfile)
      if (inpmode.eq.-1) inpmode=1

C     Convert time shifts <tst,ten> from <hh.mm> into fractional time
      call hhmm2frac(tst,frac)
      tst = frac
      call hhmm2frac(ten,frac)
      ten = frac

c     Set the time for the first data file (depending on forward/backward mode)
      if (fbflag.eq.1) then
        tstart = -tst
      else
        tstart = tst
      endif

c     Read the constant grid parameters (nx,ny,nz,xmin,xmax,ymin,ymax,pollon,pollat)
c     The negative <-fid> of the file identifier is used as a flag for parameter retrieval  
      filename = charp//dat(1)
      varname  = tvar(1)
      call input_open (fid,filename)
      call input_grid (-fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >                 tstart,pollon,pollat,rd,rd,nz,rd,rd,rd,timecheck)
      call input_close(fid)

C     Allocate memory for some meteorological arrays
      allocate(spt0(nx*ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array spt0 ***'   ! Surface pressure
      allocate(spt1(nx*ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array spt1 ***'
      allocate(p3t0(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array p3t0 ***'   ! Pressure
      allocate(p3t1(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array p3t1 ***'
      allocate(f3t0(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array p3t0 ***'   ! Lidar field
      allocate(f3t1(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array p3t1 ***'

c	  Allocate memory for output field
	  allocate(out_pos(ntim,npre),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array out_pos ***'
	  allocate(out_val(ntim,npre),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array out_val ***'
	  allocate(out_cnt(ntim,npre),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array out_cnt ***'

C     Get memory for trajectory arrays
      allocate(trainp(ntra,ntim,ncol),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra      ***' 

c     Set the flags for periodic domains
      if ( abs(xmax-xmin-360.).lt.eps ) then
         per = 1
      elseif ( abs(xmax-xmin-360.+dx).lt.eps ) then
         per = 2
      else
         per = 0
      endif

C     Set logical flag for periodic data set (hemispheric or not)
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

c     Write some status information
      print*,'---- INPUT PARAMETERS -----------------------------------'
      print*
      print*,'  Input trajectory file  : ',trim(inpfile)
      print*,'  Format of input file   : ',inpmode
      print*,'  Output netCDF    file  : ',trim(outfile)
      print*,'  Format of output file  : ',trim(outmode)
      print*,'  Forward/backward       : ',fbflag
      print*,'  #tra                   : ',ntra
      print*,'  #col                   : ',ncol
      print*,'  #tim                   : ',ntim
      print*,'  No time check          : ',trim(timecheck)
      print*,'  Interpolation mode     : ',trim(intmode)
      do i=1,ntrace0
         if (compfl(i).eq.0) then
            print*,'  Tracing field          : ',
     >                 trim(tvar(i)), fac(i), ' 0 ', tfil(i)
         else
            print*,'  Tracing field          : ',
     >                trim(tvar(i)),' : online calc not supported'
         endif
      enddo
      print*,'  Output (pmin,pmax,n)   : ',pmin,pmax,npre
      print*,'  Centering              : ',trim(centering)
      print*,'  Orientation            : ',trim(direction)
      print*,'  Coordinate Dump        : ',trim(dumpcoord)
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
      print*,'  First/last input file  : ',trim(dat(1)),
     >                                     ' ... ',
     >                                     trim(dat(numdat)) 
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

c     --------------------------------------------------------------------
c     Load the input trajectories
c     --------------------------------------------------------------------

c     Read the input trajectory file
      call ropen_tra(fid,inpfile,ntra,ntim,ncol,reftime,varsinp,inpmode)
      call read_tra (fid,trainp,ntra,ntim,ncol,inpmode)
      call close_tra(fid,inpmode)

c     Check that first four columns correspond to time,lon,lat,p
      if ( (varsinp(1).ne.'time' ).or.
     >     (varsinp(2).ne.'xpos' ).and.(varsinp(2).ne.'lon' ).or.
     >     (varsinp(3).ne.'ypos' ).and.(varsinp(3).ne.'lat' ).or.
     >     (varsinp(4).ne.'ppos' ).and.(varsinp(4).ne.'p'   ) )
     >then
         print*,' ERROR: problem with input trajectories ...'
         stop
      endif
      varsinp(1) = 'time'
      varsinp(2) = 'lon'
      varsinp(3) = 'lat'
      varsinp(4) = 'p'

c     Write some status information of the input trajectories
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

c     Check that first time is 0 - otherwise the tracing will produce
c     wrong results because later in the code only absolute times are
c     considered: <itime0   = int(abs(tfrac-tstart)/timeinc) + 1>. This 
c     will be changed in a future version.
      if ( abs( trainp(1,1,1) ).gt.eps ) then
         print*,' ERROR: First time of trajectory must be 0, i.e. '
         print*,'     correspond to the reference date. Otherwise'
         print*,'     the tracing will give wrong results... STOP'
         stop
      endif

c     If requested, open the coordinate dump file
      if ( dumpcoord.eq.'yes' ) then
         open(10,file=trim(outfile)//'.coord')
      endif

c     --------------------------------------------------------------------
c     Trace the fields (fields available on input files)
c     --------------------------------------------------------------------

      print*
      print*,'---- LIDAR FROM PRIMARY AND SECONDARY DATA FILES ------'

c     Loop over all tracing fields
      do i=1,ntrace0

c	      Skip all fields marked for online calculation
          if ( compfl(i).eq.1 ) goto 110
         
c	      Init the output fields: position and lidar field
	      do k=1,ntim
	      	do l=1,npre
	      		out_pos(k,l) = 0.
	      	    out_val(k,l) = 0.
	      	    out_cnt(k,l) = 0.
	      	 enddo
	      enddo

c         Write some status information
          print*
          print*,' Now lidaring           : ',
     >         trim(tvar(i)),compfl(i),' ',trim(tfil(i))

c         Reset flags for load manager
          iloaded0 = -1
          iloaded1 = -1

c         Reset the counter for fields outside domain
          noutside = 0

c         Loop over all times
          do j=1,ntim

c            Convert trajectory time from hh.mm to fractional time
             call hhmm2frac(trainp(1,j,1),tfrac)

c            Get the times which are needed
             itime0   = int(abs(tfrac-tstart)/timeinc) + 1
             time0    = tstart + fbflag * real(itime0-1) * timeinc
             itime1   = itime0 + 1
             time1    = time0 + fbflag * timeinc
             if ( itime1.gt.numdat ) then
                itime1 = itime0
                time1  = time0
             endif

c            Load manager: Check whether itime0 can be copied from itime1
             if ( itime0.eq.iloaded1 ) then
                f3t0     = f3t1
                p3t0     = p3t1
                spt0     = spt1
                iloaded0 = itime0
             endif

c            Load manager: Check whether itime1 can be copied from itime0
             if ( itime1.eq.iloaded0 ) then
                f3t1     = f3t0
                p3t1     = p3t0
                spt1     = spt0
                iloaded1 = itime1
             endif

c            Load manager:  Load first time (tracing variable and grid)
             if ( itime0.ne.iloaded0 ) then

                filename = tfil(i)//dat(itime0)
                call frac2hhmm(time0,tload)
                varname  = tvar(i) 
                write(*,'(a23,a20,a3,a5,f7.2)') 
     >               '    ->  loading          : ',
     >               trim(filename),'  ',trim(varname),tload
                call input_open (fid,filename)                
                call input_wind 
     >               (fid,varname,f3t0,tload,stagz,mdv,
     >               xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck) 
    
                call input_grid      
     >               (fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >                tload,pollon,pollat,p3t0,spt0,nz,ak,bk,stagz,
     >                timecheck)
                call input_close(fid)
                
                iloaded0 = itime0

             endif

c            Load manager: Load second time (tracing variable and grid)
             if ( itime1.ne.iloaded1 ) then
                
                filename = tfil(i)//dat(itime1)
                call frac2hhmm(time1,tload)
                varname  = tvar(i) 
                write(*,'(a23,a20,a3,a5,f7.2)') 
     >               '    ->  loading          : ',
     >               trim(filename),'  ',trim(varname),tload
                call input_open (fid,filename)
                call input_wind 
     >               (fid,varname,f3t1,tload,stagz,mdv,
     >               xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)      
                call input_grid      
     >               (fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >               tload,pollon,pollat,p3t1,spt1,nz,ak,bk,stagz,
     >               timecheck)
                call input_close(fid)
                
                iloaded1 = itime1
                
             endif

c            Loop over all trajectories
             do k=1,ntra
                
c               Set the trajectory position 
                x0_tra = trainp(k,j,2)                          ! Longitude
                y0_tra = trainp(k,j,3)                          ! Latitude
                p0_tra = trainp(k,j,4)                          ! Pressure

c               Get rotation angle - orient normal to trajectory
                if ( direction.eq.'normal' ) then

                   vx0 = 1.                     
                   vy0 = 0.

                   if ( j.lt.ntim ) then
                      lat =  0.5 * ( trainp(k,j,3) + trainp(k,j+1,3) )
                      vx1 = ( trainp(k,j+1,2) - trainp(k,j,2) ) * 
     >                      cos( lat * pi180 )
                      vy1 = ( trainp(k,j+1,3) - trainp(k,j,3) )
                   else
                      lat =  0.5 * ( trainp(k,j,3) + trainp(k,j-1,3) )
                      vx1 = ( trainp(k,j,2) - trainp(k,j-1,2) ) * 
     >                      cos( lat * pi180 )
                      vy1 = ( trainp(k,j,3) - trainp(k,j-1,3) )
                   endif

                   if ( vx1.gt.180  ) vx1 = vx1 - 360
                   if ( vx1.lt.-180 ) vx1 = vx1 + 360.

                   call getangle (vx0,vy0,vx1,vy1,rotation)
                   rotation = -rotation
                   
                else
                   rotation = 0.
                endif

c               Set the relative time
                call hhmm2frac(trainp(k,j,1),tfrac)
                reltpos0 = fbflag * (tfrac-time0)/timeinc    

c               Loop over pressure profile (or other positions for horizontal mode)
	        do l=1,npre

c                     Vertical
                      if ( direction.eq.'vertical' ) then
                        x0 = x0_tra
                        y0 = y0_tra
                        p0 = pmin + real(l-1)/real(npre-1) * (pmax-pmin)
                        if ( centering.eq.'yes' )then
                           p0 = p0 + trainp(k,j,4)
                        endif

c                     Longitude
                      elseif ( direction.eq.'lon' ) then
                        x0 = pmin + real(l-1)/real(npre-1) * (pmax-pmin)
                        y0 = y0_tra
                        p0 = p0_tra
                        if ( centering.eq.'yes' )then
                           x0 = x0 + x0_tra
                        endif
                         
c                     Latitude
                      elseif ( direction.eq.'lat' ) then
                        x0 = x0_tra
                        y0 = pmin + real(l-1)/real(npre-1) * (pmax-pmin)
                        p0 = p0_tra
                        if ( centering.eq.'yes' )then
                           y0 = y0 + y0_tra
                        endif

c                     Normal to trajerctory
                      elseif ( direction.eq.'normal' ) then

c                        Set the coordinate in the rotated system
                         rlat = pmin + 
     >                             real(l-1)/real(npre-1) * (pmax-pmin)
                         rlon = 0.

c                        Transform it back to geographical lon/lat
                         call getenvir_b (x0_tra,y0_tra,rotation,
     >                                              x0,y0,rlon,rlat,1)

c                        Pressure unchanged
                         p0 = p0_tra

                      endif

c                     Handle periodic boundaries in zonal direction
                      if ( (x0.gt.xmax).and.(per.ne.0) ) x0 = x0 - 360.
                      if ( (x0.lt.xmin).and.(per.ne.0) ) x0 = x0 + 360.

c                     Handle pole problems for hemispheric data (taken from caltra.f)
                      if ((hem.eq.1).and.(y0.gt.90.)) then
                         y0=180.-y0
                         x0=x0+per/2.
                      endif
                      if ((hem.eq.1).and.(y0.lt.-90.)) then
                         y0=-180.-y0
                         x0=x0+per/2.
                      endif
                      if (y0.gt.89.99) then
                         y0=89.99
                      endif     

c                 If requested, dump the lidar coordinates
                  if ( (dumpcoord.eq.'yes').and.(i.eq.1) ) then
                     write(10,'(3f10.2)') x0,y0,trainp(k,j,1)
                     write(10,'(3f10.2)') x0_tra,y0_tra,5.
                  endif

C                 Get the index where to interpolate (x0,y0,p0)
                  if ( (abs(x0-mdv).gt.eps).and.
     >                 (abs(y0-mdv).gt.eps) )
     >            then
                     call get_index4 (xind,yind,pind,x0,y0,p0,reltpos0,
     >                                p3t0,p3t1,spt0,spt1,3,
     >                                nx,ny,nz,xmin,ymin,dx,dy,mdv)
                  else
                     xind = mdv
                     yind = mdv
                     pind = mdv
                  endif

c                 If requested, apply nearest-neighbor interpolation
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

                  endif

c                 Do the interpolation: everthing is ok
                  if ( (xind.ge.1.).and.(xind.le.real(nx)).and.
     >                 (yind.ge.1.).and.(yind.le.real(ny)).and.
     >                 (pind.ge.1.).and.(pind.le.real(nz)) )
     >            then
                     f0 = int_index4(f3t0,f3t1,nx,ny,nz,
     >                               xind,yind,pind,reltpos0,mdv)

c                 Set to missing data
                  else
                     f0       = mdv
                  endif

c	              Save result to output array
                  if (abs(f0-mdv).gt.eps) then
                     out_val(j,l) = out_val(j,l) + f0 * fac(i)
                     out_cnt(j,l) = out_cnt(j,l) + 1.

	              endif

c              End loop over all pressure levels
	           enddo

c	           Save output - time index
	           ind_time = j

c                  Save output - space index for 'no centering'
	           if ( centering.eq.'no' ) then
                      if ( direction.eq.'vertical') then 
                         ind_pre  = nint( real(npre) *
     >          	       ( (p0_tra - pmin)/(pmax-pmin) ) + 1.)
                      elseif ( direction.eq.'lon') then 
                         ind_pre  = nint( real(npre) *
     >          	       ( (x0_tra - pmin)/(pmax-pmin) ) + 1.)
                      elseif ( direction.eq.'lat') then 
                         ind_pre  = nint( real(npre) *
     >          	       ( (y0_tra - pmin)/(pmax-pmin) ) + 1.)
                      endif

c                  Save output - space index for 'centering'
	           else
	           	  ind_pre  = nint( real(npre) *
     >          	       ( (0.            - pmin)/(pmax-pmin) ) + 1.)
	           endif

c                  Update the output array
	           if ( (ind_time.ge.1).and.(ind_time.le.ntim).and.
     >              (ind_pre .ge.1).and.(ind_pre .le.npre) )
     >         then
                    out_pos(ind_time,ind_pre) =
     >                          	out_pos(ind_time,ind_pre) + 1.
	           endif

c	         End loop over all trajectories
             enddo

c	      End loop over all times
          enddo

c	      Write the trajectory position to netCDF file - only once
	      if ( i.eq.1 ) then
	      	  cdfname  = outfile
	      	  varname  = 'POSITION'
	      	  longname = 'position of trajectory points'
	      	  unit     = 'none'
	      	  time     = 0.
              do k=1,npre
              	levels(k) = pmin + real(k-1)/real(npre-1) * (pmax-pmin)
              enddo
              do k=1,ntim
                 times(k) = trainp(1,k,1)
              enddo
              call writecdf2D_cf
     >            (cdfname,varname,longname,unit,out_pos,time,levels,
     >             times,npre,ntim,1,1,direction)
	      endif

c	      If no valid lidar count: set the field to missing data
          do k=1,ntim
          	do l=1,npre
          		if (abs(out_cnt(k,l)).lt.eps) then
          			out_val(k,l) = mdv
          	    endif
          	 enddo
          enddo

c	      If requested, calculate the mean of the lidar field
	      if ( outmode.eq.'mean' ) then
	      	do k=1,ntim
	      		do l=1,npre
	      			if ( (abs(out_val(k,l)-mdv).gt.eps).and.
     >                   (abs(out_cnt(k,l)    ).gt.0. ) )
     >              then
	      				out_val(k,l) = out_val(k,l) / out_cnt(k,l)
	      		    endif
	      		 enddo
	          enddo
	      endif

c	      Write the lidar field and count
	      cdfname  = outfile
	      if (outmode.eq.'sum' ) then
	         varname  = trim(tvar(i))//'_SUM'
	      elseif (outmode.eq.'mean' ) then
	         varname  = trim(tvar(i))//'_MEAN'
	      endif
	      longname = 'sum over all '//trim(tvar(i))//' profiles'
	      unit     = 'not given'
	      time     = 0.
          call writecdf2D_cf
     >            (cdfname,varname,longname,unit,out_val,time,levels,
     >             times,npre,ntim,0,1,direction)

	  	  cdfname  = outfile
	      varname  = trim(tvar(i))//'_CNT'
	      longname = 'counts of all '//trim(tvar(i))//' profiles'
	      unit     = 'not given'
	      time     = 0.
          call writecdf2D_cf
     >            (cdfname,varname,longname,unit,out_cnt,time,levels,
     >             times,npre,ntim,0,1,direction)

c         Exit point for loop over all tracing variables
 110      continue

c	   End loop over all lidar variables
       enddo


c     --------------------------------------------------------------------
c     Write output to netCDF file
c     --------------------------------------------------------------------

c     Write status information
      print*
      print*,'---- WRITE OUTPUT LIDAR FIELDS --------------------------'
      print*

c     Close coord dump file
      print*,' LIDAR written to      : ',trim(outfile)
      if ( dumpcoord.eq.'yes' ) then
        print*,' Coordinates dumped to : ',trim(outfile)//'.coord'
      endif

c     Write some status information, and end of program message
      print*  
      print*,'---- STATUS INFORMATION --------------------------------'
      print*
      print*,' ok'
      print*
      print*,'              *** END OF PROGRAM LIDAR ***'
      print*,'========================================================='


      end 


c     ********************************************************************
c     * INPUT / OUTPUT SUBROUTINES                                       *
c     ********************************************************************

c     --------------------------------------------------------------------
c     Subroutines to write 2D CF netcdf output file
c     --------------------------------------------------------------------

      subroutine writecdf2D_cf
     >          (cdfname,varname,longname,unit,arr,time,levels,times,
     >           npre,ntim,crefile,crevar,direction)

c     Create and write to the CF netcdf file <cdfname>. The variable
c     with name <varname> and with time <time> is written. The data
c     are in the two-dimensional array <arr>.  The flags <crefile> and
c     <crevar> determine whether the file and/or the variable should
c     be created.

      USE netcdf

      IMPLICIT NONE

c     Declaration of input parameters
      character*80 cdfname
      character*80 varname,longname,unit
      integer      npre,ntim
      real         arr(ntim,npre)
      real         levels(npre)
      real         times (ntim)
      real         time
      integer      crefile,crevar
      character*80 direction

c     Numerical epsilon
      real         eps
      parameter    (eps=1.e-5)

c     Local variables
      integer      ierr
      integer      ncID
      integer      LevDimId,    varLevID
      integer      TimeDimID,   varTimeID
      real         timeindex
      integer      i,j
      integer      nvars,varids(100)
      integer      ndims,dimids(100)
      real         timelist(1000)
      integer      ntimes
      integer      ind
      integer      varID

c     Quick an dirty solution for fieldname conflict
      if ( varname.eq.'time' ) varname = 'TIME'

c     Initially set error to indicate no errors.
      ierr = 0

c     ---- Create the netCDF - skip if <crefile=0> ----------------------
      if ( crefile.ne.1 ) goto 100

c     Create the file
      ierr = nf90_create(trim(cdfname), NF90_CLOBBER, ncID)

c     Define dimensions
      ierr=nf90_def_dim(ncID,'level',npre, LevDimID )
      ierr=nf90_def_dim(ncID,'time' ,ntim, TimeDimID)

c     Define space coordinate 
      ierr = nf90_def_var(ncID,'level',NF90_FLOAT,
     >     (/ LevDimID /),varLevID)
      if ( direction.eq.'vertical' ) then
        ierr = nf90_put_att(ncID, varLevID, "standard_name","level")
        ierr = nf90_put_att(ncID, varLevID, "units"        ,"hPa")
      elseif ( direction.eq.'lat' ) then
        ierr = nf90_put_att(ncID, varLevID, "standard_name","latitude")
        ierr = nf90_put_att(ncID, varLevID, "units"        ,"deg")
      elseif ( direction.eq.'lon' ) then
        ierr = nf90_put_att(ncID, varLevID, "standard_name","longitude")
        ierr = nf90_put_att(ncID, varLevID, "units"        ,"deg")
      elseif ( direction.eq.'normal' ) then
        ierr = nf90_put_att(ncID, varLevID, "standard_name","normal")
        ierr = nf90_put_att(ncID, varLevID, "units"        ,"deg")
      endif

c     Define time coordinate
      ierr = nf90_def_var(ncID,'time',NF90_FLOAT,
     >     (/ TimeDimID /), varTimeID)
      ierr = nf90_put_att(ncID, varTimeID, "long_name",    "time")
      ierr = nf90_put_att(ncID, varTimeID, "units",       "hours")

c     Write global attributes
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'Conventions', 'CF-1.0')
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'title',
     >     'pseudo-lidar from trajectory file')
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'source',
     >     'Lagranto Trajectories')
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'institution',
     >     'ETH Zurich, IACETH')

c     Check whether the definition was successful
      ierr = nf90_enddef(ncID)
      if (ierr.gt.0) then
         print*, 'An error occurred while attempting to ',
     >        'finish definition mode.'
         stop
      endif

c     Write coordinate data
      ierr = nf90_put_var(ncID,varLevID  ,levels)
      ierr = nf90_put_var(ncID,varTimeID ,times )

c     Close netCDF file
      ierr = nf90_close(ncID)

 100  continue

c     ---- Define a new variable - skip if <crevar=0> -----------------------

      if ( crevar.ne.1 ) goto 110

      print*,'Now defining ',trim(varname)

c     Open the file for read/write access
      ierr = nf90_open  (trim(cdfname), NF90_WRITE  , ncID)

c     Get the IDs for dimensions
      ierr = nf90_inq_dimid(ncID,'level', LevDimID )
      ierr = nf90_inq_dimid(ncID,'time' , TimeDimID)

c     Enter define mode
      ierr = nf90_redef(ncID)

c     Write definition and add attributes
      ierr = nf90_def_var(ncID,varname,NF90_FLOAT,
     >                   (/ TimeDimID, LevDimID /),varID)
      ierr = nf90_put_att(ncID, varID, "long_name" , longname )
      ierr = nf90_put_att(ncID, varID, "units"     , unit     )
      ierr = nf90_put_att(ncID, varID, '_FillValue', -999.99  )

c     Check whether definition was successful
      ierr = nf90_enddef(ncID)
      if (ierr.gt.0) then
         print*, 'An error occurred while attempting to ',
     >           'finish definition mode.'
         stop
      endif

c     Close netCDF file
      ierr = nf90_close(ncID)

 110  continue

c     ---- Write data --------------------------------------------------

c     Open the file for read/write access
      ierr = nf90_open  (trim(cdfname), NF90_WRITE , ncID)

c     Get the varID
      ierr = nf90_inq_varid(ncID,varname, varID )
      if (ierr.ne.0) then
         print*,'Variable ',trim(varname),' is not defined on ',
     >          trim(cdfname)
         stop
      endif

c     Write data block
      ierr = nf90_put_var(ncID,varID,arr,
     >                    start = (/ 1, 1 /),
     >                    count = (/ ntim, npre/) )

c     Check whether writing was successful
      ierr = nf90_close(ncID)
      if (ierr.ne.0) then
         write(*,*) trim(nf90_strerror(ierr))
         write(*,*) 'An error occurred while attempting to ',
     >              'close the netcdf file.'
         write(*,*) 'in clscdf_CF'
      endif

      end

c     ********************************************************************************
c     * Coordinate rotation - lidar normal to trajectory                             *
c     ********************************************************************************

c     --------------------------------------------------------------------------------
c     Backward coordinate transformation (Rotated lon/lat -> True lon/lat)
c     --------------------------------------------------------------------------------

      SUBROUTINE getenvir_b (clon,clat,rotation,
     >                       lon,lat,rlon,rlat,n)

      implicit none

c     Declaration of input and output parameters
      integer     n
      real        clon,clat,rotation
      real        lon(n), lat(n)
      real        rlon(n),rlat(n)

c     Auxiliary variables 
      real         pollon,pollat
      integer      i
      real         rlon1(n),rlat1(n)

c     Externals
      real         lmstolm,phstoph
      external     lmstolm,phstoph

c     First coordinate transformation (make the local coordinate system parallel to equator)
      pollon=-180.
      pollat=90.+rotation
      do i=1,n
         rlon1(i)=90.+lmstolm(rlat(i),rlon(i)-90.,pollat,pollon)
         rlat1(i)=phstoph(rlat(i),rlon(i)-90.,pollat,pollon)            
      enddo

c     Second coordinate transformation (make the local coordinate system parallel to equator)
      pollon=clon-180.
      if (pollon.lt.-180.) pollon=pollon+360.
      pollat=90.-clat
      do i=1,n
         lon(i)=lmstolm(rlat1(i),rlon1(i),pollat,pollon)
         lat(i)=phstoph(rlat1(i),rlon1(i),pollat,pollon)            
      enddo

      END

c     ---------------------------------------------------------------------
c     Determine the angle between two vectors
c     ---------------------------------------------------------------------

      SUBROUTINE getangle (vx1,vy1,vx2,vy2,angle)

c     Given two vectors <vx1,vy1> and <vx2,vy2>, determine the angle (in deg) 
c     between the two vectors.

      implicit none

c     Declaration of subroutine parameters
      real vx1,vy1
      real vx2,vy2
      real angle

c     Auxiliary variables and parameters
      real len1,len2,len3
      real val1,val2,val3
      real pi
      parameter (pi=3.14159265359)

      len1=sqrt(vx1*vx1+vy1*vy1)
      len2=sqrt(vx2*vx2+vy2*vy2)

      if ((len1.gt.0.).and.(len2.gt.0.)) then
         vx1=vx1/len1
         vy1=vy1/len1
         vx2=vx2/len2
         vy2=vy2/len2
         
         val1=vx1*vx2+vy1*vy2
         val2=-vy1*vx2+vx1*vy2
         
         len3=sqrt(val1*val1+val2*val2)
         
         if ( (val1.ge.0.).and.(val2.ge.0.) ) then
            val3=acos(val1/len3)
         else if ( (val1.lt.0.).and.(val2.ge.0.) ) then
            val3=pi-acos(abs(val1)/len3)
         else if ( (val1.ge.0.).and.(val2.le.0.) ) then
            val3=-acos(val1/len3)
         else if ( (val1.lt.0.).and.(val2.le.0.) ) then
            val3=-pi+acos(abs(val1)/len3)
         endif
      else
         val3=0.
      endif
      
      angle=180./pi*val3                                                                                     

      END

c     --------------------------------------------------------------------------------
c     Transformation routine: LMSTOLM and PHSTOPH from library gm2em            
c     --------------------------------------------------------------------------------

      REAL FUNCTION LMSTOLM (PHIS, LAMS, POLPHI, POLLAM)
C
C**** LMSTOLM  -   FC:BERECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE FUER
C****                 EINEN PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
C****                 IM ROTIERTEN SYSTEM. DER NORDPOL DES SYSTEMS HAT
C****                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   AUFRUF   :   LAM = LMSTOLM (PHIS, LAMS, POLPHI, POLLAM)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   BERECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE FUER
C**                EINEN PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
C**                IM ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
C**                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   VERSIONS-
C**   DATUM    :   03.05.90
C**
C**   EXTERNALS:   KEINE
C**   EINGABE-
C**   PARAMETER:   PHIS     REAL   GEOGR. BREITE DES PUNKTES IM ROT.SYS.
C**                LAMS     REAL   GEOGR. LAENGE DES PUNKTES IM ROT.SYS.
C**                POLPHI   REAL   WAHRE GEOGR. BREITE DES NORDPOLS
C**                POLLAM   REAL   WAHRE GEOGR. LAENGE DES NORDPOLS
C**   AUSGABE-
C**   PARAMETER:   WAHRE GEOGRAPHISCHE LAENGE ALS WERT DER FUNKTION
C**                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   D.MAJEWSKI
 
      REAL        LAMS,PHIS,POLPHI,POLLAM
 
      DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /
 
      ZSINPOL = SIN(ZPIR18*POLPHI)
      ZCOSPOL = COS(ZPIR18*POLPHI)
      ZLAMPOL = ZPIR18*POLLAM
      ZPHIS   = ZPIR18*PHIS
      ZLAMS   = LAMS
      IF(ZLAMS.GT.180.0) ZLAMS = ZLAMS - 360.0
      ZLAMS   = ZPIR18*ZLAMS
 
      ZARG1   = SIN(ZLAMPOL)*(- ZSINPOL*COS(ZLAMS)*COS(ZPHIS)  +
     1                          ZCOSPOL*           SIN(ZPHIS)) -
     2          COS(ZLAMPOL)*           SIN(ZLAMS)*COS(ZPHIS)
      ZARG2   = COS(ZLAMPOL)*(- ZSINPOL*COS(ZLAMS)*COS(ZPHIS)  +
     1                          ZCOSPOL*           SIN(ZPHIS)) +
     2          SIN(ZLAMPOL)*           SIN(ZLAMS)*COS(ZPHIS)
      IF (ABS(ZARG2).LT.1.E-30) THEN
        IF (ABS(ZARG1).LT.1.E-30) THEN
          LMSTOLM =   0.0
        ELSEIF (ZARG1.GT.0.) THEN
              LMSTOLAM =  90.0
            ELSE
              LMSTOLAM = -90.0
            ENDIF
      ELSE
        LMSTOLM = ZRPI18*ATAN2(ZARG1,ZARG2)
      ENDIF
 
      RETURN
      END


      REAL FUNCTION PHSTOPH (PHIS, LAMS, POLPHI, POLLAM)
C
C**** PHSTOPH  -   FC:BERECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE FUER
C****                 EINEN PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
C****                 ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
C****                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   AUFRUF   :   PHI = PHSTOPH (PHIS, LAMS, POLPHI, POLLAM)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   BERECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE FUER
C**                EINEN PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
C**                ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
C**                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   VERSIONS-
C**   DATUM    :   03.05.90
C**
C**   EXTERNALS:   KEINE
C**   EINGABE-
C**   PARAMETER:   PHIS     REAL   GEOGR. BREITE DES PUNKTES IM ROT.SYS.
C**                LAMS     REAL   GEOGR. LAENGE DES PUNKTES IM ROT.SYS.
C**                POLPHI   REAL   WAHRE GEOGR. BREITE DES NORDPOLS
C**                POLLAM   REAL   WAHRE GEOGR. LAENGE DES NORDPOLS
C**   AUSGABE-
C**   PARAMETER:   WAHRE GEOGRAPHISCHE BREITE ALS WERT DER FUNKTION
C**                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   D.MAJEWSKI
 
      REAL        LAMS,PHIS,POLPHI,POLLAM
 
      DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /
 
      SINPOL = SIN(ZPIR18*POLPHI)
      COSPOL = COS(ZPIR18*POLPHI)
      ZPHIS  = ZPIR18*PHIS
      ZLAMS  = LAMS
      IF(ZLAMS.GT.180.0) ZLAMS = ZLAMS - 360.0
      ZLAMS  = ZPIR18*ZLAMS
      ARG     = COSPOL*COS(ZPHIS)*COS(ZLAMS) + SINPOL*SIN(ZPHIS)
 
      PHSTOPH = ZRPI18*ASIN(ARG)
 
      RETURN
      END


      REAL FUNCTION LMTOLMS (PHI, LAM, POLPHI, POLLAM)
C
C%Z% Modul %M%, V%I% vom %G%, extrahiert am %H%
C
C**** LMTOLMS  -   FC:UMRECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE LAM
C****                 AUF EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
C****                 IM ROTIERTEN SYSTEM. DER NORDPOL DES SYSTEMS HAT
C****                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   AUFRUF   :   LAM = LMTOLMS (PHI, LAM, POLPHI, POLLAM)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   UMRECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE LAM AUF
C**                EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
C**                ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
C**                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   VERSIONS-
C**   DATUM    :   03.05.90
C**
C**   EXTERNALS:   KEINE
C**   EINGABE-
C**   PARAMETER:   PHI    REAL BREITE DES PUNKTES IM GEOGR. SYSTEM
C**                LAM    REAL LAENGE DES PUNKTES IM GEOGR. SYSTEM
C**                POLPHI REAL GEOGR.BREITE DES N-POLS DES ROT. SYSTEMS
C**                POLLAM REAL GEOGR.LAENGE DES N-POLS DES ROT. SYSTEMS
C**   AUSGABE-
C**   PARAMETER:   WAHRE GEOGRAPHISCHE LAENGE ALS WERT DER FUNKTION
C**                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   G. DE MORSIER
 
      REAL        LAM,PHI,POLPHI,POLLAM
 
      DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /
 
      ZSINPOL = SIN(ZPIR18*POLPHI)
      ZCOSPOL = COS(ZPIR18*POLPHI)
      ZLAMPOL =     ZPIR18*POLLAM
      ZPHI    =     ZPIR18*PHI
      ZLAM    = LAM
      IF(ZLAM.GT.180.0) ZLAM = ZLAM - 360.0
      ZLAM    = ZPIR18*ZLAM
 
      ZARG1   = - SIN(ZLAM-ZLAMPOL)*COS(ZPHI)
      ZARG2   = - ZSINPOL*COS(ZPHI)*COS(ZLAM-ZLAMPOL)+ZCOSPOL*SIN(ZPHI)
      IF (ABS(ZARG2).LT.1.E-30) THEN
        IF (ABS(ZARG1).LT.1.E-30) THEN
          LMTOLMS =   0.0
        ELSEIF (ZARG1.GT.0.) THEN
              LMTOLMS =  90.0
            ELSE
              LMTOLMS = -90.0
            ENDIF
      ELSE
        LMTOLMS = ZRPI18*ATAN2(ZARG1,ZARG2)
      ENDIF
 
      RETURN
      END


      REAL FUNCTION PHTOPHS (PHI, LAM, POLPHI, POLLAM)
C
C%Z% Modul %M%, V%I% vom %G%, extrahiert am %H%
C
C**** PHTOPHS  -   FC:UMRECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE PHI
C****                 AUF EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
C****                 IM ROTIERTEN SYSTEM. DER NORDPOL DES SYSTEMS HAT
C****                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   AUFRUF   :   PHI = PHTOPHS (PHI, LAM, POLPHI, POLLAM)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   UMRECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE PHI AUF
C**                EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
C**                ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
C**                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   VERSIONS-
C**   DATUM    :   03.05.90
C**
C**   EXTERNALS:   KEINE
C**   EINGABE-
C**   PARAMETER:   PHI    REAL BREITE DES PUNKTES IM GEOGR. SYSTEM
C**                LAM    REAL LAENGE DES PUNKTES IM GEOGR. SYSTEM
C**                POLPHI REAL GEOGR.BREITE DES N-POLS DES ROT. SYSTEMS
C**                POLLAM REAL GEOGR.LAENGE DES N-POLS DES ROT. SYSTEMS
C**   AUSGABE-
C**   PARAMETER:   ROTIERTE BREITE PHIS ALS WERT DER FUNKTION
C**                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   G. DE MORSIER
 
      REAL        LAM,PHI,POLPHI,POLLAM
 
      DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /
 
      ZSINPOL = SIN(ZPIR18*POLPHI)
      ZCOSPOL = COS(ZPIR18*POLPHI)
      ZLAMPOL = ZPIR18*POLLAM
      ZPHI    = ZPIR18*PHI
      ZLAM    = LAM
      IF(ZLAM.GT.180.0) ZLAM = ZLAM - 360.0
      ZLAM    = ZPIR18*ZLAM
      ZARG    = ZCOSPOL*COS(ZPHI)*COS(ZLAM-ZLAMPOL) + ZSINPOL*SIN(ZPHI)
 
      PHTOPHS = ZRPI18*ASIN(ZARG)
 
      RETURN
      END
      
