      PROGRAM create_startf

c     **************************************************************
c     * Create a <startfile> for <lagrangto>. It can be chosen     *
c     * whether to start from an isentropic or an isobaric         *
c     * surface. The starting points are equidistantly distributed *
c     * Michael Sprenger / Autumn 2004                             *
c     **************************************************************

      use netcdf

      implicit none

c     --------------------------------------------------------------
c     Set parameters
c     --------------------------------------------------------------

c     Maximum number of starting positions
      integer           nmax
      parameter         (nmax=4000000)

c     Maximum number of model levels
      integer           nlevmax
      parameter         (nlevmax=200)
      
c     Grid constant (distance in km corresponding to 1 deg at the equator)
      real              deltat
      parameter         (deltat=111.)

c     Mathematical constant (conversion degree -> radian)
      real              pi180
      parameter         (pi180=3.14159/180.)

c     Numerical epsilon
      real              eps
      parameter         (eps=0.00001)

c     --------------------------------------------------------------
c     Set variables
c     --------------------------------------------------------------

c     Filenames and output format
      character*80      pfile0,pfile1                  ! P filenames
      character*80      sfile0,sfile1                  ! S filenames
      character*80      ofile                          ! Output filename
      integer           oformat                        ! Output format
      real              timeshift                      ! Time shift relative to data files <*0>
      real              timeinc                        ! Time increment between input files

c     Horizontal grid
      character*80      hmode                          ! Horizontale mode
      real              lat1,lat2,lon1,lon2            ! Lat/lon boundaries
      real              ds,dlon,dlat                   ! Distance and lat/lon shifts
      character*80      hfile                          ! Filename
      integer           hn                             ! Number of entries in lat/lon list 
      real              latlist(nmax)                  ! List of latitudes
      real              lonlist(nmax)                  ! List of longitudes
      integer           pn                             ! Number of entries in lat/lon poly
      real              latpoly(500)                   ! List of polygon latitudes
      real              lonpoly(500)                   ! List of polygon longitudes
      real              loninpoly,latinpoly            ! Lon/lat inside polygon
      character*80      regionf                        ! Region file
      integer           iregion                        ! Region number
      real              xcorner(4),ycorner(4)          ! Vertices of region

c     Vertical grid
      character*80      vmode                          ! Vertical mode
      real              lev1,lev2,levlist(nmax)        ! Single levels, and list of levels
      character*80      vfile                          ! Filename
      integer           vn                             ! Number of entries

c     Unit of vertical axis
      character*80      umode                          ! Unit of vertical axis

c     Flag for 'no time check'
      character*80      timecheck                      ! Either 'no' or 'yes'

c     List of all starting positions
      integer           start_n                        ! Number of coordinates
      real              start_lat(nmax)                ! Latitudes
      real              start_lon(nmax)                ! Longitudes
      real              start_lev(nmax)                ! Levels (depending on vertical unit)
      real              start_pre(nmax)                ! Level in hPa
      integer           reftime(6)                     ! Reference time
      character*80      vars(10)                       ! Name of output fields (time,lon,lat,p)
      real,allocatable, dimension (:,:,:) :: tra       ! Trajectories (ntra,ntim,ncol)
      real              latmin,latmax
      real              lonmin,lonmax
      real              premin,premax

c     Grid description
      real              pollon,pollat                  ! Longitude/latitude of pole
      real              ak(nlevmax)                    ! Vertical layers and levels
      real              bk(nlevmax)
      real              xmin,xmax                      ! Zonal grid extension
      real              ymin,ymax                      ! Meridional grid extension
      integer           nx,ny,nz                       ! Grid dimensions
      real              dx,dy                          ! Horizontal grid resolution
      real,allocatable, dimension (:,:,:) :: pr        ! 3d pressure
      real,allocatable, dimension (:,:)   :: prs       ! surface pressure
      real,allocatable, dimension (:,:,:) :: th        ! 3d potential temperature
      real,allocatable, dimension (:,:)   :: ths       ! surface poential temperature
      real,allocatable, dimension (:,:,:) :: pv        ! 3d potential vorticity
      real,allocatable, dimension (:,:)   :: pvs       ! surface potential vorticiy
      real,allocatable, dimension (:,:,:) :: in        ! 3d 'dummy' array with vertical indices
      character*80      varname                        ! Name of input variable
      integer           fid                            ! File identifier
      real              stagz                          ! Vertical staggering
      real              mdv                            ! Missing data values
      real              tstart,tend                    ! Time on P and S file
      real              rid,rjd,rkd                    ! Real grid position

c     Auxiliary variable
      integer           i,j,k
      real              lon,lat
      real              rd
      integer           stat,flag
      real              tmp1,tmp2
      real              tfrac,frac
      real              radius,dist
      character*80      string
      character*80      selectstr
      character*80      umode_save
      character*80      maskname
      real              maskvalue
      character*2       maskoper
      integer           cdfid,varid,ierr
      integer           nmask
      integer           indx,indy
      character         ch1,ch2
      integer           len
      real,allocatable, dimension (:,:,:) :: fld0
      real,allocatable, dimension (:,:,:) :: fld1
      real,allocatable, dimension (:,:  ) :: sfc0
      real,allocatable, dimension (:,:)   :: sfc1
      real,allocatable, dimension (:,:)   :: mask

c     Externals 
      real              int_index3     ! 3d interpolation
      external          int_index3
      real              sdis           ! Speherical distance
      external          sdis
      integer           inregion       ! In/out of region
      external          inrehion

c     ------------------------------------------------------------------
c     Start of program, Read parameters
c     ------------------------------------------------------------------

c     Write start message
      print*,'========================================================='
      print*,'         *** START OF PROGRAM CREATE_STARTF ***'
      print*

c     Read parameter file
      open(10,file='create_startf.param')

c      Input P and S file
       read(10,*) pfile0,pfile1
       read(10,*) sfile0,sfile1
       read(10,*) ofile

c      Read name of region file
       read(10,*) regionf

c      Reference time
       do i=1,6
          read(10,*) reftime(i)
       enddo

c      Time shift relative to data files <pfile0,sfile0> - format (hh.mm)
       read(10,*) timeshift

c      Read timeincrement between input files
       read(10,*) timeinc

c      Parameters for horizontal grid
       read(10,*) hmode
       if ( hmode.eq.'file' ) then                  ! from file
          read(10,*) hfile
       elseif ( hmode.eq.'line' ) then              ! along a line
          read(10,*) lon1,lon2,lat1,lat2,hn
       elseif ( hmode.eq.'box.eqd' ) then           ! box: 2d equidistant
          read(10,*) lon1,lon2,lat1,lat2,ds
       elseif ( hmode.eq.'box.grid' ) then          ! box: 2d grid
          read(10,*) lon1,lon2,lat1,lat2
       elseif ( hmode.eq.'point' ) then             ! single point
          read(10,*) lon1,lat1
       elseif ( hmode.eq.'shift' ) then             ! centre + shifted
          read(10,*) lon1,lat1,dlon,dlat
       elseif ( hmode.eq.'polygon.eqd' ) then       ! polygon: 2d equidistant
          read(10,*) hfile,ds
       elseif ( hmode.eq.'polygon.grid' ) then      ! polygon: 2d grid
          read(10,*) hfile
       elseif ( hmode.eq.'circle.eqd' ) then        ! circle: 2d equidistant
          read(10,*) lon1,lat1,radius,ds 
       elseif ( hmode.eq.'circle.grid' ) then       ! circle: 2d grid
          read(10,*) lon1,lat1,radius
       elseif ( hmode.eq.'region.eqd' ) then        ! region: 2d equidistant
          read(10,*) iregion,ds 
       elseif ( hmode.eq.'region.grid' ) then       ! iregion: 2d grid
          read(10,*) iregion
       elseif ( hmode.eq.'mask.grid' ) then         ! 0/1 masked - grid points
          read(10,*) hfile,maskname
       elseif ( hmode.eq.'mask.eqd' ) then         ! 0/1 masked - euqidistant points
          read(10,*) hfile,maskname,ds
       else
          print*,' ERROR: horizontal mode not supported ',trim(hmode)
          stop
       endif

c      Parameters for vertical grid
       read(10,*) vmode
       if ( vmode.eq.'file') then                   ! from file
          read(10,*) vfile
       elseif ( vmode.eq.'level' ) then             ! single level (explicit command)
          read(10,*) lev1                         
       elseif ( vmode.eq.'list') then               ! a list
          read(10,*) vn
          read(10,*) (levlist(i),i=1,vn)
       elseif ( vmode.eq.'profile') then            ! a profile
          read(10,*) lev1,lev2,vn
       elseif ( vmode.eq.'grid') then               ! grid points
          read(10,*) lev1,lev2
       else
          print*,' ERROR: vertical mode not supported ',trim(vmode)
          stop
       endif

c      Read units of vertical axis
       read(10,*) umode
       if ( ( umode.ne.'hPa'     ).and.
     >      ( umode.ne.'hPa,agl' ).and.   
     >      ( umode.ne.'K'       ).and.
     >      ( umode.ne.'PVU'     ).and.
     >      ( umode.ne.'INDEX'   )  ) 
     > then  
          print*,' ERROR: unit not supported ',trim(umode)
          stop
       endif

c     Read selection criterion (dummy read)
      read(10,*) selectstr

c     Read flag for 'no time check'
      read(10,*) timecheck 

c     Close parameter file
      close(10)

c     Decide which output format is used (1..4: trajectory format, -1: triple list)
      call mode_tra(oformat,ofile)
      
c     Decide whether all lat/lon/lev coordaintes are read from one file
      if ( (hmode.eq.'file').and.(vmode.eq.'nil') ) then
         hmode='file3'
      elseif ( (hmode.eq.'file').and.(vmode.ne.'nil') ) then
         hmode='file2'
      endif

c     Convert timeshift (hh.mm) into a fractional time shift
      call hhmm2frac(timeshift,tfrac)
      if (tfrac.gt.0.) then
         tfrac=tfrac/timeinc
      else
         tfrac=0.
      endif

c     If a mask file is provided, no time interpolation is allowed
      if ( (hmode.eq.'mask.grid').or.(hmode.eq.'mask.eqd') ) then
         if ( abs(tfrac).gt.eps ) then
            print*,' ERROR: no intermediate times allowed for ',
     >             trim(hmode)
            stop
         endif
      endif

c     Read the region coordinates if needed
      if ( (hmode.eq.'region.eqd' ).or.
     >     (hmode.eq.'region.grid') ) then
         
         open(10,file=regionf)
          
 50       read(10,*,end=51) string
          
          if ( string(1:1).ne.'#' ) then
             call regionsplit(string,i,xcorner,ycorner)
             if ( i.eq.iregion ) goto 52
          endif

          goto 50

 51      close(10)
         
         print*,' ERROR: region ',iregion,' not found on ',trim(regionf)
         stop
 
 52      continue
          
      endif

c     Write some status information
      print*,'---- INPUT PARAMETERS -----------------------------------'
      print*
      if ( timeshift.gt.0. ) then
         print*,'  P file                     : ',trim(pfile0),
     >                                            ' ',
     >                                            trim(pfile1)
         print*,'  S file                     : ',trim(sfile0),
     >                                            ' ',
     >                                            trim(sfile1)
      else
         print*,'  P file                     : ',trim(pfile0)
         print*,'  S file                     : ',trim(sfile0)
      endif
      print*,'  Output file                : ',trim(ofile) 
      print*
      if (oformat.eq.-1) then
         print*,'  Output format              : (lon,lat,lev)-list'
      else
         print*,'  Output format              : ',oformat
      endif
      print*
      print*,'  Reference time (year)      : ',reftime(1)
      print*,'                 (month)     : ',reftime(2)
      print*,'                 (day)       : ',reftime(3)
      print*,'                 (hour)      : ',reftime(4)
      print*,'                 (min)       : ',reftime(5)
      print*,'  Time range                 : ',reftime(6)
      print*
      print*,'  Time shift                 : ',timeshift,' + ',
     >                                         trim(pfile0)
      print*,'  Region file                : ',trim(regionf)
      print*
      print*,'  hmode                      : ',trim(hmode)
      if ( hmode.eq.'file2' ) then        
        print*,'      filename [lat/lon]      : ',trim(hfile)
      elseif ( hmode.eq.'file3' ) then        
        print*,'      filename [lat/lon/lev]  : ',trim(hfile)
      elseif ( hmode.eq.'line' ) then   
        write(*,'(a30,4f10.2,i4)') 
     >         '      lon1,lon2,lat1,lat2,n   : ',lon1,lon2,lat1,lat2,hn
      elseif ( hmode.eq.'box.eqd' ) then     
        write(*,'(a30,5f10.2)') 
     >         '      lon1,lon2,lat1,lat2,ds  : ',lon1,lon2,lat1,lat2,ds
      elseif ( hmode.eq.'box.grid' ) then    
        write(*,'(a30,4f10.2)') 
     >         '      lon1,lon2,lat1,lat2     : ',lon1,lon2,lat1,lat2 
      elseif ( hmode.eq.'point' ) then   
        print*,'      lon,lat                 : ',lon1,lat1
      elseif ( hmode.eq.'shift' ) then   
        write(*,'(a30,4f10.2)') 
     >         '      lon,lat,dlon,dlat       : ',lon1,lat1,dlon,dlat
      elseif ( hmode.eq.'polygon.eqd' ) then   
        write(*,'(a30,a10,f10.2)') 
     >         '      hfile, ds               : ',trim(hfile),ds
      elseif ( hmode.eq.'polygon.grid' ) then   
        write(*,'(a30,a10)') 
     >         '      hfile                   : ',trim(hfile)
      elseif ( hmode.eq.'circle.eqd' ) then   
        write(*,'(a30,4f10.2)') 
     >         '      lonc,latc,radius, ds    : ',lon1,lat1,radius,ds
      elseif ( hmode.eq.'circle.grid' ) then   
        write(*,'(a30,3f10.2)') 
     >         '      lonc,latc,radius        : ',lon1,lat1,radius
      elseif ( hmode.eq.'region.eqd' ) then   
        write(*,'(a30,i4,1f10.2)') 
     >         '      iregion, ds             : ',iregion,ds
        write(*,'(a30,4f10.2)')
     >         '      xcorner                 : ',(xcorner(i),i=1,4)
        write(*,'(a30,4f10.2)')
     >         '      ycorner                 : ',(ycorner(i),i=1,4)
      elseif ( hmode.eq.'region.grid' ) then   
        write(*,'(a30,i4)') 
     >         '      iregion                 : ',iregion
        write(*,'(a30,4f10.2)')
     >         '      xcorner                 : ',(xcorner(i),i=1,4)
        write(*,'(a30,4f10.2)')
     >         '      ycorner                 : ',(ycorner(i),i=1,4)
      elseif ( hmode.eq.'mask.eqd' ) then
        write(*,'(a30,a15,a15,f10.2)')
     >         '      hfile,variable,ds       : ',trim(hfile),
     >                                            trim(maskname),
     >                                            ds
      elseif ( hmode.eq.'mask.grid' ) then
        write(*,'(a30,a15,a15)')
     >         '      hfile,variable          : ',trim(hfile),
     >                                            trim(maskname)
      endif

      print*
      print*,'  vmode                      : ',trim(vmode)
      if ( vmode.eq.'file') then 
         print*,'      filename               : ',trim(vfile)
      elseif ( vmode.eq.'level' ) then 
         print*,'      level                  : ',lev1                         
      elseif ( vmode.eq.'list') then 
         print*,'      n                      : ',vn
         print*,'      level(i)               : ',(levlist(i),i=1,vn)
      elseif ( vmode.eq.'profile') then 
         print*,'      lev1,lev2,n            : ',lev1,lev2,vn
      elseif ( vmode.eq.'grid') then 
         print*,'      lev1,lev2              : ',lev1,lev2
      endif

      print* 
      print*,'  umode                      : ',trim(umode)
      print*
      print*,'  time check                 : ',trim(timecheck)
      print*

c     <mask.grid> and <mask.eqd>: split variable into name, operator,value
      if ( (hmode.eq.'mask.grid').or.(hmode.eq.'mask.eqd') ) then
          len = len_trim(maskname)
          do i=1,len_trim(maskname)-1
            ch1 = maskname(i:i)
            ch2 = maskname(i+1:i+1)

            if ( (ch1.eq.'<').and.(ch2.eq.'>') ) then
               read(maskname(i+2:len),*) maskvalue
               maskname = maskname(1:i-1)
               maskoper = 'ne'
               goto 90
            elseif ( (ch1.eq.'<').and.(ch2.eq.'=') ) then
               read(maskname(i+2:len),*) maskvalue
               maskname = maskname(1:i-1)
               maskoper = 'le'
               goto 90
            elseif ( (ch1.eq.'>').and.(ch2.eq.'=') ) then
               read(maskname(i+2:len),*) maskvalue
               maskname = maskname(1:i-1)
               maskoper = 'ge'
               goto 90
            elseif ( (ch1.eq.'=').and.(ch2.eq.'=') ) then
               read(maskname(i+2:len),*) maskvalue
               maskname = maskname(1:i-1)
               maskoper = 'eq'
               goto 90
            elseif ( ch1.eq.'=' ) then
               read(maskname(i+1:len),*) maskvalue
               maskname = maskname(1:i-1)
               maskoper = 'eq'
               goto 90
            elseif ( ch1.eq.'>' ) then
               read(maskname(i+1:len),*) maskvalue
               maskname = maskname(1:i-1)
               maskoper = 'gt'
               goto 90
            elseif ( ch1.eq.'<' ) then
               read(maskname(i+1:len),*) maskvalue
               maskname = maskname(1:i-1)
               maskoper = 'lt'
               goto 90
             endif

          enddo
      endif

  90  continue

c     ------------------------------------------------------------------
c     Read grid parameters from inital files
c     ------------------------------------------------------------------

c     Get the time of the first and second data file
      tstart = -timeshift                              ! Format hh.mm
      call hhmm2frac(tstart,frac)
      frac   = frac + timeinc          
      call frac2hhmm(frac,tend)                        ! Format hh.mm                  

c     Convert timeshift (hh.mm) into a fractional time shift
      tfrac=real(int(timeshift))+
     >          100.*(timeshift-real(int(timeshift)))/60.
      if (tfrac.gt.0.) then
         tfrac=tfrac/timeinc
      else
         tfrac=0.
      endif

c     Read the constant grid parameters (nx,ny,nz,xmin,xmax,ymin,ymax,
c     pollon,pollat) The negative <-fid> of the file identifier is used 
c     as a flag for parameter retrieval
      varname  = 'U'
      nx       = 1
      ny       = 1
      nz       = 1
      call input_open (fid,pfile0)
      call input_grid (-fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >                 tstart,pollon,pollat,rd,rd,nz,rd,rd,rd,timecheck)
      call input_close(fid)

c     Check whether region coordinates are within the domain
      if ( (hmode.eq.'region.eqd' ).or.
     >     (hmode.eq.'region.grid') ) then

         do i=1,4
            if ( (xcorner(i).lt.xmin).or.
     >           (ycorner(i).lt.ymin).or.
     >           (xcorner(i).gt.xmax).or.
     >           (ycorner(i).gt.ymax) )
     >      then
               print*,' ERROR: region not included in data domain...'
               print*,'     ',trim(string)
               print*,'     ',(xcorner(j),j=1,4)
               print*,'     ',(ycorner(j),j=1,4)
               stop
            endif
               
         enddo

      endif

C     Check if the number of levels is too large
      if (nz.gt.nlevmax) goto 993

c     Allocate memory for 3d arrays: pressure, theta, pv
      allocate(pr(nx,ny,nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array pr ***'
      allocate(prs(nx,ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array prs **'
      allocate(th(nx,ny,nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array th ***'
      allocate(ths(nx,ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array ths **'
      allocate(pv(nx,ny,nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array pv ***'
      allocate(pvs(nx,ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array pvs **'
      allocate(in(nx,ny,nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array in ***'
      allocate(mask(nx,ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array mask ***'

c     Allocate memory for temporary arrays for time interpolation
      allocate(fld0(nx,ny,nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tmp0 ***'
      allocate(fld1(nx,ny,nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tmp1 ***'
      allocate(sfc0(nx,ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array sfc0 ***'
      allocate(sfc1(nx,ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array sfc1 ***'

c     ------ Index -----------------------------------------------------

c     Init the dummy array with vertical index
      do i=1,nx
         do j=1,ny
            do k=1,nz
               in(i,j,k) = real(k)
            enddo
         enddo
      enddo

c     ------ Pressure --------------------------------------------------

c     Read pressure from first data file (pfile0) on U-grid; we have to set
c     mdv explicitely, because it's not read from netCDF
      call input_open (fid,pfile0)
      varname='U'
      stagz  = -0.5
      call input_grid                                      
     >     (fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >      tstart,pollon,pollat,fld0,sfc0,nz,ak,bk,stagz,timecheck)
      mdv = -999.99 
      call input_close(fid)

c     Read or set pressure for second data file (pfile1)
      if ( timeshift.ne.0.) then
         call input_open (fid,pfile1)
         varname='U'
         call input_grid                                      
     >        (fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >        tend,pollon,pollat,fld1,sfc1,nz,ak,bk,stagz,timecheck)
         call input_close(fid)
      else
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  fld1(i,j,k) = fld0(i,j,k)
               enddo
               sfc1(i,j) = sfc0(i,j)
            enddo
         enddo
      endif

c     Time interpolation to get the final pressure field
      do i=1,nx
         do j=1,ny
            do k=1,nz
               pr(i,j,k) = (1.-tfrac) * fld0(i,j,k) +
     >                     tfrac      * fld1(i,j,k)
            enddo
            prs(i,j) = (1.-tfrac) * sfc0(i,j) +
     >                 tfrac      * sfc1(i,j)
         enddo
      enddo

c     ------ Potential temperature -------------------------------------

       if ( (umode.eq.'K').or.(umode.eq.'PVU') ) then

c         Read potential temperature from first data file <sfile0>
          call input_open (fid,sfile0)
          varname='TH'                                      ! Theta                                      
          call input_wind
     >         (fid,varname,fld0,tstart,stagz,mdv,
     >         xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)
          call input_close(fid)

c         Read or set potential temperature for second data file (sfile1)
          if ( timeshift.ne.0.) then
             call input_open (fid,sfile1)
             varname='TH' 
             call input_wind
     >            (fid,varname,fld1,tend,stagz,mdv,
     >            xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)
             call input_close(fid)
          else
             do i=1,nx
                do j=1,ny
                   do k=1,nz
                      fld1(i,j,k) = fld0(i,j,k)
                   enddo
                enddo
             enddo
          endif

c         Time interpolation to get the final potential temperature field
          do i=1,nx
             do j=1,ny
                do k=1,nz
                   th(i,j,k) = (1.-tfrac) * fld0(i,j,k) +
     >                         tfrac      * fld1(i,j,k)
                enddo
             enddo
          enddo
          
c         Set the surface potential temperature
          do i=1,nx                            
             do j=1,ny
                ths(i,j)=th(i,j,1)
             enddo
          enddo
          
       endif

c     ------ Potential vorticity -----------------------------------------

       if ( (umode.eq.'PVU') ) then

c         Read potential vorticity from first data file <sfile0>
          call input_open (fid,sfile0)
          varname='PV'            
          call input_wind
     >         (fid,varname,fld0,tstart,stagz,mdv,
     >         xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)
          call input_close(fid)
          
c         Read or set potential vorticity for second data file (sfile1)
          if ( timeshift.ne.0.) then
             call input_open (fid,sfile1)
             varname='PV' 
             call input_wind
     >            (fid,varname,fld1,tend,stagz,mdv,
     >            xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)
             call input_close(fid)
          else
             do i=1,nx
                do j=1,ny
                   do k=1,nz
                      fld1(i,j,k) = fld0(i,j,k)
                   enddo
                enddo
             enddo
          endif

c         Time interpolation to get the final potential vorticity field
          do i=1,nx
             do j=1,ny
                do k=1,nz
                   pv(i,j,k) = (1.-tfrac) * fld0(i,j,k) +
     >                         tfrac      * fld1(i,j,k)
                enddo
             enddo
          enddo

c         Set the surface potential vorticity
          do i=1,nx                          
             do j=1,ny
                pvs(i,j)=pv(i,j,1)
             enddo
          enddo
       endif

c     ------ Load 0/1 label (mask) ------------------------------------
      if ( (hmode.eq.'mask.grid').or.(hmode.eq.'mask.eqd') ) then

c         Explicit netCDF calls - no consistency check for grids!
          ierr = NF90_OPEN(hfile,nf90_nowrite, cdfid)
          IF ( ierr /= nf90_NoErr ) PRINT *,NF90_STRERROR(ierr)
          ierr = NF90_INQ_VARID(cdfid,maskname,varid)
          IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
          ierr = nf90_get_var(cdfid,varid,mask)
          IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
          ierr = NF90_CLOSE(cdfid)
          IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)

          nmask = 0
          do i=1,nx
            do j=1,ny

                if ( maskoper.eq.'eq' ) then
                   if ( abs(mask(i,j)-maskvalue).lt.eps ) then
                      nmask = nmask + 1
                      mask(i,j) = 1.
                   else
                      mask(i,j) = 0.
                   endif
                endif

                if ( maskoper.eq.'ne' ) then
                   if ( abs(mask(i,j)-maskvalue).gt.eps ) then
                      nmask = nmask + 1
                      mask(i,j) = 1.
                   else
                      mask(i,j) = 0.
                   endif
                endif

                if ( maskoper.eq.'gt' ) then
                   if ( mask(i,j).gt.maskvalue ) then
                      nmask = nmask + 1
                      mask(i,j) = 1.
                   else
                      mask(i,j) = 0.
                   endif
                endif

                if ( maskoper.eq.'lt' ) then
                   if ( mask(i,j).lt.maskvalue ) then
                      nmask = nmask + 1
                      mask(i,j) = 1.
                   else
                      mask(i,j) = 0.
                   endif
                endif

                if ( maskoper.eq.'ge' ) then
                   if ( mask(i,j).ge.maskvalue ) then
                      nmask = nmask + 1
                      mask(i,j) = 1.
                   else
                      mask(i,j) = 0.
                   endif
                endif

                if ( maskoper.eq.'le' ) then
                   if ( mask(i,j).le.maskvalue ) then
                      nmask = nmask + 1
                      mask(i,j) = 1.
                   else
                      mask(i,j) = 0.
                   endif
                endif

             enddo
          enddo

          print*
          print*,'  Mask read from file -> ',nmask,' 1-points'
          print*

      endif

c     Write some status information
      print*,'---- CONSTANT GRID PARAMETERS ---------------------------'
      print*
      print*,'  xmin,xmax        : ',xmin,xmax
      print*,'  ymin,ymax        : ',ymin,ymax
      print*,'  dx,dy            : ',dx,dy
      print*,'  pollon,pollat    : ',pollon,pollat
      print*,'  nx,ny,nz         : ',nx,ny,nz
      print*
      print*,'  Pressure loaded  : ',trim(pfile0),' ',trim(pfile1)
      if ( (umode.eq.'K').or.(umode.eq.'PVU') ) then
         print*,'  Theta loaded     : ',trim(sfile0),' ',trim(sfile1)
      endif
      if ( (umode.eq.'PVU') ) then
         print*,'  PV loaded        : ',trim(sfile0),' ',trim(sfile1)
      endif
      print*

c     ------------------------------------------------------------------
c     Determine the expanded list of starting coordinates
c     ------------------------------------------------------------------

c     Write some status information
      print*,'---- EXPAND LIST OF STARTING POSITIONS  -----------------'
      print*

c     ------ Read lat/lon/lev from <hfile> -----------------------------
      if ( hmode.eq.'file3' ) then
         start_n = 0
         open(10,file=hfile)
 100       continue
              start_n = start_n + 1
              read(10,*,end=101) start_lon(start_n),
     >                           start_lat(start_n),
     >                           start_lev(start_n)
              goto 100
 101       continue
           start_n = start_n - 1
         close(10)
         goto 400
      endif

c     ------ Get lat/lon (horizontal) coordinates ---------------------

c     Read lat/lon from <hfile> 
      if ( hmode.eq.'file2' ) then
         hn = 0
         open(10,file=hfile)
 200       continue
              hn = hn + 1
              read(10,*,end=201) lonlist(hn),
     >                           latlist(hn)
              goto 200
 201       continue
           hn = hn - 1
         close(10)
      endif

c     Get lat/lon along a line (linear in lat/lon space)
      if ( hmode.eq.'line' ) then
         do i=1,hn
            lonlist(i) = lon1 + real(i-1)/real(hn-1)*(lon2-lon1)
            latlist(i) = lat1 + real(i-1)/real(hn-1)*(lat2-lat1)
         enddo
      endif

c     Lat/lon box: equidistant
      if ( hmode.eq.'box.eqd' ) then
         hn  = 0
         lat = lat1
         do while ( lat.le.lat2 )      
           lon = lon1
           do while ( lon.le.lon2 ) 
             hn          = hn+1
             lonlist(hn) = lon
             latlist(hn) = lat
             lon         = lon+ds/(deltat*cos(pi180*lat))
           enddo
           lat = lat+ds/deltat
        enddo
      endif

c     Lat/lon box: grid
      if ( hmode.eq.'box.grid' ) then
         hn = 0
         do j=1,ny
            do i=1,nx
               lon = xmin + real(i-1) * dx
               lat = ymin + real(j-1) * dy
               if ( (lon.ge.lon1).and.(lon.le.lon2).and.
     >              (lat.ge.lat1).and.(lat.le.lat2) )
     >         then
                  hn          = hn+1
                  lonlist(hn) = lon
                  latlist(hn) = lat
               endif
            enddo
         enddo
      endif

c     Get single starting point
      if ( hmode.eq.'point' ) then
         hn          = 1
         lonlist(hn) = lon1
         latlist(hn) = lat1
      endif

c     Get shifted and central starting point
      if ( hmode.eq.'shift' ) then
         hn         = 5
         lonlist(1) = lon1
         latlist(1) = lat1
         lonlist(2) = lon1+dlon
         latlist(2) = lat1
         lonlist(3) = lon1-dlon
         latlist(3) = lat1
         lonlist(4) = lon1
         latlist(4) = lat1+dlat
         lonlist(5) = lon1
         latlist(5) = lat1-dlat
      endif

c     Lat/lon polygon: grid
      if ( hmode.eq.'polygon.grid' ) then

c        Read list of polygon coordinates
         pn = 0
         open(10,file=hfile)
           read(10,*) loninpoly,latinpoly
 210       continue
              pn = pn + 1
              read(10,*,end=211) lonpoly(pn),
     >                           latpoly(pn)

              print*,pn,lonpoly(pn),latpoly(pn)
              
              goto 210
 211       continue
           pn = pn - 1
         close(10)

c        Define the polygon boundaries
         call DefSPolyBndry(latpoly,lonpoly,pn,latinpoly,loninpoly)

c        Get the grid points inside the polygon
         hn = 0
         do j=1,ny
            do i=1,nx
               lon = xmin + real(i-1) * dx
               lat = ymin + real(j-1) * dy

               call LctPtRelBndry(lat,lon,flag)

               if ( (flag.eq.1).or.(flag.eq.2) ) then
                  hn          = hn+1
                  lonlist(hn) = lon
                  latlist(hn) = lat
               endif

            enddo
         enddo
         
      endif

c     Lat/lon polygon: equidistant
      if ( hmode.eq.'polygon.eqd' ) then

c        Read list of polygon coordinates
         pn = 0

         open(10,file=hfile)
           read(10,*) loninpoly,latinpoly
 220       continue
              pn = pn + 1
              read(10,*,end=221) lonpoly(pn),
     >                           latpoly(pn)
              goto 220
 221       continue
           pn = pn - 1
         close(10)


c        Define the polygon boundaries
         call DefSPolyBndry(latpoly,lonpoly,pn,latinpoly,loninpoly)

c        Get the grid points inside the polygon
         hn  = 0
         lat = -90.
         do while ( lat.le.90. )      
           lon = -180.
           do while ( lon.lt.180. )

              call LctPtRelBndry(lat,lon,flag)
               
               if ( (flag.eq.1).or.(flag.eq.2) ) then
                  hn          = hn+1
                  lonlist(hn) = lon
                  latlist(hn) = lat

               endif
               
               lon = lon+ds/(deltat*cos(pi180*lat))
           enddo
           lat = lat+ds/deltat

        enddo

      endif

c     Circle: equidistant
      if ( hmode.eq.'circle.eqd' ) then
         hn  = 0
         lat = ymin
         do while ( lat.le.ymax )      
           lon = xmin
           do while ( lon.le.xmax ) 
              dist = sdis(lon1,lat1,lon,lat)
              if ( dist.le.radius ) then
                 hn          = hn+1
                 lonlist(hn) = lon
                 latlist(hn) = lat
              endif
              lon = lon+ds/(deltat*cos(pi180*lat))
           enddo
           lat = lat+ds/deltat
        enddo
      endif

c     Circle: grid
      if ( hmode.eq.'circle.grid' ) then
         hn = 0
         do j=1,ny
            do i=1,nx
               lon = xmin + real(i-1) * dx
               lat = ymin + real(j-1) * dy
               dist = sdis(lon1,lat1,lon,lat)
               if ( dist.le.radius ) then
                  hn          = hn+1
                  lonlist(hn) = lon
                  latlist(hn) = lat
               endif
            enddo
         enddo
         
      endif

c     Region: equidistant
      if ( hmode.eq.'region.eqd' ) then
         hn  = 0
         lat = ymin
         do while ( lat.le.ymax )      
           lon = xmin
           do while ( lon.le.xmax ) 
              flag = inregion(lon,lat,xcorner,ycorner)
              if ( flag.eq.1 ) then
                 hn          = hn+1
                 lonlist(hn) = lon
                 latlist(hn) = lat
              endif
              lon = lon+ds/(deltat*cos(pi180*lat))
           enddo
           lat = lat+ds/deltat
        enddo
      endif

c     Region: grid
      if ( hmode.eq.'region.grid' ) then
         hn = 0
         do j=1,ny
            do i=1,nx
               lon  = xmin + real(i-1) * dx
               lat  = ymin + real(j-1) * dy
               flag = inregion(lon,lat,xcorner,ycorner)
               if ( flag.eq.1 ) then
                  hn          = hn+1
                  lonlist(hn) = lon
                  latlist(hn) = lat
               endif
            enddo
         enddo
      endif

c     mask: grid
      if ( hmode.eq.'mask.grid' ) then
         hn  = 0
         do i=1,nx
          do j=1,ny
            if ( mask(i,j).gt.0.5 ) then
               hn = hn +1
               lonlist(hn) = xmin + real(i-1) * dx
               latlist(hn) = ymin + real(j-1) * dy
            endif
          enddo
         enddo
      endif

c     mask: equidistant
      if ( hmode.eq.'mask.eqd' ) then
         hn  = 0
         lat = -90.+dy
         do while ( lat.le.ymax )
           lon = -180.
           do while ( lon.le.xmax )
              indx = nint( (lon-xmin)/dx + 1. )
              indy = nint( (lat-ymin)/dy + 1. )
              if ( indx.lt. 1 ) indx =  1
              if ( indy.lt. 1 ) indy =  1
              if ( indx.gt.nx ) indx = nx
              if ( indy.gt.ny ) indy = ny
              if ( mask(indx,indy).gt.0.5 ) then
                 hn          = hn+1
                 lonlist(hn) = lon
                 latlist(hn) = lat
              endif
              lon = lon+ds/(deltat*cos(pi180*lat))
           enddo
           lat = lat+ds/deltat
        enddo
      endif


c     ------ Get lev (vertical) coordinates -------------------------

c     Read level list from file
      if ( vmode.eq.'file' ) then
         vn = 0
         open(10,file=vfile)
 300       continue
              vn = vn + 1
              read(10,*,end=301) levlist(vn)
              goto 300
 301       continue
           vn = vn - 1
         close(10)
      endif
      
c     Get single starting level
      if ( vmode.eq.'level' ) then
         vn          = 1
         levlist(vn) = lev1
      endif
      
c     Get level profile
      if ( vmode.eq.'profile' ) then
         do i=1,vn
            levlist(i) = lev1 + real(i-1)/real(vn-1)*(lev2-lev1)
         enddo
      endif

c     Get all grid points in a layer: at the moment set the list of levels to 
c     all indices from 1 to nz; later the correct subset of indices will be chosen
      if ( vmode.eq.'grid' ) then
         vn = nz
         do i=1,vn
            levlist(i) = real(i)
         enddo
         umode_save = umode
         umode      = 'INDEX'
         
      endif

c     ------ Compile the complete list of starting positions  ------

c     Get all starting points in specified vertical coordinate system
      start_n = 0
      do i=1,vn
         do j=1,hn

            start_n = start_n + 1
            start_lon(start_n) = lonlist(j)
            start_lat(start_n) = latlist(j)
            start_lev(start_n) = levlist(i)

         enddo
      enddo

c     ------ Exit point of this section
 400  continue

c     Write status information
      print*,'  # expanded points : ', start_n
      print*
     
c     ------------------------------------------------------------------
c     Transform starting levels into pressure
c     ------------------------------------------------------------------

c     Write some status information
      print*,'---- STARTING POSITIONS ---------------------------------'
      print*

c     Vertical mode <hPa,asl> or simply <hPa>
      if ( (umode.eq.'hPa,asl').or.(umode.eq.'hPa') ) then

         do i=1,start_n
            start_pre(i) = start_lev(i)
         enddo

c     Vertical mode <hPa,agl>
      elseif ( umode.eq.'hPa,agl' ) then

         do i=1,start_n
            call get_index3(rid,rjd,rkd,start_lon(i),start_lat(i),1050.,
     >                      3,pr,prs,nx,ny,nz,xmin,ymin,dx,dy)
            tmp1 = int_index3 (prs,nx,ny,1,rid,rjd,1,mdv)
            start_pre(i) = tmp1 - start_lev(i) 
         enddo
         
c     Vertical mode <K>
      elseif ( umode.eq.'K' ) then

         do i=1,start_n
            call get_index3(rid,rjd,rkd,start_lon(i),start_lat(i),
     >                start_lev(i),1,th,ths,nx,ny,nz,xmin,ymin,dx,dy)
            tmp1 = int_index3 (pr,nx,ny,nz,rid,rjd,rkd,mdv)
            start_pre(i) = tmp1
         enddo
         
c     Vertical mode <PVU>
      elseif ( umode.eq.'PVU' ) then

         do i=1,start_n
            call get_index3(rid,rjd,rkd,start_lon(i),start_lat(i),
     >               start_lev(i),2,pv,pvs,nx,ny,nz,xmin,ymin,dx,dy)
            tmp1 = int_index3 (pr,nx,ny,nz,rid,rjd,rkd,mdv)
            start_pre(i) = tmp1
         enddo

c     Vertical mode <INDEX>
      elseif ( umode.eq.'INDEX' ) then

         do i=1,start_n
            call get_index3(rid,rjd,rkd,start_lon(i),start_lat(i),
     >               1050.,2,pv,pvs,nx,ny,nz,xmin,ymin,dx,dy)
            rkd = start_lev(i)
            tmp1 = int_index3 (pr,nx,ny,nz,rid,rjd,rkd,mdv)
            start_pre(i) = tmp1
         enddo

      endif

c     ------------------------------------------------------------------
c     Remove invalid points from the list
c     ------------------------------------------------------------------

c     Select the correct subset if <vmode=grid>: starting points outside the layer
c     will receive a <mdv> vertical pressure and will be removed
      if ( vmode.eq.'grid' ) then

         do i=1,start_n

c           Get the pressure at the grid point
            if ( ( umode_save.eq.'hPa'      ).or.
     >            (umode_save.eq.'hPa,asl') ) 
     >      then

               call get_index3(rid,rjd,rkd,start_lon(i),start_lat(i),
     >                         start_pre(i),3,pr,prs,nx,ny,nz,xmin,
     >                         ymin,dx,dy)
               tmp1 = int_index3 (pr,nx,ny,nz,rid,rjd,rkd,mdv)
               
c           Get pressure AGL at grid point
            elseif ( umode_save.eq.'hPa,agl' ) then

               call get_index3(rid,rjd,rkd,start_lon(i),
     >                         start_lat(i),start_pre(i),3,pr,prs,
     >                         nx,ny,nz,xmin,ymin,dx,dy)
               tmp1 = int_index3 (pr,nx,ny,nz,rid,rjd,rkd,mdv)
               call get_index3(rid,rjd,rkd,start_lon(i),
     >                         start_lat(i),1050.,3,pr,prs,nx,ny,
     >                         nz,xmin,ymin,dx,dy)
               tmp2 = int_index3 (prs,nx,ny,1,rid,rjd,1,mdv)
               tmp1 = tmp2 - tmp1

c           Get potential temperature at grid point
            elseif ( umode_save.eq.'K' ) then

               call get_index3(rid,rjd,rkd,start_lon(i),start_lat(i),
     >                         start_pre(i),3,pr,prs,nx,ny,nz,
     >                         xmin,ymin,dx,dy)
               tmp1 = int_index3 (th,nx,ny,nz,rid,rjd,rkd,mdv)

c           Get potential vorticity at the grid point
            elseif ( umode_save.eq.'PVU' ) then

               call get_index3(rid,rjd,rkd,start_lon(i),start_lat(i),
     >                         start_pre(i),3,pr,prs,nx,ny,nz,xmin,
     >                         ymin,dx,dy)
               tmp1 = int_index3 (pv,nx,ny,nz,rid,rjd,rkd,mdv)

c           Get vertical index at the grid point
            elseif ( umode_save.eq.'INDEX' ) then
               
              call get_index3(rid,rjd,rkd,start_lon(i),start_lat(i),
     >                        start_pre(i),3,pr,prs,nx,ny,nz,
     >                        xmin,ymin,dx,dy)
              tmp1 = int_index3 (in,nx,ny,nz,rid,rjd,rkd,mdv)

           endif

c          Remove points outside layer
           if ( ( tmp1.lt.lev1).or.(tmp1.gt.lev2) ) then
              start_pre(i) = mdv
           endif

         enddo         
         
      endif

c     Check whether the starting levels are valid (in data domain)
      do i=1,start_n

         call get_index3(rid,rjd,rkd,start_lon(i),start_lat(i),1050.,
     >                   3,pr,prs,nx,ny,nz,xmin,ymin,dx,dy)
         tmp1 = int_index3 (prs,nx,ny, 1,rid,rjd,real( 1),mdv)           ! Surface
         tmp2 = int_index3 (pr ,nx,ny,nz,rid,rjd,real(nz),mdv)           ! Top of domain

         if ( (start_pre(i).gt.tmp1).or.
     >        (start_pre(i).lt.tmp2).or.
     >        (start_lon(i).lt.xmin).or. 
     >        (start_lon(i).gt.xmax).or. 
     >        (start_lat(i).lt.ymin).or.
     >        (start_lat(i).gt.ymax) )
     >   then
            start_pre(i) = mdv
         endif

      enddo

c     Remove all starting points outside the domain
      i        = 1
      do while ( i.le.start_n )
         
         if ( abs(start_pre(i)-mdv).lt.eps ) then
            
            if ( vmode.ne.'grid') then
             print*,'  Outside ', start_lon(i),start_lat(i),start_lev(i)
            endif
             
            do j=i,start_n
               start_lon(j) = start_lon(j+1)
               start_lat(j) = start_lat(j+1)
               start_pre(j) = start_pre(j+1)
               start_lev(j) = start_lev(j+1)
            enddo
            start_n = start_n - 1

         else
         
            i = i + 1

         endif
         
      enddo

c     Write some status information
      latmin = start_lat(1)
      latmax = start_lat(1)
      lonmin = start_lon(1)
      lonmax = start_lon(1)
      premin = start_pre(1)
      premax = start_pre(1)
      do i=1,start_n
         if (start_lat(i).lt.latmin) latmin = start_lat(i)
         if (start_lat(i).gt.latmax) latmax = start_lat(i)
         if (start_lon(i).lt.lonmin) lonmin = start_lon(i)
         if (start_lon(i).gt.lonmax) lonmax = start_lon(i)
         if (start_pre(i).lt.premin) premin = start_pre(i)
         if (start_pre(i).gt.premax) premax = start_pre(i)
      enddo
      print*,'  min(lat),max(lat) : ', latmin,latmax
      print*,'  min(lon),max(lon) : ', lonmin,lonmax
      print*,'  min(pre),max(pre) : ', premin,premax
      print*
      print*,'  # starting points : ', start_n
      print*

      
c     ------------------------------------------------------------------
c     Write starting positions to output file
c     ------------------------------------------------------------------

c     Output as a trajectory file (with only one time == 0)
      if (oformat.ne.-1) then

         allocate(tra(start_n,1,5),stat=stat)

         vars(1)  ='time'
         vars(2)  ='lon'
         vars(3)  ='lat'
         vars(4)  ='p'
         vars(5)  ='level'
         call wopen_tra(fid,ofile,start_n,1,5,reftime,vars,oformat)

         do i=1,start_n
            tra(i,1,1) = 0.
            tra(i,1,2) = start_lon(i)
            tra(i,1,3) = start_lat(i)
            tra(i,1,4) = start_pre(i)
            tra(i,1,5) = start_lev(i)
         enddo
         call write_tra(fid,tra,start_n,1,5,oformat)
         
         call close_tra(fid,oformat)

c     Output as a triple list (corresponding to <startf> file)
      else
         
         fid = 10
         open(fid,file=ofile)
          do i=1,start_n
             write(fid,'(3f10.3)') start_lon(i),start_lat(i),
     >                             start_pre(i) 
          enddo
         close(fid)
         
      endif

c     Write some status information, and end of program message
      print*  
      print*,'---- STATUS INFORMATION --------------------------------'
      print*
      print*,'ok'
      print*
      print*,'       *** END OF PROGRAM CREATE_STARTF ***'
      print*,'========================================================='
        
c     ------------------------------------------------------------------
c     Exception handling
c     ------------------------------------------------------------------

      stop

 993  write(*,*) '*** ERROR: problems with array size'
      call exit(1)

      end

c     --------------------------------------------------------------------------
c     Split a region string and get corners of the domain
c     --------------------------------------------------------------------------

      subroutine regionsplit(string,iregion,xcorner,ycorner)

c     The region string comes either as <lonw,lone,lats,latn> or as <lon1,lat1,
c     lon2,lat2,lon3,lat3,lon4,lat4>: split it into ints components and get the
c     four coordinates for the region
      
      implicit none

c     Declaration of subroutine parameters
      character*80    string
      real            xcorner(4),ycorner(4)
      integer         iregion

c     Local variables
      integer         i,n
      integer         il,ir
      real            subfloat (80)
      integer         stat
      integer         len

c     ------- Split the string
      i    = 1
      n    = 0
      stat = 0
      il   = 1
      len  = len_trim(string)

 100  continue

c     Find start of a substring
      do while ( stat.eq.0 )
         if ( string(i:i).ne.' ' ) then
            stat = 1
            il   = i
         else
            i = i + 1
         endif
      enddo

c     Find end of substring
      do while ( stat.eq.1 )         
         if ( ( string(i:i).eq.' ' ) .or. ( i.eq.len ) ) then
            stat = 2
            ir   = i
         else
            i    = i + 1
         endif
      enddo

c     Convert the substring into a number
      if ( stat.eq.2 ) then
         n = n + 1
         read(string(il:ir),*) subfloat(n)
         stat = 0
      endif

      if ( i.lt.len ) goto 100


c     -------- Get the region number
      
      iregion = nint(subfloat(1))

c     -------- Get the corners of the region
      
      if ( n.eq.5 ) then     ! lonw(2),lone(3),lats(4),latn(5)

         xcorner(1) = subfloat(2)
         ycorner(1) = subfloat(4)

         xcorner(2) = subfloat(3)
         ycorner(2) = subfloat(4)
 
         xcorner(3) = subfloat(3)
         ycorner(3) = subfloat(5)
         
         xcorner(4) = subfloat(2)
         ycorner(4) = subfloat(5)
        
      elseif ( n.eq.9 ) then     ! lon1,lat1,lon2,lat2,lon3,lon4,lat4

         xcorner(1) = subfloat(2)
         ycorner(1) = subfloat(3)

         xcorner(2) = subfloat(4)
         ycorner(2) = subfloat(5)

         xcorner(3) = subfloat(6)
         ycorner(3) = subfloat(7)
         
         xcorner(4) = subfloat(8)
         ycorner(4) = subfloat(9)
 
      else
         
         print*,' ERROR: invalid region specification '
         print*,'     ',trim(string)
         stop
         
      endif
         

      end

c     --------------------------------------------------------------------------
c     Decide whether lat/lon point is in or out of region
c     --------------------------------------------------------------------------
      
      integer function inregion (lon,lat,xcorner,ycorner)
      
c     Decide whether point (lon/lat) is in the region specified by <xcorner(1..4),
c     ycorner(1..4).
      
      implicit none
      
c     Declaration of subroutine parameters
      real    lon,lat
      real    xcorner(4),ycorner(4)

c     Local variables
      integer flag
      real    xmin,xmax,ymin,ymax
      integer i

c     Reset the flag
      flag = 0

c     Set some boundaries
      xmax = xcorner(1)
      xmin = xcorner(1)
      ymax = ycorner(1)
      ymin = ycorner(1)
      do i=2,4
        if (xcorner(i).lt.xmin) xmin = xcorner(i)
        if (xcorner(i).gt.xmax) xmax = xcorner(i)
        if (ycorner(i).lt.ymin) ymin = ycorner(i)
        if (ycorner(i).gt.ymax) ymax = ycorner(i)
      enddo

c     Do the tests - set flag=1 if all tests pased
      if (lon.lt.xmin) goto 970
      if (lon.gt.xmax) goto 970
      if (lat.lt.ymin) goto 970
      if (lat.gt.ymax) goto 970
      
      if ((lon-xcorner(1))*(ycorner(2)-ycorner(1))-
     >    (lat-ycorner(1))*(xcorner(2)-xcorner(1)).gt.0.) goto 970
      if ((lon-xcorner(2))*(ycorner(3)-ycorner(2))-
     >    (lat-ycorner(2))*(xcorner(3)-xcorner(2)).gt.0.) goto 970
      if ((lon-xcorner(3))*(ycorner(4)-ycorner(3))-
     >    (lat-ycorner(3))*(xcorner(4)-xcorner(3)).gt.0.) goto 970
      if ((lon-xcorner(4))*(ycorner(1)-ycorner(4))-
     >    (lat-ycorner(4))*(xcorner(1)-xcorner(4)).gt.0.) goto 970

      flag = 1

c     Return the value
 970  continue
      
      inregion = flag
      
      return
      
      end

c     --------------------------------------------------------------------------
c     Spherical distance between lat/lon points                                                       
c     --------------------------------------------------------------------------

      real function sdis(xp,yp,xq,yq)
c
c     calculates spherical distance (in km) between two points given
c     by their spherical coordinates (xp,yp) and (xq,yq), respectively.
c
      real      re
      parameter (re=6370.)
      real      pi180
      parameter (pi180=3.14159/180.)
      real      xp,yp,xq,yq,arg

      arg=sin(pi180*yp)*sin(pi180*yq)+
     >    cos(pi180*yp)*cos(pi180*yq)*cos(pi180*(xp-xq))
      if (arg.lt.-1.) arg=-1.
      if (arg.gt.1.) arg=1.

      sdis=re*acos(arg)

      end


c     ****************************************************************
c     * Given some spherical polygon S and some point X known to be  *
c     * located inside S, these routines will determine if an arbit- *
c     * -rary point P lies inside S, outside S, or on its boundary.  *
c     * The calling program must first call DefSPolyBndry to define  *
c     * the boundary of S and the point X. Any subsequent call to    *
c     * subroutine LctPtRelBndry will determine if some point P lies *
c     * inside or outside S, or on its boundary. (Usually            *
c     * DefSPolyBndry is called once, then LctPrRelBndry is called   *
c     * many times).                                                 *
c     *                                                              * 
c     * REFERENCE:            Bevis, M. and Chatelain, J.-L. (1989)  * 
c     *                       Maflaematical Geology, vol 21.         *
c     * VERSION 1.0                                                  *
c     ****************************************************************

      Subroutine DefSPolyBndry(vlat,vlon,nv,xlat, xlon)

c     ****************************************************************
c     * This mmn entry point is used m define ~e spheric~ polygon S  *
c     * and the point X.                                             *
c     * ARGUMENTS:                                                   *
c     * vlat,vlon (sent) ... vectors containing the latitude and     * 
c     *                      longitude of each vertex of the         *
c     *                      spherical polygon S. The ith.vertex is  *
c     *                      located at [vlat(i),vlon(i)].           *
c     * nv        (sent) ... the number of vertices and sides in the *
c     *                      spherical polygon S                     *
c     * xlat,xlon (sent) ... latitude and longitude of some point X  *
c     *                      located inside S. X must not be located *
c     *                      on any great circle that includes two   *
c     *                      vertices of S.                          *
c     *                                                              *
c     * UNITS AND SIGN CONVENTION:                                   *
c     *  Latitudes and longitudes are specified in degrees.          *
c     *  Latitudes are positive to the north and negative to the     *
c     *  south.                                                      *
c     *  Longitudes are positive to the east and negative to the     *
c     *  west.                                                       *
c     *                                                              * 
c     * VERTEX ENUMERATION:                                          * 
c     * The vertices of S should be numbered sequentially around the *
c     * border of the spherical polygon. Vertex 1 lies between vertex*
c     * nv and vertex 2. Neighbouring vertices must be seperated by  *
c     * less than 180 degrees. (In order to generate a polygon side  *
c     * whose arc length equals or exceeds 180 degrees simply        *
c     * introduce an additional (pseudo)vertex). Having chosen       *
c     * vertex 1, the user may number the remaining vertices in      *
c     * either direction. However if the user wishes to use the      *
c     * subroutine SPA to determine the area of the polygon S (Bevis *
c     * & Cambareri, 1987, Math. Geol., v.19, p. 335-346) then he or *
c     * she must follow the convention whereby in moving around the  *
c     * polygon border in the direction of increasing vertex number  *
c     * clockwise bends occur at salient vertices. A vertex is       *
c     * salient if the interior angle is less than 180 degrees.      *
c     * (In the case of a convex polygon this convention implies     *
c     * that vertices are numbered in clockwise sequence).           *
c     ****************************************************************

      implicit none
      
      integer mxnv,nv

c     ----------------------------------------------------------------
c     Edit next statement to increase maximum number of vertices that 
c     may be used to define the spherical polygon S               
c     The value of parameter mxnv in subroutine LctPtRelBndry must match
c     that of parameter mxnv in this subroutine, as assigned above.
c     ----------------------------------------------------------------
      parameter (mxnv=500)

      real  vlat(nv),vlon(nv),xlat,xlon,dellon
      real  tlonv(mxnv),vlat_c(mxnv),vlon_c(mxnv),xlat_c,xlon_c
      integer i,ibndry,nv_c,ip
 
      data ibndry/0/
      
      common/spolybndry/vlat_c,vlon_c,nv_c,xlat_c,xlon_c,tlonv,ibndry

      if (nv.gt.mxnv) then
         print *,'nv exceeds maximum allowed value'
         print *,'adjust parameter mxnv in subroutine DefSPolyBndry'
         stop
      endif

      ibndry=1                  ! boundary defined at least once (flag)
      nv_c=nv                   ! copy for named common
      xlat_c=xlat               ! . . . .
      xlon_c=xlon               !

      do i=1,nv
         vlat_c(i)=vlat(i)      ! "
         vlon_c(i)=vlon(i)      !

         call TrnsfmLon(xlat,xlon,vlat(i),vlon(i),tlonv(i))

         if (i.gt.1) then
            ip=i-1
         else
            ip=nv
         endif
         
         if ((vlat(i).eq.vlat(ip)).and.(vlon(i).eq.vlon(ip))) then
            print *,'DefSPolyBndry detects user error:'
            print *,'vertices ',i,' and ',ip,' are not distinct'
            print*,'lat ',i,ip,vlat(i),vlat(ip)
            print*,'lon ',i,ip,vlon(i),vlon(ip)            
            stop
         endif

         if (tlonv(i).eq.tlonv(ip)) then
            print *,'DefSPolyBndry detects user error:'
            print *,'vertices ',i,' & ',ip,' on same gt. circle as X'
            stop
         endif

         if (vlat(i).eq.(-vlat(ip))) then
            dellon=vlon(i)-vlon(ip)
            if (dellon.gt.+180.) dellon=dellon-360.
            if (dellon.lt.-180.) dellon=dellon-360.
            if ((dellon.eq.+180.0).or.(dellon.eq.-180.0)) then
               print *,'DefSPolyBndry detects user error:'
               print *,'vertices ',i,' and ',ip,' are antipodal'
               stop
            endif
         endif
      enddo

      return
      
      end


c     ****************************************************************
 
      Subroutine LctPtRelBndry(plat,plon,location)

c     ****************************************************************

c     ****************************************************************
c     * This routine is used to see if some point P is located       *
c     * inside, outside or on the boundary of the spherical polygon  *
c     * S previously defined by a call to subroutine DefSPolyBndry.  *
c     * There is a single restriction on point P: it must not be     *
c     * antipodal to the point X defined in the call to DefSPolyBndry*
c     * (ie.P and X cannot be seperated by exactly 180 degrees).     *
c     * ARGUMENTS:                                                   *  
c     * plat,plon (sent)... the latitude and longitude of point P    *
c     * location (returned)... specifies the location of P:          *
c     *                        location=0 implies P is outside of S  *
c     *                        location=1 implies P is inside of S   *
c     *                        location=2 implies P on boundary of S *
c     *                        location=3 implies user error (P is   *
c     *                                     antipodal to X)          *
c     * UNFfS AND SIGN CONVENTION:                                   * 
c     *  Latitudes and longitudes are specified in degrees.          *
c     *  Latitudes are positive to the north and negative to the     *
c     *  south.                                                      *    
c     *  Longitudes are positive to the east and negative to the     *
c     *  west.                                                       *
c     ****************************************************************
      
      implicit none
      
      integer mxnv

c     ----------------------------------------------------------------
c     The statement below must match that in subroutine DefSPolyBndry
c     ----------------------------------------------------------------

      parameter (mxnv=500)

      real tlonv(mxnv),vlat_c(mxnv),vlon_c(mxnv),xlat_c,xlon_c
      real plat,plon,vAlat,vAlon,vBlat,vBlon,tlonA,tlonB,tlonP
      real tlon_X,tlon_P,tlon_B,dellon
      integer i,ibndry,nv_c,location,icross,ibrngAB,ibrngAP,ibrngPB
      integer ibrng_BX,ibrng_BP,istrike

      common/spolybndry/vlat_c,vlon_c,nv_c,xlat_c,xlon_c,tlonv,ibndry

      if (ibndry.eq.0) then     ! user has never defined the bndry
         print*,'Subroutine LctPtRelBndry detects user error:'
         print*,'Subroutine DefSPolyBndry must be called before'
         print*,'subroutine LctPtRelBndry can be called'
         stop
      endif

      if (plat.eq.(-xlat_c)) then
         dellon=plon-xlon_c
         if (dellon.lt.(-180.)) dellon=dellon+360.
         if (dellon.gt.+180.) dellon=dellon-360.
         if ((dellon.eq.+180.0).or.(dellon.eq.-180.)) then
            print*,'Warning: LctPtRelBndry detects case P antipodal
     >           to X'
            print*,'location of P relative to S is undetermined'
            location=3
            return
         endif
      endif 

      location=0                ! default ( P is outside S)
      icross=0                  ! initialize counter

      if ((plat.eq.xlat_c).and.(plon.eq.xlon_c)) then
         location=1
         return
      endif

      
      call TrnsfmLon (xlat_c,xlon_c,plat,plon,tlonP)

      do i=1,nv_c              ! start of loop over sides of S 

         vAlat=vlat_c(i)
         vAlon=vlon_c(i)
         tlonA=tlonv(i)

         if (i.lt.nv_c) then
            vBlat=vlat_c(i+1)
            vBlon=vlon_c(i+1)
            tlonB=tlonv(i+1)
         else
            vBlat=vlat_c(1)
            vBlon=vlon_c(1)
            tlonB=tlonv(1)
         endif
         
         istrike=0
         
         if (tlonP.eq.tlonA) then
            istrike=1
         else
            call EastOrWest(tlonA,tlonB,ibrngAB)
            call EastOrWest(tlonA,tlonP,ibrngAP)
            call EastOrWest(tlonP,tlonB,ibrngPB)
            

            if((ibrngAP.eq.ibrngAB).and.(ibrngPB.eq.ibrngAB)) istrike=1
         endif

         
         if (istrike.eq.1) then

            if ((plat.eq.vAlat).and.(plon.eq.vAlon)) then
               location=2       ! P lies on a vertex of S
               return
            endif
            call TrnsfmLon(vAlat,vAlon,xlat_c,xlon_c,tlon_X)
            call TrnsfmLon(vAlat,vAlon,vBlat,vBlon,tlon_B)
            call TrnsfmLon(vAlat,vAlon,plat,plon,tlon_P)
            
            if (tlon_P.eq.tlon_B) then
               location=2       ! P lies on side of S
               return 
            else
               call EastOrWest(tlon_B,tlon_X,ibrng_BX)
               call EastOrWest(tlon_B,tlon_P,ibrng_BP)
               if(ibrng_BX.eq.(-ibrng_BP)) icross=icross+1
            endif
            
         endif
      enddo                     ! end of loop over the sides of S


c     if the arc XP crosses the boundary S an even number of times then P
c     is in S

      if (mod(icross,2).eq.0) location=1

      return

      end


c     ****************************************************************
      
      subroutine TrnsfmLon(plat,plon,qlat,qlon,tranlon)

c     ****************************************************************
c     * This subroutine is required by subroutines DefSPolyBndry &   *
c     * LctPtRelBndry. It finds the 'longitude' of point Q in a      *
c     * geographic coordinate system for which point P acts as a     *
c     * 'north pole'. SENT: plat,plon,qlat,qlon, in degrees.         *
c     * RETURNED: tranlon, in degrees.                               *
c     ****************************************************************

      implicit none

      real pi,dtr,plat,plon,qlat,qlon,tranlon,t,b
      parameter (pi=3.141592654,dtr=pi/180.0)
 
      if (plat.eq.90.) then
         tranlon=qlon
      else
         t=sin((qlon-plon)*dtr)*cos(qlat*dtr)
         b=sin(dtr*qlat)*cos(plat*dtr)-cos(qlat*dtr)*sin(plat*dtr)
     >    *cos((qlon-plon)*dtr)
         tranlon=atan2(t,b)/dtr
      endif

      return
      end

c     ****************************************************************

      subroutine EastOrWest(clon,dlon,ibrng)

c     ****************************************************************
c     * This subroutine is required by subroutine LctPtRelBndry.     *
c     * This routine determines if in travelling the shortest path   *
c     * from point C (at longitude clon) to point D (at longitude    *
c     * dlon) one is heading east, west or neither.                  *
c     * SENT: clon,dlon; in degrees. RETURNED: ibrng                 *
c     * (1=east,-1=west, 0=neither).                                 *
c     ****************************************************************

      implicit none
      real clon,dlon,del
      integer ibrng
      del=dlon-clon
      if (del.gt.180.) del=del-360.
      if (del.lt.-180.) del=del+360.
      if ((del.gt.0.0).and.(del.ne.180.)) then
         ibrng=-1               ! (D is west of C)
      elseif ((del.lt.0.0).and.(del.ne.-180.)) then
         ibrng=+1               ! (D is east of C)
      else
         ibrng=0                ! (D north or south of C)
      endif
      return
      end
