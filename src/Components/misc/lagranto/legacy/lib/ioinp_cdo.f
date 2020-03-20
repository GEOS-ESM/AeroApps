c     ************************************************************
c
c                   Customized for GEOS-5!!!!!
c      
c     * This package provides input routines to read the wind    *
c     * and other fields from IVE necdf files. The routines are  *
c     *                                                          *
c     * 1) input_open  : to open a data file                     *
c     * 2) input_grid  : to read the grid information, including *
c     *                  the vertical levels                     *
c     * 3) input_wind  : to read the wind components             *
c     * 4) input_close : to close an input file                  *
c     *                                                          *
c     * The file is characterised by an filename <filename> and  *
c     * a file identifier <fid>. The horizontal grid is given by *
c     * <xmin,xmax,ymin,ymax,dx,dy,nx,ny> where the pole of the  *
c     * rotated grid is given by <pollon,pollat>. The vertical   *
c     * grid is characterised by the surface pressure <ps> and   *
c     * the pressure at staggered <slev> and unstaggered <ulev>  *
c     * levels. The number of levels is given by <nz>. Finally,  *
c     * the retrieval of the wind <field> with name <fieldname>  *
c     * is characterised by a <time> and a missing data value    *
c     * <mdv>.                                                   *
c     *                                                          *
c     * Author: Michael Sprenger, Autumn 2008                    *
c     ************************************************************

      subroutine get_eta ( ak, bk, nz )

      integer, intent(in) :: nz
      real, intent(out) :: ak(nz) ! mid layer, in hPa
      real, intent(out) :: bk(nz) ! mid layer, non-dimensional
      
        real a72(73), b72(73)

        data a72 / 
     .    1.0000000, 2.0000002, 3.2700005, 4.7585009, 6.6000011, 
     .    8.9345014, 11.970302, 15.949503, 21.134903, 27.852606, 
     .    36.504108, 47.580610, 61.677911, 79.513413, 101.94402, 
     .    130.05102, 165.07903, 208.49704, 262.02105, 327.64307, 
     .    407.65710, 504.68010, 621.68012, 761.98417, 929.29420, 
     .    1127.6902, 1364.3402, 1645.7103, 1979.1604, 2373.0405, 
     .    2836.7806, 3381.0007, 4017.5409, 4764.3911, 5638.7912, 
     .    6660.3412, 7851.2316, 9236.5722, 10866.302, 12783.703, 
     .    15039.303, 17693.003, 20119.201, 21686.501, 22436.301, 
     .    22389.800, 21877.598, 21214.998, 20325.898, 19309.696, 
     .    18161.897, 16960.896, 15625.996, 14290.995, 12869.594, 
     .    11895.862, 10918.171, 9936.5219, 8909.9925, 7883.4220, 
     .    7062.1982, 6436.2637, 5805.3211, 5169.6110, 4533.9010, 
     .    3898.2009, 3257.0809, 2609.2006, 1961.3106, 1313.4804, 
     .    659.37527, 4.8048257, 0.0000000 /


        data b72 / 
     .    0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
     .    0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
     .    0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
     .    0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
     .    0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
     .    0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
     .    0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
     .    0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
     .    0.0000000, 8.1754130e-09,0.0069600246,0.028010041,0.063720063, 
     .   0.11360208,0.15622409,0.20035011,0.24674112,0.29440312, 
     .   0.34338113,0.39289115,0.44374018,0.49459020,0.54630418, 
     .   0.58104151,0.61581843,0.65063492,0.68589990,0.72116594, 
     .   0.74937819,0.77063753,0.79194696,0.81330397,0.83466097, 
     .   0.85601798,0.87742898,0.89890800,0.92038701,0.94186501, 
     .   0.96340602,0.98495195, 1.0000000 /

         if ( nz .ne. 72 ) then
            print *, 'Number of vertical layers must be 72!!!!'
            stop
         end if

         do k = 1, nz
            ak(k) = 0.5 * ( a72(k) + a72(k+1) ) / 100.
            bk(k) = 0.5 * ( b72(k) + b72(k+1) )
         end do

         return
         
      end

      
      
c     ------------------------------------------------------------
c     Open input file
c     ------------------------------------------------------------

      subroutine input_open (fid,filename)

c     Open the input file with filename <filename> and return the
c     file identifier <fid> for further reference. 

      use netcdf

      implicit none

c     Declaration of subroutine parameters
      integer      fid              ! File identifier
      character*80 filename         ! Filename

c     Declaration of auxiliary variables
      integer      ierr

c     Open netcdf file
      ierr = NF90_OPEN(TRIM(filename),nf90_nowrite, fid)
      IF ( ierr /= nf90_NoErr ) PRINT *,NF90_STRERROR(ierr)

c     Exception handling
      return

      end
      
c     ------------------------------------------------------------
c     Read information about the grid
c     ------------------------------------------------------------
      
      subroutine input_grid 
     >                   (fid,fieldname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >                    time,pollon,pollat,p3,ps,nz,ak,bk,stagz,
     >                    timecheck)

c     Read grid information at <time> from file with identifier <fid>. 
c     The horizontal grid is characterized by <xmin,xmax,ymin,ymax,dx,dy>
c     with pole position at <pollon,pollat> and grid dimension <nx,ny>.
c     The 3d arrays <p3(nx,ny,nz)> gives the vertical coordinates, either
c     on the staggered or unstaggered grid (with <stagz> as the flag).
c     The surface pressure is given in <ps(nx,ny)>. If <fid> is negative, 
c     only the grid dimensions and grid parameters (xmin...pollat,nz) are 
c     determined and returned (this is needed for dynamical allocation of 
c     memory).

      use netcdf

      implicit none

c     Declaration of subroutine parameters 
      integer      fid                 ! File identifier
      real         xmin,xmax,ymin,ymax ! Domain size
      real         dx,dy               ! Horizontal resolution
      integer      nx,ny,nz            ! Grid dimensions
      real         pollon,pollat       ! Longitude and latitude of pole
      real         p3(nx,ny,nz)        ! Staggered levels
      real         ps(nx,ny)           ! Surface pressure
      real         time                ! Time of the grid information
      real         ak(nz),bk(nz)       ! Ak and Bk for layers or levels
      real         stagz               ! Vertical staggering (0 or -0.5)
      character*80 fieldname           ! Variable from which to take grid info
      character*80 timecheck           ! Either 'yes' or 'no'
      
c     Numerical and physical parameters
      real          eps                 ! Numerical epsilon
      parameter    (eps=0.001)

c     Netcdf variables
      integer      vardim(4)
      real         varmin(4),varmax(4)
      real         mdv
      real         stag(4)
      integer      ndim
      character*80 cstfile
      integer      cstid
      integer      nvars
      character*80 vars(100)
      integer        dimids (nf90_max_var_dims),dimid
      character*80   dimname(nf90_max_var_dims)
      character*80   stdname
      real,allocatable, dimension (:)     :: lon,lat,lev
      real,allocatable, dimension (:)     :: times
      real,allocatable, dimension (:,:)   :: tmp2
      real,allocatable, dimension (:,:,:) :: tmp3
      real,allocatable, dimension (:)     :: aktmp,bktmp
      character*80  units
      character*80  leveltype
      integer       nakbktmp
      integer       vertical_swap

c     Auxiliary variables
      integer      ierr       
      integer      i,j,k
      integer      isok
      real         tmp(200)
      character*80 varname
      real         rtime
      integer      varid
      integer      cdfid
      integer      stat
      real         delta
      integer      closear
      real         maxps,minps

c     ---- Read data from netCDF file as they are ---------------------

c     Set file identifier
      if (fid.lt.0) then
        cdfid = -fid
      else 
        cdfid = fid
      endif

c     Special handling if 3D pressure is
      if ( fieldname.eq.'P' ) then
         fieldname = 'U'
      endif

c     Get number of dimensions of variable -> ndim
      ierr = NF90_INQ_VARID(cdfid,fieldname,varid)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_inquire_variable(cdfid, varid, ndims  = ndim)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      if ( ndim.ne.4 ) then
          print*,' ERROR: netCDF variables need to be 4D'
          print*,'      ',trim(fieldname)
          stop
      endif

c     Get dimensions -> vardim(1:ndim),dimname(1:ndim)
      ierr = NF90_INQ_VARID(cdfid,fieldname,varid)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_inquire_variable(cdfid, varid, 
     >                                   dimids = dimids(1:ndim))
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      do i=1,ndim
           ierr = nf90_inquire_dimension(cdfid, dimids(i), 
     >                               name = dimname(i) )
           IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
           ierr = nf90_inquire_dimension(cdfid, dimids(i),len=vardim(i))
           IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      enddo

c     Get dimension of AK,BK
ccc      varname = 'nhym'
cccc      ierr = NF90_INQ_DIMID(cdfid,varname,dimid)
ccc      ierr = nf90_inquire_dimension(cdfid, dimid,len=nakbktmp)
cccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      nakbktmp = vardim(3)
      
c     Check whether the list of dimensions is OK
      if ( ( dimname(1).ne.'lon'  ).or.
     >     ( dimname(2).ne.'lat'  ).or. 
     >     ( dimname(3).ne.'lev'  ).and.( dimname(3).ne.'lev_2'  ).or.
     >     ( dimname(4).ne.'time' ) )
     >then
        print*,' ERROR: the dimensions of the variable are not correct'
        print*,'        expected -> lon / lat / lev / time'
        print*, ( trim(dimname(i))//' / ',i=1,ndim )
        stop
      endif

c     Check about the type of vertical levels
ccc      varname=dimname(3)
ccc      ierr = NF90_INQ_VARID(cdfid,varname,varid)
ccc      IF(ierr /= nf90_NoErr) PRINT *,trim(varname),NF90_STRERROR(ierr)
ccc      ierr = nf90_get_att(cdfid, varid, "standard_name", leveltype)
ccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
ccc      if ( (leveltype.ne.'hybrid_sigma_pressure').and.
ccc     >     (leveltype.ne.'air_pressure'         ) )
ccc     >then
ccc         print*,' ERROR: input netCDF data must be on hybrid-sigma'
ccc         print*,'        or air pressure levels!',trim(leveltype)

ccc          print *, 'WARNING: Assuming hybrid sigma-pressure coordinates'
         leveltype = 'hybrid_sigma_pressure'

ccc      endif

c     Allocate memory for reading arrays
      allocate(tmp2(vardim(1),vardim(2)),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tmp2     ***'
      allocate(tmp3(vardim(1),vardim(2),vardim(3)),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tmp3     ***'
      allocate(lon(vardim(1)),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array lon     ***' 
      allocate(lat(vardim(2)),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array lat     ***' 
      allocate(lev(vardim(3)),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array lev     ***'
      allocate(times(vardim(4)),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array times   ***'
      allocate(aktmp(nakbktmp),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array aktmp   ***'
      allocate(bktmp(nakbktmp),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array bktmp   ***'

c     Get domain longitudes, latitudes and levels
      varname = dimname(1)
      ierr = NF90_INQ_VARID(cdfid,varname,varid)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_get_var(cdfid,varid,lon)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      varname = dimname(2)
      ierr = NF90_INQ_VARID(cdfid,varname,varid)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_get_var(cdfid,varid,lat)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      varname = dimname(3)
      ierr = NF90_INQ_VARID(cdfid,varname,varid)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_get_var(cdfid,varid,lev)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      
c     Get ak and bk
ccc      varname='hyam'
ccc      ierr = NF90_INQ_VARID(cdfid,varname,varid)
ccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
ccc      ierr = nf90_get_var(cdfid,varid,aktmp)
ccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
ccc      varname='hybm'
ccc      ierr = NF90_INQ_VARID(cdfid,varname,varid)
ccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
ccc      ierr = nf90_get_var(cdfid,varid,bktmp)
ccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      call get_eta(aktmp, bktmp, vardim(3) )
      
      
c     Check that unit of ak is in hPa - if necessary correct it
ccc      varname='hyam'
ccc      ierr = NF90_INQ_VARID(cdfid,varname,varid)
ccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
ccc      ierr = nf90_get_att(cdfid, varid, "units", units)
ccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
ccc      if ( units.eq.'Pa' ) then
ccc         do k=1,nakbktmp
ccc            aktmp(k) = 0.01 * aktmp(k)
ccc         enddo
ccc      endif

c     Check that unit of lev is in hPa - if necessary correct it
      if ( leveltype.eq.'air_pressure' ) then
         varname='lev'
         ierr = NF90_INQ_VARID(cdfid,varname,varid)
         IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
         ierr = nf90_get_att(cdfid, varid, "units", units)
         IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
         if ( units.eq.'Pa' ) then
            do k=1,vardim(3)
               lev(k) = 0.01 * lev(k)
            enddo
         endif
      endif

c     Decide whether to swap vertical levels - highest pressure at index 1
      vertical_swap = 1
      if ( leveltype.eq.'hybrid_sigma_pressure') then
        if ( (aktmp(1) + bktmp(1) * 1000.).gt.
     >       (aktmp(2) + bktmp(2) * 1000.) )
     >  then
          vertical_swap = 0
        endif
      elseif ( leveltype.eq.'air_pressure') then
        if ( lev(1).gt.lev(2) ) then
          vertical_swap = 0
        endif
      endif
ccc      print*,' Vertical SWAP P -> ', vertical_swap

c     Get time information (check if time is correct)
      varname = 'time'
      ierr = NF90_INQ_VARID(cdfid,varname,varid)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_get_var(cdfid,varid,times)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      isok=0
      do i=1,vardim(4)
        if (abs(time-times(i)).lt.eps) then
               isok = 1
               rtime = times(i)
        elseif (timecheck.eq.'no') then
               isok = 1
               rtime = times(1)
        endif
      enddo
      if ( isok.eq.0 ) then
         print*,' ERROR: time ',rtime,' not found on netCDF file' 
         stop
      endif

c     Read surface pressure
      varname='PS'
      ierr = NF90_INQ_VARID(cdfid,varname,varid)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_get_var(cdfid,varid,tmp2)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      tmp2 = tmp2 / 100. ! from Pa to hPa (GEOS-5)

      
c     Check that surface pressure is in hPa
      maxps = -1.e19
      minps =  1.e19
      do i=1,vardim(1)
        do j=1,vardim(2)
             if (tmp2(i,j).gt.maxps) maxps = tmp2(i,j)
             if (tmp2(i,j).lt.minps) minps = tmp2(i,j)
        enddo
      enddo
      if ( (maxps.gt.1500.).or.(minps.lt.300.) ) then
         print*,' ERROR: surface pressre PS must be in hPa'
         print*,'       ',maxps,minps
         stop
      endif

c     ---- Define output of subroutine --------------------------------

c     If not full list of vertical levels, reduce AK,BK arrays
      if ( (leveltype.eq.'hybrid_sigma_pressure').and.
     >     (nakbktmp.ne.vardim(3) ) )
     >then
c         print*,' WARNING: only subset of vertical levels used...'
         do k=1,vardim(3)
            if ( vertical_swap.eq.1 ) then
               aktmp(k) = aktmp( k+nakbktmp-vardim(3) )
               bktmp(k) = bktmp( k+nakbktmp-vardim(3) )
            endif
         enddo
      endif

c     Set the grid dimensions and constants
      nx      = vardim(1)
      ny      = vardim(2)
      nz      = vardim(3)
      xmin    = lon(1)
      ymin    = lat(1)
      xmax    = lon(nx)
      ymax    = lat(ny)
      dx      = (xmax-xmin)/real(nx-1)
      dy      = (ymax-ymin)/real(ny-1)
      pollon  = 0.
      pollat  = 90.
      stagz   = 0.
      delta   = xmax-xmin-360.
      if (abs(delta+dx).lt.eps) then
          xmax    = xmax + dx
          nx      = nx + 1
          closear = 1
      else
          closear = 0
      endif

c     Save the output arrays (if fid>0) - close arrays on request
      if ( fid.gt.0 ) then

c        Calculate layer pressures
         if (leveltype.eq.'hybrid_sigma_pressure' ) then
            do i=1,vardim(1)
             do j=1,vardim(2)
                 do k=1,vardim(3)
                  tmp3(i,j,k)=aktmp(k)+bktmp(k)*tmp2(i,j)
                enddo
              enddo
          enddo
           
         elseif (leveltype.eq.'air_pressure' ) then
           do i=1,vardim(1)
              do j=1,vardim(2)
                 do k=1,vardim(3)
                  tmp3(i,j,k)=lev(k)
                 enddo
              enddo
           enddo
         endif

c        Get PS - close array on demand
         do j=1,vardim(2)
           do i=1,vardim(1)
             ps(i,j) = tmp2(i,j)
           enddo
           if (closear.eq.1) ps(vardim(1)+1,j) = ps(1,j)
         enddo

c        Get P3 - close array on demand + vertical swap
         do j=1,vardim(2)
           do k=1,vardim(3)
             do i=1,vardim(1)
               if ( vertical_swap.eq.1 ) then
                  p3(i,j,k) = tmp3(i,j,vardim(3)-k+1)
               else
                  p3(i,j,k) = tmp3(i,j,k)
               endif
             enddo
             if (closear.eq.1) p3(vardim(1)+1,j,k) = p3(1,j,k)
           enddo
         enddo

c        Get AK,BK - vertical swap on demand
         if ( leveltype.eq.'hybrid_sigma_pressure' ) then
           do k=1,vardim(3)
              if ( vertical_swap.eq.1 ) then
                 ak(k) = aktmp(vardim(3)-k+1)
                 bk(k) = bktmp(vardim(3)-k+1)
              else
                 ak(k) = aktmp(k)
                 bk(k) = bktmp(k)
              endif
           enddo
         elseif (leveltype.eq.'air_pressure' ) then
           do k=1,vardim(3)
              if ( vertical_swap.eq.1 ) then
                 ak(k) = lev(vardim(3)-k+1)
                 bk(k) = 0.
              else
                ak(k) = lev(k)
                bk(k) = 0.
              endif
           enddo
         endif

      endif


      return
      
      end

c     ------------------------------------------------------------
c     Read wind information
c     ------------------------------------------------------------

      subroutine input_wind (fid,fieldname,field,time,stagz,mdv,
     >                       xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,
     >                       timecheck)

c     Read the wind component <fieldname> from the file with identifier
c     <fid> and save it in the 3d array <field>. The vertical staggering 
c     information is provided in <stagz> and gives the reference to either
c     the layer or level field from <input_grid>. A consistency check is
c     performed to have an agreement with the grid specified by <xmin,xmax,
c     ymin,ymax,dx,dy,nx,ny,nz>.

      use netcdf

      implicit none

c     Declaration of variables and parameters
      integer      fid                 ! File identifier
      character*80 fieldname           ! Name of the wind field
      integer      nx,ny,nz            ! Dimension of fields
      real         field(nx,ny,nz)     ! 3d wind field
      real         stagz               ! Staggering in the z direction
      real         mdv                 ! Missing data flag
      real         xmin,xmax,ymin,ymax ! Domain size
      real         dx,dy               ! Horizontal resolution
      real         time                ! Time
      character*80 timecheck           ! Either 'yes' or 'no'

c     Numerical and physical parameters
      real        eps                 ! Numerical epsilon
      parameter  (eps=0.001)
      real        notimecheck         ! 'Flag' for no time check
      parameter  (notimecheck=7.26537)

c     Netcdf variables
      integer      ierr
      character*80 varname
      integer      vardim(4)
      real         varmin(4),varmax(4)
      real         stag(4)
      integer      ndim
      real         times(10)
      integer      ntimes
      character*80 cstfile
      integer      cstid
      real         aklay(200),bklay(200),aklev(200),bklev(200)
      real         ps(nx,ny)
      integer      dimids (nf90_max_var_dims)
      character*80 dimname(nf90_max_var_dims)
      integer      varid
      integer      cdfid
      real,allocatable, dimension (:)     :: lon,lat,lev
      real,allocatable, dimension (:,:)   :: tmp2
      real,allocatable, dimension (:,:,:) :: tmp3
      real,allocatable, dimension (:)     :: aktmp,bktmp
      character*80  leveltype
      integer       vertical_swap
      character*80  units
      integer       nakbktmp
      integer       dimid

c     Auxiliary variables
      integer      isok
      integer      i,j,k
      integer      nz1
      real         rtime
      integer      closear
      integer      stat
      real         delta
      integer      pressure

c     Init mdv
      mdv = -999.
      
c     Special handling if 3D pressure is
      if ( fieldname.eq.'P' ) then
         fieldname = 'U'
         pressure  = 1
      else
         pressure = 0
      endif

c     Get the number of dimensions -> ndim
      ierr = NF90_INQ_VARID(fid,fieldname,varid)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_inquire_variable(fid, varid, ndims  = ndim)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)

c     Get the dimensions of the arrays -> varid(1:ndim)
      ierr = NF90_INQ_VARID(fid,fieldname,varid)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_inquire_variable(fid, varid, 
     >                                   dimids = dimids(1:ndim))
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      do i=1,ndim
           ierr = nf90_inquire_dimension(fid, dimids(i), 
     >                               name = dimname(i) )
           IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
           ierr = nf90_inquire_dimension(fid, dimids(i),len=vardim(i))
           IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      enddo

c     Check whether the list of dimensions is OK
      if ( ( dimname(1).ne.'lon'  ).or.
     >     ( dimname(2).ne.'lat'  ).or.
     >     ( dimname(3).ne.'lev'  ).and.( dimname(3).ne.'lev_2'  ).or.
     >     ( dimname(4).ne.'time' ) )
     >then
        print*,' ERROR: the dimensions of the variable are not correct'
        print*,'        expected -> lon / lat / lev / time'
        print*, ( trim(dimname(i))//' / ',i=1,ndim )
        stop
      endif

c     Get dimension of AK,BK
ccc      varname = 'nhym'
ccc      ierr = NF90_INQ_DIMID(fid,varname,dimid)
ccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
ccc      ierr = nf90_inquire_dimension(fid, dimid,len=nakbktmp)
ccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      nakbktmp = vardim(3)
      
c     Check about the type of vertical levels
ccc      varname=dimname(3)
ccc      ierr = NF90_INQ_VARID(fid,varname,varid)
ccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
ccc      ierr = nf90_get_att(fid, varid, "standard_name", leveltype)
ccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
ccc      if ( (leveltype.ne.'hybrid_sigma_pressure').and.
ccc     >     (leveltype.ne.'air_pressure'         ) )
ccc     >then
ccc         print*,' ERROR: input netCDF data must be on hybrid-sigma'
ccc         print*,'        or air pressure levels!',trim(leveltype)
ccc   stop
ccc         print *, 'WARNING: Assuming hybrid sigma-pressure coordinates'
         leveltype = 'hybrid_sigma_pressure'
ccc      endif

c     Allocate memory for reading arrays - depending on <closear>
      allocate(tmp2(vardim(1),vardim(2)),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tmp2     ***'
      allocate(tmp3(vardim(1),vardim(2),vardim(3)),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tmp3     ***'
      allocate(lon(vardim(1)),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array lon     ***'
      allocate(lat(vardim(2)),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array lat     ***'
      allocate(lev(vardim(3)),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array lev     ***'
      allocate(aktmp(nakbktmp),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array aktmp   ***'
      allocate(bktmp(nakbktmp),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array bktmp   ***'

c     Get domain boundaries - longitude, latitude, levels
      varname = dimname(1)
      ierr = NF90_INQ_VARID(fid,varname,varid)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_get_var(fid,varid,lon)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      varname = dimname(2)
      ierr = NF90_INQ_VARID(fid,varname,varid)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_get_var(fid,varid,lat)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      varname = dimname(3)
      ierr = NF90_INQ_VARID(fid,varname,varid)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_get_var(fid,varid,lev)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      varmin(1) = lon(1)
      varmax(1) = lon( vardim(1) )
      varmin(2) = lat(1)
      varmax(2) = lat( vardim(2) )

c     Get ak and bk
ccc      varname='hyam'
ccc      ierr = NF90_INQ_VARID(fid,varname,varid)
ccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
ccc      ierr = nf90_get_var(fid,varid,aktmp)
ccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
ccc      varname='hybm'
ccc      ierr = NF90_INQ_VARID(fid,varname,varid)
ccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
ccc      ierr = nf90_get_var(fid,varid,bktmp)
ccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      call get_eta(aktmp, bktmp, vardim(3) )

c     Check that unit of ak is in hPa - if necessary correct it
ccc      varname='hyam'
ccc      ierr = NF90_INQ_VARID(fid,varname,varid)
ccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
ccc      ierr = nf90_get_att(fid, varid, "units", units)
ccc      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
ccc      if ( units.eq.'Pa' ) then
ccc         do k=1,nakbktmp
ccc            aktmp(k) = 0.01 * aktmp(k)
ccc         enddo
ccc      endif

c     Check that unit of lev is in hPa - if necessary correct it
      if ( leveltype.eq.'air_pressure' ) then
         varname='lev'
         ierr = NF90_INQ_VARID(fid,varname,varid)
         IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
         ierr = nf90_get_att(fid, varid, "units", units)
         IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
         if ( units.eq.'Pa' ) then
            do k=1,vardim(3)
               lev(k) = 0.01 * lev(k)
            enddo
         endif
      endif

c     Decide whether to swap vertical levels
      vertical_swap = 1
      if ( leveltype.eq.'hybrid_sigma_pressure') then
        if ( (aktmp(1) + bktmp(1) * 1000.).gt.
     >       (aktmp(2) + bktmp(2) * 1000.) )
     >  then
          vertical_swap = 0
        endif
      elseif ( leveltype.eq.'air_pressure') then
        if ( lev(1).gt.lev(2) ) then
          vertical_swap = 0
        endif
      endif

c     Read data /or set 3D pressure field) 
      if ( pressure.eq.0 ) then
         ierr = NF90_INQ_VARID(fid,fieldname,varid)
         IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
         ierr = nf90_get_var(fid,varid,tmp3)
         IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      else
         if (leveltype.eq.'hybrid_sigma_pressure' ) then
            do i=1,vardim(1)
              do j=1,vardim(2)
                 do k=1,vardim(3)
                  tmp3(i,j,k)=aktmp(k)+bktmp(k)*tmp2(i,j)
                 enddo
              enddo
           enddo
         elseif (leveltype.eq.'air_pressure' ) then
           do i=1,vardim(1)
              do j=1,vardim(2)
                 do k=1,vardim(3)
                  tmp3(i,j,k)=lev(k)
                 enddo
              enddo
           enddo
         endif
      endif
  
c     If the field is 2D, expand it to 3D - simple handling of 2D tracing
      if ( vardim(3).eq.1 ) then
         do i=1,vardim(1)
            do j=1,vardim(2)
               do k=1,vardim(3)
                  tmp3(i,j,k) = tmp3(i,j,1)
               enddo
            enddo
         enddo
      endif

c     Decide whether to close arrays
      delta = varmax(1)-varmin(1)-360.
      if (abs(delta+dx).lt.eps) then
          closear = 1
      else
          closear = 0
      endif

c     Save output array - close array and swap on demand
      do j=1,vardim(2)
        do k=1,vardim(3)
          do i=1,vardim(1)
             if ( vertical_swap.eq.1 ) then
                 field(i,j,k) = tmp3(i,j,vardim(3)-k+1)
             else
                 field(i,j,k) = tmp3(i,j,k)
             endif
          enddo
          if (closear.eq.1) field(vardim(1)+1,j,k) = field(1,j,k)
        enddo
      enddo
         
c     Exit point
      return
 
      end

c     ------------------------------------------------------------
c     Close input file
c     ------------------------------------------------------------

      subroutine input_close(fid)

c     Close the input file with file identifier <fid>.

      use netcdf

      implicit none

c     Declaration of subroutine parameters
      integer fid

c     Auxiliary variables
      integer ierr

c     Close file
      ierr = NF90_CLOSE(fid)
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      
      end
      
c     ------------------------------------------------------------
c     Get a list of variables on netCDF file
c     ------------------------------------------------------------

      subroutine input_getvars(fid,vnam,nvars)

c     List of variables on netCDF file

      use netcdf

      implicit none

c     Declaration of subroutine parameters
      integer      fid
      integer      nvars
      character*80 vnam(200)

c     Auxiliary variables
      integer ierr
      integer i
      integer nDims,nGlobalAtts,unlimdimid

      ierr = nf90_inquire(fid, nDims, nVars, nGlobalAtts, unlimdimid)
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      
      do i=1,nVars
         ierr = nf90_inquire_variable(fid, i, name = vnam(i))
      enddo
 
      end
