      PROGRAM lsl2rdf
      
c     ***********************************************************************
c     * Convert a trajectory file into a lat/lon/p RDF-netCDF file          *
c     * RDF = Reverse Domain Filling                                        *
c     * Michael Sprenger / Spring, summer 2010                              *
c     ***********************************************************************

      use netcdf

      implicit none
      
c     ----------------------------------------------------------------------
c     Declaration of parameters
c     ----------------------------------------------------------------------

c     Real comparison
      real       eps
      parameter (eps=0.001)

c     ----------------------------------------------------------------------
c     Declaration of variables
c     ----------------------------------------------------------------------

c     Input and output format for trajectories (see iotra.f)
      character*80                           inpfile     ! Input filename
      character*80                           outfile     ! Output filename

c     Trajectories
      integer                                ntra        ! Number of trajectories
      integer                                ntim        ! Number of times
      integer                                ncol        ! Number of columns
      real,allocatable, dimension (:,:,:) :: tra         ! Trajectories (ntra,ntim,ncol)
      character*80                           vars(100)   ! Variable names
      integer                                refdate(6)  ! Reference date

c     Output RDF-netCDF file
      character*80                           cdfname
      character*80                           varname
      character*80                           longname
      character*80                           unit
      integer                                nx,ny,nz
      real,allocatable, dimension (:) ::     lon,lat,lev
      real,allocatable, dimension (:,:,:) :: arr
      real                                   time
      integer                                crefile,crevar,cretime
      integer,allocatable, dimension (:,:,:) :: ind
      integer                                indx,indy,indz

c     Auxiliary variables
      integer                                inpmode
      integer                                stat
      integer                                fid
      integer                                i,j,k,ivar,itim
      real                                   tmp(10000)
      integer                                ok
      integer                                count
      
c     ----------------------------------------------------------------------
c     Read input data
c     ----------------------------------------------------------------------

c     Write start message
      print*,'========================================================='
      print*,'              *** START OF PROGRAM LSL2RDF ***'
      print*

c     Read parameters
      open(10,file='lsl2rdf.param')
       read(10,*) inpfile
       read(10,*) outfile
       read(10,*) ntra,ntim,ncol
      close(10)
      
c     Determine the formats
      call mode_tra(inpmode,inpfile)
      if (inpmode.eq.-1) inpmode=1

c     Allocate memory for input trajectory
      allocate(tra(ntra,ntim,ncol),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra      ***' 

c     Read inpufile
      call ropen_tra(fid,inpfile,ntra,ntim,ncol,refdate,vars,inpmode)
      call read_tra (fid,tra,ntra,ntim,ncol,inpmode)
      call close_tra(fid,inpmode)

c     ----------------------------------------------------------------------
c     Determine output grid and get mapping {output array} <-> {trajectory}
c     ----------------------------------------------------------------------

c     Get a list of longitudes at time 0
      nx = 0
      do i=1,ntra
        ok = 1
        do j=1,nx
           if ( abs(tmp(j)-tra(i,1,2)).lt.eps ) then
              ok = 0
              goto 100
           endif
        enddo
 100    if ( ok.eq.1 ) then
          nx      = nx + 1
          tmp(nx) = tra(i,1,2)
        endif
      enddo

      if ( nx.gt.1 ) then
        call sort(nx,tmp)
      endif

      allocate(lon(nx),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array lon      ***'
      do i=1,nx
        lon(i) = tmp(i)
      enddo

c     Get a list of latitudes at time 0
      ny = 0
      do i=1,ntra
        ok = 1
        do j=1,ny
           if ( abs(tmp(j)-tra(i,1,3)).lt.eps ) then
              ok = 0
           endif
        enddo
 101    if ( ok.eq.1 ) then
          ny      = ny + 1
          tmp(ny) = tra(i,1,3)
        endif
      enddo

      if ( ny.gt.1 ) then
        call sort(ny,tmp)
      endif

      allocate(lat(ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array lat      ***'
      do i=1,ny
        lat(i) = tmp(i)
      enddo

c     Get a list of pressure levels at time 0
      nz = 0
      do i=1,ntra
        ok = 1
        do j=1,nz
           if ( abs(tmp(j)-tra(i,1,4)).lt.eps ) then
              ok = 0
              goto 102
           endif
        enddo
 102    if ( ok.eq.1 ) then
          nz      = nz + 1
          tmp(nz) = tra(i,1,4)
        endif
      enddo

      if ( nz.gt.1 ) then
         call sort(nz,tmp)
      endif

      allocate(lev(nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array lev      ***'
      do i=1,nz
        lev(i) = tmp(i)
      enddo

c     Write output grid
      print*,'---- OUTPUT GRID ----------------------------------------'
      print*
      write(*,'(a10,1000f5.0)') 'lon:',(lon(i),i=1,nx)
      write(*,'(a10,1000f5.0)') 'lat:',(lat(i),i=1,ny)
      write(*,'(a10,1000f5.0)') 'lev:',(lev(i),i=1,nz)
      print*

c     Decide whether several levels are really needed
      if ( (nx*ny).eq.ntra ) then
        if ( nz.ne.1 ) then
            print*
            print*,'All levels will be mapped to index level 1'
            nz = 1
        endif
      endif

c     Allocate memory for output RDF-netCDF
      allocate(arr(nx,ny,nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array arr      ***'
      allocate(ind(nx,ny,nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array ind      ***'

c     Init the mapping array <-> trajectory
      do i=1,nx
        do j=1,ny
            do k=1,nz
                ind(i,j,k) = 0
            enddo
        enddo
      enddo

c     Now get for each grid point the right trajectory index
      do i=1,ntra

         indx = 0
         do j=1,nx
            if ( abs(lon(j)-tra(i,1,2)).lt.eps  ) then
              indx = j
              goto 200
            endif
         enddo
 200     continue

         indy = 0
         do j=1,ny
            if ( abs(lat(j)-tra(i,1,3)).lt.eps  ) then
              indy = j
              goto 201
            endif
         enddo
 201     continue

         if ( nz.gt.1 ) then
           indz = 0
           do j=1,nz
            if ( abs(lev(j)-tra(i,1,4)).lt.eps  ) then
              indz = j
              goto 202
            endif
           enddo
 202       continue
         else
           indz = 1
         endif

         if ( (indx.ne.0).and.(indy.ne.0).and.(indz.ne.0) ) then
            ind(indx,indy,indz) = i
         endif

      enddo

c     Check whether all grid points are linked to a traejctory
      count = 0
      do i=1,nx
        do j=1,ny
            do k=1,nz
                if ( ind(i,j,k).eq.0 ) count = count + 1
            enddo
        enddo
      enddo
      if ( count.ne.0 ) then
        print*,' WARNING: not all grid points in RDF domain linked'
        print*,'          to a trajectory... proceed anyway.'
      endif


c     ----------------------------------------------------------------------
c     Fill the output array
c     ----------------------------------------------------------------------

      print*,'---- MAPPING TO RDF GRID --------------------------------'
      print*

c     Loop over all variables
      crefile = 1
      do ivar=2,ncol

c        Loop over all times
         cretime = 1
         crevar  = 1
         do itim=1,ntim

c            Get the time from the first trajectory
             time = tra(1,itim,1)
             write(*,'(a20,a10,f8.2)')
     >            ' variable,time : ',trim( vars(ivar) ),time

c            Do the remapping
             do i=1,nx
                do j=1,ny
                    do k=1,nz
                       if ( ind(i,j,k).ne.0 ) then
                        arr(i,j,k) = tra( ind(i,j,k), itim, ivar )
                       endif
                    enddo
                enddo
             enddo

c            Save the array to RDF-netCDF
             call writecdf2D_cf
     >          (outfile,vars(ivar),nx,lon,ny,lat,nz,lev,
     >           arr,time,crefile,crevar,cretime)

             crefile = 0
             crevar  = 0

         enddo

      enddo

c     Write end message
      print*
      print*,'              *** END OF PROGRAM LSL2RDF ***'
      print*,'========================================================='

      end

c     ********************************************************************
c     * SUBROUTINES                                                      *
c     ********************************************************************

c     --------------------------------------------------------------------
c     Subroutines to write the CF netcdf output file
c     --------------------------------------------------------------------

      subroutine writecdf2D_cf
     >          (cdfname,varname,nx,lon,ny,lat,nz,lev,
     <           arr,time,crefile,crevar,cretime)

      USE netcdf

      IMPLICIT NONE

c     Declaration of input parameters
      character*80 cdfname
      character*80 varname,longname,unit
      integer      nx,ny,nz
      real         lon(nx),lat(ny),lev(nz)
      real         arr(nx,ny,nz)
      real         time
      integer      crefile,crevar,cretime

c     Local variables
      integer      ierr
      integer      ncID
      integer      LonDimId,    varLonID
      integer      LatDimID,    varLatID
      integer      LevDimID,    varLevID
      integer      TimeDimID,   varTimeID
      real         timeindex
      integer      i
      integer      nvars,varids(100)
      integer      ndims,dimids(100)
      real         timelist(1000)
      integer      ntimes
      integer      ind
      integer      varID
      real         xmin,xmax,ymin,ymax,zmin,zmax

c     Quick an dirty solution for fieldname conflict
      if ( varname.eq.'time' ) varname = 'TIME'

c     Initially set error to indicate no errors.
      ierr = 0

c     ---- Create the netCDF - skip if <crefile=0> ----------------------
      if ( crefile.ne.1 ) goto 100

c     Create the file
      ierr = nf90_create(trim(cdfname), NF90_CLOBBER, ncID)

c     Define dimensions
      ierr=nf90_def_dim(ncID,'longitude' ,nx            , LonDimID )
      ierr=nf90_def_dim(ncID,'latitude'  ,ny            , LatDimID )
      ierr=nf90_def_dim(ncID,'level    ' ,nz            , LevDimID )
      ierr=nf90_def_dim(ncID,'time'      ,nf90_unlimited, TimeDimID)

c     Define coordinate Variables
      ierr = nf90_def_var(ncID,'longitude',NF90_FLOAT,
     >     (/ LonDimID /),varLonID)
      ierr = nf90_put_att(ncID, varLonID, "standard_name","longitude")
      ierr = nf90_put_att(ncID, varLonID, "units"      ,"degree_east")

      ierr = nf90_def_var(ncID,'latitude',NF90_FLOAT,
     >     (/ LatDimID /),varLatID)
      ierr = nf90_put_att(ncID, varLatID, "standard_name", "latitude")
      ierr = nf90_put_att(ncID, varLatID, "units"    ,"degree_north")

      ierr = nf90_def_var(ncID,'level',NF90_FLOAT,
     >     (/ LevDimID /),varLevID)
      ierr = nf90_put_att(ncID, varLevID, "standard_name", "level")
      ierr = nf90_put_att(ncID, varLevID, "units"    ,"nil")

      ierr = nf90_def_var(ncID,'time',NF90_FLOAT,
     >     (/ TimeDimID /), varTimeID)

c     Write general global attributes
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'Conventions', 'CF-1.0')
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'title',
     >     'RDF trajectory')
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'source',
     >     'Lagranto Trajectories')
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'institution',
     >     'ETH Zurich, IACETH')

c     Write grid extension
      xmin = lon(1)
      xmax = lon(1)
      do i=2,nx
        if ( lon(i).lt.xmin ) xmin = lon(i)
        if ( lon(i).gt.xmax ) xmax = lon(i)
      enddo

      ymin = lat(1)
      ymax = lat(1)
      do i=2,ny
        if ( lat(i).lt.ymin ) ymin = lat(i)
        if ( lat(i).gt.ymax ) ymax = lat(i)
      enddo

      zmin = lev(1)
      zmax = lev(1)
      do i=2,nz
        if ( lev(i).lt.zmin ) zmin = lev(i)
        if ( lev(i).gt.zmax ) zmax = lev(i)
      enddo

      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'nx',nx )
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'ny',ny )
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'nz',nz )
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'xmin',xmin )
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'xmax',xmax )
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'ymin',ymin )
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'ymax',ymax )
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'zmin',zmin )
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'zmax',zmax )

c     Check whether the definition was successful
      ierr = nf90_enddef(ncID)
      if (ierr.gt.0) then
         print*, 'An error occurred while attempting to ',
     >        'finish definition mode.'
         stop
      endif

c     Write coordinate data
      ierr = nf90_put_var(ncID,varLonID ,lon )
      ierr = nf90_put_var(ncID,varLatID ,lat )
      ierr = nf90_put_var(ncID,varLevID ,lev )

c     Close netCDF file
      ierr = nf90_close(ncID)

 100  continue

c     ---- Define a new variable - skip if <crevar=0> -----------------------

      if ( crevar.ne.1 ) goto 110

c     Open the file for read(write access
      ierr = nf90_open  (trim(cdfname), NF90_WRITE  , ncID)

c     Get the IDs for dimensions
      ierr = nf90_inq_dimid(ncID,'longitude'     , LonDimID )
      ierr = nf90_inq_dimid(ncID,'latitude'      , LatDimID )
      ierr = nf90_inq_dimid(ncID,'level'         , LevDimID )
      ierr = nf90_inq_dimid(ncID,'time'          , TimeDimID)

c     Enter define mode
      ierr = nf90_redef(ncID)

c     Write definition and add attributes
      ierr = nf90_def_var(ncID,varname,NF90_FLOAT,
     >             (/ LonDimID, LatDimID, LevDimID, TimeDimID /),varID)
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

c     ---- Create a new time (unlimited dimension) - skip if <cretime=0> ------

      if ( cretime.ne.1 ) goto 120

c     Open the file for read/write access
      ierr = nf90_open  (trim(cdfname), NF90_WRITE, ncID)

c     Get the list of times on the netCDF file
      ierr = nf90_inq_dimid(ncID,'time', TimeDimID)
      if ( ierr.ne.0 ) then
         print*,'Time dimension is not defined on ',trim(cdfname),
     >          ' .... Stop'
         stop
      endif
      ierr = nf90_inquire_dimension(ncID, TimeDimID, len = ntimes)
      ierr = nf90_inq_varid(ncID,'time', varTimeID)
      if ( ierr.ne.0 ) then
         print*,'Variable time is not defined on ',trim(cdfname),
     >          ' ... Stop'
         stop
      endif
      ierr = nf90_get_var(ncID,varTimeID,timelist(1:ntimes))

c     Decide whether a new time must be written
      ind = 0
      do i=1,ntimes
         if ( time.eq.timelist(i) ) ind = i
      enddo

c     Extend the time list if required
      if ( ind.eq.0 ) then
         ntimes           = ntimes + 1
         timelist(ntimes) = time
         ierr = nf90_put_var(ncID,varTimeID,timelist(1:ntimes))
      endif

c     Close netCDF file
      ierr = nf90_close(ncID)

 120  continue

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

c     Get the time index
      ierr = nf90_inq_dimid(ncID,'time', TimeDimID)
      if ( ierr.ne.0 ) then
         print*,'Time dimension is not defined on ',trim(cdfname),
     >          ' .... Stop'
         stop
      endif
      ierr = nf90_inquire_dimension(ncID, TimeDimID, len = ntimes)
      ierr = nf90_inq_varid(ncID,'time', varTimeID)
      if ( ierr.ne.0 ) then
         print*,'Variable time is not defined on ',trim(cdfname),
     >          ' ... Stop'
         stop
      endif
      ierr = nf90_get_var(ncID,varTimeID,timelist(1:ntimes))
      ind = 0
      do i=1,ntimes
         if ( time.eq.timelist(i) ) ind = i
      enddo
      if (ind.eq.0) then
         print*,'Time',time,' is not defined on the netCDF file',
     >          trim(cdfname),' ... Stop'
         stop
      endif

c     Write data block
      ierr = nf90_put_var(ncID,varID,arr,
     >                    start = (/  1,  1,  1, ind /),
     >                    count = (/ nx, ny, nz,   1 /) )

c     Check whether writing was successful
      ierr = nf90_close(ncID)
      if (ierr.ne.0) then
         write(*,*) trim(nf90_strerror(ierr))
         write(*,*) 'An error occurred while attempting to ',
     >              'close the netcdf file.'
         write(*,*) 'in clscdf_CF'
      endif

      end

c     --------------------------------------------------------------------
c     Sort an array (from Muerical Recipes)
c     --------------------------------------------------------------------

      SUBROUTINE SORT(N,RA)
      DIMENSION RA(N)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END
      

      
