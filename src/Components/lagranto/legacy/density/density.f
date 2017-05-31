      PROGRAM density

      use netcdf

      implicit none

c     ---------------------------------------------------------------------
c     Declaration of variables
c     ---------------------------------------------------------------------
      
c     Parameter and working arrays
      real                                    radius
      character*80                            runit
      integer                                 nx,ny
      integer                                 nlonlat
      real                                    dlonlat
      real                                    xmin,ymin,dx,dy
      real                                    clon,clat
      integer                                 ntime,nfield,ntra
      character*80                            inpfile
      character*80                            outfile
      character*80                            mode
      real                                    param
      integer                                 opts,npts
      integer                                 step
      character*80                            gridtype
      character*80                            field
      integer                                 crefile,crevar
      real,allocatable,    dimension (:,:) :: cnt,res,fld,area
      real,allocatable,    dimension (:)   :: traj
      real,allocatable,    dimension (:)   :: olon,olat,otim,ofld
      real,allocatable,    dimension (:)   :: nlon,nlat,ntim,nfld

c     Output format
      character*80                            outformat

c     Physical and mathematical constants
      real                                    pi180
      parameter                               (pi180=3.14159/180.)
      real                                    deltay
      parameter                               (deltay=111.)
      real                                    eps
      parameter                               (eps=0.001)

c     Input trajectories (see iotra.f)
      integer                                 inpmode
      real,allocatable, dimension (:,:,:) ::  trainp     
      integer                                 reftime(6)      
      character*80                            varsinp(100)   
      integer,allocatable, dimension (:) ::   sel_flag
      character*80                            sel_file
      character*80                            sel_format

c     Auxiliary variables
      character*80                            cdfname,varname
      integer                                 i,j,k
      integer                                 stat
      integer,allocatable, dimension (:,:) :: connect0
      integer                                 connectval0
      integer,allocatable, dimension (:,:) :: connect1
      integer                                 connectval1
      integer,allocatable, dimension (:,:) :: connect2
      integer                                 connectval2
      real                                    slat
      integer                                 ipre
      real                                    addvalue
      real                                    xmax,ymax
      real ,allocatable, dimension (:)  ::    odist,ndist
      real                                    dt
      integer                                 fid
      integer                                 dynamic_grid
      real                                    ycen,xcen
      integer                                 indx,indy
      character*80                            unit
      real                                    pollon,pollat
      real                                    rlon0,rlat0,rlon,rlat
      real                                    lon,lat
      real                                    crot
      integer                                 count
      character*80                            longname, varunit
      real                                    time
      integer                                 ind
      integer                                 ifield
      real                                    hhmm,frac
      integer                                 ierr,ncID

c     External functions
      real         lmstolm,lmtolms
      real         phstoph,phtophs
      external     lmstolm,lmtolms,phstoph,phtophs
      
      real         sdis
      external     sdis

c     ---------------------------------------------------------------------
c     Preparations
c     ---------------------------------------------------------------------
 
c     Write start message
      print*,'========================================================='
      print*,'              *** START OF PROGRAM DENSITY ***'
      print*

c     Read input parameters
      open(10,file='density.param')
       read(10,*) inpfile
       read(10,*) outfile
       read(10,*) field
       read(10,*) ntime,nfield,ntra
       read(10,*) gridtype
       if ( gridtype.eq.'latlon' ) then 
          read(10,*) nx,ny,xmin,ymin,dx,dy
       elseif ( gridtype.eq.'rotated') then
          read(10,*) clon,clat,nlonlat,dlonlat
       else
          print*,' ERROR: unsupported grid type ',trim(gridtype)
          stop
       endif
       read(10,*) radius,runit
       read(10,*) mode
       read(10,*) param
       read(10,*) step
       read(10,*) sel_file
       read(10,*) sel_format
       read(10,*) crefile
       read(10,*) crevar
      close(10)

c     Get the grid parameters if <crefile=0>
      if ( crefile.eq.0 ) then

           ierr = nf90_open  (trim(outfile), NF90_NOWRITE  , ncID)

           ierr = nf90_get_att(ncID, NF90_GLOBAL, 'grid'   ,gridtype ) 
           ierr = nf90_get_att(ncID, NF90_GLOBAL, 'clon'   ,clon     )
           ierr = nf90_get_att(ncID, NF90_GLOBAL, 'clat'   ,clat     )
           ierr = nf90_get_att(ncID, NF90_GLOBAL, 'nlonlat',nlonlat  )
           ierr = nf90_get_att(ncID, NF90_GLOBAL, 'dlonlat',dlonlat  )
           ierr = nf90_get_att(ncID, NF90_GLOBAL, 'nx'     ,nx       )
           ierr = nf90_get_att(ncID, NF90_GLOBAL, 'ny'     ,ny       )
           ierr = nf90_get_att(ncID, NF90_GLOBAL, 'dx'     ,dx       )
           ierr = nf90_get_att(ncID, NF90_GLOBAL, 'dy'     ,dy       )
           ierr = nf90_get_att(ncID, NF90_GLOBAL, 'xmin'   ,xmin     )
           ierr = nf90_get_att(ncID, NF90_GLOBAL, 'ymin'   ,ymin     )

           ierr = nf90_close(ncID)

           print*,'**** GRID PARAMETERS IMPORTED ',
     >            'FROM NETCDF FILE!!!! ****'
           print*
           
      endif

c     Check for consistency
      if ( (step.ne.0).and.(mode.ne.'keep') ) then
         print*," ERROR: interpolation is only possible for all",
     >                   ' time steps... Stop'
         stop
      endif

c     Set the number of times (just code aesthetics)
      opts=ntime

c     Set grid parameters for rotated grid
      if ( gridtype.eq.'rotated' ) then
         nx   = nlonlat
         ny   = nlonlat
         dx   = dlonlat
         dy   = dlonlat
         xmin = - real(nlonlat-1)/2. * dx
         xmax = + real(nlonlat-1)/2. * dx
         ymin = - real(nlonlat-1)/2. * dy
         ymax = + real(nlonlat-1)/2. * dy
      endif
      
c     Set the flag for dynamic grid adjustment
      if ( (nx.eq.0).or.(ny.eq.0) ) then
         dynamic_grid = 1
      else
         dynamic_grid = 0
      endif

c     Print status information
      print*,'---- INPUT PARAMETERS -----------------------------------'
      print* 
      print*,'Input                : ',trim(inpfile)
      print*,'Output               : ',trim(outfile)
      print*,'Field                : ',trim(field)
      print*,'Trajectory           : ',ntime,nfield,ntra
      print*,'Grid type            : ',trim(gridtype)
      if ( dynamic_grid.eq.1 ) then
         print*,'Grid                 : dynamic (see below)'
      elseif ( gridtype.eq.'latlon' ) then
         print*,'Grid   nlon,nlat     : ',nx,ny
         print*,'       lonmin,latmin : ',xmin,ymin
         print*,'       dlon,dlat     : ',dx,dy
      elseif ( gridtype.eq.'rotated' ) then
         print*,'Grid   clon,clat     : ',clon,clat
         print*,'       nlonlat       : ',nlonlat
         print*,'       dlonlat       : ',dlonlat
      endif
      print*,'Filter radius        : ',radius,' ',trim(runit)
      print*,'Mode                 : ',trim(mode)
      if ( ( mode.eq.'time'  ).or.
     >     ( mode.eq.'space' ).or.
     >     (mode.eq.'grid' ) ) 
     >then
         print*,'Parameter            : ',param
      endif
      if ( step.eq.0 ) then
         print*,'Time step            : all'
      elseif (step.gt.0) then
         print*,'Time step            : ',step
      endif
      print*,'Selection file       : ',trim(sel_file)
      print*,'Selection format     : ',trim(sel_file)
      print*,'Flag <crefile>       : ',crefile
      print*,'Flag <crevar>        : ',crevar

c     Check whether mode is valid
      if ((mode.ne.'keep'  ).and.
     >    (mode.ne.'time'  ).and.
     >    (mode.ne.'space' ).and.
     >    (mode.ne.'grid'  ))  
     >then
         print*,' ERROR: Invalid mode ',trim(mode)
         stop
      endif

c     Allocate memory for old and new (reparameterised) trajectory
      allocate(olon(ntime),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array olon ***'
      allocate(olat(ntime),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array olat ***'
      allocate(otim(ntime),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array otim ***'
      allocate(nlon(1000*ntime),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array nlon ***'
      allocate(nlat(1000*ntime),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array nlat ***'
      allocate(ntim(1000*ntime),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array ntim ***'
      allocate(odist(ntime),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array odist ***'
      allocate(ndist(1000*ntime),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array ndist ***'
      allocate(ofld(ntime),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array ofld ***'
      allocate(nfld(1000*ntime),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array nfld ***'

c     Allocate memory for complete trajectory set
      allocate(trainp(ntra,ntime,nfield),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array trainp ***'
      allocate(sel_flag(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array sel_flag ***'

c     Allocate memory for auxiliary fields
      allocate(traj(nfield),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array traj ***'

c     Set the format of the input file
      call mode_tra(inpmode,inpfile)
      if (inpmode.eq.-1) inpmode=1

c     Read the input trajectory file
      call ropen_tra(fid,inpfile,ntra,ntime,nfield,
     >                   reftime,varsinp,inpmode)
      call read_tra (fid,trainp,ntra,ntime,nfield,inpmode)
      call close_tra(fid,inpmode)

c     Bring trajectories into grid interval
      print*,xmin,ymin,dx,dy,nx,ny
      if ( ( abs(xmin+180.).lt.eps).and.
     >     ( abs(ymin+90. ).lt.eps).and.
     >     ( abs(dx-1.    ).lt.eps).and.
     >     ( abs(dy-1.    ).lt.eps).and.
     >     ( nx.eq.360    ).and.
     >     ( ny.eq.180    ) )
     >then
        print*,'HALLO'
        do i=1,ntra
          do j=1,ntime
            if ( trainp(i,j,2).lt.-180. ) then
                trainp(i,j,2) = trainp(i,j,2) + 360.
            endif
            if ( trainp(i,j,2).gt.180. ) then
                trainp(i,j,2) = trainp(i,j,2) - 360.
            endif
         enddo
        enddo
      endif

c     Check that first four columns correspond to time,lon,lat,p
      if ( (varsinp(1).ne.'time' ).or.
     >     (varsinp(2).ne.'xpos' ).and.(varsinp(2).ne.'lon' ).or.
     >     (varsinp(3).ne.'ypos' ).and.(varsinp(3).ne.'lat' ).or.
     >     (varsinp(4).ne.'ppos' ).and.(varsinp(4).ne.'p'   ) )
     >then
         print*,' ERROR: problem with input trajectories ...'
         stop
      endif
      varsinp(1) = 'TIME'
      varsinp(2) = 'lon'
      varsinp(3) = 'lat'
      varsinp(4) = 'p'

c     Get the index of the field (if needed)
      if ( field.ne.'nil' ) then
         ifield = 0
         do i=1,nfield
            if ( varsinp(i).eq.field ) ifield = i
         enddo
         if ( ifield.eq.0 ) then
            print*,' ERROR: field ',trim(field),' not found... Stop'
            stop
         endif
      endif

c     Write some status information of the input trajectories
      print*
      print*,'---- INPUT TRAJECTORIES ---------------------------------'
      print*
      print*,' Reference time (year)  : ',reftime(1)
      print*,'                (month) : ',reftime(2)
      print*,'                (day)   : ',reftime(3)
      print*,'                (hour)  : ',reftime(4)
      print*,'                (min)   : ',reftime(5)
      print*,' Time range (min)       : ',reftime(6)
      do i=1,nfield
         if ( i.ne.ifield ) then
            print*,' Var                    :',i,trim(varsinp(i))
         else
            print*,' Var                    :',i,trim(varsinp(i)),
     >                                        '       [ gridding ]'
         endif
      enddo
      print*,' List of selected times'
      do i=1,ntime
         if ( (step.eq.0).or.(step.eq.i) ) then
            print*,'     ',i,'  -> ',trainp(1,i,1)
         endif
      enddo
      print*

c     Select flag: all trajectories are selected
      if ( sel_file.eq.'nil' ) then

         do i=1,ntra
            sel_flag(i) = 1
         enddo

c     Select flag: index file
      elseif ( sel_format.eq.'index' ) then

         do i=1,ntra
            sel_flag(i) = 0
         enddo
         
         open(10,file=sel_file)
 142      read(10,*,end=141) ind
          sel_flag(ind) = 1
          goto 142 
 141     continue
         close(10)

c     Select flag: boolean file
      elseif ( sel_format.eq.'boolean' ) then
       
         open(10,file=sel_file)
          do i=1,ntra
            read(10,*) ind
            if ( ind.eq.1 ) sel_flag(i) = ind
          enddo
         close(10)
         
      endif

c     Write status information
      if ( sel_file.eq.'nil' ) then
          print*,' Selected trajectories  : all ',ntra          
       else
          count = 0
          do i=1,ntra
             if ( sel_flag(i).eq.1 ) count = count + 1
          enddo
          print*,' #selected trajectories : ',count,
     >            ' [ ',real(count)/real(ntra) * 100.,' % ] '
       endif
       print*

c     ---------------------------------------------------------------------
c     Coordinate transformations and grid adjustment
c     ---------------------------------------------------------------------

c     Transform from lat/lon to rotated lat/lon, if requested
      if ( gridtype.eq.'rotated') then

         crot = 0.

         pollon=clon-180.
         if (pollon.lt.-180.) pollon=pollon+360.
         pollat=90.-clat
         do i=1,ntra
            do j=1,ntime
               
               if ( sel_flag(i).eq.1 ) then

c                Get lat/lon coordinates for trajectory point
                 lon = trainp(i,j,2)
                 lat = trainp(i,j,3)

c                First Rotation
                 pollon=clon-180.
                 if (pollon.lt.-180.) pollon=pollon+360.
                 pollat=90.-clat
                 rlon0=lmtolms(lat,lon,pollat,pollon)
                 rlat0=phtophs(lat,lon,pollat,pollon)            

c                Second rotation
                 pollon=-180.
                 pollat=90.+crot
                 rlon=90.+lmtolms(rlat0,rlon0-90.,pollat,pollon)
                 rlat=phtophs(rlat0,rlon0-90.,pollat,pollon)   

c                Get rotated latitude and longitude
 100             if (rlon.lt.xmin) then
                  rlon=rlon+360.
                  goto 100
                 endif
 102             if (rlon.gt.(xmin+real(nx-1)*dx)) then
                  rlon=rlon-360.
                  goto 102
                 endif

c                Set the new trajectory coordinates
                 trainp(i,j,2) = rlon
                 trainp(i,j,3) = rlat

              endif

            enddo
         enddo
            
      endif

c     Dynamic grid adjustment
      if ( dynamic_grid.eq.1 ) then

c        Get the grid parameters
         xmin =  180.
         ymin =   90.
         xmax = -180.
         ymax =  -90.

         do i=1,ntra

            if ( sel_flag(i).eq.1 ) then

              if ( step.eq.0 ) then
               do j=1,ntime
                  if ( trainp(i,j,2).lt.xmin) xmin =  trainp(i,j,2)
                  if ( trainp(i,j,2).gt.xmax) xmax =  trainp(i,j,2)
                  if ( trainp(i,j,3).lt.ymin) ymin =  trainp(i,j,3)
                  if ( trainp(i,j,3).gt.ymax) ymax =  trainp(i,j,3)
               enddo
              else
                if ( trainp(i,step,2).lt.xmin) xmin =  trainp(i,step,2)
                if ( trainp(i,step,2).gt.xmax) xmax =  trainp(i,step,2)
                if ( trainp(i,step,3).lt.ymin) ymin =  trainp(i,step,3)
                if ( trainp(i,step,3).gt.ymax) ymax =  trainp(i,step,3)
              endif
            
            endif

         enddo

c        Get first guess for "optimal" grid
         nx = 400
         ny = 400
         dx = (xmax - xmin)/real(nx-1)
         dy = (ymax - ymin)/real(ny-1)

c        Make the grid spacing equal in zonal and meridional direction
         if ( dx.gt.dy ) then
            
            dy = dx
            ny = (ymax - ymin)/dy + 1
            if (ny.lt.nx/2)              ny = nx / 2
            if ( real(ny)*dy .ge. 180. ) ny = 180./dy + 1
            ycen = 0.5* (ymin+ymax)
            ymin = ycen - 0.5 * real(ny/2) * dy
            if (ymin.le.-90.) ymin = -90.

         else
               
            dx = dy
            nx = (xmax - xmin)/dx + 1
            if (nx.lt.ny/2)              nx = ny / 2
            if ( real(nx)*dx .ge. 360. ) nx = 360./dx + 1
            xcen = 0.5* (xmin+xmax)
            xmin = xcen - 0.5 * real(nx/2) * dx
            if (xmin.le.-180.) xmin = -180.

         endif
            
c        Write information
         print*
         print*,'---- DYNAMIC GRID ADJUSTMENT',
     >          ' ----------------------------'  
         print*
         print*,'Grid   nlon,nlat     : ',nx,ny
         print*,'       lonmin,latmin : ',xmin,ymin
         print*,'       dlon,dlat     : ',dx,dy
         print*

c     Write grid information for rotated grid (if not already done
      elseif ( gridtype.eq.'rotated') then
         
         print*
         print*,'---- GRID PARAMETERS -------',
     >          ' ----------------------------'  
         print*
         print*,'Grid   nlon,nlat     : ',nx,ny
         print*,'       lonmin,latmin : ',xmin,ymin
         print*,'       dlon,dlat     : ',dx,dy
         print*
       

      endif

c     Set the grid boundaries
      xmax=xmin+real(nx-1)*dx
      ymax=ymin+real(ny-1)*dy

c     Allocate memory for output array and auxiliary gridding array 
      allocate(cnt(nx,ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array cnt  ***'
      allocate(res(nx,ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array res  ***'
      allocate(fld(nx,ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array fld  ***'
      allocate(area(nx,ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array area ***'

      allocate(connect0(nx,ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array connect0 ***'
      allocate(connect1(nx,ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array connect1 ***'
      allocate(connect2(nx,ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array connect2 ***'


c     Init the output array
      do i=1,nx
         do j=1,ny
            connect0(i,j) = 0
            connect1(i,j) = 0
            connect2(i,j) = 0
            cnt(i,j)      = 0.
            res(i,j)      = 0.
            fld(i,j)      = 0.
         enddo
      enddo  

c     ---------------------------------------------------------------------
c     Gridding
c     ---------------------------------------------------------------------

c     Write some status information 
      print*,'---- GRIDDING -------------------------------------------'
      print*

c     Loop over all entries of sampling table
      connectval0 = 0
      connectval1 = 0
      connectval2 = 0
      count       = 0

      do i=1,ntra

         if (mod(i,100).eq.0) print*,i,' of ',ntra

c        Skip all trajectories which are not selected
         if ( sel_flag(i).eq.0 ) goto 300

c        ------- Read a complete trajectory ---------------------------
         do j=1,ntime
            otim(j) = trainp(i,j,1)
            olon(j) = trainp(i,j,2)
            olat(j) = trainp(i,j,3)
            if ( field.ne.'nil' ) then
               ofld(j) =trainp(i,j,ifield)
            endif
         enddo

c        -------- Convert hh.m time into fractional time --------------
         do j=1,ntime
            hhmm    = otim(j)
            call hhmm2frac (hhmm,frac)
            otim(j) = frac
         enddo

c        -------- Interpolation ---------------------------------------

c        Keep the trajectory points as they are
         if ( ( mode.eq.'keep').and.(step.eq.0) ) then
            npts=opts
            do j=1,opts
               ntim(j)=otim(j)
               nlon(j)=olon(j)
               nlat(j)=olat(j)
               if ( field.ne.'nil' ) then
                  nfld(j)=ofld(j)
               endif
            enddo

c        Select a single time step
         elseif ( ( mode.eq.'keep').and.(step.gt.0) ) then
            npts    = 1
            ntim(1) = otim(step)
            nlon(1) = olon(step)
            nlat(1) = olat(step)
            if ( field.ne.'nil' ) then
               nfld(1) = ofld(step)
            endif

c        Perform a reparameterisation in time
         else if ( (mode.eq.'time').and.(step.eq.0) ) then

c           Get the new number of trajectory points
            npts=nint(abs(otim(opts)-otim(1))/param)+1

c           Handle date line problem
            do j=2,opts
               if ( (olon(j-1)-olon(j)).gt.180. ) then
                  olon(j) = olon(j) + 360.
               else if ( (olon(j-1)-olon(j)).lt.-180. ) then
                  olon(j) = olon(j) - 360.
               endif
            enddo
 
c           Cubic spline fitting
            call curvefit(otim,olon,opts,ntim,nlon,npts)
            call curvefit(otim,olat,opts,ntim,nlat,npts)
            if ( field.ne.'nil' ) then
               call curvefit(otim,ofld,opts,ntim,nfld,npts)
            endif

c           Reverse date line handling
            do j=1,npts
               if ( nlon(j).gt.xmax ) then
                  nlon(j) = nlon(j) -360.
               else if ( nlon(j).lt.xmin ) then
                  nlon(j) = nlon(j) +360.
               endif
            enddo

c        Perform a reparameterisation with equally spaced gridpoint
         elseif ( (mode.eq.'space').and.(step.eq.0) ) then
            
c           Calculate the distance and spacing
            odist(1) = 0.
            unit     = 'km'
            do j=2,ntime
               odist(j)=odist(j-1) + 
     >                  sdis(olon(j-1),olat(j-1),olon(j),olat(j),unit)
            enddo
            
c           Determine the new number of trajectory points
            npts=nint(odist(ntime)/param)+1
            if (npts.eq.0) then
               npts=1.
            endif
        
c           Handle date line problem
            do j=2,opts
               if ( (olon(j-1)-olon(j)).gt.180. ) then
                  olon(j) = olon(j) + 360.
               else if ( (olon(j-1)-olon(j)).lt.-180. ) then
                  olon(j) = olon(j) - 360.
               endif
            enddo
                  
c           Cubic spline fitting
            call curvefit(odist,olon,opts,ndist,nlon,npts)
            call curvefit(odist,olat,opts,ndist,nlat,npts)
            call curvefit(odist,otim,opts,ndist,ntim,npts)
            if ( field.ne.'nil' ) then
               call curvefit(odist,ofld,opts,ndist,nfld,npts)
            endif

c           Reverse date line handling
            do j=1,npts
               if ( nlon(j).gt.xmax ) then
                  nlon(j) = nlon(j) -360.
               else if ( nlon(j).lt.xmin ) then
                  nlon(j) = nlon(j) +360.
               endif
            enddo

c        Perform a reparameterisation with equally spaced gridpoint
         elseif ( (mode.eq.'grid').and.(step.eq.0) ) then
            
c           Calculate the distance and spacing
            odist(1) = 0.
            unit     = 'deg'
            do j=2,ntime
               odist(j)=odist(j-1) + 
     >                  sdis(olon(j-1),olat(j-1),olon(j),olat(j),unit)
            enddo
            
c           Determine the new number of trajectory points
            npts=nint(odist(ntime)/param)+1
            if (npts.eq.0) then
               npts=1.
            endif
        
c           Handle date line problem
            do j=2,opts
               if ( (olon(j-1)-olon(j)).gt.180. ) then
                  olon(j) = olon(j) + 360.
               else if ( (olon(j-1)-olon(j)).lt.-180. ) then
                  olon(j) = olon(j) - 360.
               endif
            enddo
                  
c           Cubic spline fitting
            call curvefit(odist,olon,opts,ndist,nlon,npts)
            call curvefit(odist,olat,opts,ndist,nlat,npts)
            call curvefit(odist,otim,opts,ndist,ntim,npts)
            if ( field.ne.'nil' ) then
               call curvefit(odist,ofld,opts,ndist,nfld,npts)
            endif

c           Reverse date line handling
            do j=1,npts
               if ( nlon(j).gt.xmax ) then
                  nlon(j) = nlon(j) -360.
               else if ( nlon(j).lt.xmin ) then
                  nlon(j) = nlon(j) +360.
               endif
            enddo

         endif

c        -------- Do the gridding -------------------------------------

c        Gridding of trajectory
         do j=1,npts

c           Check whether point is in data domain
	    if ( (nlon(j).gt.xmin).and.(nlon(j).lt.xmax).and.
     >       (nlat(j).gt.ymin).and.(nlat(j).lt.ymax))
     >      then

c              Increase counter for gridded points
               count = count + 1

c              ----------------- Gridding: simple count -----------------
               connectval0 = connectval0+1
               
               addvalue    = 1.
               
               call  gridding1
     >              (nlat(j),nlon(j),addvalue,
     >               radius,runit,connect0,connectval0,
     >               cnt,nx,ny,xmin,ymin,dx,dy)

c              ----------------- Gridding: residence time ---------------
               connectval1 = connectval1+1
               
               if ( ntime.eq.1 ) then
                  addvalue = 0.
               elseif ( j.eq.1 )  then
                  addvalue=abs(ntim(2)-ntim(1))
               else
                  addvalue=abs(ntim(j)-ntim(j-1))
               endif
               
               call  gridding1
     >              (nlat(j),nlon(j),addvalue,
     >               radius,runit,connect1,connectval1,
     >               res,nx,ny,xmin,ymin,dx,dy)


c              --------------- Gridding: field -------------------------
               if ( field.ne.'nil' ) then

                   connectval2 = connectval2+1
               
                   addvalue    = nfld(j)
               
                   call  gridding1
     >                  (nlat(j),nlon(j),addvalue,
     >                  radius,runit,connect2,connectval2,
     >                  fld,nx,ny,xmin,ymin,dx,dy)

               endif

	    endif

         enddo

c        Exit point for loop over all trajectories
 300     continue

      enddo

c     Write status information
      print*
      print*,' # gridded points       : ',count

c     ---------------------------------------------------------------------
c     Unit conversions and output to netCDF file
c     ---------------------------------------------------------------------

c     Write some status information 
      print*
      print*,'---- WRITE OUTPUT ---------------------------------------'
      print*

c     Area (in km^2)
      do i=1,nx	         
         do j=1,ny	
            slat=ymin+real(j-1)*dy
            if (abs(abs(slat)-90.).gt.eps) then
               area(i,j) = dy*dx*cos(pi180*slat)*deltay**2
            else
               area(i,j) = 0.
            endif
         enddo
      enddo

c     Normalise gridded field
      if ( field.ne.'nil' ) then
         do i=1,nx
            do j=1,ny
               if ( cnt(i,j).gt.0. ) then
                  fld(i,j) = fld(i,j) / cnt(i,j)
               endif
            enddo
         enddo
      endif

c     Set the time for the output netCDF files - if a composite is
c     calculatd, then the time is set to 
      if ( step.eq.0 ) then
         time = -999.
         print*,'   ... COMPOSITE OVER ALL TRAJECTORY TIMES (-999)'
         print*
      else
         time = trainp(1,step,1)
      endif

c     Write output to CF netCDF
      cdfname  = outfile
      
      varname  = 'COUNT'
      longname = 'trajectory counts'
      varunit  = 'counts per grid point'
      call  writecdf2D_cf (cdfname,varname,longname,varunit,gridtype,
     >       clon,clat,nlonlat,dlonlat,cnt,time,dx,dy,xmin,ymin,nx,
     >       ny,crefile,crefile,1)
      write(*,'(a8,a10,a5,a10,a10,f7.2,a2)') 
     >     '    ... ',trim(varname),' -> ',trim(cdfname),
     >     ' [ time = ',time,' ]'    

      varname  = 'RESIDENCE'
      longname = 'residence time'
      varunit  = 'hours per grid point'

      print*,'crefile = ',crefile

      call  writecdf2D_cf (cdfname,varname,longname,varunit,gridtype,
     >       clon,clat,nlonlat,dlonlat,res,time,dx,dy,xmin,ymin,nx,
     >       ny,0,crefile,1)
      write(*,'(a8,a10,a5,a10,a10,f7.2,a2)') 
     >     '    ... ',trim(varname),' -> ',trim(cdfname),
     >     ' [ time = ',time,' ]'    

      varname  = 'AREA'
      longname = 'area corresponding to grid points'
      varunit  = 'square kilometers'
      call  writecdf2D_cf (cdfname,varname,longname,varunit,gridtype,
     >       clon,clat,nlonlat,dlonlat,area,time,dx,dy,xmin,ymin,nx,
     >       ny,0,crefile,1)
      write(*,'(a8,a10,a5,a10,a10,f7.2,a2)') 
     >     '    ... ',trim(varname),' -> ',trim(cdfname),
     >     ' [ time = ',time,' ]'    
      
      if ( field.ne.'nil' ) then
         varname  = field
         longname = field
         varunit  = 'as on trajectory file'
         call  writecdf2D_cf (cdfname,varname,longname,varunit,gridtype,
     >       clon,clat,nlonlat,dlonlat,fld,time,dx,dy,xmin,ymin,nx,
     >       ny,0,crevar,1)
 
         write(*,'(a8,a10,a5,a10,a10,f7.2,a2)') 
     >        '    ... ',trim(varname),' -> ',trim(cdfname),
     >        ' [ time = ',time,' ]'    
      endif

c     Write status information
      print*
      print*,'              *** END OF PROGRAM DENSITY **'
      print*,'========================================================='

      end

c     ********************************************************************
c     * GRIDDING SUBROUTINES                                             *
c     ********************************************************************

c     ---------------------------------------------------------------------
c     Gridding of one single data point (smoothing in km, deg, gridp)
c     ---------------------------------------------------------------------

      subroutine gridding1 (lat,lon,addval,radius,unit,
     >                      connect,connectval,
     >                      out,nx,ny,xmin,ymin,dx,dy)

      implicit none

c     Declaration of subroutine parameters
      real         lat,lon
      integer      nx,ny
      real         xmin,ymin,dx,dy
      real         out(nx,ny)
      real         radius
      character*80 unit
      integer      connectval
      integer      connect(nx,ny)
      real         addval

c     Auxiliary variables
      integer   i,j,k
      integer   mu,md,nr,nl,n,m
      integer   stackx(nx*ny),stacky(nx*ny)
      integer   tab_x(nx*ny),tab_y(nx*ny)
      real      tab_r(nx*ny)
      integer   sp
      real      lat2,lon2
      real      dist,sum
      real      xmax
      integer   periodic
      integer   test

c     Numerical epsilon
      real      eps
      parameter (eps=0.01)

c     Externals
      real      sdis,weight
      external  sdis,weight

c     Check whether lat/lon point is valid
      xmax=xmin+real(nx-1)*dx
      if (lon.lt.xmin-eps) lon=lon+360.
      if (lon.gt.xmax+eps) lon=lon-360.
      if (abs(lat-90).lt.eps) lat=90.
      if (abs(lat+90).lt.eps) lat=-90.
      if ((abs(lat).gt.(90.+eps)).or.
     >    (lon.lt.xmin-eps).or.(lon.gt.xmax+eps)) then
         print*,'Invalid lat/lon point ',lat,lon
         return
      endif

c     Set flag for periodic domain
      if (abs(xmax-xmin-360.).lt.eps) then
         periodic=1
      else if (abs(xmax-xmin-360+dx).lt.eps) then
         periodic=2
      else
         periodic=0
      endif

c     Get indices of one coarse grid point within search radius
      i=nint((lon-xmin)/dx)+1
      if ((i.eq.nx).and.(periodic.eq.1)) i=1
      j=nint((lat-ymin)/dy)+1
      lat2=ymin+real(j-1)*dy
      lon2=xmin+real(i-1)*dx
      dist=sdis(lon,lat,lon2,lat2,unit)
      if (dist.gt.radius) then
         print*,'1: Search radius is too small...'
         stop
      endif

c     Get connected points
      k=0
      stackx(1)=i
      stacky(1)=j
      sp=1
      do while (sp.ne.0) 
         
c        Get an element from stack
         n=stackx(sp)
         m=stacky(sp)
         sp=sp-1
                  
c        Get distance from reference point
         lat2=ymin+real(m-1)*dy
         lon2=xmin+real(n-1)*dx
         dist=sdis(lon,lat,lon2,lat2,unit)

c        Check whether distance is smaller than search radius: connected
         if (dist.lt.radius) then

c           Make entry in filter mask
            k=k+1
            tab_x(k)=n
            tab_y(k)=m
            tab_r(k)=weight(dist,radius)

c           Mark this point as visited
            connect(n,m)=connectval
                     
c           Get coordinates of neighbouring points
            nr=n+1
            if ((nr.gt.nx)  .and.(periodic.eq.0)) nr=nx
            if ((nr.gt.nx-1).and.(periodic.eq.1)) nr=1
            if ((nr.gt.nx)  .and.(periodic.eq.2)) nr=1
            nl=n-1
            if ((nl.lt.1).and.(periodic.eq.0)) nl=1
            if ((nl.lt.1).and.(periodic.eq.1)) nl=nx-1
            if ((nl.lt.1).and.(periodic.eq.2)) nl=nx
            mu=m+1
            if (mu.gt.ny) mu=ny
            md=m-1
            if (md.lt.1) md=1

c           Update stack
            if (connect(nr,m).ne.connectval) then
               connect(nr,m)=connectval
               sp=sp+1
               stackx(sp)=nr
               stacky(sp)=m
            endif
            if (connect(nl,m).ne.connectval) then
               connect(nl,m)=connectval
               sp=sp+1
               stackx(sp)=nl
               stacky(sp)=m
            endif
            if (connect(n,mu).ne.connectval) then
               connect(n,mu)=connectval
               sp=sp+1
               stackx(sp)=n
               stacky(sp)=mu
            endif
            if (connect(n,md).ne.connectval) then
               connect(n,md)=connectval
               sp=sp+1
               stackx(sp)=n
               stacky(sp)=md
            endif
         endif
         
      end do

      if (k.ge.1) then
         sum=0.
         do i=1,k
            sum=sum+tab_r(i)
         enddo
         do i=1,k
            out(tab_x(i),tab_y(i))=out(tab_x(i),tab_y(i))+
     >                             addval*tab_r(i)/sum

            if ((tab_x(i).eq.1).and.(periodic.eq.1)) then
               out(nx,tab_y(i))=out(nx,tab_y(i))+
     >                             addval*tab_r(i)/sum
            endif
         enddo
      else
         print*,'2: Search radius is too small...'
         stop
      endif

      end


c     ----------------------------------------------------------------------
c     Get spherical distance between lat/lon points
c     ----------------------------------------------------------------------
            
      real function sdis(xp,yp,xq,yq,unit)

c     Calculates spherical distance (in km) between two points given
c     by their spherical coordinates (xp,yp) and (xq,yq), respectively.

      real         re
      parameter    (re=6370.)
      real         xp,yp,xq,yq,arg
      character*80 unit
      real         dlon

      if ( unit.eq.'km' ) then

         arg=sind(yp)*sind(yq)+cosd(yp)*cosd(yq)*cosd(xp-xq)
         if (arg.lt.-1.) arg=-1.
         if (arg.gt.1.) arg=1.
         sdis=re*acos(arg)
         
      elseif ( unit.eq.'deg' ) then

         dlon = xp-xq
         if ( dlon.gt. 180. ) dlon = dlon - 360.
         if ( dlon.lt.-180. ) dlon = dlon + 360.
         sdis = sqrt( dlon**2 + (yp-yq)**2 )

      endif
      

c     Quick and dirty trick to avoid zero distances
      if (sdis.eq.0.) sdis=0.1

      end

c     ----------------------------------------------------------------------
c     Weight function for the filter mask
c     ----------------------------------------------------------------------
 
      real function weight (r,radius)

c     Attribute to each distanc r its corresponding weight in the filter mask

      implicit none

c     Declaration of subroutine parameters
      real r
      real radius

c     Simple 0/1 mask
      if (r.lt.radius) then
         weight=exp(-r/radius)
      else
         weight=0.
      endif

      end


c     ********************************************************************
c     * REPARAMETERIZATION SUBROUTINES                                   *
c     ********************************************************************

c     -------------------------------------------------------------
c     Interpolation of the trajectory with a natural cubic spline
c     -------------------------------------------------------------

      SUBROUTINE curvefit (time,lon,n,
     >                     sptime,splon,spn)

c     Given the curve <time,lon> with <n> data points, fit a
c     cubic spline to this curve. The new curve is returned in 
c     <sptime,splon,spn> with <spn> data points. The parameter
c     <spn> specifies on entry the number of spline interpolated points
c     along the curve.
      
      implicit none

c     Declaration of subroutine parameters
      integer n
      real time(n),lon(n)
      integer spn
      real sptime(spn),splon(spn)

c     Auxiliary variables
      real y2ax(n)
      real dt
      real s
      integer i
      real order

c     Determine whether the input array is ascending or descending
      if (time(1).gt.time(n)) then
         order=-1.
      else
         order= 1.
      endif

c     Bring the time array into ascending order
      do i=1,n
         time(i)=order*time(i)
      enddo

c     Prepare the (natural) cubic spline interpolation
      call spline (time,lon,n,1.e30,1.e30,y2ax)
      dt=(time(n)-time(1))/real(spn-1)
      do i=1,spn
         sptime(i)=time(1)+real(i-1)*dt
      enddo
      
c     Do the spline interpolation
      do i=1,spn
         call splint(time,lon,y2ax,n,sptime(i),s)
         splon(i)=s
      enddo

c     Change the time arrays back
      do i=1,spn
         sptime(i)=order*sptime(i)
      enddo
      do i=1,n
         time(i)=order*time(i)
      enddo

      return
      end

c     -------------------------------------------------------------
c     Basic routines for spline interpolation (Numerical Recipes)
c     -------------------------------------------------------------

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
         print*,'bad xa input in splint'
         stop
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      END

c     ********************************************************************
c     * INPUT / OUTPUT SUBROUTINES                                       *
c     ********************************************************************


c     --------------------------------------------------------------------
c     Subroutines to write the CF netcdf output file
c     --------------------------------------------------------------------

      subroutine writecdf2D_cf 
     >          (cdfname,varname,longname,unit,gridtype,clon,clat,
     >           nlonlat,dlonlat,arr,time,dx,dy,xmin,ymin,nx,ny,
     >           crefile,crevar,cretime)

c     Create and write to the CF netcdf file <cdfname>. The variable
c     with name <varname> and with time <time> is written. The data
c     are in the two-dimensional array <arr>. The list <dx,dy,xmin,
c     ymin,nx,ny> specifies the output grid. The flags <crefile> and
c     <crevar> determine whether the file and/or the variable should
c     be created; correspondingly for the unlimited dimension <time>
c     with the flag <cretime>.

      USE netcdf

      IMPLICIT NONE

c     Declaration of input parameters
      character*80 cdfname
      character*80 varname,longname,unit
      integer      nx,ny
      real         arr(nx,ny)
      real         dx,dy,xmin,ymin
      real         time
      integer      crefile,crevar,cretime
      character*80 gridtype
      real         clon,clat
      integer      nlonlat
      real         dlonlat

c     Local variables
      integer      ierr
      integer      ncID
      integer      LonDimId,    varLonID
      integer      LatDimID,    varLatID
      integer      TimeDimID,   varTimeID
      real         longitude(nx)
      real         latitude (ny)
      real         timeindex
      integer      i
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
      ierr=nf90_def_dim(ncID,'longitude',nx            , LonDimID )
      ierr=nf90_def_dim(ncID,'latitude' ,ny            , LatDimID )
      ierr=nf90_def_dim(ncID,'time'     ,nf90_unlimited, TimeDimID)
      
c     Define coordinate Variables 
      ierr = nf90_def_var(ncID,'longitude',NF90_FLOAT,
     >     (/ LonDimID /),varLonID)
      ierr = nf90_put_att(ncID, varLonID, "standard_name","longitude")
      ierr = nf90_put_att(ncID, varLonID, "units"      ,"degree_east")
      
      ierr = nf90_def_var(ncID,'latitude',NF90_FLOAT,
     >     (/ LatDimID /),varLatID)
      ierr = nf90_put_att(ncID, varLatID, "standard_name", "latitude")
      ierr = nf90_put_att(ncID, varLatID, "units"    ,"degree_north")
      
      ierr = nf90_def_var(ncID,'time',NF90_FLOAT, 
     >     (/ TimeDimID /), varTimeID)
      ierr = nf90_put_att(ncID, varTimeID, "axis",            "T")
      ierr = nf90_put_att(ncID, varTimeID, "calendar", "standard")
      ierr = nf90_put_att(ncID, varTimeID, "long_name",    "time")
      ierr = nf90_put_att(ncID, varTimeID, "units",       "hours")
      
c     Write global attributes 
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'Conventions', 'CF-1.0')
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'title',  
     >     'Trajectory Densities')
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'source', 
     >     'Lagranto Trajectories')
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'institution', 
     >     'ETH Zurich, IACETH')
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'grid',trim(gridtype) ) 
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'clon',clon )
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'clat',clat )
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'nlonlat',nlonlat )
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'dlonlat',dlonlat )
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'nx',nx )
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'ny',ny )
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'dx',dx )
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'dy',dy )
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'xmin',xmin )
      ierr = nf90_put_att(ncID, NF90_GLOBAL, 'ymin',ymin )

c     Write coordinate data 
      do i = 1,nx+1
         longitude(i) = xmin + real(i-1) * dx
      enddo
      do i = 1,ny+1
         latitude(i)  = ymin + real(i-1) * dy
      enddo
      
c     Check whether the definition was successful
      ierr = nf90_enddef(ncID)
      if (ierr.gt.0) then
         print*, 'An error occurred while attempting to ', 
     >        'finish definition mode.'
         stop
      endif
      
c     Write coordinate data  
      ierr = nf90_put_var(ncID,varLonID ,longitude)
      ierr = nf90_put_var(ncID,varLatID ,latitude )
      
c     Close netCDF file 
      ierr = nf90_close(ncID)
      
 100  continue

c     ---- Define a new variable - skip if <crevar=0> -----------------------

      if ( crevar.ne.1 ) goto 110
      
c     Open the file for read(write access
      ierr = nf90_open  (trim(cdfname), NF90_WRITE  , ncID)

c     Get the IDs for dimensions
      ierr = nf90_inq_dimid(ncID,'longitude', LonDimID )
      ierr = nf90_inq_dimid(ncID,'latitude' , LatDimID )
      ierr = nf90_inq_dimid(ncID,'time'     , TimeDimID)

c     Enter define mode
      ierr = nf90_redef(ncID)

c     Write definition and add attributes
      ierr = nf90_def_var(ncID,varname,NF90_FLOAT,
     >                   (/ LonDimID, LatDimID, TimeDimID /),varID)
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
      print*,trim(varname),' defined on ',trim(cdfname)

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
     >                    start = (/ 1, 1, ind /), 
     >                    count = (/ nx, ny, 1 /) )

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
c     * Transformation routine: LMSTOLM and PHSTOPH from library gm2em               *
c     ********************************************************************************

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

c     ------------------------------------------------------------------
c     Compute Cos/Sin of an argument in Degree instead of Radian
c     ------------------------------------------------------------------

      real function cosd(arg)

      real,intent(IN) :: arg
      real,parameter :: grad2rad=3.1415926/180.
      cosd=cos(arg*grad2rad)
      return
      end

      real function sind(arg)

      real,intent(IN) :: arg
      real,parameter :: grad2rad=3.1415926/180.
      sind=sin(arg*grad2rad)
      return
      end


