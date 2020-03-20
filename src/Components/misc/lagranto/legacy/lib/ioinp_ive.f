c     ************************************************************
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

c     ------------------------------------------------------------
c     Open input file
c     ------------------------------------------------------------

      subroutine input_open (fid,filename)

c     Open the input file with filename <filename> and return the
c     file identifier <fid> for further reference. 

      implicit none

c     Declaration of subroutine parameters
      integer      fid              ! File identifier
      character*80 filename         ! Filename

c     Declaration of auxiliary variables
      integer      ierr

c     Open IVE netcdf file
      call cdfopn(filename,fid,ierr)
      if (ierr.ne.0) goto 900
      
c     Exception handling
      return
      
 900  print*,'Cannot open input file ',trim(filename)
      stop

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
      real         times(10)
      integer      ntimes
      real         aklay(nz),bklay(nz),aklev(nz),bklev(nz)
      integer      nvars
      character*80 vars(100)

c     Auxiliary varaibles
      integer      ierr       
      integer      i,j,k
      integer      isok
      real         tmp(200)
      character*80 varname
      real         rtime
      integer      is2d
      integer      plev

c     Init the flag for 2D variables - assume a 3D field
      is2d = 0

c     Init the flag for pressure levels (PS is not needed)
      plev = 0

c     Inquire dimensions and grid constants if <fid> is negative
      if (fid.lt.0) then

c        Get grid info for the selected variable
         if ( fieldname.eq.'PLEV' ) then
            varname = 'PS'
            stagz   = 0.
            call getdef(-fid,varname,ndim,mdv,vardim,
     >                  varmin,varmax,stag,ierr)
            if (ierr.ne.0) goto 900
            call getcfn(-fid,cstfile,ierr)
            if (ierr.ne.0) goto 903
            call cdfopn(cstfile,cstid,ierr)
            if (ierr.ne.0) goto 903
            call getlevs(cstid,vardim(3),tmp,tmp,tmp,tmp,ierr)
            if (ierr.ne.0) goto 903
            call clscdf(cstid,ierr)
            if (ierr.ne.0) goto 903
            
         elseif ( ( fieldname.eq.'PLAY' ).or.( fieldname.eq.'P' ) ) then       
            varname = 'PS'
            stagz   = -0.5
            call getdef(-fid,varname,ndim,mdv,vardim,
     >                  varmin,varmax,stag,ierr)
            if (ierr.ne.0) goto 900
            call getcfn(-fid,cstfile,ierr)
            if (ierr.ne.0) goto 903
            call cdfopn(cstfile,cstid,ierr)
            if (ierr.ne.0) goto 903
            call getlevs(cstid,vardim(3),tmp,tmp,tmp,tmp,ierr)
            if (ierr.ne.0) goto 903
            call clscdf(cstid,ierr)
            if (ierr.ne.0) goto 903

         else
            varname = fieldname
            call getdef(-fid,varname,ndim,mdv,vardim,
     >                  varmin,varmax,stag,ierr)
            if (ierr.ne.0) goto 900
            
         endif

c        Set the grid dimensions and constants - vardim(3) is taken from constants file
         nx   = vardim(1)
         ny   = vardim(2)
         nz   = vardim(3)
         xmin = varmin(1)
         ymin = varmin(2)
         xmax = varmax(1)
         ymax = varmax(2)
         dx   = (xmax-xmin)/real(nx-1)
         dy   = (ymax-ymin)/real(ny-1)

c        Get pole position - if no constants file available, set default pole
         call getcfn(-fid,cstfile,ierr)
         if (ierr.eq.0) then  
            call cdfopn(cstfile,cstid,ierr)
            if (ierr.ne.0) goto 903
            call getpole(cstid,pollon,pollat,ierr)
            if (ierr.ne.0) goto 903
            call clscdf(cstid,ierr)
            if (ierr.ne.0) goto 903
         else
            pollon = 0.
            pollat = 90.
         endif

c     Get non-constant grid parameters (surface pressure and vertical grid)
      else
         
c        Special handling for fieldname 'P.ML' -> in this case the fields
c        P and PS are available on the data file and can be read in. There
c        is no need to reconstruct it from PS,AK and BK. This mode is
c        used for model-level (2D) trajectories
         if ( fieldname.eq.'P.ML' ) then

c           Get the right time to read
            call gettimes(fid,times,ntimes,ierr)
            if (ierr.ne.0) goto 901
            isok=0
            do i=1,ntimes
               if (abs(time-times(i)).lt.eps) then
                  isok = 1
                  rtime = times(i)
               elseif (timecheck.eq.'no') then
                  isok = 1
                  rtime = times(1)
               endif
            enddo

c           Read surface pressure and 3D pressure
            varname='PS'
            call getdat(fid,varname,rtime,0,ps,ierr)
            if (ierr.ne.0) goto 904
            varname='P'
            call getdat(fid,varname,rtime,0,p3,ierr)
            if (ierr.ne.0) goto 904

c           Set MDV to 1050. - otherwise interpolation routines don't work
            do i=1,nx
              do j=1,ny
                do k=1,nz
                   if ( p3(i,j,k).lt.0. ) p3(i,j,k) = 1050.
                enddo
              enddo
            enddo

c           Don't care about other stuff - finish subroutine
            goto 800

         endif

c        Get full grid info - in particular staggering flag; set flag for 2D tracing
         if ( fieldname.eq.'PLEV' ) then
            varname = 'PS'
            stagz   = 0.
            call getdef(fid,varname,ndim,mdv,vardim,
     >                  varmin,varmax,stag,ierr)
            if (ierr.ne.0) goto 900
            call getcfn(fid,cstfile,ierr)
            if (ierr.ne.0) goto 903
            call cdfopn(cstfile,cstid,ierr)
            if (ierr.ne.0) goto 903
            call getlevs(cstid,vardim(3),tmp,tmp,tmp,tmp,ierr)
            if (ierr.ne.0) goto 903
            call clscdf(cstid,ierr)
            if (ierr.ne.0) goto 903
            
         elseif ( ( fieldname.eq.'PLAY' ).or.( fieldname.eq.'P' ) ) then       
            varname = 'PS'
            stagz   = -0.5
            call getdef(fid,varname,ndim,mdv,vardim,
     >                  varmin,varmax,stag,ierr)
            if (ierr.ne.0) goto 900
            call getcfn(fid,cstfile,ierr)
            if (ierr.ne.0) goto 903
            call cdfopn(cstfile,cstid,ierr)
            if (ierr.ne.0) goto 903
            call getlevs(cstid,vardim(3),tmp,tmp,tmp,tmp,ierr)
            if (ierr.ne.0) goto 903
            call clscdf(cstid,ierr)
            if (ierr.ne.0) goto 903

         else
            varname=fieldname
            call getdef(fid,varname,ndim,mdv,vardim,
     >                  varmin,varmax,stag,ierr)
            if (ierr.ne.0) goto 900
            if (vardim(3).eq.1) is2d = 1
         endif

c        Get time information (check if time is correct)
         call gettimes(fid,times,ntimes,ierr)
         if (ierr.ne.0) goto 901
         isok=0
         do i=1,ntimes
            if (abs(time-times(i)).lt.eps) then
               isok = 1
               rtime = times(i)
            elseif (timecheck.eq.'no') then
               isok = 1
               rtime = times(1)
            endif
         enddo
         if ( isok.eq.0) goto 905

c        If 2D tracing requested: take dummay values for PS, AKLEV,BKLEV,AKLAY,BKLAY
         if ( is2d.eq.1 ) then
            
            do i=1,nx
               do j=1,ny
                  ps(i,j) = 1050.
               enddo
            enddo
            
            do k=1,nz
               aklev(k) = 0.
               bklev(k) = real(nz-k)/real(nz-1) 
               aklay(k) = 0.
               bklay(k) = real(nz-k)/real(nz-1) 
            enddo

c        3D tracing - read PS, AKLEV,BKLEV,AKLAY;BKLAY
         else

c           Read the level coefficients from the constants file
            call getcfn(fid,cstfile,ierr)
            if (ierr.ne.0) goto 903
            call cdfopn(cstfile,cstid,ierr)
            if (ierr.ne.0) goto 903
            call getlevs(cstid,vardim(3),aklev,bklev,aklay,bklay,ierr)
            if (ierr.ne.0) goto 903
            call clscdf(cstid,ierr)
            if (ierr.ne.0) goto 903

c           Check whether PS is needed to get the 3d pressure field
            plev = 1
            do i=1,nz
              if ( (abs(stagz).lt.eps).and.(abs(bklev(i)).gt.eps) ) then
                plev = 0
              endif
              if ( (abs(stagz).gt.eps).and.(abs(bklay(i)).gt.eps) ) then
                plev = 0
              endif
            enddo

c           Read surface pressure if needed
            if ( plev.ne.1 ) then
              varname='PS'
              call getdat(fid,varname,rtime,0,ps,ierr)
              if (ierr.ne.0) goto 904
            else
              do i=1,nx
                do j=1,ny
                    ps(i,j) = 1000.
                enddo
              enddo
            endif

         endif

c        Calculate layer and level pressures
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  if ( abs(stagz).lt.eps ) then
                     p3(i,j,k)=aklev(k)+bklev(k)*ps(i,j)
                  else
                     p3(i,j,k)=aklay(k)+bklay(k)*ps(i,j)
                  endif
               enddo
            enddo
         enddo

c        Set the ak and bk for the vertical grid
         do k=1,nz
            if ( abs(stagz).lt.eps ) then
               ak(k)=aklev(k)
               bk(k)=bklev(k)
            else
               ak(k)=aklay(k)
               bk(k)=bklay(k)
            endif
         enddo

      endif
      
c     Exit point for subroutine
 800  continue
      return
      
c     Exception handling
 900  print*,'Cannot retrieve grid dimension from ',fid
      stop
 901  print*,'Cannot retrieve grid parameters from ',fid
      stop
 902  print*,'Grid inconsistency detected for ',fid
      stop
 903  print*,'Problem with level coefficients from ',trim(cstfile)
      stop
 904  print*,'Cannot read surface pressure from ',trim(cstfile)
      stop
 905  print*,'Cannot find time ',time,' on ',fid
      stop
 906  print*,'Unable to get grid info ',fid
      stop
      
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

c     Auxiliary variables
      integer      isok
      integer      i,j,k
      integer      nz1
      real         rtime

c     Read variable definition - for P, PLEV and PLAY: load also ak,bk
      if ( ( fieldname.eq.'PLEV' ).or.
     >     ( fieldname.eq.'PLAY' ).or.
     >     ( fieldname.eq.'P'    ) )
     >then
         call getcfn(fid,cstfile,ierr)
         if (ierr.ne.0) goto 905
         call cdfopn(cstfile,cstid,ierr)
         if (ierr.ne.0) goto 905
         call getlevs(cstid,nz1,aklev,bklev,aklay,bklay,ierr)
         if (ierr.ne.0) goto 905
         call clscdf(cstid,ierr)
         if (ierr.ne.0) goto 905
         varname = 'PS'
         call getdef(fid,varname,ndim,mdv,vardim,
     >               varmin,varmax,stag,ierr)
         vardim(3) = nz1
         if (ierr.ne.0) goto 900

      else
         varname = fieldname
         call getdef(fid,varname,ndim,mdv,vardim,
     >               varmin,varmax,stag,ierr)
         if (ierr.ne.0) goto 900
         stagz=stag(3)
      endif

c     Get time information (set time to first one in the file)
      call gettimes(fid,times,ntimes,ierr)
      if (ierr.ne.0) goto 902
      isok=0
      do i=1,ntimes
         if (abs(time-times(i)).lt.eps) then
            isok = 1
            rtime = times(i)
         elseif (timecheck.eq.'no') then
            isok = 1
            rtime = times(1)
         endif
      enddo
      if ( isok.eq.0 ) goto 904

c     Read  field
      if ( ( fieldname.eq.'P' ).or.(fieldname.eq.'PLAY') )  then       ! P, PLAY
         stagz   = -0.5
         varname = 'PS'
         call getdat(fid,varname,rtime,0,ps,ierr)
         if (ierr.ne.0) goto 903
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  field(i,j,k)=aklay(k)+bklay(k)*ps(i,j)
               enddo
            enddo
         enddo
         
      elseif ( fieldname.eq.'PLEV' )  then                             ! PLEV
         stagz   = 0.
         varname = 'PS'
         call getdat(fid,varname,rtime,0,ps,ierr)
         if (ierr.ne.0) goto 903
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  field(i,j,k)=aklev(k)+bklev(k)*ps(i,j)
               enddo
            enddo
         enddo

      else                                                             ! Other fields
         varname=fieldname
         call getdat(fid,varname,rtime,0,field,ierr)
         if (ierr.ne.0) goto 903

      endif

c     If the field is 2D, expand it to 3D - simple handling of 2D tracing
      if ( vardim(3).eq.1 ) then
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  field(i,j,k) = field(i,j,1)
               enddo
            enddo
         enddo
      endif
         


c     Exception handling
      return
      
 900  print*,'Cannot retrieve definition for ',trim(varname),'  ',fid
      stop
 901  print*,'Grid inconsistency detected for ',trim(varname),'  ',fid
      stop
 902  print*,'Cannot retrieve time for ',trim(varname),'  ',fid
      stop
 903  print*,'Cannot load wind component ',trim(varname),'  ',fid
      stop
 904  print*,'Cannot load time ',time,' for ',trim(varname),'  ',fid
      stop
 905  print*,'Cannot load time vertical grid AK, BK from file  ',fid
      stop
      
      end

c     ------------------------------------------------------------
c     Close input file
c     ------------------------------------------------------------

      subroutine input_close(fid)

c     Close the input file with file identifier <fid>.

      implicit none

c     Declaration of subroutine parameters
      integer fid

c     Auxiliary variables
      integer ierr

c     Close file
      call clscdf(fid,ierr)
 
      end
      
c     ------------------------------------------------------------
c     Get a list of variables on netCDF file
c     ------------------------------------------------------------

      subroutine input_getvars(fid,vnam,nvars)

c     List of variables on netCDF file

      implicit none

c     Declaration of subroutine parameters
      integer      fid
      integer      nvars
      character*80 vnam(200)

c     Auxiliary variables
      integer ierr
      
c     Get list and save      
      call getvars(fid,nvars,vnam,ierr)
 
      end
