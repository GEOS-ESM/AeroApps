      PROGRAM reformat
      
c     ***********************************************************************
c     * Change format of a trajectory file                                  *
c     * Michael Sprenger / Spring, summer 2010                              *
c     ***********************************************************************

      implicit none
      
c     ----------------------------------------------------------------------
c     Declaration of variables
c     ----------------------------------------------------------------------

c     Mode
      character*80                           mode

c     Input and output files
      character*80                           inpfile     ! Input filename
      character*80                           outfile     ! Output filename

c     Trajectories
      integer                                ntra        ! Number of trajectories
      integer                                ntim        ! Number of times
      integer                                ncol        ! Number of columns
      real,allocatable, dimension (:,:,:) :: tra         ! Trajectories (ntra,ntim,ncol)
      character*80                           vars(100)   ! Variable names
      integer                                refdate(6)  ! Reference date
      real,allocatable, dimension (:)     :: num         ! Output number

c     Auxiliary variables
      integer                                inpmode
      integer                                stat
      integer                                fid
      integer                                i,j,k
      real                                   dist
      real                                   lon0,lat0,lon1,lat1

c     Externals
      real                                   sdis
      external                               sdis

c     ----------------------------------------------------------------------
c     Preparations
c     ----------------------------------------------------------------------

c     Read parameters
      open(10,file='traj2num.param')
       read(10,*) inpfile
       read(10,*) outfile
       read(10,*) ntra,ntim,ncol
       read(10,*) mode
      close(10)
      
c     Check that a valid mode is selected
      if ( mode.eq.'boost' ) goto 10
      print*,' Unknown mode ',trim(mode)
      stop
 10   continue

c     Determine the formats
      call mode_tra(inpmode,inpfile)
      if (inpmode.eq.-1) inpmode=1

c     Allocate memory
      allocate(tra(ntra,ntim,ncol),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra      ***' 
      allocate(num(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array num      ***'

c     Read inpufile
      call ropen_tra(fid,inpfile,ntra,ntim,ncol,refdate,vars,inpmode)
      call read_tra (fid,tra,ntra,ntim,ncol,inpmode)
      call close_tra(fid,inpmode)
    
c     ----------------------------------------------------------------------
c     Mode = 'boost': get the maximum distance traveled in one time step
c     ----------------------------------------------------------------------

      if ( mode.ne.'boost') goto 100

      do i=1,ntra

        num(i) = 0.

        do j=2,ntim

c          Get spherical distance between data points
           lon0 = tra(i,j-1,2)
           lat0 = tra(i,j-1,3)
           lon1 = tra(i,j  ,2)
           lat1 = tra(i,j  ,3)
           dist = sdis( lon1,lat1,lon0,lat0 )

           if ( dist.gt.num(i) ) num(i) = dist

        enddo

      enddo

 100  continue

c     ----------------------------------------------------------------------
c     Write output file
c     ----------------------------------------------------------------------

      open(10,file=outfile)
        do i=1,ntra
            write(10,*) num(i)
        enddo
      close(10)
      
      end

c     ***********************************************************************
c     * SUBROUTINES                                                         *
c     ***********************************************************************

c     -----------------------------------------------------------------------
c     Spherical distance between lat/lon points
c     -----------------------------------------------------------------------

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
      

      
