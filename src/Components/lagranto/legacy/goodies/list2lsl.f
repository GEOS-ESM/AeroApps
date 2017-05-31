      PROGRAM list2lsl
      
c     ***********************************************************************
c     * Convert a lat/lon/p list to a trajectory file                       *
c     * Michael Sprenger / Spring, summer 2010                              *
c     ***********************************************************************

      implicit none
      
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
      real                                   timevalue   ! Time 

c     Auxiliary variables
      integer                                inpmode
      integer                                outmode
      integer                                stat
      integer                                fid
      integer                                i
      
c     ----------------------------------------------------------------------
c     Do the reformating
c     ----------------------------------------------------------------------

c     Read parameters
      open(10,file='list2lsl.param')
       read(10,*) inpfile
       read(10,*) outfile
       read(10,*) ntra
       read(10,*) (refdate(i),i=1,6)
       read(10,*) timevalue
      close(10)
      
c     Determine the formats
      call mode_tra(outmode,outfile)
      if (outmode.eq.-1) outmode=1

c     Set parameters for output file
      ntim=1
      ncol=4
      vars(1)='time'
      vars(2)='lon'
      vars(3)='lat'
      vars(4)='p'

c     Allocate memory
      allocate(tra(ntra,ntim,ncol),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra      ***' 

c     Read inpufile
      fid = 10
      open(fid,file=inpfile)
       do i=1,ntra
          tra(i,1,1) = timevalue
          read(fid,*) tra(i,1,2),tra(i,1,3),tra(i,1,4)
       enddo
      close(fid)

c     Write output file
      call wopen_tra(fid,outfile,ntra,ntim,ncol,refdate,vars,outmode)
      call write_tra(fid,tra,ntra,ntim,ncol,outmode)
      call close_tra(fid,outmode)
      
      end

      

      
