      PROGRAM reformat
      
c     ***********************************************************************
c     * Change format of a trajectory file                                  *
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
      open(10,file='reformat.param')
       read(10,*) inpfile
       read(10,*) outfile
       read(10,*) ntra,ntim,ncol
      close(10)
      
c     Determine the formats
      call mode_tra(inpmode,inpfile)
      if (inpmode.eq.-1) inpmode=1
      call mode_tra(outmode,outfile)
      if (outmode.eq.-1) outmode=1

c     Allocate memory
      allocate(tra(ntra,ntim,ncol),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra      ***' 

c     Read inpufile
      call ropen_tra(fid,inpfile,ntra,ntim,ncol,refdate,vars,inpmode)
      call read_tra (fid,tra,ntra,ntim,ncol,inpmode)
      call close_tra(fid,inpmode)
    
c     Write output file
      call wopen_tra(fid,outfile,ntra,ntim,ncol,refdate,vars,outmode)
      call write_tra(fid,tra,ntra,ntim,ncol,outmode)
      call close_tra(fid,outmode)
      
      end

      

      
