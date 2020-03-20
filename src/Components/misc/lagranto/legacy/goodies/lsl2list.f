      PROGRAM lsl2list
      
c     ***********************************************************************
c     * Convert a trajectory file into a lat/lon/p list                     *
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
      integer                                stat
      integer                                fid
      integer                                i,j
      
c     ----------------------------------------------------------------------
c     Do the reformating
c     ----------------------------------------------------------------------

c     Read parameters
      open(10,file='lsl2list.param')
       read(10,*) inpfile
       read(10,*) outfile
       read(10,*) ntra,ntim,ncol
      close(10)
      
c     Determine the formats
      call mode_tra(inpmode,inpfile)
      if (inpmode.eq.-1) inpmode=1

c     Allocate memory
      allocate(tra(ntra,ntim,ncol),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra      ***' 

c     Read inpufile
      call ropen_tra(fid,inpfile,ntra,ntim,ncol,refdate,vars,inpmode)
      call read_tra (fid,tra,ntra,ntim,ncol,inpmode)
      call close_tra(fid,inpmode)

c     Write output file inpufile
      fid = 10
      open(fid,file=outfile)
       do i=1,ntra
          do j=1,ntim
             write(fid,'(f9.2,f8.2,i6)') 
     >               tra(i,j,2),tra(i,j,3),nint(tra(i,j,4))
          enddo
       enddo
      close(fid)

      end

      

      
