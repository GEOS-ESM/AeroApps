      PROGRAM tracal
      
c     ***********************************************************************
c     * Calculations with trajectories                                      *
c     * Michael Sprenger / Spring, summer 2010                              *
c     ***********************************************************************

      use evaluate

      implicit none

c     ----------------------------------------------------------------------
c     Declaration of variables
c     ----------------------------------------------------------------------

c     Input and output format for trajectories (see iotra.f)
      character*80                           inpfile     ! Input filename
      character*80                           outfile     ! Output filename
      character*80                           expr        ! Expression for calculation

c     Trajectories
      integer                                ntra        ! Number of trajectories
      integer                                ntim        ! Number of times
      integer                                ncol        ! Number of columns
      real,allocatable, dimension (:,:,:) :: trainp      ! Trajectories (ntra,ntim,ncol  )
      real,allocatable, dimension (:,:,:) :: traout      ! Trajectories (ntra,ntim,ncol+1)
      character*80                           vars(100)   ! Variable names
      integer                                refdate(6)  ! Reference date

c     Auxiliary variables
      integer                                inpmode
      integer                                outmode
      integer                                stat
      integer                                fid
      integer                                i,j,k
      character (len=24)                     col,new
      real                                   value
      character                              ch
      integer                                ileft,iright

c     ----------------------------------------------------------------------
c     Do the reformating
c     ----------------------------------------------------------------------

c     Read parameters
      open(10,file='tracal.param')
       read(10,*) inpfile
       read(10,*) outfile
       read(10,*) expr
       read(10,*) ntra,ntim,ncol
      close(10)
      
c	  Get the name of the output field
	  ileft  = 1
	  iright = 0
	  do i=1,80
	  	if ( expr(i:i).eq.'=' ) iright = i
	  enddo
	  if ( iright.eq.0 ) then
	      vars(ncol+1) = 'CALC'
	      expr = 'CALC='//trim(expr)
	  else
	      vars(ncol+1) = trim( expr(1:iright-1) )
	  endif
	  new = vars( ncol+1 )

      print*,'inp  = ',trim(inpfile)
      print*,'out  = ',trim(outfile)
      print*,'expr = ',trim(expr),' ---> ',trim( new )

c     Determine the formats
      call mode_tra(inpmode,inpfile)
      if (inpmode.eq.-1) inpmode=1
      call mode_tra(outmode,outfile)
      if (outmode.eq.-1) outmode=1

c     Allocate memory
      allocate(trainp(ntra,ntim,ncol),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array trainp   ***'
      allocate(traout(ntra,ntim,ncol+1),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array traout   ***'

c     Read inpufile
      call ropen_tra(fid,inpfile,ntra,ntim,ncol,refdate,vars,inpmode)
      call read_tra (fid,trainp,ntra,ntim,ncol,inpmode)
      call close_tra(fid,inpmode)

c	  Copy to output trajectory
	  do i=1,ntra
	  	do j=1,ntim
	  	  do k=1,ncol
	  	     traout(i,j,k) = trainp(i,j,k)
	  	  enddo
	  	enddo
      enddo

c	  Loop over all trajectories
	  do i=1,ntra
	  	do j=1,ntim
	  		do k=1,ncol

c	          Attribute trajectory values to symbols
	          col = vars(k)
 	  		  call defparam( col, traout(i,j,k) )

	  		enddo

c	        Evaluate expression
	        call evaleqn(expr)

c	        Get the result
	        call getparam(new,value)
	        traout(i,j,ncol+1) = value

	    enddo
	  enddo
    
c     Write output file
      call wopen_tra(fid,outfile,ntra,ntim,ncol+1,refdate,vars,outmode)
      call write_tra(fid,traout,ntra,ntim,ncol+1,outmode)
      call close_tra(fid,outmode)
      
      end

      

      
