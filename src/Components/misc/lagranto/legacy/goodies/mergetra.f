      PROGRAM mege
      
c     ***********************************************************************
c     * Merge two trajectory files                                          *
c     * Michael Sprenger / Spring, summer 2010                              *
c     ***********************************************************************

      implicit none

c     ----------------------------------------------------------------------
c     Declaration of parameters
c     ----------------------------------------------------------------------

c     Numerical epsilon
      real      eps
      parameter (eps=0.0001)

c     ----------------------------------------------------------------------
c     Declaration of variables
c     ----------------------------------------------------------------------

c     Merging mode
      integer                                mode
      character*20                           datecheck

c     Input and output format for trajectories (see iotra.f)
      character*80                           inpfile1       ! Input filename 1
      character*80                           inpfile2       ! Input filename 2
      character*80                           outfile        ! Output filename
      integer                                inpmode1       ! Input format 1
      integer                                inpmode2       ! Input format 2
      integer                                outmode        ! Output format
      
c     Trajectories
      integer                                ntra_inp1       ! Number of trajectories
      integer                                ntim_inp1       ! Number of times
      integer                                ncol_inp1       ! Number of columns
      real,allocatable, dimension (:,:,:) :: tra_inp1        ! Trajectories (ntra,ntim,ncol)
      character*80                           vars_inp1(100)  ! Variable names
      real                                   time_inp1(5000) ! Times of input trajectory
      integer                                refdate1(6)     ! Reference date

      integer                                ntra_inp2       ! Number of trajectories
      integer                                ntim_inp2       ! Number of times
      integer                                ncol_inp2       ! Number of columns
      real,allocatable, dimension (:,:,:) :: tra_inp2        ! Trajectories (ntra,ntim,ncol)
      character*80                           vars_inp2(100)  ! Variable names
      real                                   time_inp2(5000) ! Times of input trajectory
      integer                                refdate2(6)     ! Reference date

      integer                                ntra_out       ! Number of trajectories
      integer                                ntim_out       ! Number of times
      integer                                ncol_out       ! Number of columns
      real,allocatable, dimension (:,:,:) :: tra_out        ! Trajectories (ntra,ntim,ncol)
      character*80                           vars_out(100)  ! Variable names
  
      real                                   time_out(10000)! Times of output trajectory
      integer                                refdate(6)     ! Reference date

c     Auxiliary variables
      integer                                i,j,k
      integer                                ind(10000)
      integer                                itr(10000)
      integer                                isok
      integer                                stat
      integer                                fid
      real                                   rswap
      integer                                iswap

c     ----------------------------------------------------------------------
c     Read and handle parameters
c     ----------------------------------------------------------------------

c     Read parameters
      open(10,file='mergetra.param')
       read(10,*) inpfile1
       read(10,*) inpfile2
       read(10,*) outfile
       read(10,*) ntra_inp1,ntim_inp1,ncol_inp1
       read(10,*) ntra_inp2,ntim_inp2,ncol_inp2
       read(10,*) datecheck
      close(10)

c     Determine the formats
      call mode_tra(inpmode1,inpfile1)
      if (inpmode1.eq.-1) inpmode1=1
      call mode_tra(inpmode2,inpfile2)
      if (inpmode2.eq.-1) inpmode2=1
      call mode_tra(outmode,outfile)
      if (outmode.eq.-1) outmode=1

c     ----------------------------------------------------------------------
c     Read input trajectories
c     ----------------------------------------------------------------------

c     Allocate memory
      allocate(tra_inp1(ntra_inp1,ntim_inp1,ncol_inp1),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra_inp1    ***' 
      allocate(tra_inp2(ntra_inp2,ntim_inp2,ncol_inp2),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra_inp2    ***' 

c     Read inpufile
      call ropen_tra(fid,inpfile1,ntra_inp1,ntim_inp1,ncol_inp1,
     >               refdate1,vars_inp1,inpmode1)
      call read_tra (fid,tra_inp1,ntra_inp1,ntim_inp1,ncol_inp1,
     >               inpmode1)
      call close_tra(fid,inpmode1)

      call ropen_tra(fid,inpfile2,ntra_inp2,ntim_inp2,ncol_inp2,
     >               refdate2,vars_inp2,inpmode2)
      call read_tra (fid,tra_inp2,ntra_inp2,ntim_inp2,ncol_inp2,
     >               inpmode2)
      call close_tra(fid,inpmode2)

c     Get the times of the trajectories
      do i=1,ntim_inp1
         time_inp1(i) = tra_inp1(1,i,1)
      enddo
      do i=1,ntim_inp2
         time_inp2(i) = tra_inp2(1,i,1)
      enddo

c     Check format 
      if (  vars_inp1(1).ne.'time' ) goto 990
      if ( (vars_inp1(2).ne.'lon').and.
     >     (vars_inp1(2).ne.'xpos') ) goto 990
      if ( (vars_inp1(3).ne.'lat').and.
     >     (vars_inp1(3).ne.'ypos') ) goto 990
      if ( (vars_inp1(4).ne.'p'  ).and.
     >     (vars_inp1(4).ne.'ppos') ) goto 990

      if (  vars_inp2(1).ne.'time' ) goto 990
      if ( (vars_inp2(2).ne.'lon').and.
     >     (vars_inp2(2).ne.'xpos') ) goto 990
      if ( (vars_inp2(3).ne.'lat').and.
     >     (vars_inp2(3).ne.'ypos') ) goto 990
      if ( (vars_inp2(4).ne.'p'  ).and.
     >     (vars_inp2(4).ne.'ppos') ) goto 990

c     Check whether the reference dates are equal
      if ( datecheck.eq.'datecheck' ) then
        do i=1,5
          if ( refdate1(i).ne.refdate2(i) ) then
            print*,' ERROR; Reference dates must be the same... Stop'
            print*,(refdate1(j),j=1,5)
            print*,(refdate2(j),j=1,5)
            stop
          endif
        enddo
      endif
            
c     ----------------------------------------------------------------------
c     Decide what to do
c     ----------------------------------------------------------------------

c     Init the mode
      mode = 0

c     ---- Check whether only the columns should be combined (mode 1) -------
      print*,'Testing for mode 1 (combine columns)'

      if ( ntim_inp1.ne.ntim_inp2 ) goto 100
      print*,'  -> ntim       ok'

      if ( ntra_inp1.ne.ntra_inp2 ) goto 100
      print*,'  -> ntra       ok'

      do i=1,ntim_inp1
         if ( time_inp1(i).ne.time_inp2(i) ) goto 100
      enddo
      print*,'  -> times      ok'

      do i=1,ntra_inp1
         do j=1,ntim_inp1
            if ( tra_inp1(i,j,1).ne.tra_inp2(i,j,1) ) goto 100
            if ( tra_inp1(i,j,2).ne.tra_inp2(i,j,2) ) goto 100
            if ( tra_inp1(i,j,3).ne.tra_inp2(i,j,3) ) goto 100
            if ( tra_inp1(i,j,4).ne.tra_inp2(i,j,4) ) goto 100
         enddo
      enddo
      print*,'  -> lon,lat,p  ok     => Mode 1 accepted'
      
      mode = 1
      goto 130


 100  continue

c     ---- Check whether second file to appended to first one (mode 2) -----
      print*,'Testing for mode 2 (append file 2 to file 1)'
      
      if ( ntim_inp1.ne.ntim_inp2 ) goto 110
      print*,'  -> ntim       ok'

      if ( ncol_inp1.ne.ncol_inp2 ) goto 110
      print*,'  -> ncol       ok'

      do i=1,ntim_inp1
         if ( time_inp1(i).ne.time_inp2(i) ) goto 110
      enddo
      print*,'  -> times      ok'

      do i=1,ncol_inp1
         if ( vars_inp1(i).ne.vars_inp2(i) ) goto 110
      enddo
      print*,'  -> vars       ok    => Mode 2 accepted'
      
      mode = 2
      goto 130

 110  continue

c     ----- Check whether to combine different times (mode 3) --------------
      print*,'Testing for mode 3 (combining times)'
      
      if ( ntra_inp1.ne.ntra_inp2 ) goto 120
      print*,'  -> ntra       ok'

      if ( ncol_inp1.ne.ncol_inp2 ) goto 120
      print*,'  -> ncol       ok'

      do i=1,ncol_inp1
         if ( vars_inp1(i).ne.vars_inp2(i) ) goto 120
      enddo
      print*,'  -> vars       ok    => Mode 3 accepted'

      mode = 3
      goto 130
      
 120  continue

c     -----  Stop if no valid merging mode could be determined -------------
      print*,' ERROR: could not determine merging mode ...'
      stop


c     Exit point for mode decision
 130  continue       

c     ----------------------------------------------------------------------
c     Merging depending on mode
c     ----------------------------------------------------------------------

c     ------ Merging of columns (mode 1) -----------------------------------
      if ( mode.ne.1 ) goto 200

c     Get the complete list of additional field
      do i=1,ncol_inp1
         ncol_out           = ncol_out + 1
         vars_out(ncol_out) = vars_inp1(i)
         ind(ncol_out)      = i 
         itr(ncol_out)      = 1
      enddo
      do i=5,ncol_inp2
         
         isok=1
         do j=1,ncol_out
            if ( vars_inp2(i).eq.vars_out(j) ) isok=0
         enddo
         
         if ( isok.eq.1 ) then
            ncol_out           = ncol_out + 1
            vars_out(ncol_out) = vars_inp2(i)
            ind(ncol_out)      = i 
            itr(ncol_out)      = 2
         endif

      enddo      

c     Allocate memory for output trajectory
      ntra_out = ntra_inp1
      ntim_out = ntim_inp1
      allocate(tra_out(ntra_out,ntim_out,ncol_out),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra_out    ***' 

c     Save the trajectory
      do i=1,ntra_out
         do j=1,ntim_out
            do k=1,ncol_out
               
               if ( itr(k).eq.1) then
                  tra_out(i,j,k) = tra_inp1(i,j,ind(k))
               else
                  tra_out(i,j,k) = tra_inp2(i,j,ind(k))
               endif

            enddo
         enddo
      enddo   

 200  continue

c     ------ Appending files (mode 2) --------------------------------------
      if ( mode.ne.2 ) goto 210

c     Allocate memory for output trajectory
      ntra_out = ntra_inp1 + ntra_inp2
      ntim_out = ntim_inp1
      ncol_out = ncol_inp1
      allocate(tra_out(ntra_out,ntim_out,ncol_out),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra_out    ***' 

c     Set the field names for the output trajectory
      do i=1,ncol_out
         vars_out(i) = vars_inp1(i)
      enddo

c     Save the trajectory
      do i=1,ntra_out
         do j=1,ntim_out
            do k=1,ncol_out
               
               if ( i.le.ntra_inp1 ) then
                  tra_out(i,j,k) = tra_inp1(i          ,j,k)
               else
                  tra_out(i,j,k) = tra_inp2(i-ntra_inp1,j,k)
               endif

            enddo
         enddo
      enddo

 210  continue

c     ------ Combining times (mode 3) --------------------------------------
      if ( mode.ne.3 ) goto 220

c     Get a list of all output times
      ntim_out = 0
      do i=1,ntim_inp1
         isok = 1
         do j=1,ntim_out
            if ( time_inp1(i).eq.time_out(j) ) isok=0
         enddo
         if (isok.eq.1 ) then
            ntim_out           = ntim_out + 1
            time_out(ntim_out) = time_inp1(i)
            itr(ntim_out)      = 1
            ind(ntim_out)      = i
         endif
      enddo
      do i=1,ntim_inp2
         isok = 1
         do j=1,ntim_out
            if ( time_inp2(i).eq.time_out(j) ) isok=0
         enddo
         if (isok.eq.1 ) then
            ntim_out           = ntim_out + 1
            time_out(ntim_out) = time_inp2(i)
            itr(ntim_out)      = 2
            ind(ntim_out)      = i
         endif
      enddo

c     Sort the times
      do i=1,ntim_out
         do j=i+1,ntim_out
            if ( time_out(j).lt.time_out(i) ) then

               rswap       = time_out(i)
               time_out(i) = time_out(j)
               time_out(j) = rswap

               iswap       = itr(i)
               itr(i)      = itr(j)
               itr(j)      = iswap

               iswap       = ind(i)
               ind(i)      = ind(j)
               ind(j)      = iswap

            endif
         enddo
      enddo

c     Allocate memory for output trajectory
      ntra_out = ntra_inp1 
      ncol_out = ncol_inp1
      allocate(tra_out(ntra_out,ntim_out,ncol_out),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra_out    ***' 

c     Set the field names for the output trajectory
      do i=1,ncol_out
         vars_out(i) = vars_inp1(i)
      enddo

c     Save the trajectory
      do i=1,ntra_out
         do j=1,ntim_out
            do k=1,ncol_out
               
               if ( itr(j).eq.1 ) then
                  tra_out(i,j,k) = tra_inp1(i,ind(j),k)
               else
                  tra_out(i,j,k) = tra_inp2(i,ind(j),k)
               endif

            enddo
         enddo
      enddo



 220  continue


c     ----------------------------------------------------------------------
c     Write the output trajectory
c     ----------------------------------------------------------------------
 
      call wopen_tra(fid,outfile,ntra_out,ntim_out,ncol_out,
     >               refdate1,vars_out,outmode)
      call write_tra(fid,tra_out,ntra_out,ntim_out,ncol_out,outmode)
      call close_tra(fid,outmode)
 
c     ----------------------------------------------------------------------
c     Error handling
c     ----------------------------------------------------------------------
      
      stop
   
 990  print*,'First columns must be <time,lon,lat,p>... Stop'
      stop
      

      end







