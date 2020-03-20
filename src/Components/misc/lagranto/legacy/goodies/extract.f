      PROGRAM extract
      
c     ***********************************************************************
c     * Extract columns of a trajectory file                                *
c     * Michael Sprenger / Spring, summer 2010                              *
c     ***********************************************************************

      implicit none

c     ----------------------------------------------------------------------
c     Declaration of parameters
c     ----------------------------------------------------------------------

c     Numerical epsiloN
      real      eps
      parameter (eps=0.0001)

c     ----------------------------------------------------------------------
c     Declaration of variables
c     ----------------------------------------------------------------------

c     Extraction mode
      character*80                           mode                             

c     Input and output format for trajectories (see iotra.f)
      character*80                           inpfile       ! Input filename
      character*80                           outfile       ! Output filename

c     Trajectories
      integer                                ntra_inp      ! Number of trajectories
      integer                                ntim_inp      ! Number of times
      integer                                ncol_inp      ! Number of columns
      real,allocatable, dimension (:,:,:) :: tra_inp       ! Trajectories (ntra,ntim,ncol)
      character*80                           vars_inp(100) ! Variable names

      integer                                ntra_out      ! Number of trajectories
      integer                                ntim_out      ! Number of times
      integer                                ncol_out      ! Number of columns
      real,allocatable, dimension (:,:,:) :: tra_out       ! Trajectories (ntra,ntim,ncol)
      integer,allocatable, dimension (:)  :: ind           ! Index for selection
      integer,allocatable, dimension (:)  :: isok          ! Index for selection
      character*80                           vars_out(100) ! Variable names
  
      real                                   time_inp(10000) ! Times of input trajectory
      real                                   time_out(10000) ! Times of output trajectory
      integer                                refdate(6)    ! Reference date
      integer                                ind_time(10000) ! Index for time selection

c     Auxiliary variables
      integer                                inpmode
      integer                                outmode
      integer                                stat
      integer                                fid
      integer                                i,j,k,n,j0,j1
      character*80                           str
      character*80                           split_str(100)
      integer                                split_n
      integer                                isstr,nvars,ileft,iright
      character*80                           vars(100)
      character                              ch
      real                                   tmp0,tmp1
      integer                                ind1
	  character*2000                         linestr
	  integer                                istr(100)
      integer                                nstr
      character*80                           strsplit(100)
      integer                                flag

c     ----------------------------------------------------------------------
c     Read and handle parameters
c     ----------------------------------------------------------------------

c     Read parameters
      open(10,file='extract.param')
       read(10,*) inpfile
       read(10,*) outfile
       read(10,*) mode
       read(10,*) ntra_inp,ntim_inp,ncol_inp
       read(10,*) str
      close(10)

c     Split the input string
      isstr   = 0
      split_n = 0
      do i=1,80
         ch = str(i:i)
         if ( (isstr.eq.0).and.(ch.ne.' ') ) then
            isstr=1
            ileft=i
         elseif ( (isstr.eq.1).and.(ch.eq.' ') ) then
            isstr              = 0
            iright             = i-1
            split_n            = split_n+1
            split_str(split_n) = str(ileft:iright)
         endif
      enddo

c     Determine the formats
      call mode_tra(inpmode,inpfile)
      if (inpmode.eq.-1) inpmode=1
      call mode_tra(outmode,outfile)
      if ( (mode.ne.'-startf').and.(outmode.eq.-1) ) then
         outmode=1
      endif

c     ----------------------------------------------------------------------
c     Read input trajectories
c     ----------------------------------------------------------------------

c     Allocate memory
      allocate(tra_inp(ntra_inp,ntim_inp,ncol_inp),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra_inp    ***' 
      allocate(ind(ntra_inp),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array ind        ***' 
      allocate(isok(ntra_inp),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array isok       ***' 

c     Read inpufile
      fid = 10
      call ropen_tra(fid,inpfile,ntra_inp,ntim_inp,ncol_inp,
     >               refdate,vars_inp,inpmode)
      call read_tra (fid,tra_inp,ntra_inp,ntim_inp,ncol_inp,inpmode)
      call close_tra(fid,inpmode)

c     Check format 
      if (  vars_inp(1).ne.'time') goto 990
      if ( (vars_inp(2).ne.'lon').and.(vars_inp(2).ne.'xpos') ) goto 990
      if ( (vars_inp(3).ne.'lat').and.(vars_inp(3).ne.'ypos') ) goto 990
      if ( (vars_inp(4).ne.'p'  ).and.(vars_inp(4).ne.'ppos') ) goto 990

c     ----------------------------------------------------------------------
c     Option -vars : Extract columns of variables
c     ----------------------------------------------------------------------

      if ( mode.ne.'-var' ) goto 100

c     Set the first for columns of the output
      ncol_out   = 4
      vars_out(1)='time'
      vars_out(2)='lon'
      vars_out(3)='lat'
      vars_out(4)='p'
      ind(1)     =1
      ind(2)     =2
      ind(3)     =3
      ind(4)     =4

c     Get final list of extraction columns (set number of columns)
      do i=1,split_n

         if (split_str(i).eq.'to') then
            do j=1,ncol_inp
               if ( vars_inp(j).eq.split_str(i-1) ) j0 = j + 1
            enddo
            do j=1,ncol_inp
               if ( vars_inp(j).eq.split_str(i+1) ) j1 = j - 1
            enddo
            do j=j0,j1
               ncol_out           = ncol_out + 1
               vars_out(ncol_out) = vars_inp(j)
            enddo
         else
            ncol_out           = ncol_out + 1
            vars_out(ncol_out) = split_str(i)
         endif

      enddo

c     Set the dimensions of the output trajectory
      ntra_out = ntra_inp
      ntim_out = ntim_inp

c     Allocate memory for output trajectory
      allocate(tra_out(ntra_out,ntim_out,ncol_out),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra_out    ***' 

c     Extract <time,lon,lat,p> columns
      do i=1,ntra_out
         do j=1,ntim_out
            tra_out(i,j,1) = tra_inp(i,j,1)
            tra_out(i,j,2) = tra_inp(i,j,2)
            tra_out(i,j,3) = tra_inp(i,j,3)
            tra_out(i,j,4) = tra_inp(i,j,4)
         enddo
      enddo

c     Get indices for new columns (1..4 are already ok)
      do i=5,ncol_out
         ind(i)=0
      enddo
      do i=5,ncol_out
         do j=1,ncol_inp
            if ( vars_inp(j).eq.vars_out(i) ) ind(i) = j
         enddo
      enddo

c     Check if all selected columns are available
      do i=1,ncol_out
         if ( ind(i).eq.0 ) then
            print*,'Invalid column in ',trim(str)
            stop
         endif
      enddo

c     Extract the column
      do i=1,ntra_out
         do j=1,ntim_out
            do k=1,ncol_out
               tra_out(i,j,k) = tra_inp(i,j,ind(k))
            enddo
         enddo
      enddo

 100  continue

c     ----------------------------------------------------------------------
c     Option -times : Extract times of trajectories
c     ----------------------------------------------------------------------

      if ( mode.ne.'-time' ) goto 110

c     Set the dimension of the output trajectory
      ntim_out = 0

c     Get the list of times for the input trajectory
      do i=1,ntim_inp
         time_inp(i) = tra_inp(1,i,1)
      enddo

c     Get final list of extraction times (set number of times)
      do i=1,split_n

         if (split_str(i).eq.'to') then
            read(split_str(i-1),*) tmp0
            do j=1,ntim_inp
               if ( time_inp(j).eq.tmp0 ) j0 = j + 1
            enddo
            read(split_str(i+1),*) tmp0
            do j=1,ntim_inp
               if ( time_inp(j).eq.tmp0 ) j1 = j - 1
            enddo
            do j=j0,j1
               ntim_out           = ntim_out + 1
               time_out(ntim_out) = time_inp(j)
            enddo
         elseif (split_str(i).eq.'first') then
            ntim_out           = ntim_out + 1
            time_out(ntim_out) = time_inp(1)
         elseif (split_str(i).eq.'last') then
            ntim_out           = ntim_out + 1
            time_out(ntim_out) = time_inp(ntim_inp)
         else
            ntim_out           = ntim_out + 1
            read(split_str(i),*) tmp0
            time_out(ntim_out) = tmp0
         endif

      enddo

c     Get the indices of the selected times
      do i=1,ntim_out
         ind_time(i) = 0
      enddo
      do i=1,ntim_out
         do j=1,ntim_inp
            if ( abs(time_out(i)-time_inp(j)).lt.eps) ind_time(i) = j
         enddo
      enddo
      do i=1,ntim_out
         if ( ind_time(i).eq.0) then
            print*,' Invalid time ',time_out(i)
            stop
         endif
      enddo
      
c     Set dimensions of output trajectory
      ntra_out = ntra_inp
      ncol_out = ncol_inp

c     Allocate memory for output trajectory
      allocate(tra_out(ntra_out,ntim_out,ncol_out),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra_out    ***' 

c     Copy the selected times to the output trajectory
      do i=1,ntra_out
         do j=1,ntim_out
            do k=1,ncol_out
               tra_out(i,j,k) = tra_inp(i,ind_time(j),k)
            enddo
         enddo
      enddo

c     Copy meta information
      do i=1,ncol_out
         vars_out(i) = vars_inp(i)
      enddo

 110  continue

c     ----------------------------------------------------------------------
c     Option -tra : Extract trajectories by number
c     ----------------------------------------------------------------------

      if ( mode.ne.'-tra' ) goto 120

c     Set the dimension of the output trajectory
      ntra_out = 0

c     Get final list of extraction times (set number of times)
      do i=1,split_n

         if (split_str(i).eq.'to') then
            read(split_str(i-1),*) tmp0
            read(split_str(i+1),*) tmp1
            do j=nint(tmp0)+1,nint(tmp1)-1
               ntra_out        = ntra_out + 1
               ind(ntra_out)   = j   
            enddo
         elseif (split_str(i).eq.'first') then
            ntra_out           = ntra_out + 1
            ind(ntra_out)      = 1
         elseif (split_str(i).eq.'last') then
            ntra_out           = ntra_out + 1
            ind(ntra_out)      = ntra_inp
         else
            ntra_out           = ntra_out + 1
            read(split_str(i),*) tmp0
            ind(ntra_out)      = nint(tmp0) 
         endif

      enddo

c     Check whether selected trajectories are ok
      do i=1,ntra_out
         if ( (ind(i).lt.1).or.(ind(i).gt.ntra_inp) ) then
            print*,'Invalid trajectory selected ',ind(i)
            stop
         endif
      enddo

c     Set dimensions of output trajectory
      ntim_out = ntim_inp
      ncol_out = ncol_inp
     
c     Allocate memory for output trajectory
      allocate(tra_out(ntra_out,ntim_out,ncol_out),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra_out    ***' 

c     Copy the selected times to the output trajectory
      do i=1,ntra_out
         do j=1,ntim_out
            do k=1,ncol_out
               tra_out(i,j,k) = tra_inp(ind(i),j,k)
            enddo
         enddo
      enddo

c     Copy meta information
      do i=1,ncol_out
         vars_out(i) = vars_inp(i)
      enddo

 120  continue

c     ----------------------------------------------------------------------
c     Option -startf : Extract starting positions for the trajectory file
c     ----------------------------------------------------------------------

      if ( mode.ne.'-startf' ) goto 130

c     Set the first for columns of the output
      ncol_out   = 4
      vars_out(1)='time'
      vars_out(2)='lon'
      vars_out(3)='lat'
      vars_out(4)='p'
      ind(1)     =1
      ind(2)     =2
      ind(3)     =3
      ind(4)     =4

c     Set dimensions of output trajectory
      ntim_out = 1
      ntra_out = ntra_inp
     
c     Allocate memory for output trajectory
      allocate(tra_out(ntra_out,ntim_out,ncol_out),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra_out    ***' 

c     Copy the selected times to the output trajectory
      do i=1,ntra_out
         do k=1,ncol_out
            tra_out(i,1,k) = tra_inp(i,1,k)
         enddo
      enddo

c     Copy meta information
      do i=1,ncol_out
         vars_out(i) = vars_inp(i)
      enddo

 130  continue

c     ----------------------------------------------------------------------
c     Option -index : Extract trajectories by index file
c     ----------------------------------------------------------------------

      if ( mode.ne.'-index' ) goto 140

c     Read the index file
      open(10,file=str)

        ntra_out = 1
 142    read(10,*,end=141) ind(ntra_out)
        ntra_out = ntra_out + 1
        goto 142 
 141    continue
        ntra_out = ntra_out - 1

      close(10)

c     Check whether selected trajectories are ok
      do i=1,ntra_out
         if ( (ind(i).lt.1).or.(ind(i).gt.ntra_inp) ) then
            print*,'Invalid trajectory selected ',ind(i)
            stop
         endif
      enddo

c     Set dimensions of output trajectory
      ntim_out = ntim_inp
      ncol_out = ncol_inp
     
c     Allocate memory for output trajectory
      allocate(tra_out(ntra_out,ntim_out,ncol_out),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra_out    ***' 

c     Copy the selected times to the output trajectory
      do i=1,ntra_out
         do j=1,ntim_out
            do k=1,ncol_out
               tra_out(i,j,k) = tra_inp(ind(i),j,k)
            enddo
         enddo
      enddo

c     Copy meta information
      do i=1,ncol_out
         vars_out(i) = vars_inp(i)
      enddo

 140  continue

c     ----------------------------------------------------------------------
c     Option -boolean : Extract trajectories by boolean file
c     ----------------------------------------------------------------------

      if ( mode.ne.'-boolean' ) goto 150

c     Read the index file
      open(10,file=str)
        ntra_out = 0
        do i=1,ntra_inp
           read(10,*) ind1
           if ( ind1.eq.1 ) then
              ntra_out      = ntra_out + 1
              ind(ntra_out) = i
           endif
        enddo
      close(10)

c     Check whether selected trajectories are ok
      do i=1,ntra_out
         if ( (ind(i).lt.1).or.(ind(i).gt.ntra_inp) ) then
            print*,'Invalid trajectory selected ',ind(i)
            stop
         endif
      enddo

c     Set dimensions of output trajectory
      ntim_out = ntim_inp
      ncol_out = ncol_inp
     
c     Allocate memory for output trajectory
      allocate(tra_out(ntra_out,ntim_out,ncol_out),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra_out    ***' 

c     Copy the selected times to the output trajectory
      do i=1,ntra_out
         do j=1,ntim_out
            do k=1,ncol_out
               tra_out(i,j,k) = tra_inp(ind(i),j,k)
            enddo
         enddo
      enddo

c     Copy meta information
      do i=1,ncol_out
         vars_out(i) = vars_inp(i)
      enddo

 150  continue

c     ----------------------------------------------------------------------
c     Option -pattern : Extract trajectories which match a regular expression
c     ----------------------------------------------------------------------

	  if ( mode.ne.'-pattern' ) goto 160

c	  All times and columns are extracted
      ncol_out = ncol_inp
      ntim_out = ntim_inp
      ntra_out = 0

c     Split the search string
	  nstr   = 0
	  ileft  = 0
	  iright = 0
      do i=1,len_trim(str)
      	if ( (str(i:i).eq.' ').and.(ileft.eq.0) ) then
      	   ileft = ileft + 1
        elseif ( (str(i:i).ne.' ').and.(ileft.eq.0) ) then
           ileft  = i
           iright = 0
        elseif ( (str(i:i).ne.' ').and.(ileft.ne.0) ) then
           iright = i
        elseif ( (str(i:i).eq.' ').and.(ileft.ne.0) ) then
           nstr           = nstr + 1
           strsplit(nstr) = str(ileft:iright)
           ileft          = 0
           iright         = 0
        endif
      enddo
      if ( (ileft.ne.0).and.(iright.ne.0) ) then
           nstr           = nstr + 1
           strsplit(nstr) = str(ileft:iright)
           ileft          = 0
      endif

c     Loop over the trajectories - check for matching pattern
	  do n=1,ntra_inp

	    ind(n) = 0
	  	do i=1,ntim_inp

           write(linestr,'(1f7.2,f9.2,f8.2,i6,100f10.3)')
     >                   (tra_inp(n,i,j),j=1,3),             ! time, lon, lat
     >                   nint(tra_inp(n,i,4)),               ! p
     >                   (tra_inp(n,i,j),j=5,ncol_inp)       ! fields

           flag = 1
           do k=1,nstr
             istr(k) = index(trim(linestr),trim(strsplit(k)))
             if ( istr(k).eq.0 ) flag = 0
           enddo
	       if ( flag.eq.1 ) ind(n) = 1

	  	enddo
	  	if ( ind(n).eq.1 ) ntra_out = ntra_out + 1

	  enddo

c     Allocate memory for output trajectory
      allocate(tra_out(ntra_out,ntim_out,ncol_out),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra_out    ***'

c     Copy the selected times to the output trajectory
      ntra_out = 0
      do i=1,ntra_inp
      	if (ind(i).eq.1) then
      	    ntra_out = ntra_out + 1
         	do j=1,ntim_out
            	do k=1,ncol_out
               		tra_out(ntra_out,j,k) = tra_inp(i,j,k)
            	enddo
         	enddo
         endif
      enddo

c     Copy meta information
      do i=1,ncol_out
         vars_out(i) = vars_inp(i)
      enddo

 160  continue

c     ----------------------------------------------------------------------
c     Option -leaving : Extract all trajectories which leave domain
c     ----------------------------------------------------------------------

      if ( mode.ne.'-leaving' ) goto 170

c     Set dimensions of output trajectory
      ntim_out = ntim_inp
      ncol_out = ncol_inp
      ntra_out = 0

c     Copy the meta data
      do i=1,ncol_out
         vars_out(i) = vars_inp(i)
      enddo

c     Determine the number of trajectories leaving domain
      do i=1,ntra_inp
         isok(i) = 1
         do j=1,ntim_inp
            if ( tra_inp(i,j,4).lt.0. ) isok(i) = 0
         enddo         
         if ( isok(i).eq.0 ) then
            ntra_out = ntra_out + 1
         endif
      enddo
     
c     Allocate memory for output trajectory
      allocate(tra_out(ntra_out,ntim_out,ncol_out),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra_out    ***' 

c     Copy the selected trajectories to the output trajectory
      ntra_out = 0
      do i=1,ntra_inp
         if ( isok(i).eq.0 ) then
            ntra_out = ntra_out + 1
            do j=1,ntim_inp
               do k=1,ncol_out
                  tra_out(ntra_out,j,k) = tra_inp(i,j,k)
               enddo
            enddo
         endif
      enddo
         
c     Copy meta information

 170  continue

c     ----------------------------------------------------------------------
c     Option -staying : Extract all trajectories which stay in domain
c     ----------------------------------------------------------------------

      if ( mode.ne.'-staying' ) goto 180

c     Set dimensions of output trajectory
      ntim_out = ntim_inp
      ncol_out = ncol_inp
      ntra_out = 0

c     Copy the meta data
      do i=1,ncol_out
         vars_out(i) = vars_inp(i)
      enddo

c     Determine the number of trajectories staying in domain
      do i=1,ntra_inp
         isok(i) = 1
         do j=1,ntim_inp
            if ( tra_inp(i,j,4).lt.0. ) isok(i) = 0
         enddo         
         if ( isok(i).eq.1 ) then
            ntra_out = ntra_out + 1
         endif
      enddo
     
c     Allocate memory for output trajectory
      allocate(tra_out(ntra_out,ntim_out,ncol_out),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra_out    ***' 

c     Copy the selected trajectories to the output trajectory
      ntra_out = 0
      do i=1,ntra_inp
         if ( isok(i).eq.1 ) then
            ntra_out = ntra_out + 1
            do j=1,ntim_inp
               do k=1,ncol_out
                  tra_out(ntra_out,j,k) = tra_inp(i,j,k)
               enddo
            enddo
         endif
      enddo
         
c     Copy meta information

 180  continue

c     ----------------------------------------------------------------------
c     Write output trajectories
c     ----------------------------------------------------------------------

c     Write output as trajectory file
      if (outmode.ge.1) then

         call wopen_tra(fid,outfile,ntra_out,ntim_out,ncol_out,
     >                  refdate,vars_out,outmode)
         call write_tra(fid,tra_out,ntra_out,ntim_out,ncol_out,outmode)
         call close_tra(fid,outmode)

c     Write output as (lon, lat, p)-list
      else

         open(10,file=outfile)
         do i=1,ntra_out
            write(10,'(3f10.2)') tra_out(i,1,2),    ! lon
     >                           tra_out(i,1,3),    ! lat
     >                           tra_out(i,1,4)     ! p
         enddo
         close(10)

      endif



!c     ----------------------------------------------------------------------
c     Error handling
c     ----------------------------------------------------------------------
      
      stop
   
 990  print*,'First columns must be <time,lon,lat,p>... Stop'
      stop
      

      end

      

      
