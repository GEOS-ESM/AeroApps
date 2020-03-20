      PROGRAM select

c     **************************************************************
c     * Select trajectories from LSL file                          *
c     * Michael Sprenger / January, February 2008                  *
c     **************************************************************

      implicit none

c     -------------------------------------------------------------- 
c     Declaration of parameters
c     --------------------------------------------------------------

c     Maximum number of columns per trajectory
      integer          maxcol
      parameter        (maxcol=100)
      
c     Maximum number of commands
      integer          maxcmd
      parameter        (maxcmd=10000)

c     -------------------------------------------------------------
c     Declaration of variables
c     --------------------------------------------------------------

c     Input and output files
      character*120     inp_lslfile                  ! Input lsl file
      character*120     out_lslfile                  ! Output lsl file
      character*120     inp_criteria                 ! Input file with criteria
      integer           inpmode                      ! Format of input file
      integer           outmode                      ! Format of output file
      character*80      outformat                    ! Trajectory/Boolean/Index of output   
      character*80      regionf                      ! Name of the regionfile

c     Trajectories
      integer          ntim                          ! Number of times 
      integer          ncol                          ! Number of columns (including time...)
      integer          ntrainp,ntraout               ! Number of trajectories
      character*80     vars(maxcol)                  ! Names of trajectory columns
      integer          basetime(6)                   ! Base time of trajectory (first line in lsl)
      real,allocatable, dimension (:,:,:) :: trainp  ! Input trajectories
      real,allocatable, dimension (:,:,:) :: traout  ! Output trajectories
      integer,allocatable,dimension (:)   :: selflag ! Flag for selection
      real,allocatable, dimension (:)     :: time    ! Times of the trajectory
      integer,allocatable,dimension (:,:) :: trigger ! Trigger column 
      integer          itrigger                      ! Column index for trigger
      character*80     addtrigger                    ! Flag whether to add trigger column

c     Command stack
      real             cmd(maxcmd)                   ! Decoded slection criterion
      integer          ncmd                          ! Number of commands

c     Common block for initialisation of polygon check
      real    tlonv(2000),vlat_c(2000),vlon_c(2000)
      real    xlat_c,xlon_c
      integer ibndry,nv_c
      data    ibndry/0/
      common  /spolybndry/vlat_c,vlon_c,nv_c,xlat_c,xlon_c,tlonv,ibndry

c     Auxiliary variables
      integer          stat                          ! Logical (state) variable
      integer          fid,fod,fcr                   ! File identifier for input and output
      integer          i,j,k                         ! Index counter
      integer          isok                          ! Flag for selected trajectory
      real             param(1000)                   ! List of parameters
      integer          nparam                        ! Number of parameters
      character*80     specialstr                    ! Name of special command
      integer          len
      integer,allocatable,dimension (:) :: trigger1  ! Trigger column 
      real,allocatable,dimension (:,:) :: trainp1    ! A single trajectory
      character        ch

c     -------------------------------------------------------------- 
c     Preparations
c     --------------------------------------------------------------

c     Write start message
      print*,'========================================================='
      print*,'              *** START OF PROGRAM SELECT ***'
      print*

c     Read parameter file
      open(10,file='select.param')
       read(10,*) inp_lslfile
       read(10,*) out_lslfile
       read(10,*) outformat
       read(10,*) inp_criteria
       read(10,*) ntrainp
       read(10,*) ntim
       read(10,*) ncol
       read(10,*) regionf
       read(10,*) addtrigger
      close(10)
      
c     Set the formats of the input and output files
      call mode_tra(inpmode,inp_lslfile)
      if (inpmode.eq.-1) inpmode=1
      call mode_tra(outmode,out_lslfile)
      if ( (outmode.eq.-1).and.(outformat.ne.'startf') ) then
         outmode=1
      endif

c     Allocate memory for a single trajectory
      allocate(trainp(ntrainp,ntim,ncol),stat=stat)
      if (stat.ne.0) stop '*** error allocating array trainp   ***'
      allocate(time(ntim),stat=stat)
      if (stat.ne.0) stop '*** error allocating array time     ***'
      allocate(selflag(ntrainp),stat=stat)
      if (stat.ne.0) stop '*** error allocating array selflag  ***'
      allocate(trigger(ntrainp,ntim),stat=stat)
      if (stat.ne.0) stop '*** error allocating array trigger  ***'
      allocate(trigger1(ntim),stat=stat)
      if (stat.ne.0) stop '*** error allocating array trigger1 ***'
      allocate(trainp1(ntim,ncol),stat=stat)
      if (stat.ne.0) stop '*** error allocating array trainp1  ***'

c     Read the input trajectory file
      call ropen_tra(fid,inp_lslfile,ntrainp,ntim,ncol,
     >               basetime,vars,inpmode)
      call read_tra (fid,trainp,ntrainp,ntim,ncol,inpmode)
      call close_tra(fid,inpmode)

c     Check that first four columns correspond to time,lon,lat,p
      if ( (vars(1).ne.'time' ).or.
     >     (vars(2).ne.'xpos' ).and.(vars(2).ne.'lon' ).or.
     >     (vars(3).ne.'ypos' ).and.(vars(3).ne.'lat' ).or.
     >     (vars(4).ne.'ppos' ).and.(vars(4).ne.'p'   ) )
     >then
         print*,' ERROR: problem with input trajectories ...'
         stop
      endif
      vars(1) = 'time'
      vars(2) = 'lon'
      vars(3) = 'lat'
      vars(4) = 'p'

c     Get the trajectory times from first trajectory
      do i=1,ntim
         time(i)=trainp(1,i,1)
      enddo

c     Init the trigger field - first check whether it is already available
      itrigger = 0
      do i=1,ncol
         if ( vars(i).eq.'TRIGGER' ) itrigger = i
      enddo

      if ( itrigger.ne.0 ) then
         do i=1,ntrainp
            do j=1,ntim
               trigger(i,j) = nint( trainp(i,j,itrigger) )
            enddo
         enddo
      else
         do i=1,ntrainp
            do j=1,ntim
               trigger(i,j) = 0
            enddo
         enddo
      endif

c     Write some info about the trajectory
      print*,'---- INPUT PARAMETERS -----------------------------------'
      write(*,*) 
      write(*,*) 'Input file    : ',trim(inp_lslfile)
      write(*,*) 'Output file   : ',trim(out_lslfile)
      write(*,*) 'Output format : ',trim(outformat)
      write(*,*) 'Criteria file : ',trim(inp_criteria)
      write(*,*) '# tra         : ',ntrainp
      write(*,*) '# time        : ',ntim
      write(*,*) '# col         : ',ncol
      write(*,*) 'Region file   : ',trim(regionf)
      write(*,*) 'Add trigger   : ',trim(addtrigger)
      
      print*
      print*,'---- INPUT TRAJECTORY FILE ------------------------------'
      print*
      write(*,'(1x,a12,i4,a10)') 'Vars       : ',1,trim(vars(1))
      do i=2,ncol
         write(*,'(1x,a12,i4,a10)') '             ',i,trim(vars(i))
      enddo
      print*
      write(*,'(1x,a12,i4,f10.2)') 'Time       : ',1,time(1)
      do i=2,ntim
         write(*,'(1x,a12,i4,f10.2)') '             ',i,time(i)
      enddo
      print*
      write(*,'(1x,a12,i4,i10)') 'Base date  : ',1,basetime(1)
      write(*,'(1x,a12,i4,i10)') '             ',2,basetime(2)
      write(*,'(1x,a12,i4,i10)') '             ',3,basetime(3)
      write(*,'(1x,a12,i4,i10)') '             ',4,basetime(4)
      write(*,'(1x,a12,i4,i10)') '             ',5,basetime(5)      
      print*
      write(*,'(1x,a12,i4,i10)') 'Time range : ',6,basetime(6)
      print*
      if ( itrigger.ne.0 ) then
         print*,'TRIGGER FIELD FOUND IN COLUMN ',itrigger
         print*
      endif

c     Read and decode the selection criterion
      fcr = 10
      open(fcr,file=inp_criteria)
      ncmd=maxcmd
      call decode(fcr,cmd,ncmd,vars,ncol,time,ntim,regionf)
      close(fcr)

      print*
      print*,'---- PSEUDO CODE FOR SELECTION --------------------------'
      print*
      call dumpcode(cmd,ncmd,vars,ncol,time,ntim)

c     -------------------------------------------------------------- 
c     Loop over all trajectories - selection
c     --------------------------------------------------------------

c     Prepare string and parameters for SPECIAL commands
      if ( cmd(1).eq.0 ) then

c        Get command string
         j   = 2
         len = nint(cmd(j))
         specialstr = ''
         do k=1,len
            j = j + 1
            specialstr = trim(specialstr)//char(nint(cmd(j)))
         enddo
         
c        Get paramters
         j      = j + 1
         nparam = nint(cmd(j))
         do k=1,nparam
            j        = j + 1
            param(k) = cmd(j)
         enddo
         
      endif
        
c     Init the counter for selected trajectories
      ntraout = 0

c     Loop over all trajectories
      do i=1,ntrainp

c       Copy a single trajectory to <trainp1> and <trigger1>
        do j=1,ntim
           do k=1,ncol
              trainp1(j,k) = trainp(i,j,k)
           enddo
           trigger1(j) = trigger(i,j)
        enddo

C	    Skip the trajectory if missing data are found for positions
        isok = 1
        do j=1,ntim
           if ( trainp1(j,4).lt.0. ) isok = 0
        enddo

c       Decide whether the trajectory is selected (handle SPECIAL commands)
        if ( isok.eq.1 ) then
           if (cmd(1).ne.0 ) then
              call select_tra (isok,cmd,ncmd,trainp1,trigger1,ntim,ncol)

           else
              call special    (isok,specialstr,trainp1,ntim,ncol,
     >                         vars,time,param,nparam)
           endif
	    endif

c       The trigger might be changed in the selection - copy it
        do j=1,ntim
            trigger(i,j) = trigger1(j)
        enddo

c       Set flag for selected trajectories
        if (isok.eq.1) then
           selflag(i)          = 1
           ntraout             = ntraout + 1
        else
           selflag(i)          = 0
        endif

      enddo

c     -------------------------------------------------------------- 
c     Write output trajectories
c     --------------------------------------------------------------

c     ------ Write output trajectories -----------------------------
      if ( outformat.eq.'trajectory' ) then


c        Define the trigger field if it is not yet defined
         if ( ( addtrigger.eq.'-trigger' ).and.(itrigger.eq.0) ) then
            ncol       = ncol + 1
            vars(ncol) = 'TRIGGER'
            itrigger   = ncol
         endif

c        Allocate memory for output trajectory
         allocate(traout(ntraout,ntim,ncol),stat=stat)
         if (stat.ne.0) stop '*** error allocating array apply   ***'

c        Set output trajectories
         j = 0
         do i=1,ntrainp
            if (selflag(i).eq.1) then
               j             = j + 1 
               traout(j,1:ntim,1:ncol) = trainp(i,1:ntim,1:ncol)
               if ( itrigger.ne.0 ) then
                  traout(j,1:ntim,itrigger) = real(trigger(i,1:ntim))
               endif
            endif
         enddo

c        Write trajectories
         call wopen_tra(fod,out_lslfile,ntraout,ntim,ncol,
     >        basetime,vars,outmode)
         call write_tra(fod,traout,ntraout,ntim,ncol,outmode)
         call close_tra(fod,outmode)
         
c     ------ Write index list -------------------------------------
      elseif ( outformat.eq.'index' ) then
         
         open(10,file=out_lslfile)
          do i=1,ntrainp
             if ( selflag(i).eq.1) write(10,*) i
          enddo
         close(10)

c     ------ Write boolean list -----------------------------------
      elseif ( outformat.eq.'boolean' ) then
         
         open(10,file=out_lslfile)
          do i=1,ntrainp
             write(10,'(i1)') selflag(i)
          enddo
         close(10)

c     ------ Write count -------------------------------------------
      elseif ( outformat.eq.'count' ) then
 
         open(10,file=out_lslfile)
          write(10,'(i7)') ntraout
         close(10)

c     ------ Write starting positions -----------------------------
      elseif ( outformat.eq.'startf' ) then

c        Allocate memory for output trajectory
         allocate(traout(ntraout,1,ncol),stat=stat)
         if (stat.ne.0) stop '*** error allocating array apply   ***'

c        Set output trajectories
         j = 0
         do i=1,ntrainp
            if (selflag(i).eq.1) then
               j             = j + 1 
               traout(j,1,:) = trainp(i,1,:)
            endif
         enddo

c        Write trajectories
         if (outmode.ne.-1) then
            
           call wopen_tra(fod,out_lslfile,ntraout,1,ncol,
     >                    basetime,vars,outmode)
           call write_tra(fod,traout,ntraout,1,ncol,outmode)
           call close_tra(fod,outmode)       

c        Output as a triple list (corresponding to <startf> file)
         else
         
           fid = 10
           open(fid,file=out_lslfile)
            do i=1,ntraout
              write(fid,'(3f10.3)') traout(i,1,2),   ! longitude
     >                            traout(i,1,3),   ! latitude
     >                            traout(i,1,4)    ! pressure
            enddo
           close(fid)

        endif
  
      endif

c     Write some status information, and end of program message
      print*  
      print*,'---- STATUS INFORMATION --------------------------------'
      print*
      print*,' # input trajectories  : ',ntrainp
      print*,' # output trajectories : ',ntraout
      print*
      print*,'              *** END OF PROGRAM SELECT ***'
      print*,'========================================================='

      stop
      
c     Exception handling
 100  stop 'select: First column in input trajectory must be <time>'
 101  stop 'select: Input trajectory file is empty'

      end

c     -------------------------------------------------------------- 
c     Dump the command list
c     --------------------------------------------------------------
      
      subroutine dumpcode(out,n,vars,nvars,times,ntimes)

c     Write the command list to screen. The command list is decoded
c     by call to <decode>

      implicit none

c     Declaration of subroutine parameters
      integer      n
      real         out(n)
      integer      nvars
      character*80 vars(nvars)
      integer      ntimes
      real         times(ntimes)

c     A single command
      character*80 cmd
      character*80 var,mode,strtim
      integer      nval
      integer      ntim
      integer      ivar,imode,icmd,itime

c     Auxiliary variables
      integer      i,j

c     Loop through the complete list
      i=0
 100  if (i.lt.n) then 

         write(*,*) '---------------------------------------'

c        Get command
         i=i+1
         icmd=nint(out(i))

c        Special handling of SPECIAL commands
         if ( icmd.eq.0 ) then

c           Write 'header' for SPECIAL command
            write(*,'(i5,f15.4,10x,a10)') i,out(i),'SPECIAL'

c           Write command string
            i = i + 1
            icmd = nint(out(i))
            write(*,'(i5,f15.4,10x,a10)') i,out(i),'LEN(CMD)'
            do j=1,icmd
               i    = i + 1
               ivar = nint(out(i))
               write(*,'(i5,f15.4,10x,a10)') i,out(i),char(ivar)
            enddo

c           Write parameters
            i    = i + 1
            nval = nint(out(i))
            write(*,'(i5,f15.4,10x,a10)') i,out(i),'#PARAMETER'
            do j=1,nval
               i=i+1
               if ( var.ne.'INPOLYGON' ) then
                  write(*,'(i5,f15.4)') i,out(i)
               else
                  write(*,'(i5,f15.4,a2)') i,out(i),char(nint(out(i)))
               endif
            enddo

c           Nothing else to do - exit
            goto 200
            
         endif

c        Set the command
         if (icmd.eq.  1) cmd='GT'
         if (icmd.eq.  2) cmd='LT'
         if (icmd.eq.  3) cmd='IN'
         if (icmd.eq.  4) cmd='OUT'
         if (icmd.eq.  5) cmd='EQ'
         if (icmd.eq.  6) cmd='TRUE'
         if (icmd.eq.  7) cmd='FALSE'
         if (icmd.eq.  8) cmd='ALL'
         if (icmd.eq.  9) cmd='ANY'
         if (icmd.eq. 10) cmd='NONE'
         if (icmd.eq. -1) cmd='BEGIN'
         if (icmd.eq. -2) cmd='END'
         if (icmd.eq. -3) cmd='AND'
         if (icmd.eq. -4) cmd='OR'
         write(*,'(i5,f15.4,10x,a10)') i,out(i),trim(cmd)
         if (icmd.lt.0) goto 100

c        Get variable
         i=i+1
         ivar=nint(out(i))
         if ( ivar.eq. -1 ) then
            var = 'DIST'
         elseif ( ivar.eq. -2 ) then
            var = 'DIST0'
         elseif ( ivar.eq. -3 ) then
            var = 'INPOLYGON'
         elseif ( ivar.eq. -4 ) then
            var = 'INBOX'
         elseif ( ivar.eq. -5 ) then
            var = 'INCIRCLE'
         elseif ( ivar.eq. -6 ) then
            var = 'INREGION'
         elseif ( ivar.eq. -7 ) then
            var = 'TRIGGER'
         elseif ( ivar.eq. -8 ) then
            var = 'VERT0'
         else
            var=vars(ivar)
         endif
         write(*,'(i5,f15.4,10x,a10)') i,out(i),trim(var)

c        Get variable mode
         i=i+1
         imode=nint(out(i))
         if (imode.eq.1 ) mode='VALUE'
         if (imode.eq.2 ) mode='MEAN'
         if (imode.eq.3 ) mode='MAX'
         if (imode.eq.4 ) mode='MIN'
         if (imode.eq.5 ) mode='VAR'
         if (imode.eq.6 ) mode='SUM'
         if (imode.eq.7 ) mode='CHANGE'
         if (imode.eq.8 ) mode='DIFF'
         if (imode.eq.9 ) mode='RANGE'
         if (imode.eq.10) mode='ABS'
         write(*,'(i5,f15.4,10x,a10)') i,out(i),trim(mode)

c        Get values
         i=i+1
         nval=nint(out(i))
         write(*,'(i5,f15.4,10x,a10)') i,out(i),'#PARAMETER'
         do j=1,nval
            i=i+1

            if ( var.ne.'INPOLYGON' ) then
               write(*,'(i5,f15.4)') i,out(i)
            else
               write(*,'(i5,f15.4,a2)') i,out(i),char(nint(out(i)))
            endif
         enddo

c        Get the number of times
         i=i+1
         ntim=nint(out(i))

c        the number of times is variable - depending on TRIGGER
         if ( ntim .eq. -993 ) then
            write(*,'(i5,f15.4,7x,7x,a15)') i,out(i),'TIMES @ TRIGGER'

c        Get the detailed list of times
         else
            write(*,'(i5,f15.4,7x,7x,a6)') i,out(i),'#TIMES'
            do j=1,ntim
               i=i+1
               write(*,'(i5,f15.4,f7.0)') i,out(i),times(nint(out(i)))
            enddo
         endif

c        Get time mode
         i=i+1
         itime=nint(out(i))
         if (itime.eq.1) strtim='ALL'
         if (itime.eq.2) strtim='ANY'         
         if (itime.eq.3) strtim='NONE'
         if (itime.lt.0) strtim='TRIGGER'

         if ( strtim.ne.'TRIGGER' ) then
            write(*,'(i5,f15.4,10x,a10)') i,out(i),trim(strtim)
         else
            write(*,'(i5,f15.4,10x,a10,a3,i3)') i,out(i),trim(strtim),
     >                                        ' ->',abs(itime)
         endif

         goto 100
         
      endif

c     Exit point
 200  continue

      write(*,*) '---------------------------------------'

      end
         

c     -------------------------------------------------------------- 
c     Read and decode a selection set
c     --------------------------------------------------------------

      subroutine decode(fid,out,n,vars,nvars,times,ntimes,regionf)

c     A selection file is opened with file id <fid> and transformed
c     into a set of commands applied to the trajectories. On input
c     <n> sets the maximum dimension of <out>, on output it gives the
c     total length of the command string. The output is a list of
c     commands with the following format:
c     
c        out(i)             = Command 
c        out(i+1)           = Column index of variable
c        out(i+2)           = Mode for variable
c        out(i+3)           = Number of parameters (n)
c        out(i+4:i+4+n)     = Parameters
c        out(i+5+n)         = Number of times 
c        out(i+6+n:i+6+n+m) = List of time indices (m)
c        out(i+7+n+m)       = Time mode
c     
c     For SPECIAL commands (to be coded in <special.f>) the format is:
c     
c        out(i)             = Length of command string (n)
c        out(i+1:i+n)       = Command string      
c        out(i+n+1)         = Number of parameters (m)
c        out(i+n+2:i+n+1+m) = List of parameters
c     
c     The following coding is used
c      
c        Command        Variable mode    Time mode
c        ----------     -------------    ---------
c        GT       1     VALUE    1       ALL      1
c        LT       2     MEAN     2       ANY      2
c        IN       3     MAX      3       NONE     3
c        OUT      4     MIN      4       TRIGGER -i (i the trigger index)
c        EQ       5     VAR      5
c        BEGIN   -1     SUM      6
c        END     -2     CHANGE   7 
c        AND     -3     DIFF     8
c        OR      -4     RANGE    9
c        TRUE     6     ABS     10
c        FALSE    7
c        SPECIAL  0
c        ALL      8 (TRIGGER)
c        ANY      9 (TRIGGER)
c        NONE    10 (TRIGGER)


c     Several "implicit variables" are supported - out(i+1):
c     
c        DIST      -1 : Path length of the trajectory
c        DIST0     -2 : Distance from starting position
c        INPOLYGON -3 : Specified polygon region
c        INBOX     -4 : Longitude/latitude rectangle
c        INCIRCLE  -5 : Within a specified radius 
c        INREGION  -6 : Within a specified rehion on the region file
c        TRIGGER   -7 : Trigger field
c        VERT0     -8 : Vertical distance from starting position
c         
c     For the special commands BEGIN, END, AND and OR, only one field
c     in <out> is used.

      implicit none

c     Declaration of subroutine parameters
      integer      fid
      integer      n
      real         out(n)
      integer      nvars
      character*80 vars(nvars)
      integer      ntimes
      real         times(ntimes)
      character*80 regionf

c     Numerical epsilon
      real         eps
      parameter    (eps=0.001)

c     A single command term
      character*80 cmd
      character*80 var,mode
      integer      nval
      real         val(1000)
      integer      ntim
      real         tim(1000)
      character*80 tmode

c     Specification of a polygon
      character*80 filename
      integer      pn                             ! Number of entries in lat/lon poly
      real         latpoly(1000)                  ! List of polygon latitudes
      real         lonpoly(1000)                  ! List of polygon longitudes
      real         loninpoly,latinpoly            ! Lon/lat inside polygon

c     Specification of a region
      real         xcorner(4)
      real         ycorner(4)
      integer      iregion
      character*80 string

c     Transformation to UPN (handling of logical operators)
      integer      nlogical
      integer      ilogical(n)
      real         tmp(n)
      integer      mlogical
      integer      isor

c     Auxiliary variables
      integer      i,j
      integer      j1,j2
      integer      flag(ntimes)
      integer      count
      integer      ok
      integer      itrigger,ttrigger

c     Common block for initialisation of polygon check
      real    tlonv(1000),vlat_c(1000),vlon_c(1000)
      real    xlat_c,xlon_c
      integer ibndry,nv_c
      common  /spolybndry/vlat_c,vlon_c,nv_c,xlat_c,xlon_c,tlonv,ibndry

c     ------ Decode single commands ---------------------

c     Reset the filename for polygons
      filename='nil'

c     Reset the counter for logical commands
      nlogical=0

c     Loop through all commands
 100  continue

c       Read next command (handle special cases)
        read(fid,*) cmd

c       Special handling of SPECIAL commands
        if ( cmd.eq.'SPECIAL' ) then

c          Set the flag for SPECIAL command
           n      = 1
           out(n) = 0.

c          Read the command string
           read(fid,*) var

c          Add the command string
           n      = n + 1
           out(n) = len_trim(var)
           do j=1,len_trim(var)
              n      = n + 1
              out(n) = ichar(var(j:j))
           enddo

c          Read the parameters
           read(fid,*) nval
           read(fid,*) (val(i),i=1,nval)

c          Add the parameters
           n = n + 1
           out(n)=real(nval)
           do i=1,nval
              n=n+1
              out(n)=val(i)
           enddo

c          Goto exit point - nothing more top do
           goto 350

        endif

c       Handle structure commands
        if ( cmd.eq.'BEGIN') then
           out(1)=-1.
           n=1
           nlogical=1
           ilogical(1)=n
           goto 200
        elseif ( cmd.eq.'AND'  ) then
           n=n+1
           out(n)=-3.
           nlogical=nlogical+1
           ilogical(nlogical)=n
           goto 200
        elseif ( cmd.eq.'OR'   ) then
           n=n+1
           out(n)=-4.
           nlogical=nlogical+1
           ilogical(nlogical)=n
           goto 200
        elseif ( cmd.eq.'END'  ) then
           n=n+1
           out(n)=-2.
           nlogical=nlogical+1
           ilogical(nlogical)=n
           goto 300
        endif

c       Read other fields associated with the command
        read(fid,*) var,mode

c       Read parameter
        if ( var.eq.'INPOLYGON' ) then 
           read(fid,*) nval
           read(fid,*) filename
           filename = trim(filename)
        else
           read(fid,*) nval
           read(fid,*) (val(i),i=1,nval)
        endif

c       Read times (on request, change to special trigger times)
        read(fid,*) ntim
        if ( ntim.eq.-993 ) then
           ttrigger = 1
           ntim     = 1
           tim(1)   = -999.
        else
           ttrigger = 0
           read(fid,*) (tim(i),i=1,ntim)
        endif
        read(fid,*) tmode

c       Bring CAPITAL "TIME,LAT,LON,P" into "time,lat,lon,p"
        if (var.eq.'TIME') var='time'
        if (var.eq.'LAT' ) var='lat'
        if (var.eq.'LON' ) var='lon'
        if (var.eq.'P'   ) var='p'

c       If the time mode is 'TRIGGER', all times of a trajectory
c       must be considered
        if ( tmode.eq.'TRIGGER' ) then
           itrigger = nint(tim(1))
           ntim     = 1
           tim(1)   = -999.
        endif

c       Special times: transform into real time
        do i=1,ntim
           if ( abs(tim(i)+996.).lt.eps ) tim(i)=times(1)
           if ( abs(tim(i)+995.).lt.eps ) tim(i)=times(ntimes)
        enddo

c       Check whether times are valid 
        ok=0
        do i=1,ntim
           if ( (abs(tim(i)+994.).gt.eps).and.
     >          (abs(tim(i)+999.).gt.eps) )
     >     then
              do j=1,ntimes
                 if ( abs(tim(i)-times(j)).lt.eps ) then
                    ok=ok+1
                 endif
              enddo
           else
              ok=ok+1
           endif
        enddo
        if (ok.ne.ntim) goto 400
        
c       Select all times which are included in the criterion
        do i=1,ntimes
           flag(i)=0
        enddo
        i=1
 150    if (i.le.ntim) then
           
c          A list of times
           if ( (abs(tim(i)+994.).lt.eps) ) then
              j1=0
              do j=1,ntimes
                 if ( abs(tim(i-1)-times(j)).lt.eps ) then
                    j1=j
                 endif
              enddo
              j2=0
              do j=1,ntimes
                 if ( abs(tim(i+1)-times(j)).lt.eps ) then
                    j2=j
                 endif
              enddo   
              if ( (j1.eq.0).or.(j2.eq.0) ) goto 400
              do j=j1,j2
                 flag(j)=i
              enddo
              i=i+1

c          Explicitly given time value
           else
              do j=1,ntimes
                 if ( abs(tim(i)-times(j)).lt.eps ) then
                    flag(j)=i
                 endif
              enddo
              
           endif

           i=i+1
           goto 150

        endif

c       Write command identifier
        n=n+1
        if (cmd.eq.'GT'    ) out(n)= 1.
        if (cmd.eq.'LT'    ) out(n)= 2.
        if (cmd.eq.'IN'    ) out(n)= 3.
        if (cmd.eq.'OUT'   ) out(n)= 4.
        if (cmd.eq.'EQ'    ) out(n)= 5.
        if (cmd.eq.'TRUE'  ) out(n)= 6.
        if (cmd.eq.'FALSE' ) out(n)= 7.
        if (cmd.eq.'ALL  ' ) out(n)= 8.
        if (cmd.eq.'ANY  ' ) out(n)= 9.
        if (cmd.eq.'NONE ' ) out(n)=10.

c       Write index for variable - force implicit trigger
        ok=0
        do j=1,nvars
           if (vars(j).eq.var) ok=j
        enddo

        if (var.eq.'TRIGGER') ok = 0
        
        if (ok.eq.0) then
           if (var.eq.'DIST') then
              ok = -1
           elseif (var.eq.'DIST0') then
              ok = -2
           elseif (var.eq.'INPOLYGON') then
              ok = -3
           elseif (var.eq.'INBOX') then
              ok = -4
           elseif (var.eq.'INCIRCLE') then
              ok = -5
           elseif (var.eq.'INREGION') then
              ok = -6
           elseif (var.eq.'TRIGGER') then
              ok = -7
           elseif (var.eq.'VERT0') then
              ok = -8
           else
             goto 400
          endif
        endif
        n=n+1
        out(n)=real(ok)

c       Write mode for variable
        ok=0
        if (mode.eq.'VALUE'  ) ok=1
        if (mode.eq.'MEAN'   ) ok=2
        if (mode.eq.'MAX'    ) ok=3
        if (mode.eq.'MIN'    ) ok=4
        if (mode.eq.'VAR'    ) ok=5
        if (mode.eq.'SUM'    ) ok=6
        if (mode.eq.'CHANGE' ) ok=7
        if (mode.eq.'DIFF'   ) ok=8
        if (mode.eq.'RANGE'  ) ok=9
        if (mode.eq.'ABS'    ) ok=10
        if (ok.eq.0) goto 400
        n=n+1
        out(n)=real(ok)

c       Write the parameter values: INPOLYGON 
        if ( var.eq.'INPOLYGON' ) then

           n      = n+1
           out(n) = len_trim(filename)
           do j=1,len_trim(filename)
              n      = n + 1
              out(n) = ichar(filename(j:j))
           enddo

c       Write parameter value: INREGION 
        elseif ( var.eq.'INREGION' ) then

           iregion = nint(val(1))

           open(fid+1,file=regionf)          

 50        read(fid+1,*,end=51) string

           if ( string(1:1).ne.'#' ) then
              call regionsplit(string,i,xcorner,ycorner)
              if ( i.eq.iregion ) goto 52
           endif

           goto 50
           
 51        close(fid+1)        

           print*,' ERROR: region ',iregion,' not found on ',
     >                                         trim(regionf)
           stop

 52        continue
           
           n      = n + 1
           out(n) = 8                   ! Number of parameters
           do i=1,4
              n      = n + 1
              out(n) = xcorner(i)
           enddo
           do i=1,4
              n      = n + 1
              out(n) = ycorner(i)
           enddo

c       Write parameter values: all other cases
        else
           n=n+1
           out(n)=real(nval)
           do i=1,nval
              n=n+1
              out(n)=val(i)
           enddo
        endif

c       Special time handling: only trigger times are cosidered
        if ( ttrigger.eq.1 ) then
           n = n+1
           out(n)=-993.

c       All times are selected
        elseif ( abs(tim(1)+999.).lt.eps ) then
           n=n+1
           out(n)=real(ntimes)
           do i=1,ntimes
              n=n+1
              out(n)=real(i)
           enddo

c       A selection of times is given 
        else
           count=0
           do i=1,ntimes
              if ( flag(i).ne.0 ) then
                 count=count+1
              endif
           enddo
           n=n+1
           out(n)=real(count)
           do i=1,count
             do j=1,ntimes
               if (flag(j).eq.i) then
                 n      = n + 1
                 out(n) = real(j)
               endif
             enddo
           enddo

        endif

c       Write the time mode
        if ( tmode.eq.'ALL') then
           n=n+1
           out(n)=1.
        elseif ( tmode.eq.'ANY') then
           n=n+1
           out(n)=2.
        elseif ( tmode.eq.'NONE') then
           n=n+1
           out(n)=3.
        elseif ( tmode.eq.'TRIGGER') then
           n=n+1
           out(n)=-real(itrigger)
        endif
           
c     End loop: handle single command
 200  continue
      goto 100

c     End loop: loop over all commands 
 300  continue

c     ------ Read polygon file, if requested -----------
      if ( filename.ne.'nil' ) then

        print*
        print*,
     >     '---- POLYGON --------------------------------------------'

           print*
           print*,'Filename = ',trim(filename)
           print*
         
c        Read list of polygon coordinates from file
         pn = 0
         open(fid+1,file=filename)
           read(fid+1,*) loninpoly,latinpoly
           print*,'Inside (lon/lat) =',loninpoly,latinpoly
           print*
 510       continue
              pn = pn + 1
              read(fid+1,*,end=511) lonpoly(pn),
     >                              latpoly(pn)

              print*,pn,lonpoly(pn),latpoly(pn)
              
              goto 510
 511       continue
           pn = pn - 1
         close(fid+1)

c        Define the polygon boundaries
         call DefSPolyBndry(latpoly,lonpoly,pn,latinpoly,loninpoly)

      endif

c     ------ Transform to UPN --------------------------

c     Check whether logical commands are ok
      mlogical=nint(out(ilogical(1)))
      if ( mlogical.ne.-1) goto 400
      mlogical=nint(out(ilogical(nlogical)))
      if ( mlogical.ne.-2) goto 400

c     No transformation necessary if only one command
      if (nlogical.eq.2) goto 350
      
c     Copy the output to temporary list
      do i=1,n
         tmp(i)=out(i)
      enddo

c     Set BEGIN statement
      n=1
      out(n)=-1.

c     Reorder commands and transform to UPN
      isor=0
      do i=1,nlogical-1

c        Get the logical command
         mlogical=nint(out(ilogical(i)))

c        Connecting OR
         if (mlogical.eq.-4) then
            if (isor.eq.1) then
               n=n+1
               out(n)=-4.
            else
               isor=1
            endif
         endif

c        Put the command onto the stack
         do j=ilogical(i)+1,ilogical(i+1)-1
            n=n+1
            out(n)=tmp(j)
         enddo

c        Connecting AND operator
         if ( mlogical.eq.-3 ) then
            n=n+1
            out(n)=-3.
         endif
         
      enddo
      
c     Set final connecting OR
      if (isor.eq.1) then
         n=n+1
         out(n)=-4.
      endif

c     Set END statement
      n=n+1
      out(n)=-2.

c     ------ Exit point ---------------------------------

 350  continue
      return

c     ----- Exception handling --------------------------
 400  print*,'Invalid selection criterion... Stop'
      stop

      end


c     -------------------------------------------------------------- 
c     Decide whether a trajectory is selected or not
c     --------------------------------------------------------------

      subroutine select_tra (select,cmd,n,tra,trigger,ntim,ncol)

c     Decide whether a single trajectory is selected (<select=1>) or
c     is not selected <select=0> according to the selection criterion
c     given in <cmd(ncmd)>. The selection criterion <cmd(ncmd)> is 
c     returned from the call to the subroutine <decode>. The trajectory
c     is given in <tra(ntim,ncol)> where <ntim> is the number of times
c     and <ncol> is the number of columns.
c     
c     Important note: the structure of <tra(ntim,ncol)> must match to the
c     call parameter <vars,nvars,times,ntimes> in subroutine <decode>.

      implicit none

c     Declaration of subroutine parameters
      integer   select
      integer   n
      real      cmd(n)
      integer   ntim,ncol
      real      tra(ntim,ncol)
      integer   trigger(ntim)

c     Numerical epsilon (for test of equality)
      real      eps
      parameter (eps=0.000001)

c     A single command and the associated field
      integer   icmd,ivar,imode,itime,nsel,nval
      integer   time(ntim)
      real      param(100)
      real      var   (ntim)
      integer   intvar(ntim)

c     Boolean values for a single time, a single command and build-up
      integer   stack(100)
      integer   nstack
      integer   istrue(ntim)
      integer   decision

c     Auxiliary variables
      integer   i,j,k
      real      tmp,mea
      integer   istack1,istack2
      real      lat0,lon0,lat1,lon1
      real      length(ntim)
      integer   flag
      real      dist
      real      xcorner(4),ycorner(4)
      integer   iparam
      character ch
      real      varmin,varmax
      real      lev0,lev1

c     Common block for initialisation of polygon check
      real    tlonv(1000),vlat_c(1000),vlon_c(1000)
      real    xlat_c,xlon_c
      integer ibndry,nv_c
      common  /spolybndry/vlat_c,vlon_c,nv_c,xlat_c,xlon_c,tlonv,ibndry

c     Externals
      real      sdis           ! Spherical distance
      external  sdis
      integer   inregion       ! Boolean flag for regions
      external  inregion

c     Reset the decision stack (with locical values)
      nstack=0

c     Loop through the complete command list
      i=0
 100  if (i.lt.n) then  

c        --- Get the command -------------------------------
         i=i+1
         icmd=nint(cmd(i))

c        --- Handle structural commands (BEGIN, END, AND, OR)

c        Handle BEGIN statement
         if ( icmd.eq.-1) then
            nstack=0
            goto 200
         endif

c        Handle END statement
         if (icmd.eq.-2) then
            goto 300
         endif

c        Handle AND statement
         if (icmd.eq.-3) then
            istack1=stack(nstack)
            nstack=nstack-1
            istack2=stack(nstack)
            if ((istack1.eq.1).and.(istack2.eq.1)) then
               stack(nstack)=1
            else
               stack(nstack)=0
            endif
            goto 200
         endif
            
c        Handle OR statement
         if (icmd.eq.-4) then
            istack1=stack(nstack)
            nstack=nstack-1
            istack2=stack(nstack)
            if ((istack1.eq.1).or.(istack2.eq.1)) then
               stack(nstack)=1
            else
               stack(nstack)=0
            endif
            goto 200
         endif

c        --- Get all command details (parameters, modes, times)

c        Get variable (<ivar> gets the column index in <tra>)
         i=i+1
         ivar=nint(cmd(i))
         
c        Get variable mode 
         i=i+1
         imode=nint(cmd(i))

c        Get parameter values
         i=i+1
         nval=nint(cmd(i))
         do j=1,nval
            i=i+1
            param(j)=cmd(i)
         enddo
         
c        Get times (<time(j)> gets the row indices of <tra>)
         i=i+1
         nsel=nint(cmd(i))
         if ( nsel .eq. -993 ) then
            nsel = 0
            do k=1,ntim
              if ( trigger(k).ne.0 ) then
                nsel       = nsel + 1
                time(nsel) = k
              endif
           enddo
         else
           do j=1,nsel
              i=i+1
              time(j)=nint(cmd(i))
           enddo
         endif

c        If no times are selected, exit with non-select status
         if ( nsel.eq.0 ) then
             stack(1) = 0
             goto 300
         endif

c        Get time mode
         i=i+1
         itime=nint(cmd(i))

c        --- Prepare field values for analysis -----------

c        Implicit variable: DIST
         if ( ivar.eq. -1 ) then
            length(1) = 0.
            do j=2,ntim
               lon0      = tra(j-1,2)
               lat0      = tra(j-1,3)
               lon1      = tra(j  ,2)
               lat1      = tra(j  ,3)
               length(j) = length(j-1) + sdis(lon0,lat0,lon1,lat1)
            enddo
            do j=1,nsel
               var(j) = length( time(j) )
            enddo

               
c        Implict variable: DIST0
         elseif ( ivar.eq. -2 ) then
            do j=1,nsel
               lon0   = tra(1      ,2)
               lat0   = tra(1      ,3)
               lon1   = tra(time(j),2)
               lat1   = tra(time(j),3)
               var(j) = sdis(lon0,lat0,lon1,lat1)
            enddo

c        Implict variable: INPOLYGON
         elseif ( ivar.eq. -3 ) then
            do j=1,nsel
               lon1   = tra(time(j),2)
               lat1   = tra(time(j),3)               
               call LctPtRelBndry(lat1,lon1,flag)
               if ( (flag.eq.1).or.(flag.eq.2) ) then
                  var(j) = 1.
               else
                  var(j) = 0.
               endif
            enddo

c        Implict variable: INBOX 
         elseif ( ivar.eq. -4 ) then
            do j=1,nsel
               lon1   = tra(time(j),2)
               lat1   = tra(time(j),3)               
               if ( ( lon1.ge.param(1) ).and.       ! lonmin
     >              ( lon1.le.param(2) ).and.       ! lonmax
     >              ( lat1.ge.param(3) ).and.       ! latmin
     >              ( lat1.le.param(4) ) )          ! latmax
     >         then
                  var(j) = 1
               else
                  var(j) = 0
               endif
            enddo

c        Implict variable: INCIRCLE (lonc=param(1),latc=param(2),radius=param(3))
         elseif ( ivar.eq. -5 ) then
            do j=1,nsel
               lon1   = tra(time(j),2)
               lat1   = tra(time(j),3)               
               
               dist = sdis( lon1,lat1,param(1),param(2) ) 

               if ( dist.le.param(3) ) then
                  var(j) = 1
               else
                  var(j) = 0
               endif
            enddo

c        Implict variable: INREGION (xcorner=param(1..4),ycorner=param(5..8) )
         elseif ( ivar.eq.-6 ) then
            
            do j=1,4
               xcorner(j) = param(j  )
               ycorner(j) = param(j+4)
            enddo

            do j=1,nsel
               lon1   = tra(time(j),2)
               lat1   = tra(time(j),3)  
               var(j) = inregion (lon1,lat1,xcorner,ycorner)               
            enddo

c        Implict variable: TRIGGER
         elseif ( ivar.eq. -7 ) then
            do j=1,nsel
               intvar(j) = trigger( time(j) ) 
            enddo

c        Implicit variable: VERT0
         elseif ( ivar.eq. -8 ) then
            do j=1,nsel
               lev0   = tra(1      ,4)
               lev1   = tra(time(j),4)
               var(j) = lev0 - lev1
            enddo

c        Explicit variable (column index <ivar>)
         else
            do j=1,nsel
               var(j) = tra(time(j),ivar)
            enddo

         endif

c        Take MEAN of the variable (mean of selected times)
         if (imode.eq.2) then
            tmp=0.
            do j=1,nsel
               tmp=tmp+var(j)
            enddo
            var(1)=tmp/real(nsel)
            nsel=1

c        Take MAX of the variable (maximum of selected times)
         elseif (imode.eq.3) then
            tmp=var(1)
            do j=2,nsel
               if (var(j).gt.tmp) tmp=var(j)
            enddo
            var(1)=tmp
            nsel=1
            
c        Take MIN of the variable (minimum of selected times)
         elseif (imode.eq.4) then
            tmp=var(1)
            do j=2,nsel
               if (var(j).lt.tmp) tmp=var(j)
            enddo
            var(1)=tmp
            nsel=1

c        Take VAR of the variable (variance over all selected times)
         elseif (imode.eq.5) then
            tmp=0.
            do j=1,nsel
               tmp=tmp+var(j)
            enddo
            mea=tmp/real(nsel)
            do j=1,nsel
              tmp=tmp+(var(j)-mea)**2
            enddo
            var(1)=1./real(nsel-1)*tmp
            nsel=1

c        Take SUM of the variable (sum over all selected times)
         elseif (imode.eq.6) then
            tmp=0.
            do j=1,nsel
               tmp=tmp+var(j)
            enddo
            var(1)=tmp
            nsel=1
            
c        Take CHANGE of the variable (absolute difference between first and last time)
         elseif (imode.eq.7) then
            var(1)=abs(var(1)-var(nsel))
            nsel=1

c        Take DIFF of the variable (first minus last time)
         elseif (imode.eq.8) then
            var(1)=var(1)-var(nsel)
            nsel=1

c        Take RANGE of the variable
         elseif (imode.eq.9) then
            varmax=var(1)
            varmin=var(1)
            do j=2,nsel
               if (var(j).gt.varmax) varmax=var(j)
               if (var(j).lt.varmin) varmin=var(j)
            enddo
            var(1) = varmax - varmin
            nsel=1

c        Take the absolute value of the variable
         elseif (imode.eq.10) then
            do j=1,nsel
               var(j)=abs( var(j) )
            enddo
         endif

c        --- Apply the operators to the single values ---

         do j=1,nsel

c           GT
            if (icmd.eq.1) then
               if (var(j).gt.param(1)) then
                  istrue(j)=1
               else
                  istrue(j)=0
               endif
               
c           LT
            elseif (icmd.eq.2) then
               if (var(j).lt.param(1)) then
                  istrue(j)=1
               else
                  istrue(j)=0
               endif

c           IN
            elseif (icmd.eq.3) then
               if ( (var(j).gt.param(1)).and.
     >              (var(j).lt.param(2)) ) 
     >         then
                  istrue(j)=1
               else
                  istrue(j)=0
               endif

c           OUT
            elseif (icmd.eq.4) then
               if ( (var(j).lt.param(1)).or.
     >              (var(j).gt.param(2)) ) 
     >         then
                  istrue(j)=1
               else
                  istrue(j)=0
               endif

c           EQ
            elseif (icmd.eq.5) then
               if (abs(var(j)-param(1)).lt.eps) then
                  istrue(j)=1
               else
                  istrue(j)=0
               endif

c           TRUE
            elseif (icmd.eq.6) then
               if (abs(var(j)).lt.eps) then
                  istrue(j)=0
               else
                  istrue(j)=1
               endif
 
c           FALSE
            elseif (icmd.eq.7) then
               if (abs(var(j)).lt.eps) then
                  istrue(j)=1
               else
                  istrue(j)=0
               endif

c           ALL
            elseif (icmd.eq.8) then
               istrue(j) = 1
               do k=1,nval
                  iparam = nint(param(k))-1
                  if (btest(intvar(j),iparam).eqv..false.) then
                     istrue(j) = 0
                  endif
               enddo
               
c           ANY
            elseif (icmd.eq.9) then
               istrue(j) = 0
               do k=1,nval
                  iparam = nint(param(k))-1
                  if (btest(intvar(j),iparam).eqv..true.) then
                     istrue(j) = 1
                  endif
               enddo

c           NONE
            elseif (icmd.eq.10) then
               istrue(j) = 1
               do k=1,nval
                  iparam = nint(param(k))-1
                  if (btest(intvar(j),iparam).eqv..true.) then
                     istrue(j) = 0
                  endif
               enddo
   
            endif

         enddo

c        --- Determine the overall boolean value ----------
            
c        ALL
         if (itime.eq.1) then
            decision=1
            do j=1,nsel
               if (istrue(j).eq.0) then
                  decision=0
                  goto 110
               endif
            enddo
 110        continue

c        ANY
         elseif (itime.eq.2) then
            decision=0
            do j=1,nsel
               if (istrue(j).eq.1) then
                  decision=1
                  goto 120
               endif
            enddo
 120        continue
           
c        NONE
         elseif (itime.eq.3) then
            decision=1
            do j=1,nsel
               if (istrue(j).eq.1) then
                  decision=0
                  goto 130
               endif
            enddo
 130        continue

c        TRIGGER
         elseif (itime.lt.0) then
            decision=1
            do j=1,nsel
               if (istrue(j).eq.1) then
                  trigger(j) = ior( trigger(j), 2**(abs(itime)-1) )
               endif
            enddo
            
         endif

c        --- Put the new boolean value onto the stack

         nstack=nstack+1
         stack(nstack)=decision

c        Exit point for loop
 200     continue
         goto 100

      endif

c     Return the decision (selected or non-selected)
 300  continue

      select=stack(1)

      end
      

c     --------------------------------------------------------------------------
c     Split a region string and get corners of the domain
c     --------------------------------------------------------------------------

      subroutine regionsplit(string,iregion,xcorner,ycorner)

c     The region string comes either as <lonw,lone,lats,latn> or as <lon1,lat1,
c     lon2,lat2,lon3,lat3,lon4,lat4>: split it into ints components and get the
c     four coordinates for the region
      
      implicit none

c     Declaration of subroutine parameters
      character*80    string
      real            xcorner(4),ycorner(4)
      integer         iregion

c     Local variables
      integer         i,n
      integer         il,ir
      real            subfloat (80)
      integer         stat
      integer         len

c     ------- Split the string
      i    = 1
      n    = 0
      stat = 0
      il   = 1
      len  = len_trim(string)

 100  continue

c     Find start of a substring
      do while ( stat.eq.0 )
         if ( string(i:i).ne.' ' ) then
            stat = 1
            il   = i
         else
            i = i + 1
         endif
      enddo

c     Find end of substring
      do while ( stat.eq.1 )         
         if ( ( string(i:i).eq.' ' ) .or. ( i.eq.len ) ) then
            stat = 2
            ir   = i
         else
            i    = i + 1
         endif
      enddo

c     Convert the substring into a number
      if ( stat.eq.2 ) then
         n = n + 1
         read(string(il:ir),*) subfloat(n)
         stat = 0
      endif

      if ( i.lt.len ) goto 100


c     -------- Get the region number
      
      iregion = nint(subfloat(1))

c     -------- Get the corners of the region
      
      if ( n.eq.5 ) then     ! lonw(2),lone(3),lats(4),latn(5)

         xcorner(1) = subfloat(2)
         ycorner(1) = subfloat(4)

         xcorner(2) = subfloat(3)
         ycorner(2) = subfloat(4)
 
         xcorner(3) = subfloat(3)
         ycorner(3) = subfloat(5)
         
         xcorner(4) = subfloat(2)
         ycorner(4) = subfloat(5)
        
      elseif ( n.eq.9 ) then     ! lon1,lat1,lon2,lat2,lon3,lon4,lat4

         xcorner(1) = subfloat(2)
         ycorner(1) = subfloat(3)

         xcorner(2) = subfloat(4)
         ycorner(2) = subfloat(5)

         xcorner(3) = subfloat(6)
         ycorner(3) = subfloat(7)
         
         xcorner(4) = subfloat(8)
         ycorner(4) = subfloat(9)
 
      else
         
         print*,' ERROR: invalid region specification '
         print*,'     ',trim(string)
         stop
         
      endif
         

      end

c     --------------------------------------------------------------------------
c     Decide whether lat/lon point is in or out of region
c     --------------------------------------------------------------------------
      
      integer function inregion (lon,lat,xcorner,ycorner)
      
c     Decide whether point (lon/lat) is in the region specified by <xcorner(1..4),
c     ycorner(1..4).
      
      implicit none
      
c     Declaration of subroutine parameters
      real    lon,lat
      real    xcorner(4),ycorner(4)

c     Local variables
      integer flag
      real    xmin,xmax,ymin,ymax
      integer i

c     Reset the flag
      flag = 0

c     Set some boundaries
      xmax = xcorner(1)
      xmin = xcorner(1)
      ymax = ycorner(1)
      ymin = ycorner(1)
      do i=2,4
        if (xcorner(i).lt.xmin) xmin = xcorner(i)
        if (xcorner(i).gt.xmax) xmax = xcorner(i)
        if (ycorner(i).lt.ymin) ymin = ycorner(i)
        if (ycorner(i).gt.ymax) ymax = ycorner(i)
      enddo

c     Do the tests - set flag=1 if all tests pased
      if (lon.lt.xmin) goto 970
      if (lon.gt.xmax) goto 970
      if (lat.lt.ymin) goto 970
      if (lat.gt.ymax) goto 970
      
      if ((lon-xcorner(1))*(ycorner(2)-ycorner(1))-
     >    (lat-ycorner(1))*(xcorner(2)-xcorner(1)).gt.0.) goto 970
      if ((lon-xcorner(2))*(ycorner(3)-ycorner(2))-
     >    (lat-ycorner(2))*(xcorner(3)-xcorner(2)).gt.0.) goto 970
      if ((lon-xcorner(3))*(ycorner(4)-ycorner(3))-
     >    (lat-ycorner(3))*(xcorner(4)-xcorner(3)).gt.0.) goto 970
      if ((lon-xcorner(4))*(ycorner(1)-ycorner(4))-
     >    (lat-ycorner(4))*(xcorner(1)-xcorner(4)).gt.0.) goto 970

      flag = 1

c     Return the value
 970  continue
      
      inregion = flag
      
      return
      
      end

   
c     --------------------------------------------------------------------------
c     Spherical distance between lat/lon points                                                       
c     --------------------------------------------------------------------------

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


c     ****************************************************************
c     * Given some spherical polygon S and some point X known to be  *
c     * located inside S, these routines will determine if an arbit- *
c     * -rary point P lies inside S, outside S, or on its boundary.  *
c     * The calling program must first call DefSPolyBndry to define  *
c     * the boundary of S and the point X. Any subsequent call to    *
c     * subroutine LctPtRelBndry will determine if some point P lies *
c     * inside or outside S, or on its boundary. (Usually            *
c     * DefSPolyBndry is called once, then LctPrRelBndry is called   *
c     * many times).                                                 *
c     *                                                              * 
c     * REFERENCE:            Bevis, M. and Chatelain, J.-L. (1989)  * 
c     *                       Maflaematical Geology, vol 21.         *
c     * VERSION 1.0                                                  *
c     ****************************************************************

      Subroutine DefSPolyBndry(vlat,vlon,nv,xlat, xlon)

c     ****************************************************************
c     * This mmn entry point is used m define ~e spheric~ polygon S  *
c     * and the point X.                                             *
c     * ARGUMENTS:                                                   *
c     * vlat,vlon (sent) ... vectors containing the latitude and     * 
c     *                      longitude of each vertex of the         *
c     *                      spherical polygon S. The ith.vertex is  *
c     *                      located at [vlat(i),vlon(i)].           *
c     * nv        (sent) ... the number of vertices and sides in the *
c     *                      spherical polygon S                     *
c     * xlat,xlon (sent) ... latitude and longitude of some point X  *
c     *                      located inside S. X must not be located *
c     *                      on any great circle that includes two   *
c     *                      vertices of S.                          *
c     *                                                              *
c     * UNITS AND SIGN CONVENTION:                                   *
c     *  Latitudes and longitudes are specified in degrees.          *
c     *  Latitudes are positive to the north and negative to the     *
c     *  south.                                                      *
c     *  Longitudes are positive to the east and negative to the     *
c     *  west.                                                       *
c     *                                                              * 
c     * VERTEX ENUMERATION:                                          * 
c     * The vertices of S should be numbered sequentially around the *
c     * border of the spherical polygon. Vertex 1 lies between vertex*
c     * nv and vertex 2. Neighbouring vertices must be seperated by  *
c     * less than 180 degrees. (In order to generate a polygon side  *
c     * whose arc length equals or exceeds 180 degrees simply        *
c     * introduce an additional (pseudo)vertex). Having chosen       *
c     * vertex 1, the user may number the remaining vertices in      *
c     * either direction. However if the user wishes to use the      *
c     * subroutine SPA to determine the area of the polygon S (Bevis *
c     * & Cambareri, 1987, Math. Geol., v.19, p. 335-346) then he or *
c     * she must follow the convention whereby in moving around the  *
c     * polygon border in the direction of increasing vertex number  *
c     * clockwise bends occur at salient vertices. A vertex is       *
c     * salient if the interior angle is less than 180 degrees.      *
c     * (In the case of a convex polygon this convention implies     *
c     * that vertices are numbered in clockwise sequence).           *
c     ****************************************************************

      implicit none
      
      integer mxnv,nv

c     ----------------------------------------------------------------
c     Edit next statement to increase maximum number of vertices that 
c     may be used to define the spherical polygon S               
c     The value of parameter mxnv in subroutine LctPtRelBndry must match
c     that of parameter mxnv in this subroutine, as assigned above.
c     ----------------------------------------------------------------
      parameter (mxnv=2000)

      real  vlat(nv),vlon(nv),xlat,xlon,dellon
      real  tlonv(mxnv),vlat_c(mxnv),vlon_c(mxnv),xlat_c,xlon_c
      integer i,ibndry,nv_c,ip
 
      data ibndry/0/
      
      common/spolybndry/vlat_c,vlon_c,nv_c,xlat_c,xlon_c,tlonv,ibndry

      if (nv.gt.mxnv) then
         print *,'nv exceeds maximum allowed value'
         print *,'adjust parameter mxnv in subroutine DefSPolyBndry'
         stop
      endif

      ibndry=1                  ! boundary defined at least once (flag)
      nv_c=nv                   ! copy for named common
      xlat_c=xlat               ! . . . .
      xlon_c=xlon               !

      do i=1,nv
         vlat_c(i)=vlat(i)      ! "
         vlon_c(i)=vlon(i)      !

         call TrnsfmLon(xlat,xlon,vlat(i),vlon(i),tlonv(i))

         if (i.gt.1) then
            ip=i-1
         else
            ip=nv
         endif
         
         if ((vlat(i).eq.vlat(ip)).and.(vlon(i).eq.vlon(ip))) then
            print *,'DefSPolyBndry detects user error:'
            print *,'vertices ',i,' and ',ip,' are not distinct'
            print*,'lat ',i,ip,vlat(i),vlat(ip)
            print*,'lon ',i,ip,vlon(i),vlon(ip)            
            stop
         endif

         if (tlonv(i).eq.tlonv(ip)) then
            print *,'DefSPolyBndry detects user error:'
            print *,'vertices ',i,' & ',ip,' on same gt. circle as X'
            stop
         endif

         if (vlat(i).eq.(-vlat(ip))) then
            dellon=vlon(i)-vlon(ip)
            if (dellon.gt.+180.) dellon=dellon-360.
            if (dellon.lt.-180.) dellon=dellon-360.
            if ((dellon.eq.+180.0).or.(dellon.eq.-180.0)) then
               print *,'DefSPolyBndry detects user error:'
               print *,'vertices ',i,' and ',ip,' are antipodal'
               stop
            endif
         endif
      enddo

      return
      
      end


c     ****************************************************************
 
      Subroutine LctPtRelBndry(plat,plon,location)

c     ****************************************************************

c     ****************************************************************
c     * This routine is used to see if some point P is located       *
c     * inside, outside or on the boundary of the spherical polygon  *
c     * S previously defined by a call to subroutine DefSPolyBndry.  *
c     * There is a single restriction on point P: it must not be     *
c     * antipodal to the point X defined in the call to DefSPolyBndry*
c     * (ie.P and X cannot be seperated by exactly 180 degrees).     *
c     * ARGUMENTS:                                                   *  
c     * plat,plon (sent)... the latitude and longitude of point P    *
c     * location (returned)... specifies the location of P:          *
c     *                        location=0 implies P is outside of S  *
c     *                        location=1 implies P is inside of S   *
c     *                        location=2 implies P on boundary of S *
c     *                        location=3 implies user error (P is   *
c     *                                     antipodal to X)          *
c     * UNFfS AND SIGN CONVENTION:                                   * 
c     *  Latitudes and longitudes are specified in degrees.          *
c     *  Latitudes are positive to the north and negative to the     *
c     *  south.                                                      *    
c     *  Longitudes are positive to the east and negative to the     *
c     *  west.                                                       *
c     ****************************************************************
      
      implicit none
      
      integer mxnv

c     ----------------------------------------------------------------
c     The statement below must match that in subroutine DefSPolyBndry
c     ----------------------------------------------------------------

      parameter (mxnv=2000)

      real tlonv(mxnv),vlat_c(mxnv),vlon_c(mxnv),xlat_c,xlon_c
      real plat,plon,vAlat,vAlon,vBlat,vBlon,tlonA,tlonB,tlonP
      real tlon_X,tlon_P,tlon_B,dellon
      integer i,ibndry,nv_c,location,icross,ibrngAB,ibrngAP,ibrngPB
      integer ibrng_BX,ibrng_BP,istrike

      common/spolybndry/vlat_c,vlon_c,nv_c,xlat_c,xlon_c,tlonv,ibndry

      if (ibndry.eq.0) then     ! user has never defined the bndry
         print*,'Subroutine LctPtRelBndry detects user error:'
         print*,'Subroutine DefSPolyBndry must be called before'
         print*,'subroutine LctPtRelBndry can be called'
         stop
      endif

      if (plat.eq.(-xlat_c)) then
         dellon=plon-xlon_c
         if (dellon.lt.(-180.)) dellon=dellon+360.
         if (dellon.gt.+180.) dellon=dellon-360.
         if ((dellon.eq.+180.0).or.(dellon.eq.-180.)) then
            print*,'Warning: LctPtRelBndry detects case P antipodal
     >           to X'
            print*,'location of P relative to S is undetermined'
            location=3
            return
         endif
      endif 

      location=0                ! default ( P is outside S)
      icross=0                  ! initialize counter

      if ((plat.eq.xlat_c).and.(plon.eq.xlon_c)) then
         location=1
         return
      endif

      
      call TrnsfmLon (xlat_c,xlon_c,plat,plon,tlonP)

      do i=1,nv_c              ! start of loop over sides of S 

         vAlat=vlat_c(i)
         vAlon=vlon_c(i)
         tlonA=tlonv(i)

         if (i.lt.nv_c) then
            vBlat=vlat_c(i+1)
            vBlon=vlon_c(i+1)
            tlonB=tlonv(i+1)
         else
            vBlat=vlat_c(1)
            vBlon=vlon_c(1)
            tlonB=tlonv(1)
         endif
         
         istrike=0
         
         if (tlonP.eq.tlonA) then
            istrike=1
         else
            call EastOrWest(tlonA,tlonB,ibrngAB)
            call EastOrWest(tlonA,tlonP,ibrngAP)
            call EastOrWest(tlonP,tlonB,ibrngPB)
            

            if((ibrngAP.eq.ibrngAB).and.(ibrngPB.eq.ibrngAB)) istrike=1
         endif

         
         if (istrike.eq.1) then

            if ((plat.eq.vAlat).and.(plon.eq.vAlon)) then
               location=2       ! P lies on a vertex of S
               return
            endif
            call TrnsfmLon(vAlat,vAlon,xlat_c,xlon_c,tlon_X)
            call TrnsfmLon(vAlat,vAlon,vBlat,vBlon,tlon_B)
            call TrnsfmLon(vAlat,vAlon,plat,plon,tlon_P)
            
            if (tlon_P.eq.tlon_B) then
               location=2       ! P lies on side of S
               return 
            else
               call EastOrWest(tlon_B,tlon_X,ibrng_BX)
               call EastOrWest(tlon_B,tlon_P,ibrng_BP)
               if(ibrng_BX.eq.(-ibrng_BP)) icross=icross+1
            endif
            
         endif
      enddo                     ! end of loop over the sides of S


c     if the arc XP crosses the boundary S an even number of times then P
c     is in S

      if (mod(icross,2).eq.0) location=1

      return

      end


c     ****************************************************************
      
      subroutine TrnsfmLon(plat,plon,qlat,qlon,tranlon)

c     ****************************************************************
c     * This subroutine is required by subroutines DefSPolyBndry &   *
c     * LctPtRelBndry. It finds the 'longitude' of point Q in a      *
c     * geographic coordinate system for which point P acts as a     *
c     * 'north pole'. SENT: plat,plon,qlat,qlon, in degrees.         *
c     * RETURNED: tranlon, in degrees.                               *
c     ****************************************************************

      implicit none

      real pi,dtr,plat,plon,qlat,qlon,tranlon,t,b
      parameter (pi=3.141592654,dtr=pi/180.0)
 
      if (plat.eq.90.) then
         tranlon=qlon
      else
         t=sin((qlon-plon)*dtr)*cos(qlat*dtr)
         b=sin(dtr*qlat)*cos(plat*dtr)-cos(qlat*dtr)*sin(plat*dtr)
     >    *cos((qlon-plon)*dtr)
         tranlon=atan2(t,b)/dtr
      endif

      return
      end

c     ****************************************************************

      subroutine EastOrWest(clon,dlon,ibrng)

c     ****************************************************************
c     * This subroutine is required by subroutine LctPtRelBndry.     *
c     * This routine determines if in travelling the shortest path   *
c     * from point C (at longitude clon) to point D (at longitude    *
c     * dlon) one is heading east, west or neither.                  *
c     * SENT: clon,dlon; in degrees. RETURNED: ibrng                 *
c     * (1=east,-1=west, 0=neither).                                 *
c     ****************************************************************

      implicit none
      real clon,dlon,del
      integer ibrng
      del=dlon-clon
      if (del.gt.180.) del=del-360.
      if (del.lt.-180.) del=del+360.
      if ((del.gt.0.0).and.(del.ne.180.)) then
         ibrng=-1               ! (D is west of C)
      elseif ((del.lt.0.0).and.(del.ne.-180.)) then
         ibrng=+1               ! (D is east of C)
      else
         ibrng=0                ! (D north or south of C)
      endif
      return
      end
