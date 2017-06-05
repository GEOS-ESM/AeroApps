
      SUBROUTINE special (flag,cmd,tra,ntim,ncol,
     >                    vars,times,param,nparam)

c     ***************************************************************************
c     *                                                                         *
c     * OUTPUT:  flag           -> 1 if trajectory is selected, 0 if not        *
c     *                                                                         *
c     * INPUT:   cmd            <- command string                               *
c     *          tra(ntim,ncol) <- single trajectory: indices time,column       *
c     *          ntim           <- number of times                              *
c     *          ncol           <- number of columns (including time,lon,lat,p) *
c     *          vars(ncol)     <- names of columns                             *
c     *          times(ntim)    <- List of times
c     *          param(nparam)  <- parameter values                             *
c     *          nparam         <- number of parameters                         *
c     *                                                                         *
c     ***************************************************************************

      implicit none
      
c     ---------------------------------------------------------------------------
c     Declaration of subroutine parameters
c     ---------------------------------------------------------------------------

      integer       flag           ! Boolean flag whether trajectory is selected
      character*80  cmd            ! Command string
      integer       ntim,ncol      ! Dimension of single trajectory
      real          tra(ntim,ncol) ! Single trajectory
      character*80  vars(ncol)     ! Name of columns
      real          times(ntim)    ! List of times
      integer       nparam         ! # parameters
      real          param(nparam)  ! List of parameters

c     ---------------------------------------------------------------------------
c     Declaration of local variables
c     ---------------------------------------------------------------------------

      integer       i
      integer       ip,i0,i1

c     --------------------------------------------------------------------------  %)
c     SPECIAL:WCB:ascent,first,last                                               %)
c         : Detect Warm Conveyor Belts (WCB); the air stream must ascend at least %)
c         : <ascent=param(1)> hPa between the two times <first=param(2)> and      %)
c         : <last=param(3)>. Note, the lowest pressure is allowed to occur at any %)
c         : time between <first> and <last>.                                      %)
c     --------------------------------------------------------------------------- %)

      if ( cmd.eq.'WCB' ) then

c        Reset the flag for selection
         flag = 0

c        Pressure is in the 4th column
         ip = 4

c        Get times
         i0 = 0
         i1 = 0
         do i=1,ntim
            if ( param(2).eq.times(i) ) i0 = i
            if ( param(3).eq.times(i) ) i1 = i
         enddo
         if ( (i0.eq.0).or.(i1.eq.0) ) then
            print*,' ERROR: invalid times in SPECIAL:WCB... Stop'
            stop
         endif

c        Check for ascent 
         do i=i0+1,i1
            if ( ( tra(i0,ip)-tra(i,ip) ) .gt. param(1) ) flag = 1
         enddo

      endif

c     ---------------------------------------------------------------------------


      end
