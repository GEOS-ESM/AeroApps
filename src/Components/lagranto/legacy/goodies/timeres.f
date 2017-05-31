      PROGRAM timeresolution
      
c     ***********************************************************************
c     * Change time resolution of of a trajectory file                      *
c     * Michael Sprenger / Winter 2010                                      *
c     ***********************************************************************

      implicit none
      
c     ----------------------------------------------------------------------
c     Declaration of variables
c     ----------------------------------------------------------------------

c     Parameters
      character*80                           inpfile     
      character*80                           outfile     
      integer                                ntra,otim,ncol
      real                                   timeres
      character*80                           unit
      character*80                           mode
      real                                   shift

c     Trajectories
      character*80                           vars(100)   
      integer                                refdate(6)  
      integer                                ntim       
      real,allocatable, dimension (:,:,:) :: trainp        
      real,allocatable, dimension (:,:,:) :: traout        
      real,allocatable, dimension (:)     :: timold,timnew
      real,allocatable, dimension (:)     :: fldold,fldnew

c     Numerical constants
      real                                   eps
      parameter                              (eps=0.001)

c     Auxiliary variables
      integer                                inpmode
      integer                                outmode
      integer                                stat
      integer                                fid
      integer                                i,j,k
      real                                   hhmm,tfrac
      real                                   range
      integer                                date0(5),date1(5)

c     ----------------------------------------------------------------------
c     Parameter handling
c     ----------------------------------------------------------------------

c     Read parameters
      open(10,file='timeres.param')
       read(10,*) inpfile
       read(10,*) outfile
       read(10,*) ntra,otim,ncol
       read(10,*) timeres
       read(10,*) unit
       read(10,*) mode
       read(10,*) shift
      close(10)

c     Change unit to hours
      if ( unit.eq.'min') then
         timeres = 1./60. * timeres
         unit    = 'h'
      endif

c     Determine the formats
      call mode_tra(inpmode,inpfile)
      if (inpmode.eq.-1) inpmode=1
      call mode_tra(outmode,outfile)
      if (outmode.eq.-1) outmode=1

c     ----------------------------------------------------------------------
c     Read input trajectory and allocate memory
c     ----------------------------------------------------------------------

c     Allocate memory for input trajectories
      allocate(trainp(ntra,otim,ncol),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array trainp   ***' 
      allocate(timold(otim),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array timold   ***' 
      allocate(fldold(otim),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array fldold   ***' 

c     Read inpufile
      call ropen_tra(fid,inpfile,ntra,otim,ncol,refdate,vars,inpmode)
      call read_tra (fid,trainp,ntra,otim,ncol,inpmode)
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

c     Convert all times from hhmm to fractional time
      do i=1,ntra
         do j=1,otim
            hhmm = trainp(i,j,1)
            call hhmm2frac(hhmm,tfrac)
            trainp(i,j,1) = tfrac
         enddo
      enddo

c     Get the time range in hours
      range = ( trainp(1,otim,1) - trainp(1,1,1) ) 

c     If timeres=0, keep the original resolution
      if ( abs(timeres).lt.eps ) then
         timeres = trainp(1,2,1) - trainp(1,1,1)
         print*,'Keeping time resolution',timeres
      endif

c     Determine the new number of times
      ntim = nint( abs( range ) / timeres ) + 1 

c     Check that the time range and new time resolution are consistent
      if ( abs( real(ntim-1) * timeres - abs(range) ).gt.eps ) then
         print*,' ERROR: time range and resolution are not compatible'
         print*,'   range              = ',range
         print*,'   (ntim-1) * timeres = ',real(ntim-1) * timeres
         stop
      endif

c     Allocate memory for output trajectories
      allocate(traout(ntra,ntim,ncol),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array trainp   ***' 
      allocate(timnew(ntim),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array timnew   ***' 
      allocate(fldnew(ntim),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array fldnew   ***' 

c     ----------------------------------------------------------------------
c     Change time resolution
c     ----------------------------------------------------------------------

c     Define the old and new times
      do i=1,otim
         timold(i) = trainp(1,i,1)
      enddo
      do i=1,ntim
         timnew(i) = timold(1) + real(i-1) * timeres
      enddo

c     Change time resolution
      do i=1,ntra
      do k=2,ncol

c        Copy old field
         do j=1,otim
            fldold(j) = trainp(i,j,k)
         enddo
            
c        Exception: Handle date line problem for longitude
         if ( k.eq.2 ) then 
            do j=2,otim
               if ( (fldold(j-1)-fldold(j)).gt.180. ) then
                  fldold(j) = fldold(j) + 360.
               else if ( (fldold(j-1)-fldold(j)).lt.-180. ) then
                  fldold(j) = fldold(j) - 360.
               endif
            enddo
         endif
         
c        Cubic spline fitting
         if ( mode.eq.'cubic' ) then
            call cubicfit (timold,fldold,otim,timnew,fldnew,ntim)
         else if (mode.eq.'linear' ) then
            call linearfit(timold,fldold,otim,timnew,fldnew,ntim)
         endif

c        Exception: Reverse date line handling for longitude
         if ( k.eq.2 ) then
            do j=1,ntim
               if ( fldnew(j).gt.180. ) then
                  fldnew(j) = fldnew(j) -360.
               else if ( fldnew(j).lt.-180. ) then
                  fldnew(j) = fldnew(j) +360.
               endif
            enddo
         endif

c        Save the new field in the output trajectory
         do j=1,ntim
            traout(i,j,1) = timnew(j)
            traout(i,j,k) = fldnew(j)
         enddo

      enddo
      enddo

c     ----------------------------------------------------------------------
c     Shift time axis - change reference date
c     ----------------------------------------------------------------------

c     Shift the reference date
      date0(1) = refdate(1)
      date0(2) = refdate(2)
      date0(3) = refdate(3)
      date0(4) = refdate(4)
      date0(5) = refdate(5)
      call newdate (date0,shift,date1)
      refdate(1) = date1(1)
      refdate(2) = date1(2)
      refdate(3) = date1(3)
      refdate(4) = date1(4)
      refdate(5) = date1(5)

c     Shift the times
      do i=1,ntra
      do j=1,ntim
        traout(i,j,1) = traout(i,j,1) - shift
      enddo
      enddo

c     ----------------------------------------------------------------------
c     Write output trajectory
c     ----------------------------------------------------------------------

c     Convert all times from fractional to hhmm time
      do i=1,ntra
         do j=1,ntim
            tfrac = traout(i,j,1)
            call frac2hhmm(tfrac,hhmm)
            traout(i,j,1) = hhmm
         enddo
      enddo
    
c     Write output file
      call wopen_tra(fid,outfile,ntra,ntim,ncol,refdate,vars,outmode)
      call write_tra(fid,traout,ntra,ntim,ncol,outmode)
      call close_tra(fid,outmode)
      
      end


c     ********************************************************************
c     * REPARAMETERIZATION SUBROUTINES                                   *
c     ********************************************************************

c     -------------------------------------------------------------
c     Interpolation of the trajectory with linear interpolation
c     -------------------------------------------------------------

      SUBROUTINE linearfit (time,lon,n,sptime,splon,spn)

c     Given the curve <time,lon> with <n> data points, fit a
c     linear fit to this curve. The new curve is returned in 
c     <sptime,splon,spn> with <spn> data points. The parameter
c     <spn> specifies on entry the number of interpolated points
c     along the curve.
      
      implicit none

c     Declaration of subroutine parameters
      integer n
      real    time(n),lon(n)
      integer spn
      real    sptime(spn),splon(spn)

c     Auxiliary variables
      real    dt
      real    s
      integer i,j,iold
      real    order

c     Determine whether the input array is ascending or descending
      if (time(1).gt.time(n)) then
         order=-1.
      else
         order= 1.
      endif

c     Bring the time array into ascending order
      do i=1,n
         time(i)=order*time(i)
      enddo

c     Prepare the linear interpolation: define the new times
      dt=(time(n)-time(1))/real(spn-1)
      do i=1,spn
         sptime(i)=time(1)+real(i-1)*dt
      enddo
      
c     Do the interpolation
      iold = 1
      do i=1,spn

c        Decide which interval of the old time series must be taken
         do j=iold,n-1
            if ( ( sptime(i).ge.time(j  ) ).and.
     >           ( sptime(i).lt.time(j+1) ) ) 
     >      then
               iold = j
               exit
            endif
         enddo
         
c        Do the linear interpolation
         splon(i) = lon(iold) + ( lon(iold+1) - lon(iold) ) * 
     >       ( sptime(i) - time(iold) ) / ( time(iold+1) - time(iold) ) 

      enddo

c     Change the time arrays back: original order
      do i=1,spn
         sptime(i)=order*sptime(i)
      enddo
      do i=1,n
         time(i)=order*time(i)
      enddo

      return
      end


c     -------------------------------------------------------------
c     Interpolation of the trajectory with a natural cubic spline
c     -------------------------------------------------------------

      SUBROUTINE cubicfit (time,lon,n,sptime,splon,spn)

c     Given the curve <time,lon> with <n> data points, fit a
c     cubic spline to this curve. The new curve is returned in 
c     <sptime,splon,spn> with <spn> data points. The parameter
c     <spn> specifies on entry the number of spline interpolated points
c     along the curve.
      
      implicit none

c     Declaration of subroutine parameters
      integer n
      real    time(n),lon(n)
      integer spn
      real    sptime(spn),splon(spn)

c     Auxiliary variables
      real    y2ax(n)
      real    dt
      real    s
      integer i
      real    order

c     Determine whether the input array is ascending or descending
      if (time(1).gt.time(n)) then
         order=-1.
      else
         order= 1.
      endif

c     Bring the time array into ascending order
      do i=1,n
         time(i)=order*time(i)
      enddo

c     Prepare the (natural) cubic spline interpolation
      call spline (time,lon,n,1.e30,1.e30,y2ax)
      dt=(time(n)-time(1))/real(spn-1)
      do i=1,spn
         sptime(i)=time(1)+real(i-1)*dt
      enddo
      
c     Do the spline interpolation
      do i=1,spn
         call splint(time,lon,y2ax,n,sptime(i),s)
         splon(i)=s
      enddo

c     Change the time arrays back
      do i=1,spn
         sptime(i)=order*sptime(i)
      enddo
      do i=1,n
         time(i)=order*time(i)
      enddo

      return
      end

c     -------------------------------------------------------------
c     Basic routines for spline interpolation (Numerical Recipes)
c     -------------------------------------------------------------

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
         print*, 'bad xa input in splint'
         stop
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      END
    

      
