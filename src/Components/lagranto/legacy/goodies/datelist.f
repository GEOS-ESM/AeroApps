      PROGRAM datelist

c     **********************************************************************
c     * Handling of date lists                                             *
c     * Michael Sprenger                                                   *
c     **********************************************************************

      implicit none

c     ---------------------------------------------------------------------
c     Declaration of variables
c     ---------------------------------------------------------------------

c     Parameters
      character*80              datefile
      character*80              mode 
      character*80              startdate
      character*80              finaldate
      character*80              refdate
      integer                   interval
      real                      randpercent

c     Date list
      integer                   ndates
      integer,allocatable,dimension(:)    :: year,month,day,hour

c     Auxiliary variables
      integer                     i
      integer                     year1,month1,day1,hour1,min1
      integer                     year2,month2,day2,hour2,min2
      character                   direction
      integer                     date1(5),date2(5)
      character*11                datestr,datestr1,datestr2
      character                   ch
      real                        diff
      real                        time
      character*80                datefile1,datefile2
      integer                     timestamp1,timestamp2
      integer                     oldstamp1 ,oldstamp2
      integer                     state

c     ---------------------------------------------------------------------
c     Preparations
c     ---------------------------------------------------------------------

c     Read parameter file
      open(10,file='datelist.param')

       read(10,*) datefile                           ! General parameters
       read(10,*) mode

       if ( mode.eq.'-create' ) then                 ! Create date list
          read(10,*) startdate
          read(10,*) finaldate
          read(10,*) interval

       elseif ( mode.eq.'-totime' ) then             ! Convert to time list
          read(10,*) refdate

       elseif ( mode.eq.'-todate' ) then             ! Convert to date list
          read(10,*) refdate
       
       elseif ( mode.eq.'-onlyin1' ) then            ! Dates only in file 1
          read(10,*) datefile1
          read(10,*) datefile2

       else                                          ! Invalid mode
          print*,' ERROR: invalid mode for datelist'
          stop
       endif
       
      close(10)

c     ---------------------------------------------------------------------
c     Create a date list (-create)
c     ---------------------------------------------------------------------

      if ( mode.ne.'-create' ) goto 100

c     Check whether interval is ok
      if ( ( interval.le.0 ).or.(interval.gt.24) ) then
         print*,'Interval must be between 1 h and 24 h... Stop'
         stop
      endif

c     Extract dates and times
      read(startdate( 1: 4),*) year1
      read(startdate( 5: 6),*) month1
      read(startdate( 7: 8),*) day1
      read(startdate(10:11),*) hour1
      read(startdate(12:13),*) min1
      
      read(finaldate( 1: 4),*) year2
      read(finaldate( 5: 6),*) month2
      read(finaldate( 7: 8),*) day2
      read(finaldate(10:11),*) hour2
      read(finaldate(12:13),*) min2

c     Get direction of the date file
      if (year2.gt.year1) then
         direction = 'f'
         goto 101
      elseif (year2.lt.year1) then
         direction = 'b'
         goto 101
      endif

      if (month2.gt.month1) then
         direction = 'f'
         goto 101
      elseif (month2.lt.month1) then
         direction = 'b'
         goto 101
      endif

      if (day2.gt.day1) then
         direction = 'f'
         goto 101
      elseif (day2.lt.day1) then
         direction = 'b'
         goto 101
      endif
      
      if (hour2.gt.hour1) then
         direction = 'f'
         goto 101
      elseif (hour2.lt.hour1) then
         direction = 'b'
         goto 101
      endif

      if (min2.gt.min1) then
         direction = 'f'
         goto 101
      elseif (min2.lt.min1) then
         direction = 'b'
         goto 101
      endif

      direction = 'f'

 101  continue

c     Set the interval step depending on the direction
      if ( direction.eq.'b' ) then
         interval = -interval
      endif

c     Save the dates in arrays
      date1(1) = year1
      date1(2) = month1
      date1(3) = day1
      date1(4) = hour1 
      date1(5) = 0

      date2(1) = year2
      date2(2) = month2
      date2(3) = day2
      date2(4) = hour2 
      date2(5) = 0


c     Get starting and ending date for the date list
      if ( direction.eq.'f' ) then

         do while ( mod(date1(4),interval) .ne. 0 )
            date1(4) = date1(4) - 1
         enddo

         if (min2.ne.0) call newdate(date2,1.,date2)

         do while ( mod(date2(4),interval) .ne. 0 )
            date2(4) = date2(4) + 1
         enddo

      else

         if (min1.ne.0) call newdate(date1,1.,date1)

         do while ( mod(date1(4),interval) .ne. 0 )
            date1(4) = date1(4) + 1
         enddo
       
         do while ( mod(date2(4),interval) .ne. 0 )
            date2(4) = date2(4) - 1
         enddo
         
      endif

c     Create and write the datefile
      if ( datefile.ne.'/dev/stdout') then
         open(10,file=datefile)
      endif

 102  continue

       call datestring(datestr,date1(1),date1(2),date1(3),date1(4) )

       if ( datefile.ne.'/dev/stdout') then
          write(10,*) datestr
       else
          write(*,*) datestr
       endif

       if ( ( date1(1).ne.date2(1) ).or.
     >      ( date1(2).ne.date2(2) ).or.
     >      ( date1(3).ne.date2(3) ).or.
     >      ( date1(4).ne.date2(4) ) ) 
     > then
          diff = real(interval)
          call newdate(date1,diff,date1)
          goto 102
       endif

      if ( datefile.ne.'/dev/stdout') then
         close(10)
      endif

 100  continue

c     ---------------------------------------------------------------------
c     Convert dates to a list of times
c     ---------------------------------------------------------------------

      if ( mode.ne.'-totime' ) goto 112

c     Extract reference date
      read(refdate( 1: 4),*) year1
      read(refdate( 5: 6),*) month1
      read(refdate( 7: 8),*) day1
      read(refdate(10:11),*) hour1
      read(refdate(12:13),*) min1
      
c     Loop through the date file
      open(10,file=datefile)

 111  read(10,*,end=110) datestr
      
c     Extract date
      read(datestr( 1: 4),*) year2
      read(datestr( 5: 6),*) month2
      read(datestr( 7: 8),*) day2
      read(datestr(10:11),*) hour2
      min1 = 0

c     Get the time difference
      date1(1) = year1
      date1(2) = month1
      date1(3) = day1
      date1(4) = hour1 
      date1(5) = 0

      date2(1) = year2
      date2(2) = month2
      date2(3) = day2
      date2(4) = hour2 
      date2(5) = 0

      call timediff(date2,date1,diff)
      

c     Write it to screen
      write(*,'(i6)') nint(diff)

      goto 111

 110  continue
      
c     Close datefile
      close(10)
      
c     Exit point      
  112 continue

c     ---------------------------------------------------------------------
c     Convert times to a list of dates 
c     ---------------------------------------------------------------------

      if ( mode.ne.'-todate' ) goto 122

c     Extract reference date
      read(refdate( 1: 4),*) year1
      read(refdate( 5: 6),*) month1
      read(refdate( 7: 8),*) day1
      read(refdate(10:11),*) hour1
      read(refdate(12:13),*) min1
      
c     Loop through the date file
      open(10,file=datefile)

 121  read(10,*,end=120) time

c     Calculate the new date
      date1(1) = year1
      date1(2) = month1
      date1(3) = day1
      date1(4) = hour1 
      date1(5) = 0
      call newdate(date1,time,date2)
      call datestring(datestr,date2(1),date2(2),date2(3),date2(4) )

c     Write it to screen
      write(*,'(a11)') trim(datestr)

      goto 121

 120  continue
      
c     Close datefile
      close(10)
      
c     Exit point      
 122  continue
   
c     ---------------------------------------------------------------------
c     Extract all dates which are only in one datefile 
c     ---------------------------------------------------------------------

      if ( mode.ne.'-onlyin1' ) goto 134
      
c     Set reference date   
      date1(1) = 1979
      date1(2) = 1
      date1(3) = 1
      date1(4) = 0
      date1(5) = 0
      
c     Open the output file
      if ( datefile.ne.'/dev/stdout') then
         open(30,file=datefile)
      endif

c     Loop through the input date files
      open(10,file=datefile1)
      open(20,file=datefile2)
    
c     Loop through both date files
      state     = 0
      oldstamp1 = 0
      oldstamp2 = 0
      
 131  if ( (state.eq.1).or.(state.eq.0) ) then 
         read(10,*,end=133) datestr1
      endif
      if ( (state.eq.2).or.(state.eq.0) ) then 
         read(20,*,end=130) datestr2
      endif

c     Get time stamp for both date strings
      if ( (state.eq.1).or.(state.eq.0) ) then 
         read(datestr1( 1: 4),*) year2
         read(datestr1( 5: 6),*) month2
         read(datestr1( 7: 8),*) day2
         read(datestr1(10:11),*) hour2
         date2(1) = year2
         date2(2) = month2
         date2(3) = day2
         date2(4) = hour2 
         date2(5) = 0
         call timediff(date2,date1,diff)
         timestamp1 =  nint( diff )
         if ( timestamp1.lt.oldstamp1 ) then
            print*,' ERROR: datelist must be ordered ',trim(datefile1)
            stop
         else
            oldstamp1 = timestamp1
         endif
      endif
      
      if ( (state.eq.1).or.(state.eq.0) ) then 
         read(datestr2( 1: 4),*) year2
         read(datestr2( 5: 6),*) month2
         read(datestr2( 7: 8),*) day2
         read(datestr2(10:11),*) hour2
         date2(1) = year2
         date2(2) = month2
         date2(3) = day2
         date2(4) = hour2 
         date2(5) = 0
         call timediff(date2,date1,diff)
         timestamp2 =  nint( diff )
         if ( timestamp2.lt.oldstamp2 ) then
            print*,' ERROR: datelist must be ordered ',trim(datefile2)
            stop
         else
            oldstamp2 = timestamp2
         endif
      endif
      
c     Write output and set new state
      if ( timestamp1.gt.timestamp2 ) then
        state = 2
c        print*,trim(datestr1)//'.'//trim(datestr2)//' >'

      else if ( timestamp1.lt.timestamp2 ) then
        state = 1

        if ( datefile.ne.'/dev/stdout') then
           write(30,*) datestr1
        else
           write(*,*) datestr1
        endif

c        print*,trim(datestr1)//'.'//trim(datestr2)//' -> out'

      else if (timestamp1.eq.timestamp2 ) then
        state = 0
c        print*,trim(datestr1)//'.'//trim(datestr2)//' ='

      endif

      goto 131

c     Exit point for parallel reading through files
 130  continue
 
c     Write remaining part of datefile 1
 132  continue
        read(10,*,end=133) datestr1
        if ( datefile.ne.'/dev/stdout') then
           write(30,*) datestr1
        else
           write(*,*) datestr1
        endif
      goto 132
 133  continue
      
c     Close datefile
      close(10)
      close(20)
      close(30)
      
c     Exit point      
 134  continue
 
c     ---------------------------------------------------------------------
c     End
c     ---------------------------------------------------------------------

      end
  
