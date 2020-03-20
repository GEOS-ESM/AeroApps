      program gettidiff
C     =================

      implicit none
 
      integer   idate(5),irefdat(5)
      real      ihdiff

      integer   iargc
      character*(80) arg
      integer	nc1,nc2,flag1,flag2
 
c     check for sufficient requested arguments
      if (iargc().ne.2) then
         print*,
     >  'USAGE: gettidiff date1 date2 (format (YY)YYMMDD_HH(MM))'
         call exit(1)
      endif
 
c     read and transform input
      call getarg(1,arg)
      call lenchar(arg,nc1)
      call checkchar(arg,'_',flag1)

      idate(5)   = 0
      irefdat(5) = 0

      if (flag1.eq.7) then
        read(arg(1:2),'(i2)',err=120)idate(1)
        read(arg(3:4),'(i2)',err=120)idate(2)
        read(arg(5:6),'(i2)',err=120)idate(3)
        read(arg(8:9),'(i2)',err=120)idate(4)
        if (nc1.eq.11) then
          read(arg(10:11),'(i2)',err=120)idate(5)
        else if (nc1.ne.9) then
          print*,
     >   'USAGE: gettidiff date1 date2 (format (YY)YYMMDD_HH(MM))'
          call exit(1)
        endif
      else if (flag1.eq.9) then
        read(arg(1:4),'(i4)',err=120)idate(1)
        read(arg(5:6),'(i2)',err=120)idate(2)
        read(arg(7:8),'(i2)',err=120)idate(3)
        read(arg(10:11),'(i2)',err=120)idate(4)
        if (nc1.eq.13) then
          read(arg(12:13),'(i2)',err=120)idate(5)
        else if (nc1.ne.11) then
          print*,
     >   'USAGE: gettidiff date1 date2 (format (YY)YYMMDD_HH(MM))'
          call exit(1)
        endif
      else
        print*,
     > 'USAGE: gettidiff date1 date2 (format (YY)YYMMDD_HH(MM))'
        call exit(1)
      endif
 
      call getarg(2,arg)
      call lenchar(arg,nc2)
      call checkchar(arg,'_',flag2)
      if (flag1.ne.flag2) then
        print*,
     > 'error: both dates must be in same format (YY)YYMMDD_HH(MM)'
        call exit(1)
      endif

      if (flag2.eq.7) then
        read(arg(1:2),'(i2)',err=120)irefdat(1)
        read(arg(3:4),'(i2)',err=120)irefdat(2)
        read(arg(5:6),'(i2)',err=120)irefdat(3)
        read(arg(8:9),'(i2)',err=120)irefdat(4)
        if (nc2.eq.11) then
          read(arg(10:11),'(i2)',err=120)irefdat(5)
        else if (nc2.ne.9) then
          print*,
     >   'USAGE: gettidiff date1 date2 (format (YY)YYMMDD_HH(MM))'
          call exit(1)
        endif
      else if (flag2.eq.9) then
        read(arg(1:4),'(i4)',err=120)irefdat(1)
        read(arg(5:6),'(i2)',err=120)irefdat(2)
        read(arg(7:8),'(i2)',err=120)irefdat(3)
        read(arg(10:11),'(i2)',err=120)irefdat(4)
        if (nc2.eq.13) then
          read(arg(12:13),'(i2)',err=120)irefdat(5)
        else if (nc2.ne.11) then
          print*,
     >   'USAGE: gettidiff date1 date2 (format (YY)YYMMDD_HH(MM))'
          call exit(1)
        endif
      else
        print*,
     >   'USAGE: gettidiff date1 date2 (format (YY)YYMMDD_HH(MM))'
        call exit(1)
      endif

      call timediff(idate,irefdat,ihdiff)

      if (int(100.*ihdiff).eq.100*int(ihdiff)) then
        write(*,*)int(ihdiff)
      else
        write(*,'(f7.2)')ihdiff
      endif

      goto 200
 
 120  write(*,*)
     >"*** error: date must be in format (YY)YYMMDD_HH(MM) ***"
 
 200  continue
      end

      subroutine lenchar(string,lstr)
C     ===============================
 
      character*(*)     string
      integer   n,lstr
 
      do n=1,len(string)
        if (string(n:n).eq."") then
          lstr=n-1
          goto 100
        endif
      enddo
 100  continue
      end

      subroutine checkchar(string,char,flag)
C     ======================================

      character*(*)     string
      character*(1)     char
      integer   n,flag

      flag=0
      do n=1,len(string)
        if (string(n:n).eq.char) then
          flag=n
          return
        endif
      enddo
      end
