!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
!
! !PROGRAM: zeit_ci.f
!
! !DESCRIPTION:  A check in program for timing processes in shell scripts
!
! !USAGE:        This is to be used as zeit_ci.x [-r fname] [-v] [name]
!                where:  -r    Time register file name; default is .zeit
!                        -v    Verbose mode: prints current wallclock time
!                        name  Name of program to be timed
            

      implicit none
      integer i,j,n
      character*8 dd
      character*10 tt
      character*5 zz 
      integer vv(8)
      character*200 name,argin(3),user_def_reg
      integer*4 iargc
      logical present,regon
      regon=.false.
c Test for the existence of .zeit.  If it's not around then create
c a new file later.  Otherwise open the old one now.
      
c     Process argument list
      n=iargc()
      select case(n)
      case(0)             !Print usage list
         call zeit_usage('zeit_ci')
         stop
      case(1)                  
         call getarg(1,argin(1))                 
         if(argin(1).eq.'-v')then
            name='main'
            inquire(file='.zeit',exist=present)
            if(present)then
               open(11,file='.zeit',status='old',position='append')
            else
               open(11,file='.zeit',status='unknown')
            endif
            call date_and_time(dd,tt,zz,vv)
            write(11,100)name,dd(1:4),dd(5:6),dd(7:8),tt(1:2),
     $           tt(3:4),tt(5:6),tt(8:10)
 100        format('zeit_ci ',a50,1x,'Date ',a4,a2,
     $           a2,' Time ',a2,1x,a2,1x,a2,1x,a3)
            call date_and_time(dd,tt,zz,vv)
            write(6,101)name,dd(1:4),dd(5:6),dd(7:8),tt(1:2),
     $           tt(3:4),tt(5:6),tt(8:10)
         else
            name=argin(1)
            inquire(file='.zeit',exist=present)
            if(present)then
               open(11,file='.zeit',status='old',position='append')
            else
               open(11,file='.zeit',status='unknown')
            endif
            call date_and_time(dd,tt,zz,vv)
            write(11,100)name,dd(1:4),dd(5:6),dd(7:8),tt(1:2),
     $           tt(3:4),tt(5:6),tt(8:10)
         endif
      case(2)                     
         call getarg(1,argin(1))
         call getarg(2,argin(2))      
         if(argin(1).eq.'-v')then
            name=argin(2)
            inquire(file='.zeit',exist=present)
            if(present)then
               open(11,file='.zeit',status='old',position='append')
            else
               open(11,file='.zeit',status='unknown')
            endif
            call date_and_time(dd,tt,zz,vv)
            write(11,100)name,dd(1:4),dd(5:6),dd(7:8),tt(1:2),
     $           tt(3:4),tt(5:6),tt(8:10)
            write(6,101)name,dd(1:4),dd(5:6),dd(7:8),tt(1:2),
     $           tt(3:4),tt(5:6),tt(8:10)
 101        format('zeit_ci ',a16,1x,'Date ',a4,'/',a2,'/',
     $           a2,' Time ',a2,
     $           ':',a2,':',a2,'.',a3)
         elseif(argin(1).eq.'-r')then
            user_def_reg=argin(2)
            name='main'
            inquire(file=user_def_reg,exist=present)
            if(present)then
               open(11,file=user_def_reg,status='old',
     $              position='append')
            else
               open(11,file=user_def_reg,status='unknown')
            endif
            call date_and_time(dd,tt,zz,vv)
            write(11,100)name,dd(1:4),dd(5:6),dd(7:8),tt(1:2),
     $           tt(3:4),tt(5:6),tt(8:10)
         else
            call zeit_usage('zeit_ci')
            stop
         endif
      case(3)
         call getarg(1,argin(1))
         call getarg(2,user_def_reg)
         call getarg(3,argin(3))
         if(argin(1).eq.'-r')then
            if(argin(3).eq.'-v')then
               name='main'
               inquire(file=user_def_reg,exist=present)
               if(present)then
                  open(11,file=user_def_reg,status='old',
     $                 position='append')
               else
                  open(11,file=user_def_reg,status='unknown')
               endif
               call date_and_time(dd,tt,zz,vv)
               write(11,100)name,dd(1:4),dd(5:6),dd(7:8),tt(1:2),
     $              tt(3:4),tt(5:6),tt(8:10)
               write(6,101)name,dd(1:4),dd(5:6),dd(7:8),tt(1:2),
     $              tt(3:4),tt(5:6),tt(8:10)
            else
               name=argin(3)
               inquire(file=user_def_reg,exist=present)
               if(present)then
                  open(11,file=user_def_reg,status='old',
     $                 position='append')
               else
                  open(11,file=user_def_reg,status='unknown')
               endif
               call date_and_time(dd,tt,zz,vv)
               write(11,100)name,dd(1:4),dd(5:6),dd(7:8),tt(1:2),
     $              tt(3:4),tt(5:6),tt(8:10)
            endif
         else
            call zeit_usage('zeit_ci')
            stop
         endif
      case(4)
         call getarg(1,argin(1))
         call getarg(2,argin(2))
         call getarg(3,argin(3))
         call getarg(4,name)
         if(argin(1).eq.'-r'.and.argin(3).eq.'-v')then
            user_def_reg=argin(2)
            call date_and_time(dd,tt,zz,vv)
            write(6,101)name,dd(1:4),dd(5:6),dd(7:8),tt(1:2),
     $           tt(3:4),tt(5:6),tt(8:10)
            inquire(file=user_def_reg,exist=present)
            if(present)then
               open(11,file=user_def_reg,status='old',
     $              position='append')
            else
               open(11,file=user_def_reg,status='unknown')
            endif
            call date_and_time(dd,tt,zz,vv)
            write(11,100)name,dd(1:4),dd(5:6),dd(7:8),tt(1:2),
     $           tt(3:4),tt(5:6),tt(8:10)
         else
            call zeit_usage('zeit_ci')
            stop
         endif
      case default
         call zeit_usage('zeit_ci')
         stop
      end select
      stop
      end

