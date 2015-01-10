c                       code modifications based on:  04/02/93 - analy.f
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sort(expid,verbose,luverb,nnobs,
     1  rlats,rlons,rlevs,kx,kt,del,sigU,sigO,sigF,tstamp,
     1  maxreg,ktmax,iregbeg,ireglen,ityplen)

c.... Sort data by regions and set region vectors iregbeg and ireglen;
c.... Sort data in each region by kt-kx-rlats-rlons-rlevs and set the
c.... array ityplen; print (optional) summary of data in all regions,
c.... print data in selected regions.

c.......................................................................
	use m_die,only : die
        use m_mpout,only: mpout
	use m_stdio,only : stderr

      character*4 myname
      parameter(myname='sort')

      character*8   expid

c.......................................................................
c.... Control parameter for output.

      logical      verbose
      integer      luverb

c.......................................................................
c.... Srorage for data items.

      integer      nnobs
      real         rlats(nnobs)
      real         rlons(nnobs)
      real         rlevs(nnobs)
      integer      kx(nnobs)
      integer      kt(nnobs)
      real         del(nnobs)
      real         sigU(nnobs)
      real         sigO(nnobs)
      real         sigF(nnobs)
      real         tstamp(nnobs)
c-    real         qcflag(nnobs)		! no longer used

      integer      maxreg
      integer      ktmax			! must == ktmax.h
      integer      iregbeg(maxreg)
      integer      ireglen(maxreg)
      integer      ityplen(ktmax,maxreg)

c	..Locals
	integer kr
	integer i
	real    rlatmin,rlatmax
	real	rlonmin,rlonmax

c.......................................................................
c.... Sort data by regions and set region vectors iregbeg and ireglen.

      call regsort( verbose, luverb,  nnobs,   kx,      kt, 
     $              rlats,   rlons,   rlevs,   del,
     &		    sigU,    sigO,    sigF,
     $              tstamp,  maxreg,  iregbeg, ireglen, ierr )

      if( ierr.ne.0 ) then
	 write(stderr,'(2a,i4)') 'sort',
     &	    ': error return from regsort(), status = ',ierr
	 call die('sort')
      endif
      write(mpout,'(2a)') myname,': returnd, regsort()'

c.......................................................................
c.... Sort data in each region by kt-kx-rlats-rlons-rlevs and set
c.... the array ityplen.


      call typsort( verbose, luverb,  nnobs,   kx,      kt, 
     $              rlats,   rlons,   rlevs,   del,
     &		    sigU,    sigO,    sigF,
     $              tstamp,  maxreg,  iregbeg, ireglen, ktmax, 
     $              ityplen                                            )
      write(mpout,'(2a)') myname,': returnd, typsort()'

c.......................................................................
c.... Print summary of data in all regions.

!      do 300 idreg = 1, maxreg
!         call prtregn( verbose, luverb,  idreg,   .false.,
!     $                 maxreg,  iregbeg, ireglen, ityplen,
!     $                 nnobs,   kx,      kt, 
!     $                 rlats,   rlons,   rlevs,
!     $                 del,     sigO,    sigF,    tstamp )
!  300 continue
!      write(mpout,'(2a)') myname,': returnd, prtregn()'
c.......................................................................
c.... Print data in selected regions.
c
c     call prtregn( verbose, luverb,       34,  .true.,
c    $              maxreg,  iregbeg, ireglen, ityplen,
c    $              nnobs,   kx,      kt, 
c    $              rlats,   rlons,   rlevs,
c    $              del,     sigO,    sigF,    tstamp )

	write(mpout,'(2a)') myname,'*regional_summary::'
	write(mpout,'(a)') '# kr    lat              lon'
	do kr=1,maxreg
	  if(ireglen(kr).gt.0) then
	    rlatmin=rlats(iregbeg(kr))
	    rlatmax=rlats(iregbeg(kr))
	    rlonmin=rlons(iregbeg(kr))
	    rlonmax=rlons(iregbeg(kr))
	    do i=1,ireglen(kr)-1
	      if(rlatmin.gt.rlats(iregbeg(kr)+i))
     &		rlatmin=rlats(iregbeg(kr)+i)
	      if(rlatmax.lt.rlats(iregbeg(kr)+i))
     &		rlatmax=rlats(iregbeg(kr)+i)

	      if(rlonmin.gt.rlons(iregbeg(kr)+i))
     &		rlonmin=rlons(iregbeg(kr)+i)
	      if(rlonmax.lt.rlons(iregbeg(kr)+i))
     &		rlonmax=rlons(iregbeg(kr)+i)
	    end do

	    write(mpout,'(i4,4f8.2)')
     &		kr,rlatmin,rlatmax,rlonmin,rlonmax
	  endif
	end do
	write(mpout,'(2a)') '::'

      return
      end

      subroutine sort_all(expid,verbose,luverb,nnobs,rlats,rlons,rlevs,
     .	kx,kt,
     .                    kid,maxreg,ktmax,iregbeg,ireglen,ityplen)
c******************************************************************************
c*                                                                            *
c*                                                                            *
c*                                                                            *
c*                                                                            *
c******************************************************************************
      use m_die,only : die

      character*4 myname
      parameter(myname='sort_all')

      character*8   expid

      real         rlats(nnobs)
      real         rlons(nnobs)
      real         rlevs(nnobs)
      integer      kx(nnobs)
      integer      kt(nnobs)
      integer      kid(nnobs)

      integer      iregbeg(maxreg)
      integer      ireglen(maxreg)
      integer      ityplen(ktmax,maxreg)

      logical      verbose

      !t1=second(xxxx)

c-----Sort data by regions and set region vectors iregbeg and ireglen.
      call regsort_all(verbose,luverb,nnobs,kx,kt,rlats,rlons,rlevs,
     .                 kid,maxreg,iregbeg,ireglen,ierr)

      if(ierr.ne.0) then
      write(*,'(2a,i4)') 'sort',
     .	': error return from regsort(), status = ',ierr
      call die('sort')
      endif
      write(*,'(2a)') myname,': returned, regsort_all()'

c-----Sort data in each region by kt-kx-rlats-rlons-rlevs and set ityplen.
      call typsort_all(verbose,luverb,nnobs,kx,kt,rlats,rlons,rlevs,
     .                 kid,maxreg,iregbeg,ireglen,ktmax,ityplen)
      write(*,'(2a)') myname,': returned, typsort_all()'

c-----Print summary of data in all regions.
      write(*,'(2a)') myname,'*regional_summary::'
      write(*,'(a)') '# kr    lat              lon'

      do kr=1,maxreg
      if(ireglen(kr).gt.0) then
      rlatmin=rlats(iregbeg(kr))
      rlatmax=rlats(iregbeg(kr))
      rlonmin=rlons(iregbeg(kr))
      rlonmax=rlons(iregbeg(kr))
      do i=1,ireglen(kr)-1
      if(rlatmin.gt.rlats(iregbeg(kr)+i)) rlatmin=rlats(iregbeg(kr)+i)
      if(rlatmax.lt.rlats(iregbeg(kr)+i)) rlatmax=rlats(iregbeg(kr)+i)
      if(rlonmin.gt.rlons(iregbeg(kr)+i)) rlonmin=rlons(iregbeg(kr)+i)
      if(rlonmax.lt.rlons(iregbeg(kr)+i)) rlonmax=rlons(iregbeg(kr)+i)
      enddo
      write(*,'(i4,4f8.2,i8)') kr,rlatmin,rlatmax,rlonmin,rlonmax,
     &                      ireglen(kr)
      endif
      enddo

      write(*,'(2a)') '::'
      !t2=second(xxxx)

      !write(*,*) 'tnet in sort_all:',t2-t1

      return
      end
