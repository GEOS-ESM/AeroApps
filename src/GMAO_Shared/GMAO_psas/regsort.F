C  02Feb95  - Jing G.	- Changed CRAY to _UNICOS for consistency and
C			  to follow the guide lines.
c                                                   02/06/93 - regsort.f
c  modification for dynamic storage on CRAY         05/28/93
c  03oct94 - Implemented CRAY specifics with IFDEFs (A. da S.) 
c  ??????? - Todling - removed reference to common.blk
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.... Sort data by regions and set region vectors iregbeg and ireglen.

      subroutine REGSORT ( verbose, luverb,  nnobs,
     $                     kx,      kt,      rlats,    rlons,    rlevs,
     $                     del,     sigU,    sigO,     sigF,     tstamp,
     $                     maxreg,  iregbeg, ireglen,  ierr         )

        Use m_Spherical_Partition, only : xyz2reg
        use m_spherical_geometry, only : ll2xyz
	use m_die, only : die
        use m_mall,only : mall_ison,mall_mci,mall_mco
c.......................................................................
c.... Argument declarations.

      logical      verbose
      integer      luverb
      integer      nnobs
      integer      kx(nnobs)
      integer      kt(nnobs)
      real         rlats(nnobs)
      real         rlons(nnobs)
      real         rlevs(nnobs)
      real         del(nnobs)
      real         sigU(nnobs)
      real         sigO(nnobs)
      real         sigF(nnobs)
      real         tstamp(nnobs)
c-    real         qcflag(*)		! no longer used
      integer      maxreg
      integer      iregbeg(maxreg)
      integer      ireglen(maxreg)
      integer      ierr

c.......................................................................
c.... Local work space

	real,    allocatable :: x(:),y(:),z(:)
	integer, allocatable :: iregn(:),iperm(:)

	character*7 myname
	parameter  (myname='regsort')

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	allocate(x(nnobs),y(nnobs),z(nnobs),iregn(nnobs),iperm(nnobs),
     &	  stat=ierr )
	if(ierr.ne.0) call die(myname,'allocate()',ierr)
	if(mall_ison()) then
	  call mall_mci(x,myname)
	  call mall_mci(y,myname)
	  call mall_mci(z,myname)
	  call mall_mci(iregn,myname)
	  call mall_mci(iperm,myname)
	endif

c.......................................................................
c.... Compute (x,y,z) coordinates of data points from (rlons,rlats).

      call LL2XYZ ( rlons, rlats, nnobs, x, y, z, ierr )

      if( ierr.ne.0 ) then
         print  *, ' regsort: error return from ll2xyz'
		if(mall_ison()) then
		  call mall_mco(x,myname)
		  call mall_mco(y,myname)
		  call mall_mco(z,myname)
		  call mall_mco(iregn,myname)
		  call mall_mco(iperm,myname)
		endif
	 deallocate(x,y,z,iregn,iperm)
         return
      endif

c.......................................................................
c.... Compute regn from (x,y,z)

      call XYZ2REG ( nnobs, x, y, z, iregn )

c-      print *, ' regsort: returned from xyz2reg, nnobs = ', nnobs

c.......................................................................
c.... Sort by region.

      call INDEXXI ( nnobs, iregn, iperm )

c-    print *, ' regsort: returned from indexxi (sort on iregn)'
      print *, ' regsort: returned from indexxi, nnobs = ', nnobs

c.......................................................................
c.... Apply sorting permutation to data arrays.

      call PERMUTI ( kx,    iperm, nnobs, kx    )
      call PERMUTI ( kt,    iperm, nnobs, kt    )
      call PERMUTR ( rlats, iperm, nnobs, rlats )
      call PERMUTR ( rlons, iperm, nnobs, rlons )
      call PERMUTR ( rlevs, iperm, nnobs, rlevs )
      call PERMUTR ( del,   iperm, nnobs, del   )
      call PERMUTR ( sigU,  iperm, nnobs, sigU  )
      call PERMUTR ( sigO,  iperm, nnobs, sigO  )
      call PERMUTR ( sigF,  iperm, nnobs, sigF  )
      call PERMUTR ( tstamp,iperm, nnobs, tstamp)
c-    call PERMUTR ( qcflag,iperm, nnobs, qcflag)

c-     print *, ' regsort: returned from all calls to permuti & permutr'

c.......................................................................
c.... Fill iregbeg and ireglen

	ireglen(1:maxreg)=0
	do n=1,nnobs
	  ireglen(iregn(n))=ireglen(iregn(n))+1
	end do

	ksofar=1
	do ireg=1,maxreg
	  iregbeg(ireg)=ksofar
	  ksofar=ksofar+ireglen(ireg)
	end do

c-     print *, ' regsort: before returned to the parent'

	if(mall_ison()) then
	  call mall_mco(x,myname)
	  call mall_mco(y,myname)
	  call mall_mco(z,myname)
	  call mall_mco(iregn,myname)
	  call mall_mco(iperm,myname)
	endif
      deallocate(x,y,z,iregn,iperm)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end

      subroutine regsort_all(verbose,luverb,nnobs,kx,kt,rlats,rlons,
     .                       rlevs,
     .                       kid,maxreg,iregbeg,ireglen,ierr)
c******************************************************************************
c*                                                                            *
c*                                                                            *
c*                                                                            *
c*                                                                            *
c******************************************************************************

      implicit none

      integer luverb
      integer nnobs
      integer maxreg
      integer ierr

      real x(nnobs)
      real y(nnobs)
      real z(nnobs)
      real rlats(nnobs)
      real rlons(nnobs)
      real rlevs(nnobs)

      integer kx(nnobs)
      integer kt(nnobs)
      integer kid(nnobs)
      integer iregn(nnobs)
      integer iperm(nnobs)
      integer iregbeg(maxreg)
      integer ireglen(maxreg)

      logical verbose

      integer n, ireg, ksofar

c-----compute (x,y,z) coordinates of data points from (rlons,rlats)
      call ll2xyz(rlons,rlats,nnobs,x,y,z,ierr)

      if(ierr.ne.0) then
      write(*,*) ' regsort: error return from ll2xyz'
      return
      endif

c-----compute regn from (x,y,z)
      call xyz2reg(nnobs,x,y,z,iregn)

c-----sort by region
      call indexxi(nnobs,iregn,iperm)
      write(*,*) ' regsort: returned from indexxi, nnobs = ', nnobs

c-----apply sorting permutation to data arrays.
      call permuti(kid   ,iperm,nnobs,kid   )
      call permuti(kx,    iperm,nnobs,kx    )
      call permuti(kt,    iperm,nnobs,kt    )
      call permutr(rlats, iperm,nnobs,rlats )
      call permutr(rlons, iperm,nnobs,rlons )
      call permutr(rlevs, iperm,nnobs,rlevs )

c-----fill iregbeg and ireglen
      ireglen(1:maxreg)=0
      do n=1,nnobs
      ireglen(iregn(n))=ireglen(iregn(n))+1
      end do

      ksofar=1
      do ireg=1,maxreg
      iregbeg(ireg)=ksofar
      ksofar=ksofar+ireglen(ireg)
      end do

      return
      end
