!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_CorUxpy 
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_CorUxpy
      use m_costs, only : tune_cost_start
      use m_costs, only : tune_cost_stop
      use m_costs, only : tuneQ

      implicit none
      private	! except

      public :: symCxpy		! y=Cx+y for a symmetric matrix
      interface symCxpy; module procedure mt_symCxpy_; end interface

#ifdef _TODO_
      public :: recCxpy		! y=Cx+y for a symmetric matrix
      interface recCxpy; module procedure mt_recCxpy_; end interface
#endif

! !REVISION HISTORY:
!       30Oct01 - Tom Clune <clune@sgi.com>
!                 . Uses tuneQ() to determine whether to do costs.
! 	23Mar99 - Jing Guo <guo@dao> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_CorUxpy'

#include "assert.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: mt_symCxpy_ - An interface reserved for multitasking
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine mt_symCxpy_(mNav,kind_C,		&
	aNav, nSeg,krNav,ktNav, nv,kx,ks,kl,	&
	nvecs,x,y, stat)

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : lsize

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize

      use m_die,      only : die,perr,assert_
      implicit none

		! Matrix block attributes

      type(AttrVect),intent(inout) :: mNav	! a list of blocks
      integer,       intent(in   ) :: kind_C	! the kind of C

		! A navigator of attribute and vector blocks

      type(Navigator)        ,intent(in) :: aNav

		! Block attributes

      integer                ,intent(in) :: nSeg
      integer,dimension(nSeg),intent(in) :: krNav
      integer,dimension(nSeg),intent(in) :: ktNav

		! Data attributes

      integer                ,intent(in) :: nv	! the vector size
      integer,dimension(  nv),intent(in) :: kx	!
      integer,dimension(  nv),intent(in) :: ks	!
      integer,dimension(  nv),intent(in) :: kl	!

		! Vectors

      integer                 ,intent(in)    :: nvecs ! no. of vectors
      real,dimension(nvecs,nv),intent(in)    :: x
      real,dimension(nvecs,nv),intent(inout) :: y

      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	23Mar99 - Jing Guo <guo@dao> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::mt_symCxpy_'

  integer :: ier
  integer :: lcMatx
  integer :: leMatx

  if(present(stat)) stat=0
	ASSERT(nSeg == lsize(aNav))

	! a place holder of OpenMP implementation
  lcMatx=1
  leMatx=lsize(mNav)

	! local: lcMatx,leMatx,ier
	! reduction: y? maybe
	! shared: all others

  call ser_symCxpy_(lcMatx,leMatx,mNav,kind_C,	&
	aNav, nSeg,krNav,ktNav, nv,kx,ks,kl,	&
	nvecs, aNav, nv, x,y, ier)

  if(ier /= 0) then
    call perr(myname_,'ser_symCxpy_()',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

end subroutine mt_symCxpy_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ser_symCxpy - y=Cx+y, with a symmectric sparse matrix C.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine ser_symCxpy_(lcm,lem,mNav,kind_cov,	&
	aNav, nSeg,krNav,ktNav, na,kx,ks,kl,		&
	nvecs,vNav, nv, x,y, istat	)

      use m_AttrVect,      only : AttrVect
      use m_AttrVect,      only : AttrVect_lsize => lsize
      use m_AttrVect,      only : AttrVect_indexIA => indexIA
      use m_AttrVect,      only : ptr_iAttr

      use m_symMatx,only : symMatx_irow
      use m_symMatx,only : symMatx_icol
      use m_symMatx,only : symMatx_ierr

      use m_Navigator,	   only : Navigator
      use m_Navigator,	   only : get
      use m_block_corUxpy, only : sCxpy
      use m_die  ,only : die,perr,assert_
      use m_stdio,only : stderr
      use m_zeit ,only : zeit_ci,zeit_co
      use m_mpout,only : mpout, mpout_log

      implicit none

		! The local navigator segment of the matrix blocks

      integer,intent(in) :: lcm	! lbound of mNav blocks
      integer,intent(in) :: lem	! ubound of mNav blocks
      type(AttrVect),intent(inout) :: mNav	! matrix blocks
      integer,       intent(in)    :: kind_cov	! correlation type

		! The navigator the attribute blocks

      type(Navigator)        ,intent(in) :: aNav

		! The block attributes

      integer                ,intent(in) :: nSeg
      integer,dimension(nSeg),intent(in) :: krNav
      integer,dimension(nSeg),intent(in) :: ktNav

		! The attributes of the vector elements

      integer,intent(in) :: na	! size of vectors
      integer,dimension(  na),intent(in) :: kx
      integer,dimension(  na),intent(in) :: ks
      integer,dimension(  na),intent(in) :: kl

		! The input and output vectors

      integer,intent(in) :: nvecs		! nvecs of x and y

      type(Navigator),intent(in) :: vNav	! Navigator of x and y

      integer,intent(in) :: nv	! size of vectors
      real,dimension(nvecs,nv),intent(in)    :: x
      real,dimension(nvecs,nv),intent(inout) :: y

      integer,intent(out):: istat	! return status

		! Ajustable-size arrays are used instead of assumed-
		! shape arrays.  It allows the arrays to be used in
		! this procedure in traditional (F77) styles of
		! apperent continue storages.  It implies that
		! duplication processes may take place when this 
		! procedure is called.

! !REVISION HISTORY:
!
! 	26Mar99 - Jing Guo <guo@dao>
!		- Developed based on mt_sym_Cxpy()
!	
! 	16Feb96 - J. Guo	- as sym_Cxpy()
! 	09Nov98 - Jing Guo <guo@thunder> -
!		. revised to utilize a m_subBlocks module for a more
!		  controllable load-balancing scheme.
!		. revised the name to mt_sym_Cxpy().  The name sym_Cxpy
!		  is now reserved for the original interface definition.
!	03Feb99 - Jing Guo <guo@thunder> -
!		. Added share(nb_mA) to the multitasking directives.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ser_symCxpy_'

  integer,parameter :: irow=symMatx_irow
  integer,parameter :: icol=symMatx_icol
  integer,parameter :: ierr=symMatx_ierr

!=======================================================================
	! Local vars.

  integer :: kri,krj	! region index
  integer :: kti,ktj	! type index

  integer :: lci,lei,lni ! location indices of i- attribute vector block
  integer :: lcj,lej,lnj ! location indices of j- attribute vector block

  integer :: lcvi,levi,lnvi ! location indices of i- vector block
  integer :: lcvj,levj,lnvj ! location indices of j- vector block

  integer :: ir,ic	! indices of row and columns
  integer :: ier	! error status

  integer :: ib
  integer,pointer,dimension(:,:) :: mAttr


!-----------------------------------------------------------------------
  istat=0 
  if(nvecs <= 0) return 
  if(lem < lcm)  return ! no assigned blocks

!-----------------------------------------------------------------------
	! Compute the partial products by blocks

    mAttr => ptr_iAttr(mNav)
    ier = 0
    do ib=lcm,lem	! for all block in the local segment

		! Refer to the block-matrix attributes

      ir=mAttr(irow,ib)		! block-row index
      ic=mAttr(icol,ib)		! block-colum index

		! Refer to the block-row attributes

      call get(aNav,ir,lc=lci,le=lei,ln=lni)
	kri=krNav(ir)
	kti=ktNav(ir)

      call get(vNav,ir,lc=lcvi,le=levi,ln=lnvi)

	ASSERT(lni==lnvi)


      if(ir==ic) then
		! This block is _on_ the block-diagonal line


        call sCxpy(kind_cov,		&
	  kri,kti,lni,			&
	  kx(  lci:lei),		&
	  ks(  lci:lei),		&
	  kl(  lci:lei),		&
	  nvecs,x(:,lcvi:levi),		&
		y(:,lcvi:levi),		&
	  ier				)


      else

#ifdef _TODO_
		! Refer to the block-column attributes

	call get(aNav,ic,lc=lcj,le=lej,ln=lnj)
		krj=krNav(ic)
		ktj=ktNav(ic)

	call get(vNav,ic,lc=lcvj,le=levj,ln=lnvj)

		ASSERT(lnj==lnvj)

		! This block is _off_ the block-diagonal line


	call xCxpy(kind_cov,			&
	  kri,kti,lni,				&
	  kx(  lci:lei),			&
	  ks(  lci:lei),			&
	  kl(  lci:lei),			&
	  krj,ktj,lnj,				&
	  kx(  lcj:lej),			&
	  ks(  lcj:lej),			&
	  kl(  lcj:lej),			&
	  nvecs,x(:,lcvi:levi),y(:,lcvj:levj),	&
		x(:,lcvj:levj),y(:,lcvi:levi),	&
	  ier					)

#endif

      endif

      mAttr(ierr,ib)=ier

    end do
!_______________________________________________________________________

	! Verify the error code

    istat=0
    do ib=lcm,lem
      ier=mAttr(ierr,ib)
      if(ier/=0) then

	ir=mAttr(irow,ib)
	ic=mAttr(icol,ib)

	if(ic==ir) then
	  call perr(myname_,'sCxpy()',ier)
	else
	  call perr(myname_,'xCxpy()',ier)
	endif

	call get(aNav,ir,lc=lci,ln=lni)
		kri=krNav(ir)
		kti=ktNav(ir)

	write(stderr,'(2a,3i5,2i6)') myname_,	&
		': ir,kr,kt,lc,ln =',ir,kri,kti,lci,lni

	if(ic==ir) then
	  write(stderr,'(2a,i4,a)') myname_,	&
		': ic,kr,kt,lc,ln = ',ic,' [...]'

	else
	  call get(aNav,ic,lc=lcj,ln=lnj)
		krj=krNav(ic)
		ktj=ktNav(ic)

	  write(stderr,'(2a,3i5,2i6)') myname_,	&
		': ic,kr,kt,lc,ln =',ic,krj,ktj,lcj,lnj
	endif

	istat=max(abs(istat),abs(ier))
      endif
    end do

	nullify(mAttr)

    if(istat/=0) then
      call perr(myname_,'sCxpy()/xCxpy()',istat)
      return
    endif
!_______________________________________________________________________

end subroutine ser_symCxpy_

#ifdef _TODO_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: mt_recCxpy_ - A recCxpy interface reserved for multitasking
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine mt_recCxpy_(mNav,kind_C,		&
	aNav_row, nSeg_row,krNav_row,ktNav_row,	&
		na_row,kx_row,ks_row,kl_row,	&
	aNav_col, nSeg_col,krNav_col,ktNav_col,	&
		na_col,kx_col,ks_col,kl_col,	&
	nvecs,x,y, stat)

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : lsize
      use m_AttrVect, only : AttrVect_indexIA => indexIA
      use m_AttrVect, only : AttrVect_indexRA => indexRA

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize

      use m_die,      only : die,perr,assert_
      use m_stdio,    only : stderr
      implicit none

	! Matrix navigator

      type(AttrVect),intent(inout) :: mNav	! block navigator
      integer,       intent(in)    :: kind_C	! matrix type

	! Row attributes

			! Row navigator

      type(Navigator),intent(in) :: aNav_row

			! block attributes

      integer,intent(in) :: nSeg_row
      integer,dimension(nSeg_row),intent(in) :: krNav_row
      integer,dimension(nSeg_row),intent(in) :: ktNav_row

			! data attributes

      integer,intent(in) :: na_row
      integer,dimension(  na_row),intent(in) :: kx_row
      integer,dimension(  na_row),intent(in) :: ks_row
      integer,dimension(  na_row),intent(in) :: kl_row

			! Column navigator

      type(Navigator),intent(in) :: aNav_col

			! block attributes

      integer,intent(in) :: nSeg_col
      integer,dimension(nSeg_col),intent(in) :: krNav_col
      integer,dimension(nSeg_col),intent(in) :: ktNav_col

			! data attributes

      integer,intent(in) :: na_col
      integer,dimension(  na_col),intent(in) :: kx_col
      integer,dimension(  na_col),intent(in) :: ks_col
      integer,dimension(  na_col),intent(in) :: kl_col

      integer,intent(in) :: nvecs
      real,dimension(nvecs,na_col),intent(in)    :: x
      real,dimension(nvecs,na_row),intent(inout) :: y

      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	23Mar99 - Jing Guo <guo@dao> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::mt_recCxpy_'

  integer :: ier
  integer :: lcMatx,leMatx

  if(present(stat)) stat=0

	ASSERT(nSeg_row == lsize(aNav_row))
	ASSERT(nSeg_col == lsize(aNav_col))

  lcMatx=1
  leMatx=lsize(mNav)

  call ser_recCxpy_(lcMatx,leMatx,mNav,kind_C,	&
	aNav_row,nSeg_row,krNav_row,ktNav_row,	&
		na_row,kx_row,ks_row,kl_row,	&
	aNav_col,nSeg_col,krNav_col,ktNav_col,	&
		na_col,kx_col,ks_col,kl_col,	&
  	nvecs,	aNav_col,na_col,x,		&
		aNav_row,na_row,y,		&
	ier)

  if(ier /= 0) then
    call perr(myname_,'ser_rec_Cxpy_()',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

end subroutine mt_recCxpy_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ser_recCxpy - y=Cx+y, with a rectangular sparse matrix C.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine ser_recCxpy_(lcm,lem,mNav,kind_C,	&
	aNav_row,nSeg_row,krNav_row,ktNav_row,		&
		na_row,kx_row,ks_row,kl_row,		&
	aNav_col,nSeg_col,krNav_col,ktNav_col,		&
		na_col,kx_col,ks_col,kl_col,		&
  	nvecs,	xNav,nx,x,yNav,ny,y,			&
	istat)

      use m_AttrVect,only : AttrVect
      use m_AttrVect,only : lsize
      use m_AttrVect,only : ptr_iAttr

      use m_recMatx, only : recMatx_irow
      use m_recMatx, only : recMatx_icol
      use m_recMatx, only : recMatx_ierr

      use m_Navigator,	   only : Navigator
      use m_Navigator,	   only : get

      use m_block_corUxpy, only : rCxpy
      use m_die,           only : perr,assert_
      use m_stdio,         only : stderr
      use m_zeit,	   only : zeit_ci,zeit_co

      implicit none

	! Matrix navigator

      integer,intent(in) :: lcm,lem	! lbound and ubound of mNav
      type(AttrVect),intent(inout) :: mNav	! block navigator
      integer,intent(in) :: kind_C	! matrix type

	! Row attributes

			! Row navigator

      type(Navigator),intent(in) :: aNav_row

			! block attributes

      integer,intent(in) :: nSeg_row
      integer,dimension(nSeg_row),intent(in) :: krNav_row
      integer,dimension(nSeg_row),intent(in) :: ktNav_row

			! data attributes

      integer,intent(in) :: na_row
      integer,dimension(  na_row),intent(in) :: kx_row
      integer,dimension(  na_row),intent(in) :: ks_row
      integer,dimension(  na_row),intent(in) :: kl_row

			! Column navigator

      type(Navigator),intent(in) :: aNav_col

			! block attributes

      integer,intent(in) :: nSeg_col
      integer,dimension(nSeg_col),intent(in) :: krNav_col
      integer,dimension(nSeg_col),intent(in) :: ktNav_col

			! data attributes

      integer,intent(in) :: na_col
      integer,dimension(  na_col),intent(in) :: kx_col
      integer,dimension(  na_col),intent(in) :: ks_col
      integer,dimension(  na_col),intent(in) :: kl_col

      integer,intent(in) :: nvecs

      type(Navigator),intent(in) :: xNav
      integer,intent(in) :: nx
      real,dimension(nvecs,nx),intent(in)    :: x

      type(Navigator),intent(in) :: yNav
      integer,intent(in) :: ny
      real,dimension(nvecs,ny),intent(inout) :: y

      integer,intent(out) :: istat

      ! Ajustable-size arrays are used instead of assumed-
      ! shape arrays.  It allows the arrays to be used in
      ! this procedure in traditional (F77) styles of
      ! apperent continue storages.  It implies that
      ! duplication processes may take place when this 
      ! procedure is called.

! !REVISION HISTORY:
!
! 	26Mar99 - Jing Guo <guo@dao>
!		- Developed based on mt_rec_Cxpy()
!	
! 	16Feb96 - J. Guo	- as rec_Cxpy()
! 	09Nov98 - Jing Guo <guo@thunder> -
!		. revised to utilize a m_subBlocks module for a more
!		  controllable load-balancing scheme.
!		. revised the name to mt_rec_Cxpy().  The name rec_Cxpy
!		  is now reserved for the original interface definition.
!	03Feb99 - Jing Guo <guo@thunder> -
!		. Added share(nb_mA) to the multitasking directives.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ser_recCxpy_'

  integer,parameter :: irow=recMatx_irow
  integer,parameter :: icol=recMatx_icol
  integer,parameter :: ierr=recMatx_ierr

!=======================================================================
	! Local vars.

  integer :: ib,ir,ic
  integer :: kri,kti,lci,lei,lni,lcvi,levi,lnvi	! y
  integer :: krj,ktj,lcj,lej,lnj,lcvj,levj,lnvj	! x

  integer,pointer,dimension(:,:) :: mAttr
  integer :: ier	! error status


!-----------------------------------------------------------------------
  if(nvecs <= 0) return

!-----------------------------------------------------------------------

	! Compute the partial products by blocks

!________________________________________

    mAttr => ptr_iAttr(mNav)
    do ib=lcm,lem

		! References to the block-matrix attributes

      ir=mAttr(irow,ib)		! block-row index
      ic=mAttr(icol,ib)		! block-colum index

		! References to the block-row attributes
	!________________________________________

      call get(aNav_row,ir,lc=lci,le=lei,ln=lni)
	kri=krNav_row(ir)
	kti=ktNav_row(ir)

      call get(yNav,ir,lc=lcvi,le=levi,ln=lnvi)

	ASSERT(lni==lnvi)

		! References to the block-column attributes
	!________________________________________

      call get(aNav_col,ic,lc=lcj,le=lej,ln=lnj)
	krj=krNav_col(ic)
	ktj=ktNav_col(ic)

      call get(xNav,ic,lc=lcvj,le=levj,ln=lnvj)

	ASSERT(lnj==lnvj)

	!________________________________________

        call rCxpy(kind_C,		&
	  kri,kti, lni,			&
		kx_row(   lci: lei),	&
		ks_row(   lci: lei),	&
		kl_row(   lci: lei),	&
	  krj,ktj, lnj,			&
		kx_col(   lcj: lej),	&
		ks_col(   lcj: lej),	&
		kl_col(   lcj: lej),	&
	  nvecs,x(:,lcvj:levj),		&
	  	y(:,lcvi:levi),		&
	  ier)


      mAttr(ierr,ib)=ier


    end do
!_______________________________________________________________________

	! Verify the error code

    istat=0
    do ib=lcm,lem
      if(mAttr(ierr,ib) /= 0) then
	ier=mAttr(ierr,ib)
	ir =mAttr(irow,ib)
	ic =mAttr(icol,ib)

	call perr(myname_,'rCxpy()',ier)

	call get(aNav_row,ir,lc=lci,ln=lni)
		kri=krNav_row(ir)
		kti=ktNav_row(ir)

	call get(aNav_col,ic,lc=lcj,ln=lnj)
		krj=krNav_col(ic)
		ktj=ktNav_col(ic)

	write(stderr,'(2a,3i4,2i5)') myname_,	&
		': ir,kr,kt,lc,ln =',ir,kri,kti,lci,lni
	write(stderr,'(2a,3i4,2i5)') myname_,	&
		': ic,kr,kt,lc,ln =',ic,krj,ktj,lcj,lnj

	istat=max(abs(istat),abs(ier))
      endif
    end do

	nullify(mAttr)

    if(istat/=0) then
      call perr(myname_,'rCxpy(*)',istat)
      return
    endif
!_______________________________________________________________________

end subroutine ser_recCxpy_
#endif

end module m_CorUxpy
