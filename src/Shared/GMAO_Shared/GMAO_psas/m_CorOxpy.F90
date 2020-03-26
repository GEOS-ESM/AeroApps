!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_CorOxpy 
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_CorOxpy
      use m_costs, only : tune_cost_start
      use m_costs, only : tune_cost_stop
      use m_costs, only : tuneQ

      implicit none
      private	! except

      public :: symCxpy		! y=Cx+y for a symmetric matrix
      public :: recCxpy		! y=Cx+y for a symmetric matrix

      interface symCxpy; module procedure mt_symCxpy_; end interface
      interface recCxpy; module procedure mt_recCxpy_; end interface

! !REVISION HISTORY:
!       14Nov01 - Tom Clune <clune@sgi.com>
!               . Distributes blocks under OpenMP according to ithd attribute.
!       30Oct01 - Tom Clune <clune@sgi.com>
!                 . Uses tuneQ() to determine whether to do costs.
! 	23Mar99 - Jing Guo <guo@dao> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_CorOxpy'

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
	aNav, nSeg,krNav,ktNav, nv,kx,qr,kl,	&
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
      real   ,dimension(3,nv),intent(in) :: qr	!
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
	aNav, nSeg,krNav,ktNav, nv,kx,qr,kl,	&
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
	aNav, nSeg,krNav,ktNav, na,kx,qr,kl,		&
	nvecs,vNav, nv, x,y, istat	)

      use m_AttrVect,      only : AttrVect
      use m_AttrVect,      only : AttrVect_lsize => lsize
      use m_AttrVect,      only : AttrVect_indexIA => indexIA
      use m_AttrVect,      only : ptr_iAttr

      use m_symMatx,only : symMatx_irow
      use m_symMatx,only : symMatx_icol
      use m_symMatx,only : symMatx_ierr
#ifdef _OPENMP
      use m_symMatx,only : symMatx_ithd
#endif

      use m_Navigator,	   only : Navigator
      use m_Navigator,	   only : get
      use m_block_corOxpy, only : sCxpy
      use m_block_corOxpy, only : xCxpy
      use m_die  ,only : perr,assert_, die
      use m_stdio,only : stderr
      use m_zeit ,only : zeit_ci,zeit_co
      use m_mpout,only : mpout
      use m_mall, only : mall_ison, mall_mci, mall_mco

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
      real   ,dimension(3,na),intent(in) :: qr
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
! 	16Feb96 - J. Guo	- as sym_Cxpy()
! 	09Nov98 - Jing Guo <guo@thunder> -
!		. revised to utilize a m_subBlocks module for a more
!		  controllable load-balancing scheme.
!		. revised the name to mt_sym_Cxpy().  The name sym_Cxpy
!		  is now reserved for the original interface definition.
!	03Feb99 - Jing Guo <guo@thunder> -
!		. Added share(nb_mA) to the multitasking directives.
! 	26Mar99 - Jing Guo <guo@dao>
!		- Developed based on mt_sym_Cxpy()
!	
!	31Oct01 - Tom Clune <clune@sgi.com>
!		. Removed array section notation for performance safety
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ser_symCxpy_'

  integer,parameter :: irow=symMatx_irow
  integer,parameter :: icol=symMatx_icol
  integer,parameter :: ierr=symMatx_ierr
#ifdef _OPENMP
  integer,parameter :: ithd=symMatx_ithd
#endif
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

#ifdef _OPENMP
  Integer, External :: OMP_GET_MAX_THREADS, OMP_GET_THREAD_NUM
  real, allocatable :: y_part(:,:,:)
  Integer :: minseg, maxseg, i_thread, n_threads
  Logical, Allocatable :: seg_thread(:,:)
#endif


!-----------------------------------------------------------------------
  istat=0 
  if(nvecs <= 0) return 
  if(lem < lcm)  return ! no assigned blocks

!-----------------------------------------------------------------------
	! Compute the partial products by blocks

    mAttr => ptr_iAttr(mNav)

#ifdef _OPENMP
        minseg = Minval(mAttr( (/irow,icol/), lcm:lem))
        maxseg = Maxval(mAttr( (/irow,icol/), lcm:lem))
        n_threads = OMP_GET_MAX_THREADS()

        Allocate(y_part(nvecs, nv, 1:n_threads-1), &
             &   seg_thread(minseg:maxseg, 1:n_threads-1), STAT = ier)
           If (ier /= 0) Call die(myname_,'allocate(y_part)',ier)
           If (mall_ison()) Then
              Call mall_mci(y_part, myname)
              Call mall_mci(seg_thread, myname)
           End If
           seg_thread = .false.
#endif

!$OMP PARALLEL &
!$OMP &     PRIVATE(ib,ir,ic,lci,lei,lni,kri,kti,lcvi,levi,lnvi) ,      &
!$OMP &     PRIVATE(lcj,lej,lnj,krj,ktj,lcvj,levj,lnvj,ier,i_thread),    &
!$OMP &     SHARED(lcm,lem,mAttr,aNav,krNav,ktNav,vNav,kind_cov,nvecs), &
!$OMP &     SHARED(qr,kx,kl,x,y,y_part,seg_thread,minseg,maxseg, n_threads)

#ifdef _OPENMP
   i_thread = OMP_GET_THREAD_NUM()
#endif

    do ib=lcm,lem	! for all block in the local segment

#ifdef _OPENMP
     If (mAttr(symMatx_ithd,ib) >= n_threads .or. mAttr(symMatx_ithd,ib) < 0 ) Then
        !$OMP MASTER
        mAttr(ierr,ib)=-9999 ! no such thread
        !$OMP END MASTER
        cycle
     End If

     If (mAttr(symMatx_ithd,ib) /= i_thread) Cycle ! Someone else should do this one
#endif

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
#ifdef _OPENMP
      If (i_thread > 0) Then	
        If (.not. seg_thread(ir,i_thread)) Then
           seg_thread(ir,i_thread) = .true.
           y_part(:,lcvi:levi,i_thread) = 0
        End If
        call sCxpy(kind_cov,		&
	  kri,kti,lni,			&
	  kx(  lci),		&
	  qr(1,lci),		&
	  kl(  lci),		&
	  nvecs,x(1,lcvi),		&
		y_part(1,lcvi,i_thread),		&
	  ier				)
      Else
#endif
        call sCxpy(kind_cov,		&
	  kri,kti,lni,			&
	  kx(  lci),		&
	  qr(1,lci),		&
	  kl(  lci),		&
	  nvecs,x(1,lcvi),		&
		y(1,lcvi),		&
	  ier				)
#ifdef _OPENMP
      End If
#endif

      else

		! Refer to the block-column attributes

	call get(aNav,ic,lc=lcj,le=lej,ln=lnj)
		krj=krNav(ic)
		ktj=ktNav(ic)

	call get(vNav,ic,lc=lcvj,le=levj,ln=lnvj)

		ASSERT(lnj==lnvj)

		! This block is _off_ the block-diagonal line


#ifdef _OPENMP
        If (i_thread > 0) Then
           If (.not. seg_thread(ir,i_thread)) Then
              seg_thread(ir,i_thread) = .true.
              y_part(:,lcvi:levi,i_thread) = 0
           End If
           If (.not. seg_thread(ic,i_thread)) Then
              seg_thread(ic,i_thread) = .true.
              y_part(:,lcvj:levj,i_thread) = 0
           End If
	call xCxpy(kind_cov,			&
&	  kri,kti,lni,				&
&	  kx(  lci),			&
&	  qr(1,lci),			&
&	  kl(  lci),			&
&	  krj,ktj,lnj,				&
&	  kx(  lcj),			&
&	  qr(1,lcj),			&
&	  kl(  lcj),			&
&	  nvecs,x(1,lcvi),y_part(1,lcvj,i_thread),	&
&		x(1,lcvj),y_part(1,lcvi,i_thread),	&
&	  ier					)
      Else
#endif
	call xCxpy(kind_cov,			&
	  kri,kti,lni,				&
	  kx(  lci),			&
	  qr(1,lci),			&
	  kl(  lci),			&
	  krj,ktj,lnj,				&
	  kx(  lcj),			&
	  qr(1,lcj),			&
	  kl(  lcj),			&
	  nvecs,x(1,lcvi),y(1,lcvj),	&
		x(1,lcvj),y(1,lcvi),	&
	  ier					)
#ifdef _OPENMP
      End If
#endif

   end if
      mAttr(ierr,ib)=ier

    end do
!     accumulate y_part into y
!     ------------------------
#ifdef _OPENMP
!$OMP BARRIER               ! ensure that partial sums are ready for combining
!$OMP DO SCHEDULE(STATIC)
  Do ir = minseg, maxseg
     Do i_thread = 1, n_threads-1
        If (seg_thread(ir,i_thread)) Then
           call get(vNav,ir,lc=lcvi,le=levi,ln=lnvi)
           y(:,lcvi:levi) = y(:,lcvi:levi) + y_part(:,lcvi:levi,i_thread)
        End If
     End Do
  End Do
!$OMP END DO

!$OMP END PARALLEL
        If (mall_ison()) Then
           Call mall_mco(y_part, myname)
           Call mall_mco(seg_thread, myname)
        End If
    Deallocate(seg_thread, y_part)
#endif
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

	write(stderr,'(2a,3i5,2i8)') myname_,	&
		': ir,kr,kt,lc,ln =',ir,kri,kti,lci,lni

	if(ic==ir) then
	  write(stderr,'(2a,i4,a)') myname_,	&
		': ic,kr,kt,lc,ln = ',ic,' [...]'

	else
	  call get(aNav,ic,lc=lcj,ln=lnj)
		krj=krNav(ic)
		ktj=ktNav(ic)

	  write(stderr,'(2a,3i5,2i8)') myname_,	&
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
		na_row,kx_row,qr_row,kl_row,	&
	aNav_col, nSeg_col,krNav_col,ktNav_col,	&
		na_col,kx_col,qr_col,kl_col,	&
	nvecs,x,y, stat)

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_lsize => lsize
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
      real   ,dimension(3,na_row),intent(in) :: qr_row
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
      real   ,dimension(3,na_col),intent(in) :: qr_col
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
  leMatx=AttrVect_lsize(mNav)

  call ser_recCxpy_(lcMatx,leMatx,mNav,kind_C,	&
	aNav_row,nSeg_row,krNav_row,ktNav_row,	&
		na_row,kx_row,qr_row,kl_row,	&
	aNav_col,nSeg_col,krNav_col,ktNav_col,	&
		na_col,kx_col,qr_col,kl_col,	&
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
		na_row,kx_row,qr_row,kl_row,		&
	aNav_col,nSeg_col,krNav_col,ktNav_col,		&
		na_col,kx_col,qr_col,kl_col,		&
  	nvecs,	xNav,nx,x,yNav,ny,y,			&
	istat)

      use m_AttrVect,only : AttrVect
      use m_AttrVect,only : lsize
      use m_AttrVect,only : ptr_iAttr

      use m_recMatx, only : recMatx_irow
      use m_recMatx, only : recMatx_icol
      use m_recMatx, only : recMatx_ierr
#ifdef _OPENMP
      use m_recMatx, only : recMatx_ithd
#endif

      use m_Navigator,	   only : Navigator
      use m_Navigator,	   only : get

      use m_block_corOxpy, only : rCxpy
      use m_die,           only : perr,assert_,die
      use m_stdio,         only : stderr
      use m_zeit,	   only : zeit_ci,zeit_co
      use m_mpout,only : mpout
      use m_mall, only : mall_ison, mall_mci, mall_mco

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
      real   ,dimension(3,na_row),intent(in) :: qr_row
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
      real   ,dimension(3,na_col),intent(in) :: qr_col
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
#ifdef _OPENMP
  integer,parameter :: ithd=recMatx_ithd
#endif

!=======================================================================
	! Local vars.

  integer :: ib,ir,ic
  integer :: kri,kti,lci,lei,lni,lcvi,levi,lnvi	! y
  integer :: krj,ktj,lcj,lej,lnj,lcvj,levj,lnvj	! x

  integer,pointer,dimension(:,:) :: mAttr
  integer :: ier	! error status


#ifdef _OPENMP
  Integer, External :: OMP_GET_MAX_THREADS, OMP_GET_THREAD_NUM
  Integer :: minseg, maxseg
  Integer :: n_threads, i_thread
  Real, Allocatable    :: y_part(:,:,:)
  Logical, Allocatable :: seg_thread(:,:)
#endif

!-----------------------------------------------------------------------
  if(nvecs <= 0) return

!-----------------------------------------------------------------------

	! Compute the partial products by blocks

!________________________________________

    mAttr => ptr_iAttr(mNav)

#ifdef _OPENMP
        minseg = Minval(mAttr( irow, lcm:lem))
        maxseg = Maxval(mAttr( irow, lcm:lem))
        n_threads = OMP_GET_MAX_THREADS()
        Allocate(y_part(nvecs, ny,1:n_threads-1), &
             &   seg_thread(minseg:maxseg, 1:n_threads-1), STAT = ier)
           If (ier /= 0) Call die(myname_,'allocate(y_part)',ier)
           seg_thread = .false.
#endif

!$OMP PARALLEL &
!$OMP &     PRIVATE(ib,ir,ic,lci,lei,lni,kri,kti,lcvi,levi,lnvi) ,      &
!$OMP &     PRIVATE(lcj,lej,lnj,krj,ktj,lcvj,levj,lnvj,ier, i_thread),    &
!$OMP &     SHARED(lcm,lem,mAttr,aNav_row,krNav_row,ktNav_row,yNav,kind_C,nvecs), &
!$OMP &     SHARED(qr_row,kx_row,kl_row,aNav_col,krNav_col,ktNav_col,xNav), &
!$OMP &     SHARED(qr_col,kx_col,kl_col,x,y,y_part,seg_thread,minseg,maxseg,n_threads)

#ifdef _OPENMP
   i_thread = OMP_GET_THREAD_NUM()
#endif

    do ib=lcm,lem
#ifdef _OPENMP
     If (mAttr(recMatx_ithd,ib) >= n_threads .or. mAttr(recMatx_ithd,ib) < 0 ) Then
        !$OMP MASTER
        mAttr(ierr,ib)=-9999 ! no such thread
        !$OMP END MASTER
        cycle
     End If

     If (mAttr(recMatx_ithd,ib) /= i_thread) Cycle ! Someone else should do this one
#endif

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

#ifdef _OPENMP
   If (i_thread > 0) Then
      If (.not. seg_thread(ir,i_thread)) Then
         seg_thread(ir,i_thread) = .true.
         y_part(1,lcvi:levi,i_thread) = 0
      End If
        call rCxpy(kind_C,		&
	  kri,kti, lni,			&
		kx_row(   lci),	&
		qr_row(1, lci),	&
		kl_row(   lci),	&
	  krj,ktj, lnj,			&
		kx_col(   lcj),	&
		qr_col(1, lcj),	&
		kl_col(   lcj),	&
	  nvecs,x(1,lcvj),		&
	  	y_part(1,lcvi,i_thread),&
	  ier)
   Else
#endif
        call rCxpy(kind_C,		&
	  kri,kti, lni,			&
		kx_row(   lci),	&
		qr_row(1, lci),	&
		kl_row(   lci),	&
	  krj,ktj, lnj,			&
		kx_col(   lcj),	&
		qr_col(1, lcj),	&
		kl_col(   lcj),	&
	  nvecs,x(1,lcvj),		&
	  	y(1,lcvi),		&
	  ier)
#ifdef _OPENMP
   End If  
#endif


      mAttr(ierr,ib)=ier

    end do
!     accumulate y_part into y
!     ------------------------
#ifdef _OPENMP
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC)
  Do ir = minseg, maxseg
     Do i_thread = 1, n_threads-1
        If (seg_thread(ir,i_thread)) Then
           call get(yNav,ir,lc=lcvi,le=levi,ln=lnvi)
           y(:,lcvi:levi) = y(:,lcvi:levi) + y_part(:,lcvi:levi,i_thread)
        End If
     End Do
  End Do
!$OMP END DO
!$OMP END PARALLEL
    Deallocate(seg_thread, y_part)
#endif
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

	write(stderr,'(2a,3i5,2i8)') myname_,	&
		': ir,kr,kt,lc,ln =',ir,kri,kti,lci,lni
	write(stderr,'(2a,3i5,2i8)') myname_,	&
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

end module m_CorOxpy
