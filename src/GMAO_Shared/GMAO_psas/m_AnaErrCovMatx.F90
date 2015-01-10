!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_AnaErrCovMatx - Apply Pa to a vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_AnaErrCovMatx
      use m_CGSolver,only : CGSolver
      use m_FcstErrCovMatx, only : FcstErrCovMatx
      implicit none
      private	! except

      public :: AnaErrCovMatx
      public :: AnaErrCovMatx_init
      public :: AnaErrCovMatx_Cx
      public :: AnaErrCovMatx_clean
      public :: AnaErrCovMatx_trPa
      public :: nbandPa

      interface AnaErrCovMatx_init; module procedure init_; end interface
      interface AnaErrCovMatx_Cx  ; module procedure	&
	opPa_; end interface
      interface AnaErrCovMatx_trPa; module procedure trPa_; end interface
      interface AnaErrCovMatx_clean; module procedure clean_; end interface
      interface nbandPa; module procedure nband_; end interface

      type AnaErrCovMatx
     	private
     	integer :: nband
        type(CGSolver) :: pcg
        type(FcstErrCovMatx) :: pfcov
        type(FcstErrCovMatx) :: hpcov
        type(FcstErrCovMatx) :: phcov
      end type AnaErrCovMatx

! !REVISION HISTORY:
! 	28Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!	28Aug01 - Todling - adapted to Pa
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_AnaErrCovMatx'

#include "assert.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize a solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(pa,ob_corObs,ob_corPhi,ob_corPsi, &
                                  ai_corPhi,ai_corPsi, &
                     root,comm,rsrc)
      use m_die,    only : die
      use m_psasrc, only : psasrc_open,psasrc_close
      use m_CGSolver,only : CGSolver_init
      use m_CGSolver,only : nband
      use m_FcstErrCovMatx, only : FcstErrCovMatx_init
      use m_CorAttrF,only : CorAttrF
      use m_CorAttrX,only : CorAttrX

      implicit none
      type(AnaErrCovMatx),intent(out) :: pa
      type(CorAttrX),intent(in) :: ob_corObs	! C^o_c & C^o_u
      type(CorAttrF),intent(in) :: ob_corPhi	! C^phi
      type(CorAttrF),intent(in) :: ob_corPsi	! C^psi, C^chi
      type(CorAttrF),intent(in) :: ai_corPhi	! C^phi
      type(CorAttrF),intent(in) :: ai_corPsi	! C^psi, C^chi
      integer,intent(in) :: root
      integer,intent(in) :: comm
      character(len=*),optional,intent(in) :: rsrc

! !REVISION HISTORY:
! 	28Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!	28Aug01 - Todling - adapted to Pa
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'

  call CGSolver_init(pa%pcg,ob_corObs,ob_corPhi,ob_corPhi, &
        root,comm,rsrc=rsrc)

  pa%nband=nband(pa%pcg)

  call FcstErrCovMatx_init(pa%pfcov,ai_corPhi,ai_corPsi,   &
     	comm,pa%nband)

  call FcstErrCovMatx_init(pa%phcov,ai_corPhi,ai_corPsi,   &
     	ob_corPhi,ob_corPsi,comm,pa%nband)

  call FcstErrCovMatx_init(pa%hpcov,ob_corPhi,ob_corPsi,   &
     	ai_corPhi,ai_corPsi,comm,pa%nband)


end subroutine init_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nband_ - top band number
!
! !DESCRIPTION:
!
! !INTERFACE:

    function nband_(pa)
      implicit none
      type(AnaErrCovMatx),intent(in) :: pa
      integer :: nband_

! !REVISION HISTORY:
!       30May00 - Jing Guo <guo@dao.gsfc.nasa.gov>
!               - initial prototype/prolog/code
!	28Aug01 - Todling - adapted to Pa
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nband_'
  nband_ = pa%nband

end function nband_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean the solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(pa)
      use m_CGSolver,only : CGSolver_clean
      use m_FcstErrCovMatx, only : FcstErrCovMatx_clean

      implicit none
      type(AnaErrCovMatx),intent(inout) :: pa

! !REVISION HISTORY:
! 	23Aug01	- Todling- Initial code, based on J.Guo's CGSolver
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  pa%nband=-1
  call CGSolver_clean(pa%pcg)

  call FcstErrCovMatx_clean(pa%phcov)
  call FcstErrCovMatx_clean(pa%hpcov)
  call FcstErrCovMatx_clean(pa%pfcov)

end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: opPa_ - Operator to apply Pa to a vector
!
! !DESCRIPTION:
!
! Applies Pa to a vector:
!
!	$P^a y = [P^f - P^f*H^T(HP^fH^T + R)^{-1}*H*P^f] y$
!
! As an option the user may obtain the result of inc(Pa) y, where
! inc(Pa) is defined as
!
!	$inc(P^a) \equiv P^f*H^T(HP^fH^T + R)^{-1}*H*P^f$
!
!
! !INTERFACE:

	subroutine opPa_(  pa,root,comm,		&
           ob_Nav,ob_krNav,ob_ktNav,			&
           lobs,ob_sigU,ob_sigC,ob_corObs,	        &
        ob_ktab,ob_jtab,ob_sigPhi,ob_corPhi,           	&
	  ai_Nav,ai_krNav,ai_ktNav,                     &
     linc,ai_ktab,ai_jtab,ai_sigPhi,ai_corPhi,		& 
             ai_xx, ai_yy, what				)

        use m_stdio
        use m_zeit
	use m_die,   only : MP_die, die

	use m_sigAdef

        use m_CorAttrF,only : CorAttrF
        use m_CorAttrX,only : CorAttrX
        use m_sigmaPhi,only : sigmaPhi

        use m_FcstErrCovMatx,only : FcstErrCovMatx_Cx
        use m_CGSolver,only : CGSolver_solve
        use m_CGSolver,only : nband

        use m_Navigator,only : Navigator

        use m_mpif90,only : MP_comm_rank
        use m_mpif90,only : MP_type

        use m_mpout,only : mpout,mpout_ison,mpout_log

	implicit none

!     ..Region-decomposition maps
!     ============================
      type(AnaErrCovMatx),intent(inout) :: pa
      integer ,intent(in) :: root
      integer ,intent(in) :: comm

      type(Navigator),intent(in) :: ob_Nav
      integer,dimension(:),intent(in) :: ob_krNav
      integer,dimension(:),intent(in) :: ob_ktNav

      integer ,intent(in) :: lobs		! local no. of observations

      real   ,dimension(:),intent(in) :: ob_sigU
      real   ,dimension(:),intent(in) :: ob_sigC
      type(CorAttrX),intent(in) :: ob_corObs

      integer,dimension(:),intent(in) :: ob_ktab
      integer,dimension(:),intent(in) :: ob_jtab

      type(sigmaPhi),intent(in) :: ob_sigPhi
      type(CorAttrF),intent(in) :: ob_corPhi
!?    type(CorAttrF),intent(in) :: ob_corPsi


      type(Navigator),intent(in) :: ai_Nav
      integer,dimension(:),intent(in) :: ai_krNav
      integer,dimension(:),intent(in) :: ai_ktNav

      integer,intent(in) :: linc
      integer,dimension(:),intent(in) :: ai_ktab
      integer,dimension(:),intent(in) :: ai_jtab
      type(sigmaPhi),intent(in) :: ai_sigPhi
      type(CorAttrF),intent(in) :: ai_corPhi


!     ..Major quantities 
!     ==================

	  integer, intent(in), optional :: what ! operator defined according to
                                                !  what=-1 ==> op(Pf) 
                                                !  what= 0 ==> op(inc(Pa)) 
                                                !  what= 1 ==> op(Pa) default

	  real,dimension(:,:),intent(in)  :: ai_xx !  input  vector
	  real,dimension(:,:),intent(out) :: ai_yy !  output vector


! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
! !REMARKS:
!
! !REVISION HISTORY:
!       28Aug01 - Todling - Initial code.
!EOP ___________________________________________________________________

!   Dynamic memory allocation

        real,    allocatable :: ob_xx(:,:)      ! workspace
        real,    allocatable :: ob_yy(:,:)      ! workspace
        real,    allocatable :: ai_zz(:,:)      ! workspace


!-----------------------------------------------------------------------
	character(len=*), parameter :: myname_=myname//'::opPa_'

!	..Locals
	integer nvecs   ! no. of vectors in ai_xx
	integer ivec	! vector index
	integer ierr	! status from subroutines
	integer loc	! location pointer
	integer do_pa   ! determine which operator to apply; default=Pa

!       Preliminary setting/checking
!       ----------------------------
        if(lobs==0)return
        if(linc==0)return

        do_pa = 1
        if ( present(what) ) then
             do_pa = what
        end if

        call zeit_ci(myname_)

!       %-------------------------------------------------%
!       | Perform following matrix vector multiplication  |
!       |               y <--- Pf*x                       |
!       %-------------------------------------------------%

          if(do_pa/=0) & 
	  call FcstErrCovMatx_Cx(pa%pfcov,comm,         	        &
	        ai_Nav,ai_krNav,ai_ktNav,                      		&
	          linc,ai_ktab,ai_jtab,ai_sigPhi,ai_corPhi,ai_corPhi,   &
	        ai_xx,ai_yy)

          if (do_pa==-1) then
              call zeit_co(myname_)
              return   ! done calculating Pf*x
          end if

!     %---------------------%
!     | Allocate workspace  |
!     %---------------------%
       nvecs = size(ai_xx,2)
       allocate ( ob_xx(lobs,nvecs), ob_yy(lobs,nvecs), stat=ierr )
        if(ierr.ne.0) &
          call MP_die(myname_,': allocate(work) error, stat =',ierr)


     
!       %-------------------------------------------------%
!       | Perform following matrix vector multiplication  |
!       |               v <--- H*Pf*x                     |
!       %-------------------------------------------------%
 
 	  call FcstErrCovMatx_Cx(pa%hpcov,comm,               		&
	        ob_Nav,ob_krNav,ob_ktNav,                       	&
	          lobs,ob_ktab,ob_jtab,ob_sigPhi,ob_corPhi,ob_corPhi,   &
	        ai_Nav,ai_krNav,ai_ktNav,                       	&
	          linc,ai_ktab,ai_jtab,ai_sigPhi,ai_corPhi,ai_corPhi,   &
	        ai_xx,ob_xx)

!       %------------------------------------------------%
!       | Perform matrix vector multiplication           |
!       |         w <--- (H*Pf*H'+R)^{-1}*v              |
!       %------------------------------------------------%
 
	  call CGSolver_solve(pa%pcg,comm,			&
		ob_Nav,ob_krNav,ob_ktNav,			&
		lobs,ob_sigU,ob_sigC,ob_corObs,			&
		ob_ktab,ob_jtab,ob_sigPhi,ob_corPhi,ob_corPhi,	&
		ob_xx,ob_yy,ierr				)

        	if(ierr/=0) call die(myname_,'CGSolver_solve()',ierr)

!       %-------------------------------------------------%
!       | Perform following matrix vector multiplication  |
!       |               z <--- Pf*H'*w                    |
!       %-------------------------------------------------%

	  allocate ( ai_zz(linc,nvecs), stat=ierr )
	        if(ierr.ne.0) &
          	   call MP_die(myname_,': allocate(ai_zz) error, stat =',ierr)

	  call FcstErrCovMatx_Cx(pa%phcov,comm,         	        &
	        ai_Nav,ai_krNav,ai_ktNav,                      		&
	          linc,ai_ktab,ai_jtab,ai_sigPhi,ai_corPhi,ai_corPhi,   &
	        ob_Nav,ob_krNav,ob_ktNav,                      	 	&
	          lobs,ob_ktab,ob_jtab,ob_sigPhi,ob_corPhi,ob_corPhi,   &
	        ob_yy,ai_zz)

!       %---------------------------%
!       | Finally update y:= y - z  |
!       %---------------------------%

          if (do_pa==1) then 
	      ai_yy = ai_yy - ai_zz
	  else if (do_pa==0) then
	      ai_yy =         ai_zz
          else
              call MP_die(myname_,': should never get here error',99)
	  endif

	  deallocate ( ai_zz, stat=ierr )
	        if(ierr.ne.0) &
          	   call MP_die(myname_,': deallocate(ai_zz) error, stat =',ierr)

      call zeit_co(myname_)

      deallocate ( ob_xx, ob_yy, stat=ierr )
	    if(ierr.ne.0) &
                call MP_die(myname_,': deallocate(work) error, stat =',ierr)

      end subroutine opPA_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: trPa_ - Estimates trace of Pa.
!
! !DESCRIPTION:
!
!  This routine uses Girard (1989) approach to estimate the trace of
!  the analysis error covariance matrix.
!
! !INTERFACE:

	subroutine trPa_(  pa,root,comm,		&
           ob_Nav,ob_krNav,ob_ktNav,			&
           lobs,ob_sigU,ob_sigC,ob_corObs,	        &
        ob_ktab,ob_jtab,ob_sigPhi,ob_corPhi,           	&
	  ai_Nav,ai_krNav,ai_ktNav,                     &
     linc,ai_ktab,ai_jtab,ai_sigPhi,ai_corPhi,		& 
      	 nst, trace, what, stat				)

        use m_stdio
        use m_zeit
	use m_die,   only : MP_die, die

	use m_sigAdef

        use m_CorAttrF,only : CorAttrF
        use m_CorAttrX,only : CorAttrX
        use m_sigmaPhi,only : sigmaPhi

        use m_FcstErrCovMatx,only : FcstErrCovMatx_Cx
        use m_CGSolver,only : CGSolver_solve
        use m_CGSolver,only : nband

        use m_Navigator,only : Navigator

        use m_mpif90,only : MP_comm_rank
        use m_mpif90,only : MP_type
        use m_mpif90,only : MP_SUM

        use m_mpout,only : mpout,mpout_ison,mpout_log

        use m_parDOT,only : parDOT

	implicit none

!     ..Region-decomposition maps
!     ============================
      type(AnaErrCovMatx),intent(inout) :: pa
      integer ,intent(in) :: root
      integer ,intent(in) :: comm

      type(Navigator),intent(in) :: ob_Nav
      integer,dimension(:),intent(in) :: ob_krNav
      integer,dimension(:),intent(in) :: ob_ktNav

      integer ,intent(in) :: lobs		! local no. of observations

      real   ,dimension(:),intent(in) :: ob_sigU
      real   ,dimension(:),intent(in) :: ob_sigC
      type(CorAttrX),intent(in) :: ob_corObs

      integer,dimension(:),intent(in) :: ob_ktab
      integer,dimension(:),intent(in) :: ob_jtab

      type(sigmaPhi),intent(in) :: ob_sigPhi
      type(CorAttrF),intent(in) :: ob_corPhi
!?    type(CorAttrF),intent(in) :: ob_corPsi


      type(Navigator),intent(in) :: ai_Nav
      integer,dimension(:),intent(in) :: ai_krNav
      integer,dimension(:),intent(in) :: ai_ktNav

      integer,intent(in) :: linc
      integer,dimension(:),intent(in) :: ai_ktab
      integer,dimension(:),intent(in) :: ai_jtab
      type(sigmaPhi),intent(in) :: ai_sigPhi
      type(CorAttrF),intent(in) :: ai_corPhi

      integer, intent(in), optional :: what ! operator defined according to
                                            !  what=-1 ==> op(Pf) 
                                            !  what= 0 ==> op(inc(Pa)) 
                                            !  what= 1 ==> op(Pa) default


!     ..Output quantities
!     ===================


	integer, intent(out) :: nst	! no. samples used in estimate of trace
	real   , intent(out) :: trace	! trace estimte
	integer, intent(out), optional :: stat	! error code


! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
! !REMARKS:
!
! !REVISION HISTORY:
!       28Aug01 - Todling - Initial code.
!EOP ___________________________________________________________________

        integer  myID

!   Dynamic memory allocation

	real,allocatable :: ai_xx(:,:)
	real,allocatable :: ai_yy(:,:)
	real,allocatable ::    xx(:)

!-----------------------------------------------------------------------
	character(len=*), parameter :: myname_=myname//'::trPa_'

        include "trPa.h"

!	..Locals
	integer :: minc		! global no. of grid points

	integer i, j, k, jj	! data indeces
	integer ivec		! vector index
	integer ierr		! status from subroutines
	integer ns_trace	! location pointer
	real	alpha
	real	beta
        real	rerr1       	! rel error in randzed trace estimation
        real    tol_rte		! tolerance for randzed trace estimation
        real    fact		! normalization factor
        real    trace0		! exact trace of innov. cov/corr
        logical trace_ok


!     Preliminary setting/checking
!     ----------------------------
      if(present(stat)) stat=0
      nst  =0
      trace=0.
      if(linc==0)return

      call MP_comm_rank(comm,myID,ierr)
        if(ierr/=0) call MP_die(myname,'MP_comm_rank()',ierr)

      call MPI_allreduce( linc, minc, 1, MP_type(linc), MP_SUM, comm, ierr )
        if(ierr.ne.0) &
            call MP_die(myname_,': MPI_allreduce(linc) error, stat =',ierr)

       ns_trace = min(float(nst0),0.001*minc) ! Guess for no. samples required to estimate trace

!     %---------------------%
!     | Allocate workspace  |
!     %---------------------%
       allocate ( ai_xx(linc,ns_trace), ai_yy(linc,ns_trace), stat=ierr )
        if(ierr.ne.0) &
          call MP_die(myname_,': allocate(work) error, stat =',ierr)
     
        call zeit_ci(myname_)

!     %------------------------------------------------%
!     | RANDOMIZED TRACE ESTIMATION                    |
!     %------------------------------------------------%
!
       if(myID==ROOT) call zufalli ( 0 )   ! Initialize random number generator

       i = 0
       fact     = 0.0
       trace_ok = .false.
       tol_rte  = trace_tol    ! Force this trace estimate to be more accurate
                               ! than that determined by the evals
       do while ( .not. trace_ok )
         i   = i + 1
         if ( (.not. trace_ok)  .and.  i > ns_itrmax ) then
             write(stdout,'(3a,/,a)') myname_, ': WARNING --- ', &
                          'Maximum number of iterations reached ', &
                          '*** Convergence not acheived during trace estimation ***'
             trace_ok = .true.
             nst = i * ns_trace   ! this is just a guess at this point
             exit
         end if
         nst = ns_trace

!        %------------------------------------------------%
!        | Generate random vector x ~ N(0,I)              |
!        %------------------------------------------------%
           call pnormalen_ ( linc, nst, ai_xx, root, comm )

!        %--------------------------------------%
!        | Apply Pa to random N(0,1) vector     |
!        |         y <--- Pa*x                  |
!        %--------------------------------------%
 
 	 call opPa_(  pa,root,comm,			&
            ob_Nav,ob_krNav,ob_ktNav,			&
            lobs,ob_sigU,ob_sigC,ob_corObs,	        &
         ob_ktab,ob_jtab,ob_sigPhi,ob_corPhi,           &
	   ai_Nav,ai_krNav,ai_ktNav,                    &
      linc,ai_ktab,ai_jtab,ai_sigPhi,ai_corPhi,		& 
              ai_xx(:,1:nst), ai_yy, what=what		)

!        %-----------------%
!        | Estimate trace  |
!        %-----------------%
         rerr1 = trace
         do j = 1, nst
            jj    = (i -1)*nst + j
            alpha = (jj-1.)/jj
            beta  =     1. /jj
            fact  = alpha * fact   +  beta * parDOT(ai_xx(:,j),ai_xx(:,j),comm)
            trace = alpha * trace  +  beta * parDOT(ai_xx(:,j),ai_yy(:,j),comm)
         end do
         if (trace<=0. .or. fact<0) then
             if(present(stat)) stat = 99
             deallocate ( ai_xx, ai_yy, stat=ierr )
	     if(ierr.ne.0) &
                call MP_die(myname_,': deallocate(work) error, stat =',ierr)
             return
         end if
         trace = minc * trace / fact

         rerr1 = abs( rerr1 - trace  ) / trace
         if(myID==ROOT) &
         write(stdout,'(a,1p,2(a,e13.6))') myname_, &
                      ' rerr1 = ', rerr1, ' trace = ', trace
         if ( rerr1 < tol_rte ) then
             if(myID==ROOT) &
             write(stdout,'(2a,i8,a)') myname_, &
                    ': Determined trace in ', i, ' iterations '
             trace_ok = .true.
             nst = i * ns_trace
         end if

       end do ! < endwhile >
       if(myID==ROOT) &
       write(stdout,'(2a,1p,e11.4,a,i10)') myname_, &
                     ': Trace estimation result: ', trace, &
                     '  number of rnd-vecs used: ', nst

       if(trace.lt.1.e-10) then
          if(myID==ROOT) &
          write(stderr,'(2a,1p,e11.4)') myname_, &
            ': potential problem, trace seems too small =',trace
          if(present(stat)) stat = 98
       endif

      call zeit_co(myname_)

      deallocate ( ai_xx, ai_yy, stat=ierr )
	    if(ierr.ne.0) &
                call MP_die(myname_,': deallocate(work) error, stat =',ierr)

      end subroutine trPa_

      subroutine pnormalen_ ( linc, nvecs, xloc, root, comm )

      use m_mpif90,only : MP_comm_rank

      use m_Collector,only : Collector,Collector_init,clean,globalSize
      use m_CollectorComm,only : scatterv
      
      use m_die,only : MP_die
      implicit none

      integer, intent(in) :: linc
      integer, intent(in) :: nvecs
      real   , dimension(:,:), intent(out) :: xloc
      integer, intent(in) :: root
      integer, intent(in) :: comm

      character(len=*), parameter :: myname_ = myname//'::pnormalen_'
      integer :: minc, myID
      integer    ierr, i
      real   , allocatable :: xx(:)

      type(Collector) :: coll

      call MP_comm_rank(comm,myID,ierr)
        if(ierr/=0) call MP_die(myname,'MP_comm_rank()',ierr)

      call Collector_init(coll,linc,comm)
      minc = globalSize(coll)

!  Allocate memory for global array
!  -------------------------------- 
      if(myID==ROOT)then
         allocate ( xx(minc), stat=ierr )
      else
         allocate ( xx(   0), stat=ierr )
      end if
        if(ierr.ne.0) &
            call MP_die(myname_,': allocate(xx) error, stat =',ierr)

!  Create nvecs-Normally distributed random vectors and distribute them
!  --------------------------------------------------------------------
      do i = 1, nvecs
         if(myID==ROOT) call normalen ( minc, xx )

         call scatterv( xx, xloc(1:linc,i), coll, root, comm, ierr )
            if(ierr/=0) call MP_die(myname_,': MPI_scatter(xx) error, stat =',ierr)
      end do

!  Release memory
!  --------------
      deallocate ( xx, stat=ierr )
         if(ierr.ne.0) &
            call MP_die(myname_,': deallocate(xx) error, stat =',ierr)

      call clean(coll)

      end subroutine pnormalen_

   end module m_AnaErrCovMatx
