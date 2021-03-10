!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_costs - matrix block costs
!
! !DESCRIPTION:
!
! !INTERFACE:

  module m_costs

      use config, only : ktmax
      implicit none

      private      ! except

      public :: costs_init      ! Cost function initialization
      public :: costs_clean     ! Cost function cleaning up
      public :: mb_cost         ! Matrix-block formation/access cost
				! function evaluation
      public :: tuneQ           ! Whether to compute tuning statistics

!      public :: MXcost_clas    !  The maximum number of cost classes.
!      public :: MXcost_type    !  The maximum number of cost types within a
                               !  given cost class

!      public :: compute_sparsity ! Compute Matrix Sparsity
      public :: tune_cost_start
      public :: tune_cost_stop

      public :: print_cost_statistics

      interface costs_init ; module procedure init_ ; end interface
      interface costs_clean; module procedure clean_; end interface

      interface mb_cost   ; module procedure &
         sym_mb_cost_,                       &
         rec_mb_cost_
      end interface

!      interface compute_sparsity   ; module procedure &
!         sym_compute_sparsity_ ,                      &
!         rec_compute_sparsity_
!      end interface

      interface tune_cost_start
	 module procedure tune_cost_start_
      end interface

      interface tune_cost_stop
	 module procedure tune_cost_stop_
      end interface

      interface print_cost_statistics
	 module procedure print_cost_statistics_
      end interface

! !REVISION HISTORY:
!        9Jul99 - J.W. Larson <jlarson@dao> - initial prototype/prolog/code
!        6Oct99 - J.W. Larson <jlarson@dao.gsfc.nasa.gov> - replaced 
!                 occurances of MP_REAL with MP_type(message) for 
!                 better support of 32 and 64-bit platforms.  This change
!                 affects the routines sym_compute_sparsity_() and
!                 rec_compute_sparsity_().
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_costs'

#include "assert.H"
#include "lapack.H"

  integer,save :: cost_class(ktmax,ktmax)

!   Relative costs.  The matrix costs(:,:) is arranged as follows
!        Rows are defined by cost classes
!        Columns define types of costs within an individual cost class.
!
!   The maximum number of cost classes is MXcost_clas:
  integer , parameter :: MXcost_clas = 4

!   The maximum number of cost types is MXcost_type:
  integer , parameter :: MXcost_type = 3

!   In the current implementation, there are four cost classes, one 
!   for each row of costs(:,:):
!          *  The first row is the cost to compute height vs. height 
!             (and thus SLP-SLP) forecast error correlations.
!          *  The second row is the cost to compute height vs. height
!             gradient (that is h vs. h_l and h vs. h_m and vice-versa)
!             forecast error correlations.  This cost also applies to 
!             the computation of PSL vs. PSL gradient correlations.
!          *  The third row is the cost to compute height gradient vs. 
!             height gradient correlations (that is h_l vs. h_l, h_l vs. 
!             h_m, et cetera).  This cost also applies to the computation
!             of PSL gradient vs. PSL gradient correlations.
!          *  The fourth row is the cost to compute univariate q vs. q
!             forecast error correlations.
!
!
! The following documentation is now incorrect.
!   In the current implementation, there are three cost types, one for
!   each column of costs(:,:):
!          *  The first column is the cost to form the block from 
!             nearest-neighbor table look-ups.
!          *  The second column is the cost to form the block from 
!             table look-ups and linear interpolation.
!          *  The third column is the memory bandwidth cost to access
!             the stored block from memory.

  real,parameter :: min_cost = 1.e-6

  real,save    :: costs(MXcost_clas,MXcost_type)
  integer,save,target :: poly_set(2,4) = reshape(source = (/ 0,0,  1,0,  0,1,  1,1 /), &
       &                                        shape = (/ 2, 4 /))
  integer,save,target :: poly_set_dia(2,3) = reshape(source = (/ 0,0,  1,0,  2,0 /), &
       &                                        shape = (/ 2, 3 /))

  real,save,target :: cost_coefs_rec(size(poly_set,2),           MXcost_clas, MXcost_type)
  real,save,target :: cost_coefs_sym(size(poly_set,2),           MXcost_clas, MXcost_type)
  real,save,target :: cost_coefs_dia(size(poly_set_dia,2), MXcost_clas, MXcost_type)

!    Debugging switch DEBUG:
  logical , parameter :: DEBUG = .false.
  logical :: tuning_is_on = .false. ! no tuning by default

  type statistics
     real, pointer :: s(:,:,:)
  end type statistics
  type (statistics), save :: stat_rec(MXcost_clas)
  type (statistics), save :: stat_sym(MXcost_clas)
  type (statistics), save :: stat_dia(MXcost_clas)

  ! o2k hardware counter stuff

#ifdef _OPENMP
  Integer, External :: OMP_GET_THREAD_NUM
  Integer, Parameter :: MAX_THREADS=100
  Integer, save :: start_count(16,0:MAX_THREADS-1)
#else
  Integer, save :: start_count
#endif
  Integer, save :: max_clock_tick
  Integer, save :: clock_rate
  Real,    save :: clock_rate_inv

  logical,save :: mb_cost_initialized=.false.
  
contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize cost data objects cost_class(:,:) 
!                    and costs(:,:).
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(comm,root)

    use config, only : ktmax
    use m_ioutil,  only : luavail, opnieee,clsieee
    use m_mpif90
    use m_mpout,only : mpout_log
    use m_die  ,only : die,MP_die,assert_
    use m_mall ,only : mall_ison,mall_mci

    implicit none
    integer,intent(in) :: comm
    integer,intent(in) :: root

! !REVISION HISTORY:
!       09Jul99 - J.W. Larson <jlarson@dao> - initial prototype/prolog/code
!       30Oct01 - T. Clune <clune@sgi.com> 
!                 . Added logical flag "tune" that is determined via environment variable.
!                 . Long ago added sophisticated cost tuning mechanisms
!EOP ___________________________________________________________________

  character (len=*) , parameter :: myname_=myname//'::init_'
  integer i, class, nn
  integer :: unit, ierr, myid
  character(len=128) :: infile='costs.in'
  character(len=255) :: env_val ! buffer for environment variables

  if (mb_cost_initialized)	&
	call die(myname_,'multiple object definition')

  mb_cost_initialized = .true.

!   Set cost classes

!   Initialize all kti, ktj classes to a default value of 4:
      cost_class(1:ktmax,1:ktmax) = 4

!   Set cost classes for sea-level puv analysis:
      cost_class(1:2,1:2) = 3
      cost_class(1:2,3) = 2
      cost_class(3,1:2) = 2
      cost_class(3,3) = 1

!   Set cost classes for upper-air huv analysis:
      cost_class(4:5,4:5) = 3
      cost_class(4:5,6) = 2
      cost_class(6,4:5) = 2
      cost_class(6,6) = 1

!     Initialize costs(:,:) to 1.0:

      costs(1:4,1:3) = 1.0
      cost_coefs_rec = 0
      cost_coefs_sym = 0

!     Cost per element to form matrix blocks from look-ups (based on 
!     the # fp ops per element):
      costs(1,1) = 11.0
      costs(2,1) = 14.5
      costs(3,1) = 22.5
      costs(4,1) = 12.0

      Call MP_comm_rank(comm,myid,ierr)

      If (myid == root) Then
	 unit = luavail()
	 call opnieee(unit,infile,"old",ierr)

	 if (ierr /= 0) then
	    call mpout_log(myname_,'Using default costs values')

	    cost_coefs_rec(1:4,1,1)=(/273266., -503., -1702.,  28. /) / 4.e8
	    cost_coefs_rec(1:4,2,1)=(/ 89044., -128.,  -662.,  29. /) / 4.e8
	    cost_coefs_rec(1:4,3,1)=(/ 48390.,  -41., -1817.,  38. /) / 4.e8
	    cost_coefs_rec(1:4,4,1)=(/273266., -503., -1702.,  28. /) / 4.e8

	    cost_coefs_sym(1:4,1,1)=	&
		(/ .94894E+05,-.78193E+03,-.92088E+03, .39731E+02/) / 4.e8
	    cost_coefs_sym(1:4,2,1)=	&
		(/ .29512E+05,-.38329E+03,-.50861E+03, .44648E+02/) / 4.e8
	    cost_coefs_sym(1:4,3,1)=	&
		(/ .40384E+05,-.36197E+03,-.64097E+03, .56142E+02/) / 4.e8
	    cost_coefs_sym(1:4,4,1)=	&
		(/ .94894E+05,-.78193E+03,-.92088E+03, .39731E+02/) / 4.e8

	    cost_coefs_dia(1:3,1,1)=	&
		(/ .62958E+06,-.83096E+04, .31284E+02/) / 4.e8
	    cost_coefs_dia(1:3,3,1)=	&
		(/ .11080E+04,-.29164E+03, .34223E+02/) / 4.e8
	    cost_coefs_dia(1:3,4,1)=	&
		(/ .62958E+06,-.83096E+04, .31284E+02/) / 4.e8

	 Else
!!!!!!!!!!!
	    read(unit) cost_coefs_rec
	    read(unit) cost_coefs_sym
	    read(unit) cost_coefs_dia
!!!!!!!!!!!
	    call clsieee(unit, ierr)
	    if (ierr /= 0) call die(myname_,	&
		'clsieee("'//trim(infile)//'")',ierr)
     
	 End If
      End If

      ! Clear stats
      Do class = 1, MXcost_clas
	 nn = size(poly_set,2)
	 Allocate(stat_sym(class)%s(nn,nn,0:1), stat = ierr)
	    if (ierr /= 0) call die(myname_,'allocate(stat%s)',ierr)
            if(mall_ison()) call mall_mci(stat_sym(class)%s,myname)
	 stat_sym(class)%s = 0
	 Allocate(stat_rec(class)%s(nn,nn,0:1), stat = ierr)
	    if (ierr /= 0) call die(myname_,'allocate(stat%s)',ierr)
            if(mall_ison()) call mall_mci(stat_rec(class)%s,myname)
	 nn = size(poly_set_dia,2)
	 Allocate(stat_dia(class)%s(nn,nn,0:1), stat = ierr)
            if (ierr /= 0) call die(myname_,'allocate(stat_dia%s)',ierr)
	    if(mall_ison()) call mall_mci(stat_dia(class)%s,myname)

	 stat_sym(class)%s = 0
	 stat_rec(class)%s = 0
	 stat_dia(class)%s = 0
      End Do

     do i = 1, MXcost_clas
        stat_rec(i)%s = 0
        stat_sym(i)%s = 0
        stat_dia(i)%s = 0
     end do
      
      call MPI_BCAST(cost_coefs_rec, size(cost_coefs_rec),	&
	MP_type(cost_coefs_rec(1,1,1)), root, comm, ierr)
		if(ierr/=0) call MP_die(myname_,'MPI_bcast(rec)',ierr)
      call MPI_BCAST(cost_coefs_sym, size(cost_coefs_sym),	&
	MP_type(cost_coefs_sym(1,1,1)), root, comm, ierr)
		if(ierr/=0) call MP_die(myname_,'MPI_bcast(sym)',ierr)
      call MPI_BCAST(cost_coefs_dia, size(cost_coefs_dia), &
	MP_type(cost_coefs_dia(1,1,1)), root, comm, ierr)
		if(ierr/=0) call MP_die(myname_,'MPI_bcast(diag)',ierr)

      if (ierr /= 0) call MP_die(myname_,'MPI_bcase()',ierr)

!     Cost per element to form matrix blocks from look-ups and 
!     linear interpolation (currently the # of fp ops per element):
      costs(1,2) = 15.0
      costs(2,2) = 22
      costs(3,2) = 40.0
      costs(4,2) = 17.0

      if (DEBUG) then
         print *, myname,":: cost_class = "
         do i=1,ktmax
            write (*,*) cost_class(i,1:ktmax)
         end do
         print *, myname,":: costs(1:4,1) = ",costs(1:4,1)
         print *, myname,":: costs(1:4,2) = ",costs(1:4,2)
         print *, myname,":: costs(1:4,3) = ",costs(1:4,3)
      endif


      call system_clock(count_rate = clock_rate, count_max = max_clock_tick)
      clock_rate_inv = 1./clock_rate ! save away
!!$      clock_rate_inv = 1./4.e8 ! 400 mhz

!  Now determine whether this run should produce new tuning statistics
!  -------------------------------------------------------------------
    Call getenv('TUNE_COSTS',env_val)
    if ( len_trim(env_val) .gt. 0 ) then
         read(env_val,*,IOSTAT=ierr) tuning_is_on
             ALWAYS_ASSERT(ierr == 0)
    else
       tuning_is_on = .false.
    end if




    end subroutine init_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean up the object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_()
      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mco
      implicit none

! !REVISION HISTORY:
! 	25Jan01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: class
  integer :: ier

  if(.not.mb_cost_initialized) call die(myname_,'undefined object')

  do class=1,MXcost_clas
	if(mall_ison()) call mall_mco(stat_rec(class)%s,myname)
    deallocate(stat_rec(class)%s,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(stat)',ier)

	if(mall_ison()) call mall_mco(stat_sym(class)%s,myname)
    deallocate(stat_sym(class)%s,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(stat)',ier)

	if(mall_ison()) call mall_mco(stat_dia(class)%s,myname)
    deallocate(stat_dia(class)%s,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(diag)',ier)
  end do

  mb_cost_initialized=.false.

end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sym_mb_cost_ - compute the formation/access cost for a 
!            given block of a symmetric matrix.
!
! !DESCRIPTION:
!            This function returns the total cost associated with the
!            matrix block represented by the 4-tuple (kri,kti,krj,ktj)
!            for the symmetric matrix.
!
! !INTERFACE:

    function sym_mb_cost_( kind_mat, kind_cor, kri,kti,lni, krj,ktj,lnj)

      use config  , only : ktmax
      use m_die   , only : die

      implicit none

      integer , intent(IN) :: kind_mat ! matrix structure flag
      integer , intent(IN) :: kind_cor ! correlation type flag

      integer , intent(IN) :: kri      ! row kr index
      integer , intent(IN) :: kti      ! row kt index
      integer , intent(IN) :: lni      ! number of rows in matrix block

      integer , intent(IN) :: krj      ! column kr index
      integer , intent(IN) :: ktj      ! column kt index
      integer , intent(IN) :: lnj      ! number of columns in matrix block

      real :: sym_mb_cost_ 

! !REVISION HISTORY:
!       09Jul99 - J.W. Larson <jlarson@dao> - initial prototype/prolog/code
!EOP ___________________________________________________________________

      character(len=*),parameter :: myname_=myname//'::sym_mb_cost_'

      integer :: class                 ! cost class for the matrix block
      integer :: ctype                 ! cost type (TBD from kind_mat)
      real    :: num_elements          ! number of elements in matrix block
      integer :: idx

      include "kind_covs.h"
      include "kind_mats.h"

      if(.not.mb_cost_initialized) call die(myname_,	&
	'object not initialized')

!     Argument checking:

      if ( (kti <= 0) .or. (kti > ktmax) )	&
         call die(myname_,"bad argument kti",kti)

      if ( lni <= 0 )	&
         call die(myname_,"bad argument lni",lni)

      if ( (ktj <= 0) .or. (ktj > ktmax) )	&
         call die(myname_,"bad argument ktj",ktj)

      if ( lnj <= 0 )	&
         call die(myname_,"bad argument lnj",lnj)

!________________________________________

!     Determine which cost class is to be used:

      select case (kind_cor)
      case (kind_covF,kind_covS,kind_covV,kind_covU,kind_covC)
         class = cost_class(kti,ktj)
      case default
	 call die(myname_,"unsupported kind_cor",kind_cor)
      end select

!     Erorr class is undefined if (cost <= 0 or cost > MXcost_clas) 

      if ((class <= 0) .or. (class > MXcost_clas))	&
        call die(myname_,"Undefined Cost class",class)

!     Determine cost type ctype from kind_mat.
!     For non-univariate cases, BLAS Cxpy is fused with the table lookups. ctype = 1
!     For univariate cases,     BLAS Cxpy is separate.                     ctype = 2
!     nvecs = 1 is assumed

      select case (kind_mat)
      case (kind_5mat, kind_Rmat)
         ctype = 1
      case (kind_Umat)
	 ctype = 2 
      case default
         call die(myname_,"unsupported kind_mat",kind_mat)
      end select

!     Use class and ctype to index the array costs(:,:), and
!     compute the total cost for this matrix block:

!     Compute the number of elements in the matrix block:

      if ((kri == krj) .and. (kti == ktj)) then
	 sym_mb_cost_ = cost_estimate(lni, lnj, cost_coefs_dia(:,class,ctype), poly_set_dia)
      else
	 sym_mb_cost_ = cost_estimate(lni, lnj, cost_coefs_sym(:,class,ctype), poly_set)
      endif
      ! ensure that cost is non-negative
      sym_mb_cost_ = max(MIN_COST,sym_mb_cost_)
!!$      sym_mb_cost_ = cost_coefs_sym(1,class,idx) + lni*cost_coefs_sym(2,class,idx) + &
!!$	   lnj * cost_coefs_sym(3,class,idx) + lni*lnj*cost_coefs_sym(4,class,idx)

!!$         num_elements = lni * (lnj+1) / 2.
!!$      sym_mb_cost_ = num_elements * costs(class,ctype)

    end function sym_mb_cost_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rec_mb_cost_ - compute the formation/access cost for a given 
!            block of a rectangular matrix.
!
! !DESCRIPTION:
!
!            This function returns the total cost associated with the
!            matrix block represented by the 4-tuple (kri,kti,krj,ktj)
!            for the rectangular matrix.
!
! !INTERFACE:

    function rec_mb_cost_( kind_cor, kri, kti, lni, krj, ktj, lnj )

      use config  , only : ktmax
      use m_die   , only : die

      implicit none

      integer , intent(IN) :: kind_cor ! correlation type flag

      integer , intent(IN) :: kri      ! row kr index
      integer , intent(IN) :: kti      ! row kt index
      integer , intent(IN) :: lni      ! number of rows in matrix block

      integer , intent(IN) :: krj      ! column kr index
      integer , intent(IN) :: ktj      ! column kt index
      integer , intent(IN) :: lnj      ! number of columns in matrix block

      real :: rec_mb_cost_ 

! !REVISION HISTORY:
!       09Jul99 - J.W. Larson <jlarson@dao> - initial prototype/prolog/code
!EOP ___________________________________________________________________

      character(len=*),parameter :: myname_=myname//'::rec_mb_cost_'

      integer :: class                 ! cost class for the matrix block
      integer :: ctype                 ! cost type (TBD)
      real    :: num_elements          ! number of elements in matrix block
      include "kind_covs.h"

      if(.not.mb_cost_initialized) call die(myname_,	&
	'object not initialized')

!     Argument checking:

      if ( (kti <= 0) .or. (kti > ktmax) )	&
         call die(myname_,"bad argument kti",kti)

      if ( lni <= 0 )	&
         call die(myname_,"bad argument lni",lni)

      if ( (ktj <= 0) .or. (ktj > ktmax) )	&
         call die(myname_,"bad argument ktj",ktj)

      if ( lnj <= 0 )	&
         call die(myname_,"bad argument lnj",lnj)

!________________________________________

!     Compute the number of elements in the matrix block:

      num_elements = lni * lnj

!     Determine which cost class is to be used:

      select case (kind_cor)
      case (kind_covF,kind_covS,kind_covV)
         class = cost_class(kti,ktj)
      case default
	 call die(myname_,"unsupported kind_cor = ",kind_cor)
      end select

!     Erorr class is undefined if (cost <= 0 or cost > MXcost_clas) 
      if ((class <= 0) .or. (class > MXcost_clas))	&
        call die(myname_,"undefined cost class = ",class)

!     Determine the cost type from kind_mat.  Currently, ctype
!     is hard-wired with value 1 (table lookups).

      ctype = 1

!     Use class and ctype to index the array costs(:,:), and
!     compute the total cost for this matrix block:

!!$      rec_mb_cost_ = num_elements * costs(class,ctype)

!!$      rec_mb_cost_ = cost_coefs_rec(1,class,ctype) + lni*cost_coefs_rec(2,class,ctype) + &
!!$                     lnj * cost_coefs_rec(3,class,ctype) + num_elements*cost_coefs_rec(4,class,ctype)
      rec_mb_cost_ = cost_estimate(lni, lnj, cost_coefs_rec(:,class,ctype), poly_set)
      ! ensure that cost is non-negative
      rec_mb_cost_ = max(MIN_COST,rec_mb_cost_)

    end function rec_mb_cost_

    Function cost_estimate(ni, nj, weights, poly_set) result(est)
      integer, intent(in) :: ni, nj
      real, intent(in) :: weights(:)
      integer, intent(in) :: poly_set(:,:)
      real :: est

      integer :: i
      real :: xi, xj, tmp

      est = 0
      xi = real(ni)
      xj = real(nj)
      tmp = 0
      do i = 1, size(poly_set,2)
	 tmp = tmp + weights(i) * xi**poly_set(1,i) * xj**poly_set(2,i)
      end do
      est = tmp
    end Function cost_estimate
	    
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sym_compute_sparsity_ - compute the sparsity of a square
!            symmetric matrix from its navigator mNav and its dimension 
!            nrow.
!
! !DESCRIPTION:
!
!            This routine processes the symmetric matrix navigator mNav 
!            to compute the local number of matrix elements lnum\_elements,
!            taking into account the fact that only the diagonal and 
!            half of the off-diagonal blocks are represented in mNav.
!            A global sum of the local element counts yields the global
!            element count gnum\_elements.  This is compared with the 
!            product of the number of rows (nrow) and columns (nrow) in
!            the matrix to yield the sparsity.  Both gnum\_elements and
!            sparsity are returned on all PE's.
!
!
! !INTERFACE:

    subroutine sym_compute_sparsity_( mNav, nrow, lnum_elements,     &
                                    gnum_elements, sparsity, root, & 
                                    comm, ierr                      )

      use m_mpif90

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_indexIA => indexIA
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_AttrVect, only : ptr_iAttr
      use m_AttrVect, only : ptr_rAttr

      use m_die   , only : die

      implicit none

      type(AttrVect) ,  intent(IN)  :: mNav     ! Matrix Navigator
      integer        ,  intent(IN)  :: nrow     ! global number of rows

      real           ,  intent(OUT) :: lnum_elements     ! local number
                                                         ! of matrix elements
      real           ,  intent(OUT) :: gnum_elements     ! global number
                                                         ! of matrix elements
      real           ,  intent(OUT) :: sparsity          ! matrix sparsity

      integer        ,  intent(IN)  :: root              ! root PE
      integer        ,  intent(IN)  :: comm              ! communicator

      integer        ,  intent(OUT) :: ierr              ! error flag

! !REVISION HISTORY:
!       19Jul99 - J.W. Larson <jlarson@dao> - initial prototype/prolog/code
!        6Oct99 - J.W. Larson <jlarson@dao.gsfc.nasa.gov> - replaced 
!                 occurances of MP_REAL with MP_type(message) for 
!                 better support of 32 and 64-bit platforms.
!EOP ___________________________________________________________________

      character (len=*),parameter :: myname_=myname//'::sym_compute_sparsity_'
      integer,pointer,dimension(:,:) :: iAttr
      real   ,pointer,dimension(:,:) :: rAttr

!  Local work variables
      integer irow, icol, inels  ! row, column, and number of elements
      integer nBlox              ! number of matrix blocks in mNav
      integer i                  ! a humble loop counter

!  Argument checks:

      if ( nrow <= 0 )	&
         call die(myname_,"bad argument nrow",nrow)

!  Determine local number of matrix blocks in mNav:

   nBlox = AttrVect_lsize(mNav)

!  Determine row and column indices in mNav:
   irow = AttrVect_indexRA(mNav,item='irow')
   icol = AttrVect_indexRA(mNav,item='icol')

!  Determine index in mNav for item nels (number of elements)

   inels = AttrVect_indexRA(mNav,item='nels')

!  Initialize local and global numbers of matrix elements:

   lnum_elements = 0.
   gnum_elements = 0.

		! Alias mNav%iAttr and mNav%rAttr

	iAttr => ptr_iAttr(mNav)
	rAttr => ptr_rAttr(mNav)

!  Compute local number of matrix elements:

   do i=1,nBlox
      if (iAttr(irow,i) == iAttr(icol,i)) then
         lnum_elements = lnum_elements + rAttr(inels,i)
      else
         lnum_elements = lnum_elements + 2. * rAttr(inels,i)
      endif
   end do

		! Unalias mNav%iAttr and mNav%rAttr

	nullify(iAttr)
	nullify(rAttr)
   
!  Compute global number of matrix elements:

   call MPI_REDUCE( lnum_elements, gnum_elements, 1,   &
                    MP_type(lnum_elements), MP_SUM,    &
                    root, MP_COMM_WORLD, ierr           )

!  Give all PE's this information:

   call MPI_BCAST( gnum_elements, 1, MP_type(gnum_elements), &
                   root, comm, ierr )

!  Compute the sparsity of the square symmetric matrix

   sparsity = gnum_elements / (float(nrow) * float(nrow))

   end subroutine sym_compute_sparsity_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rec_compute_sparsity_ - compute the sparsity of a rectangular
!            matrix.
!
! !DESCRIPTION:
!
!            This routine processes the matrix navigator mNav to 
!            compute the local number of matrix elements lnum\_elements.
!            A global sum of the local element counts yields the global
!            element count gnum\_elements.  This is compared with the 
!            product of the number of rows (nrow) and columns (ncol) in
!            the matrix to yield the sparsity.  Both gnum\_elements and
!            sparsity are returned on all PE's.
!
! !INTERFACE:

    subroutine rec_compute_sparsity_( mNav, nrow, ncol, lnum_elements, &
                                      gnum_elements, sparsity, root,   &
                                      comm, ierr)

      use m_mpif90

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_AttrVect, only : ptr_rAttr

      use m_die   , only : die

      implicit none

      type(AttrVect) ,  intent(IN)  :: mNav     ! Matrix Navigator
      integer        ,  intent(IN)  :: nrow     ! global number of rows
      integer        ,  intent(IN)  :: ncol     ! global number of columns

      real           ,  intent(OUT) :: lnum_elements     ! local number
                                                         ! of matrix elements
      real           ,  intent(OUT) :: gnum_elements     ! global number
                                                         ! of matrix elements
      real           ,  intent(OUT) :: sparsity          ! matrix sparsity

      integer        ,  intent(IN)  :: root              ! root PE
      integer        ,  intent(IN)  :: comm              ! communicator

      integer        ,  intent(OUT) :: ierr              ! error flag

! !REVISION HISTORY:
!       19Jul99 - J.W. Larson <jlarson@dao> - initial prototype/prolog/code
!        6Oct99 - J.W. Larson <jlarson@dao.gsfc.nasa.gov> - replaced 
!                 occurances of MP_REAL with MP_type(message) for 
!                 better support of 32 and 64-bit platforms.
!EOP ___________________________________________________________________

      character(len=*),parameter :: myname_=myname//'::rec_compute_sparsity_'

      real   ,pointer,dimension(:,:) :: rAttr

!  Local work variables
      integer inels              ! number of elements
      integer nBlox              ! number of matrix blocks in mNav
      integer i                  ! a humble loop counter

!  Argument checks:

      if ( nrow <= 0 )	&
         call die(myname_,"bad argument nrow",nrow)

      if ( ncol <= 0 )	&
         call die(myname_,"bad argument ncol",ncol)

!  Determine local number of matrix blocks in mNav:

   nBlox = AttrVect_lsize(mNav)

!  Determine index in mNav for item nels (number of elements)

   inels = AttrVect_indexRA(mNav,item='nels')

!  Initialize local and global numbers of matrix elements:

   lnum_elements = 0.
   gnum_elements = 0.

!  Compute local number of matrix elements:

	rAttr => ptr_rAttr(mNav)	! Alias mNav%rAttr

   do i=1,nBlox
      lnum_elements = lnum_elements + rAttr(inels,i)
   end do

	nullify(rAttr)			! Unalias mNav%rAttr

!  Compute global number of matrix elements:

   call MPI_REDUCE( lnum_elements, gnum_elements, 1,  &
                    MP_type(lnum_elements), MP_SUM,   &
                    root, MP_COMM_WORLD, ierr           )

!  Give all PE's this information:

   call MPI_BCAST( gnum_elements, 1, MP_type(gnum_elements), &
                   root, comm, ierr )

!  Compute the sparsity of the rectangular matrix

   sparsity = gnum_elements / (float(nrow) * float(ncol))

   end subroutine rec_compute_sparsity_


   subroutine tune_cost_start_()
     use m_die,only : assert_
     Implicit None

     Character(len=*), Parameter :: myname_ = myname//'::tune_cost_start_'
     integer :: i
     integer :: ier
     integer :: ithread

#ifdef _OPENMP
	ithread=OMP_GET_THREAD_NUM()
	ASSERT(ithread <= MAX_THREADS)
     call system_clock(count=start_count(1,ithread))
#else
     call system_clock(count=start_count)
#endif
!!$     call irtc(start_count)

   end subroutine tune_cost_start_

   subroutine tune_cost_stop_(kti, lni, ktj, lnj, mtype, ctype)
     use m_die, only : die, assert_
     use m_mpout, only : mpout_log, mpout
     Implicit None
!!$     integer, intent(in) :: kri, krj  ! segment types
     integer, intent(in) :: kti, ktj  ! segment types
     integer, intent(in) :: lni, lnj  ! segment sizes
     character(len=3), intent(in) :: mtype
     integer, intent(in) :: ctype

      include "kind_covs.h"
      include "kind_mats.h"

     Character(len=*), Parameter :: myname_ = myname//'::tune_cost_stop_'
     real :: block_time
     integer :: class
!RT  integer, external :: my_pe

     integer :: ier
     integer :: i, j, k
     integer :: stop_count
     integer :: idx
     real :: estimate, delta, temp
     real :: xlni, xlnj
     real, pointer :: weights(:), st(:,:,:)
     integer, pointer :: nd(:,:)
     Integer :: c0
     Integer :: scount
     integer :: ithread

     call system_clock(c0)
!!$     call irtc(c0)

#ifdef _OPENMP
	ithread=OMP_GET_THREAD_NUM()
	ASSERT(ithread <= MAX_THREADS)
	scount = start_count(1,ithread)
#else
	scount = start_count
#endif
     If (c0 < scount) Then
        block_time = (c0-scount + max_clock_tick) * clock_rate_inv
     Else
	block_time = (c0-scount) * clock_rate_inv
     End If

     ! assume a simple model
     ! t = a + b*lni + c*lnj + d*lni*lnj

     class = cost_class(kti,ktj)


     Select Case (mtype)
     Case ('dia')
	weights => cost_coefs_dia(:,class,ctype)
	nd => poly_set_dia
	st => stat_dia(class)%s
     Case ('sym')
	weights => cost_coefs_sym(:,class,ctype)
	nd => poly_set
	st => stat_sym(class)%s
     Case ('rec')
	weights => cost_coefs_rec(:,class,ctype)
	nd => poly_set
	st => stat_rec(class)%s
     Case Default
	Call die(myname_,"bad argument mtyp: "//mtype)
     End Select

     estimate = cost_estimate(lni, lnj, weights, nd)
     delta = block_time - estimate

     ! write stats for further analysis
     xlni = real(lni)
     xlnj = real(lnj)
     
     Do k = 0, 1
	do j = 1, size(nd,2)
	   do i = 1, size(nd,2)
!$OMP ATOMIC
	      st(i,j,k) =  st(i,j,k) + xlni**(nd(1,i)+nd(1,j))*xlnj**(nd(2,i)+nd(2,j)) * delta**k
	   end do
	end do
     end do

   end subroutine tune_cost_stop_

   subroutine print_cost_statistics_()
     use m_mpif90
     use m_die, only : assert_,die
     use m_mall,only : mall_ison,mall_mci,mall_mco
     use m_ioutil,  only : luavail, opnieee, clsieee, opntext, clstext
     use m_mpout, only : mpout_log, mpout
     Implicit None
     Character(len=*), Parameter :: myname_ = myname//'::print_cost_statistics_'

     integer :: m
     integer :: n
     real,    Allocatable :: mat(:,:), rhs(:)
     integer, Allocatable :: ipiv(:)
     integer :: info
     integer :: class

!RT  integer, external :: my_pe
     integer :: unit, i, j, k, idx, max_idx, ier
     integer :: root, myid, nn_sym, nn_rec, nn_dia, nn
     character(128) :: infile
     type (statistics) :: stat_tot_sym(MXcost_clas)
     type (statistics) :: stat_tot_rec(MXcost_clas)
     type (statistics) :: stat_tot_dia(MXcost_clas)
     real, pointer :: stats(:,:,:)
     real :: new_cost_coefs_rec(size(poly_set,2), MXcost_clas, MXcost_type)
     real :: new_cost_coefs_sym(size(poly_set,2), MXcost_clas, MXcost_type)
     real :: new_cost_coefs_dia(size(poly_set_dia,2), MXcost_clas, MXcost_type)
     integer, pointer :: nd(:,:)


     root = 0
     call MP_COMM_RANK(MP_COMM_WORLD,myid,ier)
     ASSERT(ier == 0)

     Call mpout_log(myname_,'accumulating statistics', myid)

     Do class = 1, MXcost_clas
	If (myid == root) Then

	   nn_dia= size(poly_set_dia,2)
	   Allocate(stat_tot_dia(class)%s(nn_dia,nn_dia,0:1), stat = ier)
	   ALWAYS_ASSERT(ier == 0)
              if(mall_ison())	call mall_mci(stat_tot_dia(class)%s,myname)

	   nn_sym = size(poly_set,2)
	   Allocate(stat_tot_sym(class)%s(nn_sym,nn_sym,0:1), stat = ier)
	   ALWAYS_ASSERT(ier == 0)
              if(mall_ison())	call mall_mci(stat_tot_sym(class)%s,myname)

	   nn_rec = size(poly_set,2)
	   Allocate(stat_tot_rec(class)%s(nn_rec,nn_rec,0:1), stat = ier)
	   ALWAYS_ASSERT(ier == 0)
              if(mall_ison())	call mall_mci(stat_tot_rec(class)%s,myname)
	else

	   Allocate(stat_tot_dia(class)%s(1,1,0:0), stat_tot_sym(class)%s(1,1,0:0), stat_tot_rec(class)%s(1,1,0:0), stat = ier) ! stub
	   ALWAYS_ASSERT(ier == 0)
	      if(mall_ison())	Then
		call mall_mci(stat_tot_dia(class)%s,myname)
		call mall_mci(stat_tot_sym(class)%s,myname)
		call mall_mci(stat_tot_rec(class)%s,myname)
	     end if
	End If

	Call MPI_REDUCE( stat_dia(class)%s, stat_tot_dia(class)%s, size(stat_dia(class)%s),   &
	     MP_type(stat_dia(class)%s(1,1,0)), MP_SUM,    &
	     root, MP_COMM_WORLD, ier           )
	ALWAYS_ASSERT(ier == 0)

	Call MPI_REDUCE( stat_sym(class)%s, stat_tot_sym(class)%s, size(stat_sym(class)%s),   &
	     MP_type(stat_sym(class)%s(1,1,0)), MP_SUM,    &
	     root, MP_COMM_WORLD, ier           )
	ALWAYS_ASSERT(ier == 0)

	Call MPI_REDUCE( stat_rec(class)%s, stat_tot_rec(class)%s, size(stat_rec(class)%s),   &
	     MP_type(stat_rec(class)%s(1,1,0)), MP_SUM,    &
	     root, MP_COMM_WORLD, ier           )
	ALWAYS_ASSERT(ier == 0)

     End Do
     
     If (myid == root) Then

 
     new_cost_coefs_rec = cost_coefs_rec
     new_cost_coefs_sym = cost_coefs_sym
     new_cost_coefs_dia = cost_coefs_dia

     Do idx = 1, 3
	Do class = 1, MXcost_clas ! loop over classes

	   Select Case (idx)
	   Case (1)
	      nd => poly_set_dia
	      stats => stat_tot_dia(class)%s
	   Case (2)
	      nd => poly_set
	      stats => stat_tot_sym(class)%s
	   Case (3)
	      nd => poly_set
	      stats => stat_tot_rec(class)%s
	   End Select
	   
	   if (stats(1,1,0) == 0) Then
	      cycle ! no stat_tots for this class
	   end if
	   
	   nn = size(nd,2)
	   Allocate(mat(nn,nn), rhs(nn), ipiv(nn), stat = ier)

	   ALWAYS_ASSERT(ier == 0)
	   
	   if(mall_ison()) then
	     call mall_mci(mat ,myname)
	     call mall_mci(rhs ,myname)
	     call mall_mci(ipiv,myname)
	   endif

	   Do j = 1, nn
	      Do i = 1, nn
		 mat(i,j)  = stats(i,j,0)
	      End Do
	      rhs(j)    = stats(1,j,1)
	   end Do
	   
	   ! solve matrix ...
	   ! 1) LU decomposition
	   ! 2) Solve (overwrite rhs)
	   
	   
	   call _SGETRF(nn, nn, mat, nn, ipiv, info)
	   ALWAYS_ASSERT(info == 0)

	   call _SGETRS('N', nn, 1, mat, nn, ipiv, rhs, nn, info)
	   ALWAYS_ASSERT(info == 0)
	   
	   write(mpout,'(2i4,a9,4(4x,e15.5))')class, idx, ' delta = ', rhs
	   Select Case(idx)
	   Case (1)
	      new_cost_coefs_dia(:,class,1) = rhs + cost_coefs_dia(:,class,1)
	      write(mpout,'(8x,a9,4(4x,e15.5))')             '   new = ', cost_coefs_dia(:,class,1)
	   Case (2)
	      new_cost_coefs_sym(:,class,1) = rhs + cost_coefs_sym(:,class,1)
	      write(mpout,'(8x,a9,4(4x,e15.5))')             '   new = ', cost_coefs_sym(:,class,1)
	   Case (3)
	      new_cost_coefs_rec(:,class,1) = rhs + cost_coefs_rec(:,class,1)
	      write(mpout,'(8x,a9,4(4x,e15.5))')             '   new = ', cost_coefs_rec(:,class,1)
	   End Select
	   
	   ASSERT(info == 0)

	   if(mall_ison()) then
	     call mall_mco(mat ,myname)
	     call mall_mco(rhs ,myname)
	     call mall_mco(ipiv,myname)
	   endif

	   Deallocate(mat, rhs, ipiv, stat=ier)
	   ALWAYS_ASSERT(ier == 0)
	End Do
     End Do
     nullify(nd)
     nullify(stats)

     unit = luavail()
     infile = 'costs.out'
     call opnieee(unit,infile,"unknown",ier)
     ASSERT(ier == 0)
	   
!!!!!!!!!!!
     write(unit) new_cost_coefs_rec
     write(unit) new_cost_coefs_sym
     write(unit) new_cost_coefs_dia

!!!!!!!!!!!
     call clsieee(unit, ier)
     ASSERT(ier == 0)


!!$     ! Write stat_tot for more detailed analysis
!!$     infile = 'cost_stat.out'
!!$     unit = luavail()
!!$     call opntext(unit,infile,"unknown",ier)
!!$     do idx = 1, max_idx
!!$	if (idx == 1) Then ! non-diagonal case
!!$	   m = 4
!!$	   n = 4
!!$	else
!!$	   m = 3
!!$	   n = 3
!!$	end if
!!$
!!$	do class = 1, 4 ! loop over classes
!!$	
!!$	   if (stat_tot(class,idx)%s(1,1,0) == 0) cycle ! no stat_tots for this class
!!$	   Do k = 0, 1
!!$	      Do j = 0, 8
!!$		 Do i = 0,8
!!$		    write(unit,'(5i3,e20.15)') class, idx, k, j, i, stat_tot(class,idx)%s(i,j,k)
!!$		 End Do
!!$	      End Do
!!$	   End Do
!!$	end do
!!$     end do
!!$     call clstext(unit, ier)

  End If
  ! clean up
     Do class = 1, MXcost_clas
	   if(mall_ison()) Then
	     call mall_mco(stat_tot_dia(class)%s,myname)
	     call mall_mco(stat_tot_sym(class)%s,myname)
	     call mall_mco(stat_tot_rec(class)%s,myname)
	  end if
	Deallocate(stat_tot_dia(class)%s, stat_tot_sym(class)%s, stat_tot_rec(class)%s, stat=ier)
	ALWAYS_ASSERT(ier == 0)
     End Do

   end subroutine print_cost_statistics_

   Logical Function tuneQ()
     tuneQ = tuning_is_on
   End Function tuneQ

 end module m_costs

