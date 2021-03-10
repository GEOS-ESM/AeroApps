!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_amat - Masking A matrix for reduced PSAS
!
! !DESCRIPTION:
! module to generate and store masking A matrix for redpsas project
! Peter Lyster lys@dao.gsfc.nasa.gov.
!
! The A matrix a an nreg x nreg matrix (were nreg is the number of regions 
! in the the PSAS that are used to aggregate sub-blocks of error covariance 
! matrices) that masks the error covariance matrices of the PSAS.  Specifically, 
! it masks M=HP^fH^T+R and P^fH^T. The A matrix is designed to have a greater 
! degree of sparsity than the error covariance matrices and thus may be used in
! this way to reduce the computational complexity of the calculation.  The 
! masking operation is designed to leave the matrices M and P^fH^T approximately 
! unchanged (at least on the 'diagonals'), while because the masking  operation 
! is of the Hadamard product type it preserves the positive semi definiteness 
! property of M -- this is important to guarantee a stable solution to 
! the PSAS equations.  Check P.Lyster's personal notes around 20000216.

! The A matrix is formed by the outer product A = BB^T (i.e., is positive
! semi definite) where B may be nreg x nreg (B^a) or nreg x nnod (B^b).
! In the software below assume nreg=maxreg==80 for the current singly-refined 
! triangle regions based on the icosahedral decomposition on the sphere. 
! "Connectivity" is an important concept here because locality is important 
! in assessing the impact of the changes that the A matrix cause -- in effect, 
! you want minimal changes to the proper influence of  nearby observations 
! on analysis points but you are willing to accept a reduction 
! of the influence of remote observations on analysis points if this is physically 
! justifiable.  The B matrix is the basic connectivity that gives rise (through BB^T) 
! the final connectivity in A.  All connectivity for the refined icosahedral regions 
! can be defined from a node-region connectivity graph, called 'bindex(6,nreg)' below.  
! For the current software there are nreg=80 regions in the singly-refined icosahedral 
! regions, and 42 associated nodes (i.e., nreg=2*node-4). For the icosahedral 
! mesh 12 of the nodes are special, they are like 'poles', and are abutted by 
! 5 regions. All other (nnodes - 12) nodes  -- there are 30 other nodes for 
! the singly-refined icosahedral regions -- are abutted by  6 regions.  The 
! bindex matrix then lists the regions, numbered according to the Pfaendtner 
! convention (DAO Office Note 96-04 at http://dao.gsfc.nasa.gov/) that abut 
! to nodes ordered implicitly 1 to 42.  This effectively defines the node 
! numbering that is used here (a new node numbering will have to be devised 
! for more than single level refinement of the icosahedral regions).

! Before explaining how the B matrices are formed, I will provide here a 
! nomenclature for describing the connectivity of the A matrices.  The diagonal 
! elements of A apply to regions that are connected to themselves in the 
! sense that all observations and analysis points in these regions are geographically 
! local.  Setting these diagonal elements to 1 is a minimal
! requirement for guaranteeing that at least the influence of local observations
! on nearby analysis points is preserved.  Each region is connected to three
! other regions by sharing a common side -- this is defined here as type (i) 
! connectivity.  Each region is connected to slightly more distant regions 
! because they share common node and a common side with the type (i) regions 
! (in a sense, the type (i) connected regions lie between them connected to 
! a side of each).  There are however two kinds of such regions -- type (ii) regions 
! share this connectivity but the shared node is one of the 12 unique nodes described 
! above, while type (iii) connectivity is where they share this connectivity but 
! their shared node is one of the (nnodes - 12).  Finally type (iv) regions 
! are connected to a particular region in that both share a common one of 
! the (nnodes - 12) nodes but don't share a side with one of the type (i) connected 
! regions, i.e., they are further removed, and in a sense, spatially mirror-reflect each 
! other 180 degrees across their common node (I suggest to put in a postscript graphic 
! here to illustrate this similar to the graphic on page 1 of my notes of 20000216).

! As described above there are two B^a and B^b corresponding the region-region and 
! region-node connectivity respectively,i.e., either type of B matrix generates a 
! legitimate but different A matrix -- its a matter of exploring the options.  
! A further superscript k (B^{a,k} and B^{b,k}) defines how far out the connectivity of 
! B is accreted in terms of layers that share vertices.  I haven't defined this 
! concept rigorously yet -- and I think this gets at the heart of the difficulty of 
! defining a truly spatially homogeneous A, which is important for ensuring
! that the application of the A matrix retains the integrity of the physical solution 
! -- but for now I will define the concepts experientially in terms of the 
! implementations below.  Further subscript j (B^{a,k}_j and B^{b,k}_j) just 
! differentiates between different B's with connectivity at the same k accretion level.  
! Note that for every B matrix there are often different possible A's depending on 
! how you value the indices within B. For the present we always set all values of 
! B to 1 (20000310).

! Hence,

! B^{a,0}_0 is the identity connectivity (the minimum possible type B^a connectivity) 
!    which gives
!  A^{a,0}_0 == B^{a,0}_0 (B^{a,0}_0)^T = the identity (note that the letter 
!    superscript of the A matrix is just set to the superscript of the B matrix 
!    from which it derives, i.e., type 'a'  in this case).

! B^{b,0}_0 is where regions are connected to their own defining nodes (the minimum 
!    possible type B^b connectivity) which gives
!  A^{b,0}_0 == B^{b,0}_0 (B^{b,0}_0)^T which has (experientially) types
!    (i), (ii), (iii), and (iv) connectivity, and AM (i.e., Ao(HP^fH^T+R) where 'o' 
!    means the Hadamard product of elements at the region-region granularity) has 
!    16.4% sparsity (actually this should be called "1 minus sparsity" ... so sue me).  
!    When the connectivity indices, i.e., the values of,  B^{b,0}_0 are all set to 1, 
!    and A^{b,0}_0 is normalized to diagonal values of 1, then the off-diagonal elements 
!    of A^{b,0}_0 with type (i) connectivity  have value 2/3 and the other off-diagonal 
!    elements of (type (ii),  (iii), and (iv) connectivity) have values 1/3.  
!    This matrix is sometimes called Amat_0 and is generated in the module subroutine 
!    init_0_() below.  This matrix is also used to generate other B matrices.  
!    Therefore (see next paragraphs), B^{a,1}_0 is formed from A^{b,0}_0 by putting 1 
!    in the positions which have value 1 and 2/3.  Similarly, B^{a,1}_3 formed from 
!    A^{b,0}_0 by putting 1 in the positions which have value 1, 2/3, and 1/3 in A^{b,0}_0. 

! B^{a,1}_0 is where each region is connected to only corresponding type (i) 
!    regions (with index value 1 by default)  which gives
!  A^{a,1}_0 == B^{a,1}_0 (B^{a,1}_0)^T which has (experientially) type (i), (ii), 
!    and (iii) connectivity, and AM has  13.6% sparsity.  When the connectivity indices of, 
!    i.e., the values of,  B^{a,1}_0 are all  set to 1, and A^{a,1}_0 is normalized 
!    to diagonal values of 1, then the off-diagonal elements of A^{a,1}_0 with 
!    type (i) connectivity have value 1/2 and the other off-diagonal elements (type (ii) and 
!    (iii) connectivity) have values 1/4.  This matrix is sometimes called Amat_1 and 
!    is generated in the module subroutine init_1_() below. 

! B^{a,1}_1 is where each region is connected to corresponding types (i) and (ii) 
! regions, and I haven't implemented that yet (20000310).

! B^{a,1}_2 is where each region is connected to corresponding types (i), (ii), 
! and (iii) regions, and I haven't implemented that yet (20000310).

! B^{a,1}_3 is where each region is connected to corresponding types (i), (ii), 
!    (iii), and (iv) regions  (with index value 1 by default) which gives
!  A^{a,1}_3 == B^{a,1}_3 (B^{a,1}_3)^T which has (experientially) the 
!    same connectivity as the PSAS's M matrix with level_for_banded_approximation set 
!    to 5 (i.e., neglect corr R > 6000 km).  With this, AM has 40.4% sparsity (i.e., the same 
!    sparsity as M, but the off-diagonal values of AM are less than those of M because the
!    off diagonal values of A are less than 1).  This matrix 
!    is sometimes called Amat_2 and is generated in the module subroutine 
!    init_2_() below.

!
! !INTERFACE:

 module m_amat

  implicit none
  private   ! except

! The A matrix structure

   public :: amat

! The sparsity specification for "level_of_band_approximation:" or
! kind_mat.  Note that they are just temporarily defined here as the
! first step to separate these parameters from bands0 and CGSolver.

   public :: maxband
   public :: seplimband

! The functions that generate/clean the A matrix

   public :: amat_init
   public :: amat_clean

   interface amat_init ; module procedure init_ ; end interface
   interface amat_clean; module procedure clean_; end interface

! !REVISION HISTORY:
!       27Sep99 - Peter Lyster <lys@dao.gsfc.nasa.gov> and Banglin Zhang: 
!                 initial prototype/prolog/code.
!        8Mar00 - Peter Lyster: revamp to initialize amat online.
!EOP ___________________________________________________________________

   character(len=*),parameter :: myname='m_amat'

   real,save,allocatable :: amat(:,:)
   logical,save :: initialized_=.false.

! This parameter is expected to be phased out.

	! MAXBAND replaces nbandcg.  This parameter was used to specify
	! the sparsity of the full matrix.  However, the sparsity of a
	! correlation matrix should be explicitly controlled by the
	! support of its models.  Therefore, this number is now used
	! only as a fixed flag for a "full" matrix.

   ! integer,parameter :: MAXBAND=5

	! Put back the user defiiable maxband (i.e. nband)

   integer,save :: maxband
   real   ,save :: seplimband

   real,dimension(4:5),parameter :: seplim=(/26.50,58.25/)

		! Defining the function labels.  See the code below for
		! the details of function implementations.

   integer,parameter :: TOP_HAT	   =0
   integer,parameter :: WITCHES_HAT=1
   integer,parameter :: COSINE_HILL=2
   integer,parameter :: DEFAULT_HAT=TOP_HAT

   logical,parameter :: CHECK_AMAT=.true.

#include "assert.H"

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - public init to initialize any of a number of amats
!
! !DESCRIPTION:
!     subroutine to call the initialization routines of different masking 
!     matrices A which are used to impose sparsity on the PSAS multivariate
!     error covariance matrices.   The decision is made by a variable,
!     amat_type.
!
! !INTERFACE:

 subroutine init_(gpart,comm,root,rsrc)

   use m_GlobalPartition,only : GlobalPartition
   use m_GlobalPartition,only : ptr_Partition

   use m_Spherical_Partition,only : NumberOfRegions

   use m_psasrc ,only : psasrc_open,psasrc_close
   use m_inpak90,only : i90_label,i90_gfloat,i90_gint

   use m_mpout  ,only : mpout_log
   use m_die    ,only : die,perr,MP_die,assert_
   use m_mall   ,only : mall_ison,mall_mci
   use m_mpif90 ,only : MP_comm_rank,MP_type

   implicit none

   type(GlobalPartition),intent(in) :: gpart
   integer,intent(in) :: comm
   integer,intent(in) :: root
   character(len=*),optional,intent(in) :: rsrc

   include "kind_mats.h"	! use only kind_5mat

! !REVISION HISTORY:
!	24Jan02	- Jing Guo
!		. Modified the interface.
!		. Moved psasrc_open()/psasrc_close() here from m_AE.
!       07Oct00 - Peter Lyster <lys@dao.gsfc.nasa.gov> Added init from
!		  PSAS resource
!EOP ___________________________________________________________________

   character(len=*),parameter :: myname_=myname//'::init_'

! Generally we should use amat_type == AMAT_BxBT which invokes an
! algorithm which directly calculates a B matrix based on a primitive
! (generating) function on a grid of refined icosahedral spherical
! triangles.  Hence this can use 20,80,320,1280,5120,... regions.
! Keep procedures for amat_type == AMAT_UNITY (default) for heritage
! purposes.

   character(len=*),parameter :: rsrcACUTOFF ="amat_cutdistance:"
   character(len=*),parameter :: rsrcFUNCTION="amat_function:"

! This is a phased out parameter.  It is defined here to issue a
! warning only.

   character(len=*),parameter ::	&
	rsrcLevelOfBanded="level_for_banded_approximation:"

   integer,parameter :: AMAT_UNITY=-1
   integer,parameter :: AMAT_BxBT =11	! Why this number?

   integer :: myID
   integer :: icount,jcount
   integer :: ier		! will this be re-initialized?
   integer :: amat_type
   real    :: amat_cutdistance 
   integer :: amat_function
   integer :: i1,i2
   integer :: nReg
   integer :: ibuf(3)

  if(initialized_) call die(myname_,'multiple definitions')
!_______________________________________________________________________

	call MP_comm_rank(comm,myID,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

if(myID==root) then

	! Loading the resource handle to initialize configuration
	! variables.

  call psasrc_open(rsrc=rsrc,stat=ier)
	if(ier/=0) call die(myname_,'psasrc_open()',ier)
!________________________________________

! Each init is responsible for its own error detection and response

	! Specify the cutoff model and the cutoff value.

  amat_type=AMAT_BxBT
  amat_cutdistance = -1.
  amat_function=DEFAULT_HAT	! default function if not specified.

	! Define amat_cutdistance, if a user has defined it.

  call i90_label(rsrcACUTOFF,ier)
  if(ier==0) then
    amat_cutdistance=i90_gfloat(ier)
    if(ier/=0) call die(myname_,'i90_gfloat("'//rsrcACUTOFF//'")',ier)
  endif

	! Redefine amat_type from amat_cutdistance, since amat_type can
	! not be explicitly define by a user.

  if(amat_cutdistance < 0.) amat_type=AMAT_UNITY

	! Specify the form of the generating function

  select case(amat_type)
  case(AMAT_BxBT)

    call i90_label(rsrcFUNCTION,ier)
    if(ier==0) then
      amat_function=i90_gint(ier)
      if(ier/=0) call die(myname_,'i90_gint("'//rsrcFUNCTION//'")',ier)
    endif
  end select

	! This parameter is yet to be phased out.  Now it is not clear
	! how this should be done.

  maxband=kind_5mat
  call i90_label(rsrcLevelOfBanded,ier)
    if(ier==0) then
      maxband=i90_gint(ier)
      if(ier/=0) call die(myname_,	&
	'i90_gint("'//rsrcLevelOfBanded//'")',ier)
    endif

	ASSERT(maxband>=lbound(seplim,1))
	ASSERT(maxband<=ubound(seplim,1))

	! Release PSASRC

  call psasrc_close(stat=ier)
	if(ier /= 0) call die(myname_,'psasrc_close()',ier)
endif

	ibuf(1)=amat_type
	ibuf(2)=amat_function
	ibuf(3)=maxband
call MPI_bcast(ibuf,size(ibuf),MP_type(ibuf),root,comm,ier)
    if(ier/=0) call MP_die(myname_,'MPI_bcast(amat_parameters)',ier)
	amat_type=ibuf(1)
	amat_function=ibuf(2)
	maxband=ibuf(3)

seplimband=seplim(maxband)
	
call MPI_bcast(amat_cutdistance,1,MP_type(amat_cutdistance),	&
	root,comm,ier)
    if(ier/=0) call MP_die(myname_,'MPI_bcast(amat_cutdistance)',ier)

initialized_=.true.
!_______________________________________________________________________
!
! Allocate the module data storage, amat(:)

  nReg=NumberOfRegions(ptr_Partition(gpart))

  allocate(amat(nReg,nReg),stat=ier)

	if(ier /= 0) call die(myname_,'allocate()',ier)
	if(mall_ison()) call mall_mci(amat,myname)

! Define the default configuration

  amat(:,:)=1.

  select case(amat_type)
  case(AMAT_BxBT)
    call init_general_(gpart,amat_cutdistance,amat_function,comm)

		! Print message only if the call is successful.

    call mpout_log(myname_,	&
	': AMAT_BxBT, using amat_cutdistance = ',amat_cutdistance)

    select case(amat_function)
    case(TOP_HAT)
      call mpout_log(myname_,'"'//rsrcFUNCTION//'" = TOP_HAT')
    case(WITCHES_HAT)
      call mpout_log(myname_,'"'//rsrcFUNCTION//'" = WITCHES_HAT')
    case(COSINE_HILL)
      call mpout_log(myname_,'"'//rsrcFUNCTION//'" = COSINE_HILL')
    case default
      call die(myname_,'Unknown generating function',amat_function)
    end select

  case default
    call mpout_log(myname_,': AMAT_UNITY, amat_cutdistance = oo')
  end select
 
 end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_general_ - initialize the A matrix version 1x
!
! !DESCRIPTION:
!     subroutine to initialize the masking matrix amat (see above description)
!     which is used to impose sparsity on the PSAS multivariate
!     error covariance matrices.  This one is the general way of generating
!     the A matrix based on specifying the support of the primitive (generating)
!     matrix B using 'cutoff' (in units of cosine of angle between regions).
!
! !INTERFACE:

 subroutine init_general_(gpart,bcutoff,which_function,comm)

   use m_mpout,only : mpout_log,mpout
   use m_mall ,only : mall_ison,mall_mci,mall_mco,mall_ci,mall_co
   use m_die  ,only : die,MP_die,assert_

   use m_GlobalPartition,only : GlobalPartition
   use m_GlobalPartition,only : ptr_Partition
   use m_GlobalPartition,only : inquire

   use m_Spherical_Partition, only : Spherical_Partition
   use m_Spherical_Partition, only : NumberOfRegions
   use m_Spherical_Partition, only : GetRegion
   use m_Spherical_Partition, only : BaseRegion
   use m_Spherical_Partition, only : ListAreas
   use m_Spherical_Partition, only : ListCenters
   use m_Spherical_Partition, only : xyz2reg
   use m_Spherical_Partition, only : Initialize
   use m_Spherical_Partition, only : clean

   use m_Spherical_Triangle , only : Spherical_Triangle, Separation
   use m_Spherical_Triangle , only : SEPANG_CIRCUMCENTER
   use m_Spherical_Triangle , only : SEPANG_CENTER_OF_MASS

   use const,only : Radius_of_Earth	! in meters.

   implicit none

   type(GlobalPartition),intent(in) :: gpart
   real   ,intent(in) :: bcutoff
   integer,intent(in) :: comm
   integer,intent(in) :: which_function

! !REVISION HISTORY:
!       16Oct00 - Peter Lyster <lys@dao.gsfc.nasa.gov> initial
!EOP ___________________________________________________________________

   character(len=*),parameter :: myname_=myname//'::init_general_'

!________________________________________
! Basic configurations:
!
! I. Area Adjustement:

   logical,parameter :: AREA_ADJUSTMENT=.true.

!  The purpose of the area adjustment is to more accurately evaluate the
!  A(x_1,x_2) matrix using an approximation to the integral
!  integral[B(x_1,y)B(x_2,y) dy] where the integral is over the surface
!  of the sphere.  Note that in a matrix form, A=BSB' is p.d. if S is
!  p.d..  In "area_adjustment", dy ~ areas(:)=diag(S).
!
! II. Method of Computing Separation Angle:

   integer,parameter :: METHOD_OF_SEP = SEPANG_CENTER_OF_MASS

! SEPANG_CENTER_OF_MASS specifies that the seperation angle between two
! spherical trangles is based on the distance between their centers of
! mass.
!
!   integer,parameter :: METHOD_OF_SEP = SEPANG_CIRCUMCENTER
!
! SEPANG_CIRCUMCENTER is another choice of the method of computing the
! separation angle between two spherical triangels.  It specifies that
! the separation angle is computed based on the distance between the
! circumcenters of the two triangles.

! Workspace:

   type (Spherical_Partition) :: sp_base
   Type (Spherical_Triangle),allocatable,dimension(:) :: tri
   integer,allocatable,dimension(:) :: kr_base
   real   ,allocatable,dimension(:,:) :: xyz
   real   ,Allocatable,dimension(:) :: areas
   real   ,Allocatable,dimension(:) :: tmpv
   real   ,allocatable,dimension(:,:) :: bmat
   real   ,allocatable,dimension(:) :: am_base

! Local variables:

   integer :: i,j,k,l,ij,n
   integer :: nReg,mReg,nr_base
   integer :: lv_base

   real :: pi,deg
   real :: arg,eps
   real :: cos_sepang,cos_cutoff
   real :: rad_sepang,rad_cutoff

   integer :: ier
!________________________________________

	! Totol number of regions defined by the given partition,
	! either refinable (compress=.false.) or not refinable
	! (compress=.true.).

    nReg=size(amat,1)
	ASSERT(nReg==size(amat,2))
!________________________________________

    if(AREA_ADJUSTMENT) call mpout_log(myname_,	&
	': Using AREA_ADJUSTMENT')

    select case(METHOD_OF_SEP)
    case(SEPANG_CENTER_OF_MASS)
      call mpout_log(myname_,	&
	': Using SEPANG_CENTER_OF_MASS for separation angles')
    case(SEPANG_CIRCUMCENTER)
      call mpout_log(myname_,	&
	': Using SEPANG_CIRCUMCENTER for separation angles')
    case default
      call die(myname_,'unknown METHOD_OF_SEP',METHOD_OF_SEP)
    end select
!________________________________________

	! Create a compressed Spherical_Partition (sp_base) at the
	! _BaseLevel_ (lv_base) for the formation of the base amat
	! (am_base).

    call inquire(gpart,baselevel=lv_base)

    call Initialize(lv_base,partition=sp_base,compress=.true.)

	! Verify if the compressed partition has the same size as
	! the size of the baselevel partition of the refinable
	! partition.

    nr_base=NumberOfRegions(partition=sp_base)

    mReg=NumberOfRegions(atlevel=lv_base,	&
	partition=ptr_Partition(gpart))

	ASSERT(mReg==nr_base)

	! Allocate workspace

    allocate(kr_base(nReg),xyz(3,nReg),			&
	tri(nr_base),areas(nr_base),tmpv(nr_base),	&
	bmat(nr_base,nr_base),am_base(nr_base*(nr_base+1)/2), stat=ier)
	if(ier/=0) call die(myname_,'allocate(bmat,...)',ier)

	if(mall_ison()) then
	  call mall_mci(xyz     ,myname_)
	  call mall_mci(kr_base ,myname_)
	  call mall_ci(size(tri),myname_)
	  call mall_mci(areas   ,myname_)
	  call mall_mci(tmpv    ,myname_)
	  call mall_mci(bmat    ,myname_)
	  call mall_mci(am_base ,myname_)
	endif

! Form the B matrix, which is used as a generating matrix for the
! reduced PSAS masking matrix A.
!________________________________________

! Prepare constants.  Note bcutoff is in kilometers.

    pi  = 4.*atan(1.)
    deg = pi/180.

! bcutoff for determining the support of the generating function for B.

    rad_cutoff = bcutoff*1000./Radius_of_Earth
    cos_cutoff = cos(rad_cutoff)

! Area Adjustement if required.

    areas(:) = 1.
    if(AREA_ADJUSTMENT) call ListAreas(areas,partition=sp_base)

! Get all base-regions

    do j=1,nr_base
      tri(j)=getRegion(j,partition=sp_base)
    end do
!________________________________________

! Compute B'

    bmat(:,:)=0.0
    do i=1,nr_base
      do j=1,nr_base

        rad_sepang = deg*Separation(tri(i),tri(j),method=METHOD_OF_SEP)

			! Although B (bmat) is expected to be
			! symmetric, expressions in this loop are all
			! assumed for defining B^T.

	if(rad_sepang < rad_cutoff) then

          select case(which_function)
          case(TOP_HAT)
            bmat(j,i)= 1.0

          case(WITCHES_HAT)
	    cos_sepang=cos(rad_sepang)
            bmat(j,i)= ( cos_sepang-cos_cutoff )/( 1.0-cos_cutoff )

          case(COSINE_HILL)

	    arg=rad_sepang/rad_cutoff
            bmat(j,i)= 0.5 * ( 1.+cos( pi * arg ) )

          end select
	endif
      enddo
    enddo
!________________________________________

! Compute the upper triangle of A=BSB', note that bmat is symmetric and
! it is B' in the mult-add expression in this code.

    do j=1,nr_base
      l=(j-1)*j/2

	  ! am_base(i,j)=sum{bmat(k,i)*bmat(k,j)*areas(k), k=1,nr_base}

      tmpv(:)=areas(:)*bmat(:,j)
      do i=1,j
	ij=l+i
        am_base(ij)=dot_product(bmat(:,i),tmpv(:))
      enddo
    enddo

! Prepare the normalization fectors: diagonal elements.

    eps=tiny(1.)

    do j=1,nr_base
      l=(j-1)*j/2+j

		! Assertion to avoid deviding by zero

	ASSERT(am_base(l)>eps)

      tmpv(j)=1./sqrt(am_base(l))
    end do

! Normalize the upper triangle of a subset matrix of A

    do j=1,nr_base
      l=(j-1)*j/2
      do i=1,j
	ij=l+i
        am_base(ij) = am_base(ij)*tmpv(i)*tmpv(j)
      enddo
    enddo
!________________________________________

! Map region IDs defined by GlobalPartition to regions defined by the
! base partition, using their CircumCenters.

    call ListCenters(xyz,partition=ptr_Partition(gpart))

    call xyz2reg(nReg,xyz(1,:),xyz(2,:),xyz(3,:),	&
	kr_base(:),partition=sp_base)

! Fill the upper triangle of A by am_base(:)

    do j=1,nReg
      l=(kr_base(j)-1)*kr_base(j)/2
      do i=1,j
	ij=l+kr_base(i)
        amat(i,j)=am_base(ij)
      enddo
    enddo

! Fill the lower triangle of A using its upper triangle.

    do j=1,nReg-1
      do i=j+1,nReg
        amat(i,j)=amat(j,i)
      enddo
    enddo
!________________________________________

    if(CHECK_AMAT) then

      n=0
      do j=1,nr_base
	l=(j-1)*j/2
	if(am_base(l+j) .le. 1.0E-6) n=n+1
	do i=1,j-1
	  ij=l+i
	  if(am_base(ij) .le. 1.0E-6) n=n+2
	end do
      end do

      n=nint(100.*n/(nr_base*nr_base))
      call mpout_log(myname_,'am_base sparsity (%) = ',n)

    endif

! Deallocate workspace

    call clean(sp_base)

	if(mall_ison()) then
	  call mall_mco(xyz     ,myname_)
	  call mall_mco(kr_base ,myname_)
	  call mall_co(size(tri),myname_)
	  call mall_mco(areas   ,myname_)
	  call mall_mco(tmpv    ,myname_)
	  call mall_mco(bmat    ,myname_)
	  call mall_mco(am_base ,myname_)
	endif

    deallocate(kr_base,xyz,tri,areas,tmpv,bmat,am_base,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(bmat,...)',ier)

 end subroutine init_general_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - deallocate the A matrix
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine clean_()
   use m_die ,only : die
   use m_mall,only : mall_ison,mall_mco
   implicit none

! !REVISION HISTORY:
!       27Sep99 - Peter Lyster <lys@dao.gsfc.nasa.gov> initial
!EOP ___________________________________________________________________

   character(len=*),parameter :: myname_=myname//'::clean_'
   integer :: ier

   if(.not.initialized_) call die(myname_,'undefined object')

	if(mall_ison()) call mall_mco(amat,myname)
     deallocate(amat,stat=ier)
	if(ier/=0) call die(myname_,'deallocate()',ier)

   initialized_=.false.
   maxband=-1

 end subroutine clean_

end module m_amat
