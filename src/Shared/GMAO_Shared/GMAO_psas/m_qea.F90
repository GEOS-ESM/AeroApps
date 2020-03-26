!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_qea - Q.E.A parameter settings
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_qea
      implicit none
      private	! except

      public :: npole		! No. of grid points at poles
      public :: eaytresh	! treshold latitude for q.e.a. grid

      public :: lea_beg		! pointers to array segments
      public :: lea_len		! sizes of array segments
      public :: ea_lon		! longitudes of q.e.a. grid points

      public :: j_north
      public :: j_south

      public :: qea_init
      public :: qea_clean

      interface qea_init;  module procedure init_;  end interface
      interface qea_clean; module procedure clean_; end interface

! !REVISION HISTORY:
! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- Changed from qea.h to m_qea(.F90)
!
! 18Jul95 - Jing G.	- Separated from gridxx.h, for Q.E.A. grid
!			  settings only.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_qea'

      integer,parameter :: npole=4	! No. grid points at poles
					! NOTE: In reality, at the
					!       poles all possible
					!       longitudes are 
					!       possible. Leave
					!       npole alone.

      real,save :: eaytresh		! treshold latitude for a
					! quasi-equal-area grid



!     Indices of quasi-equal area  (q.e.a.) and regular lat/lon grid
!_______________________________________________________________________

		! Pointer to begining q.e.a. gridpoint for a given
		! latitude    

      integer,save,dimension(:),allocatable :: lea_beg

		! How many zonal gridpoints at that longitude

      integer,save,dimension(:),allocatable :: lea_len

		! Longitudes of q.e.q. gridpoints (lat not needed right
		! now)

      real,   save,dimension(:),allocatable :: ea_lon

		! From j_south to j_north lat/lon gridpoints coincide
		! with q.e.a. grid points.

      integer,save :: j_north   
      integer,save :: j_south

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - create working data for a q.e.a. grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(idim,jdim)
      use m_die,only : die
      use m_mall,only : mall_mci,mall_ison
      implicit none
      integer,intent(in) :: idim
      integer,intent(in) :: jdim

! !REVISION HISTORY:
! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: ier

  allocate(ea_lon(idim*jdim),lea_beg(jdim),lea_len(jdim),stat=ier)
	if(ier /= 0) call die(myname_,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(ea_lon,myname)
	  call mall_mci(lea_beg,myname)
	  call mall_mci(lea_len,myname)
	endif

end subroutine init_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean the working data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_()
      use m_die,only : die
      use m_mall,only : mall_mco,mall_ison
      implicit none

! !REVISION HISTORY:
! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  integer :: ier

	if(mall_ison()) then
	  call mall_mco(lea_len,myname)
	  call mall_mco(lea_beg,myname)
	  call mall_mco(ea_lon,myname)
	endif

  deallocate(ea_lon,lea_beg,lea_len,stat=ier)
	if(ier /= 0) call die(myname_,'deallocate()',ier)

end subroutine clean_

end module m_qea
