!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_kt_CorO - correlation matrix generators for correl.Obs.Err.
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_kt_CorO
      implicit none
      private	! except

      public :: oHHcor1
      public :: oHHcorx
      public :: sHHmcxpy
      public :: xHHmcxpy
      public :: rHHmcxpy
      public :: sHH1cxpy
      public :: xHH1cxpy
      public :: rHH1cxpy

! !REVISION HISTORY:
! 	23Feb01	- Tom Clune <clune@sgi.com>
!                 Added new interfaces that avoid storing the block
!                 but rather apply elements to right/left vectors as
!                 individual elements are constructed.  Significant
!                 performance improvement.
!	29Aug00	- Jing Guo
!		- initial prototype/prolog/code
!		- combined oHHcor1.F and oHHcorx.F into this module
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_kt_CorO'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: oHHcor1 - forms a packed symmetric correlation matrix
!
! !DESCRIPTION:
!
!	oHHcor1 forms the upper triangle of a symmetric <height,height>
!	correlation matrix.  The returned matrix is in a packed form in
!	the order of ((i=1,j),j=1,len).  The correlation function is
!	separable and in the same form as in the OI of GEOS/DAS.1.2.
!
! !INTERFACE:

    subroutine oHHcor1(len,kx,qr,kl, Mtyp, corrO)
      use m_redwin,only : redwin
      use m_redwin,only : qxWtb,mxWtb
      use m_redwin,only : redwin_initialized
      use m_die,only : die
      implicit none

      integer,intent(in) :: len
      integer,dimension(  len),intent(in) :: kx
      real   ,dimension(3,len),intent(in) :: qr
      integer,dimension(  len),intent(in) :: kl

      character(len=1),intent(inout) :: Mtyp	! "I" or "U", or ..
      real,dimension(len*(len+1)/2),intent(out) :: corrO

! !REVISION HISTORY:
! 	29Aug00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- modified for a unified argument style and as a
!		  module procedure.
! 	09Nov95 - J. Guo	- rewrite the code and the prolog.  The
!		recoding 1) separated forecast error correlation matrix
!		from observation error correlation matrix;  2) used a
!		new structure in hope that it will be optimized better;
!		and 3) getting ready for the new changes (see BUGS:)
!		that may effect the results.
!
!  29Sep95  - Jing G.	- added #ifdef for dynamic dimension check.
!  18Sep95  - Jing G.	- Created multiple copies of statistical tables
!  27Mar95  - Jing G.	- Modified for text file based data tables.
!  19Jan95  - Jing G.	- Added wobs table to replace pindx2().  One
!			  could use rlevs for the same purpose to reduce
!			  the overhead, since rlevs has no real purpose
!			  in this subroutine.
!   7oct94  Meta S.  rearranged observation error correlation test
!                     to avoid index problem
!  10may94  Meta S.  change observation error correlation tabulation
!                     to use cos(dist/rade)
!  28apr94  Meta S.  change forecast error correlation tabulation
!                     to use cos(dist/rade), added prologue
!  02feb94  Meta S.  change radiosonde obs covariance criterion to specify
!                    -matching- lat/lon and alter 'coslim' definition
!                    (now included in "hHHotab.h", "hHHotab.h")
!  15dec93  Meta S.  add observation covariances, covariance tables
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::oHHcor1'

	include	"lvmax.h"

		! observation error class reference
	include "ktmax.h"
	include "kxmax.h"
	include "kxtabl.h"

		! indirect vertical correlation matrix
	include "MX_voecH.h"
	include "voecHH.h"

		! indirect horizotal correlation matrix factors
	include "MX_hoecH.h"
	include "hoecHH.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	integer ij,i,j,l
	integer ki,kj,itau,ktau
	integer kxi,kxj
	integer jhtype,jvtype
	real tau,ctau,xtau,redwin_
	real qri_x,qri_y,qri_z
	real qrj_x,qrj_y,qrj_z
	real vercor,horcor
	real Ocoslimj,qxOtb1j,qxOtb2j,OObeg2j
	logical correlated

!-----------------------------------------------------------------------
	real HH

!=======================================================================
	if(.not.redwin_initialized())	&
		call die(myname_,'redwin not initialized')

	if(len.le.0) return

	corrO(1)=1.

	do j=2,len

	  ij=j*(j-1)/2

	  kxj=kx(j) ! the true kx code

			! the type of the horizontal correlation
			! function and the type of the vertical
			! correlation function

	  jhtype=i_hoecHc(kxj)
	  jvtype=i_voecHc(kxj)

	  correlated=jhtype.gt.0.and.jhtype.le.MX_hoecH .and.	&
		     jvtype.gt.0.and.jvtype.le.MX_voecH

			! (that is horizontally correlated)

	  if(correlated) then

		! a preliminary check if there is any (other than the
		! diagonals with the same kx).

	    correlated=.false.
	    i=1
	    do while(i.lt.j.and..not.correlated)
	      kxi=kx(i)
	      correlated=kxi.eq.kxj
	      i=i+1
	    end do
	  endif

	  if(correlated) then	! at least for some elements
	    Mtyp='U'		! a packed matrix, not 'I' any more

	    qrj_x=qr(1,j)
	    qrj_y=qr(2,j)
	    qrj_z=qr(3,j)

	    kj = kl(j)

	    Ocoslimj=Ocoslim(jhtype)
	    qxOtb1j=qxOtb1(jhtype)
	    qxOtb2j=qxOtb2(jhtype)
	    OObeg2j=OObeg2(jhtype)

	    do i=1,j-1

	      l=ij+i	! packed matrix index to the elements in the
			! upper triangle of the matrix

			! cossep is the distance between i and j in
			! cos(great_circle_angular_distance)

	      corrO(l)=0.

	      kxi=kx(i)
	      if(kxi.eq.kxj) then

		 qri_x = qr(1,i)
		 qri_y = qr(2,i)
		 qri_z = qr(3,i)
		 tau = qri_x*qrj_x + qri_y*qrj_y + qri_z*qrj_z

		horcor=0.
	        if(tau.gt.Ocoslimj) then


		  ki=kl(i)

			! table indexing

		  ctau=1.-tau
#ifdef _LINEAR
		  xtau=qxOtb1j*ctau+1.
		  if(ctau.gt.OObeg2j)		&
			xtau=qxOtb2j*(ctau-OObeg2j)+nOOtb1+1.
		  itau=min(int(xtau),nHHotab-1)
		  xtau=xtau-itau
		  HH=hoecHH(ki,itau,jhtype)+hoecHH(kj,itau,jhtype)
		  HH=HH + xtau*			&
		    (hoecHH(ki,itau+1,jhtype)+	&
		     hoecHH(kj,itau+1,jhtype)-	HH)
#else
		  xtau=qxOtb1j*ctau+1.5
		  if(ctau.gt.OObeg2j)		&
			xtau=qxOtb2j*(ctau-OObeg2j)+nOOtb1+1.5

		  itau=int(xtau)	! in [1,nHHotab] is assumed

			! a non-separable correlation function form

		  HH=hoecHH(ki,itau,jhtype)+hoecHH(kj,itau,jhtype)
#endif
		  horcor = HH

		  ktau=int(qxWtb*ctau+.5)+1
		  redwin_=0.
		  if(ktau<=mxWtb) redwin_=redwin(ktau)

		  vercor=voecHH(ki,kj,jvtype)*redwin_

		  corrO(l)=horcor*vercor

		endif	! tau.gt.Ocoslimj
	      endif	! kxi.eq.kxj
	    end do		! i=1,j-1
	  else
	    do i=1,j-1
	      corrO(ij+i)=0.
	    end do
	  endif

	  corrO(ij+j)=1.	! where i=j
	end do		! j=2,len

end subroutine oHHcor1

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sHHcxpy - applies a symmetric correlation matrix
!
! !DESCRIPTION:
!
!	sHHcxpy applies a symmetric <height,height> correlation matrix
!       to right and left vectors.  The correlation function is separable
!       and in the same form as in the OI of GEOS/DAS.1.2.  Symmery is used
!       to reduce computation by roughly 50%.
!
! !INTERFACE:

    subroutine sHHmcxpy(len,kx,qr,kl, nvecs, x, y)
      use m_redwin,only : redwin
      use m_redwin,only : qxWtb,mxWtb
      use m_redwin,only : redwin_initialized
      use m_die,only : die
      implicit none

      integer,intent(in) :: len
      integer,dimension(  len),intent(in) :: kx
      real   ,dimension(3,len),intent(in) :: qr
      integer,dimension(  len),intent(in) :: kl

      integer, intent(in) :: nvecs
      real,dimension(nvecs, len), intent(in)    :: x
      real,dimension(nvecs, len), intent(inout) :: y

! !REVISION HISTORY:
! 	23Feb01	- Tom Clune <clune@sgi.com>
!                 Initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sHHmcxpy'

	include	"lvmax.h"

		! observation error class reference
	include "ktmax.h"
	include "kxmax.h"
	include "kxtabl.h"

		! indirect vertical correlation matrix
	include "MX_voecH.h"
	include "voecHH.h"

		! indirect horizotal correlation matrix factors
	include "MX_hoecH.h"
	include "hoecHH.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	integer ij,i,j,l
	integer ki,kj,itau,ktau
	integer kxi,kxj
	integer jhtype,jvtype
	real tau,ctau,xtau,redwin_
	real qri_x,qri_y,qri_z
	real qrj_x,qrj_y,qrj_z
	real vercor,horcor
	real Ocoslimj,qxOtb1j,qxOtb2j,OObeg2j
	logical correlated

!-----------------------------------------------------------------------
	real HH, corrO

!=======================================================================
	if(.not.redwin_initialized())	&
		call die(myname_,'redwin not initialized')

	if(len.le.0) return

	y(:,1) = y(:,1) + x(:,1)

	do j=2,len

	  ij=j*(j-1)/2

	  kxj=kx(j) ! the true kx code

			! the type of the horizontal correlation
			! function and the type of the vertical
			! correlation function

	  jhtype=i_hoecHc(kxj)
	  jvtype=i_voecHc(kxj)

	  correlated=jhtype.gt.0.and.jhtype.le.MX_hoecH .and.	&
		     jvtype.gt.0.and.jvtype.le.MX_voecH

			! (that is horizontally correlated)

	  if(correlated) then

		! a preliminary check if there is any (other than the
		! diagonals with the same kx).

	    correlated=.false.
	    i=1
	    do while(i.lt.j.and..not.correlated)
	      kxi=kx(i)
	      correlated=kxi.eq.kxj
	      i=i+1
	    end do
	  endif

	  if(correlated) then	! at least for some elements

	    qrj_x=qr(1,j)
	    qrj_y=qr(2,j)
	    qrj_z=qr(3,j)

	    kj = kl(j)

	    Ocoslimj=Ocoslim(jhtype)
	    qxOtb1j=qxOtb1(jhtype)
	    qxOtb2j=qxOtb2(jhtype)
	    OObeg2j=OObeg2(jhtype)

	    do i=1,j-1

	      l=ij+i	! packed matrix index to the elements in the
			! upper triangle of the matrix

			! cossep is the distance between i and j in
			! cos(great_circle_angular_distance)

	      kxi=kx(i)
	      if(kxi.eq.kxj) then

		 qri_x = qr(1,i)
		 qri_y = qr(2,i)
		 qri_z = qr(3,i)
		 tau = qri_x*qrj_x + qri_y*qrj_y + qri_z*qrj_z

		horcor=0.
	        if(tau.gt.Ocoslimj) then


		  ki=kl(i)

			! table indexing

		  ctau=1.-tau
#ifdef _LINEAR
		  xtau=qxOtb1j*ctau+1.
		  if(ctau.gt.OObeg2j)		&
			xtau=qxOtb2j*(ctau-OObeg2j)+nOOtb1+1.
		  itau=min(int(xtau),nHHotab-1)
		  xtau=xtau-itau
		  HH=hoecHH(ki,itau,jhtype)+hoecHH(kj,itau,jhtype)
		  HH=HH + xtau*			&
		    (hoecHH(ki,itau+1,jhtype)+	&
		     hoecHH(kj,itau+1,jhtype)-	HH)
#else
		  xtau=qxOtb1j*ctau+1.5
		  if(ctau.gt.OObeg2j)		&
			xtau=qxOtb2j*(ctau-OObeg2j)+nOOtb1+1.5

		  itau=int(xtau)	! in [1,nHHotab] is assumed

			! a non-separable correlation function form

		  HH=hoecHH(ki,itau,jhtype)+hoecHH(kj,itau,jhtype)
#endif
		  horcor = HH

		  ktau=int(qxWtb*ctau+.5)+1
		  redwin_=0.
		  if(ktau<=mxWtb) redwin_=redwin(ktau)

		  vercor=voecHH(ki,kj,jvtype)*redwin_

		  corrO = horcor*vercor
		  y(:,i) = y(:,i) + corrO * x(:,j)
		  y(:,j) = y(:,j) + corrO * x(:,i)

		endif	! tau.gt.Ocoslimj
	      endif	! kxi.eq.kxj
	    end do		! i=1,j-1
	  else
	     ! no contribution
	  endif
	  y(:,j) = y(:,j) + x(:,j) ! where i=j
	end do		! j=2,len

      end subroutine sHHmcxpy

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: oHHcorx - forms a <H,H> observation correlation matrix
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine oHHcorx(leni,kxi,qri,kli,lenj,kxj,qrj,klj, Mtyp,corrO)
      use m_redwin,only : redwin
      use m_redwin,only : qxWtb,mxWtb
      use m_redwin,only : redwin_initialized
      use m_die,only : die
      implicit none

      integer,intent(in) :: leni
      integer,dimension(  leni),intent(in) :: kxi
      real   ,dimension(3,leni),intent(in) :: qri
      integer,dimension(  leni),intent(in) :: kli

      integer,intent(in) :: lenj
      integer,dimension(  lenj),intent(in) :: kxj
      real   ,dimension(3,lenj),intent(in) :: qrj
      integer,dimension(  lenj),intent(in) :: klj

      character(len=1),intent(inout) :: Mtyp	! "N"ormal or "Z"ero
      real,dimension(leni,lenj),intent(out) :: corrO ! <H,H> correlation

! !REVISION HISTORY:
! 	29Aug00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		. New prolog since the old revision history text is
!		  a mess in the earlier version.
!		. Unified the argument style and modified to be a
!		  module procedure.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::oHHcorx'

	include "lvmax.h"

		! observation error class references
	include "ktmax.h"
	include "kxmax.h"
	include "kxtabl.h"

		! indirect vertical correlation matrix
	include "MX_voecH.h"
	include "voecHH.h"

		! indirect horizontal correlaiton matrix factors
	include "MX_hoecH.h"
	include "hoecHH.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	integer ij,i,j,l
	integer ki,kj,itau,ktau
	integer kx_i,kx_j
	integer jvtype,jhtype
	real tau,ctau,xtau,redwin_
	real qri_x,qri_y,qri_z
	real qrj_x,qrj_y,qrj_z
	real vercor,horcor
	real Ocoslimj,qxOtb1j,qxOtb2j,OObeg2j
	logical correlated

!-----------------------------------------------------------------------
	real HH

!=======================================================================
	if(.not.redwin_initialized())	&
		call die(myname_,'redwin not initialized')

	if(leni.le.0.or.lenj.le.0) return

	Mtyp='Z'
	do j=1,lenj

	   kx_j = kxj(j)

			! types of the horizontal correlation functions
			! types of the vertical correlation functions

	  jhtype=i_hoecHc(kx_j)
	  jvtype=i_voecHc(kx_j)

	  correlated=jhtype.gt.0.and.jhtype.le.MX_hoecH .and.	&
		     jvtype.gt.0.and.jvtype.le.MX_voecH

	  if(correlated) then

		! a preliminary check if there is any (other than the
		! diagonals with the same correlation ID (pix).

	    correlated=.false.
	    i=1
	    do while(i.le.leni.and..not.correlated)
	      kx_i=kxi(i)
	      correlated=kx_i.eq.kx_j
	      i=i+1
	    end do
	  endif

	  if(correlated) then
	    Mtyp='N'		! a "N"onzero matrix

	    qrj_x=qrj(1,j)
	    qrj_y=qrj(2,j)
	    qrj_z=qrj(3,j)


	    kj=klj(j)

	    Ocoslimj=Ocoslim(jhtype)
	    qxOtb1j=qxOtb1(jhtype)
	    qxOtb2j=qxOtb2(jhtype)
	    OObeg2j=OObeg2(jhtype)

	    do i=1,leni

			! cossep is the distance between i and j in
			! cos(great_circle_angular_distance)

	      corrO(i,j)=0.

	      kx_i=kxi(i)
	      if(kx_i.eq.kx_j) then

		 qri_x = qri(1,i)
		 qri_y = qri(2,i)
		 qri_z = qri(3,i)
		 tau = qri_x*qrj_x + qri_y*qrj_y + qri_z*qrj_z

		horcor=0.
	        if(tau.gt.Ocoslimj) then

		   ki = kli(i)

			! table indexing

		  ctau=1.-tau
#ifdef _LINEAR
		  xtau=qxOtb1j*ctau+1.
		  if(ctau.gt.OObeg2j)		&
			xtau=qxOtb2j*(ctau-OObeg2j)+nOOtb1+1.
		  itau=min(int(xtau),nHHotab-1)
		  xtau=xtau-itau
		  HH=hoecHH(ki,itau,jhtype)+hoecHH(kj,itau,jhtype)
		  HH=HH + xtau*			&
		    (hoecHH(ki,itau+1,jhtype)+	&
		     hoecHH(kj,itau+1,jhtype)-	HH)
#else
		  xtau=qxOtb1j*ctau+1.5
		  if(ctau.gt.OObeg2j)		&
			xtau=qxOtb2j*(ctau-OObeg2j)+nOOtb1+1.5

		  itau=int(xtau)	! in [1,nHHotab] is assumed

			! a non-separable correlation function form

		  HH=hoecHH(ki,itau,jhtype)+hoecHH(kj,itau,jhtype)
#endif
		  horcor = HH

		  ktau=int(qxWtb*ctau+.5)+1
		  redwin_=0.
		  if(ktau<=mxWtb) redwin_=redwin(ktau)

		  vercor=voecHH(ki,kj,jvtype)*redwin_

		  corrO(i,j)=horcor*vercor
		endif	! tau.gt.Ocoslimj
	      endif	! kx_i.eq.kx_j
	    end do		! i=1,leni

	  else
	    do i=1,leni
	      corrO(i,j)=0.
	    end do
	  endif

	end do		! j=1,lenj

end subroutine oHHcorx
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: xHHcxpy - applies a <H,H> observation correlation matrix
!
! !DESCRIPTION:
!               applies block to right and left vectors.
!
! !INTERFACE:

    subroutine xHHmcxpy(leni,kxi,qri,kli,lenj,kxj,qrj,klj, nvecs,xi,yj,xj,yi)
      use m_redwin,only : redwin
      use m_redwin,only : qxWtb,mxWtb
      use m_redwin,only : redwin_initialized
      use m_die,only : die
      implicit none

      integer,intent(in) :: leni
      integer,dimension(  leni),intent(in) :: kxi
      real   ,dimension(3,leni),intent(in) :: qri
      integer,dimension(  leni),intent(in) :: kli

      integer,intent(in) :: lenj
      integer,dimension(  lenj),intent(in) :: kxj
      real   ,dimension(3,lenj),intent(in) :: qrj
      integer,dimension(  lenj),intent(in) :: klj

      integer, intent(in) :: nvecs
      real,dimension(nvecs,leni), intent(in)    :: xi
      real,dimension(nvecs,lenj), intent(inout) :: yj
      real,dimension(nvecs,lenj), intent(in)    :: xj
      real,dimension(nvecs,leni), intent(inout) :: yi

! !REVISION HISTORY:
! 	23Feb01	- Tom Clune <clune@sgi.com>
!                 Initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::xHHmcxpy'

	include "lvmax.h"

		! observation error class references
	include "ktmax.h"
	include "kxmax.h"
	include "kxtabl.h"

		! indirect vertical correlation matrix
	include "MX_voecH.h"
	include "voecHH.h"

		! indirect horizontal correlaiton matrix factors
	include "MX_hoecH.h"
	include "hoecHH.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	integer ij,i,j,l
	integer ki,kj,itau,ktau
	integer kx_i,kx_j
	integer jvtype,jhtype
	real tau,ctau,xtau,redwin_
	real qri_x,qri_y,qri_z
	real qrj_x,qrj_y,qrj_z
	real vercor,horcor
	real Ocoslimj,qxOtb1j,qxOtb2j,OObeg2j
	logical correlated

!-----------------------------------------------------------------------
	real HH, corrO

!=======================================================================
	if(.not.redwin_initialized())	&
		call die(myname_,'redwin not initialized')

	if(leni.le.0.or.lenj.le.0) return

	do j=1,lenj

	   kx_j = kxj(j)

			! types of the horizontal correlation functions
			! types of the vertical correlation functions

	  jhtype=i_hoecHc(kx_j)
	  jvtype=i_voecHc(kx_j)

	  correlated=jhtype.gt.0.and.jhtype.le.MX_hoecH .and.	&
		     jvtype.gt.0.and.jvtype.le.MX_voecH

	  if(correlated) then

		! a preliminary check if there is any (other than the
		! diagonals with the same correlation ID (pix).

	    correlated=.false.
	    i=1
	    do while(i.le.leni.and..not.correlated)
	      kx_i=kxi(i)
	      correlated=kx_i.eq.kx_j
	      i=i+1
	    end do
	  endif

	  if(correlated) then

	    qrj_x=qrj(1,j)
	    qrj_y=qrj(2,j)
	    qrj_z=qrj(3,j)


	    kj=klj(j)

	    Ocoslimj=Ocoslim(jhtype)
	    qxOtb1j=qxOtb1(jhtype)
	    qxOtb2j=qxOtb2(jhtype)
	    OObeg2j=OObeg2(jhtype)

	    do i=1,leni

			! cossep is the distance between i and j in
			! cos(great_circle_angular_distance)

	      kx_i=kxi(i)
	      if(kx_i.eq.kx_j) then

		 qri_x = qri(1,i)
		 qri_y = qri(2,i)
		 qri_z = qri(3,i)
		 tau = qri_x*qrj_x + qri_y*qrj_y + qri_z*qrj_z

		horcor=0.
	        if(tau.gt.Ocoslimj) then

		   ki = kli(i)

			! table indexing

		  ctau=1.-tau
#ifdef _LINEAR
		  xtau=qxOtb1j*ctau+1.
		  if(ctau.gt.OObeg2j)		&
			xtau=qxOtb2j*(ctau-OObeg2j)+nOOtb1+1.
		  itau=min(int(xtau),nHHotab-1)
		  xtau=xtau-itau
		  HH=hoecHH(ki,itau,jhtype)+hoecHH(kj,itau,jhtype)
		  HH=HH + xtau*			&
		    (hoecHH(ki,itau+1,jhtype)+	&
		     hoecHH(kj,itau+1,jhtype)-	HH)
#else
		  xtau=qxOtb1j*ctau+1.5
		  if(ctau.gt.OObeg2j)		&
			xtau=qxOtb2j*(ctau-OObeg2j)+nOOtb1+1.5

		  itau=int(xtau)	! in [1,nHHotab] is assumed

			! a non-separable correlation function form

		  HH=hoecHH(ki,itau,jhtype)+hoecHH(kj,itau,jhtype)
#endif
		  horcor = HH

		  ktau=int(qxWtb*ctau+.5)+1
		  redwin_=0.
		  if(ktau<=mxWtb) redwin_=redwin(ktau)

		  vercor=voecHH(ki,kj,jvtype)*redwin_
		  
		  corrO = horcor*vercor

		  yj(:,j) = yj(:,j) + corrO * xi(:,i)
		  yi(:,i) = yi(:,i) + corrO * xj(:,j)

		endif	! tau.gt.Ocoslimj
	      endif	! kx_i.eq.kx_j
	    end do		! i=1,leni

	  else
	     ! do nothing
	  endif

	end do		! j=1,lenj

end subroutine xHHmcxpy
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rHHcxpy - applies a <H,H> observation correlation matrix
!
! !DESCRIPTION:
!                applies block to right vector
!
! !INTERFACE:

    subroutine rHHmcxpy(leni,kxi,qri,kli,lenj,kxj,qrj,klj, nvecs,x,y)
      use m_redwin,only : redwin
      use m_redwin,only : qxWtb,mxWtb
      use m_redwin,only : redwin_initialized
      use m_die,only : die
      implicit none

      integer,intent(in) :: leni
      integer,dimension(  leni),intent(in) :: kxi
      real   ,dimension(3,leni),intent(in) :: qri
      integer,dimension(  leni),intent(in) :: kli

      integer,intent(in) :: lenj
      integer,dimension(  lenj),intent(in) :: kxj
      real   ,dimension(3,lenj),intent(in) :: qrj
      integer,dimension(  lenj),intent(in) :: klj

      integer, intent(in) :: nvecs
      real,dimension(nvecs,lenj), intent(in)    :: x
      real,dimension(nvecs,leni), intent(inout) :: y

! !REVISION HISTORY:
! 	23Feb01	- Tom Clune <clune@sgi.com>
!                 Initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rHHmcxpy'

	include "lvmax.h"

		! observation error class references
	include "ktmax.h"
	include "kxmax.h"
	include "kxtabl.h"

		! indirect vertical correlation matrix
	include "MX_voecH.h"
	include "voecHH.h"

		! indirect horizontal correlaiton matrix factors
	include "MX_hoecH.h"
	include "hoecHH.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	integer ij,i,j,l
	integer ki,kj,itau,ktau
	integer kx_i,kx_j
	integer jvtype,jhtype
	real tau,ctau,xtau,redwin_
	real qri_x,qri_y,qri_z
	real qrj_x,qrj_y,qrj_z
	real vercor,horcor
	real Ocoslimj,qxOtb1j,qxOtb2j,OObeg2j
	logical correlated

!-----------------------------------------------------------------------
	real HH, corrO

!=======================================================================
	if(.not.redwin_initialized())	&
		call die(myname_,'redwin not initialized')

	if(leni.le.0.or.lenj.le.0) return

	do j=1,lenj

	   kx_j = kxj(j)

			! types of the horizontal correlation functions
			! types of the vertical correlation functions

	  jhtype=i_hoecHc(kx_j)
	  jvtype=i_voecHc(kx_j)

	  correlated=jhtype.gt.0.and.jhtype.le.MX_hoecH .and.	&
		     jvtype.gt.0.and.jvtype.le.MX_voecH

	  if(correlated) then

		! a preliminary check if there is any (other than the
		! diagonals with the same correlation ID (pix).

	    correlated=.false.
	    i=1
	    do while(i.le.leni.and..not.correlated)
	      kx_i=kxi(i)
	      correlated=kx_i.eq.kx_j
	      i=i+1
	    end do
	  endif

	  if(correlated) then

	    qrj_x=qrj(1,j)
	    qrj_y=qrj(2,j)
	    qrj_z=qrj(3,j)


	    kj=klj(j)

	    Ocoslimj=Ocoslim(jhtype)
	    qxOtb1j=qxOtb1(jhtype)
	    qxOtb2j=qxOtb2(jhtype)
	    OObeg2j=OObeg2(jhtype)

	    do i=1,leni

			! cossep is the distance between i and j in
			! cos(great_circle_angular_distance)

	      kx_i=kxi(i)
	      if(kx_i.eq.kx_j) then

		 qri_x = qri(1,i)
		 qri_y = qri(2,i)
		 qri_z = qri(3,i)
		 tau = qri_x*qrj_x + qri_y*qrj_y + qri_z*qrj_z

		horcor=0.
	        if(tau.gt.Ocoslimj) then

		   ki = kli(i)

			! table indexing

		  ctau=1.-tau
#ifdef _LINEAR
		  xtau=qxOtb1j*ctau+1.
		  if(ctau.gt.OObeg2j)		&
			xtau=qxOtb2j*(ctau-OObeg2j)+nOOtb1+1.
		  itau=min(int(xtau),nHHotab-1)
		  xtau=xtau-itau
		  HH=hoecHH(ki,itau,jhtype)+hoecHH(kj,itau,jhtype)
		  HH=HH + xtau*			&
		    (hoecHH(ki,itau+1,jhtype)+	&
		     hoecHH(kj,itau+1,jhtype)-	HH)
#else
		  xtau=qxOtb1j*ctau+1.5
		  if(ctau.gt.OObeg2j)		&
			xtau=qxOtb2j*(ctau-OObeg2j)+nOOtb1+1.5

		  itau=int(xtau)	! in [1,nHHotab] is assumed

			! a non-separable correlation function form

		  HH=hoecHH(ki,itau,jhtype)+hoecHH(kj,itau,jhtype)
#endif
		  horcor = HH

		  ktau=int(qxWtb*ctau+.5)+1
		  redwin_=0.
		  if(ktau<=mxWtb) redwin_=redwin(ktau)

		  vercor=voecHH(ki,kj,jvtype)*redwin_

		  
		  corrO = horcor*vercor
		  y(:,i) = y(:,i) + corrO * x(:,j)
		endif	! tau.gt.Ocoslimj
	      endif	! kx_i.eq.kx_j
	    end do		! i=1,leni

	  else
	     ! do nothing
	  endif

	end do		! j=1,lenj

end subroutine rHHmcxpy
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sHHcxpy - applies a symmetric correlation matrix
!
! !DESCRIPTION:
!
!	sHHcxpy applies a symmetric <height,height> correlation matrix
!       to right and left vectors.  The correlation function is separable
!       and in the same form as in the OI of GEOS/DAS.1.2.  Symmery is used
!       to reduce computation by roughly 50%.
!
! !INTERFACE:

    subroutine sHH1cxpy(len,kx,qr,kl, nvecs, x, y)
      use m_redwin,only : redwin
      use m_redwin,only : qxWtb,mxWtb
      use m_redwin,only : redwin_initialized
      use m_die,only : die
      implicit none

      integer,intent(in) :: len
      integer,dimension(  len),intent(in) :: kx
      real   ,dimension(3,len),intent(in) :: qr
      integer,dimension(  len),intent(in) :: kl

      integer, intent(in) :: nvecs
      real,dimension(len), intent(in)    :: x
      real,dimension(len), intent(inout) :: y

! !REVISION HISTORY:
! 	23Feb01	- Tom Clune <clune@sgi.com>
!                 Initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sHH1cxpy'

	include	"lvmax.h"

		! observation error class reference
	include "ktmax.h"
	include "kxmax.h"
	include "kxtabl.h"

		! indirect vertical correlation matrix
	include "MX_voecH.h"
	include "voecHH.h"

		! indirect horizotal correlation matrix factors
	include "MX_hoecH.h"
	include "hoecHH.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	integer ij,i,j,l
	integer ki,kj,itau,ktau
	integer kxi,kxj
	integer jhtype,jvtype
	real tau,ctau,xtau,redwin_
	real qri_x,qri_y,qri_z
	real qrj_x,qrj_y,qrj_z
	real vercor,horcor
	real Ocoslimj,qxOtb1j,qxOtb2j,OObeg2j
	logical correlated

!-----------------------------------------------------------------------
	real HH, corrO

!=======================================================================
	if(.not.redwin_initialized())	&
		call die(myname_,'redwin not initialized')

	if(len.le.0) return

	y(1) = y(1) + x(1)

	do j=2,len

	  ij=j*(j-1)/2

	  kxj=kx(j) ! the true kx code

			! the type of the horizontal correlation
			! function and the type of the vertical
			! correlation function

	  jhtype=i_hoecHc(kxj)
	  jvtype=i_voecHc(kxj)

	  correlated=jhtype.gt.0.and.jhtype.le.MX_hoecH .and.	&
		     jvtype.gt.0.and.jvtype.le.MX_voecH

			! (that is horizontally correlated)

	  if(correlated) then

		! a preliminary check if there is any (other than the
		! diagonals with the same kx).

	    correlated=.false.
	    i=1
	    do while(i.lt.j.and..not.correlated)
	      kxi=kx(i)
	      correlated=kxi.eq.kxj
	      i=i+1
	    end do
	  endif

	  if(correlated) then	! at least for some elements

	    qrj_x=qr(1,j)
	    qrj_y=qr(2,j)
	    qrj_z=qr(3,j)

	    kj = kl(j)

	    Ocoslimj=Ocoslim(jhtype)
	    qxOtb1j=qxOtb1(jhtype)
	    qxOtb2j=qxOtb2(jhtype)
	    OObeg2j=OObeg2(jhtype)

	    do i=1,j-1

	      l=ij+i	! packed matrix index to the elements in the
			! upper triangle of the matrix

			! cossep is the distance between i and j in
			! cos(great_circle_angular_distance)

	      kxi=kx(i)
	      if(kxi.eq.kxj) then

		 qri_x = qr(1,i)
		 qri_y = qr(2,i)
		 qri_z = qr(3,i)
		 tau = qri_x*qrj_x + qri_y*qrj_y + qri_z*qrj_z

		horcor=0.
	        if(tau.gt.Ocoslimj) then


		  ki=kl(i)

			! table indexing

		  ctau=1.-tau
#ifdef _LINEAR
		  xtau=qxOtb1j*ctau+1.
		  if(ctau.gt.OObeg2j)		&
			xtau=qxOtb2j*(ctau-OObeg2j)+nOOtb1+1.
		  itau=min(int(xtau),nHHotab-1)
		  xtau=xtau-itau
		  HH=hoecHH(ki,itau,jhtype)+hoecHH(kj,itau,jhtype)
		  HH=HH + xtau*			&
		    (hoecHH(ki,itau+1,jhtype)+	&
		     hoecHH(kj,itau+1,jhtype)-	HH)
#else
		  xtau=qxOtb1j*ctau+1.5
		  if(ctau.gt.OObeg2j)		&
			xtau=qxOtb2j*(ctau-OObeg2j)+nOOtb1+1.5

		  itau=int(xtau)	! in [1,nHHotab] is assumed

			! a non-separable correlation function form

		  HH=hoecHH(ki,itau,jhtype)+hoecHH(kj,itau,jhtype)
#endif
		  horcor = HH

		  ktau=int(qxWtb*ctau+.5)+1
		  redwin_=0.
		  if(ktau<=mxWtb) redwin_=redwin(ktau)

		  vercor=voecHH(ki,kj,jvtype)*redwin_

		  corrO = horcor*vercor
		  y(i) = y(i) + corrO * x(j)
		  y(j) = y(j) + corrO * x(i)

		endif	! tau.gt.Ocoslimj
	      endif	! kxi.eq.kxj
	    end do		! i=1,j-1
	  else
	     ! no contribution
	  endif
	  y(j) = y(j) + x(j) ! where i=j
	end do		! j=2,len

      end subroutine sHH1cxpy

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: xHHcxpy - applies a <H,H> observation correlation matrix
!
! !DESCRIPTION:
!               applies block to right and left vectors.
!
! !INTERFACE:

    subroutine xHH1cxpy(leni,kxi,qri,kli,lenj,kxj,qrj,klj, nvecs,xi,yj,xj,yi)
      use m_redwin,only : redwin
      use m_redwin,only : qxWtb,mxWtb
      use m_redwin,only : redwin_initialized
      use m_die,only : die
      implicit none

      integer,intent(in) :: leni
      integer,dimension(  leni),intent(in) :: kxi
      real   ,dimension(3,leni),intent(in) :: qri
      integer,dimension(  leni),intent(in) :: kli

      integer,intent(in) :: lenj
      integer,dimension(  lenj),intent(in) :: kxj
      real   ,dimension(3,lenj),intent(in) :: qrj
      integer,dimension(  lenj),intent(in) :: klj

      integer, intent(in) :: nvecs
      real,dimension(leni), intent(in)    :: xi
      real,dimension(lenj), intent(inout) :: yj
      real,dimension(lenj), intent(in)    :: xj
      real,dimension(leni), intent(inout) :: yi

! !REVISION HISTORY:
! 	23Feb01	- Tom Clune <clune@sgi.com>
!                 Initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::xHH1cxpy'

	include "lvmax.h"

		! observation error class references
	include "ktmax.h"
	include "kxmax.h"
	include "kxtabl.h"

		! indirect vertical correlation matrix
	include "MX_voecH.h"
	include "voecHH.h"

		! indirect horizontal correlaiton matrix factors
	include "MX_hoecH.h"
	include "hoecHH.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	integer ij,i,j,l
	integer ki,kj,itau,ktau
	integer kx_i,kx_j
	integer jvtype,jhtype
	real tau,ctau,xtau,redwin_
	real qri_x,qri_y,qri_z
	real qrj_x,qrj_y,qrj_z
	real vercor,horcor
	real Ocoslimj,qxOtb1j,qxOtb2j,OObeg2j
	logical correlated

!-----------------------------------------------------------------------
	real HH, corrO

!=======================================================================
	if(.not.redwin_initialized())	&
		call die(myname_,'redwin not initialized')

	if(leni.le.0.or.lenj.le.0) return

	do j=1,lenj

	   kx_j = kxj(j)

			! types of the horizontal correlation functions
			! types of the vertical correlation functions

	  jhtype=i_hoecHc(kx_j)
	  jvtype=i_voecHc(kx_j)

	  correlated=jhtype.gt.0.and.jhtype.le.MX_hoecH .and.	&
		     jvtype.gt.0.and.jvtype.le.MX_voecH

	  if(correlated) then

		! a preliminary check if there is any (other than the
		! diagonals with the same correlation ID (pix).

	    correlated=.false.
	    i=1
	    do while(i.le.leni.and..not.correlated)
	      kx_i=kxi(i)
	      correlated=kx_i.eq.kx_j
	      i=i+1
	    end do
	  endif

	  if(correlated) then

	    qrj_x=qrj(1,j)
	    qrj_y=qrj(2,j)
	    qrj_z=qrj(3,j)


	    kj=klj(j)

	    Ocoslimj=Ocoslim(jhtype)
	    qxOtb1j=qxOtb1(jhtype)
	    qxOtb2j=qxOtb2(jhtype)
	    OObeg2j=OObeg2(jhtype)

	    do i=1,leni

			! cossep is the distance between i and j in
			! cos(great_circle_angular_distance)

	      kx_i=kxi(i)
	      if(kx_i.eq.kx_j) then

		 qri_x = qri(1,i)
		 qri_y = qri(2,i)
		 qri_z = qri(3,i)
		 tau = qri_x*qrj_x + qri_y*qrj_y + qri_z*qrj_z

		horcor=0.
	        if(tau.gt.Ocoslimj) then

		   ki = kli(i)

			! table indexing

		  ctau=1.-tau
#ifdef _LINEAR
		  xtau=qxOtb1j*ctau+1.
		  if(ctau.gt.OObeg2j)		&
			xtau=qxOtb2j*(ctau-OObeg2j)+nOOtb1+1.
		  itau=min(int(xtau),nHHotab-1)
		  xtau=xtau-itau
		  HH=hoecHH(ki,itau,jhtype)+hoecHH(kj,itau,jhtype)
		  HH=HH + xtau*			&
		    (hoecHH(ki,itau+1,jhtype)+	&
		     hoecHH(kj,itau+1,jhtype)-	HH)
#else
		  xtau=qxOtb1j*ctau+1.5
		  if(ctau.gt.OObeg2j)		&
			xtau=qxOtb2j*(ctau-OObeg2j)+nOOtb1+1.5

		  itau=int(xtau)	! in [1,nHHotab] is assumed

			! a non-separable correlation function form

		  HH=hoecHH(ki,itau,jhtype)+hoecHH(kj,itau,jhtype)
#endif
		  horcor = HH

		  ktau=int(qxWtb*ctau+.5)+1
		  redwin_=0.
		  if(ktau<=mxWtb) redwin_=redwin(ktau)

		  vercor=voecHH(ki,kj,jvtype)*redwin_
		  
		  corrO = horcor*vercor

		  yj(j) = yj(j) + corrO * xi(i)
		  yi(i) = yi(i) + corrO * xj(j)

		endif	! tau.gt.Ocoslimj
	      endif	! kx_i.eq.kx_j
	    end do		! i=1,leni

	  else
	     ! do nothing
	  endif

	end do		! j=1,lenj

end subroutine xHH1cxpy
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rHHcxpy - applies a <H,H> observation correlation matrix
!
! !DESCRIPTION:
!                applies block to right vector
!
! !INTERFACE:

    subroutine rHH1cxpy(leni,kxi,qri,kli,lenj,kxj,qrj,klj, nvecs,x,y)
      use m_redwin,only : redwin
      use m_redwin,only : qxWtb,mxWtb
      use m_redwin,only : redwin_initialized
      use m_die,only : die
      implicit none

      integer,intent(in) :: leni
      integer,dimension(  leni),intent(in) :: kxi
      real   ,dimension(3,leni),intent(in) :: qri
      integer,dimension(  leni),intent(in) :: kli

      integer,intent(in) :: lenj
      integer,dimension(  lenj),intent(in) :: kxj
      real   ,dimension(3,lenj),intent(in) :: qrj
      integer,dimension(  lenj),intent(in) :: klj

      integer, intent(in) :: nvecs
      real,dimension(lenj), intent(in)    :: x
      real,dimension(leni), intent(inout) :: y

! !REVISION HISTORY:
! 	23Feb01	- Tom Clune <clune@sgi.com>
!                 Initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rHH1cxpy'

	include "lvmax.h"

		! observation error class references
	include "ktmax.h"
	include "kxmax.h"
	include "kxtabl.h"

		! indirect vertical correlation matrix
	include "MX_voecH.h"
	include "voecHH.h"

		! indirect horizontal correlaiton matrix factors
	include "MX_hoecH.h"
	include "hoecHH.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	integer ij,i,j,l
	integer ki,kj,itau,ktau
	integer kx_i,kx_j
	integer jvtype,jhtype
	real tau,ctau,xtau,redwin_
	real qri_x,qri_y,qri_z
	real qrj_x,qrj_y,qrj_z
	real vercor,horcor
	real Ocoslimj,qxOtb1j,qxOtb2j,OObeg2j
	logical correlated

!-----------------------------------------------------------------------
	real HH, corrO

!=======================================================================
	if(.not.redwin_initialized())	&
		call die(myname_,'redwin not initialized')

	if(leni.le.0.or.lenj.le.0) return

	do j=1,lenj

	   kx_j = kxj(j)

			! types of the horizontal correlation functions
			! types of the vertical correlation functions

	  jhtype=i_hoecHc(kx_j)
	  jvtype=i_voecHc(kx_j)

	  correlated=jhtype.gt.0.and.jhtype.le.MX_hoecH .and.	&
		     jvtype.gt.0.and.jvtype.le.MX_voecH

	  if(correlated) then

		! a preliminary check if there is any (other than the
		! diagonals with the same correlation ID (pix).

	    correlated=.false.
	    i=1
	    do while(i.le.leni.and..not.correlated)
	      kx_i=kxi(i)
	      correlated=kx_i.eq.kx_j
	      i=i+1
	    end do
	  endif

	  if(correlated) then

	    qrj_x=qrj(1,j)
	    qrj_y=qrj(2,j)
	    qrj_z=qrj(3,j)


	    kj=klj(j)

	    Ocoslimj=Ocoslim(jhtype)
	    qxOtb1j=qxOtb1(jhtype)
	    qxOtb2j=qxOtb2(jhtype)
	    OObeg2j=OObeg2(jhtype)

	    do i=1,leni

			! cossep is the distance between i and j in
			! cos(great_circle_angular_distance)

	      kx_i=kxi(i)
	      if(kx_i.eq.kx_j) then

		 qri_x = qri(1,i)
		 qri_y = qri(2,i)
		 qri_z = qri(3,i)
		 tau = qri_x*qrj_x + qri_y*qrj_y + qri_z*qrj_z

		horcor=0.
	        if(tau.gt.Ocoslimj) then

		   ki = kli(i)

			! table indexing

		  ctau=1.-tau
#ifdef _LINEAR
		  xtau=qxOtb1j*ctau+1.
		  if(ctau.gt.OObeg2j)		&
			xtau=qxOtb2j*(ctau-OObeg2j)+nOOtb1+1.
		  itau=min(int(xtau),nHHotab-1)
		  xtau=xtau-itau
		  HH=hoecHH(ki,itau,jhtype)+hoecHH(kj,itau,jhtype)
		  HH=HH + xtau*			&
		    (hoecHH(ki,itau+1,jhtype)+	&
		     hoecHH(kj,itau+1,jhtype)-	HH)
#else
		  xtau=qxOtb1j*ctau+1.5
		  if(ctau.gt.OObeg2j)		&
			xtau=qxOtb2j*(ctau-OObeg2j)+nOOtb1+1.5

		  itau=int(xtau)	! in [1,nHHotab] is assumed

			! a non-separable correlation function form

		  HH=hoecHH(ki,itau,jhtype)+hoecHH(kj,itau,jhtype)
#endif
		  horcor = HH

		  ktau=int(qxWtb*ctau+.5)+1
		  redwin_=0.
		  if(ktau<=mxWtb) redwin_=redwin(ktau)

		  vercor=voecHH(ki,kj,jvtype)*redwin_

		  
		  corrO = horcor*vercor
		  y(i) = y(i) + corrO * x(j)
		endif	! tau.gt.Ocoslimj
	      endif	! kx_i.eq.kx_j
	    end do	! i=1,leni

	  else
	     ! do nothing
	  endif

	end do		! j=1,lenj

end subroutine rHH1cxpy
end module m_kt_CorO
