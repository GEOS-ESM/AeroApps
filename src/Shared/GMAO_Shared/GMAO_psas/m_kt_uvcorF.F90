!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_kt_uvcorF - correlation matrix generators for Fcst.Err.
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_kt_uvcorF
      implicit none
      private	! except

      public :: fQQcor1
      public :: fQQcorx

      ! These versions fuse the construction with of a QQ block with
      ! the BLAS mat-vec multiply.  This significantly reduces
      ! demands on memory and thereby improves performance.
      public :: sQQmcxpy
      public :: xQQmcxpy
      public :: rQQmcxpy
      public :: sQQ1cxpy
      public :: xQQ1cxpy
      public :: rQQ1cxpy

! !REVISION HISTORY:
! 	23Feb01	- Tom Clune <clune@sgi.com>
!                 Added new interfaces that avoid storing the block
!                 but rather apply elements to right/left vectors as
!                 individual elements are constructed.  Significant
!                 performance improvement.
!	29Aug00	- Jing Guo
!		- initial prototype/prolog/code
!		- combined fQQcor1.F etc. into this module
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_kt_uvcorF'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: fQQcor1 - packed symmetrix correlation matrix, Q-Q
!
! !DESCRIPTION:
!
!	fQQcor1 forms the upper triangle of a symmetric <q,q>
!	correlation matrix.  The returned matrix is in a packed form in
!	the order of ((i=1,j),j=1,len).  The correlation function is
!	separable and in the same form as in the OI of GEOS/DAS.1.2.
!
! !INTERFACE:

    subroutine fQQcor1(len,qr,kl,corrF)
      use m_redwin,only : redwin
      use m_redwin,only : qxWtb,mxWtb
      use m_redwin,only : redwin_initialized
      use m_die,only : die
      implicit none

		! the size of the matrix
      integer,intent(in) :: len

		! the attributes of rows and columns

      real   ,dimension(3,len),intent(in) :: qr
      integer,dimension(  len),intent(in) :: kl

		! a packed matrix

      real,dimension(len*(len+1)/2),intent(out) :: corrF

! !REVISION HISTORY:
! 	29Aug00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- unified argument style and modified to be a module
!		  procedure.
!
!	05Dec95 - J. Guo
!
!		1) separated forecast error correlation matrix from
!		   observation error correlation matrix;
!		2) used a new structure in hope that it will be
!		   optimized better; and
!		3) remove some redundant computations, expecially in
!		   wind related correlation functions.
!		4) getting ready for the new changes (see BUGS:) that
!		   may effect the results.
!
!  29Sep95  - Jing G.	- added #ifdef for dynamic dimension check
!  18Sep95  - Jing G.	- Created multiple copies of statistical tables
!  27Mar95  - Jing G.	- Modified for text file based data tables.
!  19Jan95  - Jing G.	- Added wobs table to replace pindx2().  One
!			  could use rlevs for the same purpose to
!			  reduce the overhead, since rlevs has no real
!			  purpose in this subroutine.
!  25jul94  Meta S.  original subroutine based on HHcor1
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::fQQcor1'

!-----------------------------------------------------------------------
		! Indirect vertical correlation matrix
	include	"lvmax.h"
	include "vfecQQ.h"

!-----------------------------------------------------------------------
		! indirect horizontal correlation matrix factors
	include "hfecQQ.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	integer ij,i,j,l
	integer ki,kj,itau,ktau
	real tau,ctau,xtau,redwin_
	real qri_x,qri_y,qri_z
	real qrj_x,qrj_y,qrj_z
	real vercor,horcor

!-----------------------------------------------------------------------
	real QQ
!=======================================================================
	if(.not.redwin_initialized())	&
		call die(myname_,'redwin not initialized')

	if(len.le.0) return

	corrF(1)=1.

	do j=2,len

	  ij=j*(j-1)/2

	  qrj_x=qr(1,j)
	  qrj_y=qr(2,j)
	  qrj_z=qr(3,j)

	  kj=kl(j)

	  do i=1,j-1	! i<j

		! l is an array element index of a packed upper
		! triangular matrix

	    l=ij+i

		! cossep is the polarity index of the two (i
		! and j) location vectors.  It is defined as
		! (q_i,q_j) = cos(angular_separation_of_i&j)

	    qri_x = qr(1,i)
	    qri_y = qr(2,i)
	    qri_z = qr(3,i)
	    tau = qri_x*qrj_x + qri_y*qrj_y + qri_z*qrj_z


	    corrF(l)=0.
	    if(tau.gt.Qcoslim) then

	      ki=kl(i)

			! table indexing

	      ctau=1.-tau
#ifdef _LINEAR
	      xtau=qxQtb1*ctau+1.
	      if(ctau.gt.QQbeg2) xtau=qxQtb2*(ctau-QQbeg2)+nQQtb1+1.
	      itau=min(int(xtau),nQQtab-1)
	      xtau=xtau-itau
	      QQ=hfecQQ(ki,itau)+hfecQQ(kj,itau)
	      QQ=QQ + xtau*		&
     &		(hfecQQ(ki,itau+1)+hfecQQ(kj,itau+1) - QQ)
#else
	      xtau=qxQtb1*ctau+1.5
	      if(ctau.gt.QQbeg2) xtau=qxQtb2*(ctau-QQbeg2)+nQQtb1+1.5

	      itau=int(xtau)	! in [1,nQQtab] is assumed

			! a non-separable correlation function form

	      QQ=hfecQQ(ki,itau)+hfecQQ(kj,itau)
#endif
	      horcor =	QQ

	      ktau=int(qxWtb*ctau+.5)+1
	      redwin_=0.
	      if(ktau<=mxWtb) redwin_=redwin(ktau)

	      vercor=vfecQQ(ki,kj)*redwin_

	      corrF(l)=horcor*vercor
	    endif		! tau .gt. Qcoslim
	  end do		! i=1,j-1

	  l=ij+j		! when i=j
	  corrF(l)=1.
	end do			! j=2,len

end subroutine fQQcor1

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sQQcxpy - symmetrix correlation matrix, Q-Q
!
! !DESCRIPTION:
!
!	sQQcxpy applies a symmetric <q,q> correlation matrix to right
!       and left hand vectors.The correlation function is separable and
!       in the same form as in the OI of GEOS/DAS.1.2.  Symmery is used
!       to reduce computation by roughly 50%.
!
! !INTERFACE:

    subroutine sQQmcxpy(len,qr,kl,nvecs,x,y)
      use m_redwin,only : redwin
      use m_redwin,only : qxWtb,mxWtb
      use m_redwin,only : redwin_initialized
      use m_die,only : die
      implicit none

		! the size of the matrix
      integer,intent(in) :: len

		! the attributes of rows and columns

      real   ,dimension(3,len),intent(in) :: qr
      integer,dimension(  len),intent(in) :: kl

		! a packed matrix

      Integer, Intent(In) :: nvecs
      real,dimension(nvecs,len), intent(in)    :: x
      real,dimension(nvecs,len), intent(inout) :: y

! !REVISION HISTORY:
! 	23Feb01	- Tom Clune <clune@sgi.com>
!                 Initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sQQcxpy'

!-----------------------------------------------------------------------
		! Indirect vertical correlation matrix
	include	"lvmax.h"
	include "vfecQQ.h"

!-----------------------------------------------------------------------
		! indirect horizontal correlation matrix factors
	include "hfecQQ.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	integer ij,i,j,l
	integer ki,kj,itau,ktau
	real tau,ctau,xtau,redwin_
	real qri_x,qri_y,qri_z
	real qrj_x,qrj_y,qrj_z
	real vercor,horcor

!-----------------------------------------------------------------------
	real QQ, corrF
!=======================================================================
	if(.not.redwin_initialized())	&
		call die(myname_,'redwin not initialized')

	if(len.le.0) return

	y(:,1) = y(:,1) + x(:,1)

	do j=2,len

	  ij=j*(j-1)/2

	  qrj_x=qr(1,j)
	  qrj_y=qr(2,j)
	  qrj_z=qr(3,j)

	  kj=kl(j)

	  do i=1,j-1	! i<j

		! l is an array element index of a packed upper
		! triangular matrix

	    l=ij+i

		! cossep is the polarity index of the two (i
		! and j) location vectors.  It is defined as
		! (q_i,q_j) = cos(angular_separation_of_i&j)

	    qri_x = qr(1,i)
	    qri_y = qr(2,i)
	    qri_z = qr(3,i)
	    tau = qri_x*qrj_x + qri_y*qrj_y + qri_z*qrj_z


	    if(tau.gt.Qcoslim) then

	      ki=kl(i)

			! table indexing

	      ctau=1.-tau
#ifdef _LINEAR
	      xtau=qxQtb1*ctau+1.
	      if(ctau.gt.QQbeg2) xtau=qxQtb2*(ctau-QQbeg2)+nQQtb1+1.
	      itau=min(int(xtau),nQQtab-1)
	      xtau=xtau-itau
	      QQ=hfecQQ(ki,itau)+hfecQQ(kj,itau)
	      QQ=QQ + xtau*		&
     &		(hfecQQ(ki,itau+1)+hfecQQ(kj,itau+1) - QQ)
#else
	      xtau=qxQtb1*ctau+1.5
	      if(ctau.gt.QQbeg2) xtau=qxQtb2*(ctau-QQbeg2)+nQQtb1+1.5

	      itau=int(xtau)	! in [1,nQQtab] is assumed

			! a non-separable correlation function form

	      QQ=hfecQQ(ki,itau)+hfecQQ(kj,itau)
#endif
	      horcor =	QQ

	      ktau=int(qxWtb*ctau+.5)+1
	      redwin_=0.
	      if(ktau<=mxWtb) redwin_=redwin(ktau)

	      vercor=vfecQQ(ki,kj)*redwin_

	      corrF=horcor*vercor
	      y(:,i) = y(:,i) + corrF * x(:,j)
	      y(:,j) = y(:,j) + corrF * x(:,i)
	    endif		! tau .gt. Qcoslim
	  end do		! i=1,j-1

	  l=ij+j		! when i=j
	  y(:,j) = y(:,j) + x(:,j)
	end do			! j=2,len

      end subroutine sQQmcxpy

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: fQQcorx - Fcst. err. correlation matrix, Q-Q
!
! !DESCRIPTION:
!
!	fQQcorx() forms a forecast error correlation matrix for
!	a univariate moisture variable.
!
! !INTERFACE:

    subroutine fQQcorx(leni,qri,kli, lenj,qrj,klj, corrF)
      use m_redwin,only : redwin
      use m_redwin,only : qxWtb,mxWtb
      use m_redwin,only : redwin_initialized
      use m_die,only : die
      implicit none

		! row size of the matrix

      integer,intent(in) :: leni

		! attributes of the row elements

      real   ,dimension(3,leni),intent(in) :: qri
      integer,dimension(  leni),intent(in) :: kli

		! column size of the matrix

      integer,intent(in) :: lenj

		! attributes of the column elements

      real   ,dimension(3,lenj),intent(in) :: qrj
      integer,dimension(  lenj),intent(in) :: klj


		! a correlation matrix

      real,dimension(leni,lenj),intent(out) :: corrF

! !REVISION HISTORY:
! 	29Aug00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- unified argument style and modified to be a module
!		  procedure.
!
! 	09Nov95 - Jing G.	- rewrite the code and the prolog.
!
!		1) separated forecast error correlation matrix from
!		   observation error correlation matrix;
!		2) used a new structure in hope that it will be
!		   optimized better; and
!		3) remove some redundant computations, expecially in
!		   wind related correlation functions.
!		4) getting ready for the new changes (see BUGS:) that
!		   may effect the results.
!
!  29Sep95  - Jing G.	- added #ifdef for dynamic dimension check
!  18Sep95  - Jing G.	- Created multiple copies of statistical tables
!  19Jan95  - Jing G.	- Added wobs table to replace pindx2().  One
!			  could use rlevs for the same purpose to reduce
!			  the overhead, since rlevs has no real purpose
!			  in this subroutine.
!
!  07oct94  Meta S., A.daS. slight mod of obs correlation code
!  10may94  Meta S.  change observation error correlation tabulation
!                     to use cos(dist/rade), added prologue
!  28apr94  Meta S.  change forecast error correlation tabulation
!                     to use cos(dist/rade), added prologue
!  15dec93  Meta S.  add observation covariances, covariance tables
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::fQQcorx'

!-----------------------------------------------------------------------
!   ..Indirect vertical correlaton matrix

	include	"lvmax.h"
	include "vfecQQ.h"

!-----------------------------------------------------------------------
!   ..Indirect horizontal correlation matrix factors

	include "hfecQQ.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	real	qri_x	! x of qri (location vector of rows)
	real	qri_y	! y of qri
	real	qri_z	! z of qri
	real	qrj_x	! x of qri (location vector of columns)
	real	qrj_y	! y of qri
	real	qrj_z	! z of qri
	integer ij,i,j,l
	integer ki,kj,itau,ktau
	real tau,ctau,xtau,redwin_

	real vercor,horcor

!-----------------------------------------------------------------------
	real QQ
!=======================================================================
	if(.not.redwin_initialized())	&
		call die(myname_,'redwin not initialized')

	if(leni.le.0.or.lenj.le.0) return

!   ..Loop over columns

	do j=1,lenj

	  qrj_x=qrj(1,j)
	  qrj_y=qrj(2,j)
	  qrj_z=qrj(3,j)

	  kj=klj(j)

			! loop over rows

	  do i=1,leni

			! cossep is the polarity index of the two (i
			! and j) location vectors.  It is defined as
			! (q_i,q_j) = cos(angular_separation_of_i&j)

	     qri_x = qri(1,i)
	     qri_y = qri(2,i)
	     qri_z = qri(3,i)

	     tau=qri_x*qrj_x + qri_y*qrj_y + qri_z*qrj_z

	    corrF(i,j)=0.
	    if(tau.gt.Qcoslim) then

	      ki=kli(i)

			! table indexing

	      ctau=1.-tau
#ifdef _LINEAR
	      xtau=qxQtb1*ctau+1.
	      if(ctau.gt.QQbeg2) xtau=qxQtb2*(ctau-QQbeg2)+nQQtb1+1.
	      itau=min(int(xtau),nQQtab-1)
	      xtau=xtau-itau
	      QQ=hfecQQ(ki,itau)+hfecQQ(kj,itau)
	      QQ=QQ + xtau*		&
     &		(hfecQQ(ki,itau+1)+hfecQQ(kj,itau+1) - QQ)
#else
	      xtau=qxQtb1*ctau+1.5
	      if(ctau.gt.QQbeg2) xtau=qxQtb2*(ctau-QQbeg2)+nQQtb1+1.5

	      itau=int(xtau)	! in [1,nQQtab] is assumed

			! a non-separable correlation function form

	      QQ=hfecQQ(ki,itau)+hfecQQ(kj,itau)
#endif
	      horcor =	QQ

	      ktau=int(qxWtb*ctau+.5)+1
	      redwin_=0.
	      if(ktau<=mxWtb) redwin_=redwin(ktau)

	      vercor=vfecQQ(ki,kj)*redwin_

	      corrF(i,j)=horcor*vercor
	    endif		! if (tau .gt. Qcorlim)

	  end do		! i=1,leni
	end do			! j=1,lenj
end subroutine fQQcorx

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rQQcxpy - Fcst. err. correlation matrix, Q-Q
!
! !DESCRIPTION:
!
!	rQQcxpy() applies a forecast error correlation matrix for
!	a univariate moisture variable to a right vector.
!
! !INTERFACE:

    subroutine rQQmcxpy(leni,qri,kli, lenj,qrj,klj,nvecs,x,y)
      use m_redwin,only : redwin
      use m_redwin,only : qxWtb,mxWtb
      use m_redwin,only : redwin_initialized
      use m_die,only : die
      implicit none

		! row size of the matrix
      integer,intent(in) :: leni
		! attributes of the row elements
      real   ,dimension(3,leni),intent(in) :: qri
      integer,dimension(  leni),intent(in) :: kli
		! column size of the matrix
      integer,intent(in) :: lenj
		! attributes of the column elements
      real   ,dimension(3,lenj),intent(in) :: qrj
      integer,dimension(  lenj),intent(in) :: klj
		! a correlation matrix
      integer, intent(in)                 :: nvecs
      real,dimension(nvecs,lenj), intent(in)    :: x
      real,dimension(nvecs,leni), intent(inout) :: y

! !REVISION HISTORY:
! 	23Feb01	- Tom Clune <clune@sgi.com>
!                 Initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::fQQcorx_Cx'

!-----------------------------------------------------------------------
!   ..Indirect vertical correlaton matrix

	include	"lvmax.h"
	include "vfecQQ.h"

!-----------------------------------------------------------------------
!   ..Indirect horizontal correlation matrix factors

	include "hfecQQ.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	real	qri_x	! x of qri (location vector of rows)
	real	qri_y	! y of qri
	real	qri_z	! z of qri
	real	qrj_x	! x of qri (location vector of columns)
	real	qrj_y	! y of qri
	real	qrj_z	! z of qri
	integer ij,i,j,l
	integer ki,kj,itau,ktau
	real tau,ctau,xtau,redwin_

	real vercor,horcor
	real corrf

!-----------------------------------------------------------------------
	real QQ
!=======================================================================
	if(.not.redwin_initialized())	&
		call die(myname_,'redwin not initialized')

	if(leni.le.0.or.lenj.le.0) return

!   ..Loop over columns

	do j=1,lenj

	  qrj_x=qrj(1,j)
	  qrj_y=qrj(2,j)
	  qrj_z=qrj(3,j)

	  kj=klj(j)

			! loop over rows

	  do i=1,leni

			! cossep is the polarity index of the two (i
			! and j) location vectors.  It is defined as
			! (q_i,q_j) = cos(angular_separation_of_i&j)

	     qri_x = qri(1,i)
	     qri_y = qri(2,i)
	     qri_z = qri(3,i)

	     tau=qri_x*qrj_x + qri_y*qrj_y + qri_z*qrj_z

	    if(tau.gt.Qcoslim) then

	      ki=kli(i)

			! table indexing

	      ctau=1.-tau
#ifdef _LINEAR
	      xtau=qxQtb1*ctau+1.
	      if(ctau.gt.QQbeg2) xtau=qxQtb2*(ctau-QQbeg2)+nQQtb1+1.
	      itau=min(int(xtau),nQQtab-1)
	      xtau=xtau-itau
	      QQ=hfecQQ(ki,itau)+hfecQQ(kj,itau)
	      QQ=QQ + xtau*		&
     &		(hfecQQ(ki,itau+1)+hfecQQ(kj,itau+1) - QQ)
#else
	      xtau=qxQtb1*ctau+1.5
	      if(ctau.gt.QQbeg2) xtau=qxQtb2*(ctau-QQbeg2)+nQQtb1+1.5

	      itau=int(xtau)	! in [1,nQQtab] is assumed

			! a non-separable correlation function form

	      QQ=hfecQQ(ki,itau)+hfecQQ(kj,itau)
#endif
	      horcor =	QQ

	      ktau=int(qxWtb*ctau+.5)+1
	      redwin_=0.
	      if(ktau<=mxWtb) redwin_=redwin(ktau)

	      vercor=vfecQQ(ki,kj)*redwin_

	      corrF=horcor*vercor
	      y(:,i) = y(:,i) + corrF * x(:,j)

	    endif		! if (tau .gt. Qcorlim)

	  end do		! i=1,leni
	end do			! j=1,lenj
end subroutine rQQmcxpy

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: xQQcxpy - Fcst. err. correlation matrix, Q-Q
!
! !DESCRIPTION:
!
!	xQQcxpy() applies a forecast error correlation matrix for
!	a univariate moisture variable to both right and left hand vectors
!
! !INTERFACE:

    subroutine xQQmcxpy(leni,qri,kli, lenj,qrj,klj,nvecs,xi,yj,xj,yi)
      use m_redwin,only : redwin
      use m_redwin,only : qxWtb,mxWtb
      use m_redwin,only : redwin_initialized
      use m_die,only : die
      implicit none

		! row size of the matrix
      integer,intent(in) :: leni
		! attributes of the row elements
      real   ,dimension(3,leni),intent(in) :: qri
      integer,dimension(  leni),intent(in) :: kli
		! column size of the matrix
      integer,intent(in) :: lenj
		! attributes of the column elements
      real   ,dimension(3,lenj),intent(in) :: qrj
      integer,dimension(  lenj),intent(in) :: klj
		! a correlation matrix
      integer, intent(in)                 :: nvecs
      real,dimension(nvecs,leni), intent(in)    :: xi
      real,dimension(nvecs,lenj), intent(inout) :: yj
      real,dimension(nvecs,lenj), intent(in)    :: xj
      real,dimension(nvecs,leni), intent(inout) :: yi

! !REVISION HISTORY:
! 	23Feb01	- Tom Clune <clune@sgi.com>
!                 Initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::fQQcorx_xCx'

!-----------------------------------------------------------------------
!   ..Indirect vertical correlaton matrix

	include	"lvmax.h"
	include "vfecQQ.h"

!-----------------------------------------------------------------------
!   ..Indirect horizontal correlation matrix factors

	include "hfecQQ.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	real	qri_x	! x of qri (location vector of rows)
	real	qri_y	! y of qri
	real	qri_z	! z of qri
	real	qrj_x	! x of qri (location vector of columns)
	real	qrj_y	! y of qri
	real	qrj_z	! z of qri
	integer ij,i,j,l
	integer ki,kj,itau,ktau
	real tau,ctau,xtau,redwin_

	real vercor,horcor
	real corrf

!-----------------------------------------------------------------------
	real QQ
!=======================================================================
	if(.not.redwin_initialized())	&
		call die(myname_,'redwin not initialized')

	if(leni.le.0.or.lenj.le.0) return

!   ..Loop over columns

	do j=1,lenj

	  qrj_x=qrj(1,j)
	  qrj_y=qrj(2,j)
	  qrj_z=qrj(3,j)

	  kj=klj(j)

			! loop over rows

	  do i=1,leni

			! cossep is the polarity index of the two (i
			! and j) location vectors.  It is defined as
			! (q_i,q_j) = cos(angular_separation_of_i&j)

	     qri_x = qri(1,i)
	     qri_y = qri(2,i)
	     qri_z = qri(3,i)

	     tau=qri_x*qrj_x + qri_y*qrj_y + qri_z*qrj_z

	    if(tau.gt.Qcoslim) then

	      ki=kli(i)

			! table indexing

	      ctau=1.-tau
#ifdef _LINEAR
	      xtau=qxQtb1*ctau+1.
	      if(ctau.gt.QQbeg2) xtau=qxQtb2*(ctau-QQbeg2)+nQQtb1+1.
	      itau=min(int(xtau),nQQtab-1)
	      xtau=xtau-itau
	      QQ=hfecQQ(ki,itau)+hfecQQ(kj,itau)
	      QQ=QQ + xtau*		&
     &		(hfecQQ(ki,itau+1)+hfecQQ(kj,itau+1) - QQ)
#else
	      xtau=qxQtb1*ctau+1.5
	      if(ctau.gt.QQbeg2) xtau=qxQtb2*(ctau-QQbeg2)+nQQtb1+1.5

	      itau=int(xtau)	! in [1,nQQtab] is assumed

			! a non-separable correlation function form

	      QQ=hfecQQ(ki,itau)+hfecQQ(kj,itau)
#endif
	      horcor =	QQ

	      ktau=int(qxWtb*ctau+.5)+1
	      redwin_=0.
	      if(ktau<=mxWtb) redwin_=redwin(ktau)

	      vercor=vfecQQ(ki,kj)*redwin_

	      corrF=horcor*vercor
	      yj(:,j) = yj(:,j) + corrF * xi(:,i)
	      yi(:,i) = yi(:,i) + corrF * xj(:,j)

	    endif		! if (tau .gt. Qcorlim)

	  end do		! i=1,leni
	end do			! j=1,lenj
      end subroutine xQQmcxpy

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sQQcxpy - symmetrix correlation matrix, Q-Q
!
! !DESCRIPTION:
!
!	sQQcxpy applies a symmetric <q,q> correlation matrix to right
!       and left hand vectors.The correlation function is separable and
!       in the same form as in the OI of GEOS/DAS.1.2.  Symmery is used
!       to reduce computation by roughly 50%.
!
! !INTERFACE:

    subroutine sQQ1cxpy(len,qr,kl,nvecs,x,y)
      use m_redwin,only : redwin
      use m_redwin,only : qxWtb,mxWtb
      use m_redwin,only : redwin_initialized
      use m_die,only : die
      implicit none

		! the size of the matrix
      integer,intent(in) :: len

		! the attributes of rows and columns
      real   ,dimension(3,len),intent(in) :: qr
      integer,dimension(  len),intent(in) :: kl

                ! operands
      Integer, Intent(In) :: nvecs
      real,dimension(len), intent(in)    :: x
      real,dimension(len), intent(inout) :: y

! !REVISION HISTORY:
! 	23Feb01	- Tom Clune <clune@sgi.com>
!                 Initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sQQcxpy'

!-----------------------------------------------------------------------
		! Indirect vertical correlation matrix
	include	"lvmax.h"
	include "vfecQQ.h"

!-----------------------------------------------------------------------
		! indirect horizontal correlation matrix factors
	include "hfecQQ.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	integer ij,i,j,l
	integer ki,kj,itau,ktau
	real tau,ctau,xtau,redwin_
	real qri_x,qri_y,qri_z
	real qrj_x,qrj_y,qrj_z
	real vercor,horcor

!-----------------------------------------------------------------------
	real QQ, corrF
!=======================================================================
	if(.not.redwin_initialized())	&
		call die(myname_,'redwin not initialized')

	if(len.le.0) return

	y(1) = y(1) + x(1)

	do j=2,len

	  ij=j*(j-1)/2

	  qrj_x=qr(1,j)
	  qrj_y=qr(2,j)
	  qrj_z=qr(3,j)

	  kj=kl(j)

	  do i=1,j-1	! i<j

		! l is an array element index of a packed upper
		! triangular matrix

	    l=ij+i

		! cossep is the polarity index of the two (i
		! and j) location vectors.  It is defined as
		! (q_i,q_j) = cos(angular_separation_of_i&j)

	    qri_x = qr(1,i)
	    qri_y = qr(2,i)
	    qri_z = qr(3,i)
	    tau = qri_x*qrj_x + qri_y*qrj_y + qri_z*qrj_z


	    if(tau.gt.Qcoslim) then

	      ki=kl(i)

			! table indexing

	      ctau=1.-tau
#ifdef _LINEAR
	      xtau=qxQtb1*ctau+1.
	      if(ctau.gt.QQbeg2) xtau=qxQtb2*(ctau-QQbeg2)+nQQtb1+1.
	      itau=min(int(xtau),nQQtab-1)
	      xtau=xtau-itau
	      QQ=hfecQQ(ki,itau)+hfecQQ(kj,itau)
	      QQ=QQ + xtau*		&
     &		(hfecQQ(ki,itau+1)+hfecQQ(kj,itau+1) - QQ)
#else
	      xtau=qxQtb1*ctau+1.5
	      if(ctau.gt.QQbeg2) xtau=qxQtb2*(ctau-QQbeg2)+nQQtb1+1.5

	      itau=int(xtau)	! in [1,nQQtab] is assumed

			! a non-separable correlation function form

	      QQ=hfecQQ(ki,itau)+hfecQQ(kj,itau)
#endif
	      horcor =	QQ

	      ktau=int(qxWtb*ctau+.5)+1
	      redwin_=0.
	      if(ktau<=mxWtb) redwin_=redwin(ktau)

	      vercor=vfecQQ(ki,kj)*redwin_

	      corrF=horcor*vercor
	      y(i) = y(i) + corrF * x(j)
	      y(j) = y(j) + corrF * x(i)
	    endif		! tau .gt. Qcoslim
	  end do		! i=1,j-1

	  l=ij+j		! when i=j
	  y(j) = y(j) + x(j)
	end do			! j=2,len

      end subroutine sQQ1cxpy

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rQQcxpy - Fcst. err. correlation matrix, Q-Q
!
! !DESCRIPTION:
!
!	rQQcxpy() applies a forecast error correlation matrix for
!	a univariate moisture variable to a right vector.
!
! !INTERFACE:

    subroutine rQQ1cxpy(leni,qri,kli, lenj,qrj,klj,nvecs,x,y)
      use m_redwin,only : redwin
      use m_redwin,only : qxWtb,mxWtb
      use m_redwin,only : redwin_initialized
      use m_die,only : die
      implicit none

		! row size of the matrix
      integer,intent(in) :: leni
		! attributes of the row elements
      real   ,dimension(3,leni),intent(in) :: qri
      integer,dimension(  leni),intent(in) :: kli
		! column size of the matrix
      integer,intent(in) :: lenj
		! attributes of the column elements
      real   ,dimension(3,lenj),intent(in) :: qrj
      integer,dimension(  lenj),intent(in) :: klj
                ! operands
      Integer, Intent(In) :: nvecs
      real,dimension(lenj), intent(in)    :: x
      real,dimension(leni), intent(inout) :: y

! !REVISION HISTORY:
! 	23Feb01	- Tom Clune <clune@sgi.com>
!                 Initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::fQQcorx_Cx'

!-----------------------------------------------------------------------
!   ..Indirect vertical correlaton matrix

	include	"lvmax.h"
	include "vfecQQ.h"

!-----------------------------------------------------------------------
!   ..Indirect horizontal correlation matrix factors

	include "hfecQQ.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	real	qri_x	! x of qri (location vector of rows)
	real	qri_y	! y of qri
	real	qri_z	! z of qri
	real	qrj_x	! x of qri (location vector of columns)
	real	qrj_y	! y of qri
	real	qrj_z	! z of qri
	integer ij,i,j,l
	integer ki,kj,itau,ktau
	real tau,ctau,xtau,redwin_

	real vercor,horcor
	real corrf

!-----------------------------------------------------------------------
	real QQ
!=======================================================================
	if(.not.redwin_initialized())	&
		call die(myname_,'redwin not initialized')

	if(leni.le.0.or.lenj.le.0) return

!   ..Loop over columns

	do j=1,lenj

	  qrj_x=qrj(1,j)
	  qrj_y=qrj(2,j)
	  qrj_z=qrj(3,j)

	  kj=klj(j)

			! loop over rows

	  do i=1,leni

			! cossep is the polarity index of the two (i
			! and j) location vectors.  It is defined as
			! (q_i,q_j) = cos(angular_separation_of_i&j)

	     qri_x = qri(1,i)
	     qri_y = qri(2,i)
	     qri_z = qri(3,i)

	     tau=qri_x*qrj_x + qri_y*qrj_y + qri_z*qrj_z

	    if(tau.gt.Qcoslim) then

	      ki=kli(i)

			! table indexing

	      ctau=1.-tau
#ifdef _LINEAR
	      xtau=qxQtb1*ctau+1.
	      if(ctau.gt.QQbeg2) xtau=qxQtb2*(ctau-QQbeg2)+nQQtb1+1.
	      itau=min(int(xtau),nQQtab-1)
	      xtau=xtau-itau
	      QQ=hfecQQ(ki,itau)+hfecQQ(kj,itau)
	      QQ=QQ + xtau*		&
     &		(hfecQQ(ki,itau+1)+hfecQQ(kj,itau+1) - QQ)
#else
	      xtau=qxQtb1*ctau+1.5
	      if(ctau.gt.QQbeg2) xtau=qxQtb2*(ctau-QQbeg2)+nQQtb1+1.5

	      itau=int(xtau)	! in [1,nQQtab] is assumed

			! a non-separable correlation function form

	      QQ=hfecQQ(ki,itau)+hfecQQ(kj,itau)
#endif
	      horcor =	QQ

	      ktau=int(qxWtb*ctau+.5)+1
	      redwin_=0.
	      if(ktau<=mxWtb) redwin_=redwin(ktau)

	      vercor=vfecQQ(ki,kj)*redwin_

	      corrF=horcor*vercor
	      y(i) = y(i) + corrF * x(j)

	    endif		! if (tau .gt. Qcorlim)

	  end do		! i=1,leni
	end do			! j=1,lenj
      end subroutine rQQ1cxpy

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: xQQcxpy - Fcst. err. correlation matrix, Q-Q
!
! !DESCRIPTION:
!
!	xQQcxpy() applies a forecast error correlation matrix for
!	a univariate moisture variable to both right and left hand vectors
!
! !INTERFACE:

    subroutine xQQ1cxpy(leni,qri,kli, lenj,qrj,klj,nvecs,xi,yj,xj,yi)
      use m_redwin,only : redwin
      use m_redwin,only : qxWtb,mxWtb
      use m_redwin,only : redwin_initialized
      use m_die,only : die
      implicit none

		! row size of the matrix
      integer,intent(in) :: leni
		! attributes of the row elements
      real   ,dimension(3,leni),intent(in) :: qri
      integer,dimension(  leni),intent(in) :: kli
		! column size of the matrix
      integer,intent(in) :: lenj
		! attributes of the column elements
      real   ,dimension(3,lenj),intent(in) :: qrj
      integer,dimension(  lenj),intent(in) :: klj
                ! operands
      Integer, Intent(In) :: nvecs
      real,dimension(leni), intent(in)    :: xi
      real,dimension(lenj), intent(inout) :: yj
      real,dimension(lenj), intent(in)    :: xj
      real,dimension(leni), intent(inout) :: yi

! !REVISION HISTORY:
! 	23Feb01	- Tom Clune <clune@sgi.com>
!                 Initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::fQQcorx_xCx'

!-----------------------------------------------------------------------
!   ..Indirect vertical correlaton matrix

	include	"lvmax.h"
	include "vfecQQ.h"

!-----------------------------------------------------------------------
!   ..Indirect horizontal correlation matrix factors

	include "hfecQQ.h"

!-----------------------------------------------------------------------
!	..Runing indices and temporary variables etc.

	real	qri_x	! x of qri (location vector of rows)
	real	qri_y	! y of qri
	real	qri_z	! z of qri
	real	qrj_x	! x of qri (location vector of columns)
	real	qrj_y	! y of qri
	real	qrj_z	! z of qri
	integer ij,i,j,l
	integer ki,kj,itau,ktau
	real tau,ctau,xtau,redwin_

	real vercor,horcor
	real corrf

!-----------------------------------------------------------------------
	real QQ
!=======================================================================
	if(.not.redwin_initialized())	&
		call die(myname_,'redwin not initialized')

	if(leni.le.0.or.lenj.le.0) return

!   ..Loop over columns

	do j=1,lenj

	  qrj_x=qrj(1,j)
	  qrj_y=qrj(2,j)
	  qrj_z=qrj(3,j)

	  kj=klj(j)

			! loop over rows

	  do i=1,leni

			! cossep is the polarity index of the two (i
			! and j) location vectors.  It is defined as
			! (q_i,q_j) = cos(angular_separation_of_i&j)

	     qri_x = qri(1,i)
	     qri_y = qri(2,i)
	     qri_z = qri(3,i)

	     tau=qri_x*qrj_x + qri_y*qrj_y + qri_z*qrj_z

	    if(tau.gt.Qcoslim) then

	      ki=kli(i)

			! table indexing

	      ctau=1.-tau
#ifdef _LINEAR
	      xtau=qxQtb1*ctau+1.
	      if(ctau.gt.QQbeg2) xtau=qxQtb2*(ctau-QQbeg2)+nQQtb1+1.
	      itau=min(int(xtau),nQQtab-1)
	      xtau=xtau-itau
	      QQ=hfecQQ(ki,itau)+hfecQQ(kj,itau)
	      QQ=QQ + xtau*		&
     &		(hfecQQ(ki,itau+1)+hfecQQ(kj,itau+1) - QQ)
#else
	      xtau=qxQtb1*ctau+1.5
	      if(ctau.gt.QQbeg2) xtau=qxQtb2*(ctau-QQbeg2)+nQQtb1+1.5

	      itau=int(xtau)	! in [1,nQQtab] is assumed

			! a non-separable correlation function form

	      QQ=hfecQQ(ki,itau)+hfecQQ(kj,itau)
#endif
	      horcor =	QQ

	      ktau=int(qxWtb*ctau+.5)+1
	      redwin_=0.
	      if(ktau<=mxWtb) redwin_=redwin(ktau)

	      vercor=vfecQQ(ki,kj)*redwin_

	      corrF=horcor*vercor
	      yj(j) = yj(j) + corrF * xi(i)
	      yi(i) = yi(i) + corrF * xj(j)

	    endif		! if (tau .gt. Qcorlim)

	  end do		! i=1,leni
	end do			! j=1,lenj
      end subroutine xQQ1cxpy

end module m_kt_uvcorF
