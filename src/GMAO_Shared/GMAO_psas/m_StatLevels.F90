!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_StatLevels - Levels for statistics tables
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_StatLevels
      implicit none
      private	! except

      public :: nStatLevel		! The class data structure
      public :: StatLevels		! The class data structure
      public :: listStat

      interface listStat; module procedure kxstat_; end interface

! !REVISION HISTORY:
! 	06Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_StatLevels'

!  integer,parameter :: nStatLevel=25
!  real,dimension(nStatLevel),parameter :: StatLevels =	&
!	(/ 1000.,925.,850.,700.,600.,500.,400.,300.,250.,200.,	&
!	    150.,100., 70., 50., 30., 20., 10.,  7.,  5.,  2.,	&
!	      1.,  .7,  .5,  .2,  .1	/)

  integer,parameter :: nStatLevel=18
  real,dimension(nStatLevel),parameter :: StatLevels =	&
	(/ 1000.,850.,700.,500.,400.,300.,250.,200.,150.,100.,	&
	     70., 50., 30., 10.,  5.,  2.,  1.,  .4	/)

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: kxstat_ - list simple statistics
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine kxstat_(lu,kxs,kts,rlevs,vals,header)
      use m_die,only : die
      implicit none
      integer,intent(in) :: lu	! output unit
      integer,dimension(:),intent(in) :: kxs
      integer,dimension(:),intent(in) :: ktS
      real   ,dimension(:),intent(in) :: rlevs
      real   ,dimension(:),intent(in) :: vals
      character(len=*)    ,intent(in) :: header

! !REVISION HISTORY:
! 	02Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::kxstat_'
  logical,parameter :: NEAREST=.true.
  integer :: nval,nsum
  integer :: nkx,nkt
  integer :: ikx,ikt
  logical :: kxheader
  integer :: ier
  character(len=16) :: cinta,cintb
  integer,allocatable,dimension(:) :: ilev
  real   ,allocatable,dimension(:) :: wlev

  include "ktmax.h"
  include "kxmax.h"

  if(len_trim(header)==0) then
    write(lu,'(/2a)') myname_,'::'
  else
    write(lu,'(/a)') trim(header)
  endif

  nval=size(vals)

	allocate(ilev(nval),wlev(nval), stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  call slogtab(NEAREST,nStatLevel,StatLevels,nval,rlevs,ilev,wlev)

  nsum=0		! a check-sum for all
  do ikx=1,kxmax
    kxheader=.true.

    nkx=0	! a check-sum for the ikx
    do ikt=1,ktmax

      nkt=0
      select case(ikt)
      case (ktHH,ktUU,ktVV,ktQQ)

	call stat3d_(lu,ikx,ikt,vals,kxs,kts,rlevs,ilev,kxheader,nkt)

      case (ktslp,ktUs,ktVs)

	call stat2d_(lu,ikx,ikt,vals,kxs,kts,kxheader,nkt)

      end select

      nkx=nkx+nkt	! a check-sum for the ikx
    end do

    if(nkx>0) then
      write(cinta,'(i16)') ikx
      write(cintb,'(i16)') nkx
      write(lu,'(4a)') ' total (kx ',	&
	trim(adjustl(cinta)),') = ',trim(adjustl(cintb))
    endif

    nsum=nsum+nkx		! a check-sum for all
  end do

  write(cinta,'(i16)') nval
  write(cintb,'(i16)') nsum
  write(lu,'(4a)') '  total (all kx) = ',	&
	trim(adjustl(cinta)),' / ',		&
	trim(adjustl(cintb))

	deallocate(ilev,wlev,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine kxstat_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: stat3d_ - statistics for all levels if any
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine stat3d_(lu,ikx,ikt,vals,kxs,kts,rlevs,ilev,kxheader,nkt)
      implicit none
      integer,intent(in) :: lu
      integer,intent(in) :: ikx
      integer,intent(in) :: ikt
      real   ,dimension(:),intent(in) :: vals
      integer,dimension(:),intent(in) :: kxs
      integer,dimension(:),intent(in) :: kts
      real   ,dimension(:),intent(in) :: rlevs
      integer,dimension(:),intent(in) :: ilev
      logical,intent(inout) :: kxheader
      integer,intent(out) :: nkt

! !REVISION HISTORY:
! 	02Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::stat3d_'

  include "ktmax.h"
  include "kttabl.h"
  include "realvals.h"

  real   ,dimension(nStatLevel) :: avg
  real   ,dimension(nStatLevel) :: alv
  real   ,dimension(nStatLevel) :: var
  integer,dimension(nStatLevel) :: nav
  real   ,dimension(nStatLevel) :: amx
  real   ,dimension(nStatLevel) :: amn
  integer,dimension(nStatLevel) :: imx
  integer,dimension(nStatLevel) :: imn

  real    :: amag,xmag
  real    :: stdv,del
  integer :: ipw
  integer :: i,l
  character(len=4) :: cpw

  nkt=0

  avg(:)=0.
  alv(:)=0.
  nav(:)=0
  amx(:)=0.
  amn(:)=0.
  imx(:)=0
  imn(:)=0
  amag  =0.

  do i=1,size(vals)
    if(kxs(i)/=ikx .or. kts(i)/=ikt) cycle

    nkt=nkt+1
    l=ilev(i)

    avg(l)=avg(l)+vals(i)
    alv(l)=alv(l)+rlevs(i)
    nav(l)=nav(l)+1

    if(nav(l)==1) then
      amx(l)=vals(i)
      imx(l)=i
      amn(l)=vals(i)
      imn(l)=i
      amag  =abs(vals(i))
    endif

    if(vals(i) > amx(l)) then
      amx(l)=vals(i)
      imx(l)=i
    endif

    if(vals(i) < amn(l)) then
      amn(l)=vals(i)
      imn(l)=i
    endif

    if(abs(vals(i)) > amag) amag=abs(vals(i))
  end do

  if(nkt==0) return

  do l=1,nStatLevel
    if(nav(l)>0) then
      avg(l)=avg(l)/nav(l)
      alv(l)=alv(l)/nav(l)
    endif
  end do

  ipw=0
  if(amag >RSMLVAL) then
    xmag=log(amag)/log(10.)
    ipw=floor(xmag)
  endif
  amag=10.**ipw

  write(cpw,'(a,i2.2)') 'e+',ipw
  if(ipw<0) write(cpw,'(a,i2.2)') 'e-',ipw

  var(:)=0.
  do i=1,size(vals)
    if(kxs(i)/=ikx .or. kts(i)/=ikt) cycle

    l=ilev(i)

    del=vals(i)-avg(l)
    var(l)=var(l)+del*del
  end do

  do l=1,nStatLevel
    if(nav(l)>1) var(l)=var(l)/(nav(l)-1)
  end do

  if(kxheader) then
    call wtheader_(lu,ikx)
    kxheader=.false.
  endif

  do l=1,nStatLevel
    if(nav(l)<=0) cycle

    stdv=sqrt(var(l))

    write(lu,	&
      &	'(a8,x,a4,x,i4,x,f6.1,x,f6.3,x,f5.3,2(x,f6.3,x,i5))')	&
      &	ktname(ikt),cpw,nav(l),alv(l),				&
      &	avg(l)/amag,stdv/amag,					&
      &	amn(l)/amag,imn(l),amx(l)/amag,imx(l)
  end do

end subroutine stat3d_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: stat2d_ - statistics for all levels if any
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine stat2d_(lu,ikx,ikt,vals,kxs,kts,kxheader,nkt)
      implicit none
      integer,intent(in) :: lu
      integer,intent(in) :: ikx
      integer,intent(in) :: ikt
      real   ,dimension(:),intent(in) :: vals
      integer,dimension(:),intent(in) :: kxs
      integer,dimension(:),intent(in) :: kts
      logical,intent(in) :: kxheader
      integer,intent(out) :: nkt

! !REVISION HISTORY:
! 	02Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::stat2d_'

  include "ktmax.h"
  include "kttabl.h"
  include "realvals.h"

  real    :: avg
  real    :: var
  integer :: nav
  real    :: amx
  real    :: amn
  integer :: imx
  integer :: imn

  real    :: amag,xmag
  real    :: stdv,del
  integer :: ipw
  integer :: i
  character(len=4) :: cpw

  nkt=0

  avg =0.
  nav =0
  amx =0.
  amn =0.
  imx =0
  imn =0
  amag=0.

  do i=1,size(vals)
    if(kxs(i)/=ikx .or. kts(i)/=ikt) cycle

    nkt=nkt+1

    avg=avg+vals(i)
    nav=nav+1

    if(nav==1) then
      amx =vals(i)
      imx =i
      amn =vals(i)
      imn =i
      amag=abs(vals(i))
    endif

    if(vals(i) > amx) then
      amx=vals(i)
      imx=i
    endif

    if(vals(i) < amn) then
      amn=vals(i)
      imn=i
    endif

    if(abs(vals(i)) > amag) amag=abs(vals(i))
  end do

  if(nkt==0) return

  if(nav>0) avg=avg/nav

  ipw=0
  if(amag >RSMLVAL) then
    xmag=log(amag)/log(10.)
    ipw=floor(xmag)
  endif
  amag=10.**ipw

  write(cpw,'(a,i2.2)') 'E+',ipw
  if(ipw<0) write(cpw,'(a,i2.2)') 'E-',ipw

  var=0.
  do i=1,size(vals)
    if(kxs(i)/=ikx .or. kts(i)/=ikt) cycle

    del=vals(i)-avg
    var=var+del*del
  end do
  if(nav>1) var=var/(nav-1)

  if(kxheader) call wtheader_(lu,ikx)

  stdv=sqrt(var)

  write(lu,	&
      &	'(a8,x,a4,x,i4,x,a6,x,f6.3,x,f5.3,2(x,f6.3,x,i5))')	&
      &	ktname(ikt),cpw,nav,'  --- ',				&
      &	avg/amag,stdv/amag,					&
      &	amn/amag,imn,amx/amag,imx

end subroutine stat2d_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: wtheader_ - output a header
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine wtheader_(lu,ikx)
      implicit none
      integer,intent(in) :: lu
      integer,intent(in) :: ikx

! !REVISION HISTORY:
! 	02Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::wtheader_'
  include "ktmax.h"
  include "kxmax.h"
  include "kxtabl.h"

  write(lu,'(a,i4,2a)') 'kx ',ikx,', ',trim(kxdesc(ikx))
  write(lu,'(3x,a,5x,a,3x,5(a,3x))')	&
      &	'kt','10^','n','pres',		&
      &	'mean+/-stdv','mini(imin)','maxi(imax)'

end subroutine wtheader_

end module m_StatLevels
