!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_mvcorF_nsep_bldr - Build m_mvcorF_nsep for block operations
!
! !DESCRIPTION:
!   This module initialize the data in module m_mvcorF_nsep, such as,
!   matrix table vfecHH, tabulated window function and derivatives, and
!   form a matrix of length scales.
!
! !INTERFACE:

    module m_mvcorF_nsep_bldr
      implicit none
      private	! except

      public :: mvcorF_nsep_init
      public :: mvcorF_nsep_clean

      interface mvcorF_nsep_init ; module procedure	&
	set_fecHH; end interface
      interface mvcorF_nsep_clean; module procedure	&
	clean_   ; end interface

! !REVISION HISTORY:
! 	19Nov01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- adapted from set_fecHH.F90 for non-seperable multi-
!		  variate correlation operators.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_mvcorF_nsep_bldr'

  integer,parameter :: nHHtab=15000

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: set_fecHH - initialize non-sep. multi-var. err. cor. ops.
!
! !DESCRIPTION:
!
!   set_fecHH() computes the matrix table vfecHH, and calls intp_hCorb
!   to tabulate the window function and its derivatives winHH, DwinHH,
!   and DDwinHH at a single level.  The length scales and the multiplier
!   for the cross-correlation function are also tabulated. 
!
! !INTERFACE:

    subroutine set_fecHH()
      use hfecH_tbl
      use vfecH_tbl
      use m_mvcorF_nsep

      use rlev_imat,only : MXveclev,nveclev,pveclev

!use config, only   : MXpar_hc, lvmax_vc, MX_fecH

      use const    ,only : R_Earth
      use m_chars  ,only : uppercase
      use m_realkinds,only : kind_R8
      use m_mpout  ,only : mpout,mpout_flush
      use m_stdio  ,only : stdout,stderr
      use m_die    ,only : die
      implicit none

! !REVISION HISTORY:
!       21Feb2001- G. Gaspari  Designed the new routine
!       27Nov98 - J. Larson  Incorporated Tom Clune's Cache 
!                            Optimizations:
!                              1) Fused Imat Tables hfecRR and 
!                                 hfecTT into hfecRRTT.
!                              2) Switched order of IMAT tau 
!                                 and level indices.
!			       3) Added call to Routine intp_hCorb()
!				  to accomodate optimizations 1) and 2).
! 	12Jan96 - J. Guo	- programmed and added the prolog
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::set_fecHH'
  include "kind_covs.h"
  include "MX_hfecH.h"

!_______________________________________________________________________

real,parameter :: rade = .001*R_Earth

integer lv,i,k,l,km,ii,jj
real c                    !  cutoff for winHH = .5*support of winHH
real HHmx, HHinc
real(kind_R8) :: ctaus(nHHtab)

real   CL(nveclev,nveclev)   ! multiplier: sqrt(Li*Lj)/[ .5*(Li+Lj) ]
real   rLa2(nveclev,nveclev) ! rLa2 = [ .5*(Li+Lj) ]**2

!  plev_hfecH(1,1)            = 1000
!  pars_hfecH(1,1,km)         = 6000,  pars_hfecH(2,1,km)          = 524.9
!     ...                                 ...
!  plev_hfecH(1,nlev_hfecH)   = .4
!  pars_hfecH(1,nlev_hfecH,km)= 6000,  pars_hfecH(2,nlev_hfecH,km) = 7026.5
!_______________________________________________________________________

  c=1000.                    ! A minimum cutoff distance
  do km=1,MX_fecH
    do lv=1,nlev_hfecH(km)
      if(pars_hfecH(1,lv,km).gt. c) c=.5*pars_hfecH(1,lv,km)
    end do
  end do

		! Counting in all three forecast error correlation
		! models

  call mvcorF_nsep_setid(kind_covF,reset=.true.)
  call mvcorF_nsep_setid(kind_covS)
  call mvcorF_nsep_setid(kind_covV)

#ifndef NDEBUG
  write(mpout,*) 'nHHtab= ',nHHtab
  write(mpout,*) 'cutoff= 1/2 the support=  ',c
#endif

		! Allocate the memory for the data storeage

  call mvcorF_nsep_alloc(nHHtab,nveclev)
!_______________________________________________________________________

! 
!  c = 3000;        Override input for testing. (remove later)
!
! window function ==0 for z = sqrt(2*ctaus)/cw >= 2, where cw:=c/rade, i.e.,
! ctaus >= 2*cw**2
!   
  HHmx = 2*(c/rade)**2          ! CORRECT: window fcn == 0 for ctaus > HHmx 
!  HHmx = 1 - cos(2.*c/rade)   ! incorrect - the current value assumed in psas
				! JG: I think it depends on one's
				! understanding of what this HHmx and
				! HHinc are.

  HHinc=HHmx/nHHtab             ! grid spacing

  ctaus(1) = 0.0                ! first table entry
  do i = 2,nHHtab
   ctaus(i) = HHinc*(i-1)       ! use a one-scale tabulation, res is HHinc
  end do

!_______________________________________________________________________

	! Loop for all hfecH tables
do km=1,MX_fecH

	! Verify the model name.

  if(uppercase(type_hfecH(km)) /= "WIN-POWERLAW")	&
    call die(myname_,'not supported, "'//trim(type_hfecH(km))//'"')

	! Save the resolution and the range, for future extensions.

  Hcoslim(km)=1.-ctaus(nHHtab)
				! set corr. to zero for 1-ctau > Hcoslim
                                ! i.e., ctau < ctaus(nHHtab) 
  qxHtb(km)  =1./HHinc		! used in m_mvcorF_nsep_bmop procedures

!________________________________________
!
!  intp_hCorb returns:
!  ==================
!   rLa2(nveclev,nveclev): {aveLij**2}, aveLij=.5*(Li+Lj), i,j=1..nveclev 
!   CL(nveclev,nveclev)  : {sqrt(Li*Lj)/aveLij}, i,j=1..nveclev
!   winHH, DwinHH, DDwinHH: window fcn, 1st deriv., 2nd deriv
!   tabulated on ctaus uniform grid of nHHtab points.  Tabulated only on
!   the support of winHH, i.e., for 0 <= ctaus <= HHmx:=2*(c/rade)**2 
! 
    call intp_hCorb(nlev_hfecH(km),plev_hfecH(1,km),		&
	size(pars_hfecH,1),pars_hfecH(1,1,km),			&
	nveclev,pveclev,nHHtab,ctaus,c,winHH,DwinHH,DDwinHH, rLa2,CL)

    winH(:) = winHH(:)              ! used to redefine W''(1) in fDDcorx.F
    DwinH(:) = DwinHH(:)            ! .. and in fDDcor1.F (i.e. TT[itau=1] 
    winH(1) = 0.0; DwinH(1) = 0.0   ! force TT = 0. for itau=1 only in fDD*.F 

    do ii = 1,nveclev
      do jj = ii,nveclev
       CLfact(jj,ii,km) = CL(jj,ii)  ! correct form, but not implementable.
       CLfact(ii,jj,km) = CL(jj,ii)
!       CLfact(jj,ii,km) = 1. ; CLfact(ii,jj,km) = 1.  ! mult. of cross-corr
       rLav2(jj,ii,km) = rLa2(jj,ii)
       rLav2(ii,jj,km) = rLa2(jj,ii)
!    if(km .eq. 1) then 
!      print *, ii,jj,pveclev(ii),pveclev(jj),rLav2(ii,jj,km) 
!    end if
      end do
    end do

!-----------------------------------------------------------------------
! vfecHH(p_i,p_j):  Tabulate the vertical correlations on pveclev

  call intp_vCor( name_vfecH(km),				&
		  size(plev_vfecH,1),nlev_vfecH(km),		&
		  plev_vfecH(1,km),corr_vfecH(1,1,km),		&
		  size(vfecHH,1),nveclev,pveclev, vfecHH(1,1,km))

!-----------------------------------------------------------------------
! Scale vfecHH for fecHH, fecHD and fecDD, such that
! vfecHH*hfecHH, vfecHD*hfecHR, and vfecDD*hfecRR are
! normalized.  
! Normalization factors (diagonal elements).  normDD are also
! used to determine Sigmas of a D variable (error stdev.).
!
!  See PSAS II, Eqns. (50-52)

! DWp0 = P'(1) + S'(1-),   P(u) = powerlaw, S(u) = window fcn, u=1-ctaus  
!      = 1/L**2 + S'(1-),          L = powerlaw length scale

  do k=1,nveclev
    normDD(k,km) = 1./sqrt( 1./ rLav2(k,k,km) + DwinHH(1) )
  end do

! normalize: vfecHD and vfecDD

  do l=1,nveclev
    do k=1,nveclev
      vfecHD(k,l,km) = vfecHH(k,l,km) * normDD(l,km)
      vfecDD(k,l,km) = vfecHH(k,l,km) * normDD(k,km) * normDD(l,km)
    end do
  end do        ! l=

!********1********2********3********4********5********6********7*
!    if(km .eq. 1) then
!     print *,'km=  ',km
!     print *, 'set_fecHH.F90: Tabulated window corr. fun.,first 500 rows '
!     print *,' '
!    do i = 1,500
!     print '(6F12.4)', winHH(1+6*(i-1):6*i)
!    end do
!    print *,'  '
!    do i = 1,500
!     print '(6F12.4)', DwinHH(1+6*(i-1):6*i)
!    end do
!    print *,'  '
!    do i = 1,500
!     print '(6F12.4)', DDwinHH(1+6*(i-1):6*i)
!    end do
!    print *,'  '
!     print '(6F12.4)', winHH(nHHtab-5:nHHtab)
!     print '(6F12.4)', DwinHH(nHHtab-5:nHHtab)
!     print '(6F12.4)', DDwinHH(nHHtab-5:nHHtab)
!    endif
!    print *,' '; print *,' '
!**********************************************************************

end do     ! km = 1,MX_fecHH

end subroutine set_fecHH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean the module
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_()
      use m_mvcorF_nsep,only : mvcorF_nsep_dealloc
      implicit none

! !REVISION HISTORY:
! 	19Nov01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  call mvcorF_nsep_dealloc()

end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: intp_hCorb - compute horizontal correlation components
!
! !DESCRIPTION:
!
!   Tabulates the window function GASPARI_COHN and computes the length
!   scale as a function of plev, where plev is the array of all pressure
!   levels used in the analysis. 
! 
! !INTERFACE:

    subroutine intp_hCorb(nlev_in,plev_in,mpar_in,pars_in,       &
                	nlev,   plev,   nctau,ctaus,c,    	 &
			hW,DhW,DDhW,rLa2,CL                      )

use const  ,only : R_Earth	! radius of Earth
use m_realkinds,only : kind_R8
use m_stdio,only : stderr	! standard error output unit
use m_die,  only : die

implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer, intent(in)	:: nlev_in	! no. of plev_in in psas.rc file
  real,    intent(in)	:: plev_in(nlev_in)	! input pressure levels
  integer, intent(in)   :: mpar_in      ! pars_in(1,:)=2*c, pars_in(2,:)=L 
  real,    intent(in)   :: pars_in(mpar_in,nlev_in) 
  integer, intent(in)	:: nlev		! no. of plev used in analysis
  real,    intent(in)	:: plev(nlev)	! the "nveclev" levels in the analysis
  integer, intent(in)	:: nctau	! size of ctaus, the horz. grid
  real(kind_R8),intent(in) :: ctaus(nctau)	! ctaus: the two scale horz. grid
  real,    intent(in)   :: c            ! cutoff = .5*support of hW

  real,    intent(out)  :: hW(nctau)        ! window function
  real,    intent(out)  :: DhW(nctau)       ! its first derivative
  real,    intent(out)  :: DDhW(nctau)      ! its second derivative
  real,    intent(out)  :: rLa2(nlev,nlev)    ! cross-corr scales squared
  real,    intent(out)  :: CL(nlev,nlev)      ! multiplier for cross-corr

! !REVISION HISTORY:
! 	19Nov01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::intp_hCorb'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Local workspace and variables

  integer i, k, jl
  real    rL(nlev)
  real    rLjl, rLjr, aveL, rm 
  real    z, z2, tmp, cw, cw2, cw4 

  real         a0,        a2,  a3, a4, a5
  real         a0d,  a1d, a2d, a3d
  real  am1dd, a0dd, a1dd
 
  real         bm1,   b0,    b1,   b2,  b3, b4, b5
  real  bm3d,  bm1d,  b0d,   b1d,  b2d, b3d
  real  bm5dd, bm3dd, bm1dd, b0dd, b1dd

  logical,	    parameter	:: roundoff=.true.

  real,parameter :: rade = .001*R_Earth
!_______________________________________________________________________

!***********************************************************************
!
!  Tabulate the window function  hW, and its first and
!  second derivatives DhW and DDhW, respectively.
!
    a0= 1. ; a2=-5./3.; a3= 5./8.; a4= 1./2.; a5=-1./4.  ! win fn coeffs.
    a0d = 10./3; a1d= -1.875; a2d=-2.0; a3d=1.25         ! deriv. coeffs.
    am1dd = 15./8.;  a0dd= 4.; a1dd= -15./4.             ! 2nd deriv. coeffs. 

    bm1=-2./3.; b0= 4.; b1=-5.; b2= 5./3.; b3= 5./8.     ! win. coeff.: [1,2]
    b4=-1./2.; b5= 1./12.                                ! win. coeff.: [1,2]
    bm3d = -2./3.; bm1d = 5.; b0d=-10./3.; b1d= -15./8.  ! deriv. coeff.: [1,2]
    b2d = 2.; b3d = -5./12.                              ! deriv. coeff.: [1,2]
    b1dd = 5./4.; b0dd = -4.; bm1dd = 15./8.; bm3dd=5.   ! 2nd deriv co.: [1,2]
    bm5dd= -2.                                           ! 2nd deriv co.: [1,2]

     cw = c/rade; cw2 = cw*cw; cw4 = cw2*cw2             ! normalize c to S^2    
     do i = 1,nctau
      z2 = 2.*ctaus(i)/cw2; z = sqrt(z2)             ! scale to [0,2]
      tmp = z2*(a4*z2 + a2)
      if(z.le.1.) then
       hW(i)  = z2*z*(a5*z2+a3) + tmp + a0
       DhW(i) = (z*(z*(z*a3d+a2d)+a1d)+a0d )/cw2
       if(z > 0.) Then
          DDhW(i)= (a1dd*z + a0dd + am1dd/z)/cw4 
       Else
          DDhW(i) = 0
       End if
      elseif(z > 1. .and. z < 2.) then
       hW(i) = z*(z2*(b5*z2+b3)+b1) + b0 + bm1/z - tmp
       DhW(i)= ( z*(z*(z*b3d+b2d)+b1d)+b0d               & 
     &          + (bm3d/z2+bm1d)/z ) / cw2 
       DDhW(i)= (((bm5dd/z2+bm3dd)/z2+bm1dd)/z  &
     &          + b0dd + b1dd*z)/ cw4 
      else
       hW(i) = 0.;  DhW(i) = 0.;  DDhW(i) = 0.
      endif
     end do
!***********************************************************
!     call Lscale(nlev_in,plev_in,mpar_in,pars_in,nlev,plev,c,rLa2,CL)
!*****************************************************
! Assign the length scale matrix rLa2 by linear interpolation
! from the input levels plev_in to the output levels plev
!
! Assume:  plev_in(1) > plev_in(2) > ... > plev_in(nlev_in)
!
! find rL(i) = length scale at plev(i) by log10 (input pressure plev_in)-
! linear (input length scale, pars_in(2,:) ) interpolation,
! e.g., the point (rL(i),y) on the line through (rLjl,yl) & (rLjr,yr)
! where rLjl = pars_in(2,jl),       rLjr = pars_in(2,jl-1),
!         yl = log10(plev_in(jl)),  yr = log10(plev_in(jr)), y = log10(plev(i))
!
 do i = 1,nlev
  jl = 1;                  ! index of lower pressure level
  do while( plev(i) <= plev_in(jl) .and. jl < nlev_in)
   jl = jl+1;
  end do
 
  if (jl == 1) rL(i)=pars_in(2,1)      ! plev_in(1) < plev(i)
  if (jl > 1 .and. jl < nlev_in) then  ! plev_in(jl)< plev(i) <= plev_in(jl-1)
   rLjl = pars_in(2,jl); rLjr = pars_in(2,jl-1)
   rm = (rLjr-rLjl)/log10(plev_in(jl-1)/plev_in(jl))
   rL(i) = rLjl + rm*log10(plev(i)/plev_in(jl))
  end if
  if (jl == nlev_in) then
   if (plev(i) <= plev_in(jl)) then    ! plev(i) smaller than all input levels
    rL(i)=pars_in(2,nlev_in)           ! .. set to smallest scale
   else                                ! plev_in(jl)< plev(i) <= plev_in(jl-1)
    rLjl = pars_in(2,jl); rLjr = pars_in(2,jl-1)
    rm = (rLjr-rLjl)/log10(plev_in(jl-1)/plev_in(jl))
    rL(i) = rLjl + rm*log10(plev(i)/plev_in(jl))
   end if
  end if
 end do     ! outer loop, i = 1,nlev
!
!  Now assign the squared length scales for the cross-correlation function
!
  do k = 1,nlev
   do i = k,nlev
    aveL = .5*(rL(k)+rL(i)); rLa2(i,k) = (aveL/rade)**2
    CL(i,k) = sqrt(rL(k)*rL(i))/aveL
   end do
  end do

!  print *,' inside intp_hCorb' 
!  do k = 1,nlev
!   print *, plev(k), rL(k)
!  end do

end subroutine intp_hCorb
end module m_mvcorF_nsep_bldr
