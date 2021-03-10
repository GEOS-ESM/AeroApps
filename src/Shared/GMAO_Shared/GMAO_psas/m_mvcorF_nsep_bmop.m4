changequote([,])dnl
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_mvcorF_nsep_bmop - multi-var. corr. block matrix operators
!
! !DESCRIPTION:
!
! !INTERFACE:


    Module m_mvcorF_nsep_bmop
      Implicit none
      Private	! except

      ! These routines return the elements in an array
      Public :: fHHcor1
      Public :: fDDcor1
      Public :: fHHcorx
      Public :: fDDcorx
      Public :: fHDcorx

      ! These routines fuse BLAS mat-vec routines with
      ! the computation of block elements for more
      ! efficient use of cache/memory.
      ! Each of these comes in two varieties depending
      ! on whether an ensemble analysis is being performed.
      ! Additionally, except for the *HH* cases, multiple
      ! blocks are treated simultaneously to take advantage
      ! of common table lookups.  (dh/dx, dh/dy)

      Public :: sHH1cxpy, xHH1cxpy, rHH1cxpy
      Public :: xHD1cxpy, rHD1cxpy, rDH1cxpy
      Public :: sDD1cxpy, xDD1cxpy, rDD1cxpy

      Public :: sHHmcxpy, xHHmcxpy, rHHmcxpy
      Public :: xHDmcxpy, rHDmcxpy, rDHmcxpy
      Public :: sDDmcxpy, xDDmcxpy, rDDmcxpy

! !REVISION HISTORY:
!       20Nov01 - Jing Guo
!               . Merged in Greg Gaspari''s non-separable
!		  implementation.
!               . Replaced COMMON blocks with a restructured module.
!       05Mar01 - Tom Clune <clune@sgi.com>
!               . significant redesign to accomodate more interfaces while
!                 simultaneously reducing maintenance complexity
!                 "m4" now used to reuse common source text.
!               . added new routines for the case nvec=1 for performance reasons
!               . fused 4 DD blocks, 2 HD blocks, and 2 HD blocks
!       23Feb01 - Tom Clune <clune@sgi.com>
!               . significant redesign to directly call
!                 low level routines that comput elements
!                 The BLAS call is "fused" into those routines
!                 thereby eliminating the need to store the entire
!                 block.  (Improves cache usage!)
!	29Aug00	- Jing Guo
!		- initial prototype/prolog/code
!		- combined fHHcor1.F etc. into this module
!EOP ___________________________________________________________________

      Character(Len=*), Parameter :: myname='m_mvcorF_nsep_bmop'

!-----------------------------------------------------------------------


! Long section of m4 macros for constructing routines with maximal
! reuse of source.

!--------------------------------------------
! Special treatment for nvec = 1 significantly
! speeds execution.
!--------------------------------------------

!--------------------------------------------
define([M4_SELECT_KMAT], [dnl
   call mvcorF_nsep_getid(kind_cov,kmat)])


!-----------------------------------------------------
define([M4_LOOP_OVER_BLOCK], [dnl
ifelse( OPERATION, SYM_FILL,	[dnl
M4_SYM_LOOP()],
	OPERATION, ASY_FILL,	[dnl
M4_ASY_LOOP(i,j)],
	OPERATION, SCXPY,	[dnl
M4_SYM_LOOP()],
	OPERATION, XCXPY,	[dnl
M4_ASY_LOOP(i,j)],
	OPERATION, RCXPY,	[dnl
M4_ASY_LOOP(i,j)])])

!-----------------------------------------------------
define([M4_SYM_LOOP], [dnl
   l = 0

   Outer: Do j = 1, ln
M4_UNPACK_QR(j)
ifelse(BLOCK_TYPE, DD, [M4_UNPACK_QD(j)])
      kj = kl(j)

      Inner: Do i = 1, j - 1
		! l is an array element [index] of a packed upper
		! triangular matrix
	 l = l + 1
M4_UNPACK_QR(i)
ifelse(BLOCK_TYPE, DD, [M4_UNPACK_QD(i)])
	 ki = kl(i)
M4_COMPUTE_ELEMENT_IJ

      End Do Inner

      l = l + 1; i = j
M4_UNPACK_QR(i)
ifelse(BLOCK_TYPE, DD, [M4_UNPACK_QD(i)])
	 ki = kl(i)

M4_COMPUTE_DIAG_ELEMENT
   End Do Outer])

dnl !-----------------------------------------------------
dnl ! The loop order can be changed by switching the arguments
dnl ! of the macro.  i.e. changing (i,j) to (j,i) will change
dnl ! the loop order from (do i; do j) to (do j; do i).
dnl !-----------------------------------------------------
define([M4_ASY_LOOP], [dnl
   Outer: Do $1 = 1, ln$1
M4_UNPACK_QR($1)
ifelse(	BLOCK_TYPE, DD, [dnl
M4_UNPACK_QD($1)],
	BLOCK_TYPE, HD, [dnl
ifelse($1,j,[M4_UNPACK_QD($1)])],
	BLOCK_TYPE, DH, [dnl
ifelse($1,i,[M4_UNPACK_QD($1)])])
	 k$1 = kl$1($1)

      Inner: Do $2 = 1, ln$2
M4_UNPACK_QR($2)
ifelse(	BLOCK_TYPE, DD, [dnl
M4_UNPACK_QD($2)],
	BLOCK_TYPE, HD, [dnl
ifelse($2,j,[M4_UNPACK_QD($2)])],
	BLOCK_TYPE, DH, [dnl
ifelse($2,i,[M4_UNPACK_QD($2)])])
         k$2 = kl$2($2)

M4_COMPUTE_ELEMENT_IJ()
      End Do Inner

   End Do Outer])


dnl !-----------------------------------------------------
dnl ! Compute inner product of unpacked quantities
dnl !-----------------------------------------------------
define(DOT,$1_x*$2_x + $1_y*$2_y + $1_z*$2_z) dnl

dnl !-----------------------------------------------------
dnl ! This macro unpacks a 3 vector into scalar variables
dnl ! as a hint to the compiler to keep them in register.
dnl !-----------------------------------------------------
define([M4_UNPACK_QR], [dnl
ifdef([SYMMETRIC], [dnl
         qr$1_x = qr(1,$1)
	 qr$1_y = qr(2,$1)
	 qr$1_z = qr(3,$1)],[dnl
	 qr$1_x = qr$1(1,$1)
	 qr$1_y = qr$1(2,$1)
	 qr$1_z = qr$1(3,$1)])])


dnl !-----------------------------------------------------
dnl ! This macro unpacks a 3 vector into scalar variables
dnl ! as a hint to the compiler to keep them in register.
dnl !-----------------------------------------------------
define([M4_UNPACK_QD], [dnl
ifdef([SYMMETRIC], [dnl
ifdef([FILL],[dnl
         qm$1_x = qd(1,$1)
         qm$1_y = qd(2,$1)
         qm$1_z = qd(3,$1)],[dnl
         qm$1[]1_x = qd(1,$1,1)
         qm$1[]1_y = qd(2,$1,1)
         qm$1[]1_z = qd(3,$1,1)
         qm$1[]2_x = qd(1,$1,2)
         qm$1[]2_y = qd(2,$1,2)
         qm$1[]2_z = qd(3,$1,2)])],[dnl
ifdef([FILL], [dnl
         qm$1_x = qd$1(1,$1)
         qm$1_y = qd$1(2,$1)
         qm$1_z = qd$1(3,$1)],[dnl
         qm$1[]1_x = qd$1(1,$1,1)
         qm$1[]1_y = qd$1(2,$1,1)
         qm$1[]1_z = qd$1(3,$1,1)
         qm$1[]2_x = qd$1(1,$1,2)
         qm$1[]2_y = qd$1(2,$1,2)
         qm$1[]2_z = qd$1(3,$1,2)])])])

dnl !-----------------------------------------------------
dnl ! Compute "itau" which is an index into the correlation tables
dnl ! based upon "tau".
dnl !-----------------------------------------------------
define([M4_COMPUTE_ITAU], [dnl
       ctau=1.-tau
       itau=int(qxHtb(kmat)*ctau)+1	! not the nearest?!

       		! Scale vercor with redwin, a windowing factor to
		! modulate the horrizontal correlation function.  It
		! is also known as "reduced" PSAS.

       ktau=int(qxWtb*ctau+.5)+1
       redwin_=0.
       if(ktau<=mxWtb) redwin_=redwin(ktau)
       vercor=vercor*redwin_
])
dnl
dnl !----------------------------------------------
dnl ! Macro to compute the vertical correlation from
dnl ! lookup tables.
dnl !----------------------------------------------
#define vfecDH(ki,kj,kmat) vfecHD(kj,ki,kmat)
define([M4_COMPUTE_VERCOR], [ dnl
vercor = vfec[]BLOCK_TYPE (ki,kj,kmat)
])
dnl
dnl
dnl !-----------------------------------------------------
dnl ! This is a master macro for getting the coefficients
dnl ! used to compute the horizontal correlation.
dnl ! The are computed separately, since they may be reused
dnl ! for related blocks
dnl !-----------------------------------------------------
define([M4_COMPUTE_HORZ_COEF], [dnl
ifelse(BLOCK_TYPE,HH, [dnl
M4_COMPUTE_HORZ_COEF_HH],
       BLOCK_TYPE,DD, [dnl
M4_COMPUTE_HORZ_COEF_DD],
       BLOCK_TYPE,DH, [dnl
M4_COMPUTE_HORZ_COEF_DH], [dnl
M4_COMPUTE_HORZ_COEF_HD])])


dnl !-----------------------------------------------------
dnl ! HH coefficient used in HH horizontal correlation
dnl !-----------------------------------------------------
define([M4_COMPUTE_HORZ_COEF_HH], [dnl
    HH=winHH(itau)/(1.+ctau/rLav2(kj,ki,kmat))
])
dnl
dnl !-----------------------------------------------------
dnl ! HD coefficient used in HD horizontal correlation
dnl !-----------------------------------------------------
define([M4_COMPUTE_HORZ_COEF_HD],[dnl
    rLij2 = rLav2(kj,ki,kmat)
    PL = 1./(1. + ctau/rLij2)				! powerlaw in
    HR = PL*(DwinHH(itau) + winHH(itau)*PL/rLij2)	! derivatives
])
dnl
dnl !-----------------------------------------------------
dnl ! DH coefficient used in DH horizontal correlation
dnl !-----------------------------------------------------
define([M4_COMPUTE_HORZ_COEF_DH],[dnl
    rLij2 = rLav2(ki,kj,kmat)
    PL = 1./(1. + ctau/rLij2)				! powerlaw in
    HR = PL*(DwinHH(itau) + winHH(itau)*PL/rLij2)	! derivatives
])
dnl
dnl
dnl !-----------------------------------------------------
dnl ! DD coefficient used in HH horizontal correlation
dnl !-----------------------------------------------------
define([M4_COMPUTE_HORZ_COEF_DD], [dnl
    rLij2 = rLav2(kj,ki,kmat)
    PL = 1./(1. + ctau/rLij2)	! powerlaw: PL(u), u=1.-ctau
    PLa= PL/rLij2		! PL*PLa = "PL'(u)"
    RR = PL*(DwinHH(itau) + winHH(itau)*PLa)	! HR in fHDcorx.F
    TT = PL*(DDwinHH(itau) +	&
	 2.*PLa*(DwinH(itau) + PLa*winH(itau)) )
])


dnl !-----------------------------------------------------
dnl ! Having already computed the coefficients from the
dnl ! lookup tables, compute horcor and vercor and then
dnl ! the matrix element itself.  (Perhaps computing
dnl ! several elements for blocks that use the same 
dnl ! coefficients.)
dnl !-----------------------------------------------------
define([M4_COMPUTE_ELEM_FROM_COEFS], [dnl
ifelse(BLOCK_TYPE, HH,     [dnl
M4_COMPUTE_HH_ELEM],
       BLOCK_TYPE, DD,  [dnl 
ifdef([FILL], [dnl
dnl    ! only do lone block
M4_COMPUTE_DD_1_ELEM], [dnl
dnl    ! do 2 blocks at once
M4_COMPUTE_DD_2X2_ELEM])],[
ifelse(BLOCK_TYPE, HD,  [dnl 
ifelse(OPERATION, ASY_FILL,[dnl
dnl    ! only do lone block
M4_COMPUTE_HD_ELEM()],     [dnl
dnl    ! do 2 blocks at once
M4_COMPUTE_HD_ELEM(1)
M4_COMPUTE_HD_ELEM(2)])],
       BLOCK_TYPE, DH,   [dnl
ifelse(OPERATION, ASY_FILL,[dnl
dnl    ! only do lone block
M4_COMPUTE_DH_ELEM()],     [dnl
dnl    ! do 2 blocks at once
M4_COMPUTE_DH_ELEM(1)
M4_COMPUTE_DH_ELEM(2)])])])])


dnl !-----------------------------------------------------
define([M4_COMPUTE_HH_ELEM], [dnl
       horcor = CLfact(kj,ki,kmat) * HH
       element = horcor*vercor])

dnl !-----------------------------------------------------
define([M4_COMPUTE_HD_ELEM], [dnl
       dot_rm = DOT(qri,qmj$1)
       horcor = CLfact(kj,ki,kmat) * HR * dot_rm
       element$1 = horcor * vercor])

dnl !-----------------------------------------------------
define([M4_COMPUTE_DH_ELEM], [dnl
       dot_mr = DOT(qmi$1,qrj)
       horcor = CLfact(ki,kj,kmat) * HR * dot_mr
       element$1 = horcor * vercor])

dnl !-----------------------------------------------------
define([M4_COMPUTE_DD_1_ELEM], [dnl
      dot_rm = DOT(qri,qmj)
      dot_mr = DOT(qmi,qrj)
      dot_mm = DOT(qmi,qmj)
      horcor = CLfact(kj,ki,kmat) *	&
	(RR * dot_mm + TT * dot_rm * dot_mr)
      element = horcor*vercor])

dnl !-----------------------------------------------------
define([M4_COMPUTE_DD_2X2_ELEM], [dnl
      dot_rm_1 = DOT(qri,qmj1)
      dot_mr_1 = DOT(qmi1,qrj)
      dot_mm = DOT(qmi1,qmj1)
      horcor = CLfact(kj,ki,kmat) *	&
	(RR * dot_mm + TT * dot_rm_1 * dot_mr_1)
      element_11 = horcor*vercor

      dot_rm_2 = DOT(qri,qmj2)
      dot_mm = DOT(qmi1,qmj2)
      horcor = CLfact(kj,ki,kmat) *	&
	(RR * dot_mm + TT * dot_rm_2 * dot_mr_1)
      element_12 = horcor*vercor

      dot_mr_2 = DOT(qmi2,qrj)
      dot_mm = DOT(qmi2,qmj1)
      horcor = CLfact(kj,ki,kmat) *	&
	(RR * dot_mm + TT * dot_rm_1 * dot_mr_2)
      element_21 = horcor*vercor

      dot_mm = DOT(qmi2,qmj2)
      horcor = CLfact(kj,ki,kmat) *	&
	(RR * dot_mm + TT * dot_rm_2 * dot_mr_2)
      element_22 = horcor*vercor])


dnl! -----------------------------------------------------
dnl! 
dnl!-----------------------------------------------------
define([M4_COMPUTE_ELEMENT_IJ],[dnl
		! tau is the polarity [index] of two (i and j)
		! location vectors.  It is defined as
		! (q_i,q_j) = cos(angular_separation_of_i&j)
    tau = DOT(qri,qrj)
M4_COMPUTE_VERCOR
    If (tau > Hcoslim(kmat) .and. vercor /= 0.) Then
		    ! tau>cos(dlim/rade) iff ctau<ctaus(nHHtab)
M4_COMPUTE_ITAU
M4_COMPUTE_HORZ_COEF
M4_COMPUTE_ELEM_FROM_COEFS
M4_USE_ELEMENT
ifelse(OPERATION,  SYM_FILL, [dnl
    Else
       corr(l) = 0.],
       OPERATION, ASY_FILL, [dnl
    Else
       corr(i,j)=0.])
    End If])

dnl! -----------------------------------------------------
dnl! 
dnl!-----------------------------------------------------
define([M4_COMPUTE_DIAG_ELEMENT], [dnl
ifdef([FILL],[dnl
      corr(l) = 1.],[dnl
ifelse(BLOCK_TYPE, HH, [dnl
dnl ! HH case
      Y(:,j) = Y(:,j) + X(:,j)],[dnl
dnl ! DD case
		! tau is the polarity [index] of two (i and j)
		! location vectors.  It is defined as
		! (q_i,q_j) = cos(angular_separation_of_i&j)
      tau = DOT(qri,qrj)
      ki = kj

M4_COMPUTE_VERCOR
    If (tau > Hcoslim(kmat) .and. vercor /= 0.) Then
		    ! tau>cos(dlim/rade) iff ctau<ctaus(nHHtab)
M4_COMPUTE_ITAU
M4_COMPUTE_HORZ_COEF

M4_COMPUTE_ELEM_FROM_COEFS
         Y(:,j,1) = Y(:,j,1) + X(:,j,1) + element_12 * X(:,j,2)
         Y(:,j,2) = Y(:,j,2) + X(:,j,2) + element_21 * X(:,j,1)
      Else	
         Y(:,j,1) = Y(:,j,1) + X(:,j,1)
         Y(:,j,2) = Y(:,j,2) + X(:,j,2)
      End If
])])])


dnl! -----------------------------------------------------
dnl! 
dnl!-----------------------------------------------------
define([M4_USE_ELEMENT], [dnl
ifelse(OPERATION, SYM_FILL, [       corr(l) = element],
       OPERATION, ASY_FILL, [       corr(i,j) = element],[dnl
ifelse(OPERATION, SCXPY, [dnl
! SCXPY
ifelse(BLOCK_TYPE, HH, [dnl
           Y(:,i) = Y(:,i) + element * X(:,j)
	   Y(:,j) = Y(:,j) + element * X(:,i)],
       BLOCK_TYPE, DD, [dnl
           Y(:,i,1) = Y(:,i,1) + element_11 * X(:,j,1) + element_12 * X(:,j,2)
           Y(:,j,1) = Y(:,j,1) + element_11 * X(:,i,1) + element_21 * X(:,i,2)
	   Y(:,i,2) = Y(:,i,2) + element_21 * X(:,j,1) + element_22 * X(:,j,2)
           Y(:,j,2) = Y(:,j,2) + element_12 * X(:,i,1) + element_22 * X(:,i,2)])],
       OPERATION, XCXPY, [dnl
! XCXPY
ifelse(BLOCK_TYPE, HH, [dnl
           YI(:,i) = YI(:,i) + element * XJ(:,j)
	   YJ(:,j) = YJ(:,j) + element * XI(:,i)],
       BLOCK_TYPE, HD, [dnl
           YI(:,i)   = YI(:,i) + element1 * XJ(:,j,1) + element2 * XJ(:,j,2)
	   YJ(:,j,1) = YJ(:,j,1) + element1 * XI(:,i)
	   YJ(:,j,2) = YJ(:,j,2) + element2 * XI(:,i)], [dnl
dnl
	   YI(:,i,1) = YI(:,i,1) + element_11 * XJ(:,j,1) + element_12 * XJ(:,j,2)
	   YI(:,i,2) = YI(:,i,2) + element_21 * XJ(:,j,1) + element_22 * XJ(:,j,2)
	   YJ(:,j,1) = YJ(:,j,1) + element_11 * XI(:,i,1) + element_21 * XI(:,i,2)
	   YJ(:,j,2) = YJ(:,j,2) + element_12 * XI(:,i,1) + element_22 * XI(:,i,2)])],[dnl
dnl rCxpy
! RCXPY
ifelse(BLOCK_TYPE, HH, [dnl
           Y(:,i) = Y(:,i) + element * X(:,j)],
       BLOCK_TYPE, HD, [dnl
           Y(:,i) = Y(:,i) + element1 * X(:,j,1) + element2 * X(:,j,2)],[dnl
ifelse(BLOCK_TYPE, DH,[dnl
           Y(:,i,1) = Y(:,i,1) + element1 * X(:,j)
	   Y(:,i,2) = Y(:,i,2) + element2 * X(:,j)],[dnl

	   Y(:,i,1) = Y(:,i,1) + element_11 * X(:,j,1) + element_12 * X(:,j,2)
	   Y(:,i,2) = Y(:,i,2) + element_21 * X(:,j,1) + element_22 * X(:,j,2)])])])])])


dnl!-----------------------------------------------------
dnl! macro arglist:
dnl!   $1: nvecs
dnl!   $2: blocktype  HH/HD/DD
dnl!   $3: operation  SCXPY, XCXPY, RCXPY, SYM_FILL, ASY_FILL
dnl!-----------------------------------------------------
define([M4_ROUTINE_NAME],[dnl
ifelse($3,SCXPY,s$2[]ifelse($1,1,1,m)cxpy,$3,XCXPY,x[]$2[]ifelse($1,1,1,m)cxpy,[dnl
ifelse($3,RCXPY,r$2[]ifelse($1,1,1,m)cxpy,$3,ASY_FILL,f$2corx,f$2cor1)])])

dnl!-----------------------------------------------------
define([M4_ARG_LIST],[dnl
ifelse($3, SCXPY, [dnl
          ln, qr, ifelse($2, DD, [qd,]) kl, nvecs, x, y],
       $3, SYM_FILL, [dnl
          ln, qr, ifelse($2, DD, [qd,]) kl, corr], [dnl
ifelse($3, ASY_FILL, [dnl
          lni, qri, ifelse($2, DD, [qdi,], $2, DH, [qdi,]) kli, &
          lnj, qrj, ifelse($2, DD, [qdj,], $2, HD, [qdj,]) klj,	&
	  corr],
       $3, XCXPY,  [dnl
          lni, qri, ifelse($2, DD, [qdi,], $2, DH, [qdi,]) kli,	&
          lnj, qrj, ifelse($2, DD, [qdj,], $2, HD, [qdj,]) klj,	&
	  nvecs, xi, yj, xj, yi], [dnl
dnl RCXPY
          lni, qri, ifelse($2, DD, [qdi,], $2, DH, [qdi,]) kli,	&
          lnj, qrj, ifelse($2, DD, [qdj,], $2, HD, [qdj,]) klj,	&
	  nvecs, x, y])])])

dnl!-----------------------------------------------------
define([M4_QD_DECLARATION],
[dnl
  ifdef([FILL],
  [
     Real   , Intent(In)   :: qd(3, ln)],
  [
     Real   , Intent(In)   :: qd(3, ln, 2)])
])
dnl!-----------------------------------------------------
define([M4_QD_DECLARATION_IJ],
[dnl
  ifdef([FILL],
  [
     Real   , Intent(In)   :: qd[]$1(3, ln[]$1)],
  [
     Real   , Intent(In)   :: qd[]$1(3, ln[]$1, 2)])
])
dnl!-----------------------------------------------------
define([M4_SYM_ATTRIBUTES_DECLARATIONS],
[
     Integer, Intent(In)   :: ln
     Real,    Intent(In)   :: qr(3, ln) dnl
  ifelse(BLOCK_TYPE, DD, [M4_QD_DECLARATION])
     Integer, Intent(In)   :: kl(ln)
])

dnl!-----------------------------------------------------
define([M4_ASY_ATTRIBUTES_DECLARATIONS],
[
dnl -row attributes-
     Integer, Intent(In)   :: lni
     Real,    Intent(In)   :: qri(3, lni)
dnl
  ifelse(
  BLOCK_TYPE, DH, [M4_QD_DECLARATION_IJ(i)],
  BLOCK_TYPE, DD, [M4_QD_DECLARATION_IJ(i)])
     Integer, Intent(In)   :: kli(lni)
dnl
dnl -column attributes-
     Integer, Intent(In)   :: lnj
     Real,    Intent(In)   :: qrj(3, lnj)
dnl
  ifelse(
  BLOCK_TYPE, HD, [M4_QD_DECLARATION_IJ(j)],
  BLOCK_TYPE, DD, [M4_QD_DECLARATION_IJ(j)])
     Integer, Intent(In)   :: klj(lnj)
])

dnl!-----------------------------------------------------
define([M4_ARGUMENT_DECLARATION], [
     Integer, Intent(In)   :: kind_cov dnl
dnl -attributes-
  ifdef([SYMMETRIC],
[M4_SYM_ATTRIBUTES_DECLARATIONS],
[M4_ASY_ATTRIBUTES_DECLARATIONS])
dnl
dnl -block storage or rhs/lhs vectors-
  ifdef([FILL],
[M4_BLOCKSTORAGE_DECLARATIONS],
[M4_VECTORS_DECLARATIONS])
])

dnl!-----------------------------------------------------
define([M4_BLOCKSTORAGE_DECLARATIONS], [dnl
  ifdef([SYMMETRIC],
  [
     Real   , Intent(Out)  :: corr(ln*(ln+1)/2)],
  [
     Real   , Intent(Out)  :: corr(lni,lnj)])
])
dnl!-----------------------------------------------------
define([M4_VECTORS_DECLARATIONS], [
     Integer, Intent(In)   :: nvecs
  dnl
  ifelse(
  OPERATION, SCXPY,
  [dnl
    ifelse(
    BLOCK_TYPE, HH, [
     Real, Intent(In)    :: X(nvecs, ln)
     Real, Intent(InOut) :: Y(nvecs, ln)],
    BLOCK_TYPE, DD, [
     Real, Intent(In)    :: X(nvecs, ln, 2)
     Real, Intent(InOut) :: Y(nvecs, ln, 2)])
  ],
  OPERATION, RCXPY,
  [dnl
    ifelse(
    BLOCK_TYPE, HH, [
     Real, Intent(In)    :: X(nvecs, lnj)
     Real, Intent(InOut) :: Y(nvecs, lni)],
    BLOCK_TYPE, DD, [
     Real, Intent(In)    :: X(nvecs, lnj, 2)
     Real, Intent(InOut) :: Y(nvecs, lni, 2)],
    BLOCK_TYPE, HD, [
     Real, Intent(In)    :: X(nvecs, lnj, 2)
     Real, Intent(InOut) :: Y(nvecs, lni)],
    BLOCK_TYPE, DH, [
     Real, Intent(In)    :: X(nvecs, lnj)
     Real, Intent(InOut) :: Y(nvecs, lni, 2)])
  ],
  OPERATION, XCXPY,
  [dnl
    ifelse(
    BLOCK_TYPE, HH, [
     Real, Intent(In)    :: XI(nvecs, lni)
     Real, Intent(InOut) :: YJ(nvecs, lnj)
     Real, Intent(In)    :: XJ(nvecs, lnj)
     Real, Intent(InOut) :: YI(nvecs, lni)],
    BLOCK_TYPE, DD, [
     Real, Intent(In)    :: XI(nvecs, lni, 2)
     Real, Intent(InOut) :: YJ(nvecs, lnj, 2)
     Real, Intent(In)    :: XJ(nvecs, lnj, 2)
     Real, Intent(InOut) :: YI(nvecs, lni, 2)],
    BLOCK_TYPE, HD, [
     Real, Intent(In)    :: XI(nvecs, lni)
     Real, Intent(InOut) :: YJ(nvecs, lnj, 2)
     Real, Intent(In)    :: XJ(nvecs, lnj, 2)
     Real, Intent(InOut) :: YI(nvecs, lni)])
  ]) dnl --ifelse(OPERATION,..)--
]) dnl --define(M4_VECTORS_DECLARATIONS)--


dnl!-----------------------------------------------------
dnl! M4 arglist:
dnl!   $1: nvecs
dnl!   $2: blocktype  HH/HD/DD
dnl!   $3: operation  SCXPY, XCXPY, RCXPY, ASY_FILL, SYM_FILL
dnl!-----------------------------------------------------
dnl
dnl============================================================
define([M4_CXPY],[dnl
dnl============================================================
dnl
define(BLOCK_TYPE,[$2]) dnl
define(OPERATION,[$3]) dnl
ifelse($1,1, [dnl
define([X],x([shift]([$]*)))   dnl
define([Y],y([shift]([$]*)))   dnl
define([XI],xi([shift]([$]*))) dnl
define([XJ],xj([shift]([$]*))) dnl
define([YI],yi([shift]([$]*))) dnl
define([YJ],yj([shift]([$]*)))]) dnl
ifelse($3,SCXPY,[define(SYMMETRIC)],$3,SYM_FILL,[define(SYMMETRIC)])dnl
ifelse($3,SYM_FILL,[define(FILL)],$3,ASY_FILL,[define(FILL)])dnl
 
   Subroutine M4_ROUTINE_NAME($*)(kind_cov,	&
	M4_ARG_LIST($*))
     Use m_mvcorF_nsep,only : mvcorF_nsep_getid
     use m_mvcorF_nsep,only : vfecHH
     use m_mvcorF_nsep,only : vfecHD
     use m_mvcorF_nsep,only : vfecDD
     use m_mvcorF_nsep,only : Hcoslim
     use m_mvcorF_nsep,only : qxHtb
     use m_mvcorF_nsep,only : rLav2
     use m_mvcorF_nsep,only : CLfact
     use m_mvcorF_nsep,only :   winHH
     use m_mvcorF_nsep,only :  DwinHH
     use m_mvcorF_nsep,only : DDwinHH
     use m_mvcorF_nsep,only :   winH
     use m_mvcorF_nsep,only :  DwinH

     		! Reduced PSAS window function table.

     use m_redwin,only : redwin	! a table of window function values.
     use m_redwin,only : qxWtb	! the unit length of tau for the table
     use m_redwin,only : mxWtb	! the maximum index of the table.
     use m_redwin,only : redwin_initialized
     use m_die,only : die

     Implicit None

M4_ARGUMENT_DECLARATION($*)

!----------------------------------------------------------------------
! Local Variables inclusively for all cases

      Character(Len=*), Parameter :: myname_ = myname//'::M4_ROUTINE_NAME($*)'

      Integer :: kmat
      Integer :: i, j, l
      Integer :: ki, kj
      Real    :: qri_x, qri_y, qri_z
      Real    :: qrj_x, qrj_y, qrj_z
      Real    :: qmi_x, qmi_y, qmi_z
      Real    :: qmj_x, qmj_y, qmj_z
      Real    :: qmi1_x, qmi1_y, qmi1_z
      Real    :: qmj1_x, qmj1_y, qmj1_z
      Real    :: qmi2_x, qmi2_y, qmi2_z
      Real    :: qmj2_x, qmj2_y, qmj2_z
      Real    :: tau, ctau, xtau
      Integer :: itau, ktau
      Real    :: HH, HR, RR, TT
      Real    :: horcor, vercor, redwin_
      Real    :: dot_rr
      Real    :: dot_mr, dot_mr_1, dot_mr_2
      Real    :: dot_rm, dot_rm_1, dot_rm_2
      Real    :: dot_mm
      Real    :: element, element1, element2
      Real    :: element_11, element_12, element_21, element_22
      Real    :: PL,PLa,rLij2

!======================================================================
	if(.not.redwin_initialized())	&
		call die(myname_,'redwin not initialized')

ifelse($3,SCXPY,[dnl
      If (ln <= 0) Return],
       $3, SYM_FILL,[dnl
      If (ln <= 0) Return],[dnl
      If (lni <=0 .or. lnj <= 0) Return])
ifelse($1,1,,[dnl
      If (nvecs <= 0) Return])

M4_SELECT_KMAT()
M4_LOOP_OVER_BLOCK()

End Subroutine M4_ROUTINE_NAME($*)

ifdef([FILL],[undefine([FILL])])dnl
ifdef([SYMMETRIC],[undefine([SYMMETRIC])])dnl
undefine([OPERATION])  dnl
undefine([BLOCK_TYPE])  dnl
ifelse($1,1,[dnl
undefine([X])  dnl
undefine([Y])  dnl
undefine([XI]) dnl
undefine([XJ]) dnl
undefine([YI]) dnl
undefine([YJ])])
])

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  M4_CXPY(1, HH, SYM_FILL)
  M4_CXPY(1, DD, SYM_FILL)
  M4_CXPY(1, HH, ASY_FILL)
  M4_CXPY(1, HD, ASY_FILL)
  M4_CXPY(1, DD, ASY_FILL)
  
  
  M4_CXPY(1, HH, SCXPY)
  M4_CXPY(n, HH, SCXPY)
  M4_CXPY(1, DD, SCXPY)
  M4_CXPY(n, DD, SCXPY)
  
  M4_CXPY(1, HH, XCXPY)
  M4_CXPY(n, HH, XCXPY)
  M4_CXPY(1, HD, XCXPY)
  M4_CXPY(n, HD, XCXPY)
  M4_CXPY(1, DD, XCXPY)
  M4_CXPY(n, DD, XCXPY)
  
  M4_CXPY(1, HH, RCXPY)
  M4_CXPY(n, HH, RCXPY)
  M4_CXPY(1, HD, RCXPY)
  M4_CXPY(n, HD, RCXPY)
  M4_CXPY(1, DH, RCXPY)
  M4_CXPY(n, DH, RCXPY)
  M4_CXPY(1, DD, RCXPY)
  M4_CXPY(n, DD, RCXPY)


End Module m_mvcorF_nsep_bmop



