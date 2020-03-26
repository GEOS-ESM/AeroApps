changequote([,])
undefine(include)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_mvcorF_base_bmop - m.v.correlation matrix generators for Fcst.Err.
!
! !DESCRIPTION:
!
! !INTERFACE:


    Module m_mvcorF_base_bmop
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

      Character(Len=*), Parameter :: myname='m_mvcorF_base_bmop'

! Include files used for tables and such.

	Include "MX_hfecH.h"
	Include "lvmax.h"
	Include "vfecHH.h"

!-----------------------------------------------------------------------
		! indirect horizontal correlation matrix factors
#include "hfecHH.h"

!-----------------------------------------------------------------------
!   ..Options of correlation matrices

	Include "kind_covs.h"
!-----------------------------------------------------------------------


! Long section of m4 macros for constructing routines with maximal
! reuse of source.

!--------------------------------------------
! Special treatment for nvec = 1 significantly
! speeds execution.
!--------------------------------------------

!--------------------------------------------
define([M4_SELECT_KMAT], [dnl
   if(kind_cov.eq.kind_covF) then
      kmat=kmat_HGHT
   elseif(kind_cov.eq.kind_covS) then
      kmat=kmat_STRM
   elseif(kind_cov.eq.kind_covV) then
      kmat=kmat_VELP
   else
      return
   endif])


!-----------------------------------------------------
define([M4_LOOP_OVER_BLOCK], [dnl
ifdef([SYMMETRIC], [dnl
M4_SYM_LOOP],   [dnl
M4_ASY_LOOP])])

!-----------------------------------------------------
define([M4_SYM_LOOP], [dnl
   l = 1; j = 1; i = 1
ifelse(BLOCK_TYPE, DD, [M4_UNPACK_QD(j)])
M4_COMPUTE_DIAG_ELEMENT

   Outer: Do j = 2, ln
M4_UNPACK_QR(j)
ifelse(BLOCK_TYPE, DD, [M4_UNPACK_QD(j)])
      kj = kl(j)

      Inner: Do i = 1, j - 1
	 l = l + 1
	 ki = kl(i)
M4_UNPACK_QR(i)
M4_COMPUTE_ELEMENT_IJ

      End Do Inner

      l = l + 1; i = j
M4_COMPUTE_DIAG_ELEMENT

   End Do Outer])


define([M4_ASY_LOOP], [dnl
   Outer: Do j = 1, lnj
M4_UNPACK_QR(j)
ifelse(BLOCK_TYPE, DD, [M4_UNPACK_QD(j)])
ifelse(BLOCK_TYPE, HD, [M4_UNPACK_QD(j)])
      kj = klj(j)

      Inner: Do i = 1, lni
M4_UNPACK_QR(i)
         ki = kli(i)
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
#ifdef _LINEAR
       xtau=qxHtb1*ctau+1.
       if(ctau > HHbeg2) xtau=qxHtb2*(ctau-HHbeg2)+nHHtb1+1.
       itau=min(int(xtau),nHHtab-1)
       xtau=xtau-itau
#else
       xtau=qxHtb1*ctau+1.5
       if(ctau > HHbeg2) xtau=qxHtb2*(ctau-HHbeg2)+nHHtb1+1.5
       ! round off to the table grid.  Not checking
       ! is done here, but corrected results depend on
       ! itau in the range of [1,nHHtab].
       itau=int(xtau)
[#endif] ])
dnl
dnl !----------------------------------------------
dnl ! Macro to compute the vertical correlation from
dnl ! lookup tables.
dnl !----------------------------------------------
#define vfecDH(ki,kj,kmat) vfecHD(kj,ki,kmat)
define([M4_COMPUTE_VERCOR], [ dnl
vercor = weight_*vfec[]BLOCK_TYPE (ki,kj,kmat)
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
ifelse(BLOCK_TYPE,HH,[dnl
M4_COMPUTE_HORZ_COEF_HH],
       BLOCK_TYPE,DD, [dnl
M4_COMPUTE_HORZ_COEF_DD],[dnl
M4_COMPUTE_HORZ_COEF_HD])])


dnl !-----------------------------------------------------
dnl ! HH coefficient used in HH horizontal correlation
dnl !-----------------------------------------------------
define([M4_COMPUTE_HORZ_COEF_HH], [dnl
#ifdef _LINEAR
#ifdef CACHE_TABLE_OPT1
    HH=hfecHH(itau,ki,kmat)+hfecHH(itau,kj,kmat)
    HH=HH + xtau*(hfecHH(itau+1,ki,kmat)+hfecHH(itau+1,kj,kmat) - HH)
#else
    HH=hfecHH(ki,itau,kmat)+hfecHH(kj,itau,kmat)
    HH=HH + xtau*(hfecHH(ki,itau+1,kmat)+hfecHH(kj,itau+1,kmat) - HH)
#endif
#else
#ifdef CACHE_TABLE_OPT1
    HH=hfecHH(itau,ki,kmat)+hfecHH(itau,kj,kmat)
#else
    HH=hfecHH(ki,itau,kmat)+hfecHH(kj,itau,kmat)
#endif
[#endif]])
dnl
dnl !-----------------------------------------------------
dnl ! HD coefficient used in HH horizontal correlation
dnl !-----------------------------------------------------
define([M4_COMPUTE_HORZ_COEF_HD],[dnl
#ifdef _LINEAR B
#ifdef CACHE_TABLE_OPT1
    HR=HFECRRTT(1,itau,ki,kmat)+HFECRRTT(1,itau,kj,kmat)
    HR=HR + xtau*		&
	 &		(HFECRRTT(1,itau+1,ki,kmat)+	&
	 &		 HFECRRTT(1,itau+1,kj,kmat) - HR)
#else
    HR=hfecRR(ki,itau,kmat)+hfecRR(kj,itau,kmat)
    HR=HR + xtau*		&
	 &		(hfecRR(ki,itau+1,kmat)+hfecRR(kj,itau+1,kmat) - HR)
#endif
#else
			! a non-separable correlation function form
#ifdef CACHE_TABLE_OPT1
    HR=HFECRRTT(1,itau,ki,kmat)+HFECRRTT(1,itau,kj,kmat)
#else
    HR=hfecRR(ki,itau,kmat)+hfecRR(kj,itau,kmat)
#endif
[#endif]])
dnl
dnl
dnl !-----------------------------------------------------
dnl ! DD coefficient used in HH horizontal correlation
dnl !-----------------------------------------------------
define([M4_COMPUTE_HORZ_COEF_DD], [dnl
#ifdef _LINEAR
#ifdef CACHE_TABLE_OPT1
    RR=HFECRRTT(1,itau,ki,kmat)+HFECRRTT(1,itau,kj,kmat)
    RR=RR + xtau*			&
	 &(HFECRRTT(1,itau+1,ki,kmat)+	&
	 & HFECRRTT(1,itau+1,kj,kmat) - RR)
    TT=HFECRRTT(2,itau,ki,kmat)+HFECRRTT(2,itau,kj,kmat)
    TT=TT + xtau*			&
	 &(HFECRRTT(2,itau+1,ki,kmat)+	&
	 & HFECRRTT(2,itau+1,kj,kmat) - TT)
#else
    
    RR=hfecRR(ki,itau,kmat)+hfecRR(kj,itau,kmat)
    RR=RR + xtau*		&
	 &(hfecRR(ki,itau+1,kmat)+hfecRR(kj,itau+1,kmat) - RR)
    TT=hfecTT(ki,itau,kmat)+hfecTT(kj,itau,kmat)
    TT=TT + xtau*		&
	 &(hfecTT(ki,itau+1,kmat)+hfecTT(kj,itau+1,kmat) - TT)
    
#endif
#else
    ! a non-separable correlation function form
#ifdef CACHE_TABLE_OPT1
    RR=HFECRRTT(1,itau,ki,kmat)+HFECRRTT(1,itau,kj,kmat)
    TT=HFECRRTT(2,itau,ki,kmat)+HFECRRTT(2,itau,kj,kmat)
#else
    RR=hfecRR(ki,itau,kmat)+hfecRR(kj,itau,kmat)
    TT=hfecTT(ki,itau,kmat)+hfecTT(kj,itau,kmat)
#endif
[#endif]])



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
M4_UNPACK_QD(i)
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
M4_UNPACK_QD(i)
ifelse(OPERATION, ASY_FILL,[dnl
dnl    ! only do lone block
M4_COMPUTE_DH_ELEM()],     [dnl
dnl    ! do 2 blocks at once
M4_COMPUTE_DH_ELEM(1)
M4_COMPUTE_DH_ELEM(2)])])])])


dnl !-----------------------------------------------------
define([M4_COMPUTE_HH_ELEM], [dnl
       horcor = HH
       element = horcor*vercor])

dnl !-----------------------------------------------------
define([M4_COMPUTE_HD_ELEM], [dnl
       dot_rm = DOT(qri,qmj$1)
       horcor = HR * dot_rm
       element$1 = horcor * vercor])

dnl !-----------------------------------------------------
define([M4_COMPUTE_DH_ELEM], [dnl
       dot_mr = DOT(qmi$1,qrj)
       horcor = HR * dot_mr
       element$1 = horcor * vercor])

dnl !-----------------------------------------------------
define([M4_COMPUTE_DD_1_ELEM], [dnl
      dot_rm = DOT(qri,qmj)
      dot_mr = DOT(qmi,qrj)
      dot_mm = DOT(qmi,qmj)
      horcor = RR * dot_mm + TT * dot_rm * dot_mr
      element = horcor*vercor])

dnl !-----------------------------------------------------
define([M4_COMPUTE_DD_2X2_ELEM], [dnl
      dot_rm_1 = DOT(qri,qmj1)
      dot_mr_1 = DOT(qmi1,qrj)
      dot_mm = DOT(qmi1,qmj1)
      horcor = RR * dot_mm + TT * dot_rm_1 * dot_mr_1
      element_11 = horcor*vercor

      dot_rm_2 = DOT(qri,qmj2)
      dot_mm = DOT(qmi1,qmj2)
      horcor = RR * dot_mm + TT * dot_rm_2 * dot_mr_1
      element_12 = horcor*vercor

      dot_mr_2 = DOT(qmi2,qrj)
      dot_mm = DOT(qmi2,qmj1)
      horcor = RR * dot_mm + TT * dot_rm_1 * dot_mr_2
      element_21 = horcor*vercor

      dot_mm = DOT(qmi2,qmj2)
      horcor = RR * dot_mm + TT * dot_rm_2 * dot_mr_2
      element_22 = horcor*vercor])


dnl! -----------------------------------------------------
dnl! 
dnl!-----------------------------------------------------
define([M4_COMPUTE_ELEMENT_IJ],[dnl
    tau = DOT(qri,qrj)
    If (tau > Hcoslim) Then
M4_COMPUTE_ITAU
M4_COMPUTE_VERCOR
M4_COMPUTE_HORZ_COEF
M4_COMPUTE_ELEM_FROM_COEFS
  M4_USE_ELEMENT
ifelse(OPERATION,  SYM_FILL, [dnl
    Else
       corr(l) = 0],
       OPERATION, ASY_FILL, [dnl
    Else
       corr(i,j)=0])
    End If])

dnl! -----------------------------------------------------
dnl! 
dnl!-----------------------------------------------------
define([M4_COMPUTE_DIAG_ELEMENT], [dnl
ifdef([FILL],[dnl
      corr(l) = 1],[dnl
ifelse(BLOCK_TYPE, HH, [dnl
dnl ! HH case
      Y(:,j) = Y(:,j) + X(:,j)],[dnl
dnl ! DD case
M4_UNPACK_QR(j)
M4_UNPACK_QR(i)
      tau = DOT(qri,qrj)
      kj = kl(j)
      ki = kj

      If (tau > Hcoslim) Then
M4_COMPUTE_ITAU
M4_COMPUTE_VERCOR
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
kind_cov, ln, qr, ifelse($2, DD, [qd,]) kl, nvecs, x, y, weight],
       $3, SYM_FILL, [dnl
kind_cov, ln, qr, ifelse($2, DD, [qd,]) kl, corr, weight], [dnl
ifelse($3, ASY_FILL, [dnl
kind_cov, lni, qri, ifelse($2, DD, [qdi,], $2, DH, [qdi,]) kli, dnl
          lnj, qrj, ifelse($2, DD, [qdj,], $2, HD, [qdj,]) klj, corr, weight],
       $3, XCXPY,  [dnl
kind_cov, lni, qri, ifelse($2, DD, [qdi,], $2, DH, [qdi,]) kli, dnl
          lnj, qrj, ifelse($2, DD, [qdj,], $2, HD, [qdj,]) klj, nvecs, xi, yj, xj, yi, weight], [dnl
dnl RCXPY
kind_cov, lni, qri, ifelse($2, DD, [qdi,], $2, DH, [qdi,]) kli, dnl
          lnj, qrj, ifelse($2, DD, [qdj,], $2, HD, [qdj,]) klj, nvecs, x, y, weight])])])

dnl!-----------------------------------------------------
define([M4_ARGUMENT_DECLARATION],[dnl
     Integer, Intent(In) :: kind_cov
dnl ! Parameters
ifdef([SYMMETRIC],[dnl
     Integer, Intent(In)   :: ln
     Real,    Intent(In)   :: qr(3, ln)
ifelse(BLOCK_TYPE, DD, [dnl
     Real, Intent(In) :: qd(3,ln ifdef([FILL],,[,2]))])
     Integer, Intent(In)   :: kl(ln)],[dnl
dnl ! Asymmetric
     Integer, Intent(In)   :: lni
     Real,    Intent(In)   :: qri(3, lni)
ifelse(BLOCK_TYPE, HH,, BLOCK_TYPE, HD,, [
     Real, Intent(In) :: qdi(3,lni ifdef([FILL],,[,2]))])
     Integer, Intent(In)   :: kli(lni)

     Integer, Intent(In)   :: klj(lnj)
     Integer, Intent(In)   :: lnj
     Real,    Intent(In)   :: qrj(3, lnj)
ifelse(BLOCK_TYPE, HH,, BLOCK_TYPE, DH,, [dnl
     Real, Intent(In) :: qdj(3,lnj ifdef([FILL],,[,2]))])
     Integer, Intent(In)   :: klj(lnj)])
dnl rhs,lhs,block storage
ifdef([FILL], [dnl
     Real, Intent(Out) :: corr(ifdef([SYMMETRIC],[ln*(ln+1)/2],[lni,lnj]))],[dnl
     Integer, Intent(In) :: nvecs
dnl vectors x, y
ifelse(OPERATION, SCXPY, [dnl
ifelse(BLOCK_TYPE, HH, [
     Real, Intent(In)    :: X(nvecs, ln)
     Real, Intent(InOut) :: Y(nvecs, ln)],[
     Real, Intent(In)    :: X(nvecs, ln,2)
     Real, Intent(InOut) :: Y(nvecs, ln,2)])],
       OPERATION, RCXPY, [dnl
ifelse(BLOCK_TYPE, HH, [
     Real, Intent(In)    :: X(nvecs, lnj)
     Real, Intent(InOut) :: Y(nvecs, lni)],
       BLOCK_TYPE, DD, [
     Real, Intent(In)    :: X(nvecs, lnj, 2)
     Real, Intent(InOut) :: Y(nvecs, lni, 2)], [
ifelse(BLOCK_TYPE, HD,[
     Real, Intent(In)    :: X(nvecs, lnj, 2)
     Real, Intent(InOut) :: Y(nvecs, lni)], [
     Real, Intent(In)    :: X(nvecs, lnj)
     Real, Intent(InOut) :: Y(nvecs, lni, 2)])])],[
dnl XCXPY
ifelse(BLOCK_TYPE, HH, [
     Real, Intent(In)    :: XI(nvecs, lni)
     Real, Intent(InOut) :: YJ(nvecs, lnj)
     Real, Intent(In)    :: XJ(nvecs, lnj)
     Real, Intent(InOut) :: YI(nvecs, lni)],
       BLOCK_TYPE, DD, [
     Real, Intent(In)    :: XI(nvecs, lni, 2)
     Real, Intent(InOut) :: YJ(nvecs, lnj, 2)
     Real, Intent(In)    :: XJ(nvecs, lnj, 2)
     Real, Intent(InOut) :: YI(nvecs, lni, 2)],[
dnl HD
     Real, Intent(In)    :: XI(nvecs, lni)
     Real, Intent(InOut) :: YJ(nvecs, lnj, 2)
     Real, Intent(In)    :: XJ(nvecs, lnj, 2)
     Real, Intent(InOut) :: YI(nvecs, lni)])])])

     Real, Optional, Intent(in) :: weight
])




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
 
   Subroutine M4_ROUTINE_NAME($*)(M4_ARG_LIST($*))
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
      Integer :: itau
      Real    :: HH, HR, RR, TT
      Real    :: horcor, vercor
      Real    :: dot_rr
      Real    :: dot_mr, dot_mr_1, dot_mr_2
      Real    :: dot_rm, dot_rm_1, dot_rm_2
      Real    :: dot_mm
      Real    :: element, element1, element2
      Real    :: element_11, element_12, element_21, element_22
      Integer, External :: my_pe
      Real    :: weight_

!======================================================================

ifelse($3,SCXPY,[dnl
      If (ln <= 0) Return],
       $3, SYM_FILL,[dnl
      If (ln <= 0) Return],[dnl
      If (lni <=0 .or. lnj <= 0) Return])
ifelse($1,1,,[dnl
      If (nvecs <= 0) Return])

  weight_=1.
  if(present(weight)) weight_=weight

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


End Module m_mvcorF_base_bmop



