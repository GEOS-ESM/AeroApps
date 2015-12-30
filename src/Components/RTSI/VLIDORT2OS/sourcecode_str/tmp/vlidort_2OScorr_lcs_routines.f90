! ###############################################################
! #                                                             #
! #                    THE VECTOR LIDORT MODEL                  #
! #                                                             #
! #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
! #   -      --         -        -        -         -           #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Author :      Robert. J. D. Spurr                          #
! #                                                             #
! #  Address :     RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                Tel: (617) 492 1183                          #
! #                                                             #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT, 2.4RTC,  #
! #                   2.5, 2.6, 2.7                             #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #  Release Date :   August 2014    (2.7)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
! #                                                             #
! #       SURFACE-LEAVING / BRDF-SCALING      (2.7)             #
! #       TAYLOR Series / OMP THREADSAFE      (2.7)             #
! #                                                             #
! ###############################################################

!    #####################################################
!    #                                                   #
!    #   This Version of VLIDORT comes with a GNU-style  #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###########################################################
! #                                                         #
! #  2OS History :                                          #
! #                                                         #
! #     Mark 1: 2007 for OCO-Mk1, publication               #
! #     Mark 2: 2009, with BRDFs                            #
! #     Mark 3: 2013, Re-engineered Model:                  #
! #              * Including Multiple Geometry              #
! #              * Merging with Stand-alone FO code         #
! #              * New linearization for Bulk properties    #
! #     Mark 4: 2014, integration with VLIDORT              #
! #                                                         #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! # PUBLIC Subroutines in this Module                       #
! #                                                         #
! #                 Calculate_SecondOrder_LCSPlus           #
! #                 Calculate_FirstOrder_LCSPlus            #
! #                 Add_Fourier_Component_LCSPlus           #
! #                 Calculate_Multipliers_LCPlus            #
! #                 Set_avsecant_LCPlus                     #
! #                                                         #
! ###########################################################

module vlidort_2OScorr_lcs_routines

   use VLIDORT_PARS, only : fpk, zero, one, two, three, four,       &
                            six, half, quarter, deg_to_rad,         &
                            BIGEXP, TAYLOR_SMALL, Smallnum6,        &
                            MAXLAYERS, MAXMOMENTS, MAXSTREAMS,      &
                            MAXSTREAMS_P2, MAX_ATMOSWFS, MAX_SURFACEWFS

   use vlidort_2OScorr_utilities, only : Exptrans, Exptrans_L,        &
                                         Make_Trans23, Make_Trans23_P

   implicit none

public  :: Calculate_SecondOrder_LCSPlus,  &
           Calculate_FirstOrder_LCSPlus,   &
           Add_Fourier_Component_LCSPlus,  &
           Calculate_Multipliers_LCPlus,   &
           Set_avsecant_LCPlus


contains

subroutine Calculate_SecondOrder_LCSPlus &
       ( do_LC_Jacobians, do_LS_Jacobians,m,n,nstreams,nstokes,npars,nspars, & ! Input Control
         qweights,omega,xv,xa,  LC_omega,LC_xa,               & ! Input Basics
         O2_AVTrans,            LC_O2_AVTrans,                & ! Input R2 Multiplier
         O2_QVMult_d,           LC_O2_QVMult_d,               & ! Input R2 Multipliers
         O2_QAMult_d,           LC_O2_QAMult_d,               & ! Input R2 Multipliers
         O2_QAVMult_du,         LC_O2_QAVMult_du,             & ! Input R2 Multipliers
         O2_QAVMult_dd,         LC_O2_QAVMult_dd,             & ! Input R2 Multipliers 
         Ptc,Pts,Prc,Prs,       LC_Ptc,LC_Pts,LC_Prc,LC_Prs,  & ! Input scattering matrices
         R1c,R1s,R1cscal,LC_R1c,LC_R1s,LC_R1cscal,LS_R1c,LS_R1s,LS_R1cscal,  & ! Input R1 sources
         R2c,R2s,R2cscal,LC_R2c,LC_R2s,LC_R2cscal,LS_R2c,LS_R2s,LS_R2cscal )   ! Output R2 reflectances

!  Purpose: Update the second-order Reflectance field in Layer n.
!           Also updates the Bulk-property Jacobians

   implicit none

!  Inputs
!  ======

!  1. Standard
!  -----------

!  Control (Fourier #, Layer #, NSTREAMS, NSTOKES)

   integer  , intent(in) :: m, n, nstreams, nstokes

!  General inputs (quadrature weights, ssa, secants)

   real(fpk), intent(in) :: qweights(MAXSTREAMS)
   real(fpk), intent(in) :: omega, xv, xa

!  source phase matrices

   real(fpk), intent(in) :: Ptc(2,MAXSTREAMS,4,4),Pts(2,MAXSTREAMS,4,4)
   real(fpk), intent(in) :: Prc(2,MAXSTREAMS,4,4),Prs(2,MAXSTREAMS,4,4)

!  Second-order Transmittances and multipliers

   real(fpk), intent(in) :: O2_AVTrans
   real(fpk), intent(in) :: O2_QVMult_d(MAXSTREAMS)
   real(fpk), intent(in) :: O2_QAMult_d(MAXSTREAMS)
   real(fpk), intent(in) :: O2_QAVMult_dd(MAXSTREAMS)
   real(fpk), intent(in) :: O2_QAVMult_du(MAXSTREAMS)

!  2. Linearized
!  -------------

!  Linearization control

   logical  , intent(in) :: do_LC_Jacobians, do_LS_Jacobians
   integer  , intent(in) :: npars, nspars

!  Linearizated General inputs (secants, ssa)

   real(fpk), intent(in) :: LC_omega(MAX_ATMOSWFS), LC_xa(MAX_ATMOSWFS)

!  Linearized source phase matrices

   real(fpk), intent(in) :: LC_Ptc(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LC_Pts(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LC_Prc(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LC_Prs(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)

!  Linearized Second-order Transmittances and multipliers

   real(fpk), intent(in) :: LC_O2_AVTrans(MAX_ATMOSWFS)
   real(fpk), intent(in) :: LC_O2_QVMult_d(MAXSTREAMS,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LC_O2_QAMult_d(MAXSTREAMS,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LC_O2_QAVMult_dd(MAXSTREAMS,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LC_O2_QAVMult_du(MAXSTREAMS,MAX_ATMOSWFS)

!  Outputs
!  =======

!  Modified first order inputs (standard & Linearized)

   real(fpk), intent(inout) :: R1c(2,MAXSTREAMS,4,4)
   real(fpk), intent(inout) :: R1s(2,MAXSTREAMS,4,4)
   real(fpk), intent(inout) :: R1cscal(2,MAXSTREAMS)
   real(fpk), intent(inout) :: LC_R1c(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)
   real(fpk), intent(inout) :: LC_R1s(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)
   real(fpk), intent(inout) :: LC_R1cscal(2,MAXSTREAMS,MAX_ATMOSWFS)
   real(fpk), intent(inout) :: LS_R1c(2,MAXSTREAMS,4,4,MAX_SURFACEWFS)
   real(fpk), intent(inout) :: LS_R1s(2,MAXSTREAMS,4,4,MAX_SURFACEWFS)
   real(fpk), intent(inout) :: LS_R1cscal(2,MAXSTREAMS,MAX_SURFACEWFS)

!  Modified second order outputs (standard & Linearized)

   real(fpk), intent(inout) :: R2c(4),R2s(4),R2cscal
   real(fpk), intent(inout) :: LC_R2c(4,MAX_ATMOSWFS),  LC_R2s(4,MAX_ATMOSWFS),  LC_R2cscal(MAX_ATMOSWFS)
   real(fpk), intent(inout) :: LS_R2c(4,MAX_SURFACEWFS),LS_R2s(4,MAX_SURFACEWFS),LS_R2cscal(MAX_SURFACEWFS)

!  local variables
!  ===============

   integer   :: i, k, l, nd, p

   real(fpk) :: dqw, pgv, pha, dglq, dhlq, hw
   real(fpk) :: L_dqw(MAX_ATMOSWFS),  L_pgv(MAX_ATMOSWFS), L_pha(MAX_ATMOSWFS)
   real(fpk) :: L_dglq(MAX_ATMOSWFS), L_dhlq(MAX_ATMOSWFS), L_hw(MAX_ATMOSWFS)

   real(fpk) :: S1(4),S2(4),S3(4),S4(4),S1scal,S2scal
   real(fpk) :: L_S1(4,MAX_ATMOSWFS),   L_S2(4,MAX_ATMOSWFS)  ,  LS_S1(4,MAX_SURFACEWFS),   LS_S2(4,MAX_SURFACEWFS)  
   real(fpk) :: L_S3(4,MAX_ATMOSWFS),   L_S4(4,MAX_ATMOSWFS)  ,  LS_S3(4,MAX_SURFACEWFS),   LS_S4(4,MAX_SURFACEWFS)
   real(fpk) :: L_S1scal(MAX_ATMOSWFS), L_S2scal(MAX_ATMOSWFS),  LS_S1scal(MAX_SURFACEWFS), LS_S2scal(MAX_SURFACEWFS)

   real(fpk) :: V1(4),V2(4,4),V4(4),V6(4,4),V1scal,V2scal
   real(fpk) :: L_V1(4,MAX_ATMOSWFS),   L_V2(4,4,MAX_ATMOSWFS),  LS_V1(4,MAX_SURFACEWFS),   LS_V2(4,4,MAX_SURFACEWFS)
   real(fpk) :: L_V4(4,MAX_ATMOSWFS),   L_V6(4,4,MAX_ATMOSWFS),  LS_V4(4,MAX_SURFACEWFS),   LS_V6(4,4,MAX_SURFACEWFS)
   real(fpk) :: L_V1scal(MAX_ATMOSWFS), L_V2scal(MAX_ATMOSWFS),  LS_V1scal(MAX_SURFACEWFS), LS_V2scal(MAX_SURFACEWFS)

!  initial section
!  ---------------

!  Zero the  S-matrices

   S1       = zero ; S2       = zero
   S1scal   = zero ; S2scal   = zero
   S3       = zero ; S4       = zero

   L_S1     = zero ; L_S2     = zero
   L_S1scal = zero ; L_S2scal = zero
   L_S3     = zero ; L_S4     = zero

   LS_S1     = zero ; LS_S2     = zero
   LS_S1scal = zero ; LS_S2scal = zero
   LS_S3     = zero ; LS_S4     = zero

   nd     = n    !  Dummy variable, debug only

   hw = Omega * quarter
   if ( do_LC_Jacobians ) L_hw(1:npars) = LC_omega(1:npars) * quarter 

!  Enter the k-loop
!  ================

   do k = 1, nstreams

!  factors in front of Multipliers

      dglq = hw * O2_QAVMult_du(k) ; pgv = xv * O2_QVMult_d(k)
      dhlq = hw * O2_QAVMult_dd(k) ; pha = xa * O2_QAMult_d(k) 

      if ( do_LC_Jacobians ) then
        do p = 1, npars
          L_dglq(p) = L_hw(p) * O2_QAVMult_du(k) + hw * LC_O2_QAVMult_du(k,p) 
          L_dhlq(p) = L_hw(p) * O2_QAVMult_dd(k) + hw * LC_O2_QAVMult_dd(k,p) 
          L_pgv(p)  = xv * LC_O2_QVMult_d(k,p)
          L_pha(p)  = xa * LC_O2_QAMult_d(k,p) + LC_xa(p) * O2_QAMult_d(k)
        enddo
      endif

!  V matrices
!  ----------

!  Part 1: The I and Q components

      do l = 1, 2
        V1(l)     = pgv * R1c(2,k,l,1)   + dglq * Prc(2,k,l,1)
        V2(1:2,l) = pha * R1c(1,k,1:2,l) + dhlq * Prc(1,k,1:2,l)
        if ( do_LC_Jacobians ) then
          do p = 1, npars
            L_V1(l,p)     =  L_pgv(p)  * R1c(2,k,l,1)   +  pgv * LC_R1c(2,k,l,1,p)  &
                           + L_dglq(p) * Prc(2,k,l,1)   + dglq * LC_Prc(2,k,l,1,p)
            L_V2(1:2,l,p) =  L_pha(p)  * R1c(1,k,1:2,l) +  pha * LC_R1c(1,k,1:2,l,p) &
                           + L_dhlq(p) * Prc(1,k,1:2,l) + dhlq * LC_Prc(1,k,1:2,l,p)
          enddo
        endif
        if ( do_LS_Jacobians ) then
          do p = 1, nspars
            LS_V1(l,p)     = pgv * LS_R1c(2,k,l,1,p)
            LS_V2(1:2,l,p) = pha * LS_R1c(1,k,1:2,l,p)
          enddo
        endif
      enddo

      if (m .gt. 0) then
        do l = 3, 4
          V4(l)     = pgv * R1s(2,k,l,1)   + dglq * Prs(2,k,l,1)
          V6(1:2,l) = pha * R1s(1,k,1:2,l) + dhlq * Prs(1,k,1:2,l)
          if ( do_LC_Jacobians ) then
            do p = 1, npars
              L_V4(l,p)     =   L_pgv(p) * R1s(2,k,l,1) +  pgv * LC_R1s(2,k,l,1,p)  &
                             + L_dglq(p) * Prs(2,k,l,1) + dglq * LC_Prs(2,k,l,1,p)
              L_V6(1:2,l,p) =  L_pha(p)  * R1s(1,k,1:2,l) +  pha * LC_R1s(1,k,1:2,l,p) &
                             + L_dhlq(p) * Prs(1,k,1:2,l) + dhlq * LC_Prs(1,k,1:2,l,p)
            enddo
          endif
          if ( do_LS_Jacobians ) then
            do p = 1, nspars
              LS_V4(l,p)     = pgv * LS_R1s(2,k,l,1,p)
              LS_V6(1:2,l,p) = pha * LS_R1s(1,k,1:2,l,p)
            enddo
          endif
        enddo
      endif

!  Part 2: the U component

      if (m .gt. 0) then
        do i = 3, nstokes
          V2(i,3:4) = pha * R1c(1,k,i,3:4) + dhlq * Prc(1,k,i,3:4)
          V6(i,1:2) = pha * R1s(1,k,i,1:2) + dhlq * Prs(1,k,i,1:2)
          if ( do_LC_Jacobians ) then
            do p = 1, npars
              L_V2(i,3:4,p)  =  L_pha(p) * R1c(1,k,i,3:4) +  pha * LC_R1c(1,k,i,3:4,p)  &
                             + L_dhlq(p) * Prc(1,k,i,3:4) + dhlq * LC_Prc(1,k,i,3:4,p)
              L_V6(i,1:2,p) =  L_pha(p)  * R1s(1,k,i,1:2) +  pha * LC_R1s(1,k,i,1:2,p) &
                             + L_dhlq(p) * Prs(1,k,i,1:2) + dhlq * LC_Prs(1,k,i,1:2,p)
            enddo
          endif
          if ( do_LS_Jacobians ) then
            do p = 1, nspars
              LS_V2(i,3:4,p) =  pha * LS_R1c(1,k,i,3:4,p)
              LS_V6(i,1:2,p) =  pha * LS_R1s(1,k,i,1:2,p)
            enddo
          endif
        enddo
      endif

!  Part 3: the scalar calculation

      V1scal = pgv * R1cscal(2,k) + dglq * Prc(2,k,1,1)
      V2scal = pha * R1cscal(1,k) + dhlq * Prc(1,k,1,1)
      if ( do_LC_Jacobians ) then
        do p = 1, npars
          L_V1scal(p) =  L_pgv(p) * R1cscal(2,k) +  pgv * LC_R1cscal(2,k,p) &
                      + L_dglq(p) * Prc(2,k,1,1) + dglq * LC_Prc(2,k,1,1,p)
          L_V2scal(p) =  L_pha(p) * R1cscal(1,k) +  pha * LC_R1cscal(1,k,p) &
                      + L_dhlq(p) * Prc(1,k,1,1) + dhlq * LC_Prc(1,k,1,1,p)
        enddo
      endif
      if ( do_LS_Jacobians ) then
        do p = 1, nspars
          LS_V1scal(p) = pgv * LS_R1cscal(2,k,p)
          LS_V2scal(p) = pha * LS_R1cscal(1,k,p)
        enddo
      endif

!  debug

!      write(94,'(2i3,1p5e20.12)')m,n,V4(3),V4(4),V6(3,1)
!      if (m.eq.3.and.n.eq.60)pause'gronk'

!  S Matrices
!  ----------

      dqw = Omega * qweights(k)
      if ( do_LC_Jacobians ) then
        L_dqw(1:npars) = LC_Omega(1:npars) * qweights(k)
      endif

!  Part 1: I and Q components

      do i = 1, 2
         do l = 1, 2
           S1(i) = S1(i) + dqw * Ptc(1,k,i,l) * V1(l)
           S2(i) = S2(i) + dqw * Ptc(2,k,l,1) * V2(i,l)
         enddo
         if ( do_LC_Jacobians ) then
           do p = 1, npars
             do l = 1, 2
               L_S1(i,p) = L_S1(i,p) + L_dqw(p) * Ptc(1,k,i,l)      * V1(l)    + &
                                            dqw * LC_Ptc(1,k,i,l,p) * V1(l)    + &
                                            dqw * Ptc(1,k,i,l)      * L_V1(l,p)
               L_S2(i,p) = L_S2(i,p) + L_dqw(p) * Ptc(2,k,l,1)      * V2(i,l)    + &
                                            dqw * LC_Ptc(2,k,l,1,p) * V2(i,l)    + &
                                            dqw * Ptc(2,k,l,1)      * L_V2(i,l,p)
             enddo
           enddo
         endif
         if ( do_LS_Jacobians ) then
           do p = 1, nspars
             do l = 1, 2
               LS_S1(i,p) = LS_S1(i,p) + dqw * Ptc(1,k,i,l) * LS_V1(l,p)
               LS_S2(i,p) = LS_S2(i,p) + dqw * Ptc(2,k,l,1) * LS_V2(i,l,p)
             enddo
           enddo
         endif
         if (m .gt. 0) then
           do l = 3,4
             S3(i) = S3(i) + dqw * Pts(1,k,i,l) * V4(l)
             S4(i) = S4(i) + dqw * Pts(2,k,l,1) * V6(i,l)
           enddo
           if ( do_LC_Jacobians ) then
             do p = 1, npars
               do l = 3, 4
                 L_S3(i,p) = L_S3(i,p) + L_dqw(p) * Pts(1,k,i,l)      * V4(l)    + &
                                              dqw * LC_Pts(1,k,i,l,p) * V4(l)    + &
                                              dqw * Pts(1,k,i,l)      * L_V4(l,p)
                 L_S4(i,p) = L_S4(i,p) + L_dqw(p) * Pts(2,k,l,1)      * V6(i,l)    + &
                                              dqw * LC_Pts(2,k,l,1,p) * V6(i,l)    + &
                                              dqw * Pts(2,k,l,1)      * L_V6(i,l,p)
               enddo
             enddo
           endif
           if ( do_LS_Jacobians ) then
             do p = 1, nspars
               do l = 3, 4
                 LS_S3(i,p) = LS_S3(i,p) + dqw * Pts(1,k,i,l) * LS_V4(l,p)
                 LS_S4(i,p) = LS_S4(i,p) + dqw * Pts(2,k,l,1) * LS_V6(i,l,p)
               enddo
             enddo
           endif
         endif
      enddo

!  Part 2: U component

      do i = 3, nstokes
         if (m .gt. 0) then
           do l = 3, 4
             S1(i) = S1(i) + dqw * Ptc(1,k,i,l) * V4(l)
             S2(i) = S2(i) + dqw * Pts(2,k,l,1) * V2(i,l)
           enddo
           if ( do_LC_Jacobians ) then
             do p = 1, npars
               do l = 3, 4
                 L_S1(i,p) = L_S1(i,p) + L_dqw(p) * Ptc(1,k,i,l)      * V4(l)    + &
                                              dqw * LC_Ptc(1,k,i,l,p) * V4(l)    + &
                                              dqw * Ptc(1,k,i,l)      * L_V4(l,p)
                 L_S2(i,p) = L_S2(i,p) + L_dqw(p) * Pts(2,k,l,1)      * V2(i,l)    + &
                                              dqw * LC_Pts(2,k,l,1,p) * V2(i,l)    + &
                                              dqw * Pts(2,k,l,1)      * L_V2(i,l,p)
               enddo
             enddo
           endif
           if ( do_LS_Jacobians ) then
             do p = 1, nspars
               do l = 3, 4
                 LS_S1(i,p) = LS_S1(i,p) + dqw * Ptc(1,k,i,l) * LS_V4(l,p)
                 LS_S2(i,p) = LS_S2(i,p) + dqw * Pts(2,k,l,1) * LS_V2(i,l,p)
               enddo
             enddo
           endif
           do l = 1, 2
             S3(i) = S3(i) + dqw * Pts(1,k,i,l) * V1(l)
             S4(i) = S4(i) + dqw * Ptc(2,k,l,1) * V6(i,l)
           enddo
           if ( do_LC_Jacobians ) then
             do p = 1, npars
               do l = 1, 2
                 L_S3(i,p) = L_S3(i,p) + L_dqw(p) * Pts(1,k,i,l)      * V1(l)    + &
                                              dqw * LC_Pts(1,k,i,l,p) * V1(l)    + &
                                              dqw * Pts(1,k,i,l)      * L_V1(l,p)
                 L_S4(i,p) = L_S4(i,p) + L_dqw(p) * Ptc(2,k,l,1)      * V6(i,l)    + &
                                              dqw * LC_Ptc(2,k,l,1,p) * V6(i,l)    + &
                                              dqw * Ptc(2,k,l,1)      * L_V6(i,l,p)
               enddo
             enddo
           endif
           if ( do_LS_Jacobians ) then
             do p = 1, nspars
               do l = 1, 2
                 LS_S3(i,p) = LS_S3(i,p) + dqw * Pts(1,k,i,l) * LS_V1(l,p)
                 LS_S4(i,p) = LS_S4(i,p) + dqw * Ptc(2,k,l,1) * LS_V6(i,l,p)
               enddo
             enddo
           endif
         endif
      enddo

!  Part 3: scalar

      S1scal = S1scal + dqw * Ptc(1,k,1,1) * V1scal
      S2scal = S2scal + dqw * Ptc(2,k,1,1) * V2scal
      if ( do_LC_Jacobians ) then
        do p = 1, npars
          L_S1scal = L_S1scal + L_dqw(p) * Ptc(1,k,1,1) * V1scal + &
                                  dqw    * ( LC_Ptc(1,k,1,1,p) * V1scal + Ptc(1,k,1,1) * L_V1scal(p) )
          L_S2scal = L_S2scal + L_dqw(p) * Ptc(2,k,1,1) * V2scal + &
                                  dqw    * ( LC_Ptc(2,k,1,1,p) * V2scal + Ptc(2,k,1,1) * L_V2scal(p) )
        enddo
      endif
      if ( do_LS_Jacobians ) then
        do p = 1, nspars
          LS_S1scal = LS_S1scal + dqw * Ptc(1,k,1,1) * LS_V1scal(p)
          LS_S2scal = LS_S2scal + dqw * Ptc(2,k,1,1) * LS_V2scal(p)
        enddo
      endif

!  end k loop

   enddo

!  Recursion (Using the solar/LOS transmittance multiplier)
!  =========

!  Update the R2c/R2s terms, with the transmittance

   do i = 1, 2
      if ( do_LC_Jacobians ) then
         LC_R2c(i,1:npars) = R2c(i) * LC_O2_AVTRANS(1:npars) + LC_R2c(i,1:npars)*O2_AVTRANS
      endif
      if ( do_LS_Jacobians ) then
         LS_R2c(i,1:nspars) = LS_R2c(i,1:nspars)*O2_AVTRANS
      endif
      R2c(i) = R2c(i)*O2_AVTRANS
   enddo

   if (nstokes .eq. 3) then
      do i = 3, nstokes
         if ( do_LC_Jacobians ) then
            LC_R2s(i,1:npars) = R2s(i) * LC_O2_AVTRANS(1:npars) + LC_R2s(i,1:npars)*O2_AVTRANS
         endif
         if ( do_LS_Jacobians ) then
            LS_R2s(i,1:nspars) = LS_R2s(i,1:nspars)*O2_AVTRANS
         endif
         R2s(i) = R2s(i)*O2_AVTRANS
      enddo
   endif

! update the scalar corrections

   if ( do_LC_Jacobians ) then
      LC_R2cscal(1:npars) = R2cscal * LC_O2_AVTRANS(1:npars) + LC_R2cscal(1:npars)*O2_AVTRANS
   endif
   if ( do_LS_Jacobians ) then
      LS_R2cscal(1:nspars) = LS_R2cscal(1:nspars)*O2_AVTRANS
   endif
   R2cscal = R2cscal*O2_AVTRANS

   do i = 1, 2
      if ( do_LC_Jacobians ) then
         LC_R2c(i,1:npars) = LC_R2c(i,1:npars) + half * ( L_S1(i,1:npars) + L_S2(i,1:npars) )
      endif
      if ( do_LS_Jacobians ) then
         LS_R2c(i,1:nspars) = LS_R2c(i,1:nspars) + half * ( LS_S1(i,1:nspars) + LS_S2(i,1:nspars) )
      endif
      R2c(i) = R2c(i) + half * ( S1(i) + S2(i) )
   enddo

   R2cscal = R2cscal + half * ( S1scal + S2scal )
   if ( do_LC_Jacobians ) then
      LC_R2cscal(1:npars) = LC_R2cscal(1:npars) + half * ( L_S1scal(1:npars) + L_S2scal(1:npars) )
   endif
   if ( do_LS_Jacobians ) then
      LS_R2cscal(1:nspars) = LS_R2cscal(1:nspars) + half * ( LS_S1scal(1:nspars) + LS_S2scal(1:nspars) )
   endif

   if (m .gt. 0) then
      do i = 1, 2
         R2c(i) = R2c(i) + half * ( S3(i) - S4(i) )
         if ( do_LC_Jacobians ) then
            LC_R2c(i,1:npars) = LC_R2c(i,1:npars) + half * ( L_S3(i,1:npars) - L_S4(i,1:npars) )
         endif
         if ( do_LS_Jacobians ) then
            LS_R2c(i,1:nspars) = LS_R2c(i,1:nspars) + half * ( LS_S3(i,1:nspars) - LS_S4(i,1:nspars) )
         endif
      enddo
      if (nstokes .eq. 3) then
         do i = 3, nstokes
            R2s(i) = R2s(i) + half  * ( S1(i) + S2(i) - S3(i) + S4(i) )
            if ( do_LC_Jacobians ) then
               LC_R2s(i,1:npars) = LC_R2s(i,1:npars) + half * &
                 ( L_S1(i,1:npars) + L_S2(i,1:npars) - L_S3(i,1:npars) + L_S4(i,1:npars) )
            endif
            if ( do_LS_Jacobians ) then
               LS_R2s(i,1:nspars) = LS_R2s(i,1:nspars) + half * &
                 ( LS_S1(i,1:nspars) + LS_S2(i,1:nspars) - LS_S3(i,1:nspars) + LS_S4(i,1:nspars) )
            endif
         enddo
      endif
   endif

!      write(94,'(2i3,1p5e20.12)')m,n,R2s(1),R2s(2),R2s(3)
!      if (m.eq.3.and.n.eq.60)pause'gronk'

!  Finish

   return
end subroutine Calculate_SecondOrder_LCSPlus

subroutine Calculate_FirstOrder_LCSPlus &
        ( do_LC_Jacobians, do_LS_Jacobians,m,nstreams,npars,nspars, & ! Input Control
          Prc,Prs,          LC_Prc,LC_Prs,             & ! Scattering input
          O1_QVTrans,       LC_O1_QVTrans,             & ! First-order multiplier input
          O1_QVMult,        LC_O1_QVMult,              & ! First-order multiplier input
          O1_QATrans,       LC_O1_QATrans,             & ! First-order multiplier input
          O1_QAMult,        LC_O1_QAMult,              & ! First-order multiplier input
          R1c,R1s,R1cscal,LC_R1c,LC_R1s,LC_R1cscal,LS_R1c,LS_R1s,LS_R1cscal ) ! First-order output

!  Purpose: Update the First-Order Reflectance for Layer n
!           Also updates the Bulk-property Jacobians

   implicit none

!  Inputs
!  ======

!  1. Standard
!  -----------

!  Control (Fourier #, NSTREAMS, NSTOKES)

   integer  , intent(in) :: m, nstreams

!  source phase matrices

   real(fpk), intent(in) :: Prc(2,MAXSTREAMS,4,4),Prs(2,MAXSTREAMS,4,4)

!  First-order transmittances/multipliers

   real(fpk), intent(in) :: O1_QVTrans(MAXSTREAMS), O1_QATrans(MAXSTREAMS)
   real(fpk), intent(in) :: O1_QVMult (MAXSTREAMS), O1_QAMult (MAXSTREAMS)

!  2. Linearized
!  -------------

!  Linearization control

   logical  , intent(in) :: do_LC_Jacobians, do_LS_Jacobians
   integer  , intent(in) :: npars, nspars

!  source phase matrices

   real(fpk), intent(in) :: LC_Prc(2,MAXSTREAMS,4,4,MAX_ATMOSWFS),LC_Prs(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)

!  First-order transmittances/multipliers

   real(fpk), intent(in) :: LC_O1_QVTrans(MAXSTREAMS,MAX_ATMOSWFS), LC_O1_QATrans(MAXSTREAMS,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LC_O1_QVMult (MAXSTREAMS,MAX_ATMOSWFS), LC_O1_QAMult (MAXSTREAMS,MAX_ATMOSWFS)

!  Modified outputs
!  ================

   real(fpk), intent(inout) :: R1c(2,MAXSTREAMS,4,4)
   real(fpk), intent(inout) :: R1s(2,MAXSTREAMS,4,4)
   real(fpk), intent(inout) :: R1cscal(2,MAXSTREAMS)

   real(fpk), intent(inout) :: LC_R1c(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)
   real(fpk), intent(inout) :: LC_R1s(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)
   real(fpk), intent(inout) :: LC_R1cscal(2,MAXSTREAMS,MAX_ATMOSWFS)

   real(fpk), intent(inout) :: LS_R1c(2,MAXSTREAMS,4,4,MAX_SURFACEWFS)
   real(fpk), intent(inout) :: LS_R1s(2,MAXSTREAMS,4,4,MAX_SURFACEWFS)
   real(fpk), intent(inout) :: LS_R1cscal(2,MAXSTREAMS,MAX_SURFACEWFS)

!  local variables
!  ===============

   integer :: j,i,k1,k2,p

!  viewing direction + all j-quadrature directions
!  -----------------------------------------------

!     update the R1 matrices themselves

   do j = 1, nstreams
!  R1scal updates.......
      if ( do_LC_Jacobians ) then
         do p = 1, npars
            LC_R1cscal(1,j,p) =  LC_O1_QVTRANS(j,p) *    R1cscal(1,j)    &
                               +    O1_QVTRANS(j)   * LC_R1cscal(1,j,p)  &
                               + LC_O1_QVMULT(j,p)  *    Prc(1,j,1,1)    &
                               +    O1_QVMULT(j)    * LC_Prc(1,j,1,1,p)
         enddo
      endif
      if ( do_LS_Jacobians ) then
         do p = 1, nspars
            LS_R1cscal(1,j,p) =  O1_QVTRANS(j)   * LS_R1cscal(1,j,p)
         enddo
      endif
      R1cscal(1,j) = O1_QVTrans(j) * R1cscal(1,j) + O1_QVMult(j)  * Prc(1,j,1,1)
!  R1c and R1s...............
      do k1 = 1, 4
         do k2 = 1, 4
            if ( do_LC_Jacobians ) then
               do p = 1, npars
                 LC_R1c(1,j,k1,k2,p) =  LC_O1_QVTRANS(j,p) *    R1c(1,j,k1,k2)   &
                                      +    O1_QVTRANS(j)   * LC_R1c(1,j,k1,k2,p) &
                                      + LC_O1_QVMULT(j,p)  *    Prc(1,j,k1,k2)   &
                                      +    O1_QVMULT(j)    * LC_Prc(1,j,k1,k2,p) 
               enddo
            endif
            if ( do_LS_Jacobians ) then
               do p = 1, nspars
                 LS_R1c(1,j,k1,k2,p) =  O1_QVTRANS(j)   * LS_R1c(1,j,k1,k2,p)
               enddo
            endif
            R1c(1,j,k1,k2) = O1_QVTrans(j) * R1c(1,j,k1,k2) + O1_QVMult(j)  * Prc(1,j,k1,k2)
            if (m .gt. 0) then
               if ( do_LC_Jacobians ) then
                  do p = 1, npars
                    LC_R1s(1,j,k1,k2,p) =  LC_O1_QVTRANS(j,p) *    R1s(1,j,k1,k2)   &
                                         +    O1_QVTRANS(j)   * LC_R1s(1,j,k1,k2,p) &
                                         + LC_O1_QVMULT(j,p)  *    Prs(1,j,k1,k2)   &
                                         +    O1_QVMULT(j)    * LC_Prs(1,j,k1,k2,p) 
                  enddo
               endif
               if ( do_LS_Jacobians ) then
                  do p = 1, nspars
                    LS_R1s(1,j,k1,k2,p) = O1_QVTRANS(j)   * LS_R1s(1,j,k1,k2,p)
                  enddo
               endif
               R1s(1,j,k1,k2) = O1_QVTrans(j) * R1s(1,j,k1,k2) + O1_QVMult(j) * Prs(1,j,k1,k2)
            endif
         enddo
      enddo
   enddo

!  solar direction
!  ---------------

   do i = 1, nstreams
!  R1scal updates.......
      if ( do_LC_Jacobians ) then
         do p = 1, npars
            LC_R1cscal(2,i,p) =  LC_O1_QATrans(i,p) *    R1cscal(2,i)   &
                               +    O1_QATrans(i)   * LC_R1cscal(2,i,p) &
                               + LC_O1_QAMult(i,p)  *    Prc(2,i,1,1)   &
                               +    O1_QAMult(i)    * LC_Prc(2,i,1,1,p)

         enddo
      endif
      if ( do_LS_Jacobians ) then
         do p = 1, nspars
            LS_R1cscal(2,i,p) =  O1_QATrans(i)   * LS_R1cscal(2,i,p)
         enddo
      endif
      R1cscal(2,i) = O1_QATrans(i) * R1cscal(2,i) + O1_QAMult(i)  * Prc(2,i,1,1)
!  R1c and R1s...............
      do k1 = 1, 4
         do k2 = 1, 4
            if ( do_LC_Jacobians ) then
               do p = 1, npars
                  LC_R1c(2,i,k1,k2,p) =  LC_O1_QATrans(i,p) *    R1c(2,i,k1,k2)   &
                                       +    O1_QATrans(i)   * LC_R1c(2,i,k1,k2,p) &
                                       + LC_O1_QAMult(i,p)  *    Prc(2,i,k1,k2)   &
                                       +    O1_QAMult(i)    * LC_Prc(2,i,k1,k2,p)
               enddo
            endif
            if ( do_LS_Jacobians ) then
               do p = 1, nspars
                  LS_R1c(2,i,k1,k2,p) = O1_QATrans(i)   * LS_R1c(2,i,k1,k2,p)
               enddo
            endif
            R1c(2,i,k1,k2) = O1_QATrans(i) * R1c(2,i,k1,k2) + O1_QAMult(i)  * Prc(2,i,k1,k2)
            if (m .gt. 0) then
               if ( do_LC_Jacobians ) then
                  do p = 1, npars
                     LC_R1s(2,i,k1,k2,p) =  LC_O1_QATrans(i,p) *    R1s(2,i,k1,k2)   &
                                          +    O1_QATrans(i)   * LC_R1s(2,i,k1,k2,p) &
                                          + LC_O1_QAMult(i,p)  *    Prs(2,i,k1,k2)   &
                                          +    O1_QAMult(i)    * LC_Prs(2,i,k1,k2,p)
                  enddo
               endif
               if ( do_LS_Jacobians ) then
                  do p = 1, nspars
                     LS_R1s(2,i,k1,k2,p) = O1_QATrans(i)   * LS_R1s(2,i,k1,k2,p)
                  enddo
               endif
               R1s(2,i,k1,k2) = O1_QATrans(i) * R1s(2,i,k1,k2) + O1_QAMult(i)  * Prs(2,i,k1,k2)
            endif
         enddo
      enddo
   enddo

!  End

   return
end subroutine Calculate_firstOrder_LCSPlus

subroutine Calculate_Multipliers_LCPlus &
        ( MaxGeoms, nstreams, nlayers, ngeoms, npars, do_Jacobians,    & ! Inputs
          qstreams, avg_secants, geoms, opdeps, omegas,                & ! Inputs
          LC_avg_secants, L_opdeps, L_omegas,                          & ! Inputs
          O1_QVTrans, O1_QVMult, O1_QATrans, O1_QAMult,                & ! First-order output
          LC_O1_QVTrans, LC_O1_QVMult, LC_O1_QATrans, LC_O1_QAMult,    & ! First-order output
          O2_AVTrans, O2_AVMult, O2_QVMult_d, O2_QAMult_d,             & ! Second-order output
          O2_QAVMult_du, O2_QAVMult_dd,                                & ! Second-order output
          LC_O2_AVTrans, LC_O2_AVMult, LC_O2_QVMult_d, LC_O2_QAMult_d, & ! Second-order output
          LC_O2_QAVMult_du, LC_O2_QAVMult_dd  )                          ! Second-order output

!  Purpose : calculate transmittances and integrated multipliers for R1/R2 fields

!   Convention for directional output as follows

!      Q = quadrature stream direction
!      S = solar scattering direction
!      A = Solar Average-secant direction
!      V = Los viewing direction

   implicit none

!  Inputs
!  ------

!  Control integers

   integer, intent(in) :: MaxGeoms, nlayers, nstreams, ngeoms

!  Jacobian control

   logical, intent(in) :: do_Jacobians
   integer, intent(in) :: npars

!  Cosines Quadrature

   real(fpk), intent(in) :: qstreams(MAXSTREAMS)

!  Average secants

   real(fpk), intent(in) :: avg_secants   (MAXLAYERS,MaxGeoms)
   real(fpk), intent(in) :: LC_avg_secants(MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)

!  Geometries stored in array Geoms(i,j,*)
!    i = 1/2/3 = SZA/VZA/AZM, j = 1/2/3/4 = Degrees/Radians/Cosines/sines

   real(fpk), intent(in) :: geoms(3,4,MaxGeoms)

!  optical properties (these will be deltam-scaled)

   real(fpk), intent(in) :: opdeps(MAXLAYERS)
   real(fpk), intent(in) :: omegas(MAXLAYERS)
   real(fpk), intent(in) :: L_opdeps(MAXLAYERS,MAX_ATMOSWFS)
   real(fpk), intent(in) :: L_omegas(MAXLAYERS,MAX_ATMOSWFS)

!  outputs
!  -------

!  First-order Transmittances and multipliers. Only defined for Quadrature directions.

   real(fpk), intent(out) :: O1_QVTrans(MAXLAYERS,MAXSTREAMS,MaxGeoms)
   real(fpk), intent(out) :: O1_QATrans(MAXLAYERS,MAXSTREAMS,MaxGeoms)

   real(fpk), intent(out) :: O1_QVMult(MAXLAYERS,MAXSTREAMS,MaxGeoms)
   real(fpk), intent(out) :: O1_QAMult(MAXLAYERS,MAXSTREAMS,MaxGeoms)

!  Second-order Transmittances and multipliers

   real(fpk), intent(out) :: O2_AVTrans(MAXLAYERS,MaxGeoms)
   real(fpk), intent(out) :: O2_AVMult (MAXLAYERS,MaxGeoms)

   real(fpk), intent(out) :: O2_QVMult_d(MAXLAYERS,MAXSTREAMS,MaxGeoms)
   real(fpk), intent(out) :: O2_QAMult_d(MAXLAYERS,MAXSTREAMS,MaxGeoms)
   real(fpk), intent(out) :: O2_QAVMult_dd(MAXLAYERS,MAXSTREAMS,MaxGeoms)
   real(fpk), intent(out) :: O2_QAVMult_du(MAXLAYERS,MAXSTREAMS,MaxGeoms)

!  Linearized First-order Transmittances and multipliers

   real(fpk), intent(out) :: LC_O1_QVTrans(MAXLAYERS,MAXSTREAMS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk), intent(out) :: LC_O1_QATrans(MAXLAYERS,MAXSTREAMS,MAX_ATMOSWFS,MaxGeoms)

   real(fpk), intent(out) :: LC_O1_QVMult(MAXLAYERS,MAXSTREAMS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk), intent(out) :: LC_O1_QAMult(MAXLAYERS,MAXSTREAMS,MAX_ATMOSWFS,MaxGeoms)

!  Linearized Second-order Transmittances and multipliers

   real(fpk), intent(out) :: LC_O2_AVTrans(MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk), intent(out) :: LC_O2_AVMult (MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)

   real(fpk), intent(out) :: LC_O2_QVMult_d(MAXLAYERS,MAXSTREAMS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk), intent(out) :: LC_O2_QAMult_d(MAXLAYERS,MAXSTREAMS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk), intent(out) :: LC_O2_QAVMult_dd(MAXLAYERS,MAXSTREAMS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk), intent(out) :: LC_O2_QAVMult_du(MAXLAYERS,MAXSTREAMS,MAX_ATMOSWFS,MaxGeoms)

!  local variables
!  ---------------

   integer   :: k,n,L,p
   real(fpk) :: xv,xs,xsxv, xa,xk,xvxk,xsxk,deltaus,O2_QAMult_u, Help1
   real(fpk) :: sv_plus, kv_plus, kv_minus, ka_plus, ka_minus
   real(fpk) :: dxv,dxs,dxa,dxk,atrans,vtrans,ktrans,hw,term1,term2
   real(fpk) :: L_xa(MAX_ATMOSWFS),L_hw(MAX_ATMOSWFS),L_deltaus(MAX_ATMOSWFS),L_term1
   real(fpk) :: L_dxv(MAX_ATMOSWFS),L_dxs(MAX_ATMOSWFS),L_dxa(MAX_ATMOSWFS),L_dxk(MAX_ATMOSWFS)
   real(fpk) :: L_atrans(MAX_ATMOSWFS),L_vtrans(MAX_ATMOSWFS),L_ktrans(MAX_ATMOSWFS),LC_O2_QAMult_u

!  start of code
!  =============

   do L = 1, ngeoms

!  secants

      xv = one / geoms(2,3,L)      ! Los secant
      xs = one / geoms(1,3,L)      ! solar secant for scattering
      xsxv = xs * xv

!  layer loop

      do n = 1, nlayers

!  Average secant

         xa = avg_secants(n,L) ! solar average-secant for transmittance
         if ( do_Jacobians ) then
            L_xa (1:npars) = LC_avg_secants(n,1:npars,L)
         endif

!  Help variables

         deltaus = opdeps(n)
         dxv = deltaus * xv ; call ExpTrans(BigExp,dxv,vtrans)
         dxa = deltaus * xa ; call ExpTrans(BigExp,dxa,atrans)
         dxs = deltaus * xs
         hw  = omegas(n) * quarter
         if ( do_Jacobians ) then
            L_hw(1:npars)      = L_omegas(n,1:npars) * quarter
            L_deltaus(1:npars) = L_opdeps(n,1:npars)
            L_dxv(1:npars) = L_deltaus(1:npars) * xv
            L_dxs(1:npars) = L_deltaus(1:npars) * xs
            L_dxa(1:npars) = L_deltaus(1:npars) * xa + deltaus * L_xa(1:npars)
         endif

!  Solar to viewing (second Order). Tranmsittance/Multiplier

         sv_plus = one / ( xv + xa)  ; Help1 = xv * sv_plus
         O2_AVTrans(n,L) = vtrans * atrans
         O2_AVMult (n,L) = Help1 * ( one - O2_AVTrans(n,L) )
         if ( do_Jacobians ) then
            call ExpTrans_L(BigExp,npars,dxv,L_dxv,vtrans,L_vtrans)
            call ExpTrans_L(BigExp,npars,dxa,L_dxa,atrans,L_atrans)
            do p = 1, npars
               LC_O2_AVTrans(n,p,L) = L_vtrans(p) * atrans + L_atrans(p) * vtrans
               LC_O2_AVMult (n,p,L) = - Help1*LC_O2_AVTrans(n,p,L) - L_xa(p) * O2_AVMult(n,L) * sv_plus
            enddo
         endif

!  Loop over Nstreams (quadrature)
!  -------------------------------

         do k = 1, nstreams

!  Quadrature secant value

            xk  = one/qstreams(k)

!  Help variables

            dxk  = deltaus * xk ; call ExpTrans(BigExp,dxk,ktrans)
            xvxk = xk * xv ; kv_plus   = xk + xv ; kv_minus   = xk - xv
            xsxk = xk * xs ; ka_plus   = xk + xa ; ka_minus   = xk - xa
            if ( do_Jacobians ) then
               L_dxk(1:npars) = L_deltaus(1:npars) * xk
               call ExpTrans_L(BigExp,npars,dxk,L_dxk,ktrans,L_ktrans)
            endif

!  First Order: viewing direction + all quadrature directions
!    (Formerly, Facj, Chibj)

            term1 =  xvxk / kv_plus
            O1_QVTrans(n,k,L) = ktrans * vtrans ; term2 = one - O1_QVTrans(n,k,L)
            O1_QVMult (n,k,L) = term1 * hw * term2

            if ( do_Jacobians ) then
               do p = 1, npars
                  LC_O1_QVTrans(n,k,p,L) = L_ktrans(p) * vtrans + L_vtrans(p) * ktrans
                  LC_O1_QVMult (n,k,p,L) = term1 * ( L_hw(p) * term2 - hw * LC_O1_QVTrans(n,k,p,L) )
               enddo
            endif

!  First Order: Solar (Av-secant) direction + all quadrature directions
!    (Formerly, Faci, Chibi)

            term1 =  xsxk * hw / ka_plus
            O1_QATrans(n,k,L) = ktrans * atrans ; term2 = one - O1_QATrans(n,k,L) 
            O1_QAMult (n,k,L) = term1 * term2

            if ( do_Jacobians ) then
               do p = 1, npars
                  L_term1 = ( xsxk * L_hw(p) - term1 * L_xa(p) ) / ka_plus
                  LC_O1_QATrans(n,k,p,L) = L_ktrans(p) * atrans + L_atrans(p) * ktrans
                  LC_O1_QAMult (n,k,p,L) = L_term1 * term2 - term1 * LC_O1_QATrans(n,k,p,L)
               enddo
            endif

!  Second Order Multipliers (Formerly PHG, PHH, G, and H )
!     General case (No Small-number series)    
           
            O2_QVMult_d(n,k,L) = ( O2_AVTrans(n,L) - O1_QATrans(n,k,L) ) / kv_minus
            O2_QAMult_d(n,k,L) = ( O2_AVTrans(n,L) - O1_QVTrans(n,k,L) ) / ka_minus
            O2_QAMult_u        = ( one - O1_QVTrans(n,k,L) )  / kv_plus
            O2_QAVMult_du(n,k,L) = xsxk * ( O2_AVMult(n,L) - xv * O2_QVMult_d(n,k,L) ) / ka_plus
            O2_QAVMult_dd(n,k,L) = xsxk * ( O2_AVMult(n,L) - xv * O2_QAMult_u        ) / ka_minus

            if ( do_Jacobians ) then
               do p = 1, npars
                  LC_O2_QVMult_d(n,k,p,L) = ( LC_O2_AVTrans(n,p,L) - LC_O1_QATrans(n,k,p,L) ) / kv_minus
                  LC_O2_QAMult_d(n,k,p,L) = ( LC_O2_AVTrans(n,p,L) - LC_O1_QVTrans(n,k,p,L)  &
                                               + O2_QAMult_d(n,k,L) * L_xa(p) ) / ka_minus
                  LC_O2_QAMult_u          = - LC_O1_QVTrans(n,k,p,L) / kv_plus
                  LC_O2_QAVMult_du(n,k,p,L) = ( xsxk * ( LC_O2_AVMult(n,p,L) - xv * LC_O2_QVMult_d(n,k,p,L) ) &
                                                  - O2_QAVMult_du(n,k,L) * L_xa(p) ) / ka_plus
                  LC_O2_QAVMult_dd(n,k,p,L) = ( xsxk * ( LC_O2_AVMult(n,p,L) - xv * LC_O2_QAMult_u          ) &
                                                  + O2_QAVMult_dd(n,k,L) * L_xa(p) ) / ka_minus
               enddo
            endif

!  end k loop

         enddo

!  end layer loop

      enddo

!  end geometry loop

   enddo

!  Finish

   return
end subroutine Calculate_Multipliers_LCPlus

subroutine Set_avsecant_LCPlus &
     ( MaxGeoms, do_plane_parallel,do_LC_jacobians,nlayers,ngeoms,npars,geoms,&
       SunChapman,opdeps,L_opdeps, avg_secants, LC_avg_secants )

   implicit none

!  inputs and output

   integer  , intent(in)  :: MaxGeoms
   logical  , intent(in)  :: do_plane_parallel, do_LC_jacobians
   integer  , intent(in)  :: nlayers,ngeoms,npars
   real(fpk), intent(in)  :: geoms(3,4,MaxGeoms)

   real(fpk), intent(in)  :: opdeps(MAXLAYERS)
   real(fpk), intent(in)  :: L_opdeps(MAXLAYERS,MAX_ATMOSWFS)
   real(fpk), intent(in)  :: SunChapman (MAXLAYERS,MAXLAYERS,MaxGeoms)

   real(fpk), intent(out) :: avg_secants(MAXLAYERS,MaxGeoms)
   real(fpk), intent(out) :: LC_avg_secants(MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)

! local variables

   real(fpk) :: delta(MAXLAYERS), L_delta(MAXLAYERS,MAX_ATMOSWFS)
   real(fpk) :: sum, sum1, diff, lambda,xs
   integer   :: n,n1,m,L,p,nrev(MAXLAYERS)

!  Zero output

   avg_secants    = zero
   LC_avg_secants = zero

!  Plane parallel

   if ( do_plane_parallel ) then
      do L = 1, ngeoms
         xs = geoms(1,3,L) ; xs  = one / xs
         do n = 1, nlayers
            avg_secants(n,L) = xs
         enddo
      enddo
      return
   endif

!  Pseudo-spherical approximation

   do n = 1, nlayers
      n1 = nlayers + 1 - n ; nrev(n) = n1
      delta(n1) = opdeps(n)
   enddo

   do n = 1, nlayers
      n1 = nrev(n)
      do L = 1, ngeoms
         sum1 = zero
         do m = 1, n - 1
            sum1 = sum1 + SunChapman(n-1,m,L) * delta(m)
         enddo
         sum = zero
         do m = 1, n
            sum = sum + SunChapman(n,m,L) * delta(m)
         enddo
         avg_secants(n1,L) = (sum-sum1)/delta(n)
      enddo
   enddo
   if ( .not. do_LC_jacobians ) return

!  Linearized average secants

   do n = 1, nlayers
      n1 = nrev(n)
      do p = 1, npars
         L_delta(n1,p) = L_opdeps(n,p)
      enddo
   enddo

   do n = 1, nlayers
      n1 = nrev(n)
      do L = 1, ngeoms
         lambda = avg_secants(n1,L)
         do p = 1, npars
            sum1 = zero ; sum = zero
            do m = 1, n - 1
               sum1 = sum1 + SunChapman(n-1,m,L) * L_delta(m,p)
            enddo
            do m = 1, n
               sum = sum + SunChapman(n,m,L) * L_delta(m,p)
            enddo
            diff = sum - sum1
            LC_avg_secants(n1,p,L) = ( diff - L_delta(n,p)*lambda ) / delta(n)
         enddo
      enddo
   enddo

!  Finish

   return
end subroutine Set_avsecant_LCPlus

subroutine Add_Fourier_Component_LCSPlus &
        (m,nstokes,phi,do_lc_jacobians,do_ls_jacobians,npars,nspars,        &
         R2c,R2s,R2cscal,LC_R2c,LC_R2s,LC_R2cscal,LS_R2c,LS_R2s,LS_R2cscal, &
         R2,Icorr,LC_R2,LC_Icorr,LS_R2,LS_Icorr)

   implicit none

!  inputs

   logical  , intent(in) :: do_lc_jacobians,do_ls_jacobians
   integer  , intent(in) :: m,nstokes,npars,nspars
   real(fpk), intent(in) :: phi,R2c(4),R2s(4),R2cscal
   real(fpk), intent(in) :: LC_R2c(4,MAX_ATMOSWFS)  ,LS_R2c(4,MAX_SURFACEWFS)
   real(fpk), intent(in) :: LC_R2s(4,MAX_ATMOSWFS)  ,LS_R2s(4,MAX_SURFACEWFS)
   real(fpk), intent(in) :: LC_R2cscal(MAX_ATMOSWFS),LS_R2cscal(MAX_SURFACEWFS)

!  outputs

   real(fpk), intent(inout) :: R2(4),Icorr
   real(fpk), intent(inout) :: LC_R2(4,MAX_ATMOSWFS)  , LC_Icorr(MAX_ATMOSWFS)
   real(fpk), intent(inout) :: LS_R2(4,MAX_SURFACEWFS), LS_Icorr(MAX_SURFACEWFS)

!  local variables

   integer   :: ki
   real(fpk) :: cosmph,fac,sinmph,rm

!  start of code

   rm = real(m,fpk) * DEG_TO_RAD ; cosmph = cos(rm*phi)
   if (nstokes .eq. 3)   sinmph = sin(rm*phi)
   fac = two ; if (m .eq. 0) fac = one

!  I and Q components

   do ki = 1, 2
      R2(ki) = R2(ki)+fac*R2c(ki)*cosmph
      if (do_LS_Jacobians) then
         LS_R2(ki,1:nspars) = LS_R2(ki,1:nspars) + fac * LS_R2c(ki,1:nspars)*cosmph
      endif
      if (do_LC_Jacobians) then
         LC_R2(ki,1:npars)  = LC_R2(ki,1:npars)  + fac * LC_R2c(ki,1:npars)*cosmph
      endif
   enddo

!  U component

   if (nstokes .eq. 3) then
      ki = nstokes ; R2(ki) = R2(ki)+ fac*R2s(ki)*sinmph
      if (do_LS_Jacobians) then
         LS_R2(ki,1:nspars) = LS_R2(ki,1:nspars) + fac * LS_R2s(ki,1:nspars)*sinmph
      endif
      if (do_LC_Jacobians) then
         LC_R2(ki,1:npars)  = LC_R2(ki,1:npars)  + fac * LC_R2s(ki,1:npars)*sinmph
      endif
   endif

!  scalar correction

   Icorr = Icorr+fac*(R2c(1)-R2cscal)*cosmph
   if (do_LS_Jacobians) then
      LS_Icorr(1:nspars) = LS_Icorr(1:nspars) + &
           fac * (LS_R2c(1,1:nspars)-LS_R2cscal(1:nspars)) * cosmph
   endif
   if (do_LC_Jacobians) then
      LC_Icorr(1:npars)  = LC_Icorr(1:npars)  + &
           fac * (LC_R2c(1,1:npars)-LC_R2cscal(1:npars)) * cosmph
   endif

!  Finish

   return
end subroutine Add_Fourier_Component_LCSPlus

!  End module

end module vlidort_2OScorr_lcs_routines
