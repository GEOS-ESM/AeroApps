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
! #                 Calculate_SecondOrder_LPSPlus           #
! #                 Calculate_FirstOrder_LPSPlus            #
! #                 Add_Fourier_Component_LPSPlus           #
! #                 Calculate_Multipliers_LPPlus            #
! #                 Set_avsecant_LPPlus                     #
! #                                                         #
! ###########################################################

module vlidort_2OScorr_lps_routines

   use VLIDORT_PARS, only : fpk, zero, one, two, three, four,       &
                            six, half, quarter, deg_to_rad,         &
                            BIGEXP, TAYLOR_SMALL, Smallnum6,        &
                            MAXLAYERS, MAXMOMENTS, MAXSTREAMS,      &
                            MAXSTREAMS_P2, MAX_ATMOSWFS, MAX_SURFACEWFS

   use vlidort_2OScorr_utilities, only : Exptrans, Exptrans_L,        &
                                         Make_Trans23, Make_Trans23_P

   implicit none

public  :: Calculate_SecondOrder_LPSPlus,  &
           Calculate_FirstOrder_LPSPlus,   &
           Add_Fourier_Component_LPSPlus,  &
           Calculate_Multipliers_LPPlus,   &
           Set_avsecant_LPPlus


contains

subroutine Calculate_SecondOrder_LPSPlus &
       ( do_LP_Jacobians,do_LS_Jacobians,m,layer,nstreams,nstokes,nlayers,npars,nspars,  & ! Input Control        
         qweights,omega,xv,xa,      LP_omega,LP_xa,               & ! Inputs
         O2_AVTrans,                LP_O2_AVTrans,                & ! Input R2 Multiplier
         O2_QVMult_d,               LP_O2_QVMult_d,               & ! Input R2 Multipliers
         O2_QAMult_d,               LP_O2_QAMult_d,               & ! Input R2 Multipliers
         O2_QAVMult_du,             LP_O2_QAVMult_du,             & ! Input R2 Multipliers
         O2_QAVMult_dd,             LP_O2_QAVMult_dd,             & ! Input R2 Multipliers 
         Ptc,Pts,Prc,Prs,           LP_Ptc,LP_Pts,LP_Prc,LP_Prs,  & ! Inputs
         R1c,R1s,R1cscal,LP_R1c,LP_R1s,LP_R1cscal,LS_R1c,LS_R1s,LS_R1cscal,  & ! Input R1 sources
         R2c,R2s,R2cscal,LP_R2c,LP_R2s,LP_R2cscal,LS_R2c,LS_R2s,LS_R2cscal )   ! Output R2 reflectances

!  Purpose: Update the second-order Reflectance field in Layer n.
!           Also updates the Profile Jacobians

   implicit none

!  1. Standard
!  -----------

!  Control (Fourier #, Layer #, NSTREAMS, NSTOKES)

   integer  , intent(in) :: m, layer, nstreams, nstokes

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

   logical  , intent(in) :: do_LP_Jacobians, do_LS_Jacobians
   integer  , intent(in) :: npars(MAXLAYERS), nspars, nlayers

!  Linearizated General inputs (secants, ssa)

   real(fpk), intent(in) :: LP_omega(MAX_ATMOSWFS)
   real(fpk), intent(in) :: LP_xa(MAXLAYERS,MAX_ATMOSWFS)

!  Linearized source phase matrices

   real(fpk), intent(in) :: LP_Ptc(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LP_Pts(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LP_Prc(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LP_Prs(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)

!  Linearized Second-order Transmittances and multipliers

   real(fpk), intent(in) :: LP_O2_AVTrans(MAXLAYERS,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LP_O2_QVMult_d(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LP_O2_QAMult_d(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LP_O2_QAVMult_dd(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LP_O2_QAVMult_du(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Outputs
!  =======

!  Modified first order inputs (standard & Linearized)

   real(fpk), intent(inout) :: R1c(2,MAXSTREAMS,4,4)
   real(fpk), intent(inout) :: R1s(2,MAXSTREAMS,4,4)
   real(fpk), intent(inout) :: R1cscal(2,MAXSTREAMS)
   real(fpk), intent(inout) :: LP_R1c(2,MAXSTREAMS,4,4,MAXLAYERS,MAX_ATMOSWFS)
   real(fpk), intent(inout) :: LP_R1s(2,MAXSTREAMS,4,4,MAXLAYERS,MAX_ATMOSWFS)
   real(fpk), intent(inout) :: LP_R1cscal(2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
   real(fpk), intent(inout) :: LS_R1c(2,MAXSTREAMS,4,4,MAX_SURFACEWFS)
   real(fpk), intent(inout) :: LS_R1s(2,MAXSTREAMS,4,4,MAX_SURFACEWFS)
   real(fpk), intent(inout) :: LS_R1cscal(2,MAXSTREAMS,MAX_SURFACEWFS)

!  Modified second order outputs (standard & Linearized)

   real(fpk), intent(inout) :: R2c(4),R2s(4),R2cscal
   real(fpk), intent(inout) :: LP_R2c(4,MAXLAYERS,MAX_ATMOSWFS)
   real(fpk), intent(inout) :: LP_R2s(4,MAXLAYERS,MAX_ATMOSWFS)
   real(fpk), intent(inout) :: LP_R2cscal(MAXLAYERS,MAX_ATMOSWFS)
   real(fpk), intent(inout) :: LS_R2c(4,MAX_SURFACEWFS),LS_R2s(4,MAX_SURFACEWFS),LS_R2cscal(MAX_SURFACEWFS)

!  local variables
!  ===============

   integer   :: i, k, l, n, nd, p

   real(fpk) :: dqw, pgv, pha, dglq, dhlq, HW
   real(fpk) :: L_dqw(MAX_ATMOSWFS),  L_hw(MAX_ATMOSWFS)
   real(fpk) :: L_pgv(MAXLAYERS,MAX_ATMOSWFS),  L_pha(MAXLAYERS,MAX_ATMOSWFS)
   real(fpk) :: L_dglq(MAXLAYERS,MAX_ATMOSWFS), L_dhlq(MAXLAYERS,MAX_ATMOSWFS)

   real(fpk) :: S1(4),S2(4),   L_S1(4,MAXLAYERS,MAX_ATMOSWFS),   L_S2(4,MAXLAYERS,MAX_ATMOSWFS)
   real(fpk) :: S3(4),S4(4),   L_S3(4,MAXLAYERS,MAX_ATMOSWFS),   L_S4(4,MAXLAYERS,MAX_ATMOSWFS)
   real(fpk) :: S1scal,S2scal, L_S1scal(MAXLAYERS,MAX_ATMOSWFS), L_S2scal(MAXLAYERS,MAX_ATMOSWFS)
   real(fpk) :: V1(4),V2(4,4), L_V1(4,MAXLAYERS,MAX_ATMOSWFS),   L_V2(4,4,MAXLAYERS,MAX_ATMOSWFS)
   real(fpk) :: V4(4),V6(4,4), L_V4(4,MAXLAYERS,MAX_ATMOSWFS),   L_V6(4,4,MAXLAYERS,MAX_ATMOSWFS)
   real(fpk) :: V1scal,V2scal, L_V1scal(MAXLAYERS,MAX_ATMOSWFS), L_V2scal(MAXLAYERS,MAX_ATMOSWFS)

   real(fpk) :: LS_S1(4,MAX_SURFACEWFS),   LS_S2(4,MAX_SURFACEWFS)  
   real(fpk) :: LS_S3(4,MAX_SURFACEWFS),   LS_S4(4,MAX_SURFACEWFS)
   real(fpk) :: LS_S1scal(MAX_SURFACEWFS), LS_S2scal(MAX_SURFACEWFS)
   real(fpk) :: LS_V1(4,MAX_SURFACEWFS),   LS_V2(4,4,MAX_SURFACEWFS)
   real(fpk) :: LS_V4(4,MAX_SURFACEWFS),   LS_V6(4,4,MAX_SURFACEWFS)
   real(fpk) :: LS_V1scal(MAX_SURFACEWFS), LS_V2scal(MAX_SURFACEWFS)

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

   nd     = layer    !  Dummy variable, debug only

   hw = Omega * quarter
   if ( do_LP_Jacobians ) L_hw(1:npars(layer)) = LP_omega(1:npars(layer)) * quarter 

!  Enter the k-loop
!  ================

   do k = 1, nstreams

!  factors in front of Multipliers

      dglq = hw * O2_QAVMult_du(k) ; pgv = xv * O2_QVMult_d(k)
      dhlq = hw * O2_QAVMult_dd(k) ; pha = xa * O2_QAMult_d(k) 

      L_dglq = zero ; L_dhlq = zero ; L_pgv = zero ; L_pha = zero

      if ( do_LP_Jacobians ) then
        do n = 1, nlayers
          if ( n.eq.layer ) then
            do p = 1, npars(n)
               L_dglq(n,p) = L_hw(p) * O2_QAVMult_du(k) + hw * LP_O2_QAVMult_du(k,n,p) 
               L_dhlq(n,p) = L_hw(p) * O2_QAVMult_dd(k) + hw * LP_O2_QAVMult_dd(k,n,p) 
               L_pgv(n,p)  = xv * LP_O2_QVMult_d(k,n,p)
               L_pha(n,p)  = xa * LP_O2_QAMult_d(k,n,p) + LP_xa(n,p) * O2_QAMult_d(k)
            enddo
          else
            do p = 1, npars(n)
               L_dglq(n,p) = hw * LP_O2_QAVMult_du(k,n,p) 
               L_dhlq(n,p) = hw * LP_O2_QAVMult_dd(k,n,p) 
               L_pgv(n,p)  = xv * LP_O2_QVMult_d(k,n,p)
               L_pha(n,p)  = xa * LP_O2_QAMult_d(k,n,p) + LP_xa(n,p) * O2_QAMult_d(k)
            enddo
          endif
        enddo
      endif

!  V matrices
!  ----------

!  Part 1: The I and Q components

      do l = 1, 2
        V1(l)     = pgv * R1c(2,k,l,1)   + dglq * Prc(2,k,l,1)
        V2(1:2,l) = pha * R1c(1,k,1:2,l) + dhlq * Prc(1,k,1:2,l)
        if ( do_LP_Jacobians ) then
          do n = 1, nlayers
            do p = 1, npars(n)
              L_V1(l,n,p)     =  L_pgv(n,p) * R1c(2,k,l,1)   + pgv * LP_R1c(2,k,l,1,n,p)   &
                              + L_dglq(n,p) * Prc(2,k,l,1)
              L_V2(1:2,l,n,p) =  L_pha(n,p) * R1c(1,k,1:2,l) + pha * LP_R1c(1,k,1:2,l,n,p) &
                              + L_dhlq(n,p) * Prc(1,k,1:2,l) 
            enddo
            if ( n.eq.layer ) then
              do p = 1, npars(n)
                L_V1(l,n,p)     =  L_V1(l,n,p)      + dglq * LP_Prc(2,k,l,1,p)
                L_V2(1:2,l,n,p) =  L_V2(1:2,l,n,p)  + dhlq * LP_Prc(1,k,1:2,l,p)
              enddo
            endif
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
          if ( do_LP_Jacobians ) then
            do n = 1, nlayers
              do p = 1, npars(n)
                L_V4(l,n,p)     = L_pgv(n,p) * R1s(2,k,l,1)   + pgv * LP_R1s(2,k,l,1,n,p)   &
                                + L_dglq(n,p) * Prs(2,k,l,1)
                L_V6(1:2,l,n,p) = L_pha(n,p) * R1s(1,k,1:2,l) + pha * LP_R1s(1,k,1:2,l,n,p) &
                                + L_dhlq(n,p) * Prs(1,k,1:2,l) 
              enddo
              if ( n.eq.layer ) then
                do p = 1, npars(n)
                  L_V4(l,n,p)     =  L_V4(l,n,p)      + dglq * LP_Prs(2,k,l,1,p)
                  L_V6(1:2,l,n,p) =  L_V6(1:2,l,n,p)  + dhlq * LP_Prs(1,k,1:2,l,p)
                enddo
              endif
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
          if ( do_LP_Jacobians ) then
            do n = 1, nlayers
              do p = 1, npars(n)
                L_V2(i,3:4,n,p)  =  L_pha(n,p)  * R1c(1,k,i,3:4) + pha * LP_R1c(1,k,i,3:4,n,p) &
                                  + L_dhlq(n,p) * Prc(1,k,i,3:4)
                L_V6(i,1:2,n,p) =   L_pha(n,p)  * R1s(1,k,i,1:2) + pha * LP_R1s(1,k,i,1:2,n,p) &
                                  + L_dhlq(n,p) * Prs(1,k,i,1:2)
              enddo
              if ( n.eq.layer ) then
                do p = 1, npars(n)
                  L_V2(i,3:4,n,p) = L_V2(i,3:4,n,p) + dhlq * LP_Prc(1,k,i,3:4,p)
                  L_V6(i,1:2,n,p) = L_V6(i,1:2,n,p) + dhlq * LP_Prs(1,k,i,1:2,p)
                enddo
              endif
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
      if ( do_LP_Jacobians ) then
        do n = 1, nlayers
          do p = 1, npars(n)
            L_V1scal(n,p) =  L_pgv(n,p) * R1cscal(2,k) +  pgv * LP_R1cscal(2,k,n,p) &
                           + L_dglq(n,p) * Prc(2,k,1,1)
            L_V2scal(n,p) =  L_pha(n,p) * R1cscal(1,k) +  pha * LP_R1cscal(1,k,n,p) &
                           + L_dhlq(n,p) * Prc(1,k,1,1)
          enddo
          if ( n.eq.layer ) then
            do p = 1, npars(n)
              L_V1scal(n,p) = L_V1scal(n,p) + dglq * LP_Prc(2,k,1,1,p)
              L_V2scal(n,p) = L_V2scal(n,p) + dhlq * LP_Prc(1,k,1,1,p)
            enddo
          endif
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
      if ( do_LP_Jacobians ) then
        L_dqw(1:npars(layer)) = LP_Omega(1:npars(layer)) * qweights(k)
      endif

!  Part 1: I and Q components

      do i = 1, 2

        if ( do_LP_Jacobians ) then
          do n = 1, nlayers
            do p = 1, npars(n)
              do l = 1, 2
                L_S1(i,n,p) = L_S1(i,n,p) + dqw * Ptc(1,k,i,l) * L_V1(l,n,p)
                L_S2(i,n,p) = L_S2(i,n,p) + dqw * Ptc(2,k,l,1) * L_V2(i,l,n,p)
              enddo
            enddo
            if ( n.eq.layer ) then
              do p = 1, npars(n)
                do l = 1, 2
                  L_S1(i,n,p) = L_S1(i,n,p) + L_dqw(p) * Ptc(1,k,i,l)      * V1(l)   + &
                                                dqw    * LP_Ptc(1,k,i,l,p) * V1(l) 
                  L_S2(i,n,p) = L_S2(i,n,p) + L_dqw(p) * Ptc(2,k,l,1)      * V2(i,l) + &
                                                dqw    * LP_Ptc(2,k,l,1,p) * V2(i,l) 
                enddo
              enddo
            endif           
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

        do l = 1, 2
          S1(i) = S1(i) + dqw * Ptc(1,k,i,l) * V1(l)
          S2(i) = S2(i) + dqw * Ptc(2,k,l,1) * V2(i,l)
        enddo

        if (m .gt. 0) then

          if ( do_LP_Jacobians ) then
            do n = 1, nlayers
              do p = 1, npars(n)
                do l = 3, 4
                  L_S3(i,n,p) = L_S3(i,n,p) + dqw * Pts(1,k,i,l) * L_V4(l,n,p)
                  L_S4(i,n,p) = L_S4(i,n,p) + dqw * Pts(2,k,l,1) * L_V6(i,l,n,p)
                enddo
              enddo
              if ( n.eq.layer ) then
                do p = 1, npars(n)
                  do l = 3, 4
                    L_S3(i,n,p) = L_S3(i,n,p) + L_dqw(p) * Pts(1,k,i,l)      * V4(l)   + &
                                                  dqw    * LP_Pts(1,k,i,l,p) * V4(l) 
                    L_S4(i,n,p) = L_S4(i,n,p) + L_dqw(p) * Pts(2,k,l,1)      * V6(i,l) + &
                                                  dqw    * LP_Pts(2,k,l,1,p) * V6(i,l) 
                  enddo
                enddo
              endif           
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
 
          do l = 3,4
            S3(i) = S3(i) + dqw * Pts(1,k,i,l) * V4(l)
            S4(i) = S4(i) + dqw * Pts(2,k,l,1) * V6(i,l)
          enddo

        endif
      enddo

!  Part 2: U component

      do i = 3, nstokes
        if (m .gt. 0) then

          if ( do_LP_Jacobians ) then
            do n = 1, nlayers
              do p = 1, npars(n)
                do l = 3, 4
                  L_S1(i,n,p) = L_S1(i,n,p) + dqw * Ptc(1,k,i,l) * L_V4(l,n,p)
                  L_S2(i,n,p) = L_S2(i,n,p) + dqw * Pts(2,k,l,1) * L_V2(i,l,n,p)
                enddo
              enddo
              if ( n.eq.layer ) then
                do p = 1, npars(n)
                  do l = 3, 4
                    L_S1(i,n,p) = L_S1(i,n,p) + L_dqw(p) * Ptc(1,k,i,l)      * V4(l)   + &
                                                  dqw    * LP_Ptc(1,k,i,l,p) * V4(l) 
                    L_S2(i,n,p) = L_S2(i,n,p) + L_dqw(p) * Pts(2,k,l,1)      * V2(i,l) + &
                                                  dqw    * LP_Pts(2,k,l,1,p) * V2(i,l) 
                  enddo
                enddo
              endif           
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

          do l = 3, 4
            S1(i) = S1(i) + dqw * Ptc(1,k,i,l) * V4(l)
            S2(i) = S2(i) + dqw * Pts(2,k,l,1) * V2(i,l)
          enddo

          if ( do_LP_Jacobians ) then
            do n = 1, nlayers
              do p = 1, npars(n)
                do l = 1, 2
                  L_S3(i,n,p) = L_S3(i,n,p) + dqw * Pts(1,k,i,l) * L_V1(l,n,p)
                  L_S4(i,n,p) = L_S4(i,n,p) + dqw * Ptc(2,k,l,1) * L_V6(i,l,n,p)
                enddo
              enddo
              if ( n.eq.layer ) then
                do p = 1, npars(n)
                  do l = 1, 2
                    L_S3(i,n,p) = L_S3(i,n,p) + L_dqw(p) * Pts(1,k,i,l)      * V1(l)   + &
                                                  dqw    * LP_Pts(1,k,i,l,p) * V1(l) 
                    L_S4(i,n,p) = L_S4(i,n,p) + L_dqw(p) * Ptc(2,k,l,1)      * V6(i,l) + &
                                                  dqw    * LP_Ptc(2,k,l,1,p) * V6(i,l) 
                  enddo
                enddo
              endif           
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

          do l = 1, 2
            S3(i) = S3(i) + dqw * Pts(1,k,i,l) * V1(l)
            S4(i) = S4(i) + dqw * Ptc(2,k,l,1) * V6(i,l)
          enddo

        endif

!  End k-stream loop

      enddo

!  Part 3: scalar

      if ( do_LP_Jacobians ) then
        do n = 1, nlayers
          do p = 1, npars(n)
            L_S1scal(n,p) = L_S1scal(n,p) + dqw * Ptc(1,k,1,1) * L_V1scal(n,p)
            L_S2scal(n,p) = L_S2scal(n,p) + dqw * Ptc(2,k,1,1) * L_V2scal(n,p)
          enddo
          if ( n.eq.layer ) then
            do p = 1, npars(n)
              L_S1scal(n,p) = L_S1scal(n,p) + L_dqw(p) * Ptc(1,k,1,1)      * V1scal + &
                                                dqw    * LP_Ptc(1,k,1,1,p) * V1scal
              L_S2scal(n,p) = L_S2scal(n,p) + L_dqw(p) * Ptc(2,k,1,1)      * V2scal + &
                                                dqw    * LP_Ptc(2,k,1,1,p) * V2scal 
            enddo
          endif           
        enddo
      endif
      if ( do_LS_Jacobians ) then
        do p = 1, nspars
          LS_S1scal = LS_S1scal + dqw * Ptc(1,k,1,1) * LS_V1scal(p)
          LS_S2scal = LS_S2scal + dqw * Ptc(2,k,1,1) * LS_V2scal(p)
        enddo
      endif
      S1scal = S1scal + dqw * Ptc(1,k,1,1) * V1scal
      S2scal = S2scal + dqw * Ptc(2,k,1,1) * V2scal

!  end k loop

   enddo

!  Recursion (Using the solar/LOS transmittance multiplier)
!  =========

!  Update the R2c/R2s terms, with the transmittance

   do i = 1, 2
     if ( do_LP_Jacobians ) then
       do n = 1, nlayers
         do p = 1, npars(n)
           LP_R2c(i,n,p) = R2c(i) * LP_O2_AVTRANS(n,p) + LP_R2c(i,n,p)*O2_AVTRANS
         enddo
       enddo
     endif
     if ( do_LS_Jacobians ) then
       LS_R2c(i,1:nspars) = LS_R2c(i,1:nspars)*O2_AVTRANS
     endif
     R2c(i) = R2c(i)*O2_AVTRANS
   enddo

   if (nstokes .eq. 3) then
     do i = 3, nstokes
       if ( do_LP_Jacobians ) then
         do n = 1, nlayers
           do p = 1, npars(n)
             LP_R2s(i,n,p) = R2s(i) * LP_O2_AVTRANS(n,p) + LP_R2s(i,n,p)*O2_AVTRANS
           enddo
         enddo
       endif
       if ( do_LS_Jacobians ) then
         LS_R2s(i,1:nspars) = LS_R2s(i,1:nspars)*O2_AVTRANS
       endif
       R2s(i) = R2s(i)*O2_AVTRANS
     enddo
   endif

! update the scalar corrections

   if ( do_LP_Jacobians ) then
     do n = 1, nlayers
       do p = 1, npars(n)
          LP_R2cscal(n,p) = R2cscal * LP_O2_AVTRANS(n,p) + LP_R2cscal(n,p)*O2_AVTRANS
       enddo
     enddo
   endif
   if ( do_LS_Jacobians ) then
     LS_R2cscal(1:nspars) = LS_R2cscal(1:nspars)*O2_AVTRANS
   endif
   R2cscal = R2cscal*O2_AVTRANS

   do i = 1, 2
     if ( do_LP_Jacobians ) then
       do n = 1, nlayers
         do p = 1, npars(n)
           LP_R2c(i,n,p) = LP_R2c(i,n,p) + half * ( L_S1(i,n,p) + L_S2(i,n,p) )
         enddo
       enddo
     endif
     if ( do_LS_Jacobians ) then
       LS_R2c(i,1:nspars) = LS_R2c(i,1:nspars) + half * ( LS_S1(i,1:nspars) + LS_S2(i,1:nspars) )
     endif
     R2c(i) = R2c(i) + half * ( S1(i) + S2(i) )
   enddo

   R2cscal = R2cscal + half * ( S1scal + S2scal )
   if ( do_LP_Jacobians ) then
     do n = 1, nlayers
       do p = 1, npars(n)
          LP_R2cscal(n,p) = LP_R2cscal(n,p) + half * ( L_S1scal(n,p) + L_S2scal(n,p) )
       enddo
     enddo
   endif
   if ( do_LS_Jacobians ) then
     LS_R2cscal(1:nspars) = LS_R2cscal(1:nspars) + half * ( LS_S1scal(1:nspars) + LS_S2scal(1:nspars) )
   endif

   if (m .gt. 0) then
     do i = 1, 2
       R2c(i) = R2c(i) + half * ( S3(i) - S4(i) )
       if ( do_LP_Jacobians ) then
         do n = 1, nlayers
           do p = 1, npars(n)
             LP_R2c(i,n,p) = LP_R2c(i,n,p) + half * ( L_S3(i,n,p) - L_S4(i,n,p) )
           enddo
         enddo
       endif
       if ( do_LS_Jacobians ) then
         LS_R2c(i,1:nspars) = LS_R2c(i,1:nspars) + half * ( LS_S3(i,1:nspars) - LS_S4(i,1:nspars) )
       endif
     enddo

     if (nstokes .eq. 3) then
       do i = 3, nstokes
         R2s(i) = R2s(i) + half  * ( S1(i) + S2(i) - S3(i) + S4(i) )
         if ( do_LP_Jacobians ) then
           do n = 1, nlayers
             do p = 1, npars(n)
               LP_R2s(i,n,p) = LP_R2s(i,n,p) + half * &
                 ( L_S1(i,n,p) + L_S2(i,n,p) - L_S3(i,n,p) + L_S4(i,n,p) )
             enddo
           enddo
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
end subroutine Calculate_SecondOrder_LPSPlus

subroutine Calculate_FirstOrder_LPSPlus &
        ( do_LP_Jacobians, do_LS_Jacobians,m,layer,nstreams,nlayers,npars,nspars, & ! Input Control
          Prc,Prs,          LP_Prc,LP_Prs,              & ! Scattering input
          O1_QVTrans,       LP_O1_QVTrans,              & ! First-order multiplier input
          O1_QVMult,        LP_O1_QVMult,               & ! First-order multiplier input
          O1_QATrans,       LP_O1_QATrans,              & ! First-order multiplier input
          O1_QAMult,        LP_O1_QAMult,               & ! First-order multiplier input
          R1c,R1s,R1cscal,LP_R1c,LP_R1s,LP_R1cscal,LS_R1c,LS_R1s,LS_R1cscal ) ! First-order output

!  Purpose: Update the First-Order Reflectance for Layer n
!           Also updates the Profile Jacobians

   implicit none

!  Inputs
!  ======

!  1. Standard
!  -----------

!  Control (Fourier #, NSTREAMS, NSTOKES)

   integer  , intent(in) :: m, layer, nstreams

!  source phase matrices

   real(fpk), intent(in) :: Prc(2,MAXSTREAMS,4,4),Prs(2,MAXSTREAMS,4,4)

!  First-order transmittances/multipliers

   real(fpk), intent(in) :: O1_QVTrans(MAXSTREAMS), O1_QATrans(MAXSTREAMS)
   real(fpk), intent(in) :: O1_QVMult (MAXSTREAMS), O1_QAMult (MAXSTREAMS)

!  2. Linearized
!  -------------

!  Linearization control
 
   logical  , intent(in) :: do_LP_Jacobians,do_LS_Jacobians
   integer  , intent(in) :: npars(MAXLAYERS), nspars,nlayers

!  source phase matrices

   real(fpk), intent(in) :: LP_Prc(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LP_Prs(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)

!  First-order transmittances/multipliers
!     (No Cross-Layer terms for QV)

   real(fpk), intent(in) :: LP_O1_QVTrans(MAXSTREAMS,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LP_O1_QVMult (MAXSTREAMS,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LP_O1_QATrans(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
   real(fpk), intent(in) :: LP_O1_QAMult (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Modified outputs
!  ================

   real(fpk), intent(inout) :: R1c(2,MAXSTREAMS,4,4)
   real(fpk), intent(inout) :: R1s(2,MAXSTREAMS,4,4)
   real(fpk), intent(inout) :: R1cscal(2,MAXSTREAMS)

   real(fpk), intent(inout) :: LP_R1c(2,MAXSTREAMS,4,4,MAXLAYERS,MAX_ATMOSWFS)
   real(fpk), intent(inout) :: LP_R1s(2,MAXSTREAMS,4,4,MAXLAYERS,MAX_ATMOSWFS)
   real(fpk), intent(inout) :: LP_R1cscal(2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

   real(fpk), intent(inout) :: LS_R1c(2,MAXSTREAMS,4,4,MAX_SURFACEWFS)
   real(fpk), intent(inout) :: LS_R1s(2,MAXSTREAMS,4,4,MAX_SURFACEWFS)
   real(fpk), intent(inout) :: LS_R1cscal(2,MAXSTREAMS,MAX_SURFACEWFS)

!  local variables

   integer :: j,i,k1,k2,n,p

!  viewing direction + all j-quadrature directions
!  -----------------------------------------------

!     update the R1 matrices themselves

   do j = 1, nstreams
!  R1scal updates.......
      if ( do_LP_Jacobians ) then
         do n = 1, nlayers
            do p = 1, npars(n)
              LP_R1cscal(1,j,n,p) = LP_R1cscal(1,j,n,p) * O1_QVTrans(j)
            enddo
            if ( n.eq.layer ) then
               do p = 1, npars(n)
                  LP_R1cscal(1,j,n,p) = LP_R1cscal(1,j,n,p) + LP_Prc(1,j,1,1,p) *    O1_QVMULT(j)    &
                                                            +    R1cscal(1,j)   * LP_O1_QVTrans(j,p) &
                                                            +    Prc(1,j,1,1)   * LP_O1_QVMULT(j,p)
               enddo
            endif
         enddo
      endif
      if ( do_LS_Jacobians ) then
         do p = 1, nspars
            LS_R1cscal(1,j,p) =  O1_QVTRANS(j)   * LS_R1cscal(1,j,p)
         enddo
      endif
      R1cscal(1,j) = R1cscal(1,j) * O1_QVTrans(j) + Prc(1,j,1,1) * O1_QVMULT(j)
!  R1c and R1s...............
      do k1 = 1, 4
         do k2 = 1, 4
            if ( do_LP_Jacobians ) then
               do n = 1, nlayers
                  do p = 1, npars(n)
                     LP_R1c(1,j,k1,k2,n,p) = LP_R1c(1,j,k1,k2,n,p) *    O1_QVTrans(j)
                  enddo
                  if ( n.eq.layer ) then
                     do p = 1, npars(n)
                        LP_R1c(1,j,k1,k2,n,p) = LP_R1c(1,j,k1,k2,n,p) + LP_Prc(1,j,k1,k2,p) * O1_QVMULT(j)       &
                                                                      +    R1c(1,j,k1,k2)   * LP_O1_QVTrans(j,p) &
                                                                      +    Prc(1,j,k1,k2)   * LP_O1_QVMULT(j,p)  
                     enddo
                  endif
               enddo
            endif
            if ( do_LS_Jacobians ) then
               do p = 1, nspars
                 LS_R1c(1,j,k1,k2,p) =  O1_QVTRANS(j)   * LS_R1c(1,j,k1,k2,p)
               enddo
            endif
            R1c(1,j,k1,k2) = R1c(1,j,k1,k2) * O1_QVTrans(j) + Prc(1,j,k1,k2) * O1_QVMULT(j)
            if (m .gt. 0) then
               if ( do_LP_Jacobians ) then
                  do n = 1, nlayers
                     do p = 1, npars(n)
                        LP_R1s(1,j,k1,k2,n,p) = LP_R1s(1,j,k1,k2,n,p) * O1_QVTrans(j)
                     enddo
                     if ( n.eq.layer ) then
                        do p = 1, npars(n)
                           LP_R1s(1,j,k1,k2,n,p) = LP_R1s(1,j,k1,k2,n,p) + LP_Prs(1,j,k1,k2,p) * O1_QVMULT(j)       &
                                                                         +    R1s(1,j,k1,k2)   * LP_O1_QVTrans(j,p) &
                                                                         +    Prs(1,j,k1,k2)   * LP_O1_QVMULT(j,p)
                        enddo
                     endif
                  enddo
               endif
               if ( do_LS_Jacobians ) then
                  do p = 1, nspars
                    LS_R1s(1,j,k1,k2,p) = O1_QVTRANS(j)   * LS_R1s(1,j,k1,k2,p)
                  enddo
               endif
               R1s(1,j,k1,k2) = R1s(1,j,k1,k2) * O1_QVTrans(j) + Prs(1,j,k1,k2) * O1_QVMULT(j)
            endif
         enddo
      enddo
   enddo

!  solar direction
!  ---------------

   do i = 1, nstreams
!  R1scal updates.......
      if ( do_LP_Jacobians ) then
         do n = 1, nlayers
            do p = 1, npars(n)
               LP_R1cscal(2,i,n,p) = LP_R1cscal(2,i,n,p) *    O1_QATRANS(i)     &
                                    +   R1cscal(2,i)     * LP_O1_QATRANS(i,n,p)
            enddo
            if ( n.eq.layer ) then
               do p = 1, npars(n)
                  LP_R1cscal(2,i,n,p) = LP_R1cscal(2,i,n,p) + LP_Prc(2,i,1,1,p) * O1_QAMULT(i) &
                                                            +    Prc(2,i,1,1)   * LP_O1_QAMULT(i,n,p)
               enddo
            endif
         enddo
      endif
      if ( do_LS_Jacobians ) then
         do p = 1, nspars
            LS_R1cscal(2,i,p) =  O1_QATrans(i)   * LS_R1cscal(2,i,p)
         enddo
      endif
      R1cscal(2,i) = R1cscal(2,i) * O1_QATRANS(i) + Prc(2,i,1,1) * O1_QAMULT(i)
!  R1c and R1s...............
      do k1 = 1, 4
         do k2 = 1, 4
            if ( do_LP_Jacobians ) then
               do n = 1, nlayers
                  do p = 1, npars(n)
                     LP_R1c(2,i,k1,k2,n,p) = LP_R1c(2,i,k1,k2,n,p) *    O1_QATRANS(i)     &
                                            +   R1c(2,i,k1,k2)     * LP_O1_QATRANS(i,n,p) &
                                            +   Prc(2,i,k1,k2)     * LP_O1_QAMULT(i,n,p)
                  enddo
                  if ( n.eq.layer ) then
                     do p = 1, npars(n)
                        LP_R1c(2,i,k1,k2,n,p) = LP_R1c(2,i,k1,k2,n,p) + LP_Prc(2,i,k1,k2,p) * O1_QAMULT(i)
                     enddo
                  endif
               enddo
            endif
            if ( do_LS_Jacobians ) then
               do p = 1, nspars
                  LS_R1c(2,i,k1,k2,p) = O1_QATrans(i)   * LS_R1c(2,i,k1,k2,p)
               enddo
            endif
            R1c(2,i,k1,k2) = R1c(2,i,k1,k2) * O1_QATRANS(i) + Prc(2,i,k1,k2) * O1_QAMULT(i)
            if (m .gt. 0) then
               if ( do_LP_Jacobians ) then
                  do n = 1, nlayers
                     do p = 1, npars(n)
                        LP_R1s(2,i,k1,k2,n,p) = LP_R1s(2,i,k1,k2,n,p) *    O1_QATRANS(i)     &
                                               +   R1s(2,i,k1,k2)     * LP_O1_QATRANS(i,n,p) &
                                               +   Prs(2,i,k1,k2)     * LP_O1_QAMULT(i,n,p) 
                     enddo
                     if ( n.eq.layer ) then
                        do p = 1, npars(n)
                           LP_R1s(2,i,k1,k2,n,p) = LP_R1s(2,i,k1,k2,n,p) + LP_Prs(2,i,k1,k2,p) * O1_QAMULT(i)
                        enddo
                     endif
                  enddo
               endif
               if ( do_LS_Jacobians ) then
                  do p = 1, nspars
                     LS_R1s(2,i,k1,k2,p) = O1_QATrans(i)   * LS_R1s(2,i,k1,k2,p)
                  enddo
               endif
               R1s(2,i,k1,k2) = R1s(2,i,k1,k2) * O1_QATRANS(i) + Prs(2,i,k1,k2) * O1_QAMULT(i)
            endif
         enddo
      enddo
   enddo

!  End

   return
end subroutine Calculate_firstOrder_LPSPlus

subroutine Calculate_Multipliers_LPPlus &
        ( MaxGeoms, nstreams, nlayers, ngeoms, npars, do_Jacobians,    & ! Inputs
          qstreams, avg_secants, geoms, opdeps, omegas,                & ! Inputs
          LP_avg_secants, L_opdeps, L_omegas,                          & ! Inputs
          O1_QVTrans, O1_QVMult, O1_QATrans, O1_QAMult,                & ! First-order output
          LP_O1_QVTrans, LP_O1_QVMult, LP_O1_QATrans, LP_O1_QAMult,    & ! First-order output
          O2_AVTrans, O2_AVMult, O2_QVMult_d, O2_QAMult_d,             & ! Second-order output
          O2_QAVMult_du, O2_QAVMult_dd,                                & ! Second-order output
          LP_O2_AVTrans, LP_O2_AVMult, LP_O2_QVMult_d, LP_O2_QAMult_d, & ! Second-order output
          LP_O2_QAVMult_du, LP_O2_QAVMult_dd  )                          ! Second-order output

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
   integer, intent(in) :: npars(MAXLAYERS)

!  Cosines Quadrature

   real(fpk), intent(in) :: qstreams(MAXSTREAMS)

!  Average secants

   real(fpk), intent(in) :: avg_secants   (MAXLAYERS,MaxGeoms)
   real(fpk), intent(in) :: LP_avg_secants(MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)

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

   real(fpk), intent(out) :: LP_O1_QVTrans(MAXLAYERS,MAXSTREAMS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk), intent(out) :: LP_O1_QVMult (MAXLAYERS,MAXSTREAMS,MAX_ATMOSWFS,MaxGeoms)

   real(fpk), intent(out) :: LP_O1_QATrans(MAXLAYERS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk), intent(out) :: LP_O1_QAMult (MAXLAYERS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)

!  Linearized Second-order Transmittances and multipliers

   real(fpk), intent(out) :: LP_O2_AVTrans(MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk), intent(out) :: LP_O2_AVMult (MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)

   real(fpk), intent(out) :: LP_O2_QVMult_d(MAXLAYERS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk), intent(out) :: LP_O2_QAMult_d(MAXLAYERS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk), intent(out) :: LP_O2_QAVMult_dd(MAXLAYERS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)
   real(fpk), intent(out) :: LP_O2_QAVMult_du(MAXLAYERS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)

!  local variables
!  ---------------

   integer   :: k,n,L,p,m,np,mp
   real(fpk) :: xv,xs,xsxv,xa,xk,xvxk,xsxk,deltaus,O2_QAMult_u,Help1
   real(fpk) :: sv_plus, kv_plus, kv_minus, ka_plus, ka_minus
   real(fpk) :: dxv,dxs,dxa,dxk,atrans,vtrans,ktrans,hw,term1,term2
   real(fpk) :: L_xa(MAXLAYERS,MAX_ATMOSWFS),L_hw(MAX_ATMOSWFS),L_deltaus(MAX_ATMOSWFS),L_term1,L_term2
   real(fpk) :: L_dxv(MAX_ATMOSWFS),L_dxs(MAX_ATMOSWFS),L_dxa(MAXLAYERS,MAX_ATMOSWFS),L_dxk(MAX_ATMOSWFS)
   real(fpk) :: L_atrans(MAXLAYERS,MAX_ATMOSWFS),L_vtrans(MAX_ATMOSWFS),L_ktrans(MAX_ATMOSWFS),LP_O2_QAMult_u

!  start of code
!  =============

   do L = 1, ngeoms

!  secants

      xv = one / geoms(2,3,L)      ! Los secant
      xs = one / geoms(1,3,L)      ! solar secant for scattering
      xsxv = xs * xv

!  layer loop
!  ----------

      do n = 1, nlayers

!  Average secant

         xa = avg_secants(n,L) ! solar average-secant for transmittance
         if ( do_Jacobians ) then
            do m = n, nlayers
               mp = npars(m) ; L_xa (m,1:mp) = LP_avg_secants(n,m,1:mp,L)
            enddo
         endif

!  Help variables

         deltaus = opdeps(n)
         dxv = deltaus * xv ; call ExpTrans(BigExp,dxv,vtrans)
         dxa = deltaus * xa ; call ExpTrans(BigExp,dxa,atrans)
         dxs = deltaus * xs
         hw  = omegas(n) * quarter
         if ( do_Jacobians ) then
            np = npars(n) ; L_hw(1:np) = L_omegas(n,1:np) * quarter
            L_deltaus(1:np) = L_opdeps(n,1:np)
            L_dxv(1:np)     = L_deltaus(1:np) * xv
            L_dxs(1:np)     = L_deltaus(1:np) * xs
            do m = 1, nlayers
               mp = npars(m) ; L_dxa(m,1:mp) = deltaus * L_xa(m,1:mp)
               if ( m.eq.n ) L_dxa(m,1:mp) = L_dxa(m,1:mp) + L_deltaus(1:mp) * xa
            enddo
         endif

!  Solar to viewing (second Order). Tranmsittance/Multiplier

         sv_plus = one / ( xv + xa)  ; Help1 = xv * sv_plus
         O2_AVTrans(n,L) = vtrans * atrans
         O2_AVMult (n,L) = Help1 * ( one - O2_AVTrans(n,L) )
         if ( do_Jacobians ) then
            call ExpTrans_L(BigExp,npars(n),dxv,L_dxv,vtrans,L_vtrans)
            do m = 1, nlayers
               call ExpTrans_L(BigExp,npars(m),dxa,L_dxa(m,:),atrans,L_atrans(m,:))
               do p = 1, npars(m)
                  LP_O2_AVTrans(n,m,p,L) = L_atrans(m,p) * vtrans
                  LP_O2_AVMult (n,m,p,L) = - Help1 * LP_O2_AVTrans(n,m,p,L) - L_xa(m,p) * O2_AVMult(n,L) * sv_plus
               enddo
               if ( m.eq.n ) then
                  do p = 1, npars(m)
                     L_Term1 = L_vtrans(p) * atrans
                     LP_O2_AVTrans(n,m,p,L) = LP_O2_AVTrans(n,m,p,L) + L_Term1
                     LP_O2_AVMult (n,m,p,L) = LP_O2_AVMult (n,m,p,L) - Help1 * L_Term1
                  enddo
               endif
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
               L_dxk(1:np) = L_deltaus(1:np) * xk
               call ExpTrans_L(BigExp,np,dxk,L_dxk,ktrans,L_ktrans)
            endif

!  First Order: viewing direction + all quadrature directions
!    (Formerly, Facj, Chibj). THIS HAS NO CROSS-LAYER DERIVATIVES !!!!!!!!!!!!!!!!

            term1 =  xvxk / kv_plus
            O1_QVTrans(n,k,L) = ktrans * vtrans ; term2 = one - O1_QVTrans(n,k,L)
            O1_QVMult (n,k,L) = term1 * hw * term2
            if ( do_Jacobians ) then
               do p = 1, npars(n)
                  LP_O1_QVTrans(n,k,p,L) = L_ktrans(p) * vtrans + L_vtrans(p) * ktrans
                  LP_O1_QVMult (n,k,p,L) = term1 * ( L_hw(p) * term2 - hw * LP_O1_QVTrans(n,k,p,L) )
               enddo
            endif

!  First Order: Solar (Av-secant) direction + all quadrature directions
!    (Formerly, Faci, Chibi). THIS TERM HAS CROSS_LAYER DERIVATIVES

            term1 =  xsxk * hw / ka_plus
            O1_QATrans(n,k,L) = ktrans * atrans ; term2 = one - O1_QATrans(n,k,L) 
            O1_QAMult (n,k,L) = term1 * term2

            if ( do_Jacobians ) then
               do m = 1, nlayers
!               do m = n, nlayers
                  do p = 1, npars(m)
                     L_term1 = - term1 * L_xa(m,p) / ka_plus
                     LP_O1_QATrans(n,k,m,p,L) = L_atrans(m,p) * ktrans
                     LP_O1_QAMult (n,k,m,p,L) = - term1 * LP_O1_QATrans(n,k,m,p,L) - L_xa(m,p) * O1_QAMult(n,k,L) / ka_plus
                  enddo
                  if ( m.eq.n ) then
                     do p = 1, npars(m)
                        L_Term2 = L_ktrans(p) * atrans
                        L_term1 = xsxk * L_hw(p) / ka_plus
                        LP_O1_QATrans(n,k,m,p,L) = LP_O1_QATrans(n,k,m,p,L) + L_Term2
                        LP_O1_QAMult (n,k,m,p,L) = LP_O1_QAMult (n,k,m,p,L) - L_Term2 * Term1 + Term2 * L_Term1
                     enddo
                  endif
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
               do m = 1, nlayers
                  do p = 1, npars(m)
                     LP_O2_QVMult_d(n,k,m,p,L) = ( LP_O2_AVTrans(n,m,p,L) - LP_O1_QATrans(n,k,m,p,L) ) / kv_minus
                     LP_O2_QAMult_d(n,k,m,p,L) = ( LP_O2_AVTrans(n,m,p,L) + O2_QAMult_d(n,k,L) * L_xa(m,p) ) / ka_minus


                     LP_O2_QAVMult_du(n,k,m,p,L) = ( xsxk * ( LP_O2_AVMult(n,m,p,L) - xv * LP_O2_QVMult_d(n,k,m,p,L) ) &
                                                  - O2_QAVMult_du(n,k,L) * L_xa(m,p) ) / ka_plus

                     LP_O2_QAVMult_dd(n,k,m,p,L) = ( xsxk * LP_O2_AVMult(n,m,p,L) &
                                                  + O2_QAVMult_dd(n,k,L) * L_xa(m,p) ) / ka_minus
                  enddo
                  if ( m.eq.n ) then
                     do p = 1, npars(m)
                        LP_O2_QAMult_d(n,k,m,p,L) = LP_O2_QAMult_d(n,k,m,p,L) - LP_O1_QVTrans(n,k,p,L) / ka_minus
                        LP_O2_QAMult_u          = - LP_O1_QVTrans(n,k,p,L) / kv_plus
                        LP_O2_QAVMult_dd(n,k,m,p,L) = LP_O2_QAVMult_dd(n,k,m,p,L) -  xsxk * xv * LP_O2_QAMult_u / ka_minus
                    enddo
                  endif
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
end subroutine Calculate_Multipliers_LPPlus

subroutine Set_avsecant_LPPlus &
     ( MaxGeoms, do_plane_parallel, do_LP_jacobians, & ! Flag inputs
       nlayers, ngeoms, npars, geoms,                & ! Control inputs
       SunChapman, opdeps, L_opdeps,                 & ! Geometrical and Optical inputs
       avg_secants, LP_avg_secants )                   ! Output

!  The convention in 2OS is that  Layers are counted from BOA--> TOA (upwards)
!  The convention in FO/LIDORT is Layers are counted from TOA--> BOA (downwards)

!  Thus if using the LIDORT convention for the Chapman factors, one must reverse
!     the indexing on opdeps, and L_opdeps to use it properly

!  Calculations are done from Top downwards in this routine

   implicit none

!  inputs and output

   integer  , intent(in)  :: MaxGeoms
   logical  , intent(in)  :: do_plane_parallel, do_LP_jacobians
   integer  , intent(in)  :: nlayers,ngeoms,npars(MAXLAYERS)
   real(fpk), intent(in)  :: geoms(3,4,MaxGeoms)

   real(fpk), intent(in)  :: opdeps(MAXLAYERS)
   real(fpk), intent(in)  :: L_opdeps(MAXLAYERS,MAX_ATMOSWFS)
   real(fpk), intent(in)  :: SunChapman (MAXLAYERS,MAXLAYERS,MaxGeoms)

   real(fpk), intent(out) :: avg_secants(MAXLAYERS,MaxGeoms)
   real(fpk), intent(out) :: LP_avg_secants(MAXLAYERS,MAXLAYERS,MAX_ATMOSWFS,MaxGeoms)

! local variables

   real(fpk) :: delta(MAXLAYERS), L_delta(MAXLAYERS,MAX_ATMOSWFS)
   real(fpk) :: sum, sum1, lambda,xs, fac
   integer   :: n,n1,m,m1,L,p,nrev(MAXLAYERS)

!  Zero output

   avg_secants    = zero
   LP_avg_secants = zero

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
   if ( .not. do_LP_jacobians ) return

!  Linearized average secants

   do n = 1, nlayers
      n1 = nrev(n)
      do p = 1, npars(n1)
         L_delta(n1,p) = L_opdeps(n,p)
      enddo
   enddo

   do n = 1, nlayers
      n1 = nrev(n)
      do L = 1, ngeoms
         lambda = avg_secants(n1,L)
         do m = 1, n
            m1 = nrev(m)
            if ( m.eq.n ) then
               fac = (SunChapman(n,n,L)-lambda)/delta(n)
               do p = 1, npars(m1)
                  LP_avg_secants(n1,m1,p,L) = fac*L_delta(n,p)
               enddo
            else if (m .lt. n) then
               fac = (SunChapman(n,m,L)-SunChapman(n-1,m,L))/delta(n)
               do p = 1, npars(m1)
                  LP_avg_secants(n1,m1,p,L) = fac*L_delta(m,p)
               enddo
            endif
         enddo
      enddo
   enddo

!  Finish

   return
end subroutine Set_avsecant_LPPlus

subroutine Add_Fourier_Component_LPSPlus &
        (m,nstokes,phi,do_lp_jacobians,do_ls_jacobians,nlayers,npars,nspars, & ! Inputs
         R2c,R2s,R2cscal,LP_R2c,LP_R2s,LP_R2cscal,LS_R2c,LS_R2s,LS_R2cscal,  & ! Inputs
         R2,Icorr,LP_R2,LP_Icorr,LS_R2,LS_Icorr)                               ! Modified Output

   implicit none

!  inputs

   logical  , intent(in) :: do_lp_jacobians,do_ls_jacobians
   integer  , intent(in) :: m,nstokes,nlayers,npars(MAXLAYERS),nspars
   real(fpk), intent(in) :: phi,R2c(4),R2s(4),R2cscal
   real(fpk), intent(in) :: LP_R2c(4,MAXLAYERS,MAX_ATMOSWFS)  ,LS_R2c(4,MAX_SURFACEWFS)
   real(fpk), intent(in) :: LP_R2s(4,MAXLAYERS,MAX_ATMOSWFS)  ,LS_R2s(4,MAX_SURFACEWFS)
   real(fpk), intent(in) :: LP_R2cscal(MAXLAYERS,MAX_ATMOSWFS),LS_R2cscal(MAX_SURFACEWFS)

!  outputs

   real(fpk), intent(inout) :: R2(4),Icorr
   real(fpk), intent(inout) :: LP_R2(4,MAXLAYERS,MAX_ATMOSWFS),  LS_R2(4,MAX_SURFACEWFS)
   real(fpk), intent(inout) :: LP_Icorr(MAXLAYERS,MAX_ATMOSWFS), LS_Icorr(MAX_SURFACEWFS)

!  local variables

   integer   :: ki, n, nq
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
      if (do_LP_Jacobians) then
         do n = 1, nlayers
            nq = npars(n)
            LP_R2(ki,n,1:nq)  = LP_R2(ki,n,1:nq)  + fac * LP_R2c(ki,n,1:nq)*cosmph
         enddo
      endif
   enddo

!  U component

   if (nstokes .eq. 3) then
      ki = nstokes ; R2(ki) = R2(ki)+ fac*R2s(ki)*sinmph
      if (do_LS_Jacobians) then
         LS_R2(ki,1:nspars) = LS_R2(ki,1:nspars) + fac * LS_R2s(ki,1:nspars)*sinmph
      endif
      if (do_LP_Jacobians) then
         do n = 1, nlayers
            nq = npars(n)
            LP_R2(ki,n,1:nq)  = LP_R2(ki,n,1:nq)  + fac * LP_R2s(ki,n,1:nq)*sinmph
         enddo
      endif
   endif

!  scalar correction

   Icorr = Icorr+fac*(R2c(1)-R2cscal)*cosmph
   if (do_LS_Jacobians) then
      LS_Icorr(1:nspars) = LS_Icorr(1:nspars) + &
           fac * (LS_R2c(1,1:nspars)-LS_R2cscal(1:nspars)) * cosmph
   endif
   if (do_LP_Jacobians) then
      do n = 1, nlayers
         nq = npars(n)
            LP_Icorr(n,1:nq)  = LP_Icorr(n,1:nq)  + &
                  fac * (LP_R2c(1,n,1:nq)-LP_R2cscal(n,1:nq)) * cosmph
      enddo
   endif

!  Finish

   return
end subroutine Add_Fourier_Component_LPSPlus

!  End module

end module vlidort_2OScorr_lps_routines

