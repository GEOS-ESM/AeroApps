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
! #                 Calculate_SecondOrder                   #
! #                 Calculate_Zmatrix_ALT                   #
! #                 Scale_SecondOrder                       #
! #                 Calculate_Multipliers                   #
! #                 Calculate_FirstOrder                    #
! #                 Add_Fourier_Component                   #
! #                 Test_Fourier_Convergence                #
! #                 Set_avsecant                            #
! #                                                         #
! #                 [Calculate_Zmatrix]  Removed, 3/4/15    #
! #                                                         #
! ###########################################################

module vlidort_2OScorr_routines

   use VLIDORT_PARS, only : fpk, zero, one, two, three, four,       &
                            quarter, six, half, deg_to_rad,         &
                            BIGEXP, TAYLOR_SMALL, Smallnum6,        &
                            MAXLAYERS, MAXMOMENTS, MAXSTREAMS, MAXSTREAMS_p2

   use vlidort_2OScorr_utilities, only : Exptrans, Taylor_1, Taylor_3,  Make_Trans23

   implicit none

!  Note that Calculate_Zmatrix is not thread-safe, because of save statements
!       The alternative routine Calculate_Zmatrix_ALT is thread-safe

public  :: Calculate_SecondOrder,    &
           Calculate_Zmatrix_ALT,    &
           Scale_SecondOrder,        &
           Calculate_Multipliers,    &
           Calculate_FirstOrder,     &
           Add_Fourier_Component,    &
           Test_Fourier_Convergence, &
           Set_avsecant

!  Removed 3/4/15
!           Calculate_Zmatrix,    &

contains

subroutine Calculate_SecondOrder &
        ( m,n,nstreams,nstokes,               & ! Input Control
          qweights,omega,xv,xa,               & ! Input local
          Ptc,Pts,Prc,Prs,                    & ! Input Local scattering matrices
          O2_AVTrans,O2_QVMult_d,O2_QAMult_d, & ! Input R2 Multipliers
          O2_QAVMult_du,O2_QAVMult_dd,        & ! Input R2 Multipliers 
          R1c,R1s,R1cscal,                    & ! Input R1 source field
          R2c,R2s,R2cscal )                     ! Output R2 reflectance

!  Purpose: Update the second-order Reflectance field in Layer n.

   implicit none

!  inputs
!  ======

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

!  Outputs
!  =======

!  Modified first order inputs

   real(fpk), intent(inout) :: R1c(2,MAXSTREAMS,4,4)
   real(fpk), intent(inout) :: R1s(2,MAXSTREAMS,4,4)
   real(fpk), intent(inout) :: R1cscal(2,MAXSTREAMS)

!  Modified second order outputs

   real(fpk), intent(inout) :: R2c(4),R2s(4),R2cscal

!  local variables
!  ===============

   integer   :: i, k, l, nd
   real(fpk) :: dqw, pgv, pha, dglq, dhlq, hw

   real(fpk) :: S1(4),S2(4)
   real(fpk) :: S3(4),S4(4)
   real(fpk) :: S1scal,S2scal
   real(fpk) :: V1(4),V2(4,4)
   real(fpk) :: V4(4),V6(4,4)
   real(fpk) :: V1scal,V2scal

!  initial section
!  ---------------

!  Zero the  S-matrices

   S1     = zero ; S2     = zero
   S1scal = zero ; S2scal = zero
   S3     = zero ; S4     = zero

   nd     = n    !  Dummy variable, debug only

   hw = Omega * quarter

!  Enter the k-loop
!  ================

   do k = 1, nstreams

!  factors in front of Multipliers

      dglq = hw * O2_QAVMult_du(k) ; pgv = xv * O2_QVMult_d(k)
      dhlq = hw * O2_QAVMult_dd(k) ; pha = xa * O2_QAMult_d(k) 

!  V matrices
!  ----------

!  Part 1: The I and Q components

      do l = 1, 2
         V1(l)     = pgv * R1c(2,k,l,1)   + dglq * Prc(2,k,l,1)
         V2(1:2,l) = pha * R1c(1,k,1:2,l) + dhlq * Prc(1,k,1:2,l)
      enddo        
      if (m .gt. 0) then
         do l = 3, 4
            V4(l)     = pgv * R1s(2,k,l,1)   + dglq * Prs(2,k,l,1)
            V6(1:2,l) = pha * R1s(1,k,1:2,l) + dhlq * Prs(1,k,1:2,l)
         enddo
      endif

!  Part 2: the U component

      if (m .gt. 0) then
         do i = 3, nstokes
            V2(i,3:4) = pha * R1c(1,k,i,3:4) + dhlq * Prc(1,k,i,3:4)
            V6(i,1:2) = pha * R1s(1,k,i,1:2) + dhlq * Prs(1,k,i,1:2)
         enddo
      endif

!  Part 3: the scalar calculation

      V1scal = pgv * R1cscal(2,k) + dglq * Prc(2,k,1,1)
      V2scal = pha * R1cscal(1,k) + dhlq * Prc(1,k,1,1)

!  debug

!      write(97,'(2i3,1p5e20.12)')m,n,V4(3),V4(4),V6(3,1)
!      if (m.eq.3.and.n.eq.60)pause'gronk'

!  S Matrices
!  ----------

      dqw = Omega * qweights(k)

!  Part 1: I and Q components

      do i = 1, 2
         do l = 1, 2
            S1(i) = S1(i) + dqw * Ptc(1,k,i,l) * V1(l)
            S2(i) = S2(i) + dqw * Ptc(2,k,l,1) * V2(i,l)
         enddo
         if (m .gt. 0) then
            do l = 3,4
              S3(i) = S3(i) + dqw * Pts(1,k,i,l) * V4(l)
              S4(i) = S4(i) + dqw * Pts(2,k,l,1) * V6(i,l)
            enddo
         endif
      enddo

!  Part 2: U component

      do i = 3, nstokes
         if (m .gt. 0) then
            do l = 3, 4
              S1(i) = S1(i) + dqw * Ptc(1,k,i,l) * V4(l)
              S2(i) = S2(i) + dqw * Pts(2,k,l,1) * V2(i,l)
            enddo
            do l = 1, 2
              S3(i) = S3(i) + dqw * Pts(1,k,i,l) * V1(l)
              S4(i) = S4(i) + dqw * Ptc(2,k,l,1) * V6(i,l)
            enddo
         endif
      enddo

!  Part 3: scalar

      S1scal = S1scal + dqw * Ptc(1,k,1,1) * V1scal
      S2scal = S2scal + dqw * Ptc(2,k,1,1) * V2scal

!  end k loop

   enddo

!  Recursion (Using the solar/LOS transmittance multiplier)
!  =========

!  Update the R2c/R2s terms, with the transmittance

   do i = 1, 2
      R2c(i) = R2c(i)*O2_AVTrans
   enddo
   do i = 3, nstokes
      R2s(i) = R2s(i)*O2_AVTrans
   enddo

! update the scalar corrections

   R2cscal = R2cscal*O2_AVTrans
   do i = 1, 2
      R2c(i) = R2c(i) + half * ( S1(i) + S2(i) )
   enddo
   R2cscal = R2cscal + half * ( S1scal + S2scal )

! Update sine terms
  
   if (m .gt. 0) then
      do i = 1, 2
         R2c(i) = R2c(i) + half * ( S3(i) - S4(i) )
      enddo
      do i = 3, nstokes
         R2s(i) = R2s(i) + half  * ( S1(i) + S2(i) - S3(i) + S4(i) )
      enddo
   endif

!  debug
!      write(*,'(2i3,1p5e20.12)')m,n,R2cscal
!      write(94,'(2i3,1p5e20.12)')m,n,R2c(1),R2c(2),R2c(3)
!      write(94,'(2i3,1p5e20.12)')m,n,R2s(1),R2s(2),R2s(3)
!      if (m.eq.3.and.n.eq.60)pause'gronk'
!   if(n.eq.20)pause 'Feb 3'

!  Finish

   return
end subroutine Calculate_SecondOrder

subroutine Calculate_FirstOrder &
        ( m, n, nstreams, Prc, Prs,  & ! Control + scattering   input
          O1_QVTrans, O1_QVMult,     & ! First-order multiplier input
          O1_QATrans, O1_QAMult,     & ! First-order multiplier input
          R1c, R1s, R1cscal )          ! First-order output

!  Purpose: Update the First-Order Reflectance for Layer n

   implicit none

!  inputs
!  ======

!  Control and scattering inputs

   integer  , intent(in) :: m, n, nstreams
   real(fpk), intent(in) :: Prc(2,MAXSTREAMS,4,4),Prs(2,MAXSTREAMS,4,4)

!  Multiplier inputs

   real(fpk), intent(in) :: O1_QVTrans(MAXSTREAMS), O1_QATrans(MAXSTREAMS)
   real(fpk), intent(in) :: O1_QVMult (MAXSTREAMS), O1_QAMult (MAXSTREAMS)

!  Modified outputs
!  ================

   real(fpk), intent(inout) :: R1c(2,MAXSTREAMS,4,4)
   real(fpk), intent(inout) :: R1s(2,MAXSTREAMS,4,4)
   real(fpk), intent(inout) :: R1cscal(2,MAXSTREAMS)

!  local variables
!  ===============

   integer :: j,i,k1,k2

!  viewing direction + all j-quadrature directions
!     update the R1 matrices themselves

   do j = 1, nstreams
      R1cscal(1,j) =  O1_QVTrans(j) * R1cscal(1,j) &
                    + O1_QVMult(j)  * Prc(1,j,1,1)
      do k1 = 1, 4
         do k2 = 1, 4
            R1c(1,j,k1,k2) =  O1_QVTrans(j) * R1c(1,j,k1,k2) &
                            + O1_QVMult(j)  * Prc(1,j,k1,k2)
            if (m .gt. 0) then
               R1s(1,j,k1,k2) =   O1_QVTrans(j) * R1s(1,j,k1,k2) &
                                + O1_QVMult(j)  * Prs(1,j,k1,k2)
            endif
         enddo
      enddo
   enddo

!  solar direction
!  ---------------

   do i = 1, nstreams
      R1cscal(2,i) =  O1_QATrans(i) * R1cscal(2,i) &
                    + O1_QAMult(i)  * Prc(2,i,1,1)
      do k1 = 1, 4
         do k2 = 1, 4
            R1c(2,i,k1,k2) =  O1_QATrans(i) * R1c(2,i,k1,k2)  &
                            + O1_QAMult(i)  * Prc(2,i,k1,k2)
            if (m .gt. 0) then
               R1s(2,i,k1,k2) =  O1_QATrans(i) * R1s(2,i,k1,k2) &
                               + O1_QAMult(i)  * Prs(2,i,k1,k2)
            endif
         enddo
      enddo
   enddo

!  Finish

   return
end subroutine Calculate_firstOrder

!

subroutine Calculate_Zmatrix_ALT &
        ( m, nstokes, nstreams, ns1, ns2, nmoments,  & ! Inputs
          rl, r2Lp1, sql4, sqlm, binfac, coefsm,     & ! Inputs
          xmu, xsi, Plm_2all, Plm_mlp, Plm_mlm,      & ! Inputs
          Ptc,Pts,Prc,Prs )

   implicit none

   
!  inputs

   integer  , intent(in) :: m,nstreams,nstokes,nmoments

!  Coefficients

   real(fpk), intent(in) :: coefsm(0:MAXMOMENTS,6)

!  Help variables (formerly saved)

   integer  , intent(inout) :: ns1, ns2
   real(fpk), intent(inout) :: rl(0:2*MAXMOMENTS)
   real(fpk), intent(inout) :: r2lp1(0:MAXMOMENTS)
   real(fpk), intent(inout) :: sqlm(0:MAXMOMENTS),sql4(0:MAXMOMENTS)
   real(fpk), intent(inout) :: binfac
   real(fpk), intent(inout) :: xmu(MAXSTREAMS_p2) 
   real(fpk), intent(inout) :: xsi(MAXSTREAMS_p2)
   real(fpk), intent(inout) :: Plm_2all(MAXSTREAMS_p2,5)
   real(fpk), intent(inout) :: Plm_mlp(MAXSTREAMS_p2)
   real(fpk), intent(inout) :: Plm_mlm(MAXSTREAMS_p2)

!  outputs

   real(fpk), intent(out) :: Ptc(2,MAXSTREAMS,4,4),Pts(2,MAXSTREAMS,4,4), &
                             Prc(2,MAXSTREAMS,4,4),Prs(2,MAXSTREAMS,4,4)

!  local variables
!  ---------------

!  help variables

   integer   :: l,j,k1,k2,i,lold,lnew,itmp
   real(fpk) :: parity,u
   real(fpk) :: Zmpls_xv(nstreams,4,4),Zmmin_xv(nstreams,4,4)
   real(fpk) :: Zmpls_xs(nstreams,4,4),Zmmin_xs(nstreams,4,4)
   real(fpk) :: Plm(MAXSTREAMS_p2,3,2)
   real(fpk) :: DPDpl(MAXSTREAMS_p2,4),DPDmi(MAXSTREAMS_p2,4)
   real(fpk) :: DSD(4,4), SPj, f1new, f1old, tmp, f2new,f2newa,f2old

!  Zero the transmittance and reflectance matrices

   Ptc = zero ; Pts = zero ; Prc = zero ; Prs = zero

!  Zero local values

   Zmpls_xv = zero ; Zmmin_xv = zero
   Zmpls_xs = zero ; Zmmin_xs = zero

!  Start development of spherical functions

   lold = 1 ; lnew = 2
   do i = 1, ns2
      if (m .ne. 0) then
         Plm(i,1,lnew) = binfac*xsi(i)**m
      else
         Plm(i,1,lnew) = binfac
      endif
      Plm(i,1,lold)   = zero
      Plm(i,2:3,lnew) = zero ; Plm(i,2:3,lold) = zero
   enddo

!  Start coefficient loop
!  ----------------------

   parity = - one
   do l = m, nmoments
      parity = -parity

!  Initial Associated Legendre functions

      if (l .eq. max0(m,2)) then
         if (m .eq. 0) then
            do i = 1, ns2
              Plm(i,2,lnew) = Plm_2all(i,1)
              Plm(i,3,lnew) = Plm(i,2,lnew)
            enddo
         else if (m .eq. 1) then
            do i = 1, ns2
              Plm(i,2,lnew) = Plm_2all(i,2)
              Plm(i,3,lnew) = Plm_2all(i,3)
            enddo
         else if (m .eq. 2) then
            do i = 1, ns2
              Plm(i,2,lnew) = Plm_2all(i,4)
              Plm(i,3,lnew) = Plm_2all(i,5)
            enddo
         else
            do i = 1, ns2
              Plm(i,2,lnew) = Plm_mlm(i)
              Plm(i,3,lnew) = Plm_mlp(i)
            enddo
         endif
      endif
      do i = 1, ns2
         DPDpl(i,1:3) = Plm(i,1:3,lnew)
         DPDpl(i,4)   = DPDpl(i,1)
         DPDmi(i,1)   = parity*Plm(i,1,lnew)
         DPDmi(i,2)   = parity*Plm(i,3,lnew)
         DPDmi(i,3)   = parity*Plm(i,2,lnew)
         DPDmi(i,4)   = DPDmi(i,1)
      enddo

!  here is the DSD ("Greek" matrix)

      DSD (1, 1) = coefsm (l, 1)
      DSD (2, 1) = half * coefsm (l, 5)
      DSD (2, 2) = quarter * (coefsm (l, 2) + coefsm (l, 3) )
      DSD (3, 2) = quarter * (coefsm (l, 2) - coefsm (l, 3) )
      DSD (3, 1) = DSD (2, 1)
      DSD (1, 2) = DSD (2, 1)
      DSD (1, 3) = DSD (2, 1)
      DSD (2, 3) = DSD (3, 2)
      DSD (3, 3) = DSD (2, 2)

      DSD (1, 4) = zero
      DSD (2, 4) = half * coefsm (l, 6)
      DSD (3, 4) = - DSD (2, 4)
      DSD (4, 4) = coefsm (l, 4)
      DSD (4, 1) = zero
      DSD (4, 2) = - DSD (2, 4)
      DSD (4, 3) = - DSD (3, 4)

!  now make the Z matrices for scattering: Polar to View streams

      do k2 = 1, 4
         do k1 = 1, 4
            do j = 1, nstreams
               SPj = DSD(k1,k2) * DPDpl(j,k2)
               Zmpls_xv(j,k1,k2) = Zmpls_xv(j,k1,k2)+DPDpl(ns2,k1)*SPj
               Zmmin_xv(j,k1,k2) = Zmmin_xv(j,k1,k2)+DPDmi(ns2,k1)*SPj
            enddo
         enddo
      enddo

!  now make the Z matrices for scattering: Sun to Polar streams

      do k1 = 1, 4
         SPj = DSD(k1,1) * DPDpl(ns1,1)
         do i = 1, nstreams
            Zmpls_xs(i,k1,1) = Zmpls_xs(i,k1,1)+DPDpl(i,k1)*SPj
            Zmmin_xs(i,k1,1) = Zmmin_xs(i,k1,1)+DPDmi(i,k1)*SPj
         enddo
      enddo

!  Finishing point

      if (l .eq. nmoments) goto 12

!  Upgrade Plm functions

      f1new = r2lp1(l) / sqlm(l+1)
      f1old = sqlm(l)  /  sqlm(l+1)
      do i = 1, ns2
         u  = xmu(i)
         Plm(i,1,lold) = f1new*u*Plm(i,1,lnew)-f1old*Plm(i,1,lold)
      enddo

      if (l .ge. max0(m,2)) then
         tmp    = one/(rl(l)*sql4(l+1)*sqlm(l+1))
         f2new  =  r2lp1(l) * rl(l)   * rl(l+1) * tmp
         f2newa =  r2lp1(l) * two     * rl(m)   * tmp
         f2old  = rl(l+1)   * sql4(l) * sqlm(l) * tmp
         do i = 1, ns2
            u = xmu(i)
            Plm(i,2,lold) = (f2new*u+f2newa)*Plm(i,2,lnew)-f2old*Plm(i,2,lold)
            Plm(i,3,lold) = (f2new*u-f2newa)*Plm(i,3,lnew)-f2old*Plm(i,3,lold)
         enddo
      endif
      itmp = lnew ; lnew = lold ; lold = itmp

!  End coefficient loop

   enddo

12 continue

!  Transformation 2/3

   do j = 1, nstreams
      call Make_Trans23(nstreams,j,Zmmin_xv)
      call Make_Trans23(nstreams,j,Zmpls_xv)
      call Make_Trans23(nstreams,j,Zmmin_xs)
      call Make_Trans23(nstreams,j,Zmpls_xs)
   enddo

!      if ( layer.eq.30)write(*,*)m,l,Zmmin_xv(8,2,3),Zmmin_xv(8,3,2)
!      if ( m.eq.3.and.layer.eq.30)pause

!  Set results for Polar incoming - View angle outgoing

   do j = 1, nstreams
      do k1 = 1, 2
         Ptc(1,j,k1,1:2) = Zmpls_xv(j,k1,1:2)
         Prc(1,j,k1,1:2) = Zmmin_xv(j,k1,1:2)
         if (m .gt. 0) then
            Pts(1,j,k1,3:4) = -Zmpls_xv(j,k1,3:4)
            Prs(1,j,k1,3:4) = -Zmmin_xv(j,k1,3:4)
         endif
      enddo
      do k1 = 3, nstokes
         Ptc(1,j,k1,3:4) = Zmpls_xv(j,k1,3:4)
         Prc(1,j,k1,3:4) = Zmmin_xv(j,k1,3:4)
         if (m .gt. 0) then
            Pts(1,j,k1,1:2) = Zmpls_xv(j,k1,1:2)
            Prs(1,j,k1,1:2) = Zmmin_xv(j,k1,1:2)
         endif
      enddo
   enddo

!  Set results for Solar incoming - Polar angle outgoing

   do i = 1, nstreams
      Ptc(2,i,1:2,1) = Zmpls_xs(i,1:2,1)
      Prc(2,i,1:2,1) = Zmmin_xs(i,1:2,1)
      if (m .gt. 0) then
         Pts(2,i,3:4,1) = Zmpls_xs(i,3:4,1)
         Prs(2,i,3:4,1) = Zmmin_xs(i,3:4,1)
      endif
   enddo

!  Finish routine

   return
end subroutine Calculate_Zmatrix_ALT

subroutine Scale_SecondOrder &
         ( MaxGeoms,nstokes,n,L,xv,xa,  & ! I
           facL,omega1, omega2, opdep,  & ! I
           FO_Zmatrix, FO_R1saved,      & ! I
           R2c,R2cscal )                  ! Modified

   implicit none

!  Inputs
!  ------

!  Control integers

   integer, intent(in) :: MaxGeoms,nstokes, n, L

!  Optical properties and angles

   real(fpk), intent(in) :: omega1, omega2, opdep, xv, xa, facl

!  First-order results, only for Deltam-scaling

   real(fpk), intent(in) :: FO_Zmatrix(MAXLAYERS,4,MaxGeoms)
   real(fpk), intent(in) :: FO_R1saved(0:MAXLAYERS,4,MaxGeoms)

!  Modified outputs

   real(fpk), intent(inout) :: R2c(4),R2cscal

!  Local variables

   integer   :: k1
   real(fpk) :: caph,omeg,bh,bhp1,omfh,bbb1,bbbb
   real(fpk) :: acons,bcons,aafac,bbfac,R2scal_z,R2_z

!  Develop the multipliers (aafac and bbfac)

   caph  = xv + xa
   omeg  = omega1 + omega2
   bbb1  = xv * 0.25_fpk / caph / caph
   bbbb  = omeg * xa * bbb1

   bh    = opdep * caph
   bhp1  = 1.0d0 + bh
 
   acons = opdep * facL
   bcons = bbbb * ( 1.0d0 - facL * bhp1 )
 
   omfh  = omega2 * caph
   aafac = omfh * acons
   bbfac = omfh * bcons

!  Compute source terms and add layer contribution

   R2scal_z = aafac * FO_R1saved(n-1,1,L) + bbfac * FO_Zmatrix(n,1,L)
   R2cscal  = R2cscal + R2scal_z
   do k1 = 1, nstokes
      R2_z = aafac * FO_R1saved(n-1,k1,L) + bbfac * FO_Zmatrix(n,k1,L)
      R2c(k1) = R2c(k1) + R2_z
   enddo

!  Finish routine

   return
end subroutine Scale_SecondOrder

subroutine Calculate_Multipliers &
        ( MaxGeoms, nstreams, nlayers, ngeoms,             & ! Inputs
          qstreams, avg_secants, geoms, opdeps, omegas,    & ! Inputs
          O1_QVTrans, O1_QVMult, O1_QATrans,  O1_QAMult,   & ! First-order output
          O2_AVTrans, O2_AVMult, O2_QVMult_d, O2_QAMult_d, & ! Second-order output
          O2_QAVMult_du, O2_QAVMult_dd  )                    ! Second-order output

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

   integer  , intent(in) :: MaxGeoms, nlayers, nstreams, ngeoms

!  Cosines Quadrature

   real(fpk), intent(in) :: qstreams(MAXSTREAMS)

!  Average secants

   real(fpk), intent(in) :: avg_secants(MAXLAYERS,MaxGeoms)

!  Geometries stored in array Geoms(i,j,*)
!    i = 1/2/3 = SZA/VZA/AZM, j = 1/2/3/4 = Degrees/Radians/Cosines/sines

   real(fpk), intent(in) :: geoms(3,4,MaxGeoms)

!  optical properties (these will be deltam-scaled)

   real(fpk), intent(in) :: opdeps(MAXLAYERS)
   real(fpk), intent(in) :: omegas(MAXLAYERS)

!  outputs
!  -------

!  First-order Transmittances and multipliers. Only defined for Quadrature directions.

   real(fpk), intent(out) :: O1_QVTrans(MAXLAYERS,MAXSTREAMS,MaxGeoms)
   real(fpk), intent(out) :: O1_QATrans(MAXLAYERS,MAXSTREAMS,MaxGeoms)

   real(fpk), intent(out) :: O1_QVMult(MAXLAYERS,MAXSTREAMS,MaxGeoms)
   real(fpk), intent(out) :: O1_QAMult(MAXLAYERS,MAXSTREAMS,MaxGeoms)

!  Second-order Transmittances and multipliers

   real(fpk), intent(out) :: O2_AVTrans(MAXLAYERS,MaxGeoms)
   real(fpk), intent(out) :: O2_AVMult(MAXLAYERS,MaxGeoms)

   real(fpk), intent(out) :: O2_QVMult_d(MAXLAYERS,MAXSTREAMS,MaxGeoms)
   real(fpk), intent(out) :: O2_QAMult_d(MAXLAYERS,MAXSTREAMS,MaxGeoms)
   real(fpk), intent(out) :: O2_QAVMult_dd(MAXLAYERS,MAXSTREAMS,MaxGeoms)
   real(fpk), intent(out) :: O2_QAVMult_du(MAXLAYERS,MAXSTREAMS,MaxGeoms)

!  local variables
!  ---------------

   integer   :: k,n,L
   real(fpk) :: xv,xs,xa,xk,xvxk,xsxk,deltaus,O2_QAMult_u
   real(fpk) :: sv_plus, kv_plus, kv_minus, ka_plus, ka_minus
   real(fpk) :: dxv,dxs,dxa,dxk,atrans,vtrans,ktrans,hw
   real(fpk) :: arg1, arg2, arg3, T1, T3

!  start of code
!  =============

   do L = 1, ngeoms

!  secants

      xv = one / geoms(2,3,L)      ! Los secant
      xs = one / geoms(1,3,L)      ! solar secant for scattering

!  layer loop

      do n = 1, nlayers

!  Average secant

         xa = avg_secants(n,L) ! solar average-secant for transmittance

!  Help variables

         deltaus = opdeps(n)
         dxv = deltaus * xv ; call ExpTrans(BIGEXP,dxv,vtrans)
         dxa = deltaus * xa ; call ExpTrans(BIGEXP,dxa,atrans)
         dxs = deltaus * xs
         hw  = omegas(n) * quarter

!  Solar to viewing (second Order). Tranmsittance/Multiplier

         sv_plus = xv + xa
         O2_AVTrans(n,L) = vtrans * atrans
         O2_AVMult (n,L) = xv * ( one - O2_AVTrans(n,L) ) / sv_plus

!  Loop over Nstreams (quadrature)
!  -------------------------------

         do k = 1, nstreams

!  Quadrature secant value

            xk  = one/qstreams(k)

!  Help variables

            dxk  = deltaus * xk ; call ExpTrans(BIGEXP,dxk,ktrans)
            xvxk = xk * xv ; kv_plus   = xk + xv ; kv_minus   = xk - xv
            xsxk = xk * xs ; ka_plus   = xk + xa ; ka_minus   = xk - xa

!  First Order: viewing direction + all quadrature directions
!    (Formerly, Facj, Chibj)

            O1_QVTrans(n,k,L) = ktrans * vtrans
            O1_QVMult (n,k,L) = xvxk * hw * ( one - O1_QVTrans(n,k,L) ) / kv_plus

!  First Order: Solar (Av-secant) direction + all quadrature directions
!    (Formerly, Faci, Chibi)

            O1_QATrans(n,k,L) = ktrans * atrans
            O1_QAMult (n,k,L) = xsxk * hw * ( one - O1_QATrans(n,k,L) ) / ka_plus

!  Second Order Multipliers (Formerly PHG, PHH, G, and H )
!     General case (No Small-number series)    

            O2_QVMult_d(n,k,L) = atrans * ( vtrans - ktrans ) / kv_minus
            O2_QAMult_d(n,k,L) = vtrans * ( atrans - ktrans ) / ka_minus
            O2_QAMult_u        = ( one - O1_QVTrans(n,k,L) )  / kv_plus
            O2_QAVMult_du(n,k,L) = xsxk * ( O2_AVMult(n,L) - xv * O2_QVMult_d(n,k,L) ) / ka_plus
            O2_QAVMult_dd(n,k,L) = xsxk * ( O2_AVMult(n,L) - xv * O2_QAMult_u        ) / ka_minus

!   Small-number series 
!    Debugged 05 February 2013. Taylor_3 works.

            if ( kv_minus .lt. TAYLOR_SMALL ) then
               arg1 = deltaus ; call Taylor_1(kv_minus,arg1,T1)
               O2_QVMult_d(n,k,L) = deltaus * atrans * vtrans * T1
            endif
            if ( ka_minus .lt. TAYLOR_SMALL ) then
               arg1 = deltaus ; call Taylor_1(ka_minus,arg1,T1)
               arg3 = O2_AVTrans(n,L) ; arg2 = one/sv_plus  ; call Taylor_3(ka_minus,arg1,arg2,arg3,T3)
               O2_QAMult_d  (n,k,L) = deltaus * vtrans * atrans * T1
               O2_QAVMult_dd(n,k,L) = xsxk * xv * T3 / sv_plus
            endif

!  end k loop

         enddo

!  end layer loop

      enddo

!  end geometry loop

   enddo

!  Finish

   return
end subroutine Calculate_Multipliers

subroutine Add_Fourier_Component &
        (m,nstokes, phi,R2c,R2s,R2cscal,R2,Icorr)

   implicit none

!  inputs

   integer  , intent(in) :: m,nstokes
   real(fpk), intent(in) :: phi,R2c(4),R2s(4),R2cscal

!  outputs

   real(fpk), intent(inout) :: R2(4),Icorr

!  local variables

   integer   :: ki
   real(fpk) :: cosmph,fac,sinmph,rm

!  start of code

   rm = real(m,fpk) * deg_to_rad ; cosmph = cos(rm*phi)
   if (nstokes .eq. 3)   sinmph = sin(rm*phi)
   fac = two ; if (m .eq. 0) fac = one

!  I and Q components

   do ki = 1, 2
      R2(ki) = R2(ki)+fac*R2c(ki)*cosmph
   enddo

!  U component

   if (nstokes .eq. 3) then
      R2(nstokes) = R2(nstokes)+ fac*R2s(nstokes)*sinmph
   endif

!  scalar correction

   Icorr = Icorr+fac*(R2c(1)-R2cscal)*cosmph

!  Finish

   return
end subroutine Add_Fourier_Component

subroutine Test_Fourier_Convergence &
    (verbo,m,nfoumax,nstokes,epsilon,muv,mus,almost,nextm,R2s,R2c)

   implicit none

!  inputs

   integer  , intent(in)    :: m,nfoumax,nstokes
   logical  , intent(in)    :: verbo
   real(fpk), intent(in)    :: epsilon,muv,mus

!  Modified inputs

!mick fix 4/21/2015 - removed saved attribte from "almost" - added now to argument list
   logical  , intent(inout) :: almost
   logical  , intent(inout) :: nextm
   real(fpk), intent(inout) :: R2s(4),R2c(4)

!  local variables

   integer :: i,count
   logical :: vertic

!  Initialize

   if (m .eq. 0) almost = .false.
   nextm = .true.

!  Stop Fourier series if you are beyond Foumax

   if (m .ge. nfoumax) then
      if(verbo)write(*,*)'Stop Fourier after m = ', m,' (nfoumax reached)'
      nextm = .false. ; return
   endif

!  Do not examine convergence for m < 2

   if ( m .lt. 2) then
      nextm = .true. ; return
   endif

!  Examine verticality condition

   vertic = .true.
   if ( (abs(muv - one) .gt. Smallnum6) .and. &
        (abs(mus - one) .gt. Smallnum6) ) vertic = .false.
   if (vertic) then
      if(verbo)write(*,*)'Stop Fourier after m = ', m,' (Verticality)'
      nextm = .false. ; return
   endif

!  General case, examine convergence

   nextm = .false.
   count = 0
   do i = 1, 2
      if ((mus*abs(R2c(i))) .gt. epsilon) count = count+1
   enddo
   do i = 3, nstokes
      if ((mus*abs(R2s(i))) .gt. epsilon) count = count+1
   enddo
   if (count .eq. nstokes) nextm = .true.

!  Almost condition

   if ((.not. nextm) .and. (.not. almost)) then
      nextm  = .true.;  almost = .true.
   else
      almost = .not. nextm
   endif
   if (verbo .and. .not. nextm)  &
      write(*,*)'Stop Fourier after m = ', m,' (for extra points)'

!  Done

  return
end subroutine Test_Fourier_Convergence

subroutine Set_avsecant &
     ( MaxGeoms, do_plane_parallel, nlayers, ngeoms, &
       geoms, SunChapman, opdeps, avg_secants )

   implicit none

!  inputs and output

   integer  , intent(in)  :: MaxGeoms
   logical  , intent(in)  :: do_plane_parallel
   integer  , intent(in)  :: nlayers,ngeoms

   real(fpk), intent(in)  :: geoms(3,4,MaxGeoms)
   real(fpk), intent(in)  :: opdeps(MAXLAYERS)
   real(fpk), intent(in)  :: SunChapman (MAXLAYERS,MAXLAYERS,MaxGeoms)
   real(fpk), intent(out) :: avg_secants(MAXLAYERS,MaxGeoms)

! local variables

   real(fpk) :: delta(MAXLAYERS), sum, sum1, xs
   integer   :: n,n1,m,L

   if ( do_plane_parallel ) then
      do L = 1, ngeoms
         xs = geoms(1,3,L) ; xs  = one / xs
         do n = 1, nlayers
            avg_secants(n,L) = xs
         enddo
      enddo
   else
      do n = 1, nlayers
         delta(nlayers+1-n) = opdeps(n)
      enddo
      do n = 1, nlayers
         n1 = nlayers + 1 - n
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
   endif

   return
end subroutine Set_avsecant

!  finish module

end module vlidort_2OScorr_routines


