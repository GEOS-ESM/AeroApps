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
! #                 Calculate_Zmatrix_LAPlus_ALT            #
! #                                                         #
! #      [Calculate_Zmatrix_LAPlus]   removed, 3/4/15       #
! #                                                         #
! ###########################################################

module vlidort_2OScorr_la_routines

   use VLIDORT_PARS, only : fpk, zero, one, two, four,              &
                            six, half, quarter, deg_to_rad,         &
                            BIGEXP, TAYLOR_SMALL, Smallnum6,        &
                            MAXLAYERS, MAXMOMENTS,                  &
                            MAXSTREAMS, MAXSTREAMS_P2, MAX_ATMOSWFS

   use vlidort_2OScorr_utilities, only : Make_Trans23, Make_Trans23_P

!  Note that Calculate_Zmatrix_LAPlus is not thread-safe, because of save statements
!       The alternative routine Calculate_Zmatrix_LAPlus_ALT is thread-safe

   implicit none

public  :: Calculate_Zmatrix_LAPlus_ALT

!  Removed
!           Calculate_Zmatrix_LAPlus 

contains


subroutine Calculate_Zmatrix_LAPlus_ALT &
        ( do_Jacobians, m, nstokes, nstreams, ns1, ns2, nmoments,  & ! Inputs
          npars, rl, r2Lp1, sql4, sqlm, binfac, coefsm, L_coefsm,  & ! Inputs
          xmu, xsi, Plm_2all, Plm_mlp, Plm_mlm,                    & ! Inputs
          Ptc,Pts,Prc,Prs, L_Ptc,L_Pts,L_Prc,L_Prs )

   implicit none

!  inputs
!  ======

!  Control

   logical  , intent(in) :: do_Jacobians
   integer  , intent(in) :: m,nstreams,nstokes,nmoments,npars

!  Coefficients

   real(fpk), intent(in) :: coefsm  (0:MAXMOMENTS,6)
   real(fpk), intent(in) :: L_coefsm(0:MAXMOMENTS,6,MAX_ATMOSWFS)

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
!  =======

!  Basic quantities

   real(fpk), intent(out) :: Ptc(2,MAXSTREAMS,4,4),Pts(2,MAXSTREAMS,4,4), &
                             Prc(2,MAXSTREAMS,4,4),Prs(2,MAXSTREAMS,4,4)

!  Linearized quantities

   real(fpk), intent(out) :: L_Ptc(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)
   real(fpk), intent(out) :: L_Pts(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)
   real(fpk), intent(out) :: L_Prc(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)
   real(fpk), intent(out) :: L_Prs(2,MAXSTREAMS,4,4,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  help variables

   integer   :: l,j,k1,k2,i,lold,lnew,itmp, p
   real(fpk) :: parity,u
   real(fpk) :: Zmpls_xv(nstreams,4,4),Zmmin_xv(nstreams,4,4)
   real(fpk) :: Zmpls_xs(nstreams,4,4),Zmmin_xs(nstreams,4,4)
   real(fpk) :: Plm(MAXSTREAMS_p2,3,2)
   real(fpk) :: DPDpl(MAXSTREAMS_p2,4),DPDmi(MAXSTREAMS_p2,4)
   real(fpk) :: DSD(4,4), SPj, f1new, f1old, tmp, f2new,f2newa,f2old

!  new local variables for linearisation

   real(fpk) :: L_DSD(4,4,MAX_ATMOSWFS)
   real(fpk) :: L_Zmpls_xv(nstreams,4,4,npars)
   real(fpk) :: L_Zmmin_xv(nstreams,4,4,npars)
   real(fpk) :: L_Zmpls_xs(nstreams,4,4,npars)
   real(fpk) :: L_Zmmin_xs(nstreams,4,4,npars)

!  Zero the transmittance and reflectance matrices

   Ptc   = zero ;   Pts = zero ;   Prc = zero ;   Prs = zero
   L_Ptc = zero ; L_Pts = zero ; L_Prc = zero ; L_Prs = zero

!  Zero local values

   Zmpls_xv   = zero ; Zmmin_xv   = zero
   Zmpls_xs   = zero ; Zmmin_xs   = zero
   L_Zmpls_xv = zero ; L_Zmmin_xv = zero
   L_Zmpls_xs = zero ; L_Zmmin_xs = zero

!  Start development of spherical functions

   lold = 1 ; lnew = 2
   do i = 1, ns2
      if (m .ne. 0) then
         Plm(i,1,lnew) = binfac*xsi(i)**m
      else
         Plm(i,1,lnew) = binfac
      endif
      Plm(i,1,lold) = zero
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
      DSD (2, 2) = quarter * ( coefsm (l, 2) + coefsm (l, 3) )
      DSD (3, 2) = quarter * ( coefsm (l, 2) - coefsm (l, 3) )
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

!  here is the linearised DSD

      if ( do_Jacobians ) then
         do p = 1, npars
           L_DSD (1, 1, p) = L_coefsm (l, 1, p)
           L_DSD (2, 1, p) = half * L_coefsm (l, 5, p)
           L_DSD (2, 2, p) = quarter * ( L_coefsm (l,2,p) + L_coefsm (l,3,p) )
           L_DSD (3, 2, p) = quarter * ( L_coefsm (l,2,p) - L_coefsm (l,3,p) )
           L_DSD (3, 1, p) = L_DSD (2, 1, p)
           L_DSD (1, 2, p) = L_DSD (2, 1, p)
           L_DSD (1, 3, p) = L_DSD (2, 1, p)
           L_DSD (2, 3, p) = L_DSD (3, 2, p)
           L_DSD (3, 3, p) = L_DSD (2, 2, p)
           L_DSD (1, 4, p) = zero
           L_DSD (2, 4, p) = half * L_coefsm (l, 6, p)
           L_DSD (3, 4, p) = - L_DSD (2, 4, p)
           L_DSD (4, 4, p) = L_coefsm (l, 4, p)
           L_DSD (4, 1, p) = zero
           L_DSD (4, 2, p) = - L_DSD (2, 4, p)
           L_DSD (4, 3, p) = - L_DSD (3, 4, p)
         enddo
      endif

!  now make the Z matrices for scattering: Polar to View streams

      do k2 = 1, 4
         do k1 = 1, 4
            do j = 1, nstreams
               SPj = DSD(k1,k2) * DPDpl(j,k2)
               Zmpls_xv(j,k1,k2) = Zmpls_xv(j,k1,k2)+DPDpl(ns2,k1)*SPj
               Zmmin_xv(j,k1,k2) = Zmmin_xv(j,k1,k2)+DPDmi(ns2,k1)*SPj
            enddo
            if ( do_Jacobians ) then
               do j = 1, nstreams
                  do p = 1, npars
                     SPj = L_DSD(k1,k2,p)*DPDpl(j,k2)
                     L_Zmpls_xv(j,k1,k2,p) = L_Zmpls_xv(j,k1,k2,p)+DPDpl(ns2,k1)*SPj
                     L_Zmmin_xv(j,k1,k2,p) = L_Zmmin_xv(j,k1,k2,p)+DPDmi(ns2,k1)*SPj
                  enddo
               enddo
            endif
         enddo
      enddo

!  now make the Z matrices for scattering: Sun to Polar streams

      do k1 = 1, 4
         SPj = DSD(k1,1) * DPDpl(ns1,1)
         do i = 1, nstreams
            Zmpls_xs(i,k1,1) = Zmpls_xs(i,k1,1)+DPDpl(i,k1)*SPj
            Zmmin_xs(i,k1,1) = Zmmin_xs(i,k1,1)+DPDmi(i,k1)*SPj
         enddo
         if ( do_Jacobians ) then
            do i = 1, nstreams
               do p = 1, npars
                  SPj = L_DSD(k1,1,p) * DPDpl(ns1,1)
                  L_Zmpls_xs(i,k1,1,p) = L_Zmpls_xs(i,k1,1,p)+DPDpl(i,k1)*SPj
                  L_Zmmin_xs(i,k1,1,p) = L_Zmmin_xs(i,k1,1,p)+DPDmi(i,k1)*SPj
               enddo
            enddo
         endif
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
   if ( do_Jacobians ) then
      do p = 1, npars
         do j = 1, nstreams
            call Make_Trans23_P(nstreams,npars,j,p,L_Zmmin_xv)
            call Make_Trans23_P(nstreams,npars,j,p,L_Zmpls_xv)
            call Make_Trans23_P(nstreams,npars,j,p,L_Zmmin_xs)
            call Make_Trans23_P(nstreams,npars,j,p,L_Zmpls_xs)
         enddo
      enddo
   endif  

!      if ( layer.eq.30)write(*,*)m,l,Zmmin_xv(8,2,3),Zmmin_xv(8,3,2)
!      if ( m.eq.3.and.layer.eq.30)pause

!  Set results for Polar incoming - View angle outgoing

   do j = 1, nstreams
      do k1 = 1, 2
         Ptc(1,j,k1,1:2) = Zmpls_xv(j,k1,1:2)
         Prc(1,j,k1,1:2) = Zmmin_xv(j,k1,1:2)
         if ( do_Jacobians ) then
            do p = 1, npars
               L_Ptc(1,j,k1,1:2,p) = L_Zmpls_xv(j,k1,1:2,p)
               L_Prc(1,j,k1,1:2,p) = L_Zmmin_xv(j,k1,1:2,p)
            enddo
         endif
         if (m .gt. 0) then
            Pts(1,j,k1,3:4) = -Zmpls_xv(j,k1,3:4)
            Prs(1,j,k1,3:4) = -Zmmin_xv(j,k1,3:4)
            if ( do_Jacobians ) then
               do p = 1, npars
                  L_Pts(1,j,k1,3:4,p) = -L_Zmpls_xv(j,k1,3:4,p)
                  L_Prs(1,j,k1,3:4,p) = -L_Zmmin_xv(j,k1,3:4,p)
               enddo
            endif
         endif
      enddo
      do k1 = 3, nstokes
         Ptc(1,j,k1,3:4) = Zmpls_xv(j,k1,3:4)
         Prc(1,j,k1,3:4) = Zmmin_xv(j,k1,3:4)
         if ( do_Jacobians ) then
            do p = 1, npars
               L_Ptc(1,j,k1,3:4,p) = L_Zmpls_xv(j,k1,3:4,p)
               L_Prc(1,j,k1,3:4,p) = L_Zmmin_xv(j,k1,3:4,p)
            enddo
         endif
         if (m .gt. 0) then
            Pts(1,j,k1,1:2) = Zmpls_xv(j,k1,1:2)
            Prs(1,j,k1,1:2) = Zmmin_xv(j,k1,1:2)
            if ( do_Jacobians ) then
               do p = 1, npars
                  L_Pts(1,j,k1,1:2,p) = L_Zmpls_xv(j,k1,1:2,p)
                  L_Prs(1,j,k1,1:2,p) = L_Zmmin_xv(j,k1,1:2,p)
               enddo
            endif
         endif
      enddo
   enddo

!  Set results for Solar incoming - Polar angle outgoing

   do i = 1, nstreams
      Ptc(2,i,1:2,1) = Zmpls_xs(i,1:2,1)
      Prc(2,i,1:2,1) = Zmmin_xs(i,1:2,1)
      if ( do_Jacobians ) then
         do p = 1, npars
            L_Ptc(2,i,1:2,1,p) = L_Zmpls_xs(i,1:2,1,p)
            L_Prc(2,i,1:2,1,p) = L_Zmmin_xs(i,1:2,1,p)
         enddo
      endif
      if (m .gt. 0) then
         Pts(2,i,3:4,1) = Zmpls_xs(i,3:4,1)
         Prs(2,i,3:4,1) = Zmmin_xs(i,3:4,1)
         if ( do_Jacobians ) then
            do p = 1, npars
               L_Pts(2,i,3:4,1,p) = L_Zmpls_xs(i,3:4,1,p)
               L_Prs(2,i,3:4,1,p) = L_Zmmin_xs(i,3:4,1,p)
            enddo
         endif
      endif
   enddo

!  Finish routine

   return
end subroutine Calculate_Zmatrix_LAPlus_ALT

!  Finish module

end module vlidort_2OScorr_la_routines
