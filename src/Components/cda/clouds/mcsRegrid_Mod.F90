Module mcsRegrid_Mod

  use mod_consts ! to be replaced with MAPL_Constants!!!

  implicit NONE

!...................................................................................

Contains

  subroutine getEdgeCoords1 ( pe, he, km, T, q, delp, ptop )

   implicit NONE

   integer, intent(in) :: km                 ! No. vertical layers
   real,    intent(in) :: T(km)              ! Dry temperature [K]
   real,    intent(in) :: q(km)              ! Specific humidity [kg/kg]
   real,    intent(in) :: delp(km)           ! layer pressure thicjness [Pa]
   real,    intent(in) :: ptop               ! top (edge) pressure [Pa]

   real, intent(out)   :: pe(km+1)            ! pressure at edges
   real, intent(out)   :: he(km+1)            ! height (above sfc) at edges
 
!
! ------------------------------- pe(k),   he(k)
!
! ............................... pm(k),   hm(k), delp(k)
!
! ------------------------------- pe(k+1), he(k+1)
!

   integer :: k
   real :: Tv
   real :: mixr(km), delh(km)

!  Mixing ratio from specific humidity
!  -----------------------------------
   mixr = q / ( 1.0 - q )

!  Construct edge pressures
!  ------------------------
   pe(1) = ptop
   do k = 1, km
      pe(k+1) = pe(k) + delp(k)
   end do

!  Construct mid-layer pressures and layer thickness
!  -------------------------------------------------
   do k = 1, km
      Tv = T(k) * ( 1 + 0.61 * mixr(k) )
      delh(k)  = Rgas * Tv * log(pe(k+1)/pe(k)) / grav
   end do

!  Compute Geo-potential height at edges and mid-layer
!  ---------------------------------------------------
   he(km+1) = 0.0  ! set the surface to zero height
   do k = km, 1, -1
      he(k) = he(k+1) + delh(k)
   end do

 end subroutine getEdgeCoords1

!---

  subroutine getTerrainPressure1 ( ps, zs, km, pe, ze, k_off ) 

  integer, intent(in) :: km        ! number of vertical levels
  real,    intent(in) :: pe(km+1)  ! model edge pressure [Pa]
  real,    intent(in) :: ze(km+1)  ! model edge height (above sea level) [m]
  real,    intent(in) :: zs        ! terrain height [m]

  integer, intent(in), optional :: k_off ! For robustness use mean theta in the
                                         !  first 4 layers abobe the surface
  
  real, intent(out)   :: ps        ! terrain surface pressure [m]

!
! Hydrostatic reduction of model surface pressure to the real
! surface pressure at terrain height zs. The mean virtual potential
! temperature of the first k_off  model layer is used for this reduction
! via the hydrostatic relation:
!
!          d z = - (cp/g) theta d ( p**kappa)
!
! which gives the following similarity relationship:
!
!         z2 - z1       cp           zs - z1
!       ----------- = - -- theta = ----------- 
!       p2^k - p1^k      g         ps^k - p1^k
!
! or
!                                    zs - z1
!       ps^k - p1^k =  (p2^k - p1^k) -------
!                                    z2 - z1
!

   real(kind=8) :: z1, z2, pk1, pk2 ! Do calculation with real*8
   integer :: k1, k2, koff

   if (present(k_off) ) then
       koff = k_off
   else
       koff = 4
   end if

   k1 = km + 1;     k2  = k1 - koff
   z1 = ze(k1);     pk1 = pe(k1)**kappa
   z2 = ze(k2);     pk2 = pe(k2)**kappa

   ps = (pk1 + (pk2 - pk1)*(zs-z1)/(z2-z1))**(1./kappa)

  end subroutine getTerrainPressure1

!---

  subroutine interpPressure1 ( di_pe, nk, ps, k_surface, di_ze, zs, g5_pe, g5_ze, km, rc ) 

     integer, intent(in)  :: nk           ! number of DISORT layers
     integer, intent(in)  :: km           ! number of GEOS-5 layers

     real, intent(in)     :: zs           ! Terrain height
     real, intent(in)     :: di_ze(nk+1)  ! DISORT edge height
     real, intent(in)     :: g5_pe(km+1)  ! GEOS-5 edge pressure
     real, intent(in)     :: g5_ze(km+1)  ! GEOS-5 edge height
 
     real, intent(out)    :: di_pe(nk+1)  ! DISORT edge pressure
     real, intent(out)    :: ps           ! Terrain modified surface pressure
     integer, intent(out) :: k_surface    ! surface level 
     integer, intent(out) :: rc

!                          ----
      real    :: alpha, z, g5_logP(km+1), di_logP(nk+1)
      logical :: got_it
      integer :: k, m, n, k_last

      rc = 0

!     Determine the last DISORT index above terrain surface
!     -----------------------------------------------------
      do k = 1, nk+1
         k_surface = k
         if ( di_ze(k) <= zs ) exit
      end do

!     Make sure GEOS-5 top is higher than DISORT's
!     --------------------------------------------
      g5_logP = log(g5_pe)
      if ( di_ze(1) > g5_ze(1) ) then
           rc = 1
           return
      end if

!     Linear interpolation in log-P
!     -----------------------------
      n = 1
      k_loop: do k = 1, k_surface
         z = di_ze(k)
         got_it = .FALSE.
         m_loop: do m = n, km
            if ( (z<=g5_ze(m)) .AND. (z>g5_ze(m+1)) ) then
               n = m
               got_it = .TRUE.
               exit m_loop
            end if
         end do m_loop
         if ( got_it ) then
            alpha = (z - g5_ze(n+1)) / (g5_ze(n)-g5_ze(n+1))
            di_logP(k) = alpha*g5_logP(n) + (1-alpha) * g5_logP(n+1)
            k_last = k
         else
            exit k_loop
         end if
      end do k_loop

 !    Interpolated pressure levels
 !    ----------------------------
      di_pe = UNDEF
      di_pe(1:k_last) = exp(di_logP(1:k_last))

!     Hydrostatically, extrapolate below model surface
!     ------------------------------------------------
      if ( k_surface > k_last ) then
         do k = k_last+1, k_surface-1
            call getTerrainPressure1 ( di_pe(k), di_ze(k), km, g5_pe, g5_ze ) 
         end do
      end if

!     Terrain pressure at surface, always
!     -----------------------------------
      call getTerrainPressure1 ( ps, zs, km, g5_pe, g5_ze ) 

      di_pe(k_surface) = ps

    end subroutine interpPressure1

!---

    subroutine edgeT1 ( di_Te, nk, ks, di_pe, g5_pe, g5_Tm, g5_Tskin, km, rc ) 

        implicit NONE

        integer, intent(in)  :: nk   ! number of DISORT layers
        integer, intent(in)  :: ks   ! DISORT surface level
        integer, intent(in)  :: km   ! number of GEOS-5 layers
  
        real, intent(in)     :: di_pe(nk+1)
        real, intent(in)     :: g5_pe(km+1)
        real, intent(in)     :: g5_Tm(km)
        real, intent(in)     :: g5_Tskin

        real, intent(out)    :: di_Te(nk+1)
        integer, intent(out) :: rc

!                              ---

        integer :: k, iv, kord 
        real :: di_logPm(nk), g5_logPe(km+1) 

        rc = 0

!       DISORT mid-level pressure
!       -------------------------
        do k = 1, nk
           di_logPm(k) = log(0.5 * (di_pe(k) + di_pe(k+1)))
        end do
        g5_logPe = log(g5_pe)

!       Compute Edge Pressure for k = 2 thru nk using map_ppm
!       -----------------------------------------------------
        iv = 1
        kord = 3
        call map_ppm ( di_Te(2:nk), nk-1, di_logPm, g5_logPe, g5_Tm, km, iv, kord )

!       Upper boundary - isotherm extrapolation
!       ---------------------------------------
        di_Te(1) = di_Te(2)  

!       Lower boundary - skin temperature
!       ---------------------------------
        di_Te(ks) = g5_Tskin

!       UNDEFs below surface
!       --------------------
        do k = ks+1, nk+1
           di_Te(k) = UNDEF
        end do

      end subroutine edgeT1

!---

    subroutine edgeQ1 ( di_Qe, nk, ks, di_pe, g5_pe, g5_Qm, km) 

        implicit NONE

        integer, intent(in)  :: nk   ! number of DISORT layers
        integer, intent(in)  :: ks   ! DISORT surface level
        integer, intent(in)  :: km   ! number of GEOS-5 layers
  
        real, intent(in)     :: di_pe(nk+1)
        real, intent(in)     :: g5_pe(km+1)
        real, intent(in)     :: g5_Qm(km)

        real,    intent(out) :: di_Qe(nk+1)

!                              ---

        integer :: k, iv, kord 
        real :: di_logPm(nk), g5_logPe(km+1) 

!       DISORT mid-level pressure
!       -------------------------
        do k = 1, nk
           di_logPm(k) = log(0.5 * (di_pe(k) + di_pe(k+1)))
        end do
        g5_logPe = log(g5_pe)

!       Remap to edge for k = 2 thru nk using map_ppm
!       ---------------------------------------------
        iv = 0
        kord = 3
        call map_ppm ( di_Qe(2:nk), nk-1, di_logPm, g5_logPe, g5_Qm, km, iv, kord )

!       Upper boundary - simple extrapolation
!       -------------------------------------
        di_Qe(1) = di_Qe(2)

!       Lower boundary - lowest GEOS-5 layer value
!       ------------------------------------------
        di_Qe(ks) = g5_Qm(km)

!       UNDEFs below surface
!       --------------------
        do k = ks+1, nk+1
           di_Qe(k) = UNDEF
        end do

      end subroutine edgeQ1

!---
      subroutine map_ppm ( di_q, nk, di_pe, g5_pe, g5_q, km, iv, kord )

        implicit NONE
        integer, intent(in) :: nk   ! number of DISORT layers
        integer, intent(in) :: km   ! number of GEOS-5 layers
        integer, intent(in) :: iv   ! use 0 for tracers
        integer, intent(in) :: kord ! algorithm order
  
        real, intent(in)    :: di_pe(nk+1)
        real, intent(in)    :: g5_pe(km+1)
        real, intent(in)    :: g5_q(km)

        real, intent(out)   :: di_q(nk)

!                              ---

      real pe1(1,km+1)
      real pe2(1,nk+1)
      real  q1(1,1,km)
      real  q2(1,1,nk)

      ! Copy in...
      pe1(1,:) = g5_pe(:)
      pe2(1,:) = di_pe(:)
      q1(1,1,:) = g5_q(:)

      call map_ppm_ ( km, pe1, q1, &
                      nk, pe2, q2, &
                      1, 1, 1, 1, 1, 1, iv, kord, UNDEF)
      
      ! Copy out...
      di_q(:) = q2(1,1,:)  

      end subroutine map_ppm

!---
      subroutine map_ppm_ ( km,   pe1,   q1, &
                            kn,   pe2,   q2, &
                            im, i1, i2, j, jfirst, jlast, iv, kord, undef)

! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate

      implicit none

      real       r3, r23
      parameter (r3 = 1./3., r23 = 2./3.)

! Input:
      integer i1, i2
      integer im, jfirst, jlast, iv, kord
      integer km                             ! Original vertical dimension
      integer kn                             ! Target vertical dimension

      real pe1(im,km+1)
      real pe2(im,kn+1)
      real  q1(im,jfirst:jlast,km)
      real undef
! Output
      real  q2(im,jfirst:jlast,kn)

! Local
      real   dp1(i1:i2,km)
      real  q4(4,i1:i2,km)

      integer i, j, k, l, ll, k0
      real    pl, pr, qsum, delp, esl

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            q4(1,i,k) = q1(i,j,k)
         enddo
      enddo

! Compute vertical subgrid distribution
      call ppm2m( q4, dp1, km, i1, i2, iv, kord )

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) .lt. pe1(i,1)) then
!        q2(i,j,k) = q4(2,i,1)
         q2(i,j,k) = undef
         goto 555
      elseif(pe2(i,k) .ge. pe1(i,l) .and. pe2(i,k) .le. pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) .le. pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l)) &
                         *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+ &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)* &
                     (r3*(1.+pl*(1.+pl))))
              do ll=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) .gt. pe1(i,ll+1) ) then
! Whole layer..
                 qsum = qsum + dp1(i,ll)*q4(1,i,ll)
                 else
                 delp = pe2(i,k+1)-pe1(i,ll)
                  esl = delp / dp1(i,ll)
                 qsum = qsum + delp*(q4(2,i,ll)+0.5*esl* &
                       (q4(3,i,ll)-q4(2,i,ll)+q4(4,i,ll)*(1.-r23*esl)))
                 k0 = ll
                 goto 123
                 endif
              enddo
              goto 123
           endif
      endif
100   continue
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
1000  continue

      return
    end subroutine map_ppm_

      subroutine ppm2m(a4, delp, km, i1, i2, iv, kord)

! iv = 0: positive definite scalars
! iv = 1: others

      implicit none
! Input
      integer km, lmt, iv
      integer i1, i2
      integer kord
      real    delp(i1:i2,km)
      real    a4(4,i1:i2,km)

! local arrays.
      real   dc(i1:i2,km)
      real   h2(i1:i2,km)
      real delq(i1:i2,km)

      real a1, a2, c1, c2, c3, d1, d2
      real qmax, qmin, cmax, cmin
      real qm, dq, tmp

! Local scalars:
      integer i, k, km1
      real qmp
      real lac
      integer it

      km1 = km - 1
       it = i2 - i1 + 1

      do k=2,km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
            a4(4,i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo

      do k=2,km1
         do i=i1,i2
            c1 = (delp(i,k-1)+0.5*delp(i,k))/a4(4,i,k+1)
            c2 = (delp(i,k+1)+0.5*delp(i,k))/a4(4,i,k)
            tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) / &
                                    (a4(4,i,k)+delp(i,k+1))
            qmax = max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1)) - a4(1,i,k)
            qmin = a4(1,i,k) - min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))
            dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp)
         enddo
      enddo

!****6***0*********0*********0*********0*********0*********0**********72
! 4th order interpolation of the provisional cell edge value
!****6***0*********0*********0*********0*********0*********0**********72

      do k=3,km1
      do i=i1,i2
      c1 = delq(i,k-1)*delp(i,k-1) / a4(4,i,k)
      a1 = a4(4,i,k-1) / (a4(4,i,k) + delp(i,k-1))
      a2 = a4(4,i,k+1) / (a4(4,i,k) + delp(i,k))
      a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(a4(4,i,k-1)+a4(4,i,k+1)) * &
                  ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) - &
                               delp(i,k-1)*a1*dc(i,k  ) )
      enddo
      enddo

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      do i=i1,i2
      d1 = delp(i,1)
      d2 = delp(i,2)
      qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
      dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
      c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
      c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1**2)
      a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
      a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
      dc(i,1) =  a4(1,i,1) - a4(2,i,1)
! No over- and undershoot condition
      cmax = max(a4(1,i,1), a4(1,i,2))
      cmin = min(a4(1,i,1), a4(1,i,2))
      a4(2,i,2) = max(cmin,a4(2,i,2))
      a4(2,i,2) = min(cmax,a4(2,i,2))
      enddo

      if(iv .eq. 0) then
         do i=i1,i2
            a4(2,i,1) = max(0.,a4(2,i,1))
            a4(2,i,2) = max(0.,a4(2,i,2))
         enddo
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1**2)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
         dc(i,km) = a4(3,i,km) -  a4(1,i,km)
!****6***0*********0*********0*********0*********0*********0**********72
! No over- and undershoot condition
         cmax = max(a4(1,i,km), a4(1,i,km1))
         cmin = min(a4(1,i,km), a4(1,i,km1))
         a4(2,i,km) = max(cmin,a4(2,i,km))
         a4(2,i,km) = min(cmax,a4(2,i,km))
!****6***0*********0*********0*********0*********0*********0**********72
      enddo

      if(iv .eq. 0) then
      do i=i1,i2
         a4(2,i,km) = max(0.,a4(2,i,km))
         a4(3,i,km) = max(0.,a4(3,i,km))
      enddo
      endif

      do k=1,km1
         do i=i1,i2
            a4(3,i,k) = a4(2,i,k+1)
         enddo
      enddo
 
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
 
! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      if(kord .eq. 7) then
!****6***0*********0*********0*********0*********0*********0**********72
! Huynh's 2nd constraint
!****6***0*********0*********0*********0*********0*********0**********72
      do k=2, km1
         do i=i1,i2
            h2(i,k) = delq(i,k) - delq(i,k-1)
         enddo
      enddo

      do 1000 k=3, km-2
        do i=i1,i2
! Right edges
         qmp   = a4(1,i,k)                 + 2.0*delq(i,k-1)
         lac   = a4(1,i,k) + 1.5*h2(i,k-1) + 0.5*delq(i,k-1)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(3,i,k) = min(max(a4(3,i,k), qmin), qmax)
! Left  edges
         qmp   = a4(1,i,k)                 - 2.0*delq(i,k)
         lac   = a4(1,i,k) + 1.5*h2(i,k+1) - 0.5*delq(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(2,i,k) = min(max(a4(2,i,k), qmin), qmax)
! Recompute A6
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to prevent negatives
         if (iv .eq. 0) then
             call kmppm(dc(i1,k), a4(1,i1,k), it, 2)
         endif
1000  continue

      else
 
         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv .eq. 0) lmt = min(2, lmt)

      do k=3, km-2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm(dc(i1,k), a4(1,i1,k), it, lmt)
      enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      return
      end subroutine ppm2m

      subroutine kmppm(dm, a4, im, lmt)
      implicit none

      real       r12
      parameter (r12 = 1./12.)

      integer im, lmt
      real    a4(4,*)
      real    dm(*)

      integer i
      real da1, da2, a6da
      real fmin

      if ( lmt .eq. 3 ) return
! Full constraint

      if(lmt .eq. 0) then
      do i=1,im
      if(dm(i) .eq. 0.) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da .lt. -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da .gt. da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt .eq. 2) then
! Positive definite constraint

      do i=1,im
      if( abs(a4(3,i)-a4(2,i)) .lt. -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
      if( fmin .lt. 0. ) then
      if(a4(1,i).lt.a4(3,i) .and. a4(1,i).lt.a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
      elseif(a4(3,i) .gt. a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
      else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
      endif
      endif
      endif
      enddo

      elseif (lmt .eq. 1) then
              write(6,*) 'semi-monotonic not implemented!'
              stop
      endif

      return
      end subroutine kmppm

    end Module mcsRegrid_Mod
