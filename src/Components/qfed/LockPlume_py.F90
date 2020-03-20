!
! Python front end to Lock Plume Rise
!

#define R_(a) reshape(a,(/1,km/))

   subroutine plumeg(km, jm, im, u, v, t, q, delp, ptop, b_star, u_star, ztop )

      use LockPlume_Mod, only: mpbl_depth, derived
      implicit NONE

      integer, intent(in)  :: km, im, jm      ! dimensionms
      real,    intent(in)  :: u(km,jm,im)     ! zonal wind (m/s)
      real,    intent(in)  :: v(km,jm,im)     ! meridional wind (m/s)
      real,    intent(in)  :: t(km,jm,im)     ! temperature (K)
      real,    intent(in)  :: q(km,jm,im)     ! pressure
      real,    intent(in)  :: delp(km,jm,im)  ! pressure thickness (Pa)
      real,    intent(in)  :: ptop            ! pressure at top (Pa)
      real,    intent(in)  :: u_star(jm,im)   ! surface velocity scale (m/s)
      real,    intent(in)  :: b_star(jm,im)   ! buoyancy scale (m s-2)

      real,    intent(out) :: ztop(jm,im)     ! plume height (m)

!
!       Computes Lock plume top given transposed arrays on input.
!

!                             ---

      integer :: i, j

      do i = 1, im
         do j = 1, jm
            call plume(km, u(:,j,i), v(:,j,i), t(:,j,i), q(:,j,i), delp(:,j,i), &
                       ptop, b_star(j,i), u_star(j,i), ztop(j,i) )
         end do
      end do

    end subroutine plumeg

   subroutine plume(km,u, v, t, q, delp, ptop, b_star, u_star, ztop )

      use LockPlume_Mod, only: mpbl_depth, derived
      implicit NONE

      integer, intent(in)  :: km
      real,    intent(in)  :: u(km)    ! zonal wind (m/s)
      real,    intent(in)  :: v(km)    ! meridional wind (m/s)
      real,    intent(in)  :: t(km)    ! temperature (K)
      real,    intent(in)  :: q(km)    ! pressure
      real,    intent(in)  :: delp(km) ! pressure thickness (Pa)
      real,    intent(in)  :: ptop     ! pressure at top (Pa)
      real,    intent(in)  :: u_star   ! surface velocity scale (m/s)
      real,    intent(in)  :: b_star   ! buoyancy scale (m s-2)

      real,    intent(out) :: ztop     ! plume height (m)

!                             ---

      real :: z(km), p(km), ztop_(1)
      real :: tpfac, entrate, pceff
      integer :: ipbl, n, ncol=1, i=1

!     Hardwire these for now
!     ----------------------
      tpfac = 1.0 ! GEOS uses 20 by default
      entrate = 0.2/200.
      pceff = 0.5

!     Vertical coordinates
!     --------------------
      call Derived( p, z, km, T, q, delp, ptop )

!     Go for it
!     ---------
      call mpbl_depth(i, ncol, km, tpfac, entrate, pceff, &
                      R_(t), R_(q), R_(u), R_(v), R_(z), R_(p), &
                      (/b_star/), (/u_star/), ipbl, ztop_ )

      ztop = ztop_(1)
      
    end subroutine plume

