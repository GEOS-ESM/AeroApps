!
! Read GEOS-4 files with gridded T, H, etc and call the non-ESMF version
! portion of the PlumeRise. This is a stand alone driver for testing 
! this module.
!

   Program ut_PlumeRise

   use m_set_eta
   use m_STrTemplate
   use m_die

   use rconstants, only: cp, rgas, g, kappa=>rocp, p00
   use PlumeRise_Mod

   implicit NONE

   type(PlumeRise) :: plume

   integer, parameter :: im=288, jm=181, km=55

   real, dimension(im,jm   )   :: ps, phis, lwi
   real, dimension(im,jm,km)   :: THETA, Q, H, RHO, PM
   real, dimension(im,jm,km+1) :: HE, PE, PEK

   real, dimension(im,jm,N_BIOME) :: fire_areas, z1, z2

   integer :: ia, iz, ja, jz, i, j   

   real, dimension(km+1) :: ak, bk

   integer :: nymd = 20061001, nhms = 120000
   integer :: ks, k, n
   real    :: ptop, pint

   character(len=*), parameter :: phis_fn = '/home/dasilva/out/plume/phis.hdf'
   character(len=*), parameter :: diag_fn = '/home/dasilva/out/plume/diag.hdf'

   real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0
   real, parameter :: ha = 10000. ! 1 hectare = 10,000 meter^2

   logical :: sfc2top

   integer ::  is = 241, js = 91


!                            -----


!  Read gridded fields
!  -------------------
   call Read_ ( phis_fn, 'phis', nymd, 0, im, jm, 1, var2d=phis )
   call Read_ ( diag_fn, 'oro', nymd, nhms, im, jm, 1, var2d=lwi )
   call Read_ ( diag_fn, 'surfp', nymd, nhms, im, jm, 1, var2d=ps )
   call Read_ ( diag_fn, 'q', nymd, nhms, im, jm, km, var3d=Q )
   call Read_ ( diag_fn, 'hghte', nymd, nhms, im, jm, km, var3d=H )
   sfc2top = ( h(1,1,km) .gt.  h(1,1,km-1) ) 
   if ( sfc2top ) then
     HE(:,:,2:km+1) = H               ! H is workspace here
     HE(:,:,1) = phis / g
   else
     HE(:,:,1:km) = H               ! H is workspace here
     HE(:,:,km+1) = phis / g
   end if
   call Read_ ( diag_fn, 'h', nymd, nhms, im, jm, km, var3d=H )

   q  = q / ( 1. - q )  ! specific humidity => mixing ratio

!  Get eta grid coeffs
!  -------------------
   call set_eta ( km, ks, ptop, pint, ak, bk )

   sfc2top = ( h(1,1,km) .gt.  h(1,1,km-1) ) 
   if ( sfc2top ) then
        ak(1:km+1) = ak(km+1:1:-1)
        bk(1:km+1) = bk(km+1:1:-1)
   end if

!  Construct pressures
!  -------------------
   do k = 1, km + 1
      pe(:,:,k) = ak(k) + bk(k) * ps
   end do
   do k = 1, km
      pm(:,:,k) = 0.5 * ( pe(:,:,k) + pe(:,:,k+1) )
   end do

!  Compute theta (TO DO: remove virtual effects)
!  ---------------------------------------------
   pek = (pe/p00) ** kappa
   do k = 1, km
      theta(:,:,k) = - (g/cp) * (  he(:,:,k+1) -  he(:,:,k) ) &
                              / ( pek(:,:,k+1) - pek(:,:,k) ) 
        rho(:,:,k) = - (1/g)  * (  pe(:,:,k+1) -  pe(:,:,k) ) &
                              / (  he(:,:,k+1) -  he(:,:,k) ) 
   end do

!  Sample output over amazon
!  -------------------------
!           '    222.5     1.5    29.5     0.0    11.7'
   print *, '    H(m)    p (hPa) THETA(K)  rho    q(g/kg)'
   i = is; j = js
   do k = km, 1, -1
      print 10, h(i,j,k), pm(i,j,k)/100, THETA(i,j,k), &
               rho(i,j,k), 1000*q(i,j,k)
   end do
     print 10, he(i,j,1)
10 format(1x,5f8.2)

!  Set fire sizes only over land gridpoints
!  ----------------------------------------
   do n = 1, N_BIOME
      where(lwi == LAND)
            fire_areas(:,:,n) = 20. * ha  ! in meter^2
      elsewhere
            fire_areas(:,:,n) = 0.0
      end where
   end do

!  Stats
!  -----
   print *, 'phis = ', minval(phis), maxval(phis)
   print *, '  ps = ', minval(ps),   maxval(ps)
   print *, ' oro = ', minval(lwi),  maxval(lwi)
   print *, 'Theta= ', minval(THETA),    maxval(THETA)
   print *, '   q = ', minval(Q),    maxval(Q)
   print *, '   h = ', minval(H),    maxval(H)
   print *, '  he = ', minval(HE),   maxval(HE)
   print *, '  pm = ', minval(PM),   maxval(PM)
   print *, ' rho = ', minval(RHO),  maxval(RHO)
   do n = 1, N_BIOME
   print *, 'area = ', minval(fire_areas(:,:,n)),  maxval(fire_areas(:,:,n)), &
                       nint(100.*count(fire_areas(:,:,n)/=0.0)/(im*jm)), '%'
   end do


!  Calculate plume rise
!  --------------------
   ia = is; iz = is; ja = js; jz =js
   call PlumeRise_Run ( plume, ia, iz, ja, jz, km, &
                        theta(ia:iz,ja:jz,:),  rho(ia:iz,ja:jz,:), &
                        pm(ia:iz,ja:jz,:), q(ia:iz,ja:jz,:),  &
                        h(ia:iz,ja:jz,:),  he(ia:iz,ja:jz,:),  &
                        fire_areas(ia:iz,ja:jz,:),            &
                        z1(ia:iz,ja:jz,:), z2(ia:iz,ja:jz,:) )

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!-BOP
!
! !IROUTINE:  Read_ --- Reads fields from file and distribute
!
! !INTERFACE:
!
   subroutine Read_ ( filen, varn, nymd, nhms, im, jm, km, &
                      var2d, var3d, cyclic )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

  character(len=*), intent(in) :: filen   ! GFIO compatible file name
  character(len=*), intent(in) :: varn    ! variable name
  integer, intent(in)          :: nymd, nhms ! date/time

                                          ! Distributed grid info:
  integer,      intent(in)     :: im      !   global zonal dimension
  integer,      intent(in)     :: jm      !   global zonal dimension
  integer,      intent(in)     :: km      !   vertical dimension


! !OUTPUT PARAMETERS:

  real, OPTIONAL, intent(out)   :: var2d(:,:)
  real, OPTIONAL, intent(out)   :: var3d(:,:,:)
  logical, OPTIONAL, intent(in) :: cyclic ! whether time dimension is periodic

! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
!  28oct2003 da Silva  First crack.
!  03ayg2004 da Silva  Uses GetVarT for time interpolation
!  18nov2004 da Silva  Added cyclic option for climatological files.
!  30may2005 da Silva  Introduced template expansing in file names.
!
!-EOP
!-------------------------------------------------------------------------

    character(len=*), parameter ::  myname = 'Read_'
    logical :: tcyclic

    integer  :: READ_ONLY=1, nokm=0
    real, allocatable :: v2d(:,:), v3d(:,:,:)
    integer :: fid, rc, ios
    character(len=255) :: fname

!   Consistency check
!   -----------------
    if ( .not. ( present(var2d) .or. present(var3d) ) ) then
       call die ( myname, 'missing var2d or var3d' )
    else if ( present(var2d) .and. present(var3d) ) then
       call die ( myname, 'either var2d or var3d, but not both' )
    end if

    if ( present(cyclic) ) then
        tcyclic = cyclic
    else
        tcyclic = .false. ! by default time dimension is not periodic
    end if


!   Expand templates
!   ----------------
    if ( index(filen,'%') .gt. 0 ) then
         call StrTemplate ( fname, filen, xid='unknown', &
                            nymd=nymd, nhms=nhms )
    else
         fname = filen
    end if


!   Read file
!   ---------

!       print *, myname // ': Opening GFIO file ' // trim(fname)

!      Open file, get first time
!      -------------------------
       call GFIO_Open ( fname, READ_ONLY, fid, rc )
       if ( rc .ne. 0 ) then
          call die(myname,'cannot open GFIO file '//trim(fname))
       end if

!      Read global array
!      -----------------
       print *, myname // ': reading ' // trim(varn) // &
                 ' from ' // trim(fname)
       if ( present(var2d) ) then
            call GFIO_GetVarT1 ( fid, trim(varn), nymd, nhms, im, jm, nokm, 1, &
                               var2d, rc, tcyclic, fid )
            if ( rc .ne. 0 ) call die(myname,'cannot read '//trim(varn) )
       else
            call GFIO_GetVarT1 ( fid, trim(varn), nymd, nhms, im, jm, 1, km, &
                               var3d, rc, tcyclic, fid )
!           print *, 'rc = ', rc, im, jm, km
            if ( rc .ne. 0 ) call die(myname,'cannot read '//trim(varn) )
       end if

!      Close file
!      ----------
       call GFIO_Close ( fid, rc )
!       print *, myname // ': Closing GFIO file ' // trim(fname)

!   All done
!   --------
    return

end subroutine Read_

end Program ut_PlumeRise
