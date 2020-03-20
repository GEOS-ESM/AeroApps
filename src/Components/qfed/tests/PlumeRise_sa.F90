!
! Read GEOS-4 files with gridded T, H, etc and call the non-ESMF version
! portion of the PlumeRise. This is a stand alone driver for testing 
! this module.
!

   Program PlumeRise_sa

   use m_set_eta
   use m_STrTemplate
   use m_die

   use rconstants, only: cp, rgas, g, kappa=>rocp, p00
   use PlumeRise_Mod

   implicit NONE

   character(len=*), parameter ::  myname = 'PlumeRise'

   type(PlumeRise) :: plume

   integer :: im, jm, km, lm
   integer :: j_60S, i_60W, j_EQ, is, js

   real, allocatable, target, dimension(:,:)   :: ps, phis, lwi
   real, allocatable, target, dimension(:,:,:) :: THETA, Q, H, HE1, RHO, PM
   real, allocatable, target, dimension(:,:,:) :: HE, PE, PEK

   logical                             :: do_biome(N_BIOME)
   real, allocatable, target, dimension(:,:,:) :: fire_areas, z1, z2, p1, p2
   real, allocatable, target, dimension(:)     :: ak, bk

   integer :: ia, iz, ja, jz, i, j   

   integer :: ks, k, n
   real    :: ptop, pint, undef = 1.E+15

   character(len=255) :: metx_fn, dynv_fn, out_fn
   integer :: nymd, nhms

   real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0
   real, parameter :: ha = 10000. ! 1 hectare = 10,000 meter^2

   real :: fsize ! user defined fire size in ha 

   logical :: sfc2top, verbose

!  Metadata for output GFIO file
!  -----------------------------
   integer, parameter :: NVARS = 4 * N_BIOME
   character(len=255)          :: title, source, contact, levunits
   character(len=255), pointer :: vname(:), vtitle(:), vunits(:)
   real,    pointer :: lat(:), lon(:), lev(:), area(:)
   real             :: dlon, dlat
   real,    pointer :: valid_range(:,:), packing_range(:,:)
   integer, pointer :: kmvar(:)
   integer          :: timeinc = 060000
   type field
       real, pointer :: data(:,:)
   end type field
   type(field) var(NVARS)


!                            -----

!  Parse command line
!  ------------------
   call ParseCmd_()

!  Allocate memory
!  ---------------
   call DoMalloc_ ( dynv_fn )

!  Read gridded fields
!  -------------------
   call Read_ ( metx_fn, 'phis', -1, 0, im, jm, 1, var2d=phis )
   call Read_ ( dynv_fn, 'oro', nymd, nhms, im, jm, 1, var2d=lwi )
   call Read_ ( dynv_fn, 'surfp', nymd, nhms, im, jm, 1, var2d=ps )
   call Read_ ( dynv_fn, 'q', nymd, nhms, im, jm, km, var3d=Q )
   call Read_ ( dynv_fn, 'hghte', nymd, nhms, im, jm, km, var3d=HE1 )
   call Read_ ( dynv_fn, 'h', nymd, nhms, im, jm, km, var3d=H )

!  Derived quantities
!  ------------------
   call Derived_()

!  Determine fire size, mask
!  -------------------------
   call FireSize_()

!  Calculate plume rise
!  --------------------
   ia =  1; iz = im; ja =  1; jz =jm
!!!   ia = is; iz = is; ja = js; jz =js
   call PlumeRise_Run ( plume, ia, iz, ja, jz, km, &
                        theta(ia:iz,ja:jz,:),  rho(ia:iz,ja:jz,:), &
                        pm(ia:iz,ja:jz,:), q(ia:iz,ja:jz,:),  &
                        h(ia:iz,ja:jz,:),  he(ia:iz,ja:jz,:), &
                        fire_areas(ia:iz,ja:jz,:),            &
                        p1_plume = p1(ia:iz,ja:jz,:),         &
                        p2_plume = p2(ia:iz,ja:jz,:),         &
                        z1_plume = z1(ia:iz,ja:jz,:),         &
                        z2_plume = z2(ia:iz,ja:jz,:)          )

!  For now, Boreal forest uses same parameters as tropical forests
!  ---------------------------------------------------------------
   z1(:,:,BOREAL_FOREST) = z1(:,:,TROPICAL_FOREST)  
   z2(:,:,BOREAL_FOREST) = z2(:,:,TROPICAL_FOREST)  

!  Mask out ocean/antarctica values
!  --------------------------------
   do n = 1, N_BIOME
      where(lwi /= LAND )
            z1(:,:,n) = undef
            z2(:,:,n) = undef
            p1(:,:,n) = undef
            p2(:,:,n) = undef
      end where
   end do
   z1(:,1:j_60S,:) = undef
   z2(:,1:j_60S,:) = undef
   p1(:,1:j_60S,:) = undef
   p2(:,1:j_60S,:) = undef

!  Write out results
!  -----------------
   call Write_ ( out_fn, nymd, nhms )    ! writes out GFIO file

CONTAINS

!............................................................................

  subroutine Derived_()

   sfc2top = ( hE1(1,1,km) .gt.  hE1(1,1,km-1) ) 
   if ( sfc2top ) then
     HE(:,:,2:km+1) = HE1      
     HE(:,:,1) = phis / g
   else
     HE(:,:,1:km) = HE1
     HE(:,:,km+1) = phis / g
   end if

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
!srf-convert to dry theta
      theta(:,:,k) =  (1./(1+.61*q(:,:,k)))*theta(:,:,k)
			      
   end do

!  Sample output over amazon
!  -------------------------
!           '    222.5     1.5    29.5     0.0    11.7'
   print *, '    H(m)    p (hPa) THETA(K)  rho    q(g/kg)'
   i = is;  j = js
   do k = km, 1, -1
      print 10, h(i,j,k), pm(i,j,k)/100, THETA(i,j,k), &
               rho(i,j,k), 1000*q(i,j,k)
   end do
   print 10, he(i,j,1)
10 format(1x,5f8.2)

   print *, 'phis = ', minval(phis), maxval(phis)
   print *, '  ps = ', minval(ps),   maxval(ps)
   print *, ' oro = ', minval(lwi),  maxval(lwi)
   print *, 'Theta= ', minval(THETA),    maxval(THETA)
   print *, '   q = ', minval(Q),    maxval(Q)
   print *, '   h = ', minval(H),    maxval(H)
   print *, '  he = ', minval(HE),   maxval(HE)
   print *, '  pm = ', minval(PM),   maxval(PM)
   print *, ' rho = ', minval(RHO),  maxval(RHO)
   print *

  end subroutine Derived_

!............................................................................

  subroutine FireSize_()

!  Set fire sizes only over land gridpoints
!  ----------------------------------------
   do n = 1, N_BIOME
      where(lwi == LAND )
            fire_areas(:,:,n) = fsize * ha  ! in meter^2
      elsewhere
            fire_areas(:,:,n) = 0.0
      end where
      if ( .NOT. do_biome(n) ) fire_areas(:,:,n) = 0.0
   end do
   fire_areas(:,1:j_60S,:) = 0.0 ! no fires south of 60S

!  Stats
!  -----
   do n = 1, N_BIOME
   print *, 'area = ', minval(fire_areas(:,:,n)),  maxval(fire_areas(:,:,n)), &
                       nint(100.*count(fire_areas(:,:,n)/=0.0)/(im*jm)), '%'
   end do
   print *

  end subroutine FireSize_

!............................................................................

  subroutine DoMalloc_( fname )
   character(len=*), intent(in) :: fname

   character(len=*), parameter ::  myname = 'DoMalloc_'
   integer :: ios

!  Point roughly over the amazon basin
!  -----------------------------------
   j_60S = 1 + nint( (90.0 -  60.) * (jm-1) / 180. ) 
   j_EQ  = 1 + nint( (90.0 -   0.) * (jm-1) / 180. ) 
   i_60W = 1 + nint( ( 0.0 + 300.) *   im   / 360. ) 

   is = i_60W
   js = j_EQ

!  sanity check
!  ------------
   if ( im < 1 .or. jm<1 .or. km<1 ) then
      print *, 'im, jm, km = ', im, jm, km 
      call die(myname,'invalid dimensions')
   endif

   allocate ( ps(im,jm), &
              phis(im,jm), &
              lwi(im,jm), &
              THETA(im,jm,km), &
              Q(im,jm,km), &
              H(im,jm,km), &
              HE1(im,jm,km), &
              HE(im,jm,km+1), &
              PE(im,jm,km+1), &
              PEK(im,jm,km+1), &
              RHO(im,jm,km), &
              PM(im,jm,km), &
              fire_areas(im,jm,N_BIOME), &
              z1(im,jm,N_BIOME), &
              z2(im,jm,N_BIOME), &
              p1(im,jm,N_BIOME), &
              p2(im,jm,N_BIOME), &
              ak(km+1), bk(km+1), &
                                         stat=ios )

    if ( ios /= 0 ) call die(myname, 'cannot allocate memory' )

  end subroutine DoMalloc_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!-BOP
!
! !IROUTINE:  Read_ --- Reads fields from file and distribute
!
! !INTERFACE:
!
   subroutine Read_ ( filen, varn, nymd1, nhms1, im, jm, km, &
                      var2d, var3d, cyclic )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

  character(len=*), intent(in) :: filen   ! GFIO compatible file name
  character(len=*), intent(in) :: varn    ! variable name
  integer, intent(in)          :: nymd1, nhms1 ! date/time

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
    integer :: fid, rc, ios, nymd, nhms, incSecs
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
                            nymd=nymd1, nhms=nhms1 )
    else
         fname = filen
    end if


!   Read file
!   ---------

       print *, myname // ': Opening GFIO file ' // trim(fname), nymd1, nhms1

!      Open file, get first time
!      -------------------------
       call GFIO_Open ( fname, READ_ONLY, fid, rc )
       if ( rc .ne. 0 ) then
          call die(myname,'cannot open GFIO file '//trim(fname))
       end if

!      File date is negative, get it from file
!      ---------------------------------------
       if ( nymd1 < 0 ) then
          call GetBegDateTime ( fid, nymd, nhms, incSecs, rc )
          if ( rc .ne. 0 ) then
             call die(myname,'cannot determine firs date/time for GFIO file '//trim(fname))
          end if
            print *, 'Using first time/date for file ('//trim(fname)// &
                     ': ', nymd, nhms
       else
          nymd = nymd1
          nhms = nhms1
       endif

!      Read global array
!      -----------------
       print *, myname // ': reading ' // trim(varn) // &
                 ' from ' // trim(fname)
       if ( present(var2d) ) then
            call GFIO_GetVarT1 ( fid, trim(varn), nymd, nhms, im, jm, nokm, 1, &
                               var2d, rc, tcyclic, fid )
            if ( rc .ne. 0 ) call die(myname,'cannot read '//trim(varn) )
            print *, varn,': ', minval(var2d), maxval(var2d)
       else
            call GFIO_GetVarT1 ( fid, trim(varn), nymd, nhms, im, jm, 1, km, &
                               var3d, rc, tcyclic, fid )
!           print *, 'rc = ', rc, im, jm, km
            if ( rc .ne. 0 ) call die(myname,'cannot read '//trim(varn) )
           print *, varn,': ', minval(var3d), maxval(var3d)
       end if

!      Close file
!      ----------
       call GFIO_Close ( fid, rc )
!       print *, myname // ': Closing GFIO file ' // trim(fname)


!   All done
!   --------
    return 

  end subroutine Read_

!....................................................................

subroutine Write_ ( fname, nymd, nhms )    ! writes out GFIO file
  implicit NONE
  character(len=*), intent(in) :: fname        ! diagnostic file name
  integer, intent(in)          :: nymd, nhms   ! date/time
!                              ---
  integer, parameter :: READ_WRITE = 0
  logical :: fexists  
  integer :: i,j,l, fid, rc, rcs(NVARS), prec=0
  
  call meta_init_()
  call grid_()
  inquire ( file=trim(fname), exist=fexists )  
  if ( fexists ) then
      if ( verbose ) print *, ' [] Writing to ', trim(fname)
     call GFIO_Open ( fname, READ_WRITE, fid, rc ) ! open existing file
     if ( rc/=0 ) call die(myname,'cannot open GFIO file '//trim(fname) )  
  else
    if ( verbose ) print *, ' [] Creating ', trim(fname)
    call GFIO_Create ( fname, title, source, contact, undef,       &
                       im, jm, km, lon, lat, lev, levunits,         &
                       nymd, nhms, timeinc,                         &
                       nvars, vname, vtitle, vunits, kmvar,         &
                       valid_range, packing_range, prec,            &
                       fid, rc )
    if ( rc/=0 ) call die(myname,'cannot create GFIO file '//trim(fname) )
  end if
  rcs = 0
!ams  print *,' writeGFIO:NVARS ',NVARS
  do i = 1, NVARS
      if ( verbose ) print *, '       o writing ', trim(vname(i)), &
                                      ' on ', nymd, nhms
      call GFIO_PutVar ( fid, trim(vname(i)), nymd, nhms, im, jm, 0, 1, &
                         var(i)%data, rcs(i) )
  end do
  if ( any(rcs .ne. 0) ) call die(myname,'cannot write fields')
  call GFIO_Close ( fid, rc )
  call meta_clean_()

end subroutine Write_

subroutine meta_init_ ( )       ! allocates local memory
  integer rc, i, j, k

  allocate ( vname(nvars), vunits(nvars), vtitle(nvars), kmvar(nvars), &
             valid_range(2,nvars), packing_range(2,nvars),             &
             lon(im), lat(jm), lev(km), &
             stat=rc )
  if ( rc/=0 ) call die(myname,'cannot allocate mem for metadata')
  levunits = 'hPa'
  title = 'Plume Rise Estimates'
  source = 'Global Modeling and Assimilation Office, NASA/GSFC'
  contact = 'data@gmao.gsfc.nasa.gov'
  kmvar(1:nvars)  = 0  ! all 2D quantites
  valid_range = undef
  packing_range = undef
  i = 0
!

! Tropical Forest
! ---------------
  i = i + 1; if (i>NVARS) call die(myname,'too many variables')
  vname(i)  = 'z1_tf'
  vtitle(i) = 'Bottom Plume Height (Tropical Forest)'
  vunits(i) = 'meter'
  var(i)%data => z1(:,:,1)

  i = i + 1; if (i>NVARS) call die(myname,'too many variables')
  vname(i)  = 'z2_tf'
  vtitle(i) = 'Top Plume Height (Tropical Forest)'
  vunits(i) = 'meter'
  var(i)%data => z2(:,:,1)

  i = i + 1; if (i>NVARS) call die(myname,'too many variables')
  vname(i)  = 'p1_tf'
  vtitle(i) = 'Bottom Plume Pressure (Tropical Forest)'
  vunits(i) = 'hPa'
  var(i)%data => p1(:,:,1)

  i = i + 1; if (i>NVARS) call die(myname,'too many variables')
  vname(i)  = 'p2_tf'
  vtitle(i) = 'Top Plume Pressure (Tropical Forest)'
  vunits(i) = 'hPa'
  var(i)%data => p2(:,:,1)

! Boreal Forest
! -------------
  i = i + 1; if (i>NVARS) call die(myname,'too many variables')
  vname(i)  = 'z1_bf'
  vtitle(i) = 'Bottom Plume Height (Boreal Florest)'
  vunits(i) = 'meter'
  var(i)%data => z1(:,:,2)

  i = i + 1; if (i>NVARS) call die(myname,'too many variables')
  vname(i)  = 'z2_bf'
  vtitle(i) = 'Top Plume Height (Boreal Florest)'
  vunits(i) = 'meter'
  var(i)%data => z2(:,:,2)

  i = i + 1; if (i>NVARS) call die(myname,'too many variables')
  vname(i)  = 'p1_bf'
  vtitle(i) = 'Bottom Plume Pressure (Boreal Florest)'
  vunits(i) = 'hPa'
  var(i)%data => p1(:,:,2)

  i = i + 1; if (i>NVARS) call die(myname,'too many variables')
  vname(i)  = 'p2_bf'
  vtitle(i) = 'Top Plume Pressure (Boreal Florest)'
  vunits(i) = 'hPa'
  var(i)%data => p2(:,:,2)

! Savanna
! -------
  i = i + 1; if (i>NVARS) call die(myname,'too many variables')
  vname(i)  = 'z1_sv'
  vtitle(i) = 'Bottom Plume Height (Savanna)'
  vunits(i) = 'meter'
  var(i)%data => z1(:,:,3)

  i = i + 1; if (i>NVARS) call die(myname,'too many variables')
  vname(i)  = 'z2_sv'
  vtitle(i) = 'Top Plume Height (Savanna)'
  vunits(i) = 'meter'
  var(i)%data => z2(:,:,3)

  i = i + 1; if (i>NVARS) call die(myname,'too many variables')
  vname(i)  = 'p1_sv'
  vtitle(i) = 'Bottom Plume Pressure (Savanna)'
  vunits(i) = 'hPa'
  var(i)%data => p1(:,:,3)

  i = i + 1; if (i>NVARS) call die(myname,'too many variables')
  vname(i)  = 'p2_sv'
  vtitle(i) = 'Top Plume Pressure (Savanna)'
  vunits(i) = 'hPa'
  var(i)%data => p2(:,:,3)

! Grassland
! ---------
  i = i + 1; if (i>NVARS) call die(myname,'too many variables')
  vname(i)  = 'z1_gl'
  vtitle(i) = 'Bottom Plume Height (Grassland)'
  vunits(i) = 'meter'
  var(i)%data => z1(:,:,4)

  i = i + 1; if (i>NVARS) call die(myname,'too many variables')
  vname(i)  = 'z2_gl'
  vtitle(i) = 'Top Plume Height (Grassland)'
  vunits(i) = 'meter'
  var(i)%data => z2(:,:,4)

  i = i + 1; if (i>NVARS) call die(myname,'too many variables')
  vname(i)  = 'p1_gl'
  vtitle(i) = 'Bottom Plume Pressure (Grassland)'
  vunits(i) = 'hPa'
  var(i)%data => p1(:,:,4)

  i = i + 1; if (i>NVARS) call die(myname,'too many variables')
  vname(i)  = 'p2_gl'
  vtitle(i) = 'Top Plume Pressure (Grassland)'
  vunits(i) = 'hPa'
  var(i)%data => p2(:,:,4)

  if ( i /= NVARS ) then
     print *, myname//': actual/declared NVARS: ', i, nvars
     call die(myname, 'less variables than NVARS have been set' )
  end if

end subroutine meta_init_

subroutine meta_clean_()             ! de-allocates local memory
     integer rc
     deallocate ( vname, vunits, vtitle, kmvar, &
             valid_range, packing_range, lon, lat,  &
             stat=rc )
     if ( rc/=0 ) call die(myname,'cannot deallocate mem for metadata')
end subroutine meta_clean_

!..................................................................
   subroutine grid_()

   integer :: i, j
   real :: d2r, pi, factor, phi_1, phi_2

   dlon = 360.0 / im
   dlat = 180.0 / ( jm - 1)

   pi = 4.0 * atan ( 1.0 )
   d2r = pi / 180.0

   do i = 1, im
      lon(i) = 0.0 + (i-1) * dlon
   end do

   do j = 1, jm
      lat(j) = -90. + (j-1)*dlat
   end do

   do k = 1, km
      lev(k) = k
   end do

 end subroutine grid_

!....................................................................

  subroutine ParseCmd_()

     integer        iarg, argc, iargc, i
     character(len=255)  argv, token, out_tmpl

     integer :: fid, lm, nvars, ngatts, rc(4), incSecs, nymd0, nhms0
     integer :: READ_ONLY=1
 
print *
print *,'NAME'
print *,'     '//myname//' - Estimate Smoke Plume Heights'

!    Default values
!    --------------
     metx_fn = '/home/dasilva/out/plume/phis_2x25.hdf'
     dynv_fn = '/home/dasilva/out/plume/diag_2x25.hdf'
     out_tmpl = 'pr_v1.hghts.%y4%m2%d2_%h2%n2z.hdf'
     do_biome = (/ .true., .false., .true.,  .true. /)
     verbose = .false.
!!!     nymd = 20061001
!!!     nhms = 120000
     nymd = -1
     nhms = -1
     fsize = 20.

!     Parse command line 
!     ------------------
      argc =  iargc()
      if ( argc < 2 ) call usage()

      iarg = 0
      do i = 1, argc
         iarg = iarg + 1
         if ( iarg .gt. argc ) call usage()
         call GetArg ( iarg, argv )
         if ( index(argv,'-o') .gt. 0 ) then
           if ( iarg+1 .gt. argc ) call usage()
           iarg = iarg + 1
           call GetArg ( iarg, out_tmpl )
         elseif ( index(argv,'-v') .gt. 0 ) then
           verbose = .true.
         elseif ( index(argv,'-skipTF') .gt. 0 ) then
           do_biome(1) = .false.
         elseif ( index(argv,'-skipBF') .gt. 0 ) then
           do_biome(2) = .false.
         elseif ( index(argv,'-skipSV') .gt. 0 ) then
           do_biome(3) = .false.
         elseif ( index(argv,'-skipGL') .gt. 0 ) then
           do_biome(4) = .false.
         elseif ( index(argv,'-date') .gt. 0 ) then
           if ( iarg+1 .gt. argc ) call usage()
           iarg = iarg + 1
           call GetArg ( iArg, argv ); read(argv,*) nymd
         elseif ( index(argv,'-time') .gt. 0 ) then
           if ( iarg+1 .gt. argc ) call usage()
           iarg = iarg + 1
           call GetArg ( iArg, argv ); read(argv,*) nhms
         elseif ( index(argv,'-fsize') .gt. 0 ) then
           if ( iarg+1 .gt. argc ) call usage()
           iarg = iarg + 1
           call GetArg ( iArg, argv ); read(argv,*) fsize
         else
           exit
         endif
       end do

       if ( (argc-iarg+1) < 2 ) call usage()

       call GetArg ( iarg+0, metx_fn )
       call GetArg ( iarg+1, dynv_fn )

!    Get dimensions/time from diag file (assumes 1 time/file)
!    --------------------------------------------------------
     call GFIO_Open ( dynv_fn, READ_ONLY, fid, rc(1) )
     call GFIO_DimInquire ( fid, im, jm, km, lm, nvars, ngatts, rc(2))
     call GetBegDateTime ( fid, nymd0, nhms0, incSecs, rc(3) )
     if ( nymd < 0 ) nymd = nymd0
     if ( nhms < 0 ) nhms = nhms0
     if ( nymd < 0 .or. nhms < 0 ) then
        if ( lm/=1 ) then
           print *, myname//': file '//trim(dynv_fn)//' has more than 1 time step' 
           print *, myname//': using only first time ', nymd, nhms
        end if
     endif
     call GFIO_Close ( fid, rc(4) )
     if ( any(rc /= 0) ) call die(myname,'cannot query '//trim(dynv_fn))

     call StrTemplate ( out_fn, out_tmpl, nymd = nymd, nhms = nhms )

      print *
      print *, '           Dimensions: ', im, jm
      print *, '            Date/Time: ', nymd, nhms
      print *
      print *, 'Input  PHIS file name: ', trim(metx_fn)
      print *, 'Input  DIAG file name: ', trim(dynv_fn)
      print *, 'Output GFIO  template: ', trim(out_tmpl)
      print *, 'Output GFIO file name: ', trim(out_fn)
      print *, '          Biome flags: ', do_biome
      print *, '            Fire size: ', fsize
      print *

      return

!     Abnormal exit
!     ------------
999   continue
      call usage()

  end subroutine ParseCmd_

!....................................................................

end Program PlumeRise_sa

!..................................................................

   subroutine usage ()

print *
print *,'SYNOPSIS'
print *,'     PlumeRise_sa.x  [OPTIONS]  metx_fn dynv_fn'
print *
print *,'DESCRIPTION'
print *,'     Reads files with the geopotential surface height (metx_fn)'
print *,'     and thermodynamic profiles (dynv_fn) and estimate smoke'
print *,'     plume injection heights using the 1D PyroCu model of'
print *,'     Freitas et al (2006)'
print *
print *, 'OPTIONS'
print *,'     -o  out_tmpl      output HDF template file name; default is '
print *,'                       PlumeRise_v1.hghts.%y4%m2%d2_%h2%n2z.hdf'
print *,'     -date nymd        date to read from file; default is first'
print *,'                         date on file'
print *,'     -time nhms        time to read from file; default is first'
print *,'                         time on file'
print *,'     -fsize  size      fire size in hectars (ha); default is 20.'
print *,'     -skipTF           Skip alculation for TROPICAL FOREST'
print *,'     -skipBF           Skip alculation for BOREAL FOREST'
print *,'     -skipSV           Skip alculation for SAVANNA'
print *,'     -skipGL           Skip alculation for GRASSLAND'
print *,'     -v                enables verbose mode'
print *
print *

  call exit(1)

end subroutine usage

