!
!  Simple f77 wrapper for the Python interface to the Mie Calculator.
!  This version works from columns already interpolated to obs location.



subroutine getMieDims ( mieTableFile, nCh, nRh, nBin, nMom, nPol, rc)

! Retrieve dimensions of MieTable

  use Chem_MieTableMod, only: Chem_MieTableGetDims
  implicit NONE
  character(len=*), intent(in)  :: mieTableFile ! input Mie Table file name            
  integer,          intent(out) :: nCh
  integer,          intent(out) :: nRh
  integer,          intent(out) :: nBin
  integer,          intent(out) :: nMom
  integer,          intent(out) :: nPol
  integer,          intent(out) :: rc ! return code                             

  call Chem_MieTableGetDims(mieTableFile, nCh, nRh, nBin, nMom, nPol, rc)

end subroutine getMieDims

!.............................................................................

 subroutine getEdgeVars ( km, nobs, airdens, delp, ptop, &
                          pe, ze, te )
!
! Get pressure and temperature at the edge of layers, 
!

  implicit NONE

  integer,          intent(in)  :: km               ! number vertical layers
  integer,          intent(in)  :: nobs             ! number of profiles

  real,             intent(in)  :: airdens(km,nobs)
  real,             intent(in)  :: delp(km,nobs)
  real,             intent(in)  :: ptop             ! top (edge) pressure [Pa]

  real,             intent(out) :: pe(km+1,nobs)    ! edge pressure [Pa]
  real,             intent(out) :: ze(km+1,nobs)    ! edge height above sfc [m]
  real,             intent(out) :: te(km+1,nobs)    ! edge Temperature [K]

!                               ---

  integer             :: k, n
 
  real                :: alpha, pm(km), tm(km)
  
  real, parameter     :: grav = 9.80616
  real, parameter     :: Rgas = 287.
  real, parameter     :: kappa = 2.0 / 7.0

  do n = 1, nobs

        pe(1,n) = ptop
        do k = 1, km
           pe(k+1,n) = pe(k,n) + delp(k,n)
        end do

        ze(km+1,n) = 0.0 ! height above surface
        do k = km, 1, -1
           ze(k,n) = ze(k+1,n) + delp(k,n) / ( airdens(k,n) * grav ) 
        end do

        do k = 1, km
           pm(k) = ( pe(k,n) + pe(k+1,n) ) / 2.0
           tm(k) = pm(k) / ( airdens(k,n) * Rgas )
        end do

        te(1,n) = tm(1) ! isothermal at highest level
        do k = 2, km
           alpha = log(pe(k,n)/pm(k-1))/log(pm(k)/pm(k-1))
           te(k,n) = tm(k-1) + alpha * (tm(k) - tm(k-1) )
        end do
                              ! dry adiabatic
        te(km+1,n) = tm(km) * (pe(km+1,n)/pm(km))**kappa

     end do

end subroutine getEdgeVars

!...................................................................................

 subroutine getAOPscalar ( km, nobs, nch, nq, rcfile, channels, vname, verbose, &
                           qm, rh, &
                           tau, ssa, g, rc )

! Returns aod, ssa and asymmetry factor profiles.

  use Chem_MieMod
  implicit NONE

  integer,          intent(in)  :: km               ! number vertical layers
  integer,          intent(in)  :: nobs             ! number of profiles

  character(len=*), intent(in)  :: rcfile           ! resource file, e.g., Aod_EOS.rc

  integer,          intent(in)  :: nch              ! number of channels
  real,             intent(in)  :: channels(nch)

  integer,          intent(in)  :: nq               ! number of tracers
  character,        intent(in)  :: vname(nq,16)     ! variable name

  integer,          intent(in)  :: verbose

  real,             intent(in)  :: qm(km,nq,nobs)   ! (mixing ratio) * delp/g
  real,             intent(in)  :: rh(km,nobs)      ! relative humidity

  real,             intent(out) :: tau(km,nch,nobs) ! aerosol optical depth
  real,             intent(out) :: ssa(km,nch,nobs) ! single scattering albedo
  real,             intent(out) :: g(km,nch,nobs)   ! asymmetry factor
 
  integer,          intent(out) :: rc

!                               ---

  type(Chem_Mie)      :: mieTables
  real                :: idxChannel(nch) ! this should have been integer
  integer             :: idxTable
  character(len=16)   :: vname_(nq)

  integer :: iq, n, m, i, k
  real :: tau_, ssa_, g_

  rc = 0

! Deal with f2py strange handling of strings
! ------------------------------------------
  do iq = 1, nq
     do n = 1, 16
        vname_(iq)(n:n) = vname(iq,n)
     end do
  end do

! Create the Mie Tables
! ---------------------
  mieTables = Chem_MieCreate(rcfile,rc)
  if ( rc /= 0 ) then
     print *, 'Cannot create Mie tables from '//trim(rcfile)
     return
  end if

! Determine channel indices
! -------------------------
  do n = 1, nch
     idxChannel(n) = -1 ! this is really the channel index
     do m = 1, mieTables%nch
        if ( abs(channels(n) - (1.e9)*mieTables%channels(m)) < 1. ) then
           idxChannel(n) = m
           exit
         end if
      end do
   end do
   if ( any(idxChannel<0) ) then
        print *, 'Mie resource files does not set the required channel'
        print *, 'Channels requested:  ', channels
        print *, 'Channels on RC file: ', 1.e+9 * mieTables%channels
        rc = 99
        return
   end if

! Initialize output arrays to zero
! --------------------------------
  tau = 0.0
  ssa = 0.0
  g = 0.0

! Loop over aerosol species
! -------------------------
  do iq = 1,nq
     idxTable = Chem_MieQueryIdx(mieTables,vname_(iq),rc)
     if(idxTable == -1) cycle
     if ( rc/=0 ) then
        print *, 'cannot get Mie index for '//vname_(iq)
        return
     end if

     if (verbose==1) &
          print *, '[+] Adding '//trim(vname_(iq))//' contribution'

!    Loop over nobs, km, nch
!    --------------------------
     do i = 1, nobs
        do n = 1, nch
           do k =1, km
        
            call Chem_MieQuery(mieTables, idxTable, idxChannel(n), &
                               qm(k,iq,i), rh(k,i), tau=tau_, ssa=ssa_, gasym=g_)
 
                            tau(k,n,i) = tau(k,n,i) + tau_
                            ssa(k,n,i) = ssa(k,n,i) + ssa_ * tau_ 
                              g(k,n,i) =   g(k,n,i) +   g_ * ssa_ * tau_

            end do  ! end nch
         end do ! end km
      end do  ! end nobs
  end do ! end tracers

!  Normalize ssa and g
!  -------------------
   where (     tau > 0.0 ) ssa = ssa / tau
   where ( ssa*tau > 0.0 )   g =  g / ( ssa * tau ) 

   call Chem_MieDestroy(mieTables, rc)
   if ( rc /= 0 ) then
     print *, 'Cannot destroy MieTables'
     return
  end if


   end subroutine getAOPscalar

!...................................................................................

 subroutine getAOPvector ( km, nobs, nch, nq, rcfile, channels, vname, verbose, &
                           qm, rh,nMom,nPol, &
                           tau, ssa, g, pmom, rc )

! Returns aod, ssa, asymmetry factor, phase function profiles.

  use Chem_MieMod
  implicit NONE

  integer,          intent(in)  :: km               ! number vertical layers
  integer,          intent(in)  :: nobs             ! number of profiles

  character(len=*), intent(in)  :: rcfile           ! resource file, e.g., Aod_EOS.rc

  integer,          intent(in)  :: nch              ! number of channels
  real,             intent(in)  :: channels(nch)

  integer,          intent(in)  :: nMom             ! number of legender momemts 
  integer,          intent(in)  :: nPol             ! number of components of the scattering matrix

  integer,          intent(in)  :: nq               ! number of tracers
  character,        intent(in)  :: vname(nq,16)     ! variable name

  integer,          intent(in)  :: verbose

  real,             intent(in)  :: qm(km,nq,nobs)   ! (mixing ratio) * delp/g
  real,             intent(in)  :: rh(km,nobs)      ! relative humidity

  real,             intent(out) :: tau(km,nch,nobs) ! aerosol optical depth
  real,             intent(out) :: ssa(km,nch,nobs) ! single scattering albedo
  real,             intent(out) :: g(km,nch,nobs)   ! asymmetry factor
  real,             intent(out) :: pmom(km,nMom,nPol,nch,nobs) ! elements of scattering phase matrix
 
  integer,          intent(out) :: rc

! -------------------------------------------------------------------------------------
! PMOM Notes:
!
! The paramer "pmom" contains elements of the phase matrix and it is dimensioned as
! (npol,nmom,radius,rh,lambda) where
!
!    nmom = number of phase function moments in the file
!    npol = index of the phase function quantity
!    npol = 4 or 6
!
! index  moments of quantity
! -----  -------------------
! 1      P11
! 2      P12
! 3      P33
! 4      P34
! 5      P22
! 6      P44
!
! Note that for Mie based tables we only have npol = 4.  For the ellipsoid
! tables we have npol = 6 and have appended the P22 and P44 as the last
! two indices.  From symmetry, for Mie based tables we know P22 = P11 and
! P44 = P33.
!
! In this code, if nPol=6 is requested, we allways return 6 elements of the
! phase matrix, setting P22=P11 and P44=P33 if needed.
!
! -------------------------------------------------------------------------------------

!                               ---

  type(Chem_Mie)      :: mieTables
  real                :: idxChannel(nch) ! this should have been integer
  integer             :: idxTable
  character(len=16)   :: vname_(nq)

  integer             :: iq, n, m, i, k,iMom, iPol, nPol_
  real                :: tau_, ssa_, g_
  real, pointer       :: pmom_(:,:)
  logical             :: spherical_ext ! whether to set P22=P11 and P44=P11
  integer             :: iP11=1, iP12=2, iP33=3, iP34=4, iP22=5, iP44=6

  rc = 0
  
! Deal with f2py strange handling of strings
! ------------------------------------------
  do iq = 1, nq
     do n = 1, 16
        vname_(iq)(n:n) = vname(iq,n)
     end do
  end do

! Create the Mie Tables
! ---------------------
  mieTables = Chem_MieCreate(rcfile,rc)
  print *, 'mietables', mieTables%bc_optics_file
  if ( rc /= 0 ) then
     print *, 'Cannot create Mie tables from '//trim(rcfile)
     return
  end if

  if ( nMom > mieTables%nMom ) then ! mieTables%nMom is writen in Aod_EOS.rc file
     print *, 'mieTables do not have enough moments', nMom, mieTables%nMom
     rc = 99
     return
  end if

! Determine channel indices
! -------------------------
  do n = 1, nch
     idxChannel(n) = -1 ! this is really the channel index
     do m = 1, mieTables%nch
        if ( abs(channels(n) - (1.e9)*mieTables%channels(m)) < 1. ) then
           idxChannel(n) = m
           exit
         end if
      end do
   end do
   if ( any(idxChannel<0) ) then
        print *, 'Mie resource files does not set the required channel'
        print *, 'Channels requested:  ', channels
        print *, 'Channels on RC file: ', 1.e+9 * mieTables%channels
        rc = 99
        return
   end if

! Allocate memory for phase function; size determined by what kind
! table has been loaded; n_moments given in Aod_EOS.rc
! ----------------------------------------------------------------
!  print *,'nmom mietable', mieTables%nMom,mieTables%vtableUse%nMom
!  print *,'npol mietable', mieTables%nPol,mieTables%vtableUse%nPol
 

! Initialize output arrays to zero
! --------------------------------
  tau = 0.0
  ssa = 0.0
  g = 0.0
  pmom=0.0

! Loop over aerosol species
! -------------------------
  do iq = 1,nq

     idxTable = Chem_MieQueryIdx(mieTables,vname_(iq),rc)
     if(idxTable == -1) cycle
     if ( rc/=0 ) then
        print *, 'cannot get Mie index for '//vname_(iq)
        return
     end if

     if (verbose==1) &
          print *, '[+] Adding '//trim(vname_(iq))//' contribution'

!    Check number of moments on file for this species
!    ------------------------------------------------
     if ( nMom > mieTables%vtableUse%nMom ) then 
        print *, 'ERROR: mieTables do not have enough moments', nMom, mieTables%vtableUse%nMom 
        rc = 666
        return
     end if

!    Handle possible case of non-spherical dust, but spherical everything else
!    -------------------------------------------------------------------------
     if ( nPol .LE. mieTables%vtableUse%nPol ) then  ! user requests fewer polarizations
          nPol_ = nPol
          print*,'nPol mie table', mieTables%vtableUse%nPol
          spherical_ext = .FALSE.
     else if ( nPol == 6 .AND. mieTables%vtableUse%nPol == 4 ) then
          nPol_ = 4              ! file only has 4, user wants 6
          spherical_ext = .TRUE. ! special case, will set P22=P11, P44=P33
     else
          rc = 777
          print *, 'ERROR: inconsistent number of polarizations: ', nPol, mieTables%vtableUse%nPol 
          return
     end if
     print *, 'npol_ test',iq, nPol_
 
     allocate(pmom_(mieTables%nMom,nPol_),stat=rc)
     if ( rc /= 0 ) then
     print *, 'Cannot allocate memory for pmom_'
     return
     end if
!     STOP
!    Loop over nobs, km, nch
!    --------------------------
     do i = 1, nobs
        do n = 1, nch
           do k =1, km
            
            call Chem_MieQuery(mieTables, idxTable, idxChannel(n), &
                               qm(k,iq,i), rh(k,i), tau=tau_, ssa=ssa_, gasym=g_,pmom=pmom_)
 
            tau(k,n,i) = tau(k,n,i) + tau_
            ssa(k,n,i) = ssa(k,n,i) + ssa_ * tau_ 
              g(k,n,i) =   g(k,n,i) +   g_ * ssa_ * tau_
            do ipol=1, nPol_      
               do imom=1, nMom   
                  pmom(k,imom,ipol,n,i) = pmom(k,imom,ipol,n,i) &
                                        + pmom_(imom,ipol) * ssa_ * tau_  
                  
               end do
            end do

!           Special handling, spherical symmetry
!           ------------------------------------
            if ( spherical_ext ) then
               pmom(k,:,iP22,n,i) = pmom(k,:,iP11,n,i)  
               pmom(k,:,iP44,n,i) = pmom(k,:,iP33,n,i)  
            end if

            end do  ! end nch
         end do ! end km
      end do  ! end nobs
  deallocate(pmom_)
  end do ! end tracers
  print*, 'end tracer', iq
!  Normalize ssa and g
!  -------------------
   where (     tau > 0.0 ) ssa = ssa / tau
   where ( ssa*tau > 0.0 )   g =  g / ( ssa * tau ) 
   do i = 1,nobs
      do n = 1, nch
         do k = 1, km
            do ipol=1, nPol   ! normalize  Pmom 
               do imom=1, nMom 
                  if (( ssa(k,n,i)  * tau(k,n,i)  ) > 0.0 ) then
                  pmom(k,imom,ipol,n,i)  = pmom(k,imom,ipol,n,i)  / ( ssa(k,n,i)  * tau(k,n,i)  ) 
                  end if
               end do
            end do   
         end do
      end do
   end do
   
   

   call Chem_MieDestroy(mieTables, rc)
   if ( rc /= 0 ) then
     print *, 'Cannot destroy MieTables'
     return
   end if

   end subroutine getAOPvector

!...................................................................................

 subroutine getO3tauOMI ( km, nobs, nch, channels, verbose, o3, t, delp, &
                          tau, rc )

! Return the extinction profile due to gaseous ozone absorption, based on the 
! OMI channels and the extinction cross sections provided by Can Li convolved
! with OMI slit function

  implicit NONE

  integer,          intent(in)  :: km               ! number vertical layers
  integer,          intent(in)  :: nobs             ! number of profiles
  integer,          intent(in)  :: nch              ! number of channels
  real,             intent(in)  :: channels(nch)

  integer,          intent(in)  :: verbose

  real,             intent(in)  :: o3(km,nobs)     ! ozone volume mixing ratio
  real,             intent(in)  :: t(km,nobs)      ! temperature [K]
  real,             intent(in)  :: delp(km,nobs)   ! pressure thickness of level [Pa]

  real,             intent(out) :: tau(km,nch,nobs) ! aerosol optical depth
 
  integer,          intent(out) :: rc

! -------------------------------------------------------------------------------------
! Cross section notes
!
! Can Li provides O3 cross section coefficients at OMI wavelengths convolved with
! OMI slit function.  Given ozone (converted to # cm-2 in vertical level) calculate
! the absorbtion cross section (cm2) derived from the data coefficients:
!  xsec = 1.e-20 * a0 * (1. + a1*T + a2*T*T)
! where T is the temperature [deg C] and xsec is [cm2]

  real, parameter :: omichan(16) = &
          (/   308.700, 312.340, 312.610, 317.350, 317.620, &
               322.420, 325.000, 331.060, 331.340, 339.660, &
               354.000, 359.880, 360.150, 379.950, 388.000, 471.000 /)
  real, parameter :: a0(16) = &
          (/   11.9240, 6.97293, 6.63798, 3.66506, 3.72993, &
               2.08910, 1.56904, 0.689132, 0.692929, 0.121311, &
               0.00884554, 0.00656508, 0.00628016, 0.000462138, 0.000513032, 0. /)
  real, parameter :: a1(16) = &
          (/   0.00264655, 0.00329817, 0.00343651, 0.00371339, 0.00330990, &
               0.00344314, 0.00283151, 0.00402877, 0.00385602, 0.0127579, &
               0.0197702,  0.00803209, 0.0104759,  0.0186137,  0.00696032, 0. /)
  real, parameter :: a2(16) = &
          (/   2.11937e-05, 2.01685e-05, 2.10381e-05, 2.74817e-05, 2.66713e-05, &
               2.58214e-05, 2.34310e-05, 3.11755e-05, 2.57440e-05, 7.52920e-05, &
               0.000136205, 6.48508e-05, 6.97569e-05, 8.70968e-05, -3.86025e-06, 0. /)

  integer, parameter  :: nch_omi = 16

  real, parameter     :: grav = 9.80616   ! acceleration of gravity [m s-2]
  real, parameter     :: Na   = 6.022e23  ! Avogadro's numbers
  real, parameter     :: mair = 0.02897   ! air mole weight [kg mole-1]


  integer             :: idxChannel(nch)
  integer             :: iq, n, m, i, k
  real                :: xsec, t_, fac

  rc = 0
  
! Determine channel indices
! -------------------------
  do n = 1, nch
     idxChannel(n) = -1 ! this is really the channel index
     do m = 1, nch_omi
        if ( abs(channels(n) - omichan(m)) < 1. ) then
           idxChannel(n) = m
           exit
         end if
      end do
   end do
   if ( any(idxChannel<0) ) then
        print *, 'Requested channel not available'
        print *, 'Channels requested:  ', channels
        print *, 'Channels on RC file: ', omichan
        rc = 99
        return
   end if


! Initialize output arrays to zero
! --------------------------------
  tau = 0.0

! Loop over nobs, km, nch
! --------------------------
  fac = 1.e-4 * na/mair/grav
  do i = 1, nobs
     do n = 1, nch
        m = idxchannel(n)
        do k =1, km
            t_ = t(k,i)-273.15  ! degrees C
            xsec = 1.e-20 * a0(m) * (1. + a1(m)*t_ + a2(m)*t_**2)
            if(xsec < 1.e-30) xsec = 1.e-30
            tau(k,n,i) = xsec*o3(k,i)*delp(k,i)*fac
        end do  ! end nch
     end do ! end km
  end do  ! end nobs

  end subroutine getO3tauOMI

!...................................................................................

 subroutine getO3tauOMPS ( km, nobs, nch, channels, verbose, o3, t, delp, &
                           tau, rc )

! Return the extinction profile due to gaseous ozone absorption, based on the 
! OMPS channels and the extinction cross sections provided by Can Li convolved
! with OMPS slit function

  implicit NONE

  integer,          intent(in)  :: km               ! number vertical layers
  integer,          intent(in)  :: nobs             ! number of profiles
  integer,          intent(in)  :: nch              ! number of channels
  real,             intent(in)  :: channels(nch)

  integer,          intent(in)  :: verbose

  real,             intent(in)  :: o3(km,nobs)     ! ozone volume mixing ratio
  real,             intent(in)  :: t(km,nobs)      ! temperature [K]
  real,             intent(in)  :: delp(km,nobs)   ! pressure thickness of level [Pa]

  real,             intent(out) :: tau(km,nch,nobs) ! aerosol optical depth
 
  integer,          intent(out) :: rc

! -------------------------------------------------------------------------------------
! Cross section notes
!
! Can Li provides O3 cross section coefficients at OMPS wavelengths convolved with
! OMPS slit function.  Given ozone (converted to # cm-2 in vertical level) calculate
! the absorbtion cross section (cm2) derived from the data coefficients:
!  xsec = 1.e-20 * a0 * (1. + a1*T + a2*T*T)
! where T is the temperature [deg C] and xsec is [cm2]

  real, parameter :: omichan(16) = &
          (/   308.700, 312.340, 312.610, 317.350, 317.620, &
               322.420, 325.000, 331.060, 331.340, 339.660, &
               354.000, 359.880, 360.150, 379.950, 388.000, 471.000 /)
  real, parameter :: a0(16) = &
          (/   11.9184, 7.02340, 6.70705, 3.66829, 3.62842, &
               2.04284, 1.39918, 0.645718,0.617893,0.126645, &
               0.00935895, 0.00644011,  0.00620974, 0.000450954, 0.000516689, 0./)
  real, parameter :: a1(16) = &
          (/      0.00260584, 0.00327230, 0.00336150, 0.00373411, &
                  0.00357413, 0.00355121, 0.00377594, 0.00440684, &
                  0.00481047, 0.0115148,  0.0193700,  0.00988421, &
                  0.0108591,  0.0185442,  0.00692930, 0./)
  real, parameter :: a2(16) = &
          (/     2.04564e-05, 2.02593e-05, 2.09936e-05, 2.72707e-05, &
                 2.65668e-05, 2.62994e-05, 2.47051e-05, 3.04936e-05, &
                 2.96094e-05, 7.17826e-05, 0.000135241, 7.71440e-05, &
                 7.59650e-05, 8.85767e-05, -5.70736e-06, 0./)

  integer, parameter  :: nch_omi = 16

  real, parameter     :: grav = 9.80616   ! acceleration of gravity [m s-2]
  real, parameter     :: Na   = 6.022e23  ! Avogadro's numbers
  real, parameter     :: mair = 0.02897   ! air mole weight [kg mole-1]


  integer             :: idxChannel(nch)
  integer             :: iq, n, m, i, k
  real                :: xsec, t_, fac

  rc = 0
  
! Determine channel indices
! -------------------------
  do n = 1, nch
     idxChannel(n) = -1 ! this is really the channel index
     do m = 1, nch_omi
        if ( abs(channels(n) - omichan(m)) < 1. ) then
           idxChannel(n) = m
           exit
         end if
      end do
   end do
   if ( any(idxChannel<0) ) then
        print *, 'Requested channel not available'
        print *, 'Channels requested:  ', channels
        print *, 'Channels on RC file: ', omichan
        rc = 99
        return
   end if


! Initialize output arrays to zero
! --------------------------------
  tau = 0.0

! Loop over nobs, km, nch
! --------------------------
  fac = 1.e-4 * na/mair/grav
  do i = 1, nobs
     do n = 1, nch
        m = idxchannel(n)
        do k =1, km
            t_ = t(k,i)-273.15  ! degrees C
            xsec = 1.e-20 * a0(m) * (1. + a1(m)*t_ + a2(m)*t_**2)
            if(xsec < 1.e-30) xsec = 1.e-30
            tau(k,n,i) = xsec*o3(k,i)*delp(k,i)*fac
        end do  ! end nch
     end do ! end km
  end do  ! end nobs

  end subroutine getO3tauOMPS

!...................................................................................

 subroutine getSO2tauOMI ( km, nobs, nch, channels, verbose, so2, t, delp, &
                           tau, rc )

! Return the extinction profile due to gaseous ozone absorption, based on the 
! OMI channels and the extinction cross sections provided by Can Li convolved
! with OMI slit function

  implicit NONE

  integer,          intent(in)  :: km               ! number vertical layers
  integer,          intent(in)  :: nobs             ! number of profiles
  integer,          intent(in)  :: nch              ! number of channels
  real,             intent(in)  :: channels(nch)

  integer,          intent(in)  :: verbose

  real,             intent(in)  :: so2(km,nobs)    ! SO2 mass mixing ratio
  real,             intent(in)  :: t(km,nobs)      ! temperature [K]
  real,             intent(in)  :: delp(km,nobs)   ! pressure thickness of level [Pa]

  real,             intent(out) :: tau(km,nch,nobs) ! aerosol optical depth
 
  integer,          intent(out) :: rc

! -------------------------------------------------------------------------------------
! Cross section notes
!
! Can Li provides SO2 cross section coefficients at OMI wavelengths convolved with
! OMI slit function.  Given SO2 (converted to # cm-2 in vertical level) calculate
! the absorbtion cross section (cm2) derived from the data coefficients:
!  xsec = 1.e-20 * a0 * (1. + a1*T + a2*T*T)
! where T is the temperature [deg C] and xsec is [cm2]

  real, parameter :: omichan(16) = &
          (/   308.700, 312.340, 312.610, 317.350, 317.620, &
               322.420, 325.000, 331.060, 331.340, 339.660, &
               354.000, 359.880, 360.150, 379.950, 388.000, 471.000 /)
  real, parameter :: a0(16) = &
          (/   45.7172, 13.1613, 14.4654, 8.44695, 7.53393, &
               2.41046, 0.939168, 0.140702, 0.111248, 0.0117402, &
               0.00909536, 0.0143170, 0.00701996, 0.0144760, 0.0314615, 0. /)
  real, parameter :: a1(16) = &
          (/    -0.000903696, 0.000393644, 0.000263059, 0.000331801, 0.00117448, &
                 0.00199029,  0.00521772,  0.0106194,   0.0141481,   0.0347246, &
                 -0.0183250, -0.00523190, -0.00890666, -0.00396620, -0.00130793, 0. /)
  real, parameter :: a2(16) = &
          (/     3.92349e-07, -8.07379e-07, -3.42160e-06, 1.25723e-06, 4.57435e-06, &
                 1.94237e-05,  4.83911e-05,  0.000152681, 0.000205431, 0.00111805, &
                 0.000480777,  0.000464484,  0.00118211,  2.23354e-05, 3.31965e-05, 0. /)

  integer, parameter  :: nch_omi = 16

  real, parameter     :: grav = 9.80616   ! acceleration of gravity [m s-2]
  real, parameter     :: Na   = 6.022e23  ! Avogadro's numbers
  real, parameter     :: mair = 0.02897   ! air mole weight [kg mole-1]
  real, parameter     :: mso2 = 0.064     ! SO2 mole weight [kg mole-1]


  integer             :: idxChannel(nch)
  integer             :: iq, n, m, i, k
  real                :: xsec, t_, fac

  rc = 0
  
! Determine channel indices
! -------------------------
  do n = 1, nch
     idxChannel(n) = -1 ! this is really the channel index
     do m = 1, nch_omi
        if ( abs(channels(n) - omichan(m)) < 1. ) then
           idxChannel(n) = m
           exit
         end if
      end do
   end do
   if ( any(idxChannel<0) ) then
        print *, 'Requested channel not available'
        print *, 'Channels requested:  ', channels
        print *, 'Channels on RC file: ', omichan
        rc = 99
        return
   end if


! Initialize output arrays to zero
! --------------------------------
  tau = 0.0

! Loop over nobs, km, nch
! --------------------------
  fac = 1.e-4 * na/mair/grav*(mair/mso2)
  do i = 1, nobs
     do n = 1, nch
        m = idxchannel(n)
        do k =1, km
            t_ = t(k,i)-273.15  ! degrees C
            xsec = 1.e-20 * a0(m) * (1. + a1(m)*t_ + a2(m)*t_**2)
            if(xsec < 1.e-30) xsec = 1.e-30
            tau(k,n,i) = xsec*so2(k,i)*delp(k,i)*fac
        end do  ! end nch
     end do ! end km
  end do  ! end nobs

  end subroutine getSO2tauOMI

!...................................................................................

 subroutine getSO2tauOMPS ( km, nobs, nch, channels, verbose, so2, t, delp, &
                            tau, rc )

! Return the extinction profile due to gaseous ozone absorption, based on the 
! OMPS channels and the extinction cross sections provided by Can Li convolved
! with OMPS slit function

  implicit NONE

  integer,          intent(in)  :: km               ! number vertical layers
  integer,          intent(in)  :: nobs             ! number of profiles
  integer,          intent(in)  :: nch              ! number of channels
  real,             intent(in)  :: channels(nch)

  integer,          intent(in)  :: verbose

  real,             intent(in)  :: so2(km,nobs)    ! SO2 mass mixing ratio
  real,             intent(in)  :: t(km,nobs)      ! temperature [K]
  real,             intent(in)  :: delp(km,nobs)   ! pressure thickness of level [Pa]

  real,             intent(out) :: tau(km,nch,nobs) ! aerosol optical depth
 
  integer,          intent(out) :: rc

! -------------------------------------------------------------------------------------
! Cross section notes
!
! Can Li provides SO2 cross section coefficients at OMPS wavelengths convolved with
! OMPS slit function.  Given SO2 (converted to # cm-2 in vertical level) calculate
! the absorbtion cross section (cm2) derived from the data coefficients:
!  xsec = 1.e-20 * a0 * (1. + a1*T + a2*T*T)
! where T is the temperature [deg C] and xsec is [cm2]

  real, parameter :: omichan(16) = &
          (/   308.700, 312.340, 312.610, 317.350, 317.620, &
               322.420, 325.000, 331.060, 331.340, 339.660, &
               354.000, 359.880, 360.150, 379.950, 388.000, 471.000 /)
  real, parameter :: a0(16) = &
          (/   35.7355, 13.9534, 15.6521, 8.15167, 7.99223, &
               2.26554, 0.961155, 0.142323, 0.117623, 0.00994909, &
               0.0110929, 0.0144042, 0.00926030, 0.0132294, 0.0303940, 0. /)
  real, parameter :: a1(16) = &
          (/    -3.17972e-05, 0.000628254, 0.000250902, 0.000519986, &
                 0.000608687, 0.00247980,  0.00455231,  0.0119573, &
                 0.0147213,   0.0444674,  -0.0214338,  -0.00672439, &
                -0.00672764, -0.000797140,-0.000876723, 0. /)
  real, parameter :: a2(16) = &
          (/      3.07185e-06, 1.51349e-06,-1.04839e-06, 2.35145e-06, &
                  2.87820e-06, 2.12412e-05, 4.32081e-05, 0.000174339, &
                  0.000212400, 0.00124433,  0.000577062, 0.000705503, &
                  0.00123442,  0.000126878, 4.59309e-05, 0. /)

  integer, parameter  :: nch_omi = 16

  real, parameter     :: grav = 9.80616   ! acceleration of gravity [m s-2]
  real, parameter     :: Na   = 6.022e23  ! Avogadro's numbers
  real, parameter     :: mair = 0.02897   ! air mole weight [kg mole-1]
  real, parameter     :: mso2 = 0.064     ! SO2 mole weight [kg mole-1]


  integer             :: idxChannel(nch)
  integer             :: iq, n, m, i, k
  real                :: xsec, t_, fac

  rc = 0
  
! Determine channel indices
! -------------------------
  do n = 1, nch
     idxChannel(n) = -1 ! this is really the channel index
     do m = 1, nch_omi
        if ( abs(channels(n) - omichan(m)) < 1. ) then
           idxChannel(n) = m
           exit
         end if
      end do
   end do
   if ( any(idxChannel<0) ) then
        print *, 'Requested channel not available'
        print *, 'Channels requested:  ', channels
        print *, 'Channels on RC file: ', omichan
        rc = 99
        return
   end if


! Initialize output arrays to zero
! --------------------------------
  tau = 0.0

! Loop over nobs, km, nch
! --------------------------
  fac = 1.e-4 * na/mair/grav*(mair/mso2)
  do i = 1, nobs
     do n = 1, nch
        m = idxchannel(n)
        do k =1, km
            t_ = t(k,i)-273.15  ! degrees C
            xsec = 1.e-20 * a0(m) * (1. + a1(m)*t_ + a2(m)*t_**2)
            if(xsec < 1.e-30) xsec = 1.e-30
            tau(k,n,i) = xsec*so2(k,i)*delp(k,i)*fac
        end do  ! end nch
     end do ! end km
  end do  ! end nobs

  end subroutine getSO2tauOMPS

!...................................................................................

  subroutine getExt ( km, nobs, nch, nq, channels, vname, verbose, &
                     qc,qm, rh, ext, sca, bsc, absc_SFC,absc_TOA, rc )

! Returns aerosol extinction profile.

  use Chem_MieMod
  implicit NONE

  integer,          intent(in)  :: km               ! number vertical layers
  integer,          intent(in)  :: nobs             ! number of profiles

  integer,          intent(in)  :: nch              ! number of channels
  real,             intent(in)  :: channels(nch)

  integer,          intent(in)  :: nq               ! number of tracers
  character,        intent(in)  :: vname(nq,16)     ! variable name

  integer,          intent(in)  :: verbose

  real,             intent(in)  :: qc(km,nq,nobs)   ! (mixing ratio) * (air density)
  real,             intent(in)  :: qm(km,nq,nobs)   ! (mixing ratio) * delp/g
  real,             intent(in)  :: rh(km,nobs)      ! relative humidity

  real,             intent(out) :: ext(km,nch,nobs) ! total aerosol extinction
  real,             intent(out) :: sca(km,nch,nobs) ! scattering extinction
  real,             intent(out) :: bsc(km,nch,nobs) ! total aero backscatter (toa)
  real,             intent(out) :: absc_TOA(km,nch,nobs) ! attenuated aero bascatter (from toa)
  real,             intent(out) :: absc_SFC(km,nch,nobs) ! attenuated aero bascatter (from surface)

  integer,          intent(out) :: rc

!                               ---

  type(Chem_Mie)      :: mieTables
  real                :: idxChannel(nch) ! this should have been integer
  integer             :: idxTable
  character(len=16)   :: vname_(nq)

  integer :: iq, n, m, i, k, l
  real :: ext_, bsc_,ssa_,bext_,tau_,taulev
  real :: tau(km,nch,nobs)

  rc = 0

! Deal with f2py strange handling of strings
! ------------------------------------------
  do iq = 1, nq
     do n = 1, 16
        vname_(iq)(n:n) = vname(iq,n)
     end do
  end do

! Create the Mie Tables
! ---------------------
  mieTables = Chem_MieCreate('Aod_EOS.rc',rc)
  if ( rc /= 0 ) then
     print *, 'Cannot create Mie tables from Aod_EOS.rc'
     return
  end if

! Determine channel indices
! -------------------------
  do n = 1, nch
     idxChannel(n) = -1 ! this is really the channel index
     do m = 1, mieTables%nch
        if ( abs(channels(n) - (1.e9)*mieTables%channels(m)) < 1. ) then
           idxChannel(n) = m
           exit
         end if
      end do
   end do
   if ( any(idxChannel<0) ) then
        print *, 'Mie resource files does not set the required channel'
        print *, 'Channels requested:  ', channels
        print *, 'Channels on RC file: ', 1.e+9 * mieTables%channels
        rc = 99
        return
   end if

! Initialize output arrays to zero
! --------------------------------
  ext = 0.0
  sca = 0.0
  bsc = 0.0
  tau = 0.0
  absc_TOA = 0.0
  absc_SFC = 0.0

! Loop over aerosol species
! -------------------------
  do iq = 1,nq
     idxTable = Chem_MieQueryIdx(mieTables,vname_(iq),rc)
     if(idxTable == -1) cycle
     if ( rc/=0 ) then
        print *, 'cannot get Mie index for '//vname_(iq)
        return
     end if

     if (verbose==1) &
          print *, '[+] Adding '//trim(vname_(iq))//' contribution'

!    Loop over nobs, km, nch
!    --------------------------
     do i = 1, nobs
        do n = 1, nch
           do k =1, km
        
            call Chem_MieQuery(mieTables, idxTable, idxChannel(n), &
                               qc(k,iq,i), rh(k,i), tau=ext_,&
                               ssa=ssa_,bext=bext_, bbck=bsc_ )
 
                               ext(k,n,i) = ext(k,n,i) + ext_
                               sca(k,n,i) = sca(k,n,i) + ssa_ * ext_
                               bsc(k,n,i) = bsc(k,n,i) + bsc_*qc(k,iq,i)

            call Chem_MieQuery(mieTables, idxTable, idxChannel(n), &
                               qm(k,iq,i), rh(k,i), tau=tau_)
 
                               tau(k,n,i) = tau(k,n,i) + tau_            
            end do  ! end nch
         end do ! end km
      end do  ! end nobs
  end do ! end tracers

  bsc= bsc*1e03 ! in km-1 sr-1
  ext= ext*1e03 ! in km-1 
  sca= sca*1e03 ! in km-1   



!  Attenuated backscatter from space 
 
    absc_TOA(1,:,:) = bsc(1,:,:)*exp(-tau(1,:,:))
    do n = 1, nch
       do k = 2, km        
         do i = 1, nobs
           taulev = 0.
           do l = 1, k-1
              taulev = taulev + tau(l,n,i)
           enddo
           taulev = taulev + 0.5 *tau(k,n,i)
           absc_TOA(k,n,i) = bsc(k,n,i)*exp(-2.*taulev)
        enddo
      enddo
    enddo

!  Attenuated backscatter from surface 
   
    absc_SFC(km,:,:)= bsc(km,:,:)*exp(-tau(km,:,:))
    do n = 1, nch
       do k = km -1, 1, -1        
         do i = 1, nobs
           taulev = 0.
           do l = km, k+1, -1
              taulev = taulev + tau(l,n,i)
           enddo
        taulev = taulev + 0.5 * tau(k,n,i)
        absc_SFC(k,n,i) = bsc(k,n,i)*exp(-2.*taulev)
        enddo
      enddo
    enddo
 
end subroutine getExt


!..............................................................

! Read the OMI research retrieval required Fresnel Ocean Reflectivity LUT
! Function of wavelength (340, 380 nm) and nodes of SZA, VZA, VAA
  subroutine readFresnel (filename, ocean_ler)
  
  character(len=*), intent(in)             :: filename  ! LUT filename
  real, intent(out), dimension(2,16,16,16) :: ocean_ler ! LUT values

  integer                                  :: ix, iy, iz, iw
  real*4                                   :: datai

  
  open(unit=11,file=filename, form='unformatted')
  do iw = 1, 16
  do iz = 1, 16
  do iy = 1, 16
  do ix = 1, 2
   read(11) datai
   ocean_ler(ix,iy,iz,iw) = datai
  enddo
  enddo
  enddo
  enddo
  close(unit=11)

end subroutine readFresnel

