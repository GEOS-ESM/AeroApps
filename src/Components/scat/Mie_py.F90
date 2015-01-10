!
!  Simple f77 wrapper for the Python interface to the Mie Calculator.
!

Module Mie_py

CONTAINS

subroutine getMieGridded (filename, im, jm, km, nch, nMom, nPol, &
                          nymd, nhms, channels, verbose, &
                          tau, ssa, pe, ze, te, rc, &
                          g, pmom)
!
! Given a Chem Bundle file name, returns the aerosol optical depth (tau),
! single scattering albedo (ssa) and asymmetry factor (g). It requires
! that th following Chem Registry files be present in the current direcory:
!
! Chem_MieRegistry.rc --- used to decide which variable ro read from 
!                         Chem Bundle
! Aod_EOS.rc --- used to specify channels, location of Mie tables, etc.
!
  use Chem_MieMod
  use Chem_RegistryMod
  use Chem_BundleMod

  implicit NONE
  character(len=*), intent(in)  :: filename ! Chem bundle with aerosol
                                            ! mixing ratio
  integer,          intent(in)  :: im   ! number of longitudes on file
  integer,          intent(in)  :: jm   ! number of latitudes on file
  integer,          intent(in)  :: km   ! number of levels on file
  integer,          intent(in)  :: nch  ! number of channels
  integer,          intent(in)  :: nMom ! number of legender momemts 
  integer,          intent(in)  :: nPol ! number of components of the scattering matrix
  integer,          intent(in)  :: nymd ! date, eg., 20080629
  integer,          intent(in)  :: nhms ! time, eg., 120000
                                        ! wavelengths [nm] 
  real,             intent(in)  :: channels(nch) ! wavenumber [nm]
  integer,          intent(in)  :: verbose
  
  real,             intent(out) :: tau(im,jm,km,nch) ! aerosol optical depth
  real,             intent(out) :: ssa(im,jm,km,nch) ! single scattering albedo
  real, optional,   intent(out) :: g(im,jm,km,nch)   ! asymmetry factor
  real, optional,   intent(out) :: pmom(im,jm,km,nch,nMom,nPol) ! elements of scattering phase matrix
  real,             intent(out) :: pe(im,jm,km+1)    ! edge pressure [Pa]
  real,             intent(out) :: ze(im,jm,km+1)    ! edge height above sfc [m]
  real,             intent(out) :: te(im,jm,km+1)    ! edge Temperature [K]
  integer,          intent(out) :: rc ! return code
!                               ---
  type(Chem_Registry) :: chemReg
  type(Chem_Mie)      :: mieTables
  type(Chem_Bundle)   :: w_c 
  integer             :: i, j, k, m, n, iq, nymd_, nhms_, idxTable, fid, iPol, iMom
 
  real                :: alpha, tau_, ssa_, g_
  real, pointer       :: pmom_(:,:)
  real                :: idxChannel(nch) ! this should have been integer
  real                :: airdens(im,jm,km), pm(km), tm(km)
  
  integer, parameter  :: READ_ONLY=1

  real, parameter     :: grav = 9.80616
  real, parameter     :: Rgas = 287.
  real, parameter     :: kappa = 2.0 / 7.0

  

! Load the Chem Registry
! ----------------------
  chemReg = Chem_RegistryCreate(rc,'Chem_MieRegistry.rc')
  if(rc/=0) then
     print *, 'cannot create Chem Registry'
     return
  else
     if (verbose==1) &
          call Chem_RegistryPrint(chemReg)
  end if

! Create the Mie Tables
! ---------------------
  mieTables = Chem_MieCreate('Aod_EOS.rc',rc)
  if ( rc /= 0 ) then
     print *, 'Cannot create Mie tables from Aod_EOS.rc'
     return
  end if
!  print*, 'test1',mieTables%nMom
!  print*, 'test1',mieTables%nPol

  if ( nMom > mieTables%nMom ) then ! mieTables%nMom is writen in Aod_EOS.rc file
     print *, 'mieTables do not have enough moments', nMom, mieTables%nMom
     rc = 99
     return
  end if
  

! Allocate memory for phase function; size determined by what kind
! table has been loaded; n_moments given in Aod_EOS.rc
! ----------------------------------------------------------------
  allocate(pmom_(mieTables%nMom,mieTables%nPol),stat=rc)
  if ( rc /= 0 ) then
     print *, 'Cannot allocate memory for pmom_'
     return
  end if

! Read in Chem Bundle
! -------------------
  if ( nymd<0 ) then ! get last time on file
        call Chem_BundleRead(filename, nymd_, nhms_, w_c, rc, & 
                             ChemReg=chemReg, timidx=0)
  else ! get specified date/time
        nymd_ = nymd
        nhms_ = nhms
        call Chem_BundleRead(filename, nymd_, nhms_, w_c, rc, & 
                             ChemReg=chemReg )
  end if

  if (rc==0) then
     if (verbose==1) &
         print *, '[r] read Chem Bundle '//trim(filename)
  else
     print *, '[x] cannot read Chem Bundle '//trim(filename)
     return
  end if

! Read in air density for temperature calculation
! -----------------------------------------------
  call GFIO_Open ( filename, READ_ONLY, fid, rc )
  if (rc/=0) then
     print *,'cannot open '//trim(filename)
     return
  end if
  call GFIO_GetVar (fid, 'AIRDENS', nymd, nhms, im, jm, 1, km, airdens, rc )
  if (rc/=0) then
     print *,'cannot read AIRDENS from '//trim(filename)
     return
  end if
  call GFIO_Close ( fid, rc )
  if (rc/=0) then
     print *,'cannot close file '//trim(filename)
     return
  end if

  
! Check consistency of dimensions
! -------------------------------
  if ( im /= w_c%grid%im .OR. &
       jm /= w_c%grid%jm .OR. &
       km /= w_c%grid%km  ) then
     print *, 'expecting im, jm, km = ', im, jm, km
     print *, 'but found im, jm, km = ', w_c%grid%im, w_c%grid%jm, w_c%grid%km
     rc = 97
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
  if ( present(g) ) g = 0.0
  if ( present(pmom) ) pmom = 0.0
     
  
! Loop over aerosol species
! -------------------------
  do iq = 1,chemReg%nq
!     print*, 'name', chemReg%nq, chemReg%vname(iq)
     idxTable = Chem_MieQueryIdx(mieTables,chemReg%vname(iq),rc)
     if(idxTable == -1) cycle
     if ( rc/=0 ) then
        print *, 'cannot get Mie index for '//chemReg%vname(iq)
        return
     end if

     if (verbose==1) &
          print *, '[+] Adding '//trim(chemReg%vname(iq))//' contribution'
     
  
!    Loop over channel, x, y, z
!    --------------------------
     do n = 1, nch
        do k = 1, km
           do j = 1, jm
              do i = 1, im
                           
              call Chem_MieQuery(mieTables, idxTable, idxChannel(n),    &
                         w_c%qa(iq)%data3d(i,j,k)*w_c%delp(i,j,k)/grav, &
                         w_c%rh(i,j,k), tau=tau_, ssa=ssa_, gasym=g_)

                            tau(i,j,k,n)    = tau(i,j,k,n)    + tau_
                            ssa(i,j,k,n)    = ssa(i,j,k,n)    + ssa_ * tau_ 
              if ( present(g) ) then      
                            g(i,j,k,n)      =   g(i,j,k,n)    + g_ * ssa_ * tau_
              end if

              if ( present(pmom) ) then
                   call Chem_MieQuery(mieTables, idxTable, idxChannel(n),    &
                         w_c%qa(iq)%data3d(i,j,k)*w_c%delp(i,j,k)/grav, &
                         w_c%rh(i,j,k),tau=tau_, ssa=ssa_, pmom=pmom_)

                         do ipol=1, nPol      
                            do imom=1, nMom   
             pmom(i,j,k,n,imom,ipol) = pmom(i,j,k,n,imom,ipol) + pmom_(imom,ipol) * ssa_ * tau_  
                            end do
                         end do
              end if
                          
              end do ! longitudes
          end do ! latitudes
       end do ! levels
    end do ! channels
        
   end do ! aerosol tracers

! Normalize ssa and g
! -------------------
   where (     tau > 0.0 ) ssa = ssa / tau
   where ( ssa*tau > 0.0 )  g =  g / ( ssa * tau ) 


   if (present (pmom)) then
    do n = 1, nch
       do k = 1, km
          do j = 1, jm
             do i = 1, im
 
                do ipol=1, nPol   ! normalize  Pmom 
                   do imom=1, nMom 
                   if (( ssa(i,j,k,n)  * tau(i,j,k,n)  ) > 0.0 ) then
                   pmom(i,j,k,n,imom,ipol)  = pmom(i,j,k,n,imom,ipol)  / ( ssa(i,j,k,n)  * tau(i,j,k,n)  ) 
                   end if
                   end do
                end do
             
             end do
          end do
      end do
   end do
  end if

! Finally, compute mid-layer pressures
! ------------------------------------
  do j = 1, jm
     do i = 1, im

        pe(i,j,1) = w_c%grid%ptop
        do k = 1, km
           pe(i,j,k+1) = pe(i,j,k) + w_c%delp(i,j,k)
        end do

        ze(i,j,km+1) = 0.0 ! height above surface
        do k = km, 1, -1
           ze(i,j,k) = ze(i,j,k+1) + w_c%delp(i,j,k) / ( airdens(i,j,k) * grav ) 
        end do

        do k = 1, km
           pm(k) = ( pe(i,j,k) + pe(i,j,k+1) ) / 2.0
           tm(k) = pm(k) / ( airdens(i,j,k) * Rgas )
        end do

        te(i,j,1) = tm(1) ! isothermal at highest level
        do k = 2, km
           alpha = log(pe(i,j,k)/pm(k-1))/log(pm(k)/pm(k-1))
           te(i,j,k) = tm(k-1) + alpha * (tm(k) - tm(k-1) )
        end do
                              ! dry adiabatic
        te(i,j,km+1) = tm(km) * (pe(i,j,km+1)/pm(km))**kappa

     end do
  end do

! Chem bundle no longer needed
! ----------------------------
  call Chem_BundleDestroy(w_c, rc)
  if ( rc /= 0 ) then
     print *, 'Cannot destroy Chem_Bundle'
     return
  end if
  deallocate(pmom_)
  
  call Chem_MieDestroy(mieTables, rc)
  if ( rc /= 0 ) then
     print *, 'Cannot destroy MieTables'
     return
  end if


  if (verbose==1) &
       print *, '[x] All done!'

end subroutine getMieGridded
!.............................................................................

subroutine getMieObs(filename, km, nch, nMom, nPol,nobs, nymd, nhms, channels, &
                      lon, lat, verbose, tau, ssa, pe, ze, te, rc,   &
                      g, pmom )
!
! Given a Chem Bundle file name and lat/lon coordinates at observation
! locatpns, returns the aerosol optical depth (tau), single scattering 
! albedo (ssa) and asymmetry factor (g) at these locations. It requires
! that the following Chem Registry files be present in the current directory:
!
! Chem_MieRegistry.rc --- used to decide which variable ro read from 
!                         Chem Bundle
! Aod_EOS.rc --- used to specify channels, location of Mie tables, etc.
!

  implicit NONE
  character(len=*), intent(in)  :: filename ! Chem bundle with aerosol
                                            ! mixing ratio
  integer,          intent(in)  :: km   ! number of levels on file
  integer,          intent(in)  :: nch   ! number of channels
  integer,          intent(in)  :: nMom ! number of momemts legender
  integer,          intent(in)  :: nPol ! number of components of the scattering matrix
  integer,          intent(in)  :: nobs ! number of observations
  integer,          intent(in)  :: nymd ! date, eg., 20080629
  integer,          intent(in)  :: nhms ! time, eg., 120000
                                        ! wavelengths [nm] 
  real,             intent(in)  :: channels(nch)

  real,             intent(in)  :: lat(nobs)
  real,             intent(in)  :: lon(nobs)
 
 
  integer,          intent(in)  :: verbose

  real,             intent(out) :: tau(km,nch,nobs) ! aerosol optical depth
  real,             intent(out) :: ssa(km,nch,nobs) ! single scattering albedo
  real, optional,   intent(out) :: g(km,nch,nobs)   ! asymmetry factor
  real, optional,   intent(out) :: pmom(km,nch,nobs,nMom,nPol) ! phase function 
 
  real,             intent(out) :: pe(km+1,nobs)    ! edge pressure [Pa]
  real,             intent(out) :: ze(km+1,nobs)    ! edge height above sfc [m]
  real,             intent(out) :: te(km+1,nobs)    ! edge Temperature [K]

  integer,          intent(out) :: rc ! return code
!                               ---


  integer :: im, jm, km_, tm
  real, allocatable :: tau_(:,:,:,:)
  real, allocatable :: ssa_(:,:,:,:)
  real, allocatable ::   g_(:,:,:,:)
  real, allocatable :: pmom_(:,:,:,:,:,:)
  real, allocatable ::  pe_(:,:,:)
  real, allocatable ::  ze_(:,:,:)
  real, allocatable ::  te_(:,:,:)
 
  
! Get size of gridded fields
! --------------------------
  call getDims (filename, im, jm, km_, tm, rc)
  if ( rc /= 0 ) return

  if ( km /= km_ ) then
     print *, 'Inconsistent vertical dimensions: ', km, km_
     rc = 1
     return
  end if


! Allocate necessary space
! ------------------------
  allocate( tau_(im,jm,km,nch), ssa_(im,jm,km,nch), g_(im,jm,km,nch), &
            pmom_(im,jm,km,nch,nMom,nPol), pe_(im,jm,km+1), ze_(im,jm,km+1), te_(im,jm,km+1), stat=rc)
  if ( rc /= 0 ) return

! Read in gridded Mie parameters
! ------------------------------
  if (present(pmom)) then 

  call getMieGridded (filename, im, jm, km, nch, nmom, &
                      npol, nymd, nhms, channels, verbose, &
                      tau_, ssa_, pe_, ze_, te_, rc, &
                      g=g_, pmom=pmom_)

  else
  call getMieGridded (filename, im, jm, km, nch, nmom, &
                      npol, nymd, nhms, channels, verbose, &
                      tau_, ssa_, pe_, ze_, te_, rc, &
                      g=g_) 
  end if


! Next interpolate to obs locations
! ---------------------------------
  call Interpxy_ ( lon, lat, nobs, im, jm, km, nch, tau_, tau )
  call Interpxy_ ( lon, lat, nobs, im, jm, km, nch, ssa_, ssa )
  call Interpxy_ ( lon, lat, nobs, im, jm, km+1, 1,  pe_,  pe )
  call Interpxy_ ( lon, lat, nobs, im, jm, km+1, 1,  ze_,  ze )
  call Interpxy_ ( lon, lat, nobs, im, jm, km+1, 1,  te_,  te )
 
  if (present(g)) then
  call Interpxy_ ( lon, lat, nobs, im, jm, km, nch,   g_,   g )
  end if
  
  if (present(pmom)) then
  call Interpxy_bis ( lon, lat, nobs,nMom, nPol, im, jm, km, nch, pmom_, pmom )
  end if

! De-allocate memory
! ------------------
  deallocate( tau_, ssa_, g_, pmom_, pe_, ze_, te_ )

end subroutine getMieObs

!..........................................................................

subroutine getAerObs(filename, im, jm, km, nobs, nymd, nhms, &
                      lat, lon, verbose, qm, rh, pe, ze, te, rc)

!
! Given a Chem Bundle file name and lat/lon coordinates at observation
! locatpns, returns the mass mixing ratio and the relative humidity at these locations. 
! It requires that th following Chem Registry files be present in the current direcory:
!
! Chem_MieRegistry.rc --- used to decide which variable ro read from 
!                         Chem Bundle

  use Chem_BundleMod
  use Chem_RegistryMod

  implicit NONE
  character(len=*), intent(in)  :: filename ! Chem bundle with aerosol
                                            ! mixing ratio

  integer,          intent(in)  :: im   ! number of longitudes on file
  integer,          intent(in)  :: jm   ! number of latitudes on file
  integer,          intent(in)  :: km   ! number of levels on file
  integer,          intent(in)  :: nobs ! number of observations
  integer,          intent(in)  :: nymd ! date, eg., 20080629
  integer,          intent(in)  :: nhms ! time, eg., 120000
                                        ! wavelengths [nm] 

  real,             intent(in)  :: lat(nobs)
  real,             intent(in)  :: lon(nobs)
 
 
  integer,          intent(in)  :: verbose
                                                  ! 18 : number of tracers
  real,             intent(out) :: qm(18,km,nobs) ! mass mixing ratio
  real,             intent(out) :: rh(18,km,nobs) ! relative humidity

  real,             intent(out) :: pe(km+1,nobs)    ! edge pressure [Pa]
  real,             intent(out) :: ze(km+1,nobs)    ! edge height above sfc [m]
  real,             intent(out) :: te(km+1,nobs)    ! edge Temperature [K]

  integer,          intent(out) :: rc ! return code

  
!                               ---
  type(Chem_Bundle)   :: w_c
  type(Chem_Registry) :: chemReg
  integer             ::  km_, fid, iq,i,j,k,nymd_, nhms_ 
  real                ::  alpha

  real, allocatable   ::  pe_(:,:,:)
  real, allocatable   ::  ze_(:,:,:)
  real, allocatable   ::  te_(:,:,:)
  
  real                :: qm_(18,im,jm,km)
  real                :: rh_(18,im,jm,km)

  integer, parameter  :: READ_ONLY=1
  real                :: airdens(im,jm,km), pm(km), tm(km)
  real, parameter     :: grav = 9.80616
  real, parameter     :: Rgas = 287.
  real, parameter     :: kappa = 2.0 / 7.0

! Allocate necessary space
! ------------------------
  allocate( pe_(im,jm,km+1), ze_(im,jm,km+1), te_(im,jm,km+1),stat=rc )
  if ( rc /= 0 ) return

! Load the Chem Registry
! ----------------------
  chemReg = Chem_RegistryCreate(rc,'Chem_MieRegistry.rc')
  if(rc/=0) then
     print *, 'cannot create Chem Registry'
     return
  else
     if (verbose==1) &
          call Chem_RegistryPrint(chemReg)
  end if

! Get size of gridded fields
! --------------------------
  call getDims (filename, im, jm, km_, tm, rc)
  if ( rc /= 0 ) return

  if ( km /= km_ ) then
     print *, 'Inconsistent vertical dimensions: ', km, km_
     rc = 1
     return
  end if


! Read in Chem Bundle
! -------------------
  if ( nymd<0 ) then ! get last time on file
        call Chem_BundleRead(filename, nymd_, nhms_, w_c, rc, & 
                             ChemReg=chemReg, timidx=0)
  else ! get specified date/time
        nymd_ = nymd
        nhms_ = nhms
        call Chem_BundleRead(filename, nymd_, nhms_, w_c, rc, & 
                             ChemReg=chemReg )
  end if
  

  if (rc==0) then
     if (verbose==1) &
         print *, '[r] read Chem Bundle '//trim(filename)
  else
     print *, '[x] cannot read Chem Bundle '//trim(filename)
     return
  end if

! Check consistency of dimensions
! -------------------------------
  if ( im /= w_c%grid%im .OR. &
       jm /= w_c%grid%jm .OR. &
       km /= w_c%grid%km  ) then
     print *, 'expecting im, jm, km = ', im, jm, km
     print *, 'but found im, jm, km = ', w_c%grid%im, w_c%grid%jm, w_c%grid%km
     rc = 97
     return
  end if


  do iq = 1,chemReg%nq !aerosol tracers
     
     if (verbose==1) &
          print *, '[+] Adding '//trim(chemReg%vname(iq))//' contribution'
    
!    Loop over channel, x, y, z
!    --------------------------
        do k = 1, km
           do j = 1, jm
              do i = 1, im
    
              qm_(iq,i,j,k) = w_c%qa(iq)%data3d(i,j,k)*w_c%delp(i,j,k)/grav 
              rh_(iq,i,j,k) = w_c%rh(i,j,k)
              
              end do
           end do
        end do
  end do ! end tracers

! Next interpolate to obs locations
! ---------------------------------
  
  call Interpxy_aer ( lon, lat, nobs, im, jm, km, chemReg%nq,   1,  qm_, qm )
  call Interpxy_aer ( lon, lat, nobs, im, jm, km, chemReg%nq,  1,  rh_, rh )

! Read in air density for temperature calculation
! -----------------------------------------------
  call GFIO_Open ( filename, READ_ONLY, fid, rc )
  if (rc/=0) then
     print *,'cannot open '//trim(filename)
     return
  end if
  call GFIO_GetVar (fid, 'AIRDENS', nymd, nhms, im, jm, 1, km, airdens, rc )
  if (rc/=0) then
     print *,'cannot read AIRDENS from '//trim(filename)
     return
  end if
  call GFIO_Close ( fid, rc )
  if (rc/=0) then
     print *,'cannot close file '//trim(filename)
     return
  end if


!  Compute mid-layer pressures
! ------------------------------------
  do j = 1, jm
     do i = 1, im

        pe_(i,j,1) = w_c%grid%ptop
        do k = 1, km
           pe_(i,j,k+1) = pe_(i,j,k) + w_c%delp(i,j,k)
        end do

        ze_(i,j,km+1) = 0.0 ! height above surface
        do k = km, 1, -1
           ze_(i,j,k) = ze_(i,j,k+1) + w_c%delp(i,j,k) / ( airdens(i,j,k) * grav ) 
        end do

        do k = 1, km
           pm(k) = ( pe_(i,j,k) + pe_(i,j,k+1) ) / 2.0
           tm(k) = pm(k) / ( airdens(i,j,k) * Rgas )
        end do

        te_(i,j,1) = tm(1) ! isothermal at highest level
        do k = 2, km
           alpha = log(pe_(i,j,k)/pm(k-1))/log(pm(k)/pm(k-1))
           te_(i,j,k) = tm(k-1) + alpha * (tm(k) - tm(k-1) )
        end do
                              ! dry adiabatic
        te_(i,j,km+1) = tm(km) * (pe_(i,j,km+1)/pm(km))**kappa
  end do
  end do

! Next interpolate to obs locations
! ---------------------------------
  
  call Interpxy_ ( lon, lat, nobs, im, jm, km+1, 1,  pe_,  pe )
  call Interpxy_ ( lon, lat, nobs, im, jm, km+1, 1,  ze_,  ze )
  call Interpxy_ ( lon, lat, nobs, im, jm, km+1, 1,  te_,  te )
  
! De-allocate memory
! ------------------
  deallocate( pe_, ze_, te_)

! Chem bundle no longer needed
! ----------------------------
  call Chem_BundleDestroy(w_c, rc)
  if ( rc /= 0 ) then
     print *, 'Cannot destroy Chem_Bundle'
     return
  end if
  
  call Chem_RegistryDestroy(chemReg, rc)
  if ( rc /= 0 ) then
     print *, 'Cannot destroy Chemisty Registry'
     return
  end if
  
  if (verbose==1) &
       print *, '[x] All done!'

end subroutine getAerObs

!.............................................................................

 subroutine getMieObs2 (filename, km, nch, nMom, nPol, &
                        nobs, nymd, nhms, lon,lat, channels, verbose, &
                        tau, ssa, pe, ze, te, rc, g, pmom)
!
! Given a Chem Bundle file name, returns the mass mixing ratio and relative humidity . It requires
! that th following Chem Registry files be present in the current direcory:
!
! Chem_MieRegistry.rc --- used to decide which variable ro read from 
!                         Chem Bundle
! Aod_EOS.rc --- used to specify channels, location of Mie tables, etc.
!
  use Chem_MieMod
  use Chem_RegistryMod
  use Chem_BundleMod

  implicit NONE
  character(len=*), intent(in)  :: filename ! Chem bundle with aerosol
                                            ! mixing ratio
  integer,          intent(in)  :: km   ! number of levels on file
  integer,          intent(in)  :: nch  ! number of channels
  integer,          intent(in)  :: nMom ! number of legender momemts 
  integer,          intent(in)  :: nPol ! number of components of the scattering matrix
  integer,          intent(in)  :: nobs! number of components of the scattering matrix
  integer,          intent(in)  :: nymd ! date, eg., 20080629
  integer,          intent(in)  :: nhms ! time, eg., 120000
                                        ! wavelengths [nm] 
  real,             intent(in)  :: lat(nobs)
  real,             intent(in)  :: lon(nobs)
 


  real,             intent(in)  :: channels(nch) ! wavenumber [nm]
  integer,          intent(in)  :: verbose

  real,             intent(out) :: tau(km,nch,nobs) ! aerosol optical depth
  real,             intent(out) :: ssa(km,nch,nobs) ! single scattering albedo
  real, optional,   intent(out) :: g(km,nch,nobs)   ! asymmetry factor
  real, optional,   intent(out) :: pmom(km,nch,nobs,nMom,nPol) ! elements of scattering phase matrix
  real,             intent(out) :: pe(km+1,nobs)    ! edge pressure [Pa]
  real,             intent(out) :: ze(km+1,nobs)    ! edge height above sfc [m]
  real,             intent(out) :: te(km+1,nobs)    ! edge Temperature [K]
  integer,          intent(out) :: rc ! return code
!                               ---
  type(Chem_Registry) :: chemReg
  type(Chem_Mie)      :: mieTables

 
  integer             :: i, n, k, m, iq, nymd_, nhms_, idxTable, fid, iPol, iMom
  integer             :: im, jm, km_, tm
  real                :: alpha, tau_, ssa_, g_
  real, pointer       :: pmom_(:,:)
  real                :: idxChannel(nch) ! this should have been integer
  
  
  integer, parameter  :: READ_ONLY=1
  real, allocatable   :: qm(:,:,:)
  real, allocatable   :: rh(:,:,:)

! Allocate necessary space
! ------------------------
  allocate( qm(18,km,nobs), rh(18,km,nobs), stat=rc)
  if ( rc /= 0 ) return

! Load the Chem Registry
! ----------------------
  chemReg = Chem_RegistryCreate(rc,'Chem_MieRegistry.rc')
  if(rc/=0) then
     print *, 'cannot create Chem Registry'
     return
  else
     if (verbose==1) &
          call Chem_RegistryPrint(chemReg)
  end if

! Get size of gridded fields
! --------------------------
  call getDims (filename, im, jm, km_, tm, rc)
  if ( rc /= 0 ) return

  if ( km /= km_ ) then
     print *, 'Inconsistent vertical dimensions: ', km, km_
     rc = 1
     return
  end if
  

! Create the Mie Tables
! ---------------------
  mieTables = Chem_MieCreate('Aod_EOS.rc',rc)
  if ( rc /= 0 ) then
     print *, 'Cannot create Mie tables from Aod_EOS.rc'
     return
  end if
!  print*, 'test1',mieTables%nMom
!  print*, 'test1',mieTables%nPol

  if ( nMom > mieTables%nMom ) then ! mieTables%nMom is writen in Aod_EOS.rc file
     print *, 'mieTables do not have enough moments', nMom, mieTables%nMom
     rc = 99
     return
  end if
  

! Allocate memory for phase function; size determined by what kind
! table has been loaded; n_moments given in Aod_EOS.rc
! ----------------------------------------------------------------
  allocate(pmom_(mieTables%nMom,mieTables%nPol),stat=rc)
  if ( rc /= 0 ) then
     print *, 'Cannot allocate memory for pmom_'
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
  if ( present(g) ) g = 0.0
  if ( present(pmom) ) pmom = 0.0

! return the qm and rh for each obs
! ---------------------------------
     
  call getAerObs (filename, im, jm, km, nobs, nymd, nhms, lat, lon,  verbose, &
                  qm, rh, pe, ze, te,rc)
 
! Loop over aerosol species
! -------------------------
  do iq = 1,chemReg%nq
!     print*, 'name', chemReg%nq, chemReg%vname(iq)
     idxTable = Chem_MieQueryIdx(mieTables,chemReg%vname(iq),rc)
     if(idxTable == -1) cycle
     if ( rc/=0 ) then
        print *, 'cannot get Mie index for '//chemReg%vname(iq)
        return
     end if

     if (verbose==1) &
          print *, '[+] Adding '//trim(chemReg%vname(iq))//' contribution'

!    Loop over nobs, km, nch
!    --------------------------
     
      do i = 1, nobs
         do k =1, km
            do n = 1, nch
        
            call Chem_MieQuery(mieTables, idxTable, idxChannel(n), &
                            qm(iq,k,i), rh(iq,k,i), tau=tau_, ssa=ssa_, gasym=g_)
 
                            tau(k, n,i)    = tau(k,n,i)    + tau_
                            ssa(k, n,i)    = ssa(k,n,i)    + ssa_ * tau_ 
            if ( present(g) ) then      
                            g(k, n,i)       =   g(k,n,i)    + g_ * ssa_ * tau_
            end if

            if ( present(pmom) ) then
                    
              call Chem_MieQuery(mieTables, idxTable, idxChannel(n), &
                         qm(iq,k,i), rh(iq,k,i), tau=tau_, ssa=ssa_, pmom=pmom_)
                do ipol=1, nPol      
                   do imom=1, nMom   
                   pmom(k,n,i,imom,ipol) = pmom(k,n,i,imom,ipol) + pmom_(imom,ipol) * ssa_ * tau_  
                   end do
                end do
            end if  
            end do  ! end nch
         end do ! end km
      end do  ! end nobs
  end do ! end tracers

! Normalize ssa and g
! -------------------
   where (     tau > 0.0 ) ssa = ssa / tau
   where ( ssa*tau > 0.0 )  g =  g / ( ssa * tau ) 


   if (present (pmom)) then
    do i = 1, nobs
       do n = 1, nch
          do k =1, km
 
            do ipol=1, nPol   ! normalize Pmom   
               do imom=1, nMom 
                  
              if (( ssa(k,n,i)  * tau(k,n,i)  ) > 0.0 ) then
              pmom(k,n,i,imom,ipol)  = pmom(k,n,i,imom,ipol)  / ( ssa(k,n,i)  * tau(k,n,i)  ) 
              end if
               
               end do
            end do
          end do
      end do
    end do
  end if
! De-allocate memory
! ------------------
  deallocate( qm, rh, pmom_ )
   
  call Chem_MieDestroy(mieTables, rc)
  if ( rc /= 0 ) then
     print *, 'Cannot destroy MieTables'
     return
  end if

  call Chem_RegistryDestroy(chemReg, rc)
  if ( rc /= 0 ) then
     print *, 'Cannot destroy Chemisty Registry'
     return
  end if

end subroutine getMieObs2

end Module Mie_py

!..............................................................

subroutine getDims (filename, im, jm, km, tm, rc)
!
! Return dimensions on GFIO-compatible file.
!
  implicit NONE
  character(len=*), intent(in)  :: filename ! input GFIO compatible file name
  integer,          intent(out) :: im ! number of longitudes on file
  integer,          intent(out) :: jm ! number of latitudes on file
  integer,          intent(out) :: km ! number of levels on file
  integer,          intent(out) :: tm ! number of times on file
  integer,          intent(out) :: rc ! return code
!                   ---
  integer, parameter :: READ_ONLY=1
  integer :: fid, nvars, ngatts
  call GFIO_Open ( filename, READ_ONLY, fid, rc )
  if (rc/=0) then
     print *,'cannot open '//trim(filename)
     return
  end if
  call GFIO_DimInquire ( fid, im, jm, km, tm, nvars, ngatts, rc)
  if (rc/=0) then
     print *,'cannot get dimensions of file '//trim(filename)
     return
  end if
  call GFIO_Close ( fid, rc )
  if (rc/=0) then
     print *,'cannot close file '//trim(filename)
     return
  end if
end subroutine getDims

!.............................................................................
subroutine getmiedims(mieTableFile, nCh, nRh, nBin, nMom, nPol, rc)
!                                                                                            
! Retrieve dimensions of MieTable                                                            
!                                                                                            
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

end subroutine getmiedims


!.............................................................................
 

subroutine scalar(filename, km, nch, nMom, nPol, nobs, nymd, nhms, &
                  channels,lon, lat, verbose, tau, ssa, pe, ze,   &
                  te, g, rc)

  use Mie_py 
  character(len=*), intent(in)  :: filename ! Chem bundle with aerosol
                                            ! mixing ratio
  integer,          intent(in)  :: km   ! number of levels on file
  integer,          intent(in)  :: nch   ! number of channels
  integer,          intent(in)  :: nobs ! number of observations
  integer,          intent(in)  :: nMom ! number of momemts legender
  integer,          intent(in)  :: nPol ! number of components of the scattering matrix
  integer,          intent(in)  :: nymd ! date, eg., 20080629
  integer,          intent(in)  :: nhms ! time, eg., 120000
                                        ! wavelengths [nm] 
  real,             intent(in)  :: channels(nch)

  real,             intent(in)  :: lat(nobs)
  real,             intent(in)  :: lon(nobs)
 
 
  integer,          intent(in)  :: verbose

  real,             intent(out) :: tau(km,nch,nobs) ! aerosol optical depth
  real,             intent(out) :: ssa(km,nch,nobs) ! single scattering albedo
  real,             intent(out) :: g(km,nch,nobs)   ! asymmetry factor
 
  real,             intent(out) :: pe(km+1,nobs)    ! edge pressure [Pa]
  real,             intent(out) :: ze(km+1,nobs)    ! edge height above sfc [m]
  real,             intent(out) :: te(km+1,nobs)    ! edge Temperature [K]

  integer,          intent(out) :: rc ! return code
 
! qm gridded -> Mie parameters gridded -> Mie parameters for each obs 
!  call getMieObs(filename, km, nch, nMom, nPol, nobs, nymd, nhms, &
!                channels,lon, lat, verbose, tau, ssa, pe, ze, &
!                 te, rc, g = g)


! qm gridded -> qm for each obs -> Mie parameters for each obs 
  call getMieObs2(filename, km, nch, nMom, nPol, &
                  nobs, nymd, nhms, lon,lat, channels, verbose, &
                  tau, ssa, pe, ze, te, rc, g=g)

  
  

 end subroutine scalar

!.............................................................................
subroutine vector(filename, km, nch, nMom, nPol, nobs, nymd, nhms, &
                  channels,lon, lat, verbose, tau, ssa, pe, ze,   &
                  te, g, pmom, rc )
 
  use Mie_py 
  character(len=*), intent(in)  :: filename ! Chem bundle with aerosol
                                            ! mixing ratio
  integer,          intent(in)  :: km   ! number of levels on file
  integer,          intent(in)  :: nch   ! number of channels
  integer,          intent(in)  :: nobs ! number of observations
  integer,          intent(in)  :: nMom ! number of momemts legender
  integer,          intent(in)  :: nPol ! number of components of the scattering matrix
  integer,          intent(in)  :: nymd ! date, eg., 20080629
  integer,          intent(in)  :: nhms ! time, eg., 120000
                                        ! wavelengths [nm] 
  real,             intent(in)  :: channels(nch)

  real,             intent(in)  :: lat(nobs)
  real,             intent(in)  :: lon(nobs)
 
 
  integer,          intent(in)  :: verbose

  real,             intent(out) :: tau(km,nch,nobs) ! aerosol optical depth
  real,             intent(out) :: ssa(km,nch,nobs) ! single scattering albedo
  real, optional,   intent(out) :: g(km,nch,nobs)   ! asymmetry factor
  real, optional,   intent(out) :: pmom(km,nch,nobs,nMom,nPol) ! phase function 
 
  real,             intent(out) :: pe(km+1,nobs)    ! edge pressure [Pa]
  real,             intent(out) :: ze(km+1,nobs)    ! edge height above sfc [m]
  real,             intent(out) :: te(km+1,nobs)    ! edge Temperature [K]

  integer,          intent(out) :: rc ! return code

! qm gridded -> Mie parameters gridded -> Mie parameters for each obs 
! call getMieObs(filename, km, nch, nMom, nPol, nobs, nymd, nhms, &
 !               channels,lon, lat, verbose, tau, ssa, pe, ze, &
 !               te, g = g, pmom = pmom, rc=rc)

! qm gridded -> qm for each obs -> Mie parameters for each obs
 call getMieObs2(filename, km, nch, nMom, nPol, &
                 nobs, nymd, nhms, lon, lat, channels, verbose, &
                 tau, ssa, pe, ze, te, g = g, pmom = pmom, rc=rc)
 
 end subroutine vector

!.............................................................................

subroutine getTopObs (filename, nobs, lon, lat, zs, rc )
!
! Get surface topography at obs location.
!

  implicit NONE
  character(len=*), intent(in)  :: filename ! Chem bundle with aerosol
                                            ! mixing ratio
  integer,          intent(in)  :: nobs ! number of observations

  real,             intent(in)  :: lat(nobs)
  real,             intent(in)  :: lon(nobs)
 
  real,             intent(out) :: zs(nobs) ! pressure at mid-layer [Pa]

  integer,          intent(out) :: rc ! return code

!                       ---

  integer, parameter :: READ_ONLY=1
  integer :: nymd, nhms, incSecs
  integer :: im, jm, km, tm, nvars, ngatts, fid
  real, allocatable ::  zs_(:,:)

  call GFIO_Open ( filename, READ_ONLY, fid, rc )
  if (rc/=0) return

  call GFIO_DimInquire ( fid, im, jm, km, tm, nvars, ngatts, rc)
  if (rc/=0) return

  allocate(zs_(im,jm),stat=rc)
  if ( rc/=0 ) return

  call GetBegDateTime ( fid, nymd, nhms, incSecs, rc )
  if ( rc/=0 ) return

  call GFIO_GetVar (fid, 'zs', nymd, nhms, im, jm, 0, 1, zs_, rc )
  if (rc/=0) return

  call GFIO_Close ( fid, rc )
  if (rc/=0) return

  call Interpxy_ ( lon, lat, nobs, im, jm, 1, 1, zs_, zs )

  deallocate(zs_)

  end subroutine getTopObs

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Interp_--- Interpolates 3D field to observation locations
!
! !INTERFACE:
!
  subroutine Interpxy_ ( lon, lat, nobs, im, jm, km, nch, gField, &
                         oField )

! !USES:
!
      Implicit NONE

! !INPUT PARAMETERS:
!

      integer, intent(in)        :: nobs      ! Number of observations
      real,    intent(in)        :: lon(nobs) ! longitude in degrees [-180,+180]
      real,    intent(in)        :: lat(nobs) ! latitude  in degrees [-90,90]

      integer, intent(in)        :: im        ! zonal dimension
      integer, intent(in)        :: jm        ! meridional dimension
      integer, intent(in)        :: km        ! vertical dimension: 
                                              ! = 1 for 2D fields
                                              ! = km for mid-layer fields
                                              ! = km+1 for edge fields
      integer, intent(in)        :: nch        ! number of channels

                                              ! Gridded Field
      real,    intent(in)        :: gField(im,jm,km,nch) 

! !OUTPUT PARAMETERS:
!
                                              ! Interpolated profile
      real,    intent(out)       :: oField(km,nch,nobs)

! !DESCRIPTION: This routine interpolates gridded model fields to observation
!               locations. This routine implements only the horizontal
!  interpolation.
!
!  IMPORTANT:   The input lon coordinates must be in [-180, 180].
!               The input field cannot ahve any UNDEFs.
!
! !SEE ALSO:
!
!              Module m_insitu which uses the same linear interpolation algorithm.
!
!
! !REVISION HISTORY:
!
!  10feb2010  da Silva  Simplified m_interp routine for profile interpolation.
!
!EOP
!-------------------------------------------------------------------------
 
    character(len=*), parameter :: myname = 'Interpxy_'


! Local
      integer i, j, k, nob
      real    o_lon, o_lat
      real    m_dlon, m_dlat
      real    alfa, beta

      real a11(km,nch)         !W-S
      real a12(km,nch)         !W-N
      real a21(km,nch)         !E-S
      real a22(km,nch)         !E-N
      real a00(km,nch)         !temp

      integer i1, i2

!                         ------

     if ( nobs .eq. 0 ) return    ! nothing to do (keep this)

     m_dlon = float(im) / 360.
     m_dlat = float(jm-1) / 180.
     
!    Loop over observations
!    ----------------------
     do nob = 1, nobs

!       Longitude
!       ---------
        o_lon = 1. + (lon(nob)+180.) * m_dlon
        i   = min(im, int( o_lon ))
        alfa  = o_lon - i
        if(i .eq. im) then
           i1 = im
           i2 = 1
        else
           i1 = i
           i2 = i + 1
        endif
        
!       Latitude
!       --------
        o_lat = 1. + (lat(nob) + 90.) * m_dlat
        j   = min( jm-1, int( o_lat ) )
        beta  = o_lat - j
        
        a11 = gField(i1,j,  :,:)
        a21 = gField(i2,j,  :,:)
        a12 = gField(i1,j+1,:,:)
        a22 = gField(i2,j+1,:,:)
        a00 = a11 + alfa * ( a21 - a11 )

        oField(:,:,nob) = a00 + beta * ( a12 + alfa * ( a22 - a12 )  -  a00 )
     
  end do ! loop over obs
  
  return
  
end subroutine Interpxy_



!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Interp_--- Interpolates 3D field to observation locations
!
! !INTERFACE:
!
  subroutine Interpxy_bis ( lon, lat, nobs,nMom, nPol, im, jm, km, nch, gField, &
                         oField )

! !USES:
!
      Implicit NONE

! !INPUT PARAMETERS:
!

      integer, intent(in)        :: nobs      ! Number of observations
      real,    intent(in)        :: lon(nobs) ! longitude in degrees [-180,+180]
      real,    intent(in)        :: lat(nobs) ! latitude  in degrees [-90,90]

      integer, intent(in)        :: im        ! zonal dimension
      integer, intent(in)        :: jm        ! meridional dimension
      integer, intent(in)        :: km        ! vertical dimension: 
                                              ! = 1 for 2D fields
                                              ! = km for mid-layer fields
                                              ! = km+1 for edge fields
      integer, intent(in)        :: nch        ! number of channels
      integer,  intent(in)       :: nMom
      integer,  intent(in)       :: nPol

                                              ! Gridded Field
      real,    intent(in)        :: gField(im,jm,km,nch,nMom,nPol) 

     

! !OUTPUT PARAMETERS:
!
                                              ! Interpolated profile
      real,    intent(out)       :: oField(km,nch,nobs,nMom,nPol)

! !DESCRIPTION: This routine interpolates gridded model fields to observation
!               locations. This routine implements only the horizontal
!  interpolation.
!
!  IMPORTANT:   The input lon coordinates must be in [-180, 180].
!               The input field cannot ahve any UNDEFs.
!
! !SEE ALSO:
!
!              Module m_insitu which uses the same linear interpolation algorithm.
!
!
! !REVISION HISTORY:
!
!  10feb2010  da Silva  Simplified m_interp routine for profile interpolation.
!
!EOP
!-------------------------------------------------------------------------
 
    character(len=*), parameter :: myname = 'Interpxy_bis'


! Local
      integer i, j, k, nob
      real    o_lon, o_lat
      real    m_dlon, m_dlat
      real    alfa, beta

      real a11(km,nch, nMom, nPol)         !W-S
      real a12(km,nch, nMom, nPol)         !W-N
      real a21(km,nch, nMom, nPol)         !E-S
      real a22(km,nch, nMom, nPol)         !E-N
      real a00(km,nch, nMom, nPol)         !temp

      integer i1, i2

!                         ------

     if ( nobs .eq. 0 ) return    ! nothing to do (keep this)

     m_dlon = float(im) / 360.
     m_dlat = float(jm-1) / 180.
     
!    Loop over observations
!    ----------------------
     do nob = 1, nobs

!       Longitude
!       ---------
        o_lon = 1. + (lon(nob)+180.) * m_dlon
        i   = min(im, int( o_lon ))
        alfa  = o_lon - i
        if(i .eq. im) then
           i1 = im
           i2 = 1
        else
           i1 = i
           i2 = i + 1
        endif
        
!       Latitude
!       --------
        o_lat = 1. + (lat(nob) + 90.) * m_dlat
        j   = min( jm-1, int( o_lat ) )
        beta  = o_lat - j
        
        a11 = gField(i1,j,  :,:, :, :)
        a21 = gField(i2,j,  :,:, :, :)
        a12 = gField(i1,j+1,:,:, :, :)
        a22 = gField(i2,j+1,:,:, :, :)
        a00 = a11 + alfa * ( a21 - a11 )

        oField(:,:,nob,:,:) = a00 + beta * ( a12 + alfa * ( a22 - a12 )  -  a00 )
     
  end do ! loop over obs
  
  return
  
end subroutine Interpxy_bis

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Interp_--- Interpolates 3D field to observation locations
!
! !INTERFACE:
!
  subroutine Interpxy_aer ( lon, lat, nobs, im, jm, km, nq, nch, gField, &
                         oField )

! !USES:
!
      Implicit NONE

! !INPUT PARAMETERS:
!

      integer, intent(in)        :: nobs      ! Number of observations
      real,    intent(in)        :: lon(nobs) ! longitude in degrees [-180,+180]
      real,    intent(in)        :: lat(nobs) ! latitude  in degrees [-90,90]

      integer, intent(in)        :: im        ! zonal dimension
      integer, intent(in)        :: jm        ! meridional dimension
      integer, intent(in)        :: km        ! vertical dimension: 
                                              ! = 1 for 2D fields
                                              ! = km for mid-layer fields
                                              ! = km+1 for edge fields
      integer, intent(in)        :: nch        ! number of channels
      integer, intent(in)        :: nq        ! number of tracers

                                              ! Gridded Field
      real,    intent(in)        :: gField(nq,im,jm,km,nch) 

! !OUTPUT PARAMETERS:
!
                                              ! Interpolated profile
      real,    intent(out)       :: oField(nq,km,nch,nobs)

! !DESCRIPTION: This routine interpolates gridded model fields to observation
!               locations. This routine implements only the horizontal
!  interpolation.
!
!  IMPORTANT:   The input lon coordinates must be in [-180, 180].
!               The input field cannot ahve any UNDEFs.
!
! !SEE ALSO:
!
!              Module m_insitu which uses the same linear interpolation algorithm.
!
!
! !REVISION HISTORY:
!
!  10feb2010  da Silva  Simplified m_interp routine for profile interpolation.
!
!EOP
!-------------------------------------------------------------------------
 
    character(len=*), parameter :: myname = 'Interpxy_'


! Local
      integer i, j, k, nob
      real    o_lon, o_lat
      real    m_dlon, m_dlat
      real    alfa, beta

      real a11(nq,km,nch)         !W-S
      real a12(nq,km,nch)         !W-N
      real a21(nq,km,nch)         !E-S
      real a22(nq,km,nch)         !E-N
      real a00(nq,km,nch)         !temp

      integer i1, i2

!                         ------

     if ( nobs .eq. 0 ) return    ! nothing to do (keep this)

     m_dlon = float(im) / 360.
     m_dlat = float(jm-1) / 180.
     
!    Loop over observations
!    ----------------------
     do nob = 1, nobs

!       Longitude
!       ---------
        o_lon = 1. + (lon(nob)+180.) * m_dlon
        i   = min(im, int( o_lon ))
        alfa  = o_lon - i
        if(i .eq. im) then
           i1 = im
           i2 = 1
        else
           i1 = i
           i2 = i + 1
        endif
        
!       Latitude
!       --------
        o_lat = 1. + (lat(nob) + 90.) * m_dlat
        j   = min( jm-1, int( o_lat ) )
        beta  = o_lat - j
        
        a11 = gField(:,i1,j,  :,:)
        a21 = gField(:,i2,j,  :,:)
        a12 = gField(:,i1,j+1,:,:)
        a22 = gField(:,i2,j+1,:,:)
        a00 = a11 + alfa * ( a21 - a11 )

        oField(:,:,:,nob) = a00 + beta * ( a12 + alfa * ( a22 - a12 )  -  a00 )
     
  end do ! loop over obs
  
  return
  
end subroutine Interpxy_aer




