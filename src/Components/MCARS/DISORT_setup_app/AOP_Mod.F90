! AOP_Mod.F90
! Class to calculate aerosol properties

#include "ErrorHandle.h"

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

module AOP_Mod

  use ErrorHandleMod
  use Chem_MieMod
  implicit none

  type AOP
    type(Chem_Mie) :: mieTables
  end type AOP

contains

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function AOP_create (rcfile, rc) result (self)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    character(len=*), intent(in) :: rcfile  ! resource file, e.g., Aod_MODIS.rc
    integer, intent(out), optional :: rc  ! return code
    type (AOP) :: self

    __Iam__('AOP_create')

    ! Create the Mie Tables
    ! ---------------------
    self%mieTables = Chem_MieCreate(rcfile,status)
    print *, 'mieTables: ', trim(self%mieTables%bc_optics_file)
    if (status /= 0) then
       print *, 'Cannot create Mie tables from '//trim(rcfile)
       ASSERT_(.FALSE.)
    end if

    RETURN_(SUCCESS)

  end function AOP_create

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AOP_calculate (self, &
    km, nobs, wavelen_nm, nMom, nPol, nq, tracer, verbose, &
    qm, rh, tau, ssa, g, pmom, rc)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Returns aod, ssa, asymmetry factor, phase function profiles.

    implicit none
    type (AOP), intent(inout) :: self

    integer,          intent(in)  :: km               ! number vertical layers
    integer,          intent(in)  :: nobs             ! number of profiles

    real,             intent(in)  :: wavelen_nm       ! wavelength (nm)

    integer,          intent(in)  :: nMom             ! num of Legendre momemts
    integer,          intent(in)  :: nPol             ! num of comps of scattering matrix

    integer,          intent(in)  :: nq               ! number of tracers
    character(len=*), intent(in)  :: tracer(nq)       ! tracer names

    integer,          intent(in)  :: verbose

    real,             intent(in)  :: qm(km,nq,nobs)   ! (mixing ratio) * delp/g
    real,             intent(in)  :: rh(km,nobs)      ! relative humidity

    real,             intent(out) :: tau (km,nobs)           ! aerosol optical depth
    real,             intent(out) :: ssa (km,nobs)           ! single scattering albedo
    real,             intent(out) :: g   (km,nobs)           ! asymmetry factor
    real,             intent(out) :: pmom(km,nobs,nMom,nPol) ! elements of scattering phase matrix

    integer, intent(out), optional :: rc  ! return code

    ! locals
    real          :: idxChannel  ! should have been an int
    integer       :: idxTable

    integer       :: iq, m, i, k, iMom, iPol, nPol_
    real          :: tau_, ssa_, g_
    real, pointer :: pmom_(:,:)

    logical       :: spherical_ext ! whether to set P22=P11 and P44=P11
    integer,parameter :: iP11=1, iP12=2, iP33=3, iP34=4, iP22=5, iP44=6

    __Iam__('AOP_calculate')

    ! make sure mie table has enough moments
    if (nMom > self%mieTables%nMom) then
       print *, 'mieTables do not have enough moments', nMom, self%mieTables%nMom
       ASSERT_(.FALSE.)
    end if

    ! Determine channel index
    idxChannel = -1
    do m = 1, self%mieTables%nch
      if (abs(wavelen_nm - 1.e9*self%mieTables%channels(m)) < 1.) then
        idxChannel = m
        exit
      end if
    end do
    if (idxChannel < 0) then
      print *, 'Mie resource file does not set the required channel'
      print *, 'Channel requested:   ', wavelen_nm
      print *, 'Channels in RC file: ', 1.e+9 * self%mieTables%channels
      ASSERT_(.FALSE.)
    end if

    ! Initialize output arrays to zero
    tau  = 0.
    ssa  = 0.
    g    = 0.
    pmom = 0.

    ! Loop over aerosol species
    do iq = 1,nq

      idxTable = Chem_MieQueryIdx(self%mieTables,tracer(iq),status)
      if (idxTable == -1) cycle
      if (status /= 0) then
        print *, 'cannot get Mie index for '//tracer(iq)
        ASSERT_(.FALSE.)
      end if

      if (verbose.eq.1) &
        print *, '[+] Adding '//tracer(iq)//' contribution'

      ! Check number of moments on file for this species
      if (nMom > self%mieTables%vtableUse%nMom) then 
        print *, 'mieTables do not have enough moments', nMom, self%mieTables%vtableUse%nMom 
        ASSERT_(.FALSE.)
      end if

      ! Handle possible case of non-spherical dust, but spherical everything else
      if (nPol < self%mieTables%vtableUse%nPol) then
        ! user requests fewer polarizations
        nPol_ = nPol
        spherical_ext = .FALSE.
      else if (nPol == 6 .and. self%mieTables%vtableUse%nPol == 4) then
        ! file only has 4, user wants 6
        ! special case, will set P22=P11, P44=P33
        nPol_ = 4
        spherical_ext = .TRUE.
      else
        print *, 'inconsistent number of polarizations: ', nPol, self%mieTables%vtableUse%nPol 
        ASSERT_(.FALSE.)
      end if
 
      ! Loop over nobs, km
      allocate (pmom_(self%mieTables%nMom,nPol_),__STAT__)
      do i = 1, nobs
        do k =1, km
            
          call Chem_MieQuery(self%mieTables, idxTable, idxChannel, &
            qm(k,iq,i), rh(k,i), tau=tau_, ssa=ssa_, gasym=g_, pmom=pmom_)
  
          tau(k,i) = tau(k,i) + tau_
          ssa(k,i) = ssa(k,i) + tau_ * ssa_ 
            g(k,i) =   g(k,i) + tau_ * ssa_ * g_
          do iPol = 1, nPol_      
            do iMom = 1, nMom   
              pmom(k,i,iMom,iPol) = pmom(k,i,iMom,iPol) + pmom_(iMom,iPol) * ssa_ * tau_  
            end do
          end do

          ! Special handling, spherical symmetry
          if (spherical_ext) then
            pmom(k,i,:,iP22) = pmom(k,i,:,iP11)  
            pmom(k,i,:,iP44) = pmom(k,i,:,iP33)  
          end if

        end do ! end km
      end do  ! end nobs
      deallocate(pmom_,__STAT__)

    end do ! end tracers

    ! Normalize ssa and g
    where (    tau > 0.) ssa = ssa / tau
    where (ssa*tau > 0.)   g =   g / (ssa*tau) 

    ! Normalize pmom 
    do i = 1, nobs
      do k = 1, km
        if (ssa(k,i) * tau(k,i) > 0.) then
          pmom(k,i,:,:) = pmom(k,i,:,:) / (ssa(k,i) * tau(k,i)) 
        end if
      end do
    end do
   
    RETURN_(SUCCESS)
  end subroutine AOP_calculate

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AOP_destroy (self, rc)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    implicit none
    type (AOP), intent(inout) :: self
    integer, intent(out), optional :: rc
    __Iam__('AOP_destroy')
    call Chem_MieDestroy(self%mieTables, status)
    if (status /= 0) then
      print *, 'Cannot destroy MieTables'
      ASSERT_(.FALSE.)
    end if
    RETURN_(SUCCESS)
  end subroutine AOP_destroy

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AOP_PrepareFromMCS (tables, &
    cloud_file, spec_file, nq, tracer, isb, ise, ifb, ife, &
    nsc, nfc, nch, channels, km, rh, qm, rc)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    use netcdf
    use MAPL_ConstantsMod, only : MAPL_GRAV
    implicit none

    type (AOP), intent(inout) :: tables	     ! Mie tables
    character(len=*), intent(in) :: cloud_file, spec_file  ! MCS filename for ICA input

    integer, intent(in) :: nq                   ! num of tracers
    character(len=*), intent(in) :: tracer(nq)	! tracer names
    
    integer, intent(in) :: isb, ise  ! sample scans  isb:ise
    integer, intent(in) :: ifb, ife  ! sample frames ifb:ife

    integer, intent(out) :: nsc      ! number of scans  in isb:ise
    integer, intent(out) :: nfc      ! number of frames in ifb:ife

    integer, intent(out) :: nch                    ! number of channels in MCS
    real, intent(out), allocatable :: channels(:)  ! (nch) wavelength (nm)

    integer, intent(out) :: km                     ! number of layers in MCS
    real, intent(out), allocatable :: rh(:,:,:)	   ! (km,nfc,nsc)    relative humidity
    real, intent(out), allocatable :: qm(:,:,:,:)  ! (km,nq,nfc,nsc) (mix ratio)*delp/g

    integer, intent(out), optional :: rc

    real*4, parameter :: missing_r4 = 1.e15
    real,   parameter :: missing    = 1.e15

    character(len=nf90_max_name) :: name
    integer :: ncId, varId, dimId, dimIds(4)
    integer :: ns, nf, iq, n, k
    character(len=MAXSTRING) :: icaFile, speciesFile

    real, allocatable :: height(:)         ! layer midpoints (kilometers)
    real, allocatable :: delp(:,:,:)       ! pressure thickness [Pa]

    real*4, allocatable :: v1d(:), v2d(:,:), v3d(:,:,:)
    real*4 :: fillValue

    __Iam__('PhaseFunc')

    ! ~~~~~~~~~~~~~
    ! read ICA file
    ! ~~~~~~~~~~~~~

    ! basics
    write(icaFile,'(a)') cloud_file
    write(*,'(x,"reading: ",a)') trim(icaFile)
    NCCHECK( nf90_open(trim(icaFile), nf90_nowrite, ncId) )
    NCCHECK( nf90_inq_dimid(ncId, 'height', dimId) )
    NCCHECK( nf90_inquire_dimension(ncId, dimId, name, km) )
    print *, 'layers detected:', km
    NCCHECK( nf90_inq_varid(ncId, 'height', varId) )
    allocate(height(km),__STAT__)
    NCCHECK( nf90_get_var(ncId, varId, height) )
    print *, 'layer midpoint heights (kilometers):'
    print *, height
    NCCHECK( nf90_inq_dimid(ncId, 'height_edge', dimId) )
    NCCHECK( nf90_inquire_dimension(ncId, dimId, name, n) )
    if (n.ne.km+1) then
      print *, 'height_edge size inconsistency'
      ASSERT_(.FALSE.)
    end if
    NCCHECK( nf90_inq_dimid(ncId, 'nscans', dimId) )
    NCCHECK( nf90_inquire_dimension(ncId, dimId, name, ns) )
    NCCHECK( nf90_inq_dimid(ncId, 'nframes', dimId) )
    NCCHECK( nf90_inquire_dimension(ncId, dimId, name, nf) )
    write(*,'(x,i4," scans x ",i4," frames =", i8, " obs")') ns, nf, ns * nf

    ! validate sampling
    if (.not.(1.le.isb.and.isb.le.ise.and.ise.le.ns)) then
      print *, 'isb, ise, ns: ', isb, ise, ns
      print *, 'require 1 <= isb <= ise <= ns'
      ASSERT_(.FALSE.)
    end if
    nsc = ise - isb + 1
    if (.not.(1.le.ifb.and.ifb.le.ife.and.ife.le.nf)) then
      print *, 'ifb, ife, nf: ', ifb, ife, nf
      print *, 'require 1 <= ifb <= ife <= nf'
      ASSERT_(.FALSE.)
    end if
    nfc = ife - ifb + 1

    ! get RH 
    allocate(v3d(km,nfc,nsc),__STAT__)
    NCCHECK( nf90_inq_varid(ncId, 'RHMCS', varId) )
    NCCHECK( nf90_get_var(ncId, varId, v3d, start=[1,ifb,isb], count=[km,nfc,nsc]) )
    NCCHECK( nf90_get_att(ncId, varId, '_FillValue', fillValue) )
    allocate(rh(km,nfc,nsc),__STAT__)
    where (v3d.eq.fillValue)
      rh = missing
    elsewhere (v3d.le.0.)
      rh = 0.
    elsewhere
      rh = v3d
    end where
    deallocate(v3d,__STAT__)

    ! get delp
    allocate(v3d(km+1,nfc,nsc),__STAT__)
    NCCHECK( nf90_inq_varid(ncId, 'PE', varId) )
    NCCHECK( nf90_get_var(ncId, varId, v3d, start=[1,ifb,isb], count=[km+1,nfc,nsc]) )
    NCCHECK( nf90_get_att(ncId, varId, 'missing_value', fillValue) )
    where (v3d.eq.fillValue) v3d = missing_r4
    allocate(delp(km,nfc,nsc),__STAT__)
    do k = 1, km
      where ((v3d(k,:,:).ne.missing_r4).and.(v3d(k+1,:,:).ne.missing_r4))
        delp(k,:,:) = v3d(k+1,:,:) - v3d(k,:,:)
      elsewhere
        delp(k,:,:) = missing
      end where
    end do
    deallocate(v3d,__STAT__)
    NCCHECK( nf90_close(ncId) )

    ! ~~~~~~~~~~~~~~~~~
    ! read species file
    ! ~~~~~~~~~~~~~~~~~

    ! basics
    write(speciesFile,'(a)') spec_file
    write(*,'(x,"reading: ",a)') trim(speciesFile)
    NCCHECK( nf90_open(trim(speciesFile), nf90_nowrite, ncId) )
    NCCHECK( nf90_inq_dimid(ncId, 'height', dimId) )
    NCCHECK( nf90_inquire_dimension(ncId, dimId, name, n) )
    if (n.ne.km) then
      print *, 'height size inconsistency'
      ASSERT_(.FALSE.)
    end if
    NCCHECK( nf90_inq_varid(ncId, 'height', varId) )
    allocate(v1d(km),__STAT__)
    NCCHECK( nf90_get_var(ncId, varId, v1d) )
    if (any(real(v1d).ne.height)) then
      print *, 'height inconsistency'
      ASSERT_(.FALSE.)
    end if
    deallocate(v1d,__STAT__)
    deallocate(height,__STAT__)
    NCCHECK( nf90_inq_dimid(ncId, 'wavelength', dimId) )
    NCCHECK( nf90_inquire_dimension(ncId, dimId, name, nch) )
    print *, 'channels detected:', nch
    NCCHECK( nf90_inq_varid(ncId, 'wavelength', varId) )
    allocate(channels(nch),__STAT__)
    NCCHECK( nf90_get_var(ncId, varId, channels) )
    print *, 'wavelengths (nm):'
    print *, channels
    NCCHECK( nf90_inq_dimid(ncId, 'nscans', dimId) )
    NCCHECK( nf90_inquire_dimension(ncId, dimId, name, n) )
    if (n.ne.ns) then
      print *, 'number of scans inconsistency'
      ASSERT_(.FALSE.)
    end if
    NCCHECK( nf90_inq_dimid(ncId, 'nframes', dimId) )
    NCCHECK( nf90_inquire_dimension(ncId, dimId, name, n) )
    if (n.ne.nf) then
      print *, 'number of frames inconsistency'
      ASSERT_(.FALSE.)
    end if

    ! read each required tracer
    allocate(v3d(km,nfc,nsc),__STAT__)
    allocate(qm(km,nq,nfc,nsc),__STAT__)
    do iq = 1, nq
      print *, 'reading tracer ', trim(tracer(iq)), ' ...'
      NCCHECK( nf90_inq_varid(ncId, trim(tracer(iq)), varId) )
      NCCHECK( nf90_get_var(ncId, varId, v3d, start=[1,ifb,isb], count=[km,nfc,nsc]) )
      NCCHECK( nf90_get_att(ncId, varId, '_FillValue', fillValue) )
      do k = 1, km
        where (v3d(k,:,:).eq.fillValue.or.delp(k,:,:).eq.missing)
          qm(k,iq,:,:) = missing
        elsewhere (v3d(k,:,:).le.0.)
          qm(k,iq,:,:) = 0.
        elsewhere
          qm(k,iq,:,:) = v3d(k,:,:) * delp(k,:,:) / MAPL_GRAV
        end where
      end do
    end do
    deallocate(v3d,__STAT__)
    deallocate(delp,__STAT__)
    NCCHECK( nf90_close(ncId) )

    RETURN_(SUCCESS)

  end subroutine AOP_PrepareFromMCS

end module AOP_Mod
