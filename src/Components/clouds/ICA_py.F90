! ICA_py.F90
! generate ICA subcolumns of a grid column, etc.
! pre-2016  pmn  original MCS version
! 020816    pmn  horizontal dims collapsed to 1D in simRe and simTau for SCS
!                -- now mcs.py will need to ravel horizontal dims explicitly

!...........................................................................
!BOP

! !ROUTINE: genICA --- Generate sub-columns of water variables 

! !INTERFACE:

  subroutine genICA ( ncols, km, &
                      pTop, delp, T, qv, ql, qi, f, &
                      pLim, Lp, Ls, &
                      qv_g, ql_g, qi_g, rc)

    use mod_utils
    use mod_consts
    use mod_StaticGridcolumnStatisticalModel

    implicit none

! !INPUT PARAMETERS:

    integer, intent(in) :: ncols       ! no. subcols to generate
    integer, intent(in) :: km          ! no. vertical layers

    ! Gridcolumn profiles:
    real*4,  intent(in) :: pTop        ! top pressure [Pa]
    real*4,  intent(in) :: delp (km)   ! layer pressure thickness [Pa]
    real*4,  intent(in) :: T    (km)   ! layer temperature [K]
    real*4,  intent(in) :: qv   (km)   ! layer specific humidity  [kg/kg]
    real*4,  intent(in) :: ql   (km)   ! layer specific liq water [kg/kg]
    real*4,  intent(in) :: qi   (km)   ! layer specific ice water [kg/kg]
    real*4,  intent(in) :: f    (km)   ! layer cloud fraction [0-1]

    real*4,  intent(in) :: pLim        ! largest pressure [Pa] at which
                                       !   to force zero variability

    real*4,  intent(in) :: Lp          ! Vertical decorrelation "length"
                                       !   for correlation matrix [Pa]
                                       !   typical: 100.e2 Pa

    real*4,  intent(in) :: Ls          ! Decorr length for Riishojgaard's
                                       !   flow-dependent correlation []
                                       !   (turned off unless > 0)
                                       !   typical: 0.1
! !OUTPUT PARAMETERS:

    real*4,  intent(out) :: qv_g (ncols,km)  ! specific humidity  [kg/kg]
    real*4,  intent(out) :: ql_g (ncols,km)  ! specific liq water [kg/kg]
    real*4,  intent(out) :: qi_g (ncols,km)  ! specific ice water [kg/kg]

    integer, intent(out) :: rc               ! error code

! !DESCRIPTION:
!
!   Generates subcolumns of
!   \bi
!     \item qv (specific humidity),
!     \item qc (specific cloud liq water content),
!     \item qi (specific cloud ice water content).
!   \ei
!   The *\_g  are the output subcolumns. 
!
!   Notes:
!   \bn
!     \item for vapor, liquid, and ice *paths* [kg/m2],
!        just multiply by rho * dz ~ delp / grav.
!    \item for *total condensate*, use
!        qc_g  = ql_g  + qi_g 
!   \en

!EOP
!...........................................................................

    ! parameters (potential input arguments):
    ! maximum allowed total saturation ratios in different phases.
    ! For mixed phase regions, the Smax will scale between the two
    ! based on the temperature.
    real*4, parameter :: Smax_liq = 1.10, Smax_ice = 1.40

    ! locals
    integer :: k
    real*4, dimension(km) :: w_v, w_l, w_i, w_t, work
    type (GridcolumnStatisticalModel) :: gsm
    type (MoistureSubcolumnEnsemble) :: mse

! -----

    ! PC - always start with no errors
    rc = 0
    exception = .false.

    ! prepare moisture contents
    ! first convert 'specific' quantities (mass per mass air)
    ! to mixing ratios (mass per mass dry air)
    work = 1. - (qv + ql + qi)
    w_v = qv / work
    w_l = ql / work
    w_i = qi / work

    ! now convert mixing ratios to contents (per Norris et al. 2008)
    ! so that finally w_x = rho_x / rho*, where rho* is the "virtual
    ! density", such that rho* = p / (Rd T) = rhod + rhov / eps. Then
    ! w_sat has the particularly simple but exact form 
    !   w_sat = eps * esat(T) / p.
    work = 1. + w_v / eps
    w_v = w_v / work
    w_l = w_l / work
    w_i = w_i / work

    ! total water content
    w_t = w_v + w_i + w_l

    ! validate f
    if (any(f < 0. .or. f > 1.)) then
      call myError('f out of [0,1]')
      rc = 1
      return
    end if

    ! setup GSM
    call GridcolumnStatisticalModel_create (gsm, &
      pTop, pLim, Smax_liq, Smax_ice, &
      km, delp, T, w_t, w_v, f, &
      Lp, Ls)
    if (exception) then
      rc = 2
      return
    end if
        
    ! generate subcolumns
    call GridcolumnStatisticalModel_generate (gsm, ncols, mse)
    if (exception) then
      rc = 3
      return
    end if
        
    ! load outputs
    if (gsm%klm1 .ge. 1) then
      qv_g (:,1:gsm%klm1) = spread(gsm%top_svap,                     1,ncols)
      ql_g (:,1:gsm%klm1) = spread(gsm%top_scon*gsm%fliq(1:gsm%klm1),1,ncols)
      qi_g (:,1:gsm%klm1) = spread(gsm%top_scon*gsm%fice(1:gsm%klm1),1,ncols)
    end if
    do k = 1, gsm%nkv
      qv_g (:,gsm%klm1+k) =  mse%svap (:,k)
      ql_g (:,gsm%klm1+k) =  mse%scon (:,k) * gsm%fliq(gsm%klm1+k)
      qi_g (:,gsm%klm1+k) =  mse%scon (:,k) * gsm%fice(gsm%klm1+k)
    end do

  end subroutine genICA


!...........................................................................
!BOP

! !ROUTINE: genICATauRe --- Generate sub-columns of Tau and Re

! !INTERFACE:

  subroutine genICATauRe ( ncols, km, &
                           pTop, delp, T, qv, ql, qi, f, &
                           pref, u, qils, qian, &
                           pLim, Lp, Ls, &
                           taul, taui, rel, rei, rc)

    use mod_utils
    use mod_consts
    use mod_StaticGridcolumnStatisticalModel
    use mod_simTauRe

    implicit none


! !INPUT PARAMETERS:

    integer, intent(in) :: ncols        ! no. subcols to generate
    integer, intent(in) :: km           ! no. vertical layers

    ! Gridcolumn profiles
    real*4,  intent(in) :: pTop         ! top pressure [Pa]
    real*4,  intent(in) :: delp (km)    ! layer pressure thickness [Pa]
    real*4,  intent(in) :: T    (km)    ! layer temperature [K]
    real*4,  intent(in) :: qv   (km)    ! layer specific humidity  [kg/kg]
    real*4,  intent(in) :: ql   (km)    ! layer specific liq water [kg/kg]
    real*4,  intent(in) :: qi   (km)    ! layer specific ice water [kg/kg]
    real*4,  intent(in) :: f    (km)    ! layer cloud fraction [0-1]

    ! These 4 for effective radius only:
    real*4,  intent(in) :: pref (km+1)  ! reference edge pressures [Pa]
    real*4,  intent(in) :: u    (km)    ! layer zonal (east) velocity [m/s]
    real*4,  intent(in) :: qils (km)    ! largescale ice wtr content [kg/kg]
    real*4,  intent(in) :: qian (km)    ! anvil      ice wtr content [kg/kg]

    real*4,  intent(in) :: pLim         ! largest pressure [Pa] at which
                                        !   to force zero variability

    real*4,  intent(in) :: Lp           ! Vertical decorrelation "length"
                                        !   for correlation matrix [Pa]
                                        !   typical: 100.e2 Pa

    real*4,  intent(in) :: Ls           ! Decorr length for Riishojgaard's
                                        !   flow-dependent correlation []
                                        !   (turned off unless > 0)
                                        !   typical: 0.1
! !OUTPUT PARAMETERS:

    real*4,  intent(out) :: taul (ncols,km)  ! liquid opt thicknesses []
    real*4,  intent(out) :: taui (ncols,km)  ! ice    opt thicknesses []
    real*4,  intent(out) :: rel  (ncols,km)  ! liquid effective radii [m]
    real*4,  intent(out) :: rei  (ncols,km)  ! ice    effective radii [m]

    integer, intent(out) :: rc               ! error code

! !DESCRIPTION:
!
!   Generates subcolumns of
!   \bi
!     \item tau (cloud optical thickness),
!     \item re  (effective radius)
!   \ei
!   for liquid and ice phases.

!EOP
!...........................................................................

    ! parameters (potential input arguments):
    ! maximum allowed total saturation ratios in different phases.
    ! For mixed phase regions, the Smax will scale between the two
    ! based on the temperature.
    real*4, parameter :: Smax_liq = 1.10, Smax_ice = 1.40

    ! locals
    real*4, dimension(km) :: w_v, w_l, w_i, w_t, work
    type (GridcolumnStatisticalModel) :: gsm
    type (simTauRe) :: simulator

! -----

    rc = 0

    ! prepare moisture contents
    ! first convert 'specific' quantities (mass per mass air)
    !   to mixing ratios (mass per mass dry air)
    work = 1. - (qv + ql + qi)
    w_v = qv / work
    w_l = ql / work
    w_i = qi / work

    ! now convert mixing ratios to contents (per Norris et al. 2008)
    ! so that finally w_x = rho_x / rho*, where rho* is the "virtual
    ! density", such that rho* = p / (Rd T) = rhod + rhov / eps. Then
    ! w_sat has the particularly simple but exact form 
    !   w_sat = eps * esat(T) / p.
    work = 1. + w_v / eps
    w_v = w_v / work
    w_l = w_l / work
    w_i = w_i / work

    ! total water content
    w_t = w_v + w_i + w_l

    ! validate f
    if (any(f < 0. .or. f > 1.)) then
       call myError('f out of [0,1]')
       rc = 1
       return
    end if

    ! setup GSM
    call GridcolumnStatisticalModel_create (gsm, &
      pTop, pLim, Smax_liq, Smax_ice, &
      km, delp, T, w_t, w_v, f, &
      Lp, Ls)
    if (exception) then
      rc = 2
      return
    end if
        
    ! setup tauRe simulator
    call simTauRe_create (simulator, &
      gsm, pref, u, qils, qian)
    if (exception) then
      rc = 3
      return
    end if
        
    ! genewrate tau & re
    call simTauRe_generate (simulator, ncols, &
      taul, taui, rel, rei)
    if (exception) then
      rc = 4
      return
    end if
        
  end subroutine genICATauRe

!...........................................................................
!BOP

! !ROUTINE: genICA_GEOS5like --- Generate sub-columns of water variables 

! !INTERFACE:

  subroutine genICA_GEOS5like ( ncols, km, &
                                pTop, delp, T, qv, ql, qi, f_in, pref, &
                                qv_g, ql_g, qi_g, rc)

    use mod_utils
    use mod_icefrac
    use MAPL_SatVaporMod
    implicit none

! !INPUT PARAMETERS:

    integer, intent(in) :: ncols       ! no. subcols to generate
    integer, intent(in) :: km          ! no. vertical layers

    ! Gridcolumn profiles:
    real*4,  intent(in) :: pTop        ! top pressure [Pa]
    real*4,  intent(in) :: delp (km)   ! layer pressure thickness [Pa]
    real*4,  intent(in) :: T    (km)   ! layer temperature [K]
    real*4,  intent(in) :: qv   (km)   ! layer specific humidity  [kg/kg]
    real*4,  intent(in) :: ql   (km)   ! layer specific liq water [kg/kg]
    real*4,  intent(in) :: qi   (km)   ! layer specific ice water [kg/kg]
    real*4,  intent(in) :: f_in (km)   ! layer cloud fraction [0-1]
    real*4,  intent(in) :: pref (km+1) ! *reference* edge press [Pa]

! !OUTPUT PARAMETERS:

    real*4,  intent(out) :: qv_g (ncols,km)  ! specific humidity  [kg/kg]
    real*4,  intent(out) :: ql_g (ncols,km)  ! specific liq water [kg/kg]
    real*4,  intent(out) :: qi_g (ncols,km)  ! specific ice water [kg/kg]

    integer, intent(out) :: rc               ! error code

! !DESCRIPTION:
!
!   Generates subcolumns of
!   \bi
!     \item qv (specific humidity),
!     \item qc (specific cloud liq water content),
!     \item qi (specific cloud ice water content).
!   \ei
!   The *\_g  are the output subcolumns. 
!   Uses Maximum-random overlap of homogeneous clouds
!     in low, middle, and high layers, ala GEOS-5.

!EOP
!...........................................................................

    ! parameters
    real*4, parameter :: &
      PRS_LOW_MID=700.e2, PRS_MID_HIGH=400.e2

    ! locals
    integer :: k, LCLDMH, LCLDLM, band, klo(3), khi(3)
    real*4 :: pTope, pBote
    real*4, dimension(km) :: &
      pMid, fice, fliq, qSat, f, qlC, qiC, qvC, qv0
    real*4 :: U(ncols,3) 

! -----

    rc = 0

    ! Evaluate pMid
    ! -------------
    pBote = pTop
    do k = 1, km
      pTope = pBote
      pBote = pTope + delp(k)
      pMid(k) = (pBote + pTope) / 2.
    end do

    ! Evaluate qSat
    ! -------------
    call icefrac (km, T, fice, fliq)
    qSat = fliq * MAPL_EQsat(T,pMid,OverIce=.false.) + &
           fice * MAPL_EQsat(T,pMid,OverIce=.true. )

    ! Determine model level separating high-middle and low-middle clouds
    ! ------------------------------------------------------------------

    ! validate band boundaries
    if (.not.(PRS_LOW_MID > PRS_MID_HIGH)) then
      call myError('NOT: PRS_LOW_MID > PRS_MID_HIGH')
      rc = 1
      return
    end if
    if (.not.(PRS_LOW_MID < pref(km))) then
      call myError('NOT: PRS_LOW_MID < pref(km)')
      rc = 2
      return
    end if

    ! pick LCLDMH and LCLDLM
    k = 1
    do while (pref(k) < PRS_MID_HIGH)
      k=k+1
    end do
    LCLDMH = k
    do while (pref(k) < PRS_LOW_MID)
      k=k+1
    end do
    LCLDLM = k

    ! choose k indicies for each band
    klo = [1,        LCLDMH,   LCLDLM]
    khi = [LCLDMH-1, LCLDLM-1, km    ]

    ! In-cloud & clear portion moisture variables
    ! -------------------------------------------
    ! avoid nearly clear or overcast
    f = f_in
    where (f < 0.01)
      ! force clear
      f = 0.
      qlC = 0.
      qiC = 0.
      qv0 = qv
      qvC = qv
    elsewhere (f > 0.99)
      ! force overcast
      f = 1.
      qlC = ql
      qiC = qi
      qvC = qv
      qv0 = qv
    elsewhere
      ! partially cloudy ...
      ! in-cloud condensates:
      qlC = ql/f
      qiC = qi/f
      ! in-cloud & clear portion vapor:
      ! qv = f * qvC + (1-f) * qv0
      where (qv > qSat)
        ! Too much vapor for qSat.
        ! Can't have clear vapor > in-cloud value.
        ! So, force qv0 = qvC, therefore = qv.
        qvC = qv
        qv0 = qv
      elsewhere (qv < f*qSat)
        ! qSat in-cloud would give qv0 < 0.
        ! So, force qv0 = 0.
        qv0 = 0.
        qvC = qv/f
      elsewhere
        qvC = qSat
        qv0 = (qv-f*qSat)/(1.-f)
      end where
    end where

    ! generate uniform random numbers
    call random_number(U)

    ! maximum overlap within each band
    ! bands independent of each other
    do band = 1,3
      do k = klo(band), khi(band)
        where (U(:,band) < f(k))
          ql_g(:,k) = qlC(k)
          qi_g(:,k) = qiC(k)
          qv_g(:,k) = qvC(k)
        elsewhere
          ql_g(:,k) = 0.
          qi_g(:,k) = 0.
          qv_g(:,k) = qv0(k)
        end where
      end do
    end do

  end subroutine genICA_GEOS5like

!...........................................................................
!BOP

! !ROUTINE: genICA_COSP --- Generate sub-columns of water variables 

! !INTERFACE:

  subroutine genICA_COSP ( ncols, km, &
                           pTop, delp, T, qv, ql, qi, f_in, &
                           qv_g, ql_g, qi_g, rc)

    use mod_utils
    use mod_icefrac
    use MAPL_SatVaporMod
    implicit none

! !INPUT PARAMETERS:

    integer, intent(in) :: ncols       ! no. subcols to generate
    integer, intent(in) :: km          ! no. vertical layers

    ! Gridcolumn profiles:
    real*4,  intent(in) :: pTop        ! top pressure [Pa]
    real*4,  intent(in) :: delp (km)   ! layer pressure thickness [Pa]
    real*4,  intent(in) :: T    (km)   ! layer temperature [K]
    real*4,  intent(in) :: qv   (km)   ! layer specific humidity  [kg/kg]
    real*4,  intent(in) :: ql   (km)   ! layer specific liq water [kg/kg]
    real*4,  intent(in) :: qi   (km)   ! layer specific ice water [kg/kg]
    real*4,  intent(in) :: f_in (km)   ! layer cloud fraction [0-1]

! !OUTPUT PARAMETERS:

    real*4,  intent(out) :: qv_g (ncols,km)  ! specific humidity  [kg/kg]
    real*4,  intent(out) :: ql_g (ncols,km)  ! specific liq water [kg/kg]
    real*4,  intent(out) :: qi_g (ncols,km)  ! specific ice water [kg/kg]

    integer, intent(out) :: rc               ! error code

! !DESCRIPTION:
!
!   Generates subcolumns of
!   \bi
!     \item qv (specific humidity),
!     \item qc (specific cloud liq water content),
!     \item qi (specific cloud ice water content).
!   \ei
!   The *\_g  are the output subcolumns. 
!   Uses COSP (i.e., SCOPS) overlap of homogeneous clouds.

!EOP
!...........................................................................

    ! locals
    integer :: k
    real*4 :: pTope, pBote
    real*4, dimension(km) :: &
      pMid, fice, fliq, qSat, f, qlC, qiC, qvC, qv0
    logical :: cloudy(ncols,km) 

! -----

    rc = 0

    ! Evaluate pMid
    ! -------------
    pBote = pTop
    do k = 1, km
      pTope = pBote
      pBote = pTope + delp(k)
      pMid(k) = (pBote + pTope) / 2.
    end do

    ! Evaluate qSat
    ! -------------
    call icefrac (km, T, fice, fliq)
    qSat = fliq * MAPL_EQsat(T,pMid,OverIce=.false.) + &
           fice * MAPL_EQsat(T,pMid,OverIce=.true. )

    ! In-cloud & clear portion moisture variables
    ! -------------------------------------------
    ! avoid nearly clear or overcast
    f = f_in
    where (f < 0.01)
      ! force clear
      f = 0.
      qlC = 0.
      qiC = 0.
      qv0 = qv
      qvC = qv
    elsewhere (f > 0.99)
      ! force overcast
      f = 1.
      qlC = ql
      qiC = qi
      qvC = qv
      qv0 = qv
    elsewhere
      ! partially cloudy ...
      ! in-cloud condensates:
      qlC = ql/f
      qiC = qi/f
      ! in-cloud & clear portion vapor:
      ! qv = f * qvC + (1-f) * qv0
      where (qv > qSat)
        ! Too much vapor for qSat.
        ! Can't have clear vapor > in-cloud value.
        ! So, force qv0 = qvC, therefore = qv.
        qvC = qv
        qv0 = qv
      elsewhere (qv < f*qSat)
        ! qSat in-cloud would give qv0 < 0.
        ! So, force qv0 = 0.
        qv0 = 0.
        qvC = qv/f
      elsewhere
        qvC = qSat
        qv0 = (qv-f*qSat)/(1.-f)
      end where
    end where

    ! generate SCOPS overlap
    call scops(ncols,km,f,cloudy)

    ! apply SCOPS overlap
    do k = 1,km
      where (cloudy(:,k))
        ql_g(:,k) = qlC(k)
        qi_g(:,k) = qiC(k)
        qv_g(:,k) = qvC(k)
      elsewhere
        ql_g(:,k) = 0.
        qi_g(:,k) = 0.
        qv_g(:,k) = qv0(k)
      end where
    end do

  end subroutine genICA_COSP

  subroutine scops(ncols,km,f,cloudy)

!   Simplified maximum-random SCOPS generator
!     - one grid-column version
!     - no convective cloud
!   PMN 05-03-2013 / 06-04-13

!   *****************************COPYRIGHT****************************
!   (c) British Crown Copyright 2009, the Met Office.
!   All rights reserved.
!  
!   Redistribution and use in source and binary forms, with or without
!   modification, are permitted provided that the
!   following conditions are met:
!  
!       * Redistributions of source code must retain the above
!         copyright  notice, this list of conditions and the following
!         disclaimer.
!       * Redistributions in binary form must reproduce the above
!         copyright notice, this list of conditions and the following
!         disclaimer in the documentation and/or other materials
!         provided with the distribution.
!       * Neither the name of the Met Office nor the names of its
!         contributors may be used to endorse or promote products
!         derived from this software without specific prior written
!         permission.
!  
!   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
!   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!  
!   *****************************COPYRIGHT*******************************
!   *****************************COPYRIGHT*******************************
!   *****************************COPYRIGHT*******************************

    integer, intent(in) :: ncols       ! no. subcols to generate
    integer, intent(in) :: km          ! no. vertical layers
    real*4,  intent(in) :: f (km)      ! layer cloud fraction [0-1]

    logical, intent(out) :: cloudy (ncols,km)

    ! locals
    integer :: k
    real*4 :: minf, U(ncols), rank(ncols)

    ! default clear
    cloudy = .false.

    ! top to bottom
    do k = 1, km

      if (k .eq. 1) then

        ! randomized top layer
        call random_number(rank)

      else

        ! maximally overlap the common cloudiness
        ! between the current and previous level,
        ! since can only "overlap" cloud that exists
        ! at both levels. Otherwise randomly overlap.
        call random_number(U)
        minf = min(f(k-1),f(k))
        where (rank .ge. minf) &
          rank = minf + (1. - minf) * U

      end if

      ! set cloudy mask
      where (rank .lt. f(k)) cloudy(:,k) = .true.

    end do

  end subroutine scops

!...........................................................................
!BOP

! !ROUTINE: simRe --- Simulate effective radius on a set of pixels

! !INTERFACE:

  subroutine simRe ( np, km, &
                     pTop, delp, T, &
                     pref, u, qils, qian, &
                     rel, rei)

    use mod_reff, only : get_tempor, reff

    implicit none


! !INPUT PARAMETERS:

    integer, intent(in) :: np             ! no. of pixels
    integer, intent(in) :: km             ! no. vertical layers

    ! subcolumn profiles
    real*4,  intent(in) :: pTop           ! top pressure [Pa]
    real*4,  intent(in) :: delp (np,km)   ! pressure thickness [Pa]
    real*4,  intent(in) :: T    (np,km)   ! temperature [K]
    real*4,  intent(in) :: pref (km+1)    ! *reference* edge press [Pa]
    real*4,  intent(in) :: u    (np,km)   ! zonal (east) velocity [m/s]
    real*4,  intent(in) :: qils (np,km)   ! largescale ice water [kg/kg]
    real*4,  intent(in) :: qian (np,km)   ! anvil      ice water [kg/kg]

! !OUTPUT PARAMETERS:

    real*4, intent(out) :: rel  (np,km)   ! liquid effective radii [m]
    real*4, intent(out) :: rei  (np,km)   ! ice    effective radii [m]

!EOP
!...........................................................................

    ! locals
    integer :: n, k
    real*4 :: tempor, pBote, pTope
    real*4 :: pMid(km)

! -----

    ! loop over pixels
    do n = 1, np

      ! midpoint pressure
      pTope = pTop
      do k = 1, km
        pBote = pTope + delp(n,k)
        pMid(k) = (pBote + pTope) / 2.
        pTope = pBote
      end do

      ! effective radii
      call get_tempor (km, pref, u(n,:), tempor)
      do k = 1, km
        call reff (T(n,k), pMid(k)/100., &
          qils(n,k), qian(n,k), tempor, &
          rel(n,k), rei(n,k))
      end do

    end do

  end subroutine simRe


!...........................................................................
!BOP

! !ROUTINE: simTau --- Simulate cloud optical depth on a set of pixels

! !INTERFACE:

  subroutine simTau ( np, km, &
                      delp, rel, rei, ql, qi, &
                      taul, taui)

    use mod_consts, only : grav
    use mod_simTauRe, only : getftauVis

    implicit none

! !INPUT PARAMETERS:

    integer, intent(in) :: np             ! no. of pixels
    integer, intent(in) :: km             ! no. vertical layers

    real*4,  intent(in) :: delp (np,km)   ! pressure thickness [Pa]
    real*4,  intent(in) :: rel  (np,km)   ! liquid effective radii [m]
    real*4,  intent(in) :: rei  (np,km)   ! ice    effective radii [m]
    real*4,  intent(in) :: ql   (np,km)   ! spec. liquid water [kg/kg]
    real*4,  intent(in) :: qi   (np,km)   ! spec. ice    water [kg/kg]

! !OUTPUT PARAMETERS:

    real*4, intent(out) :: taul (np,km)   ! liquid opt thicknesses []
    real*4, intent(out) :: taui (np,km)   ! ice    opt thicknesses []

!EOP
!...........................................................................

    ! locals
    integer :: n
    real*4 :: delpog(km), ftauLqdVis(km), ftauIceVis(km)

! -----

    ! loop over pixels
    do n = 1, np

      ! Visible COT per condensed water path (see mod_simTauRe)
      call getftauVis (rel(n,:), rei(n,:), km, ftauLqdVis, ftauIceVis)

      ! optical depths per layer
      delpog = delp(n,:) / grav
      taul(n,:) = ql(n,:) * delpog * ftaulqdVis
      taui(n,:) = qi(n,:) * delpog * ftauiceVis

    end do

  end subroutine simTau


!...........................................................................
!BOP

! !ROUTINE: qSatMCS --- Saturation specific humidity (MCS method)

! !INTERFACE:

  subroutine qSatMCS (km, pTop, delp, T, qSat)

    use mod_icefrac
    use MAPL_SatVaporMod
    implicit none

! !INPUT PARAMETERS:

    integer, intent(in) :: km          ! no. vertical layers

    ! Gridcolumn profiles:
    real*4,  intent(in) :: pTop        ! top pressure [Pa]
    real*4,  intent(in) :: delp (km)   ! layer pressure thickness [Pa]
    real*4,  intent(in) :: T    (km)   ! layer temperature [K]

! !OUTPUT PARAMETERS:

    real*4,  intent(out) :: qSat (km)  ! sat. spec. humidity [kg/kg]

! !DESCRIPTION:
!
!   Saturation specific humidity [kg/kg] by MCS method: 
!     an icefrac weighting of MAPL Qsat values over water and ice.

!EOP
!...........................................................................

    ! locals
    integer :: k
    real*4 :: pTope, pBote
    real*4, dimension(km) :: pMid, fice, fliq

! -----

    ! Evaluate pMid
    ! -------------
    pBote = pTop
    do k = 1, km
      pTope = pBote
      pBote = pTope + delp(k)
      pMid(k) = (pBote + pTope) / 2.
    end do

    ! Evaluate qSat
    ! -------------
    call icefrac (km, T, fice, fliq)
    qSat = fliq * MAPL_EQsat(T,pMid,OverIce=.false.) + &
           fice * MAPL_EQsat(T,pMid,OverIce=.true. )

  end subroutine qSatMCS

!...........................................................................
