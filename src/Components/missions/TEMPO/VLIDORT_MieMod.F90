!
!  Simple f77 wrapper for the Python interface to the Mie Calculator.
!  This version works from columns already interpolated to obs location.
!  PC: modified from Mieobs_py.F90 - took out Mie table read.

!.............................................................................

!...................................................................................

 subroutine VLIDORT_getAOPscalar ( mieTables, km, nobs, nch, nq, channels, vname, verbose, &
                           qm, rh, &
                           tau, ssa, g, rc )

! Returns aod, ssa and asymmetry factor profiles.

  use Chem_MieMod
  implicit NONE

  type(Chem_Mie),   intent(inout)   :: mieTables
  integer,          intent(in)      :: km               ! number vertical layers
  integer,          intent(in)      :: nobs             ! number of profiles

  integer,          intent(in)      :: nch              ! number of channels
  real,             intent(in)      :: channels(nch)

  integer,          intent(in)      :: nq               ! number of tracers
  character,        intent(in)      :: vname(nq,16)     ! variable name

  integer,          intent(in)      :: verbose

  real,             intent(in)      :: qm(km,nq,nobs)   ! (mixing ratio) * delp/g
  real,             intent(in)      :: rh(km,nobs)      ! relative humidity

  real,             intent(out)     :: tau(km,nch,nobs) ! aerosol optical depth
  real,             intent(out)     :: ssa(km,nch,nobs) ! single scattering albedo
  real,             intent(out)     :: g(km,nch,nobs)   ! asymmetry factor
 
  integer,          intent(out)     :: rc

!                               ---

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

   end subroutine VLIDORT_getAOPscalar

!...................................................................................

 subroutine VLIDORT_getAOPvector ( mieTables, km, nobs, nch, nq, channels, vname, verbose, &
                           qm, rh,nMom,nPol, &
                           tau, ssa, g, pmom, rc )

! Returns aod, ssa, asymmetry factor, phase function profiles.

  use Chem_MieMod
  implicit NONE
  type(Chem_Mie),   intent(inout)   :: mieTables
  integer,          intent(in)      :: km               ! number vertical layers
  integer,          intent(in)      :: nobs             ! number of profiles

  integer,          intent(in)      :: nch              ! number of channels
  real,             intent(in)      :: channels(nch)

  integer,          intent(in)      :: nMom             ! number of legender momemts 
  integer,          intent(in)      :: nPol             ! number of components of the scattering matrix

  integer,          intent(in)      :: nq               ! number of tracers
  character,        intent(in)      :: vname(nq,16)     ! variable name

  integer,          intent(in)      :: verbose

  real,             intent(in)      :: qm(km,nq,nobs)   ! (mixing ratio) * delp/g
  real,             intent(in)      :: rh(km,nobs)      ! relative humidity

  real,             intent(out)     :: tau(km,nch,nobs) ! aerosol optical depth
  real,             intent(out)     :: ssa(km,nch,nobs) ! single scattering albedo
  real,             intent(out)     :: g(km,nch,nobs)   ! asymmetry factor
  real,             intent(out)     :: pmom(km,nch,nobs,nMom,nPol) ! elements of scattering phase matrix
 
  integer,          intent(out)     :: rc

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
          if (verbose == 1)  print*,'nPol mie table', mieTables%vtableUse%nPol
          spherical_ext = .FALSE.
     else if ( nPol == 6 .AND. mieTables%vtableUse%nPol == 4 ) then
          nPol_ = 4              ! file only has 4, user wants 6
          spherical_ext = .TRUE. ! special case, will set P22=P11, P44=P33
     else
          rc = 777
          print *, 'ERROR: inconsistent number of polarizations: ', nPol, mieTables%vtableUse%nPol 
          return
     end if
     if (verbose ==1) print *, 'npol_ test',iq, nPol_
 
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
                  pmom(k,n,i,imom,ipol) = pmom(k,n,i,imom,ipol) &
                                        + pmom_(imom,ipol) * ssa_ * tau_  
                  
               end do
            end do

!           Special handling, spherical symmetry
!           ------------------------------------
            if ( spherical_ext ) then
               pmom(k,n,i,:,iP22) = pmom(k,n,i,:,iP11)  
               pmom(k,n,i,:,iP44) = pmom(k,n,i,:,iP33)  
            end if

            end do  ! end nch
         end do ! end km
      end do  ! end nobs
  deallocate(pmom_)
  end do ! end tracers
  if (verbose == 1)  print*, 'end tracer', iq
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
                  pmom(k,n,i,imom,ipol)  = pmom(k,n,i,imom,ipol)  / ( ssa(k,n,i)  * tau(k,n,i)  ) 
                  end if
               end do
            end do   
         end do
      end do
   end do
   
   end subroutine VLIDORT_getAOPvector
