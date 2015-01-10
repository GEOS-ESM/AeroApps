!
!  Simple f77 wrapper for the Python interface to the Mie Calculator.
!  This version works from columns already interpolated to obs location.

Module AOP_Mod

  use Chem_MieMod

  type(Chem_Mie) :: mieTables

  integer :: intialized = 0

end Module AOP_Mod

!.............................................................................

  subroutine AOPcreate(rcfile,rc)
  use AOP_Mod
  implicit NONE

!
!  Load the Mie tables given rc file.
!

  character(len=*), intent(in)  :: rcfile           ! resource file, e.g., Aod_EOS.rc
  integer,          intent(out) :: rc

! Create the Mie Tables
! ---------------------
  mieTables = Chem_MieCreate(rcfile,rc)
  if ( rc /= 0 ) then
     print *, 'Cannot create Mie tables from '//trim(rcfile)
     return
  end if

  initialized = 1

  return

end subroutine AOPcreate

subroutine AOPdestroy(rc)
   use AOP_Mod
   use Chem_MieMod
   implicit NONE
   integer, intent(out) :: rc
   call Chem_MieDestroy(rc)
   if ( rc /= 0 ) then
     print *, 'Cannot destroy MieTables'
     return
  end if
  initialized = 0
end subroutine AOPdestroy

!.............................................................................

 subroutine AOPget ( which, km, nobs, nch, nq, channels, vname, qm, rh, verbose &
                     aop, rc )

! Returns aod, ssa and asymmetry factor profiles.

  use Chem_MieMod
  use AOP_Mod
  implicit NONE

  character(len=*), intent(in)  :: which            ! which AOP to compute

  integer,          intent(in)  :: km               ! number vertical layers
  integer,          intent(in)  :: nobs             ! number of profiles

  integer,          intent(in)  :: nch              ! number of channels
  real,             intent(in)  :: channels(nch)

  integer,          intent(in)  :: nq               ! number of tracers
  character,        intent(in)  :: vname(nq,16)     ! variable name

  real,             intent(in)  :: qm(km,nq,nobs)   ! (mixing ratio) * delp/g
  real,             intent(in)  :: rh(km,nobs)      ! relative humidity

  integer,          intent(in)  :: verbose

  real,             intent(out) :: aop(km,nch,nobs) ! aerosol optical property
 
  integer,          intent(out) :: rc

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
        print *, 'Mie resource files does not specify the required channel'
        print *, 'Channels requested:  ', channels
        print *, 'Channels on RC file: ', 1.e+9 * mieTables%channels
        rc = 99
        return
   end if

! Initialize output arrays to zero
! --------------------------------
  aop = 0.0

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

!    -----------------------------------------------------------------------------
!    Design consideration: 'which' should only containvariables that are additive,
!    therefore, no SSA but absorption optical depth is OK. The required normalization
!    should take place outside, in python.
!    
!    Loop over nobs, km, nch
!    --------------------------
     do i = 1, nobs
        do n = 1, nch
           do k =1, km
        
            if ( which == 'tau' ) then
                 call Chem_MieQuery(mieTables, idxTable, idxChannel(n), &
                                    qm(k,iq,i), rh(k,i), tau=aop(k,n,i) )
             else
               print *, 'AOPget: cannot handle <'//trim(what)//'>'
               rc = 1
            endif

            end do  ! end nch
         end do ! end km
      end do  ! end nobs
  end do ! end tracers

end subroutine AOPget
