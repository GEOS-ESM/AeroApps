module m_ods_obsdiags

! !REVISION HISTORY:
! 
!  22Jan2008  Todling  Initial code
!  13Feb2009  Todling  Add opt to return obs sensitivity only
!  21Jan2014  Todling  Add opt to control kind of sigo in file
!  
use m_chars,   only: lowercase
implicit none
private

public ods_obsdiags_setparam
public ods_obsdiags_getparam
public ods_obsdiags
public ods_obsdiags_init
public ods_obsdiags_clean

integer, save :: miter_        = 2
integer, save :: jiter_        = 1
logical, save :: initialized_  = .false.
logical, save :: lobsdiagsave_ = .false.
logical, save :: lobssens_     = .false.
logical, save :: ladjsigo_     = .false.
logical, save :: lreduced_     = .false.

interface ods_obsdiags_setparam
   module procedure setgsi_paramI_
   module procedure setgsi_paramL_
end interface

interface ods_obsdiags_getparam
   module procedure getgsi_paramI_
   module procedure getgsi_paramL_
end interface

interface ods_obsdiags
   module procedure ods_obsdiagsS_
   module procedure ods_obsdiagsV_
end interface

interface ods_obsdiags_init
   module procedure init_
end interface
interface ods_obsdiags_clean
   module procedure clean_
end interface


CONTAINS
      subroutine ods_obsdiagsS_ ( nlomx, tlomx, obimp, data, ioff, ii, ninfo, nobs, undef, passed )
      implicit none
      integer,intent(in)  :: ioff, ii, ninfo, nobs
      real,   intent(in)  :: undef
      real,   intent(out) :: nlomx, tlomx, obimp
      real omx
      real(4) data(ninfo,nobs)
      logical, optional, intent(out) :: passed
      logical  passed_
      integer miter, jiter, jj, idia
      real, allocatable :: muse(:)
     
        if (.not.initialized_ ) return
        if (.not.lobsdiagsave_) return

          idia = ioff
          passed_=.false.
          miter = miter_
          jiter = jiter_
          allocate(muse(miter+1))

          nlomx=undef
          tlomx=undef
          obimp=0.0
          omx  =0.0
          do jj=1,miter
            idia=idia+1
            muse(jj) = data(idia,ii)
          enddo
          do jj=1,miter+1
            idia=idia+1
            if(jj==jiter) then
               nlomx = data(idia,ii)
               if(muse(jj)>0.0) omx=data(idia,ii)
             endif
          enddo
          do jj=1,miter
            idia=idia+1
            if(jj==jiter) tlomx = data(idia,ii)
          enddo
          do jj=1,miter
            idia=idia+1
            if(jj==jiter.and.muse(jj)>0.0) then
               if (lobssens_) then
                 obimp =       data(idia,ii)
               else
                 obimp = omx * data(idia,ii)
               endif
               passed_=.true.
            endif
          enddo
          deallocate(muse)
          if(present(passed)) passed=passed_

      end subroutine ods_obsdiagsS_

      subroutine ods_obsdiagsV_ ( nlomx, tlomx, obimp, data, ioff, ii, ninfo, nobs, undef, passed )
      implicit none
      integer, intent(in)  :: ioff, ii, ninfo, nobs
      real,    intent(in)  :: undef
      real,    intent(out) :: nlomx(2), tlomx(2), obimp(2)
      logical, optional, intent(out) :: passed
      real(4) data(ninfo,nobs)
      real omx(2)
      logical  passed_
      integer miter, jiter, jj, idia
      real, allocatable :: muse(:)

     
        if (.not.initialized_ ) return
        if (.not.lobsdiagsave_) return

          idia = ioff
          passed_=.false.
          miter = miter_
          jiter = jiter_
          allocate(muse(miter+1))

          nlomx=undef
          tlomx=undef
          obimp=0.0
          omx=0.0
          do jj=1,miter
            idia=idia+1
            muse(jj) = data(idia,ii)
          enddo
          do jj=1,miter+1
            idia=idia+1
            if(jj==jiter) then
               nlomx(1) = data(idia,ii)
               if(muse(jj)>0.0) omx(1)=data(idia,ii)
            endif
            idia=idia+1
            if(jj==jiter) then
               nlomx(2) = data(idia,ii)
               if(muse(jj)>0.0) omx(2)=data(idia,ii)
             endif
          enddo
          do jj=1,miter
            idia=idia+1
            if(jj==jiter) tlomx(1) = data(idia,ii)
            idia=idia+1
            if(jj==jiter) tlomx(2) = data(idia,ii)
          enddo
          if (lobssens_) then
            do jj=1,miter
              idia=idia+1
              if(jj==jiter.and.muse(jj)>0.0) obimp(1) = data(idia,ii)
              idia=idia+1
              if(jj==jiter.and.muse(jj)>0.0) then
                 obimp(2) = data(idia,ii)
                 passed_=.true.
              endif
            enddo
          else
            do jj=1,miter
              idia=idia+1
              if(jj==jiter.and.muse(jj)>0.0) obimp(1) = omx(1) * data(idia,ii)
              idia=idia+1
              if(jj==jiter.and.muse(jj)>0.0) then
                 obimp(2) = omx(2) * data(idia,ii)
                 passed_=.true.
              endif
            enddo
          endif
          deallocate(muse)
          if(present(passed)) passed=passed_

      end subroutine ods_obsdiagsV_

      subroutine setgsi_paramI_ ( desc, value )
      implicit none
      character(len=*), intent(in) :: desc 
      integer, intent(in) :: value
      character(len=256) descript
      logical :: found
      found=.false.
      descript = lowercase(desc)
      if ( trim(descript) == 'miter' ) then
          miter_ = value
           found = .true.
      endif
      if ( trim(descript) == 'jiter' ) then
          jiter_ = value
           found = .true.
      endif
      if(.not.found)then
        print*, 'm_ods_obsdiags: option(I) ',trim(descript), ' cannot be set'
        print*, 'Aborting ...'
        call exit(1)
      endif
      end subroutine setgsi_paramI_

      subroutine setgsi_paramL_ ( desc, value )
      implicit none
      character(len=*), intent(in) :: desc 
      logical, intent(in) :: value
      character(len=256) descript
      logical :: found
      found=.false.
      descript = lowercase(desc)
      if ( trim(descript) == 'lobsdiagsave' ) then
           lobsdiagsave_ = value
           found = .true.
      endif
      if ( trim(descript) == 'lobssens'     ) then
           lobssens_     = value
           found = .true.
      endif
      if ( trim(descript) == 'ladjsigo'     ) then
           ladjsigo_     = value
           found = .true.
      endif
      if ( trim(descript) == 'reduce_diag'  ) then
           lreduced_    = value
           found = .true.
      endif
      if(.not.found)then
        print*, 'm_ods_obsdiags: option(L) ',trim(descript), ' cannot be set'
        print*, 'Aborting ...'
        call exit(1)
      endif
      end subroutine setgsi_paramL_

      subroutine getgsi_paramI_ ( desc, value )
      implicit none
      character(len=*), intent(in) :: desc 
      integer, intent(out) :: value
      character(len=256) descript
      logical :: found
      found=.false.
      descript = lowercase(desc)
      if ( trim(descript) == 'miter' ) then
           value = miter_
           found = .true.
      endif
      if ( trim(descript) == 'jiter' ) then
           value = jiter_
           found = .true.
      endif
      if(.not.found)then
        print*, 'm_ods_obsdiags: option(I) ',trim(descript), ' not found'
        print*, 'Aborting ...'
        call exit(1)
      endif
      end subroutine getgsi_paramI_

      subroutine getgsi_paramL_ ( desc, value )
      implicit none
      character(len=*), intent(in) :: desc 
      logical, intent(out) :: value
      character(len=256) descript
      logical :: found
      found=.false.
      descript = lowercase(desc)
      if ( trim(descript) == 'lobsdiagsave' ) then
           value = lobsdiagsave_
           found = .true.
      endif
      if ( trim(descript) == 'lobssens'     ) then
           value = lobssens_
           found = .true.
      endif
      if ( trim(descript) == 'ladjsigo'     ) then
           value = ladjsigo_
           found = .true.
      endif
      if ( trim(descript) == 'reduce_diag'  ) then
           value = lreduced_
           found = .true.
      endif
      if(.not.found)then
        print*, 'm_ods_obsdiags: option ',trim(descript), ' not found'
        print*, 'Aborting ...'
        call exit(1)
      endif
      end subroutine getgsi_paramL_

      subroutine init_ 
      implicit none
      initialized_  = .true.
      lobsdiagsave_ = .true.
      lobssens_     = .true.
      ladjsigo_     = .true.
      lreduced_     = .false.
      end subroutine init_

      subroutine clean_ 
      implicit none
      initialized_  = .false.
      lobsdiagsave_ = .false.
      lobssens_     = .false.
      ladjsigo_     = .false.
      lreduced_     = .false.
      end subroutine clean_
end module m_ods_obsdiags
