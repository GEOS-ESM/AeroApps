!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_odstransfer --- Implements transfers between ODS vectors
! 
! !INTERFACE:
!

    module  m_odstransfer

! !USES:
   use decompmodule, only : decomptype
   use ghostmodule, only  : ghosttype
   use m_zeit, only       : zeit_ci, zeit_co
   use parutilitiesmodule, only : parpatterntype
   implicit none

   PRIVATE
   PUBLIC  Transfertype   ! Transfer communication pattern

! !PUBLIC MEMBER FUNCTIONS:
!
   PUBLIC  Transfer_Init
   PUBLIC  Transfer_Clean
   PUBLIC  Transfer_Perform

!
! !DESCRIPTION:
! \label{MODS:Mod}
!  This module defines the observation data stream vector class.
!  It relies on the ODS library and HDF.
!
! !REVISION HISTORY: 
!  25Oct2002 Sawyer    Split off from m_ods for clarity
!
!EOP
!-------------------------------------------------------------------------

   Interface  Transfer_Init
      module procedure transfer_init_d_d
      module procedure transfer_init_d_g
      module procedure transfer_init_g_d
      module procedure transfer_init_g_g
      module procedure transferhalo_init
   end Interface

   Interface  Transfer_Clean
      module procedure transfer_clean
   end Interface

!  Transfer communication patterns
!  -------------------------------
   type  transfertype
      type(parpatterntype) pattern_r ! pattern for REAL8 arrays
      type(parpatterntype) pattern_i ! pattern for INT4 arrays
   end type transfertype


CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Transfer_init_d_d --- Create ODS vector transfer pattern
!
! !INTERFACE:
!
      subroutine Transfer_init_d_d( decomp_in,decomp_out,transfer,rc )

! !USES:
      use parutilitiesmodule, only : commglobal, INT4, REAL8,           &
                                     parpatterncreate, parpatternfree
      implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(decomptype), intent(in)     ::  decomp_in       ! In decomposition
  type(decomptype), intent(inout)  ::  decomp_out      ! Out decomposition
  type(transfertype), intent(inout)::  transfer    ! Transfer pattern
  integer,        intent(out)    ::  rc            ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - could not allocate mem

! !DESCRIPTION:
! \label{MODS:Transfer_init}
!    Builds the communication pattern for transfer between two
!    unghosted ODS vectors.  Note that the pattern is only one way;
!    to build the reverse direction this routine must be called
!    with the vectors in reverse order.
!
! !REVISION HISTORY:
!
!  17may2002   Sawyer   Initial code
!
!EOP
!-------------------------------------------------------------------------

      call parpatterncreate( commglobal, decomp_in, decomp_out,          &
                             transfer%pattern_r )
      call parpatterncreate( commglobal, decomp_in, decomp_out,          &
                             transfer%pattern_i, t=INT4 )
      rc = 0
      return

      end subroutine Transfer_init_d_d

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Transfer_init_g_d --- Create ODS vector transfer pattern
!
! !INTERFACE:
!
      subroutine Transfer_init_g_d( halo_in,decomp_out,transfer,rc )

! !USES:
      use parutilitiesmodule, only : commglobal, INT4, REAL8,           &
                                     parpatterncreate, parpatternfree
      implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(ghosttype), intent(in)      ::  halo_in         ! In decomposition
  type(decomptype), intent(inout)  ::  decomp_out      ! Out decomposition
  type(transfertype), intent(inout)::  transfer    ! Transfer pattern
  integer,        intent(out)    ::  rc            ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - could not allocate mem

! !DESCRIPTION:
! \label{MODS:Transfer_init}
!    Builds the communication pattern for transfer between one ghosted
!    and one unghosted ODS vector.  Note that the pattern is only one way;
!    to build the reverse direction this routine must be called
!    with the vectors in reverse order.
!
! !REVISION HISTORY:
!
!  17may2002   Sawyer   Initial code
!
!EOP
!-------------------------------------------------------------------------

      call parpatterncreate( commglobal, halo_in, decomp_out,          &
                             transfer%pattern_r )
      call parpatterncreate( commglobal, halo_in, decomp_out,          &
                             transfer%pattern_i, t=INT4 )
      rc = 0

      return

      end subroutine Transfer_init_g_d

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Transfer_init_d_g --- Create ODS vector transfer pattern
!
! !INTERFACE:
!
      subroutine Transfer_init_d_g( decomp_in,halo_out,transfer,rc )

! !USES:
      use parutilitiesmodule, only : commglobal, INT4, REAL8,           &
                                     parpatterncreate, parpatternfree
      implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(decomptype), intent(in)     ::  decomp_in       ! In decomposition
  type(ghosttype), intent(inout)   ::  halo_out        ! Out decomposition
  type(transfertype), intent(inout)::  transfer    ! Transfer pattern
  integer,        intent(out)    ::  rc            ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - could not allocate mem

! !DESCRIPTION:
! \label{MODS:Transfer_init}
!    Builds the communication pattern for transfer between one ghosted
!    and one unghosted ODS vector.  Note that the pattern is only one way;
!    to build the reverse direction this routine must be called
!    with the vectors in reverse order.
!
! !REVISION HISTORY:
!
!  17may2002   Sawyer   Initial code
!
!EOP
!-------------------------------------------------------------------------

      call parpatterncreate( commglobal, decomp_in, halo_out,            &
                             transfer%pattern_r )
      call parpatterncreate( commglobal, decomp_in, halo_out,            &
                             transfer%pattern_i, t=INT4 )
      rc = 0

      return

      end subroutine Transfer_init_d_g

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Transfer_init_g_g --- Create ODS vector transfer pattern
!
! !INTERFACE:
!
      subroutine Transfer_init_g_g( halo_in,halo_out,transfer,rc )

! !USES:
      use parutilitiesmodule, only : commglobal, INT4, REAL8,           &
                                     parpatterncreate, parpatternfree
      implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(ghosttype), intent(in)          ::  halo_in   ! In decomposition
  type(ghosttype), intent(inout)       ::  halo_out  ! Out decomposition
  type(transfertype), intent(inout)::  transfer  ! Transfer pattern
  integer,        intent(out)    ::  rc            ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - could not allocate mem

! !DESCRIPTION:
! \label{MODS:Transfer_init}
!    Builds the communication pattern for transfer between two
!    ghosted ODS vectors.  Note that the pattern is only one way;
!    to build the reverse direction this routine must be called
!    with the vectors in reverse order.
!
! !REVISION HISTORY:
!
!  17may2002   Sawyer   Initial code
!
!EOP
!-------------------------------------------------------------------------

      call parpatterncreate( commglobal, halo_in, halo_out,          &
                             transfer%pattern_r )
      call parpatterncreate( commglobal, halo_in, halo_out,          &
                             transfer%pattern_i, t=INT4 )
      rc = 0

      return

      end subroutine Transfer_init_g_g

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  TransferHalo_init --- Create ODS vector transfer pattern
!
! !INTERFACE:
!
      subroutine TransferHalo_Init ( halo_decomp, transfer, rc )

! !USES:
      use parutilitiesmodule, only : commglobal, INT4, REAL8,           &
                                     parpatterncreate
      use ghostmodule, only : ghostdefined
      implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(ghosttype), intent(in)    ::  halo_decomp    ! ODS descriptor
  type(transfertype), intent(inout)::  transfer ! Transfer pattern
  integer,        intent(out)    ::  rc             ! Error return code:
                                                    !  0 - all is well
                                                    !  1 - could not allocate 
                                                    !      memory

! !DESCRIPTION:
! \label{MODS:TransferHalo_init}
!    Builds the communication pattern for the update of ghost 
!    regions from the ods vector.  This is a special case of
!    {\tt ods\_transfer\_init}.  The pattern only contains that which
!    needs to be updated, not the entire contents.
!
! !REVISION HISTORY:
!
!  17Jul2002   Sawyer   Initial code
!
!EOP
!-------------------------------------------------------------------------
     if ( ghostdefined( halo_decomp ) ) then
        call parpatterncreate( commglobal, halo_decomp,transfer%pattern_r)
        call parpatterncreate( commglobal, halo_decomp,                    &
                            transfer%pattern_i, t=INT4 )
     endif
     rc = 0
     return

     end subroutine TransferHalo_Init

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Transfer_clean --- Destroys ODS vector transfer pattern
!
! !INTERFACE:
!
      subroutine Transfer_clean ( transfer, rc )

! !USES:
      use parutilitiesmodule, only : commglobal, parpatternfree
      implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(transfertype), intent(inout)::  transfer      ! Transfer pattern
  integer,        intent(out)    ::  rc            ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - could not allocate mem

! !DESCRIPTION:
! \label{MODS:Transfer_clean}
!    Frees the transfer pattern
!
! !REVISION HISTORY:
!
!  17may2002   Sawyer   Initial code
!
!EOP
!-------------------------------------------------------------------------

  call parpatternfree( commglobal, transfer%pattern_r )
  call parpatternfree( commglobal, transfer%pattern_i )
  rc = 0

  return

 end subroutine Transfer_clean


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Transfer_perform --- Perform the data communication
!
! !INTERFACE:
!
      subroutine Transfer_perform ( transfer, data_in, data_out, rc )

! !USES:
      use parutilitiesmodule, only : commglobal,                           &
                                     parbegintransfer, parendtransfer
      use m_MPodsdata, only : obs_vect
      implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(transfertype), intent(in)::  transfer      ! Transfer pattern
  type(obs_vect), intent(in)::      data_in       ! obs_vect in
  type(obs_vect), intent(inout)::   data_out      ! obs_vect out
  integer,        intent(out)::     rc           ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - problem

! !DESCRIPTION:
! \label{MODS:Transfer_perform}
!    Performs the transfer of the data only (use {\tt ods\_tranfer} 
!    for a full transfer of an ODS vector.  It is allowed that 
!    {\tt data\_in} = {\tt data\_out}.
!
! !REVISION HISTORY:
!
!  06dec2002   Sawyer   Initial code
!
!EOP
!-------------------------------------------------------------------------

  rc = 1

!
!   Distribute the individual obs vectors in ods
!   --------------------------------------------
  call parbegintransfer( transfer%pattern_r, data_in%lat, data_out%lat )
  call parendtransfer( transfer%pattern_r, data_in%lat, data_out%lat )

  call parbegintransfer( transfer%pattern_r, data_in%lon, data_out%lon )
  call parendtransfer( transfer%pattern_r, data_in%lon, data_out%lon )

  call parbegintransfer( transfer%pattern_r, data_in%lev, data_out%lev )
  call parendtransfer( transfer%pattern_r, data_in%lev, data_out%lev )

  call parbegintransfer( transfer%pattern_r, data_in%xm, data_out%xm )
  call parendtransfer( transfer%pattern_r, data_in%xm, data_out%xm )

  call parbegintransfer( transfer%pattern_r, data_in%obs, data_out%obs )
  call parendtransfer( transfer%pattern_r, data_in%obs, data_out%obs )

  call parbegintransfer( transfer%pattern_r, data_in%OmF, data_out%OmF )
  call parendtransfer( transfer%pattern_r, data_in%OmF, data_out%OmF )

  call parbegintransfer( transfer%pattern_r, data_in%OmA, data_out%OmA )
  call parendtransfer( transfer%pattern_r, data_in%OmA, data_out%OmA )

  call parbegintransfer( transfer%pattern_r, data_in%Xvec, data_out%Xvec )
  call parendtransfer( transfer%pattern_r, data_in%Xvec, data_out%Xvec )

!
! Transfer the integer quantities.  These require additions to PILGRIM
! --------------------------------------------------------------------

  call parbegintransfer( transfer%pattern_i, data_in%gid, data_out%gid )
  call parendtransfer( transfer%pattern_i, data_in%gid, data_out%gid )

  call parbegintransfer( transfer%pattern_i, data_in%kx, data_out%kx )
  call parendtransfer( transfer%pattern_i, data_in%kx, data_out%kx )

  call parbegintransfer( transfer%pattern_i, data_in%kt, data_out%kt )
  call parendtransfer( transfer%pattern_i, data_in%kt, data_out%kt )

  call parbegintransfer( transfer%pattern_i, data_in%ks, data_out%ks )
  call parendtransfer( transfer%pattern_i, data_in%ks, data_out%ks )

  call parbegintransfer( transfer%pattern_i, data_in%time, data_out%time )
  call parendtransfer( transfer%pattern_i, data_in%time, data_out%time )

  call parbegintransfer( transfer%pattern_i, data_in%qcexcl, data_out%qcexcl )
  call parendtransfer( transfer%pattern_i, data_in%qcexcl, data_out%qcexcl )

  call parbegintransfer( transfer%pattern_i, data_in%qchist, data_out%qchist )
  call parendtransfer( transfer%pattern_i, data_in%qchist, data_out%qchist )

  rc = 0

  return

 end subroutine Transfer_perform

end module m_odstransfer


