!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_SparseComm_SHMEM
!
! !DESCRIPTION:
!
! !INTERFACE:


Module m_SparseComm_SHMEM_PTR
#ifdef USE_SHMEM_PTR

#include "assert.H"
  Use m_die, only : assert_
  implicit none
  Private ! Except for
  
  Public :: SparseComm              ! public datatype
  Public :: Analyze                 ! analyze communication pattern
  Public :: PartialReduceScatter    ! implements collective operation
  Public :: PartialGather           ! implements collective operation
  Public :: DescribeComm            ! writes a description of communication pattern
  Public :: Clean                   ! deallocates SparseComm objects
  Public :: ReleaseSymmetricMemory  ! release symmetric heap memory for complete clean-up

  Interface PartialGather
     Module Procedure pgath_r2
!!$     Module Procedure pgath_i1
  End Interface


! !REVISION HISTORY:
!       29Oct01 - Tom Clune
!               - Finally added this intro section

!EOP ___________________________________________________________________

  !!! supported operations
  Integer, Parameter, Public :: GATHER         = 1
  Integer, Parameter, Public :: REDUCE_SCATTER = 2

!!! this include file is necessary for some shmem functions
  Include 'mpp/shmem.fh'

  Character(Len=*), Parameter :: myname = 'm_SparseComm'

  Integer, Parameter :: NOOP                 = 0
  Integer, Parameter :: FETCH                = 1
  Integer, Parameter :: FETCH_AND_ACCUMULATE = 2
  Integer, Parameter :: PROXY_FORWARD        = -1 ! signal in analysis, but no communication
  Integer, Parameter :: FILL_PROXY_THEN_ACCUMULATE = 4
  Integer, Parameter :: CLEAR = 8
  Integer, Parameter :: LOG2_MAX_PES   = 10 ! 1024 = maximum number of pes
  Integer, Parameter :: MAX_PES        = 2 ** LOG2_MAX_PES

  Type PacketDesc
     Integer :: displ     = -1
     Integer :: count     = -1
     Integer :: pe        = -1 ! host for the data
  End Type PacketDesc

  Type SparseComm
     Private
     Integer :: comm_id   = -1 ! Unique communicator handle
     Integer :: operation = -1 ! Reduce scatter or Gather
     Integer :: offset_in = -1
     Integer :: offset_out = -1
     Integer :: npes
     Type (PacketDesc),    Pointer :: packets(:)   => Null()
     Type (PacketDesc),    Pointer :: prefill(:)   => Null()
  End Type SparseComm


  Pointer(p_sym_buffer_i, int_symmetric_buffer)   ! (pointer, pointee)
  Pointer(p_sym_buffer_r, real_symmetric_buffer)  ! (pointer, pointee)

  Save p_sym_buffer_i
  Save p_sym_buffer_r
  Real           :: real_symmetric_buffer(1) ! Allocate with shpalloc or similar
  Integer        :: int_symmetric_buffer(1)  ! Allocate with shpalloc or similar
  Integer, Save  :: symmetric_buffer_size = 0
  Integer, Save  :: global_comm_id        = 0

Contains ! Functions and subroutines

  Subroutine Analyze(comm, seg_in, seg_out, nav, offset, operation, mp_comm)
    Use m_Navigator, only : Navigator
    Use m_Navigator, only : lsize
    Use m_Navigator, only : get
    Use m_mpout, only : mpout_log, mpout
    Use m_die, only : die
    Implicit None
    Type (SparseComm), Intent(InOut) :: comm
    Logical, Intent(In) :: seg_in(:)  ! nseg x np
    Logical, Intent(In) :: seg_out(:)
    Type (Navigator), Intent(In) :: nav  ! contains sizes of the segments
    Integer, Intent(In) :: offset
    Integer, Intent(In) :: operation
    Integer, Intent(In) :: mp_comm
  !=========================================================================

    Character(Len=*), Parameter :: myname_ = myname//'::Analyze'
    
    Pointer(p_you_have, you_have)
    Logical :: you_have(Size(seg_in)) ! Assume same size as Integer
    Pointer(p_I_have, I_have)
    Logical :: I_have(Size(seg_in)) ! Assume same size as Integer

    Logical :: transmit

    Integer :: my_id, np
    Integer :: segsize
    Integer :: p, pp
    Integer :: ier
    Integer :: n_words
    Integer :: n_packets, i_packet, n_segments, iseg, iseg1, nseg, cnt

    Integer, External :: num_pes
    Integer, External :: my_pe

!----------------------------------------------------------------------

    ASSERT(size(seg_in) == size(seg_out))
    ASSERT(CleanQ(comm))
    ASSERT(size(seg_in) == lsize(nav))

    comm%operation = operation

 if(offset<1 .or. offset>lsize(nav)) then
    comm%offset_in =0
    comm%offset_out=0
 else
    Select Case (operation)
    Case (GATHER)
       Call get(nav, offset, displ = comm%offset_in)
       comm%offset_out = 0
    Case (REDUCE_SCATTER)
       comm%offset_in  = 0
       Call get(nav, offset, displ = comm%offset_out)
    Case Default
       call die(myname_,'Illegal operation must be GATHER or REDUCE_SCATTER.',operation)
    End Select
 endif
    
    n_segments = size(seg_in)
    np = num_pes()
    Comm%npes = np
    If (np == 1) Return ! skip analysis - trivial case

    my_id = my_pe()

! Store possession flag in a global arena for easy access
! -------------------------------------------------------
    n_words = n_segments ! required
    Call ReserveSymmetricMemory(n_words_int = n_words)

    p_I_have = p_sym_buffer_i
    I_have   = seg_in
    
! Prefill info
! ------------
    n_packets = Count_Contiguous(seg_in)
    Allocate(comm%prefill(n_packets), STAT = ier)
       ALWAYS_ASSERT(ier == 0)

    transmit = .false.
    i_packet = 0
    Do iseg = 1, n_segments
       If (seg_in(iseg)) Then
          If (transmit) Then ! old packet
             comm%prefill(i_packet)%count = comm%prefill(i_packet)%count + 1
          Else ! new packet
             i_packet = i_packet + 1
             comm%prefill(i_packet)%pe    = my_id
             comm%prefill(i_packet)%displ = iseg
             comm%prefill(i_packet)%count = 1
             transmit = .true.
          End If
       Else
          transmit = .false.
       End If
    End Do

   ALWAYS_ASSERT(i_packet == n_packets)

! Convert displacement and count units from segments to points
! ------------------------------------------------------------
   Do i_packet = 1, n_packets

      iseg = comm%prefill(i_packet)%displ
      nseg = comm%prefill(i_packet)%count

      Call Get(nav, iseg, displ = comm%prefill(i_packet)%displ)         
      
      cnt = 0
      Do iseg1 = iseg, iseg + comm%prefill(i_packet)%count - 1
         Call Get(nav, iseg1, ln = segsize)
         cnt = cnt + segsize
      End Do
      comm%prefill(i_packet)%count = cnt
   End Do


!-------------------------------------------------------------------------
    Call SHMEM_BARRIER_ALL()  ! We've all put the info in the table


! Loop over processors and see who has useful data.
! -------------------------------------------------
    n_packets = 0
    Do pp = 0, np - 1 ! stagger lookups
       p = mod(pp + my_id, np)

       p_you_have = SHMEM_PTR(I_have, p)
       n_packets = n_packets + Count_Contiguous(seg_out .and. you_have(1:n_segments))

    End Do

    Allocate(comm%packets(n_packets), STAT = ier)
    ALWAYS_ASSERT(ier == 0)

! Loop over processors and fill packet descriptions
! -------------------------------------------------

    i_packet = 0
    Do pp = 0, np - 1 ! stagger lookups
       p = mod(pp + my_id, np)

       p_you_have = SHMEM_PTR(I_have, p)
       
       transmit = .false.
       Do iseg = 1, n_segments
          If (seg_out(iseg) .and. you_have(iseg)) Then
             If (transmit) Then ! already in a packet
                comm%packets(i_packet)%count = comm%packets(i_packet)%count + 1
             Else ! new packet
                i_packet = i_packet + 1
                comm%packets(i_packet)%pe    = p
                comm%packets(i_packet)%displ = iseg
                comm%packets(i_packet)%count = 1
                transmit = .true.
             End If
          Else
             transmit = .false.
          End If
       End Do
   End Do

   ALWAYS_ASSERT(i_packet == n_packets)

! Wait until all PEs are done reading from shared buffer
   Call SHMEM_BARRIER_ALL()

! Convert displacement and count units from segments to points
! ------------------------------------------------------------
   Do i_packet = 1, n_packets

      iseg = comm%packets(i_packet)%displ
      nseg = comm%packets(i_packet)%count

      Call Get(nav, iseg, displ = comm%packets(i_packet)%displ)         
      
      cnt = 0
      Do iseg1 = iseg, iseg + comm%packets(i_packet)%count - 1
         Call Get(nav, iseg1, ln = segsize)
         cnt = cnt + segsize
      End Do
      comm%packets(i_packet)%count = cnt
   End Do

 Contains
   
   Integer Function Count_Contiguous(flags)
   Logical, Intent(In) :: flags(:)

   Integer :: cnt, iseg
   Logical :: transmit

   transmit = .false.
   cnt = 0
   Do iseg = 1, Size(flags)
       If (flags(iseg)) Then
          If (transmit) Then ! old packet
          Else ! new packet
             cnt = cnt + 1
             transmit = .true.
          End If
       Else
          transmit = .false.
       End If
    End Do

    Count_Contiguous = cnt

    End Function Count_Contiguous


  End Subroutine Analyze
  !===================================================================
  Subroutine Clean(comm)
    Implicit None
    Type (SparseComm), Intent(InOut) :: comm
  !===================================================================
    Character(Len=*), Parameter :: myname_ = myname//'::Clean'
    Integer :: ier
    Integer :: round

    If (comm%npes == 1) Then ! trivial case - nothing to do
       comm%npes = -1 
       Return
    End If

    If (initQ(comm)) Then
       Deallocate(comm%packets, comm%prefill, STAT = ier)
          ALWAYS_ASSERT(ier == 0)
    End If

  End Subroutine Clean

  Logical Function initQ(comm)
    Type (SparseComm), Intent(In) :: comm
    initQ = (comm%npes >= 0)
  End Function initQ

  !==================================================
  Subroutine PartialReduceScatter(x_in, x_out, comm)
    Use m_die, only : assert_
    Use m_mpout, only : mpout, mpout_log
    Implicit None
    Real,    Intent(In)  :: x_in(:,:)  
    Real,    Intent(Out) :: x_out(:,:)
    Type (SparseComm), Intent(In) :: comm
  !==================================================

    Character(Len=*), Parameter :: myname_ = myname//'::PartialReduceScatter'

    ASSERT(comm%operation == REDUCE_SCATTER)
    If (comm%npes == 1) Then ! trivial case
       x_out = x_in
       Return
    End If

    Call PartialComm(x_in, x_out, comm)

  End Subroutine PartialReduceScatter
  !==================================================
  Subroutine pgath_r2(x_in, x_out, comm)
    Use m_die, only : assert_
    Use m_mpout, only : mpout, mpout_log
    Implicit None
    Real,    Intent(In)  :: x_in(:,:)  
    Real,    Intent(Out) :: x_out(:,:)
    Type (SparseComm), Intent(In) :: comm
  !==================================================

    Character(Len=*), Parameter :: myname_ = myname//'::PartialGather'

    ASSERT(comm%operation == GATHER)
    If (comm%npes == 1) Then ! trivial case
       x_out = x_in
       Return
    End If

    Call PartialComm(x_in, x_out, comm)

  End Subroutine Pgath_r2



  !==============================================================
  Subroutine DescribeComm(comm, label)
    Use m_mpout, only : mpout_log, mpout
    Use m_die, only : die
    Type (SparseComm), Intent(In) :: comm
    Character(Len=*), Intent(In) :: label
  !==============================================================
  End Subroutine DescribeComm

  !==========================================================
  Subroutine ReleaseSymmetricMemory()
    Use m_die, only : perr
    Use m_mall, only : mall_ci, mall_co, mall_ison
    Implicit None
  !==========================================================
    Character(Len=*), Parameter :: myname_ = myname//'::ReleaseSymmetricMemory'
    External SHPDEALLC
    Integer :: ier
    Integer, Parameter :: abort = 0 ! do not abort - return error codes instead
  !!! Error return codes
    Integer, Parameter :: address_outside_bounds   = -3
    Integer, Parameter :: block_already_free       = -4
    Integer, Parameter :: address_not_at_beginning = -5

    
    If (mall_ison()) Call mall_co(symmetric_buffer_size,myname)
        
    Call SHPDEALLC(p_sym_buffer_r,ier, abort)

    Select Case (ier)
    Case (address_outside_bounds)
       Call perr(myname_, 'address p_sym_buffer is outside bounds of symmetric heap.')
    Case (block_already_free)
       Call perr(myname_, 'symmetric buffer is already released.')
    Case (address_not_at_beginning)
       Call perr(myname_, 'p_sym_buffer is not at the beginning of the block.')
    Case Default
       ! everything is a-ok
       symmetric_buffer_size = 0
       p_sym_buffer_i = p_sym_buffer_r
    End Select

  End Subroutine ReleaseSymmetricMemory ! static

  !==================================================
  Subroutine PartialComm_r2(x_in, x_out, comm)
    Use m_die, only : assert_, die
    Use m_mpout, only : mpout, mpout_log
    Implicit None
    Type (SparseComm), Intent(In) :: comm
    Real,    Intent(In)  :: x_in(:,1+comm%offset_in:)  
    Real,    Intent(Out) :: x_out(:,1+comm%offset_out:)
  !==================================================

    Character(Len=*), Parameter :: myname_ = myname//'::PartialComm'

    Pointer(p_buffer, buffer)
    Real :: buffer(size(x_in,1),max(size(x_in,2),size(x_out,2))) ! nvecs by nobs
    Pointer(p_remote, remote_buffer)
    Real :: remote_buffer(size(x_in,1),max(size(x_in,2),size(x_out,2))) ! nvecs by nobs


    Integer :: nvecs, nobs, n_rounds, round
    Integer :: n_data, n_packets, i_packet
    Integer :: istart, istop, ii
    Integer :: ier
    Integer :: remote_pe

    Type (PacketDesc), Pointer :: packets(:) 
!------------------------------------------------------------------------------------


    nvecs = Size(x_in,1)
    nobs  = Max(size(x_in,2),Size(x_out,2))  ! one or the other will span all obs.
    ASSERT(nvecs == Size(x_out,1))

    n_data      = nvecs * nobs 
    Call ReserveSymmetricMemory(n_words_real = n_data)
    p_buffer    = p_sym_buffer_i

    n_packets = Size(comm%prefill)

    Do i_packet = 1, n_packets
       istart = 1      + comm%prefill(i_packet)%displ
       istop  = istart + comm%prefill(i_packet)%count - 1
       buffer(:,istart:istop) = x_in(:,istart:istop)
    End Do

    x_out = 0

! Synchronize once
! ----------------
    Call shmem_barrier_all() ! maybe remove later?

    n_packets = Size(comm%packets)
    Do i_packet = 1, n_packets

#ifndef NDEBUG
       Call mpout_log(myname_,'     packet',i_packet)
       Call mpout_log(myname_,'         pe',comm%packets(i_packet)%pe)
       Call mpout_log(myname_,'      count',comm%packets(i_packet)%count)
#endif

       p_remote = SHMEM_PTR(buffer, comm%packets(i_packet)%pe)

       istart = 1      + comm%packets(i_packet)%displ
       istop  = istart + comm%packets(i_packet)%count - 1
       Select Case (comm%operation)
       Case (GATHER)
          x_out(:,istart:istop) = remote_buffer(:,istart:istop)
       Case (REDUCE_SCATTER)
          x_out(:,istart:istop) = x_out(:,istart:istop) + remote_buffer(:,istart:istop)
       Case Default
          Call die(myname_,'illegal operation = ',comm%operation)
       End Select
    End Do

! Synchronize twice
! ----------------
    Call shmem_barrier_all() ! maybe remove later?

  End Subroutine PartialComm_r2

  !==========================================================
  Subroutine ReserveSymmetricMemory(n_words_int, n_words_real)
    Use m_die, only : perr
    Use m_mall, only : mall_ci, mall_co, mall_ison
    use m_mpout, only : mpout, mpout_log
    Implicit None
    Integer, Intent(In), Optional :: n_words_int
    Integer, Intent(In), Optional :: n_words_real
  !==========================================================
    
    External SHPALLOC
    Character(Len=*), Parameter :: myname_ = myname//'::ReserveSymmetricMemory'
    Integer :: n_words_tot ! IRIX assumes 4 byte words, PSAS uses 8 byte reals
    Integer :: ier
    Integer, Parameter :: abort = 0 ! do not abort - return error codes instead
  !!! Error return codes
    Integer, Parameter :: non_negative_length   = -1
    Integer, Parameter :: insufficient_memory   = -2
    Integer, Parameter :: addr_out_of_bounds    = -3
    Integer, Parameter :: block_already_free    = -4
    Integer, Parameter :: addr_not_at_beginning = -5


#ifdef sgi
    Integer, Parameter :: REAL_STORAGE_SIZE = 2
#else
    Integer, Parameter :: REAL_STORAGE_SIZE = 1
#endif 


    n_words_tot = 0
    If (Present(n_words_int)) Then
       n_words_tot = n_words_tot + n_words_int
    End If
    If (Present(n_words_real)) Then
       n_words_tot = n_words_tot + REAL_STORAGE_SIZE * n_words_real
    End If

    If (symmetric_buffer_size < n_words_tot) Then

        If (symmetric_buffer_size /= 0) Then
           If (mall_ison()) Call mall_co(symmetric_buffer_size,myname)
           Call SHPCLMOVE(p_sym_buffer_i, n_words_tot, ier, abort)
        Else
           Call SHPALLOC(p_sym_buffer_i, n_words_tot, ier, abort)
        End If

        Select Case (ier)
        Case (non_negative_length)
           Call perr(myname_, 'n_words_tot not greater than 0')
        Case (insufficient_memory)
           Call perr(myname_, 'insufficient memory available on the symmetric heap')
        Case (addr_out_of_bounds) 
           Call perr(myname_, 'addresss outside of bounds of symmetric heap')
        Case (block_already_free) 
           Call perr(myname_, 'block already free')
        Case (addr_not_at_beginning) 
           Call perr(myname_, 'address not at the beginning of a block')
        Case Default
           ! everything a-ok
           symmetric_buffer_size = n_words_tot
           p_sym_buffer_r = p_sym_buffer_i
        End Select

        If (mall_ison()) Call mall_ci(n_words_tot,myname)
        
     End If
    
  End Subroutine ReserveSymmetricMemory
  !======================================================================
  Logical Function CleanQ(comm)
    Implicit None
    Type (SparseComm), Intent(In) :: comm
  !======================================================================
    Character(Len=*), Parameter :: myname_ = myname//'::Clean'

    If (comm%npes == -1) Then
       CleanQ = .true.
    Else
       CleanQ = .not. (Associated(comm%prefill)  .or. Associated(comm%packets))
    End If

  End Function CleanQ
#endif
End Module m_SparseComm_SHMEM_PTR





