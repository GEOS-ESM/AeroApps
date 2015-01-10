!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_SparseComm_MPI_1
!
! !DESCRIPTION:
!
! !INTERFACE:


Module m_SparseComm_MPI_1
#ifdef USE_MPI_1
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
  End Interface

! !REVISION HISTORY:
!       11Feb02 - Tom Clune <clune@sgi.com>
!               . Initial MPI_1 version (based on SHMEM version)
!       29Oct01 - Tom Clune
!               . Finally added this intro section

!EOP ___________________________________________________________________

  Character(Len=*), Parameter :: myname = 'm_SparseComm_MPI_1'

  !!! supported operations
  Integer, Parameter, Public :: GATHER         = 1
  Integer, Parameter, Public :: REDUCE_SCATTER = 2

!!! this include file is necessary for some shmem functions
  !Include 'mpif.h'

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

     Integer          :: operation  = -1
     Integer          :: mpi_comm   = -1
     Integer          :: npes       = -1
     Integer          :: my_id      = -1
     Integer          :: offset_in  = -1
     Integer          :: offset_out = -1
     Type (PacketDesc), Pointer :: recv_packets(:) => Null()
     Type (PacketDesc), Pointer :: send_packets(:) => Null()
     Integer,           Pointer :: requests(:)     => Null()
  End Type SparseComm

Contains ! Functions and subroutines

  Subroutine Analyze(comm, seg_in, seg_out, nav, offset, operation, mp_comm)
    Use m_Navigator, only : Navigator
    Use m_Navigator, only : lsize
    Use m_Navigator, only : get
    Use m_die, only : die, mp_die, assert_
    Use m_mpif90, only : MP_TYPE, MP_COMM_SIZE, MP_COMM_RANK
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
    
    Logical, Allocatable :: seg_remote(:,:)

    Logical :: transmit

    Integer :: my_id, np
    Integer :: segsize
    Integer :: p, pp
    Integer :: ier
    Integer :: n_words
    Integer :: n_packets, i_packet, n_segments, iseg, iseg1, nseg, cnt, displ

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

    Call MP_comm_size(mp_comm, np, ier)
    	if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)
    Call MP_comm_rank(mp_comm, my_id, ier)
    	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)
    comm%npes     = np
    comm%my_id    = my_id
    comm%mpi_comm = mp_comm

    If (np == 1) Then
       Allocate(comm%recv_packets(0), comm%send_packets(0), comm%requests(0))
       return ! skip analysis - trivial case
    End If
! Store possession flag in a global arena for easy access
! -------------------------------------------------------
    n_words = n_segments ! required

!-------------------------------------------------------------------------
! Gather possession information
! -------------------------------------------------
   Allocate(seg_remote(n_segments,0:np-1), STAT = ier)
      if (ier /= 0) call die(myname_,'Allocate()',ier)

    Call MPI_AllGather(seg_in,   Size(seg_in),  MP_TYPE(seg_in), &
         &            seg_remote, Size(seg_in), MP_TYPE(seg_remote), mp_comm, ier)
    	if(ier/=0) call MP_die(myname_,'MPI_AllGather()',ier)
    n_packets = 0
    Do pp = 0, np - 1 ! stagger lookups
       p = mod(pp + my_id, np)
       n_packets = n_packets + Count_Contiguous(seg_out .and. seg_remote(:,p))
    End Do

    Allocate(comm%recv_packets(n_packets), STAT = ier)
    ALWAYS_ASSERT(ier == 0)

! Loop over processors and fill packet descriptions
! -------------------------------------------------

    i_packet = 0
    Do pp = 0, np - 1 ! stagger lookups
       p = mod(pp + my_id, np)

       transmit = .false.
       Do iseg = 1, n_segments
          If (seg_out(iseg) .and. seg_remote(iseg,p)) Then
             If (transmit) Then ! already in a packet
                comm%recv_packets(i_packet)%count = comm%recv_packets(i_packet)%count + 1
             Else ! new packet
                i_packet = i_packet + 1
                comm%recv_packets(i_packet)%pe    = p
                comm%recv_packets(i_packet)%displ = iseg
                comm%recv_packets(i_packet)%count = 1
                transmit = .true.
             End If
          Else
             transmit = .false.
          End If
       End Do
   End Do

   ALWAYS_ASSERT(i_packet == n_packets)

!-------------------------------------------------------------------------
! Gather possession information
! -------------------------------------------------

    Call MPI_AllGather(seg_out,  Size(seg_out), MP_TYPE(seg_out), &
         &            seg_remote, Size(seg_out), MP_TYPE(seg_remote), mp_comm, ier)
    	if(ier/=0) call MP_die(myname_,'MPI_AllGather()',ier)

    n_packets = 0
    Do pp = 0, np - 1 ! stagger lookups
       p = mod(pp + my_id, np)
       n_packets = n_packets + Count_Contiguous(seg_in .and. seg_remote(:,p))
    End Do

    Allocate(comm%send_packets(n_packets), STAT = ier)
    ALWAYS_ASSERT(ier == 0)

! Loop over processors and fill packet descriptions
! -------------------------------------------------

    i_packet = 0
    Do pp = 0, np - 1 ! stagger lookups
       p = mod(pp + my_id, np)

       transmit = .false.
       Do iseg = 1, n_segments

          If (seg_in(iseg) .and. seg_remote(iseg,p)) Then
             If (transmit) Then ! already in a packet
                comm%send_packets(i_packet)%count = comm%send_packets(i_packet)%count + 1
             Else ! new packet
                i_packet = i_packet + 1
                comm%send_packets(i_packet)%pe    = p
                comm%send_packets(i_packet)%displ = iseg
                comm%send_packets(i_packet)%count = 1
                transmit = .true.
             End If
          Else
             transmit = .false.
          End If
       End Do
   End Do

   ALWAYS_ASSERT(i_packet == n_packets)

   Deallocate(seg_remote)

! Convert displacement and count units from segments to points
! ------------------------------------------------------------

   Do i_packet = 1, size(comm%recv_packets)

      iseg = comm%recv_packets(i_packet)%displ
      nseg = comm%recv_packets(i_packet)%count


      Call Get(nav, iseg, displ = displ)
      comm%recv_packets(i_packet)%displ = displ
      
      cnt = 0
      Do iseg1 = iseg, iseg + nseg - 1
         Call Get(nav, iseg1, ln = segsize)
         cnt = cnt + segsize
      End Do
      comm%recv_packets(i_packet)%count = cnt

   End Do

! Convert displacement and count units from segments to points
! ------------------------------------------------------------

   Do i_packet = 1, size(comm%send_packets)

      iseg = comm%send_packets(i_packet)%displ
      nseg = comm%send_packets(i_packet)%count

      Call Get(nav, iseg, displ = displ)
      comm%send_packets(i_packet)%displ = displ         
      
      cnt = 0
      Do iseg1 = iseg, iseg + nseg - 1
         Call Get(nav, iseg1, ln = segsize)
         cnt = cnt + segsize
      End Do
      comm%send_packets(i_packet)%count = cnt

   End Do

! Allocate a handy array to store non-blocking requests

   n_packets = size(comm%send_packets)

   Select Case (operation)
   Case (Gather)
      n_packets = size(comm%recv_packets)
   Case (REDUCE_SCATTER)
      n_packets = size(comm%send_packets)
   End Select
   Allocate(comm%requests(n_packets), STAT = ier)
   if(ier/=0) call die(myname_,'Allocate()',ier)

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
    use m_die, only : assert_
    Implicit None
    Type (SparseComm), Intent(InOut) :: comm
  !===================================================================
    Character(Len=*), Parameter :: myname_ = myname//'::Clean'

    Integer :: ier

    If (initQ(comm)) Then
       comm%npes = -1
       Deallocate(comm%recv_packets, comm%send_packets, &
            & comm%requests, STAT = ier)
          ALWAYS_ASSERT(ier == 0)
    End If

  End Subroutine Clean

  !===================================================================
  Logical Function initQ(comm)
    Type (SparseComm), Intent(In) :: comm
    initQ = (comm%npes >= 1)
  End Function initQ

  !==================================================
  Subroutine PartialReduceScatter(x_in, x_out, comm)
    Use m_mpif90, only : MP_TYPE
    Use m_mpif90, only : MP_STATUS_SIZE
    Use m_mpif90, only : MP_MAX_ERROR_STRING
#ifdef LATER
    Use m_mpif90, only : MP_ERROR
#endif
    Use m_die, only : assert_, mp_die, die
    Use m_mpout, only : mpout, mpout_log, mpout_flush
    Implicit None
    Type (SparseComm), Intent(In) :: comm
    Real,    Intent(In)  :: x_in(:,1+comm%offset_in:)  
    Real,    Intent(Out) :: x_out(:,1+comm%offset_out:)
  !==================================================
    Character(Len=*), Parameter :: myname_ = myname//'::PartialReduceScatter'

    Integer :: nvecs, nobs
    Integer :: n_packets, i_packet
    Integer :: istart, istop, lb, ub
    Integer :: ier
    Integer :: itag, pe, icount, pe_last
    Integer, Parameter :: MAX_PACKETS = 10000
    Real, Allocatable :: buffer(:,:)

    Integer :: status(MP_STATUS_SIZE)
    Integer :: mpi_stati(MP_STATUS_SIZE,1000)

    Character(len=MP_MAX_ERROR_STRING) :: msg
    Integer :: msg_len

!------------------------------------------------------------------------------------

    ASSERT(comm%operation == REDUCE_SCATTER)
    If (comm%npes == 1) Then ! trivial case
       x_out = x_in
       Return
    End If

    nvecs = Size(x_in,1)
    nobs  = Max(size(x_in,2),Size(x_out,2))  ! one or the other will span all obs.
    ASSERT(nvecs == Size(x_out,1))

    x_out = 0

    n_packets = Size(comm%send_packets)
    pe_last = -1
    Do i_packet = 1, n_packets

       pe = comm%send_packets(i_packet)%pe
       istart = 1      + comm%send_packets(i_packet)%displ
       icount = nvecs * comm%send_packets(i_packet)%count
       If (pe == pe_last) Then
          itag = itag + 1
       Else
          itag = 1
       End If
       pe_last = pe
       
       Call mpi_isend(x_in(1,istart), icount, MP_TYPE(x_in), pe, itag, &
            & comm%mpi_comm, comm%requests(i_packet),ier)
          if(ier/=0) call MP_die(myname_,'MPI_isend()',ier)

    End Do

    lb = lbound(x_out,2)
    ub = ubound(x_out,2)
    Allocate(buffer(nvecs,lb:ub), STAT = ier)
          if(ier/=0) call die(myname_,'Allocate()',ier)

    n_packets = Size(comm%recv_packets)
    pe_last = -1
    Do i_packet = 1, n_packets
       pe = comm%recv_packets(i_packet)%pe
       istart = 1 + comm%recv_packets(i_packet)%displ
       icount = comm%recv_packets(i_packet)%count

       istop  = istart + icount - 1
       icount = icount * nvecs
       If (pe == pe_last) Then
          itag = itag + 1
       Else
          itag = 1
       End If
       pe_last = pe

       ASSERT(istart >= lb)
       ASSERT(istop <= ub)
       Call MPI_RECV(buffer(1,istart), icount, MP_TYPE(buffer), pe, itag, &
            comm%mpi_comm, status, ier)
          if(ier/=0) call MP_die(myname_,'MPI_RECV()',ier)

       x_out(:,istart:istop) = x_out(:,istart:istop) + buffer(:,istart:istop)

    End Do
    Deallocate(buffer)

    Call MPI_WAITALL(size(comm%requests), comm%requests, MPI_STATI, ier)
          if(ier /= 0) Then
#ifdef LATER
             Do i_packet = 1, Size(comm%requests)
                ier = MPI_STATI(MP_ERROR,i_packet)
                if(ier /= 0) Then
                   Call mpout_log(myname_,'MPI_WAIT  packet = ', i_packet)
                   Call mpout_log(myname_,'MPI_WAIT     ier = ', ier)
                   Call mpout_log(myname_,'MPI_WAIT npacket = ', size(comm%requests))
                   call mpout_flush()

                   msg_len = MP_MAX_ERROR_STRING
                   Call mpi_error_string(MPI_STATI(MP_ERROR,i_packet), msg, msg_len, ier)
                   call mpout_log(myname_,'MPI_error ier = ',ier)
                   call mpout_log(myname_,'MPI_error message length: ',msg_len)
                   Write(mpout,*)myname_,'    message = ',msg(1:msg_len)
                   call mpout_flush()
                   Call MP_die(myname_,'MPI_WAITALL()',ier)
                   
                End if
             End Do
#else
                   Call MP_die(myname_,'MPI_WAITALL()',ier)
#endif
          end if

  End Subroutine PartialReduceScatter

  !==================================================
  Subroutine pgath_r2(x_in, x_out, comm)
    Use m_mpif90, only : MP_TYPE
    Use m_mpif90, only : MP_STATUS_SIZE
    Use m_mpif90, only : MP_MAX_ERROR_STRING
#ifdef LATER
    Use m_mpif90, only : MP_ERROR
#endif
    Use m_die, only : assert_, mp_die, die
    Use m_mpout, only : mpout, mpout_log, mpout_flush
    Implicit None
    Type (SparseComm), Intent(In) :: comm
    Real,    Intent(In)  :: x_in (:,1+comm%offset_in:)  
    Real,    Intent(Out) :: x_out(:,1+comm%offset_out:)
  !==================================================
    Character(Len=*), Parameter :: myname_ = myname//'::PartialGather'

    Integer :: nvecs, nobs
    Integer :: n_packets, i_packet
    Integer :: istart, istop
    Integer :: ier
    Integer :: itag, pe, icount, pe_last

    Integer :: mpi_stati(MP_STATUS_SIZE,size(comm%requests))
    Character(len=MP_MAX_ERROR_STRING) :: msg
    Integer :: msg_len

!------------------------------------------------------------------------------------

    ASSERT(comm%operation == GATHER)

    If (comm%npes == 1) Then ! trivial case
       x_out = x_in
       Return
    End If

    nvecs = Size(x_in,1)
    ASSERT(nvecs == Size(x_out,1))

    nobs  = size(x_out,2)

#ifndef NDEBUG
    x_out = -Huge(1.) ! illegal value - should never actually be used
#endif

! post non-blocking receives for all packets
    n_packets = Size(comm%recv_packets)
    pe_last = -1
    Do i_packet = 1, n_packets
       pe = comm%recv_packets(i_packet)%pe
       istart = 1 + comm%recv_packets(i_packet)%displ
       icount = comm%recv_packets(i_packet)%count

       istop  = istart + icount - 1
       icount = icount * nvecs
       If (pe == pe_last) Then
          itag = itag + 1
       Else
          itag = 1
       End If
       pe_last = pe

       Call MPI_IRECV(x_out(1,istart), icount, MP_TYPE(x_out), pe, itag, &
            comm%mpi_comm, comm%requests(i_packet), ier)
           if(ier/=0) call MP_die(myname_,'MPI_IRECV()',ier)
    end Do

    n_packets = Size(comm%send_packets)
    pe_last = -1
    Do i_packet = 1, n_packets

       pe = comm%send_packets(i_packet)%pe
       istart = 1      + comm%send_packets(i_packet)%displ
       icount = nvecs * comm%send_packets(i_packet)%count
       If (pe == pe_last) Then
          itag = itag + 1
       Else
          itag = 1
       End If
       pe_last = pe
       
       Call mpi_send(x_in(1,istart), icount, MP_TYPE(x_in), pe, itag, &
            & comm%mpi_comm, ier)
          if(ier/=0) call MP_die(myname_,'MPI_send()',ier)

    End Do

    Call MPI_WAITALL(size(comm%requests), comm%requests, MPI_STATI, ier)
        if(ier/=0) Then
#ifdef LATER
           Do i_packet = 1, Size(comm%requests)
              ier = MPI_STATI(MP_ERROR,i_packet)
              if(ier /= 0) Then
                 Call mpout_log(myname_,'MPI_WAIT  packet = ', i_packet)
                 Call mpout_log(myname_,'MPI_WAIT     ier = ', ier)
                 Call mpout_log(myname_,'MPI_WAIT npacket = ', size(comm%requests))
                 write(mpout,*)'sources? ',comm%recv_packets(:)%pe
                 msg_len = MP_MAX_ERROR_STRING
                 Call mpi_error_string(MPI_STATI(MP_ERROR,i_packet), msg, msg_len, ier)
                 call mpout_log(myname_,'MPI_error ier = ',ier)
                 call mpout_log(myname_,'MPI_error message length: ',msg_len)
                 Write(mpout,*)myname_,'    message = ',msg(1:msg_len)
                 call mpout_flush()
                 Call MP_die(myname_,'MPI_WAITALL()',ier)
              End if
           End Do
#else
                 Call MP_die(myname_,'MPI_WAITALL()',ier)
#endif
        end if

  End Subroutine pgath_r2

  !==============================================================
  Subroutine DescribeComm(comm, label)
    Use m_mpout, only : mpout_log, mpout, mpout_flush
    Type (SparseComm), Intent(In) :: comm
    Character(Len=*), Intent(In)  :: label
  !==============================================================

#ifdef NDEBUG
    Write(mpout,*)'*****************************************'
    Write(mpout,*)'*** SparseComm (MPI-1) Description '//Trim(label)
    Write(mpout,*)'***    npes       = ', comm%npes
    Write(mpout,*)'***    my_id      = ', comm%my_id
    Write(mpout,*)'***    mpi comm   = ', comm%mpi_comm
    Write(mpout,*)'***    offset in  = ', comm%offset_in
    Write(mpout,*)'***    offset out = ', comm%offset_out
    Write(mpout,*)'*****************************************'
    Call mpout_flush()
#endif

  End Subroutine DescribeComm

  !==========================================================
  Subroutine ReleaseSymmetricMemory()
    Use m_die, only : perr
    Use m_mall, only : mall_ci, mall_co, mall_ison
    Implicit None
  !==========================================================
    Character(Len=*), Parameter :: myname_ = myname//'::ReleaseSymmetricMemory'

    ! BLANK STUB - unnecessary for MPI-1 version

  End Subroutine ReleaseSymmetricMemory ! static

  !==========================================================
  Subroutine ReserveSymmetricMemory(n_words_int, n_words_real)
    Use m_die, only : perr
    Use m_mall, only : mall_ci, mall_co, mall_ison
    Implicit None
    Integer, Intent(In), Optional :: n_words_int
    Integer, Intent(In), Optional :: n_words_real
  !==========================================================
    Character(Len=*), Parameter :: myname_ = myname//'::ReserveSymmetricMemory'    
    
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
       CleanQ = .not. (Associated(comm%send_packets) .or. Associated(comm%recv_packets))
    End If

  End Function CleanQ

#endif
End Module m_SparseComm_MPI_1




