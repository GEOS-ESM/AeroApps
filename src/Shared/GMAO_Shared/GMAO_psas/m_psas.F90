!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_psas - miscellaneous routines for orchestrating MPI psas
!
! !DESCRIPTION:
!
! !INTERFACE:

      module m_psas
      Implicit None
      Private

      Public :: distribute_pure_openmp
      Public :: distribute_hybrid_openmp
      Public :: psas_init
      Public :: psas_end
      Public :: HYBRID
      Public :: PURE_OMP
      Public :: PURE_MPI

      Interface distribute_hybrid_openmp
        module procedure distribute_hybrid_openmp_
      End Interface

      Interface distribute_pure_openmp
        module procedure distribute_pure_openmp_
      End Interface

!
!     Routines to initialize PSAS - only need by Jim Taft and MPI version
!
      Character(Len=*), Parameter :: myname = 'm_psas'
      Logical :: init = .false.

      Integer, Parameter :: PURE_MPI = 0
      Integer, Parameter :: HYBRID   = 1
      Integer, Parameter :: PURE_OMP = 2

    CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE: psas_init --- Initialization of Taft and MPI PSAS
!
! !INTERFACE:
!
      subroutine psas_init(I_AM_ROOT, comm, topology)
	Use m_mpif90, only : MP_COMM_WORLD, MP_INIT, MP_COMM_RANK, MP_initialized
        Use m_die, only : die, mp_die
	Integer, Optional, Intent(Out) :: I_AM_ROOT
        Integer, Optional, Intent(Out) :: comm
        Integer, Optional, Intent(In)  :: topology

! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  05mar2002  Clune <clune@sgi.com> initial implementation
!  21Mar2002  Todling added this prologue; made PURE_MPI default
!  10Sep2002  Todling turned i_am_root in to integer (root value)
!
!EOP
!-----------------------------------------------------------------------

        Character(Len=*), Parameter :: myname_ = myname//'::psas_init'
        Integer, Parameter :: ROOT = 0
        Integer :: rank
        Integer :: ier
        Integer :: topology_
        Logical :: init_mpi


!     Initialize MPI - must be done first
!     -----------------------------------
      call MP_initialized (init_mpi,ier)
      if (.not.init_mpi) then
         Call MP_INIT(ier)
            if ( ier /= 0) call mp_die ( myname_, 'MP_INIT()',ier)
      endif
      Call MP_COMM_RANK(MP_COMM_WORLD, rank, ier)
         if ( ier /= 0) call mp_die ( myname_, 'MP_COMM_RANK()',ier)	
      If (Present(comm)) comm = MP_COMM_WORLD

      topology_ = PURE_MPI

      If (Present(topology)) topology_ = topology

      If (Present(I_AM_ROOT)) I_AM_ROOT = ROOT

!      Enable dynamic assignment of the number of OpenMP threads for each MPI process.
!      -------------------------------------------------------------------------------

#ifdef USE_HYBRID
     Call OMP_SET_DYNAMIC(.true.)
     Select Case (topology_)
     Case (HYBRID)
        Call distribute_hybrid_openmp()
     Case (PURE_OMP)
        Call distribute_pure_openmp()
     Case Default
        Call die(myname_, 'Invalid topology')
     End Select
     Call OMP_SET_DYNAMIC(.FALSE.)
#endif
      end subroutine psas_init


      subroutine psas_end()
         Use m_mpif90, only : MP_Finalize
         Use m_die ,only : MP_die

         integer :: ier

         Call MP_Finalize(ier)
            if(ier/=0) call MP_die(myname,'MP_finalize()',ier)

      end subroutine psas_end


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: distribute_pure_openmp_() --- Use all CPUs for OpenMP threads on MPI root
!
! !INTERFACE:

    subroutine distribute_pure_openmp_ (maxthreads )

! !USES:
      use m_mpif90
      use m_die
      use m_mpout
      use m_zeit

    Implicit NONE

!
! !INPUT PARAMETERS:
    Integer, Optional, Intent(In) :: maxthreads
!
!
! !INPUT/OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: Sets the number of OMP threads to be the total number of available processors
!               on the mpi root process.  All other mpi processes are set to have one thread,
!               which presumably remains idle.  On SGI the root threads are explicitly assigned
!               to cpus.
!
! !REVISION HISTORY:
!
!  08Oct2001   T. Clune    Initial version
!
!EOP
!-------------------------------------------------------------------------

    character(len=*), parameter ::  myname_ = myname//'::distribute_pure_openmp_'

    Integer :: ier
    Integer :: rank, rankp
    Integer :: npes_mpi, n_cpu, cpu
    Integer :: n_threads, i_thread, j
    character(len=255) :: env_val

!   Access MPI info
!   ---------------

#if defined(USE_HYBRID) && defined(sgi)
    call zeit_ci('pure 1')
    Call MP_COMM_RANK(MP_COMM_WORLD, rank, ier)
         if ( ier /= 0) call die ( myname_, 'MP_COMM_RANK()')
    Call MP_COMM_SIZE(MP_COMM_WORLD, npes_mpi, ier)
         if ( ier /= 0) call die ( myname_, 'MP_COMM_SIZE()')

    Call getenv('N_CPU',env_val)
    if ( len_trim(env_val) .gt. 0 ) then
         read(env_val,*) n_cpu
    else
       call die('distribute_pure_openmp_','env variable N_CPU not set' )
    end if

    If (n_cpu <= 0) Call die( myname_, ' N_CPU must be positive ' )
    If (mod(n_cpu, npes_mpi) /= 0) call warn ( myname_, ' uneven PE distribution ' )

    call zeit_co('pure 1')
    call zeit_ci('pure 2')

!   Kill existing threads if already in progress
!   --------------------------------------------
    If (init) Then
#ifdef _OPENMP
       Call mp_destroy()
#endif
       Call mpi_barrier(MP_COMM_WORLD,ier)
           If (ier /= 0) Call die( myname_, ' barrier failure')
#ifdef _OPENMP
     Call runanywhere()
#endif
    Else
       init = .true.
    End If
	
    call zeit_co('pure 2')
    call zeit_ci('pure 3')

    If (rank > 0) Then

!      Set the number of threads
!      -------------------------
       n_threads = 1
#ifdef _OPENMP
     Call OMP_SET_NUM_THREADS(n_threads)
#endif

!$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i_thread, cpu), &
!$OMP &   SHARED(n_threads, rank, n_cpu, npes_mpi, mpout)

       Do i_thread = 0, n_threads - 1

!         assign this thread to the corresponding cpu <* NONPORTABLE *>
!         -------------------------------------------------------------
          cpu = rank * n_cpu / npes_mpi
          call mp_assign_to_cpu(cpu)

       End Do
!$OMP END PARALLEL DO

    End If

    call zeit_co('pure 3')
    call zeit_ci('pure 4')

    If (rank == 0) Then

!      Set the number of threads
!      -------------------------
       n_threads = n_cpu - npes_MPI + 1
       If (Present(maxthreads)) n_threads = max(1,Min(maxthreads,n_threads))

#ifdef _OPENMP
       Call OMP_SET_NUM_THREADS(n_threads)
#endif

!$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i_thread, cpu, j, rankp), SHARED(n_threads,n_cpu,npes_mpi,mpout,rank)
       Do i_thread = 0, n_threads - 1

!         assign this thread to the corresponding cpu <* NONPORTABLE *>
!         -------------------------------------------------------------
          ! skip MPI team masters?
          cpu = 0
          rankp = 1
          Do j = 0, i_thread - 1
             cpu = cpu + 1
             If (cpu == rankp*n_cpu/npes_mpi) Then ! skip
                cpu = cpu + 1
                rankp = rankp + 1
             End If
          End Do
          Call mp_assign_to_cpu(cpu)
       End Do
!$OMP END PARALLEL DO

    End If
    call zeit_co('pure 4')

#endif

   End subroutine distribute_pure_openmp_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: distribute_hybrid_openmp_() --- Use all CPUs for OpenMP threads on MPI root
!
! !INTERFACE:

    subroutine distribute_hybrid_openmp_ ( )

! !USES:
      use m_mpif90
      use m_die
      use m_mpout

    Implicit NONE

!
! !INPUT PARAMETERS:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: Distributes available cpus among all mpi processes.
!
! !REVISION HISTORY:
!
!  08Oct2001   T. Clune    Initial version
!
!EOP
!-------------------------------------------------------------------------

    character(len=*), parameter ::  myname_ = myname//'::distribute_hybrid_openmp_'

    Integer :: ier
    Integer :: rank
    Integer :: npes_mpi
    Integer :: n_threads, n_cpu, i_thread
    Integer :: cpu_top, cpu_bot, cpu
    Character(len=100) :: env_val ! buffer for reading environment variables

!   Access MPI info
!   ---------------

#if defined(sgi) && defined(USE_HYBRID)

    Call MP_COMM_RANK(MP_COMM_WORLD, rank, ier)
         if ( ier /= 0) call die ( myname_, 'MP_COMM_RANK()')
    Call MP_COMM_SIZE(MP_COMM_WORLD, npes_mpi, ier)
         if ( ier /= 0) call die ( myname_, 'MP_COMM_SIZE()')

    Call getenv('N_CPU',env_val)
    if ( Len_Trim(env_val) .gt. 0 ) then
         read(env_val,*) n_cpu
    else
       call die('distribute_hybrid_openmp_','env variable N_CPU not set' )
    end if

    If (n_cpu <= 0) Call die( myname_, ' N_CPU must be positive ' )
    If (mod(n_cpu, npes_mpi) /= 0) call warn ( myname_, ' uneven PE distribution ' )

!   Kill existing threads
!   ---------------------
    Call mp_destroy()
    Call runanywhere()

    cpu_bot = (rank*n_cpu)/npes_mpi
    cpu_top = ((rank+1)*n_cpu)/npes_mpi - 1
    n_threads = cpu_top - cpu_bot + 1
#ifdef _OPENMP
    Call OMP_SET_NUM_THREADS(n_threads)
#endif

!$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i_thread, cpu), SHARED(n_threads,cpu_bot)
    Do i_thread = 0, n_threads - 1

!      assign this thread to the corresponding cpu <* NONPORTABLE *>
!      -------------------------------------------------------------
       cpu = cpu_bot + i_thread
       call mp_assign_to_cpu(cpu)

    End Do
!$OMP END PARALLEL DO

#endif
   End subroutine distribute_hybrid_openmp_

!........................................................................................

!
! On some systems, e.g. MacOS X, one need this trick to force linking data only modules
!

     subroutine psas0()  

       use OEclass_tbl
       use hfecQ_tbl
       use vfecQ_tbl
       use rlev_imat
       use rlat_imat
       use FEalpha_imat
       use FEsigW_imat
       use FEalpha_tabl
       use hoecH_tbl
       use voecH_tbl

       call OEclass_tbl0()
       call hfecQ_tbl0()
       call vfecQ_tbl0()
       call rlev_imat0()
       call rlat_imat0()
       call FEalpha_imat0()
       call FEsigW_imat0()
       call FEalpha_tabl0()
       call hoecH_tbl0()
       call voecH_tbl0()

     end subroutine psas0

   end module m_psas

