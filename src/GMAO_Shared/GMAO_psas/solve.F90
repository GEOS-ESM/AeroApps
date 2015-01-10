Program Solve
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !Program:  solve --- invokes PSAS from within ANA as a split executable.
!
! !INTERFACE:

  Use m_GetAI, Only : GetAI
  Use m_psas,  Only : psas_init, psas_end
  Use m_psas,  Only : HYBRID
  Use m_mpif90,Only : MP_comm_rank
  Use m_mpout, Only : mpout_open, mpout_close, mpout
  Use m_mpout, Only : mpout_setflush, mpout_log
  Use m_zeit,  Only : zeit_ci, zeit_co, MWTIME
  Use m_zeit,  Only : zeit_flush, zeit_allflush
  Use m_die,   Only : MP_die
  
  Implicit None

! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
!  05mar2002  Clune <clune@sgi.com> initial implementation
!  21Mar2002  Todling ifdef'ed HYBRID initialization
!
!EOP
!-----------------------------------------------------------------------
  Character(Len=*), Parameter :: myname = 'solve.x'
  Integer :: ier
  Integer :: rank
  Integer :: root
  Integer :: comm
  Logical :: i_am_root
  

!  Initialize parallelism
#ifdef USE_HYBRID
  Call psas_init(root, comm, topology = HYBRID)
#else
  Call psas_init(root, comm)
#endif

! Inquire processor ID
! --------------------
  call MP_comm_rank ( comm,rank,ier )
    if(ier/=0) call MP_die(myname,'MP_comm_rank()',ier)
  i_am_root = (rank == root)

!  Call mpout_open(pfix='psas',mask=127)
  Call mpout_setflush(flush = .true.)
  Call mpout_log(myname,'Starting PSAS Solver... ')

     Call zeit_ci('solve')  
  Call GetAI(comm, root)
     Call zeit_co('solve')

!  call mpout_close()


!  Call mpout_open(pfix='zeit_psas',mask=127)
  If  (I_am_root)  call zeit_flush(lu=mpout,umask=MWTIME)
  call zeit_allflush(comm=comm,root=root,lu=mpout,	umask=MWTIME)

  Call mpout_log(myname,'Done with PSAS Solver... ')

  call mpout_close()
  call psas_end()

End Program Solve


