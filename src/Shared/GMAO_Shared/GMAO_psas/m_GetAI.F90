!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_GetAI --- a handshake between OpenMP ANA and MPI-PSAS
!
! !INTERFACE:
!

   MODULE  m_GetAI

! !USES:

     use m_FileResolv, Only: ExecResolv
     Use m_die,      Only : MP_die
     Use m_die,      Only : die

     Implicit None
!
! !PUBLIC MEMBER FUNCTIONS:
!

  Public :: GetAI
  Public :: GetAI_External ! ana.x interface
  Public :: GetAI_Internal ! solve.x interface
  Public :: GetAI_setAOD_sigF
  Public :: get_sigf_W

  Interface GetAI
     Module Procedure GetAI_External
     Module Procedure GetAI_Internal
  End Interface!

  Interface GetAI_setAOD_sigF; Module Procedure setAOD_sigF; End Interface
  Interface get_sigf_W; Module Procedure get_sigf_W_; End Interface

! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
!  03mar2002  Clune     Initial implementation.
!  04Mar2003  Todling   Move moisture-related mods from m_ana here:
!                       wired-in ktWW and ktTPW (this is bad)
!
!EOP
!-------------------------------------------------------------------------

  Character(Len=*), Parameter :: myname = 'm_GetAI'

  Character(Len=*), Parameter :: file_to_psas   = 'ana-psas.bin'
  Character(Len=*), Parameter :: file_from_psas = 'psas-ana.bin'

  integer, external :: system

! The following should really come from m_odsmeta
! -----------------------------------------------
  integer, parameter  :: ktWW    =  7
  integer, parameter  :: ktTPW   = 18
  integer, parameter  :: ktSkinT = 38
  integer, parameter  :: ktAOD   = 45

! Analysis level assigned to TPW observations
! -------------------------------------------
  real, parameter     :: levTPW = 850.      ! (hPa)

Contains

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE: GetAI_Internal --- Serial code invoked from ANA to drive psas
!
! !INTERFACE:
!
  Subroutine GetAI_Internal(comm, root)
! !USES:
!
    
    Use m_AISolver, Only : Solve
    Use m_sigFi_lookups, Only : sigFi_remove
    Use m_sigFi_lookups, Only : MIXRE
    Use m_mpif90,   Only : MP_comm_rank
    Use m_mpif90,   Only : MP_type
    Use m_ioutil,   Only : luavail, opnieee, clsieee
    Use m_zeit,     Only : zeit_ci, zeit_co
    

    Integer, Intent(In) :: comm
    Integer, Intent(In) :: root
!
! !DESCRIPTION: Loads user options from resource file.
!
! !REVISION HISTORY:
!
!  03mar2002  Clune     Initial implementation.
!  22mar2002  Todling   Small mod in error msg to distinguishing them
!  10sep2002  Todling   Passing root value instead of who's root
!
!-------------------------------------------------------------------------

    Character(Len=*), Parameter :: myname_ = myname // '::GetAI_Internal'

    Integer :: ier
    Integer :: rank
    Integer :: unit
    Logical :: ts_ana, aod_ana
    Logical :: i_am_root

    Integer              :: ninc
    Real,    Allocatable :: ainc(:), alat(:), alon(:), alev(:)
    Integer, Allocatable :: akts(:)
    Integer              :: nobs
    Real,    Allocatable :: xvec(:), lat(:), lon(:), tlev(:), ttim(:)
    Integer, Allocatable :: kx(:), ks(:), tkt(:)
    Real,    Allocatable :: txmso(:), tomf(:)

    ! 0) Initialize communications/parallelism
    ! 1) Read data
    ! 1.1) Open the input file
    ! 1,2) Read Size information
    ! 1.3) Allocate arrays
    ! 1.4) Read array values
    ! 1.5) Close file

    ! 2) Invoke Solve in PSAS

    ! 3) Write Data
    ! 3.1) Open output file
    ! 3.2) Write size information
    ! 3.3) Write array values
    ! 3.4) Close file
    ! 3.5) Deallocate arrays

    ! 4) Close parallelism

!   Inquire processor ID
!   --------------------
    call MP_comm_rank ( comm,rank,ier )
      if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)
    i_am_root = ( rank==root )

! Read data
!   Open the input file
!   -------------------
    If (i_am_root) Then
       unit = luavail()
       Call opnieee(unit, Trim(file_to_psas), 'old', ier)
       If (ier /= 0) Call die(myname_,'opnieee("'//Trim(file_to_psas)//'")', ier)
    End If
    
!   Read Size information
!   ---------------------
    If (i_am_root) Then
       Read(unit) ninc, nobs
    Else
       ninc = 0
       nobs = 0
    ENd If

!   Allocate arrays
!   ---------------
    Allocate(ainc(ninc), alat(ninc), alon(ninc), alev(ninc), akts(ninc), STAT = ier)
       If (ier /= 0) Call die(myname_,'allocate() A', ier)
    Allocate(xvec(nobs), lat(nobs), lon(nobs), tlev(nobs), ttim(nobs), kx(nobs), &
           ks(nobs), tkt(nobs), txmso(nobs), tomf(nobs), STAT = ier)
       If (ier /= 0) Call die(myname_,'allocate() B', ier)
    
!   Read array values
!   -----------------
    If (i_am_root) Then
       Read(unit) alat, alon, alev, akts
       Read(unit) lat, lon, tlev, ttim, kx, ks, tkt, txmso, tomf
    End If
       
!   Close file
!   ----------
    If (i_am_root) Then
       Call clsieee(unit, ier, status = 'DELETE')
          If (ier /= 0) Call die(myname_,'clsieee("'//Trim(file_to_psas)//'")', ier)
    End If

!   Invoke Solve in PSAS
!   --------------------
!*****************************************************************************
    If (i_am_root) Then
       ts_ana = .false.
       aod_ana = .false.
       if ( any(tkt(1:nobs)==ktSkinT) ) then
            ts_ana = .true.
            if( any(tkt(1:nobs)/=ktSkinT) ) &
                call MP_die(myname_,'SkinT-kt being analyzed w/ other kt',98)
       end if
       if ( any(tkt(1:nobs)==ktAOD) ) then
            aod_ana = .true.
            if( any(tkt(1:nobs)/=ktAOD) ) &
                call MP_die(myname_,'AOD-kt being analyzed w/ other kt',98)
       else
            
       end if

       where ( tkt .eq. ktSkinT )
               tkt = ktWW   ! pretend this is single-level mixing ratio
       end where

       where ( tkt .eq. ktAOD )
               tkt = ktWW   ! pretend this is single-level mixing ratio
       end where

    End If

    call MPI_bcast(ts_ana, 1,MP_type(ts_ana),ROOT,comm,ier)
    call MPI_bcast(aod_ana,1,MP_type(ts_ana),ROOT,comm,ier)

    if (i_am_root) Then
       if ( aod_ana ) then
          print *, myname_//': doing AOD analysis'
       else           
          print *, myname_//': NOT doing AOD analysis, kt =', &
                   minval(tkt), maxval(tkt)
       end if
    end If


!   Define moisture/Ts sigF for PSAS
!   --------------------------------
    if ( ts_ana ) then
        call Ts_Sigf_ ( i_am_root )
    else if ( aod_ana ) then
        call setAOD_SigF ( i_am_root ) ! avoids reading of sigF file
    else
        call Moist_Sigf_ ( i_am_root )
    end if

                                                                              Call zeit_ci('psas')
    Call Solve( ninc, ainc, alat, alon, alev, akts,   &
                nobs, xvec, lat, lon, tlev,           &
                kx, ks, tkt, txmso, tomf,             &
                comm=comm, root=root )
                                                                             Call zeit_co('psas')

!   Clear definition of PSAS MIXR forecast error stats
!   --------------------------------------------------
    call sigFi_remove ( MIXRE )  

!*****************************************************************************

! Write Data
!   Open output file
!   ----------------
    If (i_am_root) Then
       unit = luavail()
       Call opnieee(unit, Trim(file_from_psas), 'unknown', ier)
       If (ier /= 0) Call die(myname_,'opnieee("'//Trim(file_from_psas)//'")', ier)
    End If

!   Write array values
!   ------------------
    If (i_am_root) Then
       print *, 'GetAIinternal:  ainc =', minval(ainc), maxval(ainc)
       Write(unit) ainc, xvec
    End If

!   Close file
!   ---------
    If (i_am_root) Then
       Call clsieee(unit, ier)
          If (ier /= 0) Call die(myname_,'clsieee("'//Trim(file_from_psas)//'")', ier)
    End If

    ! 3.5) Deallocate arrays
    Deallocate(ainc, alat, alon, alev, akts, STAT = ier)
       If (ier /= 0) Call die(myname_,'deallocate() A', ier)
    Deallocate(xvec, lat, lon, tlev, ttim, kx, &
           ks, tkt, txmso, tomf, STAT = ier)
       If (ier /= 0) Call die(myname_,'deallocate() B', ier)

  End Subroutine GetAI_Internal

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE: GetAI_External --- Serial code invoked from ANA to drive psas
!
! !INTERFACE:
!

  Subroutine GetAI_External(ninc, ainc, alat, alon, alev, akts, &
                            nobs, xvec, lat, lon, tlev, ttim, kx, ks, tkt, txmso, tomf)

! !USES:
!
    Use m_ioutil,  Only : luavail, opnieee, clsieee
    Use m_zeit,    Only : zeit_ci, zeit_co

    Integer, Intent(In)  :: ninc
    Real,    Intent(InOut) :: ainc(ninc)
    Real,    Intent(In)  :: alat(ninc), alon(ninc), alev(ninc)
    Integer, Intent(In)  :: akts(ninc)

    Integer, Intent(In) :: nobs
    Real,    Intent(InOut) :: xvec(nobs)
    Real,    Intent(In) :: lat(nobs), lon(nobs), tlev(nobs)
    Real,    Intent(In) :: ttim(nobs)
    Integer, Intent(In) :: kx(nobs), ks(nobs), tkt(nobs)
    Real,    Intent(In) :: tomf(nobs), txmso(nobs)

!
! !DESCRIPTION: Loads user options from resource file.
!
! !REVISION HISTORY:
!
!  03mar2002  Clune     Initial implementation.
!  22mar2002  Todling   Added getenv and related stuff.
!  14jun2002  da Silva  Introduced PSAS_NUM_MPI env variable
!  10sep2002  Todling   Bug fix; xvec was being treated as (In)put
!  13dec2002  Todling   Added env var name MPIRUN (halem fix)
!
!-------------------------------------------------------------------------

    Character(Len=*), Parameter :: myname_ = myname // '::GetAI_External'
    Integer :: ier
    Integer :: unit
    Integer :: npes
    Integer :: nsmp
    Character(len=8)   cnpes
    Character(len=255) mpirun
    Character(len=255) solve_dot_x
    Character(len=255) command

    ! 1) Write Data
    ! 1.1) Open output file
    ! 1.2) Write size information
    ! 1.3) Write array values
    ! 1.4) Close file

    ! 2) Invoke solve.x

    ! 3) Read data
    ! 3.1) Open the input file
    ! 3.3) Read array values
    ! 3.4) Close file

    Call zeit_ci('I/O ovrhd')

! Write Data
!   Open output file
!   ----------------
    unit = luavail()
    Call opnieee(unit, Trim(file_to_psas), 'new', ier)
       If (ier /= 0) Call die(myname_,'opnieee("'//Trim(file_to_psas)//'")', ier)

    print *, 'GetAI: ninc, nobs = ', ninc, nobs
    print *, 'GetAI:  alat =', minval(alat), maxval(alat)
    print *, 'GetAI:  alon =', minval(alon), maxval(alon)
    print *, 'GetAI:   lat =', minval(lat), maxval(lat)
    print *, 'GetAI:   lon =', minval(lon), maxval(lon)
    print *, 'GetAI:  tlev =', minval(tlev), maxval(tlev)
    print *, 'GetAI:  ttim =', minval(ttim), maxval(ttim)
    print *, 'GetAI:    kx =', minval(kx), maxval(kx)
    print *, 'GetAI:    ks =', minval(ks), maxval(ks)
    print *, 'GetAI:   tkt =', minval(tkt), maxval(tkt)
    print *, 'GetAI: txmso =', minval(txmso), maxval(txmso)
    print *, 'GetAI:  tomf =', minval(tomf), maxval(tomf)

!   Write size information
!   ----------------------
    Write(unit) ninc, nobs

!   Write array values
!   ------------------
    Write(unit) alat, alon, alev, akts
    Write(unit) lat, lon, tlev, ttim, kx, ks, tkt, txmso, tomf

!   Close file
!   ----------
    Call clsieee(unit, ier)
       If (ier /= 0) Call die(myname_,'clsieee("'//Trim(file_to_psas)//'")', ier)
    
!   Determine number of PEs and location of solve.x
!   -----------------------------------------------
    call getenv ( 'PSAS_NUM_MPI',cnpes)
    read(cnpes,'(i6)',iostat=ier) npes
    if ( ier /= 0 ) call die(myname_,&
         'env variable PSAS_NUM_MPI invalid or not found') 
    if(npes < 1) Call die(myname_,'invalid PSAS_NUM_MPI: ', npes)

    mpirun = ' '
    call getenv('MPIRUN',mpirun)             ! Unix binding
    !if(mpirun.eq.' ') Call die(myname_,'env var MPIRUN not set', 99)

!   Invoke 'solve.x'
!   ----------------
    call ExecResolv ( 'solve.x', solve_dot_x ) ! find solve.x in path
    write(command,'(3a)') trim(mpirun), ' ', trim(solve_dot_x)
    write(6,'(3a)') myname_, ': ', trim(command)
    call flush(6)
    ier = System ( trim(command) )
    if ( ier /= 0 ) then
       call die ( myname_, 'cannot run ' // trim(command) )
    end if

! Read data
!   Open the input file
!   -------------------
    unit = luavail()
    Call opnieee(unit, Trim(file_from_psas), 'old', ier)
       If (ier /= 0) Call die(myname_,'opnieee("'//Trim(file_from_psas)//'")', ier)

!   Read array values
!   -----------------
    Read(unit) ainc, xvec

    print *, 'GetAI:  ainc =', minval(ainc), maxval(ainc)

!   Close file
!   ----------
    Call clsieee(unit, ier, STATUS='DELETE')
       If (ier /= 0) Call die(myname_,'clsieee("'//Trim(file_from_psas)//'")', ier)

! Done
! ----
    Call zeit_co('I/O ovrhd')
    Return

  End Subroutine GetAI_External

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Moist_Sigf_ --- Generate sigf for pseudo rh
!
! !INTERFACE:

    subroutine Moist_Sigf_ ( i_am_root )

! !USES:

   use  m_sigFi_lookups, only : sigFi_store ! Sigf redefinition for PSAS
   use  m_sigFi_lookups, only : mixre_defined
   use  m_sigFi_lookups, only : MIXRE

    Implicit NONE

!
! !INPUT PARAMETERS:
!
   logical, intent(in) :: i_am_root
!
! !OUTPUT PARAMETERS:
!

!
! !DESCRIPTION: Creates 3D gridded pseudo-rh fcst error stdv table
!               and passes it to PSAS
!
! !REVISION HISTORY:
!
!  27sep2000  Dee        First version (hardwired table)
!  07nov2000  Dee        Check if already defined
!  21nov2000  Dee        Use get_sigf_W
!  03mar2003  Todling    Rename get_sigf_W to get_sigf_W_
!
!EOP
!-------------------------------------------------------------------------

    character(len=*), parameter ::  myname_ = myname//':Moist_Sigf_'

!   3D array for sigf (for now, depends on pressure only)
!   -----------------

    real :: sigf(2,2,10)  ! need only four lon/lat pairs per level

!   table definition:
!   ----------------
    real, parameter :: p(10) = (/ 1040.1, 1000.0, 925.0, 850.0, 700.0, &
                                   500.0,  400.0, 300.0, 100.0, 0.01 /)

    integer :: i, j, k
    real    :: s(10)

!   Check if already defined
!   ------------------------
    if ( mixre_defined ) return

!   Get fcst error stdv
!   -------------------
    call get_sigf_W_ ( s, p, 10 )
    do i = 1, 2
       do j = 1, 2
          do k = 1, 10
               sigf(i,j,k) = s(k)
          end do
       end do
    end do

!   tell PSAS to use values in sigf to override resource
!   ----------------------------------------------------
    call sigFi_store ( MIXRE, 2, 2, 10, p, sigf )

    if(i_am_root) &
    write(*,'(a)') myname_//': Defined moisture fcst error std dev for PSAS'

  end subroutine Moist_Sigf_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: setAOD_sigF_ --- Generate sigf for all control analysis variables
!
! !INTERFACE:

    subroutine setAOD_sigF ( i_am_root )

! !USES:

   use  m_sigFi_lookups, only : sigFi_store ! Sigf redefinition for PSAS
   use  m_sigFi_lookups, only : mixre_defined, hghte_defined
   use  m_sigFi_lookups, only : MIXRE, HGHTE

    Implicit NONE

!
! !INPUT PARAMETERS:
!
   logical, intent(in) :: i_am_root
!
! !OUTPUT PARAMETERS:
!

!
! !DESCRIPTION: Creates 3D gridded fcst error stdv table
!               and passes it to PSAS
!
! !REVISION HISTORY:
!
!  17feb2010  da Silva  Adapted for AOD analysis
!
!EOP
!-------------------------------------------------------------------------

    character(len=*), parameter ::  myname_ = myname//':setAOD_sigF'

!   3D array for sigf (for now, depends on pressure only)
!   -----------------

    real :: sigf(2,2,4)  ! need only four lon/lat pairs per level

!   table definition:
!   ----------------
    real :: p(4) = (/ 870., 660., 550., 470. /)
    real :: s(4) = (/ 0.45, 0.45, 0.45, 0.45 /)

    integer :: i, j, k


!   Check if already defined
!   ------------------------
    if ( mixre_defined .or. hghte_defined ) return

!   Get fcst error stdv
!   -------------------
    do i = 1, 2
       do j = 1, 2
          do k = 1, 4
               sigf(i,j,k) = s(k)
          end do
       end do
    end do

!   tell PSAS to use values in sigf to override resource
!   ----------------------------------------------------
    call sigFi_store ( HGHTE, 2, 2, 4, p, sigf ) ! not needed for AOD analsis
    call sigFi_store ( MIXRE, 2, 2, 4, p, sigf )

    if(i_am_root) &
    write(*,*) &
      myname_//': Defined fcst error std dev for PSAS, sigF = ', &
      s(:)

  end subroutine SetAOD_sigF


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: get_sigf_W_ --- Get sigf for pseudo rh
!
! !INTERFACE:

    subroutine get_sigf_W_ ( sigf, p, n )

! !USES:

    Implicit NONE

!
! !INPUT PARAMETERS:
!
    integer, intent(in)        :: n
    real,    intent(in)        :: p(n)      !  pressure in hPa

!
! !OUTPUT PARAMETERS:
!
    real,    intent(out)       :: sigf(n)   !  fcst error std dev

!
! !DESCRIPTION: Returns sigf for pseudo rh (moisture analysis variable)
!
!
! !REVISION HISTORY:
!
!  20nov2000  Dee        First version (hardwired table)
!  19jan2001  Dee        Extended table to higher level
!  08feb2001  Dee        Adjusted sigf values
!  28mar2002  Dee        Make sigf zero at 70hPa and higher
!
!EOP
!-------------------------------------------------------------------------

    character(len=*), parameter ::  myname_ = myname//':get_sigf_W_'

!   table definition:
!   ----------------
    real, parameter :: ptab(9) = (/ 1000.0, 925.0, 850.0, 700.0, 500.0, 400.0, 300.0, 100.0, 70.0 /)
    real, parameter :: stab(9) = (/   0.13,  0.13,  0.14,  0.15,  0.16,  0.16,  0.14,  0.08, 0.00 /)

    integer         :: i, k, ntab


!   table lookup with linear-in-log(p) interpolation
!   ------------------------------------------------
    ntab = size(ptab)
    do i = 1, n
       sigf(i) = stab(ntab)
       if ( p(i) .gt. ptab(1)) then
            sigf(i) = stab(1)
            cycle
       end if
       do k = 1, ntab-1
          if ( ptab(k) .ge. p(i) .and. p(i) .gt. ptab(k+1) ) then
               sigf(i) = stab(k) + (stab(k+1)-stab(k)) *      &
                        (log(ptab(k)) - log(p(i))) / (log(ptab(k)) - log(ptab(k+1)))
               exit
          end if
       end do
    end do


  end subroutine get_sigf_W_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Ts_Sigf_ --- Generate sigf for $T_s$ analysis
!
! !INTERFACE:

    subroutine Ts_Sigf_ ( i_am_root )

   use  m_sigFi_lookups, only : sigFi_store ! Sigf redefinition for PSAS
   use  m_sigFi_lookups, only : mixre_defined
   use  m_sigFi_lookups, only : MIXRE

    Implicit NONE

! !INPUT PARAMETERS:

   logical, intent(in) :: i_am_root

!
! !DESCRIPTION: Creates 3D gridded pseudo-rh fcst error stdv table
!               and passes it to PSAS
!
! !REVISION HISTORY:
!
!  27sep2000  Dee        First version (hardwired table)
!  07nov2000  Dee        Check if already defined
!  21nov2000  Dee        Use get_sigf_W
!  17dec2001  da Silva   Adapted from Moist_sigf
!  11Sep2002  da Silva & Radakovich  Changed sigF based on sigO=0.11 and 
!                                    K=0.866
!  28Mar2003  da Silva   Brought weights down to K=0.75
!  20Apr2003  da Silva   Latitude dependent weights: K=0.7  |lat|<  60,
!                                                    K=0.5  |lat|>= 60.
!
!EOP
!-------------------------------------------------------------------------

    character(len=*), parameter ::  myname_ = myname//':Ts_Sigf_'


!   3D array for sigf (for now, depends on pressure only)
!   -----------------------------------------------------
    integer, parameter :: jdim = 7
    real :: sigf(2,jdim,2)  ! need only four lon/lat pairs per level

!   table definition:
!   ----------------
    real, parameter :: p(2) = (/ 1000.0, 850.0 /)
    real, parameter :: s1(2) = (/ 0.17,  0.17 /)  ! K=0.70
    real, parameter :: s2(2) = (/ 0.11,  0.11 /)  ! K=0.50

    integer :: i, j, k
    real :: dlat, ylat

!   Check if already defined
!   ------------------------
    if ( mixre_defined ) return

!   Get fcst error stdv
!   -------------------
    dlat = 180.0 / ( jdim - 1 )
    do i = 1, 2
       do j = 1, jdim
          ylat = -90.0 + (j-1) * dlat
          if ( abs(ylat) .lt. 60.0 ) then
               do k = 1, 2
                  sigf(i,j,k) = s1(k)
               end do
          else
               do k = 1, 2
                  sigf(i,j,k) = s2(k)
               end do
          end if
       end do
    end do

!   tell PSAS to use values in sigf to override resource
!   ----------------------------------------------------
    call sigFi_store ( MIXRE, 2, jdim, 2, p, sigf )

    if(i_am_root) &
    write(*,'(a)') myname_//': Defined Ts fcst error std dev for PSAS'

  end subroutine Ts_Sigf_

End Module m_GetAI
