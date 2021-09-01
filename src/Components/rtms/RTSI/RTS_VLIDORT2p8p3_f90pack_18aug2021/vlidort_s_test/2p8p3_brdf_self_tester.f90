
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p3                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Authors :     Robert. J. D. Spurr (1)                      #
! #                Matt Christi                                 #
! #                                                             #
! #  Address (1) : RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VLIDORT_2p8p3                             #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released        August 2014                 #
! #      2.8   F90, released        May    2017                 #
! #      2.8.1 F90, released        August 2019                 # 
! #      2.8.2 F90, limited release May    2020                 # 
! #                                                             #
! #  Features Summary of Recent VLIDORT Versions:               #
! #  -------------------------------------------                #
! #                                                             #
! #      NEW: TOTAL COLUMN JACOBIANS         (2.4)              #
! #      NEW: BPDF Land-surface KERNELS      (2.4R)             #
! #      NEW: Thermal Emission Treatment     (2.4RT)            #
! #      Consolidated BRDF treatment         (2.4RTC)           #
! #      f77/f90 Release                     (2.5)              #
! #      External SS / New I/O Structures    (2.6)              #
! #                                                             #
! #      SURFACE-LEAVING / BRDF-SCALING      (2.7)              #
! #      TAYLOR Series / OMP THREADSAFE      (2.7)              #
! #      New Water-Leaving Treatment         (2.8)              #
! #      LBBF & BRDF-Telescoping, enabled    (2.8)              #
! #      Several Performance Enhancements    (2.8)              #
! #      Water-leaving coupled code          (2.8.1)            #
! #      Planetary problem, media properties (2.8.1)            #
! #      Doublet geometry post-processing    (2.8.2)            #
! #      Reduction zeroing, dynamic memory   (2.8.2)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.3, released 31 March 2021.                            #
! #     ==> Green's function RT solutions (Nstokes = 1 or 3)    #
! #     ==> Sphericity Corrections using MS source terms        #
! #     ==> BRDF upgrades, including new snow reflectance       #
! #     ==> SLEAVE Upgrades, extended water-leaving treatment   #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.3 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2021.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p3 ( Version 2.8.3 )            #
! #                                                                 #
! # VLIDORT_2p8p3 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p3 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p3  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

      PROGRAM BRDF_Self_Tester

!  Pars file

      USE VLIDORT_PARS_m

!  Module files for VBRDF

      USE VBRDF_SUP_AUX_m, Only : VBRDF_READ_ERROR
      USE VBRDF_SUP_MOD_m
      USE VBRDF_LINSUP_MOD_m

!  Implicit none

      IMPLICIT NONE

!  VBRDF supplement file inputs status structure

      TYPE(VBRDF_Input_Exception_Handling)   :: VBRDF_Sup_InputStatus

!  VBRDF supplement input structure

      TYPE(VBRDF_Sup_Inputs)                 :: VBRDF_Sup_In

!  VBRDF supplement output structure

      TYPE(VBRDF_Sup_Outputs)                :: VBRDF_Sup_Out, VBRDF_Sup_OutP

      TYPE(VBRDF_Output_Exception_Handling)  :: VBRDF_Sup_OutputStatus

!  VBRDF supplement linearized input structure

      TYPE(VBRDF_LinSup_Inputs)              :: VBRDF_LinSup_In

!  VBRDF supplement linearized output structure

      TYPE(VBRDF_LinSup_Outputs)             :: VBRDF_LinSup_Out

!  Local Variables
!  ===============

!  Help variables

      INTEGER            :: BRDF_IDX, I, IA, IB, J, M, MQ(16), NS, NSSQ, NTASKS, &
                            Q, Q1, Q2, T, TEST, UM, W
      INTEGER            :: NBEAMS, NSTREAMS, N_USER_STREAMS, N_USER_RELAZMS, NUMGEO

      INTEGER            :: FIRST_BRDF_IDX, LAST_BRDF_IDX

      REAL(fpk)          :: DX, EPS, EPSFAC, ANJac(MAXSTOKES), FDJac(MAXSTOKES)

      CHARACTER(LEN=1)   :: C1,CNS
      CHARACTER(LEN=2)   :: CBRDF
      CHARACTER(LEN=40)  :: FORM1,FORM2
      CHARACTER(LEN=50)  :: VBRDF_INPUT_FILE
      CHARACTER(LEN=120) :: TRACE

      INTEGER            :: BS_NMOMENTS_INPUT
      LOGICAL            :: DO_DEBUG_RESTORATION, GCM
      LOGICAL            :: DO_SOLAR_SOURCES, DO_USER_STREAMS


      write(*,'(/1x,a)') 'Starting BRDF_Self_Tester'

!  BASELINE CALCULATION
!  ********************

!  Start VBRDF test loop

      !General set
      FIRST_BRDF_IDX = 1
      LAST_BRDF_IDX  = MAXBRDF_IDX

      !Specific subsets
      !FIRST_BRDF_IDX = 19
      !LAST_BRDF_IDX  = 19

      DO BRDF_IDX = FIRST_BRDF_IDX, LAST_BRDF_IDX

      !Cox-Munk type only set
      !DO TEST = 1,5

      !IF (TEST == 1) THEN
      !  BRDF_IDX = 9
      !ELSEIF (TEST == 2) THEN
      !  BRDF_IDX = 10
      !ELSEIF (TEST == 3) THEN
      !  BRDF_IDX = 11
      !ELSEIF (TEST == 4) THEN
      !  BRDF_IDX = 15
      !ELSEIF (TEST == 5) THEN
      !  BRDF_IDX = 16
      !ENDIF

      !BPDF/Mod Fresnel type only set
      !DO TEST = 1,4

      !IF (TEST == 1) THEN
      !  BRDF_IDX = 12
      !ELSEIF (TEST == 2) THEN
      !  BRDF_IDX = 13
      !ELSEIF (TEST == 3) THEN
      !  BRDF_IDX = 14
      !ELSEIF (TEST == 4) THEN
      !  BRDF_IDX = 18
      !ENDIF

!  VBRDF SUPPLEMENT CALCULATION
!  ============================

!  Get the VBRDF inputs

      VBRDF_INPUT_FILE = 'vlidort_s_test/VBRDF_ReadInput_self_test.cfg'

      CALL VBRDF_LIN_INPUTMASTER ( VBRDF_INPUT_FILE, & ! Input
        VBRDF_Sup_In,         & ! Outputs
        VBRDF_LinSup_In,      & ! Outputs
        VBRDF_Sup_InputStatus ) ! Outputs

      IF ( VBRDF_Sup_InputStatus%BS_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VBRDF_READ_ERROR ( 'VBRDF_ReadInput.log', VBRDF_Sup_InputStatus )

!  Modify some VBRDF config file inputs based on the current VBRDF being tested

      CALL GET_VBRDF_PLUS(BRDF_IDX,VBRDF_Sup_In,VBRDF_LinSup_In)

      write(*,'(/1x,a)')   '************************'
      write(*,'(1x,2a)')   'Doing VBRDF ',trim(adjustl(VBRDF_Sup_In%BS_BRDF_NAMES(1)))
      write(*,'(1x,a,i2)') 'BRDF index = ',BRDF_IDX

!  Do not want debug restoration

      DO_DEBUG_RESTORATION = .false.

!  A normal calculation will require

      BS_NMOMENTS_INPUT = 2 * VBRDF_Sup_In%BS_NSTREAMS - 1

!  Linearized VBRDF call
!    The output will now be used for the Main VLIDORT calculation

      write(*,'(1x,a,i1)') 'Doing baseline test'

      CALL VBRDF_LIN_MAINMASTER ( &
        DO_DEBUG_RESTORATION,  & ! Inputs
        BS_NMOMENTS_INPUT,     & ! Inputs
        VBRDF_Sup_In,          & ! Inputs
        VBRDF_LinSup_In,       & ! Inputs
        VBRDF_Sup_Out,         & ! Outputs
        VBRDF_LinSup_Out,      & ! Outputs
        VBRDF_Sup_OutputStatus ) ! Output Status

      if ( VBRDF_Sup_OutputStatus%BS_STATUS_OUTPUT .ne. vlidort_success ) THEN
         TRACE = 'VBRDF_LIN_MAINMASTER failed' ; go to 678
      endif

!  mick note 9/19/2017 - adjusted indices on some of the BRDF output arrays to promote
!                        consistency in their use
!                      - added this table for clarity
!
!                        Index     Task for BRDF output
!                        ----- ----------------------------
!                         ib   incident solar zenith angle
!                         i    incident quadrature angle
!                         j    reflected quadrature angle
!                         um   reflected user zenith angle
!                         ia   reflected user azimuth angle

!  Write normal VBRDF supplement results (as a check)

      GCM = (VBRDF_Sup_In%BS_which_brdf(1) .eq. 10 .or. &
             VBRDF_Sup_In%BS_which_brdf(2) .eq. 10 .or. &
             VBRDF_Sup_In%BS_which_brdf(3) .eq. 10)

!     if ( GCM ) then
!        NSSQ = VBRDF_Sup_In%BS_NSTOKES * VBRDF_Sup_In%BS_NSTOKES
!     else
!        NSSQ = 1
!     endif

      ns   = VBRDF_Sup_In%BS_NSTOKES
      nssq = ns * ns
      write(cns,'(i1)') ns

      nbeams           = VBRDF_Sup_In%BS_nbeams
      nstreams         = VBRDF_Sup_In%BS_nstreams
      n_user_streams   = VBRDF_Sup_In%BS_n_user_streams
      n_user_relazms   = VBRDF_Sup_In%BS_n_user_relazms

      DO_USER_STREAMS  = VBRDF_Sup_In%BS_DO_USER_STREAMS
      DO_SOLAR_SOURCES = VBRDF_Sup_In%BS_DO_SOLAR_SOURCES

!  Define output mask (for baseline files)

      MQ = 0
      IF ( ns.eq.1 ) then
        MQ(1) = 1
      ELSE IF ( ns.eq.2 ) then
        MQ(1) = 1 ; MQ(2) = 2
        MQ(3) = 5 ; MQ(4) = 6
      ELSE IF ( ns.eq.3 ) then
        MQ(1) = 1 ; MQ(2) = 2  ; MQ(3) = 3
        MQ(4) = 5 ; MQ(5) = 6  ; MQ(6) = 7
        MQ(7) = 9 ; MQ(8) = 10 ; MQ(9) = 11
      ELSE IF ( ns.eq.4 ) then
        DO q = 1, nssq
           MQ(q) = q
        ENDDO
      ENDIF

!  Write BASELINE regular results
!  ==============================

      IF (BRDF_IDX == FIRST_BRDF_IDX) THEN
        !write(cbrdf,'(i2.2)') brdf_idx
        !open(36, file = 'vlidort_s_test/results_brdfplus_selftest.brdf' // cbrdf // &
        !                '_ref_ns' // cns // '_v3' , &
        !         status = 'unknown')
        open(36, file = 'vlidort_s_test/results_brdf_self_tester.res',&
                 status = 'unknown')

        write(36,'(/a/)') ' Number of Beams/Quad-streams/User-streams/' // &
                          'User-relazms/nstokessq/GCM flag'
        write(36,'(5i5,2x,L2)') &
          nbeams,         nstreams, &
          n_user_streams, n_user_relazms, &
          NSSQ, GCM
      ENDIF

      write(36,'(/1x,a)')         '****************************************'
      write(36,'(1x,2a,3x,a,i2)') 'VBRDF Name: ',trim(adjustl(VBRDF_Sup_In%BS_BRDF_NAMES(1))),&
                                  'VBRDF Index: ',BRDF_IDX

      IF ( DO_SOLAR_SOURCES .and. DO_USER_STREAMS ) THEN
        write(36,'(/T12,a//T3,a,2x,a/)') 'Exact BRDF, all geometries',&
                                         'sza vza azm','Exact BRDF values'
        numgeo = nbeams*n_user_streams*n_user_relazms
        do ib = 1, nbeams
          do um = 1, n_user_streams
            do ia = 1, n_user_relazms
              do q = 1, ns
                q2 = q*ns ; q1 = q2 - ns + 1
                if ( q .eq. 1) then
                  write(36,'(3i4,1p4e16.7)') ib,um,ia,VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(mq(q1):mq(q2),um,ia,ib)
                else
                  write(36,'(12x,1p4e16.7)') VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(mq(q1):mq(q2),um,ia,ib)
                end if
              enddo
              if ((ns > 1) .and. (ib*um*ia < numgeo)) write(36,*)
            enddo
          enddo
        enddo
        write(36,*)
      ENDIF

      IF ( .not.VBRDF_Sup_In%BS_DO_DIRECTBOUNCE_ONLY ) THEN

      DO M = 0, BS_NMOMENTS_INPUT

        write(36,'(/T12,a,i3//T17,a,2x,a/)') 'Fourier component # ',M,&
                                             'str str','BRDF values'
        numgeo = nstreams**2
        do i = 1, nstreams
          do j = 1, nstreams
            do q = 1, ns
              q2 = q*ns ; q1 = q2 - ns + 1
              if ( q .eq. 1) then
                write(36,'(1x,a,2i4,1p4e16.7)')'Stream-Stream',i,j,VBRDF_Sup_Out%BS_BRDF_F(M,mq(q1):mq(q2),j,i)
              else
                write(36,'(22x,1p4e16.7)') VBRDF_Sup_Out%BS_BRDF_F(M,mq(q1):mq(q2),j,i)
              end if
            enddo
            if ((ns > 1) .and. (i*j < numgeo)) write(36,*)
          enddo
        enddo

        IF ( DO_SOLAR_SOURCES ) THEN
          write(36,'(/T17,a,2x,a/)') 'sza str','BRDF values'
          numgeo = nbeams*nstreams
          do ib = 1, nbeams
            do j = 1, nstreams
              do q = 1, ns
                q2 = q*ns ; q1 = q2 - ns + 1
                if ( q .eq. 1) then
                  write(36,'(1x,a,2i4,1p4e16.7)')'Solars-Stream',ib,j,VBRDF_Sup_Out%BS_BRDF_F_0(M,mq(q1):mq(q2),j,ib)
                else
                  write(36,'(22x,1p4e16.7)') VBRDF_Sup_Out%BS_BRDF_F_0(M,mq(q1):mq(q2),j,ib)
                end if
              enddo
              if ((ns > 1) .and. (ib*j < numgeo)) write(36,*)
            enddo
          enddo
        ENDIF

        IF ( DO_USER_STREAMS ) THEN
          write(36,'(/T17,a,2x,a/)') 'str vza','BRDF values'
          numgeo = nstreams*n_user_streams
          do i = 1, nstreams
            do um = 1, n_user_streams
              do q = 1, ns
                q2 = q*ns ; q1 = q2 - ns + 1
                if ( q .eq. 1) then
                  write(36,'(1x,a,2i4,1p4e16.7)')'Stream-Users ',i,um,VBRDF_Sup_Out%BS_USER_BRDF_F(M,mq(q1):mq(q2),um,i)
                else
                  write(36,'(22x,1p4e16.7)') VBRDF_Sup_Out%BS_USER_BRDF_F(M,mq(q1):mq(q2),um,i)
                end if
              enddo
              if ((ns > 1) .and. (i*um < numgeo)) write(36,*)
            enddo
          enddo
        ENDIF

        IF ( DO_SOLAR_SOURCES .and. DO_USER_STREAMS ) THEN
          write(36,'(/T17,a,2x,a/)') 'sza vza','BRDF values'
          numgeo = nbeams*n_user_streams
          do ib = 1, nbeams
            do um = 1, n_user_streams
              do q = 1, ns
                q2 = q*ns ; q1 = q2 - ns + 1
                if ( q .eq. 1) then
                  write(36,'(1x,a,2i4,1p4e16.7)')'Solars-Users ',ib,um,VBRDF_Sup_Out%BS_USER_BRDF_F_0(M,mq(q1):mq(q2),um,ib)
                else
                  write(36,'(22x,1p4e16.7)') VBRDF_Sup_Out%BS_USER_BRDF_F_0(M,mq(q1):mq(q2),um,ib)
                end if
              enddo
              if ((ns > 1) .and. (ib*um < numgeo)) write(36,*)
            enddo
          enddo
        ENDIF

        write(36,*)

!  End Fourier componenet loop

      ENDDO

      ENDIF

!  Close regular reflectance file when all VBRDFs finished

      IF (BRDF_IDX == LAST_BRDF_IDX) close(36)

!  FINITE DIFFERENCE SECTION
!  *************************

!  Define perturbation for FD computations

      eps = 1.0d-04
      epsfac = 1.0d0 + eps

!  Open linearized reflectance file for current VBRDF & write metadata

      IF (BRDF_IDX == FIRST_BRDF_IDX) THEN
        !open(37, file = 'vlidort_s_test/results_brdfplus_selftest.brdf' // cbrdf // &
        !                '_wfs_ns' // cns // '_v3' , &
        !         status = 'unknown')
        open(37, file = 'vlidort_s_test/results_brdf_self_tester.wfs',&
                 status = 'unknown')

        write(37,'(/a/)') ' Number of SurfaceWFs/Beams/Quad-streams/' // &
                          'User-streams/User-relazms/nstokessq/GCM flag'
        write(37,'(6i5,2x,l2)') &
          VBRDF_LinSup_In%BS_n_surface_wfs, nbeams, &
          nstreams,         n_user_streams, &
          n_user_relazms,   NSSQ, GCM
      ENDIF

      write(37,'(/1x,a)')         '****************************************'
      write(37,'(1x,2a,3x,a,i2)') 'VBRDF Name: ',trim(adjustl(VBRDF_Sup_In%BS_BRDF_NAMES(1))),&
                                  'VBRDF Index: ',BRDF_IDX

!  Number of perturbation tasks

      IF ( BRDF_IDX .eq. 10 ) THEN
        !Giss Cox-Munk Surface (Vector):
        !  can only do analytic derivative wrt brdf factor & wind speed (par #1)
        ntasks = 2
      ELSEIF ( BRDF_IDX .eq. 11 ) THEN
        !Giss Cox-Munk CRI Surface (Vector):
        !  there are currently no analytic derivatives built for this kernel
        ntasks = 0
      ELSEIF ( BRDF_IDX .eq. 15 .or. BRDF_IDX .eq. 16 ) THEN
        !For New Cox-Munk Surface (Scalar) or New Giss Cox-Munk Surface (Vector):
        !  can only do analytic derivative wrt wind speed (par #1)
        ntasks = 1
      ELSEIF ( BRDF_IDX .eq. 19 ) THEN
        !Snow Surface (Scalar):
        !  can only do analytic derivative wrt brdf factor and then diameter of
        !  snow-grains and pollution contamination (par #1 & par #2)
        ntasks = 3
      ELSE
        !All others
        ntasks = 1 + VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)
      ENDIF

      IF (ntasks .gt. 0) THEN

!  Start task loop

      DO t = 1, ntasks

!  Progress

        write(*,'(1x,a,i1)') 'Doing FD test # ',t
        write(c1,'(i1)') t
        write(37,'(/1x,a,i1/1x,a)') 'For surface Jacobian # ',t,&
                                    '------------------------'

!  Read the VBRDF information again

        VBRDF_INPUT_FILE = 'vlidort_s_test/VBRDF_ReadInput_self_test.cfg'

        CALL VBRDF_INPUTMASTER ( VBRDF_INPUT_FILE, &
          VBRDF_Sup_In, &
          VBRDF_Sup_InputStatus )

        IF ( VBRDF_Sup_InputStatus%BS_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
          CALL VBRDF_READ_ERROR ( 'VBRDF_ReadInput.log', VBRDF_Sup_InputStatus )

!  Modify some VBRDF config file inputs based on the current VBRDF being tested

        CALL GET_VBRDF(BRDF_IDX,VBRDF_Sup_In)

!  Do perturbation for 1-kernel VBRDF

        IF ( BRDF_IDX .eq. 15 .or. BRDF_IDX .eq. 16 ) THEN
          !For New Cox-Munk Surface (Scalar) or New Giss Cox-Munk Surface (Vector):
          !  can only do analytic derivative wrt wind speed (par #1)
          IF ( t .eq. 1 ) THEN
            DX = VBRDF_Sup_In%BS_WINDSPEED * eps
            VBRDF_Sup_In%BS_WINDSPEED = VBRDF_Sup_In%BS_WINDSPEED * epsfac
          ENDIF
        ELSE
          !All others
          IF ( t .eq. 1 ) THEN
            DX = VBRDF_Sup_In%BS_BRDF_FACTORS(1) * eps
            VBRDF_Sup_In%BS_BRDF_FACTORS(1) = VBRDF_Sup_In%BS_BRDF_FACTORS(1) * epsfac
          ELSEIF ( t .ge. 2 .and. t .le. ntasks ) THEN
            DX = VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,t-1) * eps
            VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,t-1) = VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,t-1) * epsfac
          ENDIF

          !Special adjustment for snow kernel - have to scale perturbation of M parameter
          !  (par #2) by 1.0e-8 for computing FD Jacobian (task) #3
          IF ( ( BRDF_IDX .eq. 19 ) .and. ( t .eq. 3 ) ) DX = DX * 1.0e-8_fpk
        ENDIF

!  VBRDF call

        CALL VBRDF_MAINMASTER ( &
          DO_DEBUG_RESTORATION,  & ! Inputs
          BS_NMOMENTS_INPUT,     & ! Inputs
          VBRDF_Sup_In,          & ! Inputs
          VBRDF_Sup_OutP,        & ! Outputs
          VBRDF_Sup_OutputStatus ) ! Output Status

        IF ( VBRDF_Sup_OutputStatus%BS_STATUS_OUTPUT .ne. vlidort_success ) THEN
           TRACE = 'VBRDF_MAINMASTER failed, perturbation # '//c1 ; go to 678
        ENDIF

!  Write BASELINE analytic and FD linearized results
!  =================================================

        IF ( DO_SOLAR_SOURCES .and. DO_USER_STREAMS ) THEN

          form1 = '(3i4,1p'//cns//'e16.7,2x,1p'//cns//'e16.7)'
          form2 = '(12x,1p'//cns//'e16.7,2x,1p'//cns//'e16.7)'

          write(37,'(/T12,a//T3,a,2x,a/)') 'Exact BRDF, all geometries',&
                                           'sza vza azm','Exact BRDF values'
          numgeo = nbeams*n_user_streams*n_user_relazms
          do ib = 1, nbeams
            do um = 1, n_user_streams
              do ia = 1, n_user_relazms
                do q = 1, ns
                  q2 = q*ns ; q1 = q2 - ns + 1
                  ANJac(1:ns) = VBRDF_LinSup_Out%BS_LS_DBOUNCE_BRDFUNC(t,mq(q1):mq(q2),um,ia,ib)
                  FDJac(1:ns) = ( VBRDF_Sup_OutP%BS_DBOUNCE_BRDFUNC(mq(q1):mq(q2),um,ia,ib) &
                                 - VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(mq(q1):mq(q2),um,ia,ib) ) / DX
                  if ( q .eq. 1) then
                    write(37,form1) ib,um,ia,ANJac(1:ns),FDJac(1:ns)
                  else
                    write(37,form2) ANJac(1:ns),FDJac(1:ns)
                  end if
                enddo
                if ((ns > 1) .and. (ib*um*ia < numgeo)) write(37,*)
              enddo
            enddo
          enddo
          write(37,*)
        ENDIF

        IF ( .not.VBRDF_Sup_In%BS_DO_DIRECTBOUNCE_ONLY ) THEN

        form1 = '(1x,a,2i4,1p'//cns//'e16.7,2x,1p'//cns//'e16.7)'
        form2 =      '(22x,1p'//cns//'e16.7,2x,1p'//cns//'e16.7)'

        DO M = 0, BS_NMOMENTS_INPUT

          write(37,'(/T12,a,i3//T16,a,2x,a/)') 'Fourier component # ',M,&
                                               'str str','BRDF values'
          numgeo = nstreams**2
          do i = 1, nstreams
            do j = 1, nstreams
              do q = 1, ns
                q2 = q*ns ; q1 = q2 - ns + 1
                ANJac(1:ns) = VBRDF_LinSup_Out%BS_LS_BRDF_F(t,M,mq(q1):mq(q2),j,i)
                FDJac(1:ns) = ( VBRDF_Sup_OutP%BS_BRDF_F(M,mq(q1):mq(q2),j,i) &
                               - VBRDF_Sup_Out%BS_BRDF_F(M,mq(q1):mq(q2),j,i) ) / DX
                if ( q .eq. 1) then
                  write(37,form1)'Stream-Stream',i,j,ANJac(1:ns),FDJac(1:ns)
                else
                  write(37,form2) ANJac(1:ns),FDJac(1:ns)
                end if
              enddo
              if ((ns > 1) .and. (i*j < numgeo)) write(37,*)
            enddo
          enddo

          IF ( DO_SOLAR_SOURCES ) THEN
            write(37,'(/T16,a,2x,a/)') 'sza str','BRDF values'
            numgeo = nbeams*nstreams
            do ib = 1, nbeams
              do j = 1, nstreams
                do q = 1, ns
                  q2 = q*ns ; q1 = q2 - ns + 1
                  ANJac(1:ns) = VBRDF_LinSup_Out%BS_LS_BRDF_F_0(t,M,mq(q1):mq(q2),j,ib)
                  FDJac(1:ns) = ( VBRDF_Sup_OutP%BS_BRDF_F_0(M,mq(q1):mq(q2),j,ib) &
                                 - VBRDF_Sup_Out%BS_BRDF_F_0(M,mq(q1):mq(q2),j,ib) ) / DX
                  if ( q .eq. 1) then
                    write(37,form1)'Solars-Stream',ib,j,ANJac(1:ns),FDJac(1:ns)
                  else
                    write(37,form2) ANJac(1:ns),FDJac(1:ns)
                  end if
                enddo
                if ((ns > 1) .and. (ib*j < numgeo)) write(37,*)
              enddo
            enddo
          ENDIF

          IF ( DO_USER_STREAMS ) THEN
            write(37,'(/T16,a,2x,a/)') 'str vza','BRDF values'
            numgeo = nstreams*n_user_streams
            do i = 1, nstreams
              do um = 1, n_user_streams
                do q = 1, ns
                  q2 = q*ns ; q1 = q2 - ns + 1
                  ANJac(1:ns) = VBRDF_LinSup_Out%BS_LS_USER_BRDF_F(t,M,mq(q1):mq(q2),um,i)
                  FDJac(1:ns) = ( VBRDF_Sup_OutP%BS_USER_BRDF_F(M,mq(q1):mq(q2),um,i) &
                                 - VBRDF_Sup_Out%BS_USER_BRDF_F(M,mq(q1):mq(q2),um,i) ) / DX
                  if ( q .eq. 1) then
                    write(37,form1)'Stream-Users ',i,um,ANJac(1:ns),FDJac(1:ns)
                  else
                    write(37,form2) ANJac(1:ns),FDJac(1:ns)
                  end if
                enddo
                if ((ns > 1) .and. (i*um < numgeo)) write(37,*)
              enddo
            enddo
          ENDIF

          IF ( DO_SOLAR_SOURCES .and. DO_USER_STREAMS ) THEN
            write(37,'(/T16,a,2x,a/)') 'sza vza','BRDF values'
            numgeo = nbeams*n_user_streams
            do ib = 1, nbeams
              do um = 1, n_user_streams
                do q = 1, ns
                  q2 = q*ns ; q1 = q2 - ns + 1
                  ANJac(1:ns) = VBRDF_LinSup_Out%BS_LS_USER_BRDF_F_0(t,M,mq(q1):mq(q2),um,ib)
                  FDJac(1:ns) = ( VBRDF_Sup_OutP%BS_USER_BRDF_F_0(M,mq(q1):mq(q2),um,ib) &
                                 - VBRDF_Sup_Out%BS_USER_BRDF_F_0(M,mq(q1):mq(q2),um,ib) ) / DX
                  if ( q .eq. 1) then
                    write(37,form1)'Solars-Users ',ib,um,ANJac(1:ns),FDJac(1:ns)
                  else
                    write(37,form2) ANJac(1:ns),FDJac(1:ns)
                  end if
                enddo
                if ((ns > 1) .and. (ib*um < numgeo)) write(37,*)
              enddo
            enddo
          ENDIF

          write(37,*)

!  End Fourier component loop

        ENDDO

        ENDIF

!  End task loop

      ENDDO

      ELSE

        write(37,'(/1x,a/)') '(No Surface Jacobians for this BRDF kernel)'

      ENDIF

!  Close linearized reflectance file when all VBRDFs finished

      IF (BRDF_IDX == LAST_BRDF_IDX) close(37)

!  End VBRDF test loop

      ENDDO

      write(*,'(/1x,a)')  '************************'
      write(*,'(/1x,a/)') 'BRDF_Self_Tester done'

!  Normal stop

      stop

!  VBRDF error finish

678   continue

      write(*,*)
      write(*,'(1x,a)') trim(TRACE)
      write(*,*)'Number of error messages from VBRDF calculation = ',&
                 VBRDF_Sup_OutputStatus%BS_NOUTPUTMESSAGES
      do i = 1, VBRDF_Sup_OutputStatus%BS_NOUTPUTMESSAGES
         write(*,'(a,i2,a,a)')' Message # ', i, ': ',&
           adjustl(trim(VBRDF_Sup_OutputStatus%BS_OUTPUTMESSAGES(i)))
      enddo
      write(*,*)

      stop
      END PROGRAM BRDF_Self_Tester

!*******************************************************************************
      SUBROUTINE GET_VBRDF(BRDF_IDX,VBRDF_Sup_In)

      USE VLIDORT_PARS_m

      USE VBRDF_SUP_KERNELS_m, Only: VBRDF_Water_RefracIndex
      USE VBRDF_SUP_MOD_m

      IMPLICIT NONE

!  Inputs

      INTEGER, INTENT(IN) :: BRDF_IDX

!  InOut

      TYPE(VBRDF_Sup_Inputs), INTENT(INOUT) :: VBRDF_Sup_In

!  Local variables

      LOGICAL   :: Do_Shadow_Effect
      REAL(fpk) :: Refrac_R, Refrac_I, Salinity, SigmaSq, SurfAlb, &
                   Wavelength, WindDir, WindSpeed

!  Max VBRDF index check

      if (MAXBRDF_IDX > 19) then
        write(*,*)
        write(*,*) 'Error: subroutine GET_VBRDF in BRDF_Self_Tester currently set'
        write(*,*) 'up for VBRDFs 1-19.  If necessary, add additional VBRDFs to GET_VBRDF'
        stop
      endif

!  Misc inputs

      !For Lambertian surface
      SurfAlb = 0.1_fpk      

      !For water-related surfaces
      WindSpeed  = 15.0_fpk  !in m/s
      WindDir    = 180.0_fpk   !in degrees
      !WindDir    = 70.0_fpk  !in degrees
      SigmaSq    = 0.003_fpk + 0.00512_fpk*WindSpeed !sigma**2

      Wavelength = 0.5_fpk !in um   --> Refrac_I on order of 1.0e-9
      !Wavelength = 3.0_fpk !in um   --> Refrac_I on order of 1.0e-1
      Salinity   = 34.3_fpk  !in ppt

      !Do_Shadow_Effect = .FALSE.
      Do_Shadow_Effect = .TRUE.

      CALL VBRDF_Water_RefracIndex ( Wavelength, Salinity, Refrac_R, Refrac_I )
      Refrac_I = Zero
      !write(*,*) 'Refrac_R = ',Refrac_R
      !write(*,*) 'Refrac_I = ',Refrac_I

!  Define BRDF

      !BASIC SURFACE INPUTS
      VBRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG = .FALSE.

      VBRDF_Sup_In%BS_N_BRDF_KERNELS = 1
      VBRDF_Sup_In%BS_WHICH_BRDF(1)  = BRDF_IDX

      VBRDF_Sup_In%BS_BRDF_NAMES        = ''
      VBRDF_Sup_In%BS_BRDF_FACTORS      = 0.0_fpk
      VBRDF_Sup_In%BS_N_BRDF_PARAMETERS = 0
      VBRDF_Sup_In%BS_BRDF_PARAMETERS   = 0.0_fpk

      !BRDF INPUTS FOR A ...
      IF (BRDF_IDX == 1) THEN
        !... LAMBERTIAN SURFACE USING GENERAL BRDF CODE (ISOTROPIC: SCALAR)

        VBRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG = .TRUE.

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Lambertian'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = SurfAlb

      ELSE IF (BRDF_IDX == 2) THEN
        !... ROSS-THIN SURFACE (TREES WITH SMALL LEAVES OVER A FLAT SURFACE: SCALAR)

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Ross-thin '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk

      ELSE IF (BRDF_IDX == 3) THEN
        !... ROSS-THICK SURFACE (TREES WITH LARGE LEAVES OVER A FLAT SURFACE: SCALAR)

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Ross-thick'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk

      ELSE IF (BRDF_IDX == 4) THEN
        !... LI-SPARSE SURFACE (SPARSE TREES: SCALAR)

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Li-sparse '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 2
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = 1.0_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = 2.0_fpk

      ELSE IF (BRDF_IDX == 5) THEN
        !... LI-DENSE SURFACE (DENSE TREES: SCALAR)

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Li-dense  '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 2
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = 2.5_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = 2.0_fpk

      ELSE IF (BRDF_IDX == 6) THEN
        !... HAPKE SURFACE (???: SCALAR)

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Hapke     '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 3
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = 0.6_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = 0.06_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,3)      = 1.0_fpk

      ELSE IF (BRDF_IDX == 7) THEN
        !... ROUJEAN SURFACE (RANDOM RECTANGULAR BLOCKS: SCALAR)

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Roujean   '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk

      ELSE IF (BRDF_IDX == 8) THEN
        !... RAHMAN SURFACE (???: SCALAR)

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Rahman    '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 3
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      =  0.1_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = -0.5_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,3)      =  0.5_fpk
 
      ELSE IF (BRDF_IDX == 9) THEN          
        !... COX-MUNK SURFACE (SUNGLINT ON WATER: SCALAR)
 
        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Cox-Munk  '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 1.0_fpk !0.1_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 2
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = SigmaSq
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = Refrac_R**2

        !Note: Shadowing here determined indirectly by the "DO_SHADOW_EFFECT"
        !      flag through the 3rd parameter
        VBRDF_Sup_In%BS_DO_SHADOW_EFFECT = Do_Shadow_Effect

      ELSE IF (BRDF_IDX == 10) THEN
        !... GISS COX-MUNK SURFACE (SUNGLINT ON WATER: VECTOR)
        !    --> Note difference of parameters #1 & #2 between this kernel & kernel 9

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'GissCoxMnk'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 1.0_fpk !0.1_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 2
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = Half*SigmaSq 
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = Refrac_R

        !Note: Shadowing here determined indirectly by the "DO_SHADOW_EFFECT"
        !      flag through the 3rd parameter
        VBRDF_Sup_In%BS_DO_SHADOW_EFFECT = Do_Shadow_Effect

      ELSE IF (BRDF_IDX == 11) THEN
        !... GISS COX-MUNK COMPLEX REFRACTIVE INDEX ("CRI") SURFACE (SUNGLINT ON WATER: VECTOR)
        !    --> Note difference of parameters #1 & #2 between this kernel & kernel 9

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'GCMcomplex'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 1.0_fpk !0.1_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 3
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = Half*SigmaSq
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = Refrac_R
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,3)      = Refrac_I

        !Note: Shadowing here determined directly by the "DO_SHADOW_EFFECT" flag
        VBRDF_Sup_In%BS_DO_SHADOW_EFFECT = Do_Shadow_Effect

      ELSE IF (BRDF_IDX == 12) THEN
        !... BPDF-SOIL SURFACE (SOIL: VECTOR)
        
        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'BPDF-Soil '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 1.0_fpk !0.99_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 1
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = Refrac_R

      ELSE IF (BRDF_IDX == 13) THEN
        !... BPDF-VEGN SURFACE (VEGETATION: VECTOR)
        
        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'BPDF-Vegn '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 1
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = Refrac_R

      ELSE IF (BRDF_IDX == 14) THEN          
        !... BPDF-NDVI SURFACE (NDVI-PARAMETERIZED: VECTOR)
        
        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'BPDF-NDVI '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 3
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = Refrac_R
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = 0.2_fpk !NDVI
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,3)      = 1.0_fpk

      ELSE IF (BRDF_IDX == 15) THEN          
        !... NEW COX-MUNK SURFACE (SUNGLINT ON WATER: SCALAR)

        VBRDF_Sup_In%BS_DO_NewCMGLINT             = .TRUE.
        VBRDF_Sup_In%BS_DO_NewGCMGLINT            = .FALSE.

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'NewCMGlint'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 1.0_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 2
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = WindSpeed
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = Salinity

        VBRDF_Sup_In%BS_WINDSPEED                 = WindSpeed
        VBRDF_Sup_In%BS_WINDDIR                   = WindDir
        VBRDF_Sup_In%BS_WAVELENGTH                = Wavelength
        VBRDF_Sup_In%BS_SALINITY                  = Salinity

        VBRDF_Sup_In%BS_DO_GlintShadow            = Do_Shadow_Effect
        VBRDF_Sup_In%BS_DO_FoamOption             = .FALSE. !.TRUE.
        VBRDF_Sup_In%BS_DO_FacetIsotropy          = .TRUE. !.FALSE.

      ELSE IF (BRDF_IDX == 16) THEN          
        !... NEW GISS COX-MUNK SURFACE (SUNGLINT ON WATER: VECTOR)

        VBRDF_Sup_In%BS_DO_NewCMGLINT             = .FALSE.
        VBRDF_Sup_In%BS_DO_NewGCMGLINT            = .TRUE.

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'NewGCMGlit'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 1.0_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 2
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = WindSpeed
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = Salinity

        VBRDF_Sup_In%BS_WINDSPEED                 = WindSpeed
        VBRDF_Sup_In%BS_WINDDIR                   = WindDir
        VBRDF_Sup_In%BS_WAVELENGTH                = Wavelength
        VBRDF_Sup_In%BS_SALINITY                  = Salinity

        VBRDF_Sup_In%BS_DO_GlintShadow            = Do_Shadow_Effect
        VBRDF_Sup_In%BS_DO_FoamOption             = .FALSE. !.TRUE.
        VBRDF_Sup_In%BS_DO_FacetIsotropy          = .TRUE. !.FALSE.

      ELSE IF (BRDF_IDX == 17) THEN  
        !... ROSS-THICK HOTSPOT SURFACE (TREES WITH LARGE LEAVES OVER A FLAT SURFACE: SCALAR)

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'RtkHotSpot'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk

      ELSE IF (BRDF_IDX == 18) THEN          
        !... MODIFIED FRESNEL SURFACE (JUST WATER?: VECTOR)
        
        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'ModFresnel'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 4
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = Refrac_R
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = SigmaSq
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,3)      = 1.0_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,4)      = 0.5_fpk

      ELSE IF (BRDF_IDX == 19) THEN          
        !... SNOW SURFACE (SNOW: SCALAR)
        
        VBRDF_Sup_In%BS_WSA_VALUE = 0.75_fpk
        VBRDF_Sup_In%BS_BSA_VALUE = 0.75_fpk

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'SnowModel'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 1.0_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 3
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = 3.6_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = 5.5_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,3)      = 1.02_fpk
      END IF

      !Note! Upon adding an additional BRDF to subroutines GET_VBRDF and GET_VBRDF_PLUS,
      !      update the limit in "Max VBRDF index check" near the top of each subroutine
      !      to avoid an error due to this safety check.

      END SUBROUTINE GET_VBRDF

!*******************************************************************************
      SUBROUTINE GET_VBRDF_PLUS(BRDF_IDX,VBRDF_Sup_In,VBRDF_LinSup_In)

      USE VLIDORT_PARS_m

      USE VBRDF_SUP_KERNELS_m, Only: VBRDF_Water_RefracIndex
      USE VBRDF_SUP_MOD_m
      USE VBRDF_LINSUP_MOD_m

      IMPLICIT NONE

!  Inputs

      INTEGER, INTENT(IN) :: BRDF_IDX

!  InOut

      TYPE(VBRDF_Sup_Inputs), INTENT(INOUT) :: VBRDF_Sup_In

      TYPE(VBRDF_LinSup_Inputs), INTENT(INOUT) :: VBRDF_LinSup_In

!  Local variables

      LOGICAL   :: Do_Shadow_Effect
      REAL(fpk) :: Refrac_R, Refrac_I, Salinity, SigmaSq, SurfAlb, &
                   Wavelength, WindDir, WindSpeed

!  Max VBRDF index check

      if (MAXBRDF_IDX > 19) then
        write(*,*)
        write(*,*) 'Error: subroutine GET_VBRDF_PLUS in BRDF_Self_Tester currently set'
        write(*,*) 'up for VBRDFs 1-19.  If necessary, add additional VBRDFs to GET_VBRDF_PLUS'
        stop
      endif

!  Misc inputs

      !For Lambertian surface
      SurfAlb = 0.1_fpk      

      !For water-related surfaces
      WindSpeed  = 15.0_fpk  !in m/s
      WindDir    = 180.0_fpk   !in degrees
      !WindDir    = 70.0_fpk  !in degrees
      SigmaSq    = 0.003_fpk + 0.00512_fpk*WindSpeed !sigma**2

      Wavelength = 0.5_fpk !in um   --> Refrac_I on order of 1.0e-9
      !Wavelength = 3.0_fpk !in um   --> Refrac_I on order of 1.0e-1
      Salinity   = 34.3_fpk  !in ppt

      !Do_Shadow_Effect = .FALSE.
      Do_Shadow_Effect = .TRUE.

      CALL VBRDF_Water_RefracIndex ( Wavelength, Salinity, Refrac_R, Refrac_I )
      Refrac_I = Zero
      !write(*,*) 'Refrac_R = ',Refrac_R
      !write(*,*) 'Refrac_I = ',Refrac_I

!  Define BRDF

      !BASIC SURFACE INPUTS
      VBRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG = .FALSE.

      VBRDF_Sup_In%BS_N_BRDF_KERNELS = 1
      VBRDF_Sup_In%BS_WHICH_BRDF(1)  = BRDF_IDX

      VBRDF_Sup_In%BS_BRDF_NAMES        = ''
      VBRDF_Sup_In%BS_BRDF_FACTORS      = 0.0_fpk
      VBRDF_Sup_In%BS_N_BRDF_PARAMETERS = 0
      VBRDF_Sup_In%BS_BRDF_PARAMETERS   = 0.0_fpk

      !LINEARIZED SURFACE INPUTS
      VBRDF_LinSup_In%BS_N_SURFACE_WFS       = 0
      VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS = 0
      VBRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS = 0

      VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS = .FALSE.
      VBRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS = .FALSE.

      VBRDF_LinSup_In%BS_DO_KPARAMS_DERIVS    = .FALSE.

      VBRDF_LinSup_In%BS_DO_BSAVALUE_WF  = .FALSE.
      VBRDF_LinSup_In%BS_DO_WSAVALUE_WF  = .FALSE.
      VBRDF_LinSup_In%BS_DO_WINDSPEED_WF = .FALSE.

      !BRDF INPUTS FOR A ...
      IF (BRDF_IDX == 1) THEN
        !... LAMBERTIAN SURFACE USING GENERAL BRDF CODE (ISOTROPIC: SCALAR)

        VBRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG = .TRUE.

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Lambertian'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = SurfAlb

        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 1
        VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = 1

        VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS(1) = .TRUE.

      ELSE IF (BRDF_IDX == 2) THEN
        !... ROSS-THIN SURFACE (TREES WITH SMALL LEAVES OVER A FLAT SURFACE: SCALAR)

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Ross-thin '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk

        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 1
        VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = 1

        VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS(1) = .TRUE.

      ELSE IF (BRDF_IDX == 3) THEN
        !... ROSS-THICK SURFACE (TREES WITH LARGE LEAVES OVER A FLAT SURFACE: SCALAR)

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Ross-thick'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk

        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 1
        VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = 1

        VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS(1) = .TRUE.

      ELSE IF (BRDF_IDX == 4) THEN
        !... LI-SPARSE SURFACE (SPARSE TREES: SCALAR)

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Li-sparse '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 2
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = 1.0_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = 2.0_fpk

        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 3
        VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = 1
        VBRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS    = 2

        VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS(1) = .TRUE.
        VBRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS(1,1:2) = .TRUE.

        VBRDF_LinSup_In%BS_DO_KPARAMS_DERIVS(1)   = .TRUE.

      ELSE IF (BRDF_IDX == 5) THEN
        !... LI-DENSE SURFACE (DENSE TREES: SCALAR)

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Li-dense  '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 2
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = 2.5_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = 2.0_fpk

        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 3
        VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = 1
        VBRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS    = 2

        VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS(1) = .TRUE.
        VBRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS(1,1:2) = .TRUE.

        VBRDF_LinSup_In%BS_DO_KPARAMS_DERIVS(1)   = .TRUE.

      ELSE IF (BRDF_IDX == 6) THEN
        !... HAPKE SURFACE (???: SCALAR)

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Hapke     '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 3
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = 0.6_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = 0.06_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,3)      = 1.0_fpk

        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 4
        VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = 1
        VBRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS    = 3

        VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS(1) = .TRUE.
        VBRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS(1,1:3) = .TRUE.

        VBRDF_LinSup_In%BS_DO_KPARAMS_DERIVS(1)   = .TRUE.

      ELSE IF (BRDF_IDX == 7) THEN
        !... ROUJEAN SURFACE (RANDOM RECTANGULAR BLOCKS: SCALAR)

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Roujean   '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk

        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 1
        VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = 1

        VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS(1) = .TRUE.

      ELSE IF (BRDF_IDX == 8) THEN
        !... RAHMAN SURFACE (???: SCALAR)

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Rahman    '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 3
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      =  0.1_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = -0.5_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,3)      =  0.5_fpk
 
        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 4
        VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = 1
        VBRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS    = 3

        VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS(1) = .TRUE.
        VBRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS(1,1:3) = .TRUE.

        VBRDF_LinSup_In%BS_DO_KPARAMS_DERIVS(1)   = .TRUE.

      ELSE IF (BRDF_IDX == 9) THEN          
        !... COX-MUNK SURFACE (SUNGLINT ON WATER: SCALAR)
 
        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'Cox-Munk  '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 1.0_fpk !0.1_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 2
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = SigmaSq
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = Refrac_R**2

        !Note: Shadowing here determined indirectly by the "DO_SHADOW_EFFECT"
        !      flag through the 3rd parameter ( --> no Jacobian for 3rd parameter)
        VBRDF_Sup_In%BS_DO_SHADOW_EFFECT = Do_Shadow_Effect

        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 3
        VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = 1
        VBRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS    = 2

        VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS(1) = .TRUE.
        VBRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS(1,1:2) = .TRUE.

        VBRDF_LinSup_In%BS_DO_KPARAMS_DERIVS(1)   = .TRUE.

      ELSE IF (BRDF_IDX == 10) THEN
        !... GISS COX-MUNK SURFACE (SUNGLINT ON WATER: VECTOR)
        !    --> Note difference of parameters #1 & #2 between this kernel & kernel 9

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'GissCoxMnk'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 1.0_fpk !0.1_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 2
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = Half*SigmaSq 
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = Refrac_R

        !Note: Shadowing here determined indirectly by the "DO_SHADOW_EFFECT"
        !      flag through the 3rd parameter ( --> no Jacobian for 3rd parameter)
        VBRDF_Sup_In%BS_DO_SHADOW_EFFECT = Do_Shadow_Effect

        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 2
        VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = 1
        VBRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS    = 1

        VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS(1) = .TRUE.
        VBRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS(1,1) = .TRUE.

        VBRDF_LinSup_In%BS_DO_KPARAMS_DERIVS(1)   = .TRUE.

      ELSE IF (BRDF_IDX == 11) THEN
        !... GISS COX-MUNK COMPLEX REFRACTIVE INDEX ("CRI") SURFACE (SUNGLINT ON WATER: VECTOR)
        !    --> Note difference of parameters #1 & #2 between this kernel & kernel 9

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'GCMcomplex'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 1.0_fpk !0.1_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 3
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = Half*SigmaSq
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = Refrac_R
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,3)      = Refrac_I

        !Note: Shadowing here determined directly by the "DO_SHADOW_EFFECT" flag
        VBRDF_Sup_In%BS_DO_SHADOW_EFFECT = Do_Shadow_Effect

        !Note: No analytic surface Jacobians with this kernel

      ELSE IF (BRDF_IDX == 12) THEN
        !... BPDF-SOIL SURFACE (SOIL: VECTOR)
        
        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'BPDF-Soil '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 1.0_fpk !0.99_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 1
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = Refrac_R

        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 2
        VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = 1
        VBRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS    = 1

        VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS(1) = .TRUE.
        VBRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS(1,1) = .TRUE.

        VBRDF_LinSup_In%BS_DO_KPARAMS_DERIVS(1)   = .TRUE.

      ELSE IF (BRDF_IDX == 13) THEN
        !... BPDF-VEGN SURFACE (VEGETATION: VECTOR)
        
        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'BPDF-Vegn '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 1
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = Refrac_R

        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 2
        VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = 1
        VBRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS    = 1

        VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS(1) = .TRUE.
        VBRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS(1,1) = .TRUE.

        VBRDF_LinSup_In%BS_DO_KPARAMS_DERIVS(1)   = .TRUE.

      ELSE IF (BRDF_IDX == 14) THEN          
        !... BPDF-NDVI SURFACE (NDVI-PARAMETERIZED: VECTOR)
        
        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'BPDF-NDVI '
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 3
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = Refrac_R
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = 0.2_fpk !NDVI
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,3)      = 1.0_fpk

        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 4
        VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = 1
        VBRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS    = 3

        VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS(1) = .TRUE.
        VBRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS(1,1:3) = .TRUE.

        VBRDF_LinSup_In%BS_DO_KPARAMS_DERIVS(1)   = .TRUE.

      ELSE IF (BRDF_IDX == 15) THEN          
        !... NEW COX-MUNK SURFACE (SUNGLINT ON WATER: SCALAR)

        VBRDF_Sup_In%BS_DO_NewCMGLINT             = .TRUE.
        VBRDF_Sup_In%BS_DO_NewGCMGLINT            = .FALSE.

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'NewCMGlint'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 1.0_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 2
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = WindSpeed
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = Salinity

        VBRDF_Sup_In%BS_WINDSPEED                 = WindSpeed
        VBRDF_Sup_In%BS_WINDDIR                   = WindDir
        VBRDF_Sup_In%BS_WAVELENGTH                = Wavelength
        VBRDF_Sup_In%BS_SALINITY                  = Salinity

        VBRDF_Sup_In%BS_DO_GlintShadow            = Do_Shadow_Effect
        VBRDF_Sup_In%BS_DO_FoamOption             = .FALSE. !.TRUE.
        VBRDF_Sup_In%BS_DO_FacetIsotropy          = .TRUE. !.FALSE.

        !Note: Only surface Jacobian wrt wind speed (par #1) available here

        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 1
        VBRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS    = 1

        VBRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS(1,1) = .TRUE.
        VBRDF_LinSup_In%BS_DO_KPARAMS_DERIVS(1)   = .TRUE.

      ELSE IF (BRDF_IDX == 16) THEN          
        !... NEW GISS COX-MUNK SURFACE (SUNGLINT ON WATER: VECTOR)

        VBRDF_Sup_In%BS_DO_NewCMGLINT             = .FALSE.
        VBRDF_Sup_In%BS_DO_NewGCMGLINT            = .TRUE.

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'NewGCMGlit'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 1.0_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 2
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = WindSpeed
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = Salinity

        VBRDF_Sup_In%BS_WINDSPEED                 = WindSpeed
        VBRDF_Sup_In%BS_WINDDIR                   = WindDir
        VBRDF_Sup_In%BS_WAVELENGTH                = Wavelength
        VBRDF_Sup_In%BS_SALINITY                  = Salinity

        VBRDF_Sup_In%BS_DO_GlintShadow            = Do_Shadow_Effect
        VBRDF_Sup_In%BS_DO_FoamOption             = .FALSE. !.TRUE.
        VBRDF_Sup_In%BS_DO_FacetIsotropy          = .TRUE. !.FALSE.

        !Note: Only surface Jacobian wrt wind speed (par #1) available here

        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 1
        VBRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS    = 1

        VBRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS(1,1) = .TRUE.
        VBRDF_LinSup_In%BS_DO_KPARAMS_DERIVS(1)   = .TRUE.

      ELSE IF (BRDF_IDX == 17) THEN  
        !... ROSS-THICK HOTSPOT SURFACE (TREES WITH LARGE LEAVES OVER A FLAT SURFACE: SCALAR)

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'RtkHotSpot'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk

        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 1
        VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = 1

        VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS(1) = .TRUE.

      ELSE IF (BRDF_IDX == 18) THEN          
        !... MODIFIED FRESNEL SURFACE (JUST WATER?: VECTOR)
        
        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'ModFresnel'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 0.99_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 4
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = Refrac_R
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = SigmaSq
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,3)      = 1.0_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,4)      = 0.5_fpk

        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 5
        VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = 1
        VBRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS    = 4

        VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS(1) = .TRUE.
        VBRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS(1,1:4) = .TRUE.

        VBRDF_LinSup_In%BS_DO_KPARAMS_DERIVS(1)   = .TRUE.

      ELSE IF (BRDF_IDX == 19) THEN          
        !... SNOW SURFACE (SNOW: SCALAR)

        VBRDF_Sup_In%BS_WSA_VALUE = 0.75_fpk
        VBRDF_Sup_In%BS_BSA_VALUE = 0.75_fpk

        VBRDF_Sup_In%BS_BRDF_NAMES(1)             = 'SnowModel'
        VBRDF_Sup_In%BS_BRDF_FACTORS(1)           = 1.0_fpk
        VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)      = 3
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1)      = 3.6_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2)      = 5.5_fpk
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,3)      = 1.02_fpk

        VBRDF_LinSup_In%BS_N_SURFACE_WFS          = 3
        VBRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = 1
        VBRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS    = 2

        VBRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS(1) = .TRUE.
        VBRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS(1,1:2) = .TRUE.

        VBRDF_LinSup_In%BS_DO_KPARAMS_DERIVS(1)   = .TRUE.

      END IF

      !Note! Upon adding an additional BRDF to subroutines GET_VBRDF and GET_VBRDF_PLUS,
      !      update the limit in "Max VBRDF index check" near the top of each subroutine
      !      to avoid an error due to this safety check.

      END SUBROUTINE GET_VBRDF_PLUS

!*******************************************************************************
!*******************************************************************************

