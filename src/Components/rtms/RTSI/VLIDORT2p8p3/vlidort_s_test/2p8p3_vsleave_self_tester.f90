
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

      program VSLEAVE_Self_Tester

!  New test of VSLEAVE supplement, Version 2.8
!  --------------------------------------------

!   This is designed to produce multigeometrical output for
!      1. VSLEAVE Exact (direct-bounce) reflectance
!      2. diffuse-field (mutliply-scattered) VSLEAVE Fourier components
        
!  Upgrade for Version 2.8.1, August 2019
!  ---------------------------------------

!  Module files for VLIDORT. Strict usage, Version 2.8 upwards

      USE VLIDORT_PARS_m
      USE VSLEAVE_SUP_AUX_m, Only : VSLEAVE_READ_ERROR
      USE VSLEAVE_SUP_MOD_m
      USE VSLEAVE_LINSUP_MOD_m

!  Implicit none

      IMPLICIT NONE

!  VSLEAVE supplement file inputs status structure

      TYPE(VSLEAVE_Input_Exception_Handling)   :: VSLEAVE_Sup_InputStatus

!  VSLEAVE supplement input structures

      TYPE(VSLEAVE_Sup_Inputs)                 :: VSLEAVE_Sup_In
      TYPE(VSLEAVE_LinSup_Inputs)              :: VSLEAVE_LinSup_In

!  VSLEAVE supplement output structures

      TYPE(VSLEAVE_Sup_Outputs)                :: VSLEAVE_Sup_Out
      TYPE(VSLEAVE_LinSup_Outputs)             :: VSLEAVE_LinSup_Out 
      TYPE(VSLEAVE_Output_Exception_Handling)  :: VSLEAVE_Sup_OutputStatus

!  Saved Results
!  =============

      INTEGER, parameter :: MAXJACS = 7

      REAL(fpk) :: ISLEAVE_EX    ( MAXSTOKES, MAXBEAMS, 0:MAXJACS )
      REAL(fpk) :: LS_ISLEAVE_EX ( MAXSTOKES, MAXBEAMS, MAXJACS )

      REAL(fpk) :: ESLEAVE_EX    ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS, 0:MAXJACS )
      REAL(fpk) :: LS_ESLEAVE_EX ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS, MAXJACS )

      REAL(fpk) :: MSLEAVE_SD    ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS, 0:MAXJACS )
      REAL(fpk) :: LS_MSLEAVE_SD ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS, MAXJACS )

      REAL(fpk) :: MSLEAVE_SU    ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS, 0:MAXJACS )
      REAL(fpk) :: LS_MSLEAVE_SU ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS, MAXJACS )

!  Fluorescence data
!  =================

      DOUBLE PRECISION, DIMENSION(3,2) :: FL_DataGAUSSIANS

      !FL_DataGAUSSIANS(1,1) =   1.445d0
      !FL_DataGAUSSIANS(2,1) = 736.800d0
      !FL_DataGAUSSIANS(3,1) =  21.200d0
      !FL_DataGAUSSIANS(1,2) =   0.868d0
      !FL_DataGAUSSIANS(2,2) = 685.200d0
      !FL_DataGAUSSIANS(3,2) =   9.550d0

      DATA FL_DataGAUSSIANS / 1.445d0, 736.8d0, 21.2d0, 0.868d0, 685.2d0, 9.55d0 /

!  Local Variables
!  ===============

      INTEGER           :: test
      INTEGER           :: i, ia, ib, g, j, jw, k, m, n, o1, q, t, um
      INTEGER           :: nstokes, nbeams, nusers, nazims, ndisos, nmoms, njacs, ntests
      CHARACTER(LEN=1)  :: c1, cnj
      CHARACTER(LEN=5)  :: surf
      CHARACTER(LEN=30) :: jacform, jacheader ( maxjacs )
      CHARACTER(LEN=80) :: infile, outfile
      CHARACTER(LEN=99) :: trace

      REAL(fpk)         :: eps, epsfac, FL_Amplitude755_SAVE, CHLORCONC_SAVE, WINDSPEED_SAVE
      REAL(fpk)         :: DELPARS ( maxjacs )

!  Start test loop

      ntests = 2

      DO test=1,ntests
      !DO test=1,1
      !DO test=2,2

        write(*,*) 'doing test ',test

!  Define input file

        IF ( test == 1 ) THEN
          surf = 'land'
        ELSE
          surf = 'water'
        ENDIF
        infile = 'vlidort_s_test/2p8p3_vsleave_self_tester_' // trim(surf) // '.cfg'
        write(*,*) 'trim(infile) = |',trim(infile),'|'

!  Get the VSLEAVE inputs

        CALL VSLEAVE_LIN_INPUTMASTER ( &
          infile,                      & ! Input
          VSLEAVE_Sup_In,              & ! Outputs
          VSLEAVE_LinSup_In,           & ! Outputs
          VSLEAVE_Sup_InputStatus )      ! Outputs

        IF ( VSLEAVE_Sup_InputStatus%SL_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
          CALL VSLEAVE_READ_ERROR ( '2p8p3_vsleave_self_tester.log', VSLEAVE_Sup_InputStatus )

!  Define some inputs not handled by VSLEAVE_LIN_INPUTMASTER

        VSLEAVE_Sup_In%SL_VSLEAVE_DATAPATH = 'vlidort_s_test/data'

!  Setup test

        if (test == 1) then
          !Land (fluorescence) test
          if ( .not. VSLEAVE_Sup_In%SL_FL_DO_DataGaussian ) then
            VSLEAVE_Sup_In%SL_FL_InputGAUSSIANS(1:3,1) = FL_DataGAUSSIANS(1:3,1)
            VSLEAVE_Sup_In%SL_FL_InputGAUSSIANS(1:3,2) = FL_DataGAUSSIANS(1:3,2)
          endif
        endif

!  Setup Proxies

        nstokes = VSLEAVE_Sup_In%SL_nstokes
        nbeams  = VSLEAVE_Sup_In%SL_nbeams
        nusers  = VSLEAVE_Sup_In%SL_n_user_streams
        nazims  = VSLEAVE_Sup_In%SL_n_user_relazms
        ndisos  = VSLEAVE_Sup_In%SL_nstreams
        nmoms   = 0!2*ndisos - 1

        njacs  = VSLEAVE_LinSup_In%SL_N_SLEAVE_WFS

        IF ( test == 1 ) THEN
          !Land test
          eps = 1.0d-05
        ELSE
          !Water test
          eps = 1.0d-03
        ENDIF
        epsfac = 1.0d0 + eps

!  VSLEAVE baseline call

        CALL VSLEAVE_LIN_MAINMASTER ( &
          VSLEAVE_Sup_In,          & ! Inputs
          VSLEAVE_LinSup_In,       & ! Inputs
          VSLEAVE_Sup_Out,         & ! Outputs
          VSLEAVE_LinSup_Out,      & ! Outputs
          VSLEAVE_Sup_OutputStatus ) ! Output Status

        if ( VSLEAVE_Sup_OutputStatus%SL_STATUS_OUTPUT.ne.VLIDORT_SUCCESS ) THEN
           TRACE = 'VSLEAVE_LIN_MAINMASTER failed' ; go to 678
        endif

!  Save baseline

        q = 0
        ISLEAVE_EX(1:nstokes,1:nbeams,q) = &
          VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1:nstokes,1:nbeams)
        ESLEAVE_EX(1:nstokes,1:nusers,1:nazims,1:nbeams,q) = &
          VSLEAVE_Sup_Out%SL_SLTERM_USERANGLES(1:nstokes,1:nusers,1:nazims,1:nbeams)
        MSLEAVE_SD(0:nmoms,1:nstokes,1:ndisos,1:nbeams,q)  = &
          VSLEAVE_Sup_Out%SL_SLTERM_F_0(0:nmoms,1:nstokes,1:ndisos,1:nbeams)
        MSLEAVE_SU(0:nmoms,1:nstokes,1:nusers,1:nbeams,q)  = &
          VSLEAVE_Sup_Out%SL_USER_SLTERM_F_0(0:nmoms,1:nstokes,1:nusers,1:nbeams)

        do q = 1, njacs
          LS_ISLEAVE_EX(1:nstokes,1:nbeams,q) = &
            VSLEAVE_LinSup_Out%SL_LS_SLTERM_ISOTROPIC(q,1:nstokes,1:nbeams)
          LS_ESLEAVE_EX(1:nstokes,1:nusers,1:nazims,1:nbeams,q) = &
            VSLEAVE_LinSup_Out%SL_LS_SLTERM_USERANGLES(q,1:nstokes,1:nusers,1:nazims,1:nbeams)
          LS_MSLEAVE_SD(0:nmoms,1:nstokes,1:ndisos,1:nbeams,q)  = &
            VSLEAVE_LinSup_Out%SL_LS_SLTERM_F_0(q,0:nmoms,1:nstokes,1:ndisos,1:nbeams)
          LS_MSLEAVE_SU(0:nmoms,1:nstokes,1:nusers,1:nbeams,q)  = &
            VSLEAVE_LinSup_Out%SL_LS_USER_SLTERM_F_0(q,0:nmoms,1:nstokes,1:nusers,1:nbeams)
        end do

!  FD calcs

        if (test == 1) then
          !Land (fluorescence) test

          !(1) Fluorescence amplitude @ 755nm 
          q = 1
          FL_Amplitude755_SAVE              = VSLEAVE_Sup_In%SL_FL_Amplitude755
          VSLEAVE_Sup_In%SL_FL_Amplitude755 = FL_Amplitude755_SAVE * epsfac
          DELPARS(q)                        = FL_Amplitude755_SAVE * eps

          CALL VSLEAVE_MAINMASTER (  &
            VSLEAVE_Sup_In,          & ! Inputs
            VSLEAVE_Sup_Out,         & ! Outputs
            VSLEAVE_Sup_OutputStatus ) ! Output Status

          if ( VSLEAVE_Sup_OutputStatus%SL_STATUS_OUTPUT.ne.VLIDORT_SUCCESS ) THEN
             write(c1,'(i1)')q ; TRACE = 'VSLEAVE_MAINMASTER failed, perturbation # '//c1 ; go to 678
          endif

          ISLEAVE_EX(1:nstokes,1:nbeams,q) = &
            VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1:nstokes,1:nbeams)
          ESLEAVE_EX(1:nstokes,1:nusers,1:nazims,1:nbeams,q) = &
            VSLEAVE_Sup_Out%SL_SLTERM_USERANGLES(1:nstokes,1:nusers,1:nazims,1:nbeams)
          MSLEAVE_SD(0:nmoms,1:nstokes,1:ndisos,1:nbeams,q)  = &
            VSLEAVE_Sup_Out%SL_SLTERM_F_0(0:nmoms,1:nstokes,1:ndisos,1:nbeams)
          MSLEAVE_SU(0:nmoms,1:nstokes,1:nusers,1:nbeams,q)  = &
            VSLEAVE_Sup_Out%SL_USER_SLTERM_F_0(0:nmoms,1:nstokes,1:nusers,1:nbeams)

          VSLEAVE_Sup_In%SL_FL_Amplitude755 = FL_Amplitude755_SAVE

          !(2) Fluorescence GAUSSIAN parameters
          do q = 2, njacs
            if ( q.lt.5 ) k = q-1 ; if ( q.lt.5 ) g=1
            if ( q.ge.5 ) k = q-4 ; if ( q.ge.5 ) g=2
            VSLEAVE_Sup_In%SL_FL_InputGAUSSIANS(k,g) = FL_DataGAUSSIANS(k,g) * epsfac
            DELPARS(q)                               = FL_DataGAUSSIANS(k,g) * eps

            CALL VSLEAVE_MAINMASTER ( &
              VSLEAVE_Sup_In,          & ! Inputs
              VSLEAVE_Sup_Out,         & ! Outputs
              VSLEAVE_Sup_OutputStatus ) ! Output Status

            if ( VSLEAVE_Sup_OutputStatus%SL_STATUS_OUTPUT.ne.VLIDORT_SUCCESS ) THEN
               write(c1,'(i1)')q ; TRACE = 'VSLEAVE_MAINMASTER failed, perturbation # '//c1 ; go to 678
            endif

            !sleave_fd(jw,q) = ( sleave_pert(jw) - sleave_save(jw) ) / eps / FL_DataGAUSSIANS(k,g)

            ISLEAVE_EX(1:nstokes,1:nbeams,q) = &
              VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1:nstokes,1:nbeams)
            ESLEAVE_EX(1:nstokes,1:nusers,1:nazims,1:nbeams,q) = &
              VSLEAVE_Sup_Out%SL_SLTERM_USERANGLES(1:nstokes,1:nusers,1:nazims,1:nbeams)
            MSLEAVE_SD(0:nmoms,1:nstokes,1:ndisos,1:nbeams,q)  = &
              VSLEAVE_Sup_Out%SL_SLTERM_F_0(0:nmoms,1:nstokes,1:ndisos,1:nbeams)
            MSLEAVE_SU(0:nmoms,1:nstokes,1:nusers,1:nbeams,q)  = &
              VSLEAVE_Sup_Out%SL_USER_SLTERM_F_0(0:nmoms,1:nstokes,1:nusers,1:nbeams)

            VSLEAVE_Sup_In%SL_FL_InputGAUSSIANS(k,g) = FL_DataGAUSSIANS(k,g)

          enddo
        else
          !Water test

          !(1) Chlorophyll concentration 
          q = 1
          CHLORCONC_SAVE              = VSLEAVE_Sup_In%SL_CHLORCONC
          VSLEAVE_Sup_In%SL_CHLORCONC = CHLORCONC_SAVE * epsfac
          DELPARS(q)                  = CHLORCONC_SAVE * eps

          CALL VSLEAVE_MAINMASTER ( &
            VSLEAVE_Sup_In,          & ! Inputs
            VSLEAVE_Sup_Out,         & ! Outputs
            VSLEAVE_Sup_OutputStatus ) ! Output Status

          if ( VSLEAVE_Sup_OutputStatus%SL_STATUS_OUTPUT.ne.VLIDORT_SUCCESS ) THEN
             write(c1,'(i1)')q ; TRACE = 'VSLEAVE_MAINMASTER failed, perturbation # '//c1 ; go to 678
          endif

          ISLEAVE_EX(1:nstokes,1:nbeams,q) = &
            VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1:nstokes,1:nbeams)
          ESLEAVE_EX(1:nstokes,1:nusers,1:nazims,1:nbeams,q) = &
            VSLEAVE_Sup_Out%SL_SLTERM_USERANGLES(1:nstokes,1:nusers,1:nazims,1:nbeams)
          MSLEAVE_SD(0:nmoms,1:nstokes,1:ndisos,1:nbeams,q)  = &
            VSLEAVE_Sup_Out%SL_SLTERM_F_0(0:nmoms,1:nstokes,1:ndisos,1:nbeams)
          MSLEAVE_SU(0:nmoms,1:nstokes,1:nusers,1:nbeams,q)  = &
            VSLEAVE_Sup_Out%SL_USER_SLTERM_F_0(0:nmoms,1:nstokes,1:nusers,1:nbeams)

          VSLEAVE_Sup_In%SL_CHLORCONC = CHLORCONC_SAVE

          !(2) Surface wind speed
          q = 2
          WINDSPEED_SAVE              = VSLEAVE_Sup_In%SL_WINDSPEED
          VSLEAVE_Sup_In%SL_WINDSPEED = WINDSPEED_SAVE * epsfac
          DELPARS(q)                  = WINDSPEED_SAVE * eps

          CALL VSLEAVE_MAINMASTER (  &
            VSLEAVE_Sup_In,          & ! Inputs
            VSLEAVE_Sup_Out,         & ! Outputs
            VSLEAVE_Sup_OutputStatus ) ! Output Status

          if ( VSLEAVE_Sup_OutputStatus%SL_STATUS_OUTPUT.ne.VLIDORT_SUCCESS ) THEN
             write(c1,'(i1)')q ; TRACE = 'VSLEAVE_MAINMASTER failed, perturbation # '//c1 ; go to 678
          endif

          ISLEAVE_EX(1:nstokes,1:nbeams,q) = &
            VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1:nstokes,1:nbeams)
          ESLEAVE_EX(1:nstokes,1:nusers,1:nazims,1:nbeams,q) = &
            VSLEAVE_Sup_Out%SL_SLTERM_USERANGLES(1:nstokes,1:nusers,1:nazims,1:nbeams)
          MSLEAVE_SD(0:nmoms,1:nstokes,1:ndisos,1:nbeams,q)  = &
            VSLEAVE_Sup_Out%SL_SLTERM_F_0(0:nmoms,1:nstokes,1:ndisos,1:nbeams)
          MSLEAVE_SU(0:nmoms,1:nstokes,1:nusers,1:nbeams,q)  = &
            VSLEAVE_Sup_Out%SL_USER_SLTERM_F_0(0:nmoms,1:nstokes,1:nusers,1:nbeams)

          VSLEAVE_Sup_In%SL_WINDSPEED = WINDSPEED_SAVE

        endif

!  Define & open output file

        outfile = 'vlidort_s_test/results_vsleave_self_test_' // trim(surf) // '.all'
        open(36,file=trim(outfile),status='unknown')

!  Write results

        write(36,*)
        write(36,'(A)') '****************************************'
        if (test == 1) then
          write(36,'(T4,A)')'Land surface-leaving (Fluorescence)'
        else
          write(36,'(T9,A)')'Water surface-leaving'
        endif
        write(36,'(A)') '****************************************'

        do o1 = 1, nstokes

        write(36,'(/T10,A,i1)')'For Stokes: ',o1

!  Basic output

        write(36,'(/T10,A)')'STANDARD QUANTITIES'

        write(36,'(/A//T15,A/)')'Isotropic SLEAVE, all geometries','sza     Iso SLEAVE value'
        do ib = 1, nbeams
          write(36,'(A,T15,i2,5x,es15.6)')'Iso',ib,ISLEAVE_EX(o1,ib,0)
        enddo

        if (test == 2) then
          !Water test
          write(36,'(/A//T15,A/)')'Exact SLEAVE, all geometries','sza vza azm    Exact SLEAVE value'
          do ib = 1, nbeams
            do um = 1, nusers
              do ia = 1, nazims
                write(36,'(A,T15,3(i2,2x),es16.6)')'Exact',ib,um,ia,ESLEAVE_EX(o1,um,ia,ib,0)
              enddo
            enddo
          enddo

          write(36,'(/A)')'Fourier Components:'
          DO M = 0, NMOMS
            write(36,'(/T15,A,i2//T15,A/)')'Fourier component # ',M,'sza str    SLEAVE value'
            do ib = 1, nbeams
              do j = 1, ndisos
                write(36,'(A,T15,2(i2,2x),es14.6)')'Solars-Stream',ib,j,MSLEAVE_SD(M,o1,j,ib,0)
              enddo
            enddo

            write(36,'(/T16,A/)')'sza vza    SLEAVE value'
            do ib = 1, nbeams
              do um = 1, nusers
                write(36,'(A,T15,2(i2,2x),es14.6)')'Solars-Users ',ib,um,MSLEAVE_SU(M,o1,um,ib,0)
              enddo
            enddo
          ENDDO
        endif

!  Jacobian output

        write(36,'(//T10,A)')'WEIGHTING FUNCTIONS (ANALYTIC + FINITE_DIFFERENCE)'

        write(cnj,'(i1)') njacs
        jacheader(1:njacs) = '       Analytic    Finite-diff'

        jacform='(/A//T15,A,' // cnj // 'A/)'
        write(36,jacform)'Isotropic SLEAVE, all geometries',&
                         'sza',jacheader(1:njacs)

        do ib = 1, nbeams
          write(36,'(A,T15,i2,2x,7(2x,2es14.6))')'Iso',ib,&
            ( LS_ISLEAVE_EX(o1,ib,q),&
              (ISLEAVE_EX(o1,ib,q)-ISLEAVE_EX(o1,ib,0))/DELPARS(q),q=1,njacs )
        enddo

        if (test == 2) then
          write(36,jacform)'Exact SLEAVE, all geometries',&
                           'sza vza azm',jacheader(1:njacs)
          do ib = 1, nbeams
            do um = 1, nusers
              do ia = 1, nazims
                write(36,'(A,T15,3(i2,2x),7(2x,2es14.6))')'Exact',ib,um,ia,&
                  ( LS_ESLEAVE_EX(o1,um,ia,ib,q),&
                    (ESLEAVE_EX(o1,um,ia,ib,q)-ESLEAVE_EX(o1,um,ia,ib,0))/DELPARS(q),q=1,njacs )
              enddo
            enddo
          enddo

          write(36,'(/A)')'Fourier Components:'
          DO M = 0, NMOMS
            jacform='(/T15,A,i2//T15,A,' // cnj // 'A/)'
            write(36,jacform)'Fourier component # ',M,&
                             'sza str',jacheader(1:njacs)
            do ib = 1, nbeams
              do j = 1, ndisos
                write(36,'(A,T15,2(i2,2x),7(2x,2es14.6))')'Solars-Stream',ib,j,&
                  ( LS_MSLEAVE_SD(M,o1,j,ib,q),&
                    (MSLEAVE_SD(M,o1,j,ib,q)-MSLEAVE_SD(M,o1,j,ib,0))/DELPARS(q),q=1,njacs )
              enddo
            enddo

            jacform='(/T15,A,' // cnj // 'A/)'
            write(36,jacform)'sza vza',jacheader(1:njacs)
            do ib = 1, nbeams
              do um = 1, nusers
                 write(36,'(A,T15,2(i2,2x),7(2x,2es14.6))')'Solars-Users ',ib,um,&
                   ( LS_MSLEAVE_SU(M,o1,um,ib,q),&
                     (MSLEAVE_SU(M,o1,um,ib,q)-MSLEAVE_SU(M,o1,um,ib,0))/DELPARS(q),q=1,njacs )
              enddo
            enddo
          ENDDO
        endif

        enddo

!  Close output file

        close(36)

!  End test loop

      enddo

 98   format(3I4,1p10e16.7)
 99   format(I4,1p10e16.7)

      write(*,*)
      write(*,*) 'Main program finished successfully'

      stop

!  Error finish

678   continue

      write(*,*)
      write(*,'(1x,a)') trim(TRACE)
      write(*,*)'Number of error messages from VSLEAVE calculation = ',&
        VSLEAVE_Sup_OutputStatus%SL_NOUTPUTMESSAGES
      do i = 1, VSLEAVE_Sup_OutputStatus%SL_NOUTPUTMESSAGES
         write(*,'(a,i2,a,a)')' Message # ', i, ': ',&
           adjustl(trim(VSLEAVE_Sup_OutputStatus%SL_OUTPUTMESSAGES(i)))
      enddo

!  Finish

end program VSLEAVE_Self_Tester

