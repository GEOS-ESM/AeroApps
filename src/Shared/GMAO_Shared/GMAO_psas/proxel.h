!*************************************************************************
!                                                                        *
!                  NASA/Goddard Space Flight Center                      *
!                     Laboratory for Atmospheres                         *
!                      Data Assimilation Office                          *
!                            Code 910.3                                  *
!                                                                        *
!            PHYSICAL-SPACE STATISTICAL ANALYSIS SYSTEM (PSAS)           *
!                                                                        *
!*************************************************************************

!.........................................................................
! !BOP
!
! !FILE: proxel.h
! 
! !ROUTINE: n/a
! 
! !DESCRIPTION: 
!
!     Include file for Proximity Elimination package (proxel.f).
!
! !CALLING SEQUENCE: n/a
!
! !INPUT PARAMETERS: n/a
!
! !OUTPUT PARAMETERS: n/a
!
! !SEE ALSO: proxelim.f
!
! !SYSTEM ROUTINES: none
!
! !FILES USED: none
!
! !WRITTEN BY: Arlindo da Silva, 26sep94
! 
! !REVISION HISTORY:
!
! !EOP
!.........................................................................
!

!     Black-list data structure
!     -------------------------
      parameter ( nlist = 1000 )         ! max No. obs in list
      integer      idx_lst(nlist)       ! Index in master list
      integer      kx_lst(nlist)        ! observation data source
      integer      kt_lst(nlist)        ! observation type
      integer      ireg_lst(nlist)      ! region index
      logical      kl_lst(nlist)        ! logical flag
      real         rlats_lst(nlist)     ! latitude of obs in deg
      real         rlons_lst(nlist)     ! longitude of obs in deg
      real         rlevs_lst(nlist)     ! levels in mb
      real         del_lst(nlist)       ! innovations
      real         sigOc_lst(nlist)      ! obs. error std
      real         sigOu_lst(nlist)      ! obs. error std
      real         sigF_lst(nlist)      ! forecast error std.
      real         x_lst(nlist)         ! cartesian x on unity sphere
      real         y_lst(nlist)         ! cartesian y on unity sphere
      real         z_lst(nlist)         ! cartesian z on unity sphere

      common / rproxdat / rlats_lst, rlons_lst,
     .                   rlevs_lst, del_lst,
     .			 sigOc_list,sigOu_lst, sigF_lst,
     .                   x_lst, y_lst, z_lst,
     .                   ilist, idx_lst,
     .                   kx_lst, kt_lst, kl_lst 
 

!     Threshold log-p distance to superob observattions
!     -------------------------------------------------
!cc      parameter ( DELLNP = 0.1 )

!     Threshold horizontal distance to superob 
!      observations (in km)
!     --------------------------------------------------
!cc      parameter ( Rkm = 80. )                  ! units: km
      parameter ( Rerth = 6.37 E 3 )            ! Radius earth in km
!cc      parameter ( Rtresh = Rkm / Rerth )        ! normalized R
!cc      parameter ( Rtresh2 = Rtresh * Rtresh )  
      common / proxprm / Rtresh2, dellnp

!     Flag to print superobs
!     ----------------------
      logical debug
      common / proxlog / debug

!     A small number (smaller than longitude "resolution")
!     ----------------------------------------------------
      parameter ( eps = 0.00001 )

!     A big integer
!     -------------
      parameter ( IBIG = 2**30 )

!     Counter of eliminate obs. by data source (kx)
!     ---------------------------------------------
      integer kxelim(kxmax)

!     Separation limit (in degrees) for regions to be included
!      in search
!     --------------------------------------------------------
      parameter ( SEPLIM = 26.5 )
      Real, Parameter :: cosseplim = 0.894934

!     ------------------------------------------------
!     NOTE: kxrank is initialized in routine 'kxname0'
!     ------------------------------------------------

