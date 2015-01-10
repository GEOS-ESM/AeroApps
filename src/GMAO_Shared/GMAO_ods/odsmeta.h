!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !INCLUDE: odsmeta.h --- Defines the ODS metadata conventions
!
! !INTERFACE:

!  include 'odsmeta.h'

! !DESCRIPTION:
!
!  Defines ODS metadata conventions. To be included by SPMD and serial
!          versiosn of m_odsmeta.
!
! !REVISION HISTORY:
!
!  01oct2005  da Silva  From m_odsmeta.f
!  31Jan2005  Todling   Added varQC exclusion flag
!  17May2006  Todling   Spelled out kt names
!  14Dec2006  Todling   Move sat instrument def from scan to here
!  08Apr2008  Sienkiewicz  define OMI in satellite instrument list
!  18Apr2008  Ravi/RT   Updated KtSurfAll fields; bug fix in kt=45/46
!  17Mar2009  Meunier   Add position as observation type (lat/lon)
!  28Apr2009  Todling   Add IASI
!  05Apr2009  Todling   Add refractivity and bending angle
!  20Apr2010  Todling   Increase kxmax to 1024
!  30May2012  Todling   ATMS/CrIS; distinguish GPS ref/bend; consistent w/ kt_list.rc
!  24Oct2012  Sienkiewicz  Take out obsolete ssmis references (las,uas,etc)
!                          add 'seviri', change SSMIS to same idsats as SSMI
!
!EOP
!-------------------------------------------------------------------------
!BOC

!                             -------------------------
!                             Data Type (kt) Parameters
!                             --------------------------
!
!       10Apr93 - Jing G.       - Initial code (kt_max.h)
!       24Mar95 - Jing G.       - New version for text based data table.
!       09Mar98 - Todling       - Added all kt definition from ON-95.01
!       10Apr98 - Todling       - Updated kt according to Y.Kondratyeva
!       29Dec99 - Todling       - Embedded in ODSmeta module.
!
        integer, parameter  :: ktmax   = 128

!.......................................................................

                                                ! Surface variables
        integer, parameter  :: ktus    =  1     !   sea level zonal wind (m/s)
        integer, parameter  :: ktvs    =  2     !   sea level meridional wind (m/s)
        integer, parameter  :: ktslp   =  3     !   sea level pressure (hPa)
                                                ! Upper-air variables
        integer, parameter  :: ktuu    =  4     !   zonal wind (m/s)
        integer, parameter  :: ktvv    =  5     !   meridional wind (m/s)
        integer, parameter  :: ktHH    =  6     !   geopotential height (m)
        integer, parameter  :: ktww    =  7     !   water vapor mixing ratio (g/Kg);
                                                !     CAUTION: used to be ktqq
                                                !     won't work w/ PSAS
        integer, parameter  :: ktTT    =  8     !   temperature (K)
        integer, parameter  :: ktTd    =  9     !   dew-point temperature (K)
        integer, parameter  :: ktrh    = 10     !   relative humidity (%)
        integer, parameter  :: ktqq    = 11     !   specific humidity (g/Kg); 
                                                !     CAUTION: see ktww above
                                                ! More surface variables
        integer, parameter  :: ktus10  = 12     !   10 meter zonal wind (m/s)
        integer, parameter  :: ktTs10  = 13     !   10 meter meridional wind (K)
        integer, parameter  :: ktTds   = 14     !   10 meter dew-point temperature (K)
        integer, parameter  :: ktrhs   = 15     !   10 meter relative humidity (%)
        integer, parameter  :: ktqs10  = 16     !   10 meter specific humidity (g/Kg)

        integer, parameter  :: ktpr    = 17     !   precipitation rate (mm/day) 
        integer, parameter  :: kttpw   = 18     !   total precipitable water (mm)
        integer, parameter  :: kttlw   = 19     !   total cloud liquid water (mm)
        integer, parameter  :: ktfcc   = 20     !   fractional cloud cover (%)
                                                ! Chemistry
        integer, parameter  :: kttco3  = 21     !   total column ozone (Dobson)
        integer, parameter  :: kto3    = 22     !   layer ozone (Dobson)
        integer, parameter  :: ktuthk  = 23     !   height thickness (m)

                                                ! Various
        integer, parameter  :: ktspd2m   = 24   !   2 meter wind speed (m/s)
        integer, parameter  :: ktmxspd2m = 25   !   2 meter max wind speed (m/s)
        integer, parameter  :: ktwgust2m = 26   !   2 meter max wind gust (m/s)
        integer, parameter  :: ktt2m     = 27   !   2 meter temperature (K)
        integer, parameter  :: ktmxt2m   = 28   !   2 meter maximum temperature (K)
        integer, parameter  :: ktmmnt2m  = 29   !   2 meter minimum temperature (K)
        integer, parameter  :: ktdewt2m  = 30   !   2 meter dew-point temperature
        integer, parameter  :: ktrh2m    = 31   !   2m relative humidity (%)
        integer, parameter  :: ktsphu2m  = 32   !   2 meter specific humidity (g/Kg)
        integer, parameter  :: ktps2m    = 33   !   2 meter (surf) pressure (hPa)
        integer, parameter  :: ktvis     = 34   !   visibility (km)
        integer, parameter  :: ktsdepth  = 35   !   snow depth (mm)
        integer, parameter  :: ktwcond   = 36   !   weather condition (none)
        integer, parameter  :: ktth_UppA = 37   !   potential temp (ps=1000hPa) (K)
        integer, parameter  :: ktskint   = 38   !   skin temperature (K)
        integer, parameter  :: ktSST     = 39   !   sea surface temperature (K)
        integer, parameter  :: ktTb      = 40   !   brightness temperature (K)
        integer, parameter  :: ktTv      = 44   !   layer-mean virtual temperature (K)

        integer, parameter :: ktLogAOD = 43     !   log-transformed AOD
        integer, parameter :: ktAOD  = 45       !   average Aerosol Optical Depth (none)
        integer, parameter :: ktANGE = 46       !   Angstrom Exponent
        integer, parameter :: ktAAOD = 49       !   Aerosol Absorption Optical Depth 
        integer, parameter :: ktSOLZ = 54       !   solar zenith  (degress)
        integer, parameter :: ktSOLA = 55       !   solar azimuth (degrees)
        integer, parameter :: ktSENZ = 56       !   sensor zenith (degress)
        integer, parameter :: ktSENA = 57       !   sensor azimuth (degress)
        integer, parameter :: ktREFL = 63       !   mean reflectance (none)
        integer, parameter :: ktN_RF = 65       !   counter for optical depth (none)
        integer, parameter :: ktpr2  = 86       !   GSI precip ln(1+rain rate) (mm/hr)
        integer, parameter :: kto3mx = 87       !   3d ozone mixing ratio (ppmv)
        integer, parameter :: ktGPSr = 88       !   gps (refractivity)
        integer, parameter :: ktGPSb = 89       !   gps (bending angle)
        integer, parameter :: ktlat  = 92       !   position-type obs: latitude
        integer, parameter :: ktlon  = 93       !   position-type obs: longitude
        integer, parameter :: ktdw   = 94       !   doppler wind lidar

        integer, dimension(28), parameter :: ktSurfAll = (/ ktus, ktvs, ktslp,
     &                                                      ktus10, ktTs10, ktTds, ktrhs, ktqs10,
     &                                                      ktspd2m, ktmxspd2m, ktwgust2m, ktt2m,
     &                                                      ktmxt2m, ktmmnt2m, ktdewt2m, ktrh2m,
     &                                                      ktsphu2m, ktps2m, ktskint, ktpr, kttpw,
     &                                                      ktSST, ktpr2, 
     &                                                      ktANGE, ktSOLZ, ktSOLA, ktSENZ, ktSENA /)

        integer, dimension(13), parameter :: ktUppaAll = (/ ktuu, ktvv, ktHH, ktww,
     &                                                      ktTT, ktTd, ktrh, ktqq,
     &                                                      ktuthk, ktth_UppA,
     &                                                      ktTb, ktTv, ktdw /)

        integer, dimension( 3), parameter :: ktChemAll = (/ kttco3, kto3, kto3mx /)

        integer, dimension( 8), parameter :: ktOthrAll = (/ kttlw, ktfcc, ktvis,
     &                                                      ktsdepth, ktwcond, ktgpsr, ktgpsb,
     &                                                      ktdw /)

!.

!                             ---------------------------
!                             Data Source (kx) Parameters
!                             ---------------------------
        integer, parameter ::   kxmax = 1024
        integer, parameter ::   kxmod = 1000   ! the smallest 10's power > kxmax 

!                             ------------------
!                             QC Flag parameters
!                             ------------------
!
!     History flags
!     -------------

      integer, parameter :: H_PRE_SUSP   =   1  ! unspecified preprocessing history flag
      integer, parameter :: H_UNDERG     =   3  ! suspect due to underground check

                                                ! NCEP quality marks that we treat as suspect
      integer, parameter :: H_NCEP1      =   4  !   NCEP CQC - modified value
      integer, parameter :: H_NCEP3      =   5  !   NCEP CQC - suspect value
      integer, parameter :: H_NCEP89     =   6  !   NCEP PREVENTS - P suspect
                                                !     (underground or sfcP too high)
                                                !      or suspect specific humid.
      integer, parameter :: H_NCEP14     =   7  !   NCEP SDM purged
      integer, parameter :: H_NCEP_OTHER =   8  !   NCEP missing or unknown QM
                                                !     (QM != 0,1,2,3,8,9,13,14,15)
                                                !  Suspect mark derived from
                                                !     pressure value (obs maybe OK,
                                                !     pressure is suspect)
      integer, parameter :: H_NCEP_PRES1 =   9  ! NCEP CQC modified P value
      integer, parameter :: H_NCEP_PRES3 =  10  ! NCEP CQC suspect P value
      integer, parameter :: H_NCEP_PRES14=  11  ! NCEP SDM purged P
                                                !  Suspect marks for moisture
                                                !  vars calcfrom suspect input values

      integer, parameter :: H_SUSP_TEMP    = 12 ! moisture from suspect temp
      integer, parameter :: H_SUSP_DEWTEMP = 13 ! moisture from suspect dewpt

      integer, parameter :: H_BACKG   = 17      ! background check

      integer, parameter :: H_YELLOW  = 20      ! obs marked as suspect by "Yellow List"
      integer, parameter :: H_SIMUL   = 21      ! obs confidence level less than 1
      integer, parameter :: H_BEAUFORT = 22     ! Beaufort corrected winds
      integer, parameter :: H_SUSP_LOCATION = 23  ! Suspect location winds


!     Exclusion flags
!     ---------------
      integer, parameter :: X_PRE_BAD     =   1 ! preprocessing exclusion flag
      integer, parameter :: X_BAD_LOC     =   2 ! bad location flag
      integer, parameter :: X_UNDERG      =   3 ! underground flag
      integer, parameter :: X_OBS_FILL    =   4 ! observation fill flag
      integer, parameter :: X_FCS_FILL    =   5 ! forecast at obs location fill flag
      integer, parameter :: X_TOO_HIGH    =   6 ! obs above certain desired level
      integer, parameter :: X_PASSIVE     =   7 ! passive data type exclusion flag
      integer, parameter :: X_TIME_ACT    =   8 ! data outside active time window
      integer, parameter :: X_NOT_ANAVAR  =   9 ! not an analysis variable


                                                ! NCEP quality marks treated as exclusion
      integer, parameter :: X_NCEP13      =  10 ! NCEP CQC bad observation
      integer, parameter :: X_NCEP15      =  11 ! NCEP PREPDATA bad observation
      integer, parameter :: X_NCEP_PRES13 =  12 ! NCEP CQC bad pressure
      integer, parameter :: X_NCEP_PRES15 =  13 ! NCEP PREPDATA bad pressure

                                                ! DAO preprocessing quality marks
      integer, parameter :: X_RANGE       =  14 ! DAO range check failed
      integer, parameter :: X_DUP         =  15 ! DAO duplicate obs. (>1 in 6 hr)
      integer, parameter :: X_HYDRO       =  16 ! DAO failed hydrostatic check
                                                !  Rejection marks for moisture values
                                                !    calculated from rejected/dubious input


      integer, parameter :: X_BUDDY       = 17  ! buddy check
      integer, parameter :: X_WIND        = 18  ! wind check
      integer, parameter :: X_NOSTATS     = 19  ! no error statistics
      integer, parameter :: X_PROFILE     = 20  ! profile check
      integer, parameter :: X_BACKG       = 21  ! background check
      integer, parameter :: X_BADXM       = 22  ! bad xm for scaling

      integer, parameter :: X_BADLAYER    = 25  ! layer too thin or too thick

      integer, parameter :: X_BAD_TEMP    = 28  ! moisture from bad temperature
      integer, parameter :: X_THIN        = 29  ! observation excl. by thinner

      integer, parameter :: X_UNPHYSICAL  = 30  ! obs with unphysical value
      integer, parameter :: X_RED         = 31  ! obs excluded by "Red List"
      integer, parameter :: X_SIMUL       = 32  ! obs could not be simulated
      integer, parameter :: X_PSAS        = 33  ! excluded by PSAS
      integer, parameter :: X_RANGE_ELEV  = 34  ! Failed elevation limit check
      integer, parameter :: X_NCEP_PURGE  = 35  ! NCEP SDM purged (non-raob)
      integer, parameter :: X_NCEP_NLNQC  = 36  ! Rejected by GSI non-linear QC 
      integer, parameter :: X_ODSMATCH    = 37  ! ODS-match could not find a match

!     Descriptions of history flags:
!     -----------------------------

                                ! number of history flags in use
      integer, parameter :: nqcHmax = 23

      character(len=32), parameter :: qcHnames(nqcHmax)=(/
     1                 'unspecified preprocessing flag  ',
     2                 '                                ',
     3                 'gcm slightly underground        ',
     4                 'NCEP CQC - modified value       ',
     5                 'NCEP CQC - suspect value        ',
     6                 'NCEP PREVENTS undergnd/bad P, q ',
     7                 'NCEP SDM purged                 ',
     8                 'NCEP unknown or missing QM      ',
     9                 'NCEP CQC modified P             ',
     +                 'NCEP CQC suspect P value        ',
     1                 'NCEP SDM purged P               ',
     2                 'moisture from suspect temp.     ',
     3                 'moisture from suspect dewpt.    ',
     4                 '                                ',
     5                 '                                ',
     6                 '                                ',
     7                 'outlier wrt background          ',
     8                 '                                ',
     9                 '                                ',
     +                 'marked suspect by Yellow list   ',
     1                 'obs simulated w/ confidence < 1 ',
     2                 'Beaufort Corrected Winds        ',
     3                 'Suspect location                '     /)


!     Descriptions of exclusion flags:
!     -------------------------------

                                ! number of exclusion flags in use
      integer, parameter :: nqcXmax = 37

      character(len=33), parameter :: qcXnames(nqcXmax)=(/
     1                 'unspecified preprocessing flag  ',
     2                 'impossible location             ',
     3                 'gcm deep underground            ',
     4                 'observation value undefined     ',
     5                 'forecast value undefined        ',
     6                 'observation level too high      ',
     7                 'passive data type               ',
     8                 'outside active time window      ',
     9                 'not an analysis variable        ',
     +                 'NCEP CQC bad observation        ',
     1                 'NCEP PREPDATA bad observation   ',
     2                 'NCEP CQC bad pressure           ',
     3                 'NCEP PREPDATA bad pressure      ',
     4                 'DAO range check failed          ',
     5                 'DAO duplicate obs. (>1 in 6 hr) ',
     6                 'DAO failed hydrostatic check    ',
     7                 'failed buddy check              ',
     8                 'incomplete wind vector          ',
     9                 'improper error statistics       ',
     +                 'incomplete vertical profile     ',
     1                 'extreme outlier wrt background  ',
     2                 'bad xm value for omf scaling    ',
     3                 '                                ',
     4                 '                                ',
     5                 'layer too thin or too thick     ',
     6                 '                                ',
     7                 '                                ',
     8                 'moisture from bad temp.         ',
     9                 'observation removed by thinner  ',
     +                 'unphysical value                ',
     1                 'excluded by "Red List"          ',
     2                 'obs could not be simulated      ',
     3                 'excluded by PSAS                ',
     4                 'Failed elevation limit check    ',
     5                 'NCEP SDM purged (non-raob)      ',
     6                 'Rejected by GSI non-linear QC   ',
     7                 'ODSmatch could not find match   '/)


      integer, parameter :: nsats = 38
      character(len=*), parameter :: sats(nsats)=(/
     .                 'hirs2           ', 'hirs3           ', 'hirs4           ',
     .                 'msu             ', 'ssu             ', 'sndr            ',
     .                 'sndrd1          ', 'sndrd2          ', 'sndrd3          ',
     .                 'sndrd4          ', 'amsua           ', 'amsub           ',
     .                 'mhs             ', 'airs            ', 'hsb             ',
     .                 'goes_img        ', 'avhrr           ', 'avhrr_navy      ',
     .                 'ssmi            ', 'amsre_low       ', 'amsre_mid       ',
     .                 'amsre_hig       ', 'ssmis           ', 'seviri          ',
     .                 'sbuv2           ', 'omi             ', 'iasi            ',
     .                 'atms            ', 'cris            ', 'omieff          ',
     .                 'o3lev           ', 'tomseff         ', 'gome            ',
     .                 'mls             ', 'mls20           ', 'mls22           ',
     .                 'mls30           ', 'mls55           '/)
! note: numbers below were made up for MHS, and SSU
! note: CRIS and ATMS numbers assigned at will
! note: omieff number assigned as omi
! note: mlsoz (o3lev) assigned were made up to 9999
      integer, parameter :: idsats(nsats)=(/
     .                 0                 , 0                 , 0                 ,
     .                 200               , 350               , 45                ,
     .                 50                , 60                , 70                ,
     .                 80                , 300               , 400               ,
     .                 800               , 49                , 450               ,
     .                 250               , 600               , 650               ,
     .                 700               , 547               , 548               ,
     .                 549               , 700               , 980               ,
     .                 450               , 449               , 850               ,
     .                 900               , 950               , 449               ,
     .                 304               , 440               , 445               ,
     .                 310               , 315               , 320               ,
     .                 325               , 330               /)

      integer, parameter :: npcp = 4
      character(len=*), parameter :: pcpt(npcp)=(/
     .                 'pcp_ssmi        ', 'pcp_tmi         ', 'pcp_amsu        ',
     .                 'pcp_stage3      '/)
                                                                                                                                             


!EOC




