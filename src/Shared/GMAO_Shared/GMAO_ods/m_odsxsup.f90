!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3   !
!------------------------------------------------------------------------
!BOP
!
! !MODULE: m_odsxsup -- Extra routines supporting external applications
!
! !INTERFACE:
!

      module m_odsxsup

! !USES:
      use m_ods

      implicit none

! !PUBLIC TYPES:
!
      PRIVATE
      PUBLIC getodsmeta
      PUBLIC ncepQCXval

!
! !DESCRIPTION:
!  Module to contain routines to interface with PREPBUFR
!
!
! !REVISION HISTORY:
!   28 July  2003    Sienkiewicz   convert subroutines to module
!    5 Aug   2003    Sienkiewicz   added interface code
!   28 Apr   2004    Todling       extracted from m_prep
!   19 Nov   2004    Dee           modified getodsmeta for NCEP report types
!
!EOP
!-------------------------------------------------------------------------

  CONTAINS

!-------------------------------------------------------------------------
!      NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3   !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  ncepQCXval() --- fill in values for NCEP QCexcl marks
!
! !INTERFACE:
!
      subroutine ncepQCXval( ods )

! !USES

      implicit NONE
!
! !INPUT/OUTPUT PARAMETERS:
!
      type(ods_vect), intent(inout) ::  ods ! ODS vector

!
! !DESCRIPTION: Fill in entries for NCEP QCexcl marks used by this program
!   Reference:
!    http://www.emc.ncep.noaa.gov/mmb/papers/keyser/prepbufr.doc/table_7.htm
!
!
! !REVISION HISTORY:
!  13 Mar 2002  Sienkiewicz   Initial version
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter ::  myname = 'ncepqcxval'

! default 32 chars             ....+....1....+....2....+....3..
      ods%meta%qcx_names    = ''
      ods%meta%qcx_names(4) = 'OIQC - prev QM was 0            '
      ods%meta%qcx_names(5) = 'OIQC - prev QM was 1            '
      ods%meta%qcx_names(6) = 'OIQC - prev QM was 2            '
      ods%meta%qcx_names(7) = 'OIQC - prev QM was 3            '
      ods%meta%qcx_names(8) = 'PREVENT                         '
      ods%meta%qcx_names(9) = 'PREVENT                         '
      ods%meta%qcx_names(10)= 'ACQC superob; PROFL failed chks '
      ods%meta%qcx_names(11)= 'PROFL no pass median/shear chk  '
      ods%meta%qcx_names(12)= 'Reject /  PROFL failed shear chk'
      ods%meta%qcx_names(13)= 'Fail Auto QC; PROFL faild median'
      ods%meta%qcx_names(14)= 'SDM purge flag                  '
      ods%meta%qcx_names(15)= 'PREPRO flag for non-use by analy'
      ods%meta%qcx_names(16)= 'Unknown NCEP QM                 '
      ods%meta%qcx_names(17)= 'Missing NCEP QM                 '

      return

      end subroutine ncepQCXval

!-------------------------------------------------------------------------
!      NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3   !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  getodsmeta() --- fill in ods metadata
!
! !INTERFACE:
!
      subroutine getodsmeta( ods )

! !USES:
!
      implicit NONE

!
! !INPUT/OUTPUT PARAMETERS:
!
      type(ods_vect), intent(inout) ::  ods ! ODS vector

!
! !DESCRIPTION: Fill in ODS metadata (KX, KT tables)
!
! !REVISION HISTORY:
!  11 Mar 2002  Sienkiewicz   Initial version - incorporates modified
!                              rdkxtbl and rdkttbl from PSAS
!
!  24 Jul 2003  Sienkiewicz   Rewrite with hard-coded values for KT
!                             hardcoded NCEP 'report types' for KX
!   reference http://www.emc.ncep.noaa.gov/mmb/papers/keyser/prepbufr.doc/table_2.htm
!                             Mapping based on 'merger' between mass and
!                             wind types.
!  19 Nov 2004  Dee           Defined kx_names according to NCEP usage
!                             (see notes below)
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter ::  myname = 'getodsmeta'

      character*60 errstring

      integer iret

!     integer, parameter :: nkxmax = 287
!     character(len=150) kxnames(nkxmax)

      ods%meta%kx_meta=''
      ods%meta%kx_names=''

      ods%meta%kx_names(111) = 'SYNTHETIC TROPICAL CYCLONE STORM CENTER PSFC, Q'
      ods%meta%kx_names(112) = 'PSEUDO MEAN SEA-LEVEL PRESSURE AT TROPICAL CYCLONE STORM CENTER'
      ods%meta%kx_names(120) = 'RAWINSONDE VIRTUAL TEMPERATURE, SPECIFIC HUMIDITY, STATION PRESSURE'
      ods%meta%kx_names(122) = 'CLASS SOUNDING VIRTUAL TEMPERATURE, SPECIFIC HUMIDITY, STATION PRESSURE'
      ods%meta%kx_names(126) = 'RASS (FROM NPN OR MAP NETWORK)'
      ods%meta%kx_names(130) = 'AIREP AND PIREP AIRCRAFT SENSIBLE TEMPERATURE'
      ods%meta%kx_names(131) = 'AMDAR AIRCRAFT SENSIBLE TEMPERATURE'
      ods%meta%kx_names(132) = 'FLIGHT-LEVEL RECONNAISSANCE AND PROFILE DROPSONDE VIRTUAL TEMPERATURE, '// &
                               'SPECIFIC HUMIDITY, STATION PRESSURE'
      ods%meta%kx_names(133) = 'MDCARS AIRCRAFT SENSIBLE TEMPERATURE (SPECIFIC HUMIDITY FLAGGED FOR NON-USE BY ANALYSIS)'
      ods%meta%kx_names(134) = 'TAMDAR AIRCRAFT SENSIBLE TEMPERATURE, Q'
      ods%meta%kx_names(135) = 'CANADIAN AMDAR AIRCRAFT SENSIBLE TEMPERATURE'
      ods%meta%kx_names(150) = 'SSM/I SUPEROBED FNOC RAIN RATE (DMSP-13, DMSP-15)'
      ods%meta%kx_names(151) = 'NESDIS SFOV CLOUD TOP PRESSURE AND TEMPERATURE, CLOUD AMOUNT (GOES)'
      ods%meta%kx_names(152) = 'SSM/I SUPEROBED NEURAL NET 3 TOTAL PRECIPITABLE WATER RETRIEVALS (DMSP-13, DMSP-15)'
      ods%meta%kx_names(153) = 'GPS-INTEGRATED PRECIPITABLE WATER (GPS-IPW)'
      ods%meta%kx_names(156) = 'NESDIS LAYER PRECIPITABLE WATER RETRIEVALS OVER LAND - CLEAR (GOES)'
      ods%meta%kx_names(157) = 'NESDIS LAYER PRECIPITABLE WATER RETRIEVALS OVER LAND - CLOUD CORRECTED (GOES)'
      ods%meta%kx_names(158) = 'NESDIS LAYER PRECIPITABLE WATER RETRIEVALS OVER OCEAN - CLEAR (GOES)'
      ods%meta%kx_names(159) = 'NESDIS LAYER PRECIPITABLE WATER RETRIEVALS OVER OCEAN - CLOUD CORRECTED (GOES)'
      ods%meta%kx_names(164) = 'NESDIS RADIANCES OVER LAND - CLEAR (GOES)'
      ods%meta%kx_names(165) = 'NESDIS RADIANCES OVER LAND - CLOUD CORRECTED (GOES)'
      ods%meta%kx_names(174) = 'NESDIS RADIANCES OVER OCEAN - CLEAR (GOES)'
      ods%meta%kx_names(175) = 'NESDIS RADIANCES OVER OCEAN - CLOUD CORRECTED (GOES)'
      ods%meta%kx_names(180) = 'SURFACE MARINE (SHIP, BUOY, C-MAN) VIRTUAL TEMPERATURE, SPECIFIC HUMIDITY, '// &
                               'STATION PRESSURE (STATION PRESSURE REPORTED)'
      ods%meta%kx_names(181) = 'SURFACE LAND SYNOPTIC AND METAR STATION PRESSURE, SPECIFIC HUMIDITY (TEMPERATURE '// &
                               'NOT USED BY ANALYSIS) (STATION PRESSURE REPORTED)'
      ods%meta%kx_names(182) = 'SPLASH LEVEL VIRTUAL TEMPERATURE, SPECIFIC HUMIDITY, STATION PRESSURE (OVER OCEAN ONLY)'
      ods%meta%kx_names(183) = 'SURFACE MARINE (SHIP, BUOY, C-MAN), LAND SYNOPTIC AND METAR VIRTUAL TEMPERATURE, '// &
                               'SPECIFIC HUMIDITY, STATION PRESSURE (STATION PRESSURE NOT REPORTED)'
      ods%meta%kx_names(187) = 'SURFACE METAR VIRTUAL TEMPERATURE, SPECIFIC HUMIDITY, STATION PRESSURE (STATION '// &
                               'PRESSURE NOT REPORTED)'
      ods%meta%kx_names(190) = 'OPC/NOS POINT MEAN SEA-LEVEL PRESSURE BOGUS'
      ods%meta%kx_names(191) = 'AUSTRALIAN PAOB MEAN SEA-LEVEL PRESSURE BOGUS'
      ods%meta%kx_names(199) = 'DRIFTING BUOY SEA-LEVEL PRESSURE AND TEMPERATURE'
      ods%meta%kx_names(210) = 'SYNTHETIC TROPICAL CYCLONE U, V'
      ods%meta%kx_names(211) = 'TRMM RAIN RATE'
      ods%meta%kx_names(220) = 'RAWINSONDE U, V'
      ods%meta%kx_names(221) = 'PIBAL U, V'
      ods%meta%kx_names(222) = 'CLASS SOUNDING U, V'
      ods%meta%kx_names(223) = 'NOAA PROFILER NETWORK (NPN) WIND PROFILER U, V'
      ods%meta%kx_names(224) = 'NEXRAD VERTICAL AZIMUTH DISPLAY (VAD) U, V'
      ods%meta%kx_names(225) = 'NEXRAD RADIAL U, V'
      ods%meta%kx_names(228) = 'JAPANESE (JMA) WIND PROFILER U, V'
      ods%meta%kx_names(229) = 'WIND PROFILER FROM PILOT (PIBAL) BULLETINS U, V'
      ods%meta%kx_names(230) = 'AIREP AND PIREP AIRCRAFT U, V'
      ods%meta%kx_names(231) = 'AMDAR AIRCRAFT U, V'
      ods%meta%kx_names(232) = 'FLIGHT-LEVEL RECONNAISSANCE AND PROFILE DROPSONDE U, V'
      ods%meta%kx_names(233) = 'MDCARS AIRCRAFT U, V'
      ods%meta%kx_names(234) = 'TAMDAR AIRCRAFT U, V'
      ods%meta%kx_names(235) = 'CANADIAN AMDAR AIRCRAFT U, V'
      ods%meta%kx_names(241) = 'INDIA IR AND VISIBLE CLOUD DRIFT U, V (INSAT-2E)'
      ods%meta%kx_names(242) = 'JMA IR AND VISIBLE CLOUD DRIFT U, V AT LEVELS BELOW 850 MB (GMS, MTSAT)'
      ods%meta%kx_names(243) = 'EUMETSAT IR AND VISIBLE CLOUD DRIFT U, V AT LEVELS BELOW 850 MB (METEOSAT)'
      ods%meta%kx_names(244) = 'AVHRR/POES IR CLOUD DRIFT (NOAA, METOP) U,V'
      ods%meta%kx_names(245) = 'NESDIS IR CLOUD DRIFT U, V (GOES)'
      ods%meta%kx_names(246) = 'NESDIS IMAGER WATER VAPOR CLOUD U, V AT CLOUD TOP (GOES)'
      ods%meta%kx_names(247) = 'NESDIS IMAGER WATER VAPOR CLOUD U, V - DEEP LAYER (GOES)'
      ods%meta%kx_names(248) = 'NESDIS SOUNDER WATER VAPOR CLOUD U, V AT CLOUD TOP (GOES)'
      ods%meta%kx_names(249) = 'NESDIS SOUNDER WATER VAPOR CLOUD U, V - DEEP LAYER (GOES)'
      ods%meta%kx_names(250) = 'JMA WATER VAPOR CLOUD U, V (GMS, MTSAT)'
      ods%meta%kx_names(251) = 'NESDIS VISIBLE CLOUD DRIFT U, V (GOES)'
      ods%meta%kx_names(252) = 'JMA IR AND VISIBLE CLOUD DRIFT U, V AT LEVELS ABOVE 850 MB (GMS, MTSAT)'
      ods%meta%kx_names(253) = 'EUMETSAT IR AND VISIBLE CLOUD DRIFT U, V AT LEVELS ABOVE 850 MB (METEOSAT)'
      ods%meta%kx_names(254) = 'EUMETSAT WATER VAPOR CLOUD U, V (METEOSAT)'
      ods%meta%kx_names(255) = 'NESDIS PICTURE TRIPLET CLOUD U, V (GOES)'
      ods%meta%kx_names(256) = 'INDIA WATER VAPOR CLOUD U, V (INSAT-2E)'
      ods%meta%kx_names(257) = 'NASA/MODIS POES IR CLOUD-DRIFT U, V'
      ods%meta%kx_names(258) = 'NASA/MODIS POES WATER VAPOR IMAGER CLOUD TOP U, V'
      ods%meta%kx_names(259) = 'NASA/MODIS POES WATER VAPOR IMAGER DEEP LAYER U, V'
      ods%meta%kx_names(264) = 'SSM/I SUPEROBED RAIN RATE'
      ods%meta%kx_names(280) = 'SURFACE MARINE (SHIP, BUOY, C-MAN) U, V (STATION PRESSURE REPORTED)'
      ods%meta%kx_names(281) = 'SURFACE LAND SYNOPTIC AND METAR U, V (STATION PRESSURE REPORTED)'
      ods%meta%kx_names(282) = 'ATLAS BUOY U, V'
      ods%meta%kx_names(283) = 'SSM/I SUPEROBED NEURAL NET 3 WIND SPEED (DIRECTION SET TO ZERO, SPEED ASSIMILATED '// &
                               'DIRECTLY) (DMSP-13, DMSP-15)'
      ods%meta%kx_names(284) = 'SURFACE MARINE (SHIP, BUOY, C-MAN), LAND SYNOPTIC AND METAR U, V (STATION PRESSURE NOT REPORTED)'
      ods%meta%kx_names(285) = 'QUIKSCAT SCATTEROMETER U, V'
      ods%meta%kx_names(286) = 'ERS-2 SCATTEROMETER U, V'
      ods%meta%kx_names(287) = 'SURFACE METAR U, V (STATION PRESSURE NOT REPORTED)'
      ods%meta%kx_names(289) = 'SUPEROBED WINDSAT SCATTEROMETER U,V OVER OCEAN'
      ods%meta%kx_names(290) = 'NON-SUPEROBED ASCAT SCATTEROMETER U,V OVER OCEAN'
      ods%meta%kx_names(291) = 'NON-SUPEROBED OSCAT SCATTEROMETER U,V OVER OCEAN'
      ods%meta%kx_names(299) = 'DRIFTING BUOY U,V'


      ods%meta%kt_names=''
      ods%meta%kt_units=''
      ods%meta%kt_names (  1 ) = 'Surface (10m) zonal wind '
      ods%meta%kt_units (  1 ) = 'm/sec '
      ods%meta%kt_names (  2 ) = 'Surface (10m) meridional wind '
      ods%meta%kt_units (  2 ) = 'm/sec '
      ods%meta%kt_names (  3 ) = 'Sea level pressure '
      ods%meta%kt_units (  3 ) = 'hPa '
      ods%meta%kt_names (  4 ) = 'Upper-air zonal wind '
      ods%meta%kt_units (  4 ) = 'm/sec '
      ods%meta%kt_names (  5 ) = 'Upper-air meridional wind '
      ods%meta%kt_units (  5 ) = 'm/sec '
      ods%meta%kt_names (  6 ) = 'Upper-air geopotential height '
      ods%meta%kt_units (  6 ) = 'm '
      ods%meta%kt_names (  7 ) = 'Upper-air water vapor mixing ratio '
      ods%meta%kt_units (  7 ) = 'g/kg '
      ods%meta%kt_names (  8 ) = 'Upper-air temperature '
      ods%meta%kt_units (  8 ) = 'Kelvin '
      ods%meta%kt_names (  9 ) = 'Upper-air dew-point temperature '
      ods%meta%kt_units (  9 ) = 'Kelvin '
      ods%meta%kt_names ( 10 ) = 'Upper-air relative humidity '
      ods%meta%kt_units ( 10 ) = '% '
      ods%meta%kt_names ( 11 ) = 'Upper-air specific humidity '
      ods%meta%kt_units ( 11 ) = 'g/kg '
      ods%meta%kt_names ( 12 ) = 'Surface (10m) wind speed '
      ods%meta%kt_units ( 12 ) = 'm/sec'
      ods%meta%kt_names ( 13 ) = 'Surface (10m) temperature '
      ods%meta%kt_units ( 13 ) = 'Kelvin '
      ods%meta%kt_names ( 14 ) = 'Surface (10m) dew-point temperature '
      ods%meta%kt_units ( 14 ) = 'Kelvin '
      ods%meta%kt_names ( 15 ) = 'Surface (10m) relative humidity '
      ods%meta%kt_units ( 15 ) = '% '
      ods%meta%kt_names ( 16 ) = 'Surface (10m) specific humidity '
      ods%meta%kt_units ( 16 ) = 'g/kg '
      ods%meta%kt_names ( 17 ) = 'Precipitation rate '
      ods%meta%kt_units ( 17 ) = 'mm/day '
      ods%meta%kt_names ( 18 ) = 'Total precipitable water '
      ods%meta%kt_units ( 18 ) = 'mm '
      ods%meta%kt_names ( 19 ) = 'Total cloud liquid water '
      ods%meta%kt_units ( 19 ) = 'mm '
      ods%meta%kt_names ( 20 ) = 'Fractional cloud cover '
      ods%meta%kt_units ( 20 ) = '% '
      ods%meta%kt_names ( 21 ) = 'Total column ozone'
      ods%meta%kt_units ( 21 ) = 'Dobson '
      ods%meta%kt_names ( 22 ) = 'Ozone '
      ods%meta%kt_units ( 22 ) = 'Dobson '
      ods%meta%kt_names ( 23 ) = 'Upper-air thickness'
      ods%meta%kt_units ( 23 ) = 'm'
      ods%meta%kt_names ( 33 ) = 'Surface pressure'
      ods%meta%kt_units ( 33 ) = 'hPa'
      ods%meta%kt_names ( 39 ) = 'Sea Surface Temperature'
      ods%meta%kt_units ( 40 ) = 'K'
      ods%meta%kt_names ( 40 ) = 'Brightness Temperature'
      ods%meta%kt_units ( 40 ) = 'K'
      ods%meta%kt_names ( 41 ) = 'Mean layer virtual temperature'
      ods%meta%kt_units ( 41 ) = 'K'
      ods%meta%kt_names ( 44 ) = 'Virtual temperature (single lvl)'
      ods%meta%kt_units ( 44 ) = 'K'
      ods%meta%kt_names ( 86 ) = 'GSI rain rate log(1+rate)'
      ods%meta%kt_units ( 86 ) = 'log (mm/hr)'
      ods%meta%kt_names ( 88 ) = 'Refractivity'
      ods%meta%kt_units ( 88 ) = 'N-units'
      ods%meta%kt_names ( 89 ) = 'Bending Angle'
      ods%meta%kt_units ( 89 ) = 'NONE'

      return

      end subroutine getodsmeta


!$$$Notes:  Issues to be discussed
!$$$
!$$$ 1    PREPBUFR data is written in separate mass and wind reports.  Do we
!$$$      want to merge the mass and wind reports back together (give mass and
!$$$      wind obs from the same report the same 'ks').  Currently the routine
!$$$      is written to do this.  NOTE: Steve Bloom says yes
!$$$
!$$$ 2    NCEP report types do not map one-to-one to our KX table - how do we
!$$$      want to handle this?  It would be possible to:
!$$$        (a) Use NCEP report types - but then must merge mass/wind types if
!$$$            we merge mass and wind reports.
!$$$  or    (b) Try to define an approximate mapping from NCEP report types to
!$$$            our KX definitions.
!$$$ 3    May need addtional clarification of data to be extracted from PREPBUFR.
!$$$         e.g   Spec. hum or convert to mixing ratio?
!$$$         Currently write most as NCEP variables - "converted" rain rate

!$$$  Did not convert virtual temperature - maybe want to make new KT?
!$$$  New - convert virtual temperature to sensible temperature when Q available
!
!$$$  19 Nov 2004 (Dick): 
!$$$      Decided to implement 2(a): i.e. use NCEP report types in KX table, but 
!$$$      keeping our own KT table (with additions as needed). 
!
!  how to handle surface data?
!    SLP  -  the slp is observation,  what value for level? decided not to put
!    sfcP    -  sfcP could be observation with ELV as level?   ???
!  other vars-  use sfcP as level (if available) or ELV or what...
!     decided to use sfcP - as this is used in NCEP SSI as the level.
!
!
      end module m_odsxsup
