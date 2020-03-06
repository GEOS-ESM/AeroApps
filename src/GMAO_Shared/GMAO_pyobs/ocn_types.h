c
c  NRL Coupled Ocean Data Assimilation (NCODA) Data Types
c
c   0 = All Data Combined
c   1 = Bathy Temperatures (C)
c   2 = NOAA14 Day GAC SST (C)
c   3 = SHIP Engine Room Intake (C)
c   4 = Fixed BUOY Temperature (C)
c   5 = Drifting BUOY (C)
c   6 = NOAA14 Night GAC SST (C)
c   7 = NOAA14 Relaxed Day GAC SST (C)
c   8 = SSM/I F11 Ice (%)
c   9 = SSM/I F13 Ice (%)
c  10 = SSM/I F14 Ice (%)
c  11 = Supplemental Ice (%)
c  12 = Topex (M)
c  13 = ERS2 (M)
c  14 = GFO (M)
c  15 = MODAS Temperature (C)
c  16 = GDEM 3D Climatology (C)
c  17 = GOES8 Day SST (C)
c  18 = GOES8 Night SST (C)
c  19 = Direct Method Temperature (C)
c  20 = TESAC Temperature (C)
c  21 = SHIP Bucket (C)
c  22 = SHIP Hull Sensor (C)
c  23 = CMAN SST (C)
c  24 = NOAA15 Day GAC SST (C)
c  25 = NOAA15 Night GAC SST (C)
c  26 = NOAA15 Relaxed Day GAC SST (C)
c  27 = Mechanical BT (C)
c  28 = Hydrocast BT (C)
c  29 = SSM/I F15 Ice (%)
c  30 = In Situ Sea Surface Height Anomaly (M)
c  31 = SSM/I Ice Shelf
c  32 = TESAC Salinity (PSU)
c  33 = MODAS Salinity (PSU)
c  34 = TRACK OB Temperature (C)
c  35 = TRACK OB Salinty (PSU)
c  36 = Argo Float Temperature (C)
c  37 = Argo Float Salinity (PSU)
c  38 = Supplemental MODAS Temperature (C)
c  39 = Supplemental MODAS Salinity (PSU)
c  40 = Supplemental Sea Surface Height Anomaly (M)
c  41 = Freezing Sea Water SST (C)
c  42 = SST super ob (C)
c  43 = NOAA16 Day GAC SST (C)
c  44 = NOAA16 Night GAC SST (C)
c  45 = NOAA16 Relaxed Day GAC SST (C)
c  46 = SST Derived Surface Salinity (PSU)
c  47 = GOES10 Day SST (C)
c  48 = GOES10 Night SST (C)
c  49 = Direct Method Salinity (PSU)
c  50 = Extended Temperatures (C)
c  51 = Extended Salinity (PSU)
c  52 = Fixed BUOY Salinity (PSU)
c  53 = Jason-1 SSH (M)
c  54 = Drifting BUOY Salinity (PSU)
c  55 = Envisat SSH (M)
c  56 = NOAA17 Day GAC SST (C)
c  57 = NOAA17 Night GAC SST (C)
c  58 = NOAA17 Relaxed Day GAC SST (C)
c  59 = NOAA16 Day LAC SST (C)
c  60 = NOAA16 Night LAC SST (C)
c  61 = NOAA17 Day LAC SST (C)
c  62 = NOAA17 Night LAC SST (C)
c  63 = ERS2 Significant Wave Height (M)
c  64 = Topex Significant Wave Height (M)
c  65 = GFO Significant Wave Height (M)
c  66 = Jason-1 Significant Wave Height (M)
c  67 = Envisat Significant Wave Height (M)
c  68 = Fixed BUOY Significant Wave Height (M)
c  69 = AMSRE Day Microwave SST (C)
c  70 = GOES12 Day SST (C)
c  71 = GOES12 Night SST (C)
c  72 = AMSRE Night Microwave SST (C)
c  73 = TRMM microwave SST (C)
c  74 = AATSR (Envisat) Day SST (C)
c  75 = AATSR (Envisat) Night SST (C)
c  76 = Topex Interleaved SSH (M)
c  77 = AMSR2 Night Microwave SST (C)
c  78 = SSM/I Sea Ice Super Observations
c  79 = SST Super Observations
c  80 = SSMIS F17 Ice (%)
c  81 = Altimeter SSH Super Observations
c  82 = Altimeter SWH Super Observations
c  83 = Near Shore Ice (%)
c  84 = Aircraft Sea Surface Temperature (C)
c  85 = HF Radar U Velocity Component
c  86 = HF Radar V Velocity Component
c  87 = U Velocity
c  88 = V Velocity
c  89 = Glider Absolute U Velocity Component
c  90 = SSMIS F18 Ice (%)
c  91 = Jason-2 SSH (M)
c  92 = Glider Absolute V Velocity Component
c  93 = Surface Drifter U Velocity Component
c  94 = NOAA18 Day GAC SST (C)
c  95 = NOAA18 Night GAC SST (C)
c  96 = NOAA18 Relaxed Day GAC SST (C)
c  97 = NOAA18 Day LAC SST (C)
c  98 = NOAA18 Night LAC SST (C)
c  99 = MSG Day SST (C)
c 100 = MSG Night SST (C)
c 101 = MSG Day/Night SST (C)
c 102 = Glider Temperature Profiles (C)
c 103 = Glider Salinity Profiles (PSU)
c 104 = Surface Drifter V Velocity Component
c 105 = ADCP U Velocity Component
c 106 = HYCOM Layer Pressure (db)
c 107 = GOES11 Day SST (C)
c 108 = GOES11 Night SST (C)
c 109 = SSMIS F16 Ice (%)
c 110 = METOP-A Day GAC SST (C)
c 111 = METOP-A Night GAC SST (C)
c 112 = METOP-A Relaxed Day GAC SST (C)
c 113 = METOP-A Day LAC SST (C)
c 114 = METOP-A Night LAC SST (C)
c 115 = METOP-B Day GAC SST (C)
c 116 = METOP-B Night GAC SST (C)
c 117 = METOP-B Relaxed Day GAC SST (C)
c 118 = METOP-B Day LAC SST (C)
c 119 = METOP-B Night LAC SST (C)
c 120 = METOP-C Day GAC SST (C)
c 121 = METOP-C Night GAC SST (C)
c 122 = METOP-C Relaxed Day GAC SST (C)
c 123 = METOP-C Day LAC SST (C)
c 124 = METOP-C Night LAC SST (C)
c 125 = Jason-2 SWH (M)
c 126 = Not Used
c 127 = Jason-1 Interleaved SSH (M)
c 128 = NOAA19 Day GAC SST (C)
c 129 = NOAA19 Night GAC SST (C)
c 130 = NOAA19 Relaxed Day GAC SST (C)
c 131 = NOAA19 Day LAC SST (C)
c 132 = NOAA19 Night LAC SST (C)
c 133 = Anmimal Borne Temperature (C)
c 134 = Animal Borne Salinity (PSU)
c 135 = GOES13 Day SST (C)
c 136 = GOES13 Night SST (C)
c 137 = ADCP V Velocity Component
c 138 = AMSR2 Sea Ice (%)
c 139 = MTSAT-2 Day SST (C)
c 140 = MTSAT-2 Night SST (C)
c 141 = Expanded Real Profile Temperature (C)
c 142 = Expanded Real Profile Salinity (PSU)
c 143 = Expanded Synthetic Profile Temperature (C)
c 144 = Expanded Synthetic Profile Salinity (PSU)
c 145 = Envisat Low Orbit SSH (M)
c 146 = Physical SST (C)
c 147 = MODIS (Aqua) Ice Surface Temperature (C)
c 148 = MODIS (Terra) Ice Surface Temperature (C)
c 149 = Ice Surface Temp Super Observations
c 150 = METOP-A Ice Surface Temperature (C)
c 151 = WindSat Day Microwave SST (C)
c 152 = WindSat Night Microwave SST (C)
c 153 = NPP VIIRS Day SST (C)
c 154 = NPP VIIRS Night SST (C)
c 155 = NPP VIIRS Relaxed Day SST (C)
c 156 = NPP VIIRS Ice Surface Temperature (C)
c 157 = Argo Trajectory U Velocity (m/s)
c 158 = Argo Trajectory V Velocity (m/s)
c 159 = GOES15 Day SST (C)
c 160 = GOES15 Night SST (C)
c 161 = Altimeter SSH Cross-Track Geostrophic U Velocity (m/s)
c 162 = Altimeter SSH Cross-Track Geostrophic V Velocity (m/s)
c 163 = Wave Glider Temperature (C)
c 164 = Wave Glider Salinity (PSU)
c 165 = WindSat Sea Ice (%)
c 166 = Cryosat-2 Significant Wave Height (M)
c 167 = Cryosat-2 Sea Surface Height (M)
c 168 = Jason-1 Geodetic Orbit Sea Surface Height (M)
c 169 = AMSR2 Day Microwave SST (C)
c 170 = Aquarius Salinity (PSU)
c 171 = ISOP Temperature (C)
c 172 = ISOP Salinity (PSU)
c 173 = SSS Super Observations
c 174 = SMOS Salinity (PSU)
c 175 = GOES14 Day SST (C)
c 176 = GOES14 Night SST (C)
c 177 = Expanded U Velocity Profile (m/s)
c 178 = Expanded V Velocity Profile (m/s)
c 179 = Altika Sea Surface Height (M)
c 180 = Altika Significant Wave Height (M)
c 181 = HIMAWARI-8 Day SST (C)
c 182 = HIMAWARI-8 Night SST (C)
c 183 = Expanded ISOP Profile Temperature (C)
c 184 = Expanded ISOP Profile Salinity (PSU)
c 185 = SSMIS F19 Ice (%)
c 186 = Expendable Conductivity Temperature (C)
c 187 = Expendable Conductivity Salinity (PSU)
c 188 = ALAMO Float Temperature (C)
c 189 = ALAMO Float Salinity (PSU)
c 190 = ASTER Terra Day SST (C)
c 191 = ASTER Terra Night SST (C)
c 192 = VIIRS Suomi/NPP Ice (%)
c 193 = Jason-2 Interleaved SSH (M)
c 194 = Jason-2 Interleaved SWH (M)
c 195 = Jason-3 SSH (M)
c 196 = Jason-3 SWH (M)
c 197 = Sentinel-3A SSH (M)
c 198 = Sentinel-3A SWH (M)
c 199 = Sentinel-3B SSH (M)
c 200 = Sentinel-3B SWH (M)
c 201 = Altika Drifting Phase Sea Surface Height (M)
c 202 = Altika Drifting Phase Significant Wave Height (M)
c
c     Name                   Description
c   ---------         ---------------------------
c   MX_TYPES          maximum number data types
c   data_lbl          data type labels
c
      integer    MX_TYPES
      parameter (MX_TYPES = 202)
c
      character data_lbl (0:MX_TYPES) * 20
      data      data_lbl(0)  / '   All Data Combined' /
      data      data_lbl(1)  / '       eXpendable BT' /
      data      data_lbl(2)  / '     N14 GAC Day SST' /
      data      data_lbl(3)  / '            ERI SHIP' /
      data      data_lbl(4)  / '     Fixed BUOY Temp' /
      data      data_lbl(5)  / '  Drifting BUOY Temp' /
      data      data_lbl(6)  / '   N14 GAC Night SST' /
      data      data_lbl(7)  / '  N14 GAC RlxDay SST' /
      data      data_lbl(8)  / '       SSM/I F11 Ice' /
      data      data_lbl(9)  / '       SSM/I F13 Ice' /
      data      data_lbl(10) / '       SSM/I F14 Ice' /
      data      data_lbl(11) / '    Supplemental Ice' /
      data      data_lbl(12) / '           Topex SSH' /
      data      data_lbl(13) / '            ERS2 SSH' /
      data      data_lbl(14) / '             GFO SSH' /
      data      data_lbl(15) / '   MODAS Temperature' /
      data      data_lbl(16) / '    GDEM Climate SST' /
      data      data_lbl(17) / '       GOES8 Day SST' /
      data      data_lbl(18) / '     GOES8 Night SST' /
      data      data_lbl(19) / '  Direct Temperature' /
      data      data_lbl(20) / '   TESAC Temperature' /
      data      data_lbl(21) / '         Bucket SHIP' /
      data      data_lbl(22) / '    Hull Sensor SHIP' /
      data      data_lbl(23) / '            CMAN SST' /
      data      data_lbl(24) / '     N15 GAC Day SST' /
      data      data_lbl(25) / '   N15 GAC Night SST' /
      data      data_lbl(26) / '  N15 GAC RlxDay SST' /
      data      data_lbl(27) / '       Mechanical BT' /
      data      data_lbl(28) / '        Hydrocast BT' /
      data      data_lbl(29) / '       SSM/I F15 Ice' /
      data      data_lbl(30) / '         In Situ SSH' /
      data      data_lbl(31) / '     SSM/I Shelf Ice' /
      data      data_lbl(32) / '      TESAC Salinity' /
      data      data_lbl(33) / '      MODAS Salinity' /
      data      data_lbl(34) / '       TRACK OB Temp' /
      data      data_lbl(35) / '       TRACK OB Salt' /
      data      data_lbl(36) / '     Argo Float Temp' /
      data      data_lbl(37) / '     Argo Float Salt' /
      data      data_lbl(38) / '    MODAS Suppl Temp' /
      data      data_lbl(39) / '    MODAS Suppl Salt' /
      data      data_lbl(40) / '    Supplemental SSH' /
      data      data_lbl(41) / '         Sea Ice SST' /
      data      data_lbl(42) / '                 SST' /
      data      data_lbl(43) / '     N16 GAC Day SST' /
      data      data_lbl(44) / '   N16 GAC Night SST' /
      data      data_lbl(45) / '  N16 GAC RlxDay SST' /
      data      data_lbl(46) / '        SST Sfc Salt' /
      data      data_lbl(47) / '      GOES10 Day SST' /
      data      data_lbl(48) / '    GOES10 Night SST' /
      data      data_lbl(49) / '     Direct Salinity' /
      data      data_lbl(50) / '       Extended Temp' /
      data      data_lbl(51) / '       Extended Salt' /
      data      data_lbl(52) / '     Fixed BUOY Salt' /
      data      data_lbl(53) / '         Jason-1 SSH' /
      data      data_lbl(54) / '  Drifting BUOY Salt' /
      data      data_lbl(55) / '         Envisat SSH' /
      data      data_lbl(56) / '     N17 GAC Day SST' /
      data      data_lbl(57) / '   N17 GAC Night SST' /
      data      data_lbl(58) / '  N17 GAC RlxDay SST' /
      data      data_lbl(59) / '     N16 LAC Day SST' /
      data      data_lbl(60) / '   N16 LAC Night SST' /
      data      data_lbl(61) / '     N17 LAC Day SST' /
      data      data_lbl(62) / '   N17 LAC Night SST' /
      data      data_lbl(63) / '            ERS2 SWH' /
      data      data_lbl(64) / '           Topex SWH' /
      data      data_lbl(65) / '             GFO SWH' /
      data      data_lbl(66) / '         Jason-1 SWH' /
      data      data_lbl(67) / '         Envisat SWH' /
      data      data_lbl(68) / '      Fixed BUOY SWH' /
      data      data_lbl(69) / '       AMSRE Day SST' /
      data      data_lbl(70) / '      GOES12 Day SST' /
      data      data_lbl(71) / '    GOES12 Night SST' /
      data      data_lbl(72) / '     AMSRE Night SST' /
      data      data_lbl(73) / '         TRMM MW SST' /
      data      data_lbl(74) / '       AATSR Day SST' /
      data      data_lbl(75) / '     AATSR Night SST' /
      data      data_lbl(76) / '   Topex Intrlvd SSH' /
      data      data_lbl(77) / '     AMSR2 Night SST' /
      data      data_lbl(78) / '             Sea Ice' /
      data      data_lbl(79) / '                 SST' /
      data      data_lbl(80) / '       SSMIS F17 Ice' /
      data      data_lbl(81) / '                 SSH' /
      data      data_lbl(82) / '                 SWH' /
      data      data_lbl(83) / '      Near Shore Ice' /
      data      data_lbl(84) / '        Aircraft SST' /
      data      data_lbl(85) / '      HF Radar U Vel' /
      data      data_lbl(86) / '      HF Radar V Vel' /
      data      data_lbl(87) / '          U Velocity' /
      data      data_lbl(88) / '          V Velocity' /
      data      data_lbl(89) / '        Glider U Vel' /
      data      data_lbl(90) / '       SSMIS F18 Ice' /
      data      data_lbl(91) / '         Jason-2 SSH' /
      data      data_lbl(92) / '        Glider V Vel' /
      data      data_lbl(93) / '       Drifter U Vel' /
      data      data_lbl(94) / '     N18 GAC Day SST' /
      data      data_lbl(95) / '   N18 GAC Night SST' /
      data      data_lbl(96) / '  N18 GAC RlxDay SST' /
      data      data_lbl(97) / '     N18 LAC Day SST' /
      data      data_lbl(98) / '   N18 LAC Night SST' /
      data      data_lbl(99) / '         MSG Day SST' /
      data      data_lbl(100) / '       MSG Night SST' /
      data      data_lbl(101) / '   MSG Day/Night SST' /
      data      data_lbl(102) / '         Glider Temp' /
      data      data_lbl(103) / '         Glider Salt' /
      data      data_lbl(104) / '       Drifter V Vel' /
      data      data_lbl(105) / '          ADCP U Vel' /
      data      data_lbl(106) / '      Layer Pressure' /
      data      data_lbl(107) / '      GOES11 Day SST' /
      data      data_lbl(108) / '    GOES11 Night SST' /
      data      data_lbl(109) / '       SSMIS F16 Ice' /
      data      data_lbl(110) / '   MET-A GAC Day SST' /
      data      data_lbl(111) / ' MET-A GAC Night SST' /
      data      data_lbl(112) / 'MET-A GAC RlxDay SST' /
      data      data_lbl(113) / '   MET-A LAC Day SST' /
      data      data_lbl(114) / ' MET-A LAC Night SST' /
      data      data_lbl(115) / '   MET-B GAC Day SST' /
      data      data_lbl(116) / ' MET-B GAC Night SST' /
      data      data_lbl(117) / 'MET-B GAC RlxDay SST' /
      data      data_lbl(118) / '   MET-B LAC Day SST' /
      data      data_lbl(119) / ' MET-B LAC Night SST' /
      data      data_lbl(120) / '   MET-C GAC Day SST' /
      data      data_lbl(121) / ' MET-C GAC Night SST' /
      data      data_lbl(122) / 'MET-C GAC RlxDay SST' /
      data      data_lbl(123) / '   MET-C LAC Day SST' /
      data      data_lbl(124) / ' MET-C LAC Night SST' /
      data      data_lbl(125) / '         Jason-2 SWH' /
      data      data_lbl(126) / '            Not Used' /
      data      data_lbl(127) / ' Jason-1 Intrlvd SSH' /
      data      data_lbl(128) / '     N19 GAC Day SST' /
      data      data_lbl(129) / '   N19 GAC Night SST' /
      data      data_lbl(130) / '  N19 GAC RlxDay SST' /
      data      data_lbl(131) / '     N19 LAC Day SST' /
      data      data_lbl(132) / '   N19 LAC Night SST' /
      data      data_lbl(133) / '   Animal Borne Temp' /
      data      data_lbl(134) / '   Animal Borne Salt' /
      data      data_lbl(135) / '      GOES13 Day SST' /
      data      data_lbl(136) / '    GOES13 Night SST' /
      data      data_lbl(137) / '          ADCP V Vel' /
      data      data_lbl(138) / '       AMSR2 Sea Ice' /
      data      data_lbl(139) / '     MTSAT-2 Day SST' /
      data      data_lbl(140) / '   MTSAT-2 Night SST' /
      data      data_lbl(141) / '  Expanded Real Temp' /
      data      data_lbl(142) / '  Expanded Real Salt' /
      data      data_lbl(143) / '   Expanded Syn Temp' /
      data      data_lbl(144) / '   Expanded Syn Salt' /
      data      data_lbl(145) / ' Envisat Low-Orb SSH' /
      data      data_lbl(146) / '        Physical SST' /
      data      data_lbl(147) / ' MODIS Aqua Ice Temp' /
      data      data_lbl(148) / 'MODIS Terra Ice Temp' /
      data      data_lbl(149) / '            Ice Temp' /
      data      data_lbl(150) / '    METOP-A Ice Temp' /
      data      data_lbl(151) / '     WindSat Day SST' /
      data      data_lbl(152) / '   WindSat Night SST' /
      data      data_lbl(153) / '   NPP VIIRS Day SST' /
      data      data_lbl(154) / ' NPP VIIRS Night SST' /
      data      data_lbl(155) / 'NPP VIIRS RlxDay SST' /
      data      data_lbl(156) / '  NPP VIIRS Ice Temp' /
      data      data_lbl(157) / '  Argo Traject U Vel' /
      data      data_lbl(158) / '  Argo Traject V Vel' /
      data      data_lbl(159) / '      GOES15 Day SST' /
      data      data_lbl(160) / '    GOES15 Night SST' /
      data      data_lbl(161) / '           SSH U Vel' /
      data      data_lbl(162) / '           SSH V Vel' /
      data      data_lbl(163) / '    Wave Glider Temp' /
      data      data_lbl(164) / '    Wave Glider Salt' /
      data      data_lbl(165) / '     WindSat Sea Ice' /
      data      data_lbl(166) / '       CryoSat-2 SWH' /
      data      data_lbl(167) / '       CryoSat-2 SSH' /
      data      data_lbl(168) / ' Jason1 Geodetic SSH' /
      data      data_lbl(169) / '       AMSR2 Day SST' /
      data      data_lbl(170) / '       Aquarius Salt' /
      data      data_lbl(171) / '    ISOP Temperature' /
      data      data_lbl(172) / '       ISOP Salinity' /
      data      data_lbl(173) / '                 SSS' /
      data      data_lbl(174) / '           SMOS Salt' /
      data      data_lbl(175) / '      GOES14 Day SST' /
      data      data_lbl(176) / '    GOES14 Night SST' /
      data      data_lbl(177) / ' Expanded U Velocity' /
      data      data_lbl(178) / ' Expanded V Velocity' /
      data      data_lbl(179) / '          Altika SSH' /
      data      data_lbl(180) / '          Altika SWH' /
      data      data_lbl(181) / '   HIMAWARI8 Day SST' /
      data      data_lbl(182) / ' HIMAWARI8 Night SST' /
      data      data_lbl(183) / '  Expanded ISOP Temp' /
      data      data_lbl(184) / '  Expanded ISOP Salt' /
      data      data_lbl(185) / '       SSMIS F19 Ice' /
      data      data_lbl(186) / ' eXpendable CTD Temp' /
      data      data_lbl(187) / ' eXpendable CTD Salt' /
      data      data_lbl(188) / '    ALAMO Float Temp' /
      data      data_lbl(189) / '    ALAMO Float Salt' /
      data      data_lbl(190) / '       ASTER Day SST' /
      data      data_lbl(191) / '     ASTER Night SST' /
      data      data_lbl(192) / '       VIIRS NPP Ice' /
      data      data_lbl(193) / ' Jason-2 Intrlvd SSH' /
      data      data_lbl(194) / ' Jason-2 Intrlvd SWH' /
      data      data_lbl(195) / '         Jason-3 SSH' /
      data      data_lbl(196) / '         Jason-3 SWH' /
      data      data_lbl(197) / '     Sentinel-3A SSH' /
      data      data_lbl(198) / '     Sentinel-3A SWH' /
      data      data_lbl(199) / '     Sentinel-3B SSH' /
      data      data_lbl(200) / '     Sentinel-3B SWH' /
      data      data_lbl(201) / '       Altika DP SSH' /
      data      data_lbl(202) / '       Altika DP SWH' /
