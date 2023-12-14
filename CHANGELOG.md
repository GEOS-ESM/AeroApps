# ChangeLog

## [Unreleased]

### Added

- ASIA-AQ mission directory initialized with FLUID instance from FIREX-AQ
- pyabc module added in GMAO_aeropybs. the snket class is needed in pyabc for loading ffnet .net files
- VIIRS and MODIS NNR angstrom exponent predictions are now saved in gridded output files (Level 3)
- VIIRS and MODIS NNR ods files of AOD can be generated from angstrom exponent fit coefficients
- VIIRS NNR ods file generator now takes nsyn (number of synoptic times) as an input.
- VIIRS & MODIS NNR codes have an improved cloud filter and analomaly detection algorithm. There parameters for this can be inputs in the pcf file.

### Fixed

- bug in NNR cloud contamination check. Was not setting suspected too cloudy pixels to MISSING.

### Changed

- Update GMAOpyobs to v1.0.4
- Update GMAOpyobs to v1.0.5
- Update GMAOpyobs to v1.0.6

## [v2.0.2] - 2023-05-25

### Added

- VIIRS NNR training code - new VIIRS giant reader
- New code in GAAS_App to generate NNR ODS files from VIIRS obs

### Fixed

- fix in NNR code for new syntax of the Kfold generator in sklearn
- fix in NNR testing code to protect against cases where VIIRS standard product may not retrieve multiple wavelengths 
- fix pyods __ini__.py to use relative imports. needed for python3
- fix for check if AEfitb is None is mxd04_nnr and vx04_nnr. Using == gives error when AEfitb is initialized to an array
### Changed

- add VIIRS aerosol products to GMAO_ods kx list
- update GMAOpyobs to v1.0.3 

## [v2.0.1] - 2023-05-17

### Added


### Fixed

- fix in trj_sampler for ICARTT files that have missing values in the location/time
- fix in py_ods that was trying to do element wise comparisons on tuples and lists

### Changed

- updated GMAOpyobs to v1.0.2
- updates to NNR training code to handle angstrom exponent targets
- some other minor changes to NNR training code to be python3 compliant

## [v2.0.0] - 2023-05-17

### Added

- Added changelog enforcer

### Changed

- Updated to use latest components (matching GEOSgcm v10.25.0)
- Wholesale conversion to python3. Changes to all python scripts, and CMakelists with f2py compilation
- the components.yaml now imports the GMAOpyobs repository, which replaces the legacy GMAO_pyobs directory
- The GMAO_pyobs3 directory is no longer needed
- A new directory called GMAO_aeropyobs has been added to share. This includes python utilities that require MAPL or GFIO that were previously in GMAO_pyobs. 
- There may still be issues related to byte-to-string conversions that require fixing, but will need to be fixed as they are encountered.
- Added GMAOpyobs v1.0.1 to components.yaml
- Updated MAPL to v2.39.9

## [v1.0.0] - 2022-06-17

### Changed

- Converted to mepo, cmake, and v3 giant files

### Added

- initial support for mepo
- initial treatment for cmake
- legacy MAPL_Base/Python/MAPL into latest MAPL on develop branch

### Fixed

- Moved src/GMAO_Shared to src/Shared/GMAO_Shared
- obs_aod codes now work with v3 giant files

## [0.99.0] - 2020-03-23

Initial GitHub release.

### Added
  - Open source license
  - basic changelog (this file)
  
### Changed
  - removed some CVS legacy directories.
  
