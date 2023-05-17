# ChangeLog

## [Unreleased]

### Added

- Added changelog enforcer

### Fixed

### Changed [2023-03]

- Wholesale conversion to python3. Changes to all python scripts, and CMakelists with f2py compilation
- the components.yaml now imports the GMAOpyobs repository, which replaces the legacy GMAO_pyobs directory
- The GMAO_pyobs3 directory is no longer needed
- A new directory called GMAO_aeropyobs has been added to share. This includes python utilities that require MAPL or GFIO that were previously in GMAO_pyobs. 
- There may still be issues related to byte-to-string conversions that require fixing, but will need to be fixed as they are encountered.
- Added GMAOpyobs v1.0.1 to components.yaml
- Updated MAPL to v2.39.9

### Changed

- Updated to use latest components (matching GEOSgcm v10.25.0)

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
  
