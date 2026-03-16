# We only want to build and use one VLIDORT library

set(VLIDORT_DEFAULT 2p8p3)
set(VLIDORT ${VLIDORT_DEFAULT} CACHE STRING "VLIDORT Library to link to")
string(CONCAT VLIDORT_LIBRARY "VLIDORT" ${VLIDORT})

add_subdirectory (cda)
add_subdirectory (misc)
add_subdirectory (missions)
#add_subdirectory (qfed)
add_subdirectory (rtms)
