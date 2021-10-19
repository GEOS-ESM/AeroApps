if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${VLIDORT_LIBRARY}" AND IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${VLIDORT_LIBRARY}")
  message(STATUS "Building with VLIDORT library ${VLIDORT_LIBRARY}")
else()
  message(FATAL_ERROR "VLIDORT library ${VLIDORT_LIBRARY} does not exist")
endif()

# NOTE: VLIDORT is an old version that is needed for the OMI_ target (needs VLIDORT_LambMod)
#       and VLIDORT90 is an old version required by VLIDORT_cloud_ target (needs VLIDORT_ScatMod_cd)
#       While these are "VLIDORT Library" it is not what we consider for this variable
esma_add_subdirectories (2STREAM LIDORT3p7 VLIDORT VLIDORT90 ${VLIDORT_LIBRARY})
