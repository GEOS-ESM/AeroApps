esma_set_this ()

set (srcs 
    SOLAR.F90)


esma_add_library (${this}
  SRCS ${srcs}
  )

ecbuild_add_executable (
  TARGET test_solar.x
  SOURCES test_solar.F90
  LIBS ${this}
  )
