esma_set_this()

file (GLOB srcs legacy/lib/*.[Ffc] legacy/lib/*.[fF]90)
list (APPEND srcs
  legacy/trace/calvar.f legacy/select/special.f
  legacy/goodies/evalmod.f90 legacy/goodies/precmod.f90 legacy/goodies/stringmod.f90
  )


esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES chem_base MAPL mpeu gfio_r4
  PUBLIC_LIBS ${MKL_LIBRARIES}
  PUBLIC_INCLUDES ${MKL_INCLUDE_DIRS}
  )
