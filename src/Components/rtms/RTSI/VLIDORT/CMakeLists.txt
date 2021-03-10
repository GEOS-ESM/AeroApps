esma_set_this ()

set (src_dirs sourcecode sourcecode_str LAPACK_directory)
set (srcs VLIDORT_LambMod.F)
set_property(SOURCE VLIDORT_LambMod.F PROPERTY COMPILE_FLAGS ${EXTENDED_SOURCE})

foreach (dir ${src_dirs})
  file (GLOB srcs_ ${dir}/*.[Ffc] ${dir}/*.[fF]90)
  list (APPEND srcs ${srcs_})
endforeach()

esma_add_library (${this}
  SRCS ${srcs}
  PUBLIC_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/includes
  )
