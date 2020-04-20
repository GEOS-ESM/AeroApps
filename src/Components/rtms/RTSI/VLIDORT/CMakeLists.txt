esma_set_this ()

set (src_dirs sourcecode sourcecode_str LAPACK_directory)
set (srcs)

foreach (dir ${src_dirs})
  file (GLOB srcs_ ${dir}/*.[Ffc] ${dir}/*.[fF]90)
  list (APPEND srcs ${srcs_})
endforeach()

esma_add_library (${this}
  SRCS ${srcs}
  )
