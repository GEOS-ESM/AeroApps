esma_set_this()

add_f2py_module (duem_
  SOURCES duem_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

# these are executable scripts
set (PYSCRIPTS
    dee_gaplot-xy.py
    dee_mm.py
    deep.py
    m2_anom.py
    m2_mm.py
    )

install (
    PROGRAMS  ${PYSCRIPTS}
    DESTINATION bin/${this}
    )

install (
   DIRECTORY jobs
   DESTINATION bin/${this}
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )

# these are non-executaables
set (RCFILES
    DU_GridComp.rc
    )

set (OTHERFILES
    ginoux1.gsf
    ginoux5.gsf
    ginoux.udxt
    go.gs
    )

install (
    FILES ${RCFILES} ${OTHERFILES}
    DESTINATION bin/${this}
    )
