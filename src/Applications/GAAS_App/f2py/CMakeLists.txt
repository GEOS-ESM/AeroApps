#
# Cmake rules for creating python modules with f2py and 
#

# The OVERRIDE tells esma_set_this() to name the library as
# what we set here rather than the name of the directory

esma_set_this (OVERRIDE GAAS_App.f2py)

# cmake requirements
# ------------------
find_package(F2PY3 REQUIRED)

set (srcs             
     Geometry_icosahedron.F90
     Icosahedron_Mod.F90
)

# Library for the Fortran module
# ------------------------------
esma_add_library(${this}
    SRCS            ${srcs} 
        INCLUDES    ${CMAKE_CURRENT_BINARY_DIR}
                    ${CMAKE_BINARY_DIR}/lib
                    ${include_${this}}
)
#target_link_libraries(${this} PRIVATE OpenMP::OpenMP_Fortran)
#target_link_libraries(${this})

# f2py for buddycheck
# --------------------
#esma_add_f2py3_module(BuddyCheck_
#    SOURCES        BuddyCheck_py.F90
#    DESTINATION    lib/Python/
#    LIBRARIES      ${this} GMAO_psas GMAO_mpeu
#    INCLUDEDIRS    ${CMAKE_CURRENT_BINARY_DIR}
#                   ${CMAKE_BINARY_DIR}/lib
#                   ${esma_include}/GMAO_psas
#                   ${esma_include}/GMAO_mpeu
#                   ${include_${this}}
#    DEPENDS        ${this}		   
#        USE_OPENMP
#    USE_MPI
#)

#add_dependencies(BuddyCheck_ ${this} GMAO_mpeu GMAO_psas)

# f2py for buddycheck
# --------------------
esma_add_f2py3_module(Icosahedron_
	SOURCES        Icosahedron_py.F90
    DESTINATION    lib/Python/
    LIBRARIES      ${this} 
    INCLUDEDIRS    ${CMAKE_CURRENT_BINARY_DIR}
                   ${CMAKE_BINARY_DIR}/lib
                   ${include_${this}}
    DEPENDS        ${this}		   
    USE_MPI
    USE_OPENMP
    )
add_dependencies(Icosahedron_ ${this})
