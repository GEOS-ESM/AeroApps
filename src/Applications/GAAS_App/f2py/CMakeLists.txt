#
# Cmake rules for creating python modules with f2py and 
#

# The OVERRIDE tells esma_set_this() to name the library as
# what we set here rather than the name of the directory

esma_set_this (OVERRIDE GAAS_App.f2py)

# cmake requirements
# ------------------
find_package(F2PY3 REQUIRED)

# Library for the Fortran module
# ------------------------------
esma_add_library(${this}
    SRCS            BuddyCheck_Mod.F90
    DEPENDENCIES    GMAO_mpeu GMAO_psas
    INCLUDES        ${CMAKE_CURRENT_BINARY_DIR}
                    ${CMAKE_BINARY_DIR}/lib
                    ${esma_include}/GMAO_mpeu
                    ${esma_include}/GMAO_psas

)
target_link_libraries(${this} PRIVATE OpenMP::OpenMP_Fortran)

# f2py for buddycheck
# --------------------
esma_add_f2py3_module(BuddyCheck_
    SOURCES        BuddyCheck_py.F90
    DESTINATION    lib/Python/
    LIBRARIES      ${this} GMAO_mpeu GMAO_psas
    INCLUDEDIRS    ${CMAKE_CURRENT_BINARY_DIR}
                   ${CMAKE_BINARY_DIR}/lib
                   ${esma_include}/GMAO_mpeu
                   ${esma_include}/GMAO_psas
                   ${include_${this}}
    USE_OPENMP
    USE_MPI
)

#add_dependencies(BuddyCheck_ ${this} GMAO_mpeu GMAO_psas)

