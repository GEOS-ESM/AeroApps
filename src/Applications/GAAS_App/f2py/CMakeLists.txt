#
# Cmake rules for creating python modules with f2py and 
#

# The OVERRIDE tells esma_set_this() to name the library as
# what we set here rather than the name of the directory

#esma_set_this ()

# cmake requirements
# ------------------
find_package(F2PY3 REQUIRED)

# Library for the Fortran module
# ------------------------------
esma_add_library(BuddyCheck_Mod
    SRCS            BuddyCheck_Mod.F90
    DEPENDENCIES    GMAO_mpeu GMAO_psas
    INCLUDES        ${CMAKE_CURRENT_BINARY_DIR}
                    ${CMAKE_BINARY_DIR}/lib
                    ${esma_include}/GMAO_mpeu
                    ${esma_include}/GMAO_psas
    USE_OPENMP
)
# f2py for buddycheck
# --------------------
esma_add_f2py3_module(BuddyCheck_
    SOURCES        BuddyCheck_py.F90
    DESTINATION    lib/Python
    LIBRARIES      BuddyCheck_Mod GMAO_mpeu GMAO_psas
    INCLUDEDIRS    ${CMAKE_CURRENT_BINARY_DIR}
                   ${CMAKE_BINARY_DIR}/lib
                   ${esma_include}/GMAO_mpeu
                   ${esma_include}/GMAO_psas
		   ${include_buddycheck_mod}
    USE_OPENMP
    USE_MPI
)

add_dependencies(BuddyCheck_ BuddyCheck_Mod GMAO_mpeu GMAO_psas)

