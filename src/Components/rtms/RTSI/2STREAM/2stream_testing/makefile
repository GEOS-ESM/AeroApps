LIDORT_COMPILE = ifort -c -O3 -zero
LAPACK_COMPILE = ifort -c -O3 -zero
TWOSTREAM_COMPILE = ifort -c -O3 -zero
F77 = ifort

# Link definition
#################

LINK.f = $(F77)

# Path variables
################

ROOT_PATH = ..

SPATH_S = $(ROOT_PATH)/sourcecode/
IPATH_S = $(ROOT_PATH)/includes/
2PATH_S = $(ROOT_PATH)/2stream_sourcecode/
OBJ     = $(ROOT_PATH)/OBJECTS

#  OBJECT MODULES
#################

# LIDORT modules in directory sourcecode

#  Masters set

OBJECTS_LIDORT_MASTERS   = $(OBJ)/lidort_master.o      \
                           $(OBJ)/lidort_fourier_master.o

OBJECTS_LIDORT_L_MASTERS = $(OBJ)/lidort_l_master.o    \
                           $(OBJ)/lidort_l_fourier_master.o

#  Basic set for Radiances

OBJECTS_LIDORT_SOLUTIONS = $(OBJ)/lidort_solutions.o   \
                           $(OBJ)/lidort_multipliers.o \
                           $(OBJ)/lidort_bvproblem.o

OBJECTS_LIDORT_INTENSITY = $(OBJ)/lidort_intensity.o   \
                           $(OBJ)/lidort_corrections.o \
                           $(OBJ)/lidort_converge.o

OBJECTS_LIDORT_THERMAL   = $(OBJ)/lidort_thermalsup.o

OBJECTS_LIDORT_SETUPS    = $(OBJ)/lidort_miscsetups.o \
                           $(OBJ)/lidort_inputs.o

OBJECTS_LIDORT_BRDF      = $(OBJ)/lidort_brdf.o       \
                           $(OBJ)/lidort_brdf_kernels.o

OBJECTS_LIDORT_AUX       = $(OBJ)/lidort_aux.o        \
                           $(OBJ)/lidort_writemodules.o

# linearized set

OBJECTS_LIDORT_L_SOLUTIONS = $(OBJ)/lidort_l_solutions.o   \
                             $(OBJ)/lidort_l_multipliers.o \
                             $(OBJ)/lidort_l_bvproblem.o

OBJECTS_LIDORT_L_JACOBIANS = $(OBJ)/lidort_l_wfatmos.o     \
                             $(OBJ)/lidort_l_wfsurface.o   \
                             $(OBJ)/lidort_l_corrections.o \
                             $(OBJ)/lidort_l_converge.o

OBJECTS_LIDORT_L_SETUPS    = $(OBJ)/lidort_l_miscsetups.o  \
                             $(OBJ)/lidort_l_inputs.o

OBJECTS_LIDORT_L_THERMAL   = $(OBJ)/lidort_l_thermalsup.o

OBJECTS_LIDORT_L_BRDF      = $(OBJ)/lidort_l_brdf.o        \
                             $(OBJ)/lidort_l_brdf_kernels.o

OBJECTS_LIDORT_L_AUX       = $(OBJ)/lidort_l_writemodules.o

# 2stream modules in directory 2stream

OBJECTS_TWOSTREAM_I  = $(OBJ)/2stream_masters_basic.o

OBJECTS_TWOSTREAM_M  = $(OBJ)/2stream_solutions.o   \
                     $(OBJ)/2stream_intensity.o     \
                     $(OBJ)/2stream_inputs.o        \
                     $(OBJ)/2stream_bvproblem.o     \
                     $(OBJ)/2stream_miscsetups.o

OBJECTS_TWOSTREAM_MA  = $(OBJ)/2stream_solutions.o   \
                     $(OBJ)/2stream_intensity.o     \
                     $(OBJ)/2stream_inputs.o        \
                     $(OBJ)/2stream_bvproblem.o     \
                     $(OBJ)/2stream_aux_lapack.o    \
                     $(OBJ)/2stream_miscsetups.o

OBJECTS_TWOSTREAM_J  = $(OBJ)/2stream_l_masters_basic.o

OBJECTS_TWOSTREAM_L  = $(OBJ)/2stream_l_solutions.o   \
                     $(OBJ)/2stream_l_jacobians.o     \
                     $(OBJ)/2stream_l_bvproblem.o     \
                     $(OBJ)/2stream_l_miscsetups.o

# environment & interface modules external to the main package

OBJECTS_E1  = $(OBJ)/test_lidort_2stream.o

OBJECTS_E3  = $(OBJ)/test_2sj.o
OBJECTS_E3A  = $(OBJ)/test_2sjl.o

#OBJECTS_00  = $(OBJ)/read_prop.o

# EXECUTABLE
############

APPLICATION_1           = test_lidort_2stream.exe
APPLICATION_3           = test_2sj.exe
APPLICATION_3A          = test_2sjl.exe

#all: $(APPLICATION_1)
#all: $(APPLICATION_3) $(APPLICATION_3A)
all: $(APPLICATION_1) $(APPLICATION_3) $(APPLICATION_3A)

$(APPLICATION_1):   $(OBJECTS_E1)               \
                  $(OBJECTS_TWOSTREAM_I)        \
                  $(OBJECTS_TWOSTREAM_M)        \
                  $(OBJECTS_LIDORT_MASTERS)     \
                  $(OBJECTS_LIDORT_SOLUTIONS)   \
                  $(OBJECTS_LIDORT_INTENSITY)   \
                  $(OBJECTS_LIDORT_THERMAL)     \
                  $(OBJECTS_LIDORT_SETUPS)      \
                  $(OBJECTS_LIDORT_BRDF)        \
                  $(OBJECTS_LIDORT_AUX)
	$(LINK.f) -o $(APPLICATION_1)             \
                  $(OBJECTS_E1)                 \
                  $(OBJECTS_TWOSTREAM_I)        \
                  $(OBJECTS_TWOSTREAM_M)        \
                  $(OBJECTS_LIDORT_MASTERS)     \
                  $(OBJECTS_LIDORT_SOLUTIONS)   \
                  $(OBJECTS_LIDORT_INTENSITY)   \
                  $(OBJECTS_LIDORT_THERMAL)     \
                  $(OBJECTS_LIDORT_SETUPS)      \
                  $(OBJECTS_LIDORT_BRDF)        \
                  $(OBJECTS_LIDORT_AUX)         \
                  $(FLIBRARIES)

$(APPLICATION_3):   $(OBJECTS_E3)               \
                  $(OBJECTS_TWOSTREAM_I)        \
                  $(OBJECTS_TWOSTREAM_J)        \
                  $(OBJECTS_TWOSTREAM_MA)       \
                  $(OBJECTS_TWOSTREAM_L)
	$(LINK.f) -o $(APPLICATION_3)             \
                  $(OBJECTS_E3)                 \
                  $(OBJECTS_TWOSTREAM_I)        \
                  $(OBJECTS_TWOSTREAM_J)        \
                  $(OBJECTS_TWOSTREAM_MA)       \
                  $(OBJECTS_TWOSTREAM_L)        \
                  $(FLIBRARIES)

$(APPLICATION_3A):   $(OBJECTS_E3A)               \
                  $(OBJECTS_LIDORT_MASTERS)     \
                  $(OBJECTS_LIDORT_SOLUTIONS)   \
                  $(OBJECTS_LIDORT_INTENSITY)   \
                  $(OBJECTS_LIDORT_THERMAL)     \
                  $(OBJECTS_LIDORT_SETUPS)      \
                  $(OBJECTS_LIDORT_BRDF)        \
                  $(OBJECTS_LIDORT_AUX)         \
                  $(OBJECTS_LIDORT_L_MASTERS)   \
                  $(OBJECTS_LIDORT_L_SOLUTIONS) \
                  $(OBJECTS_LIDORT_L_JACOBIANS) \
                  $(OBJECTS_LIDORT_L_THERMAL)   \
                  $(OBJECTS_LIDORT_L_SETUPS)    \
                  $(OBJECTS_LIDORT_L_BRDF)      \
                  $(OBJECTS_LIDORT_L_AUX)
	$(LINK.f) -o $(APPLICATION_3A)             \
                  $(OBJECTS_E3A)                 \
                  $(OBJECTS_LIDORT_MASTERS)     \
                  $(OBJECTS_LIDORT_SOLUTIONS)   \
                  $(OBJECTS_LIDORT_INTENSITY)   \
                  $(OBJECTS_LIDORT_THERMAL)     \
                  $(OBJECTS_LIDORT_SETUPS)      \
                  $(OBJECTS_LIDORT_BRDF)        \
                  $(OBJECTS_LIDORT_AUX)         \
                  $(OBJECTS_LIDORT_L_MASTERS)   \
                  $(OBJECTS_LIDORT_L_SOLUTIONS) \
                  $(OBJECTS_LIDORT_L_JACOBIANS) \
                  $(OBJECTS_LIDORT_L_THERMAL)   \
                  $(OBJECTS_LIDORT_L_SETUPS)    \
                  $(OBJECTS_LIDORT_L_BRDF)      \
                  $(OBJECTS_LIDORT_L_AUX)       \
                  $(FLIBRARIES)

# COMPILATIONS
##############

distclean:
	rm -f *.exe $(OBJ)/*.o

#------------------------------------------- Environment modules

$(OBJ)/test_lidort_2stream.o: test_lidort_2stream.f \
	 $(IPATH_S)LIDORT.PARS               \
	 $(IPATH_S)LIDORT_INPUTS.VARS        \
	 $(IPATH_S)LIDORT_BOOKKEEP.VARS      \
	 $(IPATH_S)LIDORT_RESULTS.VARS
	$(LIDORT_COMPILE) test_lidort_2stream.f -o $(OBJ)/test_lidort_2stream.o

$(OBJ)/read_prop.o: read_prop.f
	$(LIDORT_COMPILE) read_prop.f -o $(OBJ)/read_prop.o

$(OBJ)/test_2sj.o: test_2sj.f
	$(LIDORT_COMPILE) test_2sj.f -o $(OBJ)/test_2sj.o

$(OBJ)/test_2sjl.o: test_2sjl.f \
	 $(IPATH_S)LIDORT.PARS               \
	 $(IPATH_S)LIDORT_INPUTS.VARS        \
	 $(IPATH_S)LIDORT_L_INPUTS.VARS      \
	 $(IPATH_S)LIDORT_BOOKKEEP.VARS      \
	 $(IPATH_S)LIDORT_RESULTS.VARS       \
	 $(IPATH_S)LIDORT_L_RESULTS.VARS
	$(LIDORT_COMPILE) test_2sjl.f -o $(OBJ)/test_2sjl.o

# ---------------------------------------------------
#----------------------LIDORT L_master modules ------
# ---------------------------------------------------

$(OBJ)/lidort_l_master.o: $(SPATH_S)lidort_l_master.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_L_INPUTS.VARS       \
         $(IPATH_S)LIDORT_RESULTS.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_l_master.f -o $(OBJ)/lidort_l_master.o

$(OBJ)/lidort_l_fourier_master.o: $(SPATH_S)lidort_l_fourier_master.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_L_INPUTS.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_THERMALSUP.VARS     \
         $(IPATH_S)LIDORT_SOLUTION.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_l_fourier_master.f -o $(OBJ)/lidort_l_fourier_master.o

# ---------------------------------------------------
#----------------------LIDORT master modules --------
# ---------------------------------------------------

$(OBJ)/lidort_master.o: $(SPATH_S)lidort_master.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_RESULTS.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_master.f -o $(OBJ)/lidort_master.o

$(OBJ)/lidort_fourier_master.o: $(SPATH_S)lidort_fourier_master.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_THERMALSUP.VARS     \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_S)LIDORT_RESULTS.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_fourier_master.f -o $(OBJ)/lidort_fourier_master.o

# ---------------------------------------------------
#----------------------LIDORT Radiance modules ------
# ---------------------------------------------------

#--------------------------- LIDORT solutions modules

$(OBJ)/lidort_solutions.o: $(SPATH_S)lidort_solutions.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_SOLUTION.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_solutions.f -o $(OBJ)/lidort_solutions.o

$(OBJ)/lidort_multipliers.o: $(SPATH_S)lidort_multipliers.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_SOLUTION.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_multipliers.f -o $(OBJ)/lidort_multipliers.o

$(OBJ)/lidort_bvproblem.o: $(SPATH_S)lidort_bvproblem.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_REFLECTANCE.VARS    \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_SOLUTION.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_bvproblem.f -o $(OBJ)/lidort_bvproblem.o

#------------------------------ LIDORT Intensity modules

$(OBJ)/lidort_intensity.o: $(SPATH_S)lidort_intensity.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_THERMALSUP.VARS     \
         $(IPATH_S)LIDORT_REFLECTANCE.VARS    \
         $(IPATH_S)LIDORT_RESULTS.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_intensity.f -o $(OBJ)/lidort_intensity.o

$(OBJ)/lidort_converge.o: $(SPATH_S)lidort_converge.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_SINGSCAT.VARS       \
         $(IPATH_S)LIDORT_REFLECTANCE.VARS    \
         $(IPATH_S)LIDORT_RESULTS.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_converge.f -o $(OBJ)/lidort_converge.o

$(OBJ)/lidort_corrections.o: $(SPATH_S)lidort_corrections.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_RESULTS.VARS        \
         $(IPATH_S)LIDORT_REFLECTANCE.VARS    \
         $(IPATH_S)LIDORT_SINGSCAT.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_corrections.f -o $(OBJ)/lidort_corrections.o

#-------------------------- LIDORT Setups modules

$(OBJ)/lidort_miscsetups.o: $(SPATH_S)lidort_miscsetups.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_REFLECTANCE.VARS    \
         $(IPATH_S)LIDORT_SETUPS.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_miscsetups.f -o $(OBJ)/lidort_miscsetups.o

$(OBJ)/lidort_inputs.o: $(SPATH_S)lidort_inputs.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_inputs.f -o $(OBJ)/lidort_inputs.o

#---------------------------- LIDORT Thermal supplement

$(OBJ)/lidort_thermalsup.o: $(SPATH_S)lidort_thermalsup.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_THERMALSUP.VARS     \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_thermalsup.f -o $(OBJ)/lidort_thermalsup.o

#---------------------------- LIDORT BRDF modules

$(OBJ)/lidort_brdf.o: $(SPATH_S)lidort_brdf.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_REFLECTANCE.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_brdf.f -o $(OBJ)/lidort_brdf.o

$(OBJ)/lidort_brdf_kernels.o: $(SPATH_S)lidort_brdf_kernels.f    \
         $(IPATH_S)LIDORT.PARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_brdf_kernels.f -o $(OBJ)/lidort_brdf_kernels.o

#------------------------------- LIDORT Auxiliary modules

$(OBJ)/lidort_aux.o: $(SPATH_S)lidort_aux.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS
	$(LAPACK_COMPILE) $(SPATH_S)lidort_aux.f -o $(OBJ)/lidort_aux.o

$(OBJ)/lidort_writemodules.o: $(SPATH_S)lidort_writemodules.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_RESULTS.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_writemodules.f -o $(OBJ)/lidort_writemodules.o

# ---------------------------------------------------
#----------------------LIDORT Jacobian modules ------
# ---------------------------------------------------

#---------------LIDORT linearized solutions modules

$(OBJ)/lidort_l_solutions.o: $(SPATH_S)lidort_l_solutions.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_S)LIDORT_L_INPUTS.VARS       \
         $(IPATH_S)LIDORT_L_SETUPS.VARS       \
         $(IPATH_S)LIDORT_L_SOLUTION.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_l_solutions.f -o $(OBJ)/lidort_l_solutions.o

$(OBJ)/lidort_l_multipliers.o: $(SPATH_S)lidort_l_multipliers.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_S)LIDORT_L_INPUTS.VARS       \
         $(IPATH_S)LIDORT_L_SETUPS.VARS       \
         $(IPATH_S)LIDORT_L_SOLUTION.VARS     \
         $(IPATH_S)LIDORT_L_MULTIPLIERS.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_l_multipliers.f -o $(OBJ)/lidort_l_multipliers.o

$(OBJ)/lidort_l_bvproblem.o: $(SPATH_S)lidort_l_bvproblem.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_REFLECTANCE.VARS    \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_S)LIDORT_L_INPUTS.VARS       \
         $(IPATH_S)LIDORT_L_SETUPS.VARS       \
         $(IPATH_S)LIDORT_L_SOLUTION.VARS     \
         $(IPATH_S)LIDORT_L_THERMALSUP.VARS   \
         $(IPATH_S)LIDORT_L_MULTIPLIERS.VARS  \
         $(IPATH_S)LIDORT_L_REFLECTANCE.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_l_bvproblem.f -o $(OBJ)/lidort_l_bvproblem.o

#---------------------------- LIDORT Linearized Thermal supplement

$(OBJ)/lidort_l_thermalsup.o: $(SPATH_S)lidort_l_thermalsup.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_THERMALSUP.VARS     \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_L_INPUTS.VARS       \
         $(IPATH_S)LIDORT_L_SETUPS.VARS       \
         $(IPATH_S)LIDORT_L_SOLUTION.VARS     \
         $(IPATH_S)LIDORT_L_THERMALSUP.VARS   \
         $(IPATH_S)LIDORT_L_MULTIPLIERS.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_l_thermalsup.f -o $(OBJ)/lidort_l_thermalsup.o

#------------------------- LIDORT weighting function modules

$(OBJ)/lidort_l_wfatmos.o: $(SPATH_S)lidort_l_wfatmos.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_S)LIDORT_THERMALSUP.VARS     \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_REFLECTANCE.VARS    \
         $(IPATH_S)LIDORT_RESULTS.VARS        \
         $(IPATH_S)LIDORT_L_SETUPS.VARS       \
         $(IPATH_S)LIDORT_L_SOLUTION.VARS     \
         $(IPATH_S)LIDORT_L_MULTIPLIERS.VARS  \
         $(IPATH_S)LIDORT_L_THERMALSUP.VARS   \
         $(IPATH_S)LIDORT_L_REFLECTANCE.VARS  \
         $(IPATH_S)LIDORT_L_RESULTS.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_l_wfatmos.f -o $(OBJ)/lidort_l_wfatmos.o

$(OBJ)/lidort_l_wfsurface.o: $(SPATH_S)lidort_l_wfsurface.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_REFLECTANCE.VARS    \
         $(IPATH_S)LIDORT_RESULTS.VARS        \
         $(IPATH_S)LIDORT_L_SOLUTION.VARS     \
         $(IPATH_S)LIDORT_L_REFLECTANCE.VARS  \
         $(IPATH_S)LIDORT_L_RESULTS.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_l_wfsurface.f -o $(OBJ)/lidort_l_wfsurface.o

$(OBJ)/lidort_l_converge.o: $(SPATH_S)lidort_l_converge.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SINGSCAT.VARS       \
         $(IPATH_S)LIDORT_REFLECTANCE.VARS    \
         $(IPATH_S)LIDORT_RESULTS.VARS        \
         $(IPATH_S)LIDORT_L_INPUTS.VARS       \
         $(IPATH_S)LIDORT_L_SINGSCAT.VARS     \
         $(IPATH_S)LIDORT_L_REFLECTANCE.VARS  \
         $(IPATH_S)LIDORT_L_RESULTS.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_l_converge.f -o $(OBJ)/lidort_l_converge.o

$(OBJ)/lidort_l_corrections.o: $(SPATH_S)lidort_l_corrections.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_SINGSCAT.VARS       \
         $(IPATH_S)LIDORT_RESULTS.VARS        \
         $(IPATH_S)LIDORT_L_INPUTS.VARS       \
         $(IPATH_S)LIDORT_L_SINGSCAT.VARS     \
         $(IPATH_S)LIDORT_L_SETUPS.VARS       \
         $(IPATH_S)LIDORT_L_SOLUTION.VARS     \
         $(IPATH_S)LIDORT_L_MULTIPLIERS.VARS  \
         $(IPATH_S)LIDORT_L_REFLECTANCE.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_l_corrections.f -o $(OBJ)/lidort_l_corrections.o

#--------------------- LIDORT Linearized Setups modules

$(OBJ)/lidort_l_miscsetups.o: $(SPATH_S)lidort_l_miscsetups.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_L_INPUTS.VARS       \
         $(IPATH_S)LIDORT_L_SETUPS.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_l_miscsetups.f -o $(OBJ)/lidort_l_miscsetups.o

$(OBJ)/lidort_l_inputs.o: $(SPATH_S)lidort_l_inputs.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_L_INPUTS.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_l_inputs.f -o $(OBJ)/lidort_l_inputs.o

#---------------------- LIDORT Linearized BRDF modules

$(OBJ)/lidort_l_brdf.o: $(SPATH_S)lidort_l_brdf.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_REFLECTANCE.VARS    \
         $(IPATH_S)LIDORT_L_INPUTS.VARS       \
         $(IPATH_S)LIDORT_L_REFLECTANCE.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_l_brdf.f -o $(OBJ)/lidort_l_brdf.o

$(OBJ)/lidort_l_brdf_kernels.o: $(SPATH_S)lidort_l_brdf_kernels.f    \
         $(IPATH_S)LIDORT.PARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_l_brdf_kernels.f -o $(OBJ)/lidort_l_brdf_kernels.o

#----------------------- LIDORT Linearized Auxiliary modules

$(OBJ)/lidort_l_writemodules.o: $(SPATH_S)lidort_l_writemodules.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_INPUTS.VARS         \
         $(IPATH_S)LIDORT_BOOKKEEP.VARS       \
         $(IPATH_S)LIDORT_RESULTS.VARS        \
         $(IPATH_S)LIDORT_L_RESULTS.VARS
	$(LIDORT_COMPILE) $(SPATH_S)lidort_l_writemodules.f -o $(OBJ)/lidort_l_writemodules.o


#-------------------------------------
#   TWO_STREAM MODULES
#-------------------------------------


#_I

$(OBJ)/2stream_masters_basic.o: $(2PATH_S)2stream_masters_basic.f
	$(TWOSTREAM_COMPILE) $(2PATH_S)2stream_masters_basic.f -o $(OBJ)/2stream_masters_basic.o

#_J

$(OBJ)/2stream_l_masters_basic.o: $(2PATH_S)2stream_l_masters_basic.f
	$(TWOSTREAM_COMPILE) $(2PATH_S)2stream_l_masters_basic.f -o $(OBJ)/2stream_l_masters_basic.o

#_M

$(OBJ)/2stream_bvproblem.o: $(2PATH_S)2stream_bvproblem.f
	$(TWOSTREAM_COMPILE) $(2PATH_S)2stream_bvproblem.f -o $(OBJ)/2stream_bvproblem.o

$(OBJ)/2stream_miscsetups.o: $(2PATH_S)2stream_miscsetups.f
	$(TWOSTREAM_COMPILE) $(2PATH_S)2stream_miscsetups.f -o $(OBJ)/2stream_miscsetups.o

$(OBJ)/2stream_solutions.o: $(2PATH_S)2stream_solutions.f
	$(TWOSTREAM_COMPILE) $(2PATH_S)2stream_solutions.f -o $(OBJ)/2stream_solutions.o

$(OBJ)/2stream_aux_lapack.o: $(2PATH_S)2stream_aux_lapack.f
	$(TWOSTREAM_COMPILE) $(2PATH_S)2stream_aux_lapack.f -o $(OBJ)/2stream_aux_lapack.o

$(OBJ)/2stream_intensity.o: $(2PATH_S)2stream_intensity.f
	$(TWOSTREAM_COMPILE) $(2PATH_S)2stream_intensity.f -o $(OBJ)/2stream_intensity.o

$(OBJ)/2stream_inputs.o: $(2PATH_S)2stream_inputs.f
	$(TWOSTREAM_COMPILE) $(2PATH_S)2stream_inputs.f -o $(OBJ)/2stream_inputs.o

#_L

$(OBJ)/2stream_l_bvproblem.o: $(2PATH_S)2stream_l_bvproblem.f
	$(TWOSTREAM_COMPILE) $(2PATH_S)2stream_l_bvproblem.f -o $(OBJ)/2stream_l_bvproblem.o

$(OBJ)/2stream_l_miscsetups.o: $(2PATH_S)2stream_l_miscsetups.f
	$(TWOSTREAM_COMPILE) $(2PATH_S)2stream_l_miscsetups.f -o $(OBJ)/2stream_l_miscsetups.o

$(OBJ)/2stream_l_solutions.o: $(2PATH_S)2stream_l_solutions.f
	$(TWOSTREAM_COMPILE) $(2PATH_S)2stream_l_solutions.f -o $(OBJ)/2stream_l_solutions.o

$(OBJ)/2stream_l_jacobians.o: $(2PATH_S)2stream_l_jacobians.f
	$(TWOSTREAM_COMPILE) $(2PATH_S)2stream_l_jacobians.f -o $(OBJ)/2stream_l_jacobians.o

