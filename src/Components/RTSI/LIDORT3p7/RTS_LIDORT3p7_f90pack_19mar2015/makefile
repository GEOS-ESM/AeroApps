#
# This file serves as the head Makefile for LIDORT and its "tester" programs
#

#
# UPGRADE TO LIDORT VERSION 3.7, May 2014
#

#
# Define some variables
#

UTIL_PATH = util
SUP_PATH = sup
LID_DEF_PATH = lidort_def
LID_MAIN_PATH = lidort_main
LID_TEST_PATH = lidort_test

MOD_PATH = mod
OBJ_PATH = obj

MOD_FILES = $(MOD_PATH)/*.mod
OBJ_FILES = $(OBJ_PATH)/*.o


#
# Define default shell make will use
#

SHELL = /bin/bash


#
# Define FORTRAN90 compiler to use (can be defined here, but usually defined on makefile command line)
#

#  Intel
#FC = ifort

#  gfortran
#FC = gfortran

#  g95
#FC = g95

#  NAG
#FC = f95


#
# Define FORTRAN90 compiler flags
#

FFLAGS = -c

# Additional flags for Intel
ifeq ($(FC), ifort)
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -module $(MOD_PATH)
	FFLAGS_DEBUG = -g -warn all -check all -traceback
	FFLAGS_OPENMP = -openmp
	FFLAGS_OPT = -O3
endif

# Additional flags for gfortran
ifeq ($(FC), gfortran)
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -J$(MOD_PATH)
	FFLAGS_DEBUG  = -g -C -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
	FFLAGS_Catch  =       -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
	FFLAGS_OPT    = -O3
	FFLAGS_OPENMP = -fopenmp -frecursive
endif

# Additional flags for g95
ifeq ($(FC), g95)
#      older g95
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -fmod=$(MOD_PATH)
#	FFLAGS_DEBUG = -g -Wall -fbounds-check
	FFLAGS_DEBUG = -g -Wall -Wuninitialized -fbounds-check
	FFLAGS_OPT = -O3
#      g95 v0.92
#	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -fmod=$(MOD_PATH)
#	FFLAGS_DEBUG = -g -Wall -fbounds-check -ftrace=full
endif

# Additional flags for NAG
#ifeq ($(FC), f95)
ifeq ($(FC), $(filter $(FC), f95 nagfor))
	FFLAGS := $(FFLAGS) -mdir $(MOD_PATH) -I$(MOD_PATH)
	FFLAGS_DEBUG = -w=obs -w=unused -C=array -C=undefined -gline
#	FFLAGS_DEBUG = -w=obs -w=unused -g -C=all -C=undefined -gline -mtrace=all -nan
	FFLAGS_OPT = -O3
        FFLAGS_OPENMP = -openmp
endif

#  for compilation detailed
#	FFLAGS := $(FFLAGS) $(FFLAGS_Catch)

# For debug build, use "make DEBUG=t"
ifeq ($(DEBUG), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_DEBUG)
endif

# For optimized build, use "make OPT=t"
ifeq ($(OPT), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_OPT)
endif

# For parallel build using OpenMP, use "make OPENMP=t"
ifeq ($(OPENMP), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_OPENMP)
endif

.SUFFIXES:

#
# Define list of source files
# (Note: ordering is important because of dependencies)
#

BASE_SOURCES =
SOURCES =
L_SOURCES =
LPS_SOURCES =
LCS_SOURCES =

BASE_SOURCES +=   \
   $(LID_DEF_PATH)/lidort_pars.f90

SOURCES +=   \
   $(BASE_SOURCES) \
   $(LID_DEF_PATH)/lidort_inputs_def.f90        \
   $(LID_DEF_PATH)/lidort_sup_brdf_def.f90      \
   $(LID_DEF_PATH)/lidort_sup_ss_def.f90        \
   $(LID_DEF_PATH)/lidort_sup_sleave_def.f90    \
   $(LID_DEF_PATH)/lidort_sup_def.f90           \
   $(LID_DEF_PATH)/lidort_outputs_def.f90       \
   $(LID_DEF_PATH)/lidort_io_defs.f90           \
   $(LID_MAIN_PATH)/lidort_aux.f90              \
   $(LID_MAIN_PATH)/lidort_Taylor.f90           \
   $(LID_MAIN_PATH)/lidort_getplanck.f90	\
   $(LID_MAIN_PATH)/lidort_geometry.f90         \
   $(LID_MAIN_PATH)/lidort_inputs.f90           \
   $(LID_MAIN_PATH)/lidort_miscsetups.f90       \
   $(LID_MAIN_PATH)/lidort_corrections.f90      \
   $(LID_MAIN_PATH)/lidort_thermalsup.f90       \
   $(LID_MAIN_PATH)/lidort_solutions.f90        \
   $(LID_MAIN_PATH)/lidort_bvproblem.f90        \
   $(LID_MAIN_PATH)/lidort_PostProcessing.f90   \
   $(LID_MAIN_PATH)/lidort_intensity.f90        \
   $(LID_MAIN_PATH)/lidort_writemodules.f90     \
   $(LID_MAIN_PATH)/lidort_masters.f90

L_SOURCES += \
   $(LID_DEF_PATH)/lidort_lin_inputs_def.f90	\
   $(LID_DEF_PATH)/lidort_lin_sup_brdf_def.f90  \
   $(LID_DEF_PATH)/lidort_lin_sup_ss_def.f90	\
   $(LID_DEF_PATH)/lidort_lin_sup_sleave_def.f90\
   $(LID_DEF_PATH)/lidort_lin_sup_def.f90	\
   $(LID_DEF_PATH)/lidort_lin_outputs_def.f90	\
   $(LID_DEF_PATH)/lidort_lin_io_defs.f90	\
   $(LID_MAIN_PATH)/lidort_l_inputs.f90         \
   $(LID_MAIN_PATH)/lidort_la_miscsetups.f90    \
   $(LID_MAIN_PATH)/lidort_lpc_solutions.f90    \
   $(LID_MAIN_PATH)/lidort_ls_corrections.f90   \
   $(LID_MAIN_PATH)/lidort_l_thermalsup.f90     \
   $(LID_MAIN_PATH)/lidort_ls_wfsurface.f90     \
   $(LID_MAIN_PATH)/lidort_ls_wfsleave.f90      \
   $(LID_MAIN_PATH)/lidort_lbbf_jacobians.f90   \
   $(LID_MAIN_PATH)/lidort_l_writemodules.f90

LPS_SOURCES += \
   $(SOURCES) \
   $(L_SOURCES) \
   $(LID_MAIN_PATH)/lidort_lp_miscsetups.f90      \
   $(LID_MAIN_PATH)/lidort_lp_corrections.f90     \
   $(LID_MAIN_PATH)/lidort_lp_bvproblem.f90       \
   $(LID_MAIN_PATH)/lidort_lp_PostProcessing.f90  \
   $(LID_MAIN_PATH)/lidort_lp_wfatmos.f90         \
   $(LID_MAIN_PATH)/lidort_lps_masters.f90

LCS_SOURCES += \
   $(SOURCES) \
   $(L_SOURCES) \
   $(LID_MAIN_PATH)/lidort_lc_miscsetups.f90      \
   $(LID_MAIN_PATH)/lidort_lc_corrections.f90     \
   $(LID_MAIN_PATH)/lidort_lc_bvproblem.f90       \
   $(LID_MAIN_PATH)/lidort_lc_PostProcessing.f90  \
   $(LID_MAIN_PATH)/lidort_lc_wfatmos.f90         \
   $(LID_MAIN_PATH)/lidort_lcs_masters.f90

# (Include supplement source files)
include $(SUP_PATH)/makefile.sup

# Main tests
# ----------

#  tests 1-3

SOURCES_SOLAR = $(SOURCES) +  \
   $(LID_TEST_PATH)/3p7T_solar_tester.f90

SOURCES_SOLAR_LPS = $(LPS_SOURCES) +  \
   $(LID_TEST_PATH)/3p7T_solar_lps_tester.f90

SOURCES_SOLAR_LCS = $(LCS_SOURCES) +  \
   $(LID_TEST_PATH)/3p7T_solar_lcs_tester.f90

#  tests 4-6

SOURCES_THERMAL = $(SOURCES) +  \
   $(BRDF_SUP_SOURCES) \
   $(LID_TEST_PATH)/3p7T_thermal_tester.f90

SOURCES_THERMAL_LPS = $(LPS_SOURCES) +  \
   $(BRDF_LINSUP_SOURCES) \
   $(LID_TEST_PATH)/3p7T_thermal_lps_tester.f90

SOURCES_THERMAL_LCS = $(LCS_SOURCES) +  \
   $(BRDF_LINSUP_SOURCES) \
   $(LID_TEST_PATH)/3p7T_thermal_lcs_tester.f90

#  tests 8-9

SOURCES_BRDF_SELF = $(BASE_SOURCES) + \
   $(BRDF_SUP_SOURCES)              \
   $(BRDF_LINSUP_SOURCES)           \
   $(LID_TEST_PATH)/3p7T_brdf_self_tester.f90

SOURCES_BRDF_FULL = $(LCS_SOURCES) +  \
   $(SLEAVE_SUP_SOURCES) \
   $(BRDF_LINSUP_SOURCES) \
   $(LID_TEST_PATH)/lidort_3p7_sup_accessories.f90 \
   $(LID_TEST_PATH)/3p7T_brdf_full_tester.f90

#  tests 13-15 (OMP)

SOURCES_SOLAR_OMP = $(SOURCES) +  \
   $(LID_TEST_PATH)/3p7T_solar_OMP_tester.f90

SOURCES_SOLAR_LPS_OMP = $(LPS_SOURCES) +  \
   $(LID_TEST_PATH)/3p7T_solar_lps_OMP_tester.f90

SOURCES_SOLAR_LCS_OMP = $(LCS_SOURCES) +  \
   $(LID_TEST_PATH)/3p7T_solar_lcs_OMP_tester.f90


#  SPECIAL TESTS (Need to use dedicated PARS files first)
#  -------------

# Tests 10-12

SOURCES_SLEAVE_SELF = $(BASE_SOURCES) + \
   $(SLEAVE_SUP_SOURCES)              \
   $(SLEAVE_LINSUP_SOURCES)           \
   $(LID_TEST_PATH)/3p7T_sleave_self_tester.f90

SOURCES_SLEAVE_FULL = $(LCS_SOURCES) +  \
   $(SLEAVE_LINSUP_SOURCES) \
   $(BRDF_SUP_SOURCES) \
   $(LID_TEST_PATH)/lidort_3p7_sup_accessories.f90 \
   $(LID_TEST_PATH)/3p7T_sleave_full_tester.f90

SOURCES_THERMAL_LBBF = $(LPS_SOURCES) +  \
   $(BRDF_LINSUP_SOURCES) \
   $(LID_TEST_PATH)/3p7T_thermal_lbbf_tester.f90

# Utilities

SOURCES_UTIL_LID =
SOURCES_UTIL_LID +=   \
   $(UTIL_PATH)/lidort_diff.f90

SOURCES_UTIL_SLV =
SOURCES_UTIL_SLV +=   \
   $(UTIL_PATH)/sleave_diff.f90

#
# Define pattern rules for creating object files:
#

#.SUFFIXES:

# For lidort main source files
$(OBJ_PATH)/%.o : $(LID_DEF_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(LID_MAIN_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(LID_TEST_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@

# For utility source files
$(OBJ_PATH)/%.o : $(UTIL_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@


#
# Define object files
#  Main tests
#  ----------

#  tests 1-3

F90SOURCES_SOLAR := $(notdir $(filter %.f90, $(SOURCES_SOLAR)))
F90OBJECTS_SOLAR := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_SOLAR)))

F90SOURCES_SOLAR_LPS := $(notdir $(filter %.f90, $(SOURCES_SOLAR_LPS)))
F90OBJECTS_SOLAR_LPS := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_SOLAR_LPS)))

F90SOURCES_SOLAR_LCS := $(notdir $(filter %.f90, $(SOURCES_SOLAR_LCS)))
F90OBJECTS_SOLAR_LCS := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_SOLAR_LCS)))

#  tests 4-6

F90SOURCES_THERMAL := $(notdir $(filter %.f90, $(SOURCES_THERMAL)))
F90OBJECTS_THERMAL := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_THERMAL)))

F90SOURCES_THERMAL_LPS := $(notdir $(filter %.f90, $(SOURCES_THERMAL_LPS)))
F90OBJECTS_THERMAL_LPS := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_THERMAL_LPS)))

F90SOURCES_THERMAL_LCS := $(notdir $(filter %.f90, $(SOURCES_THERMAL_LCS)))
F90OBJECTS_THERMAL_LCS := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_THERMAL_LCS)))

#  tests 8-9

F90SOURCES_BRDF_SELF := $(notdir $(filter %.f90, $(SOURCES_BRDF_SELF)))
F90OBJECTS_BRDF_SELF := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_BRDF_SELF)))

F90SOURCES_BRDF_FULL := $(notdir $(filter %.f90, $(SOURCES_BRDF_FULL)))
F90OBJECTS_BRDF_FULL := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_BRDF_FULL)))

#  tests 13-15 (OMP)

F90SOURCES_SOLAR_OMP := $(notdir $(filter %.f90, $(SOURCES_SOLAR_OMP)))
F90OBJECTS_SOLAR_OMP := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_SOLAR_OMP)))

F90SOURCES_SOLAR_LPS_OMP := $(notdir $(filter %.f90, $(SOURCES_SOLAR_LPS_OMP)))
F90OBJECTS_SOLAR_LPS_OMP := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_SOLAR_LPS_OMP)))

F90SOURCES_SOLAR_LCS_OMP := $(notdir $(filter %.f90, $(SOURCES_SOLAR_LCS_OMP)))
F90OBJECTS_SOLAR_LCS_OMP := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_SOLAR_LCS_OMP)))

#  Special tests
#  -------------

#  tests 10-12

F90SOURCES_SLEAVE_SELF := $(notdir $(filter %.f90, $(SOURCES_SLEAVE_SELF)))
F90OBJECTS_SLEAVE_SELF := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_SLEAVE_SELF)))

F90SOURCES_SLEAVE_FULL := $(notdir $(filter %.f90, $(SOURCES_SLEAVE_FULL)))
F90OBJECTS_SLEAVE_FULL := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_SLEAVE_FULL)))

F90SOURCES_THERMAL_LBBF := $(notdir $(filter %.f90, $(SOURCES_THERMAL_LBBF)))
F90OBJECTS_THERMAL_LBBF := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_THERMAL_LBBF)))

#  Utilities

F90SOURCES_UTIL_LID := $(notdir $(filter %.f90, $(SOURCES_UTIL_LID)))
F90OBJECTS_UTIL_LID := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_UTIL_LID)))

F90SOURCES_UTIL_SLV := $(notdir $(filter %.f90, $(SOURCES_UTIL_SLV)))
F90OBJECTS_UTIL_SLV := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_UTIL_SLV)))

#
# Define desired targets
#

main: solar   \
      thermal \
      brdf

# Tests 1-3

solar: 3p7T_solar_tester.exe \
     3p7T_solar_lps_tester.exe \
     3p7T_solar_lcs_tester.exe

# Tests 4-6

thermal: 3p7T_thermal_tester.exe \
     3p7T_thermal_lps_tester.exe \
     3p7T_thermal_lcs_tester.exe

# Test 7 ( PLACEHOLDER )

# Tests 8-9

brdf: 3p7T_brdf_self_tester.exe \
      3p7T_brdf_full_tester.exe

# Tests 10-12 ( SEE BELOW )

# Tests 13-15

OMP: 3p7T_solar_OMP_tester.exe \
     3p7T_solar_lps_OMP_tester.exe \
     3p7T_solar_lcs_OMP_tester.exe

# Tests 16-18 ( SEE BELOW )

#  selected tests with OBSERVATIONAL GEOMETRY option

obsgeo: 3p7T_solar_lps_tester.exe \
        3p7T_solar_lcs_tester.exe \
        3p7T_brdf_full_tester.exe

obsgeo_sleave: 3p7T_sleave_full_tester.exe

# Special Tests
# -------------

# Tests 10-12 not included in "main" or "obsgeo" (require different "lidort_pars.f90" settings)
# Tests 16-18 not included in "main" for NO REAL REASON

# Tests 10-11

sleave: 3p7T_sleave_self_tester.exe \
        3p7T_sleave_full_tester.exe

# Test 12

lbbf:   3p7T_thermal_lbbf_tester.exe

# Tests 16-18

Taylor: 3p7T_Taylor_tester.exe     \
        3p7T_Taylor_lcs_tester.exe \
        3p7T_Taylor_lps_tester.exe

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#  Main tests
#  ----------

#  Tests 1-3

3p7T_solar_tester.exe: $(F90OBJECTS_SOLAR)
	$(FC) $^ -o $@
3p7T_solar_lps_tester.exe: $(F90OBJECTS_SOLAR_LPS)
	$(FC) $^ -o $@
3p7T_solar_lcs_tester.exe: $(F90OBJECTS_SOLAR_LCS)
	$(FC) $^ -o $@

#  Tests 4-6

3p7T_thermal_tester.exe: $(F90OBJECTS_THERMAL)
	$(FC) $^ -o $@
3p7T_thermal_lps_tester.exe: $(F90OBJECTS_THERMAL_LPS)
	$(FC) $^ -o $@
3p7T_thermal_lcs_tester.exe: $(F90OBJECTS_THERMAL_LCS)
	$(FC) $^ -o $@

#  Tests 8-9

3p7T_brdf_self_tester.exe: $(F90OBJECTS_BRDF_SELF)
	$(FC) $^ -o $@
3p7T_brdf_full_tester.exe: $(F90OBJECTS_BRDF_FULL)
	$(FC) $^ -o $@

#  Tests 13-15

3p7T_solar_OMP_tester.exe: $(F90OBJECTS_SOLAR_OMP)
	$(FC) $(FFLAGS_OPENMP) $^ -o $@
3p7T_solar_lps_OMP_tester.exe: $(F90OBJECTS_SOLAR_LPS_OMP)
	$(FC) $(FFLAGS_OPENMP) $^ -o $@
3p7T_solar_lcs_OMP_tester.exe: $(F90OBJECTS_SOLAR_LCS_OMP)
	$(FC) $(FFLAGS_OPENMP) $^ -o $@

#3p7T_solar_OMP_tester.exe: $(F90OBJECTS_SOLAR_OMP)
#	$(FC) $^ -o $@
#3p7T_solar_lps_OMP_tester.exe: $(F90OBJECTS_SOLAR_LPS_OMP)
#	$(FC) $^ -o $@
#3p7T_solar_lcs_OMP_tester.exe: $(F90OBJECTS_SOLAR_LCS_OMP)
#	$(FC) $^ -o $@

#  specials
#  --------

# Tests 10-12

3p7T_sleave_self_tester.exe: $(F90OBJECTS_SLEAVE_SELF)
	$(FC) $^ -o $@
3p7T_sleave_full_tester.exe: $(F90OBJECTS_SLEAVE_FULL)
	$(FC) $^ -o $@
3p7T_thermal_lbbf_tester.exe: $(F90OBJECTS_THERMAL_LBBF)
	$(FC) $^ -o $@

# Diff utilities

lidort_diff: $(F90OBJECTS_UTIL_LID)
	$(FC) $^ -o $@

sleave_diff: $(F90OBJECTS_UTIL_SLV)
	$(FC) $^ -o $@

# Make clean

.PHONY: clean
clean:
	rm -f *.o $(OBJ_FILES) *.mod $(MOD_FILES) *.log *.exe

