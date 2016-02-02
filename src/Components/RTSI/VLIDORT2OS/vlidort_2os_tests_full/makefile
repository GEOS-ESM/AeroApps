#
# This file serves as the head Makefile for VLIDORT and its "tester" programs
#


#
# Define some variables
#

UTIL_PATH      = ../util
VSUP_PATH      = ../vsup
VLID_DEF_PATH  = ../vlidort_def
VLID_MAIN_PATH = ../vlidort_main
VLID_TEST_PATH = .
FO_MAIN_PATH   = ../fo_main

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
	FFLAGS_OPT = -O3
        FFLAGS_OPENMP = -openmp
endif

# Additional flags for gfortran
ifeq ($(FC), gfortran)
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -J$(MOD_PATH)

#	FFLAGS_DEBUG = -g -Wall -fbounds-check
#	FFLAGS_DEBUG = -g -Wall -fbounds-check -fbacktrace
	FFLAGS_DEBUG = -g -C -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
#	FFLAGS_DEBUG = -fstack-usage
	FFLAGS_OPT = -O3
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
#	FFLAGS_DEBUG = -g -C=all -C=undefined -gline -mtrace=all -nan
#	FFLAGS_DEBUG = -w=obs -w=unused -g -C=all -C=undefined -gline -mtrace=all -nan
	FFLAGS_OPT = -O3
        FFLAGS_OPENMP = -openmp
endif

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
BASE_SOURCES +=   \
   $(VLID_DEF_PATH)/vlidort_pars.f90

SOURCES_V2OS = 
SOURCES_V2OS += \
   $(BASE_SOURCES) \
   $(VLID_MAIN_PATH)/lapack_tools.f90		\
   $(VLID_DEF_PATH)/vlidort_inputs_def.f90	\
   $(VLID_DEF_PATH)/vlidort_sup_brdf_def.f90	\
   $(VLID_DEF_PATH)/vlidort_sup_ss_def.f90	\
   $(VLID_DEF_PATH)/vlidort_sup_sleave_def.f90	\
   $(VLID_DEF_PATH)/vlidort_sup_def.f90		\
   $(VLID_DEF_PATH)/vlidort_outputs_def.f90	\
   $(VLID_DEF_PATH)/vlidort_io_defs.f90		\
   $(VLID_DEF_PATH)/vlidort_work_def.f90	\
   $(VLID_MAIN_PATH)/vlidort_aux.f90		\
   $(VLID_MAIN_PATH)/vlidort_getplanck.f90	\
   $(VLID_MAIN_PATH)/vlidort_geometry.f90       \
   $(VLID_MAIN_PATH)/vlidort_Taylor.f90         \
   $(VLID_MAIN_PATH)/vlidort_inputs.f90		\
   $(VLID_MAIN_PATH)/vlidort_miscsetups.f90	\
   $(VLID_MAIN_PATH)/vlidort_multipliers.f90	\
   $(VLID_MAIN_PATH)/vlidort_corrections.f90	\
   $(VLID_MAIN_PATH)/vlidort_thermalsup.f90	\
   $(VLID_MAIN_PATH)/vlidort_solutions.f90	\
   $(VLID_MAIN_PATH)/vlidort_bvproblem.f90	\
   $(VLID_MAIN_PATH)/vlidort_intensity.f90	\
   $(VLID_MAIN_PATH)/vlidort_writemodules.f90	\
   $(VLID_MAIN_PATH)/vlidort_pack.f90		\
   $(VLID_MAIN_PATH)/vlidort_unpack.f90		\
   $(VLID_MAIN_PATH)/vlidort_2OScorr_aux.f90       \
   $(VLID_MAIN_PATH)/vlidort_2OScorr_utilities.f90 \
   $(VLID_MAIN_PATH)/vlidort_2OScorr_routines.f90  \
   $(VLID_MAIN_PATH)/vlidort_2OScorr_master.f90    \
   $(VLID_MAIN_PATH)/vlidort_masters.f90

L_SOURCES_V2OS =
L_SOURCES_V2OS += \
   $(VLID_DEF_PATH)/vlidort_lin_inputs_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_sup_brdf_def.f90  \
   $(VLID_DEF_PATH)/vlidort_lin_sup_ss_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_sup_sleave_def.f90\
   $(VLID_DEF_PATH)/vlidort_lin_sup_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_outputs_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_io_defs.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_work_def.f90	  \
   $(VLID_MAIN_PATH)/vlidort_l_inputs.f90	  \
   $(VLID_MAIN_PATH)/vlidort_la_miscsetups.f90	  \
   $(VLID_MAIN_PATH)/vlidort_la_corrections.f90	  \
   $(VLID_MAIN_PATH)/vlidort_ls_corrections.f90   \
   $(VLID_MAIN_PATH)/vlidort_l_thermalsup.f90	  \
   $(VLID_MAIN_PATH)/vlidort_lpc_solutions.f90	  \
   $(VLID_MAIN_PATH)/vlidort_lpc_bvproblem.f90	  \
   $(VLID_MAIN_PATH)/vlidort_lbbf_jacobians.f90   \
   $(VLID_MAIN_PATH)/vlidort_ls_wfsurface.f90	  \
   $(VLID_MAIN_PATH)/vlidort_ls_wfsleave.f90	  \
   $(VLID_MAIN_PATH)/vlidort_l_writemodules.f90   \
   $(VLID_MAIN_PATH)/vlidort_l_pack.f90		  \
   $(VLID_MAIN_PATH)/vlidort_l_unpack.f90         \
   $(VLID_MAIN_PATH)/vlidort_2OScorr_la_routines.f90

LPS_SOURCES_V2OS =
LPS_SOURCES_V2OS += \
   $(VLID_MAIN_PATH)/vlidort_lp_miscsetups.f90    \
   $(VLID_MAIN_PATH)/vlidort_lp_corrections.f90   \
   $(VLID_MAIN_PATH)/vlidort_lp_solutions.f90	  \
   $(VLID_MAIN_PATH)/vlidort_lp_bvproblem.f90     \
   $(VLID_MAIN_PATH)/vlidort_lp_wfatmos.f90       \
   $(VLID_MAIN_PATH)/vlidort_lp_pack.f90	  \
   $(VLID_MAIN_PATH)/vlidort_lp_unpack.f90	  \
   $(VLID_MAIN_PATH)/vlidort_2OScorr_lps_routines.f90 \
   $(VLID_MAIN_PATH)/vlidort_2OScorr_lps_master.f90   \
   $(VLID_MAIN_PATH)/vlidort_lps_masters.f90

LCS_SOURCES_V2OS =
LCS_SOURCES_V2OS += \
   $(VLID_MAIN_PATH)/vlidort_lc_miscsetups.f90    \
   $(VLID_MAIN_PATH)/vlidort_lc_corrections.f90   \
   $(VLID_MAIN_PATH)/vlidort_lc_solutions.f90	  \
   $(VLID_MAIN_PATH)/vlidort_lc_bvproblem.f90     \
   $(VLID_MAIN_PATH)/vlidort_lc_wfatmos.f90       \
   $(VLID_MAIN_PATH)/vlidort_lc_pack.f90	  \
   $(VLID_MAIN_PATH)/vlidort_lc_unpack.f90	  \
   $(VLID_MAIN_PATH)/vlidort_2OScorr_lcs_routines.f90 \
   $(VLID_MAIN_PATH)/vlidort_2OScorr_lcs_master.f90   \
   $(VLID_MAIN_PATH)/vlidort_lcs_masters.f90

# (Include vector supplement source files)
include $(VSUP_PATH)/makefile.vsup

# (Include first-order source files)
include $(FO_MAIN_PATH)/makefile.fo

# Main vector tests

SOURCES_V2OS_FULL_TEST = $(FO_SOURCES_Vector) + \
   $(SOURCES_V2OS) \
   $(VLID_TEST_PATH)/V2OS_full_tester.f90

SOURCES_V2OS_FULL_LPCSTEST = $(FO_SOURCES_L_Vector) + \
   $(SOURCES_V2OS)     \
   $(L_SOURCES_V2OS)   \
   $(LPS_SOURCES_V2OS) \
   $(LCS_SOURCES_V2OS) \
   $(VLID_TEST_PATH)/V2OS_full_lpcs_tester.f90

SOURCES_V2OS_FULL_OMPTEST =  $(FO_SOURCES_Vector) + \
   $(SOURCES_V2OS) \
   $(VLID_TEST_PATH)/V2OS_full_OMPtester.f90

SOURCES_V2OS_FULL_LPCSOMPTEST = $(FO_SOURCES_L_Vector) + \
   $(SOURCES_V2OS)     \
   $(L_SOURCES_V2OS)   \
   $(LPS_SOURCES_V2OS) \
   $(LCS_SOURCES_V2OS) \
   $(VLID_TEST_PATH)/V2OS_full_lpcs_OMPtester.f90

#
# Define pattern rules for creating object files:
#

#.SUFFIXES:

# For vlidort main source files
$(OBJ_PATH)/%.o : $(VLID_DEF_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_MAIN_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_TEST_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@

#
# Define object files
#

F90SOURCES_V2OS_FULL_TEST := $(notdir $(filter %.f90, $(SOURCES_V2OS_FULL_TEST)))
F90OBJECTS_V2OS_FULL_TEST := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_V2OS_FULL_TEST)))

F90SOURCES_V2OS_FULL_LPCSTEST := $(notdir $(filter %.f90, $(SOURCES_V2OS_FULL_LPCSTEST)))
F90OBJECTS_V2OS_FULL_LPCSTEST := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_V2OS_FULL_LPCSTEST)))

F90SOURCES_V2OS_FULL_OMPTEST := $(notdir $(filter %.f90, $(SOURCES_V2OS_FULL_OMPTEST)))
F90OBJECTS_V2OS_FULL_OMPTEST := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_V2OS_FULL_OMPTEST)))

F90SOURCES_V2OS_FULL_LPCSOMPTEST := $(notdir $(filter %.f90, $(SOURCES_V2OS_FULL_LPCSOMPTEST)))
F90OBJECTS_V2OS_FULL_LPCSOMPTEST := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_V2OS_FULL_LPCSOMPTEST)))

#
# Define desired targets
#

main  :  V2OS_full_tester.exe     V2OS_full_lpcs_tester.exe
openmp:  V2OS_full_OMPtester.exe  V2OS_full_lpcs_OMPtester.exe 

V2OS_full_tester.exe: $(F90OBJECTS_V2OS_FULL_TEST)
	$(FC) $^ -o $@
V2OS_full_lpcs_tester.exe: $(F90OBJECTS_V2OS_FULL_LPCSTEST)
	$(FC) $^ -o $@

V2OS_full_OMPtester.exe: $(F90OBJECTS_V2OS_FULL_OMPTEST)
	$(FC) $(FFLAGS_OPENMP) $^ -o $@
V2OS_full_lpcs_OMPtester.exe: $(F90OBJECTS_V2OS_FULL_LPCSOMPTEST)
	$(FC) $(FFLAGS_OPENMP) $^ -o $@
#V2OS_full_OMPtester.exe: $(F90OBJECTS_V2OS_FULL_OMPTEST)
#	$(FC) $^ -o $@
#V2OS_full_lpcs_OMPtester.exe: $(F90OBJECTS_V2OS_FULL_LPCSOMPTEST)
#	$(FC) $^ -o $@

.PHONY: clean
clean:
	rm -f *.o $(OBJ_FILES) *.mod $(MOD_FILES) *.log *.exe

