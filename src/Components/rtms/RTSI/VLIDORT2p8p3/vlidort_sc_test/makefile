#
# This file serves as the head Makefile for VLIDORT and its "tester" programs
#

# FOR VLIDORT VERSION 2.8.3, 31 March 2021
#     R. Spurr. RT Solutions Inc.

#
# Define some variables
#

UTIL_PATH       = util
SUP_PATH        = vsup
FO_MAIN_PATH    = vlidort_focode
VLID_DEF_PATH   = vlidort_def
VLID_MAIN_PATH  = vlidort_main
VLID_MAIN_PATH1 = vlidort_main/regular
VLID_MAIN_PATH2 = vlidort_main/linearized
VLID_TEST_PATH  = vlidort_sc_test

VLID_INTFC_PATH = $(VLID_TEST_PATH)/sphercorr_intfc

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
#FFLAGS = -c -mcmodel=large


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
#	FFLAGS_DEBUG = -w=obs -w=unused -C=array -C=undefined -gline
	FFLAGS_DEBUG = -w=obs -w=unused -C=all -C=undefined -gline
#	FFLAGS_DEBUG = -g -C=all -C=undefined -gline -mtrace=all -nan
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
   $(VLID_DEF_PATH)/vlidort_pars.f90

#BASE_SOURCES +=   \
#   $(VLID_MAIN_PATH)/vlidort_testpts.f90

#  1/31/21. Some changes for Version 2.8.3

SOURCES +=   \
   $(BASE_SOURCES) \
   $(VLID_MAIN_PATH1)/lapack_tools.f90		\
   $(VLID_DEF_PATH)/vlidort_inputs_def.f90	\
   $(VLID_DEF_PATH)/vlidort_sup_brdf_def.f90	\
   $(VLID_DEF_PATH)/vlidort_sup_ss_def.f90	\
   $(VLID_DEF_PATH)/vlidort_sup_sleave_def.f90	\
   $(VLID_DEF_PATH)/vlidort_sup_def.f90		\
   $(VLID_DEF_PATH)/vlidort_outputs_def.f90	\
   $(VLID_DEF_PATH)/vlidort_work_def.f90	\
   $(VLID_DEF_PATH)/vlidort_io_defs.f90		\
   $(VLID_MAIN_PATH1)/vlidort_aux.f90		\
   $(VLID_MAIN_PATH1)/vlidort_geometry.f90      \
   $(VLID_MAIN_PATH1)/vlidort_Taylor.f90        \
   $(VLID_MAIN_PATH1)/vlidort_pack.f90	 	\
   $(VLID_MAIN_PATH1)/vlidort_unpack.f90	\
   $(VLID_MAIN_PATH1)/vlidort_inputs.f90	\
   $(VLID_MAIN_PATH1)/vlidort_miscsetups.f90	\
   $(VLID_MAIN_PATH1)/vlidort_PostProcessing.f90\
   $(VLID_MAIN_PATH1)/vlidort_vfo_interface.f90	\
   $(VLID_MAIN_PATH1)/vlidort_thermalsup.f90	\
   $(VLID_MAIN_PATH1)/vlidort_solutions.f90	\
   $(VLID_MAIN_PATH1)/vlidort_bvproblem.f90	\
   $(VLID_MAIN_PATH1)/vlidort_intensity.f90	\
   $(VLID_MAIN_PATH1)/vlidort_converge.f90	\
   $(VLID_MAIN_PATH1)/vlidort_writemodules.f90	\
   $(VLID_MAIN_PATH1)/vlidort_mediaprops.f90    \
   $(VLID_MAIN_PATH1)/vlidort_masters_V2p8p3.f90

L_SOURCES += \
   $(VLID_DEF_PATH)/vlidort_lin_inputs_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_sup_brdf_def.f90  \
   $(VLID_DEF_PATH)/vlidort_lin_sup_ss_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_sup_sleave_def.f90\
   $(VLID_DEF_PATH)/vlidort_lin_sup_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_outputs_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_work_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_io_defs.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_l_inputs.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_l_pack.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_l_unpack.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_la_miscsetups.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_l_thermalsup.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_lpc_solutions.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_lpc_bvproblem.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_lbbf_jacobians_vector.f90 \
   $(VLID_MAIN_PATH2)/vlidort_ls_wfsurface.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_ls_wfsleave.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_l_writemodules.f90

LPS_SOURCES += \
   $(VLID_MAIN_PATH2)/vlidort_lp_miscsetups.f90    \
   $(VLID_MAIN_PATH2)/vlidort_lp_pack.f90	   \
   $(VLID_MAIN_PATH2)/vlidort_lp_unpack.f90	   \
   $(VLID_MAIN_PATH2)/vlidort_lp_solutions.f90	   \
   $(VLID_MAIN_PATH2)/vlidort_lp_bvproblem.f90     \
   $(VLID_MAIN_PATH2)/vlidort_lp_mediaprops.f90    \
   $(VLID_MAIN_PATH2)/vlidort_lp_PostProcessing.f90\
   $(VLID_MAIN_PATH2)/vlidort_lp_wfatmos.f90       \
   $(VLID_MAIN_PATH2)/vlidort_lps_converge.f90	   \
   $(VLID_MAIN_PATH2)/vlidort_vfo_lps_interface.f90\
   $(VLID_MAIN_PATH2)/vlidort_lps_masters_V2p8p3.f90

LCS_SOURCES += \
   $(VLID_MAIN_PATH2)/vlidort_lc_miscsetups.f90    \
   $(VLID_MAIN_PATH2)/vlidort_lc_pack.f90	   \
   $(VLID_MAIN_PATH2)/vlidort_lc_unpack.f90	   \
   $(VLID_MAIN_PATH2)/vlidort_lc_solutions.f90	   \
   $(VLID_MAIN_PATH2)/vlidort_lc_bvproblem.f90     \
   $(VLID_MAIN_PATH2)/vlidort_lc_mediaprops.f90    \
   $(VLID_MAIN_PATH2)/vlidort_lc_PostProcessing.f90\
   $(VLID_MAIN_PATH2)/vlidort_lc_wfatmos.f90       \
   $(VLID_MAIN_PATH2)/vlidort_lcs_converge.f90	   \
   $(VLID_MAIN_PATH2)/vlidort_vfo_lcs_interface.f90\
   $(VLID_MAIN_PATH2)/vlidort_lcs_masters_V2p8p3.f90

# (Include supplement source files)
#include $(SUP_PATH)/makefile.sup

# (Include first-order source files)
include $(FO_MAIN_PATH)/makefile.fo_vlidort

# Main spherical correction tests

SOURCES_ALL = $(FO_SOURCES_Vector) + \
   $(SOURCES) \
   $(VLID_INTFC_PATH)/SpherCorr_Geometry_Conversions.f90 \
   $(VLID_INTFC_PATH)/VLIDORT_SpherCorr_Routines.f90 \
   $(VLID_TEST_PATH)/VLIDORT_SpherCorr_tester_V3.f90

SOURCES_LCS_ALL = $(FO_SOURCES_L_Vector) + \
   $(SOURCES)      \
   $(L_SOURCES)    \
   $(LCS_SOURCES)  \
   $(VLID_INTFC_PATH)/SpherCorr_Geometry_Conversions.f90 \
   $(VLID_INTFC_PATH)/VLIDORT_LCS_SpherCorr_Routines.f90 \
   $(VLID_TEST_PATH)/VLIDORT_LCS_SpherCorr_tester_V2.f90

SOURCES_LPS_ALL = $(FO_SOURCES_L_Vector) + \
   $(SOURCES)      \
   $(L_SOURCES)    \
   $(LPS_SOURCES)  \
   $(VLID_INTFC_PATH)/SpherCorr_Geometry_Conversions.f90 \
   $(VLID_INTFC_PATH)/VLIDORT_LPS_SpherCorr_Routines.f90 \
   $(VLID_TEST_PATH)/VLIDORT_LPS_SpherCorr_tester_V2.f90

# Utilities

SOURCES_UTIL =
SOURCES_UTIL +=   \
   $(UTIL_PATH)/vlidort_diff.f90

#
# Define pattern rules for creating object files:
#

# For vlidort main source files
$(OBJ_PATH)/%.o : $(VLID_DEF_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_MAIN_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_MAIN_PATH1)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_MAIN_PATH2)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_INTFC_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_TEST_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@

# For utility source files
$(OBJ_PATH)/%.o : $(UTIL_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@

#
# Define object files
#

F90SOURCES_ALL := $(notdir $(filter %.f90, $(SOURCES_ALL)))
F90OBJECTS_ALL := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_ALL)))

F90SOURCES_LCS_ALL := $(notdir $(filter %.f90, $(SOURCES_LCS_ALL)))
F90OBJECTS_LCS_ALL := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_LCS_ALL)))

F90SOURCES_LPS_ALL := $(notdir $(filter %.f90, $(SOURCES_LPS_ALL)))
F90OBJECTS_LPS_ALL := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_LPS_ALL)))

#
# Define desired targets
#

Reg: V2p8p3_SpherCorr_tester_V3.exe
LCS: V2p8p3_LCS_SpherCorr_tester_V2.exe
LPS: V2p8p3_LPS_SpherCorr_tester_V2.exe

V2p8p3_SpherCorr_tester_V3.exe: $(F90OBJECTS_ALL)
	$(FC) $^ -o $@
V2p8p3_LCS_SpherCorr_tester_V2.exe: $(F90OBJECTS_LCS_ALL)
	$(FC) $^ -o $@
V2p8p3_LPS_SpherCorr_tester_V2.exe: $(F90OBJECTS_LPS_ALL)
	$(FC) $^ -o $@

vlidort_diff: $(F90OBJECTS_UTIL)
	$(FC) $^ -o $@

.PHONY: clean
clean:
	rm -f *.o $(OBJ_FILES) *.mod $(MOD_FILES) *.log *.exe

