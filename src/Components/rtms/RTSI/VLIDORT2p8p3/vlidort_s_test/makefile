#
# This file serves as the head Makefile for VLIDORT and its "tester" programs
#

# UPGRADE TO VLIDORT VERSION 2.7, May   2014
# UPGRADE TO VLIDORT VERSION 2.8, March 2017, together with FO code Version 1.5
# UPGRADE TO VLIDORT VERSION 2.8.1 AUGUST 2019
#      4/9/19.   Removal of TRANSFLUX MASTER Routine
#      4/19/19.  Introduction of mediaprops routines (3 in all)
#      5/1/19.   Supplement and FO codes no compilation changes
#      6/3/19.   Introduce 2 new special tests (Planetary/LWCoupling)
# UPGRADE TO VLIDORT VERSION 2.8.3, 31 March 2021

#
# Define some variables
#

UTIL_PATH       = util
VSUP_PATH       = vsup
FO_MAIN_PATH    = vlidort_focode
VLID_DEF_PATH   = vlidort_def
VLID_MAIN_PATH  = vlidort_main
VLID_MAIN_PATH1 = vlidort_main/regular
VLID_MAIN_PATH2 = vlidort_main/linearized
VLID_TEST_PATH  = vlidort_s_test

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

SOURCES +=   \
   $(BASE_SOURCES) \
   $(VLID_MAIN_PATH1)/lapack_tools.f90		 \
   $(VLID_DEF_PATH)/vlidort_inputs_def.f90	 \
   $(VLID_DEF_PATH)/vlidort_sup_brdf_def.f90	 \
   $(VLID_DEF_PATH)/vlidort_sup_ss_def.f90	 \
   $(VLID_DEF_PATH)/vlidort_sup_sleave_def.f90	 \
   $(VLID_DEF_PATH)/vlidort_sup_def.f90		 \
   $(VLID_DEF_PATH)/vlidort_outputs_def.f90	 \
   $(VLID_DEF_PATH)/vlidort_io_defs.f90		 \
   $(VLID_DEF_PATH)/vlidort_work_def.f90	 \
   $(VLID_MAIN_PATH1)/vlidort_aux.f90		 \
   $(VLID_MAIN_PATH1)/vlidort_getplanck.f90	 \
   $(VLID_MAIN_PATH1)/vlidort_geometry.f90       \
   $(VLID_MAIN_PATH1)/vlidort_Taylor.f90         \
   $(VLID_MAIN_PATH1)/vlidort_inputs.f90	 \
   $(VLID_MAIN_PATH1)/vlidort_miscsetups.f90	 \
   $(VLID_MAIN_PATH1)/vlidort_vfo_interface.f90	 \
   $(VLID_MAIN_PATH1)/vlidort_thermalsup.f90	 \
   $(VLID_MAIN_PATH1)/vlidort_solutions.f90	 \
   $(VLID_MAIN_PATH1)/vlidort_bvproblem.f90	 \
   $(VLID_MAIN_PATH1)/vlidort_PostProcessing.f90 \
   $(VLID_MAIN_PATH1)/vlidort_intensity.f90	 \
   $(VLID_MAIN_PATH1)/vlidort_converge.f90	 \
   $(VLID_MAIN_PATH1)/vlidort_writemodules.f90	 \
   $(VLID_MAIN_PATH1)/vlidort_pack.f90		 \
   $(VLID_MAIN_PATH1)/vlidort_unpack.f90	 \
   $(VLID_MAIN_PATH1)/vlidort_mediaprops.f90     \
   $(VLID_MAIN_PATH1)/vlidort_masters_V2p8p3.f90

L_SOURCES += \
   $(VLID_DEF_PATH)/vlidort_lin_inputs_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_sup_brdf_def.f90  \
   $(VLID_DEF_PATH)/vlidort_lin_sup_ss_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_sup_sleave_def.f90\
   $(VLID_DEF_PATH)/vlidort_lin_sup_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_outputs_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_io_defs.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_work_def.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_l_inputs.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_l_writemodules.f90  \
   $(VLID_MAIN_PATH2)/vlidort_la_miscsetups.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_l_pack.f90          \
   $(VLID_MAIN_PATH2)/vlidort_l_unpack.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_l_thermalsup.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_lpc_solutions.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_lpc_bvproblem.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_lbbf_jacobians_vector.f90 \
   $(VLID_MAIN_PATH2)/vlidort_ls_wfsurface.f90	  \
   $(VLID_MAIN_PATH2)/vlidort_ls_wfsleave.f90

LPS_SOURCES += \
   $(VLID_MAIN_PATH2)/vlidort_lp_miscsetups.f90      \
   $(VLID_MAIN_PATH2)/vlidort_vfo_lps_interface.f90  \
   $(VLID_MAIN_PATH2)/vlidort_lp_pack.f90	     \
   $(VLID_MAIN_PATH2)/vlidort_lp_unpack.f90	     \
   $(VLID_MAIN_PATH2)/vlidort_lp_mediaprops.f90      \
   $(VLID_MAIN_PATH2)/vlidort_lp_solutions.f90	     \
   $(VLID_MAIN_PATH2)/vlidort_lp_bvproblem.f90       \
   $(VLID_MAIN_PATH2)/vlidort_lp_PostProcessing.f90  \
   $(VLID_MAIN_PATH2)/vlidort_lp_wfatmos.f90         \
   $(VLID_MAIN_PATH2)/vlidort_lps_converge.f90	     \
   $(VLID_MAIN_PATH2)/vlidort_lps_masters_V2p8p3.f90

LCS_SOURCES += \
   $(VLID_MAIN_PATH2)/vlidort_lc_miscsetups.f90      \
   $(VLID_MAIN_PATH2)/vlidort_vfo_lcs_interface.f90  \
   $(VLID_MAIN_PATH2)/vlidort_lc_pack.f90	     \
   $(VLID_MAIN_PATH2)/vlidort_lc_unpack.f90	     \
   $(VLID_MAIN_PATH2)/vlidort_lc_mediaprops.f90      \
   $(VLID_MAIN_PATH2)/vlidort_lc_solutions.f90	     \
   $(VLID_MAIN_PATH2)/vlidort_lc_bvproblem.f90       \
   $(VLID_MAIN_PATH2)/vlidort_lc_PostProcessing.f90  \
   $(VLID_MAIN_PATH2)/vlidort_lc_wfatmos.f90         \
   $(VLID_MAIN_PATH2)/vlidort_lcs_converge.f90	     \
   $(VLID_MAIN_PATH2)/vlidort_lcs_masters_V2p8p3.f90

# (Include vector supplement source files)
include $(VSUP_PATH)/makefile.vsup

# (Include first-order source files)
include $(FO_MAIN_PATH)/makefile.fo_vlidort

# Main scalar tests

SOURCES_SOLAR = $(FO_SOURCES_Vector) + \
   $(SOURCES) \
   $(VLID_TEST_PATH)/2p8p3_solar_tester.f90

SOURCES_THERMAL = $(FO_SOURCES_Vector) + \
   $(SOURCES) \
   $(VBRDF_SUP_SOURCES)   \
   $(VBRDF_SUP_ACCESS_SOURCES) \
   $(VLID_TEST_PATH)/2p8p3_thermal_tester.f90

SOURCES_BRDF_SELF = $(VBRDF_LINSUP_SOURCES) + \
   $(VLID_TEST_PATH)/2p8p3_brdf_self_tester.f90

SOURCES_SOLAR_LPCS = $(FO_SOURCES_L_Vector) + \
   $(SOURCES) \
   $(L_SOURCES) \
   $(LPS_SOURCES) \
   $(LCS_SOURCES) \
   $(VLID_TEST_PATH)/2p8p3_solar_lpcs_tester.f90

SOURCES_THERMAL_LPCS = $(FO_SOURCES_L_Vector) + \
   $(SOURCES) \
   $(L_SOURCES) \
   $(LPS_SOURCES) \
   $(LCS_SOURCES) \
   $(VLID_TEST_PATH)/2p8p3_thermal_lpcs_tester.f90

SOURCES_BRDFPLUS = $(FO_SOURCES_L_Vector) + \
   $(SOURCES) \
   $(L_SOURCES) \
   $(LCS_SOURCES) \
   $(VBRDF_LINSUP_SOURCES) \
   $(VBRDF_LINSUP_ACCESS_SOURCES) \
   $(VLID_TEST_PATH)/2p8p3_brdfplus_tester.f90

# Special scalar tests

SOURCES_VFZMAT = $(FO_SOURCES_Vector) + \
   $(SOURCES) \
   $(VFZMAT_SUP_SOURCES)   \
   $(VLID_TEST_PATH)/2p8p3_vfzmat_tester.f90

SOURCES_VSLEAVE_SELF = $(VSLEAVE_LINSUP_SOURCES) + \
   $(VLID_TEST_PATH)/2p8p3_vsleave_self_tester.f90

#  Special scalar test for Planetary problem, Version 2.8.1. 03 June 2019

SOURCES_FULL_PLANETARY = $(FO_SOURCES_L_Vector) + \
   $(SOURCES) \
   $(L_SOURCES)           \
   $(LCS_SOURCES)         \
   $(LPS_SOURCES)         \
   $(VLID_TEST_PATH)/2p8p3_Planetary_tester.f90

#  Special test for the LW Coupling

SOURCES_LWCoupling = $(FO_SOURCES_Vector) + \
   $(SOURCES)           \
   $(VSLEAVE_SUP_SOURCES) \
   $(VSLEAVE_SUP_ACCESS_SOURCES) \
   $(VLID_TEST_PATH)/2p8p3_LWCoupling_tester.f90

# Utilities

SOURCES_UTIL =
SOURCES_UTIL +=   \
   $(UTIL_PATH)/vlidort_diff.f90

#
# Define pattern rules for creating object files:
#

#.SUFFIXES:

# For vlidort main source files
$(OBJ_PATH)/%.o : $(VLID_DEF_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_MAIN_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_MAIN_PATH1)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_MAIN_PATH2)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_TEST_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@

# For utility source files
$(OBJ_PATH)/%.o : $(UTIL_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@

#
# Define object files
#

F90SOURCES_SOLAR := $(notdir $(filter %.f90, $(SOURCES_SOLAR)))
F90OBJECTS_SOLAR := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_SOLAR)))

F90SOURCES_THERMAL := $(notdir $(filter %.f90, $(SOURCES_THERMAL)))
F90OBJECTS_THERMAL := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_THERMAL)))

F90SOURCES_BRDF_SELF := $(notdir $(filter %.f90, $(SOURCES_BRDF_SELF)))
F90OBJECTS_BRDF_SELF := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_BRDF_SELF)))

F90SOURCES_SOLAR_LPCS := $(notdir $(filter %.f90, $(SOURCES_SOLAR_LPCS)))
F90OBJECTS_SOLAR_LPCS := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_SOLAR_LPCS)))

F90SOURCES_THERMAL_LPCS := $(notdir $(filter %.f90, $(SOURCES_THERMAL_LPCS)))
F90OBJECTS_THERMAL_LPCS := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_THERMAL_LPCS)))

F90SOURCES_BRDFPLUS := $(notdir $(filter %.f90, $(SOURCES_BRDFPLUS)))
F90OBJECTS_BRDFPLUS := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_BRDFPLUS)))


F90SOURCES_VFZMAT := $(notdir $(filter %.f90, $(SOURCES_VFZMAT)))
F90OBJECTS_VFZMAT := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_VFZMAT)))

F90SOURCES_VSLEAVE_SELF := $(notdir $(filter %.f90, $(SOURCES_VSLEAVE_SELF)))
F90OBJECTS_VSLEAVE_SELF := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_VSLEAVE_SELF)))


F90SOURCES_FULL_PLANETARY := $(notdir $(filter %.f90, $(SOURCES_FULL_PLANETARY)))
F90OBJECTS_FULL_PLANETARY := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_FULL_PLANETARY)))

F90SOURCES_LWCoupling := $(notdir $(filter %.f90, $(SOURCES_LWCoupling)))
F90OBJECTS_LWCoupling := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_LWCoupling)))


F90SOURCES_UTIL := $(notdir $(filter %.f90, $(SOURCES_UTIL)))
F90OBJECTS_UTIL := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_UTIL)))


#
# Define desired targets
#

main: solar   \
      thermal \
      brdf

solar: s2p8p3_solar_tester.exe \
       s2p8p3_solar_lpcs_tester.exe

thermal: s2p8p3_thermal_tester.exe \
         s2p8p3_thermal_lpcs_tester.exe

brdf: s2p8p3_brdf_self_tester.exe \
      s2p8p3_brdfplus_tester.exe

Planetary: s2p8p3_Planetary_tester.exe

LWCoupling: s2p8p3_LWCoupling_tester.exe

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

s2p8p3_solar_tester.exe: $(F90OBJECTS_SOLAR)
	$(FC) $^ -o $@
s2p8p3_thermal_tester.exe: $(F90OBJECTS_THERMAL)
	$(FC) $^ -o $@
s2p8p3_brdf_self_tester.exe: $(F90OBJECTS_BRDF_SELF)
	$(FC) $^ -o $@
s2p8p3_solar_lpcs_tester.exe: $(F90OBJECTS_SOLAR_LPCS)
	$(FC) $^ -o $@
s2p8p3_thermal_lpcs_tester.exe: $(F90OBJECTS_THERMAL_LPCS)
	$(FC) $^ -o $@
s2p8p3_brdfplus_tester.exe: $(F90OBJECTS_BRDFPLUS)
	$(FC) $^ -o $@
s2p8p3_vfzmat_tester.exe: $(F90OBJECTS_VFZMAT)
	$(FC) $^ -o $@

# Note: These tests not included in "all" because they require different
#       "vlidort_pars.f90" settings

s2p8p3_vsleave_self_tester.exe: $(F90OBJECTS_VSLEAVE_SELF)
	$(FC) $^ -o $@
s2p8p3_Planetary_tester.exe: $(F90OBJECTS_FULL_PLANETARY)
	$(FC) $^ -o $@
s2p8p3_LWCoupling_tester.exe : $(F90OBJECTS_LWCoupling)
	$(FC) $^ -o $@


vlidort_diff: $(F90OBJECTS_UTIL)
	$(FC) $^ -o $@

.PHONY: clean
clean:
	rm -f *.o $(OBJ_FILES) *.mod $(MOD_FILES) *.log *.exe

