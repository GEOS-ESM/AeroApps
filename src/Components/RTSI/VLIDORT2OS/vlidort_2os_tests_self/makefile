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

SOURCES_2OS = 
SOURCES_2OS += \
   $(BASE_SOURCES) \
   $(VLID_DEF_PATH)/vlidort_inputs_def.f90	\
   $(VLID_DEF_PATH)/vlidort_sup_brdf_def.f90	\
   $(VLID_DEF_PATH)/vlidort_sup_ss_def.f90	\
   $(VLID_DEF_PATH)/vlidort_sup_sleave_def.f90	\
   $(VLID_DEF_PATH)/vlidort_sup_def.f90		\
   $(VLID_DEF_PATH)/vlidort_outputs_def.f90	\
   $(VLID_DEF_PATH)/vlidort_io_defs.f90		\
   $(VLID_MAIN_PATH)/vlidort_aux.f90		\
   $(VLID_MAIN_PATH)/vlidort_2OScorr_aux.f90       \
   $(VLID_MAIN_PATH)/vlidort_2OScorr_utilities.f90 \
   $(VLID_MAIN_PATH)/vlidort_2OScorr_routines.f90  \
   $(VLID_MAIN_PATH)/vlidort_2OScorr_master.f90    \

LA_SOURCES_2OS += \
   $(VLID_DEF_PATH)/vlidort_lin_inputs_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_sup_brdf_def.f90  \
   $(VLID_DEF_PATH)/vlidort_lin_sup_ss_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_sup_sleave_def.f90\
   $(VLID_DEF_PATH)/vlidort_lin_sup_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_outputs_def.f90	  \
   $(VLID_DEF_PATH)/vlidort_lin_io_defs.f90	  \
   $(VLID_MAIN_PATH)/vlidort_2OScorr_la_routines.f90

LPCS_SOURCES_2OS+= \
   $(VLID_MAIN_PATH)/vlidort_2OScorr_lcs_routines.f90 \
   $(VLID_MAIN_PATH)/vlidort_2OScorr_lps_routines.f90 \
   $(VLID_MAIN_PATH)/vlidort_2OScorr_lps_master.f90   \
   $(VLID_MAIN_PATH)/vlidort_2OScorr_lcs_master.f90 

# Main vector tests

SOURCES_2OS_TEST = \
   $(SOURCES_2OS) \
   $(VLID_TEST_PATH)/TwoOS_self_tester.f90

SOURCES_2OS_LPCSTEST = \
   $(SOURCES_2OS) \
   $(LA_SOURCES_2OS) \
   $(LPCS_SOURCES_2OS) \
   $(VLID_TEST_PATH)/TwoOS_self_lpcs_tester.f90

SOURCES_2OS_OMPTEST = \
   $(SOURCES_2OS) \
   $(VLID_TEST_PATH)/TwoOS_self_OMPtester.f90

SOURCES_2OS_LPCSOMPTEST = \
   $(SOURCES_2OS) \
   $(LA_SOURCES_2OS) \
   $(LPCS_SOURCES_2OS) \
   $(VLID_TEST_PATH)/TwoOS_self_lpcs_OMPtester.f90

# Utilities

SOURCES_UTIL =
SOURCES_UTIL +=   \
   $(UTIL_PATH)/TwoOS_diff.f90

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

# For utility source files
$(OBJ_PATH)/%.o : $(UTIL_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@

#
# Define object files
#

F90SOURCES_2OS_TEST := $(notdir $(filter %.f90, $(SOURCES_2OS_TEST)))
F90OBJECTS_2OS_TEST := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_2OS_TEST)))

F90SOURCES_2OS_LPCSTEST := $(notdir $(filter %.f90, $(SOURCES_2OS_LPCSTEST)))
F90OBJECTS_2OS_LPCSTEST := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_2OS_LPCSTEST)))

F90SOURCES_2OS_OMPTEST := $(notdir $(filter %.f90, $(SOURCES_2OS_OMPTEST)))
F90OBJECTS_2OS_OMPTEST := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_2OS_OMPTEST)))

F90SOURCES_2OS_LPCSOMPTEST := $(notdir $(filter %.f90, $(SOURCES_2OS_LPCSOMPTEST)))
F90OBJECTS_2OS_LPCSOMPTEST := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_2OS_LPCSOMPTEST)))

F90SOURCES_UTIL := $(notdir $(filter %.f90, $(SOURCES_UTIL)))
F90OBJECTS_UTIL := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_UTIL)))

#
# Define desired targets
#

main  :  TwoOS_self_tester.exe    TwoOS_self_lpcs_tester.exe    
openmp:  TwoOS_self_OMPtester.exe TwoOS_self_lpcs_OMPtester.exe 

TwoOS_self_tester.exe: $(F90OBJECTS_2OS_TEST)
	$(FC) $^ -o $@
TwoOS_self_lpcs_tester.exe: $(F90OBJECTS_2OS_LPCSTEST)
	$(FC) $^ -o $@

TwoOS_self_OMPtester.exe: $(F90OBJECTS_2OS_OMPTEST)
	$(FC) $(FFLAGS_OPENMP) $^ -o $@
TwoOS_self_lpcs_OMPtester.exe: $(F90OBJECTS_2OS_LPCSOMPTEST)
	$(FC) $(FFLAGS_OPENMP) $^ -o $@
#TwoOS_self_OMPtester.exe: $(F90OBJECTS_2OS_OMPTEST)
#	$(FC) $^ -o $@
#TwoOS_self_lpcs_OMPtester.exe: $(F90OBJECTS_2OS_LPCSOMPTEST)
#	$(FC) $^ -o $@

TwoOS_diff: $(F90OBJECTS_UTIL)
	$(FC) $^ -o $@

.PHONY: clean
clean:
	rm -f *.o $(OBJ_FILES) *.mod $(MOD_FILES) *.log *.exe

