
VLID_DEF_PATH   = ../../vlidort_def
VLID_TEST_PATH  = .

MOD_PATH = mod
OBJ_PATH = obj

MOD_FILES = $(MOD_PATH)/*.mod
OBJ_FILES = $(OBJ_PATH)/*.o

# Define default shell make will use

SHELL = /bin/bash

# Define FORTRAN90 compiler to use (can be defined here, but usually defined on makefile command line)
#  gfortran
FC = gfortran
FFLAGS = -c

# Additional flags for gfortran
ifeq ($(FC), gfortran)
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -J$(MOD_PATH)
#	FFLAGS_DEBUG = -g -Wall -fbounds-check
#	FFLAGS_DEBUG = -g -Wall -fbounds-check -fbacktrace
	FFLAGS_DEBUG = -g -C -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
	FFLAGS_OPT = -O3
        FFLAGS_OPENMP = -fopenmp -frecursive
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


# Define list of source files
# (Note: ordering is important because of dependencies)


BASE_SOURCES =
BASE_SOURCES +=   \
   $(VLID_DEF_PATH)/vlidort_pars.f90

SOURCES =
SOURCES +=   \
   $(BASE_SOURCES) \
   $(VLID_TEST_PATH)/vsleave_findpar.f90         \
   $(VLID_TEST_PATH)/vsleave_sup_inputs_def.f90  \
   $(VLID_TEST_PATH)/vsleave_sup_outputs_def.f90 \
   $(VLID_TEST_PATH)/vsleave_sup_aux.f90         \
   $(VLID_TEST_PATH)/vsleave_sup_routines.f90    \
   $(VLID_TEST_PATH)/vsleave_sup_masters.f90     \
   $(VLID_TEST_PATH)/vsleave_sup_mod.f90

SOURCES_SELF = $(SOURCES) + \
   $(VLID_TEST_PATH)/vsleave_self_tester_1.f90

# Define pattern rules for creating object files:


# For vlidort main source files
$(OBJ_PATH)/%.o : $(VLID_DEF_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_TEST_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@

# Define object files

F90SOURCES_SLEAVE_SELF := $(notdir $(filter %.f90, $(SOURCES_SELF)))
F90OBJECTS_SLEAVE_SELF := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_SLEAVE_SELF)))

# Define desired targets

All: vsleave_self_tester_1.exe
vsleave_self_tester_1.exe: $(F90OBJECTS_SLEAVE_SELF)
	$(FC) $^ -o $@

.PHONY: clean
clean:
	rm -f *.o $(OBJ_FILES) *.mod $(MOD_FILES) *.log *.exe

