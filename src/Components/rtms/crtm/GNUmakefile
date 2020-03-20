#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 3mar2004  Zaslavsky  Initial imlementation.
# 20Oct2004  da Silva  Standardization
# 16Mar2007  Kokron    Disable default O3 optimization; add LOOP_VECT optimization
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../../..
endif
this_make := $(lastword $(MAKEFILE_LIST))

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

test:
	@ echo "$(this_make):$@: ARCH = $(ARCH)"
	@ echo "$(this_make):$@: MACH = $(MACH)"
	@ echo "$(this_make):$@: SITE = $(SITE)"
	@ echo "$(this_make):$@: NODE = $(NODE)"
	@ echo "$(this_make):$@:   FC = $(FC)"
	@ echo "$(this_make):$@: _make_this_lib_ = $(_make_this_lib_)"

THIS = $(shell basename `pwd`)
LIB  = lib$(THIS).a
.PHONY: lib

_endianness_ := Native
_endianness_ := Big_Endian
_endianness_ := Little_Endian

# Choose a make-entry for the given Fortran compiler (FC) on this OS, from the
# targets listed in this GNUmakefile.  This is a major change from previous OS
# based dicision making, to a compiler based dicision making.
#

ifeq ($(FC), )
  FC := $(error FC (Fortran Compiler) has not been configured in $(THIS), ARCH=$(ARCH))
endif

# Check FC to name the proper target.  Note: FC=ifort and FC=mpiifort are treated the same.
ifeq ($(findstring ifort,$(FC)),ifort)
  _make_this_lib_ := intel
else ifeq ($(findstring mpif90,$(FC)),mpif90)
  _make_this_lib_ := intel

else ifeq ($(findstring    xlf95,$(FC)),xlf95)
  _make_this_lib_ := ibm

else ifeq ($(findstring gfortran,$(FC)),gfortran)
  _make_this_lib_ := $(error target not configured in $(THIS), FC=$(FC), ARCH=$(ARCH))

else ifeq ($(findstring    pgf95,$(FC)),pgf95)
  _make_this_lib_ := pgi

else ifeq ($(findstring     lf95,$(FC)),lf95)
  _make_this_lib_ := $(error target not configured in $(THIS), FC=$(FC), ARCH=$(ARCH))

else ifeq ($(findstring      g95,$(FC)),g95)
  _make_this_lib_ := g95

else
#else ifeq ($(findstring      f95,$(FC)),f95)
#NAME_SUNOS = f95
#else ifeq ($(findstring      f90,$(FC)),f90)
#NAME_IRIX64 = f90
#NAME_HPUX = f90
#NAME_ABSOFT = f90
  _make_this_lib_ := $(error target not configured in $(THIS), FC=$(FC), ARCH=$(ARCH))
endif

ifeq ($(_endianness_),Big_Endian)
  _make_this_lib_ := $(_make_this_lib_)_big-endian
else ifeq ($(_endianness_),Little_Endian)
  _make_this_lib_ := $(_make_this_lib_)_little-endian
endif

ifeq ($(BOPT),g)
_make_this_lib_ := $(_make_this_lib_)_debug
endif

esma_install install: $(_make_this_lib_) $(BIN)
esma_install install:
	@ echo ">>> $@ \"$(LIB)\" etc. with \"$(my_build)\" build <<<"
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -p *.a         $(ESMALIB)
	$(CP) -p *.mod       $(ESMAINC)/$(THIS)

esma_clean clean:
	$(RM) *~ *.[aox] *.[Mm][Oo][Dd]

esma_distclean distclean:
	$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]

esma_doc doc:
	@echo "Target $@ not implemented yet in `pwd`"


esma_help help:
	@echo "Standard ESMA targets:"
	@echo "% make esma_install    (builds and install under ESMADIR)"
	@echo "% make esma_clean      (removes deliverables: *.[aox], etc)"
	@echo "% make esma_distclean  (leaves in the same state as cvs co)"
	@echo "% make esma_doc        (generates PDF, installs under ESMADIR)"
	@echo "% make esma_help       (this message)"
	@echo "Environment:"
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE)"
	@echo "        FREAL = $(FREAL)"

#                  --------------------
#                  User Defined Targets
#                  --------------------

include make.filelist
SRCS := $(SRC_FILES)

OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS)))

list:
	@ for d in $(DEPS); do echo $$d; done

ifeq ("$(FOPT)","-O3")
   FOPT += $(LOOP_VECT)
endif
FREAL = $(FREAL4)

MOD_DIRS = . $(INC_MPEU)
FFLAGS =
USER_FFLAGS = -IPF_fp_relaxed -ip $(BIG_ENDIAN) $(BYTERECLEN)
  
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 

vpath % $(MOD_DIRS)

lib: $(_make_this_lib_)

%.x : $(LIB) %.o 
	$(LD) $(LDFLAGS) -o $@ $*.o $(LIB) $(LIB_SYS)

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#==============================================================================
#
# Makefile for JCSDA CRTM Build
#
#==============================================================================

# -------------
# Define macros
# -------------
include ./make.macros
# -------------
# Override ESMA compilation rules
# -------------
include ./make.rules
.SUFFIXES: .d

# This makefile
MAKE_FILE = GNUmakefile

# The library name
LIBRARY_NAME = $(LIB)

# Below are OS based targets.  They are not used, but for information only.
sunos:
	$(MAKE) -f $(MAKE_FILE) library $(SUNOS_FLAGS)

aix:
	$(MAKE) -f $(MAKE_FILE) ibm

irix64:
	$(MAKE) -f $(MAKE_FILE) library $(IRIX64_FLAGS)

linux darwin:
	$(MAKE) -f $(MAKE_FILE) library $(LINUX_FLAGS)

# IBM AIX Compiler
# ----------------
ibm_debug:
	$(MAKE) -f $(MAKE_FILE) library $(AIX_FLAGS_DEBUG)

ibm:
	$(MAKE) -f $(MAKE_FILE) library $(AIX_FLAGS_PROD)


# HP-UX Compiler
# --------------

hp_debug:
	$(MAKE) -f $(MAKE_FILE) library $(HPUX_FLAGS_DEBUG)

hp:
	$(MAKE) -f $(MAKE_FILE) library $(HPUX_FLAGS_PROD)


# SGI MIPSpro Compiler
# --------------------
sgi_debug:
	$(MAKE) -f $(MAKE_FILE) library $(IRIX64_FLAGS_DEBUG)

sgi:
	$(MAKE) -f $(MAKE_FILE) library $(IRIX64_FLAGS_PROD)


# SunOS Compiler
# --------------
sun_debug:
	$(MAKE) -f $(MAKE_FILE) library $(SUNOS_FLAGS_DEBUG)

sun:
	$(MAKE) -f $(MAKE_FILE) library $(SUNOS_FLAGS_PROD)


# Linux compilers
# ---------------
# === Entries of using Intel fortran compiler, with "ifort ... convert=little_endian"
intel_little-endian_debug:
	$(MAKE) -f $(MAKE_FILE) library $(LINUX_FLAGS_INTEL_LittleEndian_DEBUG)

intel_little-endian:
	$(MAKE) -f $(MAKE_FILE) library $(LINUX_FLAGS_INTEL_LittleEndian_PROD)

# === Entries of using Intel fortran compiler, with "ifort ... convert=big_endian"
intel_big-endian_debug:
	$(MAKE) -f $(MAKE_FILE) library $(LINUX_FLAGS_INTEL_BigEndian_DEBUG)

intel_big-endian:
	$(MAKE) -f $(MAKE_FILE) library $(LINUX_FLAGS_INTEL_BigEndian_PROD)

# === Entries of using Intel fortran compiler, with "ifort ... convert=native"
intel_debug:
	$(MAKE) -f $(MAKE_FILE) library $(LINUX_FLAGS_INTEL_DEBUG)

intel:
	$(MAKE) -f $(MAKE_FILE) library $(LINUX_FLAGS_INTEL_PROD)

##
lahey_debug:
	$(MAKE) -f $(MAKE_FILE) library $(LINUX_FLAGS_LAHEY_DEBUG)

lahey:
	$(MAKE) -f $(MAKE_FILE) library $(LINUX_FLAGS_LAHEY_PROD)

##
pgi_debug:
	$(MAKE) -f $(MAKE_FILE) library $(LINUX_FLAGS_PGI_DEBUG)

pgi:
	$(MAKE) -f $(MAKE_FILE) library $(LINUX_FLAGS_PGI_PROD)

##
g95_debug:
	$(MAKE) -f $(MAKE_FILE) library $(LINUX_FLAGS_G95_DEBUG)

g95:
	$(MAKE) -f $(MAKE_FILE) library $(LINUX_FLAGS_G95_PROD)

##
gfortran_debug:
	$(MAKE) -f $(MAKE_FILE) library $(LINUX_FLAGS_GFORTRAN_DEBUG)

gfortran:
	$(MAKE) -f $(MAKE_FILE) library $(LINUX_FLAGS_GFORTRAN_PROD)

##
absoft_debug:
	$(MAKE) -f $(MAKE_FILE) library $(LINUX_FLAGS_ABSOFT_DEBUG)

absoft:
	$(MAKE) -f $(MAKE_FILE) library $(LINUX_FLAGS_ABSOFT_PROD)

#
# Build the library and include files.  The first target (library) is phony,
# while the second target ($(LIB)) is real.  The real target allows the checking
# of the target age against the ages of its prerequisits ($(OBJS)).
#

library: $(LIB)

$(LIB): $(OBJS)
	@ echo ">>> build \"$@\" with \"$(my_build)\" build <<<"
	$(REMOVE) $(LIBRARY_NAME)
	$(ARCHIVER) $(ARCHIVER_FLAGS) $(LIBRARY_NAME) $(OBJS)

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
