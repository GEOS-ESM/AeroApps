#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 3mar2004  Zaslavsky  Initial imlementation.
# 20Oct2004  da Silva  Standardization
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------


THIS := $(shell basename `pwd`)
LIB   = lib$(THIS).a

PY_MODS = 
PY_EXTS = VLIDORT_OMI_.$(F2PYEXT)  VLIDORT_BRDF_.$(F2PYEXT) 
PY_SCRP = 

SEDCMD = s^__ESMADIR__^\$(ESMADIR)^;s^__BASEDIR__^\$(BASEDIR)^

esma_install install python_install: $(DEPS) $(BINS) $(CSHS) $(PY_MODS) $(PY_EXTS) 
	$(MKDIR) $(ESMALIB) $(ESMALIB)/Python $(ESMAETC) \
                 $(ESMAINC)/$(THIS) $(ESMABIN)

esma_clean clean:
	$(RM) -rf *~ *.[aox] *.[Mm][Oo][Dd] *.pyc *.$(F2PYEXT) 

esma_distclean distclean: clean
	$(RM) -rf *.d __cache__

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

SRCS  = #
CSRCS = #
OBJS := $(addsuffix .o, $(basename $(SRCS))) \
	$(addsuffix .o, $(basename $(CSRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS))) \
        $(addsuffix .d, $(basename $(BINS)))

FREAL = 



vpath % $(INC_DIRS) 

THIS_GFIO = GMAO_gfio_r4

INC_VLIDORT90 = $(ESMADIR)/Linux/include/VLIDORT90
LIB_VLIDORT90 = $(ESMADIR)/Linux/lib/libVLIDORT90.a


INC_DIRS = . $(INC_MPEU) $(INC_CHEM_BASE) $(INC_VLIDORT90)
USER_FMODS = $(foreach dir,$(INC_DIRS),$(M)$(dir)) 
USER_FINCS = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 

LIBS = $(LIB_CHEM_BASE) $(LIB_MAPL_BASE) $(LIB_MAPL_BASE_STUBS) \
       $(LIB_MPEU) $(LIB_GFIO) $(LIB_ESMF) $(LIB_MPI) $(LIB_SDF) $(LIB_SYS) 

$(LIB) lib : $(DEPS) $(OBJS)
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)


%.x : $(LIB) %.F90
	$(FC) $(LDFLAGS) -o $@ $(D)TEST $*.F90 $(LIBS)

ut_%.x : $(LIB) tests/ut_%.F90
	$(FC) $(LDFLAGS) -o $@ $(D)TEST tests/ut_$*.F90 $(LIBS)

# Python extension modules
# ------------------------
python: $(PY_EXTS)


Mie_.$(F2PYEXT): Mie_py.F90 
	$(F2PY) -c -m Mie_ $(M). $(M)$(INC_CHEM_BASE) \
                 Mie_py.F90 $(LIBS) \
                 only: getdims getmiedims getmiegridded getaerobs\
                 getmieobs getmieobs2 scalar vector gettopobs :

MieObs_.$(F2PYEXT): MieObs_py.F90 
	$(F2PY) -c -m MieObs_ $(M). $(M)$(INC_CHEM_BASE) \
                 MieObs_py.F90 $(LIBS) 

OMI_.$(F2PYEXT): OMI_py.F90 
	$(F2PY) -c -m OMI_ $(M). $(M)$(INC_CHEM_BASE) \
	         $(I)$(INC_VLIDORT) \
                 OMI_py.F90 $(LIB_VLIDORT) $(LIBS) \
                 only: scalar vector :

VLIDORT_BRDF_.$(F2PYEXT): VLIDORT_BRDF_py.F90 
	$(FC) -c $(M). $(M)$(INC_CHEM_BASE) \
	         $(I)$(INC_VLIDORT90) \
                VLIDORT_BRDF_py.F90 $(LIB_VLIDORT90) $(LIBS) 

VLIDORT_OMI_.$(F2PYEXT): VLIDORT_OMI_py.F90 
	$(FC) -c $(M). $(M)$(INC_CHEM_BASE) \
	         $(I)$(INC_VLIDORT90) \
                 VLIDORT_OMI_py.F90 $(LIB_VLIDORT90) $(LIBS) \

VLIDORT_cloud_.$(F2PYEXT): VLIDORT_cloud_py.F90 
	$(F2PY) -c -m VLIDORT_cloud_ $(M). $(M)$(INC_CHEM_BASE) \
	         $(I)$(INC_VLIDORT90) \
                 VLIDORT_cloud_py.F90 $(LIB_VLIDORT90) $(LIBS) \


qsat_.$(F2PYEXT): qsat_py.F90 
	$(F2PY) -c -m qsat_ $(M). $(M)$(INC_MAPL_BASE) \
                 qsat_py.F90 $(LIBS) \

binObs_.$(F2PYEXT): binObs_py.F
	$(F2PY) -c -m binObs_ $(M). \
                 binObs_py.F

test_mieobs.$(F2PYEXT): test_mieobs.F90
	$(F2PY) -c -m test_mieobs $(M). \
                 test_mieobs.F90

%_.$(F2PYEXT): %_py.F90
	$(F2PY) -c -m $*_ $(M). $<

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.
