#
# Makefile for ESMA components.
#
# REVISION HISTORY:
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

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------


THIS := $(shell basename `pwd`)
LIB   = lib$(THIS).a

PY_MODS = sixsv.py
PY_EXTS = sixsv_.$(F2PYEXT)
PY_SCRP = 

esma_install install python_install: $(DEPS) $(LIB) $(BINS) $(CSHS) $(PY_MODS) $(PY_EXTS) 
	$(MKDIR) $(ESMALIB) $(ESMALIB)/Python $(ESMAETC) \
                 $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -pr $(PY_MODS) $(PY_EXTS)    $(ESMALIB)/Python 
#	$(CP) -p $(PY_SCRP)                $(ESMABIN)

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

SRCS  = $(wildcard *.F)
CSRCS = #
OBJS := $(addsuffix .o, $(basename $(SRCS))) \
	$(addsuffix .o, $(basename $(CSRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS))) \
        $(addsuffix .d, $(basename $(BINS)))

FREAL = 

INC_DIRS = . $(INC_MPEU) $(INC_CHEM_BASE)
USER_FMODS = $(foreach dir,$(INC_DIRS),$(M)$(dir)) 
USER_FINCS = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 

vpath % $(INC_DIRS) 

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

sixsv_.$(F2PYEXT): sixsv_py.F90
	$(F2PY) -c -m sixsv_ $(M). sixsv_py.F90 $(LIB)

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.
