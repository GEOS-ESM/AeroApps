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
LIB   = # lib$(THIS).a

BINS = basic_1.x

esma_install install: $(DEPS) $(BINS)
	@echo "Nothing to install here. Running test..."
	./basic_1.x
	cat vlidort_v2p4rt_testing.mfout_1
#	$(MKDIR) $(ESMALIB) $(ESMALIB)/Python $(ESMAETC) \
#                 $(ESMAINC)/$(THIS) $(ESMABIN)
#	$(CP) -p $(BINS) $(CSHS)           $(ESMABIN)
#	$(CP) -pr $(PY_MODS) $(PY_EXTS)    $(ESMALIB)/Python 
#	$(CP) -p *.a                       $(ESMALIB)
#	$(CP) -p *.rc                      $(ESMAETC)
#	$(CP) -p *.mod                     $(ESMAINC)/$(THIS)

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

FREAL = $(FREAL4)

INC_DIRS = . $(INC_VLIDORT)
USER_FFLAGS = $(BIG_ENDIAN) \
              $(foreach dir,$(INC_DIRS),$(M)$(dir)) 

vpath % $(INC_DIRS) 

LIBS = $(LIB_VLIDORT)

%.x : $(LIB) %.f
	$(FC) $(LDFLAGS) -o $@ $*.f $(LIBS)

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.
