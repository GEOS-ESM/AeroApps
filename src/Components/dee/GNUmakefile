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

PY_EXTS = duem_.$(F2PYEXT)

esma_install install python_install: $(DEPS) $(PY_EXTS) 
	$(MKDIR) $(ESMALIB) $(ESMALIB)/Python $(ESMAETC) \
                 $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -pr $(PY_EXTS)               $(ESMALIB)/Python 

esma_clean clean:
	$(RM) -rf *~ *.[aox] *.[Mm][Oo][Dd] *.pyc *.$(F2PYEXT) 

esma_distclean distclean: clean
	$(RM) -rf *.d __cache__ *.$(F2PYEXT) *.pyc

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


# Python extension modules
# ------------------------
python: $(PY_EXTS)

duem_.$(F2PYEXT): duem_py.F90 
	$(F2PY) -c -m duem_  $(M). duem_py.F90

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
