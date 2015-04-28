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
PY_EXTS = GeoAngles_.$(F2PYEXT)
PY_SCRP = geo_warp.py geo_sampler.py

esma_install install python_install: $(DEPS) $(BINS) $(CSHS) $(PY_MODS) $(PY_EXTS) 
	$(MKDIR) $(ESMALIB) $(ESMALIB)/Python $(ESMAETC) \
                 $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -pr $(PY_MODS) $(PY_EXTS)    $(ESMALIB)/Python 
	$(CP) -p  $(PY_SCRP)               $(ESMABIN)

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

INC_DIRS = . $(INC_MAPL_BASE) $(INC_MPI) $(INC_NETCDF) $(INC_CHEM_BASE) $(INC_MPEU) $(INC_VLIDORT90)
USER_FMODS = $(foreach dir,$(INC_DIRS),$(M)$(dir)) 
USER_FINCS = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 

vpath % $(INC_DIRS) 

THIS_GFIO = GMAO_gfio_r4

# Python extension modules
# ------------------------
python: $(PY_EXTS)

%.$(F2PYEXT): %py.F90 
	$(F2PY) -c -m $* $(M). $< 

%.x: %.o sys_tracker.o
	$(CC) -c sys_tracker.c
	$(FC) $(FOPTS) $^ -o $@  $(LIB_CHEM_BASE) $(LIB_MAPL_BASE) $(LIB_MAPL_BASE_STUBS) \
	$(LIB_MPEU) $(LIB_GFIO) $(LIB_ESMF) $(LIB_MPI) $(LIB_SDF) $(LIB_SYS) \
	$(LIB_NETCDF) 


# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.
