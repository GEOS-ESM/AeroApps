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

PY_EXTS = pyobs/binObs_.$(F2PYEXT) pyobs/IGBP_.$(F2PYEXT) pyobs/sgp4_.$(F2PYEXT)
PY_SCRP = stn_sampler.py trj_sampler.py ext_sampler.py trj_xsect.py

esma_install install python_install: $(DEPS) $(LIB) $(PY_EXTS) 
	$(MKDIR) $(ESMALIB) $(ESMALIB)/Python $(ESMAETC) \
                 $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -pr pyobs             $(ESMALIB)/Python 
	$(CP) -pr $(PY_SCRP)        $(ESMABIN)
	$(RM) -rf $(ESMALIB)/Python/pyobs/CVS 
	$(CP) -p $(LIB)                    $(ESMALIB)
	$(CP) -p *.mod                     $(ESMAINC)/$(THIS)

esma_clean clean:
	$(RM) -rf *~ *.[aox] *.[Mm][Oo][Dd] pyobs/*~ pyobs/*.pyc *.$(F2PYEXT) pyobs/*.$(F2PYEXT) 

esma_distclean distclean: clean
	$(RM) -rf *.nc *.nc4 *.d __cache__

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

SRCS  = VegType_Mod.F90 glint_mod.F90  sgp4_mod.F90  TLE_mod.F90
CSRCS = VegType_io.c

OBJS := $(addsuffix .o, $(basename $(SRCS))) \
	$(addsuffix .o, $(basename $(CSRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS))) \
        $(addsuffix .d, $(basename $(BINS)))

FREAL = 

INC_DIRS = . $(INC_MPEU) $(INC_CHEM_BASE)
USER_FMODS = $(foreach dir,$(INC_DIRS),$(M)$(dir)) 
USER_FINCS = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 

vpath % $(INC_DIRS) 

THIS_GFIO = GMAO_gfio_r4

LIBS = $(LIB_CHEM_BASE) $(LIB_MAPL_BASE) $(LIB_EU) $(LIB_GFIO) \
       $(LIB_ESMF) $(LIB_MPI) $(LIB_SDF) $(LIB_SYS) 

$(LIB) lib : $(DEPS) $(OBJS)
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)


%.x : $(LIB) %.F90
	$(FC) $(LDFLAGS) -o $@ $(D)TEST $*.F90 $(LIBS)

# Python extension modules
# ------------------------
python: $(PY_EXTS)

pyobs/binObs_.$(F2PYEXT): binObs_py.F
	$(F2PY) -c -m binObs_ $(M). \
                 binObs_py.F
	$(CP) -p binObs_.$(F2PYEXT) pyobs

pyobs/IGBP_.$(F2PYEXT): IGBP_py.F90 VegType_Mod.o VegType_io.o
	$(F2PY) -c -m IGBP_ $(M). \
                 IGBP_py.F90 VegType_Mod.o VegType_io.o \
                 only: getsimpleveg getdetailedveg : 
	$(CP) -p IGBP_.$(F2PYEXT) pyobs/

pyobs/%.$(F2PYEXT): %py.F90 
	$(F2PY) -c -m $* $(M). $< $(LIB) 
	$(CP) -p $*.$(F2PYEXT) pyobs/

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
