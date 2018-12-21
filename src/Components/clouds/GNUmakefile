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

PY_EXTS = ICA_.$(F2PYEXT) mcsRegrid_.$(F2PYEXT) Aero_.$(F2PYEXT)

esma_install install python_install: $(DEPS) $(LIB) $(PY_EXTS) 
	$(MKDIR) $(ESMALIB) $(ESMALIB)/Python $(ESMAETC) \
                 $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -pr $(PY_EXTS)  *.py         $(ESMALIB)/Python 
	$(CP) -p $(LIB)                    $(ESMALIB)
	$(CP) -p *.mod                     $(ESMAINC)/$(THIS)

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

#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS = mcsRegrid_Mod.F90 \
       mod_mkl.F90 \
       mod_StaticGridcolumnStatisticalModel.F90 \
       mod_GaussianCopula.F90 \
       mod_consts.F90 \
       mod_eswi.F90 \
       mod_icefrac.F90 \
       mod_reff.F90 \
       mod_simTauRe.F90 \
       mod_triangle.F90 \
       mod_utils.F90

CSRCS = 
OBJS := $(addsuffix .o, $(basename $(SRCS))) \
	$(addsuffix .o, $(basename $(CSRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS))) \
        $(addsuffix .d, $(basename $(BINS)))

FREAL = 

INC_DIRS = . $(INC_MPEU) $(INC_CHEM_BASE) $(INC_MAPL_BASE) $(FPATH)
USER_FMODS = $(foreach dir,$(INC_DIRS),$(M)$(dir)) 
USER_FINCS = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
 
LIB_SCI_PY = -L$(MKLPATH) -lmkl_rt

vpath % $(INC_DIRS) 

THIS_GFIO = GMAO_gfio_r4

LIBS = $(LIB_CHEM_BASE) $(LIB_MAPL_BASE) $(LIB_EU) $(LIB_GFIO) \
       $(LIB_ESMF) $(LIB_MPI) $(LIB_SDF) $(LIB_SYS) 

# f2py cannot handle a -Wl flag used by Open MPI. But that option is not
# necessary for the code to compile or run. So if it is part of LIB_MPI
# we just strip it from LIBS
BADPATTERN := -Wl,--export-dynamic
LIBS_F2PY1 = $(subst $(BADPATTERN),,$(LIBS))
BADPATTERN2 := -Wl,-commons,use_dylibs
LIBS_F2PY2 = $(subst $(BADPATTERN2),,$(LIBS_F2PY1))
BADPATTERN3 := -pthread
LIBS_F2PY3 = $(subst $(BADPATTERN3),,$(LIBS_F2PY2))
BADPATTERN4 := -Wl,-flat_namespace
LIBS_F2PY = $(subst $(BADPATTERN4),,$(LIBS_F2PY3))

$(LIB) lib : $(DEPS) $(OBJS)
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)


%.x : $(LIB) %.F90
	$(FC) $(LDFLAGS) -o $@ $(D)TEST $*.F90 $(LIBS)

# Python extension modules
# ------------------------
python: $(PY_EXTS)

Aero_.$(F2PYEXT): Aero_py.F90 
	$(F2PY) -c -m Aero_  $(M). $(M)$(INC_MAPL_BASE) Aero_py.F90 $(LIB) $(LIBS_F2PY) $(LIB_SCI_PY)

ICA_.$(F2PYEXT): ICA_py.F90 
	$(F2PY) -c -m ICA_  $(M). $(M)$(INC_MAPL_BASE) ICA_py.F90 $(LIB) $(LIBS_F2PY) $(LIB_SCI_PY)

mcsRegrid_.$(F2PYEXT): mcsRegrid_py.F90 
	$(F2PY) -c -m mcsRegrid_  $(M). mcsRegrid_py.F90 $(LIB) $(LIB_SCI_PY)

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
