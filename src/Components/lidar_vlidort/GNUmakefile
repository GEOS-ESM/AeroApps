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


PY_EXTS = LidarAngles_.$(F2PYEXT) VLIDORT_.$(F2PYEXT)
PY_SCRP =  

esma_install install python_install: $(DEPS) $(LIB) $(BINS) $(PY_EXTS) 
	$(MKDIR) $(ESMALIB) $(ESMALIB)/Python $(ESMAETC) \
                 $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -pr $(PY_EXTS)               $(ESMALIB)/Python 
#	$(CP) -p  $(BINS) $(PY_SCRP)       $(ESMABIN)

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
INC_VLIDORT90  = $(ESMAINC)/VLIDORT90
LIB_VLIDORT90  = $(ESMALIB)/libVLIDORT90.a

INC_VLIDORT2OS = $(ESMAINC)/VLIDORT2OS
LIB_VLIDORT2OS = $(ESMALIB)/libVLIDORT2OS.a

INC_LIDORT3p7  = $(ESMAINC)/LIDORT3p7
LIB_LIDORT3p7  = $(ESMALIB)/libLIDORT3p7.a

INC_GEOVLIDORT = ../geo_vlidort
LIB_GEOVLIDORT = ../geo_vlidort/libgeo_vlidort.a

INC_DIRS = . $(INC_MAPL_BASE) $(INC_ESMF) $(INC_MPI) $(INC_NETCDF) $(INC_CHEM_BASE) $(INC_VLIDORT2OS) $(INC_LIDORT3p7) $(INC_GEOVLIDORT)) /usr/include
USER_FMODS = $(foreach dir,$(INC_DIRS),$(M)$(dir)) 
USER_FINCS = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 

vpath % $(INC_DIRS) 

FREAL = $(FREAL4)
THIS_CFIO = MAPL_cfio_r4
THIS_GFIO = GMAO_gfio_r4

LIBS = $(LIB_GEOVLIDORT) $(LIB_VLIDORT2OS) $(LIB_LIDORT3p7) $(LIB_CHEM_BASE) \
       $(LIB_MAPL_BASE) $(LIB_MAPL_BASE_STUBS) \
       $(LIB_MAPL_BASE_STUBS2) \
       $(LIB_MPEU) $(LIB_CFIO) $(LIB_GFIO) $(LIB_NETCDF) \
       $(LIB_ESMF) $(LIB_SDF) $(LIB_MPI) $(LIB_SYS)

#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS =  #
OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS))) 
 
$(LIB) lib : $(DEPS) $(OBJS)
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

%.x : $(LIB) %.F90
	$(FC) $(LDFLAGS) -o $@ $(D)TEST $*.F90 $(LIBS)	


# Python extension modules
# ------------------------
python: $(PY_EXTS)

VLIDORT_.$(F2PYEXT): VLIDORT_py.F90 
	$(F2PY) -c -m VLIDORT_ $(M). $(M)$(INC_GEOVLIDORT) \
                 VLIDORT_py.F90 $(LIBS)                  
# 	$(CP) -p VLIDORT_.$(F2PYEXT) $(ESMALIB)/Python 	

%.$(F2PYEXT): %py.F90 
	$(F2PY) -c -m $* $(M). $< 


# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.
