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
BINS := 

PY_EXTS = stnAngles_.$(F2PYEXT) VLIDORT_STN_.$(F2PYEXT)
PY_SCRP = g5nr_stn_sampler.py mcd12c_stn_sampler.py mcd43c_stn_sampler.py run_g5nr_stn_sampler_invariant.py \
			run_g5nr_stn_sampler.py ext_stn_vlidort_sampler.py
PY_LIBS = directsun_vlidort.py aeronet_vlidort.py stn_vlidort_aux.py pystn

esma_install install python_install: $(DEPS) $(LIB) $(BINS) $(PY_EXTS) 
	$(MKDIR) $(ESMALIB) $(ESMALIB)/Python $(ESMAETC) \
                 $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -pr $(PY_EXTS)            $(ESMALIB)/Python                 
	$(CP) -pr $(PY_LIBS)            $(ESMALIB)/Python
	$(CP) -p  $(BINS) $(PY_SCRP)    $(ESMABIN)

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
INC_VLIDORT2OS = $(ESMAINC)/VLIDORT2OS
LIB_VLIDORT2OS = $(ESMALIB)/libVLIDORT2OS.a

INC_GEOVLIDORT = ../geo_vlidort
LIB_GEOVLIDORT = ../geo_vlidort/libgeo_vlidort.a

INC_DIRS = . $(INC_MAPL_BASE) $(INC_ESMF) $(INC_MPI) $(INC_NETCDF) $(INC_CHEM_BASE) $(INC_VLIDORT2OS) $(INC_GEOVLIDORT) /usr/include
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


#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS =  #
OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS))) 
 
$(LIB) lib : $(DEPS) $(OBJS)
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

XLIBS =
ifeq ($(wildcard /etc/os-release/.*),)
   XLIBS = -L/usr/lib64 -lssl -lcrypto
endif

# Python extension modules
# ------------------------
python: $(PY_EXTS)

VLIDORT_STN_.$(F2PYEXT): VLIDORT_STN_py.F90 
	$(F2PY) -c -m VLIDORT_STN_ $(M). $(M)$(INC_GEOVLIDORT) \
                 VLIDORT_STN_py.F90 $(LIBS_F2PY) $(XLIBS)
# 	$(CP) -p VLIDORT_.$(F2PYEXT) $(ESMALIB)/Python 	

stnAngles_.$(F2PYEXT): stnAngles_py.F90 
	$(F2PY) -c -m stnAngles_ $(M). \
                 stnAngles_py.F90 $(LIBS_F2PY) $(XLIBS)

%.$(F2PYEXT): %py.F90 
	$(F2PY) -c -m $* $(M). $< 

#python python_install: $(PY_EXTS)

#stnAngles_.$(F2PYEXT): stnAngles_py.F90 
#	$(F2PY) -c -m stnAngles_  $(M). $(M)$(INC_MAPL_BASE) stnAngles_py.F90 $(LIB) $(LIBS_F2PY) $(XLIBS)

#VLIDORT_STN_.$(F2PYEXT): VLIDORT_STN_py.F90 
#	$(F2PY) -c -m VLIDORT_STN_  $(M). $(M)../geo_vlidort $(M)$(INC_MAPL_BASE) VLIDORT_STN_py.F90 $(LIB) $(LIBS_F2PY) $(XLIBS)


# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.
