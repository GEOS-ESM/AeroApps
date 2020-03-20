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
BINS  = qfed_l2b.py qfed_l2c.py qfed_l3a.py qfed_l3b.py 
CSHS  = 

PY_MODS = qfed 
PY_EXTS = qfed/PlumeRise_.$(F2PYEXT) \
          qfed/LockPlume_.$(F2PYEXT) \
          qfed/crtmmodis_.$(F2PYEXT) 
#         qfed/tausfc_.$(F2PYEXT)

esma_install install: $(DEPS) $(BINS) $(CSHS) $(PY_EXTS) 
	$(MKDIR) $(ESMALIB) $(ESMALIB)/Python $(ESMAETC) \
                 $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -p $(BINS) $(CSHS)           $(ESMABIN)
	$(CP) -pr $(PY_MODS)               $(ESMALIB)/Python 
#	$(CP) -p *.a                       $(ESMALIB)
#	$(CP) -p *.rc                      $(ESMAETC)
#	$(CP) -p *.mod                     $(ESMAINC)/$(THIS)

#	gmake -f GNUmakefile cvstag_sub

esma_clean clean:
	$(RM) -rf *~ qfed/*~ *.[aox] *.[Mm][Oo][Dd] qfed/*.pyc qfed/*.$(F2PYEXT) 

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

SRCS  = PlumeRise_ConstMod.F90 PlumeRise_Mod.F90 LockPlume_Mod.F90 qsat_Mod.F90
CSRCS = 
OBJS := $(addsuffix .o, $(basename $(SRCS))) \
	$(addsuffix .o, $(basename $(CSRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS))) \
        $(addsuffix .d, $(basename $(BINS)))

FREAL = 

LIB_PYOBS = $(ESMALIB)/libGMAO_pyobs.a
INC_PYOBS = $(ESMAINC)/GMAO_pyobs

INC_DIRS = . $(INC_MPEU) $(INC_HERMES) $(INC_PYOBS) $(INC_CRTM)
USER_FFLAGS = $(BIG_ENDIAN) \
              $(foreach dir,$(INC_DIRS),$(M)$(dir)) 

vpath % $(INC_DIRS) ./tests

$(LIB) lib : $(DEPS) $(OBJS)
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)


%.x : $(LIB) %.F90
	$(FC) $(LDFLAGS) -o $@ $(D)TEST $*.F90 $(LIB) \
              $(LIB_EU) $(LIB_HERMES) $(LIB_GFIO) $(LIB_SDF)

ut_%.x : $(LIB) tests/ut_%.F90
	$(FC) $(LDFLAGS) -o $@ $(D)TEST tests/ut_$*.F90 $(LIB) \
              $(LIB_EU) $(LIB_HERMES) $(LIB_GFIO) $(LIB_SDF)

%.xxxx : $(LIB) %.F90
	$(FC) $(LDFLAGS) -o $@ $(D)TEST $*.F90 $(LIB) \
              $(LIB_EU) $(LIB_HERMES) $(LIB_GFIO) $(LIB_SDF)

# Python extension modules
# ------------------------

qfed/PlumeRise_.$(F2PYEXT): PlumeRise_py.F90 PlumeRise_Mod.o 
	$(F2PY) -c -m PlumeRise_ $(M). $(M)$(INC_PYOBS)\
                 PlumeRise_py.F90 PlumeRise_Mod.o $(LIB_PYOBS) \
		 only: plume plumevmd getvmd biome kpbl : 
	/bin/mv PlumeRise_.$(F2PYEXT) qfed/

qfed/LockPlume_.$(F2PYEXT): LockPlume_py.F90 LockPlume_Mod.o qsat_Mod.o 
	$(F2PY) -c -m LockPlume_ $(M). \
                 LockPlume_py.F90 LockPlume_Mod.o qsat_Mod.o 
	/bin/mv LockPlume_.$(F2PYEXT) qfed/

qfed/tausfc_.$(F2PYEXT): tausfc_py.F90
	$(F2PY) -c -m tausfc_ $(M). $(M)$(INC_RADTRANS)\
                 tausfc_py.F90 $(LIB_RADTRANS) \
	$(LIB_GLATOVS) $(LIB_FASTEM) \
	$(LIB_SARTA) $(LIB_OPTRAN) $(LIB_MIT) \
        $(LIB_IDL) $(LIB_HFFP) $(LIB_IRUTIL) \
        $(LIB_EU) $(LIB_HERMES) $(LIB_GFIO) $(LIB_SDF)
	/bin/mv tausfc_.$(F2PYEXT) qfed/

qfed/crtmmodis_.$(F2PYEXT): crtmmodis_py.F90 
	$(F2PY) -c -m crtmmodis_ $(M). $(M)$(INC_CRTM) \
                 crtmmodis_py.F90  $(LIB_CRTM)
	/bin/mv crtmmodis_.$(F2PYEXT) qfed/


# CVS tag name substitution 
# ------------------------
SEDCMD = s/@CVSTAG/$(CVSTAG)/g
SRCS_CVSTAG = emissions.py mxd14_l3.py PlumeRise.py
cvstag_sub := $(addprefix $(ESMALIB)/Python/qfed/, $(SRCS_CVSTAG))

CVSTAG := $(if $(wildcard $(ESMAETC)/CVSTAG),$(shell cat $(ESMAETC)/CVSTAG),$(shell echo 'CVSTAG_UNKNOWN'))

cvstag_sub:
	@for src in $($@); do  \
	   ( $(SED) -i "$(SEDCMD)" $$src ) \
	done

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.
