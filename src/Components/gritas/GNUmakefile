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
BINS  = gritas.x grisas.x

INC_MPEU = $(ESMAINC)/GMAO_mpeu
INC_PSAS = $(ESMAINC)/psas
INC_HERMES = $(ESMAINC)/GMAO_hermes
INC_GFIO = $(ESMAINC)/GMAO_gfio

LIB_ODS  = $(ESMALIB)/libGMAO_ods.a
LIB_MPEU = $(ESMALIB)/libGMAO_mpeu.a
LIB_EU   = $(ESMALIB)/libGMAO_eu.a
LIB_PSAS = $(ESMALIB)/libpsas.a
LIB_HERMES = $(ESMALIB)/libGMAO_hermes.a
LIB_GFIO = $(ESMALIB)/libGMAO_gfio.a

esma_install install: $(DEPS) $(LIB) $(BINS)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -p *.a         $(ESMALIB)
	$(CP) -p *.rc        $(ESMAETC)
	$(CP) -p *.mod       $(ESMAINC)/$(THIS)
	$(CP) -p $(BINS)     $(ESMABIN)

esma_clean esma_distclean clean distclean:
	$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd] *.x

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

SRCS = etc.f gritas_gfio.f gritas_grads.f gritas_kxkt.f m_gritas_grids.f rdel.f

OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS))) \
        $(addsuffix .d, $(basename $(BINS)))

FOPT = $(FOPT3)
FPE  =

INC_DIRS = . $(INC_MPEU) $(INC_ODS) $(INC_PSAS) $(INC_HERMES)
USER_FFLAGS = $(BIG_ENDIAN) \
              $(foreach dir,$(INC_DIRS),$(M)$(dir)) 

vpath % $(INC_DIRS)

$(LIB) lib : $(DEPS) $(OBJS)
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

%.x : $(LIB) %.o
	$(FC) $(LDFLAGS) -o $@ $*.o $(LIB) \
                                     $(LIB_HERMES) \
              $(LIB_ODS) $(LIB_EU) $(LIB_GFIO) $(LIB_SDF)

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.
