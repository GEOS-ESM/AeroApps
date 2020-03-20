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
LIB  := lib$(THIS).a

CORE = create_startf.x caltra.x trace.x select.x density.x lidar.x

GOODIES = traj2num.x lsl2rdf.x changet.x extract.x getmima.x gettidiff.x getvars.x \
          list2lsl.x lsl2list.x mergetra.x newtime.x reformat.x timeres.x trainfo.x \
          difference.x datelist.x tracal.x

BINS = $(CORE) $(GOODIES)

SCRIPTS = $(foreach script, \
                    caltra changet create_startf create_startf.pl \
                    lagrantohelp trainfo, \
                    scripts/$(script) )

esma_install install python_install: $(DEPS) $(BINS) $(SCRIPTS)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -p $(LIB)                   $(ESMALIB)
	$(CP) -p $(BINS) $(SCRIPTS)       $(ESMABIN)

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

SRC_DIRS = legacy/lib
BIN_DIRS = $(foreach dir,$(basename $(CORE)), legacy/$(dir)) legacy/goodies

SRCS := $(foreach dir,$(SRC_DIRS), \
        $(wildcard $(dir)/*.[fFc]) $(wildcard $(dir)/*.[fF]90) ) \
         legacy/trace/calvar.f legacy/select/special.f \
         legacy/goodies/evalmod.f90 legacy/goodies/precmod.f90 legacy/goodies/stringmod.f90
OBJS := $(notdir $(addsuffix .o, $(basename $(SRCS))) )
DEPS := $(notdir $(addsuffix .d, $(basename $(SRCS))) )

FREAL = $(FREAL4)

INC_DIRS = . $(INC_NETCDF)
USER_FMODS = $(foreach dir,$(INC_DIRS),$(M)$(dir)) 
USER_FINCS = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 

vpath % $(INC_DIRS) $(SRC_DIRS) $(BIN_DIRS)

THIS_GFIO = GMAO_gfio_r4

$(LIB) lib : $(OBJS)
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(LIB)

%.x : $(LIB) %.o
	$(FC) $(LDFLAGS) -o $@  $*.o $(LIB) $(LIB_NETCDF) $(LIB_MPI) $(LIB_SYS)

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.
