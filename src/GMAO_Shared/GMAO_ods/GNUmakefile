#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 19Jan2004  Zaslavsky  Initial imlementation.
# 20Oct2004  da Silva   Standardization
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

THIS  := $(shell basename `pwd`)
FREAL = $(FREAL4)

LIB  = lib$(THIS).a

BIN_DOTX_SRC = dconv2ods.f  ods_blank.f90  ods_maker.f90  ods_sample.f  \
               ods_scan.f   ods_tally.f90  saber2ods.f ods_synX2Y.f90
BIN_NODX_SRC = odslist.f  odsmatch.f  odsselect.f  odsshuffle.f odsstats.f

BIN_DOTX = $(addsuffix .x, $(basename $(BIN_DOTX_SRC)))
BIN_NODX = $(basename $(BIN_NODX_SRC))
BINS     = $(BIN_DOTX) $(BIN_NODX)

esma_install install: $(LIB) $(BINS) python_install
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -p *.a    $(ESMALIB)
	$(CP) -p *.[Mm][Oo][Dd] $(ESMAINC)/$(THIS)
	$(CP) -p $(BINS)        $(ESMABIN)
	$(CP) -p unique_instruments.csh $(ESMABIN)
	$(CP) -p odslist.odt $(ESMAETC)
	$(CP) -p *.rc        $(ESMAETC)
	$(CP) -p obimp_summary.pl $(ESMABIN)
	$(SED) -e "s^@DASPERL^$(PERL)^" < diag2ods > $(ESMABIN)/diag2ods
	chmod 755 $(ESMABIN)/diag2ods
	chmod 755 $(ESMABIN)/obimp_summary.pl
	chmod 755 $(ESMABIN)/unique_instruments.csh

esma_clean clean: python_clean
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd] $(BINS)

esma_distclean distclean: python_clean
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd] $(BINS)

esma_doc doc:
	@echo "Target $@ not implemented yet in `pwd`"

esma_help:
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
	@echo "         OBJS= $(OBJS) "

help:
	@ echo "Available targets:"
	@ echo "  make lib             (to create libods.a)"
	@ echo "  make lib MODOBJS=    (to omit f90 modules)"
	@ echo "  make test            (to test ODS library)"
	@ echo "  make doc             (to create documentation)"
	@ echo "  make all"  
	@ echo "  make export"  
	@ echo "  make clean"  
	@ echo "  make distclean"
	@ echo "  make depend"  


INC_PILGRIM = $(ESMAINC)/pilgrim
LIB_PILGRIM = $(ESMALIB)/libpilgrim.a

MODSRC = m_ods.f90 m_roms.f90 m_odsmeta.f m_ods_structure.f m_ods_obsdiags.F90 \
         m_saber.f m_odsxsup.f90 m_Sndx.F90 m_ODS_Sample.f m_mraob.f90 

ifneq ( $(wildcard $(LIB_PILGRIM)), $(null) ) # for MP ods
   MODSRC += m_MPodsmeta.F m_MPodsdata.F90  m_odstransfer.f90 m_MPods.F90  
endif

ODSINC = ods_hdf.h  ods_worksp.h  ods_stdio.h

ODSSRC = \
ods_cal2time.f  ods_cal2min.f\
ods_julian.f    ods_ncp_itoi4.f  ods_rget.f\
ods_caldat.f    ods_namemonth.f  ods_ncp_itor4.f  ods_scalei.f\
ods_case.f      ods_ncagtc.f     ods_ncp_rtoi1.f  ods_scaleirev.f\
ods_cget.f      ods_ncagti.f     ods_ncp_rtoi2.f  ods_scaler.f\
ods_checki.f    ods_ncagtr.f     ods_ncp_rtoi4.f  ods_scalerrev.f\
ods_checkr.f    ods_ncapti.f     ods_ncp_rtor4.f  ods_setparmc.f\
ods_coards.f    ods_ncaptr.f     ods_ncvgti.f     ods_setparmi.f\
ods_defwsp.f    ods_nccreate.f   ods_ncvgtr.f     ods_setparmr.f\
ods_file.f      ods_ncg_i1toi.f  ods_ncvpti.f     ods_strsearch.f\
ods_getattp.f   ods_ncg_i1tor.f  ods_ncvptr.f     ods_strsize.f\
ods_getc.f      ods_putc.f\
ods_ncvptc.f    ods_ncvgtc.f\
ods_geti.f      ods_ncg_i2toi.f  ods_ncvwsp.f     ods_time2cal.f\
ods_min2cal.f   ods_mraob.f\
ods_getlist.f   ods_ncg_i2tor.f  ods_nget.f       ods_type.f\
ods_getp.f      ods_ncg_i4toi.f  ods_nval.f       ods_updatep.f\
ods_getr.f      ods_ncg_i4tor.f  ods_parm.f       ods_valwsp.f\
ods_handle.f    ods_ncg_r4toi.f  ods_puti.f       ods_vartype.f\
ods_iget.f      ods_ncg_r4tor.f  ods_putlist.f    ods_verindex.f\
ods_imonth.f    ods_ncopen.f     ods_putr.f       ods_vertag.f\
ods_ncp_itoi1.f ods_readp.f      ods_writep.f\
ods_julhr.f     ods_ncp_itoi2.f  ods_resetp.f     odsnxtime.f\
ods_eget.f      ods_dcget.f      ods_dconv.f      ods_dcscan.f\
ods_drad.f      ods_dsbuv.f      ods_dpcp.f       ods_rmdupl.f\
ods_do3lev.f


DOCSRC = $(MODSRC) \
         ods_file.f   ods_puti.f   ods_putr.f   ods_putc.f \
                      ods_geti.f   ods_getr.f   ods_getc.f \
                      ods_iget.f   ods_rget.f   ods_cget.f \
                      ods_eget.f   ods_nget.f \
         ods_getlist.f \
         ods_julian.f ods_caldat.f ods_time2cal.f ods_cal2time.f \
         ods_min2cal.f \
         ods_cal2min.f \
         ods_setparmi.f \
         ods_setparmr.f \
         ods_setparmc.f \
         ods_parm.f

#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS = $(MODSRC) $(ODSSRC)
OBJS := $(addsuffix .o, $(basename $(SRCS)))

ALL_SRCS = $(SRCS) $(BIN_DOTX_SRC) $(BIN_NODX_SRC)
DEPS := $(addsuffix .d, $(basename $(ALL_SRCS)))

FOPT = $(FOPT3)
NCGEN = $(BASEBIN)/ncgen
INC_MPEU = $(ESMAINC)/GMAO_mpeu
LIB_MPEU = $(ESMALIB)/libGMAO_eu.a
USER_FFLAGS = $(BIG_ENDIAN) $(MPFLAG) $(M). $(I)$(INC_SDF) $(M)$(INC_MPEU) \
              $(M)$(INC_PILGRIM)

MOD_DIRS = . $(INC_SDF) $(INC_MPEU)

vpath % $(MOD_DIRS)

$(LIB) lib : $(DEPS) $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

%.x : $(LIB) $(LIB_MPEU) %.o
	$(LD) $(LDFLAGS) -o $@ $*.o $(LIB) $(LIB_MPEU) $(LIB_SDF) $(LIB_SYS)

odslist odsselect odsshuffle odsmatch odsstats: $(LIB) $(LIB_MPEU) \
                              odslist.o odsselect.o odsshuffle.o odsmatch.o odsstats.o
	$(LD) $(LDFLAGS) -o $@ $@.o $(LIB) $(LIB_MPEU) $(LIB_SDF) $(LIB_SCI) $(LIB_SYS)

sample1.nc : sample1.cdl
	$(NCGEN) -o sample1.nc sample1.cdl

sample2.nc : sample2.cdl
	$(NCGEN) -o sample2.nc sample2.cdl

test: ods_test.x sample1.nc sample2.nc
	./ods_test.x > ods_test.out
	diff ods_test.out ods_test.dao
	@ echo '        ***** ODS test successful ****'


#                  --------------------
#                        Python
#                  --------------------

python : pyods/odsmeta.py pyods/pyods_.$(F2PYEXT)
	@echo "Python modules created; 'make python_install' to install them".

pyods/odsmeta.py : odsmeta.h
	@grep 'integer, parameter' odsmeta.h | \
        sed -e 's/integer, parameter//g' -e 's/:://g' -e 's/\!/   #/g' -e 's/^[ ]*//g' | \
        grep -v idsats \
        > pyods/odsmeta.py

pyods/pyods_.$(F2PYEXT): pyods_.F90 $(LIB) 
	$(F2PY) -c -m pyods_ $(M). pyods_.F90 $(LIB) $(LIB_MPEU) $(LIB_SDF) $(LIB_SYS)
	mv pyods_.$(F2PYEXT) pyods

python_clean:
	$(RM) -rf dist build pyods/pyods_.so pyods/odsmeta.py pyods/*.pyc pyods/*~

python_install: pyods/pyods_.$(F2PYEXT) pyods/odsmeta.py
	$(MKDIR) $(ESMALIB)/Python
	$(CP) -pr pyods $(ESMALIB)/Python

python_bdist:
	env BASELIB=$(BASELIB) ESMALIB=$(ESMALIB) ODSLIB=$(PWD) python setup.py bdist

python_system_install: pyods/pyods_.$(F2PYEXT)
	sudo env BASELIB=$(BASELIB) ESMALIB=$(ESMALIB) ODSLIB=$(PWD) python setup.py install

#                  --------------------
#                      Dependencies
#                  --------------------

#include Makefile.depend

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros
#.
