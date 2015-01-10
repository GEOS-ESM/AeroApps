#!/bin/make
#-----------------------------------------------------------------------
#       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
#-----------------------------------------------------------------------
# $Header: /cvsroot/esma/esma/src/Shared/GMAO_psas/Makefile,v 1.4 2011-02-26 05:13:52 adasilva Exp $
#
# !ROUTINE: Makefile
#
# !DESCRIPTION:
#	A makefile for PSAS
#
# !CALLING SEQUENCE: make
#
# !BUGS:
#
# !SEE ALSO: Makefile.conf Makefile.IRIX Makefile.OSF1 Makefile.UNICOS
#
# !REVISION HISTORY: (`cvs log <>' for more)
#
# 	05-Aug-97 - Jing Guo - Added the switch to reproduce etest14
#
# 	01Sep95 - J. Guo	- added the prolog
#				- Rename from Makefile
#	30Dec94 - Jing G.	-
#			- With new correlation functions calls for sea-
#			  level variables
#			- With parallelization codes in some of files
#			- New PSAS_OBJS and PSAS_INCL listings
# file: Makefile: last change: 02Nov94 (J. Guo)
# file: Makefile: last change: 03oct94 (A. da S.)
#
#.......................................................................

# Possible values of DEFS include:
#	TRACE_CALLS

SHELL		=sh

# Macro _E14 is a switch to ensure the code reproduce a zero-deff
# result as the version of psas-20nov96, when an "etest14" resource
# file is used.  The only difference is in the source file hcorfuns.F90,
# where a different expresion is used to convert an arc-distance to a
# cordal distance in subroutine gaspari-cohn().  Check the source file
# for detail.  It can be down through the command,
#
#	make libpsas DEFS=-D_ETEST14_
#

DEFS		= -DTRACE_CALLS
DEFS		=

MAKE		= make
MAKEFILE	= Makefile
MAKECONF	= Makefile.conf

CPP		= /lib/cpp
CPPFLAGS	= -P

F77		= f77
F77FLAGS	= $(FFLAGS)
F77PP		= $(FPP)
F77PPFLAGS	= $(FPPFLAGS)

F90		= f90
F90FLAGS	= $(FFLAGS)
F90PP		= $(FPP)
F90PPFLAGS	= $(FPPFLAGS)

FC		= $(F90)
FFLAGS		= $(F90FLAGS)
FPP		= $(CPP)
FPPFLAGS	= $(CPPFLAGS)

LD		= $(F90)
LDFLAGS		=

LN		= ln -s
RM		= rm -f
CP		= cp -p

ZIP		= compress
ZIPTAG		= .Z
ZCAT		= zcat
TAR		= tar
AR		= ar cq

include	$(MAKECONF)
#sinclude $(MAKECONF)

.SUFFIXES:
.SUFFIXES: .H .h .F90 .f90 .F77 .f77 .F .f .o

.f77.o:
	@$(LN) $*.f77 $*.f
	$(F77) -c $(F77FLAGS) $*.f
	@$(RM) $*.f

.F77.o:
	@$(RM) $*.f
	$(F77PP) $(F77PPFLAGS) $(DEFS) $*.F77 $*.f
	$(F77) -c $(F77FLAGS) $*.f
	@$(RM) $*.f

.f.o:
	$(FC) -c $(FFLAGS) $*.f

.F.o:
	@$(RM) $*.f
	$(FPP) $(FPPFLAGS) $(DEFS) $*.F $*.f
	$(FC) -c $(FFLAGS) $*.f
	@$(RM) $*.f

.f90.o:
	@if [ -f $*.[Ff] ]; then \
	   echo "make: WARNING! coexistent $*.[Ff] and $*.f90" ; fi
	$(F90) -c $(F90FLAGS) $*.f90

.F90.o:
	@ if [ -f $*.[Ff] ]; then \
	   echo "make: WARNING! coexistent $*.[Ff] and $*.F90" ; fi
	@ $(RM) $*.f90
	$(F90PP) $(F90PPFLAGS) $(DEFS) $*.F90 $*.f90
	$(F90) -c $(F90FLAGS) $*.f90
	@ $(RM) $*.f90

.H.h:
	$(CPP) $(CPPFLAGS) $(DEFS) -Dsys`uname -s` $*.H $*.h

help:
	@ echo "make psas.x to create the program excutable"

conf:	$(MAKECONF)

$(MAKECONF):
	@ CONF=$@.`uname -n`; \
	  if [ ! -r $${CONF} ]; then \
	    CONF=$@.`uname -s`; fi; \
	  $(RM) $${CONF}
	  $(LN) $${CONF} $@

#=======================================================================
# Targets specific macros

LIBpsas		= libpsas.a

LIBpsas_MODS	= \
 FEalpha_imat.o FEalpha_tabl.o     config.o      const.o     cordriv.o\
    cormats.o  hfecH_tbl.o  hfecQ_tbl.o  hoecH_tbl.o OEclass_tbl.o\
  rlat_imat.o  vfecH_tbl.o  vfecQ_tbl.o  voecH_tbl.o      Cprodm.o\
 FEsigW_tabl.o FEsigW_imat.o   hcorfuns.o  m_inpak90.o  m_chars.o \
 m_grads.o

LIBpsas_OBJS	= \
    aj_Alf.o     aj_Bet.o     aj_Gam.o     bands0.o     cgnorm.o\
    conjgr.o     Cprod1.o     Cprodx.o   diagcorD.o   diagcorF.o\
  diagcorM.o   diagcorO.o   diagcorU.o    dupelim.o    fDDcor1.o\
   fDDcorx.o    fHDcorx.o    fHHcor1.o    fHHcorx.o    fQQcor1.o\
   fQQcorx.o    gCprodx.o     gdstat.o     Gea2ll.o   getAIall.o\
   getAinc.o     gridxx.o       hunt.o imat_alpha.o\
   index3r.o    indexxi.o    indexxr.o    indexxs.o    iniainc.o\
    initkl.o   initRSRC.o              intp_ctaus.o  intp_hCor.o\
 intp_sigO.o  intp_vCor.o   inxSlist.o      ipick.o    krname0.o\
   ktname0.o    kxname0.o   listvals.o    ll2qvec.o     ll2xyz.o\
   llboxes.o     lnblnk.o     lstins.o    luavail.o     lvstat.o\
 merg_lats.o merg_plevs.o      mvPHx.o     mv_Alf.o     mv_Bet.o\
    mv_Gam.o    obssmry.o     obstat.o   offdcorF.o    oHHcor1.o\
   opntext.o      op_Mx.o    pardisp.o    permuti.o\
   permutl.o    permutr.o     proxel.o    prtdata.o    prtregn.o\
  psasexit.o   psasrcbd.o      qtrig.o     qtrig0.o     qvmvdv.o\
    rdpars.o    rdkttbl.o    rdkxtbl.o   rdlevels.o     \
   rdoetbl.o    rdvctbl.o   rec_Cxpy.o    regsort.o   restrict.o\
    sepang.o     setbox.o    seticos.o   setpindx.o     setpix.o\
 set_fecHH.o  set_fecQQ.o set_FEhCor.o set_FEvCor.o  set_oecHH.o\
set_OEclas.o set_OEhCor.o set_OEvCor.o    slintab.o    slogtab.o\
      smex.o     smexsh.o    solve4x.o       sort.o     sparse.o\
    spline.o   sym_Cxpy.o                tabIlist.o   tabRlist.o\
  tabSlist.o    tofront.o    typsort.o    uHHcor1.o     vindex.o\
   xyz2reg.o      zeits.o    opnieee.o     amatch.o   offdcorO.o\
   oHHcorx.o	  op_Pf.o tabl_FEalpha.o  intp_sigF.o   savevect.o\
  savemat1.o   savematx.o  indexGvec.o     \
  dervsigF_slW.o dervsigF_slD.o \
  dervsigF_upW.o dervsigF_upD.o \
   getsigF.o     dervsigF_slp.o \
  mark_nsig.o   getivec.o mv_diag.o imat_sigW.o tabl_FEsigW.o\
  getAImix.o	getAIzuv.o  getAIpuv.o   mark_reg.o   makedelv.o\
  stdv_FEqc.o  stdv_OEqc.o

LIBpsas_INCL	= \
     bands.h     gridxx.h     hfecHH.h     hfecQQ.h             \
    hoecHH.h    icosdef.h  kind_covs.h  kind_mats.h     krname.h\
     ktmax.h     kttabl.h   ktwanted.h      kxmax.h     kxtabl.h\
    latmax.h    levtabl.h      lvmax.h     maxreg.h    maxrlen.h\
  MXpar_hc.h     mxpass.h     mxvecs.h   MX_hoecH.h   MX_voecH.h\
   nobsmax.h     proxel.h     psasrc.h        qea.h       rade.h\
  realvals.h      stdio.H    trigtab.h     vfecHH.h     vfecQQ.h\
    voecHH.h   pres4slp.h   MX_hfecH.h

#.......................................................................

$(LIBpsas): $(LIBpsas_DEPS)
	$(RM) $(LIBpsas)
	@ echo $(AR) $(LIBpsas) $(LIBpsas_DEPS)
	@ $(AR) $(LIBpsas) $(LIBpsas_DEPS)

#.......................................................................

PSAS_OBJS	= \
 aio_grads.o    getdel2.o  psasshell.o wrtxvec.o

PSAS_INCL	= \
       aio.h      bands.h     gridxx.h      ktmax.h\
  ktwanted.h      kxmax.h     maxreg.h    nobsmax.h     psasrc.h\
     stdio.H

#.......................................................................

# Making a list of source files, except *.h
# =========================================

ALL_SRCS	= $(PSAS_OBJS) $(LIBpsas_OBJS)

ALL_INCL	= $(PSAS_INCL) $(LIBpsas_INCL)

ALL_MODS	= $(LIBpsas_MODS)

ALL_MAKE	=       Makefile              Makefile.conf.eramus\
  Makefile.conf.molotov Makefile.conf.mintz   Makefile.conf.UNICOS\
  getfix.sh	getAImix.sed	getAIpuv.sed	getAIzuv.sed

ALL_DATA	= psas.rc

all_srcs:
	@ for i in $(ALL_SRCS); do \
	    m=`basename $$i .o` ;\
	    for x in .F90 .f90 .F .f .F77 .f77; do \
	      if [ -r $$m$$x ]; then echo $$m$$x; continue 2; fi ;\
	    done; echo "$@: not found, $$i" 1>&2 ;\
	  done | sort -fu

all_incl:
	@ for i in $(ALL_INCL); do \
	    if [ -r $$i ]; then echo $$i ;\
	    else echo "$@: not found, $$i" 1>&2 ; fi ;\
	  done | sort -fu

all_mods:
	@ for i in $(ALL_MODS); do \
	    m=`basename $$i .o` ;\
	    for x in .F90 .f90 .F .f; do \
	      if [ -r $$m$$x ]; then echo $$m$$x; continue 2; fi ;\
	    done; echo "$@: not found, $$i" 1>&2 ;\
	  done | sort -fu


all_make:
	@ for i in $(ALL_MAKE); do \
	    if [ -r $$i ]; then echo $$i ;\
	    else echo "$@: not found, $$i" 1>&2 ; fi ;\
	  done | sort -fu

all_data:
	@ for i in $(ALL_DATA); do \
	    if [ -r $$i ]; then echo $$i ;\
	    else echo "$@: not found, $$i" 1>&2 ; fi ;\
	  done | sort -fu

psas_srcs:
	@ for i in $(PSAS_OBJS); do \
	    m=`basename $$i .o` ;\
	    for x in .F90 .f90 .F .f .F77 .f77; do \
	      if [ -r $$m$$x ]; then echo $$m$$x; continue 2; fi ;\
	    done; echo "$@: not found, $$i" 1>&2 ;\
	  done | sort -fu

libpsas_srcs:
	@ for i in $(LIBpsas_OBJS); do \
	    m=`basename $$i .o` ;\
	    for x in .F90 .f90 .F .f .F77 .f77; do \
	      if [ -r $$m$$x ]; then echo $$m$$x; continue 2; fi ;\
	    done; echo "$@: not found, $$i" 1>&2 ;\
	  done | sort -fu

psas_incl:
	@ for i in $(PSAS_INCL); do \
	    if [ -r $$i ]; then echo $$i ;\
	    else echo "$@: not found, $$i" 1>&2 ; fi ;\
	  done | sort -fu

libpsas_incl:
	@ for i in $(LIBpsas_INCL); do \
	    if [ -r $$i ]; then echo $$i ;\
	    else echo "$@: not found, $$i" 1>&2 ; fi ;\
	  done | sort -fu

#			PSAS Driver
#			===========

NAME	= psas
MODULE	= psasv1

$(NAME).x: $(LIBpsas) $(PSAS_OBJS)
	$(LD) -o $@ $(LDFLAGS) $(PSAS_OBJS) $(LIBpsas) $(LIBsys)

#			 PSAS Tar.Z
#			============

Arch: arch
	@ Arch=$(NAME)-`date "+%d%h%y" | tr "[A-Z]" "[a-z]"` ;\
	  mv $(NAME).tar.Z $${Arch}.tar.Z

arch: do-it-anyway
	@ echo "Making $(NAME).tar.Z ..."
	@ $(RM) $(NAME).tar.Z
	@ $(TAR) -cvf $(NAME).tar \
		`$(MAKE) -f $(MAKEFILE) all_make` \
		`$(MAKE) -f $(MAKEFILE) all_mods` \
		`$(MAKE) -f $(MAKEFILE) all_srcs` \
		`$(MAKE) -f $(MAKEFILE) all_incl` \
		`$(MAKE) -f $(MAKEFILE) all_data`
	@ $(ZIP) -v $(NAME).tar

unar: $(NAME).tar.Z
	@ echo "Untaring $(NAME).tar.Z ..."
	@ $(ZCAT) $(NAME).tar.Z | $(TAR) -xvf -

clean::
	$(RM) stdio.h *.[lo] $(LIBpsas) $(NAME).x core *.mod
	$(RM) getAImix.F getAIpuv.F getAIzuv.F

do-it-anyway:

commit:
	cvs history -m $(MODULE) >History
	@ Arch=$(MODULE)-`date "+%Hz%d%h%y" | tr "[A-Z]" "[a-z]"` ;\
	  cvs commit ;\
	  echo Tagging $(MODULE) with $${Arch} ;\
	  cvs rtag $${Arch} $(MODULE)

#  fixlib replace the existing libpsas.a with a bugfixed version.

fixlib:
	@ Arch=`date "+%d%h%y" | tr "[A-Z]" "[a-z]"` ;\
	  sh ML.sh -r $${Arch}

#  newlib mv the existing libpsas.a to olde/, and replace it with a
#  newer version.

newlib:
	@ Arch=`date "+%d%h%y" | tr "[A-Z]" "[a-z]"` ;\
	  sh ML.sh $${Arch}

rebuild: unar clean $(LIBpsas) $(NAME).x
#.......................................................................
# file dependency list
# ====================

mark_reg.o: stdio.h

# Special files (if any)
# =============

getAIpuv.F: getAIall.F getAIpuv.sed
	sh ./getfix.sh $@

getAIzuv.F: getAIall.F getAIzuv.sed
	sh ./getfix.sh $@

getAImix.F: getAIall.F getAImix.sed
	sh ./getfix.sh $@

config.o: MX_hfecH.h MX_hoecH.h MX_voecH.h MXpar_hc.h kind_covs.h kind_mats.h ktmax.h kxmax.h latmax.h lvmax.h pres4slp.h rade.h stdio.h
Cprodm.o: cordriv.o
FEalpha_imat.o: config.o
FEalpha_tabl.o: config.o
FEsigW_imat.o: config.o
FEsigW_tabl.o: config.o
hfecH_tbl.o: config.o
hfecQ_tbl.o: config.o
hoecH_tbl.o: config.o
OEclass_tbl.o: config.o
rlat_imat.o: config.o
vfecH_tbl.o: config.o
vfecQ_tbl.o: config.o
voecH_tbl.o: config.o

# regular files
# =============

aio_grads.o: aio.h gridxx.h ktmax.h ktwanted.h maxreg.h stdio.h m_inpak90.o
getdel2.o: ktmax.h kxmax.h stdio.h m_inpak90.o
psasshell.o: aio.h bands.h gridxx.h nobsmax.h psasrc.h stdio.h m_inpak90.o

aj_Alf.o: FEalpha_imat.o config.o
aj_Bet.o: config.o
aj_Gam.o: config.o
bands0.o: bands.h stdio.h m_inpak90.o
cgnorm.o: stdio.h realvals.h
conjgr.o: bands.h kind_covs.h kind_mats.h mxpass.h realvals.h stdio.h
Cprod1.o: config.o cordriv.o kind_covs.h
Cprodx.o: config.o cordriv.o kind_covs.h
dervsigF_slD.o: config.o const.o levtabl.h vfecHH.h
dervsigF_slp.o: const.o
dervsigF_slW.o: FEalpha_imat.o config.o const.o levtabl.h rlat_imat.o vfecHH.h
dervsigF_upD.o: config.o const.o levtabl.h vfecHH.h
dervsigF_upW.o: FEalpha_imat.o config.o const.o levtabl.h rlat_imat.o vfecHH.h
diagcorD.o: FEalpha_imat.o config.o cormats.o
diagcorF.o: config.o cormats.o
diagcorM.o: config.o FEsigW_imat.o
diagcorO.o: config.o cormats.o
diagcorU.o: config.o cormats.o
dupelim.o: stdio.h
fDDcor1.o: MX_hfecH.h hfecHH.h kind_covs.h lvmax.h vfecHH.h
fDDcorx.o: MX_hfecH.h hfecHH.h kind_covs.h lvmax.h stdio.h vfecHH.h
fHDcorx.o: MX_hfecH.h hfecHH.h kind_covs.h lvmax.h vfecHH.h
fHHcor1.o: MX_hfecH.h hfecHH.h kind_covs.h lvmax.h vfecHH.h
fHHcorx.o: MX_hfecH.h hfecHH.h kind_covs.h lvmax.h vfecHH.h
fQQcor1.o: hfecQQ.h lvmax.h vfecQQ.h
fQQcorx.o: hfecQQ.h lvmax.h vfecQQ.h
gCprodx.o: config.o cordriv.o kind_covs.h
Gea2ll.o: ktmax.h
getAIall.o: OEclass_tbl.o bands.h config.o ktmax.h ktwanted.h kxmax.h levtabl.h lvmax.h maxreg.h psasrc.h rlat_imat.o m_inpak90.o
getAImix.o: OEclass_tbl.o bands.h config.o ktmax.h ktwanted.h kxmax.h levtabl.h lvmax.h maxreg.h psasrc.h rlat_imat.o m_inpak90.o
getAIpuv.o: OEclass_tbl.o bands.h config.o ktmax.h ktwanted.h kxmax.h levtabl.h lvmax.h maxreg.h psasrc.h rlat_imat.o m_inpak90.o
getAIzuv.o: OEclass_tbl.o bands.h config.o ktmax.h ktwanted.h kxmax.h levtabl.h lvmax.h maxreg.h psasrc.h rlat_imat.o m_inpak90.o
getAinc.o: bands.h gridxx.h ktmax.h levtabl.h lvmax.h maxreg.h stdio.h config.o
getivec.o: config.o
getsigF.o: m_grads.o gridxx.h ktmax.h ktwanted.h pres4slp.h stdio.h
m_grads.o: stdio.h m_inpak90.o m_chars.o
gridxx.o: gridxx.h maxreg.h qea.h stdio.h m_inpak90.o
hcorfuns.o: const.o
hunt.o: stdio.h
imat_alpha.o: FEalpha_imat.o FEalpha_tabl.o config.o const.o levtabl.h rlat_imat.o
imat_sigW.o: FEsigW_imat.o FEsigW_tabl.o config.o const.o levtabl.h vfecHH.h rlat_imat.o
indexGvec.o: config.o
iniainc.o: ktmax.h ktwanted.h stdio.h m_inpak90.o
initkl.o: ktmax.h kttabl.h ktwanted.h kxmax.h kxtabl.h
initRSRC.o: config.o
intp_ctaus.o: config.o
intp_hCor.o: config.o hcorfuns.o
intp_sigF.o: config.o const.o levtabl.h vfecHH.h
intp_sigO.o: OEclass_tbl.o config.o kxtabl.h
intp_vCor.o: config.o
ipick.o: icosdef.h maxreg.h
krname0.o: krname.h maxreg.h
ktname0.o: ktmax.h kttabl.h stdio.h
kxname0.o: ktmax.h kxmax.h kxtabl.h stdio.h
ll2xyz.o: stdio.h
luavail.o: stdio.h
merg_lats.o: config.o
merg_plevs.o: config.o
makedelv.o: bands.h ktmax.h levtabl.h lvmax.h rlat_imat.o stdio.h kind_covs.h
mvPHx.o: config.o ktmax.h kttabl.h levtabl.h lvmax.h rlat_imat.o stdio.h
mv_Alf.o: FEalpha_imat.o config.o
mv_Bet.o: config.o
mv_diag.o: config.o
mv_Gam.o: config.o
obssmry.o: ktmax.h kttabl.h kxmax.h kxtabl.h
obstat.o: ktmax.h kttabl.h realvals.h
offdcorF.o: cormats.o ktmax.h stdio.h
offdcorO.o: cormats.o ktmax.h
oHHcor1.o: MX_hoecH.h MX_voecH.h hoecHH.h ktmax.h kxmax.h kxtabl.h lvmax.h voecHH.h
oHHcorx.o: MX_hoecH.h MX_voecH.h hoecHH.h ktmax.h kxmax.h kxtabl.h lvmax.h voecHH.h
op_Mx.o: FEsigW_imat.o config.o levtabl.h
op_Pf.o: FEsigW_imat.o config.o levtabl.h
proxel.o: ktmax.h kttabl.h kxmax.h kxtabl.h proxel.h stdio.h m_inpak90.o
prtdata.o: ktmax.h kttabl.h kxmax.h kxtabl.h
prtregn.o: ktmax.h kttabl.h kxmax.h kxtabl.h
psasexit.o: stdio.h
psasrcbd.o: psasrc.h
qtrig.o: trigtab.h
qtrig0.o: trigtab.h
rdkttbl.o: realvals.h stdio.h m_inpak90.o
rdkxtbl.o: realvals.h stdio.h m_inpak90.o
rdlevels.o: realvals.h stdio.h m_inpak90.o
rdoetbl.o: realvals.h stdio.h m_inpak90.o
rdpars.o: config.o m_inpak90.o
rdvctbl.o: realvals.h stdio.h m_inpak90.o
rec_Cxpy.o: Cprodm.o config.o
regsort.o: stdio.h
restrict.o: ktmax.h kxmax.h stdio.h
sepang.o: icosdef.h maxreg.h
setbox.o: ktmax.h kxmax.h stdio.h m_inpak90.o
seticos.o: icosdef.h maxreg.h
setpix.o: kxmax.h
set_fecHH.o: config.o hfecHH.h hfecH_tbl.o levtabl.h vfecHH.h vfecH_tbl.o
set_fecQQ.o: config.o hfecQQ.h hfecQ_tbl.o levtabl.h vfecQQ.h vfecQ_tbl.o
set_FEhCor.o: config.o hfecH_tbl.o hfecQ_tbl.o hcorfuns.o
set_FEvCor.o: config.o vfecH_tbl.o vfecQ_tbl.o
set_oecHH.o: config.o hoecHH.h hoecH_tbl.o levtabl.h voecHH.h voecH_tbl.o
set_OEclas.o: OEclass_tbl.o config.o kttabl.h kxtabl.h
set_OEhCor.o: OEclass_tbl.o config.o hoecH_tbl.o kxtabl.h hcorfuns.o
set_OEvCor.o: OEclass_tbl.o config.o kxtabl.h voecH_tbl.o
solve4x.o: bands.h ktmax.h levtabl.h lvmax.h rlat_imat.o stdio.h
sort.o: stdio.h
sparse.o: bands.h config.o kttabl.h
sym_Cxpy.o: Cprodm.o config.o
tabIlist.o: config.o
tabl_FEalpha.o: FEalpha_tabl.o config.o
tabl_FEsigW.o: FEsigW_tabl.o config.o
tabRlist.o: config.o
tabSlist.o: config.o
tofront.o: stdio.h
typsort.o: stdio.h
uHHcor1.o: MX_voecH.h ktmax.h kxmax.h kxtabl.h lvmax.h voecHH.h
xyz2reg.o: icosdef.h maxreg.h
stdv_FEqc.o: config.o levtabl.h psasrc.h rlat_imat.o
stdv_OEqc.o: config.o m_inpak90.o psasrc.h
