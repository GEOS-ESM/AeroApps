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


THIS = $(shell basename `pwd` )
LIB  = lib$(THIS).a

esma_install install: $(LIB) solve.x
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -p *.a         $(ESMALIB)
	$(CP) -p psas.rc.sample        $(ESMAETC)/psas.rc
	$(CP) -p *.mod       $(ESMAINC)/$(THIS)
	$(CP) -p solve.x $(ESMABIN)

esma_clean esma_distclean clean:
	$(RM) *~ *___.F90 *.[aoxd] *.[Mm][Oo][Dd] *.x

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

# Need to pull these out of the library for linking on Darwin
SRCS_special = OEclass_tbl.f90\
hfecH_tbl.f90		\
hfecQ_tbl.f90		\
hoecH_tbl.f90		\
rlat_imat.f90		\
rlev_imat.F90		\
vfecQ_tbl.f90		\
voecH_tbl.f90           \
FEsigW_imat.f90		\
FEsigW_tabl.f90		\
FEalpha_imat.f90	\
FEalpha_tabl.f90	

SRCS	= $(SRCS_special) \
config.f90		\
const.f90		\
indexxs.f90		\
inxSlist.f90		\
ll2qvec.f90		\
rdpars.f90		\
rdvctbl.f90		\
set_FEhCor.f90		\
set_FEvCor.f90		\
set_OEvCor.f90		\
setbox.f90		\
setpix.f90		\
slintab.f90		\
slogtab.f90		\
tabSlist.f90		\
tabl_FEalpha.f90	\
tabl_FEsigW.f90		\
vfecH_tbl.f90		\
PSAS_Sigs.F90		\
hcorfuns.F90		\
imat_alpha.F90		\
imat_sigW.F90		\
intp_ctaus.F90		\
intp_hCor.F90		\
intp_sigO.F90		\
intp_vCor.F90		\
m_AE.F90		\
m_AISolver.F90		\
m_AIGrid.F90		\
m_Attributes.F90	\
m_AttributesComm.F90	\
m_AttributesMAN.F90	\
m_AttrVect.F90		\
m_CGSolver.F90		\
m_ClusterDistribution.F90	\
m_Collector.F90		\
m_CollectorComm.F90	\
m_CorAttrF.F90		\
m_CorAttrX.F90		\
m_CorFxpy.F90		\
m_CorOxpy.F90		\
m_CorUxpy.F90		\
m_CorMatxF.F90		\
m_CorMatxO.F90		\
m_CorMatxU.F90		\
m_Checks.F90		\
m_Dictionary.F90	\
m_Dictionary80.F90	\
m_DictionaryTable.F90	\
m_Distribution.F90	\
m_DistributionComm.F90	\
m_EAGrid.F90		\
m_ErrCovModels.F90	\
m_FEsigFi_tabl.F90	\
m_FcstErrCov.F90	\
m_FcstErrCovMatx.F90	\
m_FcstErrCovModels.F90	\
m_GetAI.F90             \
m_GlobalPartition.F90	\
m_GrADSfiles.F90	\
m_HilbertPartMatx.F90   \
m_Icosahedron.F90       \
m_InnovCovMatx.F90	\
m_LLAGrid.F90		\
m_LLGridDistribution.F90\
m_LVGrid.F90		\
m_LevCache.F90		\
m_MultiAccessNavigator.F90	\
m_Navigator.F90		\
m_ObsErrCovMatx.F90	\
m_PartMatx.F90		\
m_Partitioner.F90	\
m_Regioner.F90		\
m_RegionIterator.F90    \
m_redwin.F90		\
m_SharedDictionary.F90	\
m_Sndx.F90		\
m_SparseComm.F90        \
m_SparseComm_SHMEM_PTR.F90    \
m_SparseComm_MPI_1.F90        \
m_Spherical_Partition.F90 \
m_Spherical_Triangle.F90 \
m_StatLevels.F90	\
m_Transposer.F90	\
m_TransposerComm.F90	\
m_aGrid.F90		\
m_ai_Operators.F90	\
m_block_corD.F90	\
m_block_corF.F90	\
m_block_corO.F90	\
m_block_corU.F90	\
m_block_corFxpy.F90	\
m_block_corOxpy.F90	\
m_block_corUxpy.F90	\
m_block_storage.F90     \
m_boxes.F90		\
m_cache.F90		\
m_compPartMatx.F90	\
m_costs.F90		\
m_diagcorM.F90		\
m_elemPart.F90		\
m_gammaChi.F90		\
m_gammaPhi.F90		\
m_gammaPsi.F90		\
m_geometry.F90          \
m_gdstat.F90		\
m_intpAP.F90		\
m_intpAPelem.F90	\
m_intpAPmiss.F90	\
m_keyChain.F90		\
m_kt_corF.F90		\
m_kt_corO.F90		\
m_kt_corU.F90		\
m_kt_uvcorF.F90		\
m_ktList.F90		\
m_llgrad.F90		\
m_lookup.F90		\
m_mvcorF_bldr.F90	\
m_mvcorF_blox.F90	\
m_mvcorF_matx.F90	\
m_mvcorF_nsep.F90	\
m_mvcorF_nsep_bldr.F90	\
m_normCor.F90		\
m_ob_Operators.F90	\
m_output.F90		\
m_parDOT.F90		\
m_psas.F90              \
m_psasrc.F90		\
m_psasversion.F90	\
m_qea.F90		\
m_recMatx.F90		\
m_rectPartMatx.F90	\
m_recurPart.F90		\
m_rootedAIGrid.F90	\
m_scatPartMatx.F90	\
m_sigDCWindErr.F90	\
m_sigFi_lookups.F90	\
m_sigma.F90		\
m_sigmaChi.F90		\
m_sigmaPhi.F90		\
m_sigmaPsi.F90		\
m_simplePart.F90	\
m_sparse.F90            \
m_spherical_geometry.F90\
m_stdv_FE.F90		\
m_subBlocks.F90		\
m_superObs.F90		\
m_symMatx.F90		\
m_wGrADS.F90		\
m_xOp_sigFi.F90		\
m_xRSRC_sigFi.F90	\
m_xTab_lats.F90		\
m_xTab_levs.F90		\
m_xTab_sigFi.F90	\
m_xtenlev.F90		\
merg_lats.F90		\
merg_plevs.F90		\
set_OEclas.F90		\
set_OEhCor.F90		\
set_fecQQ.F90		\
set_oecHH.F90		\
stdv_FEqc.F90		\
stdv_OEqc.F90		\
tabIlist.F90		\
tabRlist.F90		\
m_PHxSolver.F90         \
amatch.f		\
bands0.f		\
cgnorm.f		\
dupelim.f		\
gdstat.f		\
gridxx.f		\
hunt.f			\
index3r.f		\
indexxi.f		\
indexxr.f		\
initkl.f		\
ipick.f			\
ktname0.f		\
kxname0.f		\
listvals.f		\
lstins.f		\
lvstat.f		\
obssmry.f		\
obstat.f		\
permuti.f		\
permutl.f		\
permutr.f		\
qtrig0.f                \
rdkttbl.f		\
rdkxtbl.f		\
rdlevels.f		\
rdoetbl.f		\
seticos.f               \
smex.f			\
sort.f	                \
spline.f		\
vindex.f		\
wrtxvec.f		\
xyz2reg.f               \
getdel2.F		\
ll2xyz.F		\
regsort.F		\
restrict.F		\
savemat1.F		\
savematx.F		\
savevect.F		\
tofront.F		\
typsort.F               \
m_mvcorF_nsep_bmop___.F90

# m_mvcorF_nsep_bmop.m4

SRCS_c	 = mp_assign_to_cpu.c

OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS)))

OBJS_special := $(addsuffix .o, $(basename $(SRCS_special)))

FOPT = $(FOPT3)
FPE =

INC_DIRS = . $(INC_MPEU)
MOD_DIRS = . $(INC_MPEU)
USER_FFLAGS = $(BIG_ENDIAN) 
USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FDEFS  = $(D)NDEBUG $(D)LANGUAGE_FORTRAN $(D)TUNE_COSTS $(D)USE_MPI_1 $(D)_DLAPACK_


vpath % $(INC_DIRS) # $(MOD_DIRS)

#LIB_SCI = -llapack -lblas
LIB_SCI = -L$(SHARE)/dasilva/lib -llapack -lblas

$(LIB) lib : $(DEPS) $(OBJS)
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

solve.x : $(LIB) solve.o $(LIB_MPEU)
	$(LD) $(LDFLAGS) -o solve.x solve.o $(OBJS_special) $(LIB) \
              $(LIB_MPEU) $(LIB_MPI) $(LIB_SCI) $(LIB_SYS)

%___.F90 : %.m4
	$(M4) $*.m4 > $@

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif


  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
