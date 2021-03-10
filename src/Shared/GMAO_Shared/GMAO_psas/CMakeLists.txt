esma_set_this ()

set (SRCS_special
  OEclass_tbl.f90
  hfecH_tbl.f90
  hfecQ_tbl.f90
  hoecH_tbl.f90
  rlat_imat.f90
  rlev_imat.F90
  vfecQ_tbl.f90
  voecH_tbl.f90
  FEsigW_imat.f90
  FEsigW_tabl.f90
  FEalpha_imat.f90
  FEalpha_tabl.f90
  )

set (SRCS
  ${SRCS_special}
  config.f90
  const.f90
  indexxs.f90
  inxSlist.f90
  ll2qvec.f90
  rdpars.f90
  rdvctbl.f90
  set_FEhCor.f90
  set_FEvCor.f90
  set_OEvCor.f90
  setbox.f90
  setpix.f90
  slintab.f90
  slogtab.f90
  tabSlist.f90
  tabl_FEalpha.f90
  tabl_FEsigW.f90
  vfecH_tbl.f90
  PSAS_Sigs.F90
  hcorfuns.F90
  imat_alpha.F90
  imat_sigW.F90
  intp_ctaus.F90
  intp_hCor.F90
  intp_sigO.F90
  intp_vCor.F90
  m_AE.F90
  m_AISolver.F90
  m_AIGrid.F90
  m_Attributes.F90
  m_AttributesComm.F90
  m_AttributesMAN.F90
  m_AttrVect.F90
  m_CGSolver.F90
  m_ClusterDistribution.F90
  m_Collector.F90
  m_CollectorComm.F90
  m_CorAttrF.F90
  m_CorAttrX.F90
  m_CorFxpy.F90
  m_CorOxpy.F90
  m_CorUxpy.F90
  m_CorMatxF.F90
  m_CorMatxO.F90
  m_CorMatxU.F90
  m_Checks.F90
  m_Dictionary.F90
  m_Dictionary80.F90
  m_DictionaryTable.F90
  m_Distribution.F90
  m_DistributionComm.F90
  m_EAGrid.F90
  m_ErrCovModels.F90
  m_FEsigFi_tabl.F90
  m_FcstErrCov.F90
  m_FcstErrCovMatx.F90
  m_FcstErrCovModels.F90
  m_GetAI.F90
  m_GlobalPartition.F90
  m_GrADSfiles.F90
  m_HilbertPartMatx.F90
  m_Icosahedron.F90
  m_InnovCovMatx.F90
  m_LLAGrid.F90
  m_LLGridDistribution.F90
  m_LVGrid.F90
  m_LevCache.F90
  m_MultiAccessNavigator.F90
  m_Navigator.F90
  m_ObsErrCovMatx.F90
  m_PartMatx.F90
  m_Partitioner.F90
  m_Regioner.F90
  m_RegionIterator.F90
  m_redwin.F90
  m_SharedDictionary.F90
  m_Sndx.F90
  m_SparseComm.F90
  m_SparseComm_SHMEM_PTR.F90
  m_SparseComm_MPI_1.F90
  m_Spherical_Partition.F90
  m_Spherical_Triangle.F90
  m_StatLevels.F90
  m_Transposer.F90
  m_TransposerComm.F90
  m_aGrid.F90
  m_ai_Operators.F90
  m_block_corD.F90
  m_block_corF.F90
  m_block_corO.F90
  m_block_corU.F90
  m_block_corFxpy.F90
  m_block_corOxpy.F90
  m_block_corUxpy.F90
  m_block_storage.F90
  m_boxes.F90
  m_cache.F90
  m_compPartMatx.F90
  m_costs.F90
  m_diagcorM.F90
  m_elemPart.F90
  m_gammaChi.F90
  m_gammaPhi.F90
  m_gammaPsi.F90
  m_geometry.F90
  m_gdstat.F90
  m_intpAP.F90
  m_intpAPelem.F90
  m_intpAPmiss.F90
  m_keyChain.F90
  m_kt_corF.F90
  m_kt_corO.F90
  m_kt_corU.F90
  m_kt_uvcorF.F90
  m_ktList.F90
  m_llgrad.F90
  m_lookup.F90
  m_mvcorF_bldr.F90
  m_mvcorF_blox.F90
  m_mvcorF_matx.F90
  m_mvcorF_nsep.F90
  m_mvcorF_nsep_bldr.F90
  m_normCor.F90
  m_ob_Operators.F90
  m_output.F90
  m_psas.F90
  m_psasrc.F90
  m_psasversion.F90
  m_qea.F90
  m_recMatx.F90
  m_rectPartMatx.F90
  m_recurPart.F90
  m_rootedAIGrid.F90
  m_scatPartMatx.F90
  m_sigDCWindErr.F90
  m_sigFi_lookups.F90
  m_sigma.F90
  m_sigmaChi.F90
  m_sigmaPhi.F90
  m_sigmaPsi.F90
  m_simplePart.F90
  m_sparse.F90
  m_spherical_geometry.F90
  m_stdv_FE.F90
  m_subBlocks.F90
  m_superObs.F90
  m_symMatx.F90
  m_wGrADS.F90
  m_xOp_sigFi.F90
  m_xRSRC_sigFi.F90
  m_xTab_lats.F90
  m_xTab_levs.F90
  m_xTab_sigFi.F90
  m_xtenlev.F90
  merg_lats.F90
  merg_plevs.F90
  set_OEclas.F90
  set_OEhCor.F90
  set_fecQQ.F90
  set_oecHH.F90
  stdv_FEqc.F90
  stdv_OEqc.F90
  tabIlist.F90
  tabRlist.F90
  m_PHxSolver.F90
  amatch.f
  bands0.f
  cgnorm.f
  dupelim.f
  gdstat.f
  gridxx.f
  hunt.f
  index3r.f
  indexxi.f
  indexxr.f
  initkl.f
  ipick.f
  ktname0.f
  kxname0.f
  listvals.f
  lstins.f
  lvstat.f
  obssmry.f
  obstat.f
  permuti.f
  permutl.f
  permutr.f
  qtrig0.f
  rdkttbl.f
  rdkxtbl.f
  rdlevels.f
  rdoetbl.f
  seticos.f
  smex.f
  sort.f
  spline.f
  vindex.f
  wrtxvec.f
  xyz2reg.f
  getdel2.F
  ll2xyz.F
  regsort.F
  restrict.F
  savemat1.F
  savematx.F
  savevect.F
  tofront.F
  typsort.F
  m_mvcorF_nsep_bmop___.F90
  mp_assign_to_cpu.c
  )

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   set( CMAKE_Fortran_FLAGS_RELEASE "${FOPT3} ${BIG_ENDIAN} ${EXTENDED_SOURCE} ${ALIGNCOM}")
endif ()

set (M4 m4)

add_custom_command (OUTPUT m_mvcorF_nsep_bmop___.F90
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND ${M4} ${CMAKE_CURRENT_SOURCE_DIR}/m_mvcorF_nsep_bmop.m4 > m_mvcorF_nsep_bmop___.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES GMAO_mpeu ${MKL_LIBRARIES}
  )

target_compile_definitions (${this} PUBLIC NDEBUG LANGUAGE_FORTRAN TUNE_COSTS USE_MPI_1 _DLAPACK_)

string(REPLACE " " ";" tmp ${FREAL8})
foreach(flag ${tmp})
  target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach()



ecbuild_add_executable (TARGET solve.x SOURCES solve.F90 LIBS ${this})

string(REPLACE " " ";" tmp ${FREAL8})
foreach(flag ${tmp})
  target_compile_options (solve.x PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach()
