esma_set_this()


# these are executable scripts
set (PYSCRIPTS
    imager_swath.py
    imager_sampler.py
    run_imager_sampler.py
    g5nr_imager_sampler.py
    )

set (SCRIPTS
    setup_env
    run_imager_sampler.j
    )

install (
    PROGRAMS  ${PYSCRIPTS} ${SCRIPTS}
    DESTINATION bin/missions/${this}
    )

# these are non-executaables
set (RCFILES
    Aod_EOS.rc
    Chem_MieRegistry.rc
    sbg.pcf
    ssd650.pcf
    imager_files.pcf
    imager_sampler.pcf
    )

install (
    FILES ${RCFILES}
    DESTINATION bin/missions/${this}
    )

#install (
#   DIRECTORY rc
#   DESTINATION bin/missions/${this}
#   USE_SOURCE_PERMISSIONS
#   MESSAGE_NEVER
#   )

install (
   DIRECTORY tle
   DESTINATION bin/missions/${this}
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )

install (
   DIRECTORY sampling
   DESTINATION bin/missions/${this}
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )

