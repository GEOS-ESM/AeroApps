esma_set_this()


# these are executable scripts
set (PYSCRIPTS
    polarimeter_swath.py
    run_polarimeter_swath.py
    swath2tile_lb.py
    plot_swath2tile_lb.py
    )

set (SCRIPTS
    setup_env
    run_polarimeter_swath.j
    )

install (
    PROGRAMS  ${PYSCRIPTS} ${SCRIPTS}
    DESTINATION bin/missions/${this}
    )

# these are non-executaables
set (RCFILES
    polar07.pcf
    ssd475.pcf
    swath_files.pcf
    tile_files.pcf
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

#install (
#   DIRECTORY sampling
#   DESTINATION bin/missions/${this}
#   USE_SOURCE_PERMISSIONS
#   MESSAGE_NEVER
#   )

