
set(SCRIPTS ctl2reference.py ctl2reference_run.j)

install (
    PROGRAMS  ${SCRIPTS}
    DESTINATION bin
    )

set(CONFIGS ctl2reference.yaml)

install (
    FILES ${CONFIGS}
    DESTINATION etc
    )
