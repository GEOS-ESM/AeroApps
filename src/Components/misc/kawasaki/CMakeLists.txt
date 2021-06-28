esma_set_this ()


# these are executable scripts
set (PYSCRIPTS
    m2_anom.py
    m2_mm.py
    mk_npz.py
    plot_ave.py
    plot_full.py
    plot_kd.py
    )


install (
    PROGRAMS  ${PYSCRIPTS}
    DESTINATION bin/${this}
    )

# these are non-executaables
set (PYFILES
    npz.py
    )

set (OTHERFILES
    japan.gs
    japan_population.xlsx
    )

install (
    FILES ${PYFILES} ${OTHERFILES}
    DESTINATION bin/${this}
    )

