esma_set_this ()

# these are executable scripts
set (PYSCRIPTS 
    dust_utils.py
    plot_utils.py
    carma_utils.py
    __init__.py
    )

install (
    PROGRAMS  ${PYSCRIPTS}
    DESTINATION lib/Python/${this}
    )
