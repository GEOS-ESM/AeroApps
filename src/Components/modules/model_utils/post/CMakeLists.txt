esma_set_this ()

# these are executable scripts
set (PYSCRIPTS 
    blend.py
    monthly_mod.py
    monthly_obs.py
    __init__.py
    )

install (
    PROGRAMS  ${PYSCRIPTS}
    DESTINATION lib/Python/${this}
    )
