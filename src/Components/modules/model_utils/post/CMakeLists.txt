esma_set_this ()

# these are executable scripts
set (PYSCRIPTS 
    blend.py
    monthly_mod.py
    monthly_obs.py
    )

install (
    PROGRAMS  ${PYSCRIPTS}
    DESTINATION lib/Python/model_utils/${this}
    )
