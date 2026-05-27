esma_set_this ()

# these are executable scripts
set (PYSCRIPTS 
    budget.py
    __init__.py
    )

install (
    PROGRAMS  ${PYSCRIPTS}
    DESTINATION lib/Python/${this}
    )
