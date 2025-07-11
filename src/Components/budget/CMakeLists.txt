esma_set_this ()

# these are executable scripts
set (PYSCRIPTS 
    budget.py
    )

install (
    PROGRAMS  ${PYSCRIPTS}
    DESTINATION lib/Python
    )
