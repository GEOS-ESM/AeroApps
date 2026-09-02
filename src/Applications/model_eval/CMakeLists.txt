
set(SCRIPTS model_eval_setup.py run_eval.py run_plots.py)

install (
    PROGRAMS  ${SCRIPTS}
    DESTINATION bin
    )

