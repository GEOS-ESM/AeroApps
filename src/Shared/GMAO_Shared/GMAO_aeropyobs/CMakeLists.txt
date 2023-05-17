file (GLOB python_files pyobs/*.py)
install (FILES ${python_files} DESTINATION lib/Python/pyobs)

# Python executables
# ------------------
install(PROGRAMS ext_sampler.py  stn_sampler.py  stn_xsect.py  trj_sampler.py  trj_xsect.py DESTINATION bin)

