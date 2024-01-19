install(DIRECTORY fluid/ DESTINATION bin/fluid)
install(DIRECTORY ts_plotting DESTINATION bin/missions/ASIAAQ)

# this is executable
install (
    PROGRAMS  ts_plotting/plot_stn_timeseries.py
    DESTINATION bin/missions/ASIAAQ/ts_plotting
    )
