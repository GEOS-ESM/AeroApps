firex_plot.py 20190820 0 ../models.rc ../models/GEOSFP.yml 

for model in GEOSFP
qsub firex_catalog_geos.j
qsub firex_catalog_cams.j
qsub firex_catalog_geoscf.j
qsub firex_catalog_camchem.j
qsub firex_catalog_raqms.j
qsub firex_catalog_wrfchem.j
qsub firex_catalog_hrrr.j
qsub firex_catalog_uclawrf.j
qsub firex_catalog_arqi.j
#qsub firex_catalog_uiowawrf.j
qsub firex_catalog_all.j

##qsub firex_catalog_arqi.j

