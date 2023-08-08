wxmap.py --theme `pwd` --fcst_dt 20180517T000000 --time_dt 20180526T230000 --field aerosol --level 0 --region atlantic --lights_off --oname aerosol.png

wxmap.py --theme `pwd` --fcst_dt 20180517T000000 --time_dt 20180526T230000 --field olr --level 0 --region natl --lights_off --oname olr.png

wxmap.py --theme `pwd` --fcst_dt 20180517T000000 --time_dt 20180521T140000 --field epv --level 0 --region nps --lights_off --oname epv.png

wxmap.py --theme `pwd` --fcst_dt 20180517T000000 --time_dt 20180521T140000 --field tpw --level 0 --region ortho --lights_off --oname tpw01.png

wxmap.py --theme `pwd` --fcst_dt 20180517T000000 --time_dt 20180521T140000 --field tpw --level 0 --region australia --lights_off --oname tpw02.png

wxmap.py --theme `pwd` --fcst_dt 20180517T000000 --time_dt 20180521T140000 --field vort --level 500 --region pac --lights_off --oname vorticity.png

wxmap.py --theme `pwd` --fcst_dt 20180517T000000 --time_dt 20180521T140000 --field humid --level 0 --region nam --lights_off --oname humidity.png

wxmap.py --theme `pwd` --time_dt 19930313T180000 --field precip --level 0 --region usa --lights_off --stream MERRA2 --oname precip.png

wxmap.py --theme `pwd` --time_dt 19930313T180000 --field ptype --level 0 --region usa --lights_off --stream MERRA2 --oname precip_type.png

wxmap.py --theme `pwd` --time_dt 19930313T180000 --field wspd --level 850 --region usa --lights_off --stream MERRA2 --oname wspd.png
