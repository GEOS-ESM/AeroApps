*re-initialze everything
'reinit'

expid_ = x0034

showmem_ = 1
showinc_ = 0
showspread_ = 1
doplot_ = 0
dotend_ = 0

t0_ = 1
tmax_ = 15

nfiles_ = 0

if (dotend_ > 0)
  cint_ = 5
else
  cint_ = 3
* cint_ = 4
endif
*tmax_ = 34

if ( expid_ = "x0034_rt" | expid_ = "x0034" )
 'open ./DDF/'expid_'.inst3_3d_wxme_Np_ens.ctl'
  nfiles_ = nfiles_ + 1
 'xdfopen ./DDF/'expid_'.asm_inst3d_met_p.ddf'
  nfiles_ = nfiles_ + 1
  region_ = tpac
  stormname = "TPAC"
  region_ = tatl
  stormname = "TATL"
* region_ = eaus
* stormname = "WPAC"
* region_ = glob
* stormname = "GLO"
* region_ = trop
* stormname = "TRO"
* region_ = spol
* stormname = "SPOLE"
  date_ = 06z01aug2018
  if (dotend_ > 0)
     cint_ = 5
  else
     cint_ = 4
  endif
  t0_ = 20
  tmax_ = 21
endif

if (region_ = glob )
  'set lat  -90 90'
  'set lon -180  180'
endif
if (region_ = panw )
  'set lat  34 60'
  'set lon -170 -115'
endif
if (region_ = eusa )
  'set lat  10 40'
  'set lon -100 -50'
endif
if (region_ = wpac )
 'set lat  5 35'
 'set lon 100 160'
endif
if (region_ = tatl )
 'set lat  12 42'
*'set lon -90 -40'
 'set lon -110 -20'
endif
if (region_ = tpac )
 'set lat  -35 35'
 'set lon 100 180'
endif
if (region_ = trop )
 'set lat  -30 30'
 'set lon -180 180'
endif
if (region_ = eaus )
 'set lat -40 -10'
 'set lon 140 180'
endif
if (region_ = spac )
 'set lat -60 -10'
 'set lon 50 180'
endif
if (region_ = spol )
 'set mproj sps'
 'set lat -90 -40'
 'set lon -180 180'
endif
if (region_ = npol )
 'set mproj nps'
 'set lat  40  90'
 'set lon -180 180'
endif

* open land-mask file (as last file)
'open /home/dao_ops/GEOSadas-CURRENT/GEOSadas/src/GMAO_Shared/GEOS_Util/plots/grads_util/lwmask1440721.tabl'
nfiles_ = nfiles_ + 1

if ( expid_ = "x0027" | expid_ = "x0027B" | expid_ = "x0027C" ) 
  figname_ = "NotRecentered_spaghetti"
else
  figname_ = stormname"_spaghetti"
endif

var_ = slp
lev_ = 1000
clev1_ = 5400
clev2_ = 5700
clev1_ = 245
clev2_ = 255
clev1_ = 1017
clev2_ = 1017

if(dotend_ > 0)
  factor1_ = 4
  factor2_ = 4
else
  factor1_ = 1
  factor2_ = 1
endif
if (var_ = slp)
  if(dotend_ > 0)
    factor1_ = 0.04
    factor2_ = 0.04
* typical of central ...
*   factor1_ = 4
*   factor2_ = 4
  else
    factor1_ = 0.01
*   factor2_ = 0.01
*   factor1_ = 1
    factor2_ = 1
  endif
endif

memb_ = 1
meme_ = 32

*'set display color white'
'set display color black'
'c'
'set vpage 0.5 11 0.5 8.0'
'set grid off'
*'set mproj lambert'

* Set domain, level, and time (these should be passed as options)
* ---------------------------------------------------------------
*'set mpdset hires'
*'set mproj nps'
*'set lat  10 40'
*'set lon -100 -50'
*'set mproj orthogr'
*'set lev 500'
'set time 'date_


nt = t0_
'set t 't0_
while (nt < tmax_ )
 'c'
 'set grads off'
 'set dfile 1'
  if ( var_ != "slp" )
    'set lev 'lev_
  endif
 'set e 1'

'set gxout shaded'

* Display colored map
* -------------------
if (dotend_ = 0 )
*'set mpdset mres'
*'set poli on'
'set xlopts 1 5 .15'
'set ylopts 1 5 .15'
'set rgb 98 139 69 19'
'set map  98 1 1'
'set datawarn off'
'set rgb 92 255 227 171'
'set ccols 92 0'
'set clevs 0.0025'
'set dfile 'nfiles_
'd lwmask.'nfiles_'(t=1,z=1)'
endif

'set gxout contour'
'set dfile 1'
if ( var_ != "slp" )
  'set lev 'lev_
endif

* Ensemble mean
* -------------
  if(dotend_ > 0)
    'define ensmean=ave('var_'(lev='lev_'),e='memb_',e='meme_')-ave('var_'(lev='lev_',t-1),e='memb_',e='meme_')'
  else
    'define ensmean=ave('var_'(lev='lev_'),e='memb_',e='meme_')'
  endif

* Ensemble spread
* ---------------
  rms = 'pow('var_'-ensmean,2)'
  variance = 'ave('rms',e='memb_',e='meme_')'
  'define spread=sqrt('variance')'

* plot spread
* -----------
'set gxout shaded'
'set datawarn off'
'set rgb 40 127 255 0'
'set rgb 41 0   205 0'
'set rgb 42 0 139 0'
'set rgb 43 16 78 139'
'set rgb 44 30 144 255'
'set rgb 45 0 178 238'
'set rgb 46 0 238 238'
'set rgb 47 137 104 205'
'set rgb 48 145 44 238'
'set rgb 49 139 0 139'
'set rgb 50 139 0 0'
'set rgb 51 205 0 0'
'set rgb 52 238 64 0'
'set rgb 53 255 127 0'
'set rgb 54 205 133 0'
'set rgb 55 255 215 0'
'set rgb 56 238 238 0'
'set rgb 57 255 255 0'

if ( showinc_ = 1 )
  'set clevs 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4 4.8 5.2 5.6 6.0'
  'set ccols  0 40 41 42 43  44  45  46  47  48 49 50 51 52 53 54 '
  'display (ps.4 - ps.3)/100'
else
   if (showspread_ = 1)
    'set clevs 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6'
    'set ccols  0 40 41 42 43  44  45  46  47  48 49 50 51 52 53 '
    'define spr = spread * 'factor1_
    'display maskout(spr, abs(spr) -0.4)'
    'run cbarn.gs 0.8 0 5.5 0.3'
  endif
endif

* Spaghetti-like plot
* -------------------
'set gxout contour'
'set clab off'
e1=memb_
e2=meme_
if (showmem_ = 1)
   while (e1<=e2)
   'set e 'e1
   if (clev1_ = clev2_)
     'set ccolor '15
     'set cint 'cint_
      if ( showinc_ = 0 )
        'd 'factor1_'*'var_
      endif
   else
     'set ccolor 2'
     'set clevs 'clev1_
     'd 'factor1_'*'var_
     'set ccolor 4'
     'set clevs 'clev2_
     'd 'factor1_'*'var_
   endif
   e1=e1+1
   endwhile
endif


* Central Analysis/Assimilation
* -----------------------------
if ( showinc_ = 0 & nfiles_ > 2 )
   'set dfile 2'
    if ( var_ != "slp" )
      'set lev 'lev_
    endif
   if (dotend_ > 0)
      'set gxout shaded'
      'set rbrange -30 30'
      'set black -'cint_' 'cint_
   else
      'set ccolor 2'
      'set ccolor 13'
   endif
   'set cthick 6'
   'set e 1'
   'set cint 'cint_
   'set clab on'
   if (dotend_ = 1 | dotend_ = 3)
      'd 'factor2_'*('var_'.2-'var_'.2(t-1))'
      'cbar'
   else
      'd 'factor2_'*'var_'.2'
   endif
endif

* Central Background
* ------------------
if ( showinc_ = 0 &nfiles_ > 3 )
  'set dfile 3'
  if ( var_ != "slp" )
    'set lev 'lev_
  endif
  'set cthick 6'
  'set ccolor 9'
  'set e 1'
  'set cint 'cint_
  'set clab on'
  'd slp.3'
endif

* Plot ensemble mean
* ------------------
if ( showinc_ = 0 )
   if(dotend_ = 2)
     'set gxout shaded'
     'set rbrange -30 30'
     'set black -'cint_' 'cint_
   else
     'set gxout contour'
     'set ccolor 1'
   endif
   'set cthick 6'
   if (clev1_ = clev2_)
    'set cint 'cint_
   else
   'set clevs ' clev1_ ' ' clev2_
   endif
   'set clab off'
   if (dotend_ != 1)
     'd 'factor1_'*ensmean'
      if (dotend_ != 0)
         'cbar'
      endif
   endif
endif

* Figure title
* ------------
'q time'
fdate=subwrd(result,3)
fday=subwrd(result,6)
*'draw title Central (aqua) EnMean & Sprd (shaded) 'stormname ' 'fday' 'fdate
if (dotend_ > 0 )
  if (dotend_ = 1)
    'draw title 'var_' (mb/dy) Central 'stormname ' 'fday' 'fdate
  endif
  if (dotend_ = 2)
    'draw title 'var_' (mb/dy) EnsMean 'stormname ' 'fday' 'fdate
  endif
  if (dotend_ = 3)
    'draw title 'var_' (mb/dy) Central & EnsMean (cnt) 'stormname ' 'fday' 'fdate
  endif
else
* 'draw title Central (aqua) Ens-Mean & Members 'stormname ' 'fday' 'fdate
  'draw title  'fday' 'fdate
endif


if (doplot_ = 1)
  'enable print fig.gx'
  'print        fig.gx'
  'disable print'
* '!/gpfsm/dnb52/projects/p10/dasilva/opengrads/Contents/gxyat fig.gx -x 1100 -y 850'
* '!/gpfsm/dnb52/projects/p10/dasilva/opengrads/Contents/gxyat fig.gx -x 2048 -y 1024'
  '!/gpfsm/dnb52/projects/p10/dasilva/opengrads/Contents/gxyat -r fig.gx -x 1440 -y 720'
   if (nt < 10 )
     '!/bin/mv fig.png 'expid_'_'var_'_'figname_'00'nt'.png'
   endif
   if ( nt > 9 & nt < 100)
     '!/bin/mv fig.png 'expid_'_'var_'_'figname_'0'nt'.png'
   endif
   if ( nt > 99 )
     '!/bin/mv fig.png 'expid_'_'var_'_'figname_''nt'.png'
   endif
endif

  nt = nt + 1
  'set t 'nt
endwhile
