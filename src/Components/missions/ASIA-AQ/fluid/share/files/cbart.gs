*
*  Script to plot a colorbar
*
*  The script will assume a colorbar is wanted even if there is 
*  not room -- it will plot on the side or the bottom if there is
*  room in either place, otherwise it will plot along the bottom and
*  overlay labels there if any.  This can be dealt with via 
*  the 'set parea' command.  In version 2 the default parea will
*  be changed, but we want to guarantee upward compatibility in
*  sub-releases.
*
*
*	modifications by mike fiorino 940614 and arlindo da silva 160821
*
*	- the extreme colors are plotted as triangles
*	- the colors are boxed in white
*	- input arguments in during a run execution:
* 
*	run cbarn sf vert xmid ymid skip
*
*	sf   - scale the whole bar 1.0 = original 0.5 half the size, etc.
*	vert - 0 FORCES a horizontal bar = 1 a vertical bar
*	xmid - the x position on the virtual page the center the bar
*	ymid - the x position on the virtual page the center the bar
*       skip = how many labels to skip
*
*	if 
*	as in the original algorithm
*  

function colorbar (args)

sf=subwrd(args,1)
vert=subwrd(args,2)
xmid=subwrd(args,3)
ymid=subwrd(args,4)
skip=subwrd(args,5)

* For consistency, replace * with ''
if (sf='*'); sf=''; endif
if (vert='*'); vert=''; endif
if (xmid='*'); xmid=''; endif
if (ymid='*'); ymid=''; endif
if (skip='*'); skip=''; endif

if(skip=''); skip=1; endif

if(sf='');sf=1.0;endif

*
*  Check shading information
*
  'query shades'
  shdinfo = result
  if (subwrd(shdinfo,1)='None') 
    say 'Cannot plot color bar: No shading information'
    return
  endif

* 
*  Get plot size info
*
  'query gxinfo'
  rec2 = sublin(result,2)
  rec3 = sublin(result,3)
  rec4 = sublin(result,4)
  xsiz = subwrd(rec2,4)
  ysiz = subwrd(rec2,6)
  ylo = subwrd(rec4,4)
  xhi = subwrd(rec3,6)
  xd = xsiz - xhi

  ylolim=0.6*sf
  xdlim1=1.0*sf
  xdlim2=1.5*sf  
  barsf=0.8*sf
  yoffset=0.15*sf
  stroff=0.05*sf
  strxsiz=0.12*sf
  strysiz=0.13*sf
*
*  Decide if horizontal or vertical color bar
*  and set up constants.
*
  if (ylo<ylolim & xd<xdlim1) 
    say "Not enough room in plot for a colorbar"
    return
  endif
  cnum = subwrd(shdinfo,5)
*
*	logic for setting the bar orientation with user overides
*
  if (ylo<ylolim | xd>xdlim1)
    vchk = 1
    if(vert = 0) ; vchk = 0 ; endif
  else
    vchk = 0
    if(vert = 1) ; vchk = 1 ; endif
  endif
*
*	vertical bar
*

  if (vchk = 1 )

    xwid = 0.2*sf
    ywid = 0.5*sf

*   if(xmid = '') ; xmid = xhi+xd/2 ; endif
*   if(xmid = '') ; xmid = xhi+xd/4 ; endif
    if(xmid = '') ; xmid = xhi + 0.1 + xwid/2 ; endif
    
    xl = xmid-xwid/2
    xr = xl + xwid
    xleft  = xl
    xright = xr
    if (ywid*cnum > ysiz*barsf) 
      ywid = ysiz*barsf/cnum
    endif
    if(ymid = '') ; ymid = ysiz/2 ; endif
    yb = ymid - ywid*cnum/2
    ybot = yb
    ytop = ymid + ywid*cnum/2
    'set string 1 l 5'
    vert = 1

  else

*
*	horizontal bar
*

    ywid = 0.4
    xwid = 0.8

    if(ymid = '') ; ymid = ylo/2-ywid/2 ; endif
    yt   = ymid + yoffset
    yb   = ymid
    ytop = yt
    ybot = yb
    if(xmid = '') ; xmid = xsiz/2 ; endif
    if (xwid*cnum > xsiz*barsf)
      xwid = xsiz*barsf/cnum
    endif
    xl     = xmid - xwid*cnum/2
    xleft  = xl
    xright = xmid + xwid*cnum/2
    'set string 1 tc 5'
    vert = 0
  endif


*
*  Plot colorbar
*


  'set strsiz 'strxsiz' 'strysiz
  num = 0
  snum = 0
  while (num<cnum) 
    rec = sublin(shdinfo,num+2)
    col = subwrd(rec,1)
    hi = subwrd(rec,3)
    if (vert) 
      yt = yb + ywid
    else 
      xr = xl + xwid
    endif

    xtick = (xr-xl) / 3.0
    ytick = (yt-yb) / 3.0

    if(vert = 1)
      'set line 'col
      'draw recf 'xl' 'yb' 'xr' 'yt
      'set line 1 1 5'
      'draw line 'xr' 'yb' 'xr-xtick' 'yb
      'draw line 'xr' 'yt' 'xr-xtick' 'yt
    else
      'set line 'col
      'draw recf 'xl' 'yb' 'xr' 'yt
      'set line 1 1 5'
      'draw line 'xl' 'yb' 'xl' 'yb+ytick
      'draw line 'xr' 'yb' 'xr' 'yb+ytick
    endif

*   Put numbers under each segment of the color key
    if (num < cnum-1)
      if ( num=snum )
        if (vert) 
          xp=xr+stroff
          'draw string 'xp' 'yt' 'hi
        else
          yp=yb-stroff
         'draw string 'xr' 'yp' 'hi
        endif
        snum = snum + skip
      endif
    endif

*   Reset variables for next loop execution
    if (vert) 
      yb = yt
    else
      xl = xr
    endif
    num = num + 1

  endwhile

  'set line 1 1 5'
* 'draw rec 'xlo' 'ylo' 'xhi' 'yhi
  'draw rec 'xleft' 'ybot' 'xright' 'ytop

return
