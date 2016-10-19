"""
   This module contains auxiliary functions used by the MODIS Collection 6 Neural Net Retrieval.

   P. Castellanos, 2016.

"""

import os, sys
from   matplotlib.pyplot    import  cm, imshow, plot, figure
from   matplotlib.pyplot    import  xlabel, ylabel, title, grid, savefig, legend
import matplotlib.pyplot    as      plt
from   matplotlib.ticker    import  MultipleLocator
import matplotlib.patches   as      mpatches
from   numpy                import  c_ as cat
from   numpy                import  random, sort, pi, load, cos, log, std, exp
from   numpy                import  reshape, arange, ones, zeros, interp, sqrt
from   numpy                import  meshgrid, concatenate, squeeze
import numpy                as      np
import itertools
from   sklearn.linear_model import LinearRegression
# ------


#---------------------------------------------------------------------

def boxplot_imshow(data,plottype,blocks,masterlist,title,filename,
                   vseps=None, yrange=None,ylabel=None):
    nvars,ncomb = blocks.shape

    fig = plt.figure()
    ax  = plt.subplot(211)
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)

    if plottype is 'box':
      bp = plt.boxplot(data,showfliers=False,showbox=True,whis='range',
                whiskerprops={'linestyle':'-'})
      plt.setp(bp['boxes'], color='black')
      plt.setp(bp['whiskers'], color='black')
      plt.plot([0,ncomb+0.5],[np.median(data[:,0]),np.median(data[:,0])],color='b',ls='--',zorder=5)
    elif plottype is 'scatter':
      scat = np.mean(data,axis=0)
      plt.plot(np.arange(ncomb)+1,scat,'rD')
      plt.plot([0,ncomb+0.5],[np.mean(data[:,0]),np.mean(data[:,0])],color='b',ls='--',zorder=5,markersize=5)
    elif plottype is 'errorbar':
      scat = np.mean(data,axis=0)
      yerr_max = np.abs(np.max(data,axis=0)-scat)
      yerr_min = np.abs(np.min(data,axis=0)-scat)
      plt.errorbar(np.arange(ncomb)+1,scat,yerr=[yerr_min,yerr_max], ls='none',marker='D',color='r',ecolor='k',markersize=5)   
      plt.plot([0,ncomb+0.5],[np.mean(data[:,0]),np.mean(data[:,0])],color='b',ls='--',zorder=5)


    ax.set_xlim(0.5,ncomb+0.5)
    ax.set_xticks(np.arange(ncomb)+0.5)
    ax.set_xticklabels([])    

    if yrange is not None:
        ax.set_ylim(yrange)

    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)

    # Minor Y Tick Marks    
    # ylim = ax.get_ylim()
    # dy   = (ylim[1] - ylim[0])/(len(yticks)-1)
    # minorLocator = MultipleLocator(0.5*dy)
    # ax.yaxis.set_minor_locator(minorLocator)

    if ylabel is not None:
      ax.set_ylabel(ylabel,fontsize=14)
    # Make vertical lines to separate bins of number of inputvars
    if vseps is not None:
        for v in vseps:
            ax.plot([v,v],np.array(ax.get_ylim()),'k-')

    plt.title(title)
    #ax.minorticks_on()
    plt.grid(True,axis='y',which='both',color='0.5',linestyle='-')
    ax.set_axisbelow(True)  #Grid lines go to the back.
    plt.tick_params(
        axis='y',          # changes apply to the y-axis
        which='major',     # major ticks are affected
        direction='out',
        right='off') 


    axblocks = plt.subplot(212)
    plt.imshow(blocks,interpolation='none',aspect='auto')
    axblocks.set_yticks(np.arange(nvars))
    axblocks.set_yticklabels(masterlist)
    axblocks.set_yticks(np.arange(nvars)+0.5,minor=True)
    axblocks.set_xticks(np.arange(ncomb)+0.5,minor=True)
    # plt.draw()  # this is needed because get_window_extent needs a renderer to work
    # yax = axblocks.get_yaxis()
    # # find the maximum width of the label on the major ticks
    # pad = max(T.label.get_window_extent().width for T in yax.majorTicks)
    # yax.set_tick_params(pad=pad)

    plt.tick_params(
        axis='both',          # changes apply to both
        which='major',     # major ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        left='off',
        right='off',
        labelbottom='off') # labels along the bottom edge are off 
    plt.grid(True,which='minor',color='0.5',linestyle='-')
    axblocks.set_axisbelow(True)  #Grid lines go to the back.
    xlim = axblocks.get_xlim()
    if vseps is not None:
        for v in vseps:
            axblocks.plot([v-1,v-1],np.array(axblocks.get_ylim()),'k-')
    axblocks.set_xlim(xlim)
    
    plt.tight_layout()
    plt.subplots_adjust(right=0.99,hspace=0.001)
    fig.set_size_inches(10.5, 5.5)   
    plt.savefig(filename, transparent='true',dpi=300)
    #plt.show()
    plt.close(fig)

# ---
def SummarizeCombinations(mxd,Input_nnr,yrange=None,sortname='rmse'):
    """
    Create summary plot
    """
    outdir = mxd.outdir

    # Flatten Input_nnr into one list
    # -------------------------------
    Input = ()
    for i in Input_nnr:
      if type(i) is list:
        Input = Input + tuple(i)
      else:
        Input = Input + (i,)

    input_nnr = list(Input)

    expid      = '.'.join(input_nnr)
    ident      = mxd.ident
    nvars      = len(input_nnr)
    ngroups    = len(Input_nnr)
    ncomb      = len(mxd.comblist)
    blocks     = np.zeros([nvars,ncomb])
    masterlist = np.array(tuple(input_nnr))

    for c,comb in enumerate(mxd.comblist):
        blocks[:,c] = np.array([v == masterlist for v in comb]).sum(axis=0)

    blocks = np.insert(blocks,0,np.zeros(nvars),axis=1)  #for the original MODIS data

    newlist = []
    for var in masterlist:
        if mxd.sat == 'Terra':
          try:
            newlist.append(dict(MODVARNAMES,**VARNAMES)[var])
          except:
            newlist.append(var)
        else:
          try:
            newlist.append(dict(MYDVARNAMES,**VARNAMES)[var])
          except:
            new.append(var)

    masterlist = np.array(newlist)

    #nblocks = blocks.sum(axis=0)
    nblocks = [len(group) for group in mxd.combgroups]
    nblocks.insert(0,0)

    print "MASTERLIST", masterlist

    #--------------
    # Default sort by mean RMSE
    #------------------
    sortvalue = mxd.nnr.__dict__[sortname]
    sortvalue  = np.insert(sortvalue,0,np.array(mxd.orig.__dict__[sortname][:,0]),axis=1)
    sortmetric = np.median(sortvalue, axis=0)
    isort = np.empty(0,dtype='int')
    vseps = np.empty(0,dtype='int')
    for i in np.sort(np.unique(nblocks)):
        istart  = np.where(nblocks == i)[0].min()
        vseps   = np.append(vseps,istart+0.5)
        isort   = np.append(isort,istart + np.argsort(sortmetric[nblocks == i]))
    blocks = blocks[:,isort]

    for i in np.arange(nvars):
        blocks[i,blocks[i,:] == 1] = i+1

    blocks = np.ma.masked_array(blocks)
    blocks.mask = [blocks == 0]


    def getplotdata(mxd,varname):
      vardata = mxd.nnr.__dict__[varname]
      vardata = np.insert(vardata,0,np.array(mxd.orig.__dict__[varname][:,0]),axis=1)
      vardata = vardata[:,isort]
      
      if vardata.shape[0] == 1:
        vardata = np.append(vardata,vardata,axis=0)

      return vardata

    if mxd.nnr.slope.shape[0] == 1:
      plottype = 'scatter'
    elif mxd.nnr.slope.shape[0] >= 5:
      plottype = 'box'
    else:
      plottype = 'errorbar'


    boxplot_imshow(getplotdata(mxd,'slope'),plottype,
                   blocks,masterlist,
                   'Slope',
                   '{}/Slope.{}.{}.png'.format(outdir,expid,ident),
                   vseps = vseps,
                   yrange=[0,1])

    boxplot_imshow(getplotdata(mxd,'R'),plottype,
                   blocks,masterlist,
                   'R',
                   '{}/R.{}.{}.png'.format(outdir,expid,ident),
                   vseps = vseps,
                   yrange=[0,1])

    boxplot_imshow(getplotdata(mxd,'intercept'),plottype,
                   blocks,masterlist,
                   'Intercept',
                   '{}/Intercept.{}.{}.png'.format(outdir,expid,ident),
                   vseps = vseps,
                   yrange=yrange)

    boxplot_imshow(getplotdata(mxd,'rmse'),plottype,
                   blocks,masterlist,
                   'RMSE',
                   '{}/RMSE.{}.{}.png'.format(outdir,expid,ident),
                   vseps = vseps,
                   yrange=yrange)

    boxplot_imshow(getplotdata(mxd,'me'),plottype,
                   blocks,masterlist,
                   'Mean Bias',
                   '{}/ME.{}.{}.png'.format(outdir,expid,ident),
                   vseps = vseps,
                   yrange=yrange)

    boxplot_imshow(getplotdata(mxd,'mae'),plottype,
                   blocks,masterlist,
                   'Mean Absolute Error',
                   '{}/MAE.{}.{}.png'.format(outdir,expid,ident),
                   vseps = vseps,
                   yrange=yrange)    

#--------------------------------------------------------------------------------------

def get_Iquartiles(mxd,I=None):


  if I is None:
    Irange     = arange(mxd.nobs)
    I          = [Irange[mxd.iValid]]

  I1 = []
  I2 = []
  I3 = []
  I4 = []
  for iTest in I:
    targets  = mxd.getTargets(iTest)
    if len(targets.shape) > 1:
      targets = targets[:,0]

    p25  = np.percentile(targets,25)
    p50  = np.percentile(targets,50)
    p75  = np.percentile(targets,75)
    p100 = targets.max()

    I1.append(iTest[targets <= p25])
    I2.append(iTest[(targets>p25) & (targets<=p50)])
    I3.append(iTest[(targets>p50) & (targets<=p75)])
    I4.append(iTest[(targets>p75)])

  return I1, I2, I3, I4
#---------------------------------------------------------------------
def get_Ispecies(mxd,I=None):

  if I is None:
    Irange     = arange(mxd.nobs)
    I          = [Irange[mxd.iValid]]

  Ifdu = []
  Ifss = []
  Ifcc = []
  Ifsu = []
  for iTest in I:
    fdu  = mxd.fdu[iTest].squeeze()
    fss  = mxd.fss[iTest].squeeze()
    fcc  = mxd.fcc[iTest].squeeze()
    fsu  = mxd.fsu[iTest].squeeze()
    
    Ifdu.append(iTest[fdu > 0.5])
    Ifss.append(iTest[fss > 0.5])
    Ifcc.append(iTest[fcc > 0.5])
    Ifsu.append(iTest[fsu > 0.5])

  return Ifdu, Ifss, Ifcc, Ifsu

#---------------------------------------------------------------------
def get_ImRef(mxd,refName,refmin,refmax,I=None):

  if I is None:
    Irange     = arange(mxd.nobs)
    I          = [Irange[mxd.iValid]]


  ImRef = []
  for iTest in I:
    mRef  = mxd.__dict__[refName][iTest].squeeze()

    if refmax is None:
      ImRef.append(iTest[(mRef >= refmin)])
    else:
      ImRef.append(iTest[(mRef >= refmin) & (mRef < refmax)])

  return ImRef

#---------------------------------------------------------------------

def make_plots(mxd,expid,ident,I=None):  
  outdir = mxd.outdir
  if I is None:
    I = ones(mxd.lon.shape).astype(bool)
  # Plot KDE of corrected AOD
  # -------------------------
  # mxd.plotKDE(I=I,figfile=expid+"."+ident+"_kde-"+mxd.Target[0][1:]+"-corrected.png")
  targets  = mxd.getTargets(I).squeeze()
  results = mxd.eval(I).squeeze()
  _plotKDE(targets,results,y_label='NNR')
  title("Log("+mxd.Target[0][1:]+"+0.01)- "+ident)
  savefig(outdir+"/"+expid+"."+ident+"_kde-"+mxd.Target[0][1:]+'-corrected.png')

  # Plot KDE of uncorrected AOD
  # ---------------------------   
  original = log(mxd.mTau550[I]+0.01)
  _plotKDE(targets,original,y_label='Original MODIS')
  title("Log("+mxd.Target[0][1:]+"+0.01)- "+ident)
  savefig(outdir+"/"+expid+"."+ident+"_kde-"+mxd.Target[0][1:]+'.png')

  # Scatter diagram for testing
  # ---------------------------
  mxd.plotScat(I=I,figfile=outdir+"/"+expid+"."+ident+"_scat-"+mxd.Target[0][1:]+'.png')


#---------------------------------------------------------------------
def make_error_pdfs(mxd,Input,expid,ident,K=None,I=None,Title=None,netfileRoot=None,
                    emin=-1.5,emax=2.5):  
  outdir = mxd.outdir
  if I is None:
    I = mxd.iValid

  # Plot PDF of Error
  # -------------------------
  if K is None:
    targets  = [mxd.getTargets(I)]
    if len(targets[0].shape) > 1:
      targets[0] = targets[0][:,0]

    results  = [mxd.eval(I)[:,0]]
    original = [mxd.mTau550[I]]

    if mxd.laod:
      original[0] = log(original[0] + 0.01)

    mod04RMSE = rmse(original,targets)
    nnrRMSE   = rmse(results,targets)
  else:
    targets   = []
    original  = []
    results   = []
    mod04RMSE = []
    nnrRMSE   = []
    for k,iTest in enumerate(I):
      # Irange     = arange(mxd.nobs)
      # iValid     = Irange[mxd.iValid]
      # mxd.iTest  = iValid[iTest]
      mxd.iTest    = iTest

      targets.append(mxd.getTargets(mxd.iTest))
      if len(targets[k].shape) > 1:
            targets[k] = targets[k][:,0]

      original.append(mxd.mTau550[mxd.iTest])

      inputs = mxd.getInputs(mxd.iTest,Input=Input)

      knet = mxd.loadnet(netfileRoot+'.k={}_Tau.net'.format(str(k+1)))
      out = knet(inputs)[:,0]
      results.append(out)

      if mxd.laod:
        original[k] = log(original[k] + 0.01)

      mod04RMSE.append(rmse(original[k],targets[k]))
      nnrRMSE.append(rmse(results[k],targets[k]))

    print 'mod04RMSE',mod04RMSE
    print 'nnrRMSE',nnrRMSE
    mod04RMSE = np.mean(mod04RMSE)
    nnrRMSE   = np.mean(nnrRMSE)

  eorig = []
  ecorr = []
  for o,t,r in zip(original,targets,results):
    eorig.append(o - t)
    ecorr.append(r - t)

  if emax is None:
    emax  = np.array([[e.max() for e in eorig],[e.max() for e in ecorr]]).max()
  if emin is None:
    emin  = np.array([[e.min() for e in eorig],[e.min() for e in ecorr]]).min()

  nbins        = 100
  corrected = []
  orig      = []
  if K is None:
    cc, x = np.histogram(ecorr[0],bins=np.linspace(emin,emax,nbins+1),density=True)
    oo, x = np.histogram(eorig[0],bins=np.linspace(emin,emax,nbins+1),denstiry=True)

    corrected.append(c)
    orig.append(o)
  else:
    for k,iTest in enumerate(I):
      cc, x = np.histogram(ecorr[k],bins=np.linspace(emin,emax,nbins+1),density=True)
      oo, x = np.histogram(eorig[k],bins=np.linspace(emin,emax,nbins+1),density=True)

      corrected.append(cc)
      orig.append(oo)

    corrected = np.array(corrected).T
    orig      = np.array(orig).T

  xcen = x[:-1] + 0.5*(x[1:] - x[:-1])
  if K is not None:
    xcen = np.tile(xcen,(K,1)).T

  fig = plt.figure()
  ax  = plt.subplot(111)  
  ax.plot(xcen,orig,color='k')
  ax.plot(xcen,corrected,color='r')
  ax.set_xlim(emin,emax)
  ax.set_ylim(0,3.5)
  orig_patch      = mpatches.Patch(color='k', label='MOD04 RMSE={:1.2F}'.format(mod04RMSE))
  corrected_patch = mpatches.Patch(color='r', label='NNR RMSE={:1.2F}'.format(nnrRMSE) )
  ax.legend(handles=[orig_patch,corrected_patch])
  plt.grid(True, which='major',axis='x',color='0.50',linestyle='-')

  if Title is None:
    title("Error Log("+mxd.Target[0][1:]+"+0.01)")
  else:
    title(Title)
  savefig(outdir+"/error_pdf-"+expid+"."+ident+"-"+mxd.Target[0][1:]+'.png')  
  plt.close(fig)
#---------------------------------------------------------------------  

def TestStats(mxd,K,C):
    if K is None:
      k = 0
    else:
      k = K

    if C is None:
      c = 0
    else:
      c = C

    # regression[0,2] = slope, intercept, r-value
    out, reg = mxd.test(iprint=False)
    out = out[:,0]

    mxd.nnr.slope[k,c]     = reg[0][0]
    mxd.nnr.intercept[k,c] = reg[0][1]
    mxd.nnr.R[k,c]         = reg[0][2]

    targets  = mxd.getTargets(mxd.iTest)
    if len(targets.shape) > 1:
      targets = targets[:,0]
    original = mxd.mTau550[mxd.iTest]

    if mxd.laod:
      original = log(original + 0.01)

    mxd.nnr.rmse[k,c] = rmse(out,targets)
    mxd.nnr.mae[k,c]  = mae(out,targets)
    mxd.nnr.me[k,c]   = me(out,targets)
 
    lm = LinearRegression()
    targets.shape = targets.shape + (1,)
    lm.fit(targets,original)
    mxd.orig.slope[k,c]     = lm.coef_[0]
    mxd.orig.intercept[k,c] = lm.intercept_
    mxd.orig.R[k,c]         = sqrt(lm.score(targets,original))

    targets  = targets.squeeze()
    mxd.orig.rmse[k,c] = rmse(original,targets)
    mxd.orig.mae[k,c]  = mae(original,targets)
    mxd.orig.me[k,c]   = me(original,targets)

    
# ---
def rmse(predictions, targets):
    return sqrt((np.square(predictions - targets)).mean())
# ---
def mae(predictions, targets):
    return np.abs(predictions-targets).mean()
# ---
def me(predictions, targets):
    return (predictions-targets).mean()    
#---------------------------------------------------------------------
def SummaryPDFs(mxdx):
    K = mxdx.K
    if K is None:
      I = [mxdx.iTest]
    else:
      I = []
      Irange     = arange(mxdx.nobs)
      iValid     = Irange[mxdx.iValid]
      for iTrain, iTest in mxdx.kf:          
        I.append(iValid[iTest])

    I1, I2, I3, I4 = get_Iquartiles(mxdx,I=I)
    Ifdu, Ifss, Ifcc, Ifsu = get_Ispecies(mxdx,I=I)

    Ifna = zeros(mxdx.lon.shape).astype(bool)
    for f in (mxdx.fdu,mxdx.fss,mxdx.fcc,mxdx.fsu):
      J = f>0.5                 # all obs for which species dominate
      Ifna[J] = True

    if K is None:
      Ifna   = [Irange[~Ifna]]
    else:
      Ifna   = [Irange[~Ifna]]*K

    for c,Input in enumerate(mxdx.comblist):
      for invars in itertools.permutations(Input):
        netfileRoot = mxdx.outdir+"/"+'.'.join(invars)+"."+mxdx.ident
        filelist = glob(netfileRoot+'*.net')
        if len(filelist > 0):
          Input = invars
          break

      if len(filelist == 0):
        print '{} not found.  Need to train this combinatin of inputs'.format(netfileRoot+'*.net')
        raise
        
      make_error_pdfs(mxdx,Input,'.'.join(Input),mxdx.ident,K=K,I=I,netfileRoot=netfileRoot)

      make_error_pdfs(mxdx,Input,'.'.join(Input)+'.q25',mxdx.ident,K=K,I=I1,
                     Title="Q1 Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot)
      make_error_pdfs(mxdx,Input,'.'.join(Input)+'.q50',mxdx.ident,K=K,I=I2,
                      Title="Q2 Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot)
      make_error_pdfs(mxdx,Input,'.'.join(Input)+'.q75',mxdx.ident,K=K,I=I3,
                      Title="Q3 Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot)
      make_error_pdfs(mxdx,Input,'.'.join(Input)+'.q100',mxdx.ident,K=K,I=I4,
                      Title="Q4 Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot)

      make_error_pdfs(mxdx,Input,'.'.join(Input)+'.fdu',mxdx.ident,K=K,I=Ifdu,
                     Title="Dust Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot)
      make_error_pdfs(mxdx,Input,'.'.join(Input)+'.fss',mxdx.ident,K=K,I=Ifss,
                      Title="Sea Salt Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot)
      make_error_pdfs(mxdx,Input,'.'.join(Input)+'.fcc',mxdx.ident,K=K,I=Ifcc,
                      Title="BC+OC Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot)
      make_error_pdfs(mxdx,Input,'.'.join(Input)+'.fsu',mxdx.ident,K=K,I=Ifsu,
                      Title="Sulfate Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot)   

      make_error_pdfs(mxdx,Input,'.'.join(Input)+'.fna',mxdx.ident,K=K,I=Ifna,
                      Title="No Dominant Species Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot)                               

    if mxdx.surface == 'ocean':
      varname  = 'mRef870'
      ImRefLow = get_ImRef(mxdx,varname,0,0.02,I=I)
      ImRefMed = get_ImRef(mxdx,varname,0.02,0.05,I=I)      
      ImRefHigh = get_ImRef(mxdx,varname,0.05,None,I=I)
    elif mxdx.surface == 'land':
      varname  = 'mRef870'
      ImRefLow = get_ImRef(mxdx,varname,0,0.2,I=I)
      ImRefMed = get_ImRef(mxdx,varname,0.2,0.30,I=I)      
      ImRefHigh = get_ImRef(mxdx,varname,0.30,None,I=I)
    elif mxdx.surface == 'dbl':
      varname  = 'mRef660'
      ImRefLow = get_ImRef(mxdx,'ScatteringAngle',-1.0,-0.8,I=I)
      ImRefMed = get_ImRef(mxdx,'ScatteringAngle',-0.8,-0.5,I=I)      
      ImRefHigh = get_ImRef(mxdx,'ScatteringAngle',-0.50,None,I=I)      

    for c,Input in enumerate(mxdx.comblist):
      netfileRoot = mxdx.outdir+"/"+'.'.join(Input)+"."+mxdx.ident
      make_error_pdfs(mxdx,Input,'.'.join(Input)+'.mRef870',mxdx.ident,K=K,I=ImRefLow,
                     Title="Low Reflectance Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot)

      make_error_pdfs(mxdx,Input,'.'.join(Input)+'.mRef870',mxdx.ident,K=K,I=ImRefMed,
                     Title="Med Reflectance Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot)        

      make_error_pdfs(mxdx,Input,'.'.join(Input)+'.mRef870',mxdx.ident,K=K,I=ImRefHigh,
                     Title="High Reflectance Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot)            