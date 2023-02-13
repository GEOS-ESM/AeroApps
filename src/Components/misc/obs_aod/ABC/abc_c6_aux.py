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
from   numpy                import  reshape, arange, ones, zeros, interp
from   numpy                import  meshgrid, concatenate, squeeze
import numpy                as      np
import itertools
from   sklearn.linear_model import LinearRegression
from   glob                 import glob
from   scipy                import stats
from   nn                   import _plotKDE


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
        right=False) 


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
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        left=False,
        right=False,
        labelbottom=False) # labels along the bottom edge are off 
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
    plt.savefig(filename, transparent=True,dpi=300)
    #plt.show()
    plt.close(fig)

# ---
def SummarizeCombinations(mxd,Input_nnr,yrange=None,sortname='rmse'):
    """
    Create summary plot
    """
    outdir = mxd.plotdir
    if not os.path.exists(outdir):
      os.makedirs(outdir)
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
    retrieval  = mxd.retrieval
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
            newlist.append(var)

    masterlist = np.array(newlist)

    #nblocks = blocks.sum(axis=0)
    nblocks = [len(group) for group in mxd.combgroups]
    nblocks.insert(0,0)

    print "MASTERLIST", masterlist

    #--------------
    # Default sort by mean RMSE of first target
    #------------------
    sortvalue = mxd.nnr.__dict__[sortname][:,:,0]
    sortvalue  = np.insert(sortvalue,0,np.array(mxd.orig.__dict__[sortname][:,0,0]),axis=1)
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


    def getplotdata(mxd,varname,iTarget):
      vardata = mxd.nnr.__dict__[varname][:,:,iTarget]
      vardata = np.insert(vardata,0,np.array(mxd.orig.__dict__[varname][:,0,iTarget]),axis=1)
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

    for t in range(mxd.nTarget):
      name = mxd.Target[t][1:]
      boxplot_imshow(getplotdata(mxd,'slope',t),plottype,
                   blocks,masterlist,
                   'Slope',
                   '{}/Slope.{}.{}.{}.png'.format(outdir,retrieval,ident,name),
                   vseps = vseps,
                   yrange=[0,1])

      boxplot_imshow(getplotdata(mxd,'R',t),plottype,
                   blocks,masterlist,
                   'R',
                   '{}/R.{}.{}.{}.png'.format(outdir,retrieval,ident,name),
                   vseps = vseps,
                   yrange=[0,1])

      boxplot_imshow(getplotdata(mxd,'intercept',t),plottype,
                   blocks,masterlist,
                   'Intercept',
                   '{}/Intercept.{}.{}.{}.png'.format(outdir,retrieval,ident,name),
                   vseps = vseps,
                   yrange=yrange)

      boxplot_imshow(getplotdata(mxd,'rmse',t),plottype,
                   blocks,masterlist,
                   'RMSE',
                   '{}/RMSE.{}.{}.{}.png'.format(outdir,retrieval,ident,name),
                   vseps = vseps,
                   yrange=yrange)

      boxplot_imshow(getplotdata(mxd,'me',t),plottype,
                   blocks,masterlist,
                   'Mean Bias',
                   '{}/ME.{}.{}.{}.png'.format(outdir,retrieval,ident,name),
                   vseps = vseps,
                   yrange=yrange)

      boxplot_imshow(getplotdata(mxd,'mae',t),plottype,
                   blocks,masterlist,
                   'Mean Absolute Error',
                   '{}/MAE.{}.{}.{}.png'.format(outdir,retrieval,ident,name),
                   vseps = vseps,
                   yrange=yrange)    

#--------------------------------------------------------------------------------------

def get_Iquartiles(mxd,I=None,var=None):


  if I is None:
    Irange     = arange(mxd.nobs)
    I          = [Irange[mxd.iValid]]

  I1 = []
  I2 = []
  I3 = []
  I4 = []
  for iTest in I:
    if var is None:
      targets  = mxd.getTargets(iTest)
    else:
      targets = mxd.__dict__[var][iTest].squeeze()
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
  outdir = mxd.plotdir
  if not os.path.exists(outdir):
    os.makedirs(outdir)  
  if I is None:
    I = ones(mxd.lon.shape).astype(bool)
  # Plot KDE of corrected AOD
  # -------------------------
  # mxd.plotKDE(I=I,figfile=expid+"."+ident+"_kde-"+mxd.Target[0][1:]+"-corrected.png")
  targets  = mxd.getTargets(I)
  results = mxd.eval(I)
  # loop through targets
  for t in range(mxd.nTarget):
      fig = _plotKDE(targets[:,t],results[:,t],y_label='NNR')
      title("Log("+mxd.Target[t][1:]+"+0.01)- "+ident)
      savefig(outdir+"/"+expid+"."+ident+"_kde-"+mxd.Target[t][1:]+'-corrected.png')
      plt.close(fig)

  # Plot KDE of uncorrected AOD
  # ---------------------------  
  # loop through targets
  # plot if there's a corresponding MODIS retrieval
  for t in range(mxd.nTarget):
      name = 'm'+mxd.Target[t][1:]
      if name in mxd.__dict__:
          original = log(mxd.__dict__[name][I]+0.01)
          fig = _plotKDE(targets[:,t],original,y_label='Original MODIS')
          title("Log("+mxd.Target[t][1:]+"+0.01)- "+ident)
          savefig(outdir+"/"+expid+"."+ident+"_kde-"+mxd.Target[t][1:]+'.png')
          plt.close(fig)

          # Scatter diagram for testing
          # ---------------------------
          mxd.plotScat(iTarget=t,I=I,figfile=outdir+"/"+expid+"."+ident+"_scat-"+mxd.Target[t][1:]+'.png')

  # If more than one target
  # Plot KDE of Angstrom Exponent
  # Implicitly requires AOD 550
  # -------------------------------
  if mxd.nTarget > 1:
      bins = np.arange(0., 3., 0.05 )
      # AE from Corrected AOD
      # ----------------------
      refname = 'Tau550'
      refwav  = 550.
      # find index of reference
      for t in range(mxd.nTarget):
          name = mxd.Target[t][1:]
          if name == refname:
              reft = np.exp(targets[:,t]) # keep the + 0.01 to handle negatives in MODIS data
              refr = np.exp(results[:,t]) # keep the + 0.01 to handle negatives in MODIS data
      # calculate AE with respect to 550
      for t in range(mxd.nTarget):
          name = mxd.Target[t][1:]
          if name != refname:
              print 't,wav',t,name[3:]
              wav = float(name[3:])
              tt = np.exp(targets[:,t]) # keep the + 0.01 to handle negatives in MODIS data
              rr = np.exp(results[:,t]) # keep the + 0.01 to handle negatives in MODIS data
              AEt = -1.*np.log(reft/tt)/np.log(refwav/wav)
              AEr = -1.*np.log(refr/rr)/np.log(refwav/wav)

              fig = _plotKDE(AEt,AEr,y_label='NNR',x_bins=bins,y_bins=bins)
              title("AE 550/"+name[4:])
              savefig(outdir+"/"+expid+"."+ident+"_kde-AE"+name[3:]+'-corrected.png')
              plt.close(fig)
      
      # AE from uncorrected AOD
      # ------------------------
      # find index of reference
      refname = 'm' + refname
      for t in range(mxd.nTarget):
          name = 'm'+mxd.Target[t][1:]
          if name in mxd.__dict__:
              if name == refname:
                  refo = mxd.__dict__[name][I] + 0.01 # add 0.01 to handle negatives
      # calculate AE with respect to 550
      for t in range(mxd.nTarget):
          name = 'm'+mxd.Target[t][1:]
          if name in mxd.__dict__:
              if name != refname:
                  print 'orig t,wav',t,name[4:]
                  wav = float(name[4:])
                  oo = mxd.__dict__[name][I] + 0.01 # add 0.01 to handle negatives
                  tt = np.exp(targets[:,t]) # keep + 0.01 to handle negatives
                  AEo = -1.*np.log(refo/oo)/np.log(refwav/wav)
                  AEt = -1.*np.log(reft/tt)/np.log(refwav/wav)
                  fig = _plotKDE(AEt,AEo,y_label='Original MODIS',x_bins=bins,y_bins=bins)
                  title("AE 550/"+name[3:])
                  savefig(outdir+"/"+expid+"."+ident+"_kde-AE"+name[4:]+'.png')
                  plt.close(fig)

#---------------------------------------------------------------------

def make_plots_angstrom(mxd,expid,ident,I=None):
  outdir = mxd.plotdir
  if not os.path.exists(outdir):
    os.makedirs(outdir)
  if I is None:
    I = ones(mxd.lon.shape).astype(bool)

  # get the AOD from the predicted angstrom exponent
  # ------------------------------------------------
  targets  = mxd.getTargets(I)
  results = mxd.eval(I)
  base_wav = mxd.AE_base_wav
  for t in range(mxd.nTarget):
      tname = mxd.Target[t]
      if 'Tau' in tname:
          base_name = tname
          base_tau_t  = np.exp(targets[:,t]) - 0.01
          base_tau_r  = np.exp(results[:,t]) - 0.01

  for t in range(mxd.nTarget):
      tname = mxd.Target[t]
      if 'AE' in tname:
          wav = float(tname.split('AE')[-1])
          tau = base_tau_t*np.exp(-1.*np.log(wav/base_wav)*targets[:,t])
          targets[:,t] = np.log(tau + 0.01)

          tau = base_tau_r*np.exp(-1.*np.log(wav/base_wav)*results[:,t])
          results[:,t] = np.log(tau + 0.01)

  # Plot KDE of corrected AOD
  # -------------------------
  # loop through targets
  for t in range(mxd.nTarget):
      fig = _plotKDE(targets[:,t],results[:,t],y_label='NNR')
      title("Log("+mxd.Target[t][1:]+"+0.01)- "+ident)
      tname = mxd.Target[t]
      if 'AE' in tname:
          wavs = tname.split('AE')[-1]
          name = 'Tau'+wavs
      else:
          name = tname[1:]
      savefig(outdir+"/"+expid+"."+ident+"_kde-"+name+'-corrected.png')
      plt.close(fig)

  # Plot KDE of uncorrected AOD
  # ---------------------------
  # loop through targets
  # plot if there's a corresponding MODIS retrieval
  for t in range(mxd.nTarget):
      tname = mxd.Target[t]
      if 'AE' in tname:
          wavs = tname.split('AE')[-1]
          name = 'mTau'+ wavs
      else:
          name = 'm'+mxd.Target[t][1:]
      if name in mxd.__dict__:
          original = log(mxd.__dict__[name][I]+0.01)
          fig = _plotKDE(targets[:,t],original,y_label='Original MODIS')
          title("Log("+mxd.Target[t][1:]+"+0.01)- "+ident)
          savefig(outdir+"/"+expid+"."+ident+"_kde-"+name[1:]+'.png')
          plt.close(fig)

          if name == 'mTau550':
              # Scatter diagram for testing
              # ---------------------------
              mxd.plotScat(iTarget=t,I=I,figfile=outdir+"/"+expid+"."+ident+"_scat-"+name+'.png')

  # If more than one target
  # Plot KDE of Angstrom Exponent
  # Implicitly requires AOD 550
  # -------------------------------
  if mxd.nTarget > 1:
      bins = np.arange(0., 3., 0.05 )
      # AE from Corrected AOD
      # ----------------------
      refname = 'Tau550'
      refwav  = 550.
      # find index of reference
      for t in range(mxd.nTarget):
          name = mxd.Target[t][1:]
          if name == refname:
              reft = np.exp(targets[:,t]) # keep the + 0.01 to handle negatives in MODIS data
              refr = np.exp(results[:,t]) # keep the + 0.01 to handle negatives in MODIS data
      # calculate AE with respect to 550
      for t in range(mxd.nTarget):
          tname = mxd.Target[t]
          if 'AE' in tname:
              wavs = tname.split('AE')[-1]
              name = 'Tau'+wavs
          else:         
              name = mxd.Target[t][1:]
          if name != refname:
              print 't,wav',t,name[3:]
              wav = float(name[3:])
              tt = np.exp(targets[:,t]) # keep the + 0.01 to handle negatives in MODIS data
              rr = np.exp(results[:,t]) # keep the + 0.01 to handle negatives in MODIS data
              AEt = -1.*np.log(reft/tt)/np.log(refwav/wav)
              AEr = -1.*np.log(refr/rr)/np.log(refwav/wav)

              fig = _plotKDE(AEt,AEr,y_label='NNR',x_bins=bins,y_bins=bins)
              title("AE 550/"+name[4:])
              savefig(outdir+"/"+expid+"."+ident+"_kde-AE"+wavs+'-corrected.png')
              plt.close(fig)

      # AE from uncorrected AOD
      # ------------------------
      # find index of reference
      refname = 'm' + refname
      for t in range(mxd.nTarget):
          name = 'm'+mxd.Target[t][1:]
          if name in mxd.__dict__:
              if name == refname:
                  refo = mxd.__dict__[name][I] + 0.01 # add 0.01 to handle negatives
      # calculate AE with respect to 550
      for t in range(mxd.nTarget):
          tname = mxd.Target[t]
          if 'AE' in tname:
              wavs = tname.split('AE')[-1]
              name = 'mTau'+wavs
          else:
              name = 'm'+mxd.Target[t][1:]
          if name in mxd.__dict__:
              if name != refname:
                  print 'orig t,wav',t,name[4:]
                  wav = float(name[4:])
                  oo = mxd.__dict__[name][I] + 0.01 # add 0.01 to handle negatives
                  tt = np.exp(targets[:,t]) # keep + 0.01 to handle negatives
                  AEo = -1.*np.log(refo/oo)/np.log(refwav/wav)
                  AEt = -1.*np.log(reft/tt)/np.log(refwav/wav)
                  fig = _plotKDE(AEt,AEo,y_label='Original MODIS',x_bins=bins,y_bins=bins)
                  title("AE 550/"+name[3:])
                  savefig(outdir+"/"+expid+"."+ident+"_kde-AE"+name[4:]+'.png')
                  plt.close(fig)

#---------------------------------------------------------------------
def make_error_pdfs(mxd,Input,expid,ident,K=None,I=None,Title=None,netfileRoot=None,
                    emin=-1.5,emax=2.5):  
  outdir = mxd.plotdir + '/{}/'.format('.'.join(Input))
  if not os.path.exists(outdir):
    os.makedirs(outdir)

  if I is None:
    I = [mxd.iValid]

  # Plot PDF of Error
  # -------------------------
  # loop through targets
  # plot if there's a corresponding MODIS retrieval 
  for t in range(mxd.nTarget):
    tname = 'm'+mxd.Target[t][1:]
    if 'Tau' not in tname:
        wavs = tname.split('AE')[-1]
        wav = float(wavs)
        name = 'mTau'+wavs
    else:
        name = tname
    if name in mxd.__dict__:
      if K is None:
        targets  = mxd.getTargets(I[0])[:,t]

        inputs = mxd.getInputs(I[0],Input=Input)
        knet = mxd.loadnet(netfileRoot+'_Tau.net')
        results = knet(inputs)[:,t]

        if 'Tau' not in tname:
            base_tau_t = np.exp(mxd.getTargets(I[0])[:,mxd.AE_base_wav_i]) - 0.01
            tau = base_tau_t*np.exp(-1.*np.log(wav/mxd.AE_base_wav)*targets)
            targets = np.log(tau + 0.01)

            base_tau_r = np.exp(knet(inputs)[:,mxd.AE_base_wav_i]) - 0.01
            tau = base_tau_r*np.exp(-1.*np.log(wav/mxd.AE_base_wav)*results)
            results = np.log(tau + 0.01)

        original = mxd.__dict__[name][I[0]]

        if mxd.laod:
          original = log(original + 0.01)

        mod04RMSE = rmse(original,targets)
        nnrRMSE   = rmse(results,targets)

        targets = [targets]
        results = [results]
        original = [original]
      else:
        targets   = []
        original  = []
        results   = []
        mod04RMSE = []
        nnrRMSE   = []
        for k,iTest in enumerate(I):
          mxd.iTest    = iTest

          targets.append(mxd.getTargets(mxd.iTest)[:,t])

          inputs = mxd.getInputs(mxd.iTest,Input=Input)

          knet = mxd.loadnet(netfileRoot+'.k={}_Tau.net'.format(str(k+1)))
          out = knet(inputs)[:,t]
          results.append(out)

          if 'Tau' not in tname:
              base_tau_t = np.exp(mxd.getTargets(mxd.iTest)[:,mxd.AE_base_wav_i]) - 0.01
              tau = base_tau_t*np.exp(-1.*np.log(wav/mxd.AE_base_wav)*targets[k])
              targets[k] = np.log(tau + 0.01)

              base_tau_r = np.exp(knet(inputs)[:,mxd.base_wav_i]) - 0.01
              tau = base_tau_r*np.exp(-1.*np.log(wav/mxd.AE_base_wav)*results[k])
              results[k] = np.log(tau + 0.01)

          original.append(mxd.__dict__[name][mxd.iTest])

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
      for o,tar,r in zip(original,targets,results):
        eorig.append(o - tar)
        ecorr.append(r - tar)

      if emax is None:
        emax  = np.array([[e.max() for e in eorig],[e.max() for e in ecorr]]).max()
      if emin is None:
        emin  = np.array([[e.min() for e in eorig],[e.min() for e in ecorr]]).min()

      nbins        = 100
      corrected = []
      orig      = []
      x = np.linspace(emin,emax,nbins+1)
      xcen = x[:-1] + 0.5*(x[1:] - x[:-1])
      if K is None:
        # cc, x = np.histogram(ecorr[0],bins=np.linspace(emin,emax,nbins+1),density=True)
        # oo, x = np.histogram(eorig[0],bins=np.linspace(emin,emax,nbins+1),denstiry=True)

        kernel = stats.gaussian_kde(ecorr[0])
        cc     = kernel(xcen)

        kernel = stats.gaussian_kde(eorig[0])
        oo     = kernel(xcen)    

        corrected.append(cc)
        orig.append(oo)
      else:
        for k,iTest in enumerate(I):
          # cc, x = np.histogram(ecorr[k],bins=np.linspace(emin,emax,nbins+1),density=True)
          # oo, x = np.histogram(eorig[k],bins=np.linspace(emin,emax,nbins+1),density=True)

          kernel = stats.gaussian_kde(ecorr[k])
          cc     = kernel(xcen)

          kernel = stats.gaussian_kde(eorig[k])
          oo     = kernel(xcen)          

          corrected.append(cc)
          orig.append(oo)

      corrected = np.array(corrected).T
      orig      = np.array(orig).T


      if K is not None:
         xcen = np.tile(xcen,(K,1)).T

      fig = plt.figure()
      ax  = plt.subplot(111)  
      ax.plot(xcen,orig,color='k')
      ax.plot(xcen,corrected,color='r')
      ax.set_xlim(emin,emax)

      if mxd.surface == 'ocean':
        ax.set_ylim(0,3.5)
      else:
        ax.set_ylim(0,2.0)
      orig_patch      = mpatches.Patch(color='k', label='MOD04 RMSE={:1.2F}'.format(mod04RMSE))
      corrected_patch = mpatches.Patch(color='r', label='NNR RMSE={:1.2F}'.format(nnrRMSE) )
      ax.legend(handles=[orig_patch,corrected_patch])
      plt.grid(True, which='major',axis='both',color='0.50',linestyle='-')

      if Title is None:
        title("Error Log("+name[1:]+"+0.01)")
      else:
        title(Title)
      savefig(outdir+"/error_pdf-"+expid+"."+ident+"-"+name[1:]+'.png')  
      plt.close(fig)
#---------------------------------------------------------------------  
#---------------------------------------------------------------------
def make_error_pdfs_int(mxd,Input,expid,ident,K=None,I=None,Title=None,netfileRoot=None,
                    emin=-1.5,emax=2.5):  
  outdir = mxd.plotdir + '/{}/'.format('.'.join(Input))
  if not os.path.exists(outdir):
    os.makedirs(outdir)

  if I is None:
    I = [mxd.iValid]

  # Plot PDF of Error
  # -------------------------
  if K is None:
    targets  = [mxd.getTargets(I[0])]
    if len(targets[0].shape) > 1:
      targets[0] = targets[0][:,0]

    # results    = [mxd.eval(I)[:,0]]
    inputs = mxd.getInputs(I[0],Input=Input)
    knet = mxd.loadnet(netfileRoot)
    out = knet(inputs)[:,0]
    results = [out]


    original   = [mxd.mTau550[I[0]]]
    dboriginal = [mxd.dbmTau550[I[0]]]    

    if mxd.laod:
      original[0]   = log(original[0] + 0.01)
      dboriginal[0] = log(dboriginal[0] + 0.01)

    mod04RMSE   = rmse(original[0],targets[0])
    dbmod04RMSE = rmse(dboriginal[0],targets[0])
    nnrRMSE     = rmse(results[0],targets[0])
  else:
    targets     = []
    original    = []
    dboriginal    = []
    results     = []
    mod04RMSE   = []
    dbmod04RMSE = []
    nnrRMSE     = []
    for k,iTest in enumerate(I):
      # Irange     = arange(mxd.nobs)
      # iValid     = Irange[mxd.iValid]
      # mxd.iTest  = iValid[iTest]
      mxd.iTest    = iTest

      targets.append(mxd.getTargets(mxd.iTest))
      if len(targets[k].shape) > 1:
            targets[k] = targets[k][:,0]

      original.append(mxd.mTau550[mxd.iTest])
      dboriginal.append(mxd.dbmTau550[mxd.iTest])

      inputs = mxd.getInputs(mxd.iTest,Input=Input)

      knet = mxd.loadnet(netfileRoot+'.k={}_Tau.net'.format(str(k+1)))
      out = knet(inputs)[:,0]
      results.append(out)

      if mxd.laod:
        original[k]   = log(original[k] + 0.01)
        dboriginal[k] = log(dboriginal[k] + 0.01)

      mod04RMSE.append(rmse(original[k],targets[k]))
      dbmod04RMSE.append(rmse(dboriginal[k],targets[k]))
      nnrRMSE.append(rmse(results[k],targets[k]))

    print 'mod04RMSE',mod04RMSE
    print 'dbmod04RMSE',dbmod04RMSE
    print 'nnrRMSE',nnrRMSE
    mod04RMSE   = np.mean(mod04RMSE)
    dbmod04RMSE = np.mean(dbmod04RMSE)
    nnrRMSE     = np.mean(nnrRMSE)

  eorig   = []
  edborig = []
  ecorr   = []
  for o,dbo,t,r in zip(original,dboriginal,targets,results):
    eorig.append(o - t)
    edborig.append(dbo - t)
    ecorr.append(r - t)

  if emax is None:
    emax  = np.array([[e.max() for e in eorig],[e.max() for e in edborig],[e.max() for e in ecorr]]).max()
  if emin is None:
    emin  = np.array([[e.min() for e in eorig],[e.min() for e in edborig],[e.min() for e in ecorr]]).min()

  nbins        = 100
  corrected   = []
  orig        = []
  dborig      = []
  x = np.linspace(emin,emax,nbins+1)
  xcen = x[:-1] + 0.5*(x[1:] - x[:-1])  
  if K is None:
    # cc, x   = np.histogram(ecorr[0],bins=np.linspace(emin,emax,nbins+1),density=True)
    # oo, x   = np.histogram(eorig[0],bins=np.linspace(emin,emax,nbins+1),denstiry=True)
    # dboo, x = np.histogram(edborig[0],bins=np.linspace(emin,emax,nbins+1),denstiry=True)

    kernel = stats.gaussian_kde(ecorr[0])
    cc     = kernel(xcen)

    kernel = stats.gaussian_kde(eorig[0])
    oo     = kernel(xcen)    

    kernel = stats.gaussian_kde(edborig[0])
    dboo   = kernel(xcen)    

    corrected.append(cc)
    orig.append(oo)
    dborig.append(dboo)
  else:
    for k,iTest in enumerate(I):
      # cc, x   = np.histogram(ecorr[k],bins=np.linspace(emin,emax,nbins+1),density=True)
      # oo, x   = np.histogram(eorig[k],bins=np.linspace(emin,emax,nbins+1),density=True)
      # dboo, x = np.histogram(edborig[k],bins=np.linspace(emin,emax,nbins+1),density=True)

      kernel = stats.gaussian_kde(ecorr[k])
      cc     = kernel(xcen)

      kernel = stats.gaussian_kde(eorig[k])
      oo     = kernel(xcen)    

      kernel = stats.gaussian_kde(edborig[k])
      dboo   = kernel(xcen)    

      corrected.append(cc)
      orig.append(oo)
      dborig.append(dboo)

  corrected   = np.array(corrected).T
  orig        = np.array(orig).T
  dborig      = np.array(dborig).T


  if K is not None:
    xcen = np.tile(xcen,(K,1)).T

  fig = plt.figure()
  ax  = plt.subplot(111)  
  ax.plot(xcen,orig,color='k')
  ax.plot(xcen,dborig,color='g')
  ax.plot(xcen,corrected,color='r')
  ax.set_xlim(emin,emax)
  if mxd.surface == 'ocean':
    ax.set_ylim(0,3.5)
  else:
    ax.set_ylim(0,2.0)
  orig_patch        = mpatches.Patch(color='k', label='MOD04 DT RMSE={:1.2F}'.format(mod04RMSE))
  dborig_patch      = mpatches.Patch(color='g', label='MOD04 DB RMSE={:1.2F}'.format(dbmod04RMSE))  
  corrected_patch   = mpatches.Patch(color='r', label='NNR RMSE={:1.2F}'.format(nnrRMSE) )
  ax.legend(handles=[orig_patch,dborig_patch,corrected_patch])
  plt.grid(True, which='major',axis='both',color='0.50',linestyle='-')

  if Title is None:
    title("Error Log("+mxd.Target[0][1:]+"+0.01)")
  else:
    title(Title)
  savefig(outdir+"/error_int_pdf-"+expid+"."+ident+"-"+mxd.Target[0][1:]+'.png')  
  plt.close(fig)
#---------------------------------------------------------------------  
#---------------------------------------------------------------------
def make_error_pdfs_dbdt(mxd,mxd2,Input,expid,ident,K=None,I=None,Title=None,
                         netfileRoot=None,netfileRoot2=None,Input2=None,
                         emin=-1.5,emax=2.5):  
  outdir = mxd.plotdir + '/{}/'.format('.'.join(Input))
  if not os.path.exists(outdir):
    os.makedirs(outdir)

  if I is None:
    I = [mxd.iValid]

  # Plot PDF of Error
  # -------------------------
  if K is None:
    targets  = [mxd.getTargets(I[0])]
    if len(targets[0].shape) > 1:
      targets[0] = targets[0][:,0]

    # results    = [mxd.eval(I)[:,0]]
    # results2    = [mxd2.eval(I)[:,0]]
    inputs = mxd.getInputs(I[0],Input=Input)
    knet = mxd.loadnet(netfileRoot)
    out = knet(inputs)[:,0]
    results = [out]

    inputs = mxd2.getInputs(I[0],Input=Input2)
    knet = mxd2.loadnet(netfileRoot2)
    out = knet(inputs)[:,0]
    results2 = [out]

    original   = [mxd.mTau550[I[0]]]
    dboriginal = [mxd.dbmTau550[I[0]]]    

    if mxd.laod:
      original[0]   = log(original[0] + 0.01)
      dboriginal[0] = log(dboriginal[0] + 0.01)

    mod04RMSE   = rmse(original[0],targets[0])
    dbmod04RMSE = rmse(dboriginal[0],targets[0])
    nnrRMSE     = rmse(results[0],targets[0])
    nnrRMSE2     = rmse(results2[0],targets[0])
  else:
    targets     = []
    original    = []
    dboriginal    = []
    results     = []
    results2     = []
    mod04RMSE   = []
    dbmod04RMSE = []
    nnrRMSE     = []
    nnrRMSE2     = []
    for k,iTest in enumerate(I):
      # Irange     = arange(mxd.nobs)
      # iValid     = Irange[mxd.iValid]
      # mxd.iTest  = iValid[iTest]
      mxd.iTest    = iTest

      targets.append(mxd.getTargets(mxd.iTest))
      if len(targets[k].shape) > 1:
            targets[k] = targets[k][:,0]

      original.append(mxd.mTau550[mxd.iTest])
      dboriginal.append(mxd.dbmTau550[mxd.iTest])

      inputs = mxd.getInputs(mxd.iTest,Input=Input)

      knet = mxd.loadnet(netfileRoot+'.k={}_Tau.net'.format(str(k+1)))
      out = knet(inputs)[:,0]
      results.append(out)

      inputs = mxd2.getInputs(mxd.iTest,Input=Input2)
      knet = mxd2.loadnet(netfileRoot2+'.k={}_Tau.net'.format(str(k+1)))
      out = knet(inputs)[:,0]
      results2.append(out)

      if mxd.laod:
        original[k]   = log(original[k] + 0.01)
        dboriginal[k] = log(dboriginal[k] + 0.01)

      mod04RMSE.append(rmse(original[k],targets[k]))
      dbmod04RMSE.append(rmse(dboriginal[k],targets[k]))
      nnrRMSE.append(rmse(results[k],targets[k]))
      nnrRMSE2.append(rmse(results2[k],targets[k]))

    print 'mod04RMSE',mod04RMSE
    print 'dbmod04RMSE',dbmod04RMSE
    print 'nnrRMSE',nnrRMSE
    print 'nnrRMSE2',nnrRMSE2
    mod04RMSE   = np.mean(mod04RMSE)
    dbmod04RMSE = np.mean(dbmod04RMSE)
    nnrRMSE     = np.mean(nnrRMSE)
    nnrRMSE2     = np.mean(nnrRMSE2)

  eorig   = []
  edborig = []
  ecorr   = []
  ecorr2   = []
  for o,dbo,t,r,r2 in zip(original,dboriginal,targets,results,results2):
    eorig.append(o - t)
    edborig.append(dbo - t)
    ecorr.append(r - t)
    ecorr2.append(r2 - t)

  if emax is None:
    emax  = np.array([[e.max() for e in eorig],[e.max() for e in edborig],[e.max() for e in ecorr],[e.max() for e in ecorr2]]).max()
  if emin is None:
    emin  = np.array([[e.min() for e in eorig],[e.min() for e in edborig],[e.min() for e in ecorr],[e.min() for e in ecorr2]]).min()

  nbins        = 100
  corrected   = []
  corrected2   = []
  orig        = []
  dborig      = []
  x = np.linspace(emin,emax,nbins+1)
  xcen = x[:-1] + 0.5*(x[1:] - x[:-1])    
  if K is None:
    # cc, x   = np.histogram(ecorr[0],bins=np.linspace(emin,emax,nbins+1),density=True)
    # cc2, x   = np.histogram(ecorr2[0],bins=np.linspace(emin,emax,nbins+1),density=True)
    # oo, x   = np.histogram(eorig[0],bins=np.linspace(emin,emax,nbins+1),denstiry=True)
    # dboo, x = np.histogram(edborig[0],bins=np.linspace(emin,emax,nbins+1),denstiry=True)

    kernel = stats.gaussian_kde(ecorr[0])
    cc     = kernel(xcen)

    kernel = stats.gaussian_kde(ecorr2[0])
    cc2    = kernel(xcen)

    kernel = stats.gaussian_kde(eorig[0])
    oo     = kernel(xcen)

    kernel = stats.gaussian_kde(edborig[0])
    dboo   = kernel(xcen)

    corrected.append(cc)
    corrected2.append(cc2)
    orig.append(oo)
    dborig.append(dboo)
  else:
    for k,iTest in enumerate(I):
      # cc, x   = np.histogram(ecorr[k],bins=np.linspace(emin,emax,nbins+1),density=True)
      # cc2, x   = np.histogram(ecorr2[k],bins=np.linspace(emin,emax,nbins+1),density=True)
      # oo, x   = np.histogram(eorig[k],bins=np.linspace(emin,emax,nbins+1),density=True)
      # dboo, x = np.histogram(edborig[k],bins=np.linspace(emin,emax,nbins+1),density=True)

      kernel = stats.gaussian_kde(ecorr[k])
      cc     = kernel(xcen)

      kernel = stats.gaussian_kde(ecorr2[k])
      cc2    = kernel(xcen)

      kernel = stats.gaussian_kde(eorig[k])
      oo     = kernel(xcen)

      kernel = stats.gaussian_kde(edborig[k])
      dboo   = kernel(xcen)

      corrected.append(cc)
      corrected2.append(cc2)
      orig.append(oo)
      dborig.append(dboo)

  corrected   = np.array(corrected).T
  corrected2   = np.array(corrected2).T
  orig        = np.array(orig).T
  dborig      = np.array(dborig).T


  if K is not None:
    xcen = np.tile(xcen,(K,1)).T

  fig = plt.figure()
  ax  = plt.subplot(111)  
  ax.plot(xcen,orig,color='k')
  ax.plot(xcen,dborig,color='g')
  ax.plot(xcen,corrected,color='r')
  ax.plot(xcen,corrected2,color='c')  
  ax.set_xlim(emin,emax)
  if mxd.surface == 'ocean':
    ax.set_ylim(0,3.5)
  else:
    ax.set_ylim(0,2.0)
  orig_patch        = mpatches.Patch(color='k', label='MOD04 DT RMSE={:1.2F}'.format(mod04RMSE))
  dborig_patch      = mpatches.Patch(color='g', label='MOD04 DB RMSE={:1.2F}'.format(dbmod04RMSE))  
  corrected_patch   = mpatches.Patch(color='r', label='9-Ch NNR RMSE={:1.2F}'.format(nnrRMSE) )
  corrected_patch2   = mpatches.Patch(color='c', label='3-Ch NNR RMSE={:1.2F}'.format(nnrRMSE2) )  
  ax.legend(handles=[orig_patch,dborig_patch,corrected_patch,corrected_patch2])
  plt.grid(True, which='major',axis='both',color='0.50',linestyle='-')

  if Title is None:
    title("Error Log("+mxd.Target[0][1:]+"+0.01)")
  else:
    title(Title)
  savefig(outdir+"/error_int_pdf-"+expid+"."+ident+"-"+mxd.Target[0][1:]+'.png')  
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

    # len(regression) = nTarget
    # regression[*][0:2] = slope, intercept, r-value
    # out.shape = [ntestobs,nTarget]
    out, reg = mxd.test(iprint=False)

    reg = np.array(reg)
    mxd.nnr.slope[k,c,:]     = reg[:,0]
    mxd.nnr.intercept[k,c,:] = reg[:,1]
    mxd.nnr.R[k,c,:]         = reg[:,2]

    targets  = mxd.getTargets(mxd.iTest)

    # get other NNR STATS
    mxd.nnr.rmse[k,c,:] = rmse(out,targets)
    mxd.nnr.mae[k,c,:]  = mae(out,targets)
    mxd.nnr.me[k,c,:]   = me(out,targets)

    # loop through targets
    # get corresponding original MODIS retrieval
    for t in range(mxd.nTarget):
        name = 'm'+mxd.Target[t][1:]
        if name in mxd.__dict__:
            original = mxd.__dict__[name][mxd.iTest]
            if mxd.laod:
                original = log(original + 0.01)

            # get original MODIS stats 
            lm = LinearRegression()
            tt = targets[:,t]
            tt.shape = tt.shape + (1,)
            lm.fit(tt,original)
            mxd.orig.slope[k,c,t]     = lm.coef_[0]
            mxd.orig.intercept[k,c,t] = lm.intercept_
            mxd.orig.R[k,c,t]         = np.sqrt(lm.score(tt,original))

            tt  = tt.squeeze()
            mxd.orig.rmse[k,c,t] = rmse(original,tt)
            mxd.orig.mae[k,c,t]  = mae(original,tt)
            mxd.orig.me[k,c,t]   = me(original,tt)

    
# ---
def rmse(predictions, targets):
    return np.sqrt((np.square(predictions - targets)).mean(axis=0))
# ---
def mae(predictions, targets):
    return np.abs(predictions-targets).mean(axis=0)
# ---
def me(predictions, targets):
    return (predictions-targets).mean(axis=0)    

# ----
def MAKE_ERROR_PDFS(mxdx,Input,expid,ident,K=None,I=None,Title=None,
                    netfileRoot=None,netfileRoot2=None,Input2=None,
                    emin=-1.5,emax=2.5,doInt=False,mxdx2=None):

  if mxdx2 is None:
    if not doInt:
      make_error_pdfs(mxdx,Input,expid,ident,K=K,I=I,Title=Title,netfileRoot=netfileRoot,
                      emin=emin,emax=emax)
    else:
      make_error_pdfs_int(mxdx,Input,expid,ident,K=K,I=I,Title=Title,netfileRoot=netfileRoot,
                      emin=emin,emax=emax)    
  else:
    make_error_pdfs_dbdt(mxdx,mxdx2,Input,expid,ident,K=K,I=I,Title=Title,
                    netfileRoot=netfileRoot,netfileRoot2=netfileRoot2,Input2=Input2,
                    emin=emin,emax=emax)        

#---------------------------------------------------------------------
def SummaryPDFs(mxdx,mxdx2=None,varnames=['mRef870','mSre470'],doInt=False):
    K = mxdx.K
    Irange     = arange(mxdx.nobs)
    iValid     = Irange[mxdx.iValid]
    if K is None:
      I = [Irange]
    else:
      I = []
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
      if len(mxdx.comblist) == 1:
        invars = mxdx.comblist[0]
        netfileRoot = mxdx.outdir+"/"+'.'.join(invars)
      else:
        for invars in itertools.permutations(Input):
          netfileRoot = mxdx.outdir+"/"+'.'.join(invars)
          filelist = glob(netfileRoot+'*.net')
          if len(filelist) > 0:
            Input = invars
            break

        if len(filelist) == 0:
          raise Exception('{} not found.  Need to train this combinatin of inputs'.format(netfileRoot+'*.net'))

      # Only works if combinations are in the same order
      Input2 = None
      netfileRoot2 = None
      if mxdx2 is not None:
        if len(mxdx.comblist) == 1:
          invars = mxdx2.comblist[0]
          netfileRoot2 = mxdx2.outdir+"/"+'.'.join(invars)
          Input2 = mxdx2.comblist[0]
        else:        
          Input2 = mxdx2.comblist[c]
          for invars in itertools.permutations(Input2):
            netfileRoot2 = mxdx2.outdir+"/"+'.'.join(invars)
            filelist = glob(netfileRoot2+'*.net')
            if len(filelist) > 0:
              Input2 = invars
              break

          if len(filelist) == 0:
            raise Exception('{} not found.  Need to train this combinatin of inputs'.format(netfileRoot2+'*.net'))


      MAKE_ERROR_PDFS(mxdx,Input,'AllTest',mxdx.ident,K=K,I=I,
                     netfileRoot=netfileRoot,
                     netfileRoot2=netfileRoot2,
                     Input2=Input2,
                     doInt=doInt,
                     mxdx2=mxdx2)

      MAKE_ERROR_PDFS(mxdx,Input,'q25',mxdx.ident,K=K,I=I1,
                     Title="Q1 Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot,
                     netfileRoot2=netfileRoot2,
                     Input2=Input2,
                     doInt=doInt,
                     mxdx2=mxdx2)
      MAKE_ERROR_PDFS(mxdx,Input,'q50',mxdx.ident,K=K,I=I2,
                      Title="Q2 Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot,
                     netfileRoot2=netfileRoot2,
                     Input2=Input2,
                     doInt=doInt,
                     mxdx2=mxdx2)
      MAKE_ERROR_PDFS(mxdx,Input,'q75',mxdx.ident,K=K,I=I3,
                      Title="Q3 Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot,
                     netfileRoot2=netfileRoot2,
                     Input2=Input2,                     
                     doInt=doInt,
                     mxdx2=mxdx2)
      MAKE_ERROR_PDFS(mxdx,Input,'q100',mxdx.ident,K=K,I=I4,
                      Title="Q4 Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot,
                     netfileRoot2=netfileRoot2,
                     Input2=Input2,
                     doInt=doInt,
                     mxdx2=mxdx2)

      MAKE_ERROR_PDFS(mxdx,Input,'fdu',mxdx.ident,K=K,I=Ifdu,
                     Title="Dust Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot,
                     netfileRoot2=netfileRoot2,
                     Input2=Input2,
                     doInt=doInt,
                     mxdx2=mxdx2)
      MAKE_ERROR_PDFS(mxdx,Input,'fss',mxdx.ident,K=K,I=Ifss,
                      Title="Sea Salt Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot,
                     netfileRoot2=netfileRoot2,
                     Input2=Input2,
                     doInt=doInt,
                     mxdx2=mxdx2)
      MAKE_ERROR_PDFS(mxdx,Input,'fcc',mxdx.ident,K=K,I=Ifcc,
                      Title="BC+OC Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot,
                     netfileRoot2=netfileRoot2,
                     Input2=Input2,
                     doInt=doInt,
                     mxdx2=mxdx2)
      MAKE_ERROR_PDFS(mxdx,Input,'fsu',mxdx.ident,K=K,I=Ifsu,
                      Title="Sulfate Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot,
                     netfileRoot2=netfileRoot2,
                     Input2=Input2,
                     doInt=doInt,
                     mxdx2=mxdx2)   

      MAKE_ERROR_PDFS(mxdx,Input,'fna',mxdx.ident,K=K,I=Ifna,
                      Title="No Dominant Species Error Log("+mxdx.Target[0][1:]+"+0.01)",
                     netfileRoot=netfileRoot,
                     netfileRoot2=netfileRoot2,
                     Input2=Input2,
                     doInt=doInt,
                     mxdx2=mxdx2)                            

      if varnames is not None:
        for varname in varnames:

          I1, I2, I3, I4 = get_Iquartiles(mxdx,I=I,var=varname)

          # for c,Input in enumerate(mxdx.comblist):
          #   netfileRoot = mxdx.outdir+"/"+'.'.join(Input)
          MAKE_ERROR_PDFS(mxdx,Input,'{}Q1'.format(varname),mxdx.ident,K=K,I=I1,
                         Title="Q1 {} Error Log({}+0.01)".format(varname,mxdx.Target[0][1:]),
                         netfileRoot=netfileRoot,
                         netfileRoot2=netfileRoot2,
                         Input2=Input2,                       
                         doInt=doInt,
                         mxdx2=mxdx2)

          MAKE_ERROR_PDFS(mxdx,Input,'{}Q2'.format(varname),mxdx.ident,K=K,I=I2,
                         Title="Q2 {} Error Log({}+0.01)".format(varname,mxdx.Target[0][1:]),
                         netfileRoot=netfileRoot,
                         netfileRoot2=netfileRoot2,
                         Input2=Input2,     
                         doInt=doInt,
                         mxdx2=mxdx2)      

          MAKE_ERROR_PDFS(mxdx,Input,'{}Q3'.format(varname),mxdx.ident,K=K,I=I3,
                         Title="Q3 {} Error Log({}+0.01)".format(varname,mxdx.Target[0][1:]),
                         netfileRoot=netfileRoot,
                         netfileRoot2=netfileRoot2,
                         Input2=Input2,     
                         doInt=doInt,
                         mxdx2=mxdx2)            

          MAKE_ERROR_PDFS(mxdx,Input,'{}Q4'.format(varname),mxdx.ident,K=K,I=I4,
                         Title="Q4 {} Error Log({}+0.01)".format(varname,mxdx.Target[0][1:]),
                         netfileRoot=netfileRoot,
                         netfileRoot2=netfileRoot2,
                         Input2=Input2,     
                         doInt=doInt,
                         mxdx2=mxdx2)                    
