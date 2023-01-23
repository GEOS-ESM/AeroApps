"""
   Generic Neural Net Functionality.

   Arlindo da Silva, June 2015.

"""

import os
import pyobs.sknet as nn

from   matplotlib.pyplot import  cm, imshow, plot, figure
from   matplotlib.pyplot import  xlabel, ylabel, title, grid, savefig, legend
from   numpy             import  c_ as cat
from   numpy             import  random, sort, pi, load, cos, log, std, exp
from   numpy             import  reshape, arange, ones, zeros, interp
from   numpy             import  meshgrid, concatenate
import numpy             as      np
from   matplotlib        import  ticker
from   scipy             import  stats, mgrid
from   sklearn.model_selection import KFold
#..............................................................
class aodFormat(ticker.Formatter):
    def __call__(self,x,pos=None):
        y = exp(x)-0.01
        return '%4.2f'%y
#..............................................................

class NN(object):

    def train (self,Input=None,Target=None,nHidden=200,maxfun=2550,biases=True,
               topology=None,bounds=[-100,100], **kwargs):
        """
        Train the Neural Net, using a maximum of *maxfun* iterations.
        On input,
            Input   ---  string list with the name of the predictors;
                         if specified dataset specific default is
                         redefined.
            Target  ---  string list with the name of the targets;
                         if specified dataset specific default is
                         redefined.
           nHidden  ---  number of hidden nodes
           maxfun   ---  max number of iterations
           biases   ---  whether to include bias nodes
         topology   ---  Network topology; default is (nInput,nHidden,nTarget)
         
         Returns:
            Nothing.
        """
            
        # Possibly redefine Input/Targets
        # -------------------------------
        if Input != None:
            self.Input = Input
        if Target != None:
            self.Target = Target

        # Instantiate Neural Net
        # ----------------------
        if topology==None:
            topology = (len(self.Input), nHidden,len(self.Target))
        #self.net = nn.ffnet(nn.mlgraph(topology,biases=biases))
        self.net = nn.SKNET(nn.mlgraph(topology,biases=biases))

        # Add these attributes to net so that later on
        # we now how to apply it to regular MODIS data
        # --------------------------------------------
        self.net.InputNames = self.Input
        self.net.TargetNames = self.Target
        self.net.laod = self.laod
        if self.surface == 'ocean':
            self.net.Wind = self.Wind

        # Indices for training set
        # ------------------------
        try:
            iTrain = self.iTrain
        except AttributeError:
            iTrain = self.iValid # good QC marks
            
        # Prepare inputs and targets
        # --------------------------
        inputs  = self.getInputs(iTrain)
        targets = self.getTargets(iTrain) 

        # Train
        # -----
        bounds = [bounds]*self.net.get_params()['conec'].shape[0]
        if self.verbose>0:
            print "Starting training with %s inputs and %s targets"\
                  %(str(inputs.shape),str(targets.shape))
        self.net.train_tnc(inputs,targets, maxfun=maxfun,bounds=bounds,**kwargs)
#        self.net.train_bfgs(inputs,targets, maxfun=maxfun)


    def test(self,iprint=1,fname=None):

        # Indices for training set
        # ------------------------
        try:
            iTest = self.iTest
        except AttributeError:
            iTest = self.iValid
            
        # Prepare inputs and targets
        # --------------------------
        inputs  = self.getInputs(iTest)
        targets = self.getTargets(iTest) 

        return self.net.test(inputs,targets,iprint=iprint,filename=fname)
        
    def eval(self,I=None):
        if I is None: I = self.iValid
        inputs = self.getInputs(I)
        return self.net(inputs)

    __call__ = eval

    def derivative(self,I=None):
        if I is None: I = self.iValid
        return self.net.derivative(self.getInputs(I))
    
    def savenet(self,fname):
        nn.savenet(self.net,fname)

    def loadnet(self,fname):
        return nn.loadnet(fname)    

    def exportnet(self,fname):
        nn.exportnet(self.net,fname)
        
    def split (self,fTrain=0.9):
        """
        Splits the input dataset in training and testing subsets. No data is
        actually moved only attributes with indices iTrain/iTest are created;
        only data with an iValid Q/C flag is considered. On input, *fTrain* is
        the fraction of the dataset to be used for training.
        Returns: (nothing)
        """
        n = self.lon.size
        nTrain = int(fTrain * n)
        random.seed(32768) # so that we get the same permutation
        i = random.permutation(n)
        iValid = self.iValid[i]
        self.iTrain = i[0:nTrain][iValid[0:nTrain]] # Keep only good obs
        self.iTest  = i[nTrain:][iValid[nTrain:]]   # Keep only good obs

    def kfold(self,K=3):
        """
        Splits input dataset into K training and testing subsets.
        Only data with an iValid Q/C flag is considered.
        """
        n = self.lon.size
        self.kf = KFold(np.sum(self.iValid), n_folds=K, shuffle=True, random_state=n)


    def getInputs(self,I,Input=None):
        """
        Given a set of indices *I*, returns the corresponding
        inputs for a neural net evaluation.
        Returns: inputs
        """
        if self.verbose:
            print " "
            print "       Feature          Min      Max"
            print "  ------------------  -------  -------"
        if Input==None:
            Input = self.Input
        inputs = self.__dict__[Input[0]][I]
        if self.verbose:
            print "%20s %8.4f %8.4f"%(Input[0],inputs.min(),inputs.max())
        for var in Input[1:]:
            q = self.__dict__[var][I]
            inputs = cat[inputs,q]
            if self.verbose:
                print "%20s %8.4f %8.4f"%(var,q.min(),q.max())
        if self.verbose:
            print "  ------------------  -------  -------"
            print ""

        if len(inputs.shape) == 1:
            inputs.shape = (inputs.shape[0],1)            
        return inputs
    
    def getTargets(self,I):
        """
        Given a set of indices *I*, return the corresponding
        targets for a neural net evaluation:
        Returns: tagets
        """
        var = self.Target[0]
        if self.laod and ('Tau' in var):
            targets = log(self.__dict__[var][I] + 0.01)
        else:
            targets = self.__dict__[var][I]

        for var in self.Target[1:]:
            if self.laod and ('Tau' in var):
                targets = cat[targets,log(self.__dict__[var][I] + 0.01)]
            else:
                targets = cat[targets,self.__dict__[var][I]]

        return targets
 
    def plotKDE(self,bins=None,I=None,figfile=None,
                x_label='AERONET'):
        """
        Plot Target vs Model using a 2D Kernel Density Estime.
        """
        if I is None: I = self.iValid # All data by default
        results = self.eval(I)
        targets = self.getTargets(I)
        if self.laod:
            formatter = aodFormat()
        else:
            formatter = None
        if bins == None:
            if self.laod:
                bins = arange(-5., 1., 0.1 )
            else:
                bins = arange(0., 0.6, 0.01 )
        x_bins = bins
        y_bins = bins
        if len(targets.shape) == 1:
            x_values = targets
            y_values = results.squeeze()
        else:
            x_values = targets[:,0]            
            y_values = results[:,0]
        _plotKDE(x_values,y_values,x_bins,y_bins,y_label='NNR',
                 formatter=formatter,x_label=x_label)        
        title("Log("+self.Target[0][1:]+"+0.01) - "+self.ident)
        if figfile != None:
            savefig(figfile)
            
    def plotScat(self,iTarget=0,bins=None,I=None,figfile=None):
        """
        Plot Target vs Model as a scatter plot
        """
        if I is None: I = self.iTest # Testing data by default
        results = self.eval(I)[:,iTarget]
        targets = self.getTargets(I)
        if self.nTarget > 1:
            targets = targets[:,iTarget]
        original = log(self.__dict__['m'+self.Target[iTarget][1:]][I] + 0.01)
        if bins == None:
            bins = arange(-5., 1., 0.1 )

        figure()
        plot(targets,original,'bo',label='Original')
        plot(targets,results,'ro',label='Corrected')
        legend(loc='upper left')
        plot(bins,bins,'k')
        grid()
        xlabel('AERONET')
        ylabel('MODIS')
        title("Log("+self.Target[iTarget][1:]+"+0.01) - "+self.ident)
        if figfile != None:
            savefig(figfile)
#---------------------------------------------------------------------------------
def _cat2 (X, Y):
    """
    Given 2 arrays of same shape, returns array of shape (2,N),
    where N = X.size = Y.size
    """
    xy = concatenate((X.ravel(),Y.ravel())) # shape is (N+N)
    return reshape(xy,(2,X.size))         # shape is (2,N)
    
def _plotKDE(x_values,y_values,x_bins=None,y_bins=None,
             x_label='AERONET', y_label='MODIS',formatter=None):
        """
        Plot Target vs Model using a 2D Kernel Density Estimate.
        """

        if x_bins is None: x_bins = arange(-5., 1., 0.1 )
        if y_bins is None: y_bins = x_bins

        Nx = len(x_bins)
        Ny = len(y_bins)

        print "Evaluating 2D kernel on grid with (Nx,Ny)=(%d,%d) ..."%(Nx,Ny)
        kernel = stats.kde.gaussian_kde(_cat2(x_values,y_values))
        X, Y = meshgrid(x_bins,y_bins)   # each has shape (Ny,Nx)
        Z = kernel(_cat2(X,Y))           # shape is (Ny*Nx)
        Z = reshape(Z,X.shape)

        fig = figure()
        # ax = fig.add_axes([0.1,0.1,0.75,0.75])
        ax = fig.add_axes([0.1,0.1,0.75,0.75])
        if formatter != None:
            ax.xaxis.set_major_formatter(formatter)
            ax.yaxis.set_major_formatter(formatter)
        imshow(Z, cmap=cm.gist_earth_r, origin='lower', 
               extent=(x_bins[0],x_bins[-1],y_bins[0],y_bins[-1]) )
        plot([x_bins[0],x_bins[-1]],[y_bins[0],y_bins[-1]],'k')
        xlabel(x_label)
        ylabel(y_label)
        grid()

        return fig
