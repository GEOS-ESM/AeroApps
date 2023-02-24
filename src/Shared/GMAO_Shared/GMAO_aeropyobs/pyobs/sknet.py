"""

  Simple wrapper around FFNET to make it look like a sklearn estimator.
  It also provides sklearn style cross_validation

"""

from ffnet                    import *
from sklearn.base             import BaseEstimator
from sklearn.model_selection import cross_val_score

class SKNET(ffnet,BaseEstimator):

    def predict(self,X):
        """
        Evaluate the Feed-Forward Neural Network.
        """        
        y_ = self(X)
        return y_.ravel()

    def fit(self,X,y,**kwopts):
        """
        Train the Feed-Forward Neural Net.
        """
        self.train_tnc(X,y,**kwopts)

    def score(self,X,y):
        """
        Returns the coefficient of determination R^2 of the prediction.
       
        The coefficient R^2 is defined as (1 - u/v), where u is the
        regression sum of squares ((y - y_pred) ** 2).sum() and v is the
        residual sum of squares ((y_true - y_true.mean()) ** 2).sum().
        Best possible score is 1.0, lower values are worse.
        """
        y_ = self.predict(X)
        u = ((y_-y)**2).sum()
        v = ((y-y.mean())**2).sum()
        return (1 - u/v)

    def cross_validate(self,X,y,**kwopts):
        """
        Return cross-validation scores. See cross_val_score() for
        optional parameters.
            scores = net.cross_validate(X,y,cv=10)
        """
        scores = cross_val_score(self, X, y, **kwopts)
        return scores        
