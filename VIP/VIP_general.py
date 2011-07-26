#VIP Calculations for PLSR

import numpy as np
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt


def all_scatters(X,names,casename):
    '''
    scatter plots of all variables against one another in the X matrix
    '''
    from scat_plotter import scat_plot
    nv = len(names)
    k=1
    for i in np.arange(nv):
        for j in np.arange(nv):
            k+=1
            scat_plot(k,X[:,i],X[:,j],names[i] + ' ' + names[j],names[i],names[j],
                      casename + "_" + names[i] + "_" + names[j] + '.png',False)

def r_cor(X,y):
    '''
    emulation of the Pearsons correlation function from R
    X is a n x p matrix where n is the number of measurements
                              and p is the number of variables
    y is an m x 1 vector
    '''
    r = []
    for i in X.T:
        r.append(pearsonr(i,y)[0])
    r = np.array(r)
    r = np.atleast_2d(r)
    return r

def auto_cov(X):
    '''
    calculate the covariance of X where X is an n x p matrix and
    p is the number or variables while n is the number of measurements
    '''
    p = X.shape[1]
    # subtract off the means of X and divide by sd
    for i in np.arange(p):
        cmean = np.mean(X[:,i])
        csd   = np.std(X[:,i])
        X[:,i] = (X[:,i] - cmean) / csd
    CV = np.dot(X.T,X)
    return CV

def VIP_calc(X,y,W):
    '''
    X is the coefficient matrix
    y is the data (response variable) vector
    W is the loadings matrix - from PLSR or PLS
    '''
    # define dimensions
    p = X.shape[1] # number of variates
    H = W.shape[1] # number of components
    # calculate Pearson's R correlation
    cor2 = r_cor(X,y)**2
    

    # complete VIP calculation
    VIP = np.zeros((p,H))
    VIP[:,0] = W[:,0]**2
    Rd = []
    if H>1:
        for h in np.arange(1,H+1):
            Rd = cor2[:,0:h]
            VIP[:,h-1] = np.dot(Rd,np.transpose(W[:,0:h]*W[:,0:h]))/np.sum(Rd)
    VIP = np.sqrt(p * VIP)
    return VIP