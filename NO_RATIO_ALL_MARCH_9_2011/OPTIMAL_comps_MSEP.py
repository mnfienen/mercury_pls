'''
# function to perform bootstrap MSEP calculation for
# optimization of number of components
'''

import numpy as np
import pickle
import random
import matplotlib.pyplot as plt
import os


# get the current directory name
currd = os.getcwd().split('/')[-1]


infiles = ['cv_data.pkl','hg_data.pkl']

invars = list()

for fn in infiles:
    ifp = open(fn,'rb')
    invars.append(pickle.load(ifp))

cv_data = invars[0].squeeze()
hg_data = invars[1].squeeze()

# cv data is the array of predicted cv results (for hg) 
# now calcualte the predicted sum of squared errors

num_comps_all = cv_data.shape[1]
num_samples   = cv_data.shape[0]

PRESS = list()
for i in np.arange(num_comps_all):
    tmpv = cv_data[:,i] - hg_data
    PRESS.append(np.dot(tmpv,tmpv))
    
# now, find the number of components with lowest value of PRESS
PRESS_opt = np.argmin(PRESS)

# cache all squared errors to randomly pull from below
cv_squared_error = (cv_data[:,PRESS_opt] - hg_data)**2

# now, we use an iterative, resamping with replacement approach to calcualte the 
# standard deviation of PRESS
PRESS_std = list()
num_outer = 1000
num_inner = 100
randomfun = random.random

for outer in np.arange(num_outer):
    print str(outer+1) + 'th outer iteration of ' + str(num_outer)
    PRESS_bootstrap = list()
    for inner in np.arange(num_inner):        
        randinds = np.zeros((num_samples,1)).astype(int)
        for i in np.arange(num_samples):
            randinds[i] = round(randomfun()*(num_samples-1))
        PRESS_bootstrap.append(np.sum(cv_squared_error[randinds]))
        PRESS_std.append(np.std(PRESS_bootstrap))
med_stddev = np.median(PRESS_std)

fig = plt.figure()
# make error bars vectors
errs = np.zeros_like(PRESS)
errs[PRESS_opt] = med_stddev
#plot results
plt.hold = True
plt.plot([1,num_comps_all+1],[PRESS[PRESS_opt]+med_stddev,PRESS[PRESS_opt]+med_stddev],':')
plt.errorbar(np.arange(num_comps_all)+1,PRESS,marker='.',c='b',ms=12.0,yerr = errs)
plt.xlim((0,num_comps_all+1))
plt.xlabel('number of components used')
plt.ylabel('Prediction Sum of Squared Errors')

plt.savefig(currd + '_PRESS_figure.eps')


# print out some results
ofp = open(currd + '_Components_optimization.dat','w')
ofp.write('\n'*2 + '*'*20 + '\n'*2)
ofp.write('output from OPTIMAL_comps_MSEP.py\n')
ofp.write('Number of outer iterations --> %d\n' %(num_outer))
ofp.write('Number of inner iterations --> %d\n' %(num_inner))

ofp.write('PRESS Values\n')
for i,cp in enumerate(PRESS):
    ofp.write('%d components%15.4f\n' %(i+1,cp))
ofp.write('*'*20 + '\n')
ofp.write('Mean STD of PRESS --> %10.7f\n' %(med_stddev))

bestindy = np.nonzero(PRESS < min(PRESS)+med_stddev)[0]
numopt = np.min(bestindy) + 1
ofp.write('Number of optimal components --> %d\n' %(numopt))
ofp.write('\nFor details see ' + currd + '_PRESS_figure.eps\n')
ofp.write('*'*20 + '\n'*2)
ofp.close()
