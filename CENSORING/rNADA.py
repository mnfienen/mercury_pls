#/usr/bin/python

# Code to use the R package NADA from Helsel's textbook (2005)
# Data must be in a CSV names #inroots[i]#.csv
# One column is a unique ID with the name in variable called UNID
# Another column with the data should be named #inroots[i]#
# Final column called CEN has the flags for censoring TRUE is censored, FALSE is not
#
# N.B. --> the values should be sorted smallest to largest by the data 
# a m!ke@usgs joint
#
#

import numpy as np
import rpy2
import rpy2.robjects as rob
import matplotlib.mlab as mlab
import matplotlib.pyplot as pyplot
import rpy2.robjects.numpy2ri

# root names for all the files containing data to estimate values for NDs for 
inroots = ['D_MEHG','P_MEHG','UNF_MEHG','P_HG','SULFATE']
UNID = 'QW_SEQ'
r = rob.r
r.library('NADA')

for cr in inroots:
    indat = np.genfromtxt(cr + '.csv',names=True,dtype=None,delimiter=',')
    sampID = indat[UNID]
    
    vals = r['as.vector'](indat[cr])
    cen  = r['as.vector'](indat['CEN'].astype('bool'))
    rob.globalenv['vals'] = vals
    rob.globalenv['cen']  = cen
    
    # run the censoring fill-in algorithm 
    modros = r('mod <- ros(vals,cen)')
    
    # make a numpy array of the results
    outmodros = np.array(r['as.data.frame'](modros))
    
    ofp = open(cr + '_out.dat','w')
    ofp.write('{0:s} {1:s} {2:s} {3:s}\n'.format(UNID.rjust(12),
                                                 'OBS'.rjust(12),
                                                 'CENSORED'.rjust(12),
                                                 'MODELED'.rjust(12)))
    for i in np.arange(len(sampID)):
        ofp.write('{0:12d} {1:12.6f} {2:s} {3:12.6f}\n'.format(sampID[i],
                                                              outmodros[0,i],
                                                              str(outmodros[1,i].astype('bool')).rjust(12),
                                                              outmodros[3,i]))
    ofp.close()