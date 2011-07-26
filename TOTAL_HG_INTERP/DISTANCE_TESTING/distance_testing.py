import numpy as np
from haversine import distance as distH

indat = np.genfromtxt('calib_data_all_forDPK_wHUCs_DD_alldataresults.csv',
                      delimiter=',',
                      names=True,
                      dtype=None)
# degrees style
DD_X = np.atleast_2d(indat['X_DD_IN']).T
DD_Y = np.atleast_2d(indat['Y_DD_IN']).T

ref_DD = np.tile((DD_X[0][0],DD_Y[0][0]),[len(DD_X),1])

all_DD = np.hstack((DD_X,DD_Y))

d_DD = distH(ref_DD,all_DD)

# meters style
M_X = np.atleast_2d(indat['X_METER_IN']).T/1000.0
M_Y = np.atleast_2d(indat['Y_METER_IN']).T/1000.0

ref_M = np.tile((M_X[0][0],M_Y[0][0]),[len(M_X),1])

all_M = np.hstack((M_X,M_Y))

d_M = np.sqrt((ref_M[:,0]-all_M[:,0])**2+((ref_M[:,1]-all_M[:,1]))**2)

RPD = np.abs(2*(d_M-d_DD)/(d_M+d_DD)) * 100
ofp = open('compare_dist.out','w')
ofp.write('%16s%16s%16s\n' %('DD_STYLE','METER_STYLE','RPD'))
for i in np.arange(len(RPD)):
    ofp.write('%16.8e%16.8e%15.2f%%\n' %(d_DD[i],d_M[i],RPD[i]))
ofp.close()