'''
# function to perform bootstrap MSEP calculation for
# optimization of number of components
'''

import numpy as np
import pickle
import random
import matplotlib.pyplot as plt
import os


ifp = open('LOADINGS.pkl','rb')
loads = (pickle.load(ifp))
nr,nc = loads.shape
inds = np.arange(nc)

plt.figure()
plt.hold = True
width = 0.2
colors = ['r','g','b','k']
complabs = list()
bars = list()
for i in inds:
    complabs.append('Comp' + str(i+1))
    pltinds = inds + width*i
    btmp = (plt.bar(pltinds,loads[i,:],width,color=colors[i]))
    bars.append(btmp[0])
plt.legend(bars,('SULFATE','TOC','pH','Wetland %'))
plt.xticks(inds+width*2,complabs)
plt.ylabel('Loadings')
plt.savefig('loadings.png')
plt.savefig('loadings.eps')
    