#
# Scatter Plots - preliminary


import matplotlib               # import all the plotting matplotlib stuff
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as cols

###########
# General Scatter Plotter
###########
def mypredplot(pnum,xvals,yvals,title,outfile,xylim):
    plt.figure(num=pnum)
    plt.hold(True)
    plt.plot(xvals,yvals,'ko',markerfacecolor='w')
    plt.xlim(xylim)
    plt.ylim(xylim)
    plt.plot([min(xylim),max(xylim)],[min(xylim),max(xylim)],'k')
    plt.xlabel('Measured')
    plt.ylabel('Predicted')
    plt.title(title)
    plt.hold(False)
    plt.savefig(outfile)
