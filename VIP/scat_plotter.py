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
def scat_plot(pnum,xvals,yvals,title,xname,yname,outfile,showme):
    plt.figure(num=pnum)
    plt.plot(xvals,yvals,'bd')
    plt.xlabel(xname)
    plt.ylabel(yname)
    plt.title(title)
    plt.savefig(outfile)
    if showme:
        plt.show()