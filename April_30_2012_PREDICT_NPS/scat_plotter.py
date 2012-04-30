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
def scat_plot(pnum,xvals,yvals,title,xname,outfile,RAT,showme):
    plt.figure(num=pnum)
    plt.hold(True)
    plt.plot(xvals,yvals,'bd')
    plt.xlabel(xname)
    if RAT:
        plt.ylabel('Methyl Hg/Total Hg')
    else:
        plt.ylabel('Methyl Hg')
    plt.title(title)
    plt.hold(False)
    plt.savefig(outfile)
    if showme:
        plt.show()