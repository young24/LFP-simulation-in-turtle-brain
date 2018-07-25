# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 08:05:02 2016
plot potentials on MEA with a 11x11 figure.
@author: young
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

def find_extremum(listName):
    minValue = np.min(listName)
    maxValue = np.max(listName)
    if (abs(minValue)>abs(maxValue)):
        return minValue
    else:
        return maxValue
if __name__ == '__main__':
    # Four axes, returned as a 2-d array
    t = np.linspace(0,50,401)
    LFPonMEA = np.load('simulation_results/V6N2_L2_80_L1_230/LFPonMEA.npy')
    LFPonMEA = -LFPonMEA # inhibitory neuron
    maxAmp = find_extremum(LFPonMEA)
    f, ax = plt.subplots(11,11,sharex=True,sharey=True)
    for x in range(0,11):
        for y in range(0,11):          
            ax[y, x].plot(t,LFPonMEA[x,y,:])
            ax[y, x].axis('off') 
            if x==0 and y==0:
                ax[y, x].axis('on')
                ax[y, x].spines['top'].set_visible(False)
                ax[y, x].spines['right'].set_visible(False)
                ax[y, x].get_xaxis().tick_bottom()
                ax[y, x].get_yaxis().tick_left()
                # Find at most 101 ticks on the y-axis at 'nice' locations
                ax[y, x].yaxis.set_ticks(np.arange(maxAmp,maxAmp+1,1))
                #ax[y, x].ylabel('$\mu$V')
                ax[y, x].xaxis.set_visible(False)
                

                
    plt.show()
    