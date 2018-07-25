# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 08:05:02 2016
plot potentials on MEA with a 11x11 figure.
@author: young
"""

import numpy as np
import pylab as plt


if __name__ == '__main__':
    # Four axes, returned as a 2-d array
    t = np.linspace(0,50,401)
    phi = np.load('../sim_results/phi.npy')
    f, axarr = plt.subplots(11,11,sharex=True,sharey=True)
    tstart = 0
    tend = 80 
    
    
    for x in range(0,11):
        for y in range(0,11):
            #phiData = phi[t]
            axarr[y, x].plot(t,phi[11*x+y],color = 'black')
            axarr[y, x].axis('off')
            
    plt.show()
            