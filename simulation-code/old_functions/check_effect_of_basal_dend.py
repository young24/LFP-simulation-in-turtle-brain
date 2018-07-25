# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 10:18:19 2016
check the effect of basal dendrites
@author: superuser
"""


import numpy as np
import pylab as plt

def find_extremum(arr):
    if np.abs(np.min(arr)) > np.abs(np.max(arr)):
        extremum = np.min(arr)
    else: 
        extremum = np.max(arr)
    return extremum
        
    

if __name__ == '__main__':
    # Four axes, returned as a 2-d array
    t = np.load('plots/check_dend_effect/check_data/all_dendrites/tvec.npy')
    phi_test = np.load('plots/check_dend_effect/check_data/delete_subTree/phi.npy')
    phi_refer = np.load('plots/check_dend_effect/check_data/all_dendrites/phi.npy')
    #f, axarr = plt.subplots(11,11,sharex=True,sharey=True)
    #tstart = 0
    #tend = 80 
    #t = t[tstart:tend]
#    plt.tick_params(
#        axis='x',          # changes apply to the x-axis
#        which='both',      # both major and minor ticks are affected
#        bottom='off',      # ticks along the bottom edge are off
#        top='off',         # ticks along the top edge are off
#        labelbottom='off') # labels along the bottom edge are off
    ratio = np.zeros((11,11))
    for x in range(0,11):
        for y in range(0,11):
            #phiData = phi[t]
            extremum_refer = find_extremum(phi_refer[11*x+y])
            extremum_test = find_extremum(phi_test[11*x+y])          
            ratio[x,y] = extremum_test/extremum_refer
            #axarr[y, x].plot(t,ratio[11*x+y])
            #axarr[y, x].axis('off')
    np.save('ratio.npy',ratio)   
    #plt.show()


