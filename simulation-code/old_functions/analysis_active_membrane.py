# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 10:44:07 2016
calc_half_width
@author: superuser
input:
1. imem or vmem
2. tvec
"""
import numpy as np
import pylab as plt

def get_delay_of_peak(imem):
    startTime = 10. # ms
    factor = 0.125 # 1 step corresponds to 0.125 ms in imem
    idxMax = np.argmax(np.abs(imem),1)
    delay = idxMax * factor
    delay = delay - startTime
    return delay
    
def plot_weighted_morphology(pos, weight, whichView = 'sideView', title = 'weighted morphology'):
    cm = plt.cm.get_cmap('RdYlBu')
    #f,ax = plt.subplots(2,sharex=True)
    if whichView == 'sideView':
        x = pos[:,0]
        y = pos[:,2]
    elif whichView == 'topView':
        x = pos[:,0]
        y = pos[:,1]
    else:
        errmsg = 'please enter sideView/topView for the argument"whichView"'
        raise Exception(errmsg)
              
    sc = plt.scatter(x,y,c=weight,cmap=cm,lw=0)
    plt.xlabel('x [$\mu m$]')
    plt.ylabel('y [$\mu m$]')
    plt.title(title)
    plt.colorbar(sc)
    plt.show()
    
def interpolate_delay(order,pos,z,delay):
    f=(pos-z[order])/(z[order+1]-z[order])
    delayValue = delay[order]*f+delay[order+1]*(1-f)
    return delayValue
    

if __name__ == '__main__':
    folderPath = 'sim_results'
    imem = np.load(folderPath+'/imem.npy')
    imem = np.abs(imem)
    max_current = np.max(imem, axis=1)
    halfWidth = []
    for i in range(0,len(max_current)):
        halfWidth.append(len(np.argwhere(imem[i]>max_current[i]/2)))
    halfWidth = np.array(halfWidth)
    np.save(folderPath+'/halfWidth.npy',halfWidth)