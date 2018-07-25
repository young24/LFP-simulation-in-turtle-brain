# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 14:04:30 2016
find min and max in a list of unqual length arrays
@author: young
"""

import numpy as np
import pylab as plt


def find_extremum_in_list(listName,idx):
    minValue = np.min(listName[idx])
    maxValue = np.max(listName[idx])
    for i in range(1,len(listName)):
        if np.min(listName[i]) < minValue:
            minValue = np.min(listName[i])
        if np.max(listName[i]) > maxValue:
            maxValue = np.max(listName[i])
    if (abs(minValue)>abs(maxValue)):
        return minValue
    else:
        return maxValue
    
    
if __name__ == '__main__':
    dirName = 'averaged_result_for_real_neuron'
    idx = 0
    extremeList = []
    
    for z in range(20,301,10):
        fileName = dirName+'/'+'phi_z'+str(z)+'.npy'
        phi = np.load(fileName)
        extremum = find_extremum_in_list(phi,idx)       
        extremeList.append(extremum)
    
    # normalization
#    minValue = min(extremeList)
#    maxValue = max(extremeList)
#    extremeList = [(n-minValue)/(maxValue-minValue) for n in extremeList]
        
    plt.plot(range(20,301,10),extremeList)    
    plt.xlabel('relative synaptic z postion $\mu$m')
    plt.ylabel('minimal potential $\mu$V')
    plt.savefig('averaged_real_extreme_point_versus_z_pos', bbox_inches='tight')
    
    
            