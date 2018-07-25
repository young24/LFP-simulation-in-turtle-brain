# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 10:54:24 2016
average_along_z_for_kernel
After getting the data from get_z_profile_for_kernel, 
one still need to average the raw data along z for z precision of 10 microns.
@author: superuser
"""
import numpy as np
import os
def find_index_for_each_increment(synPosition):
    'To get the first index for each increment in a 2D array(x,y,sortedZ)'
    increment = 10
    numSyn = np.shape(synPosition)[0]
    start = int(synPosition[0,2])/(increment/2) *(increment/2) # to make the final value at x10.
    # end = (int(synPosition[numSyn-1,2])/increment + 1)*increment
    # numIndex = (end-start)/increment
    index = []
    comparedNum = start
    for i in range(0,numSyn):
        if (synPosition[i,2]>comparedNum):
            index.append(i)
            comparedNum += increment
    index.append(numSyn-1)
    return index
    
if __name__ == '__main__':
    loadFolder = '../Data/Kernel/sim_zProfileofKernel_120_f6000/'
    saveFolder = 'Kernel_somaZ120/'
    if not os.path.isdir(saveFolder):
                os.mkdir(saveFolder)
    
    increment = 10
    start = int(synPosition[0,2])/(increment/2) *(increment/2)
    index = find_index_for_each_increment(synPosition)
    for i in range(0,len(index)-1):
        phi = 0
        for j in range(index[i],index[i+1]):
            loadFile = 'phiMean'+str(j)+'.npy'
            temp = np.load(loadFolder+loadFile)
            phi += temp
        phi = phi/(index[i+1]-index[i])
        saveFileName = 'phiMeanSynZ'+str(start+increment/2)
        np.save(saveFolder+saveFileName,phi)
        start += increment
        
        
        
    