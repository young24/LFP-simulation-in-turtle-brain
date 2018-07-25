# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 11:47:18 2016
make_kernel_into_3D_array
@author: superuser
"""
import numpy as np



if __name__ == '__main__':
    increment = 10
    
    for start in range(-110,230,10):
        fileName = '../Data/Kernel/Kernel_somaZ120/phiMeanSynZ_'+str(start)+'to'+str(start+increment)+'.npy'
        if start == -110:
            kernelData = np.load(fileName)
        else: 
            temp = np.load(fileName)
            kernelData = np.dstack((kernelData,temp))
            
    np.save('kernelData.npy',kernelData)
