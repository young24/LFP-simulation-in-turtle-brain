# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 13:46:08 2016
Main_Sim_with_Kernel
@author: Young
"""
import os
from os.path import join
import time
import multiprocessing
import itertools
import numpy as np
from scipy.interpolate import RegularGridInterpolator

    
def calc_LFP(t,folder):
    print(t) # show the progress
    xLen = 11
    yLen = 11
    # lengthMEA = 500
    
    zShift = 100 # z shift between stimulated neuron and cell layer
    zResolution = 10
    MEA_to_soma = 2 * zResolution # x10 um
    
    folderName = 'test'
    kernelPotential = np.load(folderName+'/KernelPotential_soma_above_MEA_z'+str(MEA_to_soma)+'.npy')
    coordinates = np.load(folderName+'/coordinates.npy')
    axonSyn = np.load('../Data/Python/Pyramidal_axonSyn.npy')
    
    rMax = np.max(coordinates[1]) # rResolution = 30
    zMin,zMax = np.min(coordinates[0]),np.max(coordinates[0])
    
    LFPonMEA = np.zeros((xLen,yLen))
    data = kernelPotential[:,:,t]
    PHI = RegularGridInterpolator(coordinates, data) # (z.r)
    
    interval = 100
    for x_idx in range(0,xLen):
        for y_idx in range(0,yLen):
            sumLFP = 0
            for pos in axonSyn:
                x = (x_idx-(xLen-1)/2)*interval-pos[0]
                y = (y_idx-(yLen-1)/2)*interval-pos[1]
                r = np.sqrt(x**2+y**2)
                if (r <= rMax and zMin<=pos[2]+zShift<=zMax):                  
                    sumLFP += PHI([pos[2]+zShift,r]) # (z,r)
                    
            LFPonMEA[x_idx,y_idx] = sumLFP
    if not os.path.isdir(folder):
        os.mkdir(folder)
    np.save(join(folder, 'LFPonMEAt'+str(t)+'.npy'),LFPonMEA)
    
def make_files_together(folder,xLen=11,yLen=11,timeSteps=401):
    'stack different time files into a single file'
    LFPonMEA = np.zeros((xLen,yLen,timeSteps))
    for t in range(0,timeSteps):
        filePath = folder+'/LFPonMEAt'+str(t)+'.npy'
        LFPonMEA[:,:,t] = np.load(filePath)
        os.remove(filePath)
    
    return LFPonMEA

def func_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return calc_LFP(*a_b)
        
if __name__ == '__main__':
    start = time.time()
    multiprocessing.freeze_support()
    pool = multiprocessing.Pool(processes=32)
    
    
    folderPath = 'LFP_delete_basal_dend'
    if not os.path.isdir(folderPath):
        os.mkdir(folderPath)
    t = range(0,401)
    pool.map(func_star, itertools.izip(t, itertools.repeat(folderPath)))
    
    pool.close()
    pool.join()
    
    LFPonMEA = make_files_together(folderPath,xLen=11,yLen=11)
    np.save(folderPath+'/LFPonMEA.npy',LFPonMEA)
        
    end = time.time()
    print(end-start)

    
    
