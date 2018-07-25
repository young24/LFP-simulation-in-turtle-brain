# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 13:46:08 2016
Main_Sim_with_Kernel
@author: superuser
"""
import os
from os.path import join
import time
import multiprocessing
import numpy as np
from scipy.interpolate import RegularGridInterpolator

def make_2D_to_3D(data,xLen,yLen):
    'make linear xy index into 2d index'
    data3D = np.zeros((xLen,yLen,np.shape(data)[1]))
    for x in range(0,xLen):
        for y in range(0,yLen):
            data3D[x,y,:] = data[x*yLen+y,:]
    return data3D
    
def calc_LFP(t):
    print(t) # show the progress
    xLen = 11
    yLen = 11
    lengthMEA = 500
    zMin = -110
    zMax = 220
    zShift = 20 # z shift between stimulated neuron and cell layer
    x = np.linspace(-lengthMEA,lengthMEA,xLen)
    y = np.linspace(-lengthMEA,lengthMEA,yLen)
    z = np.linspace(zMin,zMax,34)
    
    kernelData = np.load('../Data/Python/kernelData_soma_z120.npy')
    axonSyn = np.load('../Data/Python/axonSyn.npy')
    
    LFPonMEA = np.zeros((xLen,yLen))
    data = kernelData[:,t,:]
    data3D = make_2D_to_3D(data,xLen,yLen)
    LFP = RegularGridInterpolator((x, y, z), data3D)
    
    interval = 100
    for x_idx in range(0,xLen):
        for y_idx in range(0,yLen):
            sumLFP = 0
            for pos in axonSyn:
                if (-lengthMEA<=((x_idx-(xLen-1)/2)*interval-pos[0])<=lengthMEA and
                -lengthMEA<=((y_idx-(xLen-1)/2)*interval-pos[1])<=lengthMEA and
                zMin<=pos[2]-zShift<=zMax):    
                    sumLFP += LFP([(x_idx-(xLen-1)/2)*interval-pos[0],
                                   (y_idx-(yLen-1)/2)*interval-pos[1],pos[2]-zShift])
            LFPonMEA[x_idx,y_idx] = sumLFP
    folder = 'LFPonMEA'
    if not os.path.isdir(folder):
        os.mkdir(folder)
    np.save(join(folder, 'LFPonMEAt'+str(t)+'.npy'),LFPonMEA)
    
def make_files_together(xLen,yLen):
    'stack different time files into a single file'
    LFPonMEA = np.zeros((xLen,yLen,401))
    for t in range(0,401):
        LFPonMEA[:,:,t] = np.load('LFPonMEA/LFPonMEAt'+str(t)+'.npy')
    
    return LFPonMEA
        

if __name__ == '__main__':
    start = time.time()
    pool = multiprocessing.Pool(processes=4)
    
    t = range(0,401)
    pool.map(calc_LFP, t)
    
    pool.close()
    pool.join()
    
    xLen = 11 # keep consistent with before ones
    yLen = 11
    LFPonMEA = make_files_together(xLen,yLen)
    np.save('LFPonMEA.npy',LFPonMEA)
    
    
    end = time.time()
    print(end-start)

    
    