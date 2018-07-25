# -*- coding: utf-8 -*-
"""
Created on Wed May 11 11:22:55 2016
convert_npy_data_to_mat
@author: young
"""
import numpy as np
import scipy.io as sio
import sys


def convert_data(fileFolder,fileName):
    ' convert from .npy to .mat'
    data = np.load(fileFolder+'/'+fileName+'.npy')
    sio.savemat(fileFolder+'/'+fileName+'.mat',mdict={fileName: data}) 
    print('succeeded')

if __name__ == '__main__':
    fileFolder,fileName = sys.argv[1],sys.argv[2]
    convert_data(fileFolder,fileName)
    
