# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 09:12:00 2016
Average the spatial influence of the morphology
@author: young
"""
import os
from os.path import join
import pylab as plt
import numpy as np

def average_data_in_dif_folder(dirName,dataName,z_pos,numRotation,idxSection):
    yfileName = dirName+'_0/'+dataName+str(z_pos)+'.npy'
    y=np.load(yfileName)
    y=y[idxSection]
    for i in range(1,numRotation,1):
        yfileName = dirName+'_'+str(i)+'/'+dataName+str(z_pos)+'.npy'
        temp_y = np.load(yfileName)
        temp_y = temp_y[idxSection]
        y = y + temp_y
    y = y/numRotation
    return y
    
if __name__ == '__main__':

    numRotation = 8
    
    dirName = 'sim_results'
    xfileName = dirName+'_0/'+'tvec.npy'
    x=np.load(xfileName)
    
    fig = plt.figure(figsize=[10, 10])
    # Extracellular_potential
    ax1 = plt.subplot(311, ylabel='$\mu$V',
                      title='Extracellular\npotential')
    # share x only
    ax2 = plt.subplot(312, sharex=ax1, ylabel='mV',
                      title='Membrane potential')
    ax3 = plt.subplot(313, sharex=ax1, xlabel='ms', ylabel='nA',
                      title='Return currents')
    
    legendList = []
    for z_pos in range(20,301,10):
        legendList.append('z:'+str(z_pos)+'$\mu$m')
        
        dataName='phi_z'
        y1 = average_data_in_dif_folder(dirName,dataName,z_pos,numRotation,idxSection=2) # centre electrode
        ax1.plot(x, y1)
        dataName='phi_z'+str(z_pos)+'.npy'
        np.save(join('averaged_result_for_real_neuron', dataName), y1)
        
        dataName='vmem_z'
        y2 = average_data_in_dif_folder(dirName,dataName,z_pos,numRotation,idxSection=0) # soma index = 0
        ax2.plot(x, y2)
        dataName='vmem_z'+str(z_pos)+'.npy'
        np.save(join('averaged_result_for_real_neuron', dataName), y2)
        
        dataName='imem_z'
        y3 = average_data_in_dif_folder(dirName,dataName,z_pos,numRotation,idxSection=0)
        ax3.plot(x, y3)
        dataName='imem_z'+str(z_pos)+'.npy'
        np.save(join('averaged_result_for_real_neuron', dataName), y3)
    
    plt.legend(legendList, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig('averaged_z_profile', bbox_inches='tight')

    