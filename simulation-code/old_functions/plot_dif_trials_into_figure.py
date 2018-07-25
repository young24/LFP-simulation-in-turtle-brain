# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 09:46:20 2016
plot different trials into a figure.

@author: young
"""
import pylab as plt
import numpy as np


####
dirName = 'sim_results'
xfileName = dirName+'/'+'tvec.npy'
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
    
    y1fileName = dirName+'/'+'phi_z'+str(z_pos)+'.npy'
    y1=np.load(y1fileName)
    y1=y1[2]
    ax1.plot(x, y1)
    y2fileName = dirName+'/'+'vmem_z'+str(z_pos)+'.npy'
    y2=np.load(y2fileName)
    y2=y2[0] # soma index = 0
    ax2.plot(x, y2)
    y3fileName = dirName+'/'+'imem_z'+str(z_pos)+'.npy'
    y3=np.load(y3fileName)
    y3=y3[0]
    ax3.plot(x, y3)

plt.legend(legendList, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig('simplest_z_profile', bbox_inches='tight')