# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 15:50:56 2016

plot_dif_radius_ex_poten_at_the_same_z
@author: superuser
"""

import pylab as plt
import numpy as np


####
dirName = 'sim_results'
xfileName = dirName+'/'+'tvec.npy'
x=np.load(xfileName)

posFileName = dirName+'/'+'syn_position.npy'
syn_position = np.load(posFileName)
r = np.zeros((syn_position.shape[0],1))
for i in range(0,syn_position.shape[0]):
    r[i]=np.sqrt(syn_position[i][0]**2+syn_position[i][1]**2)
sortedIndex = np.argsort(r,axis=0)
fig = plt.figure(figsize=[10, 10])
rel_dist = 10
legendList = []
for index in sortedIndex:
    legendList.append('r:'+str(int(r[index[0]][0]))+'$\mu$m')
    yfileName = dirName+'/'+'phi_z'+str(rel_dist)+'_'+str(index[0])+'.npy'
    y = np.load(yfileName)
    y = y[2]
    plt.plot(x, y)

plt.xlabel('ms')
plt.ylabel('$\mu$V')
plt.legend(legendList, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig('dif_radius_at_same_z.png', bbox_inches='tight')

