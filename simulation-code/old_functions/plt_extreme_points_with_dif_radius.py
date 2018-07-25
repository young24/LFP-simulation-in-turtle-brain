# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 16:36:41 2016
plt_extreme_points_with_dif_radius
@author: superuser
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
    rel_dist = 151
    dirName = 'sim_results'
    xfileName = dirName+'/'+'tvec.npy'
    x=np.load(xfileName)
    
    posFileName = dirName+'/'+'syn_position_z_'+str(rel_dist)+'.npy'
    syn_position = np.load(posFileName)
    r = np.zeros((syn_position.shape[0],1))
    for i in range(0,syn_position.shape[0]):
        r[i]=np.sqrt(syn_position[i][0]**2+syn_position[i][1]**2)
    sortedIndex = np.argsort(r,axis=0)
    
    fig = plt.figure(figsize=[10, 10])
    
    radiusList = []
    extremeList = []
    for index in sortedIndex:
        yfileName = dirName+'/'+'phi_z'+str(rel_dist)+'_'+str(index[0])+'.npy'
        y = np.load(yfileName)
        extremum = find_extremum_in_list(y,2)
        extremeList.append(extremum)
        radiusList.append(r[index[0]])
        
        
    # normalization
#    minValue = min(extremeList)
#    maxValue = max(extremeList)
#    extremeList = [(n-minValue)/(maxValue-minValue) for n in extremeList]
        
    plt.plot(radiusList,extremeList)    
    plt.xlabel('radius to soma $\mu$m')
    plt.ylabel('extreme potential point $\mu$V')
    plt.savefig('z_'+str(rel_dist)+'_extreme_point_versus_radius', bbox_inches='tight')
    
    
            