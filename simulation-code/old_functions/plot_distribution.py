# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 15:45:08 2016
plot distribution of position by histgram
@author: young
"""

import numpy as np
import pylab as plt

listName = z
numberArray = []
for i in range(0,len(listName)):
    numberArray = np.append(numberArray,listName[i])

plt.hist(numberArray,50)
plt.title('histgram of z distribution of all segments')
plt.xlabel('z $\mu$m')
plt.ylabel('number of segments')
plt.savefig('z_distribution.png')