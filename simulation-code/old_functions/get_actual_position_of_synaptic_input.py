# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 13:35:38 2016
get_actual_position_of_synaptic_input
@author: young
"""


import numpy as np
import neuron
import LFPy
# TODO: Check if kCSD is usable


if __name__ == '__main__':
    
    #cell_z_pos = 150
    cell_parameters = {
        'morphology': 'morphology_files/Import3D_mod.hoc',
        'v_init': -65,
        'passive': False,
        'nsegs_method': 'lambda_f',
        'lambda_f': 400,  # Higher number means higher spatial resolution of the cell
        'timeres_NEURON': 2**-3,  # Should be a power of 2
        'timeres_python': 2**-3,
        'tstartms': 0,
        'tstopms': 50,
        'pt3d': True,
        'verbose': True,
    }
    cell = LFPy.Cell(**cell_parameters)

    # Specify the position and rotation of the cell
    #cell.set_pos(zpos=cell_z_pos)
    cell.set_rotation(x=np.pi)
    
    actualSynPos = [0, 0, 0]
    for syn_pos in synPosition: # default value 150  # z-position of cell, relative to MEA plane at 0um 
        syn_idx = cell.get_closest_idx(x=syn_pos[0], y=syn_pos[1], z=syn_pos[2],section='dend')
        act_xPos = cell.xmid[syn_idx]
        act_yPos = cell.ymid[syn_idx]
        act_zPos = cell.zmid[syn_idx]
        temp = [act_xPos, act_yPos, act_zPos]
        actualSynPos = np.vstack((actualSynPos,temp))
        
    actualSynPos=np.delete(actualSynPos, 0, 0) 
    np.save('actaulSynPos.npy',actualSynPos)
    