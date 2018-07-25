# -*- coding: utf-8 -*-
"""
Created on Tue May 10 12:01:20 2016
plot_morphology
@author: superuser
"""
import numpy as np
import pylab as plt
import neuron
import LFPy
def make_cell(morphologyPath, cell_z_pos, upsideDown=False):

    cell_parameters = {
        'morphology': morphologyPath,
        'v_init': -65,
        'passive': True,
        'nsegs_method': 'lambda_f',
        'lambda_f': 100,  # Higher number means higher spatial resolution of the cell
        'timeres_NEURON': 2**-3,  # Should be a power of 2
        'timeres_python': 2**-3,
        'tstartms': 0,
        'tstopms': 50,
    }
    cell = LFPy.Cell(**cell_parameters)
    # Specify the position and rotation of the cell
    cell.set_pos(zpos=cell_z_pos)
    if upsideDown:
        cell.set_rotation(y=np.pi)

    plot_morph = True
    if plot_morph:
        ax2 = plt.subplot(1, 2, 1, aspect=1, xlabel='x [$\mu m$]', ylabel='y [$\mu m$]', axisbg='c')
        ax4 = plt.subplot(1, 2, 2, aspect=1, xlabel='x [$\mu m$]', ylabel='z [$\mu m$]', axisbg='c', sharex=ax2)

        for idx in range(cell.totnsegs):  # Argsort to have important compartments in front (visible)
            if idx == 0:
                ax2.plot(cell.xmid[0], cell.ymid[0], 'o', ms=15, mec='none', c='k')
                ax4.plot(cell.xmid[0], cell.zmid[0], 'o', ms=15, mec='none', c='k')
            else:
                #  Plot soma as ball
                ax2.plot([cell.xstart[idx], cell.xend[idx]], [cell.ystart[idx], cell.yend[idx]], c='k', lw=1.5)
                ax4.plot([cell.xstart[idx], cell.xend[idx]], [cell.zstart[idx], cell.zend[idx]], c='k', lw=1.5)
        
        plt.show()

if __name__ == '__main__':
    morphologyPath = 'morphology_files/V6N2.hoc'
    cell_z_pos = 300
    make_cell(morphologyPath, cell_z_pos, upsideDown=False)
    