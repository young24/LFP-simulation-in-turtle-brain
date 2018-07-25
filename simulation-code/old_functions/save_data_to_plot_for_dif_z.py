# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 08:05:02 2016

@author: young
"""
import time

import os
from os.path import join
import numpy as np
import pylab as plt
import neuron
import LFPy
import MoI
# TODO: Check if kCSD is usable

class Electrode:
    """ Class to make the electrodes.
    """
    def __init__(self, slice_thickness=300., elec_radius=15., elec_x=[0], elec_y=[0]):
        self.slice_thickness = slice_thickness
        self.elec_radius = elec_radius
        self.elec_x = np.array(elec_x, dtype=float)
        self.elec_y = np.array(elec_y, dtype=float)
        self.elec_z = 0
        self.num_elecs = len(self.elec_x)
        self.elec_clr = [plt.cm.rainbow(1./(self.num_elecs - 1) * idx) for idx in range(self.num_elecs)]


class ReturnCurrents:
    """ Class to investigate the distribution of return currents following synaptic input, and calculated the
    extracellular potential at a microelectrode array (MEA) plane.
    """
    def __init__(self, syn_pos, ext_sim_dict, cell_z_pos, rel_dist, MEA, simulate_cell=True, calculate_mapping=True):

        self.morphology_name = 'morphology_files/Import3D_mod.hoc'  # Could be .hoc or .swc and perhaps some other formats
        self.folder = 'sim_results'
        self.syn_pos = syn_pos  # [x, y, z] position of synaptic input
        self.rel_dist = rel_dist
        self.syn_weight = 0.005  # Strength of synaptic input
        self.input_spike_train = np.array([10.])  # Set time(s) of synaptic input
        self.simulate_cell = simulate_cell
        self.calculate_mapping = calculate_mapping
        cell, synapse = self.make_cell(ext_sim_dict, cell_z_pos, MEA)
        #self.plot_results(cell, synapse, MEA)

    def make_cell(self, ext_sim_dict, cell_z_pos, MEA):

        cell_parameters = {
            'morphology': self.morphology_name,
            'v_init': -65,
            'passive': False,
            'nsegs_method': 'lambda_f',
            'lambda_f': 400,  # Higher number means higher spatial resolution of the cell
            'timeres_NEURON': 2**-3,  # Should be a power of 2
            'timeres_python': 2**-3,
            'tstartms': 0,
            'tstopms': 50,
            'verbose': True,
            'pt3d': True,
        }
        cell = LFPy.Cell(**cell_parameters)

        # Setting the passive parameters of the cell:
        for sec in cell.allseclist:
            sec.Ra = 100  # Ohm cm
            sec.cm = 1  # uF / cm2
            sec.insert('pas')
            sec.g_pas = 1. / 30000
            sec.e_pas = -65

        # Specify the position and rotation of the cell
        cell.set_pos(zpos=cell_z_pos)
        #cell.set_rotation(x=np.pi,z=0*np.pi)

        plot_morph = True # the original value is True
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
            [ax4.plot(MEA.elec_x[idx], MEA.elec_z, 's', c=MEA.elec_clr[idx], zorder=10) for idx in range(MEA.num_elecs)]
            [ax2.plot(MEA.elec_x[idx], MEA.elec_y[idx], 's', c=MEA.elec_clr[idx], zorder=10) for idx in range(MEA.num_elecs)]

            ax4.plot([-950, 950], [MEA.slice_thickness, MEA.slice_thickness], color='k')
            ax4.plot([-950, 950], [MEA.elec_z, MEA.elec_z], 'k')
            ax4.axis([-1000, 1000, -20, MEA.slice_thickness + 20])
            plt.show()


        syn_idx = cell.get_closest_idx(x=self.syn_pos[0], y=self.syn_pos[1], z=self.syn_pos[2], section='dend')
        synapse_parameters = {
            'idx': syn_idx,
            'e': 0., #  Change to -90 for inhibitory input, and 0 for excitatory
            'syntype': 'Exp2Syn',
            'tau1': 1.,
            'tau2': 1.,
            'weight': self.syn_weight,
            'record_current': True,
        }
        synapse = LFPy.Synapse(cell, **synapse_parameters)
        synapse.set_spike_times(self.input_spike_train)

        # If no parameters are changed, we can just load results from file
        if self.simulate_cell:
            print "Simulating cell"
            cell.simulate(rec_imem=True, rec_vmem=True, rec_isyn=True)
            if not os.path.isdir(self.folder):
                os.mkdir(self.folder)
            print cell.imem.shape
            imemName = 'imem_z'+str(self.rel_dist)+'.npy' 
            np.save(join(self.folder, imemName), cell.imem)
            np.save(join(self.folder, 'tvec.npy'), cell.tvec)
            vmemName = 'vmem_z'+str(self.rel_dist)+'.npy' 
            np.save(join(self.folder, vmemName), cell.vmem)
            np.save(join(self.folder, 'syn_i.npy'), synapse.i)
            if os.path.isfile(join(self.folder, 'mapping_normal_saline.npy')) and not self.calculate_mapping:
                mapping = np.load(join(self.folder, 'mapping_normal_saline.npy'))
                MEA.phi = np.dot(mapping, cell.imem)
        else:
            cell.imem = np.load(join(self.folder, 'imem.npy'))
            cell.vmem = np.load(join(self.folder, 'vmem.npy'))
            cell.tvec = np.load(join(self.folder, 'tvec.npy'))
            synapse.i = np.load(join(self.folder, 'syn_i.npy'))
            MEA.phi = np.load(join(self.folder, 'phi.npy'))
        if self.calculate_mapping:
            self.make_mapping(cell, MEA, ext_sim_dict, cell_z_pos)

        return cell, synapse

    def make_mapping(self, cell, MEA, ext_sim_dict, cell_z_pos):
        moi_normal_saline = {
            'sigma_G': 0.0,  # Below electrode
            'sigma_T': 0.3,  # Tissue
            'sigma_S': 1.5,  # Saline
            'h': MEA.slice_thickness,
            'steps': 20,
            }

        moi_normal_saline = MoI.MoI(**moi_normal_saline) # Modified by wenbin

        mapping_normal_saline = moi_normal_saline.make_mapping_cython(
            ext_sim_dict, xmid=cell.xmid, ymid=cell.ymid, zmid=cell.zmid, xstart=cell.xstart,
            ystart=cell.ystart, zstart=cell.zstart, xend=cell.xend, yend=cell.yend, zend=cell.zend)
        #mappingName = 'mapping_normal_saline'+str(self.rel_dist)+'.npy'
        #np.save(join(self.folder, mappingName), mapping_normal_saline)
        if hasattr(cell, 'imem'):
            MEA.phi = 1000 * np.dot(mapping_normal_saline, cell.imem)
            phiName = 'phi_z'+str(self.rel_dist)+'.npy' 
            np.save(join(self.folder, phiName), MEA.phi)

if __name__ == '__main__':
    start = time.time()
    elec_params = {'slice_thickness': 451.,
                   'elec_radius': 1.,
                   'elec_x': [-50, 50, 0, -50, 50],   # x-position of all electrodes
                   'elec_y': [50, 50, 0, -50, -50]    # y-position of all electrodes
                   }

    folder = join('sim_results')
    MEA = Electrode(**elec_params)
    ext_sim_dict = {'use_line_source': False,
                    'n_elecs': MEA.num_elecs,
                    #'n_avrg_points': 100,
                    #'elec_radius': MEA.elec_radius,
                    'moi_steps': 20,
                    'elec_x': MEA.elec_x,
                    'elec_y': MEA.elec_y,
                    'elec_z': 0,
                    'slice_thickness': MEA.slice_thickness,
                    'include_elec': False,
                    'neural_input': '.',
                    }
    cell_z_pos = 150
    
    for rel_dist in range(20,301,10): # default value 150  # z-position of cell, relative to MEA plane at 0 um 
        synapse_pos = np.array([0, 0, cell_z_pos+rel_dist])  # Synapse is inserted at closest cellular position to specified point (x,y,z)
        rc = ReturnCurrents(synapse_pos, ext_sim_dict, cell_z_pos, rel_dist, MEA,
                        calculate_mapping=True, simulate_cell=True)
                        
    end = time.time()
    print(end-start)


