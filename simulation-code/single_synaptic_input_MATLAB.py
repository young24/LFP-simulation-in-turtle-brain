# -*- coding: utf-8 -*-
"""
Created on Wed May 11 10:30:35 2016
single_synaptic_input

This function is used to generate overview of single synaptic input,
including morpholohy plot, weighted morphology plot with return current, 
extracellular potential, input current, return current
@author: young
"""

import os
from os.path import join
import numpy as np
import scipy.io as sio
import pylab as plt
import neuron
import LFPy
import MoI
# TODO: Check if kCSD is usable
   
class Electrode:
    """ Class to make the electrodes.
    """
    def __init__(self, slice_thickness=300., elec_radius=15., elec_x=[0], elec_y=[0], elec_z=[0]):
        self.slice_thickness = slice_thickness
        self.elec_radius = elec_radius
        self.elec_x = np.array(elec_x, dtype=float)
        self.elec_y = np.array(elec_y, dtype=float)
        self.elec_z = np.array(elec_z, dtype=float)  # THIS IS AN ARRAY NOW!
        self.num_elecs = len(self.elec_x)
        self.elec_clr = [plt.cm.rainbow(1./(self.num_elecs - 1) * idx) for idx in range(self.num_elecs)]


class ReturnCurrents:
    """ Class to investigate the distribution of return currents following synaptic input, and calculated the
    extracellular potential at a microelectrode array (MEA) plane.
    """
    def __init__(self, morphologyPath, resultFolder, syn_pos, syn_weight, inputSign, isPassive, soma_z, upsideDown, simulate_cell=True, calculate_mapping=True):
        self.morphology_name = morphologyPath  # Could be .hoc or .swc and perhaps some other formats
        self.fig_name = 'testFig'
        self.folder = resultFolder #'sim_results'
        self.syn_pos = syn_pos  # [x, y, z] position of synaptic input
        self.syn_weight = syn_weight  # nA Strength of synaptic input
        self.input_spike_train = np.array([10.])  # Set time(s) of synaptic input
        self.simulate_cell = simulate_cell
        self.upsideDown = upsideDown # if the morphology is upsideDown
        self.passiveMechanism = isPassive
        self.inputSign = inputSign # 0 for excitatory input, -90 for inhibitory
        self.calculate_mapping = calculate_mapping
        cell, synapse, MEA = self.make_cell(soma_z)
        self.plot_results(cell, synapse, MEA)

    def make_cell(self, soma_z):
        
        neuron.h.celsius = 24 # Celsius, experiment temperature
        cell_parameters = {
            'morphology': morphologyPath,
            'v_init': -65,
            'passive': self.passiveMechanism,
            'nsegs_method': 'fixed_length',
            'max_nsegs_length': 10,
            'timeres_NEURON': 2**-3,  # Should be a power of 2
            'timeres_python': 2**-3,
            'tstartms': 0,
            'tstopms': 50,
            'pt3d': True,
            'verbose': True, # information about nsegs, rotation, and other information about the cell
        }   
        cell = LFPy.Cell(**cell_parameters)  
        # Specify the position and rotation of the cell
        if self.upsideDown:
            cell.set_rotation(y=np.pi)
        cell.set_pos(zpos=soma_z)
        cell = self.set_membrane_mechanism(cell, isPassive=self.passiveMechanism)
        
        elec_x, elec_y, elec_z = self.make_MEA(xStart=-600,xEnd=601,xResolution=300, 
         yStart=-600,yEnd=601,yResolution=300, 
         zStart=0,zEnd=141,zResolution=100)
        elec_params = {'slice_thickness': soma_z+np.max(cell.zend), # keep cell contact with the saline layer on the top
                       'elec_radius': 1.,
                       'elec_x': elec_x,
                       'elec_y': elec_y,
                       'elec_z': elec_z,
                       }
        MEA = Electrode(**elec_params)
        ext_sim_dict = {'use_line_source': False,
                        'n_elecs': MEA.num_elecs,
                        #'n_avrg_points': 100,
                        #'elec_radius': MEA.elec_radius,
                        'moi_steps': 20,
                        'elec_x': MEA.elec_x,
                        'elec_y': MEA.elec_y,
                        'elec_z': MEA.elec_z,
                        'slice_thickness': MEA.slice_thickness,
                        'include_elec': False,
                        'neural_input': '.',
                        }

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
            [ax4.plot(MEA.elec_x[idx], MEA.elec_z[idx], 's', c=MEA.elec_clr[idx], zorder=10) for idx in range(MEA.num_elecs)]  # MEA.elec_z IS NOW AN ARRAY
            [ax2.plot(MEA.elec_x[idx], MEA.elec_y[idx], 's', c=MEA.elec_clr[idx], zorder=10) for idx in range(MEA.num_elecs)]

            ax4.plot([-950, 950], [MEA.slice_thickness, MEA.slice_thickness], color='k')
            ax4.plot([-950, 950], [0, 0], 'k')  # PLOTTING MEA LINE
            ax4.axis([-1000, 1000, -20, MEA.slice_thickness + 20])
            plt.show()

        syn_idx = cell.get_closest_idx(x=self.syn_pos[0], y=self.syn_pos[1], z=self.syn_pos[2],section='dend')
        print('The syn_idx is:'+str(syn_idx))
                
        synapse_parameters = {
            'idx': syn_idx,
            'e': self.inputSign, #  Change to -90 for inhibitory input, and 0 for excitatory
            'syntype': 'Exp2Syn',
            'tau1': 2,
            'tau2': 10,
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
            np.save(join(self.folder, 'imem.npy'), cell.imem)
            np.save(join(self.folder, 'tvec.npy'), cell.tvec)
            np.save(join(self.folder, 'vmem.npy'), cell.vmem)
            np.save(join(self.folder, 'syn_i.npy'), synapse.i)
		# save as .mat files
            sio.savemat(join(self.folder, 'imem.mat'), mdict={'imem':cell.imem})
            sio.savemat(join(self.folder, 'tvec.mat'), mdict={'tvec':cell.tvec})
            sio.savemat(join(self.folder, 'vmem.mat'), mdict={'vmem':cell.vmem})
            sio.savemat(join(self.folder, 'syn_i.mat'), mdict={'syn_i':synapse.i})
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
            self.make_mapping(cell, MEA, ext_sim_dict)

        return cell, synapse, MEA
        
    def set_membrane_mechanism(self, cell, isPassive = True):
        if isPassive == True:
            for i,sec in enumerate(cell.allseclist):
                sec.insert('pas')
                sec.Ra = 75  # Ohm cm
                sec.g_pas = 1. / 30000
                sec.e_pas = -65
                if i == 0:            
                    sec.cm = 1.5  # uF / cm2
                else:                 
                    sec.cm = 3  # uF / cm2
        else:
            for i,sec in enumerate(cell.allseclist):
                sec.insert('hh')
                sec.Ra = 75  # Ohm cm Axial resistance
                sec.gnabar_hh = 0.083 # S/cm2
                sec.gkbar_hh = 0.03 # S/cm2
                sec.gl_hh = 0.0003 # S/cm2
                sec.ena = 55 # mV
                sec.ek = -72 # mV
                sec.el_hh = -49.3 # mV
                if i == 0:     
                    sec.cm = 1.5  # uF / cm2                
                else: 
                    if sec.name()[0:4] == 'axon':
                        sec.diam = 0.5 # microns
                    sec.cm = 3  # uF / cm2
                    
        return cell
            
    def make_MEA(self,xStart=-600,xEnd=601,xResolution=300, 
         yStart=-600,yEnd=601,yResolution=300, 
         zStart=0,zEnd=141,zResolution=100):
        ' set x,y sites, maybe z of MEA '
        elec_x, elec_y, elec_z = [], [], []   
        for z in range(zStart,zEnd,zResolution):
            for x in range(xStart,xEnd,xResolution):
                for y in range(yStart,yEnd,yResolution):    
                    elec_x.append(x), elec_y.append(y), elec_z.append(z)
                
        return elec_x,elec_y,elec_z

    def make_mapping(self, cell, MEA, ext_sim_dict):
        moi_normal_saline = {
            'sigma_G': 0.0,  # Below electrode
            'sigma_T': 0.3,  # Tissue
            'sigma_S': 1.5,  # Saline
            'h': MEA.slice_thickness,
            'steps': 20,
            }

        moi_normal_saline = MoI.MoI(**moi_normal_saline)

        mapping_normal_saline = moi_normal_saline.make_mapping_free_z(  # NOW USING METHOD THAT ALLOWS FREE Z-POSITION
            ext_sim_dict, xmid=cell.xmid, ymid=cell.ymid, zmid=cell.zmid, xstart=cell.xstart,
            ystart=cell.ystart, zstart=cell.zstart, xend=cell.xend, yend=cell.yend, zend=cell.zend)

        np.save(join(self.folder, 'mapping_normal_saline.npy'), mapping_normal_saline)
        if hasattr(cell, 'imem'):
            MEA.phi = 1000 * np.dot(mapping_normal_saline, cell.imem)
            np.save(join(self.folder, 'phi.npy'), MEA.phi)
            sio.savemat(join(self.folder, 'phi.mat'), mdict={'phi':MEA.phi})


    def plot_results(self, cell, synapse, MEA):

        time_window = [8, 18]
        soma_idx = 0
        syn_idx = synapse.idx

        cell_plot_idxs = [soma_idx, syn_idx]
        cell_plot_colors = {soma_idx: 'b', syn_idx: 'g'}
        num_cols = 5

        fig = plt.figure(figsize=[15, 5])
        plt.subplots_adjust(hspace=0.6, wspace=0.5, right=0.99, left=0.03, top=0.9)
        ax1 = fig.add_axes([0.05, 0.5, 0.3, 0.4], aspect=1, frameon=False, xticks=[], yticks=[])
        ax3 = fig.add_axes([0.05, 0.1, 0.3, 0.4], aspect=1, frameon=False, xticks=[], yticks=[])
        ax2 = plt.subplot(2, num_cols, 3, aspect=1, xlabel='x [$\mu m$]', ylabel='y [$\mu m$]', axisbg='c',
                          title='Morphology with color weighted\nby return current amplitude')
        ax4 = plt.subplot(2, num_cols, 8, aspect=1, xlabel='x [$\mu m$]', ylabel='z [$\mu m$]', axisbg='c', sharex=ax2)
        ax_ec = fig.add_subplot(1, num_cols, 4, xlim=time_window, xlabel='ms', ylabel='$\mu$V',
                                title='Extracellular\npotential')
        ax_v = plt.subplot(3, num_cols, 5, title='Membrane potential', ylabel='mV', ylim=[-70, -50], xlim=time_window)
        ax_ic = plt.subplot(3, num_cols, 10, title='Input current', ylabel='nA', xlim=time_window)
        ax_rc = plt.subplot(3, num_cols, 15, title='Return currents', ylabel='nA', xlabel='ms', xlim=time_window)

        l_elec, l_soma, l_syn = self._plot_recording_set_up(cell, ax1, ax3, MEA, soma_idx, syn_idx, cell_plot_colors)
        self._plot_current_weighted_morphology(cell, ax4, ax2, soma_idx)
        self._plot_intracellular_data(cell, cell_plot_colors, syn_idx, soma_idx, cell_plot_idxs, ax_v, ax_ic, ax_rc)
        self._plot_ec_potentials(ax_ec, MEA, cell)

        fig.legend([l_soma, l_syn, l_elec], ["Soma", "Synapse", "MEA electrode"],
                   frameon=False, numpoints=1, ncol=3, loc=3)
        simplify_axes([ax_v, ax_ec, ax_ic, ax_rc])

        mark_subplots([ax1, ax3, ax2, ax4, ax_ec, ax_v, ax_ic, ax_rc], ypos=1.1, xpos=-0.15)

        plt.savefig('%s.png' % self.fig_name)
        #plt.show()

    def _plot_intracellular_data(self, cell, cell_plot_colors, syn_idx, soma_idx, cell_plot_idxs, ax_v, ax_ic, ax_rc):

        [ax_v.plot(cell.tvec, cell.vmem[idx, :], c=cell_plot_colors[idx], lw=2) for idx in cell_plot_idxs]
        for idx in xrange(cell.totnsegs):
            ax = ax_ic if idx == syn_idx else ax_rc
            if idx == syn_idx:
                c = 'g'
                lw = 2
            elif idx == soma_idx:
                c = 'b'
                lw = 2
            else:
                c = 'k'
                lw = 1
            ax.plot(cell.tvec, cell.imem[idx], c=c, lw=lw)


    def _plot_current_weighted_morphology(self, cell, ax_side, ax_top, soma_idx):

        # Find maximum current amplitude in each compartment:
        max_current = np.max(np.abs(cell.imem), axis=1)
        # We normalize color to strongest return current, i.e. next strongest current:
        normalize_const = np.sort(max_current)[-2]
        # Find current normalized color of each cellular compartment:
        idx_clr = [plt.cm.gray_r(1./normalize_const * max_current[idx]) for idx in range(len(max_current))]

        for idx in np.argsort(max_current):  # Argsort to have important compartments in front (visible)
            if idx == soma_idx:
                ax_top.plot(cell.xmid[soma_idx], cell.ymid[soma_idx], 'o', ms=15, mec='none', c=idx_clr[soma_idx])
                ax_side.plot(cell.xmid[soma_idx], cell.zmid[soma_idx], 'o', ms=15, mec='none', c=idx_clr[soma_idx])
            else:
                #  Plot soma as ball
                ax_top.plot([cell.xstart[idx], cell.xend[idx]], [cell.ystart[idx], cell.yend[idx]],
                            c=idx_clr[idx], lw=1.5)
                ax_side.plot([cell.xstart[idx], cell.xend[idx]], [cell.zstart[idx], cell.zend[idx]],
                             c=idx_clr[idx], lw=1.5)

    def _plot_ec_potentials(self, ax, MEA, cell):
        for elec in range(MEA.num_elecs):
            ax.plot(cell.tvec, MEA.phi[elec], lw=2, c=MEA.elec_clr[elec])

    def _plot_recording_set_up(self, cell, ax_neur, side_ax, MEA, soma_idx, syn_idx, cell_plot_colors):

        for comp in xrange(len(cell.xmid)):
            if comp == 0:
                ax_neur.scatter(cell.xmid[comp], cell.ymid[comp], s=cell.diam[comp],
                                edgecolor='none', color='gray', zorder=1)
            else:
                ax_neur.plot([cell.xstart[comp], cell.xend[comp]],
                             [cell.ystart[comp], cell.yend[comp]],
                             lw=cell.diam[comp]/2, color='gray', zorder=1)

        for comp in xrange(len(cell.xmid)):
            if comp == 0:
                side_ax.scatter(cell.xmid[comp], cell.zmid[comp], s=cell.diam[comp],
                                edgecolor='none', color='gray', zorder=1)
            else:
                side_ax.plot([cell.xstart[comp], cell.xend[comp]],
                             [cell.zstart[comp], cell.zend[comp]],
                             lw=cell.diam[comp]/2, color='gray', zorder=1)

        [side_ax.plot(MEA.elec_x[idx], MEA.elec_z[idx], 's', c=MEA.elec_clr[idx], zorder=10) for idx in range(MEA.num_elecs)]  # MEA.elec_z IS NOW AN ARRAY
        [ax_neur.plot(MEA.elec_x[idx], MEA.elec_y[idx], 's', c=MEA.elec_clr[idx], zorder=10) for idx in range(MEA.num_elecs)]

        l_elec, = ax_neur.plot(MEA.elec_x[0], MEA.elec_y[0], 's', c=MEA.elec_clr[0], zorder=0)

        l_soma, = ax_neur.plot(cell.xmid[soma_idx], cell.ymid[soma_idx], '*', c=cell_plot_colors[soma_idx], ms=15)
        l_syn, = ax_neur.plot(cell.xmid[syn_idx], cell.ymid[syn_idx], '*', c=cell_plot_colors[syn_idx], ms=15)

        side_ax.plot(cell.xmid[soma_idx], cell.zmid[soma_idx], '*', c=cell_plot_colors[soma_idx], ms=15)
        side_ax.plot(cell.xmid[syn_idx], cell.zmid[syn_idx], '*', c=cell_plot_colors[syn_idx], ms=15)

        side_ax.plot(cell.xmid[soma_idx], cell.zmid[soma_idx], '*', c=cell_plot_colors[soma_idx], ms=15)
        side_ax.plot(cell.xmid[syn_idx], cell.zmid[syn_idx], '*', c=cell_plot_colors[syn_idx], ms=15)

        ax_neur.plot(cell.xmid[soma_idx], cell.ymid[soma_idx], '*', c=cell_plot_colors[soma_idx], ms=15)
        ax_neur.plot(cell.xmid[syn_idx], cell.ymid[syn_idx], '*', c=cell_plot_colors[syn_idx], ms=15)

        ax_neur.arrow(-220, -100, 30, 0, lw=1, head_width=7, color='k', clip_on=False)
        ax_neur.arrow(-220, -100, 0, 30, lw=1, head_width=7, color='k', clip_on=False)
        ax_neur.text(-170, -100, 'x', size=10, ha='center', va='center', clip_on=False)
        ax_neur.text(-220, -50, 'y', size=10, ha='center', va='center', clip_on=False)

        side_ax.arrow(-1200, -50, 80, 0, lw=1, head_width=30, color='k', clip_on=False)
        side_ax.text(-1000, -50, 'x', size=10, ha='center', va='center')
        side_ax.arrow(-1200, -50, 0, 100, lw=1, head_width=30, color='k', clip_on=False)
        side_ax.text(-1200, 180, 'z', size=10, ha='center', va='center')

        side_ax.plot([-950, 950], [MEA.slice_thickness, MEA.slice_thickness], color='k')
        side_ax.plot([-950, 950], [0, 0], 'k')  # PLOTTING BOTTOM OF MEA
        side_ax.axis([-1000, 1000, -20, MEA.slice_thickness + 20])
        side_ax.plot([1000, 1000], [0, MEA.slice_thickness],
                     color='k', lw=2, clip_on=False)
        side_ax.text(1050, MEA.slice_thickness/2, '%g $\mu m$' % MEA.slice_thickness, size=8, va='center')

        side_ax.text(900, -50, 'MEA', va='top', ha='right')
        side_ax.text(900, MEA.slice_thickness, 'Saline', va='bottom', ha='right')
        side_ax.text(900, MEA.slice_thickness-50, 'Tissue', va='top', ha='right')
        return l_elec, l_soma, l_syn


def simplify_axes(axes):
    """
    :param axes: The axes object or list that is to be simplified. Right and top axis line is removed
    :return:
    """
    if not type(axes) is list:
        axes = [axes]

    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()


def mark_subplots(axes, letters='ABCDEFGHIJKLMNOPQRSTUVWXYZ', xpos=-0.12, ypos=1.15):
    """ Marks subplots in axes (should be list or axes object) with capital letters
    """
    if not type(axes) is list:
        axes = [axes]

    for idx, ax in enumerate(axes):
        ax.text(xpos, ypos, letters[idx].capitalize(),
                horizontalalignment='center',
                verticalalignment='center',
                fontweight='demibold',
                fontsize=12,
                transform=ax.transAxes)

if __name__ == '__main__':
    inputArgument = sio.loadmat('inputFor_single_synaptic_input_MATLAB.mat')
    
    morphologyPath = inputArgument['morphologyPath'][0].encode('ascii','ignore')# 'morphology_files/reconstructed/W66N1.hoc'
     
    resultFolder = inputArgument['resultFolder'][0].encode('ascii','ignore')# 'sim_results'
    soma_z = inputArgument['soma_z'][0][0]# 150  # z-position of cell, relative to MEA plane at 0 um
    syn_pos = inputArgument['syn_pos'][0]# np.array([0, 0, soma_z+250])  # Synapse is inserted at closest cellular position to specified point (x,y,z)
    syn_weight = inputArgument['syn_weight'][0][0]# 0.05 # nA; for active membrane: 0.05; for passive membrane: 0.005
    inputSign = inputArgument['inputSign'][0][0]#0 # 0 for excitatory synapse, -90 for inhibitory
    isPassive = inputArgument['isPassive'][0][0]# False
    upsideDown = inputArgument['upsideDown'][0][0]# True # if morphology is upsideDown
    rc = ReturnCurrents(morphologyPath, resultFolder, syn_pos, syn_weight, inputSign, isPassive, soma_z, upsideDown,
                        calculate_mapping=True, simulate_cell=True)

