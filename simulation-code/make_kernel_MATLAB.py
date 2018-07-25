# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 10:33:02 2016
make_kernel

Basic Assumptions:
1. There's only 1 synapse between a axon-site and a dendrite-site. (sparse enough)
2. Every rotation angle has the same chance. So we use the same possibility to sum up all angles. 
The axis is the z axis go trough the synapse 
 

Input parameters:
1. morphology file
2. resultFolder

Output:
1. phiKernel_soma_above_MEA_z.mat # cell in MATLAB, 
2. zCoordinate.mat
3. rCoordinate.mat


Main Steps:
1. Get synapses on dendrites
2. Sort synapses by z
3. Insert synapses one by one
4. Set sample sites across the entire neuron, the open side of V oritated to the neuron
denser close to the soma and sparse far away from the soma.
5. Get potential on MEA
6. Average it radically symmetrically. (need technical support)
7. Sum up and average different synapses for the same z zone.
8. Get the kernel.

@author: Young
"""


import time
import multiprocessing
import itertools
import sys
import os
from os.path import join
import numpy as np
import scipy.io as sio
import neuron
import LFPy
import MoI


class Electrode:
    """ Class to make the electrodes."""
    def __init__(self, slice_thickness=300., elec_radius=15., elec_x=[0], elec_y=[0], elec_z=[0]):
        self.slice_thickness = slice_thickness
        self.elec_radius = elec_radius
        self.elec_x = np.array(elec_x, dtype=int)
        self.elec_y = np.array(elec_y, dtype=int)
        self.elec_z = np.array(elec_z, dtype=int)  # THIS IS AN ARRAY NOW!
        self.num_elecs = len(self.elec_x)
        
def get_inserted_synPos(morphologyPath,soma_z=150,upsideDown=True):
    'synPos is relative to soma, in other words, soma_z is 0 in this case'
    cell=make_cell(morphologyPath,soma_z,upsideDown)
    synPos = get_syn_pos(cell)
    return synPos

def get_syn_pos(cell):
    'get syn pos from the neuron, sorted by z'
    idxDend = cell.get_idx(section='dend')  
    pos = np.array([cell.xmid,cell.ymid,cell.zmid],dtype=int).T
    pos = pos[idxDend]
    pos = pos[pos[:,2].argsort()] # sorted by z
    
    return pos
    
def func_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return simulate_cell(*a_b)
    
def simulate_cell(index,morphologyPath,upsideDown,resultFolder):# Main function
    # simulation parameters
    neuron.h.celsius = 24 # Celsius, experiment temperature
    inputSign = 0 # 0 for excitatory input, -90 for inhibitory
    soma_z = 150 # um
    passiveMechanism = True
    syn_tau1 = 1. # ms
    syn_tau2 = 1. # ms
    syn_weight = 0.005 # nA, for passive membrane
    syn_input_train = np.array([10.]) # the timing of onset 
    
    cell = make_cell(morphologyPath,soma_z,upsideDown) # Build up the cell and set in the right position and orientation  
    cell = set_membrane_mechanism(cell, isPassive = passiveMechanism) # choose passive/active membrane
    
    elec_x, elec_y, elec_z = make_MEA(xStart=-600,xEnd=600,xResolution=30,
             yStart=-600,yEnd=600,yResolution=30,
             zStart=10,zEnd=30,zResolution=10)
    elec_params = {'slice_thickness': soma_z+np.max(cell.zend), # keep cell contact with the saline layer on the top
                   'elec_radius': 1., # um
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
    # insert the synapse
    synPos = get_syn_pos(cell)
    synIdx = cell.get_closest_idx(x=synPos[index,0], y=synPos[index,1], z=synPos[index,2],section='dend')  # Torbjorn
    synapse_parameters = {
        'idx': synIdx,
        'e': inputSign, #  Change to -90 for inhibitory input, and 0 for excitatory
        'syntype': 'Exp2Syn',
        'tau1': syn_tau1,
        'tau2': syn_tau2,
        'weight': syn_weight, # 0.005 for passive membrane
        'record_current': True,
    }
    synapse = LFPy.Synapse(cell, **synapse_parameters)
    synapse.set_spike_times(syn_input_train)  # Set time(s) of synaptic input
        
    #print "Simulating cell"
    cell.simulate(rec_imem=True, rec_vmem=True, rec_isyn=True)   
    phi = make_mapping(cell, MEA, ext_sim_dict)
    phiProfile = average_phi_radically(phi,MEA,[synPos[index,0],synPos[index,1]],resultFolder,rResolution=30)
    np.save(join(resultFolder, 'phiProfile'+str(index)+'.npy'),phiProfile)
    
def make_cell(morphologyPath,soma_z,upsideDown):
                   
    cell_parameters = {
            'morphology': morphologyPath,
            'v_init': -65, # mV
            'passive': False,
            'nsegs_method': 'fixed_length',
            'max_nsegs_length': 10,  # um, max segment length for method 'fixed_length'
            'timeres_NEURON': 2**-3,  # Should be a power of 2
            'timeres_python': 2**-3,
            'tstartms': 0,  # ms
            'tstopms': 50, # ms
            'pt3d': True,
            'verbose': True, # information about nsegs, rotation, and other information about the cell
        }  
    cell = LFPy.Cell(**cell_parameters)  
    # Specify the position and rotation of the cell
    if upsideDown:
        cell.set_rotation(y=np.pi)
    cell.set_pos(zpos=soma_z)
    
    return cell

def set_membrane_mechanism(cell, isPassive = True):
    if isPassive == True:
        for i,sec in enumerate(cell.allseclist):
            sec.insert('pas')
            sec.Ra = 75  # Ohm cm
            sec.g_pas = 1. / 30000 # S/cm2
            sec.e_pas = -65 # mV
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

def make_MEA(xStart=-600,xEnd=601,xResolution=30, 
             yStart=-600,yEnd=601,yResolution=30, 
             zStart=0,zEnd=141,zResolution=10):
    ' set x,y sites, maybe z of MEA '
    elec_x, elec_y, elec_z = [], [], []   
    for z in range(zStart,zEnd,zResolution):
        for x in range(xStart,xEnd,xResolution):
            for y in range(yStart,yEnd,yResolution):    
                elec_x.append(x), elec_y.append(y), elec_z.append(z)
                
    return elec_x,elec_y,elec_z

def make_mapping(cell, MEA, ext_sim_dict):
    moi_normal_saline = {
        'sigma_G': 0.0,  # Below electrode
        'sigma_T': 0.3,  # Tissue
        'sigma_S': 1.5,  # Saline
        'h': MEA.slice_thickness,
        'steps': 20,
        }
    moi_normal_saline = MoI.MoI(**moi_normal_saline)
    mapping_normal_saline = moi_normal_saline.make_mapping_free_z(
        ext_sim_dict, xmid=cell.xmid, ymid=cell.ymid, zmid=cell.zmid, xstart=cell.xstart,
        ystart=cell.ystart, zstart=cell.zstart, xend=cell.xend, yend=cell.yend, zend=cell.zend)

    if hasattr(cell, 'imem'):
        MEA.phi = 1000 * np.dot(mapping_normal_saline, cell.imem)
    return MEA.phi
    
def average_phi_radically(phi,MEA,center,resultFolder,rResolution=30): # syn position is the centerPos
    'average phi radically to get the perfect phi along r'
    x,y,z = MEA.elec_x,MEA.elec_y,MEA.elec_z
    idx = np.argwhere(z==z[0])
    r = np.sqrt((x[idx] - center[0])**2 + (y[idx] - center[1])**2)
    rRange = range(0,int(np.max(r))+rResolution,rResolution)
    binCount = np.histogram(r, bins=rRange)[0]
    idxGroup = [range(binCount[0:i].sum(),binCount[0:i+1].sum()) for i in range(0,len(binCount))]
    phiProfile = np.zeros((len(z)/len(idx),len(binCount),np.shape(phi)[1])) # (z,r,t)
    for i in range(0,np.shape(phiProfile)[0]):
        for j, value in enumerate(idxGroup):
            if value: # empty list will be false
                phiProfile[i,j,:] = phi[np.array(value,dtype=int)+i*len(idx)].sum(axis=0)/len(value)
      
    return phiProfile # (z,r,t)

def sum_over_same_z_zone(synZ,resultFolder,soma_to_bottom_MEA=150,zResolution = 10,rResolution=30):
    'sum over synapses in the same z zone and save seperate file for different z zones'
    zRange = range(int(synZ[0]),int(synZ[-1]),zResolution) # maybe desert some synapses
    binCount = np.histogram(synZ, bins=zRange)[0]
    idxGroup = [range(binCount[0:i].sum(),binCount[0:i+1].sum()) for i in range(0,len(binCount))]
    
    phiList = [np.load(resultFolder+'/phiProfile'+str(i)+'.npy') for i in range(0,sum(binCount))]
    minLen = np.min([np.shape(i)[1] for i in phiList])
    phiList = [i[:,0:minLen,:] for i in phiList]
    phiKernel = []
    for indices in idxGroup:
        if not indices:
            phiKernel.append(np.zeros(np.shape(phiList[0])))
        else:
            phiSum = sum([phiList[i] for i in indices])
            phiKernel.append(phiSum)

    for idx in range(0,np.shape(phiKernel[0])[0]):
        kernelPotential = np.array([i[idx] for i in phiKernel])
        matFilePath = resultFolder+'/phiKernel_soma_above_MEA_z'+str(soma_to_bottom_MEA-idx*zResolution)+'.mat'
        sio.savemat(matFilePath,mdict={'kernelPotential':kernelPotential})# (z,r,t)
        
    zCoordinate = np.array(zRange[0:-1])+zResolution/2 
    matFileName = resultFolder+'/zCoordinate.mat'
    sio.savemat(matFileName,mdict={'zCoordinate':zCoordinate})
    
    rCoordinate = np.array(range(0,minLen*rResolution,rResolution)) # no shift to make sure include r=0
    matFileName = resultFolder+'/rCoordinate.mat'
    sio.savemat(matFileName,mdict={'rCoordinate':rCoordinate})
              
        
if __name__ == '__main__':
    start=time.time()
    multiprocessing.freeze_support()
    
    morphologyPath = sys.argv[1] #'morphology_files/stick.hoc'
    upsideDown = sys.argv[2]
    if upsideDown == 'True':
        upsideDown = True
    elif upsideDown == 'False':
        upsideDown = False
    else:
        raise Exception('The second input should be "True" or "False"')
    resultFolder = sys.argv[3] #'test'
    if not os.path.isdir(resultFolder): 
        os.mkdir(resultFolder)
    # build the cell
    synPos = get_inserted_synPos(morphologyPath,soma_z=0,upsideDown) # to avoid some trouble in following function   
    # insert the synapse, parallel computing 
    pool = multiprocessing.Pool(processes=32)
    
    index=range(0,np.shape(synPos)[0]) # syn number    
    pool.map(func_star, itertools.izip(index, itertools.repeat(morphologyPath),itertools.repeat(upsideDown),itertools.repeat(resultFolder)))
    
    pool.close()
    pool.join()
    
    sum_over_same_z_zone(synPos[:,2],resultFolder,soma_to_bottom_MEA=150)
    
    end=time.time()
    print(end-start)