% single_synaptic_input_MATLAB
% start-script to start single_synaptic_input_MATLAB.py
morphologyPath = 'morphology_files/reconstructed/W66N1.hoc';
resultFolder = 'sim_results';
soma_z = 150;  % z-position of cell, relative to MEA plane at 0 um
syn_pos = [0, 0, soma_z+250]; %Synapse is inserted at closest cellular position to specified point (x,y,z)
syn_weight = 0.05; %nA; for active membrane: 0.05; for passive membrane: 0.005
inputSign = 0; %0 for excitatory synapse, -90 for inhibitory
isPassive = 0; % passive membrane: True; active membrane: False
upsideDown = 1; % if morphology is upsideDown

%evalStr='MoIKCSD(ele_pos, pots, gdx=0.05, gdy=0.05, xmin=-2.0, xmax=2.0,ymin=-2.0, ymax= 2.0)';
save('inputFor_single_synaptic_input_MATLAB.mat','morphologyPath','resultFolder',...
'soma_z','syn_pos','syn_weight','inputSign','isPassive','upsideDown');

commandStr = ['python',' ','single_synaptic_input_MATLAB.py'];
[status, commandOut]=system(commandStr);
