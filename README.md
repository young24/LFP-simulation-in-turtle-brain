# LFP-simulation-in-turtle-brain
This project was my bachelor thesis project, however, I have never backed them up until now since I did not use GitHub at that time.

With the help from Mark Shein-Idelson and Torbjorn, I was building a compartmental model to reproduce the local field potential (LFP) pattern in the Layer 3 (L3) with inputs from L1, L2 and L3. We were using an electrophysioligical toolkit called [LFPy](https://github.com/LFPy/LFPy) which relies on the [NEURON](http://www.neuron.yale.edu/neuron/) simulator. 

We assumed that the main input layer is the L2 since it has so many neurons, and directly projects to the L3. Further, we assumed that the L2 is filled with pyramidal neurons with identical typical morphologies which were measured from the experiment. To make the network balanced, we also took account the inhibitory inputs from the L1 and the L3.

Finally, we successfully reproduced some electrophysiological traces which matched well with the experimental data. However, due to the lack of knowledge in reptile visual cortex, and the lack of time, I did not make further contributions on the reverse problem which was using measured LFP patterns to determine neuronal morphologies.

The further development and results are in the paper: _[Large-scale mapping of cortical synaptic projections with extracellular electrode arrays](https://www.ncbi.nlm.nih.gov/pubmed/28805794)_
