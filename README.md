# LFP-simulation-in-turtle-brain
Code to simulate the local field potential in turtle's visual cortex.
This project was my bachelor thesis project, however, I have never backed them up until now since I did not use GitHub at that time.
I was using an electrophysioligical toolkit [LFPy](https://github.com/LFPy/LFPy) which relies on the [NEURON](http://www.neuron.yale.edu/neuron/) simulator. 
I was building an electrophysioligical model to reproduce the LFP pattern on the L3 with inputs from L1 (inhibitory kernel), L2 (pyramidal kernel) and L3 (inhibitory kernel). 

Assumptions:
- Assume that the L2 is isotropical filled with a typical pyramidal neuron.
- Make a kernel of the potential field around the soma of specific neurons.
- Then use the interpolation method to superimpose all influence of all synapses.

Significance: 
- If succeed, which means we can use the mechanism underlying the simulation to explain the formation of LFP pattern
- Predicting power: try to implement the active membrane mechanism (selective mechanism?), to include some high filtering effect
- Inverse problem: try different morphologies to get different LFP patterns; inversely use LFP patterns to determine neuronal morphology.


Some memory are blurred to me because first it is over 2 years, second, to be honest, I knew little about neuroscience at that time, especially for the local field potential. But that experience was precious to me since it was the first time I tasted the good science. From that time on, I am getting the sense of what a good research should be like. Gilles Laurent as my advisor impressed me a lot for his sharpness in science and his insistence in quantification.
