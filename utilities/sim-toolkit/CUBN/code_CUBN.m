% CODE_CUBN  Simulate a recurrent network with excitatory and inhibitory
%            Leaky Integrate-and-Fire (LIF) neurons with current-based
%            synapses. The network has random connectivity and the input 
%            to each neuron is given by a Poisson process with a 
%            time-varying rate identical for all the neurons. 

%
%   [E2EI,I2EI,eFR,iFR] = CODE_CUBN(NET,INPUT2E,INPUT2I,SEED1,SEED2),
%   where:  
%   NET: external structure with the network's parameters 
%   (see e.g. parameters_CUBN script).
%   INPUT2E and INPUT2I are length M vectors with the rate of the external 
%   input (signal+noise) on excitatory (INPUT2E) and inhibitory (INPUT2I)
%   neurons in each time step of the simulation, Dt, 
%   in units of [(spikes/ms)*Dt] = [spikes/time_step]. Indeed Dt is a field
%   of the NET structure whose value is in units of [ms/time_step].
%   Note that the lengths of the vectors INPUT2E and INPUT2I have to be 
%   equal and gives the simulation duration in units of time steps Dt.
%   SEED1: seed for the network configuration. It must be a positive 
%   integer. If seed1=0, its value will be ignored 
%   (i.e., seed1=0 corresponds to not passing seed1)
%   SEED2: seed for the Poisson variate generator, which determines
%   the external input to each neuron. It must be a positive integer.
%  
%  
%   Returns
%  
%   E2EI and I2EI are M length vectors with the sum of the (E2EI) AMPA 
%   currents (both recurrent AMPA and external AMPA) and (I2EI) GABA 
%   currents entering all the excitatory neurons
%   in each time step of the simulation, in units of [pA]. 
%   eFR and iFR are M length vectors with the number of spikes fired in
%   each time step of the simulation from all the (eFR)excitatory and 
%   (iFR) inhibitory neurons.
%
%  
%   Note that in our model LFP = eRm*(I2EI-E2EI), where eRm is the membrane
%   resistance of the excitatory neurons.

