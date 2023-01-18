% Conventions for variable names:
% - prefix "N" stands for "number of..."
%
% - prefix "e" stands for "excitatory"
% - prefix "i" stands for "inhibitory"
% - prefix "x" stands for "external"
%
% - prefix "e2e" stands for "excitatory to excitatory"
% - prefix "e2i" stands for "excitatory to inhibitory"
% - prefix "x2e" stands for "external to excitatory"
% - prefix "e2i" stands for "external to inhibitory"
%
% - "T" stands for greek letter "tau"



% THE NETWORK

% Time resolution [ms]
net_CUBN.Dt = 0.05;


% NETWORK PROPERTIES -----------------------------------------------------

% Number of excitatory and inhibitory neurons, eNnrn and iNnrn, and total number
% of neurons, totNnrn:
net_CUBN.eNnrn = 4000;
net_CUBN.iNnrn = 1000;

% Connectivity, i.e., connection probability, p:
net_CUBN.p = 0.2; 

% Membrane time constant, [ms]
net_CUBN.eTm = 20; 
net_CUBN.iTm = 10; 

% Leak membrane potential, [mV]
net_CUBN.V_leaky = -70.; 

% Membrane potential firing threshold, [mV]
net_CUBN.Vthr = -52; 
% the neuron fires following the following sequence:
% 1. the neuron potential is reset to a value Vres, [mV]: 
net_CUBN.eVres = -59; 
net_CUBN.iVres = -59; 
% 2. the neuron cannot fire again for a refractory period, Trp, equal to
%    2ms for E neurons and 1ms for I neurons, [ms]:
net_CUBN.eTrp = 2; 
net_CUBN.iTrp = 1; 
% 3. all post-synaptic neurons receive a spike with a delay, Tl, equal to
%    of 1ms after the time of threshold crossing, [ms]:
net_CUBN.eTl = 1;
net_CUBN.iTl = 1; 

% Rise and deacy times, Tr and Td [ms]:
% - of E => E synaptic currents:
net_CUBN.e2eTr = 0.4;
net_CUBN.e2eTd = 2.; 
% - of E => I synaptic currents:
net_CUBN.e2iTr = 0.2; 
net_CUBN.e2iTd = 1.; 
% - of I synaptic currents: rise and decay times are assumed idependent of
%   the type of neuron they act on (excitatory or inhibitory)
net_CUBN.iTr = 0.25; 
net_CUBN.iTd = 5.; 

% Synaptic reversal potentials, [mV]
net_CUBN.VsynAMPA = 0; 
net_CUBN.VsynGABA = -80;

% Synaptic efficacies [pA]
% - on inhibitory interneurons:
net_CUBN.i2iJ = 54; 
net_CUBN.e2iJ = -14; 
net_CUBN.x2iJ = -19; 
% - on pyramidal neurons:
net_CUBN.i2eJ = 42.5; 
net_CUBN.e2eJ = -10.5; 
net_CUBN.x2eJ = -13.75; 

% Membrane resistances, [GOhm]
net_CUBN.eRm = 0.04; 
net_CUBN.iRm = 0.05; 

