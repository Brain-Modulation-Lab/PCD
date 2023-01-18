function cfg =  simulation_settings(sim_i,n_sim,n_sources)
%% Settings for hybrid simulation (Realistic Scenario)
% 
% Run this script to set the parameters for the Realistic Simulation
% Scenario. 

% Here we show an example on a given trial saved at /Data folder
% ------------------------------------------

% REFERENCES
% [1] S.Cavallari, S. Panzeri and A.Mazzoni (2014) Comparison of the
% dynamics of neural interactions between current-based and conductance-based 
% integrate-and-fire recurrent networks, Frontiers in Neural Circuits
% 8:12. doi: 10.3389/fncir.2014.00012

% We used the following implementation: 
% https://senselab.med.yale.edu/ModelDB/showmodel.cshtml?model=152539&file=/CavallariEtAl2014__7_2015/#tabs-1

%
%author: mvissani
%last version: Jan 2023

% see also 
%--------------------------------------------------------------------------


% create config struct
cfg = [];

% parameter of simulation
cfg.sim.simulation_length = 9000; % units: (ms)
cfg.sim.fs_sim = 1/0.05*1000;
cfg.sim.npoints = 0.3*cfg.sim.fs_sim; % for spectrum
cfg.sim.fs_new = 1000;
cfg.sim.n_sim = n_sim; % n trials

% parameter of sources
cfg.src.n_sources = n_sources; % because with the artefact it becomes 100


% constant input
intensity = 10 + 4*rand(1,cfg.src.n_sources);
cfg.src.ext.mus = 5000 + 50*randn(1,cfg.src.n_sources);
cfg.src.ext.FWHM = 1200 + 150*randn(1,cfg.src.n_sources);
cfg.src.ext.intensity = intensity;

% put stimulus

% periodic modulation
A = max(intensity)/3*randn(1,cfg.src.n_sources);
cfg.src.ext.phi = 2*pi*rand(1,cfg.src.n_sources);
cfg.src.ext.fl = 10 + 1.5*randn(1,cfg.src.n_sources);
% OU tau
cfg.src.ext.taun = 0.16;

% set random seeds
cfg.sim.SEED_OU = sim_i*12000 + [1:cfg.src.n_sources];
cfg.sim.SEED_connections = cfg.sim.SEED_OU + cfg.src.n_sources;
cfg.sim.SEED_poisson = cfg.sim.SEED_connections + cfg.src.n_sources;