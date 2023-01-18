%% Run a demo for the simulation of the Realistic Scenario
% Neural sources (D) [following [1]] and a speech-produced audio (A) [2] (in Data folder) are
% mixed to create noisy neural signal (H = D*A)
%
% The function code_CUBN.c is a mex source code. You have to compile this routine in Matlab to generate the 
% mex file (e.g.: code_CUBN.mexw64). Note that you have to include the functions ran1.c and gasdev.c in the 
% compiling instruction in the Matlab workspace, in the following way: 
% mex code_CUBN.c ran1.c gasdev.c 
%
% References:
% [1] S.Cavallari, S. Panzeri and A.Mazzoni (2014) Comparison of the
% dynamics of neural interactions between current-based and conductance-based 
% integrate-and-fire recurrent networks, Frontiers in Neural Circuits
% 8:12. doi: 10.3389/fncir.2014.00012
% [2] Bush A, Chrabaszcz A, Peterson V, Saravanan V, Dastolfo-Hromack C, 
% Lipski WJ, Richardson RM. Differentiation of speech-induced 
% artifacts from physiological high gamma activity in intracranial 
% recordings. Neuroimage. 2022 Apr 15;250:118962. 

% see also CUBN folder in utilities (& https://senselab.med.yale.edu/ModelDB
% /showmodel.cshtml?model=152539&file=/CavallariEtAl2014__7_2015/LIF_COBN/
% for more information)

%
%author: mvissani
%last version: Jan 2023


%% Step 0: Set up the different folders
% include the folder (and subfolders) utilies in the path
addpath(genpath(['.' filesep 'utilities' filesep]))

% set random generator seed
randn('state',42);

% load exemplary audio
PATH_DATA = './Data';

%% Step 1: Simulate the neural signals (H = D*A)

% settings of simulation
n_sim = 100; % number of trials
n_sources = 99; % number of sources
n_channels = 60; % 100 number of channels

% get audio from ./Data folder
[audio, epochs] = get_audio(PATH_DATA);

% generate n_sources neural sources with n_sim trials of activity
[D, cfg_all] = generate_sources(n_sim, n_sources);

% parameters for mixing (A: noisy matrix & Agt: noiseless matrix)
A = generate_mixingmatrix(n_sources,n_channels);

% define (theoretical) SNR levels
SNR = -5:5; % -2:10
SNR_t = datasample(SNR,n_sim);

% mix neural sources and audio signals to create the (noisy) neural signals
H = create_mixedsources(D, audio, A, epochs, n_sim, n_channels, SNR_t);


