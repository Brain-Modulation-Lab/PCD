%% "Run the inter-trial phase consistency analysis"
% Raw this script to compute the inter-trial phase consistency (ITPC) [1]
% following the framework proposed by Bush et al [2]. This metric
% quantifies the contamination of the speech-induced artifact in the neural
% signal. The signal is considered contaminated if ITPC >= 3.08, where 3.08
% is the permutation-based significance level.

% Here we show an example on a given trial saved at /Data folder
% ------------------------------------------

% REFERENCES
% [1] Cohen, M. X. Analyzing Neural Time Series Data: Theory and Practice. (MIT Press, 2014). 
% [2] Bush, Alan, et al. "Differentiation of speech-induced artifacts from 
% physiological high gamma activity in intracranial recordings." 
% bioRxiv (2021).

%
%author: mvissani
%last version: Jan 2023

% see also get_coherence.m
%--------------------------------------------------------------------------
addpath(genpath(['.' filesep 'utilities' filesep]))

% load exemplary data
load('./Data/Data.mat')
% This data contains a trial of a given patient data from the dataset 
% described in [1]. Patients were instructed to repeat consonant-vowel 
% syllable triplets that were played on their earphones. iEEG recordings 
% were band-pass filtered between 2 and 250 Hz and notch filtered at 60Hz
% and its 3 first harmonycs. We use use a narrow epoch around the produced
% audio to fit PCD. Once PCD is fit, and the linear transformation matrices
% are learned, we can apply the model in a wider epoch. We follow this
% procedure here.

Xe = Data.X_tofit; %is an array of n_channels x n_sampling points from 
% which the PCD will be learned
Ze = Data.z; %is the recorded audio, used to drive the estimatation of 
% the artifactual source. These two signals should have the same number of
% data points
X_toclean_e = Data.X_toclean; %data to be denoised. Here a wider extrated
% epoch that contains Xe.
F0 = Data.F0; % Annotation of F0. Here we have one annotation per produced
% syllable.
sf = Data.sf; %sampling frequency. Here 1000 Hz.
ch_names = Data.ch_names; %channel names. Here we have data recorded from 
%ECoG and DBS electrodes. 

nchans = numel(ch_names); % number of channels

% get coherence
CoH = nan(1, nchans);
for ch = 1 : nchans
    CoH(ch) = get_coherence(X_toclean_e(ch,:),Ze);  
end

% get abs(CoH) to calculate the ITPC
ITPC = abs(CoH);

% get Real(CoH) and Im(CoH) to plot on the complex space the coherence
ReCoH = real(CoH);
ImagCoH = imag(CoH);

% plot the distribution of the ITPC values
figure('defaultAxesFontSize',12)
histogram(ITPC,floor(sqrt(nchans)))
xlabel(' ITCP [a.u.] ')
ylabel(' Counts # ')
box off

% plot the Real(CoH) and Im(CoH) on the polar plot
figure('defaultAxesFontSize',12)
scatter(ReCoH,ImagCoH,15,'k','filled')
xlabel(' Re[ITPC]')
ylabel(' Imag[ITPC]')
xline(0,'--')
yline(0,'--')
