%% "Demo: comparative analysis PCD, CAR and ICA in denosing the vibration
% artifact"

% Here we show an example on a given trial saved at /Data folder
% ------------------------------------------

% REFERENCES
% [1] Bush, Alan, et al. "Differentiation of speech-induced artifacts from 
% physiological high gamma activity in intracranial recordings." NeuroImage
% 250 (2022): 118962.
%author: vpeterson
%last version: Jan 2023

% see also apply_PCD.m, fit_PCD.m, fit_ICA.m
%--------------------------------------------------------------------------
addpath(genpath('./utilities/'))
filesep = '\';

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
ze = Data.z; %is the recorded audio, used to drive the estimatation of 
% the artifactual source. These two signals should have the same number of
% data points
X_toclean_e = Data.X_toclean; %data to be denoised. Here a wider extrated
% epoch that contains Xe.
F0 = Data.F0; % Annotation of F0. Here we have one annotation per produced
% syllable.
sf = Data.sf; %sampling frequency. Here 1000 Hz.
ch_names = Data.ch_names; %channel names. Here we have data recorded from 
%ECoG and DBS electrodes. 
%% RUN PCD

% calculate the power spectrum of the audio (z)
npoints = 0.2*sf;
[pzz,f] = pwelch(detrend(ze,0),hanning(npoints),[],2^nextpow2(npoints)*5,sf);
% Find the peak around F0

% We have annotation of the F0 for every syllable. We use the mean across 
% the tripplet
fo_ = mean(F0);

% select only gamma band
idx_gamma = find(f>50&f<250);
pzz_gamma = pzz(idx_gamma);
freq_gamma = f(idx_gamma);

% find the peak at F0
idx_peak = find(freq_gamma>0.5*fo_);
pzz_peak = pzz_gamma(idx_peak);
freq_peak = freq_gamma(idx_peak);

% define starting points to define the center and bandwith of the vibration
% artifact frequency band
Xpk = fo_;
Wpk = 40;

% here we fit using Gaussians
[yhat, parameter] = fitGaussian(Xpk, Wpk, freq_peak, pzz_peak,1);
parameter = abs(parameter);
BWs = sqrt(2*log(2))*ceil(parameter(2)); % FWHM

%% Define SSD bands

% define the bands
fl = round(parameter(1) - ceil(BWs/2));
fh = round(parameter(1) + ceil(BWs/2));
bands = [fl,fh; 4,240;  fl-1,fh+1];

%check if bands definitions fullfile SSD requirements
signal_band = bands(1,:); % signal bandpass band
noise_bp_band = bands(2,:); % noise bandpass band
noise_bs_band = bands(3,:); % noise bandstop band
if not( noise_bs_band(1) < signal_band(1) && ...
    noise_bp_band(1) < noise_bs_band(1) && ...
    signal_band(2) < noise_bs_band(2) && ...
    noise_bs_band(2) < noise_bp_band(2) )
% manual def 
BWs = 40;
fl = round(fo_ - ceil(BWs/2));
fh = round(fo_ + ceil(BWs/2));
bands = [fl,fh; 4,240;  fl-1,fh+1];
end
%sometims estimated fo is greater than 240 hz, and thus we need 
%to redefine the bands
signal_band = bands(1,:); % signal bandpass band
noise_bp_band = bands(2,:); % noise bandpass band
noise_bs_band = bands(3,:); % noise bandstop band
if not( noise_bs_band(1) < signal_band(1) && ...
    noise_bp_band(1) < noise_bs_band(1) && ...
    signal_band(2) < noise_bs_band(2) && ...
    noise_bs_band(2) < noise_bp_band(2) )
% manual def 
fo_ = 120;
BWs = 40;
fl = round(fo_ - ceil(BWs/2));
fh = round(fo_ + ceil(BWs/2));
bands = [fl,fh; 4,240;  fl-1,fh+1];
end

fprintf('Running PCD')  

nc = size(Xe, 1);

%% fit PCD
ssd_components ='PR';
[W_pcd, A_pcd, params, X_ssd, X_pcd] = fit_PCD(Xe, ze, bands, sf, 'ssd_n_components',ssd_components, 'verbose', 1);
%% apply PCD
n_comp ='diff';        
[X_denoised_PCD, Proj_matrix_pcd, idx_keep_pcd] = apply_PCD(W_pcd,A_pcd,X_toclean_e,params.vlen,n_comp, 1);

%% RUN CAR
X_denoised_CAR = apply_CAR(Xe, 'spatialfilter');

%% RUN ICA
pca_components = 0.99;
[W_ica, A_ica, params, X_pca, X_ica] = fit_ICA(Xe, ze, 'pca_n_components',pca_components, 'verbose', 1);
[score, idx] = select_ica_components(ze, X_ica,sf, fo_, 'False');
n_comp ='diff';        
[X_denoised_ICA, Proj_matrix_ica, idx_keep_ica] = apply_ICA(W_ica,A_ica,X_toclean_e,score, idx,n_comp, 1);
%% Let's have a look at the estimated source in for ICA and CAR
% we take electrode 149, same as Figure 5 in [2].
Xe_49 = Xe(49,:);

time = 0:1/sf:length(ze)/sf-1/sf;
ax = figure('units','normalized','outerposition',[0 0 1 1], "renderer","painters");
% time-frequency
nwin = 100; nvlp = 80; 
fint = 50:2:250;  
[Bz,f2,T] = spectrogram(Xe_49, hamming(nwin), nvlp,fint,sf, 'power');
Bz = 20*log10(abs(Bz));
B_idx = Bz < max(max(Bz))-40;
Bz(B_idx) = 0;


subplot(2,4,1),imagesc(T, f2, Bz);
axis xy;
title('RAW');
colorbar
xlabel('Time [s]');ylabel('Frequency [Hz]')

%% denoised by CAR 
Xe_49_CAR = X_denoised_CAR(49,:); 


% time-frequency
nwin = 100; nvlp = 80; 
fint = 50:2:250;  
[Bz,f2,T] = spectrogram(Xe_49_CAR, hamming(nwin), nvlp,fint,sf, 'power');
Bz = 20*log10(abs(Bz));
B_idx = Bz < max(max(Bz))-40;
Bz(B_idx) = 0;

subplot(2,4,2),imagesc(T, f2, Bz);
axis xy;
title('CAR')
colorbar
xlabel('Time [s]');ylabel('Frequency [Hz]')
%% denoised by ICA 
Xe_49_ICA = X_denoised_ICA(49,:); 

% time-frequency
nwin = 100; nvlp = 80; 
fint = 50:2:250;  
[Bz,f2,T] = spectrogram(Xe_49_ICA, hamming(nwin), nvlp,fint,sf, 'power');
Bz = 20*log10(abs(Bz));
B_idx = Bz < max(max(Bz))-40;
Bz(B_idx) = 0;

subplot(2,4,3),imagesc(T, f2, Bz);
axis xy;
title('ICA')

colorbar
xlabel('Time [s]');ylabel('Frequency [Hz]')
%% denoised by PCD 
Xe_49_PCD = X_denoised_PCD(49,:); 

% time-frequency
nwin = 100; nvlp = 80; 
fint = 50:2:250;  
[Bz,f2,T] = spectrogram(Xe_49_PCD, hamming(nwin), nvlp,fint,sf, 'power');
Bz = 20*log10(abs(Bz));
B_idx = Bz < max(max(Bz))-40;
Bz(B_idx) = 0;

subplot(2,4,4),imagesc(T, f2, Bz);
axis xy;
title('PCD')
colorbar
xlabel('Time [s]');ylabel('Frequency [Hz]')