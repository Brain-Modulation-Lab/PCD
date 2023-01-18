%% "Run the phase-coupling decomposition method"
% Raw data, which is contaminated with the vibration artifact [1] is
% clean via a phase coupling decompoisition method (PCD). The PCD method
% is based on spatial filterings approaches, namely SSD [2] and PCO  [3]
% (see fit_PCD.m function on the utilities folder). 
% Data is clean in a trial basis. For fitting the model narrow epochs are
% used. For applying the model, and thus, saving the clean data, a broader
% epoch is defined (6 sec of buffer from the beginning and the end of the
% produced epoch). 

% Here we show an example on a given trial saved at /Data folder
% ------------------------------------------

% REFERENCES
% [1] Bush, Alan, et al. "Differentiation of speech-induced artifacts from 
% physiological high gamma activity in intracranial recordings." NeuroImage
% 250 (2022): 118962.
% [2] Nikulin, Vadim V., Guido Nolte, and Gabriel Curio.
% "A novel method for reliable and fast extraction of neuronal EEG/MEG 
% oscillations on the basis of spatio-spectral decomposition." NeuroImage 
% 55.4 (2011): 1528-1535.
% [3] Waterstraat, Gunnar, Gabriel Curio, and Vadim V. Nikulin. 
% "On optimal spatial filtering for the detection of phase coupling in 
% multivariate neural recordings." NeuroImage 157 (2017): 331-340.
%
%author: vpeterson
%last version: Jan 2023

% see also apply_PCD.m, fit_PCD.m
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
%% Estimate the vibration artifact bandwith

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
[X_clean, Proj_matrix, idx_keep] = apply_PCD(W_pcd,A_pcd,X_toclean_e,params.vlen,n_comp, 1);

%% plot audio and first retrieve source
time = 0:1/sf:length(ze)/sf-1/sf;
% plot true audio source

ax = figure('units','normalized','outerposition',[0 0 1 1], "renderer","painters");
% time
subplot(2,3,1),plot(time, ze, 'color', rgb('Blue'), 'linewidth',1);
box off
xlabel(" Time [s] ")
ylabel( " Amplitude ")
axis tight

% frequency
[pzz,f] = pwelch(detrend(ze,0),hanning(npoints),[],2^nextpow2(npoints)*5,sf);
subplot(2,3,2), semilogy(f,pzz,'linestyle','-','color',rgb('MediumBlue'),'linewidth',1.4)
grid on
grid minor
xlabel(" Frequency [Hz] ")
ylabel( " Power [dB] ")
axis tight
title("True artifact source")
xlim([20 250])
box off

% time-frequency
nwin = 100; nvlp = 80; 
fint = 2:2:250;  
[Bz,f2,T] = spectrogram(ze, hamming(nwin), nvlp,fint,sf, 'power');
Bz = 20*log10(abs(Bz));
B_idx = Bz < max(max(Bz))-60;
Bz(B_idx) = 0;

subplot(2,3,3),imagesc(T, f2, Bz);
axis xy;

colorbar
xlabel('Time [s]');ylabel('Frequency [Hz]')

%% estimated source
ze_est = X_pcd(1,:);
% time
subplot(2,3,4), plot(time, ze_est, 'color', rgb('DarkSlateGray'), 'linewidth',1);
box off
xlabel(" Time [s] ")
ylabel( " Amplitude ")
axis tight

% frequency
[pzz_est,f] = pwelch(detrend(ze_est,0),hanning(npoints),[],2^nextpow2(npoints)*5,sf);
subplot(2,3,5),semilogy(f,pzz_est,'linestyle','-','color',rgb('DarkSlateGray'),'linewidth',1.4)
grid on
grid minor
xlabel(" Frequency [Hz] ")
ylabel( " Power [dB] ")
axis tight
title("Estimated artifact source")
xlim([20 250])
box off

% time-frequency
nwin = 100; nvlp = 80; 
fint = 2:2:250;  
[Bz,f2,T] = spectrogram(ze_est, hamming(nwin), nvlp,fint,sf, 'power');
Bz = 20*log10(abs(Bz));
B_idx = Bz < max(max(Bz))-60;
Bz(B_idx) = 0;

subplot(2,3,6),imagesc(T, f2, Bz);
axis xy;
colorbar
xlabel('Time [s]');ylabel('Frequency [Hz]')