%% "Run the inter-trial phase consistency analysis"
% Run this script to compute the inter-trial phase consistency (ITPC) [1]
% following the framework proposed by Bush et al [2]. This metric
% quantifies the contamination of the speech-induced artifact in the neural
% signal. The signal is considered contaminated if ITPC >= 3.08, where 3.08
% is the permutation-based significance level.

% Here we show an example on a given trial saved at /Data folder
% ------------------------------------------

% REFERENCES
% [1] Cohen, M. X. Analyzing Neural Time Series Data: Theory and Practice. (MIT Press, 2014).
% [2] Bush A, Chrabaszcz A, Peterson V, Saravanan V, Dastolfo-Hromack C,
% Lipski WJ, Richardson RM. Differentiation of speech-induced
% artifacts from physiological high gamma activity in intracranial
% recordings. Neuroimage. 2022 Apr 15;250:118962.

%
%author: mvissani
%last version: Jan 2023

% see also get_coherence.m
%--------------------------------------------------------------------------
addpath(genpath(['.' filesep 'utilities' filesep]))

%% Step 0: Loading data in Data
% load exemplary data
load('./Data/DataTrials.mat')
% This data contains 20 trials of 2 electrodes of a given patient data from the dataset
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
X_toclean_e = Data.X_toclean; %data to be dennChansChansnoised. Here a wider extrated
% epoch that contains Xe.
F0 = Data.F0; % Annotation of F0. Here we have one annotation per produced
% syllable.
sf = Data.sf; %sampling frequency. Here 1000 Hz.
ch_names = Data.ch_names; %channel names. Here we have data recorded from
%ECoG and DBS electrodes.

nChans = numel(DataTrials{1,1}.ch_names); % number of channels (e.g., 2)
nTrials = numel(DataTrials); % number of trials (e.g., 20)

%% Step 1: Calculate ITPC before applying PCD

% extract ITPC
CoH = nan(nTrials,nChans);
for trial_i = 1 : nTrials
    X_toclean_ = DataTrials{1,trial_i}.X_toclean; % data to clean
    Ze_ = DataTrials{1,trial_i}.z; % speech-produced audio
    CoH(trial_i,:) = get_coherence(X_toclean_, Ze_);
end
[ReCoH, ImCoH, ITPC] = CoH2ITPC(CoH);

%% Step2: Estimate the vibration artifact bandwith and apply SSD


for trial_i = 1 : nTrials
    % get data & metadata in each trial
    X_toclean_ = DataTrials{1,trial_i}.X_toclean; % data to clean
    Ze_ = DataTrials{1,trial_i}.z; % speech-produced audio
    Xe_ = DataTrials{1,trial_i}.Xtofit; % data to clean used to fit
    sf = DataTrials{1,trial_i}.sf; % sample frequency
    % We have annotation of the F0 for every syllable. We use the mean across
    % the tripplet
    F0 = mean(DataTrials{1,trial_i}.F0,'omitnan'); % fundamental frequency

    % calculate the power spectrum of the audio (z)
    npoints = 0.2*sf;
    [pzz,f] = pwelch(detrend(Ze_,0),hanning(npoints),[],2^nextpow2(npoints)*5,sf);
    % Find the peak around F0

    % select only gamma band
    idx_gamma = find(f>50&f<250);
    pzz_gamma = pzz(idx_gamma);
    freq_gamma = f(idx_gamma);

    % find the peak at F0
    idx_peak = find(freq_gamma>0.5*F0);
    pzz_peak = pzz_gamma(idx_peak);
    freq_peak = freq_gamma(idx_peak);

    % define starting points to define the center and bandwith of the vibration
    % artifact frequency band
    Xpk = F0;
    Wpk = 40;

    % here we fit using Gaussians
    [yhat, parameter] = fitGaussian(Xpk, Wpk, freq_peak, pzz_peak,1);
    parameter = abs(parameter);
    BWs = sqrt(2*log(2))*ceil(parameter(2)); % FWHM
    % Define SSD bands

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
        fl = round(F0 - ceil(BWs/2));
        fh = round(F0 + ceil(BWs/2));
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

    % fit PCD
    ssd_components ='PR';
    [W_pcd, A_pcd, params, X_ssd, X_pcd] = fit_PCD(Xe_, Ze_, bands, sf, 'ssd_n_components',ssd_components, 'verbose', 1);
    % apply PCD
    n_comp ='diff';
    [X_clean, Proj_matrix, idx_keep] = apply_PCD(W_pcd,A_pcd,X_toclean_,params.vlen,n_comp, 1);
    
    % save in same struct results
    DataTrials{1,trial_i}.X_cleaned = X_clean;
    DataTrials{1,trial_i}.Proj_matrix = Proj_matrix;
    DataTrials{1,trial_i}.idx_keep = idx_keep;
end

%% Step 2: run ITPC on clean data

% extractITPC
CoH = nan(nTrials,nChans);
for trial_i = 1 : nTrials
    X_cleaned_ = DataTrials{1,trial_i}.X_cleaned; % data to clean
    Ze_ = DataTrials{1,trial_i}.z; % speech-produced audio
    CoH(trial_i,:) = get_coherence(X_cleaned_, Ze_);
end
[ReCoH_after, ImCoH_after, ITPC_after] = CoH2ITPC(CoH);

%%
% plot the distribution of the ITPC values
figure('defaultAxesFontSize',14)
tiledlayout(1,3)
nexttile(1,[1 2])
plot(repmat([1; 2],1, nChans), [ITPC;  ITPC_after],'-ok','linewidth',1.5)
xticks([1 2])
xticklabels({'Before PCD','After PCD'})
ylabel(' ITPC ')
xlim([.5 2.5])
box off
nexttile(3,[1 1])
scatter(.2*randn(1,nChans) + 1 , ITPC - ITPC_after,35,'k','filled')
ylabel(' ITPC reduction (Before - After)')
xlim([.5 1.5])
box off

% plot the Real(CoH) and Im(CoH) on the polar plot
figure('defaultAxesFontSize',14)
quiver(ReCoH, ImCoH,ReCoH_after - ReCoH, ImCoH_after - ImCoH ,'r','LineWidth',1.5)
% hold on
% scatter(ReCoH_after,ImCoH_after,15,'r','filled')
xlabel(' Re[ITPC]')
ylabel(' Imag[ITPC]')
xline(0,'--')
yline(0,'--')
axis tight
box off


