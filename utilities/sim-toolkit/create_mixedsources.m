function H = create_mixedsources(D,audio,A, epochs, n_trialsmax, n_channels, SNR_t)
% apply the mixing matrix A to D (neural sources) and the audio.
% Mixing operation: X = D*A
% 
% Output:
%       H: struct that stores all the mixing information
%       fields:
%              - time [1 x Nsamples]: array of timestamps
%              label {1 x n_channels}: cell array of electrode labels.
%              - trial {1 x n_trialsmax}: cell array of [n_channels x
%              Nsamples] (noisy, A mixing matrix) neural signal. 
%              - lfp {1 x n_trialsmax}: cell array of [n_channels x
%              Nsamples] (noiseless, A_gt mixing matrix) neural signal. 
%              - sources {1 x n_trialsmax}: cell array of [(n_sources + 1) x
%              Nsamples] neural & audio sources.
%              - SNR_e (1 x n_trialsmax):  array of empirical SNRs. They do
%              need be approximately equal to SNR_t (theoretical SNR).

% create noiseless matrix
A_gt = zeros(size(A));
A_gt(1:end-1,:) = A(1:end-1,:);

% epoch audio signal
cfg = [];
cfg.epoch = epochs(1 : n_trialsmax,:);
audio_s = bml_redefinetrial(cfg,audio);

% set up fieldtrip struct
H = struct();
H.time = audio_s.time;
H.label = cell(1,n_channels);

for chan_i = 1: n_channels
    H.label{chan_i} = ['Chan',num2str(chan_i)];
end

% mix signals
SNR_e = [];
for trial_i = 1 : n_trialsmax
    idx_prod = dsearchn(audio.time{1,trial_i}', [epochs.prod_onset(trial_i); epochs.prod_offset(trial_i)]);
    pow_audio = bandpower(audio_s.trial{1,trial_i}(idx_prod(1) : idx_prod(2)),1000,[70 180]);
    pow_gamma = bandpower(D.trial{1,trial_i}(:,idx_prod(1) : idx_prod(2)),1000,[70 180]);
    pow_gamma = mean(pow_gamma);
    k = pow_audio/pow_gamma*10^(-SNR_t(trial_i)/10);
    SNR_e = [SNR_e  10*log10(pow_audio/(k*pow_gamma))];

    sources = [sqrt(k)*D.trial{1,trial_i}'  audio_s.trial{1,trial_i}'];
    % perform source mixing
    H.trial{1,trial_i} = sources*A ;
    H.lfp{1,trial_i} = sources*A_gt ;

    H.sources{1,trial_i} = sources;
end
H.SNR_e = SNR_e;
%H.SNR_t = SNR_t;
H.A = A;
H.A_gt = A_gt;
H.epochs = epochs(1:n_trialsmax,:);