function [ReCoH, ImCoH, ITPC] = CoH2ITPC(CoH)
% Calculate IPTC normalizing the coherence values across trials
% Inputs:
%       CoH (trials x nchans): array of coherence values


ntrials = size(CoH,1);

if ntrials > 1
    mRNC = mean(real(CoH),1,'omitnan');
    mINC = mean(imag(CoH),1,'omitnan');
    seRNC = se(real(CoH));
    seINC = se(imag(CoH));
    ReCoH = mRNC./(sqrt(seRNC + seINC)/2);
    ImCoH = mINC./(sqrt(seRNC + seINC)/2);
    ITPC = sqrt(mRNC.^2 + mINC.^2)./sqrt(seRNC.^2 + seINC.^2);
else
    % no normalization is possible; just report Re[CoH_trials] and
    % Im[Coh_trials]
    ReCoH = real(CoH);
    ImCoH = imag(CoH);
    ITPC = abs(CoH);
end

end
% helper function for standard error of mean
function y = se(x)
    y = std(x,[],1,'omitnan')./sqrt(sum(~isnan(x),1));
end