function [unmixing, mixing, params, X_pca, X_icapca] = fit_ICA(X, z, varargin)
% Implements ICA for artifact rejection. This version uses PCA to
% perform a dimentionality reduction step before running ICA.
%  
%     Input:
%     ------
%     z - (array, floats > 0) - amplitudes [,n_samples]
%     X - (array, complex) - raw data
%         [n_channels, n_samples]
%   
%    'pca_n_components' - (int > 0) - determines the number of components
%                       for the dim_reducction based on ssd
%                       defaults to 20.
%     'verbose' -  boolean
%     Output:
%     -------
%     unmixing
%     mixing
%     vlen - numpy array - the length of the mean vector for each filter
%     X_pca
%     X_ica
%    
%%
if not(isempty(X))
    [n_channels, n_samples] = size(X);
end
if not(length(z') == n_samples)
    error('X and z must have the same number of epochs!!!');
end
%% params
opt= propertylist2struct(varargin{:});
opt= set_defaults(opt, ...
    'verbose', 1, ...
    'pca_n_components', 40);
% define dic of params used
params = [];

if opt.verbose > 0
    fprintf('\n--- Begin ICA analysis ---\n')
end


%% some preprocessing
% check Nan
z(isnan(z))=0;
X(isnan(X))=0;
% save data before Hilbert transform
Xr = X;
% normalize the target funtion to have zero mean and unit-variance and
% average over all epochs
z = (z-mean(z(:)))./std(z(:));
% pre-whiten data
% %center data
% X = X - repmat(mean(X), size(X,1), 1);

X = prepare_data(X');
% params.W_white = M;
% apply PCA
if opt.verbose > 0
    fprintf('ruuning PCA\n')
end
[W_pca,X_pca,lambda_pca, ~, explained_var] = pca(X,  'Economy', false);

if floor(opt.pca_n_components) ~= ceil(opt.pca_n_components)
    %is a decimal accounting for the threshold value
    mont_curve = cumsum(explained_var)/max(cumsum(explained_var));
    [~, pca_n_components] = min(abs(mont_curve-opt.pca_n_components));
elseif strcmp(opt.pca_n_components, 'PR')
    % participation ratio
    PR = (sum(explained_var)^2)/(sum(explained_var.^2));
    pca_n_components = round(PR);
else
    %is an integer, accounting for the number of components to keep
    pca_n_components = opt.pca_n_components;
end
% save ssd parameters
params.pca_n_components = pca_n_components;
params.lambda_pca = lambda_pca;
params.W_pca = W_pca;
%this is the low-dim raw signal in which we want to compute ICA
picks_pca = W_pca(:,1:pca_n_components); %columns are filters
% %apply SSD filters to the data
% X_ssd = (real(X)' * W_ssd)';
% X = X_pca(:, 1:pca_n_components)
X = X * picks_pca;
% X = X';
params.X_pca = X_pca;
% params.Xpre = X;
params.picks_pca = picks_pca;
%get the final number of filters
num = min([n_channels, size(X,2)]);

if opt.verbose > 0
    fprintf('Using the following number of dimensions = ')
    fprintf('%d ', num);
    fprintf('\n')
end

%% ica
% Runs the Picard algorithm for ICA.

[X_ica, W_ica] = picard(X', 'whiten', false);

% save ICA params
params.W_ica = W_ica;
params.X_ica = X_ica;

%% post processing
% Note that the spatial filters are now defined in PCA space and have 
% to be mapped to sensor space to be applicable to the original data.
W = picks_pca * W_ica;
% % project back to original (un-whitened) channel space
% W = real(M) * W; 

%% apply PCD filters to the data
if nargout > 4
    X_icapca = W'*Xr;
end
%% calcule the unmixing and mixing matrices
unmixing = eye(n_channels, n_channels);
unmixing(1:pca_n_components, 1:pca_n_components) = W_ica;
unmixing = W_pca*unmixing; % columns are filters
% % project back to original (un-whitened) channel space
% unmixing = real(M) * unmixing;  
% calculate the mixing matrix
mixing = pinv(unmixing);

end
%% prepare_data
function X_prewhite = prepare_data(X)
% sanity checks and hilbert transform
if isempty(X)
    error('X must be non-empty!');
end
% now size(X) = [time, channels]!
X_prewhite = zscore(X);

end


