function [unmixing, mixing, params, X_ssd, X_pcd] = fit_PCD(X, z, bands, sampling_freq, varargin)
% Implements Phase coupling descomposition(PCD), a method for extracting
% vibration artifact induced by speech audio. This version uses SSD to
% perform a dimentionality reduction step before the denoising via the
% phase-coupling optimization (PCO) method.
% 
% The PCO implementation is based on the PCO-MEET Python Library
% https://github.com/neurophysics/meet
% The SSD implementation is based on
% https://github.com/svendaehne/matlab_SPoC
% author: Victoria Peterson
%--------------------------------------------------------------------------
% Phase Amplitude Coupling Decomposition.
% --------------------------------------
%     It maximizes the length of the "mean vector" and
%     returns the filter coefficients w and a denoised version of the
%     artifactual input data.
% 
%     Input:
%     ------
%     z - (array, floats > 0) - amplitudes [,n_samples]
%     X - (array, real) - signal to be denoised,
%         [n_channels, n_samples]
%     bands -  as in ssd 
%             3 x 2 matrix with the cut-off frequencies. 
%             First row: cut-off frequencies for band-pass of the to be extracted 
%             oscillations.
%             Second row: cut-off frequencies for the lowest and highest 
%             frequencies defining flanking intervals.
%             Third row: cut-off frequencies for the band-stop filtering of 
%             the central frequency process.
%    sampling_freq - for computing the filtering
%    'filter_order' - in a butterworth filtering
%    'bestof' - (int > 0) - the number of restarts for the optimization of the
%                        individual filters. The best filter over all
%                        these restarts with random initializations is
%                        chosen, defaults to 15.
%    'ssd_n_components' - (int > 0) - determines the number of components
%                       for the dim_reducction based on ssd
%                       defaults to 20.
%     'verbose' -  boolean
%     Output:
%     -------
%     unmixing
%     mixing
%     vlen - numpy array - the length of the mean vector for each filter
%     X_ssd
%     X_pcd
%     W - numpy array - the filters for Y, each filter is in a column of Wy
%                        (if num==1: Wy is 1d)
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
    'filter_order',5,...
    'verbose', 1, ...
    'bestof', 15, ... 
    'ssd_n_components', 40);
% define dic of params used
params = [];

if opt.verbose > 0
    fprintf('\n--- Begin PCD analysis ---\n')
end

%for the first filter
function [obj, grad] = objfun_x(w)
    [obj, grad] = PCD_obj(w, z, X);
end
%for subsequent filters
function [obj, grad] = objfun_bx(w)
    [obj, grad] = PCD_obj(w, z, bx);
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

% whiten the data
[~, X, M] = prepare_data(X);
params.W_white = M;
% apply SSD
if opt.verbose > 0
    fprintf('ruuning ssd\n')
end
[W_ssd, ~, lambda_ssd, ~]  = ssd(real(X)', bands, sampling_freq, opt.filter_order,[]); 

if floor(opt.ssd_n_components) ~= ceil(opt.ssd_n_components)
    %is a decimal accounting for the threshold value
    mont_curve = cumsum(lambda_ssd)/max(cumsum(lambda_ssd));
    [~, ssd_n_components] = min(abs(mont_curve-opt.ssd_n_components));
elseif strcmp(opt.ssd_n_components, 'PR')
    % participation ratio
    PR = (sum(lambda_ssd)^2)/(sum(lambda_ssd.^2));
    ssd_n_components = round(PR);
else
    %is an integer, accounting for the number of components to keep
    ssd_n_components = opt.ssd_n_components;
end
% save ssd parameters
params.ssd_n_components = ssd_n_components;
params.lambda_ssd = lambda_ssd;
params.W_ssd = W_ssd;
%this is the low-dim raw signal in which we want to compute PCD
picks_ssd = W_ssd(:,1:ssd_n_components); %columns are filters
%apply SSD filters to the data
X_ssd = (real(X)' * W_ssd)';

X = X' * picks_ssd;
X = X';
params.X_ssd = X_ssd;
params.Xcomplex = X;
params.picks_ssd = picks_ssd;
%get the final number of filters
num = min([n_channels, size(X,1)]);

if opt.verbose > 0
    fprintf('Using the following number of dimensions = ')
    fprintf('%d ', num);
    fprintf('\n')
end

% define the options of the minimizer functions
minoptions = optimoptions('fminunc', ...
    'Algorithm', 'trust-region', ...
    'HessUpdate', 'bfgs', ...
    'GradObj','on', ...
    'Display', 'off', ...
    'UseParallel',true);
% start optimization
    for i = 1:num
        if opt.verbose > 0
            fprintf('ruuning component ')
            fprintf('%d ', i)
            fprintf('\n')
        end
        if i == 1
            % get first filter
            % get best parameters and function values of each run
            all_best = 0;
            for k = 1:opt.bestof
                w0 = (rand(num,1)*2 - 1);
                [x_best, best_f] = fminunc(@objfun_x, w0, ...
                    minoptions);
                if best_f < all_best
                    all_best = best_f;
                    all_best_x = x_best;
                end
            end
            %save results
            vlen = all_best;
            filt = all_best_x;
        else
            %get consecutive pairs of filters
            %project data into null space of previous filters
            %this is done by getting the right eigenvectors of the filter
            %maxtrix corresponding to vanishing eigenvalues
            [~, ~, Vx] = svd(filt');
            bx = Vx(1:end,i:end)'*X;
            % get best parameters and function values of each run
            all_best = 0;
            for k = 1:opt.bestof
                w0 = (rand(num - i + 1,1)*2 - 1);
                [x_best, best_f] = fminunc(@objfun_bx, ...
                    w0, minoptions);
                if best_f < all_best
                    all_best = best_f;
                    all_best_x = x_best;
                end
            end
            % save results
            vlen = cat(2, vlen, all_best);
            filt = cat(2, filt, Vx(1:end,i:end)*all_best_x); %filters are in columns
        end
        if opt.verbose > 0
            fprintf('vlen found ')
            fprintf('%f ', vlen(end))
            fprintf('\n')
        end
    end
vlen = -1 * vlen;
% for some reason, the vlen are not decreasing sorted, as they should
[vlen,sort_idx] = sort(vlen, 'descend');
W_pcd=filt(:,sort_idx);
% save PCO params
params.vlen = vlen;
params.W_pco = W_pcd;
%% post processing
% Note that the spatial filters are now defined in SSD space and have 
% to be mapped to sensor space to be applicable to the original data.
W = picks_ssd * W_pcd;
% project back to original (un-whitened) channel space
W = real(M) * W; 

%% apply PCD filters to the data
if nargout > 4
    X_pcd = W'*Xr;
end
%% calcule the unmixing and mixing matrices
unmixing = eye(n_channels, n_channels);
unmixing(1:ssd_n_components, 1:ssd_n_components) = W_pcd;
unmixing = W_ssd*unmixing; % columns are filters
% project back to original (un-whitened) channel space
unmixing = real(M) * unmixing;  
% calculate the mixing matrix
mixing = pinv(unmixing);

end
%% prepare_data
function [Cxx, X_white, M] = prepare_data(X)
% sanity checks and hilbert transform
if not(isempty(X))
    if all(isreal(X))
        X = hilbert(X');
    end
else
    error('X must be non-empty!');
end
% now size(X) = [time, channels]!
X = X - repmat(mean(X), size(X,1), 1);
Cxx = cov(real(X));

%rk = rank(Cxx);
% whiten the data
[V, D] = eig(Cxx);
[ev_sorted, sort_idx] = sort(diag(D), 'descend');
V = V(:, sort_idx);
D = diag(ev_sorted);

M = V * diag(diag(D).^-0.5); % whitening filters are in the columns!!!
X_white = X*M;
% go back to time x channnels
X_white = X_white';

% check if whitening is ok

% if  int16(diag(cov(real(X'))))~=1
%      warning('X seems not to be whitenned!');
% end

end
%% Objetive function
% The following function defines the objective function
function [vlen, grad] = PCD_obj(w, z, x)
    % Calculation of mean vector length 
    x_filtered = w'*x;
    phase = angle(x_filtered);
    % calculate the result of the objective function
    z_norm = (z - mean(z))/std(z, 1);
   
    sum1_no_square = mean(z_norm.*cos(phase));
    sum2_no_square = mean(z_norm.*sin(phase));
    % multiply with sign, i.e. -1 if function should be minimized
    % or 1 if it should be maximized
    vlen = -1*sqrt(sum1_no_square^2 + sum2_no_square^2);
    
    if nargout > 1 % calculate gradient
        % Partial derivative of phase
        phase_dwi = bsxfun(@rdivide, ...
            bsxfun(@times, -real(x), imag(x_filtered)) + ...
            bsxfun(@times, imag(x), real(x_filtered)), ...
            real(x_filtered).^2 + imag(x_filtered).^2);
%         aux1 = -real(x) .* imag(x_filtered) + imag(x) .* real(x_filtered);
%         aux2 = real(x_filtered).^2 + imag(x_filtered).^2;
%         phase_dwi = aux1./aux2;
        % Derivative of first summand
        sum1_d = 2*sum1_no_square*mean(bsxfun(@times, phase_dwi, ...
            -z_norm.*sin(phase)), 2);
        % Derivative of second summand
        sum2_d = 2*sum2_no_square*mean(bsxfun(@times, phase_dwi, ...
            z_norm.*cos(phase)), 2);
        % Derivative of sum
        sum_d = sum1_d + sum2_d;
        % Derivative including square root
        grad = sum_d/(2*vlen);
    end

end