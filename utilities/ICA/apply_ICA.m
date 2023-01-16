function [X_denoised, proj_mat, idx_keep] = apply_ICA(W,A,X,score,idx_, n_remove, verbose)
[n_channels, ~] = size(X);%% low-rank approximation 

if isempty(n_remove)
    idx_keep = idx_;
    n_components  = [];
elseif floor(n_remove) ~= ceil(n_remove)
    %is a decimal accounting for the threshold value
    n_components =  sum(score>n_remove);
    idx_aux = setdiff(1:length(score),1:n_components);
    idx_keep= idx_(idx_aux);
elseif  strcmp(n_remove, 'PR')
    % modified version of the participation ratio
    PR = (sum(score))/(max(score));
    n_components = round(PR);
    idx_aux = setdiff(1:length(score),1:n_components);
    idx_keep= idx_(idx_aux);

elseif  strcmp(n_remove, 'diff')
    % find the elbow
    score_d = diff([1, score]);
    n_components = find(abs(score_d)< 100*min(abs(score_d)), 1);
    idx_aux = setdiff(1:length(score),1:n_components);
    idx_keep= idx_(idx_aux);

else
    % array with the components to be removed
    n_components = n_remove;
    idx_aux = setdiff(1:length(score),n_components);
    idx_keep= idx_(idx_aux);


end
idx_keep = [sort(idx_keep), length(score)+1:n_channels];
if verbose
    disp('Cleaning data\n')
    fprintf('Removing the following number of component = ')
    fprintf('%d ', n_components);
    fprintf('\n')
end
% calcule the projection matrix
proj_mat = A(idx_keep,:)'*W(:, idx_keep)';
X_denoised = proj_mat * X;
end