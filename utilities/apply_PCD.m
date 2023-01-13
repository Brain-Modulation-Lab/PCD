function [X_denoised, proj_mat, idx_keep] = apply_PCD(W,A,X,vlen,n_remove, verbose)
[n_channels, ~] = size(X);
%% low-rank approximation 
if isempty(n_remove)
    idx_keep = 1:length(vlen);
    n_components  = [];
elseif floor(n_remove) ~= ceil(n_remove)
    %is a decimal accounting for the threshold value
    n_components =  sum(vlen>n_remove);
    idx_keep = setdiff(1:length(vlen),1:n_components);
elseif  strcmp(n_remove, 'PR')
    % modified version of the participation ratio
    PR = (sum(vlen))/(max(vlen));
    n_components = round(PR);
    idx_keep = setdiff(1:length(vlen),1:n_components);
elseif  strcmp(n_remove, 'diff')
    % find the elbow
    vlen_d = diff([1, vlen]);
    n_components = find(abs(vlen_d)< 100*min(abs(vlen_d)), 1);
    idx_keep = setdiff(1:length(vlen),1:n_components);
else
    % array with the components to be removed
    n_components = n_remove;
    idx_keep = setdiff(1:length(vlen),n_components);

end

idx_keep = [idx_keep, length(vlen)+1:n_channels];
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