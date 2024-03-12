function X_CAR = apply_CAR(X, method)



% X is nchan x nsamples matrix
n_ch = size(X,1);

if strcmp(method,'spatialfilter')
% apply CAR by the spatial filter definition
    J = -ones(n_ch, n_ch)./n_ch;
    W_CAR = J;
    d = logical(eye(n_ch));

    W_CAR(d) = W_CAR(d) + 1;
    X_CAR = W_CAR*X;
    %X_CAR = X_CAR';
else
    X_CAR = X - mean(X,'omitnan');
end
end