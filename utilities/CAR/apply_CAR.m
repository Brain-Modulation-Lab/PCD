function X_CAR = apply_CAR(X, method)
[~, n_ch] = size(X);
if strcmp(method,'spatialfilter')
% apply CAR by the spatial filter definition
    J = -ones(n_ch, n_ch)./n_ch;
    W_CAR = J;
    d = logical(eye(n_ch));

    W_CAR(d) = 1;
    X_CAR = W_CAR*X';
    X_CAR = X_CAR';
else
    X_CAR = X;
    idx_ch = 1:1:n_ch;
    for c=1:n_ch
        idx_car = idx_ch~=c;

        X_CAR(:,c) = X(:,c) - mean(X,2);
    end
end
end