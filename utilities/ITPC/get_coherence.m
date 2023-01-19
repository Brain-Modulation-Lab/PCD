function NC = get_coherence(x1,x2)
% Implement the formula to calculate the ITPC [1].
% Inputs:
%       x1: [nchans1 x nsamples] array of neural signals
%       x2: [nchans2 x nsamples] array of targets
% Output:
%       NC: [nchans2 x nchans1] array of coherence between neural sifnals
%       and targets

% [1] Cohen, M. X. Analyzing Neural Time Series Data: Theory and Practice.
% (MIT Press, 2014). 

%author: mvissani
%last version: Aug 2022

if size(x1,2) ~= size(x2,2)
  l=min(size(x1,2),size(x2,2));
  x1=x1(:,1:l);
  x2=x2(:,1:l);
end
nchans1 = size(x1,1);
nchans2 = size(x2,1);
nsamples = size(x1,2);

if nchans1 >= nsamples
    x1 = x1';
end

if nchans2 >= nsamples
    x1 = x2';
end

% demean signals (dim = 2 ensures consistency if nchan1 or nchans2 is 1)
x1 = x1 - mean(x1,2,'omitnan');
x2 = x2 - mean(x2,2,'omitnan');

% calculate coherence values
NC = nan(nchans2, nchans1);
for ch2 = 1 : nchans2
    for ch1 = 1 : nchans1
    C = [dot(hilbert(x2(ch2,:)),x1(ch1,:)),norm(x1(ch1,:)),norm(x2(ch2,:))];
    % normalize value of coherence
    NC(ch2,ch1) = C(1)/(C(2)*C(3));
    end
end
end
