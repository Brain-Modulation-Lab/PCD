function coherence_norm = get_coherence(x1,x2)
% Implement the formula to calculate the ITPC [1].

% [1] Cohen, M. X. Analyzing Neural Time Series Data: Theory and Practice.
% (MIT Press, 2014). 

%author: mvissani
%last version: Aug 2022

if size(x1,2) ~= size(x2,2)
  l=min(size(x1,2),size(x2,2));
  x1=x1(1,1:l);
  x2=x2(1,1:l);
end
x1 = x1 - mean(x1);
x2 = x2 - mean(x2);
% calculate coherence values
coherence = [dot(hilbert(x2),x1),norm(x1),norm(x2)];
% normalize value of coherence
coherence_norm = coherence(1)/(coherence(2)*coherence(3));
