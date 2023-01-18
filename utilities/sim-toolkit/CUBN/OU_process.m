function n = OU_process(simulLen, Dt, tau, sigma_2, seed)
%OU_PROCESS Exact numerical solution of the Ornstein-Uhlenbeck (OU) process
%   modified from the original code of Daniel Charlebois
%
%   n = OU_process(simulLen, Dt, tau, sigma_2, seed),
%   where
%   simulLen = number of samples of the OU process in output
%   Dt = time step of the OU process in output. Units of (ms)
%   tau = relaxation time of the OU process, see eq. 2 of Cavallari et al
%   2014. Units of (ms)
%   note that Dt and tau must have the same units of time (not necessarily ms)
%   sigma_2 = it is the variance(i.e. sigma^2) of the OU process in output
%   (see eq. 2 of Cavallari et al 2014). The units of sigma^2 set the units
%   of the OU process in output. Units of [(spikes/ms)/Dt]
%   seed = seed for the random number generator (integer number)
%
%   The OU process equation can be written as (Gillespie 1996):
%   dn/dt = -n/tau + sqrt(c) * eta(t) 
%   and the parameter c (diffusion constant) corresponds to (2*sigma^2/tau) 
%   with sigma e tau defined as in the eq. 2 of Cavallari et al 2014

c=2*sigma_2/tau;

n=zeros(simulLen,1);
randn('state',seed);

n(1) = 0;

for i=2:simulLen
   r1 = randn;
   n(i) = n(i-1)*exp(-Dt/tau) + sqrt((c*tau*0.5)*(1-(exp(-Dt/tau))^2))*r1;
end
