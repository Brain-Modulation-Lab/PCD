function cfg = generate_sources_dyn(cfg)

% initialize variables
lfp_all = [];
pows_all = [];
ext_sign_all = [];

% create sources
for ii = 1 : cfg.src.n_sources
    fprintf(" Creating source number  %d...\n ",ii);
    % create network
    parameters_CUBN; % to generate the structure net_CUBN with all the parameters of the network
    M = cfg.sim.simulation_length/net_CUBN.Dt; % length of the simulation in time steps
    % create external input
    external_noise = OU_process(M, net_CUBN.Dt, 16, cfg.src.ext.taun*net_CUBN.Dt, cfg.sim.SEED_OU(ii)); %0.16
    t = 0 : net_CUBN.Dt : (M-1)*net_CUBN.Dt;
    external_signal = my_gaussian(t,cfg.src.ext.mus(ii),cfg.src.ext.FWHM(ii),cfg.src.ext.intensity(ii))'*net_CUBN.Dt;
    %external_signal = net_CUBN.Dt*(cfg.src.ext.intensity(ii) + cfg.src.ext.A(ii)*sin(2*pi*cfg.src.ext.fl(ii)*t'/1000 + cfg.src.ext.phi(ii))); % units; (spikes/ms)/cell
    ext_sign_all = [ext_sign_all external_signal];
    INPUT2E = external_signal + external_noise;
    INPUT2E = INPUT2E.*(INPUT2E>0); % to avoid negative spikes incoming
    INPUT2I = INPUT2E;
    INPUT2I = INPUT2I.*(INPUT2I>0); % to avoid negative spikes incoming
    % launch the core of the model
    [E2EI,I2EI,eFR,iFR] = code_CUBN(net_CUBN, INPUT2E, INPUT2I, cfg.sim.SEED_connections(ii), cfg.sim.SEED_poisson(ii));
    % build LFP proxy
    lfp_temp =(I2EI-E2EI)*net_CUBN.eRm;
    lfp_all = [lfp_all lfp_temp];
    % compute LFP spectrum
    [pxx,f] = pwelch(detrend(lfp_temp,0),hanning(cfg.sim.npoints),[],2^nextpow2(cfg.sim.npoints)*5,cfg.sim.fs_sim);
    pows_all = [pows_all pxx];
end

cfg.sim.out.pows = pows_all;
cfg.sim.out.ext_sign = ext_sign_all;

[p,q] = rat(cfg.sim.fs_new / cfg.sim.fs_sim);
lfp_all = resample(lfp_all,p,q);
cfg.sim.out.lfp = detrend(lfp_all,0);

cfg.sim.out.f = f;
cfg.sim.out.n_samples = size(cfg.sim.out.lfp,1);
cfg.sim.out.t = (0:cfg.sim.out.n_samples-1)/cfg.sim.fs_new;