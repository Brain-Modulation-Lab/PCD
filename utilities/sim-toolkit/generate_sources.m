function [D, cfg_all] = generate_sources(n_sim, n_sources)

% initialize variables
cfg_all = cell(1,n_sim);
D = struct();

for chan_i = 1 : n_sources 
    D.label{chan_i,1} = strcat("chan_",num2str(chan_i));
    D.fsample = 1000;
end

for sim_i = 1 : n_sim
    disp(" ##################################################### ")
    disp(" New simulation running ")
    fprintf(" Running simulation %d/%d... \n ",sim_i, n_sim);

    cfg = simulation_settings(sim_i,n_sim,n_sources);
    % create sources
    cfg = generate_sources_dyn(cfg);
    
    % create ad-hoc time vector
    cfg.sim.out.t_rel = cfg.sim.out.t - 2;

    % save results
    disp(" Saving results... ")
    %cfg.sim.sname = strcat(SURR_FOLDER,"_",num2str(sim_i),".mat");
    cfg_all{sim_i} = cfg;
    D.trial{sim_i} = cfg.sim.out.lfp';
    D.time{sim_i} = cfg.sim.out.t_rel;
    
    disp(" End simulation running ")
    disp(" ##################################################### ")

end