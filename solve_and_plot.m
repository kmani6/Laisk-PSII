[ts,ys, Fs, FvFm_sim] = calc_Species_concs(xopt,... Set of parameters. This only includes the independent variables as described by the third column in Y and Constants files
                    n_trains, n_flashes, flash_duration, flash_interval, train_interval, ... Experimental parameters
                    kidcs, PSIidcs, ... all indices needed in to calculate FvFm and prepare the variables
                    tablek, tabley,... information on the k and y variables
                    kconst, rate_inds, S, species, knames, species_idcs); % model specific variables  
 
                
                
f_s_states = plot_S_states(species, ys, ts);
f_fvfm = figure; hold on;
plot(1:length(FvFm_sim), FvFm_sim, '.-');
scatter(1:length(FvFm_exp), FvFm_exp);

f_PQ = plot_pq_redox_state(species, ys, ts);
f_PC = plot_pc_redox_state(species, ys, ts);
f_NADP = plot_NAD_redox_state(species, ys, ts);
f_H = plot_H_species(species, ys, ts);
f_fd = plot_fd_redox_state(species, ys, ts);
foo = 1;