function sqerror = calc_sqerror1(x,... Set of parameters. This only includes the independent variables as described by the third column in Y and Constants files
                    n_trains, n_flashes, flash_duration, flash_interval, train_interval, ... Experimental parameters
                    Fluorescence_k_idcs, Fluorescence_y_inds,... Indeces to calculate Fluorescence
                    kidcs, PSIidcs, ... all indices needed in to calculate FvFm and prepare the variables
                    tablek, tabley,... information on the k and y variables
                    kconst, rate_inds, S, species, knames, species_idcs,... model specific variables
                    FvFm_exp,yidcs, ATPpar,kf1indcs, kf2indcs)
                
try
                
[FvFm_sim, grads] = calc_FvFm1(x,... Set of parameters. This only includes the independent variables as described by the third column in Y and Constants files
                    n_trains, n_flashes, flash_duration, flash_interval, train_interval, ... Experimental parameters
                    Fluorescence_k_idcs, Fluorescence_y_inds,... Indeces to calculate Fluorescence
                    kidcs, PSIidcs, ... all indices needed in to calculate FvFm and prepare the variables
                    tablek, tabley,... information on the k and y variables
                    kconst, rate_inds, S, species, knames, species_idcs, yidcs,ATPpar,kf1indcs, kf2indcs); % model specific variables
                                     
sqerror = norm(FvFm_exp-FvFm_sim)^2/numel(FvFm_exp); %+ 1e-10* norm(grads)/numel(grads);

catch err
  disp(getReport(err));
    
    sqerror = nan;
    
end
