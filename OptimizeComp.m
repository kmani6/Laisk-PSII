function[] = OptimizeComp(analysis_name,datestring)
%Used after OptimizeXopt to generate 4 figures. Will analyze non-optimal
%xerr vs. fvals, generate a histogram for nonoptimal fval, generate a
%histogram for nonoptimal xerr, and take optimal values and scatter plot
%fvfmexp and fvfmopt to determine if fits are similar. Will also iterate
%through every "get_ms_paramater run" to decompose xopt vector into y
%and k values and add extra column in y and k excel sheets. 

cd results
cd(analysis_name)
cd(datestring)
files = dir;
filenames = {files.name};
filenames = filenames(3:end); 

optimalfile = filenames(end);
nonoptfiles = filenames(1:end-1);
fvalbest = inf;
FvFmErr = zeros(length(nonoptfiles),1);
XErr = zeros(length(nonoptfiles),1); 

for i = 1:(length(nonoptfiles))
    load(filenames{i}, 'xopt', 'fval', 'exitflag', 'output',...optimization parameters
        'n_trains', 'n_flashes', 'flash_duration', 'flash_interval', 'train_interval', ... Experimental parameters
        'kidcs', 'PSIidcs',... all indices needed in to calculate FvFm and prepare the variables
        'tablek', 'tabley', 'Rknames',... information on the k and y variables
        'kconst', 'rate_inds', 'S', 'species', 'knames', 'species_idcs',...
        'FvFm_exp', 'FvFm_sim_opt')
    if fval < fvalbest
        fvalbest = fval;
        optimalfile = filenames{i};
        
    end
    y_std = tabley.base_val;
    tabley.xopt = zeros(length(y_std),1);
    xoptindy = find(tabley.independent);
    
    for j = 1:length(xoptindy)
        tabley.xopt(xoptindy(j)) = xopt(j);
    end
    
    k_std = tablek.base_val;
    tablek.xopt = zeros(length(k_std),1);
    xoptindk = find(tablek.independent);
    
    for k = 1:length(xoptindk)
        tablek.xopt(xoptindk(k)) = xopt(k);
    end
    
    kr = k_std(xoptindk);
    yr = y_std(xoptindy);
    
    xknown = [yr;kr];
    XErr(i) = (norm(xopt-xknown))^2;
    FvFmErr(i) = fval;
     
end

figure 
histogram(FvFmErr)
hold on

figure 
histogram(XErr);
hold on

figure 
Color = 'green';
S = scatter(XErr,FvFmErr,[],Color);
hold on
set(S, 'SizeData', 8)

load(optimalfile, 'xopt', 'fval', 'exitflag', 'output',...optimization parameters
    'n_trains', 'n_flashes', 'flash_duration', 'flash_interval', 'train_interval', ... Experimental parameters
    'kidcs', 'PSIidcs',... all indices needed in to calculate FvFm and prepare the variables
    'tablek', 'tabley', 'Rknames',... information on the k and y variables
    'kconst', 'rate_inds', 'S', 'species', 'knames', 'species_idcs',...
    'FvFm_exp', 'FvFm_sim_opt')

y_std = tabley.base_val;
tabley.xopt = zeros(length(y_std),1);
xoptindy = find(tabley.independent);

for l = 1:length(xoptindy)
    tabley.xopt(xoptindy(l)) = xopt(l);
end

k_std = tablek.base_val;
tablek.xopt = zeros(length(k_std),1);
xoptindk = find(tablek.independent);

for s = 1:length(xoptindk)
    tablek.xopt(xoptindk(s)) = xopt(s);
end

figure 
Color = 'green';
S = scatter(1:length(FvFm_exp),FvFm_exp,[],Color);
hold on
set(S, 'SizeData', 8)

Color = 'blue';
S2 = scatter(1:length(FvFm_sim_opt),FvFm_sim_opt,[],Color);
set(S2, 'SizeData', 8)
hold on

end


