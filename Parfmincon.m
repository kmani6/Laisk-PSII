function[] = Parfmincon(analysis_name)

ts = {};
ys = {};
Fs = {};
FvFm = {};
species = {};
FvFmErr = {};
XErr = {};

parfor i = 1:4
    
    [ts{i},ys{i},Fs{i}, FvFm{i}, species{i}, FvFmErr{i}, XErr{i}] = get_optimal_parameters(analysis_name)
    
end
end
